"""Maps Nirvana Position objects to VCF record dicts."""

import re
from dataclasses import replace
from typing import Any, Callable, Dict, List, Optional, Tuple

from .constants import CSQ_FIELDS
from .models import (
    NirvanaHeader,
    PopulationFrequency,
    Position,
    Sample,
    Variant,
)


_ESCAPE_RE = re.compile(r'[% ;=,]')
_ESCAPE_MAP = {'%': '%25', ' ': '%20', ';': '%3B', '=': '%3D', ',': '%2C'}


def _escape_info_value(value: str) -> str:
    """Percent-encode special characters in VCF INFO values."""
    return _ESCAPE_RE.sub(lambda m: _ESCAPE_MAP[m.group()], value)


def _fmt_float(value: float) -> str:
    """Format a float with 6 significant figures."""
    return f"{value:.6g}"


def _get_variant_for_allele(
    variants: Optional[List[Variant]], alt_allele: str
) -> Optional[Variant]:
    """Find the variant matching a specific alt allele."""
    if not variants:
        return None
    for v in variants:
        if v.alt_allele == alt_allele:
            return v
    return None


def _per_allele_values(
    position: Position, variant_map: Dict[str, Optional["Variant"]], extractor, formatter=str
) -> Optional[str]:
    """Build a comma-separated per-allele value string for Number=A INFO fields.

    variant_map: pre-built {alt_allele: Variant} for O(1) lookups.
    extractor: callable that takes a Variant and returns a value or None.
    Returns None if all values are missing.
    """
    values = []
    all_missing = True
    for alt in position.alt_alleles:
        variant = variant_map.get(alt)
        val = extractor(variant) if variant else None
        if val is not None:
            values.append(formatter(val))
            all_missing = False
        else:
            values.append(".")
    if all_missing:
        return None
    return ",".join(values)


def normalize_alleles(
    pos: int, ref: str, alts: List[str],
    chrom: Optional[str] = None,
    ref_fetcher: Optional[Callable[[str, int, int], str]] = None,
) -> Tuple[int, str, List[str]]:
    """Trim shared prefix/suffix from REF and all ALTs to minimal VCF representation.

    Phase 1: Right-trim shared suffix (all alleles must share the trailing base).
    Phase 2: Left-trim shared prefix (keeping at least 1 base), adjusting POS.
    Phase 3: When ``ref_fetcher`` is supplied, left-shift indels through repeats
        (matches ``bcftools norm`` behavior). ``ref_fetcher(chrom, start, end)``
        uses 0-based half-open coordinates.

    Symbolic alleles (<DEL>, etc.), reference-only (ALT='.'), and spanning
    deletions ('*') are left untouched.
    """
    if not alts:
        return pos, ref, alts

    # Identify which alts are "real" (participate in normalization)
    real_indices = [
        i for i, a in enumerate(alts)
        if a not in (".", "*") and not a.startswith("<")
    ]
    if not real_indices:
        return pos, ref, alts

    new_ref = ref
    new_alts = list(alts)

    # Phase 1: Right-trim while all alleles end with same base and all len >= 2
    while True:
        if len(new_ref) < 2 or any(len(new_alts[i]) < 2 for i in real_indices):
            break
        last = new_ref[-1]
        if any(new_alts[i][-1] != last for i in real_indices):
            break
        new_ref = new_ref[:-1]
        for idx in real_indices:
            new_alts[idx] = new_alts[idx][:-1]

    # Phase 2: Left-trim while all alleles start with same base and all len >= 2
    while True:
        if len(new_ref) < 2 or any(len(new_alts[i]) < 2 for i in real_indices):
            break
        first = new_ref[0]
        if any(new_alts[i][0] != first for i in real_indices):
            break
        new_ref = new_ref[1:]
        for idx in real_indices:
            new_alts[idx] = new_alts[idx][1:]
        pos += 1

    # Phase 3: reference-based left-shift (bcftools-norm parity).
    # Algorithm (Tan, Abecasis, Kang 2015): repeatedly extend by one reference
    # base on the left, then right-trim if the new trailing base is shared by
    # all alleles. Stops when the trailing bases diverge — i.e., the variant
    # is fully left-aligned. Pure SNVs and MNVs naturally stop after one step.
    # Skip non-variants (REF == every real ALT) — they'd shift indefinitely.
    if (
        ref_fetcher is not None and chrom is not None
        and any(new_alts[i] != new_ref for i in real_indices)
    ):
        while pos > 1:
            prev = ref_fetcher(chrom, pos - 2, pos - 1)
            if not prev:
                break
            ext_ref = prev + new_ref
            ext_alts = list(new_alts)
            for idx in real_indices:
                ext_alts[idx] = prev + new_alts[idx]
            last = ext_ref[-1]
            if any(ext_alts[i][-1] != last for i in real_indices):
                break
            new_ref = ext_ref[:-1]
            for idx in real_indices:
                new_alts[idx] = ext_alts[idx][:-1]
            pos -= 1

    return pos, new_ref, new_alts


def _remap_genotype(gt: str, target_allele_index: int) -> str:
    """Remap a multi-allelic genotype to biallelic for one ALT.

    In the original multi-allelic GT:
      0 = REF, 1 = ALT1, 2 = ALT2, ...

    In the decomposed biallelic GT for ALT at 1-based index *target_allele_index*:
      0 → 0  (REF stays REF)
      target → 1  (this ALT becomes the sole ALT)
      other → .  (other ALTs become missing)
      . → .  (missing stays missing)

    Preserves phasing: "/" stays "/", "|" stays "|".
    """
    sep = "|" if "|" in gt else "/"
    alleles = gt.split(sep)
    remapped = []
    for a in alleles:
        if a == ".":
            remapped.append(".")
        else:
            idx = int(a)
            if idx == 0:
                remapped.append("0")
            elif idx == target_allele_index:
                remapped.append("1")
            else:
                remapped.append(".")
    return sep.join(remapped)


def _decompose_sample(
    sample: Sample, target_allele_index: int, num_alts: int,
) -> Sample:
    """Create a decomposed sample for one ALT allele from a multi-allelic sample.

    target_allele_index is 1-based (1 for first ALT, 2 for second, etc.).
    """
    new_gt = _remap_genotype(sample.genotype, target_allele_index)

    # AD (Number=R): keep [ref_depth, target_alt_depth]
    new_ad = None
    if sample.allele_depths is not None and len(sample.allele_depths) > target_allele_index:
        new_ad = [sample.allele_depths[0], sample.allele_depths[target_allele_index]]

    # VF (Number=A): keep [target_alt_frequency]
    new_vf = None
    vf_idx = target_allele_index - 1  # 0-based index into ALT-only list
    if sample.variant_frequencies is not None and len(sample.variant_frequencies) > vf_idx:
        new_vf = [sample.variant_frequencies[vf_idx]]

    # SR (Number=R): keep [ref_count, target_alt_count]
    new_sr = None
    if sample.split_read_counts is not None and len(sample.split_read_counts) > target_allele_index:
        new_sr = [sample.split_read_counts[0], sample.split_read_counts[target_allele_index]]

    # PR (Number=R): keep [ref_count, target_alt_count]
    new_pr = None
    if sample.paired_end_read_counts is not None and len(sample.paired_end_read_counts) > target_allele_index:
        new_pr = [sample.paired_end_read_counts[0], sample.paired_end_read_counts[target_allele_index]]

    return replace(
        sample,
        genotype=new_gt,
        allele_depths=new_ad,
        variant_frequencies=new_vf,
        split_read_counts=new_sr,
        paired_end_read_counts=new_pr,
    )


def decompose_position(position: Position) -> List[Position]:
    """Decompose a multi-allelic Position into biallelic Positions.

    If the Position has 0 or 1 ALT alleles, returns [position] unchanged.
    Otherwise, returns one Position per ALT allele, each with:
      - alt_alleles: [single_alt]
      - variants: [matching_variant] or None
      - samples: decomposed with remapped GT, sliced AD/VF
      - All position-level fields (CHROM, POS, REF, QUAL, FILTER) copied
    """
    if len(position.alt_alleles) <= 1:
        return [position]

    num_alts = len(position.alt_alleles)
    variant_map = {v.alt_allele: v for v in position.variants} if position.variants else {}
    result = []

    for alt_idx_0based, alt_allele in enumerate(position.alt_alleles):
        alt_idx_1based = alt_idx_0based + 1

        matching_variant = variant_map.get(alt_allele)
        new_variants = [matching_variant] if matching_variant else None

        new_samples = None
        if position.samples:
            new_samples = [
                _decompose_sample(s, alt_idx_1based, num_alts)
                for s in position.samples
            ]

        new_pos = replace(
            position,
            alt_alleles=[alt_allele],
            variants=new_variants,
            samples=new_samples,
        )
        result.append(new_pos)

    return result


def build_csq_string(variant: Variant, alt_allele_map: Optional[Dict[str, str]] = None) -> Optional[str]:
    """Build VEP-style CSQ string for one variant (all transcripts).

    Returns pipe-delimited transcript annotations separated by commas
    for multiple transcripts. When alt_allele_map is provided, the CSQ
    Allele field uses the normalized allele string.
    """
    if not variant.transcripts:
        return None

    csq_allele = variant.alt_allele
    if alt_allele_map:
        csq_allele = alt_allele_map.get(variant.alt_allele, variant.alt_allele)

    csq_parts = []
    for t in variant.transcripts:
        fields = [
            csq_allele,                                           # Allele
            "&".join(t.consequence) if t.consequence else "",     # Consequence
            t.hgnc or "",                                         # SYMBOL
            t.gene_id or "",                                      # Gene
            "Transcript",                                         # Feature_type
            t.transcript or "",                                   # Feature
            t.bio_type or "",                                     # BIOTYPE
            t.exons or "",                                        # EXON
            t.introns or "",                                      # INTRON
            t.hgvsc or "",                                        # HGVSc
            t.hgvsp or "",                                        # HGVSp
            t.cdna_pos or "",                                     # cDNA_position
            t.cds_pos or "",                                      # CDS_position
            t.protein_pos or "",                                  # Protein_position
            t.amino_acids or "",                                  # Amino_acids
            t.codons or "",                                       # Codons
            "YES" if t.is_canonical else "",                      # CANONICAL
            _format_prediction_score(t.poly_phen_prediction, t.poly_phen_score),  # PolyPhen
            _format_prediction_score(t.sift_prediction, t.sift_score),        # SIFT
        ]
        csq_parts.append("|".join(fields))

    return ",".join(csq_parts)


def _format_prediction_score(prediction: Optional[str], score: Optional[float]) -> str:
    """Format a prediction tool result (PolyPhen, SIFT) as 'prediction(score)'."""
    if prediction is None:
        return ""
    result = prediction
    if score is not None:
        result += f"({score})"
    return result


def _popfreq_extractor(source_attr: str, freq_attr: str):
    """Return an extractor for a population frequency field on a Variant."""
    def extract(v):
        source = getattr(v, source_attr, None)
        return getattr(source, freq_attr, None) if source else None
    return extract


_PHYLOP_EXTRACTOR = lambda v: v.phylop_score
_DANN_EXTRACTOR = lambda v: v.dann_score
_GERP_EXTRACTOR = lambda v: v.gerp_score
_REVEL_EXTRACTOR = lambda v: v.revel_score
_SVTYPE_EXTRACTOR = lambda v: v.variant_type if v.is_structural_variant else None

_POPFREQ_FIELDS = [
    ("gnomAD_AF",     _popfreq_extractor("gnomad", "all_af"), _fmt_float),
    ("gnomAD_AC",     _popfreq_extractor("gnomad", "all_ac"), str),
    ("gnomAD_AN",     _popfreq_extractor("gnomad", "all_an"), str),
    ("gnomAD_AFR_AF", _popfreq_extractor("gnomad", "afr_af"), _fmt_float),
    ("gnomAD_AMR_AF", _popfreq_extractor("gnomad", "amr_af"), _fmt_float),
    ("gnomAD_EUR_AF", _popfreq_extractor("gnomad", "nfe_af"), _fmt_float),
    ("gnomAD_EAS_AF", _popfreq_extractor("gnomad", "eas_af"), _fmt_float),
    ("gnomAD_SAS_AF", _popfreq_extractor("gnomad", "sas_af"), _fmt_float),
    ("oneKG_AF",      _popfreq_extractor("one_kg", "all_af"), _fmt_float),
    ("oneKG_AFR_AF",  _popfreq_extractor("one_kg", "afr_af"), _fmt_float),
    ("oneKG_AMR_AF",  _popfreq_extractor("one_kg", "amr_af"), _fmt_float),
    ("oneKG_EUR_AF",  _popfreq_extractor("one_kg", "eur_af"), _fmt_float),
    ("oneKG_EAS_AF",  _popfreq_extractor("one_kg", "eas_af"), _fmt_float),
    ("oneKG_SAS_AF",  _popfreq_extractor("one_kg", "sas_af"), _fmt_float),
    ("TOPMed_AF",     _popfreq_extractor("topmed", "all_af"), _fmt_float),
]


def build_info_field(
    position: Position, csq_only: bool = False,
    alt_allele_map: Optional[Dict[str, str]] = None,
) -> str:
    """Build the INFO column string from position + variant annotations."""
    parts = []

    # Build variant lookup dict once for O(1) per-allele lookups throughout
    variant_map: Dict[str, Optional[Variant]] = (
        {v.alt_allele: v for v in position.variants} if position.variants else {}
    )

    if not csq_only:
        # Position-level fields
        if position.cytogenetic_band:
            parts.append(f"CytoBand={_escape_info_value(position.cytogenetic_band)}")

        # SV fields (position-level)
        if position.sv_end is not None:
            parts.append(f"SVEND={position.sv_end}")
        if position.sv_length is not None:
            parts.append(f"SVLEN={position.sv_length}")
        if position.ci_pos is not None:
            parts.append(f"CIPOS={','.join(str(x) for x in position.ci_pos)}")
        if position.ci_end is not None:
            parts.append(f"CIEND={','.join(str(x) for x in position.ci_end)}")

        # Per-allele fields from variants
        _add_per_allele(parts, position, variant_map, "phyloP", _PHYLOP_EXTRACTOR, _fmt_float)
        _add_per_allele(parts, position, variant_map, "DANN", _DANN_EXTRACTOR, _fmt_float)
        _add_per_allele(parts, position, variant_map, "GERP", _GERP_EXTRACTOR, _fmt_float)
        _add_per_allele(parts, position, variant_map, "REVEL", _REVEL_EXTRACTOR, _fmt_float)

        # SVTYPE per-allele
        _add_per_allele(parts, position, variant_map, "SVTYPE", _SVTYPE_EXTRACTOR, str)

        # Population frequencies (gnomAD, 1000 Genomes, TOPMed)
        for info_key, extractor, formatter in _POPFREQ_FIELDS:
            _add_per_allele(parts, position, variant_map, info_key, extractor, formatter)

        # ClinVar (from first variant that has it — variable number)
        _add_clinvar_info(parts, position)

        # SpliceAI (from first variant that has it)
        _add_splice_ai_info(parts, position)

    # CSQ (always included)
    csq_strings = []
    if position.variants:
        for variant in position.variants:
            csq = build_csq_string(variant, alt_allele_map=alt_allele_map)
            if csq:
                csq_strings.append(csq)
    if csq_strings:
        parts.append(f"CSQ={','.join(csq_strings)}")

    return ";".join(parts) if parts else "."


def _add_per_allele(
    parts: list, position: Position, variant_map: Dict[str, Optional[Variant]],
    key: str, extractor, formatter
) -> None:
    """Add a per-allele INFO field if any allele has a value."""
    val = _per_allele_values(position, variant_map, extractor, formatter)
    if val is not None:
        parts.append(f"{key}={val}")


def _add_clinvar_info(parts: list, position: Position) -> None:
    """Add ClinVar INFO fields, collecting from all variants."""
    if not position.variants:
        return

    sigs = []
    ids = []
    revstats = []
    for variant in position.variants:
        if not variant.clinvar:
            continue
        for cv in variant.clinvar:
            if cv.significance:
                sigs.extend(cv.significance)
            if cv.id:
                ids.append(cv.id)
            if cv.review_status:
                revstats.append(cv.review_status)

    if ids:
        parts.append(f"CLINVAR_ID={_escape_info_value('&'.join(ids))}")
    if sigs:
        parts.append(f"CLINVAR_SIG={_escape_info_value('&'.join(sigs))}")
    if revstats:
        parts.append(f"CLINVAR_REVSTAT={_escape_info_value('&'.join(revstats))}")


_SPLICE_AI_FIELDS = [
    ("SpliceAI_AG_SCORE", "acceptor_gain_score", _fmt_float),
    ("SpliceAI_AG_DIST", "acceptor_gain_distance", str),
    ("SpliceAI_AL_SCORE", "acceptor_loss_score", _fmt_float),
    ("SpliceAI_AL_DIST", "acceptor_loss_distance", str),
    ("SpliceAI_DG_SCORE", "donor_gain_score", _fmt_float),
    ("SpliceAI_DG_DIST", "donor_gain_distance", str),
    ("SpliceAI_DL_SCORE", "donor_loss_score", _fmt_float),
    ("SpliceAI_DL_DIST", "donor_loss_distance", str),
]


def _add_splice_ai_info(parts: list, position: Position) -> None:
    """Add SpliceAI INFO fields from variants."""
    if not position.variants:
        return

    for variant in position.variants:
        if not variant.splice_ai:
            continue
        for sai in variant.splice_ai:
            for info_key, attr_name, formatter in _SPLICE_AI_FIELDS:
                val = getattr(sai, attr_name)
                if val is not None:
                    parts.append(f"{info_key}={formatter(val)}")


def build_sample_columns(
    samples: Optional[List[Sample]],
) -> Tuple[str, List[str]]:
    """Build FORMAT string and per-sample value strings.

    Returns (format_string, [sample1_values, sample2_values, ...]).
    FORMAT fields are included dynamically based on which fields have data.
    GT is always first.
    """
    if not samples:
        return "", []

    # Determine which FORMAT fields are present across all samples
    format_keys = ["GT"]  # GT always first

    for key, extractor in _FORMAT_EXTRACTORS:
        for sample in samples:
            val = extractor(sample)
            if val is not None:
                format_keys.append(key)
                break

    # Build per-sample value strings
    sample_strings = []
    for sample in samples:
        values = []
        for key in format_keys:
            if key == "GT":
                values.append(sample.genotype)
            else:
                extractor = _FORMAT_EXTRACTOR_DICT[key]
                val = extractor(sample)
                if val is None:
                    values.append(".")
                else:
                    values.append(val)
        sample_strings.append(":".join(values))

    return ":".join(format_keys), sample_strings


def _scalar_fmt(attr: str, formatter=str):
    """Create extractor for a scalar sample attribute."""
    def extract(s):
        val = getattr(s, attr)
        return formatter(val) if val is not None else None
    return extract


def _list_fmt(attr: str, formatter=str):
    """Create extractor for a list sample attribute (comma-joined)."""
    def extract(s):
        val = getattr(s, attr)
        return ",".join(formatter(x) for x in val) if val else None
    return extract


def _ft_extractor(s) -> Optional[str]:
    if s.failed_filter is None:
        return None
    return "FAIL" if s.failed_filter else "PASS"


def _dn_extractor(s) -> Optional[str]:
    if s.is_de_novo is None:
        return None
    return "true" if s.is_de_novo else "false"


# Pre-built once at module load; avoids re-creating closures on every position.
_FORMAT_EXTRACTORS: List[Tuple[str, Any]] = [
    ("DP", _scalar_fmt("total_depth")),
    ("GQ", _scalar_fmt("genotype_quality")),
    ("AD", _list_fmt("allele_depths")),
    ("VF", _list_fmt("variant_frequencies", _fmt_float)),
    ("CN", _scalar_fmt("copy_number")),
    ("FT", _ft_extractor),
    ("DN", _dn_extractor),
    ("DQ", _scalar_fmt("de_novo_quality", _fmt_float)),
    ("SR", _list_fmt("split_read_counts")),
    ("PR", _list_fmt("paired_end_read_counts")),
    ("SQ", _scalar_fmt("somatic_quality", _fmt_float)),
]
_FORMAT_EXTRACTOR_DICT: Dict[str, Any] = dict(_FORMAT_EXTRACTORS)


def _get_format_extractors() -> List[Tuple[str, Any]]:
    """Return list of (FORMAT_key, extractor_function) pairs."""
    return _FORMAT_EXTRACTORS


def map_position_to_vcf_record(
    position: Position,
    header: NirvanaHeader,
    csq_only: bool = False,
    include_samples: bool = True,
    normalize: bool = False,
    ref_fetcher: Optional[Callable[[str, int, int], str]] = None,
) -> Dict[str, Any]:
    """Convert a Position into a VCF record dict.

    Returns dict with keys: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO,
    and optionally FORMAT, samples.

    When normalize=True, REF/ALT alleles are trimmed to minimal VCF
    representation (shared prefix/suffix removed, POS adjusted). When
    ``ref_fetcher`` is also provided, indels are additionally left-shifted
    through repeats (matches ``bcftools norm`` behavior).
    """
    # Apply allele normalization if requested
    vcf_pos = position.position
    vcf_ref = position.ref_allele
    vcf_alts = list(position.alt_alleles)
    alt_allele_map = None

    if normalize and vcf_alts:
        vcf_pos, vcf_ref, vcf_alts = normalize_alleles(
            vcf_pos, vcf_ref, vcf_alts,
            chrom=position.chromosome, ref_fetcher=ref_fetcher,
        )
        if vcf_alts != list(position.alt_alleles):
            alt_allele_map = dict(zip(position.alt_alleles, vcf_alts))

    # ID: collect dbsnp from all variants, deduplicate preserving order
    ids = list(dict.fromkeys(
        rsid
        for v in (position.variants or [])
        if v.dbsnp
        for rsid in v.dbsnp
    ))
    id_str = ";".join(ids) if ids else "."

    # ALT
    alt_str = ",".join(vcf_alts) if vcf_alts else "."

    # QUAL
    if position.quality is not None:
        qual_str = _fmt_float(position.quality) if position.quality != int(position.quality) else str(int(position.quality))
    else:
        qual_str = "."

    # FILTER
    if position.filters is None:
        filter_str = "."
    elif position.filters == [] or position.filters == ["."]:
        filter_str = "."
    elif position.filters == ["PASS"]:
        filter_str = "PASS"
    else:
        filter_str = ";".join(position.filters)

    # INFO
    info_str = build_info_field(
        position, csq_only=csq_only, alt_allele_map=alt_allele_map,
    )

    record = {
        "CHROM": position.chromosome,
        "POS": vcf_pos,
        "ID": id_str,
        "REF": vcf_ref,
        "ALT": alt_str,
        "QUAL": qual_str,
        "FILTER": filter_str,
        "INFO": info_str,
    }

    # FORMAT + sample columns
    if include_samples and position.samples:
        format_str, sample_strs = build_sample_columns(position.samples)
        record["FORMAT"] = format_str
        record["samples"] = sample_strs
    else:
        record["FORMAT"] = ""
        record["samples"] = []

    return record
