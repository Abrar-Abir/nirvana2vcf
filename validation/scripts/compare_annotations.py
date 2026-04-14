#!/usr/bin/env python3
"""Compare annotation values between Nirvana JSON and nirvana2vcf VCF output.

Self-consistency check: both sides derive from the same Nirvana JSON file,
so the expectation is 100% concordance. Any mismatch is a nirvana2vcf bug.

Parses Nirvana JSON independently (does NOT import nirvana2vcf modules).
Streams both files position-by-position in lockstep.

Usage:
    python3 compare_annotations.py \
        --json data/output/HiSeq.10000.json.gz \
        --vcf data/output/phase2/HiSeq.10000.default.vcf \
        --mode default \
        --output results/phase4/default_annotation_concordance.txt \
        [--verbose]
"""

import argparse
import gzip
import re
import sys
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import orjson


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# VCF INFO value escaping (same order as mapper.py)
_ESCAPE_RE = re.compile(r'[% ;=,]')
_ESCAPE_MAP = {'%': '%25', ' ': '%20', ';': '%3B', '=': '%3D', ',': '%2C'}

# CSQ subfield names (19 fields, matching nirvana2vcf/constants.py)
CSQ_FIELDS = [
    "Allele", "Consequence", "SYMBOL", "Gene", "Feature_type", "Feature",
    "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position",
    "CDS_position", "Protein_position", "Amino_acids", "Codons",
    "CANONICAL", "PolyPhen", "SIFT",
]

# Population frequency field definitions
# (vcf_key, json_source_key, json_freq_key, format_type)
POPFREQ_FIELDS = [
    ("gnomAD_AF", "gnomad", "allAf", "float"),
    ("gnomAD_AC", "gnomad", "allAc", "str"),
    ("gnomAD_AN", "gnomad", "allAn", "str"),
    ("gnomAD_AFR_AF", "gnomad", "afrAf", "float"),
    ("gnomAD_AMR_AF", "gnomad", "amrAf", "float"),
    ("gnomAD_EUR_AF", "gnomad", "nfeAf", "float"),  # EUR = NFE
    ("gnomAD_EAS_AF", "gnomad", "easAf", "float"),
    ("gnomAD_SAS_AF", "gnomad", "sasAf", "float"),
    ("oneKG_AF", "oneKg", "allAf", "float"),
    ("oneKG_AFR_AF", "oneKg", "afrAf", "float"),
    ("oneKG_AMR_AF", "oneKg", "amrAf", "float"),
    ("oneKG_EUR_AF", "oneKg", "eurAf", "float"),
    ("oneKG_EAS_AF", "oneKg", "easAf", "float"),
    ("oneKG_SAS_AF", "oneKg", "sasAf", "float"),
    ("TOPMed_AF", "topmed", "allAf", "float"),
]

# Per-allele pathogenicity score fields
# (vcf_key, json_key, format_type)
SCORE_FIELDS = [
    ("phyloP", "phylopScore", "float"),
    ("DANN", "dannScore", "float"),
    ("GERP", "gerpScore", "float"),
]

# SpliceAI field definitions
# (vcf_key, json_key, format_type)
SPLICE_AI_FIELDS = [
    ("SpliceAI_AG_SCORE", "acceptorGainScore", "float"),
    ("SpliceAI_AG_DIST", "acceptorGainDistance", "str"),
    ("SpliceAI_AL_SCORE", "acceptorLossScore", "float"),
    ("SpliceAI_AL_DIST", "acceptorLossDistance", "str"),
    ("SpliceAI_DG_SCORE", "donorGainScore", "float"),
    ("SpliceAI_DG_DIST", "donorGainDistance", "str"),
    ("SpliceAI_DL_SCORE", "donorLossScore", "float"),
    ("SpliceAI_DL_DIST", "donorLossDistance", "str"),
]

# Ordered list of all non-CSQ INFO fields to compare.
# Matches the order mapper.py emits them in build_info_field().
ALL_INFO_FIELDS = (
    # Position-level
    ["CytoBand", "SVEND", "SVLEN", "CIPOS", "CIEND"]
    # Per-allele scores
    + ["phyloP", "DANN", "GERP", "REVEL"]
    # Per-allele structural
    + ["SVTYPE"]
    # Population frequencies
    + [f[0] for f in POPFREQ_FIELDS]
    # ClinVar
    + ["CLINVAR_ID", "CLINVAR_SIG", "CLINVAR_REVSTAT"]
    # SpliceAI
    + [f[0] for f in SPLICE_AI_FIELDS]
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _escape_info_value(value: str) -> str:
    """Percent-encode special characters in VCF INFO values."""
    return _ESCAPE_RE.sub(lambda m: _ESCAPE_MAP[m.group()], value)


def _fmt_float(value) -> str:
    """Format a float with 6 significant figures (.6g)."""
    return f"{float(value):.6g}"


def _format_prediction_score(prediction, score) -> str:
    """Format PolyPhen/SIFT as 'prediction(score)'."""
    if prediction is None:
        return ""
    result = str(prediction)
    if score is not None:
        result += f"({score})"
    return result


# ---------------------------------------------------------------------------
# Allele normalization (replicated from mapper.py)
# ---------------------------------------------------------------------------

def normalize_alleles(
    pos: int, ref: str, alts: List[str],
) -> Tuple[int, str, List[str]]:
    """Trim shared prefix/suffix from REF and all ALTs to minimal VCF
    representation, adjusting POS.

    Phase 1: Right-trim shared suffix.
    Phase 2: Left-trim shared prefix (keeping >= 1 base), adjusting POS.
    Symbolic alleles (<DEL>, etc.), '.', and '*' are untouched.
    """
    if not alts:
        return pos, ref, alts

    real_indices = [
        i for i, a in enumerate(alts)
        if a not in (".", "*") and not a.startswith("<")
    ]
    if not real_indices:
        return pos, ref, alts

    new_ref = ref
    new_alts = list(alts)

    # Phase 1: Right-trim
    while True:
        if len(new_ref) < 2 or any(len(new_alts[i]) < 2 for i in real_indices):
            break
        last = new_ref[-1]
        if any(new_alts[i][-1] != last for i in real_indices):
            break
        new_ref = new_ref[:-1]
        for idx in real_indices:
            new_alts[idx] = new_alts[idx][:-1]

    # Phase 2: Left-trim
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

    return pos, new_ref, new_alts


# ---------------------------------------------------------------------------
# FieldStats dataclass
# ---------------------------------------------------------------------------

@dataclass
class FieldStats:
    """Per-field comparison statistics."""
    name: str
    concordant: int = 0
    discordant: int = 0
    both_missing: int = 0
    only_in_json: int = 0
    only_in_vcf: int = 0
    coverage: int = 0  # positions where JSON had this field
    examples: list = field(default_factory=list)
    max_examples: int = 20

    def add_example(self, pos_key: str, json_val, vcf_val):
        if len(self.examples) < self.max_examples:
            self.examples.append((pos_key, json_val, vcf_val))

    @property
    def total_compared(self) -> int:
        return (self.concordant + self.discordant
                + self.only_in_json + self.only_in_vcf)

    @property
    def is_pass(self) -> bool:
        return (self.discordant == 0
                and self.only_in_json == 0
                and self.only_in_vcf == 0)


# ---------------------------------------------------------------------------
# JSON streaming iterator
# ---------------------------------------------------------------------------

def iter_json_positions(path: str):
    """Yield position dicts from Nirvana JSON (line-based streaming).

    Nirvana format:
      Line 1:    {"header":{...},"positions":[
      Lines 2-N: one position JSON object per line (trailing comma)
      End:       ],"genes":[...]}
    """
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as f:
        # Skip header line (line 1)
        next(f)

        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            # End of positions array
            if stripped.startswith(']'):
                break
            # Strip trailing comma
            if stripped.endswith(','):
                stripped = stripped[:-1]
            if not stripped:
                continue
            try:
                yield orjson.loads(stripped)
            except Exception:
                continue


# ---------------------------------------------------------------------------
# VCF streaming iterator & parsing
# ---------------------------------------------------------------------------

def iter_vcf_lines(path: str):
    """Yield column lists for each data line in a VCF."""
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            yield line.rstrip('\n').split('\t')


def parse_info(info_str: str) -> Dict[str, str]:
    """Parse VCF INFO column into {key: value_string} dict."""
    if info_str == '.':
        return {}
    result = {}
    for part in info_str.split(';'):
        if '=' in part:
            k, v = part.split('=', 1)
            result[k] = v
        else:
            result[part] = ""  # flag key
    return result


# ---------------------------------------------------------------------------
# Expected value builders (from JSON, replicating mapper.py logic)
# ---------------------------------------------------------------------------

def _get_json_revel(variant: dict) -> Optional[float]:
    """Extract REVEL score, handling dict {'score': X} or raw float."""
    revel = variant.get('revel')
    if revel is None:
        return None
    if isinstance(revel, dict):
        return revel.get('score')
    return revel


def _build_per_allele(alt_alleles, variant_map, extractor, formatter):
    """Build comma-separated per-allele value string (Number=A).

    Returns None if all alleles are missing.
    """
    values = []
    all_missing = True
    for alt in alt_alleles:
        variant = variant_map.get(alt)
        val = extractor(variant) if variant else None
        if val is not None:
            values.append(formatter(val))
            all_missing = False
        else:
            values.append('.')
    if all_missing:
        return None
    return ','.join(values)


def build_expected_info(pos_dict: dict, csq_only: bool = False) -> dict:
    """Build expected VCF INFO key-value pairs from JSON position dict.

    Returns {info_key: expected_value_string} matching what mapper.py
    would produce (excluding CSQ, handled separately).
    """
    result = {}

    if csq_only:
        return result

    variants = pos_dict.get('variants', [])
    alt_alleles = pos_dict.get('altAlleles', [])

    # --- Position-level fields ---
    cyto = pos_dict.get('cytogeneticBand')
    if cyto:
        result['CytoBand'] = _escape_info_value(cyto)

    sv_end = pos_dict.get('svEnd')
    if sv_end is not None:
        result['SVEND'] = str(sv_end)
    sv_length = pos_dict.get('svLength')
    if sv_length is not None:
        result['SVLEN'] = str(sv_length)
    ci_pos = pos_dict.get('ciPos')
    if ci_pos is not None:
        result['CIPOS'] = ','.join(str(x) for x in ci_pos)
    ci_end = pos_dict.get('ciEnd')
    if ci_end is not None:
        result['CIEND'] = ','.join(str(x) for x in ci_end)

    # --- Build variant lookup by alt allele ---
    variant_map = {}
    for v in variants:
        alt = v.get('altAllele')
        if alt:
            variant_map[alt] = v

    # --- Per-allele pathogenicity scores ---
    for vcf_key, json_key, fmt_type in SCORE_FIELDS:
        formatter = _fmt_float if fmt_type == "float" else str
        val = _build_per_allele(
            alt_alleles, variant_map,
            lambda v, k=json_key: v.get(k),
            formatter,
        )
        if val is not None:
            result[vcf_key] = val

    # REVEL (special handling for dict/float)
    val = _build_per_allele(
        alt_alleles, variant_map,
        lambda v: _get_json_revel(v),
        _fmt_float,
    )
    if val is not None:
        result['REVEL'] = val

    # SVTYPE (per-allele, only for structural variants)
    val = _build_per_allele(
        alt_alleles, variant_map,
        lambda v: v.get('variantType') if v.get('isStructuralVariant') else None,
        str,
    )
    if val is not None:
        result['SVTYPE'] = val

    # --- Per-allele population frequencies ---
    for vcf_key, json_source, json_freq_key, fmt_type in POPFREQ_FIELDS:
        formatter = _fmt_float if fmt_type == "float" else str

        def _extractor(v, src=json_source, key=json_freq_key):
            source = v.get(src)
            return source.get(key) if source else None

        val = _build_per_allele(alt_alleles, variant_map, _extractor, formatter)
        if val is not None:
            result[vcf_key] = val

    # --- ClinVar (collected from all variants) ---
    sigs = []
    ids = []
    revstats = []
    for variant in variants:
        clinvar_list = variant.get('clinvar')
        if not clinvar_list:
            continue
        for cv in clinvar_list:
            if cv.get('significance'):
                sigs.extend(cv['significance'])
            if cv.get('id'):
                ids.append(cv['id'])
            if cv.get('reviewStatus'):
                revstats.append(cv['reviewStatus'])

    if ids:
        result['CLINVAR_ID'] = _escape_info_value('&'.join(ids))
    if sigs:
        result['CLINVAR_SIG'] = _escape_info_value('&'.join(sigs))
    if revstats:
        result['CLINVAR_REVSTAT'] = _escape_info_value('&'.join(revstats))

    # --- SpliceAI (from all variants) ---
    for variant in variants:
        splice_ai_list = variant.get('spliceAI')
        if not splice_ai_list:
            continue
        for sai in splice_ai_list:
            for vcf_key, json_attr, fmt_type in SPLICE_AI_FIELDS:
                formatter = _fmt_float if fmt_type == "float" else str
                val = sai.get(json_attr)
                if val is not None:
                    result[vcf_key] = formatter(val)

    return result


def build_expected_id(pos_dict: dict) -> Optional[str]:
    """Build expected VCF ID column from JSON position dict.

    Returns None if no dbSNP entries.
    """
    variants = pos_dict.get('variants', [])
    # Deduplicate preserving order (same as mapper.py's dict.fromkeys)
    seen = {}
    for v in variants:
        dbsnp_list = v.get('dbsnp')
        if dbsnp_list:
            for rsid in dbsnp_list:
                if rsid not in seen:
                    seen[rsid] = True
    return ';'.join(seen.keys()) if seen else None


def build_expected_csq(
    pos_dict: dict, alt_allele_map: Optional[Dict[str, str]] = None,
) -> Optional[str]:
    """Build expected CSQ value from JSON position dict.

    Replicates mapper.py build_csq_string() + build_info_field() CSQ logic.
    """
    variants = pos_dict.get('variants', [])
    if not variants:
        return None

    csq_parts = []
    for variant in variants:
        transcripts = variant.get('transcripts')
        if not transcripts:
            continue

        csq_allele = variant.get('altAllele', '')
        if alt_allele_map:
            csq_allele = alt_allele_map.get(csq_allele, csq_allele)

        for t in transcripts:
            consequence = t.get('consequence', [])
            fields = [
                csq_allele,
                '&'.join(consequence) if consequence else '',
                t.get('hgnc') or '',
                t.get('geneId') or '',
                'Transcript',
                t.get('transcript') or '',
                t.get('bioType') or '',
                t.get('exons') or '',
                t.get('introns') or '',
                t.get('hgvsc') or '',
                t.get('hgvsp') or '',
                t.get('cdnaPos') or '',
                t.get('cdsPos') or '',
                t.get('proteinPos') or '',
                t.get('aminoAcids') or '',
                t.get('codons') or '',
                'YES' if t.get('isCanonical') else '',
                _format_prediction_score(
                    t.get('polyPhenPrediction'), t.get('polyPhenScore')),
                _format_prediction_score(
                    t.get('siftPrediction'), t.get('siftScore')),
            ]
            csq_parts.append('|'.join(fields))

    return ','.join(csq_parts) if csq_parts else None


# ---------------------------------------------------------------------------
# Comparison engine
# ---------------------------------------------------------------------------

def _compare_field(stats: FieldStats, pos_key: str, expected, actual):
    """Compare one field: expected (from JSON) vs actual (from VCF).

    Both are strings or None. Exact string match.
    """
    if expected is None and actual is None:
        stats.both_missing += 1
        return
    if expected is not None and actual is None:
        stats.only_in_json += 1
        stats.add_example(pos_key, expected, "(missing)")
        return
    if expected is None and actual is not None:
        stats.only_in_vcf += 1
        stats.add_example(pos_key, "(missing)", actual)
        return
    # Both present
    if expected == actual:
        stats.concordant += 1
    else:
        stats.discordant += 1
        stats.add_example(pos_key, expected, actual)


def compare_position(
    json_pos: dict,
    vcf_cols: list,
    mode: str,
    do_normalize: bool,
    stats: Dict[str, FieldStats],
):
    """Compare all annotation fields for one position.

    json_pos: raw JSON position dict
    vcf_cols: VCF tab-split columns
    """
    chrom = json_pos.get('chromosome', '')
    pos = json_pos.get('position', 0)
    ref = json_pos.get('refAllele', '')
    alts = json_pos.get('altAlleles', [])

    pos_key = f"{chrom}:{pos} {ref}>{','.join(alts)}"

    csq_only = (mode == 'csq_only')

    # Determine alt allele map for normalization
    alt_allele_map = None
    if do_normalize and alts:
        _norm_pos, _norm_ref, norm_alts = normalize_alleles(
            pos, ref, list(alts))
        if norm_alts != list(alts):
            alt_allele_map = dict(zip(alts, norm_alts))

    # --- Sanity check: chromosome alignment ---
    vcf_chrom = vcf_cols[0] if vcf_cols else ''
    if vcf_chrom != chrom:
        print(f"WARNING: Chromosome mismatch: JSON={chrom} VCF={vcf_chrom} "
              f"at {pos_key}", file=sys.stderr)

    # === INFO fields ===
    vcf_info = parse_info(vcf_cols[7]) if len(vcf_cols) > 7 else {}

    if not csq_only:
        expected_info = build_expected_info(json_pos, csq_only=False)

        for field_name in ALL_INFO_FIELDS:
            if field_name not in stats:
                continue
            expected_val = expected_info.get(field_name)
            actual_val = vcf_info.get(field_name)
            _compare_field(stats[field_name], pos_key, expected_val, actual_val)
            if expected_val is not None:
                stats[field_name].coverage += 1

    # === ID column (dbSNP) ===
    expected_id = build_expected_id(json_pos)
    actual_id = vcf_cols[2] if len(vcf_cols) > 2 else '.'
    actual_id_norm = None if actual_id == '.' else actual_id
    _compare_field(stats['ID'], pos_key, expected_id, actual_id_norm)
    if expected_id is not None:
        stats['ID'].coverage += 1

    # === CSQ field ===
    expected_csq = build_expected_csq(json_pos, alt_allele_map=alt_allele_map)
    actual_csq = vcf_info.get('CSQ')
    _compare_field(stats['CSQ'], pos_key, expected_csq, actual_csq)
    if expected_csq is not None:
        stats['CSQ'].coverage += 1


# ---------------------------------------------------------------------------
# CSQ drill-down (per-subfield comparison for discordant CSQ)
# ---------------------------------------------------------------------------

def drill_down_csq(expected_csq: str, actual_csq: str) -> List[str]:
    """Compare CSQ strings transcript-by-transcript, subfield-by-subfield.

    Returns list of human-readable difference descriptions.
    """
    diffs = []

    exp_transcripts = expected_csq.split(',')
    act_transcripts = actual_csq.split(',')

    if len(exp_transcripts) != len(act_transcripts):
        diffs.append(
            f"Transcript count: expected {len(exp_transcripts)}, "
            f"actual {len(act_transcripts)}")
        return diffs

    # Build lookup by (Allele, Feature) for matching
    def _make_key(transcript_str):
        fields = transcript_str.split('|')
        allele = fields[0] if len(fields) > 0 else ''
        feature = fields[5] if len(fields) > 5 else ''
        return (allele, feature)

    exp_by_key = {}
    for t in exp_transcripts:
        key = _make_key(t)
        exp_by_key.setdefault(key, []).append(t)

    act_by_key = {}
    for t in act_transcripts:
        key = _make_key(t)
        act_by_key.setdefault(key, []).append(t)

    all_keys = set(exp_by_key.keys()) | set(act_by_key.keys())
    for key in sorted(all_keys):
        exp_list = exp_by_key.get(key, [])
        act_list = act_by_key.get(key, [])

        if not exp_list:
            diffs.append(f"  Transcript {key}: only in VCF")
            continue
        if not act_list:
            diffs.append(f"  Transcript {key}: only in JSON")
            continue

        for exp_t, act_t in zip(exp_list, act_list):
            if exp_t == act_t:
                continue
            exp_fields = exp_t.split('|')
            act_fields = act_t.split('|')
            for idx, field_name in enumerate(CSQ_FIELDS):
                exp_f = exp_fields[idx] if idx < len(exp_fields) else ''
                act_f = act_fields[idx] if idx < len(act_fields) else ''
                if exp_f != act_f:
                    diffs.append(
                        f"  {key[1]} {field_name}: "
                        f"expected={exp_f!r} actual={act_f!r}")

    return diffs


# ---------------------------------------------------------------------------
# Report writer
# ---------------------------------------------------------------------------

def _pct(num: int, denom: int) -> str:
    if denom == 0:
        return "N/A"
    return f"{100.0 * num / denom:.2f}%"


def _trunc(s, maxlen=80):
    s = str(s)
    return s[:maxlen - 3] + "..." if len(s) > maxlen else s


def write_report(
    stats: Dict[str, FieldStats],
    field_names: List[str],
    positions_compared: int,
    json_path: str,
    vcf_path: str,
    mode: str,
    elapsed: float,
    out=None,
):
    """Write the comparison report."""
    if out is None:
        out = sys.stdout

    p = lambda *a, **kw: print(*a, file=out, **kw)

    p("=" * 75)
    p("ANNOTATION CONCORDANCE REPORT (Phase 4)")
    p("=" * 75)
    p(f"JSON file:  {json_path}")
    p(f"VCF file:   {vcf_path}")
    p(f"Mode:       {mode}")
    p(f"Positions:  {positions_compared:,}")
    p(f"Elapsed:    {elapsed:.1f}s")
    p()

    # --- Per-field results table ---
    p("FIELD RESULTS")
    p("-" * 75)
    header = (f"{'Field':<25s} {'Coverage':>8s} {'Concord':>8s} "
              f"{'Discord':>8s} {'JSON':>6s} {'VCF':>6s} {'Result':>8s}")
    p(header)
    p("-" * 75)

    all_pass = True
    total_concordant = 0
    total_discordant = 0
    fields_with_data = 0

    for name in field_names:
        st = stats[name]
        if st.total_compared == 0:
            p(f"{name:<25s} {'0':>8s} {'-':>8s} {'-':>8s} "
              f"{'-':>6s} {'-':>6s} {'N/A':>8s}")
            continue

        result = "PASS" if st.is_pass else "FAIL"
        if not st.is_pass:
            all_pass = False

        total_concordant += st.concordant
        total_discordant += st.discordant
        fields_with_data += 1

        p(f"{name:<25s} {st.coverage:>8,} {st.concordant:>8,} "
          f"{st.discordant:>8,} {st.only_in_json:>6,} "
          f"{st.only_in_vcf:>6,} {result:>8s}")

    p("-" * 75)
    p()

    # --- Coverage summary ---
    p("COVERAGE SUMMARY")
    p("-" * 75)
    # Group by category
    categories = [
        ("Population frequencies",
         [f[0] for f in POPFREQ_FIELDS]),
        ("Pathogenicity scores",
         ["phyloP", "DANN", "GERP", "REVEL"]),
        ("ClinVar",
         ["CLINVAR_ID", "CLINVAR_SIG", "CLINVAR_REVSTAT"]),
        ("SpliceAI",
         [f[0] for f in SPLICE_AI_FIELDS]),
        ("Structural variants",
         ["SVTYPE", "SVEND", "SVLEN", "CIPOS", "CIEND"]),
        ("Position-level",
         ["CytoBand"]),
        ("Other",
         ["ID", "CSQ"]),
    ]
    for cat_name, cat_fields in categories:
        available = [f for f in cat_fields if f in stats]
        if not available:
            continue
        max_cov = max(stats[f].coverage for f in available)
        if max_cov == 0:
            p(f"  {cat_name:<30s} no data")
        else:
            items = []
            for f in available:
                cov = stats[f].coverage
                if cov > 0:
                    items.append(f"{f}={cov:,}")
            p(f"  {cat_name:<30s} {', '.join(items)}")
    p()

    # --- Discordant examples ---
    has_examples = any(stats[n].examples for n in field_names
                       if stats[n].discordant > 0
                       or stats[n].only_in_json > 0
                       or stats[n].only_in_vcf > 0)
    if has_examples:
        p("=" * 75)
        p("DISCORDANT EXAMPLES")
        p("=" * 75)
        for name in field_names:
            st = stats[name]
            if not st.examples:
                continue
            if st.is_pass:
                continue
            p(f"\n{name} ({st.discordant} discordant, "
              f"{st.only_in_json} JSON-only, {st.only_in_vcf} VCF-only):")
            for pos_key, json_val, vcf_val in st.examples:
                p(f"  {pos_key}")
                p(f"    expected: {_trunc(json_val, 100)}")
                p(f"    actual:   {_trunc(vcf_val, 100)}")
                # CSQ drill-down for discordant CSQ entries
                if (name == 'CSQ'
                        and json_val != "(missing)"
                        and vcf_val != "(missing)"):
                    drilldown = drill_down_csq(str(json_val), str(vcf_val))
                    for dd in drilldown[:10]:
                        p(f"    {dd}")
        p()

    # --- Overall verdict ---
    p("=" * 75)
    if all_pass:
        p(f"RESULT: PASS  (100% concordance, {total_concordant:,} "
          f"field values verified across {fields_with_data} fields)")
    else:
        fail_fields = [n for n in field_names
                       if stats[n].total_compared > 0 and not stats[n].is_pass]
        p(f"RESULT: FAIL  ({total_discordant:,} discordant values "
          f"in {len(fail_fields)} field(s): {', '.join(fail_fields)})")
    p("=" * 75)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare annotation values between Nirvana JSON and "
                    "nirvana2vcf VCF output (self-consistency check).",
    )
    parser.add_argument(
        "--json", required=True,
        help="Path to Nirvana JSON file (.json or .json.gz)",
    )
    parser.add_argument(
        "--vcf", required=True,
        help="Path to nirvana2vcf VCF output",
    )
    parser.add_argument(
        "--mode", required=True,
        choices=["default", "raw", "no_samples", "csq_only"],
        help="Conversion mode (determines normalization and field presence)",
    )
    parser.add_argument(
        "--output",
        help="Write report to file (default: stdout)",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Print progress to stderr",
    )
    args = parser.parse_args()

    do_normalize = args.mode in ('default', 'no_samples', 'csq_only')
    csq_only = args.mode == 'csq_only'

    # Initialize field stats
    if csq_only:
        field_names = ['ID', 'CSQ']
    else:
        field_names = list(ALL_INFO_FIELDS) + ['ID', 'CSQ']

    all_stats = {name: FieldStats(name=name) for name in field_names}

    # Stream both files in lockstep
    t0 = time.time()
    json_iter = iter_json_positions(args.json)
    vcf_iter = iter_vcf_lines(args.vcf)

    positions_compared = 0

    for json_pos in json_iter:
        try:
            vcf_cols = next(vcf_iter)
        except StopIteration:
            print(f"WARNING: VCF ended before JSON at position "
                  f"{positions_compared + 1}", file=sys.stderr)
            break

        compare_position(
            json_pos, vcf_cols, args.mode, do_normalize, all_stats)
        positions_compared += 1

        if args.verbose and positions_compared % 2000 == 0:
            elapsed = time.time() - t0
            print(f"  Compared {positions_compared:,} positions "
                  f"({elapsed:.1f}s)", file=sys.stderr)

    # Check for extra VCF lines
    remaining_vcf = 0
    for _ in vcf_iter:
        remaining_vcf += 1
    if remaining_vcf > 0:
        print(f"WARNING: VCF has {remaining_vcf} extra lines after "
              f"JSON ended", file=sys.stderr)

    elapsed = time.time() - t0
    if args.verbose:
        print(f"Done in {elapsed:.1f}s ({positions_compared:,} positions)",
              file=sys.stderr)

    # Write report
    out_fh = open(args.output, 'w') if args.output else None
    try:
        write_report(
            all_stats, field_names, positions_compared,
            args.json, args.vcf, args.mode, elapsed, out=out_fh)
    finally:
        if out_fh:
            out_fh.close()

    # Print output path if redirected
    if args.output:
        print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
