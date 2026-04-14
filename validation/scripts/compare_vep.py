#!/usr/bin/env python3
"""Compare nirvana2vcf CSQ against native VEP CSQ output.

Cross-tool comparison: nirvana2vcf derives CSQ from Nirvana JSON, VEP produces
CSQ independently.  Exact matching is NOT expected — different transcript
versions, database builds, and algorithm choices cause legitimate divergence.
The goal is characterisation, not pass/fail.

Usage:
    python3 compare_vep.py \
        --nirvana2vcf data/output/phase2/HiSeq.10000.default.vcf \
        --vep     data/phase5/vep_output.vcf \
        --output  results/phase5/vep_comparison.txt \
        [--verbose]
"""

import argparse
import re
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Subfields we compare (subset of the union of both CSQ definitions).
COMPARE_FIELDS = [
    "Consequence",
    "SYMBOL",
    "BIOTYPE",
    "CANONICAL",
    "SIFT",
    "PolyPhen",
    "HGVSc",
    "HGVSp",
    "EXON",
    "INTRON",
]


# ---------------------------------------------------------------------------
# FieldStats dataclass (matches compare_annotations.py pattern)
# ---------------------------------------------------------------------------

@dataclass
class FieldStats:
    name: str
    concordant: int = 0
    discordant: int = 0
    both_missing: int = 0
    only_in_j2v: int = 0
    only_in_vep: int = 0
    examples: list = field(default_factory=list)
    max_examples: int = 20

    def add_example(self, pos_key: str, j2v_val, vep_val):
        if len(self.examples) < self.max_examples:
            self.examples.append((pos_key, j2v_val, vep_val))

    @property
    def total(self) -> int:
        return (self.concordant + self.discordant
                + self.only_in_j2v + self.only_in_vep)

    @property
    def pct(self) -> str:
        if self.total == 0:
            return "N/A"
        return f"{100.0 * self.concordant / self.total:.1f}%"


# ---------------------------------------------------------------------------
# VCF / CSQ parsing helpers
# ---------------------------------------------------------------------------

def parse_csq_format(header_lines: List[str], csq_id: str = "CSQ") -> List[str]:
    """Extract CSQ subfield names from the ##INFO=<ID=CSQ,...> header line."""
    for line in header_lines:
        if line.startswith(f"##INFO=<ID={csq_id},"):
            m = re.search(r'Format:\s*([^"]+)', line)
            if m:
                return m.group(1).strip().split('|')
    return []


def read_vcf(path: str) -> Tuple[List[str], List[List[str]]]:
    """Read a VCF file, returning (header_lines, data_rows).

    Each data_row is a tab-split list of column values.
    """
    headers = []
    rows = []
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                headers.append(line.rstrip('\n'))
            else:
                rows.append(line.rstrip('\n').split('\t'))
    return headers, rows


def parse_info(info_str: str) -> Dict[str, str]:
    if info_str == '.':
        return {}
    result = {}
    for part in info_str.split(';'):
        if '=' in part:
            k, v = part.split('=', 1)
            result[k] = v
        else:
            result[part] = ""
    return result


def parse_csq_entries(
    info: Dict[str, str],
    field_names: List[str],
    csq_id: str = "CSQ",
) -> List[Dict[str, str]]:
    """Parse CSQ value into list of subfield dicts."""
    raw = info.get(csq_id)
    if not raw:
        return []
    entries = []
    for entry_str in raw.split(','):
        parts = entry_str.split('|')
        d = {}
        for i, name in enumerate(field_names):
            d[name] = parts[i] if i < len(parts) else ''
        entries.append(d)
    return entries


# ---------------------------------------------------------------------------
# Variant key helpers
# ---------------------------------------------------------------------------

def _strip_chr(chrom: str) -> str:
    """Normalise contig name — remove 'chr' prefix for matching."""
    if chrom.startswith('chr'):
        return chrom[3:]
    return chrom


def _variant_key(cols: List[str]) -> Tuple[str, str, str, str]:
    """(bare_chrom, POS, REF, ALT) from VCF columns."""
    chrom = _strip_chr(cols[0])
    pos = cols[1]
    ref = cols[3]
    alt = cols[4]
    return (chrom, pos, ref, alt)


# ---------------------------------------------------------------------------
# Transcript matching helpers
# ---------------------------------------------------------------------------

def _base_transcript(tid: str) -> str:
    """Strip version suffix: ENST00000456328.2 → ENST00000456328."""
    return tid.split('.')[0] if '.' in tid else tid


def _extract_prediction(value: str) -> str:
    """Extract categorical prediction from 'prediction(score)' format."""
    if not value:
        return ''
    idx = value.find('(')
    label = value[:idx] if idx > 0 else value
    return label.strip().lower().replace(' ', '_')


# Pseudogene subtypes Ensembl/VEP distinguish that Nirvana collapses to "pseudogene".
_BIOTYPE_COLLAPSE = {
    'transcribed_pseudogene': 'pseudogene',
    'transcribed_processed_pseudogene': 'pseudogene',
    'transcribed_unprocessed_pseudogene': 'pseudogene',
    'transcribed_unitary_pseudogene': 'pseudogene',
    'translated_processed_pseudogene': 'pseudogene',
    'translated_unprocessed_pseudogene': 'pseudogene',
    'unprocessed_pseudogene': 'pseudogene',
    'processed_pseudogene': 'pseudogene',
    'polymorphic_pseudogene': 'pseudogene',
    'unitary_pseudogene': 'pseudogene',
    'IG_pseudogene': 'pseudogene',
    'TR_pseudogene': 'pseudogene',
    # Ensembl retaxonomised lincRNA/antisense into `lncRNA` after Nirvana's
    # bundled Ensembl 91 cache was frozen; Nirvana still emits `misc_RNA` for
    # many of these loci.  Same kind of version drift as the pseudogene split.
    'misc_RNA': 'lncRNA',
}


def _normalize_biotype(b: str) -> str:
    return _BIOTYPE_COLLAPSE.get(b, b)


# SO terms added after Nirvana's bundled ontology was frozen — ignored when present.
_VEP_NEWER_TERMS = {'splice_polypyrimidine_tract_variant'}


def _strip_transcript_prefix(hgvs: str) -> str:
    """Strip transcript prefix from HGVS notation: NM_001.2:c.123A>G → c.123A>G."""
    if ':' in hgvs:
        return hgvs.split(':', 1)[1]
    return hgvs


_HGVSP_INNER_PAREN = re.compile(r'^p\.\((.+)\)$')


_RC_TABLE = str.maketrans('ACGTacgt', 'TGCAtgca')


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


# HGVS `n.` (non-coding) component regexes, used to detect Nirvana's
# genomic-orientation alleles vs VEP's transcript-strand HGVS spec form.
# Only `n.` is affected: `c.` HGVSc is always transcript-oriented in both tools.
_HGVS_N_SUB = re.compile(r'^n\.([\d+\-*_]+)([ACGT]+)>([ACGT]+)$')
_HGVS_N_DEL = re.compile(r'^n\.([\d+\-*_]+)del([ACGT]+)$')
_HGVS_N_DUP = re.compile(r'^n\.([\d+\-*_]+)dup([ACGT]+)$')
_HGVS_N_INS = re.compile(r'^n\.([\d+\-*_]+)ins([ACGT]+)$')
_HGVS_N_DELINS = re.compile(r'^n\.([\d+\-*_]+)delins([ACGT]+)$')

# (regex, tuple of capture-group indices that hold allele sequences)
_HGVS_N_PATTERNS = [
    (_HGVS_N_SUB, (2, 3)),
    (_HGVS_N_DELINS, (2,)),
    (_HGVS_N_DEL, (2,)),
    (_HGVS_N_DUP, (2,)),
    (_HGVS_N_INS, (2,)),
]


def _hgvsc_strand_equivalent(a: str, b: str) -> bool:
    """True if two HGVSc strings agree up to transcript-strand reverse-complement.

    Nirvana emits genomic alleles inside HGVSc for reverse-strand non-coding
    transcripts, while VEP reverse-complements to transcript strand per HGVS.
    When the coord position is identical but the alleles are exact complements,
    that alone is diagnostic — no strand/FASTA lookup required.
    """
    for regex, allele_groups in _HGVS_N_PATTERNS:
        ma = regex.match(a)
        mb = regex.match(b)
        if not (ma and mb):
            continue
        if ma.group(1) != mb.group(1):
            return False
        return all(
            _reverse_complement(ma.group(g)) == mb.group(g)
            for g in allele_groups
        )
    return False


def _normalize_hgvsp(hgvs: str) -> str:
    """Normalise HGVSp to bare `p.<change>` form.

    Nirvana emits HGVSp as a composite HGVSc with the protein change wrapped in
    doubly-parenthesised form, e.g. `NM_001.1:c.743=(p.(Asp248=))`.  VEP emits
    the plain protein-ID form with URL-encoded `=`, e.g. `NP_001.1:p.Asp248%3D`.
    Both decode to `p.Asp248=` here.
    """
    if not hgvs:
        return ''
    # URL-decode the synonymous `%3D` VEP emits for `=`.
    value = hgvs.replace('%3D', '=')
    # Extract the `p.` portion wherever it sits.
    idx = value.find('p.')
    if idx < 0:
        return ''
    protein = value[idx:]
    # Trim any trailing wrapper parens Nirvana's composite leaves behind,
    # e.g. `p.(Asp248=))` → `p.(Asp248=)`.
    while protein.endswith(')') and protein.count('(') < protein.count(')'):
        protein = protein[:-1]
    # Strip the inner wrapping parens around the change itself.
    m = _HGVSP_INNER_PAREN.match(protein)
    if m:
        protein = 'p.' + m.group(1)
    return protein


# ---------------------------------------------------------------------------
# Comparison logic
# ---------------------------------------------------------------------------

def _load_hgnc_aliases(path: str) -> Dict[str, str]:
    """Build `any_symbol -> canonical_approved_symbol` map from HGNC TSV.

    Reads the `symbol`, `prev_symbol`, and `alias_symbol` columns from
    `hgnc_complete_set.txt`.  The latter two are pipe-delimited.  The approved
    `symbol` maps to itself and each of its former/alias names also maps to it.
    In rare collisions (a name was previously used by multiple genes), the
    first approved symbol wins — this is good enough for version-drift
    normalisation.
    """
    aliases: Dict[str, str] = {}
    with open(path, 'r', encoding='utf-8') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        try:
            i_sym = header.index('symbol')
            i_prev = header.index('prev_symbol')
            i_alias = header.index('alias_symbol')
        except ValueError:
            return aliases
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) <= max(i_sym, i_prev, i_alias):
                continue
            approved = cols[i_sym].strip()
            if not approved:
                continue
            aliases.setdefault(approved, approved)
            for field_idx in (i_prev, i_alias):
                raw = cols[field_idx].strip().strip('"')
                if not raw:
                    continue
                for name in raw.split('|'):
                    name = name.strip().strip('"')
                    if name:
                        aliases.setdefault(name, approved)
    return aliases


def _canonical_symbol(symbol: str, aliases: Optional[Dict[str, str]]) -> str:
    """Resolve symbol through HGNC alias map, or return as-is if unresolved."""
    if not symbol or not aliases:
        return symbol
    return aliases.get(symbol, symbol)


def _compare_consequence(j2v_val: str, vep_val: str) -> bool:
    """Set comparison — both can have & delimited multiple SO terms."""
    if not j2v_val and not vep_val:
        return True
    j2v_set = set(j2v_val.split('&')) if j2v_val else set()
    vep_set = set(vep_val.split('&')) if vep_val else set()
    return (j2v_set - _VEP_NEWER_TERMS) == (vep_set - _VEP_NEWER_TERMS)


def _compare_field(
    field_name: str,
    j2v_val: str,
    vep_val: str,
    hgnc_aliases: Optional[Dict[str, str]] = None,
) -> bool:
    """Compare a single CSQ subfield between nirvana2vcf and VEP."""
    if field_name == "Consequence":
        return _compare_consequence(j2v_val, vep_val)

    if field_name in ("SIFT", "PolyPhen"):
        return _extract_prediction(j2v_val) == _extract_prediction(vep_val)

    if field_name == "HGVSc":
        a = _strip_transcript_prefix(j2v_val)
        b = _strip_transcript_prefix(vep_val)
        if a == b:
            return True
        return _hgvsc_strand_equivalent(a, b)

    if field_name == "HGVSp":
        return _normalize_hgvsp(j2v_val) == _normalize_hgvsp(vep_val)

    if field_name == "CANONICAL":
        # vep_val here is VEP's MANE_SELECT transcript ID (set by the caller),
        # not VEP's CANONICAL flag — see the main compare() loop.
        return (j2v_val == 'YES') == bool(vep_val)

    if field_name == "BIOTYPE":
        return _normalize_biotype(j2v_val) == _normalize_biotype(vep_val)

    if field_name == "SYMBOL":
        return (_canonical_symbol(j2v_val, hgnc_aliases)
                == _canonical_symbol(vep_val, hgnc_aliases))

    # Exact string match for EXON, INTRON
    return j2v_val == vep_val


# ---------------------------------------------------------------------------
# Main comparison engine
# ---------------------------------------------------------------------------

def compare(
    j2v_path: str,
    vep_path: str,
    hgnc_path: Optional[str] = None,
    verbose: bool = False,
) -> Tuple[Dict[str, FieldStats], dict]:
    """Run the full comparison. Returns (per_field_stats, summary_dict)."""

    # HGNC alias map (optional) — falls back to exact SYMBOL match if absent.
    hgnc_aliases: Optional[Dict[str, str]] = None
    if hgnc_path:
        try:
            hgnc_aliases = _load_hgnc_aliases(hgnc_path)
            if verbose:
                print(f"Loaded {len(hgnc_aliases):,} HGNC alias entries "
                      f"from {hgnc_path}", file=sys.stderr)
        except OSError as e:
            print(f"WARN: could not read HGNC file {hgnc_path}: {e} "
                  f"— SYMBOL compared by exact match", file=sys.stderr)

    # Read both VCFs
    j2v_headers, j2v_rows = read_vcf(j2v_path)
    vep_headers, vep_rows = read_vcf(vep_path)

    j2v_csq_fields = parse_csq_format(j2v_headers, "CSQ")
    vep_csq_fields = parse_csq_format(vep_headers, "CSQ")

    if not j2v_csq_fields:
        print("ERROR: No CSQ header found in nirvana2vcf VCF", file=sys.stderr)
        sys.exit(1)
    if not vep_csq_fields:
        print("ERROR: No CSQ header found in VEP VCF", file=sys.stderr)
        sys.exit(1)

    if verbose:
        print(f"nirvana2vcf CSQ fields ({len(j2v_csq_fields)}): "
              f"{', '.join(j2v_csq_fields)}", file=sys.stderr)
        print(f"VEP CSQ fields ({len(vep_csq_fields)}): "
              f"{', '.join(vep_csq_fields)}", file=sys.stderr)

    # Build variant → CSQ entry lookup for both
    # Key: (bare_chrom, POS, REF, ALT)
    j2v_by_var: Dict[tuple, List[Dict[str, str]]] = defaultdict(list)
    for row in j2v_rows:
        if len(row) < 8:
            continue
        key = _variant_key(row)
        info = parse_info(row[7])
        entries = parse_csq_entries(info, j2v_csq_fields)
        j2v_by_var[key].extend(entries)

    vep_by_var: Dict[tuple, List[Dict[str, str]]] = defaultdict(list)
    for row in vep_rows:
        if len(row) < 8:
            continue
        key = _variant_key(row)
        info = parse_info(row[7])
        entries = parse_csq_entries(info, vep_csq_fields)
        vep_by_var[key].extend(entries)

    # Variant-level match statistics
    j2v_keys = set(j2v_by_var.keys())
    vep_keys = set(vep_by_var.keys())
    shared_keys = j2v_keys & vep_keys

    summary = {
        'j2v_variants': len(j2v_keys),
        'vep_variants': len(vep_keys),
        'shared_variants': len(shared_keys),
        'j2v_only_variants': len(j2v_keys - vep_keys),
        'vep_only_variants': len(vep_keys - j2v_keys),
        'total_tx_matched': 0,
        'total_tx_j2v_only': 0,
        'total_tx_vep_only': 0,
    }

    # Per-field stats
    stats = {name: FieldStats(name=name) for name in COMPARE_FIELDS}

    # Compare transcript-by-transcript within shared variants
    for vkey in sorted(shared_keys):
        j2v_entries = j2v_by_var[vkey]
        vep_entries = vep_by_var[vkey]

        pos_label = f"{vkey[0]}:{vkey[1]} {vkey[2]}>{vkey[3]}"

        # Index by base transcript ID
        j2v_by_tx: Dict[str, Dict[str, str]] = {}
        for e in j2v_entries:
            feature = e.get('Feature', '')
            if feature:
                base = _base_transcript(feature)
                j2v_by_tx[base] = e

        vep_by_tx: Dict[str, Dict[str, str]] = {}
        for e in vep_entries:
            feature = e.get('Feature', '')
            if feature:
                base = _base_transcript(feature)
                vep_by_tx[base] = e

        j2v_tx_set = set(j2v_by_tx.keys())
        vep_tx_set = set(vep_by_tx.keys())
        matched_tx = j2v_tx_set & vep_tx_set

        summary['total_tx_matched'] += len(matched_tx)
        summary['total_tx_j2v_only'] += len(j2v_tx_set - vep_tx_set)
        summary['total_tx_vep_only'] += len(vep_tx_set - j2v_tx_set)

        for tx_base in matched_tx:
            j2v_e = j2v_by_tx[tx_base]
            vep_e = vep_by_tx[tx_base]
            tx_label = f"{pos_label} {tx_base}"

            for fname in COMPARE_FIELDS:
                st = stats[fname]
                j2v_val = j2v_e.get(fname, '')
                if fname == "CANONICAL":
                    # Reframe CANONICAL as a cross-tool preferred-transcript
                    # axis: Nirvana's CANONICAL flag vs VEP's MANE_SELECT.
                    # VEP's own CANONICAL uses a different selector from
                    # Nirvana's, so the direct flag-to-flag compare is noise.
                    vep_val = vep_e.get('MANE_SELECT', '')
                else:
                    vep_val = vep_e.get(fname, '')

                if not j2v_val and not vep_val:
                    st.both_missing += 1
                    continue
                if j2v_val and not vep_val:
                    st.only_in_j2v += 1
                    st.add_example(tx_label, j2v_val, "(missing)")
                    continue
                if not j2v_val and vep_val:
                    st.only_in_vep += 1
                    st.add_example(tx_label, "(missing)", vep_val)
                    continue

                if _compare_field(fname, j2v_val, vep_val, hgnc_aliases):
                    st.concordant += 1
                else:
                    st.discordant += 1
                    st.add_example(tx_label, j2v_val, vep_val)

    return stats, summary


# ---------------------------------------------------------------------------
# Report writer
# ---------------------------------------------------------------------------

def _pct(num: int, denom: int) -> str:
    if denom == 0:
        return "N/A"
    return f"{100.0 * num / denom:.1f}%"


def _trunc(s, maxlen=100):
    s = str(s)
    return s[:maxlen - 3] + "..." if len(s) > maxlen else s


def write_report(
    stats: Dict[str, FieldStats],
    summary: dict,
    j2v_path: str,
    vep_path: str,
    elapsed: float,
    out=None,
):
    if out is None:
        out = sys.stdout

    p = lambda *a, **kw: print(*a, file=out, **kw)

    p("=" * 75)
    p("CROSS-TOOL COMPARISON: nirvana2vcf vs VEP (Phase 5a)")
    p("=" * 75)
    p(f"nirvana2vcf VCF: {j2v_path}")
    p(f"VEP VCF:      {vep_path}")
    p(f"Elapsed:      {elapsed:.1f}s")
    p()

    # --- Variant match statistics ---
    p("VARIANT MATCHING")
    p("-" * 75)
    p(f"  nirvana2vcf variants:   {summary['j2v_variants']:>8,}")
    p(f"  VEP variants:        {summary['vep_variants']:>8,}")
    p(f"  Shared variants:     {summary['shared_variants']:>8,}")
    p(f"  nirvana2vcf only:       {summary['j2v_only_variants']:>8,}")
    p(f"  VEP only:            {summary['vep_only_variants']:>8,}")
    p()

    # --- Transcript match statistics ---
    total_tx = (summary['total_tx_matched']
                + summary['total_tx_j2v_only']
                + summary['total_tx_vep_only'])
    p("TRANSCRIPT MATCHING (within shared variants)")
    p("-" * 75)
    p(f"  Matched (by base ID): {summary['total_tx_matched']:>8,}")
    p(f"  nirvana2vcf only:        {summary['total_tx_j2v_only']:>8,}")
    p(f"  VEP only:             {summary['total_tx_vep_only']:>8,}")
    p(f"  Total:                {total_tx:>8,}")
    p()

    # --- Per-field concordance table ---
    p("PER-FIELD CONCORDANCE (matched transcripts)")
    p("-" * 75)
    header = (f"{'Field':<18s} {'Agree':>8s} {'Disagree':>8s} "
              f"{'J2V-only':>8s} {'VEP-only':>8s} {'Total':>8s} {'Pct':>8s}")
    p(header)
    p("-" * 75)
    for fname in COMPARE_FIELDS:
        st = stats[fname]
        if st.total == 0 and st.both_missing == 0:
            p(f"{fname:<18s} {'-':>8s} {'-':>8s} {'-':>8s} "
              f"{'-':>8s} {'-':>8s} {'N/A':>8s}")
            continue
        p(f"{fname:<18s} {st.concordant:>8,} {st.discordant:>8,} "
          f"{st.only_in_j2v:>8,} {st.only_in_vep:>8,} "
          f"{st.total:>8,} {st.pct:>8s}")
    p("-" * 75)
    p()

    # --- Notes on expected differences ---
    p("NOTES ON EXPECTED DIFFERENCES")
    p("-" * 75)
    p("  - PolyPhen: Nirvana uses HVAR model, VEP uses HumDiv → systematic")
    p("    categorical differences are expected (informational only)")
    p("  - SIFT: Score thresholds may differ between Nirvana and VEP versions")
    p("  - Transcript versions: Nirvana bundles Ensembl 91; VEP cache may be")
    p("    newer → transcript sets differ, matched by base ID (no version)")
    p("  - HGVSc: transcript prefix stripped; for `n.` (non-coding) HGVSc on")
    p("    reverse-strand transcripts, allele reverse-complement is tolerated")
    p("    (Nirvana emits genomic orientation; VEP emits transcript strand)")
    p("  - HGVSp: normalised to bare `p.<change>` with URL-decoded `=`")
    p("  - SYMBOL: resolved through HGNC prev/alias map when --hgnc supplied,")
    p("    so Ensembl-91 symbols match their current approved HGNC names")
    p("  - CANONICAL: Nirvana `CANONICAL=YES` is compared against VEP's")
    p("    `MANE_SELECT` (non-empty) — a cross-tool stable preferred-")
    p("    transcript axis; VEP's own CANONICAL uses a different selector")
    p()

    # --- Discordant examples ---
    has_examples = any(st.examples for st in stats.values()
                       if st.discordant > 0
                       or st.only_in_j2v > 0
                       or st.only_in_vep > 0)
    if has_examples:
        p("=" * 75)
        p("DISCORDANT EXAMPLES (first 10 per field)")
        p("=" * 75)
        for fname in COMPARE_FIELDS:
            st = stats[fname]
            if not st.examples:
                continue
            p(f"\n{fname} ({st.discordant} disagree, "
              f"{st.only_in_j2v} J2V-only, {st.only_in_vep} VEP-only):")
            for pos_key, j2v_val, vep_val in st.examples[:10]:
                p(f"  {pos_key}")
                p(f"    nirvana2vcf: {_trunc(j2v_val)}")
                p(f"    VEP:      {_trunc(vep_val)}")
        p()

    # --- CONCORDANCE lines for shell script parsing ---
    p("=" * 75)
    p("CONCORDANCE SUMMARY")
    p("=" * 75)
    for fname in COMPARE_FIELDS:
        st = stats[fname]
        total = st.total
        pct_val = f"{100.0 * st.concordant / total:.1f}" if total > 0 else "N/A"
        p(f"CONCORDANCE: {fname} {st.concordant} {total} {pct_val}")
    p("=" * 75)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare nirvana2vcf CSQ against native VEP CSQ output.",
    )
    parser.add_argument(
        "--nirvana2vcf", required=True,
        help="Path to nirvana2vcf VCF output",
    )
    parser.add_argument(
        "--vep", required=True,
        help="Path to VEP-annotated VCF output",
    )
    parser.add_argument(
        "--output",
        help="Write report to file (default: stdout)",
    )
    parser.add_argument(
        "--hgnc",
        help=("Path to HGNC complete_set.txt for SYMBOL alias resolution "
              "(optional — falls back to exact match if absent)"),
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Print progress to stderr",
    )
    args = parser.parse_args()

    t0 = time.time()
    stats, summary = compare(
        args.nirvana2vcf, args.vep,
        hgnc_path=args.hgnc,
        verbose=args.verbose,
    )
    elapsed = time.time() - t0

    out_fh = open(args.output, 'w') if args.output else None
    try:
        write_report(stats, summary, args.nirvana2vcf, args.vep, elapsed, out=out_fh)
    finally:
        if out_fh:
            out_fh.close()

    if args.output:
        print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
