#!/usr/bin/env python3
"""Compare nirvana2vcf CSQ/INFO against SnpEff ANN + SnpSift ClinVar annotations.

Cross-tool comparison: nirvana2vcf derives annotations from Nirvana JSON; SnpEff
and SnpSift annotate independently.  Exact matching is NOT expected — different
transcript databases, ClinVar versions, and SO term granularity cause legitimate
divergence.  The goal is characterisation, not pass/fail.

Usage:
    python3 compare_snpeff.py \
        --nirvana2vcf data/output/phase2/HiSeq.10000.default.vcf \
        --snpeff   data/phase5/snpeff_clinvar.vcf \
        --output   results/phase5/snpeff_comparison.txt \
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

# SnpEff ANN subfield positions (standard 16-field format)
ANN_ALLELE = 0
ANN_ANNOTATION = 1       # SO consequence term(s), & delimited
ANN_IMPACT = 2            # HIGH / MODERATE / LOW / MODIFIER
ANN_GENE_NAME = 3
ANN_GENE_ID = 4
ANN_FEATURE_TYPE = 5
ANN_FEATURE_ID = 6        # transcript ID
ANN_TX_BIOTYPE = 7
ANN_RANK = 8              # exon/intron rank
ANN_HGVS_C = 9
ANN_HGVS_P = 10
ANN_CDNA_POS = 11
ANN_CDS_POS = 12
ANN_PROTEIN_POS = 13
ANN_DISTANCE = 14
ANN_ERRORS = 15

# SO term normalisation: SnpEff → Nirvana/VEP canonical form
SO_NORMALISE = {
    "intergenic_region": "intergenic_variant",
}

# Fields we compare
COMPARE_FIELDS = [
    "Consequence",
    "Gene_symbol",
    "ClinVar_sig",
    "ClinVar_revstat",
]


# ---------------------------------------------------------------------------
# FieldStats dataclass
# ---------------------------------------------------------------------------

@dataclass
class FieldStats:
    name: str
    concordant: int = 0
    discordant: int = 0
    both_missing: int = 0
    only_in_j2v: int = 0
    only_in_snpeff: int = 0
    examples: list = field(default_factory=list)
    max_examples: int = 20

    def add_example(self, pos_key: str, j2v_val, snpeff_val):
        if len(self.examples) < self.max_examples:
            self.examples.append((pos_key, j2v_val, snpeff_val))

    @property
    def total(self) -> int:
        return (self.concordant + self.discordant
                + self.only_in_j2v + self.only_in_snpeff)

    @property
    def pct(self) -> str:
        if self.total == 0:
            return "N/A"
        return f"{100.0 * self.concordant / self.total:.1f}%"


# ---------------------------------------------------------------------------
# VCF parsing helpers
# ---------------------------------------------------------------------------

def parse_csq_format(header_lines: List[str]) -> List[str]:
    """Extract CSQ subfield names from nirvana2vcf header."""
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ,"):
            m = re.search(r'Format:\s*([^"]+)', line)
            if m:
                return m.group(1).strip().split('|')
    return []


def read_vcf(path: str) -> Tuple[List[str], List[List[str]]]:
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


def _strip_chr(chrom: str) -> str:
    if chrom.startswith('chr'):
        return chrom[3:]
    return chrom


def _variant_key(cols: List[str]) -> Tuple[str, str, str, str]:
    return (_strip_chr(cols[0]), cols[1], cols[3], cols[4])


# ---------------------------------------------------------------------------
# nirvana2vcf CSQ parsing
# ---------------------------------------------------------------------------

def parse_j2v_csq_entries(
    info: Dict[str, str],
    field_names: List[str],
) -> List[Dict[str, str]]:
    raw = info.get('CSQ')
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
# SnpEff ANN parsing
# ---------------------------------------------------------------------------

def parse_ann_entries(info: Dict[str, str]) -> List[Dict[str, str]]:
    """Parse SnpEff ANN field into list of dicts with normalised keys."""
    raw = info.get('ANN')
    if not raw:
        return []
    entries = []
    for entry_str in raw.split(','):
        parts = entry_str.split('|')
        if len(parts) < 8:
            continue
        entries.append({
            'Allele': parts[ANN_ALLELE].strip(),
            'Consequence': parts[ANN_ANNOTATION].strip(),
            'Impact': parts[ANN_IMPACT].strip(),
            'Gene_Name': parts[ANN_GENE_NAME].strip(),
            'Feature_ID': parts[ANN_FEATURE_ID].strip() if len(parts) > ANN_FEATURE_ID else '',
        })
    return entries


# ---------------------------------------------------------------------------
# ClinVar normalisation
# ---------------------------------------------------------------------------

def _normalise_clinvar_sig(raw: str) -> str:
    """Normalise ClinVar significance for comparison.

    Lowercase, replace _ with space, sort set of terms.
    Handles both nirvana2vcf %-encoded and SnpSift pipe-delimited formats.
    """
    if not raw:
        return ''
    # nirvana2vcf percent-encodes: %20 for space, %2C for comma, %26 for &
    s = raw.replace('%20', ' ').replace('%2C', ',').replace('%3B', ';')
    s = s.replace('%26', '&').replace('%25', '%').replace('%3D', '=')
    # SnpSift may use / or , or | as delimiter between terms
    s = s.replace('/', '&').replace('|', '&').replace(',', '&')
    # Also handle & from nirvana2vcf joined values
    terms = [t.strip().lower().replace('_', ' ') for t in s.split('&') if t.strip()]
    return '&'.join(sorted(set(terms)))


def _normalise_clinvar_revstat(raw: str) -> str:
    if not raw:
        return ''
    s = raw.replace('%20', ' ').replace('%2C', ',').replace('%3B', ';')
    s = s.replace('%26', '&').replace('%25', '%').replace('%3D', '=')
    return s.strip().lower().replace('_', ' ')


# ---------------------------------------------------------------------------
# SO term normalisation
# ---------------------------------------------------------------------------

def _normalise_consequences(raw: str) -> set:
    """Normalise & delimited SO terms to a canonical set.

    Drops `non_coding_transcript_variant` when a more specific co-term is
    present on the same transcript: SnpEff's GRCh37.75 does not append it as
    a modifier to intronic / splice-region / UTR calls on non-coding
    transcripts, while Nirvana does. When it is the *only* term we keep it
    (that is a legitimate standalone annotation).
    """
    if not raw:
        return set()
    terms = set()
    for t in raw.split('&'):
        t = t.strip()
        t = SO_NORMALISE.get(t, t)
        if t:
            terms.add(t)
    if len(terms) > 1 and 'non_coding_transcript_variant' in terms:
        terms.discard('non_coding_transcript_variant')
    return terms


def _base_transcript(tid: str) -> str:
    return tid.split('.')[0] if '.' in tid else tid


# ---------------------------------------------------------------------------
# HGNC alias resolution (duplicated from compare_vep.py to keep the script
# self-contained; see compare_vep.py:_load_hgnc_aliases for notes)
# ---------------------------------------------------------------------------

def _load_hgnc_aliases(path: str) -> Dict[str, str]:
    """Build `any_symbol -> canonical_approved_symbol` map from HGNC TSV."""
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
    if not symbol or not aliases:
        return symbol
    return aliases.get(symbol, symbol)


# ---------------------------------------------------------------------------
# Comparison engine
# ---------------------------------------------------------------------------

def compare(
    j2v_path: str,
    snpeff_path: str,
    hgnc_path: Optional[str] = None,
    verbose: bool = False,
) -> Tuple[Dict[str, FieldStats], dict]:

    hgnc_aliases: Optional[Dict[str, str]] = None
    if hgnc_path:
        try:
            hgnc_aliases = _load_hgnc_aliases(hgnc_path)
            if verbose:
                print(f"Loaded {len(hgnc_aliases):,} HGNC alias entries "
                      f"from {hgnc_path}", file=sys.stderr)
        except OSError as e:
            print(f"WARN: could not read HGNC file {hgnc_path}: {e} "
                  f"— Gene_symbol compared by exact match", file=sys.stderr)

    j2v_headers, j2v_rows = read_vcf(j2v_path)
    snp_headers, snp_rows = read_vcf(snpeff_path)

    j2v_csq_fields = parse_csq_format(j2v_headers)
    if not j2v_csq_fields:
        print("ERROR: No CSQ header found in nirvana2vcf VCF", file=sys.stderr)
        sys.exit(1)

    # Build variant lookup
    # nirvana2vcf: variant → {CSQ entries, INFO fields}
    j2v_by_var: Dict[tuple, dict] = {}
    for row in j2v_rows:
        if len(row) < 8:
            continue
        key = _variant_key(row)
        info = parse_info(row[7])
        csq = parse_j2v_csq_entries(info, j2v_csq_fields)
        j2v_by_var[key] = {'info': info, 'csq': csq}

    # SnpEff: variant → {ANN entries, INFO fields}
    snp_by_var: Dict[tuple, dict] = {}
    for row in snp_rows:
        if len(row) < 8:
            continue
        key = _variant_key(row)
        info = parse_info(row[7])
        ann = parse_ann_entries(info)
        snp_by_var[key] = {'info': info, 'ann': ann}

    j2v_keys = set(j2v_by_var.keys())
    snp_keys = set(snp_by_var.keys())
    shared_keys = j2v_keys & snp_keys

    summary = {
        'j2v_variants': len(j2v_keys),
        'snpeff_variants': len(snp_keys),
        'shared_variants': len(shared_keys),
        'j2v_only_variants': len(j2v_keys - snp_keys),
        'snpeff_only_variants': len(snp_keys - j2v_keys),
        'consequence_tx_matched': 0,
        'consequence_tx_j2v_only': 0,
        'consequence_tx_snpeff_only': 0,
        'clinvar_compared': 0,
    }

    stats = {name: FieldStats(name=name) for name in COMPARE_FIELDS}

    for vkey in sorted(shared_keys):
        j2v_data = j2v_by_var[vkey]
        snp_data = snp_by_var[vkey]
        pos_label = f"{vkey[0]}:{vkey[1]} {vkey[2]}>{vkey[3]}"

        # --- Consequence & Gene symbol (transcript-level) ---
        # Index nirvana2vcf CSQ by base transcript ID
        j2v_by_tx: Dict[str, Dict[str, str]] = {}
        for e in j2v_data['csq']:
            feat = e.get('Feature', '')
            if feat:
                j2v_by_tx[_base_transcript(feat)] = e

        # Index SnpEff ANN by base transcript ID
        snp_by_tx: Dict[str, Dict[str, str]] = {}
        for e in snp_data['ann']:
            feat = e.get('Feature_ID', '')
            if feat:
                snp_by_tx[_base_transcript(feat)] = e

        j2v_tx_set = set(j2v_by_tx.keys())
        snp_tx_set = set(snp_by_tx.keys())
        matched_tx = j2v_tx_set & snp_tx_set

        summary['consequence_tx_matched'] += len(matched_tx)
        summary['consequence_tx_j2v_only'] += len(j2v_tx_set - snp_tx_set)
        summary['consequence_tx_snpeff_only'] += len(snp_tx_set - j2v_tx_set)

        for tx_base in matched_tx:
            j2v_e = j2v_by_tx[tx_base]
            snp_e = snp_by_tx[tx_base]
            tx_label = f"{pos_label} {tx_base}"

            # Consequence
            j2v_csq_set = _normalise_consequences(j2v_e.get('Consequence', ''))
            snp_csq_set = _normalise_consequences(snp_e.get('Consequence', ''))
            st = stats['Consequence']
            if not j2v_csq_set and not snp_csq_set:
                st.both_missing += 1
            elif j2v_csq_set and not snp_csq_set:
                st.only_in_j2v += 1
                st.add_example(tx_label, '&'.join(sorted(j2v_csq_set)), "(missing)")
            elif not j2v_csq_set and snp_csq_set:
                st.only_in_snpeff += 1
                st.add_example(tx_label, "(missing)", '&'.join(sorted(snp_csq_set)))
            elif j2v_csq_set == snp_csq_set:
                st.concordant += 1
            else:
                st.discordant += 1
                st.add_example(tx_label,
                               '&'.join(sorted(j2v_csq_set)),
                               '&'.join(sorted(snp_csq_set)))

            # Gene symbol — resolve both sides through HGNC prev/alias map
            # before comparing (same version-drift handling as compare_vep.py).
            j2v_sym = j2v_e.get('SYMBOL', '')
            snp_sym = snp_e.get('Gene_Name', '')
            st = stats['Gene_symbol']
            if not j2v_sym and not snp_sym:
                st.both_missing += 1
            elif j2v_sym and not snp_sym:
                st.only_in_j2v += 1
                st.add_example(tx_label, j2v_sym, "(missing)")
            elif not j2v_sym and snp_sym:
                st.only_in_snpeff += 1
                st.add_example(tx_label, "(missing)", snp_sym)
            elif (_canonical_symbol(j2v_sym, hgnc_aliases)
                  == _canonical_symbol(snp_sym, hgnc_aliases)):
                st.concordant += 1
            else:
                st.discordant += 1
                st.add_example(tx_label, j2v_sym, snp_sym)

        # --- ClinVar (variant-level) ---
        j2v_info = j2v_data['info']
        snp_info = snp_data['info']

        j2v_sig = _normalise_clinvar_sig(j2v_info.get('CLINVAR_SIG', ''))
        snp_sig = _normalise_clinvar_sig(snp_info.get('CLINVAR_CLNSIG', ''))

        st = stats['ClinVar_sig']
        if not j2v_sig and not snp_sig:
            st.both_missing += 1
        elif j2v_sig and not snp_sig:
            st.only_in_j2v += 1
            st.add_example(pos_label, j2v_sig, "(missing)")
        elif not j2v_sig and snp_sig:
            st.only_in_snpeff += 1
            st.add_example(pos_label, "(missing)", snp_sig)
        elif j2v_sig == snp_sig:
            st.concordant += 1
        else:
            st.discordant += 1
            st.add_example(pos_label, j2v_sig, snp_sig)

        j2v_rev = _normalise_clinvar_revstat(j2v_info.get('CLINVAR_REVSTAT', ''))
        snp_rev = _normalise_clinvar_revstat(snp_info.get('CLINVAR_CLNREVSTAT', ''))

        st = stats['ClinVar_revstat']
        if not j2v_rev and not snp_rev:
            st.both_missing += 1
        elif j2v_rev and not snp_rev:
            st.only_in_j2v += 1
            st.add_example(pos_label, j2v_rev, "(missing)")
        elif not j2v_rev and snp_rev:
            st.only_in_snpeff += 1
            st.add_example(pos_label, "(missing)", snp_rev)
        elif j2v_rev == snp_rev:
            st.concordant += 1
        else:
            st.discordant += 1
            st.add_example(pos_label, j2v_rev, snp_rev)

        if j2v_sig or snp_sig:
            summary['clinvar_compared'] += 1

    return stats, summary


# ---------------------------------------------------------------------------
# Report writer
# ---------------------------------------------------------------------------

def _trunc(s, maxlen=100):
    s = str(s)
    return s[:maxlen - 3] + "..." if len(s) > maxlen else s


def write_report(
    stats: Dict[str, FieldStats],
    summary: dict,
    j2v_path: str,
    snpeff_path: str,
    elapsed: float,
    out=None,
):
    if out is None:
        out = sys.stdout

    p = lambda *a, **kw: print(*a, file=out, **kw)

    p("=" * 75)
    p("CROSS-TOOL COMPARISON: nirvana2vcf vs SnpEff/SnpSift (Phase 5b)")
    p("=" * 75)
    p(f"nirvana2vcf VCF:   {j2v_path}")
    p(f"SnpEff VCF:     {snpeff_path}")
    p(f"Elapsed:        {elapsed:.1f}s")
    p()

    # --- Variant matching ---
    p("VARIANT MATCHING")
    p("-" * 75)
    p(f"  nirvana2vcf variants:   {summary['j2v_variants']:>8,}")
    p(f"  SnpEff variants:     {summary['snpeff_variants']:>8,}")
    p(f"  Shared variants:     {summary['shared_variants']:>8,}")
    p(f"  nirvana2vcf only:       {summary['j2v_only_variants']:>8,}")
    p(f"  SnpEff only:         {summary['snpeff_only_variants']:>8,}")
    p()

    # --- Transcript matching ---
    total_tx = (summary['consequence_tx_matched']
                + summary['consequence_tx_j2v_only']
                + summary['consequence_tx_snpeff_only'])
    p("TRANSCRIPT MATCHING (within shared variants)")
    p("-" * 75)
    p(f"  Matched (by base ID): {summary['consequence_tx_matched']:>8,}")
    p(f"  nirvana2vcf only:        {summary['consequence_tx_j2v_only']:>8,}")
    p(f"  SnpEff only:          {summary['consequence_tx_snpeff_only']:>8,}")
    p(f"  Total:                {total_tx:>8,}")
    p()

    # --- Per-field concordance ---
    p("PER-FIELD CONCORDANCE")
    p("-" * 75)
    header = (f"{'Field':<20s} {'Agree':>8s} {'Disagree':>8s} "
              f"{'J2V-only':>8s} {'SnpE-only':>9s} {'Total':>8s} {'Pct':>8s}")
    p(header)
    p("-" * 75)
    for fname in COMPARE_FIELDS:
        st = stats[fname]
        if st.total == 0 and st.both_missing == 0:
            p(f"{fname:<20s} {'-':>8s} {'-':>8s} {'-':>8s} "
              f"{'-':>9s} {'-':>8s} {'N/A':>8s}")
            continue
        p(f"{fname:<20s} {st.concordant:>8,} {st.discordant:>8,} "
          f"{st.only_in_j2v:>8,} {st.only_in_snpeff:>9,} "
          f"{st.total:>8,} {st.pct:>8s}")
    p("-" * 75)
    p()

    # --- Notes ---
    p("NOTES ON EXPECTED DIFFERENCES")
    p("-" * 75)
    p("  - Consequence: SnpEff uses 'intergenic_region' where Nirvana uses")
    p("    'intergenic_variant' — normalised before comparison")
    p("  - Consequence: `non_coding_transcript_variant` is dropped when it")
    p("    co-occurs with a more specific SO term on the same transcript —")
    p("    SnpEff omits this modifier, Nirvana appends it. Standalone")
    p("    `non_coding_transcript_variant` calls are preserved.")
    p("  - ClinVar: Version drift between Nirvana's bundled ClinVar and the")
    p("    downloaded ClinVar VCF may cause classification differences")
    p("  - Gene symbol: resolved through HGNC prev_symbol/alias_symbol map")
    p("    when --hgnc supplied, else compared by exact match (expected")
    p("    ~3-5K version-drift mismatches without the map).")
    p("  - Transcript sets: SnpEff GRCh37.75 vs Nirvana Ensembl 91 — different")
    p("    transcript databases cause different annotation scopes")
    p()

    # --- Discordant examples ---
    has_examples = any(st.examples for st in stats.values()
                       if st.discordant > 0
                       or st.only_in_j2v > 0
                       or st.only_in_snpeff > 0)
    if has_examples:
        p("=" * 75)
        p("DISCORDANT EXAMPLES (first 10 per field)")
        p("=" * 75)
        for fname in COMPARE_FIELDS:
            st = stats[fname]
            if not st.examples:
                continue
            p(f"\n{fname} ({st.discordant} disagree, "
              f"{st.only_in_j2v} J2V-only, {st.only_in_snpeff} SnpEff-only):")
            for pos_key, j2v_val, snpeff_val in st.examples[:10]:
                p(f"  {pos_key}")
                p(f"    nirvana2vcf: {_trunc(j2v_val)}")
                p(f"    SnpEff:   {_trunc(snpeff_val)}")
        p()

    # --- CONCORDANCE lines ---
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
        description="Compare nirvana2vcf CSQ/INFO against SnpEff ANN + SnpSift ClinVar.",
    )
    parser.add_argument(
        "--nirvana2vcf", required=True,
        help="Path to nirvana2vcf VCF output",
    )
    parser.add_argument(
        "--snpeff", required=True,
        help="Path to SnpEff/SnpSift annotated VCF",
    )
    parser.add_argument(
        "--output",
        help="Write report to file (default: stdout)",
    )
    parser.add_argument(
        "--hgnc",
        help=("Path to HGNC complete_set.txt for Gene_symbol alias resolution "
              "(optional — falls back to exact match if absent)"),
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Print progress to stderr",
    )
    args = parser.parse_args()

    t0 = time.time()
    stats, summary = compare(
        args.nirvana2vcf, args.snpeff,
        hgnc_path=args.hgnc,
        verbose=args.verbose,
    )
    elapsed = time.time() - t0

    out_fh = open(args.output, 'w') if args.output else None
    try:
        write_report(stats, summary, args.nirvana2vcf, args.snpeff, elapsed, out=out_fh)
    finally:
        if out_fh:
            out_fh.close()

    if args.output:
        print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
