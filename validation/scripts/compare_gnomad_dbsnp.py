#!/usr/bin/env python3
"""Compare nirvana2vcf gnomAD AF and dbSNP rsIDs against bcftools-direct annotation.

Cross-tool comparison: nirvana2vcf derives gnomAD/dbSNP from Nirvana's bundled data;
bcftools annotates directly from source VCFs.  Nirvana bundles gnomAD v2.1 while
the direct download is v2.1.1 — expect >99% match with minor patch differences.

Usage:
    python3 compare_gnomad_dbsnp.py \
        --nirvana2vcf data/output/phase2/HiSeq.10000.default.vcf \
        --annotated data/phase5/j2v_gnomad_annotated.vcf \
        --output results/phase5/gnomad_dbsnp_comparison.txt \
        [--verbose]
"""

import argparse
import math
import sys
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


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
    only_in_direct: int = 0
    examples: list = field(default_factory=list)
    max_examples: int = 20

    def add_example(self, pos_key: str, j2v_val, direct_val):
        if len(self.examples) < self.max_examples:
            self.examples.append((pos_key, j2v_val, direct_val))

    @property
    def total(self) -> int:
        return (self.concordant + self.discordant
                + self.only_in_j2v + self.only_in_direct)

    @property
    def pct(self) -> str:
        if self.total == 0:
            return "N/A"
        return f"{100.0 * self.concordant / self.total:.1f}%"


# ---------------------------------------------------------------------------
# VCF parsing helpers
# ---------------------------------------------------------------------------

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


def read_vcf_records(path: str) -> Dict[tuple, dict]:
    """Read VCF, returning dict keyed by (bare_chrom, POS, REF, ALT) → {info, id}."""
    records = {}
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8:
                continue
            key = (_strip_chr(cols[0]), cols[1], cols[3], cols[4])
            info = parse_info(cols[7])
            records[key] = {
                'info': info,
                'id': cols[2] if cols[2] != '.' else '',
            }
    return records


# ---------------------------------------------------------------------------
# Numeric comparison helpers
# ---------------------------------------------------------------------------

def _parse_af(val: str) -> Optional[float]:
    """Parse an AF value, returning None for missing/unparseable."""
    if not val or val == '.':
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def pearson_r(xs: List[float], ys: List[float]) -> Optional[float]:
    """Compute Pearson correlation coefficient. Returns None if < 2 pairs."""
    n = len(xs)
    if n < 2:
        return None
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if den_x == 0 or den_y == 0:
        return None
    return num / (den_x * den_y)


# ---------------------------------------------------------------------------
# Comparison engine
# ---------------------------------------------------------------------------

AF_TOLERANCE = 0.001


def compare(
    j2v_path: str,
    annotated_path: str,
    verbose: bool = False,
) -> Tuple[Dict[str, FieldStats], dict]:
    """Compare gnomAD AF and dbSNP rsIDs.

    The annotated VCF is the nirvana2vcf output re-annotated by bcftools with
    GNOMAD_DIRECT_AF from the gnomAD source VCF.
    """

    j2v_records = read_vcf_records(j2v_path)
    ann_records = read_vcf_records(annotated_path)

    j2v_keys = set(j2v_records.keys())
    ann_keys = set(ann_records.keys())
    shared_keys = j2v_keys & ann_keys

    stats = {
        'gnomAD_AF': FieldStats(name='gnomAD_AF'),
        'dbSNP_rsID': FieldStats(name='dbSNP_rsID'),
    }

    # For Pearson correlation on AF
    af_pairs_x: List[float] = []
    af_pairs_y: List[float] = []
    max_abs_diff: float = 0.0
    within_tol: int = 0
    outside_tol: int = 0

    for vkey in sorted(shared_keys):
        j2v_info = j2v_records[vkey]['info']
        ann_info = ann_records[vkey]['info']
        pos_label = f"{vkey[0]}:{vkey[1]} {vkey[2]}>{vkey[3]}"

        # --- gnomAD AF ---
        j2v_af_str = j2v_info.get('gnomAD_AF', '')
        direct_af_str = ann_info.get('GNOMAD_DIRECT_AF', '')

        j2v_af = _parse_af(j2v_af_str)
        direct_af = _parse_af(direct_af_str)

        st = stats['gnomAD_AF']
        if j2v_af is None and direct_af is None:
            st.both_missing += 1
        elif j2v_af is not None and direct_af is None:
            st.only_in_j2v += 1
            st.add_example(pos_label, j2v_af_str, "(missing)")
        elif j2v_af is None and direct_af is not None:
            st.only_in_direct += 1
            st.add_example(pos_label, "(missing)", direct_af_str)
        else:
            diff = abs(j2v_af - direct_af)
            if diff > max_abs_diff:
                max_abs_diff = diff
            af_pairs_x.append(j2v_af)
            af_pairs_y.append(direct_af)
            if diff <= AF_TOLERANCE:
                st.concordant += 1
                within_tol += 1
            else:
                st.discordant += 1
                outside_tol += 1
                st.add_example(pos_label, j2v_af_str,
                               f"{direct_af_str} (diff={diff:.6f})")

        # --- dbSNP rsID ---
        j2v_id = j2v_records[vkey]['id']
        ann_id = ann_records[vkey]['id']

        st = stats['dbSNP_rsID']
        if not j2v_id and not ann_id:
            st.both_missing += 1
        elif j2v_id and not ann_id:
            st.only_in_j2v += 1
            st.add_example(pos_label, j2v_id, "(missing)")
        elif not j2v_id and ann_id:
            st.only_in_direct += 1
            st.add_example(pos_label, "(missing)", ann_id)
        elif j2v_id == ann_id:
            st.concordant += 1
        else:
            # Check set overlap (multiple rsIDs separated by ;)
            j2v_set = set(j2v_id.split(';'))
            ann_set = set(ann_id.split(';'))
            if j2v_set == ann_set:
                st.concordant += 1
            else:
                st.discordant += 1
                st.add_example(pos_label, j2v_id, ann_id)

    r_value = pearson_r(af_pairs_x, af_pairs_y)

    summary = {
        'j2v_variants': len(j2v_keys),
        'annotated_variants': len(ann_keys),
        'shared_variants': len(shared_keys),
        'j2v_only_variants': len(j2v_keys - ann_keys),
        'annotated_only_variants': len(ann_keys - j2v_keys),
        'af_pairs': len(af_pairs_x),
        'af_within_tol': within_tol,
        'af_outside_tol': outside_tol,
        'af_max_diff': max_abs_diff,
        'af_pearson_r': r_value,
        'af_tolerance': AF_TOLERANCE,
    }

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
    annotated_path: str,
    elapsed: float,
    out=None,
):
    if out is None:
        out = sys.stdout

    p = lambda *a, **kw: print(*a, file=out, **kw)

    p("=" * 75)
    p("CROSS-TOOL COMPARISON: nirvana2vcf vs bcftools-direct (Phase 5c)")
    p("=" * 75)
    p(f"nirvana2vcf VCF:       {j2v_path}")
    p(f"bcftools-annotated: {annotated_path}")
    p(f"Elapsed:            {elapsed:.1f}s")
    p()

    # --- Variant matching ---
    p("VARIANT MATCHING")
    p("-" * 75)
    p(f"  nirvana2vcf variants:   {summary['j2v_variants']:>8,}")
    p(f"  Annotated variants:  {summary['annotated_variants']:>8,}")
    p(f"  Shared variants:     {summary['shared_variants']:>8,}")
    p(f"  nirvana2vcf only:       {summary['j2v_only_variants']:>8,}")
    p(f"  Annotated only:      {summary['annotated_only_variants']:>8,}")
    p()

    # --- gnomAD AF comparison ---
    st_af = stats['gnomAD_AF']
    p("gnomAD AF COMPARISON")
    p("-" * 75)
    p(f"  Tolerance:             {summary['af_tolerance']}")
    p(f"  Both present:          {summary['af_pairs']:>8,}")
    p(f"  Within tolerance:      {summary['af_within_tol']:>8,}")
    p(f"  Outside tolerance:     {summary['af_outside_tol']:>8,}")
    p(f"  Max absolute diff:     {summary['af_max_diff']:.6f}")
    r = summary['af_pearson_r']
    r_str = f"{r:.6f}" if r is not None else "N/A"
    p(f"  Pearson correlation:   {r_str}")
    p(f"  nirvana2vcf only:         {st_af.only_in_j2v:>8,}")
    p(f"  Direct only:           {st_af.only_in_direct:>8,}")
    p(f"  Both missing:          {st_af.both_missing:>8,}")
    p()

    # --- dbSNP rsID comparison ---
    st_id = stats['dbSNP_rsID']
    p("dbSNP rsID COMPARISON")
    p("-" * 75)
    p(f"  Agree:         {st_id.concordant:>8,}")
    p(f"  Disagree:      {st_id.discordant:>8,}")
    p(f"  nirvana2vcf only: {st_id.only_in_j2v:>8,}")
    p(f"  Direct only:   {st_id.only_in_direct:>8,}")
    p(f"  Both missing:  {st_id.both_missing:>8,}")
    p(f"  Total compared:{st_id.total:>8,}")
    p(f"  Concordance:   {st_id.pct:>8s}")
    p()

    # --- Notes ---
    p("NOTES ON EXPECTED DIFFERENCES")
    p("-" * 75)
    p("  - Nirvana bundles gnomAD v2.1; direct download is v2.1.1 (minor patch)")
    p("  - Multi-allelic sites: bcftools may report different AF for multi-allelic")
    p("    variants depending on decomposition / join handling")
    p("  - dbSNP: version drift between Nirvana's bundled dbSNP and the downloaded")
    p("    version may add or remove rsIDs")
    p()

    # --- Discordant examples ---
    has_examples = any(st.examples for st in stats.values()
                       if st.discordant > 0
                       or st.only_in_j2v > 0
                       or st.only_in_direct > 0)
    if has_examples:
        p("=" * 75)
        p("DISCORDANT EXAMPLES (first 10 per field)")
        p("=" * 75)
        for fname in ['gnomAD_AF', 'dbSNP_rsID']:
            st = stats[fname]
            if not st.examples:
                continue
            p(f"\n{fname} ({st.discordant} disagree, "
              f"{st.only_in_j2v} J2V-only, {st.only_in_direct} direct-only):")
            for pos_key, j2v_val, direct_val in st.examples[:10]:
                p(f"  {pos_key}")
                p(f"    nirvana2vcf: {_trunc(j2v_val)}")
                p(f"    direct:   {_trunc(direct_val)}")
        p()

    # --- CONCORDANCE lines ---
    p("=" * 75)
    p("CONCORDANCE SUMMARY")
    p("=" * 75)
    for fname in ['gnomAD_AF', 'dbSNP_rsID']:
        st = stats[fname]
        total = st.total
        pct_val = f"{100.0 * st.concordant / total:.1f}" if total > 0 else "N/A"
        p(f"CONCORDANCE: {fname} {st.concordant} {total} {pct_val}")
    r = summary['af_pearson_r']
    r_str = f"{r:.6f}" if r is not None else "N/A"
    p(f"CONCORDANCE: gnomAD_AF_correlation {r_str}")
    p("=" * 75)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare nirvana2vcf gnomAD AF and dbSNP rsIDs against "
                    "bcftools-direct annotation from source VCFs.",
    )
    parser.add_argument(
        "--nirvana2vcf", required=True,
        help="Path to nirvana2vcf VCF output",
    )
    parser.add_argument(
        "--annotated", required=True,
        help="Path to bcftools-annotated VCF with GNOMAD_DIRECT_AF",
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

    t0 = time.time()
    stats, summary = compare(args.nirvana2vcf, args.annotated, verbose=args.verbose)
    elapsed = time.time() - t0

    out_fh = open(args.output, 'w') if args.output else None
    try:
        write_report(stats, summary, args.nirvana2vcf, args.annotated, elapsed,
                     out=out_fh)
    finally:
        if out_fh:
            out_fh.close()

    if args.output:
        print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
