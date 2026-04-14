#!/usr/bin/env python3
"""Compare variant positions between two VCF files.

Extracts (CHROM, POS, REF, ALT) tuples from each VCF, expanding
multi-allelic rows into per-allele tuples, and reports set overlap.

Memory-efficient: processes one chromosome at a time using a merge
approach on two coordinate-sorted VCFs.

Usage:
    python3 scripts/compare_positions.py \
        --nirvana2vcf /tmp/nirvana2vcf_chrY.vcf \
        --reference /tmp/ref_extracts/chrY.vcf \
        [--max-examples 20]
"""

import argparse
import gzip
import sys
import time
from collections import defaultdict


# ---------------------------------------------------------------------------
# VCF reading helpers
# ---------------------------------------------------------------------------

def open_vcf(path):
    """Return a line iterator for a VCF file (plain or gzipped)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def iter_chrom_groups(path, label="file"):
    """Yield (chrom, set_of_(pos,ref,alt)) for each chromosome in a sorted VCF.

    Only one chromosome's data is in memory at a time.
    Also returns row count via a final yield of (None, row_count).
    """
    total_rows = 0
    current_chrom = None
    current_set = set()

    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            total_rows += 1

            fields = line.rstrip("\n").split("\t", 5)
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt_field = fields[4]

            if chrom != current_chrom:
                if current_chrom is not None:
                    print(f"  [{label}] {current_chrom}: {len(current_set):,} tuples",
                          file=sys.stderr)
                    yield (current_chrom, current_set)
                    del current_set
                current_chrom = chrom
                current_set = set()

            if alt_field == ".":
                continue

            for alt in alt_field.split(","):
                alt = alt.strip()
                if alt in (".", "*"):
                    continue
                current_set.add((pos, ref, alt))

    if current_chrom is not None:
        print(f"  [{label}] {current_chrom}: {len(current_set):,} tuples",
              file=sys.stderr)
        yield (current_chrom, current_set)

    yield (None, total_rows)


# ---------------------------------------------------------------------------
# Per-chromosome comparison with merge
# ---------------------------------------------------------------------------

def _analyze_chrom_discrepancies(chrom, only_j2v, only_ref, max_per_chrom=100):
    """Categorize discrepancies for one chromosome with bounded examples."""
    j2v_positions = defaultdict(list)
    for t in only_j2v:
        j2v_positions[t[0]].append(t)

    ref_positions = defaultdict(list)
    for t in only_ref:
        ref_positions[t[0]].append(t)

    shared_positions = set(j2v_positions.keys()) & set(ref_positions.keys())

    pos_matched_j2v = []
    pos_matched_ref = []
    for pos in sorted(shared_positions)[:max_per_chrom]:
        pos_matched_j2v.extend((chrom, *x) for x in j2v_positions[pos])
        pos_matched_ref.extend((chrom, *x) for x in ref_positions[pos])

    truly_unique_j2v = [(chrom, *t) for t in sorted(only_j2v)
                        if t[0] not in shared_positions][:max_per_chrom]
    truly_unique_ref = [(chrom, *t) for t in sorted(only_ref)
                        if t[0] not in shared_positions][:max_per_chrom]

    return {
        "shared_positions": len(shared_positions),
        "pos_matched_j2v": pos_matched_j2v,
        "pos_matched_ref": pos_matched_ref,
        "truly_unique_j2v": truly_unique_j2v,
        "truly_unique_ref": truly_unique_ref,
    }


def compare_per_chrom(j2v_path, ref_path, max_examples):
    """Compare two VCFs using a merge approach: one chromosome at a time.

    Reads nirvana2vcf fully (one chrom at a time), buffers reference chromosomes
    that haven't been matched yet, and compares when both sides have data.
    """
    t0 = time.time()

    # We read nirvana2vcf first (usually smaller), one chrom at a time.
    # For each nirvana2vcf chrom, we advance the reference iterator to find
    # the matching chrom. Reference chroms seen before the match are buffered
    # temporarily (just their names, since we can't compare without nirvana2vcf data).

    # Strategy: read both files completely but one chrom at a time.
    # Since both are sorted, we do a merge.

    # Accumulators
    total_shared = 0
    total_only_j2v = 0
    total_only_ref = 0
    j2v_rows = 0
    ref_rows = 0
    j2v_total_tuples = 0
    ref_total_tuples = 0
    chrom_results = []
    all_pos_matched_j2v = []
    all_pos_matched_ref = []
    all_truly_unique_j2v = []
    all_truly_unique_ref = []

    # Read both files as chrom-group iterators
    print("Reading nirvana2vcf and reference by chromosome (merge)...",
          file=sys.stderr)

    j2v_iter = iter_chrom_groups(j2v_path, "nirvana2vcf")
    ref_iter = iter_chrom_groups(ref_path, "reference")

    # Build dicts of {chrom: set} from both iterators, but only buffer
    # as needed. Since we process sequentially, we keep a small buffer.
    j2v_buffer = {}  # chroms read from j2v but not yet compared
    ref_buffer = {}  # chroms read from ref but not yet compared

    j2v_done = False
    ref_done = False

    def next_j2v():
        nonlocal j2v_rows, j2v_done
        try:
            chrom, data = next(j2v_iter)
            if chrom is None:
                j2v_rows = data
                j2v_done = True
                return None
            return (chrom, data)
        except StopIteration:
            j2v_done = True
            return None

    def next_ref():
        nonlocal ref_rows, ref_done
        try:
            chrom, data = next(ref_iter)
            if chrom is None:
                ref_rows = data
                ref_done = True
                return None
            return (chrom, data)
        except StopIteration:
            ref_done = True
            return None

    def compare_chrom(chrom, j2v_set, ref_set):
        nonlocal total_shared, total_only_j2v, total_only_ref
        nonlocal j2v_total_tuples, ref_total_tuples

        j2v_total_tuples += len(j2v_set)
        ref_total_tuples += len(ref_set)

        # Range-filter reference to nirvana2vcf range
        if j2v_set:
            j2v_min = min(t[0] for t in j2v_set)
            j2v_max = max(t[0] for t in j2v_set)
            ref_set = {t for t in ref_set if j2v_min <= t[0] <= j2v_max}

        shared = j2v_set & ref_set
        only_j = j2v_set - ref_set
        only_r = ref_set - j2v_set

        total_shared += len(shared)
        total_only_j2v += len(only_j)
        total_only_ref += len(only_r)

        if only_j or only_r:
            analysis = _analyze_chrom_discrepancies(chrom, only_j, only_r)
            all_pos_matched_j2v.extend(analysis["pos_matched_j2v"][:50])
            all_pos_matched_ref.extend(analysis["pos_matched_ref"][:50])
            all_truly_unique_j2v.extend(analysis["truly_unique_j2v"][:50])
            all_truly_unique_ref.extend(analysis["truly_unique_ref"][:50])

        total = len(shared) + len(only_j) + len(only_r)
        if total > 0:
            pct = 100.0 * len(shared) / total
            print(f"  >> {chrom}: shared={len(shared):,} only_j2v={len(only_j):,} "
                  f"only_ref={len(only_r):,} ({pct:.2f}%)", file=sys.stderr)

        chrom_results.append({
            "chrom": chrom,
            "j2v": len(j2v_set),
            "ref": len(ref_set),
            "shared": len(shared),
            "only_j2v": len(only_j),
            "only_ref": len(only_r),
        })

        # Free memory
        del shared, only_j, only_r

    # Main merge loop: read chroms from both files and compare when matched
    # We read all chroms from both files. When a chrom appears in both,
    # we compare immediately. Otherwise we track it as only_j2v or only_ref.
    # To avoid buffering too much data, we use a simple approach:
    # read j2v one chrom at a time, for each chrom, drain ref until we find it.

    # Actually, let's use a simpler approach: read all chroms from j2v,
    # but only keep one at a time. Then read all chroms from ref and compare.
    # The issue is memory... Let's do a true merge.

    # Read first entries from both
    j2v_entry = next_j2v()
    ref_entry = next_ref()

    while j2v_entry is not None or ref_entry is not None:
        if j2v_entry is not None and ref_entry is not None:
            j_chrom, j_set = j2v_entry
            r_chrom, r_set = ref_entry

            if j_chrom == r_chrom:
                # Same chrom — compare
                compare_chrom(j_chrom, j_set, r_set)
                del j_set, r_set
                # Also check if this chrom was buffered on either side
                j2v_entry = next_j2v()
                ref_entry = next_ref()
            elif j_chrom < r_chrom:
                # j2v has a chrom that ref hasn't reached yet.
                # Check if it's in the ref buffer
                if j_chrom in ref_buffer:
                    compare_chrom(j_chrom, j_set, ref_buffer.pop(j_chrom))
                else:
                    # Buffer j2v chrom, advance j2v
                    j2v_buffer[j_chrom] = j_set
                j2v_entry = next_j2v()
            else:
                # ref has a chrom that j2v hasn't reached yet
                if r_chrom in j2v_buffer:
                    compare_chrom(r_chrom, j2v_buffer.pop(r_chrom), r_set)
                else:
                    ref_buffer[r_chrom] = r_set
                ref_entry = next_ref()

        elif j2v_entry is not None:
            j_chrom, j_set = j2v_entry
            if j_chrom in ref_buffer:
                compare_chrom(j_chrom, j_set, ref_buffer.pop(j_chrom))
            else:
                # j2v-only chrom
                compare_chrom(j_chrom, j_set, set())
            j2v_entry = next_j2v()

        else:  # ref_entry is not None
            r_chrom, r_set = ref_entry
            if r_chrom in j2v_buffer:
                compare_chrom(r_chrom, j2v_buffer.pop(r_chrom), r_set)
            else:
                compare_chrom(r_chrom, set(), r_set)
            ref_entry = next_ref()

    # Handle any remaining buffered chroms
    for chrom in sorted(set(j2v_buffer.keys()) | set(ref_buffer.keys())):
        j_set = j2v_buffer.pop(chrom, set())
        r_set = ref_buffer.pop(chrom, set())
        compare_chrom(chrom, j_set, r_set)

    # Drain remaining items from iterators (to get row counts)
    if not j2v_done:
        for chrom, data in j2v_iter:
            if chrom is None:
                j2v_rows = data
    if not ref_done:
        for chrom, data in ref_iter:
            if chrom is None:
                ref_rows = data

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.1f}s", file=sys.stderr)

    # --- Print report ---
    _print_report(
        j2v_path=j2v_path,
        ref_path=ref_path,
        j2v_rows=j2v_rows,
        ref_rows=ref_rows,
        j2v_total=j2v_total_tuples,
        ref_total=ref_total_tuples,
        total_shared=total_shared,
        total_only_j2v=total_only_j2v,
        total_only_ref=total_only_ref,
        chrom_results=chrom_results,
        pos_matched_j2v=all_pos_matched_j2v,
        pos_matched_ref=all_pos_matched_ref,
        truly_unique_j2v=all_truly_unique_j2v,
        truly_unique_ref=all_truly_unique_ref,
        max_examples=max_examples,
    )


# ---------------------------------------------------------------------------
# Report formatting
# ---------------------------------------------------------------------------

def pct(num, denom):
    if denom == 0:
        return "N/A"
    return f"{100.0 * num / denom:.2f}%"


def fmt_tuple(t):
    return f"  {t[0]}\t{t[1]}\t{t[2]}\t{t[3]}"


def _print_report(j2v_path, ref_path, j2v_rows, ref_rows,
                  j2v_total, ref_total, total_shared,
                  total_only_j2v, total_only_ref,
                  chrom_results, pos_matched_j2v, pos_matched_ref,
                  truly_unique_j2v, truly_unique_ref, max_examples):
    union_size = total_shared + total_only_j2v + total_only_ref

    print("=" * 60)
    print("POSITION-LEVEL CONCORDANCE REPORT")
    print("=" * 60)
    print(f"nirvana2vcf file:  {j2v_path}")
    print(f"Reference file: {ref_path}")
    print()

    # --- Variant counts ---
    print("VARIANT COUNTS")
    print(f"  nirvana2vcf VCF rows:        {j2v_rows:>12,}")
    print(f"  Reference VCF rows:       {ref_rows:>12,}")
    print(f"  nirvana2vcf variant tuples:  {j2v_total:>12,}  "
          f"(after multi-allelic expansion)")
    print(f"  Reference variant tuples: {ref_total:>12,}  "
          f"(after multi-allelic expansion)")
    print()

    # --- Set comparison ---
    print("SET COMPARISON")
    print(f"  Shared (in both):    {total_shared:>12,}  ({pct(total_shared, union_size)})")
    print(f"  Only in nirvana2vcf:    {total_only_j2v:>12,}  ({pct(total_only_j2v, j2v_total)})")
    print(f"  Only in reference:   {total_only_ref:>12,}  ({pct(total_only_ref, ref_total)})")
    if union_size > 0:
        jaccard = total_shared / union_size
        print(f"  Jaccard similarity:  {jaccard:>12.6f}")
    print()

    # --- Per-chromosome breakdown ---
    print("PER-CHROMOSOME BREAKDOWN")
    print(f"  {'CHROM':<8} {'nirvana2vcf':>12} {'reference':>12} {'shared':>12} "
          f"{'only_j2v':>12} {'only_ref':>12} {'overlap%':>10}")
    for cr in chrom_results:
        total = cr["shared"] + cr["only_j2v"] + cr["only_ref"]
        ov = f"{100.0 * cr['shared'] / total:.2f}%" if total > 0 else "N/A"
        print(f"  {cr['chrom']:<8} {cr['j2v']:>12,} {cr['ref']:>12,} "
              f"{cr['shared']:>12,} {cr['only_j2v']:>12,} {cr['only_ref']:>12,} {ov:>10}")
    print()

    # --- Discrepancy examples ---
    if pos_matched_j2v:
        n = min(max_examples, len(pos_matched_j2v))
        print(f"EXAMPLES: Position-matched, allele-different (first {n})")
        seen = set()
        count = 0
        for t in pos_matched_j2v:
            key = (t[0], t[1])
            if key in seen:
                continue
            seen.add(key)
            j2v_at_pos = [x for x in pos_matched_j2v if (x[0], x[1]) == key]
            ref_at_pos = [x for x in pos_matched_ref if (x[0], x[1]) == key]
            print(f"  POS {key[0]}:{key[1]}")
            for x in j2v_at_pos:
                print(f"    nirvana2vcf:  REF={x[2]}  ALT={x[3]}")
            for x in ref_at_pos:
                print(f"    reference: REF={x[2]}  ALT={x[3]}")
            count += 1
            if count >= max_examples:
                break
        print()

    if truly_unique_j2v:
        n = min(max_examples, len(truly_unique_j2v))
        print(f"EXAMPLES: Truly unique to nirvana2vcf (first {n})")
        for t in truly_unique_j2v[:n]:
            print(fmt_tuple(t))
        print()

    if truly_unique_ref:
        n = min(max_examples, len(truly_unique_ref))
        print(f"EXAMPLES: Truly unique to reference (first {n})")
        for t in truly_unique_ref[:n]:
            print(fmt_tuple(t))
        print()

    # --- Pass/fail summary ---
    overlap_pct = 100.0 * total_shared / j2v_total if j2v_total else 0
    print("-" * 60)
    if overlap_pct >= 99.0:
        print(f"RESULT: PASS  (overlap = {overlap_pct:.2f}% of nirvana2vcf tuples)")
    else:
        print(f"RESULT: REVIEW NEEDED  (overlap = {overlap_pct:.2f}% "
              f"of nirvana2vcf tuples, target >= 99%)")
    print("-" * 60)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare variant positions between two VCF files.",
    )
    parser.add_argument(
        "--nirvana2vcf", required=True,
        help="Path to nirvana2vcf output VCF (plain or .gz)",
    )
    parser.add_argument(
        "--reference", required=True,
        help="Path to reference VCF extract (plain or .gz)",
    )
    parser.add_argument(
        "--max-examples", type=int, default=20,
        help="Max examples to print per discrepancy category (default: 20)",
    )
    args = parser.parse_args()

    compare_per_chrom(args.nirvana2vcf, args.reference, args.max_examples)


if __name__ == "__main__":
    main()
