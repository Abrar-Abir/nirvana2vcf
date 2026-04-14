# Phase 3 Completion Report — Round-Trip Position Concordance

**Date:** 2026-04-13
**Status:** Done

---

## Objective

Verify that the same (CHROM, POS, REF, ALT) variant tuples survive the VCF -> Nirvana JSON -> nirvana2vcf VCF round-trip. This is a self-consistency check: since nirvana2vcf converts the Nirvana JSON that was produced from the original VCF, every variant position should be preserved.

---

## Environment

| Component | Version / Detail |
|---|---|
| macOS | Darwin 25.4.0, Apple Silicon (arm64) |
| nirvana2vcf | Development install (`pip install -e .`) |
| Python | 3.9 (system) |
| bcftools | 1.23.1 (htslib 1.23.1) |

---

## Input Data

| File | Source | Site Tuples |
|---|---|---|
| `HiSeq.10000.vcf.gz` | Original Nirvana test VCF | 9,965 |
| `HiSeq.10000.{mode}.vcf` | Phase 2 nirvana2vcf outputs | 9,965 each |

---

## Dataset Characteristics

- **9,965 SNPs** on chr1, single sample (NA12878)
- **No indels**, no multi-allelic sites
- All mode VCFs produce byte-identical position columns (verified in Phase 2)
- Normalization and decomposition have no effect on this all-SNP dataset

---

## Methodology

Two independent comparison methods were used for each mode:

### Method 1: bcftools + comm

1. Extract (CHROM, POS, REF, ALT) tuples from both the original VCF and each nirvana2vcf mode VCF using `bcftools query`
2. Sort both files with `LC_ALL=C sort` for locale-independent ordering
3. Use `comm` to compute set intersection and differences:
   - `comm -12`: shared tuples
   - `comm -23`: only in original
   - `comm -13`: only in nirvana2vcf

### Method 2: Python (compare_positions.py)

A memory-efficient Python script that:
1. Reads both VCFs one chromosome at a time using a merge approach
2. Expands multi-allelic rows into per-allele (POS, REF, ALT) tuples
3. Computes set overlap per chromosome and reports Jaccard similarity
4. Categorizes discrepancies into position-matched/allele-different and truly unique

---

## Results

### bcftools + comm Concordance

| Mode | nirvana2vcf Tuples | Shared | Only Original | Only nirvana2vcf | Result |
|---|---|---|---|---|---|
| default | 9,965 | 9,965 | 0 | 0 | PASS |
| raw | 9,965 | 9,965 | 0 | 0 | PASS |
| decomposed | 9,965 | 9,965 | 0 | 0 | PASS |
| no_samples | 9,965 | 9,965 | 0 | 0 | PASS |

### Python Concordance (compare_positions.py)

| Mode | nirvana2vcf Rows | Reference Rows | nirvana2vcf Tuples | Reference Tuples | Shared | Only nirvana2vcf | Only Reference | Jaccard | Result |
|---|---|---|---|---|---|---|---|---|---|
| default | 9,965 | 9,965 | 9,965 | 9,965 | 9,965 | 0 | 0 | 1.000000 | PASS |
| raw | 9,965 | 9,965 | 9,965 | 9,965 | 9,965 | 0 | 0 | 1.000000 | PASS |
| decomposed | 9,965 | 9,965 | 9,965 | 9,965 | 9,965 | 0 | 0 | 1.000000 | PASS |
| no_samples | 9,965 | 9,965 | 9,965 | 9,965 | 9,965 | 0 | 0 | 1.000000 | PASS |

All modes achieve **100% concordance** — Jaccard similarity = 1.000000, zero discrepancies.

---

## Modes Skipped

| Mode | Flags | Reason |
|---|---|---|
| csq_only | `--csq-only` | Same positions as default (only affects INFO field content) |
| raw + decomposed | `--no-normalize --decompose` | Identical to raw for all-SNP data (no multi-allelic sites to split, no alleles to normalize) |

---

## Artifacts

| File | Path | Description |
|---|---|---|
| Original sites | `results/phase3/original_sites.tsv` | Extracted (CHROM, POS, REF, ALT) from original VCF |
| Mode sites | `results/phase3/{mode}_sites.tsv` | Extracted sites per mode (4 files) |
| bcftools report | `results/phase3/bcftools_concordance.txt` | comm-based set comparison results |
| Python reports | `results/phase3/{mode}_position_concordance.txt` | Full concordance reports (4 files) |
| Python logs | `results/phase3/{mode}_position_concordance.log` | Per-chromosome progress logs (4 files) |
| Summary | `results/phase3/summary.txt` | Combined pass/fail table |

---

## Script

| Script | Purpose |
|---|---|
| `scripts/compare_positions.py` | Memory-efficient per-chromosome VCF position comparison |
| `scripts/run_phase3.sh` | Idempotent Phase 3 automation — runs both comparison methods |

Safe to re-run. Overwrites results on each run.

---

## Conclusion

**Phase 3: PASS**

- **bcftools + comm**: 9,965/9,965 shared tuples across all 4 modes, zero only-original, zero only-nirvana2vcf
- **Python comparison**: Jaccard similarity = 1.000000 across all 4 modes, zero discrepancies
- **Two independent methods agree**: every variant position in the original VCF survives the Nirvana JSON -> nirvana2vcf round-trip

The result is expected for this all-SNP dataset: normalization and decomposition have no effect on simple SNPs, so all modes produce identical position columns. Datasets with indels and multi-allelic sites (Phases 4+, or larger public datasets) would exercise normalization and decomposition more thoroughly.

**Next:** Phase 4 — Round-trip annotation concordance (verify annotation values are faithfully converted)
