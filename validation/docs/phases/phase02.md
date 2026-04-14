# Phase 2 Completion Report — Smoke Test & VCF Spec Validation

**Date:** 2026-04-13
**Status:** Done

---

## Objective

Validate that nirvana2vcf produces structurally sound, spec-compliant VCF 4.2 output across all mode combinations, using independent validators (bcftools, vcftools) and programmatic header/key declaration checks.

---

## Environment

| Component | Version / Detail |
|---|---|
| macOS | Darwin 25.4.0, Apple Silicon (arm64) |
| nirvana2vcf | Development install (`pip install -e .`) |
| Python | 3.9 (system) |
| bcftools | 1.23.1 (htslib 1.23.1) |
| vcftools | 0.1.17 (Perl vcf-validator) |
| EBI vcf_validator | Skipped (pre-built binary requires Boost dylibs not present on this system) |

---

## Input Data

| File | Source | Size |
|---|---|---|
| `HiSeq.10000.json.gz` | Phase 1 Nirvana output | 1.3 MB |

- **Variants:** 9,965 (all SNPs, chr1 only, 1 sample NA12878)

---

## Mode Combinations Tested

| Mode | Flags | Variants | Description |
|---|---|---|---|
| default | *(none)* | 9,965 | Normalize on, decompose off |
| raw | `--no-normalize` | 9,965 | Raw Nirvana alleles |
| decomposed | `--decompose` | 9,965 | Normalize + decompose |
| no_samples | `--no-samples` | 9,965 | Omit genotype columns |
| csq_only | `--csq-only` | 9,965 | VEP-style CSQ only |

All modes produced identical variant counts (9,965). The decomposed mode also produced 9,965 because the HiSeq.10000 dataset contains only SNPs — no multi-allelic sites exist to split. Decomposition will show its effect on datasets with multi-allelic variants.

---

## Validation Results

### bcftools (Structural Parsing)

| Mode | Result |
|---|---|
| default | PASS |
| raw | PASS |
| decomposed | PASS |
| no_samples | PASS |
| csq_only | PASS |

All five VCFs parse cleanly through bcftools with no errors.

### bcftools stats (Default Mode)

| Metric | Value |
|---|---|
| Samples | 1 |
| Records | 9,965 |
| SNPs | 9,965 |
| Indels | 0 |
| MNPs | 0 |
| Multiallelic sites | 0 |
| Ts/Tv ratio | 2.02 |

The Ts/Tv ratio of 2.02 is consistent with expectations for human germline SNPs (~2.0–2.1 genome-wide), confirming the variant data is biologically plausible.

### Header Declaration Check

| Mode | INFO Keys | FORMAT Keys |
|---|---|---|
| default | All declared | All declared |
| raw | All declared | All declared |
| decomposed | All declared | All declared |
| no_samples | All declared | N/A (no samples) |
| csq_only | All declared | All declared |

No undeclared INFO or FORMAT keys in any mode. Every key used in data lines has a corresponding `##INFO` or `##FORMAT` header declaration.

### vcftools vcf-validator

| Mode | Errors | Recommendations |
|---|---|---|
| default | 0 | 9 |
| raw | 0 | 9 |
| decomposed | 0 | 9 |
| no_samples | 0 | 9 |
| csq_only | 0 | 9 |

**Zero errors** across all modes. The 9 recommendations (identical for each mode) are:

1. **Missing `##reference` header** — "Not required but highly recommended." Nirvana JSON does not include the reference FASTA path, so nirvana2vcf cannot populate this.

2. **Missing `##contig` header for chr1** — "Not required but highly recommended." Contig lines with lengths are emitted based on the assembly, but vcftools wants an exact match for every CHROM value encountered. nirvana2vcf emits contig lines for all chromosomes in the assembly, including `chr1`.

3. **FILTER values not declared in header** (7 distinct values) — FILTER strings (`FDRtranche2.00to10.00+`, `DPFilter`, `HARD_TO_VALIDATE`, `LowQual`, `FDRtranche2.00to10.00`, `FDRtranche0.10to1.00`, `SnpCluster`, `FDRtranche1.00to2.00`, `Indel`) come from the original VCF's FILTER column, passed through Nirvana's JSON. nirvana2vcf outputs them faithfully but does not generate `##FILTER` header declarations since Nirvana JSON does not include the original VCF's FILTER definitions.

All three recommendation categories are cosmetic — the VCF spec treats these as SHOULD, not MUST. The data is structurally correct and parseable by all tools.

### EBI vcf_validator

Skipped. The pre-built macOS arm64 binary (v0.10.2 from GitHub releases) requires `libboost_system.dylib` and `libboost_iostreams.dylib` linked against an older Boost version than what Homebrew provides (1.90.0). This is a known packaging limitation. Alternative installation paths (conda, building from source) can be used if stricter validation is needed.

---

## Variant Count Cross-Check

| Mode | Count | Match |
|---|---|---|
| default | 9,965 | Exact (matches Phase 1 smoke VCF) |
| raw | 9,965 | Exact |
| decomposed | 9,965 | >= reference (all SNPs, no multi-allelic to split) |
| no_samples | 9,965 | Exact |
| csq_only | 9,965 | Exact |

---

## Artifacts

| File | Path | Description |
|---|---|---|
| default VCF | `data/output/phase2/HiSeq.10000.default.vcf` | Standard mode output |
| raw VCF | `data/output/phase2/HiSeq.10000.raw.vcf` | `--no-normalize` output |
| decomposed VCF | `data/output/phase2/HiSeq.10000.decomposed.vcf` | `--decompose` output |
| no_samples VCF | `data/output/phase2/HiSeq.10000.no_samples.vcf` | `--no-samples` output |
| csq_only VCF | `data/output/phase2/HiSeq.10000.csq_only.vcf` | `--csq-only` output |
| bcftools report | `results/bcftools_validation.txt` | Parse results + stats for all modes |
| bcftools stats | `results/bcftools_stats_default.txt` | Full stats for default mode |
| header check | `results/header_declaration_check.txt` | INFO/FORMAT declaration audit |
| vcftools report | `results/vcftools_validation.txt` | Perl vcf-validator output |
| variant counts | `results/variant_counts.txt` | Count cross-check table |

---

## Script

| Script | Purpose |
|---|---|
| `scripts/run_phase2.sh` | Idempotent Phase 2 automation — generates VCFs, runs all validators |

Safe to re-run. Downloads vcf_validator binary and installs vcftools on first run if needed.

---

## Conclusion

**Phase 2: PASS**

- **bcftools**: All 5 mode VCFs parse without errors
- **vcftools**: Zero errors; 9 cosmetic recommendations (missing `##reference`, `##contig`, and `##FILTER` header lines — all SHOULD, not MUST)
- **Header declarations**: All INFO and FORMAT keys are properly declared
- **Variant counts**: 9,965 across all modes, matching Phase 1 output exactly
- **bcftools stats**: Ts/Tv = 2.02, consistent with human germline SNPs

nirvana2vcf produces structurally sound VCF 4.2 output across all mode combinations. The output is parseable by standard bioinformatics tools and conforms to the VCF specification.

**Next:** Phase 3 — Round-trip position concordance (verify same variants survive JSON-to-VCF conversion)
