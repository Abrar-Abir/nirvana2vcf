# Phase 4 Completion Report — Round-Trip Annotation Concordance

**Date:** 2026-04-13
**Status:** Done

---

## Objective

Verify that annotation values are faithfully converted from Nirvana JSON to VCF. This is a self-consistency check: both the JSON source and the VCF output derive from the same Nirvana run, so the expectation is **100% concordance**. Any mismatch is a nirvana2vcf bug.

---

## Environment

| Component | Version / Detail |
|---|---|
| macOS | Darwin 25.4.0, Apple Silicon (arm64) |
| nirvana2vcf | Development install (`pip install -e .`) |
| Python | 3.9 (system) |
| orjson | Fast JSON parsing (same as nirvana2vcf) |

---

## Input Data

| File | Source | Positions |
|---|---|---|
| `HiSeq.10000.json.gz` | Nirvana-annotated JSON (Phase 1) | 9,965 |
| `HiSeq.10000.{mode}.vcf` | Phase 2 nirvana2vcf outputs | 9,965 each |

---

## Methodology

A standalone Python script (`compare_annotations.py`) independently parses the Nirvana JSON and compares expected annotation values against actual VCF output, field by field.

### Approach

1. Stream Nirvana JSON and VCF in lockstep (1:1 position mapping)
2. For each position, reconstruct expected VCF INFO values from raw JSON using the same formatting rules as `mapper.py` (`.6g` floats, `str` ints, percent-escaping, per-allele ordering)
3. Parse actual VCF INFO values
4. Compare expected vs actual for each field — **exact string match** (no tolerance)

### Fields Compared

| Category | Fields | Count |
|---|---|---|
| Population frequencies | gnomAD (AF, AC, AN, AFR, AMR, EUR/NFE, EAS, SAS), oneKG (AF, AFR, AMR, EUR, EAS, SAS), TOPMed AF | 15 |
| Pathogenicity scores | phyloP, DANN, GERP, REVEL | 4 |
| ClinVar | CLINVAR_ID, CLINVAR_SIG, CLINVAR_REVSTAT | 3 |
| SpliceAI | AG/AL/DG/DL SCORE and DIST | 8 |
| Structural variants | SVTYPE, SVEND, SVLEN, CIPOS, CIEND | 5 |
| Position-level | CytoBand | 1 |
| dbSNP | ID column | 1 |
| CSQ | Full VEP-style CSQ string (19 pipe-delimited subfields per transcript) | 1 |
| **Total** | | **38** |

---

## Results

### Per-Mode Summary

| Mode | Fields Tested | Values Verified | Discordant | Result |
|---|---|---|---|---|
| default | 25 | 45,661 | 0 | **PASS** |
| raw | 25 | 45,661 | 0 | **PASS** |
| no_samples | 25 | 45,661 | 0 | **PASS** |
| csq_only | 2 | 7,428 | 0 | **PASS** |

**Cross-mode consistency confirmed:** default, raw, and no_samples produce identical concordance numbers (45,661 values, 25 fields with data). This is expected — annotation values are independent of normalization and sample inclusion.

### Per-Field Coverage (default mode)

| Annotation | Positions with Data | Concordant | Discordant | Result |
|---|---|---|---|---|
| CytoBand | 9,965 | 9,965 | 0 | PASS |
| phyloP | 9,473 | 9,473 | 0 | PASS |
| DANN | 7,212 | 7,212 | 0 | PASS |
| GERP | 9,965 | 9,965 | 0 | PASS |
| REVEL | 49 | 49 | 0 | PASS |
| gnomAD_AF | 162 | 162 | 0 | PASS |
| gnomAD_AC | 162 | 162 | 0 | PASS |
| gnomAD_AN | 162 | 162 | 0 | PASS |
| gnomAD_AFR_AF | 162 | 162 | 0 | PASS |
| gnomAD_AMR_AF | 162 | 162 | 0 | PASS |
| gnomAD_EUR_AF | 162 | 162 | 0 | PASS |
| gnomAD_EAS_AF | 162 | 162 | 0 | PASS |
| gnomAD_SAS_AF | 19 | 19 | 0 | PASS |
| oneKG_AF | 41 | 41 | 0 | PASS |
| oneKG_AFR_AF | 41 | 41 | 0 | PASS |
| oneKG_AMR_AF | 41 | 41 | 0 | PASS |
| oneKG_EUR_AF | 41 | 41 | 0 | PASS |
| oneKG_EAS_AF | 41 | 41 | 0 | PASS |
| oneKG_SAS_AF | 41 | 41 | 0 | PASS |
| TOPMed_AF | 98 | 98 | 0 | PASS |
| CLINVAR_ID | 24 | 24 | 0 | PASS |
| CLINVAR_SIG | 24 | 24 | 0 | PASS |
| CLINVAR_REVSTAT | 24 | 24 | 0 | PASS |
| ID (dbSNP) | 610 | 610 | 0 | PASS |
| CSQ | 6,818 | 6,818 | 0 | PASS |

### Fields with No Data (expected)

| Field | Reason |
|---|---|
| SpliceAI (8 fields) | No SpliceAI annotations in HiSeq.10000 dataset |
| SVTYPE, SVEND, SVLEN, CIPOS, CIEND | No structural variants (all-SNP dataset) |

### Coverage vs Expected

| Annotation | Expected (approx.) | Actual | Match |
|---|---|---|---|
| CytoBand | (not estimated) | 9,965 | - |
| phyloP | ~8,962 | 9,473 | Close |
| DANN | (not estimated) | 7,212 | - |
| GERP | (not estimated) | 9,965 | - |
| CSQ transcripts | ~6,818 | 6,818 | Exact |
| dbSNP rsID | ~610 | 610 | Exact |
| gnomAD | ~162 | 162 | Exact |
| TOPMed | ~98 | 98 | Exact |
| REVEL | ~49 | 49 | Exact |
| oneKG | ~41 | 41 | Exact |
| ClinVar | ~24 | 24 | Exact |
| SpliceAI | 0 | 0 | Exact |

---

## Modes Skipped

| Mode | Flags | Reason |
|---|---|---|
| decomposed | `--decompose` | Identical to default for all-SNP data (no multi-allelic sites to split) |

---

## Artifacts

| File | Path | Description |
|---|---|---|
| Concordance reports | `results/phase4/{mode}_annotation_concordance.txt` | Per-field results (4 files) |
| Comparison logs | `results/phase4/{mode}_annotation_concordance.log` | Progress logs (4 files) |
| Summary | `results/phase4/summary.txt` | Combined pass/fail table with coverage |

---

## Scripts

| Script | Purpose |
|---|---|
| `scripts/compare_annotations.py` | Standalone annotation comparison — streams JSON & VCF in lockstep, exact string match on all fields |
| `scripts/run_phase4.sh` | Idempotent Phase 4 automation — runs comparison for all modes, generates summary |

Safe to re-run. Overwrites results on each run.

---

## Conclusion

**Phase 4: PASS**

- **100% concordance** across all 25 fields with data, all 4 modes tested
- **45,661 individual field values verified** per mode (default/raw/no_samples), 7,428 for csq_only
- **Zero discordant examples** — every annotation value in the VCF exactly matches the source JSON
- **Coverage matches expectations** — gnomAD, oneKG, TOPMed, ClinVar, REVEL, dbSNP, and CSQ counts match the dataset characteristics

The self-consistency check confirms that nirvana2vcf faithfully converts all annotation values from Nirvana JSON to VCF format, including:
- Per-allele population frequencies (gnomAD, 1000 Genomes, TOPMed)
- Pathogenicity scores (phyloP, DANN, GERP, REVEL) with `.6g` float formatting
- ClinVar fields with percent-encoding
- VEP-style CSQ strings with all 19 subfields (Consequence, SYMBOL, PolyPhen, SIFT, HGVS, etc.)
- dbSNP rsIDs with deduplication

**Next:** Phase 5 — Cross-tool comparison against VEP, SnpEff, and bcftools (independent validation)
