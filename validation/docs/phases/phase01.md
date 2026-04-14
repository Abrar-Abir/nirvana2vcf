# Phase 1 Completion Report

**Date:** 2026-04-13
**Status:** Done

---

## Objective

Generate Nirvana JSON from public data and verify that nirvana2vcf can consume it end-to-end. This establishes the real-world input data for all subsequent validation phases.

---

## Environment

| Component | Version / Detail |
|---|---|
| macOS | Darwin 25.4.0, Apple Silicon (arm64) |
| Nirvana | v3.18.1 (latest release tag) |
| .NET Runtime | 6.0 |
| Genome Assembly | GRCh37 |
| nirvana2vcf | Development install (`pip install -e .`) |
| Python | 3.9 (system) |

---

## Input Data

| File | Source | Size |
|---|---|---|
| `HiSeq.10000.vcf.gz` | [Nirvana documentation](https://illumina.github.io/NirvanaDocumentation/files/HiSeq.10000.vcf.gz) | 468 KB |

- **Variants:** 9,965 (all SNPs, chr1 only)
- **Samples:** 1 (NA12878)
- **Assembly:** GRCh37

This is Illumina's bundled test file, freely redistributable. It provides a compact, well-characterized dataset covering a single chromosome with dense annotations.

---

## Pipeline Execution

### Step 1: Clone & Build Nirvana

Cloned from `https://github.com/Illumina/Nirvana.git`, checked out tag `v3.18.1`, built with `dotnet build -c Release`.

Two source-level issues required patching before Nirvana could run:

**Banner crash fix** (see [nirvana_banner.md](../fixes/nirvana_banner.md)):
The .NET SDK appends the full git SHA to `InformationalVersion`, making the version string exceed the hardcoded 75-char banner width. `run_phase1.sh` patches `CommandLineUtilities.cs` to truncate at `+` before display.

**BlockCompression arm64 fix** (see [nirvana_block_compression.md](../fixes/nirvana_block_compression.md)):
Nirvana ships a precompiled x86_64-only `libBlockCompression.dylib`. On Apple Silicon this fails at load time. `build_block_compression.sh` compiles an arm64 replacement from a minimal C shim that implements `get_library_id()`, `bgzf_compress()`, `bgzf_decompress()` and statically links libzstd for the ZSTD functions.

Three bugs were encountered and fixed during BlockCompression development:

| Bug | Symptom | Fix |
|---|---|---|
| Magic value mismatch | `get_library_id()` returned `0xCAFEFACE` (-889259314) but Nirvana expects `0xCEFAFECA` (-822411574) | Use decimal literal `return -822411574` |
| Missing gzip magic byte | Output files started with `\x1f\x08` instead of `\x1f\x8b` | Correct header byte 2 from `0x08` to `0x8b` |
| Missing compression method byte | Gzip CM field was absent; FLG byte occupied the CM position | Insert `0x08` (deflate) as byte 3 of BGZF header |

### Step 2: Download Annotation Data

Downloaded GRCh37 annotation sources using Nirvana's Downloader:

```
dotnet Downloader.dll --ga GRCh37 --out nirvana_data/
```

- **Total size:** 29 GB
- **Location:** `validation/data/nirvana_data/`
- **Cache prefix:** `Cache/GRCh37/Both` (Nirvana's `-c` flag takes a path prefix, not a directory)

### Step 3: Run Nirvana Annotation

```
dotnet Nirvana.dll \
    -c nirvana_data/Cache/GRCh37/Both \
    --sd nirvana_data/SupplementaryAnnotation/GRCh37 \
    -r nirvana_data/References/Homo_sapiens.GRCh37.Nirvana.dat \
    -i HiSeq.10000.vcf.gz \
    -o output/HiSeq.10000
```

**Nirvana output:**

```
Initialization                                         Time     Positions/s
---------------------------------------------------------------------------
Cache                                               00:00:01.2
SA Position Scan                                    00:00:00.0      474,298

Reference                                Preload    Annotation   Variants/s
---------------------------------------------------------------------------
chr1                                    00:00:00.1  00:00:00.8       11,252

Time: 00:00:03.8
```

**Output file:** `HiSeq.10000.json.gz` (1.3 MB, BGZF-compressed)

### Step 4: nirvana2vcf Smoke Test

```
nirvana2vcf -i output/HiSeq.10000.json.gz -o output/HiSeq.10000.smoke.vcf
```

Completed without errors.

---

## Output Summary

| Metric | Value |
|---|---|
| Input variants | 9,965 |
| Output variants | 9,965 |
| Variant recovery | **100%** |
| Output file size | 8.5 MB |
| VCF header lines | 77 |
| INFO declarations | 37 |
| FORMAT declarations | 12 |
| Samples | 1 (NA12878) |

### Annotation Coverage

| Annotation | Variants with data | Coverage |
|---|---|---|
| CSQ (transcript consequences) | 6,818 | 68.4% |
| phyloP (conservation) | 9,473 | 95.1% |
| gnomAD AF | 162 | 1.6% |
| REVEL | 49 | 0.5% |
| ClinVar | 24 | 0.2% |
| SpliceAI | 8 | 0.1% |

Low gnomAD/ClinVar/SpliceAI coverage is expected for this dataset -- the HiSeq.10000 test file is chr1 only and contains predominantly intergenic/intronic SNPs. Annotation density will increase with larger, clinically enriched datasets in later validation phases.

### Sample Data Lines

```
#CHROM  POS     ID  REF  ALT  QUAL  FILTER                    INFO
chr1    109     .   A    T    0     FDRtranche2.00to10.00+    CytoBand=1p36.33;GERP=0
chr1    147     .   C    A    0     FDRtranche2.00to10.00+    CytoBand=1p36.33;GERP=0
chr1    177     .   A    C    0     FDRtranche2.00to10.00+    CytoBand=1p36.33;GERP=0
```

---

## Artifacts

| File | Path | Description |
|---|---|---|
| Input VCF | `data/HiSeq.10000.vcf.gz` | Nirvana's public test file (468 KB) |
| Nirvana JSON | `data/output/HiSeq.10000.json.gz` | Annotated output (1.3 MB) |
| Nirvana index | `data/output/HiSeq.10000.json.gz.jsi` | JSON search index (2.5 KB) |
| Smoke test VCF | `data/output/HiSeq.10000.smoke.vcf` | nirvana2vcf output (8.5 MB) |
| Annotation data | `data/nirvana_data/` | GRCh37 cache + supplementary (29 GB, .gitignored) |
| Nirvana source | `data/Nirvana/` | Patched build (.gitignored) |

---

## Scripts

| Script | Purpose |
|---|---|
| `scripts/run_phase1.sh` | Idempotent end-to-end Phase 1 automation |
| `scripts/build_block_compression.sh` | Builds arm64 `libBlockCompression.dylib` for Apple Silicon |

Both scripts are safe to re-run. Each step checks for its output artifact before executing, so completed steps are skipped on subsequent runs.

---

## Issues Encountered

### 1. Cache path mismatch

**Error:** `The transcript cache file GRCh37.transcripts.ndb does not exist`

**Cause:** Script used `-c Cache/Both/GRCh37` but Nirvana's Downloader places files at `Cache/GRCh37/Both.*.ndb`. The `-c` flag expects a prefix (`Cache/GRCh37/Both`), not a directory.

**Fix:** `CACHE_PREFIX="$NIRVANA_DATA/Cache/$GENOME_ASSEMBLY/Both"`

### 2. BlockCompression x86_64/arm64 mismatch

**Error:** `Unable to find the block GZip compression library (BlockCompression)`

**Cause:** Nirvana ships only a precompiled x86_64 `libBlockCompression.dylib`. The .NET runtime on Apple Silicon cannot load it.

**Fix:** Created `build_block_compression.sh` to compile a C replacement targeting arm64. See [nirvana_block_compression.md](../fixes/nirvana_block_compression.md) for full details.

### 3. BGZF header bugs in C shim

**Error:** `gzip.BadGzipFile: Not a gzipped file` / `Unknown compression method`

**Cause:** The BGZF header template in the C shim had two errors: wrong gzip magic byte (`0x08` instead of `0x8b`) and missing compression method byte (`0x08` for deflate). Nirvana wrote output files with invalid gzip block headers that Python's gzip module rejected.

**Fix:** Corrected the 18-byte BGZF header to match the gzip spec: `1f 8b 08 04 ...` (magic, CM=deflate, FLG=FEXTRA, ...).

---

## Conclusion

Phase 1 is complete. The full pipeline -- clone Nirvana, build from source, download annotation data, annotate a public VCF, and convert with nirvana2vcf -- runs end-to-end on Apple Silicon. All 9,965 input variants survive the round trip with annotations populated across all expected categories.

The Nirvana JSON output at `data/output/HiSeq.10000.json.gz` and the smoke test VCF at `data/output/HiSeq.10000.smoke.vcf` are ready for Phase 2 (VCF spec validation) and Phase 3 (position concordance).
