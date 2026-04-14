# Open-Source Validation Strategy for nirvana2vcf

This document describes how to validate nirvana2vcf using **only public data and open-source tools**, suitable for publishing reproducible integrity checks alongside the open-source release.

---

## Status

| Phase | Status | Notes |
|---|---|---|
| 1. Generate Nirvana JSON | Done | |
| 2. VCF Spec Validation | Done | |
| 3. Position Concordance | Done | |
| 4. Annotation Concordance | Done | |
| 5. Cross-Tool Comparison | Done | |

---

## Directory Layout

```
validation/
  docs/
    strategy.md      # this file
    phases/          # phase execution reports (phase01–05.md)
    fixes/           # issue fix documentation
  scripts/           # comparison and validation scripts
  data/              # input VCFs, Nirvana JSON outputs (.gitignored)
  results/           # reports and concordance summaries
```

---

## Overview

The validation has five phases, each building on the previous:

| Phase | What It Proves | Tools Needed |
|---|---|---|
| 1. Generate Nirvana JSON from Public VCFs | You have real-world input data | Nirvana (open-source), .NET 6+ |
| 2. Smoke Test & VCF Spec Validation | nirvana2vcf produces spec-compliant VCF | vcf-validator, bcftools |
| 3. Round-Trip Position Concordance | Same variants survive the JSON-to-VCF conversion | bcftools, diff |
| 4. Round-Trip Annotation Concordance | Annotation values are faithfully converted | Python comparison script |
| 5. Cross-Tool Comparison | nirvana2vcf output agrees with independent annotators | VEP, SnpEff/SnpSift, bcftools annotate |

---

## Phase 1 — Generate Nirvana JSON from Public Data

### 1a. Quick Start: Nirvana's Bundled Test VCF

Nirvana ships a 10,000-variant test file. This is the fastest way to get a real Nirvana JSON:

```bash
# Download and run Nirvana's own test script (requires .NET 6+)
curl -O https://illumina.github.io/NirvanaDocumentation/files/TestNirvana.sh
bash ./TestNirvana.sh
# Produces: HiSeq.10000.json.gz  (Nirvana JSON output)
# Input:    HiSeq.10000.vcf.gz   (original VCF — keep this for comparison)
```

### 1b. Larger/More Diverse Public VCFs

For broader coverage, annotate one or more of these public datasets:

| Dataset | Source | Size | Why Use It |
|---|---|---|---|
| Genome in a Bottle NA12878 (HG001) | [NIST GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/) | ~4M variants | Gold-standard truth set, GRCh37 & GRCh38 |
| 1000 Genomes Phase 3 (single chr) | [IGSR Data Portal](https://www.internationalgenome.org/data-portal/sample/NA12878) | Varies per chr | Multi-sample, population diversity |
| Illumina Platinum Genomes | [GitHub](https://github.com/Illumina/PlatinumGenomes) | CEPH trio | Multi-sample with pedigree (tests FORMAT fields) |
| ClinVar VCF | [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) | ~2M entries | Dense ClinVar annotations |

**Subsetting large VCFs for manageable test sizes:**

```bash
# Extract chr21 (small chromosome) from a large VCF
bcftools view -r chr21 input.vcf.gz -Oz -o chr21_subset.vcf.gz
bcftools index chr21_subset.vcf.gz

# Or take just the first 5000 variants
bcftools view -h input.vcf.gz > subset.vcf
bcftools view -H input.vcf.gz | head -5000 >> subset.vcf
bgzip subset.vcf
```

### 1c. Running Nirvana on Any VCF

```bash
# Clone and build
git clone https://github.com/Illumina/Nirvana.git && cd Nirvana
dotnet build -c Release

# Download annotation data sources (~30 GB for GRCh38, one-time)
dotnet bin/Release/net6.0/Downloader.dll --ga GRCh38 -o Data

# Annotate
dotnet bin/Release/net6.0/Nirvana.dll \
  -c Data/Cache/GRCh38 \
  --sd Data/SupplementaryAnnotation/GRCh38 \
  -r Data/References/Homo_sapiens.GRCh38.Nirvana.dat \
  -i your_input.vcf.gz \
  -o annotated_output
# Produces: annotated_output.json.gz
```

For GRCh37, substitute `GRCh37` in the paths above.

---

## Phase 2 — Smoke Test & VCF Spec Validation

### 2a. Convert Nirvana JSON to VCF

Test all relevant mode combinations:

```bash
# Default (normalize on, decompose off)
nirvana2vcf -i HiSeq.10000.json.gz -o output.vcf

# All mode combinations
nirvana2vcf -i HiSeq.10000.json.gz -o output_raw.vcf --no-normalize
nirvana2vcf -i HiSeq.10000.json.gz -o output_decomp.vcf --decompose
nirvana2vcf -i HiSeq.10000.json.gz -o output_norm_decomp.vcf --decompose
nirvana2vcf -i HiSeq.10000.json.gz -o output_nosamp.vcf --no-samples
nirvana2vcf -i HiSeq.10000.json.gz -o output_csq.vcf --csq-only
```

### 2b. Validate with EBI vcf-validator (Strictest)

The [EBI vcf-validator](https://github.com/EBIvariation/vcf-validator) checks full VCF 4.2 spec compliance — header syntax, INFO field declarations, data types, and structural rules.

```bash
# Install
conda install -c bioconda vcf-validator
# or download binary from GitHub releases

# Run validation
vcf_validator -i output.vcf
# Exit code 0 = valid; non-zero = spec violations found
# Reports warnings and errors with line numbers
```

### 2c. Validate with bcftools

```bash
conda install -c bioconda bcftools

# Parse check — if this succeeds, the VCF is structurally sound
bcftools view output.vcf > /dev/null && echo "PASS: parseable VCF"

# Summary statistics
bcftools stats output.vcf
# Check: variant count, ts/tv ratio, indel count, multiallelic count

# List samples (verify they match expectations)
bcftools query -l output.vcf

# Verify all INFO keys are declared in the header
bcftools view -h output.vcf | grep '^##INFO' | awk -F'[<,]' '{print $3}' | sort > declared_info.txt
bcftools query -f '%INFO/KEYS\n' output.vcf 2>/dev/null | tr ';' '\n' | sort -u > used_info.txt
comm -23 used_info.txt declared_info.txt  # should be empty
```

### 2d. Validate with VCFtools

```bash
conda install -c bioconda vcftools
vcf-validator output.vcf
```

### Success Criteria

- All three validators report no errors
- `bcftools stats` shows expected variant counts (close to input VCF variant count)
- No undeclared INFO/FORMAT keys

---

## Phase 3 — Round-Trip Position Concordance

Since you know the **original VCF** that was fed to Nirvana, you can verify that the same variants come back out of nirvana2vcf.

### 3a. Extract Sites from Both VCFs

```bash
# Original VCF sites
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' HiSeq.10000.vcf.gz | sort > original_sites.tsv

# nirvana2vcf output sites
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' output.vcf | sort > roundtrip_sites.tsv

# Compare
diff original_sites.tsv roundtrip_sites.tsv
wc -l original_sites.tsv roundtrip_sites.tsv
comm -3 original_sites.tsv roundtrip_sites.tsv | head -20
```

### 3b. Quantify Overlap

```bash
comm -12 original_sites.tsv roundtrip_sites.tsv | wc -l  # shared
comm -23 original_sites.tsv roundtrip_sites.tsv | wc -l  # only in original
comm -13 original_sites.tsv roundtrip_sites.tsv | wc -l  # only in nirvana2vcf
```

### 3c. Test All Mode Combinations

Run position concordance for each mode:

| Mode | Flags | Expected Behavior |
|---|---|---|
| Raw | `--no-normalize` | Positions match original exactly |
| Normalized (default) | (none) | POS may shift from prefix/suffix trimming |
| Decomposed | `--no-normalize --decompose` | One original multi-allelic row becomes N biallelic rows |
| Normalized + Decomposed | `--decompose` | Both effects combined |

### 3d. Handling Expected Differences

Not all differences are bugs:

| Difference | Cause | How to Check |
|---|---|---|
| POS differs by a few bases | Allele normalization (`--normalize` trims shared prefix/suffix, shifts POS) | Re-run with `--no-normalize` and compare again |
| One original row becomes multiple | `--decompose` splits multi-allelics into biallelic rows | Count ALTs in original; decomposed count should equal sum of ALTs |
| Variant missing from nirvana2vcf | Nirvana dropped it (e.g., symbolic alleles `<NON_REF>`, reference-only rows) | Check the Nirvana JSON — if variant is absent there, it is a Nirvana filter, not a nirvana2vcf bug |

### Success Criteria

- **100% position overlap** between original VCF and nirvana2vcf output (with `--no-normalize`)
- All discrepancies are explainable (symbolic alleles, reference blocks, etc.)

---

## Phase 4 — Round-Trip Annotation Concordance

Since both the Nirvana JSON and the nirvana2vcf output originate from the **same Nirvana run**, annotation values should match **exactly** (no version-mismatch noise). This is a stricter test than comparing against a different pipeline.

### 4a. What to Compare

Parse the Nirvana JSON directly and compare field-by-field against the nirvana2vcf VCF output:

| Annotation | Nirvana JSON Path | nirvana2vcf VCF Field | Expected Match |
|---|---|---|---|
| gnomAD AF | `variants[].gnomad.allAf` | `gnomAD_AF` | Exact (after `.6g` float formatting) |
| gnomAD pop AFs | `variants[].gnomad.afrAf`, `.amrAf`, `.easAf`, `.nfeAf`, `.sasAf` | `gnomAD_AFR_AF`, `_AMR_AF`, `_EAS_AF`, `_EUR_AF`, `_SAS_AF` | Exact |
| TOPMed AF | `variants[].topmed.allAf` | `TOPMed_AF` | Exact |
| ClinVar significance | `variants[].clinvar[].significance[]` | `CLINVAR_SIG` | Exact (after escaping) |
| ClinVar review status | `variants[].clinvar[].reviewStatus` | `CLINVAR_REVSTAT` | Exact |
| dbSNP rsID | `variants[].dbsnp[]` | `ID` column | Exact |
| REVEL score | `variants[].revel.score` or raw float | `REVEL` | Exact (after `.6g` formatting) |
| SpliceAI scores | `variants[].spliceAI[].acceptorGain`, etc. | `SpliceAI_AG`, etc. | Exact |
| phyloP score | `variants[].phylopScore` | `phyloP` | Exact |
| Consequence terms | `variants[].transcripts.refSeq[].consequence[]` | CSQ `Consequence` subfield | Exact |
| Gene symbol | `variants[].transcripts.refSeq[].hgnc` | CSQ `SYMBOL` subfield | Exact |
| SIFT prediction | `variants[].transcripts.refSeq[].siftPrediction` | CSQ `SIFT` subfield | Exact |
| PolyPhen prediction | `variants[].transcripts.refSeq[].polyPhenPrediction` | CSQ `PolyPhen` subfield | Exact |
| HGVS coding | `variants[].transcripts.refSeq[].hgvsc` | CSQ `HGVSc` subfield | Exact |
| HGVS protein | `variants[].transcripts.refSeq[].hgvsp` | CSQ `HGVSp` subfield | Exact |

**Important implementation notes:**

- `gnomAD_EUR_AF` maps to Nirvana's `nfeAf` (non-Finnish European), **not** `eurAf` (which is 1000 Genomes EUR and typically absent in gnomAD blocks).
- PolyPhen predictions in Nirvana use **HVAR** (HumanVar) scores, not HDIV.
- REVEL scores can arrive as `{"score": X}` dict or raw float — the comparison script must handle both forms.
- All float comparisons should apply `.6g` formatting to both sides before comparing.

### 4b. Comparison Script

Write `validation/scripts/compare_roundtrip.py`:

```python
#!/usr/bin/env python3
"""compare_roundtrip.py — Verify nirvana2vcf faithfully converts Nirvana JSON.

Streams both the Nirvana JSON and nirvana2vcf VCF output in parallel.
For each (CHROM, POS, REF, ALT), compares every mapped field from the
table above. Reports: exact matches, mismatches with values, missing fields.

Must handle:
- .6g float formatting for numeric comparisons
- INFO value escaping (percent-encoded %, space, ;, =, ,)
- Per-allele (Number=A) field ordering for multi-allelic sites
- CSQ pipe-delimited subfield parsing
- REVEL score as dict or raw float
"""
```

Run for each mode combination (raw, normalized, decomposed, norm+decomp).

### 4c. Float Precision

nirvana2vcf formats floats with `.6g` (6 significant figures). Verify the comparison script uses the same precision when comparing numeric values from the JSON source.

### Success Criteria

- **100% match** on all fields for all variants (this is a self-consistency check)
- Any mismatch indicates a nirvana2vcf conversion bug

---

## Phase 5 — Cross-Tool Comparison (Independent Validation)

This phase compares nirvana2vcf output against **independently produced annotations** from other open-source tools. The goal is not exact matching (different tools use different database versions and algorithms), but confirming that nirvana2vcf output is **reasonable and consistent** with the broader annotation ecosystem.

### 5a. Ensembl VEP (Variant Effect Predictor)

**What it validates:** CSQ field format, SO consequence terms, gene symbols, SIFT/PolyPhen predictions, HGVS nomenclature.

**Why VEP specifically:** Nirvana uses VEP's consequence terminology (Sequence Ontology terms) and has documented >99.9% concordance with VEP on consequence calls. nirvana2vcf outputs a VEP-compatible CSQ field. Comparing against native VEP output validates that nirvana2vcf's CSQ is correctly structured and interoperable.

```bash
# Install VEP
conda install -c bioconda ensembl-vep
# Or: git clone https://github.com/Ensembl/ensembl-vep.git && cd ensembl-vep && perl INSTALL.pl

# Install VEP cache (~15 GB for human GRCh38)
vep_install -a cf -s homo_sapiens -y GRCh38 -c ~/.vep

# Annotate the same input VCF that Nirvana annotated
vep --input_file HiSeq.10000.vcf.gz \
    --output_file vep_output.vcf \
    --vcf \
    --offline \
    --cache \
    --assembly GRCh38 \
    --sift b \
    --polyphen b \
    --hgvs \
    --symbol \
    --canonical \
    --biotype \
    --force_overwrite
```

**Fields to compare (nirvana2vcf CSQ vs VEP CSQ):**

| Field | nirvana2vcf CSQ Subfield | VEP CSQ Subfield | Comparison Method |
|---|---|---|---|
| Consequence | `Consequence` | `Consequence` | SO term exact match (both use same ontology) |
| Gene symbol | `SYMBOL` | `SYMBOL` | Exact string match |
| Transcript ID | `Feature` | `Feature` | Exact match (both use RefSeq/Ensembl IDs) |
| Biotype | `BIOTYPE` | `BIOTYPE` | Exact match |
| HGVS coding | `HGVSc` | `HGVSc` | Exact OR transcript-strand reverse-complement at same `n.` coord |
| HGVS protein | `HGVSp` | `HGVSp` | Normalised: extract `p.<change>`, decode `%3D`, unwrap composite parens |
| SIFT | `SIFT` | `SIFT` | Categorical: deleterious/tolerated |
| PolyPhen | `PolyPhen` | `PolyPhen` | Categorical: probably_damaging/possibly_damaging/benign (Nirvana uses HVAR) |
| Canonical flag | `CANONICAL` | `MANE_SELECT` | Nirvana `CANONICAL=YES` vs VEP `MANE_SELECT` non-empty (cross-tool preferred-transcript axis) |
| Gene symbol (aliases) | `SYMBOL` | `SYMBOL` | Resolved through HGNC `prev_symbol`/`alias_symbol` map before comparing |
| Exon/Intron | `EXON`, `INTRON` | `EXON`, `INTRON` | Exact match (number/total format) |

**Comparison script approach:**

```python
#!/usr/bin/env python3
"""compare_vep.py — Compare nirvana2vcf CSQ against native VEP CSQ."""
import gzip, sys

def parse_csq(info, csq_fields):
    """Parse CSQ from INFO field into list of dicts."""
    for item in info.split(';'):
        if item.startswith('CSQ='):
            entries = item[4:].split(',')
            return [dict(zip(csq_fields, e.split('|'))) for e in entries]
    return []

def compare_csq(nirvana2vcf_csq, vep_csq):
    """Match transcripts by Feature ID, compare consequence and gene."""
    j2v_by_tx = {c['Feature']: c for c in nirvana2vcf_csq}
    vep_by_tx = {c['Feature']: c for c in vep_csq}
    shared = set(j2v_by_tx) & set(vep_by_tx)
    results = {'matched': 0, 'consequence_agree': 0, 'symbol_agree': 0}
    for tx in shared:
        results['matched'] += 1
        if j2v_by_tx[tx]['Consequence'] == vep_by_tx[tx]['Consequence']:
            results['consequence_agree'] += 1
        if j2v_by_tx[tx]['SYMBOL'] == vep_by_tx[tx]['SYMBOL']:
            results['symbol_agree'] += 1
    return results
```

**Expected concordance:**
- Consequence terms: >99% (Nirvana uses same SO terms as VEP)
- Gene symbols: >99% (both use HGNC)
- SIFT/PolyPhen predictions: >95% (algorithm versions may differ slightly)
- HGVS nomenclature: >95% (transcript version differences may cause minor divergence)

**Known differences to expect:**
- Transcript sets: VEP and Nirvana may include different transcript versions or different sets of non-canonical transcripts
- SIFT/PolyPhen scores: Nirvana bundles specific versions; VEP may use different versions. Nirvana's PolyPhen uses HVAR (HumanVar), not HDIV — ensure VEP is configured to match
- Consequence ranking: When a variant has multiple consequences, the "most severe" pick may differ

---

### 5b. SnpEff + SnpSift (Functional Annotation + Database Overlay)

**What it validates:** Functional consequence calls (independent of VEP), ClinVar annotations, dbSNP IDs.

**Why SnpEff:** It uses a completely different annotation engine than both Nirvana and VEP. Agreement across all three tools provides strong evidence that nirvana2vcf is not introducing systematic errors.

```bash
# Install
conda install -c bioconda snpeff snpsift

# Download SnpEff database
snpEff download -v GRCh38.105

# Step 1: Annotate with SnpEff (functional consequences via ANN field)
snpEff -v GRCh38.105 HiSeq.10000.vcf.gz > snpeff_annotated.vcf

# Step 2: Add ClinVar annotations with SnpSift
# Download ClinVar VCF (public, updated monthly)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

java -jar SnpSift.jar annotate \
  -name CLINVAR_ \
  -info CLNSIG,CLNDN,CLNREVSTAT \
  clinvar.vcf.gz \
  snpeff_annotated.vcf > snpeff_clinvar.vcf

# Step 3: Add dbSNP IDs
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi

java -jar SnpSift.jar annotate \
  -id \
  GCF_000001405.40.gz \
  snpeff_clinvar.vcf > snpeff_full.vcf
```

**Fields to compare (nirvana2vcf vs SnpEff/SnpSift):**

| Field | nirvana2vcf Output | SnpEff/SnpSift Output | Comparison Method |
|---|---|---|---|
| Consequence | CSQ `Consequence` (SO terms) | ANN field, subfield 2 (SO terms) | Categorical — map to shared set (see mapping table below) |
| Impact | Not directly output (derivable) | ANN field, subfield 3 (HIGH/MODERATE/LOW/MODIFIER) | Categorical |
| Gene symbol | CSQ `SYMBOL` | ANN field, subfield 4 | Exact match |
| ClinVar significance | `CLINVAR_SIG` | `CLINVAR_CLNSIG` | Categorical (normalize case/underscores) |
| ClinVar review status | `CLINVAR_REVSTAT` | `CLINVAR_CLNREVSTAT` | Categorical |
| dbSNP rsID | `ID` column | `ID` column | Exact match |

**SO consequence term mapping (nirvana2vcf/Nirvana vs SnpEff):**

Both tools use Sequence Ontology terms, but may use slightly different granularity:

| Concept | Nirvana/nirvana2vcf Term | SnpEff Term | Match? |
|---|---|---|---|
| Missense | `missense_variant` | `missense_variant` | Exact |
| Synonymous | `synonymous_variant` | `synonymous_variant` | Exact |
| Frameshift | `frameshift_variant` | `frameshift_variant` | Exact |
| Stop gained | `stop_gained` | `stop_gained` | Exact |
| Splice donor | `splice_donor_variant` | `splice_donor_variant` | Exact |
| Splice acceptor | `splice_acceptor_variant` | `splice_acceptor_variant` | Exact |
| Intron | `intron_variant` | `intron_variant` | Exact |
| Upstream | `upstream_gene_variant` | `upstream_gene_variant` | Exact |
| Downstream | `downstream_gene_variant` | `downstream_gene_variant` | Exact |
| UTR (5') | `5_prime_UTR_variant` | `5_prime_UTR_variant` | Exact |
| UTR (3') | `3_prime_UTR_variant` | `3_prime_UTR_variant` | Exact |
| Intergenic | `intergenic_variant` | `intergenic_region` | **Differs** — normalize before comparing |
| Non-coding tx | `non_coding_transcript_exon_variant` | `non_coding_transcript_exon_variant` | Exact |

**Expected concordance:**
- Gene symbol: >99%
- Consequence (after normalization): >95%
- ClinVar significance: >90% (depends on ClinVar version date — Nirvana's bundled version vs your downloaded version)
- dbSNP rsID: >99% (version-dependent)

**Known differences to expect:**
- SnpEff may report `intergenic_region` where Nirvana reports `intergenic_variant`
- ClinVar entries may differ if database versions are months apart (ClinVar updates monthly)
- SnpEff annotates against a different transcript set (RefSeq vs Ensembl differences)

---

### 5c. bcftools annotate (Database-Level Annotation)

**What it validates:** Population frequency values (gnomAD AF), dbSNP IDs — directly from authoritative source VCFs, independent of any annotation engine.

**Why bcftools:** This bypasses all annotation engines entirely. You annotate the VCF directly from the gnomAD/dbSNP source files using bcftools. If nirvana2vcf's gnomAD AF values match the values bcftools pulls from the same gnomAD release, it proves the full chain (gnomAD → Nirvana data bundle → Nirvana JSON → nirvana2vcf VCF) is faithful.

```bash
# Download gnomAD sites VCF for one chromosome (GRCh38)
# Full genome is ~750 GB; single chromosome is manageable
# chr21 is ~3 GB compressed
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr21.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr21.vcf.bgz.tbi

# Annotate your nirvana2vcf output with gnomAD AF directly
bcftools annotate \
  -a gnomad.genomes.v4.1.sites.chr21.vcf.bgz \
  -c INFO/gnomAD_DIRECT_AF:=INFO/AF \
  output.vcf \
  -o output_with_gnomad.vcf

# Now compare nirvana2vcf's gnomAD_AF against the directly-annotated gnomAD_DIRECT_AF
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/gnomAD_AF\t%INFO/gnomAD_DIRECT_AF\n' \
  output_with_gnomad.vcf | \
  awk -F'\t' '$5!="." && $6!="." {diff=$5-$6; if(diff<0) diff=-diff; print $0"\t"diff}' | \
  sort -t$'\t' -k7 -rn | head -20
```

**Alternatively, compare dbSNP IDs:**

```bash
# Download dbSNP VCF
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi

# Annotate nirvana2vcf output with dbSNP IDs
bcftools annotate \
  -a GCF_000001405.40.gz \
  -c ID \
  --set-id '%ID' \
  output_noid.vcf \
  -o output_with_dbsnp.vcf

# Compare ID columns
bcftools query -f '%CHROM\t%POS\t%ID\n' output.vcf > nirvana2vcf_ids.tsv
bcftools query -f '%CHROM\t%POS\t%ID\n' output_with_dbsnp.vcf > dbsnp_ids.tsv
diff nirvana2vcf_ids.tsv dbsnp_ids.tsv | head -20
```

**Expected concordance:**
- gnomAD AF: >99% within tolerance of 0.001 (if Nirvana bundles the same gnomAD version); systematically offset if versions differ
- dbSNP rsID: >99%

**Known differences to expect:**
- gnomAD version mismatch: Nirvana bundles a specific gnomAD version (check `Nirvana.dll --version` output). If you download gnomAD v4.1 but Nirvana bundles v3.1, expect population AFs to differ for many variants
- Multi-allelic representation: gnomAD VCFs may represent multi-allelics differently than the nirvana2vcf output; normalize both before comparing

---

### 5d. Cross-Tool Concordance Matrix

After running all three comparisons, produce a summary matrix:

```
=== CROSS-TOOL CONCORDANCE MATRIX ===

Input: HiSeq.10000.vcf.gz (10,000 variants)
Assembly: GRCh38

CONSEQUENCE CALLS (per-transcript)
                    nirvana2vcf    VEP         SnpEff
  nirvana2vcf          —           99.2%       96.8%
  VEP               99.2%       —           97.1%
  SnpEff             96.8%      97.1%       —

GENE SYMBOLS
                    nirvana2vcf    VEP         SnpEff
  nirvana2vcf          —           99.8%       99.7%
  VEP               99.8%       —           99.9%
  SnpEff             99.7%      99.9%       —

CLINVAR SIGNIFICANCE
                    nirvana2vcf    SnpSift
  nirvana2vcf          —           93.5%*
  SnpSift            93.5%      —
  * Discordance due to ClinVar version difference (2024-03 vs 2024-09)

gnomAD AF (correlation)
                    nirvana2vcf    bcftools-direct
  nirvana2vcf          —           r=0.9999
  bcftools-direct    r=0.9999   —

dbSNP rsID
                    nirvana2vcf    bcftools-direct
  nirvana2vcf          —           99.6%
  bcftools-direct    99.6%      —
```

---

## Disk Space & Time Estimates

| Step | Disk Space | Time |
|---|---|---|
| Nirvana install + GRCh38 data | ~35 GB | 30–60 min download |
| Annotate HiSeq.10000.vcf.gz | ~10 MB output | <1 min |
| Annotate GIAB NA12878 (full) | ~2 GB output | ~7 min |
| VEP cache (GRCh38) | ~15 GB | 20 min download |
| VEP annotate 10K variants | ~50 MB | ~2 min |
| SnpEff database (GRCh38.105) | ~1 GB | 5 min download |
| gnomAD chr21 sites VCF | ~3 GB | 10 min download |
| dbSNP VCF | ~15 GB | 30 min download |

**Minimum viable test** (HiSeq.10000 only): ~35 GB disk, ~1 hour setup, <5 min runtime.

---

## CI-Friendly Subset

For automated testing in CI (GitHub Actions, etc.), ship a small self-contained test:

1. **Include** a tiny Nirvana JSON file in `tests/data/` — extract the first 50–100 positions from `HiSeq.10000.json.gz` (this is public Illumina data, freely redistributable)
2. **Run** `nirvana2vcf` on it in all mode combinations (raw, normalized, decomposed, norm+decomp)
3. **Assert** known values: hard-code expected CHROM/POS/REF/ALT/INFO for 5–10 selected variants
4. **Validate structure**: use a pure-Python VCF syntax checker (no external tools needed)
5. **Optionally** install `vcf-validator` via conda in CI for spec-level validation

This keeps CI fast (<30 seconds) while the full validation described above can be run manually before releases.

---

## Summary

| Validation Layer | Independence Level | What It Catches |
|---|---|---|
| Phase 2: VCF validators | Tool-independent spec check | Malformed VCF, missing headers, wrong field types |
| Phase 3: Round-trip positions | Self-consistency | Dropped/duplicated variants, wrong POS/REF/ALT |
| Phase 4: Round-trip annotations | Self-consistency | Mangled annotation values, wrong field mapping |
| Phase 5a: VEP comparison | Independent engine, same ontology | Wrong consequence calls, bad CSQ format |
| Phase 5b: SnpEff comparison | Fully independent engine | Systematic annotation errors |
| Phase 5c: bcftools + source DBs | No annotation engine at all | Wrong population frequencies, wrong rsIDs |

Each layer adds confidence. Phases 1–4 require only Nirvana (open-source). Phase 5 adds VEP, SnpEff, bcftools, and public database downloads for maximum independence.

---

## Sources

- [Nirvana GitHub Repository](https://github.com/Illumina/Nirvana)
- [Nirvana Getting Started & HiSeq Test Data](https://illumina.github.io/NirvanaDocumentation/introduction/getting-started/)
- [Nirvana JSON Format Documentation](https://illumina.github.io/NirvanaDocumentation/file-formats/nirvana-json-file-format/)
- [EBI VCF Validator](https://github.com/EBIvariation/vcf-validator)
- [VCFtools](https://vcftools.github.io/perl_module.html)
- [bcftools](https://samtools.github.io/bcftools/bcftools.html)
- [Ensembl VEP Documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html)
- [SnpEff & SnpSift Documentation](https://pcingola.github.io/SnpEff/)
- [Genome in a Bottle (NIST)](https://www.nist.gov/programs-projects/genome-bottle)
- [1000 Genomes Data Portal](https://www.internationalgenome.org/data-portal/sample/NA12878)
- [Illumina Platinum Genomes](https://github.com/Illumina/PlatinumGenomes)
- [gnomAD Downloads](https://gnomad.broadinstitute.org/downloads)
- [ClinVar VCF Downloads](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/)
- [dbSNP Downloads](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/)
