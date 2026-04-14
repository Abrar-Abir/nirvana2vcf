# gnomAD Annotation Fixes (Phase 5c)

## Problem 1: Malformed FORMAT Fields in gnomAD VCF

After fixing VCF preprocessing issues (see `vcf_preprocessing.md`), `bcftools annotate` failed reading the gnomAD v2.1.1 annotation source:

```
[W::vcf_parse_format_dict2] FORMAT '-' at 1:874678 is not defined in the header
[E::vcf_parse_format_dict2] Could not add dummy header for FORMAT
    'damaging(0.56)|tolerated(0.58),T|missense_variant|SAMD11|...' at 1:876680
Error: VCF parse error
```

Two distinct malformations in the gnomAD genomes VCF:

1. **Literal `-` in FORMAT column** (position 874678) — the genotype field separator is a bare hyphen instead of valid FORMAT tags.
2. **CSQ annotation string in FORMAT column** (position 876680) — an entire VEP CSQ entry (`damaging(0.56)|tolerated(0.58),T|missense_variant|...`) appears where FORMAT tags should be. This contains parentheses, pipes, and commas that make it impossible for bcftools to create even a dummy header entry.

`--force` bypasses the first error but not the second — bcftools cannot construct a valid header entry for the CSQ-contaminated FORMAT string and aborts.

### Root Cause

These are data quality issues in the gnomAD v2.1.1 genomes chr1 VCF (`gnomad.genomes.r2.1.1.sites.1.vcf.bgz`). The file is distributed by the gnomAD project and cannot be fixed upstream. The malformations are in FORMAT/genotype columns which are irrelevant to the annotation operation (we only need `INFO/AF`).

### Fix Applied

Added a sites-only extraction step (Step 2b in `run_phase5.sh`) that strips all FORMAT and sample columns using awk, bypassing bcftools' FORMAT parser entirely:

```bash
gzip -dc "$GNOMAD_SUBSET" | awk '
    /^##FORMAT/ { next }
    /^#CHROM/ {
        for (i = 1; i <= 8; i++) printf "%s%s", $i, (i < 8 ? "\t" : "\n")
        next
    }
    /^#/ { print; next }
    {
        for (i = 1; i <= 8; i++) printf "%s%s", $i, (i < 8 ? "\t" : "\n")
    }
' | bgzip -c > "$GNOMAD_SITES"
tabix -p vcf "$GNOMAD_SITES"
```

The resulting sites-only VCF contains only columns 1–8 (CHROM through INFO), with `##FORMAT` header lines removed. `bcftools annotate` then reads this clean file without any `--force` flag.

---

## Problem 2: 38 GB Streaming Merge Performance

The initial implementation annotated the nirvana2vcf VCF directly against the full gnomAD genomes chr1 file:

```bash
bcftools annotate -a gnomad.genomes.r2.1.1.sites.1.vcf.bgz \
    -c "INFO/GNOMAD_DIRECT_AF:=INFO/AF" j2v_bare.vcf.gz -o annotated.vcf
```

This command ran for 10+ minutes without completing and was killed by timeout.

### Root Cause

`bcftools annotate -a` performs a **streaming merge** — it reads both the input VCF and the annotation VCF in coordinate order, advancing through both files simultaneously. Even though only ~10,000 variants need annotation, bcftools must read through the dense gnomAD genomes file from the first input position to the last.

The gnomAD genomes chr1 VCF is **38 GB** (bgzipped) and contains sites at nearly every genomic position. The input variants span positions 109 to 5,235,136 on chr1. At the observed throughput (~87K positions/minute), annotating would take approximately 1 hour.

The tabix index exists but `bcftools annotate` uses streaming merge, not random access, for the `-a` annotation source.

### Fix Applied

Replaced the single-pass annotation with a three-step pipeline:

**Step 1: Create BED regions** from nirvana2vcf variant positions:
```bash
bcftools query -f '%CHROM\t%POS\n' j2v_bare.vcf.gz \
    | awk '{print $1"\t"$2-1"\t"$2}' > j2v_regions.bed
```

**Step 2: Extract matching gnomAD sites** using tabix random access:
```bash
bcftools view --force -R j2v_regions.bed gnomad.vcf.bgz -Oz -o gnomad_subset.vcf.gz
```

`bcftools view -R` uses the tabix index to jump directly to each region, extracting only the ~2,700 gnomAD sites that overlap our 9,965 variant positions. The `--force` flag is needed here to handle the FORMAT malformations during extraction (they're in the FORMAT column which bcftools reads but we later discard). Output: ~5 MB (vs 38 GB original). Runtime: ~5 minutes.

**Step 2b: Create sites-only subset** (strips malformed FORMAT columns — see Problem 1 above).

**Step 3: Annotate against the small clean subset**:
```bash
bcftools annotate -a gnomad_subset.sites.vcf.gz \
    -c "INFO/GNOMAD_DIRECT_AF:=INFO/AF" j2v_bare.vcf.gz -o annotated.vcf
```

Runtime: < 1 second.

Each intermediate file is cached independently with `[ ! -f "$FILE" ]` guards, so re-runs skip completed steps.

## Verification

After both fixes, the full pipeline completes successfully:

```
gnomAD subset: 2685 sites
bcftools annotation complete: j2v_gnomad_annotated.vcf
```

The comparison report shows:
- **gnomAD AF**: 149/149 concordant where both present (tolerance ≤ 0.001), Pearson r = 1.000000
- **dbSNP rsID**: 610/610 = 100.0% concordance

## Affected File

`validation/scripts/run_phase5.sh` — Section 4 (Sub-Phase 5c: bcftools Direct Annotation), Steps 1–3.
