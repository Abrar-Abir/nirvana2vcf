# VCF Preprocessing Fixes (Phase 5)

## Problem 1: Missing ##FILTER Header Declarations

After fixing the contig rename issue (see `bcftools_contig_rename.md`), `bcftools annotate` still fails when reading the nirvana2vcf output VCF:

```
[W::vcf_parse_filter] FILTER 'FDRtranche2.00to10.00+' is not defined in the header
Encountered an error, cannot proceed.
```

Even with `--force`, bcftools fails during output formatting:

```
[E::vcf_format] Invalid BCF, the FILTER tag id=51 at 1:250 not present in the header
[main_vcfannotate] Error: failed to write to j2v_gnomad_annotated.vcf
```

bcftools internally assigns numeric IDs to FILTER values from header declarations. Without `##FILTER` lines, those IDs are unresolvable during output, causing a fatal write error that `--force` cannot bypass.

### Root Cause

nirvana2vcf preserves FILTER values from the source VCF data lines (DPFilter, HARD_TO_VALIDATE, FDRtranche*, Indel, LowQual, SnpCluster) but does not emit `##FILTER=<ID=...,Description=...>` header lines. The output VCF only has `##FILTER=<ID=PASS,...>` from the VCF spec default.

The HiSeq.10000 dataset uses GATK UnifiedGenotyper VQSR tranches and has 9 distinct FILTER values in data lines, none declared in the header.

### Fix Applied

Added `inject_filter_headers()` shell function in `run_phase5.sh`. It does a single-pass scan that buffers the entire file in awk (acceptable for ~10K-line VCFs), collects unique FILTER values from column 7, and injects `##FILTER` declarations before the `#CHROM` line:

```bash
inject_filter_headers() {
    local infile="$1" outfile="$2"
    # ... decompression logic ...
    awk -F'\t' '
        /^##/ { headers[++nh] = $0; next }
        /^#CHROM/ { chrom_line = $0; next }
        {
            n = split($7, arr, ";")
            for (i = 1; i <= n; i++) {
                if (arr[i] != "." && arr[i] != "PASS" && !(arr[i] in seen))
                    seen[arr[i]] = 1
            }
            data[++nd] = $0
        }
        END {
            for (i = 1; i <= nh; i++) print headers[i]
            for (f in seen) print "##FILTER=<ID=" f ",Description=\"Imported from source VCF\">"
            print chrom_line
            for (i = 1; i <= nd; i++) print data[i]
        }
    ' | bgzip -c > "$outfile"
    tabix -p vcf "$outfile"
}
```

This runs after `strip_chr_prefix()` and before `bcftools annotate`, producing a VCF that bcftools can parse cleanly without `--force`.

---

## Problem 2: awk Field Separator Corrupts Lines with Spaces

After injecting FILTER headers, `bcftools annotate` hit a new error:

```
[W::vcf_parse_format_dict2] FORMAT '-' at 1:874678 is not defined in the header
```

The FORMAT column at position 874678 was `-` instead of `GT:DP:GQ:AD:VF`. Investigation showed the line had **13 tab-delimited fields instead of 10** — extra tabs had been injected.

### Root Cause

Both `strip_chr_prefix()` and `inject_filter_headers()` used awk with the **default field separator** (whitespace). VCF data lines are tab-delimited, but the CSQ INFO field can contain literal spaces in annotation values. At position 874678, the CSQ contained:

```
tolerated - low confidence
```

With the default FS, awk split this into three additional fields on the spaces. When `strip_chr_prefix()` modified `$1` (stripping `chr`), awk reconstructed `$0` using `OFS="\t"`, joining all fields (including the spurious space-split ones) with tabs. The resulting line had 13 fields, with `FORMAT` shifted to `-` and the sample data split across extra columns.

### Why Only Some Lines Were Affected

Most CSQ values (gene symbols, HGVS notation, SO terms) contain no spaces. Only a few variants had prediction strings like `"tolerated - low confidence"` or similar free-text annotations. The corruption was silent — no error until bcftools tried to parse the FORMAT column.

### Fix Applied

Added explicit `-F'\t'` to both awk invocations so awk splits only on tabs:

```bash
# Before (broken):
awk '...' OFS="\t"

# After (fixed):
awk -F'\t' '
    BEGIN { OFS="\t" }
    ...
'
```

### Verification

```bash
# Before fix: line at 874678 has 13 fields, FORMAT is '-'
$ gzip -dc j2v_bare.vcf.gz | awk -F'\t' '$2 == 874678 {print NF, $9}'
13 -

# After fix: line has 10 fields, FORMAT is correct
$ gzip -dc j2v_bare.vcf.gz | awk -F'\t' '$2 == 874678 {print NF, $9}'
10 GT:DP:GQ:AD:VF
```

## Affected Locations

| Function | File | Purpose |
|----------|------|---------|
| `strip_chr_prefix()` | `scripts/run_phase5.sh` | Renames `chr1`→`1` in contig headers and data lines |
| `inject_filter_headers()` | `scripts/run_phase5.sh` | Injects `##FILTER` declarations from data-line values |

Both functions are also used by `bcftools_contig_rename.md`'s `strip_chr_prefix()` — that fix doc predates the FS correction and should be read in conjunction with this one.
