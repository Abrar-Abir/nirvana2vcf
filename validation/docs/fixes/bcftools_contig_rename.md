# bcftools Contig Rename Fix (Phase 5)

## Problem

When `run_phase5.sh` prepares bare-chrom VCFs (renaming `chr1` to `1` for compatibility with GRCh37 resources like gnomAD and SnpEff), the command fails immediately:

```
[W::vcf_parse] Contig 'chr1' is not defined in the header. (Quick workaround: index the file with tabix.)
Encountered an error, cannot proceed. Please check the error output above.
If feeling adventurous, use the --force option. (At your own risk!)
```

This blocks all three Phase 5 sub-phases (VEP, SnpEff, bcftools) since the bare-chrom original VCF is a shared prerequisite.

## Root Cause

The original plan used `bcftools annotate --rename-chrs` to convert `chr1` to `1`:

```bash
bcftools annotate --rename-chrs chr_to_bare.txt original.vcf.gz -Oz -o bare.vcf.gz
```

This command requires contig names to be declared in `##contig` header lines. The HiSeq.10000 original VCF (VCFv4.1, produced by GATK UnifiedGenotyper) has **zero `##contig` lines** — it defines `##INFO`, `##FORMAT`, and `##FILTER` headers but never declares contigs. bcftools treats undeclared contigs in data lines as an error and refuses to proceed.

The nirvana2vcf output VCF does have `##contig` lines (25 of them, added by nirvana2vcf from Nirvana's reference metadata), so bcftools would succeed on that file — but the original VCF fails first and halts the script under `set -euo pipefail`.

## Fix Applied

Replaced both `bcftools annotate --rename-chrs` calls with a `strip_chr_prefix()` shell function that uses awk for text-based renaming:

```bash
strip_chr_prefix() {
    local infile="$1" outfile="$2"
    if [[ "$infile" == *.gz || "$infile" == *.bgz ]]; then
        gzip -dc "$infile"
    else
        cat "$infile"
    fi | awk '
        /^##contig=<ID=chr/ { sub(/ID=chr/, "ID="); print; next }
        /^#/ { print; next }
        { sub(/^chr/, "", $1); print }
    ' OFS="\t" | bgzip -c > "$outfile"
    tabix -p vcf "$outfile"
}
```

The function:

1. **Auto-detects compressed vs uncompressed input** — checks file extension for `.gz`/`.bgz` and uses `gzip -dc` or `cat` accordingly. This matters because the original VCF is bgzipped (`.vcf.gz`) but the nirvana2vcf output is plain text (`.vcf`).

2. **Renames `##contig` headers** — lines matching `##contig=<ID=chr...` have the `chr` prefix stripped from the contig ID. This keeps the header consistent with the data lines for downstream tools.

3. **Renames data-line CHROM column** — non-header lines have `chr` stripped from the first field. Uses `sub(/^chr/, "", $1)` which only affects a leading `chr` prefix.

4. **Produces bgzipped + tabix-indexed output** — pipes through `bgzip -c` and runs `tabix -p vcf` so the result is ready for bcftools annotation and SnpEff.

The now-unused `chr_to_bare.txt` mapping file and its `CHR_RENAME` variable were removed from the script.

### Why not `--force`?

bcftools suggests `--force` in its error message, but this flag suppresses all parsing errors indiscriminately. Using it risks silently corrupting output if there are other header issues. The awk approach is safe and explicit.

## Affected Locations

Two call sites in `run_phase5.sh` were updated:

| Section | Input file | Purpose |
|---------|-----------|---------|
| Section 1 (Prerequisites) | `HiSeq.10000.vcf.gz` (original, compressed, no `##contig`) | Bare-chrom VCF for SnpEff input |
| Section 4 (bcftools) | `HiSeq.10000.default.vcf` (nirvana2vcf output, uncompressed, has `##contig`) | Bare-chrom VCF for gnomAD AF comparison |

## Verification

Tested both input types:

```
# Compressed input (original VCF, no ##contig headers)
$ strip_chr_prefix HiSeq.10000.vcf.gz /tmp/test.vcf.gz
$ zcat /tmp/test.vcf.gz | grep -v '^#' | head -1 | cut -f1
1

# Uncompressed input (nirvana2vcf VCF, has ##contig headers)
$ strip_chr_prefix HiSeq.10000.default.vcf /tmp/test.vcf.gz
$ zcat /tmp/test.vcf.gz | grep '##contig' | head -1
##contig=<ID=1,length=249250621>
$ zcat /tmp/test.vcf.gz | grep -v '^#' | head -1 | cut -f1
1
```
