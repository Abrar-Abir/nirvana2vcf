# VEP & SnpEff Installation Fixes (Phase 5a/5b)

## Problem 1: VEP segfault on macOS ARM (NDBM_File load order)

When running VEP on macOS ARM (Apple Silicon) after installing via `conda install -c bioconda ensembl-vep`, every invocation crashes with a segfault:

```
$ vep --help
[1]    12345 segmentation fault  vep --help
$ echo $?
139
```

No output, no error message. Exit code 139 = SIGSEGV.

## Root Cause

The crash originates in Perl's `AnyDBM_File` module when loaded implicitly by BioPerl's `Bio::DB::IndexedBase`. The module chain is:

```
vep -> Bio::EnsEMBL::VEP::Runner
     -> Bio::EnsEMBL::VEP::BaseVEP
     -> Bio::EnsEMBL::Variation::Utils::FastaSequence
     -> Bio::DB::Fasta
     -> Bio::DB::IndexedBase
     -> AnyDBM_File
     -> NDBM_File (compiled .bundle)
     -> SEGV
```

`AnyDBM_File` dispatches to `NDBM_File` (the first available DBM backend). The `NDBM_File.bundle` compiled extension crashes during initialization when loaded as part of BioPerl's deep `use` chain. The same module loads fine when imported first, before BioPerl:

```perl
# Crashes:
$ perl -e 'use Bio::DB::IndexedBase'    # exit 139

# Works:
$ perl -e 'use AnyDBM_File; use Bio::DB::IndexedBase; print "OK\n"'
# OK
```

This appears to be a load-order-dependent initialization bug in the conda `NDBM_File` build for `osx-arm64`. Pre-loading `AnyDBM_File` (which triggers `NDBM_File` loading) before BioPerl's deep module graph avoids the crash.

### Why `sitecustomize.pl` doesn't help

The obvious fix — placing `use AnyDBM_File;` in `sitecustomize.pl` — fails because the conda Perl 5.32 build was compiled without `-Dusesitecustomize`:

```
$ perl -V | grep sitecustomize
# (no output)
```

## Fix Applied

Replace the `vep` symlink (`bin/vep -> ../share/ensembl-vep-115.2-1/vep`) with a wrapper script that pre-loads `AnyDBM_File`:

```bash
#!/bin/bash
exec perl -MAnyDBM_File /opt/homebrew/Caskroom/miniforge/base/envs/biotools/share/ensembl-vep-115.2-1/vep "$@"
```

The same wrapper is applied to `vep_install`. The original symlinks are preserved as `.symlink.bak` files.

### Affected packages

| Package | Version | Build | Channel |
|---|---|---|---|
| perl | 5.32.1 | 7_h4614cfb_perl5 | conda-forge |
| perl-bioperl-core | 1.7.8 | pl5321hdfd78af_1 | bioconda |
| ensembl-vep | 115.2 | pl5321h2a3209d_1 | bioconda |

---

## Problem 2: SnpEff OutOfMemoryError loading GRCh37.75

SnpEff crashes with `OutOfMemoryError` when loading the GRCh37.75 database:

```
java.lang.OutOfMemoryError: Java heap space
    at org.snpeff.interval.Intron.createSpliceSiteRegionStart(Intron.java:136)
    at org.snpeff.interval.Transcript.createSpliceSites(Transcript.java:733)
    at org.snpeff.interval.Genes.createSpliceSites(Genes.java:129)
    at org.snpeff.snpEffect.SnpEffectPredictor.createGenomicRegions(SnpEffectPredictor.java:185)
    at org.snpeff.snpEffect.SnpEffectPredictor.buildForest(SnpEffectPredictor.java:134)
```

## Root Cause

The conda `snpEff` wrapper script defaults to `-Xms512m -Xmx1g`. The GRCh37.75 human genome database requires more heap to build the interval forest for splice site regions.

## Fix Applied

Pass `-Xmx4g` to snpEff in `run_phase5.sh`:

```bash
# Before:
snpEff -v GRCh37.75 "$BARE_ORIGINAL"

# After:
snpEff -Xmx4g -v GRCh37.75 "$BARE_ORIGINAL"
```

The `-Xmx4g` flag is passed through the conda wrapper to the JVM.

---

## Problem 3: VEP requires FASTA for indexed cache

VEP fails with two sequential errors when using the indexed (BAM-edited) cache:

```
ERROR: Cannot generate HGVS coordinates (--hgvs and --hgvsg) in offline mode without a FASTA file
```

After removing `--hgvs`:

```
ERROR: Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file
```

## Root Cause

The indexed VEP cache (downloaded from `indexed_vep_cache/` on the Ensembl FTP) is "BAM-edited" — it stores transcript sequences relative to a reference. VEP auto-enables `--use_transcript_ref` for BAM-edited caches, which requires a FASTA file. Without a downloaded reference FASTA (~3 GB compressed), both `--hgvs` and transcript-ref mode fail.

## Fix Applied

Two changes to the VEP invocation in `run_phase5.sh`:

1. **Removed `--hgvs`** — HGVSc/HGVSp fields show as N/A in the comparison (nirvana2vcf-only). The key comparison fields (Consequence, SYMBOL, BIOTYPE, SIFT, PolyPhen, EXON, INTRON) work without FASTA.

2. **Added `--use_given_ref`** — Overrides the auto-enabled `--use_transcript_ref`, telling VEP to use the alleles as given in the input VCF rather than looking them up in the reference.

```bash
vep \
    --input_file "$ORIGINAL_VCF" \
    --output_file "$VEP_OUTPUT" \
    --vcf --offline --cache --assembly GRCh37 \
    --sift b --polyphen b \
    --symbol --canonical --biotype --numbers \
    --merged \
    --use_given_ref \       # override BAM-edited cache auto-ref
    --force_overwrite
```

---

## Problem 4: VEP cache path detection

`run_phase5.sh` checked for the cache at `~/.vep/homo_sapiens/` but the merged cache extracts to `~/.vep/homo_sapiens_merged/`.

## Fix Applied

Added a fallback check in `run_phase5.sh`:

```bash
VEP_CACHE_DIR="${HOME}/.vep/homo_sapiens"
if [ -d "${HOME}/.vep/homo_sapiens_merged" ] && [ ! -d "$VEP_CACHE_DIR" ]; then
    VEP_CACHE_DIR="${HOME}/.vep/homo_sapiens_merged"
fi
```

---

## Summary of Script Changes

All fixes are in `validation/scripts/run_phase5.sh`:

| Line | Change | Problem |
|---|---|---|
| VEP invocation | Removed `--hgvs`, added `--use_given_ref` | FASTA requirement |
| SnpEff invocation | Added `-Xmx4g` | Java heap OOM |
| Cache detection | Added `homo_sapiens_merged/` fallback | Merged cache path |

External to the script (conda environment):

| File | Change | Problem |
|---|---|---|
| `$CONDA_PREFIX/bin/vep` | Symlink replaced with wrapper | Perl NDBM_File segfault |
| `$CONDA_PREFIX/bin/vep_install` | Symlink replaced with wrapper | Same segfault |
