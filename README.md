# nirvana2vcf

Convert [Nirvana/Illumina Connected Annotations](https://illumina.github.io/NirvanaDocumentation/) JSON output to VCF 4.2 format.

## Features

- Pure Python, zero external dependencies at runtime
- Streaming pipeline — processes one position at a time, no full-file load into memory
- Reads `.json` and `.json.gz` input
- Supports GRCh37 and GRCh38 assemblies (auto-detected from Nirvana header)
- Allele normalization — trims shared prefix/suffix to minimal VCF representation (enabled by default)
- Multi-allelic decomposition — splits multi-allelic sites into biallelic rows, like `bcftools norm -m-` (`--decompose`)
- VEP-style CSQ field with per-transcript annotations
- Annotations: gnomAD, ClinVar, SpliceAI, REVEL, DANN, GERP, phyloP, 1000 Genomes, TOPMed

## Installation

```bash
pip install -e .
```

## Usage

```bash
# Basic conversion
nirvana2vcf -i input.json.gz -o output.vcf

# Output to stdout (pipe to bcftools, etc.)
nirvana2vcf -i input.json.gz | bcftools view -f PASS

# VEP-style CSQ only (no flat INFO fields)
nirvana2vcf -i input.json -o output.vcf --csq-only

# Omit sample/genotype columns
nirvana2vcf -i input.json.gz -o output.vcf --no-samples

# Override genome assembly (instead of auto-detecting from header)
nirvana2vcf -i input.json.gz -o output.vcf --assembly GRCh37

# Disable allele normalization (keep raw Nirvana alleles)
nirvana2vcf -i input.json.gz -o output.vcf --no-normalize

# Decompose multi-allelic sites into biallelic rows (normalization is applied by default)
nirvana2vcf -i input.json.gz -o output.vcf --decompose

# Decompose without normalization (keep raw Nirvana alleles)
nirvana2vcf -i input.json.gz -o output.vcf --decompose --no-normalize
```

## Mode Reference

### Allele Normalization (`--normalize` / `--no-normalize`)

Enabled by default. Trims shared prefix and suffix bases from REF and ALT alleles to produce the
minimal VCF representation, adjusting POS accordingly. This matches what tools like `bcftools norm`
and `vt normalize` do for left-trimming (though without reference-based left-alignment).

Nirvana sometimes emits redundant flanking bases — for example, when representing an insertion
or deletion relative to a longer context sequence.

**Example — SNV emitted with flanking context bases:**

```
# Raw Nirvana output (--no-normalize):
#   Nirvana emitted C→T substitution with flanking A…GT context:
chr1  1000  .  ACGT  ATGT  .  .  ...

# Normalized (default):
#   Phase 1 (right-trim): strip shared T → ACGT→ACG, ATGT→ATG
#                          strip shared G → ACG→AC,  ATG→AT
#                          AC[-1]=C ≠ AT[-1]=T → stop
#   Phase 2 (left-trim):  strip shared A → AC→C, AT→T, POS advances to 1001
#                          len(C)=1 → stop
chr1  1001  .  C  T  .  .  ...
```

**Example — deletion with right-anchor padding:**

```
# Raw Nirvana output (--no-normalize):
#   Deletion of A, represented with trailing GT context:
chr1  1000  .  ACGT  CGT  .  .  ...

# Normalized (default):
#   Phase 1 (right-trim): strip T → ACGT→ACG, CGT→CG
#                          strip G → ACG→AC,  CG→C
#                          len(C)=1 → stop
#   Phase 2 (left-trim):  len(C)=1 → stop (REF=AC, ALT=C is already minimal)
chr1  1000  .  AC  C  .  .  ...
```

**Example — SNV that needs no trimming:**

```
# Both modes produce the same output for a clean SNV:
chr7  117548628  .  A  G  .  .  ...
```

Symbolic alleles (`<DEL>`, `<DUP>`, etc.), reference-only ALTs (`.`), and spanning deletions (`*`)
are never modified by normalization.

---

### Multi-allelic Decomposition (`--decompose`)

Disabled by default. When enabled, a position with multiple ALT alleles is split into one VCF row
per ALT allele (like `bcftools norm -m-`). Variant annotations and per-allele INFO fields are scoped
to each row. Sample genotypes, allele depths (AD), and variant frequencies (VF) are remapped:

- GT: alleles matching this ALT → `1`; other ALTs → `.` (missing); REF stays `0`
- AD: `[ref_depth, this_alt_depth]`
- VF: `[this_alt_frequency]`

**Example — tri-allelic site:**

```
# Without --decompose (one row, multi-allelic):
chr1  925952  .  GCACA  ACACA,G  .  .  gnomAD_AF=0.001,0.0005  GT:AD  1/2:10,5,3

# With --decompose (two rows, biallelic):
chr1  925952  .  GCACA  ACACA  .  .  gnomAD_AF=0.001  GT:AD  1/.:10,5
chr1  925952  .  GCACA  G      .  .  gnomAD_AF=0.0005 GT:AD  ./1:10,3
```

When combined with `--normalize` (the default), decomposition runs first, then each biallelic row
is normalized independently. This means rows may end up with different POS values if their alleles
trim differently.

---

### VEP-style CSQ (`--csq-only`)

By default, nirvana2vcf writes flat INFO fields for each annotation type (`gnomAD_AF`, `ClinVar_SIG`,
`SpliceAI_DS_AG`, etc.). With `--csq-only`, all transcript-level annotations are packed into a
single `CSQ` INFO field in the same pipe-delimited format used by Ensembl VEP, and the flat fields
are omitted.

Use `--csq-only` when your downstream tool (e.g. a variant database loader) expects VEP-annotated
VCFs and parses the CSQ field directly.

```
# Default flat fields:
INFO=gnomAD_AF=0.0032;ClinVar_SIG=Pathogenic;SpliceAI_DS_AG=0.85;REVEL=0.92

# With --csq-only:
INFO=CSQ=ENST00000357654|missense_variant|MODERATE|BRCA1|...|0.92|...
```

The CSQ format string is written to the VCF header (`##INFO=<ID=CSQ,...,Format="Allele|...">`).

---

### Omitting Samples (`--no-samples`)

By default, all samples from the Nirvana JSON are written as genotype columns in the VCF. Use
`--no-samples` to produce a sites-only VCF with no FORMAT or sample columns — useful for annotation
databases or tools that do not expect genotype data.

```
# Default (with samples):
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO                FORMAT  SAMPLE1  SAMPLE2
chr1    925952  .  A    G    .     .       gnomAD_AF=0.001     GT:DP   0/1:30   0/0:25

# With --no-samples:
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
chr1    925952  .  A    G    .     .       gnomAD_AF=0.001
```

---

### Assembly Override (`--assembly`)

Assembly is auto-detected from the Nirvana JSON header (the `genomeAssembly` field). Use
`--assembly GRCh37` or `--assembly GRCh38` to override this — for example, if the header is
missing or wrong.

The assembly affects contig `##contig` header lines and chromosome naming:
- GRCh38: `chr1`, `chr2`, ..., `chrX`, `chrY`, `chrM`
- GRCh37: `1`, `2`, ..., `X`, `Y`, `MT` (no `chr` prefix)

## Development

```bash
pip install -e ".[dev]"
python3 -m pytest -v         # run all tests
python3 -m pytest -v -k csq  # run tests matching keyword
```

## Architecture

Streaming pipeline: **parse** → **map** → **write**

- `nirvana2vcf/parser.py` — Streams Nirvana's line-based JSON format, yielding `(NirvanaHeader, Position)` tuples
- `nirvana2vcf/mapper.py` — Transforms positions into VCF record dicts (per-allele fields, CSQ, INFO escaping, normalization, decomposition)
- `nirvana2vcf/vcf_writer.py` — Writes VCF 4.2 plain text
- `nirvana2vcf/models.py` — Dataclass contracts between parser and mapper
- `nirvana2vcf/constants.py` — VCF header definitions, contig maps, CSQ field names

### Nirvana JSON format

Nirvana emits a line-based streaming format (not standard JSON):
- Line 1: `{"header":{...},"positions":[`
- Lines 2–N: one position object per line (comma-separated)
- Last line: `],"genes":[...]}`

Each position line maps directly to one VCF row (or multiple rows when using `--decompose`).

## Notes on field provenance

### `gnomAD_EUR_AF` is gnomAD European non-Finnish (NFE)

The `gnomAD_EUR_AF` INFO field is sourced from Nirvana's `gnomad.nfeAf`
(European non-Finnish) — not `gnomad.eurAf`, which is a 1000 Genomes field
that never appears inside Nirvana's `gnomad` block. The INFO header
description reflects this (`"gnomAD allele frequency (European non-Finnish)"`).

The name `gnomAD_EUR_AF` (rather than `gnomAD_NFE_AF`) is kept deliberately
for stability with downstream consumers that expect an `EUR` suffix. When
comparing against another VCF, map this field against the other source's
`NFE` column, not an `EUR` column.

### `gnomAD_*` fields reflect Nirvana's bundled gnomAD release

Nirvana's `gnomad` block combines gnomAD genomes + exomes from the
version bundled with your Nirvana data files (v2.1 for the release used
in [validation/docs/phases/phase05.md](validation/docs/phases/phase05.md)).
When cross-checking against a standalone gnomAD VCF, match both the
version and the genomes/exomes scope — otherwise concordance numbers
will be misleading. See [phase05.md](validation/docs/phases/phase05.md)
for a worked example.

### PolyPhen predictions are HVAR, not HDIV

Nirvana's `polyPhenPrediction` is trained on HumanVar (HVAR), not
HumanDiv (HDIV). When comparing `PolyPhen` (from the CSQ field) against
dbNSFP, use `dbNSFP_Polyphen2_HVAR_pred`, not `_HDIV_pred`. HDIV leans
more damaging and HVAR leans more benign, so mismatching the flavor
produces a strongly asymmetric `B` vs `D/P` disagreement pattern.

## Validation

Open-source validation using only public data and open-source tools
lives under [validation/](validation/). See
[validation/docs/strategy.md](validation/docs/strategy.md) for the
five-phase strategy and [validation/results/](validation/results/) for
concordance reports against `bcftools`, VEP, and SnpEff on Nirvana's
bundled 10,000-variant HiSeq test VCF.

## Contributing

Issues and pull requests are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md)
for development setup and guidelines.

## License

MIT — see [LICENSE](LICENSE).
