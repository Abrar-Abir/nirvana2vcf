"""VCF 4.2 plain-text output writer."""

from typing import Any, Dict, IO, List, Optional

from .constants import (
    CSQ_FIELDS,
    FORMAT_DEFINITIONS,
    GRCH37_CONTIGS,
    GRCH38_CONTIGS,
    INFO_DEFINITIONS,
)
from .models import NirvanaHeader


def get_contig_header_lines(genome_assembly: str) -> List[str]:
    """Return ##contig=<...> lines for the given assembly."""
    if genome_assembly in ("GRCh38", "hg38"):
        contigs = GRCH38_CONTIGS
    elif genome_assembly in ("GRCh37", "hg19"):
        contigs = GRCH37_CONTIGS
    else:
        contigs = GRCH38_CONTIGS  # default

    lines = []
    for name, length in contigs.items():
        lines.append(f"##contig=<ID={name},length={length}>")
    return lines


def write_vcf_header(
    out: IO[str],
    header: NirvanaHeader,
    assembly: str,
    sample_names: List[str],
    csq_only: bool = False,
) -> None:
    """Write all VCF header lines.

    Writes ##fileformat, ##source, ##INFO, ##FORMAT, ##contig,
    and the #CHROM column header line.
    """
    # File format
    out.write("##fileformat=VCFv4.2\n")

    # Source
    out.write(f"##source=nirvana2vcf (from {header.annotator})\n")

    # INFO definitions
    if csq_only:
        # Only CSQ definition
        for info_def in INFO_DEFINITIONS:
            if "ID=CSQ," in info_def:
                out.write(info_def + "\n")
                break
    else:
        for info_def in INFO_DEFINITIONS:
            out.write(info_def + "\n")

    # FORMAT definitions (only if samples present)
    if sample_names:
        for fmt_def in FORMAT_DEFINITIONS:
            out.write(fmt_def + "\n")

    # Contig definitions
    for contig_line in get_contig_header_lines(assembly):
        out.write(contig_line + "\n")

    # Column header line
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if sample_names:
        columns.append("FORMAT")
        columns.extend(sample_names)
    out.write("\t".join(columns) + "\n")


def write_vcf_record(
    out: IO[str],
    record: Dict[str, Any],
    include_samples: bool = True,
) -> None:
    """Write a single VCF data line."""
    fields = [
        str(record["CHROM"]),
        str(record["POS"]),
        str(record["ID"]),
        str(record["REF"]),
        str(record["ALT"]),
        str(record["QUAL"]),
        str(record["FILTER"]),
        str(record["INFO"]),
    ]

    if include_samples and record.get("FORMAT"):
        fields.append(record["FORMAT"])
        fields.extend(record.get("samples", []))

    out.write("\t".join(fields) + "\n")
