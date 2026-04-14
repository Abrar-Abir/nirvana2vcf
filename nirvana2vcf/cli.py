"""CLI entry point for nirvana2vcf."""

import argparse
import os
import sys

from . import __version__
from .mapper import decompose_position, map_position_to_vcf_record
from .parser import stream_positions
from .vcf_writer import write_vcf_header, write_vcf_record


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="nirvana2vcf",
        description="Convert Nirvana/Illumina Connected Annotations JSON to VCF format.",
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input Nirvana JSON file (.json or .json.gz)",
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output VCF file (default: stdout)",
    )
    parser.add_argument(
        "--csq-only", action="store_true",
        help="Output only VEP-style CSQ in INFO field (no flat annotation fields)",
    )
    parser.add_argument(
        "--no-samples", action="store_true",
        help="Omit sample/genotype columns from output",
    )
    parser.add_argument(
        "--assembly", choices=["GRCh37", "GRCh38", "auto"], default="auto",
        help="Genome assembly for contig headers (default: auto from header)",
    )
    parser.add_argument(
        "--normalize", action="store_true", default=True, dest="normalize",
        help="Normalize allele representations to minimal VCF convention (default: enabled)",
    )
    parser.add_argument(
        "--no-normalize", action="store_false", dest="normalize",
        help="Disable allele normalization (output raw Nirvana alleles)",
    )
    parser.add_argument(
        "--decompose", action="store_true", default=False, dest="decompose",
        help="Decompose multi-allelic sites into biallelic rows (like bcftools norm -m-)",
    )
    parser.add_argument(
        "--no-decompose", action="store_false", dest="decompose",
        help="Keep multi-allelic sites as single rows (default)",
    )
    parser.add_argument(
        "--version", action="version",
        version=f"%(prog)s {__version__}",
    )

    args = parser.parse_args(argv)

    # Validate input exists
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Open output
    out = None
    try:
        if args.output:
            out = open(args.output, "w", encoding="utf-8")
        else:
            out = sys.stdout

        include_samples = not args.no_samples
        header_written = False
        for nirvana_header, position in stream_positions(args.input):
            if not header_written:
                assembly = (
                    nirvana_header.genome_assembly
                    if args.assembly == "auto"
                    else args.assembly
                )
                sample_names = (
                    nirvana_header.samples if include_samples else []
                )
                write_vcf_header(
                    out, nirvana_header, assembly, sample_names,
                    csq_only=args.csq_only,
                )
                header_written = True

            positions_to_map = (
                decompose_position(position) if args.decompose
                else [position]
            )
            for pos in positions_to_map:
                record = map_position_to_vcf_record(
                    pos, nirvana_header,
                    csq_only=args.csq_only,
                    include_samples=include_samples,
                    normalize=args.normalize,
                )
                write_vcf_record(
                    out, record,
                    include_samples=include_samples,
                )
    finally:
        if out is not None and out is not sys.stdout:
            out.close()
