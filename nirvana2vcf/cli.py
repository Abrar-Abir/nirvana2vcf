"""CLI entry point for nirvana2vcf."""

import argparse
import io
import os
import sys

from . import __version__
from ._optional import require_pysam
from .mapper import decompose_position, map_position_to_vcf_record
from .parser import stream_positions
from .vcf_writer import write_vcf_header, write_vcf_record


PROGRESS_INTERVAL = 10000


def _is_bgzip_path(path):
    return path.endswith(".vcf.gz") or path.endswith(".vcf.bgz")


def _open_output(path):
    """Return (text_writer, close_fn, is_bgzip). Plain text if path is None or .vcf."""
    if path is None:
        return sys.stdout, (lambda: None), False
    if _is_bgzip_path(path):
        pysam = require_pysam("bgzip output (.vcf.gz)")
        bgz = pysam.BGZFile(path, "wb")
        text = io.TextIOWrapper(bgz, encoding="utf-8", write_through=True)

        def _close():
            text.flush()
            text.close()
            bgz.close()

        return text, _close, True
    fh = open(path, "w", encoding="utf-8")
    return fh, fh.close, False


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
        help="Output VCF file (default: stdout). Ending with .vcf.gz or .vcf.bgz writes bgzipped output (requires pysam).",
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
        "--reference", "-r", default=None,
        help="Reference FASTA (must be indexed; .fai alongside). Enables bcftools-norm-equivalent left-shifting of indels through repeats. Requires pysam.",
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
        "--tabix", action="store_true",
        help="After writing, build a tabix index (.tbi). Requires bgzipped output and pysam.",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help=f"Print progress to stderr every {PROGRESS_INTERVAL} positions.",
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

    # Validate --tabix implies bgzip output
    if args.tabix and not (args.output and _is_bgzip_path(args.output)):
        print(
            "Error: --tabix requires bgzipped output (-o ending in .vcf.gz or .vcf.bgz)",
            file=sys.stderr,
        )
        sys.exit(1)

    # Optional reference for left-shift normalization
    ref_fetcher = None
    fasta = None
    if args.reference:
        if not os.path.exists(args.reference):
            print(f"Error: Reference FASTA not found: {args.reference}", file=sys.stderr)
            sys.exit(1)
        pysam = require_pysam("--reference left-alignment")
        fasta = pysam.FastaFile(args.reference)

        def ref_fetcher(chrom, start, end, _fa=fasta):
            try:
                return _fa.fetch(chrom, start, end)
            except (KeyError, ValueError):
                # Try with/without 'chr' prefix
                alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
                try:
                    return _fa.fetch(alt, start, end)
                except (KeyError, ValueError):
                    return ""

    # Open output
    out, close_out, _ = _open_output(args.output)

    include_samples = not args.no_samples
    header_written = False
    n_positions = 0
    n_rows = 0
    try:
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
                    ref_fetcher=ref_fetcher,
                )
                write_vcf_record(
                    out, record,
                    include_samples=include_samples,
                )
                n_rows += 1

            n_positions += 1
            if args.verbose and n_positions % PROGRESS_INTERVAL == 0:
                print(
                    f"[nirvana2vcf] {n_positions} positions processed",
                    file=sys.stderr,
                )
    finally:
        close_out()
        if fasta is not None:
            fasta.close()

    if args.verbose:
        print(
            f"[nirvana2vcf] done — {n_positions} positions, {n_rows} VCF rows",
            file=sys.stderr,
        )

    if args.tabix:
        pysam = require_pysam("--tabix indexing")
        pysam.tabix_index(args.output, preset="vcf", force=True)
