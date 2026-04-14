"""End-to-end integration tests: JSON in → VCF out."""

import gzip
import io
import json

import pytest

from nirvana2vcf.cli import main
from nirvana2vcf.mapper import decompose_position, map_position_to_vcf_record
from nirvana2vcf.parser import parse_position_line, stream_positions
from nirvana2vcf.vcf_writer import write_vcf_header, write_vcf_record
from tests.conftest import (
    make_test_header as _make_header,
    MINIMAL_HEADER,
    TWO_SAMPLE_HEADER,
    MINIMAL_POSITION_SNV,
    MULTI_ALLELIC_POSITION,
    SV_POSITION,
    MISSING_QUAL_POSITION,
    TWO_SAMPLE_POSITION,
    INSERTION_EXTRA_CONTEXT,
    SNV_LONG_CONTEXT,
    DELETION_EXTRA_SUFFIX,
    THREE_ALT_POSITION,
    build_nirvana_json_string,
)


def _full_pipeline(positions, header_dict=None, **kwargs):
    """Run the full parse→map→write pipeline and return the VCF string."""
    header_dict = header_dict or MINIMAL_HEADER
    header = _make_header(header_dict)

    decompose = kwargs.pop("decompose", False)
    # Separate kwargs for header vs mapper (normalize is mapper-only)
    header_kwargs = {k: v for k, v in kwargs.items() if k != "normalize"}
    out = io.StringIO()
    write_vcf_header(out, header, header.genome_assembly, header.samples, **header_kwargs)

    for pos_dict in positions:
        pos = parse_position_line(json.dumps(pos_dict))
        positions_to_map = decompose_position(pos) if decompose else [pos]
        for p in positions_to_map:
            record = map_position_to_vcf_record(p, header, **kwargs)
            write_vcf_record(out, record)

    return out.getvalue()


class TestEndToEndMinimalSNV:
    def test_header_present(self):
        """Complete conversion has proper VCF header."""
        vcf = _full_pipeline([MINIMAL_POSITION_SNV])

        assert vcf.startswith("##fileformat=VCFv4.2\n")
        assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE001\n" in vcf

    def test_single_data_line(self):
        """Single position produces exactly 1 data line."""
        vcf = _full_pipeline([MINIMAL_POSITION_SNV])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]

        assert len(data_lines) == 1

    def test_data_line_fields(self):
        """Data line has correct CHROM, POS, ID, REF, ALT, QUAL, FILTER."""
        vcf = _full_pipeline([MINIMAL_POSITION_SNV])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        assert fields[0] == "chr1"          # CHROM
        assert fields[1] == "12345"         # POS
        assert fields[2] == "rs12345"       # ID
        assert fields[3] == "A"             # REF
        assert fields[4] == "T"             # ALT
        assert fields[5] == "200"           # QUAL
        assert fields[6] == "PASS"          # FILTER

    def test_info_contains_annotations(self):
        """INFO field contains expected annotation keys."""
        vcf = _full_pipeline([MINIMAL_POSITION_SNV])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        info = data_lines[0].split("\t")[7]

        assert "phyloP=3.5" in info
        assert "gnomAD_AF=0.00012" in info
        assert "oneKG_AF=0.0002" in info
        assert "CLINVAR_SIG=pathogenic" in info
        assert "CSQ=" in info

    def test_sample_column_present(self):
        """Sample column has genotype data."""
        vcf = _full_pipeline([MINIMAL_POSITION_SNV])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        assert len(fields) == 10  # 8 fixed + FORMAT + 1 sample
        assert fields[8].startswith("GT")
        assert fields[9].startswith("0/1")


class TestEndToEndMultiAllelic:
    def test_comma_separated_alt(self):
        """Multi-allelic produces comma-separated ALT."""
        vcf = _full_pipeline([MULTI_ALLELIC_POSITION])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        assert fields[4] == "A,GT"

    def test_per_allele_info(self):
        """Per-allele INFO values are comma-separated."""
        vcf = _full_pipeline([MULTI_ALLELIC_POSITION])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        info = data_lines[0].split("\t")[7]

        assert "gnomAD_AF=0.05,0.001" in info
        assert "phyloP=1.2,-0.5" in info


class TestEndToEndMultiplePositions:
    def test_three_positions(self):
        """Three positions produce three data lines in order."""
        vcf = _full_pipeline(
            [MINIMAL_POSITION_SNV, MULTI_ALLELIC_POSITION, SV_POSITION]
        )
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]

        assert len(data_lines) == 3
        assert data_lines[0].startswith("chr1\t12345")
        assert data_lines[1].startswith("chr2\t48010488")
        assert data_lines[2].startswith("chr7\t55000000")


class TestEndToEndSV:
    def test_sv_info_fields(self):
        """SV position has SVEND, SVLEN, CIPOS, CIEND in INFO."""
        vcf = _full_pipeline([SV_POSITION])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        info = data_lines[0].split("\t")[7]

        assert "SVEND=55050000" in info
        assert "SVLEN=-50000" in info
        assert "CIPOS=-50,50" in info
        assert "CIEND=-100,100" in info


class TestEndToEndMissingFields:
    def test_missing_qual_filter_id(self):
        """Position with missing quality/filter/dbsnp outputs dots."""
        vcf = _full_pipeline([MISSING_QUAL_POSITION])
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        assert fields[2] == "."     # ID
        assert fields[5] == "."     # QUAL
        assert fields[6] == "."     # FILTER


class TestEndToEndTwoSamples:
    def test_two_sample_columns(self):
        """Two samples produce two sample columns."""
        vcf = _full_pipeline(
            [TWO_SAMPLE_POSITION],
            header_dict=TWO_SAMPLE_HEADER,
        )
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        # 8 fixed + FORMAT + 2 samples = 11
        assert len(fields) == 11


class TestEndToEndOutputParseable:
    def test_all_lines_valid(self):
        """Every output line is either a header or has correct tab count."""
        vcf = _full_pipeline(
            [MINIMAL_POSITION_SNV, MULTI_ALLELIC_POSITION, SV_POSITION, MISSING_QUAL_POSITION]
        )
        lines = vcf.strip().split("\n")
        num_samples = len(MINIMAL_HEADER["samples"])

        for line in lines:
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                cols = line.split("\t")
                assert cols[0] == "#CHROM"
            else:
                cols = line.split("\t")
                # Should be 8 fixed + FORMAT + num_samples
                expected = 8 + 1 + num_samples
                assert len(cols) == expected, f"Bad column count: {len(cols)} in line: {line[:80]}"


class TestEndToEndViaFile:
    def test_file_round_trip(self, minimal_snv_json_gz, tmp_path):
        """Full file-based round trip: .json.gz → .vcf"""
        out_path = str(tmp_path / "output.vcf")
        main(["-i", minimal_snv_json_gz, "-o", out_path])

        with open(out_path) as f:
            content = f.read()

        assert content.startswith("##fileformat=VCFv4.2")
        data_lines = [l for l in content.strip().split("\n") if not l.startswith("#")]
        assert len(data_lines) == 1
        assert data_lines[0].startswith("chr1\t12345\trs12345\tA\tT")


class TestEndToEndNormalization:
    """Integration tests for allele normalization through the full pipeline."""

    def test_insertion_normalized_pipeline(self):
        """Full pipeline normalizes insertion with extra context."""
        vcf = _full_pipeline([INSERTION_EXTRA_CONTEXT], normalize=True)
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        assert fields[1] == "2783262"  # POS unchanged
        assert fields[3] == "C"       # REF trimmed from CA
        assert fields[4] == "CA"      # ALT trimmed from CAA

    def test_snv_context_normalized_pipeline(self):
        """Full pipeline normalizes SNV in long context."""
        vcf = _full_pipeline([SNV_LONG_CONTEXT], normalize=True)
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        assert fields[1] == "10800"
        assert fields[3] == "A"
        assert fields[4] == "T"

    def test_no_normalize_flag_preserves(self):
        """--no-normalize preserves raw Nirvana alleles."""
        vcf = _full_pipeline([INSERTION_EXTRA_CONTEXT], normalize=False)
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        fields = data_lines[0].split("\t")

        assert fields[3] == "CA"
        assert fields[4] == "CAA"

    def test_snv_identical_with_or_without_normalize(self):
        """SNV output is identical regardless of normalize flag."""
        vcf_norm = _full_pipeline([MINIMAL_POSITION_SNV], normalize=True)
        vcf_raw = _full_pipeline([MINIMAL_POSITION_SNV], normalize=False)

        # Data lines should be identical
        data_norm = [l for l in vcf_norm.strip().split("\n") if not l.startswith("#")]
        data_raw = [l for l in vcf_raw.strip().split("\n") if not l.startswith("#")]
        assert data_norm == data_raw

    def test_normalize_via_cli(self, tmp_path):
        """CLI --no-normalize flag works end-to-end."""
        json_str = build_nirvana_json_string(
            MINIMAL_HEADER, [INSERTION_EXTRA_CONTEXT]
        )
        json_path = tmp_path / "test.json.gz"
        with gzip.open(json_path, "wt", encoding="utf-8") as f:
            f.write(json_str)

        # Default (normalize=True)
        out_norm = str(tmp_path / "norm.vcf")
        main(["-i", str(json_path), "-o", out_norm])
        with open(out_norm) as f:
            norm_lines = [l for l in f.read().strip().split("\n") if not l.startswith("#")]

        # With --no-normalize
        out_raw = str(tmp_path / "raw.vcf")
        main(["-i", str(json_path), "-o", out_raw, "--no-normalize"])
        with open(out_raw) as f:
            raw_lines = [l for l in f.read().strip().split("\n") if not l.startswith("#")]

        norm_fields = norm_lines[0].split("\t")
        raw_fields = raw_lines[0].split("\t")

        assert norm_fields[3] == "C"    # Normalized REF
        assert norm_fields[4] == "CA"   # Normalized ALT
        assert raw_fields[3] == "CA"    # Original REF
        assert raw_fields[4] == "CAA"   # Original ALT


class TestEndToEndDecomposition:
    """Integration tests for multi-allelic decomposition through the full pipeline."""

    def test_biallelic_unchanged(self):
        """Single-ALT position unchanged with decompose."""
        vcf = _full_pipeline([MINIMAL_POSITION_SNV], decompose=True)
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        assert len(data_lines) == 1
        assert data_lines[0].split("\t")[4] == "T"

    def test_multi_allelic_decomposed_to_two_rows(self):
        """Two-ALT position produces two VCF rows with decompose."""
        vcf = _full_pipeline([MULTI_ALLELIC_POSITION], decompose=True)
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        assert len(data_lines) == 2
        assert data_lines[0].split("\t")[4] == "A"
        assert data_lines[1].split("\t")[4] == "GT"

    def test_decomposed_info_single_values(self):
        """Per-allele INFO fields have single values (no commas) after decomposition."""
        vcf = _full_pipeline([MULTI_ALLELIC_POSITION], decompose=True)
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        info1 = data_lines[0].split("\t")[7]
        info2 = data_lines[1].split("\t")[7]

        assert "gnomAD_AF=0.05" in info1
        assert "gnomAD_AF=0.001" in info2

    def test_decomposed_csq_scoped(self):
        """CSQ only contains transcripts for the relevant ALT."""
        vcf = _full_pipeline([MULTI_ALLELIC_POSITION], decompose=True)
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        info1 = data_lines[0].split("\t")[7]
        info2 = data_lines[1].split("\t")[7]

        assert "synonymous_variant" in info1
        assert "frameshift_variant" not in info1
        assert "frameshift_variant" in info2
        assert "synonymous_variant" not in info2

    def test_decomposed_gt_in_output(self):
        """Decomposed sample GT is remapped in output."""
        vcf = _full_pipeline([MULTI_ALLELIC_POSITION], decompose=True)
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]

        # Original GT: 1/2
        sample1 = data_lines[0].split("\t")[9]
        sample2 = data_lines[1].split("\t")[9]
        assert sample1.startswith("1/.")
        assert sample2.startswith("./1")

    def test_decompose_plus_normalize(self):
        """Decomposition + normalization produces correct biallelic rows."""
        vcf = _full_pipeline([MULTI_ALLELIC_POSITION], decompose=True, normalize=True)
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        assert len(data_lines) == 2

    def test_without_decompose_unchanged(self):
        """Without decompose flag, multi-allelic row stays as one."""
        vcf = _full_pipeline([MULTI_ALLELIC_POSITION], decompose=False)
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        assert len(data_lines) == 1
        assert data_lines[0].split("\t")[4] == "A,GT"

    def test_decompose_preserves_total_row_count(self):
        """Mixed input: total output rows = sum of per-position ALT counts."""
        vcf = _full_pipeline(
            [MINIMAL_POSITION_SNV, MULTI_ALLELIC_POSITION, SV_POSITION],
            decompose=True,
        )
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        # 1 (SNV) + 2 (multi-allelic) + 1 (SV) = 4
        assert len(data_lines) == 4

    def test_three_alt_decomposed(self):
        """Three-ALT position produces three VCF rows."""
        vcf = _full_pipeline([THREE_ALT_POSITION], decompose=True)
        data_lines = [l for l in vcf.strip().split("\n") if not l.startswith("#")]
        assert len(data_lines) == 3
        assert data_lines[0].split("\t")[4] == "T"
        assert data_lines[1].split("\t")[4] == "C"
        assert data_lines[2].split("\t")[4] == "G"

    def test_decompose_via_cli(self, tmp_path):
        """CLI --decompose flag works end-to-end."""
        json_str = build_nirvana_json_string(
            MINIMAL_HEADER, [MULTI_ALLELIC_POSITION]
        )
        json_path = tmp_path / "test.json.gz"
        with gzip.open(json_path, "wt", encoding="utf-8") as f:
            f.write(json_str)

        # Default (no decompose)
        out_default = str(tmp_path / "default.vcf")
        main(["-i", str(json_path), "-o", out_default])
        with open(out_default) as f:
            default_lines = [l for l in f.read().strip().split("\n") if not l.startswith("#")]

        # With --decompose
        out_decompose = str(tmp_path / "decompose.vcf")
        main(["-i", str(json_path), "-o", out_decompose, "--decompose"])
        with open(out_decompose) as f:
            decompose_lines = [l for l in f.read().strip().split("\n") if not l.startswith("#")]

        assert len(default_lines) == 1      # Multi-allelic: one row
        assert len(decompose_lines) == 2    # Decomposed: two rows
        assert decompose_lines[0].split("\t")[4] == "A"
        assert decompose_lines[1].split("\t")[4] == "GT"

    def test_decompose_plus_normalize_via_cli(self, tmp_path):
        """CLI --decompose + --normalize works end-to-end."""
        json_str = build_nirvana_json_string(
            MINIMAL_HEADER, [MULTI_ALLELIC_POSITION]
        )
        json_path = tmp_path / "test.json.gz"
        with gzip.open(json_path, "wt", encoding="utf-8") as f:
            f.write(json_str)

        out_path = str(tmp_path / "output.vcf")
        main(["-i", str(json_path), "-o", out_path, "--decompose"])

        with open(out_path) as f:
            data_lines = [l for l in f.read().strip().split("\n") if not l.startswith("#")]

        assert len(data_lines) == 2
