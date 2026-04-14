"""Tests for the CLI interface."""

import subprocess
import sys

import pytest

from nirvana2vcf.cli import main


class TestCLI:
    def test_basic_invocation(self, minimal_snv_json_gz, tmp_path):
        """nirvana2vcf -i input.json.gz -o output.vcf produces output file."""
        out_path = str(tmp_path / "output.vcf")
        main(["-i", minimal_snv_json_gz, "-o", out_path])

        with open(out_path) as f:
            content = f.read()

        assert content.startswith("##fileformat=VCFv4.2")
        assert "#CHROM" in content
        assert "chr1\t12345" in content

    def test_output_to_stdout(self, minimal_snv_json_gz, capsys):
        """nirvana2vcf -i input.json.gz (no -o) writes to stdout."""
        main(["-i", minimal_snv_json_gz])
        captured = capsys.readouterr()

        assert "##fileformat=VCFv4.2" in captured.out
        assert "chr1\t12345" in captured.out

    def test_missing_input_file(self):
        """nirvana2vcf -i nonexistent.json.gz exits with error."""
        with pytest.raises(SystemExit) as exc_info:
            main(["-i", "/nonexistent/path.json.gz"])
        assert exc_info.value.code == 1

    def test_csq_only_flag(self, minimal_snv_json_gz, tmp_path):
        """--csq-only produces INFO with only CSQ, no flat fields."""
        out_path = str(tmp_path / "output.vcf")
        main(["-i", minimal_snv_json_gz, "-o", out_path, "--csq-only"])

        with open(out_path) as f:
            content = f.read()

        # Header should only have CSQ INFO definition
        assert "##INFO=<ID=CSQ," in content
        assert "##INFO=<ID=phyloP," not in content

        # Data line INFO should start with CSQ=
        data_lines = [l for l in content.strip().split("\n") if not l.startswith("#")]
        assert len(data_lines) == 1
        fields = data_lines[0].split("\t")
        assert fields[7].startswith("CSQ=")

    def test_no_samples_flag(self, minimal_snv_json_gz, tmp_path):
        """--no-samples omits FORMAT and sample columns."""
        out_path = str(tmp_path / "output.vcf")
        main(["-i", minimal_snv_json_gz, "-o", out_path, "--no-samples"])

        with open(out_path) as f:
            content = f.read()

        # No FORMAT definitions in header
        assert "##FORMAT=" not in content

        # Column header should not have FORMAT or SAMPLE001
        chrom_lines = [l for l in content.split("\n") if l.startswith("#CHROM")]
        assert len(chrom_lines) == 1
        cols = chrom_lines[0].split("\t")
        assert len(cols) == 8

        # Data lines should have 8 columns
        data_lines = [l for l in content.strip().split("\n") if not l.startswith("#")]
        for line in data_lines:
            assert len(line.split("\t")) == 8

    def test_version_flag(self, capsys):
        """--version prints version and exits."""
        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])
        assert exc_info.value.code == 0

    def test_assembly_override(self, minimal_snv_json_gz, tmp_path):
        """--assembly GRCh37 uses GRCh37 contigs."""
        out_path = str(tmp_path / "output.vcf")
        main(["-i", minimal_snv_json_gz, "-o", out_path, "--assembly", "GRCh37"])

        with open(out_path) as f:
            content = f.read()

        # GRCh37 uses "1" not "chr1" for contig names
        assert "##contig=<ID=1," in content
        assert "##contig=<ID=chr1," not in content

    def test_multi_position_file(self, multi_position_json_gz, tmp_path):
        """Multiple positions produce multiple data lines."""
        out_path = str(tmp_path / "output.vcf")
        main(["-i", multi_position_json_gz, "-o", out_path])

        with open(out_path) as f:
            content = f.read()

        data_lines = [l for l in content.strip().split("\n") if not l.startswith("#")]
        assert len(data_lines) == 3

    def test_uncompressed_input(self, minimal_snv_json_file, tmp_path):
        """Uncompressed .json input works."""
        out_path = str(tmp_path / "output.vcf")
        main(["-i", minimal_snv_json_file, "-o", out_path])

        with open(out_path) as f:
            content = f.read()

        assert "chr1\t12345" in content
