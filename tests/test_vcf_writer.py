"""Tests for the VCF output writer."""

import io

import pytest

from nirvana2vcf.vcf_writer import (
    get_contig_header_lines,
    write_vcf_header,
    write_vcf_record,
)
from nirvana2vcf.mapper import map_position_to_vcf_record
from tests.conftest import (
    make_test_header as _make_header,
    parse_pos as _parse_pos,
    MINIMAL_POSITION_SNV,
    MULTI_ALLELIC_POSITION,
    SV_POSITION,
)


class TestVCFHeader:
    def test_starts_with_fileformat(self):
        """First line is ##fileformat=VCFv4.2"""
        out = io.StringIO()
        write_vcf_header(out, _make_header(), "GRCh38", ["SAMPLE001"])
        lines = out.getvalue().split("\n")

        assert lines[0] == "##fileformat=VCFv4.2"

    def test_source_line(self):
        """Source line includes annotator info."""
        out = io.StringIO()
        write_vcf_header(out, _make_header(), "GRCh38", ["SAMPLE001"])
        lines = out.getvalue().split("\n")

        assert any("##source=" in l and "Nirvana" in l for l in lines)

    def test_contains_info_definitions(self):
        """All INFO definitions are present in header."""
        out = io.StringIO()
        write_vcf_header(out, _make_header(), "GRCh38", ["SAMPLE001"])
        content = out.getvalue()

        assert "##INFO=<ID=CSQ," in content
        assert "##INFO=<ID=phyloP," in content
        assert "##INFO=<ID=gnomAD_AF," in content
        assert "##INFO=<ID=CLINVAR_SIG," in content
        assert "##INFO=<ID=SVEND," in content

    def test_contains_format_definitions(self):
        """FORMAT definitions are present when samples exist."""
        out = io.StringIO()
        write_vcf_header(out, _make_header(), "GRCh38", ["SAMPLE001"])
        content = out.getvalue()

        assert "##FORMAT=<ID=GT," in content
        assert "##FORMAT=<ID=DP," in content
        assert "##FORMAT=<ID=GQ," in content
        assert "##FORMAT=<ID=AD," in content

    def test_no_format_without_samples(self):
        """FORMAT definitions are omitted when no samples."""
        out = io.StringIO()
        write_vcf_header(out, _make_header(), "GRCh38", [])
        content = out.getvalue()

        assert "##FORMAT=" not in content

    def test_contig_lines_grch38(self):
        """GRCh38 assembly produces chr1-chrM contig lines."""
        out = io.StringIO()
        write_vcf_header(out, _make_header(), "GRCh38", [])
        content = out.getvalue()

        assert "##contig=<ID=chr1,length=248956422>" in content
        assert "##contig=<ID=chrX," in content
        assert "##contig=<ID=chrM," in content

    def test_contig_lines_grch37(self):
        """GRCh37 assembly produces 1-MT contig lines (no chr prefix)."""
        lines = get_contig_header_lines("GRCh37")

        assert any("ID=1," in l for l in lines)
        assert any("ID=X," in l for l in lines)
        assert any("ID=MT," in l for l in lines)
        # Should not have chr prefix
        assert not any("ID=chr" in l for l in lines)

    def test_column_header_with_samples(self):
        """#CHROM line includes FORMAT and sample names."""
        out = io.StringIO()
        write_vcf_header(out, _make_header(), "GRCh38", ["SAMPLE001"])
        lines = out.getvalue().strip().split("\n")
        chrom_line = [l for l in lines if l.startswith("#CHROM")][0]
        cols = chrom_line.split("\t")

        assert cols[0] == "#CHROM"
        assert cols[1] == "POS"
        assert cols[8] == "FORMAT"
        assert cols[9] == "SAMPLE001"

    def test_column_header_without_samples(self):
        """#CHROM line has 8 columns when no samples."""
        out = io.StringIO()
        write_vcf_header(out, _make_header(), "GRCh38", [])
        lines = out.getvalue().strip().split("\n")
        chrom_line = [l for l in lines if l.startswith("#CHROM")][0]
        cols = chrom_line.split("\t")

        assert len(cols) == 8
        assert cols[-1] == "INFO"

    def test_csq_only_header(self):
        """csq_only mode only includes CSQ INFO definition."""
        out = io.StringIO()
        write_vcf_header(out, _make_header(), "GRCh38", [], csq_only=True)
        content = out.getvalue()

        assert "##INFO=<ID=CSQ," in content
        assert "##INFO=<ID=phyloP," not in content
        assert "##INFO=<ID=gnomAD_AF," not in content


class TestWriteRecord:
    def test_tab_separated(self):
        """Data line is tab-separated."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header())
        out = io.StringIO()
        write_vcf_record(out, record)

        line = out.getvalue().strip()
        fields = line.split("\t")
        assert len(fields) >= 8

    def test_correct_field_values(self):
        """Data line has correct basic field values."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header())
        out = io.StringIO()
        write_vcf_record(out, record)

        fields = out.getvalue().strip().split("\t")
        assert fields[0] == "chr1"
        assert fields[1] == "12345"
        assert fields[2] == "rs12345"
        assert fields[3] == "A"
        assert fields[4] == "T"
        assert fields[5] == "200"
        assert fields[6] == "PASS"

    def test_record_with_samples(self):
        """Data line includes FORMAT and sample columns."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header())
        out = io.StringIO()
        write_vcf_record(out, record)

        fields = out.getvalue().strip().split("\t")
        # Should have CHROM-INFO (8) + FORMAT (1) + sample (1) = 10
        assert len(fields) == 10
        assert "GT" in fields[8]

    def test_record_without_samples(self):
        """Data line without samples has 8 fields."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header(), include_samples=False)
        out = io.StringIO()
        write_vcf_record(out, record, include_samples=False)

        fields = out.getvalue().strip().split("\t")
        assert len(fields) == 8

    def test_write_to_stringio(self):
        """Writer correctly writes to a StringIO object."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header())
        out = io.StringIO()

        write_vcf_record(out, record)

        result = out.getvalue()
        assert result.endswith("\n")
        assert "chr1" in result
