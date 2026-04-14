"""Tests for the streaming Nirvana JSON parser."""

import gzip
import json
import io

import pytest

from nirvana2vcf.parser import (
    parse_header,
    parse_position_line,
    parse_sample,
    parse_variant,
    stream_positions,
)
from tests.conftest import (
    MINIMAL_HEADER,
    MINIMAL_POSITION_SNV,
    MULTI_ALLELIC_POSITION,
    MISSING_QUAL_POSITION,
    SV_POSITION,
    EMPTY_SAMPLE_POSITION,
    SPLICE_AI_POSITION,
    build_nirvana_json_string,
)


class TestParseHeader:
    def test_extracts_all_fields(self):
        """Header parsing extracts annotator, assembly, samples, etc."""
        header_json = json.dumps({"header": MINIMAL_HEADER})
        first_line = header_json[:-1] + ',"positions":['

        header = parse_header(first_line)

        assert header.annotator == "Nirvana 3.22.0"
        assert header.genome_assembly == "GRCh38"
        assert header.schema_version == 6
        assert header.data_version == "91.26.50"
        assert header.samples == ["SAMPLE001"]
        assert header.creation_time == "2024-06-15 10:30:00"
        assert len(header.data_sources) == 1
        assert header.data_sources[0]["name"] == "VEP"

    def test_raises_on_invalid_first_line(self):
        """Invalid first line raises ValueError."""
        with pytest.raises(ValueError, match="positions"):
            parse_header('{"header": {"annotator": "test"}}')


class TestParsePositionLine:
    def test_parse_snv_position(self):
        """Single SNV position parses correctly."""
        line = json.dumps(MINIMAL_POSITION_SNV) + ","
        pos = parse_position_line(line)

        assert pos is not None
        assert pos.chromosome == "chr1"
        assert pos.position == 12345
        assert pos.ref_allele == "A"
        assert pos.alt_alleles == ["T"]
        assert pos.quality == 200.0
        assert pos.filters == ["PASS"]
        assert pos.cytogenetic_band == "1p36.33"

    def test_parse_multi_allelic(self):
        """Position with two alt alleles produces two Variant objects."""
        line = json.dumps(MULTI_ALLELIC_POSITION)
        pos = parse_position_line(line)

        assert pos is not None
        assert pos.alt_alleles == ["A", "GT"]
        assert len(pos.variants) == 2
        assert pos.variants[0].alt_allele == "A"
        assert pos.variants[1].alt_allele == "GT"

    def test_parse_missing_quality(self):
        """Position without quality/filters produces None values."""
        line = json.dumps(MISSING_QUAL_POSITION)
        pos = parse_position_line(line)

        assert pos is not None
        assert pos.quality is None
        assert pos.filters is None

    def test_parse_sv_position(self):
        """Structural variant position parses SV fields."""
        line = json.dumps(SV_POSITION)
        pos = parse_position_line(line)

        assert pos is not None
        assert pos.sv_end == 55050000
        assert pos.sv_length == -50000
        assert pos.ci_pos == [-50, 50]
        assert pos.ci_end == [-100, 100]

    def test_returns_none_at_genes_boundary(self):
        """Line starting with ] returns None (end of positions)."""
        assert parse_position_line('],"genes":[') is None
        assert parse_position_line("]}") is None
        assert parse_position_line("") is None

    def test_strips_trailing_comma(self):
        """Lines with trailing comma are parsed correctly."""
        line_with_comma = json.dumps(MINIMAL_POSITION_SNV) + ","
        line_without_comma = json.dumps(MINIMAL_POSITION_SNV)

        pos1 = parse_position_line(line_with_comma)
        pos2 = parse_position_line(line_without_comma)

        assert pos1.chromosome == pos2.chromosome
        assert pos1.position == pos2.position


class TestParseVariant:
    def test_parse_variant_annotations(self):
        """Variant with full annotations parses all nested structures."""
        variant_data = MINIMAL_POSITION_SNV["variants"][0]
        variant = parse_variant(variant_data)

        assert variant.vid == "1-12345-A-T"
        assert variant.variant_type == "SNV"
        assert variant.phylop_score == 3.5
        assert variant.dbsnp == ["rs12345"]
        assert variant.hgvsg == "NC_000001.11:g.12345A>T"

    def test_parse_variant_gnomad(self):
        """gnomAD population frequency is parsed."""
        variant_data = MINIMAL_POSITION_SNV["variants"][0]
        variant = parse_variant(variant_data)

        assert variant.gnomad is not None
        assert variant.gnomad.all_af == 0.00012
        assert variant.gnomad.all_ac == 15
        assert variant.gnomad.all_an == 125000

    def test_parse_variant_onekg(self):
        """1000 Genomes frequency is parsed."""
        variant_data = MINIMAL_POSITION_SNV["variants"][0]
        variant = parse_variant(variant_data)

        assert variant.one_kg is not None
        assert variant.one_kg.all_af == 0.0002

    def test_parse_variant_clinvar(self):
        """ClinVar entries are parsed."""
        variant_data = MINIMAL_POSITION_SNV["variants"][0]
        variant = parse_variant(variant_data)

        assert variant.clinvar is not None
        assert len(variant.clinvar) == 1
        assert variant.clinvar[0].id == "RCV000012345"
        assert variant.clinvar[0].significance == ["pathogenic"]
        assert variant.clinvar[0].phenotypes == ["Some disease"]

    def test_parse_variant_transcripts(self):
        """Transcript annotations are parsed."""
        variant_data = MINIMAL_POSITION_SNV["variants"][0]
        variant = parse_variant(variant_data)

        assert variant.transcripts is not None
        assert len(variant.transcripts) == 1
        t = variant.transcripts[0]
        assert t.transcript == "NM_001234.5"
        assert t.source == "RefSeq"
        assert t.hgnc == "GENE1"
        assert t.consequence == ["missense_variant"]
        assert t.amino_acids == "M/L"
        assert t.hgvsc == "NM_001234.5:c.50A>T"
        assert t.is_canonical is True

    def test_parse_variant_splice_ai(self):
        """SpliceAI entries are parsed."""
        variant_data = SPLICE_AI_POSITION["variants"][0]
        variant = parse_variant(variant_data)

        assert variant.splice_ai is not None
        assert len(variant.splice_ai) == 1
        sai = variant.splice_ai[0]
        assert sai.hgnc == "GENE2"
        assert sai.donor_gain_score == 0.8
        assert sai.donor_gain_distance == -2

    def test_parse_variant_no_annotations(self):
        """Variant with no annotations has None for optional fields."""
        variant_data = MISSING_QUAL_POSITION["variants"][0]
        variant = parse_variant(variant_data)

        assert variant.gnomad is None
        assert variant.clinvar is None
        assert variant.transcripts is None
        assert variant.dbsnp is None


class TestParseSample:
    def test_parse_full_sample(self):
        """Sample with all genotype fields parses correctly."""
        sample_data = MINIMAL_POSITION_SNV["samples"][0]
        sample = parse_sample(sample_data)

        assert sample.genotype == "0/1"
        assert sample.total_depth == 40
        assert sample.genotype_quality == 99
        assert sample.variant_frequencies == [0.5]
        assert sample.allele_depths == [20, 20]

    def test_parse_empty_sample(self):
        """Sample with isEmpty=True parses correctly."""
        sample_data = EMPTY_SAMPLE_POSITION["samples"][0]
        sample = parse_sample(sample_data)

        assert sample.is_empty is True
        assert sample.genotype == "./."

    def test_parse_sv_sample(self):
        """SV sample with split/paired-end reads parses correctly."""
        sample_data = SV_POSITION["samples"][0]
        sample = parse_sample(sample_data)

        assert sample.paired_end_read_counts == [15, 10]
        assert sample.split_read_counts == [12, 8]


class TestStreamPositions:
    def test_stream_from_json_file(self, minimal_snv_json_file):
        """One position in uncompressed file produces one Position."""
        results = list(stream_positions(minimal_snv_json_file))

        assert len(results) == 1
        header, pos = results[0]
        assert header.genome_assembly == "GRCh38"
        assert pos.chromosome == "chr1"
        assert pos.position == 12345

    def test_stream_from_gz_file(self, minimal_snv_json_gz):
        """Gzipped input produces identical result to uncompressed."""
        results = list(stream_positions(minimal_snv_json_gz))

        assert len(results) == 1
        header, pos = results[0]
        assert header.genome_assembly == "GRCh38"
        assert pos.chromosome == "chr1"

    def test_stream_multiple_positions(self, tmp_path, mixed_positions_json_string):
        """Three positions produce three Position objects in order."""
        path = tmp_path / "multi.json"
        with open(path, "w") as f:
            f.write(mixed_positions_json_string)

        results = list(stream_positions(str(path)))

        assert len(results) == 3
        assert results[0][1].chromosome == "chr1"
        assert results[1][1].chromosome == "chr2"
        assert results[2][1].chromosome == "chr7"

    def test_header_same_object_for_all_positions(self, tmp_path, mixed_positions_json_string):
        """Header object is the same instance for every yielded position."""
        path = tmp_path / "multi.json"
        with open(path, "w") as f:
            f.write(mixed_positions_json_string)

        results = list(stream_positions(str(path)))
        headers = [r[0] for r in results]

        assert all(h is headers[0] for h in headers)

    def test_stops_at_genes_section(self, tmp_path):
        """Parser stops at genes section and does not parse gene objects."""
        genes = [{"name": "BRCA1", "hgncId": 1100}]
        content = build_nirvana_json_string(
            MINIMAL_HEADER, [MINIMAL_POSITION_SNV], gene_dicts=genes
        )
        path = tmp_path / "with_genes.json"
        with open(path, "w") as f:
            f.write(content)

        results = list(stream_positions(str(path)))

        # Should get exactly 1 position, not try to parse genes
        assert len(results) == 1
