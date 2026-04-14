"""Edge case tests for unusual inputs and boundary conditions."""

import pytest

from nirvana2vcf.mapper import (
    build_info_field,
    build_sample_columns,
    map_position_to_vcf_record,
)
from nirvana2vcf.parser import parse_variant, parse_sample
from nirvana2vcf.models import Position, Variant, Sample
from tests.conftest import (
    make_test_header as _make_header,
    parse_pos as _parse_pos,
    REFERENCE_ONLY_POSITION,
    EMPTY_SAMPLE_POSITION,
)


class TestReferenceOnlyPosition:
    def test_reference_only_produces_vcf_line(self):
        """Position with altAlleles=['.'] outputs ALT='.'"""
        pos = _parse_pos(REFERENCE_ONLY_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["ALT"] == "."

    def test_reference_only_info_dot(self):
        """Reference-only position has INFO='.'"""
        pos = _parse_pos(REFERENCE_ONLY_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["INFO"] == "."


class TestNoVariantsArray:
    def test_position_without_variants(self):
        """Position without 'variants' key still produces a valid record."""
        pos_data = {
            "chromosome": "chr4",
            "position": 5000,
            "refAllele": "A",
            "altAlleles": ["C"],
            "quality": 30.0,
            "filters": ["PASS"],
            "samples": [{"genotype": "0/1", "totalDepth": 20}],
        }
        pos = _parse_pos(pos_data)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["CHROM"] == "chr4"
        assert record["POS"] == 5000
        assert record["ID"] == "."
        # INFO should be "." since no variants with annotations
        assert record["INFO"] == "."


class TestBooleanOmission:
    def test_absent_is_de_novo_means_false(self):
        """isDeNovo absent from sample means None, not error."""
        sample_data = {"genotype": "0/1", "totalDepth": 25}
        sample = parse_sample(sample_data)

        assert sample.is_de_novo is None
        assert sample.failed_filter is None

    def test_absent_boolean_in_variant(self):
        """isStructuralVariant absent means None."""
        variant_data = {
            "vid": "1-100-A-T",
            "chromosome": "chr1",
            "begin": 100,
            "end": 100,
            "refAllele": "A",
            "altAllele": "T",
            "variantType": "SNV",
        }
        variant = parse_variant(variant_data)

        assert variant.is_structural_variant is None
        assert variant.is_decomposed_variant is None


class TestLongAnnotationStrings:
    def test_clinvar_phenotype_with_special_chars(self):
        """ClinVar phenotype with semicolons and commas is escaped."""
        pos_data = {
            "chromosome": "chr17",
            "position": 41245466,
            "refAllele": "G",
            "altAlleles": ["A"],
            "quality": 500.0,
            "filters": ["PASS"],
            "samples": [{"genotype": "0/1"}],
            "variants": [
                {
                    "vid": "17-41245466-G-A",
                    "chromosome": "chr17",
                    "begin": 41245466,
                    "end": 41245466,
                    "refAllele": "G",
                    "altAllele": "A",
                    "variantType": "SNV",
                    "clinvar": [
                        {
                            "id": "RCV000999",
                            "significance": ["pathogenic"],
                            "reviewStatus": "criteria provided, multiple submitters, no conflicts",
                            "phenotypes": ["Hereditary breast cancer"],
                        }
                    ],
                }
            ],
        }
        pos = _parse_pos(pos_data)
        info = build_info_field(pos)

        # Review status should have spaces and commas escaped
        assert "CLINVAR_REVSTAT=" in info
        assert "%20" in info  # space
        assert "%2C" in info  # comma


class TestEmptySamplePlaceholder:
    def test_empty_sample_genotype(self):
        """isEmpty sample outputs ./. for GT."""
        pos = _parse_pos(EMPTY_SAMPLE_POSITION)
        format_str, samples = build_sample_columns(pos.samples)

        assert samples[0].startswith("./.")

    def test_empty_sample_other_fields_dot(self):
        """isEmpty sample with no depth outputs '.' for DP."""
        pos = _parse_pos(EMPTY_SAMPLE_POSITION)
        format_str, samples = build_sample_columns(pos.samples)

        # With only GT and DP available, DP should be "."
        parts = format_str.split(":")
        sample_parts = samples[0].split(":")
        if "DP" in parts:
            dp_idx = parts.index("DP")
            assert sample_parts[dp_idx] == "."


class TestMultipleClinVarEntries:
    def test_multiple_clinvar_concatenated(self):
        """Multiple ClinVar entries are encoded as &-separated in INFO."""
        pos_data = {
            "chromosome": "chr13",
            "position": 32914438,
            "refAllele": "T",
            "altAlleles": ["C"],
            "quality": 400.0,
            "filters": ["PASS"],
            "samples": [{"genotype": "0/1"}],
            "variants": [
                {
                    "vid": "13-32914438-T-C",
                    "chromosome": "chr13",
                    "begin": 32914438,
                    "end": 32914438,
                    "refAllele": "T",
                    "altAllele": "C",
                    "variantType": "SNV",
                    "clinvar": [
                        {
                            "id": "RCV000001",
                            "significance": ["pathogenic"],
                            "reviewStatus": "practice guideline",
                        },
                        {
                            "id": "RCV000002",
                            "significance": ["likely pathogenic"],
                            "reviewStatus": "criteria provided, single submitter",
                        },
                    ],
                }
            ],
        }
        pos = _parse_pos(pos_data)
        info = build_info_field(pos)

        assert "CLINVAR_ID=RCV000001&RCV000002" in info
        assert "CLINVAR_SIG=pathogenic&likely%20pathogenic" in info


class TestMixedSampleData:
    def test_sample_with_only_genotype(self):
        """Sample with only genotype has minimal FORMAT."""
        pos_data = {
            "chromosome": "chr1",
            "position": 100,
            "refAllele": "A",
            "altAlleles": ["T"],
            "samples": [{"genotype": "0/1"}],
            "variants": [
                {
                    "vid": "1-100-A-T",
                    "chromosome": "chr1",
                    "begin": 100,
                    "end": 100,
                    "refAllele": "A",
                    "altAllele": "T",
                    "variantType": "SNV",
                }
            ],
        }
        pos = _parse_pos(pos_data)
        format_str, samples = build_sample_columns(pos.samples)

        assert format_str == "GT"
        assert samples[0] == "0/1"


class TestVariantWithRevel:
    def test_revel_as_dict(self):
        """REVEL score from dict format is parsed."""
        variant_data = {
            "vid": "1-100-A-T",
            "chromosome": "chr1",
            "begin": 100,
            "end": 100,
            "refAllele": "A",
            "altAllele": "T",
            "variantType": "SNV",
            "revel": {"score": 0.85},
        }
        variant = parse_variant(variant_data)
        assert variant.revel_score == 0.85

    def test_revel_in_info(self):
        """REVEL score appears in INFO."""
        pos_data = {
            "chromosome": "chr1",
            "position": 100,
            "refAllele": "A",
            "altAlleles": ["T"],
            "samples": [{"genotype": "0/1"}],
            "variants": [
                {
                    "vid": "1-100-A-T",
                    "chromosome": "chr1",
                    "begin": 100,
                    "end": 100,
                    "refAllele": "A",
                    "altAllele": "T",
                    "variantType": "SNV",
                    "revel": {"score": 0.85},
                }
            ],
        }
        pos = _parse_pos(pos_data)
        info = build_info_field(pos)

        assert "REVEL=0.85" in info


class TestPhasedGenotype:
    def test_phased_genotype_preserved(self):
        """Phased genotype (|) is preserved in output."""
        sample_data = {"genotype": "0|1", "totalDepth": 30}
        sample = parse_sample(sample_data)

        assert sample.genotype == "0|1"


class TestDeNovoSample:
    def test_de_novo_in_format(self):
        """De novo sample includes DN and DQ in FORMAT."""
        pos_data = {
            "chromosome": "chr1",
            "position": 200,
            "refAllele": "G",
            "altAlleles": ["A"],
            "samples": [
                {
                    "genotype": "0/1",
                    "totalDepth": 50,
                    "isDeNovo": True,
                    "deNovoQuality": 95.5,
                }
            ],
            "variants": [
                {
                    "vid": "1-200-G-A",
                    "chromosome": "chr1",
                    "begin": 200,
                    "end": 200,
                    "refAllele": "G",
                    "altAllele": "A",
                    "variantType": "SNV",
                }
            ],
        }
        pos = _parse_pos(pos_data)
        format_str, samples = build_sample_columns(pos.samples)

        assert "DN" in format_str
        assert "DQ" in format_str
        parts = format_str.split(":")
        sample_parts = samples[0].split(":")
        dn_idx = parts.index("DN")
        assert sample_parts[dn_idx] == "true"
