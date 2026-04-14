"""Shared test fixtures with all test data embedded as Python dicts."""

import gzip
import json

import pytest


# ============================================================
# Header
# ============================================================

MINIMAL_HEADER = {
    "annotator": "Nirvana 3.22.0",
    "creationTime": "2024-06-15 10:30:00",
    "genomeAssembly": "GRCh38",
    "schemaVersion": 6,
    "dataVersion": "91.26.50",
    "dataSources": [
        {
            "name": "VEP",
            "version": "91",
            "description": "Ensembl",
            "releaseDate": "2024-01-01",
        }
    ],
    "samples": ["SAMPLE001"],
}

TWO_SAMPLE_HEADER = {
    **MINIMAL_HEADER,
    "samples": ["SAMPLE001", "SAMPLE002"],
}

# ============================================================
# Positions
# ============================================================

# 1. Minimal SNV — single sample, full annotations
MINIMAL_POSITION_SNV = {
    "chromosome": "chr1",
    "position": 12345,
    "refAllele": "A",
    "altAlleles": ["T"],
    "quality": 200.0,
    "filters": ["PASS"],
    "cytogeneticBand": "1p36.33",
    "samples": [
        {
            "genotype": "0/1",
            "variantFrequencies": [0.5],
            "totalDepth": 40,
            "genotypeQuality": 99,
            "alleleDepths": [20, 20],
        }
    ],
    "variants": [
        {
            "vid": "1-12345-A-T",
            "chromosome": "chr1",
            "begin": 12345,
            "end": 12345,
            "refAllele": "A",
            "altAllele": "T",
            "variantType": "SNV",
            "hgvsg": "NC_000001.11:g.12345A>T",
            "phylopScore": 3.5,
            "dbsnp": ["rs12345"],
            "transcripts": [
                {
                    "transcript": "NM_001234.5",
                    "source": "RefSeq",
                    "bioType": "protein_coding",
                    "geneId": "1234",
                    "hgnc": "GENE1",
                    "consequence": ["missense_variant"],
                    "codons": "Atg/Ttg",
                    "aminoAcids": "M/L",
                    "cdnaPos": "100",
                    "cdsPos": "50",
                    "proteinPos": "17",
                    "hgvsc": "NM_001234.5:c.50A>T",
                    "hgvsp": "NP_001234.1:p.(Met17Leu)",
                    "isCanonical": True,
                }
            ],
            "clinvar": [
                {
                    "id": "RCV000012345",
                    "variationId": "12345",
                    "reviewStatus": "criteria provided, single submitter",
                    "significance": ["pathogenic"],
                    "phenotypes": ["Some disease"],
                    "lastUpdatedDate": "2024-01-15",
                }
            ],
            "gnomad": {
                "allAf": 0.00012,
                "allAc": 15,
                "allAn": 125000,
            },
            "oneKg": {
                "allAf": 0.0002,
                "allAn": 5008,
            },
        }
    ],
}

# 2. Multi-allelic — two alt alleles
MULTI_ALLELIC_POSITION = {
    "chromosome": "chr2",
    "position": 48010488,
    "refAllele": "G",
    "altAlleles": ["A", "GT"],
    "quality": 150.0,
    "filters": ["PASS"],
    "samples": [
        {
            "genotype": "1/2",
            "variantFrequencies": [0.333, 0.5],
            "totalDepth": 57,
            "genotypeQuality": 85,
            "alleleDepths": [10, 19, 28],
        }
    ],
    "variants": [
        {
            "vid": "2-48010488-G-A",
            "chromosome": "chr2",
            "begin": 48010488,
            "end": 48010488,
            "refAllele": "G",
            "altAllele": "A",
            "variantType": "SNV",
            "phylopScore": 1.2,
            "dbsnp": ["rs9999"],
            "transcripts": [
                {
                    "transcript": "NM_005678.3",
                    "source": "RefSeq",
                    "bioType": "protein_coding",
                    "geneId": "5678",
                    "hgnc": "MSH6",
                    "consequence": ["synonymous_variant"],
                    "isCanonical": True,
                }
            ],
            "gnomad": {"allAf": 0.05},
        },
        {
            "vid": "2-48010488-G-GT",
            "chromosome": "chr2",
            "begin": 48010488,
            "end": 48010488,
            "refAllele": "G",
            "altAllele": "GT",
            "variantType": "insertion",
            "phylopScore": -0.5,
            "transcripts": [
                {
                    "transcript": "NM_005678.3",
                    "source": "RefSeq",
                    "bioType": "protein_coding",
                    "geneId": "5678",
                    "hgnc": "MSH6",
                    "consequence": ["frameshift_variant"],
                    "isCanonical": True,
                }
            ],
            "gnomad": {"allAf": 0.001},
        },
    ],
}

# 3. Reference-only position
REFERENCE_ONLY_POSITION = {
    "chromosome": "chr3",
    "position": 100000,
    "refAllele": "C",
    "altAlleles": ["."],
    "quality": 50.0,
    "filters": [],
    "samples": [
        {
            "genotype": "0/0",
            "totalDepth": 30,
            "genotypeQuality": 60,
        }
    ],
}

# 4. Structural variant
SV_POSITION = {
    "chromosome": "chr7",
    "position": 55000000,
    "refAllele": "N",
    "altAlleles": ["<DEL>"],
    "quality": 999.0,
    "filters": ["PASS"],
    "svEnd": 55050000,
    "svLength": -50000,
    "ciPos": [-50, 50],
    "ciEnd": [-100, 100],
    "samples": [
        {
            "genotype": "0/1",
            "totalDepth": 25,
            "pairedEndReadCounts": [15, 10],
            "splitReadCounts": [12, 8],
        }
    ],
    "variants": [
        {
            "vid": "7-55000000-N-<DEL>",
            "chromosome": "chr7",
            "begin": 55000000,
            "end": 55050000,
            "refAllele": "N",
            "altAllele": "<DEL>",
            "variantType": "deletion",
            "isStructuralVariant": True,
        }
    ],
}

# 5. Missing quality and no filter
MISSING_QUAL_POSITION = {
    "chromosome": "chrX",
    "position": 500000,
    "refAllele": "GCA",
    "altAlleles": ["G"],
    "samples": [
        {
            "genotype": "1/1",
            "totalDepth": 10,
        }
    ],
    "variants": [
        {
            "vid": "X-500000-GCA-G",
            "chromosome": "chrX",
            "begin": 500000,
            "end": 500002,
            "refAllele": "GCA",
            "altAllele": "G",
            "variantType": "deletion",
        }
    ],
}

# 6. Two samples
TWO_SAMPLE_POSITION = {
    "chromosome": "chr1",
    "position": 99999,
    "refAllele": "C",
    "altAlleles": ["T"],
    "quality": 300.0,
    "filters": ["PASS"],
    "samples": [
        {
            "genotype": "0/1",
            "totalDepth": 50,
            "genotypeQuality": 99,
            "variantFrequencies": [0.45],
            "alleleDepths": [28, 22],
        },
        {
            "genotype": "0/0",
            "totalDepth": 60,
            "genotypeQuality": 99,
            "variantFrequencies": [0.0],
            "alleleDepths": [60, 0],
        },
    ],
    "variants": [
        {
            "vid": "1-99999-C-T",
            "chromosome": "chr1",
            "begin": 99999,
            "end": 99999,
            "refAllele": "C",
            "altAllele": "T",
            "variantType": "SNV",
            "dbsnp": ["rs55555"],
            "gnomad": {"allAf": 0.1},
        }
    ],
}

# 7. Empty sample
EMPTY_SAMPLE_POSITION = {
    "chromosome": "chr5",
    "position": 12000,
    "refAllele": "A",
    "altAlleles": ["G"],
    "quality": 100.0,
    "filters": ["PASS"],
    "samples": [
        {
            "isEmpty": True,
            "genotype": "./.",
        }
    ],
    "variants": [
        {
            "vid": "5-12000-A-G",
            "chromosome": "chr5",
            "begin": 12000,
            "end": 12000,
            "refAllele": "A",
            "altAllele": "G",
            "variantType": "SNV",
        }
    ],
}

# 8. Position with failed filters
FAILED_FILTER_POSITION = {
    "chromosome": "chr10",
    "position": 77777,
    "refAllele": "T",
    "altAlleles": ["C"],
    "quality": 15.0,
    "filters": ["LowQual", "LowDP"],
    "samples": [
        {
            "genotype": "0/1",
            "totalDepth": 5,
        }
    ],
    "variants": [
        {
            "vid": "10-77777-T-C",
            "chromosome": "chr10",
            "begin": 77777,
            "end": 77777,
            "refAllele": "T",
            "altAllele": "C",
            "variantType": "SNV",
        }
    ],
}

# 9. Position with SpliceAI annotations
SPLICE_AI_POSITION = {
    "chromosome": "chr11",
    "position": 50000,
    "refAllele": "C",
    "altAlleles": ["T"],
    "quality": 250.0,
    "filters": ["PASS"],
    "samples": [
        {
            "genotype": "0/1",
            "totalDepth": 35,
        }
    ],
    "variants": [
        {
            "vid": "11-50000-C-T",
            "chromosome": "chr11",
            "begin": 50000,
            "end": 50000,
            "refAllele": "C",
            "altAllele": "T",
            "variantType": "SNV",
            "spliceAI": [
                {
                    "hgnc": "GENE2",
                    "acceptorGainScore": 0.1,
                    "acceptorGainDistance": -5,
                    "acceptorLossScore": 0.0,
                    "acceptorLossDistance": 10,
                    "donorGainScore": 0.8,
                    "donorGainDistance": -2,
                    "donorLossScore": 0.3,
                    "donorLossDistance": 15,
                }
            ],
        }
    ],
}


# 10. Insertion with extra context base (normalization test)
INSERTION_EXTRA_CONTEXT = {
    "chromosome": "chrY",
    "position": 2783262,
    "refAllele": "CA",
    "altAlleles": ["CAA"],
    "quality": 100.0,
    "filters": ["PASS"],
    "samples": [
        {
            "genotype": "0/1",
            "totalDepth": 30,
        }
    ],
    "variants": [
        {
            "vid": "Y-2783262-CA-CAA",
            "chromosome": "chrY",
            "begin": 2783262,
            "end": 2783263,
            "refAllele": "CA",
            "altAllele": "CAA",
            "variantType": "insertion",
            "gnomad": {"allAf": 0.001},
            "transcripts": [
                {
                    "transcript": "NM_TEST.1",
                    "source": "RefSeq",
                    "bioType": "protein_coding",
                    "geneId": "TEST1",
                    "hgnc": "TESTGENE",
                    "consequence": ["intron_variant"],
                    "isCanonical": True,
                }
            ],
        }
    ],
}

# 11. SNV embedded in long allele context (normalization test)
SNV_LONG_CONTEXT = {
    "chromosome": "chr1",
    "position": 10800,
    "refAllele": "ACACATGCTAGCGCGTCGGGGTG",
    "altAlleles": ["TCACATGCTAGCGCGTCGGGGTG"],
    "quality": 200.0,
    "filters": ["PASS"],
    "samples": [
        {
            "genotype": "0/1",
            "totalDepth": 40,
        }
    ],
    "variants": [
        {
            "vid": "1-10800-A-T",
            "chromosome": "chr1",
            "begin": 10800,
            "end": 10823,
            "refAllele": "ACACATGCTAGCGCGTCGGGGTG",
            "altAllele": "TCACATGCTAGCGCGTCGGGGTG",
            "variantType": "SNV",
            "phylopScore": 2.1,
            "transcripts": [
                {
                    "transcript": "NM_TEST.2",
                    "source": "RefSeq",
                    "bioType": "protein_coding",
                    "geneId": "TEST2",
                    "hgnc": "TESTGENE2",
                    "consequence": ["missense_variant"],
                    "isCanonical": True,
                }
            ],
        }
    ],
}

# 12. Deletion with extra trailing context (normalization test)
DELETION_EXTRA_SUFFIX = {
    "chromosome": "chr3",
    "position": 50000,
    "refAllele": "CATG",
    "altAlleles": ["CG"],
    "quality": 150.0,
    "filters": ["PASS"],
    "samples": [
        {
            "genotype": "1/1",
            "totalDepth": 20,
        }
    ],
    "variants": [
        {
            "vid": "3-50000-CATG-CG",
            "chromosome": "chr3",
            "begin": 50000,
            "end": 50004,
            "refAllele": "CATG",
            "altAllele": "CG",
            "variantType": "deletion",
            "gnomad": {"allAf": 0.01},
        }
    ],
}

# 13. Three-ALT multi-allelic (decomposition test)
THREE_ALT_POSITION = {
    "chromosome": "chr1",
    "position": 100,
    "refAllele": "A",
    "altAlleles": ["T", "C", "G"],
    "quality": 300.0,
    "filters": ["PASS"],
    "samples": [
        {
            "genotype": "1/3",
            "variantFrequencies": [0.3, 0.0, 0.5],
            "totalDepth": 60,
            "alleleDepths": [12, 18, 0, 30],
        }
    ],
    "variants": [
        {
            "vid": "1-100-A-T",
            "chromosome": "chr1",
            "begin": 100,
            "end": 100,
            "refAllele": "A",
            "altAllele": "T",
            "variantType": "SNV",
            "dbsnp": ["rs111"],
            "gnomad": {"allAf": 0.1},
            "transcripts": [
                {
                    "transcript": "NM_T1.1",
                    "source": "RefSeq",
                    "bioType": "protein_coding",
                    "geneId": "T1",
                    "hgnc": "GENE_T",
                    "consequence": ["missense_variant"],
                    "isCanonical": True,
                }
            ],
        },
        {
            "vid": "1-100-A-C",
            "chromosome": "chr1",
            "begin": 100,
            "end": 100,
            "refAllele": "A",
            "altAllele": "C",
            "variantType": "SNV",
            "gnomad": {"allAf": 0.02},
            "transcripts": [
                {
                    "transcript": "NM_T1.1",
                    "source": "RefSeq",
                    "bioType": "protein_coding",
                    "geneId": "T1",
                    "hgnc": "GENE_T",
                    "consequence": ["synonymous_variant"],
                    "isCanonical": True,
                }
            ],
        },
        {
            "vid": "1-100-A-G",
            "chromosome": "chr1",
            "begin": 100,
            "end": 100,
            "refAllele": "A",
            "altAllele": "G",
            "variantType": "SNV",
            "dbsnp": ["rs333"],
            "gnomad": {"allAf": 0.001},
            "transcripts": [
                {
                    "transcript": "NM_T1.1",
                    "source": "RefSeq",
                    "bioType": "protein_coding",
                    "geneId": "T1",
                    "hgnc": "GENE_T",
                    "consequence": ["stop_gained"],
                    "isCanonical": True,
                }
            ],
        },
    ],
}


# ============================================================
# Helpers
# ============================================================


def make_test_header(h=None):
    """Build a NirvanaHeader from a header dict (defaults to MINIMAL_HEADER)."""
    from nirvana2vcf.models import NirvanaHeader

    h = h or MINIMAL_HEADER
    return NirvanaHeader(
        annotator=h["annotator"],
        creation_time=h["creationTime"],
        genome_assembly=h["genomeAssembly"],
        schema_version=h["schemaVersion"],
        data_version=h["dataVersion"],
        data_sources=h["dataSources"],
        samples=h["samples"],
    )


def parse_pos(position_dict):
    """Parse a position dict into a Position object (test convenience)."""
    from nirvana2vcf.parser import parse_position_line

    return parse_position_line(json.dumps(position_dict))


def build_nirvana_json_string(header_dict, position_dicts, gene_dicts=None):
    """Build a complete Nirvana JSON string in the line-based streaming format.

    Format:
      Line 1: {"header":{...},"positions":[
      Lines 2-N: one position JSON object per line, trailing comma except last
      End: ],"genes":[...]}
    """
    if gene_dicts is None:
        gene_dicts = []

    # Line 1: header + open positions array
    header_json = json.dumps({"header": header_dict})
    # Remove closing }, append ,"positions":[
    line1 = header_json[:-1] + ',"positions":[\n'

    lines = [line1]

    for i, pos in enumerate(position_dicts):
        pos_str = json.dumps(pos)
        if i < len(position_dicts) - 1:
            pos_str += ","
        lines.append(pos_str + "\n")

    # Close positions, open genes
    lines.append('],"genes":[\n')

    for i, gene in enumerate(gene_dicts):
        gene_str = json.dumps(gene)
        if i < len(gene_dicts) - 1:
            gene_str += ","
        lines.append(gene_str + "\n")

    # Close genes and root
    lines.append("]}\n")

    return "".join(lines)


# ============================================================
# Pytest fixtures
# ============================================================


@pytest.fixture
def minimal_snv_json_string():
    return build_nirvana_json_string(MINIMAL_HEADER, [MINIMAL_POSITION_SNV])


@pytest.fixture
def multi_allelic_json_string():
    return build_nirvana_json_string(MINIMAL_HEADER, [MULTI_ALLELIC_POSITION])


@pytest.fixture
def mixed_positions_json_string():
    return build_nirvana_json_string(
        MINIMAL_HEADER,
        [MINIMAL_POSITION_SNV, MULTI_ALLELIC_POSITION, SV_POSITION],
    )


@pytest.fixture
def two_sample_json_string():
    return build_nirvana_json_string(
        TWO_SAMPLE_HEADER,
        [TWO_SAMPLE_POSITION],
    )


@pytest.fixture
def minimal_snv_json_gz(tmp_path, minimal_snv_json_string):
    """Write a gzipped Nirvana JSON file."""
    path = tmp_path / "test.json.gz"
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write(minimal_snv_json_string)
    return str(path)


@pytest.fixture
def minimal_snv_json_file(tmp_path, minimal_snv_json_string):
    """Write an uncompressed Nirvana JSON file."""
    path = tmp_path / "test.json"
    with open(path, "w") as f:
        f.write(minimal_snv_json_string)
    return str(path)


@pytest.fixture
def multi_position_json_gz(tmp_path, mixed_positions_json_string):
    """Write a gzipped multi-position Nirvana JSON file."""
    path = tmp_path / "multi.json.gz"
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write(mixed_positions_json_string)
    return str(path)


@pytest.fixture
def all_positions_json_gz(tmp_path):
    """Write a gzipped JSON with all position types."""
    content = build_nirvana_json_string(
        TWO_SAMPLE_HEADER,
        [
            MINIMAL_POSITION_SNV,
            MULTI_ALLELIC_POSITION,
            SV_POSITION,
            MISSING_QUAL_POSITION,
            TWO_SAMPLE_POSITION,
            EMPTY_SAMPLE_POSITION,
            FAILED_FILTER_POSITION,
        ],
    )
    path = tmp_path / "all.json.gz"
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write(content)
    return str(path)
