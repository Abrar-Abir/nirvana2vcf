"""Real-data integration tests against the HiSeq.10000 Nirvana dataset."""

import gzip
import io
import os

import pytest

from nirvana2vcf.cli import main


pytestmark = [pytest.mark.slow, pytest.mark.real_data]


def _parse_vcf_variants(vcf_path):
    """Return set of (CHROM, POS, REF, ALT) tuples from a VCF file."""
    variants = set()
    open_fn = gzip.open if str(vcf_path).endswith(".gz") else open
    with open_fn(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, _id, ref, alt_field = parts[0], parts[1], parts[2], parts[3], parts[4]
            for alt in alt_field.split(","):
                variants.add((chrom, pos, ref, alt))
    return variants


def _count_data_lines(vcf_path):
    count = 0
    with open(vcf_path) as f:
        for line in f:
            if not line.startswith("#"):
                count += 1
    return count


class TestHiSeqRealData:
    def test_hiseq_runs_without_error(self, hiseq_pair, tmp_path):
        _vcf_path, json_path = hiseq_pair
        out_vcf = str(tmp_path / "out.vcf")
        main(["-i", json_path, "-o", out_vcf])
        assert os.path.exists(out_vcf)
        assert os.path.getsize(out_vcf) > 0

    def test_hiseq_output_is_valid_vcf(self, hiseq_pair, tmp_path):
        _vcf_path, json_path = hiseq_pair
        out_vcf = str(tmp_path / "out.vcf")
        main(["-i", json_path, "-o", out_vcf])

        with open(out_vcf) as f:
            lines = f.readlines()

        assert lines[0].startswith("##fileformat=")
        chrom_line = next(l for l in lines if l.startswith("#CHROM"))
        assert chrom_line.startswith("#CHROM\tPOS\tID\tREF\tALT")

        data_lines = [l for l in lines if not l.startswith("#")]
        assert data_lines, "No data lines in output VCF"
        first = data_lines[0].split("\t")
        assert len(first) >= 8
        assert first[1].isdigit()
        assert first[3]  # REF
        assert first[4]  # ALT

    def test_hiseq_position_concordance(self, hiseq_pair, tmp_path):
        vcf_path, json_path = hiseq_pair
        out_vcf = str(tmp_path / "out.vcf")
        main(["-i", json_path, "-o", out_vcf])

        src = _parse_vcf_variants(vcf_path)
        out = _parse_vcf_variants(out_vcf)

        intersection = src & out
        union = src | out
        jaccard = len(intersection) / len(union) if union else 1.0
        assert jaccard >= 0.99, (
            f"Jaccard {jaccard:.4f} < 0.99  "
            f"src={len(src)} out={len(out)} "
            f"intersection={len(intersection)} union={len(union)}"
        )

    def test_hiseq_variant_count_in_range(self, hiseq_pair, tmp_path):
        _vcf_path, json_path = hiseq_pair
        out_vcf = str(tmp_path / "out.vcf")
        main(["-i", json_path, "-o", out_vcf])
        count = _count_data_lines(out_vcf)
        assert 9_000 <= count <= 11_000, f"Unexpected variant count: {count}"

    def test_hiseq_decompose_mode(self, hiseq_pair, tmp_path):
        _vcf_path, json_path = hiseq_pair
        out_normal = str(tmp_path / "normal.vcf")
        out_decomp = str(tmp_path / "decomp.vcf")
        main(["-i", json_path, "-o", out_normal])
        main(["-i", json_path, "-o", out_decomp, "--decompose"])
        normal_count = _count_data_lines(out_normal)
        decomp_count = _count_data_lines(out_decomp)
        assert decomp_count >= normal_count, (
            f"Decompose produced fewer rows: {decomp_count} < {normal_count}"
        )
