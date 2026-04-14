"""Tests for the JSON-to-VCF field mapping engine."""

import pytest

from nirvana2vcf.mapper import (
    build_csq_string,
    build_info_field,
    build_sample_columns,
    decompose_position,
    map_position_to_vcf_record,
    normalize_alleles,
    _escape_info_value,
    _remap_genotype,
)
from nirvana2vcf.models import Position
from nirvana2vcf.parser import parse_variant, parse_sample
from tests.conftest import (
    make_test_header as _make_header,
    parse_pos as _parse_pos,
    MINIMAL_POSITION_SNV,
    MULTI_ALLELIC_POSITION,
    REFERENCE_ONLY_POSITION,
    SV_POSITION,
    MISSING_QUAL_POSITION,
    TWO_SAMPLE_POSITION,
    EMPTY_SAMPLE_POSITION,
    FAILED_FILTER_POSITION,
    SPLICE_AI_POSITION,
    INSERTION_EXTRA_CONTEXT,
    SNV_LONG_CONTEXT,
    DELETION_EXTRA_SUFFIX,
    THREE_ALT_POSITION,
)


class TestBasicMapping:
    def test_snv_chrom_pos_ref_alt(self):
        """Minimal SNV maps to correct CHROM/POS/REF/ALT."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["CHROM"] == "chr1"
        assert record["POS"] == 12345
        assert record["REF"] == "A"
        assert record["ALT"] == "T"

    def test_snv_id_from_dbsnp(self):
        """dbsnp array maps to ID column."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["ID"] == "rs12345"

    def test_missing_id_produces_dot(self):
        """Position with no dbsnp outputs ID='.'"""
        pos = _parse_pos(MISSING_QUAL_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["ID"] == "."

    def test_multi_allelic_alt(self):
        """Two alt alleles produce comma-separated ALT."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["ALT"] == "A,GT"

    def test_multi_allelic_id_dedup(self):
        """dbsnp from multiple variants are deduplicated."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        # Only first variant has rs9999, second has none
        assert record["ID"] == "rs9999"


class TestQualFilter:
    def test_quality_present(self):
        """quality=200.0 outputs QUAL=200."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["QUAL"] == "200"

    def test_quality_missing(self):
        """quality=None outputs QUAL='.'"""
        pos = _parse_pos(MISSING_QUAL_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["QUAL"] == "."

    def test_filter_pass(self):
        """filters=["PASS"] outputs FILTER=PASS."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["FILTER"] == "PASS"

    def test_filter_multiple(self):
        """filters=["LowQual","LowDP"] outputs semicolon-joined."""
        pos = _parse_pos(FAILED_FILTER_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        assert record["FILTER"] == "LowQual;LowDP"

    def test_filter_empty(self):
        """filters=[] outputs FILTER='.'"""
        pos = _parse_pos(REFERENCE_ONLY_POSITION)
        record = map_position_to_vcf_record(pos, _make_header())

        # Empty filter list should produce "."
        assert record["FILTER"] == "."


class TestInfoField:
    def test_info_gnomad_af(self):
        """gnomAD allele frequency appears in INFO."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        assert "gnomAD_AF=0.00012" in info

    def test_info_gnomad_ac_an(self):
        """gnomAD AC and AN appear in INFO."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        assert "gnomAD_AC=15" in info
        assert "gnomAD_AN=125000" in info

    def test_info_gnomad_eur_af_from_nfe(self):
        """gnomAD_EUR_AF is sourced from Nirvana's nfeAf (not eurAf).

        Nirvana's gnomad block populates nfeAf; eurAf is a 1000 Genomes
        field that never appears under gnomad. The header description
        says "European non-Finnish" to reflect this.
        """
        position = {
            "chromosome": "chr1",
            "position": 100,
            "refAllele": "A",
            "altAlleles": ["T"],
            "variants": [
                {
                    "vid": "1-100-A-T",
                    "chromosome": "chr1",
                    "begin": 100,
                    "end": 100,
                    "refAllele": "A",
                    "altAllele": "T",
                    "variantType": "SNV",
                    "gnomad": {
                        "allAf": 0.1,
                        "nfeAf": 0.25,
                        "eurAf": 0.99,
                    },
                }
            ],
        }
        pos = _parse_pos(position)
        info = build_info_field(pos)

        assert "gnomAD_EUR_AF=0.25" in info
        assert "gnomAD_EUR_AF=0.99" not in info

    def test_info_phylop(self):
        """phyloP score appears in INFO."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        assert "phyloP=3.5" in info

    def test_info_onekg(self):
        """1000 Genomes AF appears in INFO."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        assert "oneKG_AF=0.0002" in info

    def test_info_multi_allelic_per_alt(self):
        """Per-allele fields output comma-separated values in alt order."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        info = build_info_field(pos)

        assert "gnomAD_AF=0.05,0.001" in info

    def test_info_multi_allelic_phylop(self):
        """Per-allele phyloP with values for both alleles."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        info = build_info_field(pos)

        assert "phyloP=1.2,-0.5" in info

    def test_info_clinvar_significance(self):
        """ClinVar significance appears in INFO."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        assert "CLINVAR_SIG=pathogenic" in info
        assert "CLINVAR_ID=RCV000012345" in info

    def test_info_clinvar_review_status_escaped(self):
        """ClinVar review status with special chars is escaped."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        # "criteria provided, single submitter" has a comma and space
        assert "CLINVAR_REVSTAT=criteria%20provided%2C%20single%20submitter" in info

    def test_info_sv_fields(self):
        """SV position outputs SVEND, SVLEN, CIPOS, CIEND."""
        pos = _parse_pos(SV_POSITION)
        info = build_info_field(pos)

        assert "SVEND=55050000" in info
        assert "SVLEN=-50000" in info
        assert "CIPOS=-50,50" in info
        assert "CIEND=-100,100" in info

    def test_info_cytoband(self):
        """CytoBand appears in INFO."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos)

        assert "CytoBand=1p36.33" in info

    def test_info_no_annotations(self):
        """Position with no annotations outputs '.' INFO."""
        pos = _parse_pos(REFERENCE_ONLY_POSITION)
        info = build_info_field(pos)

        assert info == "."

    def test_info_splice_ai(self):
        """SpliceAI scores appear in INFO."""
        pos = _parse_pos(SPLICE_AI_POSITION)
        info = build_info_field(pos)

        assert "SpliceAI_DG_SCORE=0.8" in info
        assert "SpliceAI_DG_DIST=-2" in info
        assert "SpliceAI_AG_SCORE=0.1" in info

    def test_csq_only_mode(self):
        """csq_only=True outputs only CSQ, no flat INFO fields."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        info = build_info_field(pos, csq_only=True)

        assert info.startswith("CSQ=")
        assert "gnomAD_AF" not in info
        assert "phyloP=" not in info
        assert "CytoBand" not in info


class TestCSQString:
    def test_csq_format(self):
        """CSQ field contains pipe-delimited transcript annotation."""
        variant_data = MINIMAL_POSITION_SNV["variants"][0]
        variant = parse_variant(variant_data)
        csq = build_csq_string(variant)

        assert csq is not None
        parts = csq.split("|")
        assert parts[0] == "T"               # Allele
        assert parts[1] == "missense_variant" # Consequence
        assert parts[2] == "GENE1"           # SYMBOL
        assert parts[3] == "1234"            # Gene
        assert parts[4] == "Transcript"      # Feature_type
        assert parts[5] == "NM_001234.5"     # Feature
        assert parts[6] == "protein_coding"  # BIOTYPE

    def test_csq_canonical(self):
        """Canonical transcript is marked YES."""
        variant_data = MINIMAL_POSITION_SNV["variants"][0]
        variant = parse_variant(variant_data)
        csq = build_csq_string(variant)

        parts = csq.split("|")
        assert parts[16] == "YES"  # CANONICAL

    def test_csq_no_transcripts(self):
        """Variant without transcripts returns None."""
        variant_data = MISSING_QUAL_POSITION["variants"][0]
        variant = parse_variant(variant_data)
        csq = build_csq_string(variant)

        assert csq is None


class TestSampleColumns:
    def test_format_and_sample(self):
        """Sample data maps to FORMAT string and per-sample values."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        format_str, samples = build_sample_columns(pos.samples)

        assert "GT" in format_str
        assert "DP" in format_str
        assert "GQ" in format_str
        assert "AD" in format_str
        assert "VF" in format_str
        assert len(samples) == 1
        # Sample values: 0/1:40:99:20,20:0.5
        parts = samples[0].split(":")
        assert parts[0] == "0/1"   # GT
        assert "40" in parts       # DP
        assert "99" in parts       # GQ

    def test_two_samples(self):
        """Two samples produce two sample columns."""
        pos = _parse_pos(TWO_SAMPLE_POSITION)
        format_str, samples = build_sample_columns(pos.samples)

        assert len(samples) == 2
        # First sample: 0/1, second: 0/0
        assert samples[0].startswith("0/1")
        assert samples[1].startswith("0/0")

    def test_empty_sample(self):
        """isEmpty sample outputs ./. for GT."""
        pos = _parse_pos(EMPTY_SAMPLE_POSITION)
        format_str, samples = build_sample_columns(pos.samples)

        assert len(samples) == 1
        assert samples[0].startswith("./.")

    def test_sv_sample_format(self):
        """SV sample includes SR and PR in FORMAT."""
        pos = _parse_pos(SV_POSITION)
        format_str, samples = build_sample_columns(pos.samples)

        assert "SR" in format_str
        assert "PR" in format_str

    def test_no_samples(self):
        """No samples returns empty strings."""
        format_str, samples = build_sample_columns(None)

        assert format_str == ""
        assert samples == []


class TestEscaping:
    def test_escape_semicolon(self):
        assert _escape_info_value("a;b") == "a%3Bb"

    def test_escape_equals(self):
        assert _escape_info_value("a=b") == "a%3Db"

    def test_escape_space(self):
        assert _escape_info_value("a b") == "a%20b"

    def test_escape_comma(self):
        assert _escape_info_value("a,b") == "a%2Cb"


class TestNormalizeAlleles:
    """Unit tests for the normalize_alleles() function."""

    def test_snv_unchanged(self):
        """SNVs are already minimal — no change."""
        assert normalize_alleles(100, "A", ["T"]) == (100, "A", ["T"])

    def test_already_minimal_insertion(self):
        """Minimal insertion REF=C ALT=CA — no change."""
        assert normalize_alleles(100, "C", ["CA"]) == (100, "C", ["CA"])

    def test_already_minimal_deletion(self):
        """Minimal deletion REF=GCA ALT=G — no change."""
        assert normalize_alleles(500000, "GCA", ["G"]) == (500000, "GCA", ["G"])

    def test_insertion_extra_suffix(self):
        """Real-data example: REF=CA ALT=CAA → REF=C ALT=CA (trim trailing A)."""
        assert normalize_alleles(2783262, "CA", ["CAA"]) == (2783262, "C", ["CA"])

    def test_complex_indel_long_context(self):
        """Real-data example: long shared suffix trimmed to minimal."""
        pos, ref, alts = normalize_alleles(
            10816, "CGGGGTGGAG", ["CAGGGGTGGAG"]
        )
        assert pos == 10816
        assert ref == "C"
        assert alts == ["CA"]

    def test_snv_in_long_context(self):
        """Real-data example: SNV buried in 23-base context."""
        pos, ref, alts = normalize_alleles(
            10800,
            "ACACATGCTAGCGCGTCGGGGTG",
            ["TCACATGCTAGCGCGTCGGGGTG"],
        )
        assert pos == 10800
        assert ref == "A"
        assert alts == ["T"]

    def test_deletion_with_shared_suffix(self):
        """Deletion with shared trailing base: REF=CATG ALT=CG → REF=CAT ALT=C."""
        pos, ref, alts = normalize_alleles(50000, "CATG", ["CG"])
        assert pos == 50000
        assert ref == "CAT"
        assert alts == ["C"]

    def test_multi_allelic_shared_suffix(self):
        """Multi-allelic: shared suffix trimmed across all alts."""
        pos, ref, alts = normalize_alleles(100, "ATG", ["AG", "AATG"])
        # Shared suffix 'G': ATG→AT, AG→A, AATG→AAT
        # Then shared prefix 'A': AT→T, A is len 1 so stop
        # Actually: after right-trim: AT, A, AAT — 'A' has len 1, left-trim stops
        assert pos == 100
        assert ref == "AT"
        assert alts == ["A", "AAT"]

    def test_symbolic_allele_skipped(self):
        """Symbolic alleles like <DEL> are not trimmed."""
        assert normalize_alleles(100, "N", ["<DEL>"]) == (100, "N", ["<DEL>"])

    def test_reference_only_dot_skipped(self):
        """ALT='.' (reference-only) is not trimmed."""
        assert normalize_alleles(100, "A", ["."]) == (100, "A", ["."])

    def test_spanning_deletion_skipped(self):
        """ALT='*' (spanning deletion) is not trimmed."""
        assert normalize_alleles(100, "A", ["*"]) == (100, "A", ["*"])

    def test_mixed_real_and_symbolic(self):
        """Mixed real and symbolic alts: only real alts participate in trim."""
        pos, ref, alts = normalize_alleles(100, "ATG", ["AG", "<DEL>"])
        # Only 'AG' is real. Shared suffix 'G': ATG→AT, AG→A. Left-trim 'A': len 1, stop.
        assert pos == 100
        assert ref == "AT"
        assert alts == ["A", "<DEL>"]

    def test_empty_alts(self):
        """No alts — returns unchanged."""
        assert normalize_alleles(100, "A", []) == (100, "A", [])

    def test_mnv_no_shared(self):
        """MNV with no shared prefix or suffix — unchanged."""
        assert normalize_alleles(100, "AT", ["GC"]) == (100, "AT", ["GC"])

    def test_left_trim_adjusts_pos(self):
        """Left-trimming shared prefix increments POS."""
        # REF=AATC ALT=AAGC — right-trim 'C': AAT/AAG — left-trim 'AA': T/G, POS+2
        pos, ref, alts = normalize_alleles(100, "AATC", ["AAGC"])
        assert pos == 102
        assert ref == "T"
        assert alts == ["G"]


class TestNormalizeAllelesWithReference:
    """Phase 3 reference-based left-shifting (bcftools-norm parity)."""

    @staticmethod
    def _make_fetcher(ref_seq, chrom_name="chr1"):
        """Return a 0-based half-open fetcher backed by an in-memory string."""
        def fetcher(chrom, start, end):
            if chrom != chrom_name:
                return ""
            if start < 0 or end > len(ref_seq):
                return ""
            return ref_seq[start:end]
        return fetcher

    def test_homopolymer_insertion_left_shifts(self):
        # ref: A A A A A A A (1-based positions 1..7)
        # raw: insertion of A at pos 7 (REF=A, ALT=AA) → should shift to pos 1
        fetcher = self._make_fetcher("AAAAAAA")
        pos, ref, alts = normalize_alleles(
            7, "A", ["AA"], chrom="chr1", ref_fetcher=fetcher,
        )
        assert pos == 1
        assert ref == "A"
        assert alts == ["AA"]

    def test_homopolymer_deletion_left_shifts(self):
        # ref: A A A A A A A. Deletion of one A at pos 6 (REF=AA, ALT=A) → pos 1
        fetcher = self._make_fetcher("AAAAAAA")
        pos, ref, alts = normalize_alleles(
            6, "AA", ["A"], chrom="chr1", ref_fetcher=fetcher,
        )
        assert pos == 1
        assert ref == "AA"
        assert alts == ["A"]

    def test_str_dinucleotide_left_shifts(self):
        # ref: G A G A G A (positions 1..6). Insertion of "GA" at pos 6.
        # Right-most VCF: REF=A (pos 6), ALT=AGA. Left-most: REF=G (pos 1), ALT=GAG.
        fetcher = self._make_fetcher("GAGAGA")
        pos, ref, alts = normalize_alleles(
            6, "A", ["AGA"], chrom="chr1", ref_fetcher=fetcher,
        )
        assert pos == 1
        assert ref == "G"
        assert alts == ["GAG"]

    def test_already_left_aligned_no_change(self):
        # ref: T T A C G T. Insertion of G at pos 3 (REF=A, ALT=AG).
        # The inserted base G is not present anywhere in the surrounding
        # left context, so the trailing-base check fails immediately.
        fetcher = self._make_fetcher("TTACGT")
        pos, ref, alts = normalize_alleles(
            3, "A", ["AG"], chrom="chr1", ref_fetcher=fetcher,
        )
        assert pos == 3
        assert ref == "A"
        assert alts == ["AG"]

    def test_snv_unchanged(self):
        # SNVs cannot be shifted: the trailing-base check fails immediately.
        fetcher = self._make_fetcher("ACGTACGT")
        pos, ref, alts = normalize_alleles(
            5, "A", ["G"], chrom="chr1", ref_fetcher=fetcher,
        )
        assert pos == 5
        assert ref == "A"
        assert alts == ["G"]

    def test_left_edge_no_underflow(self):
        # Insertion at pos 1 in homopolymer — already at left edge, must not crash.
        fetcher = self._make_fetcher("AAAAA")
        pos, ref, alts = normalize_alleles(
            1, "A", ["AA"], chrom="chr1", ref_fetcher=fetcher,
        )
        assert pos == 1
        assert ref == "A"
        assert alts == ["AA"]

    def test_symbolic_allele_skipped(self):
        # Symbolic ALTs are excluded from real_indices → Phase 3 has nothing to do.
        fetcher = self._make_fetcher("AAAAA")
        pos, ref, alts = normalize_alleles(
            3, "N", ["<DEL>"], chrom="chr1", ref_fetcher=fetcher,
        )
        assert pos == 3
        assert ref == "N"
        assert alts == ["<DEL>"]

    def test_off_chromosome_fetch_stops_loop(self):
        # If the fetcher returns "" (e.g. unknown contig), Phase 3 stops cleanly.
        def fetcher(chrom, start, end):
            return ""
        pos, ref, alts = normalize_alleles(
            100, "A", ["AA"], chrom="chrUnknown", ref_fetcher=fetcher,
        )
        assert pos == 100
        assert ref == "A"
        assert alts == ["AA"]

    def test_no_ref_fetcher_falls_back_to_trim_only(self):
        # Without ref_fetcher, behavior matches the prior trim-only path.
        pos, ref, alts = normalize_alleles(7, "A", ["AA"])
        assert pos == 7
        assert ref == "A"
        assert alts == ["AA"]

    def test_multi_allelic_left_shift(self):
        # ref: A A A A A. Two ALTs that both share trailing A with REF.
        # REF=A (pos 5), ALTs=[AA, AAA] → should shift while both trailing A's match.
        fetcher = self._make_fetcher("AAAAA")
        pos, ref, alts = normalize_alleles(
            5, "A", ["AA", "AAA"], chrom="chr1", ref_fetcher=fetcher,
        )
        assert pos == 1
        assert ref == "A"
        assert alts == ["AA", "AAA"]


class TestNormalizationInRecord:
    """Tests for normalization applied through map_position_to_vcf_record()."""

    def test_insertion_normalized_ref_alt(self):
        """Insertion with extra context: REF/ALT trimmed in record."""
        pos = _parse_pos(INSERTION_EXTRA_CONTEXT)
        record = map_position_to_vcf_record(pos, _make_header(), normalize=True)

        assert record["REF"] == "C"
        assert record["ALT"] == "CA"
        assert record["POS"] == 2783262

    def test_snv_long_context_normalized(self):
        """SNV in long context: trimmed to single-base REF/ALT."""
        pos = _parse_pos(SNV_LONG_CONTEXT)
        record = map_position_to_vcf_record(pos, _make_header(), normalize=True)

        assert record["REF"] == "A"
        assert record["ALT"] == "T"
        assert record["POS"] == 10800

    def test_deletion_extra_suffix_normalized(self):
        """Deletion with extra trailing base: suffix trimmed."""
        pos = _parse_pos(DELETION_EXTRA_SUFFIX)
        record = map_position_to_vcf_record(pos, _make_header(), normalize=True)

        assert record["REF"] == "CAT"
        assert record["ALT"] == "C"
        assert record["POS"] == 50000

    def test_no_normalize_preserves_original(self):
        """normalize=False preserves raw Nirvana alleles."""
        pos = _parse_pos(INSERTION_EXTRA_CONTEXT)
        record = map_position_to_vcf_record(pos, _make_header(), normalize=False)

        assert record["REF"] == "CA"
        assert record["ALT"] == "CAA"

    def test_snv_unchanged_with_normalize(self):
        """Already-minimal SNV is unchanged by normalization."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        record_norm = map_position_to_vcf_record(pos, _make_header(), normalize=True)
        record_raw = map_position_to_vcf_record(pos, _make_header(), normalize=False)

        assert record_norm["REF"] == record_raw["REF"]
        assert record_norm["ALT"] == record_raw["ALT"]
        assert record_norm["POS"] == record_raw["POS"]

    def test_normalize_csq_allele_updated(self):
        """CSQ Allele field uses normalized alt allele."""
        pos = _parse_pos(INSERTION_EXTRA_CONTEXT)
        record = map_position_to_vcf_record(pos, _make_header(), normalize=True)

        # CSQ starts with the Allele field (first pipe-delimited value)
        csq = record["INFO"].split("CSQ=")[1]
        csq_allele = csq.split("|")[0]
        assert csq_allele == "CA"  # Normalized from "CAA"

    def test_normalize_csq_allele_snv_context(self):
        """CSQ Allele field uses normalized alt for SNV-in-context."""
        pos = _parse_pos(SNV_LONG_CONTEXT)
        record = map_position_to_vcf_record(pos, _make_header(), normalize=True)

        csq = record["INFO"].split("CSQ=")[1]
        csq_allele = csq.split("|")[0]
        assert csq_allele == "T"  # Normalized from "TCACATGCTAGCGCGTCGGGGTG"

    def test_normalize_per_allele_info_preserved(self):
        """Per-allele INFO values remain correct after normalization."""
        pos = _parse_pos(INSERTION_EXTRA_CONTEXT)
        record = map_position_to_vcf_record(pos, _make_header(), normalize=True)

        assert "gnomAD_AF=0.001" in record["INFO"]


class TestRemapGenotype:
    """Unit tests for _remap_genotype()."""

    def test_het_first_alt(self):
        assert _remap_genotype("0/1", 1) == "0/1"

    def test_het_second_alt(self):
        assert _remap_genotype("0/2", 2) == "0/1"

    def test_het_second_alt_from_first(self):
        """0/2 seen from ALT1's row → 0/."""
        assert _remap_genotype("0/2", 1) == "0/."

    def test_compound_het_first(self):
        assert _remap_genotype("1/2", 1) == "1/."

    def test_compound_het_second(self):
        assert _remap_genotype("1/2", 2) == "./1"

    def test_hom_first_alt(self):
        assert _remap_genotype("1/1", 1) == "1/1"

    def test_hom_second_alt_becomes_missing(self):
        assert _remap_genotype("2/2", 1) == "./."

    def test_hom_second_alt(self):
        assert _remap_genotype("2/2", 2) == "1/1"

    def test_ref_ref_unchanged(self):
        assert _remap_genotype("0/0", 1) == "0/0"

    def test_missing_unchanged(self):
        assert _remap_genotype("./.", 1) == "./."

    def test_phased_preserved(self):
        assert _remap_genotype("0|1", 1) == "0|1"

    def test_phased_compound_het(self):
        assert _remap_genotype("1|2", 2) == ".|1"

    def test_three_alts_first(self):
        assert _remap_genotype("1/3", 1) == "1/."

    def test_three_alts_third(self):
        assert _remap_genotype("1/3", 3) == "./1"

    def test_three_alts_middle(self):
        assert _remap_genotype("2/3", 2) == "1/."

    def test_haploid_match(self):
        assert _remap_genotype("1", 1) == "1"

    def test_haploid_no_match(self):
        assert _remap_genotype("2", 1) == "."

    def test_haploid_ref(self):
        assert _remap_genotype("0", 1) == "0"


class TestDecomposePosition:
    """Unit tests for decompose_position()."""

    def test_biallelic_passthrough(self):
        """Single-ALT position returns unchanged."""
        pos = _parse_pos(MINIMAL_POSITION_SNV)
        result = decompose_position(pos)
        assert len(result) == 1
        assert result[0] is pos

    def test_multi_allelic_produces_two(self):
        """Two-ALT position decomposes into two Positions."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        result = decompose_position(pos)
        assert len(result) == 2

    def test_decomposed_alt_alleles(self):
        """Each decomposed Position has a single ALT."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        result = decompose_position(pos)
        assert result[0].alt_alleles == ["A"]
        assert result[1].alt_alleles == ["GT"]

    def test_decomposed_variants_scoped(self):
        """Each decomposed Position has only the matching variant."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        result = decompose_position(pos)
        assert len(result[0].variants) == 1
        assert result[0].variants[0].alt_allele == "A"
        assert len(result[1].variants) == 1
        assert result[1].variants[0].alt_allele == "GT"

    def test_decomposed_gt_remapped(self):
        """Genotype is remapped for each decomposed Position."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        result = decompose_position(pos)
        # Original GT: 1/2
        assert result[0].samples[0].genotype == "1/."
        assert result[1].samples[0].genotype == "./1"

    def test_decomposed_allele_depths(self):
        """AD is sliced to [ref, target_alt] for each Position."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        result = decompose_position(pos)
        # Original AD: [10, 19, 28]
        assert result[0].samples[0].allele_depths == [10, 19]
        assert result[1].samples[0].allele_depths == [10, 28]

    def test_decomposed_variant_frequencies(self):
        """VF is sliced to [target_alt_vf] for each Position."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        result = decompose_position(pos)
        # Original VF: [0.333, 0.5]
        assert result[0].samples[0].variant_frequencies == [0.333]
        assert result[1].samples[0].variant_frequencies == [0.5]

    def test_decomposed_position_level_fields_copied(self):
        """Position-level fields are duplicated to each decomposed row."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        result = decompose_position(pos)
        for r in result:
            assert r.chromosome == "chr2"
            assert r.position == 48010488
            assert r.ref_allele == "G"
            assert r.quality == 150.0
            assert r.filters == ["PASS"]

    def test_three_alts_produces_three(self):
        """Three-ALT position decomposes into three Positions."""
        pos = _parse_pos(THREE_ALT_POSITION)
        result = decompose_position(pos)
        assert len(result) == 3
        assert [r.alt_alleles[0] for r in result] == ["T", "C", "G"]

    def test_three_alts_gt_remapped(self):
        """Three-ALT genotype remapped correctly."""
        pos = _parse_pos(THREE_ALT_POSITION)
        result = decompose_position(pos)
        # Original GT: 1/3
        assert result[0].samples[0].genotype == "1/."   # ALT=T (target=1)
        assert result[1].samples[0].genotype == "./."   # ALT=C (target=2)
        assert result[2].samples[0].genotype == "./1"   # ALT=G (target=3)

    def test_three_alts_ad_sliced(self):
        """Three-ALT AD sliced to [ref, target_alt]."""
        pos = _parse_pos(THREE_ALT_POSITION)
        result = decompose_position(pos)
        # Original AD: [12, 18, 0, 30]
        assert result[0].samples[0].allele_depths == [12, 18]
        assert result[1].samples[0].allele_depths == [12, 0]
        assert result[2].samples[0].allele_depths == [12, 30]

    def test_three_alts_vf_sliced(self):
        """Three-ALT VF sliced correctly."""
        pos = _parse_pos(THREE_ALT_POSITION)
        result = decompose_position(pos)
        # Original VF: [0.3, 0.0, 0.5]
        assert result[0].samples[0].variant_frequencies == [0.3]
        assert result[1].samples[0].variant_frequencies == [0.0]
        assert result[2].samples[0].variant_frequencies == [0.5]

    def test_empty_alts_passthrough(self):
        """Position with no ALTs returns unchanged."""
        pos = Position(chromosome="chr1", position=100, ref_allele="A", alt_alleles=[])
        result = decompose_position(pos)
        assert len(result) == 1

    def test_reference_only_passthrough(self):
        """Reference-only position (ALT='.') returns unchanged."""
        pos = _parse_pos(REFERENCE_ONLY_POSITION)
        result = decompose_position(pos)
        assert len(result) == 1

    def test_no_samples_decomposed(self):
        """Position without samples decomposes correctly."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        # Replace samples with None via a new Position
        pos_no_samples = Position(
            chromosome=pos.chromosome, position=pos.position,
            ref_allele=pos.ref_allele, alt_alleles=pos.alt_alleles,
            variants=pos.variants, samples=None,
        )
        result = decompose_position(pos_no_samples)
        assert len(result) == 2
        assert all(r.samples is None for r in result)

    def test_decomposed_total_depth_preserved(self):
        """Scalar sample fields like DP are copied unchanged."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        result = decompose_position(pos)
        for r in result:
            assert r.samples[0].total_depth == 57

    def test_decomposed_info_single_values(self):
        """Per-allele INFO fields have single values after decomposition + mapping."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        result = decompose_position(pos)

        rec0 = map_position_to_vcf_record(result[0], _make_header())
        rec1 = map_position_to_vcf_record(result[1], _make_header())

        assert "gnomAD_AF=0.05" in rec0["INFO"]
        assert "gnomAD_AF=0.001" in rec1["INFO"]
        # No commas in gnomAD_AF
        for rec in [rec0, rec1]:
            for part in rec["INFO"].split(";"):
                if part.startswith("gnomAD_AF="):
                    assert "," not in part

    def test_decomposed_dbsnp_scoped(self):
        """dbSNP IDs scoped to matching variant after decomposition."""
        pos = _parse_pos(MULTI_ALLELIC_POSITION)
        result = decompose_position(pos)

        rec0 = map_position_to_vcf_record(result[0], _make_header())
        rec1 = map_position_to_vcf_record(result[1], _make_header())

        assert rec0["ID"] == "rs9999"  # ALT=A variant has rs9999
        assert rec1["ID"] == "."       # ALT=GT variant has no dbsnp
