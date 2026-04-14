"""Streaming parser for Nirvana JSON format.

Nirvana JSON uses a line-based streaming format:
  Line 1:    {"header":{...},"positions":[
  Lines 2-N: one position JSON object per line (trailing comma except last)
  End:       ],"genes":[...]}
"""

import gzip
from typing import Generator, Iterator, List, Optional, Tuple

from orjson import loads as _json_loads

from .models import (
    ClinVarEntry,
    NirvanaHeader,
    PopulationFrequency,
    Position,
    Sample,
    SpliceAIEntry,
    TranscriptAnnotation,
    Variant,
)


def open_nirvana_file(path: str) -> Iterator[str]:
    """Open a .json.gz or .json file and yield lines."""
    if path.endswith(".gz"):
        f = gzip.open(path, "rt", encoding="utf-8")
    else:
        f = open(path, "r", encoding="utf-8")
    try:
        yield from f
    finally:
        f.close()


def parse_header(first_line: str) -> NirvanaHeader:
    """Extract NirvanaHeader from the first line of the file.

    The first line looks like: {"header":{...},"positions":[
    We need to extract the header dict by finding the ,"positions":[ suffix.
    """
    line = first_line.strip()

    # Find the ,"positions":[ at the end and remove it
    positions_marker = ',"positions":['
    idx = line.rfind(positions_marker)
    if idx == -1:
        raise ValueError(
            "First line does not contain '\"positions\":[' marker. "
            "Is this a valid Nirvana JSON file?"
        )

    # Take everything before the marker, close the JSON object
    header_json = line[:idx] + "}"
    data = _json_loads(header_json)

    h = data.get("header", data)
    return NirvanaHeader(
        annotator=h.get("annotator", ""),
        creation_time=h.get("creationTime", ""),
        genome_assembly=h.get("genomeAssembly", ""),
        schema_version=h.get("schemaVersion", 0),
        data_version=h.get("dataVersion", ""),
        data_sources=h.get("dataSources", []),
        samples=h.get("samples", []),
    )


def _parse_population_freq(data: Optional[dict]) -> Optional[PopulationFrequency]:
    """Parse a population frequency dict into PopulationFrequency."""
    if not data:
        return None
    return PopulationFrequency(
        all_af=data.get("allAf"),
        all_ac=data.get("allAc"),
        all_an=data.get("allAn"),
        all_hc=data.get("allHc"),
        afr_af=data.get("afrAf"),
        amr_af=data.get("amrAf"),
        eur_af=data.get("eurAf"),
        eas_af=data.get("easAf"),
        sas_af=data.get("sasAf"),
        fin_af=data.get("finAf"),
        nfe_af=data.get("nfeAf"),
        asj_af=data.get("asjAf"),
        oth_af=data.get("othAf"),
        male_af=data.get("maleAf"),
        female_af=data.get("femaleAf"),
        failed_filter=data.get("failedFilter"),
    )


def _parse_splice_ai(data_list: Optional[list]) -> Optional[List[SpliceAIEntry]]:
    """Parse spliceAI array into list of SpliceAIEntry."""
    if not data_list:
        return None
    result = []
    for d in data_list:
        result.append(
            SpliceAIEntry(
                hgnc=d.get("hgnc"),
                acceptor_gain_score=d.get("acceptorGainScore"),
                acceptor_gain_distance=d.get("acceptorGainDistance"),
                acceptor_loss_score=d.get("acceptorLossScore"),
                acceptor_loss_distance=d.get("acceptorLossDistance"),
                donor_gain_score=d.get("donorGainScore"),
                donor_gain_distance=d.get("donorGainDistance"),
                donor_loss_score=d.get("donorLossScore"),
                donor_loss_distance=d.get("donorLossDistance"),
            )
        )
    return result


def _parse_clinvar(data_list: Optional[list]) -> Optional[List[ClinVarEntry]]:
    """Parse clinvar array into list of ClinVarEntry."""
    if not data_list:
        return None
    result = []
    for d in data_list:
        result.append(
            ClinVarEntry(
                id=d.get("id", ""),
                variation_id=d.get("variationId"),
                review_status=d.get("reviewStatus"),
                significance=d.get("significance"),
                allele_origins=d.get("alleleOrigins"),
                phenotypes=d.get("phenotypes"),
                med_gen_ids=d.get("medGenIds"),
                omim_ids=d.get("omimIds"),
                orphanet_ids=d.get("orphanetIds"),
                pub_med_ids=d.get("pubMedIds"),
                is_allele_specific=d.get("isAlleleSpecific"),
                last_updated_date=d.get("lastUpdatedDate"),
                ref_allele=d.get("refAllele"),
                alt_allele=d.get("altAllele"),
            )
        )
    return result


def _parse_transcript(data: dict) -> TranscriptAnnotation:
    """Parse a single transcript dict into TranscriptAnnotation."""
    return TranscriptAnnotation(
        transcript=data.get("transcript", ""),
        source=data.get("source", ""),
        bio_type=data.get("bioType"),
        gene_id=data.get("geneId"),
        hgnc=data.get("hgnc"),
        protein_id=data.get("proteinId"),
        consequence=data.get("consequence"),
        codons=data.get("codons"),
        amino_acids=data.get("aminoAcids"),
        cdna_pos=data.get("cdnaPos"),
        cds_pos=data.get("cdsPos"),
        protein_pos=data.get("proteinPos"),
        exons=data.get("exons"),
        introns=data.get("introns"),
        hgvsc=data.get("hgvsc"),
        hgvsp=data.get("hgvsp"),
        is_canonical=data.get("isCanonical"),
        poly_phen_score=data.get("polyPhenScore"),
        poly_phen_prediction=data.get("polyPhenPrediction"),
        sift_score=data.get("siftScore"),
        sift_prediction=data.get("siftPrediction"),
        complete_overlap=data.get("completeOverlap"),
    )


def parse_variant(data: dict) -> Variant:
    """Parse a single variant dict into Variant."""
    transcripts = None
    if "transcripts" in data:
        transcripts = [_parse_transcript(t) for t in data["transcripts"]]

    # Parse revel - can be a dict with "score" or just a float
    revel_score = None
    revel_data = data.get("revel")
    if isinstance(revel_data, dict):
        revel_score = revel_data.get("score")
    elif isinstance(revel_data, (int, float)):
        revel_score = float(revel_data)

    return Variant(
        vid=data.get("vid", ""),
        chromosome=data.get("chromosome", ""),
        begin=data.get("begin", 0),
        end=data.get("end", 0),
        ref_allele=data.get("refAllele", ""),
        alt_allele=data.get("altAllele", ""),
        variant_type=data.get("variantType", ""),
        hgvsg=data.get("hgvsg"),
        phylop_score=data.get("phylopScore"),
        dann_score=data.get("dannScore"),
        gerp_score=data.get("gerpScore"),
        is_reference_minor_allele=data.get("isReferenceMinorAllele"),
        is_structural_variant=data.get("isStructuralVariant"),
        in_low_complexity_region=data.get("inLowComplexityRegion"),
        is_decomposed_variant=data.get("isDecomposedVariant"),
        is_recomposed_variant=data.get("isRecomposedVariant"),
        dbsnp=data.get("dbsnp"),
        transcripts=transcripts,
        regulatory_regions=data.get("regulatoryRegions"),
        clinvar=_parse_clinvar(data.get("clinvar")),
        one_kg=_parse_population_freq(data.get("oneKg")),
        gnomad=_parse_population_freq(data.get("gnomad")),
        topmed=_parse_population_freq(data.get("topmed")),
        gme_variome=_parse_population_freq(data.get("gmeVariome")),
        splice_ai=_parse_splice_ai(data.get("spliceAI")),
        revel_score=revel_score,
        raw=data,
    )


def parse_sample(data: dict) -> Sample:
    """Parse a single sample dict into Sample."""
    return Sample(
        genotype=data.get("genotype", "./."),
        is_empty=data.get("isEmpty"),
        total_depth=data.get("totalDepth"),
        genotype_quality=data.get("genotypeQuality"),
        variant_frequencies=data.get("variantFrequencies"),
        allele_depths=data.get("alleleDepths"),
        copy_number=data.get("copyNumber"),
        minor_haplotype_copy_number=data.get("minorHaplotypeCopyNumber"),
        repeat_unit_counts=data.get("repeatUnitCounts"),
        failed_filter=data.get("failedFilter"),
        is_de_novo=data.get("isDeNovo"),
        de_novo_quality=data.get("deNovoQuality"),
        split_read_counts=data.get("splitReadCounts"),
        paired_end_read_counts=data.get("pairedEndReadCounts"),
        heteroplasmypercentile=data.get("heteroplasmyPercentile"),
        loss_of_heterozygosity=data.get("lossOfHeterozygosity"),
        somatic_quality=data.get("somaticQuality"),
        artifact_adjusted_quality_score=data.get("artifactAdjustedQualityScore"),
        likelihood_ratio_quality_score=data.get("likelihoodRatioQualityScore"),
        bin_count=data.get("binCount"),
    )


def parse_position_line(line: str) -> Optional[Position]:
    """Parse a single position line into a Position object.

    Returns None if the line marks the end of positions (genes section start).
    Strips trailing comma before parsing.
    """
    stripped = line.strip()

    # End of positions section
    if not stripped or stripped.startswith("]"):
        return None

    # Strip trailing comma
    if stripped.endswith(","):
        stripped = stripped[:-1]

    data = _json_loads(stripped)

    variants = None
    if "variants" in data:
        variants = [parse_variant(v) for v in data["variants"]]

    samples = None
    if "samples" in data:
        samples = [parse_sample(s) for s in data["samples"]]

    return Position(
        chromosome=data.get("chromosome", ""),
        position=data.get("position", 0),
        ref_allele=data.get("refAllele", ""),
        alt_alleles=data.get("altAlleles", []),
        quality=data.get("quality"),
        filters=data.get("filters"),
        cytogenetic_band=data.get("cytogeneticBand"),
        repeat_unit=data.get("repeatUnit"),
        ref_repeat_count=data.get("refRepeatCount"),
        sv_end=data.get("svEnd"),
        sv_length=data.get("svLength"),
        ci_pos=data.get("ciPos"),
        ci_end=data.get("ciEnd"),
        strand_bias=data.get("strandBias"),
        joint_somatic_normal_quality=data.get("jointSomaticNormalQuality"),
        variants=variants,
        samples=samples,
        clingen=data.get("clingen"),
        clingen_dosage_sensitivity_map=data.get("clingenDosageSensitivityMap"),
        raw=data,
    )


def stream_positions(path: str) -> Generator[Tuple[NirvanaHeader, Position], None, None]:
    """Stream positions from a Nirvana JSON file.

    Yields (header, position) tuples. The header is the same object for
    every position — it is yielded with each position for convenience.
    """
    lines = open_nirvana_file(path)
    header = None

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        # First non-empty line is the header line
        if header is None:
            header = parse_header(stripped)
            continue

        # Try to parse as a position
        position = parse_position_line(stripped)
        if position is None:
            # Reached end of positions (genes section or closing)
            break

        yield header, position
