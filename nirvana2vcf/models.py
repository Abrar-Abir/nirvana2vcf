from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any


@dataclass
class NirvanaHeader:
    annotator: str = ""
    creation_time: str = ""
    genome_assembly: str = ""
    schema_version: int = 0
    data_version: str = ""
    data_sources: List[Dict[str, str]] = field(default_factory=list)
    samples: List[str] = field(default_factory=list)


@dataclass
class TranscriptAnnotation:
    transcript: str = ""
    source: str = ""  # "RefSeq" or "Ensembl"
    bio_type: Optional[str] = None
    gene_id: Optional[str] = None
    hgnc: Optional[str] = None
    protein_id: Optional[str] = None
    consequence: Optional[List[str]] = None
    codons: Optional[str] = None
    amino_acids: Optional[str] = None
    cdna_pos: Optional[str] = None
    cds_pos: Optional[str] = None
    protein_pos: Optional[str] = None
    exons: Optional[str] = None
    introns: Optional[str] = None
    hgvsc: Optional[str] = None
    hgvsp: Optional[str] = None
    is_canonical: Optional[bool] = None
    poly_phen_score: Optional[float] = None
    poly_phen_prediction: Optional[str] = None
    sift_score: Optional[float] = None
    sift_prediction: Optional[str] = None
    complete_overlap: Optional[bool] = None


@dataclass
class ClinVarEntry:
    id: str = ""
    variation_id: Optional[str] = None
    review_status: Optional[str] = None
    significance: Optional[List[str]] = None
    allele_origins: Optional[List[str]] = None
    phenotypes: Optional[List[str]] = None
    med_gen_ids: Optional[List[str]] = None
    omim_ids: Optional[List[str]] = None
    orphanet_ids: Optional[List[str]] = None
    pub_med_ids: Optional[List[str]] = None
    is_allele_specific: Optional[bool] = None
    last_updated_date: Optional[str] = None
    ref_allele: Optional[str] = None
    alt_allele: Optional[str] = None


@dataclass
class PopulationFrequency:
    all_af: Optional[float] = None
    all_ac: Optional[int] = None
    all_an: Optional[int] = None
    all_hc: Optional[int] = None
    afr_af: Optional[float] = None
    amr_af: Optional[float] = None
    eur_af: Optional[float] = None
    eas_af: Optional[float] = None
    sas_af: Optional[float] = None
    fin_af: Optional[float] = None
    nfe_af: Optional[float] = None
    asj_af: Optional[float] = None
    oth_af: Optional[float] = None
    male_af: Optional[float] = None
    female_af: Optional[float] = None
    failed_filter: Optional[bool] = None


@dataclass
class SpliceAIEntry:
    hgnc: Optional[str] = None
    acceptor_gain_score: Optional[float] = None
    acceptor_gain_distance: Optional[int] = None
    acceptor_loss_score: Optional[float] = None
    acceptor_loss_distance: Optional[int] = None
    donor_gain_score: Optional[float] = None
    donor_gain_distance: Optional[int] = None
    donor_loss_score: Optional[float] = None
    donor_loss_distance: Optional[int] = None


@dataclass
class Variant:
    vid: str = ""
    chromosome: str = ""
    begin: int = 0
    end: int = 0
    ref_allele: str = ""
    alt_allele: str = ""
    variant_type: str = ""
    hgvsg: Optional[str] = None
    phylop_score: Optional[float] = None
    dann_score: Optional[float] = None
    gerp_score: Optional[float] = None
    is_reference_minor_allele: Optional[bool] = None
    is_structural_variant: Optional[bool] = None
    in_low_complexity_region: Optional[bool] = None
    is_decomposed_variant: Optional[bool] = None
    is_recomposed_variant: Optional[bool] = None
    dbsnp: Optional[List[str]] = None
    transcripts: Optional[List[TranscriptAnnotation]] = None
    regulatory_regions: Optional[List[Dict[str, Any]]] = None
    clinvar: Optional[List[ClinVarEntry]] = None
    one_kg: Optional[PopulationFrequency] = None
    gnomad: Optional[PopulationFrequency] = None
    topmed: Optional[PopulationFrequency] = None
    gme_variome: Optional[PopulationFrequency] = None
    splice_ai: Optional[List[SpliceAIEntry]] = None
    revel_score: Optional[float] = None
    raw: Optional[Dict[str, Any]] = None


@dataclass
class Sample:
    genotype: str = "./."
    is_empty: Optional[bool] = None
    total_depth: Optional[int] = None
    genotype_quality: Optional[int] = None
    variant_frequencies: Optional[List[float]] = None
    allele_depths: Optional[List[int]] = None
    copy_number: Optional[int] = None
    minor_haplotype_copy_number: Optional[int] = None
    repeat_unit_counts: Optional[List[int]] = None
    failed_filter: Optional[bool] = None
    is_de_novo: Optional[bool] = None
    de_novo_quality: Optional[float] = None
    split_read_counts: Optional[List[int]] = None
    paired_end_read_counts: Optional[List[int]] = None
    heteroplasmypercentile: Optional[List[float]] = None
    loss_of_heterozygosity: Optional[bool] = None
    somatic_quality: Optional[float] = None
    artifact_adjusted_quality_score: Optional[float] = None
    likelihood_ratio_quality_score: Optional[float] = None
    bin_count: Optional[int] = None


@dataclass
class Position:
    chromosome: str = ""
    position: int = 0
    ref_allele: str = ""
    alt_alleles: List[str] = field(default_factory=list)
    quality: Optional[float] = None
    filters: Optional[List[str]] = None
    cytogenetic_band: Optional[str] = None
    repeat_unit: Optional[str] = None
    ref_repeat_count: Optional[int] = None
    sv_end: Optional[int] = None
    sv_length: Optional[int] = None
    ci_pos: Optional[List[int]] = None
    ci_end: Optional[List[int]] = None
    strand_bias: Optional[float] = None
    joint_somatic_normal_quality: Optional[int] = None
    variants: Optional[List[Variant]] = None
    samples: Optional[List[Sample]] = None
    clingen: Optional[List[Dict[str, Any]]] = None
    clingen_dosage_sensitivity_map: Optional[List[Dict[str, Any]]] = None
    raw: Optional[Dict[str, Any]] = None
