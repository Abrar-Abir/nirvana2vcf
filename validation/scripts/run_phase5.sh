#!/usr/bin/env bash
# =============================================================================
# Phase 5: Cross-Tool Comparison
#
# Compares nirvana2vcf output against annotations produced by independent tools
# (VEP, SnpEff/SnpSift, bcftools annotate).  The goal is characterisation —
# confirming nirvana2vcf output is reasonable and consistent with the broader
# annotation ecosystem — not exact matching (different tools use different
# database versions).
#
# Prerequisites:
#   - Phase 2 complete (default VCF in data/output/phase2/)
#   - Nirvana JSON available (HiSeq.10000.json.gz in data/output/)
#   - Original VCF available (HiSeq.10000.vcf.gz in data/)
#   - python3, bcftools, bgzip, tabix available
#   - Optional: vep, snpEff, SnpSift (sub-phases skip gracefully if absent)
#
# Environment variable overrides:
#   SKIP_VEP=true       Skip VEP sub-phase
#   SKIP_SNPEFF=true    Skip SnpEff sub-phase
#   SKIP_BCFTOOLS=true  Skip bcftools sub-phase
#
# Usage:
#   bash validation/scripts/run_phase5.sh
#
# Idempotent — safe to re-run. Overwrites results.
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Section 0: Configuration
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATA_DIR="$PROJECT_ROOT/validation/data"
OUTPUT_DIR="$DATA_DIR/output"
RESULTS_DIR="$PROJECT_ROOT/validation/results"
STRATEGY_MD="$PROJECT_ROOT/validation/docs/strategy.md"

PHASE5_DATA="$DATA_DIR/phase5"
PHASE5_RESULTS="$RESULTS_DIR/phase5"
REFS_DIR="$PROJECT_ROOT/validation/refs"

# HGNC complete-set path shared by 5a (VEP) and 5b (SnpEff) comparators.
# The file is downloaded in 5a; 5b reuses whatever 5a left on disk, and both
# comparators fall back to exact symbol match if the file is absent.
HGNC_FILE="$REFS_DIR/hgnc_complete_set.txt"

# Input files
J2V_VCF="$OUTPUT_DIR/phase2/HiSeq.10000.default.vcf"
ORIGINAL_VCF="$DATA_DIR/HiSeq.10000.vcf.gz"
NIRVANA_JSON="$OUTPUT_DIR/HiSeq.10000.json.gz"

# Comparison scripts
COMPARE_VEP="$SCRIPT_DIR/compare_vep.py"
COMPARE_SNPEFF="$SCRIPT_DIR/compare_snpeff.py"
COMPARE_GNOMAD="$SCRIPT_DIR/compare_gnomad_dbsnp.py"

# Skip flags (default: false)
SKIP_VEP="${SKIP_VEP:-false}"
SKIP_SNPEFF="${SKIP_SNPEFF:-false}"
SKIP_BCFTOOLS="${SKIP_BCFTOOLS:-false}"

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log()   { echo -e "${BLUE}[Phase5]${NC} $*"; }
ok()    { echo -e "${GREEN}[  OK  ]${NC} $*"; }
warn()  { echo -e "${YELLOW}[ WARN ]${NC} $*"; }
err()   { echo -e "${RED}[ERROR ]${NC} $*" >&2; }

# Track overall results
ERRORS=0
WARNINGS=0

record_pass() { ok "$1"; }
record_warn() { warn "$1"; WARNINGS=$((WARNINGS + 1)); }
record_fail() { err "$1"; ERRORS=$((ERRORS + 1)); }

# Sub-phase results
VEP_RESULT="SKIPPED"
SNPEFF_RESULT="SKIPPED"
BCFTOOLS_RESULT="SKIPPED"

# ---------------------------------------------------------------------------
# Section 1: Prerequisites
# ---------------------------------------------------------------------------

log "Checking prerequisites..."

# nirvana2vcf VCF (Phase 2 default output)
if [ ! -f "$J2V_VCF" ]; then
    err "nirvana2vcf VCF not found: $J2V_VCF"
    err "Run Phase 2 first: bash validation/scripts/run_phase2.sh"
    exit 1
fi
ok "nirvana2vcf VCF found: $J2V_VCF"

# Original VCF
if [ ! -f "$ORIGINAL_VCF" ]; then
    err "Original VCF not found: $ORIGINAL_VCF"
    err "Run Phase 1 first: bash validation/scripts/run_phase1.sh"
    exit 1
fi
ok "Original VCF found: $ORIGINAL_VCF"

# Nirvana JSON
if [ ! -f "$NIRVANA_JSON" ]; then
    err "Nirvana JSON not found: $NIRVANA_JSON"
    exit 1
fi
ok "Nirvana JSON found: $NIRVANA_JSON"

# python3
if ! command -v python3 &>/dev/null; then
    err "python3 not found"
    exit 1
fi
ok "python3 available"

# bcftools
if ! command -v bcftools &>/dev/null; then
    err "bcftools not found (required for contig renaming and gnomAD annotation)"
    exit 1
fi
ok "bcftools available"

# bgzip / tabix
if ! command -v bgzip &>/dev/null || ! command -v tabix &>/dev/null; then
    err "bgzip/tabix not found"
    exit 1
fi
ok "bgzip/tabix available"

# Comparison scripts
for script in "$COMPARE_VEP" "$COMPARE_SNPEFF" "$COMPARE_GNOMAD"; do
    if [ ! -f "$script" ]; then
        err "Comparison script not found: $script"
        exit 1
    fi
done
ok "All comparison scripts found"

# Create working directories
mkdir -p "$PHASE5_DATA" "$PHASE5_RESULTS" "$REFS_DIR"

# --- Helper: text-based chr-to-bare rename (works even without ##contig) ---
strip_chr_prefix() {
    # Rename chr1→1 etc. in both ##contig headers and data lines via awk.
    # Handles VCFs that lack ##contig lines (where bcftools --rename-chrs fails).
    # Accepts both compressed (.gz/.bgz) and uncompressed input.
    local infile="$1" outfile="$2"
    if [[ "$infile" == *.gz || "$infile" == *.bgz ]]; then
        gzip -dc "$infile"
    else
        cat "$infile"
    fi | awk -F'\t' '
        BEGIN { OFS="\t" }
        /^##contig=<ID=chr/ { sub(/ID=chr/, "ID="); print; next }
        /^#/ { print; next }
        { sub(/^chr/, "", $1); print }
    ' | bgzip -c > "$outfile"
    tabix -p vcf "$outfile"
}

# --- Helper: inject missing ##FILTER headers ---
# nirvana2vcf preserves FILTER values from the source VCF but doesn't emit ##FILTER
# header lines.  bcftools refuses to parse such files (even with --force it segfaults
# during output formatting).  This function scans data lines for FILTER values and
# injects proper ##FILTER declarations before the #CHROM line.
inject_filter_headers() {
    local infile="$1" outfile="$2"
    if [[ "$infile" == *.gz || "$infile" == *.bgz ]]; then
        gzip -dc "$infile"
    else
        cat "$infile"
    fi | awk -F'\t' '
        /^##/ { headers[++nh] = $0; next }
        /^#CHROM/ { chrom_line = $0; next }
        {
            n = split($7, arr, ";")
            for (i = 1; i <= n; i++) {
                if (arr[i] != "." && arr[i] != "PASS" && !(arr[i] in seen)) {
                    seen[arr[i]] = 1
                }
            }
            data[++nd] = $0
        }
        END {
            for (i = 1; i <= nh; i++) print headers[i]
            for (f in seen) {
                print "##FILTER=<ID=" f ",Description=\"Imported from source VCF\">"
            }
            print chrom_line
            for (i = 1; i <= nd; i++) print data[i]
        }
    ' | bgzip -c > "$outfile"
    tabix -p vcf "$outfile"
}

# --- Prepare bare-chrom version of original VCF ---
BARE_ORIGINAL="$PHASE5_DATA/HiSeq.10000.bare.vcf.gz"
if [ ! -f "$BARE_ORIGINAL" ]; then
    log "Creating bare-chrom original VCF..."
    strip_chr_prefix "$ORIGINAL_VCF" "$BARE_ORIGINAL"
    ok "Bare-chrom original VCF: $BARE_ORIGINAL"
else
    ok "Bare-chrom original VCF already exists"
fi

# ---------------------------------------------------------------------------
# Section 2: Sub-Phase 5a — VEP
# ---------------------------------------------------------------------------

log ""
log "========================================="
log "  Sub-Phase 5a: VEP Comparison"
log "========================================="

if [ "$SKIP_VEP" = "true" ]; then
    record_warn "VEP sub-phase skipped (SKIP_VEP=true)"
else
    # Check vep command
    if ! command -v vep &>/dev/null; then
        record_warn "VEP not found — skipping sub-phase 5a"
        record_warn "Install with: conda install -c bioconda ensembl-vep"
        VEP_RESULT="SKIPPED (not installed)"
    else
        # Check VEP cache (look for homo_sapiens or homo_sapiens_merged)
        VEP_CACHE_DIR="${HOME}/.vep/homo_sapiens"
        if [ -d "${HOME}/.vep/homo_sapiens_merged" ] && [ ! -d "$VEP_CACHE_DIR" ]; then
            VEP_CACHE_DIR="${HOME}/.vep/homo_sapiens_merged"
        fi
        if [ ! -d "$VEP_CACHE_DIR" ]; then
            record_warn "VEP cache not found at ~/.vep/homo_sapiens{,_merged} — skipping"
            VEP_RESULT="SKIPPED (no cache)"
        else
            VEP_OUTPUT="$PHASE5_DATA/vep_output.vcf"

            # --- Reference FASTA for HGVS output ---
            # Offline VEP needs a bgzipped+indexed FASTA to compute HGVSc/HGVSp.
            # Without it, those fields come out empty. Tolerate missing tools —
            # if samtools or the download fails, skip HGVS rather than the phase.
            REF_FASTA="$PHASE5_DATA/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
            VEP_HGVS_ARGS=()
            FASTA_FRESHLY_DOWNLOADED=false
            if [ ! -f "$REF_FASTA" ]; then
                if ! command -v samtools &>/dev/null; then
                    record_warn "samtools not found — skipping HGVS (HGVSc/HGVSp will be empty)"
                else
                    log "Downloading Ensembl GRCh37 primary-assembly FASTA (~900 MB)..."
                    REF_FASTA_TMP="$PHASE5_DATA/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz.tmp"
                    if curl -C - -L -f -o "$REF_FASTA_TMP" \
                        "https://ftp.ensembl.org/pub/grch37/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz" \
                        2> "$PHASE5_RESULTS/fasta_download.log"; then
                        # Ensembl ships as regular gzip; VEP requires bgzip.
                        log "Recompressing FASTA with bgzip..."
                        if gzip -dc "$REF_FASTA_TMP" | bgzip -c > "$REF_FASTA" \
                            && samtools faidx "$REF_FASTA" 2>> "$PHASE5_RESULTS/fasta_download.log"; then
                            rm -f "$REF_FASTA_TMP"
                            FASTA_FRESHLY_DOWNLOADED=true
                            ok "Reference FASTA ready: $REF_FASTA"
                        else
                            record_warn "FASTA bgzip/index failed — skipping HGVS"
                            rm -f "$REF_FASTA" "$REF_FASTA_TMP" "${REF_FASTA}.fai" "${REF_FASTA}.gzi"
                        fi
                    else
                        record_warn "FASTA download failed — skipping HGVS"
                        rm -f "$REF_FASTA_TMP"
                    fi
                fi
            else
                ok "Reference FASTA already exists"
            fi

            if [ -f "$REF_FASTA" ] && [ -f "${REF_FASTA}.fai" ]; then
                VEP_HGVS_ARGS=(--hgvs --fasta "$REF_FASTA")
                # Cached VEP output predates the FASTA — force re-run.
                if [ "$FASTA_FRESHLY_DOWNLOADED" = "true" ] && [ -f "$VEP_OUTPUT" ]; then
                    log "Removing stale VEP output (new FASTA available)..."
                    rm -f "$VEP_OUTPUT"
                fi
            fi

            # --- HGNC alias map for SYMBOL normalisation ---
            # Nirvana bundles Ensembl-91 gene names; VEP 115 uses current HGNC
            # symbols.  Pulling prev_symbol/alias_symbol from HGNC lets the
            # comparator resolve the ~5K version-drift name mismatches.
            if [ ! -f "$HGNC_FILE" ]; then
                log "Downloading HGNC complete set (~15 MB)..."
                HGNC_TMP="$HGNC_FILE.tmp"
                if curl -C - -L -f -o "$HGNC_TMP" \
                    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt" \
                    2> "$PHASE5_RESULTS/hgnc_download.log"; then
                    mv "$HGNC_TMP" "$HGNC_FILE"
                    ok "HGNC complete set: $HGNC_FILE"
                else
                    record_warn "HGNC download failed — SYMBOL will use exact match"
                    rm -f "$HGNC_TMP"
                fi
            else
                ok "HGNC complete set already exists"
            fi

            # --mane adds MANE_SELECT / MANE_PLUS_CLINICAL CSQ columns; the
            # comparator uses MANE_SELECT as the cross-tool canonical axis.
            # Stale VEP output without MANE columns must be re-run.
            if [ -f "$VEP_OUTPUT" ] && ! grep -q "MANE_SELECT" "$VEP_OUTPUT"; then
                log "Removing stale VEP output (missing MANE columns)..."
                rm -f "$VEP_OUTPUT"
            fi

            if [ ! -f "$VEP_OUTPUT" ]; then
                log "Running VEP on original VCF..."
                vep \
                    --input_file "$ORIGINAL_VCF" \
                    --output_file "$VEP_OUTPUT" \
                    --vcf \
                    --offline \
                    --cache \
                    --assembly GRCh37 \
                    --sift b \
                    --polyphen b \
                    --symbol \
                    --canonical \
                    --mane \
                    --biotype \
                    --numbers \
                    --merged \
                    --use_given_ref \
                    --force_overwrite \
                    ${VEP_HGVS_ARGS[@]+"${VEP_HGVS_ARGS[@]}"} \
                    2> "$PHASE5_RESULTS/vep_run.log"
                ok "VEP annotation complete: $VEP_OUTPUT"
            else
                ok "VEP output already exists: $VEP_OUTPUT"
            fi

            # Run comparison
            VEP_REPORT="$PHASE5_RESULTS/vep_comparison.txt"
            log "Running VEP comparison..."
            HGNC_ARGS=()
            if [ -f "$HGNC_FILE" ]; then
                HGNC_ARGS=(--hgnc "$HGNC_FILE")
            fi
            python3 "$COMPARE_VEP" \
                --nirvana2vcf "$J2V_VCF" \
                --vep "$VEP_OUTPUT" \
                --output "$VEP_REPORT" \
                ${HGNC_ARGS[@]+"${HGNC_ARGS[@]}"} \
                --verbose \
                2> "$PHASE5_RESULTS/vep_comparison.log"

            if grep -q "CONCORDANCE:" "$VEP_REPORT"; then
                record_pass "VEP comparison complete"
                VEP_RESULT="DONE"
            else
                record_fail "VEP comparison produced no concordance data"
                VEP_RESULT="FAIL"
            fi
        fi
    fi
fi

# ---------------------------------------------------------------------------
# Section 3: Sub-Phase 5b — SnpEff
# ---------------------------------------------------------------------------

log ""
log "========================================="
log "  Sub-Phase 5b: SnpEff Comparison"
log "========================================="

if [ "$SKIP_SNPEFF" = "true" ]; then
    record_warn "SnpEff sub-phase skipped (SKIP_SNPEFF=true)"
else
    # Check snpEff command
    if ! command -v snpEff &>/dev/null; then
        record_warn "snpEff not found — skipping sub-phase 5b"
        record_warn "Install with: conda install -c bioconda snpeff snpsift"
        SNPEFF_RESULT="SKIPPED (not installed)"
    else
        SNPEFF_ANNOTATED="$PHASE5_DATA/snpeff_annotated.vcf"
        SNPEFF_CLINVAR_VCF="$PHASE5_DATA/snpeff_clinvar.vcf"
        CLINVAR_GRCh37="$PHASE5_DATA/clinvar_GRCh37.vcf.gz"

        # Download SnpEff GRCh37.75 database if needed
        if ! snpEff databases 2>/dev/null | grep -q "GRCh37.75"; then
            log "Downloading SnpEff GRCh37.75 database..."
            snpEff download -v GRCh37.75 2> "$PHASE5_RESULTS/snpeff_download.log" || true
        fi

        # Run SnpEff on bare-chrom VCF
        if [ ! -f "$SNPEFF_ANNOTATED" ]; then
            log "Running SnpEff on bare-chrom VCF..."
            snpEff -Xmx4g -v GRCh37.75 "$BARE_ORIGINAL" \
                > "$SNPEFF_ANNOTATED" \
                2> "$PHASE5_RESULTS/snpeff_run.log"
            ok "SnpEff annotation complete: $SNPEFF_ANNOTATED"
        else
            ok "SnpEff output already exists: $SNPEFF_ANNOTATED"
        fi

        # Download ClinVar GRCh37 if needed
        if [ ! -f "$CLINVAR_GRCh37" ]; then
            log "Downloading ClinVar GRCh37 VCF..."
            curl -C - -L -o "$CLINVAR_GRCh37" \
                "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz" \
                2> "$PHASE5_RESULTS/clinvar_download.log"
            # Download index if available
            curl -C - -L -o "${CLINVAR_GRCh37}.tbi" \
                "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi" \
                2>> "$PHASE5_RESULTS/clinvar_download.log" || true
            ok "ClinVar GRCh37 downloaded"
        else
            ok "ClinVar GRCh37 already exists"
        fi

        # Run SnpSift ClinVar overlay
        if [ ! -f "$SNPEFF_CLINVAR_VCF" ]; then
            if command -v SnpSift &>/dev/null; then
                log "Running SnpSift ClinVar annotation..."
                SnpSift annotate \
                    -name CLINVAR_ \
                    -info CLNSIG,CLNDN,CLNREVSTAT \
                    "$CLINVAR_GRCh37" \
                    "$SNPEFF_ANNOTATED" \
                    > "$SNPEFF_CLINVAR_VCF" \
                    2> "$PHASE5_RESULTS/snpsift_run.log"
                ok "SnpSift ClinVar annotation complete"
            else
                warn "SnpSift not found — using SnpEff output without ClinVar overlay"
                cp "$SNPEFF_ANNOTATED" "$SNPEFF_CLINVAR_VCF"
            fi
        else
            ok "SnpEff+ClinVar output already exists"
        fi

        # Run comparison — reuse the HGNC file downloaded earlier in 5a
        # (same alias-resolution story: Nirvana's Ensembl-91 symbols vs
        # SnpEff's Ensembl-75 symbols).  Falls back to exact match if 5a
        # was skipped or the download failed.
        SNPEFF_REPORT="$PHASE5_RESULTS/snpeff_comparison.txt"
        log "Running SnpEff comparison..."
        HGNC_ARGS=()
        if [ -f "$HGNC_FILE" ]; then
            HGNC_ARGS=(--hgnc "$HGNC_FILE")
        fi
        python3 "$COMPARE_SNPEFF" \
            --nirvana2vcf "$J2V_VCF" \
            --snpeff "$SNPEFF_CLINVAR_VCF" \
            --output "$SNPEFF_REPORT" \
            ${HGNC_ARGS[@]+"${HGNC_ARGS[@]}"} \
            --verbose \
            2> "$PHASE5_RESULTS/snpeff_comparison.log"

        if grep -q "CONCORDANCE:" "$SNPEFF_REPORT"; then
            record_pass "SnpEff comparison complete"
            SNPEFF_RESULT="DONE"
        else
            record_fail "SnpEff comparison produced no concordance data"
            SNPEFF_RESULT="FAIL"
        fi
    fi
fi

# ---------------------------------------------------------------------------
# Section 4: Sub-Phase 5c — bcftools direct (gnomAD / dbSNP)
# ---------------------------------------------------------------------------

log ""
log "========================================="
log "  Sub-Phase 5c: bcftools Direct Annotation"
log "========================================="

if [ "$SKIP_BCFTOOLS" = "true" ]; then
    record_warn "bcftools sub-phase skipped (SKIP_BCFTOOLS=true)"
else
    GNOMAD_VCF="$PHASE5_DATA/gnomad.genomes.r2.1.1.sites.1.vcf.bgz"
    J2V_BARE="$PHASE5_DATA/j2v_bare.vcf.gz"
    J2V_ANNOTATED="$PHASE5_DATA/j2v_gnomad_annotated.vcf"

    # Download gnomAD v2.1.1 chr1 GRCh37 (~5 GB) if needed
    if [ ! -f "$GNOMAD_VCF" ]; then
        log "Downloading gnomAD v2.1.1 chr1 GRCh37 genomes VCF (~5 GB)..."
        log "This is a large download and may take a while."
        curl -C - -L -o "$GNOMAD_VCF" \
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.1.vcf.bgz" \
            2> "$PHASE5_RESULTS/gnomad_download.log"
        ok "gnomAD VCF downloaded"
    else
        ok "gnomAD VCF already exists"
    fi

    # Download gnomAD index if needed
    if [ ! -f "${GNOMAD_VCF}.tbi" ]; then
        log "Downloading gnomAD index..."
        curl -C - -L -o "${GNOMAD_VCF}.tbi" \
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.1.vcf.bgz.tbi" \
            2>> "$PHASE5_RESULTS/gnomad_download.log"
        ok "gnomAD index downloaded"
    else
        ok "gnomAD index already exists"
    fi

    # Prepare bare-chrom version of nirvana2vcf output with FILTER headers injected.
    # Two steps: strip_chr_prefix (chr1→1) then inject_filter_headers (add ##FILTER).
    # bcftools annotate requires both: matching contig names AND declared FILTER IDs.
    if [ ! -f "$J2V_BARE" ]; then
        log "Creating bare-chrom nirvana2vcf VCF..."
        J2V_BARE_TMP="${J2V_BARE%.vcf.gz}.tmp.vcf.gz"
        strip_chr_prefix "$J2V_VCF" "$J2V_BARE_TMP"
        inject_filter_headers "$J2V_BARE_TMP" "$J2V_BARE"
        rm -f "$J2V_BARE_TMP" "${J2V_BARE_TMP}.tbi"
        ok "Bare-chrom nirvana2vcf VCF: $J2V_BARE"
    else
        ok "Bare-chrom nirvana2vcf VCF already exists"
    fi

    # Annotate with bcftools — two-pass approach for performance.
    # The gnomAD chr1 genomes VCF is ~38 GB.  A direct `bcftools annotate -a` does a
    # streaming merge through the entire file even for only ~10K query positions.
    # Instead: (1) extract matching gnomAD sites via tabix random access, then
    # (2) annotate against the small subset.
    #
    # --force: the gnomAD v2.1.1 genomes VCF has malformed FORMAT fields (literal '-')
    # at some sites.  Safe because we only extract INFO/AF.
    GNOMAD_SUBSET="$PHASE5_DATA/gnomad_subset.vcf.gz"
    J2V_REGIONS="$PHASE5_DATA/j2v_regions.bed"

    # Step 1: Create BED regions from j2v variant positions
    if [ ! -f "$J2V_REGIONS" ]; then
        log "Creating regions file from nirvana2vcf variant positions..."
        bcftools query -f '%CHROM\t%POS\n' "$J2V_BARE" \
            | awk '{print $1"\t"$2-1"\t"$2}' > "$J2V_REGIONS"
        ok "Regions file: $(wc -l < "$J2V_REGIONS") positions"
    else
        ok "Regions file already exists"
    fi

    # Step 2: Extract matching gnomAD sites (fast tabix random access)
    if [ ! -f "$GNOMAD_SUBSET" ]; then
        log "Extracting matching gnomAD sites (tabix indexed)..."
        bcftools view --force -R "$J2V_REGIONS" "$GNOMAD_VCF" -Oz \
            -o "$GNOMAD_SUBSET" \
            2> "$PHASE5_RESULTS/bcftools_extract.log"
        tabix -p vcf "$GNOMAD_SUBSET"
        SUBSET_COUNT=$(bcftools view -H "$GNOMAD_SUBSET" 2>/dev/null | wc -l | tr -d ' ')
        ok "gnomAD subset: $SUBSET_COUNT sites"
    else
        ok "gnomAD subset already exists"
    fi

    # Step 2b: Create sites-only gnomAD subset (strip FORMAT/sample columns).
    # The gnomAD v2.1.1 genomes VCF has malformed FORMAT fields (literal '-' at
    # some sites, and CSQ annotation strings misplaced into FORMAT at others).
    # bcftools cannot parse these even with --force.  Since we only need INFO/AF,
    # strip FORMAT/genotype columns via awk to bypass bcftools' FORMAT parser.
    GNOMAD_SITES="${GNOMAD_SUBSET%.vcf.gz}.sites.vcf.gz"
    if [ ! -f "$GNOMAD_SITES" ]; then
        log "Creating sites-only gnomAD subset..."
        gzip -dc "$GNOMAD_SUBSET" | awk '
            /^##FORMAT/ { next }
            /^#CHROM/ {
                for (i = 1; i <= 8; i++) printf "%s%s", $i, (i < 8 ? "\t" : "\n")
                next
            }
            /^#/ { print; next }
            {
                for (i = 1; i <= 8; i++) printf "%s%s", $i, (i < 8 ? "\t" : "\n")
            }
        ' | bgzip -c > "$GNOMAD_SITES"
        tabix -p vcf "$GNOMAD_SITES"
        ok "Sites-only gnomAD subset: $GNOMAD_SITES"
    else
        ok "Sites-only gnomAD subset already exists"
    fi

    # Step 3: Annotate j2v VCF with the clean sites-only gnomAD subset
    if [ ! -f "$J2V_ANNOTATED" ]; then
        log "Annotating nirvana2vcf VCF with gnomAD AF..."
        bcftools annotate \
            -a "$GNOMAD_SITES" \
            -c "INFO/GNOMAD_DIRECT_AF:=INFO/AF" \
            "$J2V_BARE" \
            -o "$J2V_ANNOTATED" \
            2> "$PHASE5_RESULTS/bcftools_annotate.log"
        ok "bcftools annotation complete: $J2V_ANNOTATED"
    else
        ok "bcftools-annotated VCF already exists"
    fi

    # Run comparison
    BCFTOOLS_REPORT="$PHASE5_RESULTS/gnomad_dbsnp_comparison.txt"
    log "Running gnomAD/dbSNP comparison..."
    python3 "$COMPARE_GNOMAD" \
        --nirvana2vcf "$J2V_VCF" \
        --annotated "$J2V_ANNOTATED" \
        --output "$BCFTOOLS_REPORT" \
        --verbose \
        2> "$PHASE5_RESULTS/gnomad_dbsnp_comparison.log"

    if grep -q "CONCORDANCE:" "$BCFTOOLS_REPORT"; then
        record_pass "gnomAD/dbSNP comparison complete"
        BCFTOOLS_RESULT="DONE"
    else
        record_fail "gnomAD/dbSNP comparison produced no concordance data"
        BCFTOOLS_RESULT="FAIL"
    fi
fi

# ---------------------------------------------------------------------------
# Section 5: Concordance Matrix
# ---------------------------------------------------------------------------

log ""
log "========================================="
log "  Concordance Matrix"
log "========================================="

MATRIX="$PHASE5_RESULTS/concordance_matrix.txt"
{
    echo "============================================="
    echo "  CROSS-TOOL CONCORDANCE MATRIX"
    echo "============================================="
    echo ""
    echo "Input: HiSeq.10000.vcf.gz (9,965 SNPs, chr1, GRCh37)"
    echo ""

    # Collect CONCORDANCE lines from all sub-phase reports
    for report in "$PHASE5_RESULTS"/vep_comparison.txt \
                  "$PHASE5_RESULTS"/snpeff_comparison.txt \
                  "$PHASE5_RESULTS"/gnomad_dbsnp_comparison.txt; do
        if [ -f "$report" ]; then
            basename_report=$(basename "$report" .txt)
            echo "--- $basename_report ---"
            grep "^CONCORDANCE:" "$report" || echo "  (no data)"
            echo ""
        fi
    done

    echo "============================================="
} > "$MATRIX"
ok "Concordance matrix: $MATRIX"

# ---------------------------------------------------------------------------
# Section 6: Summary Report
# ---------------------------------------------------------------------------

log "Generating summary report..."

SUMMARY="$PHASE5_RESULTS/summary.txt"
{
    echo "============================================="
    echo "  Phase 5: Cross-Tool Comparison Summary"
    echo "============================================="
    echo ""
    echo "nirvana2vcf VCF: $J2V_VCF"
    echo "Original VCF: $ORIGINAL_VCF"
    echo ""
    printf "%-20s %-15s\n" "Sub-Phase" "Result"
    printf "%-20s %-15s\n" "---------" "------"
    printf "%-20s %-15s\n" "5a: VEP" "$VEP_RESULT"
    printf "%-20s %-15s\n" "5b: SnpEff" "$SNPEFF_RESULT"
    printf "%-20s %-15s\n" "5c: bcftools" "$BCFTOOLS_RESULT"
    echo ""

    # Print concordance data from available reports
    for report in "$PHASE5_RESULTS"/vep_comparison.txt \
                  "$PHASE5_RESULTS"/snpeff_comparison.txt \
                  "$PHASE5_RESULTS"/gnomad_dbsnp_comparison.txt; do
        if [ -f "$report" ]; then
            echo "--- $(basename "$report" .txt) ---"
            grep "^CONCORDANCE:" "$report" | while read -r line; do
                echo "  $line"
            done
            echo ""
        fi
    done

    # Overall verdict
    COMPLETED=0
    SKIPPED_COUNT=0
    FAILED=0
    for result in "$VEP_RESULT" "$SNPEFF_RESULT" "$BCFTOOLS_RESULT"; do
        case "$result" in
            DONE) COMPLETED=$((COMPLETED + 1)) ;;
            FAIL) FAILED=$((FAILED + 1)) ;;
            *)    SKIPPED_COUNT=$((SKIPPED_COUNT + 1)) ;;
        esac
    done

    if [ "$FAILED" -gt 0 ]; then
        echo "Overall: FAIL ($COMPLETED completed, $SKIPPED_COUNT skipped, $FAILED failed)"
    elif [ "$COMPLETED" -eq 3 ]; then
        echo "Overall: PASS ($COMPLETED/3 sub-phases completed)"
    elif [ "$COMPLETED" -gt 0 ]; then
        echo "Overall: PASS (partial) ($COMPLETED/3 sub-phases completed, $SKIPPED_COUNT skipped)"
    else
        echo "Overall: SKIPPED (all sub-phases skipped — install VEP/SnpEff or set SKIP_*=false)"
    fi
    echo "============================================="
} > "$SUMMARY"

ok "Summary report: $SUMMARY"

# ---------------------------------------------------------------------------
# Section 7: Update strategy.md
# ---------------------------------------------------------------------------

if grep -q '| 5\. Cross-Tool Comparison | Not started |' "$STRATEGY_MD" 2>/dev/null; then
    log "Updating Phase 5 status in strategy.md..."
    sed -i '' 's/| 5\. Cross-Tool Comparison | Not started |/| 5. Cross-Tool Comparison | Done |/' "$STRATEGY_MD"
    ok "strategy.md updated"
else
    ok "strategy.md already up to date"
fi

# ---------------------------------------------------------------------------
# Section 8: Console Summary
# ---------------------------------------------------------------------------

echo ""
echo "============================================="
echo "  Phase 5 Complete"
echo "============================================="
echo ""
echo "  Sub-phase results:"
printf "    %-20s %-15s\n" "Sub-Phase" "Result"
printf "    %-20s %-15s\n" "---------" "------"
printf "    %-20s %-15s\n" "5a: VEP" "$VEP_RESULT"
printf "    %-20s %-15s\n" "5b: SnpEff" "$SNPEFF_RESULT"
printf "    %-20s %-15s\n" "5c: bcftools" "$BCFTOOLS_RESULT"
echo ""
echo "  Reports:"
for f in "$PHASE5_RESULTS"/*.txt; do
    [ -f "$f" ] && echo "    $f"
done
echo ""
if [ "$ERRORS" -eq 0 ]; then
    echo -e "  ${GREEN}Result: Phase 5 complete${NC} ($WARNINGS warnings, $ERRORS errors)"
else
    echo -e "  ${RED}Result: Phase 5 had errors${NC} ($WARNINGS warnings, $ERRORS errors)"
fi
echo "============================================="
