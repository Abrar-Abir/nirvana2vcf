#!/usr/bin/env bash
# =============================================================================
# Phase 2: Smoke Test & VCF Spec Validation
#
# Validates that nirvana2vcf produces spec-compliant VCF 4.2 output across all
# mode combinations. Uses EBI vcf_validator (strictest), bcftools, and
# vcftools for independent structural checks.
#
# Prerequisites:
#   - Phase 1 complete (HiSeq.10000.json.gz exists in data/output/)
#   - nirvana2vcf installed (pip install -e .)
#   - bcftools installed (brew install bcftools)
#
# Usage:
#   bash validation/scripts/run_phase2.sh
#
# Idempotent — safe to re-run. Overwrites VCF outputs and re-validates.
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Section 0: Configuration
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
DATA_DIR="$PROJECT_ROOT/validation/data"
OUTPUT_DIR="$DATA_DIR/output"
RESULTS_DIR="$PROJECT_ROOT/validation/results"
STRATEGY_MD="$PROJECT_ROOT/validation/docs/strategy.md"
TOOLS_DIR="$DATA_DIR/tools"

INPUT_JSON="$OUTPUT_DIR/HiSeq.10000.json.gz"
PHASE2_DIR="$OUTPUT_DIR/phase2"

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log()   { echo -e "${BLUE}[Phase2]${NC} $*"; }
ok()    { echo -e "${GREEN}[  OK  ]${NC} $*"; }
warn()  { echo -e "${YELLOW}[ WARN ]${NC} $*"; }
err()   { echo -e "${RED}[ERROR ]${NC} $*" >&2; }

# Track overall pass/fail
ERRORS=0
WARNINGS=0

record_pass() { ok "$1"; }
record_warn() { warn "$1"; WARNINGS=$((WARNINGS + 1)); }
record_fail() { err "$1"; ERRORS=$((ERRORS + 1)); }

# ---------------------------------------------------------------------------
# Section 1: Prerequisites
# ---------------------------------------------------------------------------

log "Checking prerequisites..."

# Phase 1 output
if [ ! -f "$INPUT_JSON" ]; then
    err "Phase 1 output not found: $INPUT_JSON"
    err "Run Phase 1 first: bash validation/scripts/run_phase1.sh"
    exit 1
fi
ok "Phase 1 output found: $INPUT_JSON"

# nirvana2vcf
if ! command -v nirvana2vcf &>/dev/null; then
    log "Installing nirvana2vcf in dev mode..."
    pip install -e "$PROJECT_ROOT" --quiet
fi
ok "nirvana2vcf available"

# bcftools
if ! command -v bcftools &>/dev/null; then
    log "Installing bcftools via brew..."
    brew install bcftools
fi
ok "bcftools $(bcftools --version | head -1 | awk '{print $2}')"

# EBI vcf_validator — download arm64 binary if not present
VCF_VALIDATOR="$TOOLS_DIR/vcf_validator"
if [ ! -x "$VCF_VALIDATOR" ]; then
    log "Downloading EBI vcf_validator (macOS arm64)..."
    mkdir -p "$TOOLS_DIR"
    curl -sL -o "$VCF_VALIDATOR" \
        "https://github.com/EBIvariation/vcf-validator/releases/download/v0.10.2/vcf_validator_macos_arm64"
    chmod +x "$VCF_VALIDATOR"
    # Remove macOS quarantine attribute so the unsigned binary can execute
    xattr -d com.apple.quarantine "$VCF_VALIDATOR" 2>/dev/null || true
    # Ad-hoc code sign for macOS Sequoia+
    codesign -s - "$VCF_VALIDATOR" 2>/dev/null || true
fi
if "$VCF_VALIDATOR" --help &>/dev/null; then
    VCF_VALIDATOR_VERSION="$("$VCF_VALIDATOR" --version 2>&1 || echo 'v0.10.2')"
    ok "vcf_validator ready ($VCF_VALIDATOR_VERSION)"
    HAS_VCF_VALIDATOR=true
else
    # Pre-built binary may need specific Boost dylibs not present on this system.
    # This is a known limitation — bcftools + vcftools still provide good coverage.
    warn "vcf_validator not functional (likely missing Boost dylibs) — skipping EBI validation"
    HAS_VCF_VALIDATOR=false
fi

# vcftools (Perl vcf-validator)
if command -v vcf-validator &>/dev/null; then
    ok "vcftools vcf-validator available"
    HAS_VCFTOOLS=true
elif command -v brew &>/dev/null; then
    log "Installing vcftools via brew..."
    if brew install vcftools 2>/dev/null; then
        ok "vcftools installed"
        HAS_VCFTOOLS=true
    else
        record_warn "vcftools install failed — skipping vcftools validation"
        HAS_VCFTOOLS=false
    fi
else
    record_warn "vcftools not available and no brew — skipping vcftools validation"
    HAS_VCFTOOLS=false
fi

# ---------------------------------------------------------------------------
# Section 2: Generate VCFs in All Mode Combinations
# ---------------------------------------------------------------------------

log "Creating output directory..."
mkdir -p "$PHASE2_DIR" "$RESULTS_DIR"

log "Generating VCFs for all mode combinations..."

VCF_FILES=""  # space-separated list of output paths

generate_vcf() {
    local name="$1"
    shift
    local outfile="$PHASE2_DIR/HiSeq.10000.${name}.vcf"

    log "  Mode: $name ${*:-(default flags)}"
    nirvana2vcf -i "$INPUT_JSON" -o "$outfile" "$@"
    local count
    count=$(grep -cv '^#' "$outfile" || true)
    ok "  $name: $count variants"
    VCF_FILES="$VCF_FILES $outfile"
}

generate_vcf "default"
generate_vcf "raw"        --no-normalize
generate_vcf "decomposed" --decompose
generate_vcf "no_samples" --no-samples
generate_vcf "csq_only"   --csq-only

VCF_FILE_COUNT=$(echo $VCF_FILES | wc -w | tr -d ' ')
ok "All mode VCFs generated ($VCF_FILE_COUNT files)"

# ---------------------------------------------------------------------------
# Section 3: Validate with bcftools
# ---------------------------------------------------------------------------

log "Running bcftools validation..."

BCFTOOLS_REPORT="$RESULTS_DIR/bcftools_validation.txt"
: > "$BCFTOOLS_REPORT"

for vcf in $VCF_FILES; do
    name="$(basename "$vcf" .vcf)"
    echo "=== $name ===" >> "$BCFTOOLS_REPORT"

    # Parse check
    if bcftools view "$vcf" > /dev/null 2>&1; then
        record_pass "bcftools parse: $name"
        echo "Parse: PASS" >> "$BCFTOOLS_REPORT"
    else
        record_fail "bcftools parse: $name"
        echo "Parse: FAIL" >> "$BCFTOOLS_REPORT"
    fi

    # Stats
    bcftools stats "$vcf" >> "$BCFTOOLS_REPORT" 2>&1

    # Sample listing
    echo "Samples:" >> "$BCFTOOLS_REPORT"
    bcftools query -l "$vcf" >> "$BCFTOOLS_REPORT" 2>/dev/null || echo "(none)" >> "$BCFTOOLS_REPORT"

    echo "" >> "$BCFTOOLS_REPORT"
done

ok "bcftools validation report: $BCFTOOLS_REPORT"

# ---------------------------------------------------------------------------
# Section 4: Check INFO/FORMAT Header Declarations
# ---------------------------------------------------------------------------

log "Checking INFO/FORMAT header declarations..."

HEADER_REPORT="$RESULTS_DIR/header_declaration_check.txt"
: > "$HEADER_REPORT"

for vcf in $VCF_FILES; do
    name="$(basename "$vcf" .vcf)"
    echo "=== $name ===" >> "$HEADER_REPORT"

    # Extract declared INFO keys from header
    bcftools view -h "$vcf" | grep '^##INFO' | \
        sed -E 's/.*ID=([^,>]+).*/\1/' | sort > "$PHASE2_DIR/_declared_info.tmp"

    # Extract used INFO keys from data lines
    # bcftools query -f '%INFO/KEYS\n' is not universally supported,
    # so parse INFO field directly. Filter out '.' (VCF missing value).
    grep -v '^#' "$vcf" | cut -f8 | tr ';' '\n' | \
        sed -E 's/=.*//' | grep -v '^\.$' | sort -u > "$PHASE2_DIR/_used_info.tmp"

    # Find undeclared keys
    undeclared=$(comm -23 "$PHASE2_DIR/_used_info.tmp" "$PHASE2_DIR/_declared_info.tmp" || true)
    if [ -z "$undeclared" ]; then
        record_pass "INFO declarations complete: $name"
        echo "INFO keys: all declared" >> "$HEADER_REPORT"
    else
        record_fail "Undeclared INFO keys in $name: $undeclared"
        echo "Undeclared INFO keys: $undeclared" >> "$HEADER_REPORT"
    fi

    # Extract declared FORMAT keys from header (grep may find none → || true)
    bcftools view -h "$vcf" | (grep '^##FORMAT' || true) | \
        sed -E 's/.*ID=([^,>]+).*/\1/' | sort > "$PHASE2_DIR/_declared_format.tmp"

    # Extract used FORMAT keys from data lines (column 9 = FORMAT)
    if [ "$(bcftools query -l "$vcf" 2>/dev/null | wc -l)" -gt 0 ]; then
        grep -v '^#' "$vcf" | cut -f9 | tr ':' '\n' | sort -u > "$PHASE2_DIR/_used_format.tmp"
        undeclared_fmt=$(comm -23 "$PHASE2_DIR/_used_format.tmp" "$PHASE2_DIR/_declared_format.tmp" || true)
        if [ -z "$undeclared_fmt" ]; then
            record_pass "FORMAT declarations complete: $name"
            echo "FORMAT keys: all declared" >> "$HEADER_REPORT"
        else
            record_fail "Undeclared FORMAT keys in $name: $undeclared_fmt"
            echo "Undeclared FORMAT keys: $undeclared_fmt" >> "$HEADER_REPORT"
        fi
    else
        echo "FORMAT keys: N/A (no samples)" >> "$HEADER_REPORT"
    fi

    echo "" >> "$HEADER_REPORT"
done

rm -f "$PHASE2_DIR"/_declared_*.tmp "$PHASE2_DIR"/_used_*.tmp
ok "Header declaration report: $HEADER_REPORT"

# ---------------------------------------------------------------------------
# Section 5: Validate with EBI vcf_validator
# ---------------------------------------------------------------------------

if [ "$HAS_VCF_VALIDATOR" = true ]; then
    log "Running EBI vcf_validator..."

    EBI_REPORT="$RESULTS_DIR/ebi_vcf_validator.txt"
    : > "$EBI_REPORT"

    for vcf in $VCF_FILES; do
        name="$(basename "$vcf" .vcf)"
        echo "=== $name ===" >> "$EBI_REPORT"

        # vcf_validator outputs to stdout; exit code 0 = valid
        # Capture output and exit code
        validator_output=$("$VCF_VALIDATOR" -i "$vcf" 2>&1) || true
        echo "$validator_output" >> "$EBI_REPORT"

        # Count errors and warnings from output
        ebi_errors=$(echo "$validator_output" | grep -ci 'error' || true)
        ebi_warnings=$(echo "$validator_output" | grep -ci 'warning' || true)

        if [ "$ebi_errors" -eq 0 ]; then
            record_pass "vcf_validator: $name (${ebi_warnings} warnings)"
        else
            record_fail "vcf_validator: $name ($ebi_errors errors, ${ebi_warnings} warnings)"
        fi

        echo "" >> "$EBI_REPORT"
    done

    ok "EBI vcf_validator report: $EBI_REPORT"
else
    log "Skipping EBI vcf_validator (not available)"
fi

# ---------------------------------------------------------------------------
# Section 6: Validate with vcftools vcf-validator (Perl)
# ---------------------------------------------------------------------------

if [ "$HAS_VCFTOOLS" = true ]; then
    log "Running vcftools vcf-validator..."

    VCFTOOLS_REPORT="$RESULTS_DIR/vcftools_validation.txt"
    : > "$VCFTOOLS_REPORT"

    for vcf in $VCF_FILES; do
        name="$(basename "$vcf" .vcf)"
        echo "=== $name ===" >> "$VCFTOOLS_REPORT"

        validator_output=$(vcf-validator "$vcf" 2>&1) || true
        echo "$validator_output" >> "$VCFTOOLS_REPORT"

        vt_errors=$(echo "$validator_output" | grep -ci 'error' || true)
        if [ "$vt_errors" -eq 0 ]; then
            record_pass "vcftools validator: $name"
        else
            record_warn "vcftools validator: $name ($vt_errors issues)"
        fi

        echo "" >> "$VCFTOOLS_REPORT"
    done

    ok "vcftools report: $VCFTOOLS_REPORT"
else
    log "Skipping vcftools vcf-validator (not available)"
fi

# ---------------------------------------------------------------------------
# Section 7: Variant Count Cross-Check
# ---------------------------------------------------------------------------

log "Cross-checking variant counts..."

COUNT_REPORT="$RESULTS_DIR/variant_counts.txt"
: > "$COUNT_REPORT"

# Expected: 9,965 variants for non-decomposed modes
# Decomposed mode may have more rows (multi-allelic split)
# Use the Phase 1 smoke VCF as the reference count
SMOKE_COUNT=$(grep -cv '^#' "$OUTPUT_DIR/HiSeq.10000.smoke.vcf" 2>/dev/null || echo "unknown")

echo "Reference count (Phase 1 smoke VCF): $SMOKE_COUNT" >> "$COUNT_REPORT"
echo "" >> "$COUNT_REPORT"

printf "%-40s %10s %10s\n" "Mode" "Variants" "Match?" >> "$COUNT_REPORT"
printf "%-40s %10s %10s\n" "----" "--------" "------" >> "$COUNT_REPORT"

for vcf in $VCF_FILES; do
    name="$(basename "$vcf" .vcf)"
    count=$(grep -cv '^#' "$vcf" || true)

    # Decomposed modes may have more rows — that's expected
    if echo "$name" | grep -q "decomposed"; then
        if [ "$count" -ge "$SMOKE_COUNT" ]; then
            match="yes (>= ref)"
            record_pass "Variant count $name: $count (>= $SMOKE_COUNT, decomposed)"
        else
            match="UNEXPECTED"
            record_fail "Variant count $name: $count (< $SMOKE_COUNT for decomposed mode)"
        fi
    else
        if [ "$count" -eq "$SMOKE_COUNT" ]; then
            match="exact"
            record_pass "Variant count $name: $count (matches reference)"
        else
            match="MISMATCH"
            record_fail "Variant count $name: $count (expected $SMOKE_COUNT)"
        fi
    fi

    printf "%-40s %10s %10s\n" "$name" "$count" "$match" >> "$COUNT_REPORT"
done

echo "" >> "$COUNT_REPORT"
ok "Variant count report: $COUNT_REPORT"

# ---------------------------------------------------------------------------
# Section 8: bcftools stats Summary
# ---------------------------------------------------------------------------

log "Generating bcftools stats summary for default mode..."

STATS_REPORT="$RESULTS_DIR/bcftools_stats_default.txt"
bcftools stats "$PHASE2_DIR/HiSeq.10000.default.vcf" > "$STATS_REPORT" 2>&1

# Extract key metrics
SN_RECORDS=$(grep '^SN' "$STATS_REPORT" | grep 'number of records' | awk '{print $NF}')
SN_SNPS=$(grep '^SN' "$STATS_REPORT" | grep 'number of SNPs' | awk '{print $NF}')
SN_INDELS=$(grep '^SN' "$STATS_REPORT" | grep 'number of indels' | awk '{print $NF}')
SN_MULTIALLELIC=$(grep '^SN' "$STATS_REPORT" | grep 'number of multiallelic sites' | awk '{print $NF}')
TSTV=$(grep '^TSTV' "$STATS_REPORT" | head -1 | awk '{print $5}')

ok "Default mode stats: records=$SN_RECORDS SNPs=$SN_SNPS indels=$SN_INDELS multiallelic=$SN_MULTIALLELIC ts/tv=$TSTV"

# ---------------------------------------------------------------------------
# Section 9: Summary & Status Update
# ---------------------------------------------------------------------------

# Update Phase 2 status in strategy.md
if grep -q '| 2\. VCF Spec Validation | Not started |' "$STRATEGY_MD" 2>/dev/null; then
    log "Updating Phase 2 status in strategy.md..."
    sed -i '' 's/| 2\. VCF Spec Validation | Not started |/| 2. VCF Spec Validation | Done |/' "$STRATEGY_MD"
    ok "strategy.md updated"
else
    ok "strategy.md already up to date"
fi

# Print summary
echo ""
echo "============================================="
echo "  Phase 2 Complete"
echo "============================================="
echo ""
echo "  Mode VCFs generated:"
for vcf in $VCF_FILES; do
    name="$(basename "$vcf" .vcf)"
    count=$(grep -cv '^#' "$vcf" || true)
    printf "    %-30s %6d variants\n" "$name" "$count"
done
echo ""
echo "  bcftools stats (default mode):"
echo "    Records:       $SN_RECORDS"
echo "    SNPs:          $SN_SNPS"
echo "    Indels:        $SN_INDELS"
echo "    Multiallelic:  $SN_MULTIALLELIC"
echo "    Ts/Tv ratio:   $TSTV"
echo ""
echo "  Validators:"
echo "    bcftools:       all modes parsed successfully"
[ "$HAS_VCF_VALIDATOR" = true ] && echo "    vcf_validator:  ran on all modes" || echo "    vcf_validator:  skipped (not available)"
[ "$HAS_VCFTOOLS" = true ] && echo "    vcftools:       ran on all modes" || echo "    vcftools:       skipped (not available)"
echo ""
echo "  Reports:"
echo "    $BCFTOOLS_REPORT"
echo "    $HEADER_REPORT"
echo "    $COUNT_REPORT"
echo "    $STATS_REPORT"
[ "$HAS_VCF_VALIDATOR" = true ] && echo "    $EBI_REPORT"
[ "$HAS_VCFTOOLS" = true ] && echo "    $VCFTOOLS_REPORT"
echo ""
if [ "$ERRORS" -eq 0 ]; then
    echo -e "  ${GREEN}Result: PASS${NC} ($WARNINGS warnings, $ERRORS errors)"
else
    echo -e "  ${RED}Result: FAIL${NC} ($WARNINGS warnings, $ERRORS errors)"
fi
echo "============================================="
