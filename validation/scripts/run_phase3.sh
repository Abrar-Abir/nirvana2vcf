#!/usr/bin/env bash
# =============================================================================
# Phase 3: Round-Trip Position Concordance
#
# Verifies that the same (CHROM, POS, REF, ALT) variant tuples survive the
# VCF -> Nirvana JSON -> nirvana2vcf VCF round-trip. Uses both a quick bcftools +
# comm comparison and a full Python comparison (compare_positions.py).
#
# Prerequisites:
#   - Phase 2 complete (mode VCFs exist in data/output/phase2/)
#   - Original VCF available (HiSeq.10000.vcf.gz in data/)
#   - bcftools installed
#   - python3 available
#
# Usage:
#   bash validation/scripts/run_phase3.sh
#
# Idempotent — safe to re-run. Overwrites results.
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

PHASE2_DIR="$OUTPUT_DIR/phase2"
PHASE3_RESULTS="$RESULTS_DIR/phase3"
ORIGINAL_VCF="$DATA_DIR/HiSeq.10000.vcf.gz"
COMPARE_SCRIPT="$SCRIPT_DIR/compare_positions.py"

# Modes to test (skip csq_only — same positions as default;
# skip --no-normalize --decompose — identical to raw for all-SNP data)
MODES=("default" "raw" "decomposed" "no_samples")

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log()   { echo -e "${BLUE}[Phase3]${NC} $*"; }
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

# Original VCF
if [ ! -f "$ORIGINAL_VCF" ]; then
    err "Original VCF not found: $ORIGINAL_VCF"
    err "Run Phase 1 first: bash validation/scripts/run_phase1.sh"
    exit 1
fi
ok "Original VCF found: $ORIGINAL_VCF"

# Phase 2 mode VCFs
for mode in "${MODES[@]}"; do
    vcf="$PHASE2_DIR/HiSeq.10000.${mode}.vcf"
    if [ ! -f "$vcf" ]; then
        err "Phase 2 VCF not found: $vcf"
        err "Run Phase 2 first: bash validation/scripts/run_phase2.sh"
        exit 1
    fi
done
ok "All Phase 2 mode VCFs found (${#MODES[@]} modes)"

# bcftools
if ! command -v bcftools &>/dev/null; then
    err "bcftools not found. Install with: brew install bcftools"
    exit 1
fi
ok "bcftools $(bcftools --version | head -1 | awk '{print $2}')"

# python3
if ! command -v python3 &>/dev/null; then
    err "python3 not found"
    exit 1
fi
ok "python3 available"

# compare_positions.py
if [ ! -f "$COMPARE_SCRIPT" ]; then
    err "compare_positions.py not found: $COMPARE_SCRIPT"
    exit 1
fi
ok "compare_positions.py found"

# ---------------------------------------------------------------------------
# Section 2: Quick bcftools + comm Comparison
# ---------------------------------------------------------------------------

log "Running bcftools + comm position comparison..."

mkdir -p "$PHASE3_RESULTS"

BCFTOOLS_REPORT="$PHASE3_RESULTS/bcftools_concordance.txt"
: > "$BCFTOOLS_REPORT"

# Extract original sites once
log "  Extracting original VCF sites..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$ORIGINAL_VCF" | LC_ALL=C sort > "$PHASE3_RESULTS/original_sites.tsv"
ORIGINAL_COUNT=$(wc -l < "$PHASE3_RESULTS/original_sites.tsv" | tr -d ' ')
ok "  Original: $ORIGINAL_COUNT site tuples"

echo "Original VCF: $ORIGINAL_VCF" >> "$BCFTOOLS_REPORT"
echo "Original site tuples: $ORIGINAL_COUNT" >> "$BCFTOOLS_REPORT"
echo "" >> "$BCFTOOLS_REPORT"

# Per-mode comparison
declare -a BCFTOOLS_RESULTS=()

for mode in "${MODES[@]}"; do
    vcf="$PHASE2_DIR/HiSeq.10000.${mode}.vcf"
    log "  Mode: $mode"

    # Extract sites
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$vcf" | LC_ALL=C sort > "$PHASE3_RESULTS/${mode}_sites.tsv"
    MODE_COUNT=$(wc -l < "$PHASE3_RESULTS/${mode}_sites.tsv" | tr -d ' ')

    # Set comparisons
    SHARED=$(comm -12 "$PHASE3_RESULTS/original_sites.tsv" "$PHASE3_RESULTS/${mode}_sites.tsv" | wc -l | tr -d ' ')
    ONLY_ORIG=$(comm -23 "$PHASE3_RESULTS/original_sites.tsv" "$PHASE3_RESULTS/${mode}_sites.tsv" | wc -l | tr -d ' ')
    ONLY_J2V=$(comm -13 "$PHASE3_RESULTS/original_sites.tsv" "$PHASE3_RESULTS/${mode}_sites.tsv" | wc -l | tr -d ' ')

    echo "=== $mode ===" >> "$BCFTOOLS_REPORT"
    echo "  nirvana2vcf site tuples: $MODE_COUNT" >> "$BCFTOOLS_REPORT"
    echo "  Shared:        $SHARED" >> "$BCFTOOLS_REPORT"
    echo "  Only original: $ONLY_ORIG" >> "$BCFTOOLS_REPORT"
    echo "  Only nirvana2vcf: $ONLY_J2V" >> "$BCFTOOLS_REPORT"

    if [ "$SHARED" -eq "$ORIGINAL_COUNT" ] && [ "$ONLY_ORIG" -eq 0 ] && [ "$ONLY_J2V" -eq 0 ]; then
        record_pass "bcftools concordance: $mode ($SHARED/$ORIGINAL_COUNT shared, 0 only-original, 0 only-j2v)"
        echo "  Result: PASS" >> "$BCFTOOLS_REPORT"
        BCFTOOLS_RESULTS+=("PASS")
    else
        record_fail "bcftools concordance: $mode ($SHARED/$ORIGINAL_COUNT shared, $ONLY_ORIG only-original, $ONLY_J2V only-j2v)"
        echo "  Result: FAIL" >> "$BCFTOOLS_REPORT"
        BCFTOOLS_RESULTS+=("FAIL")
    fi
    echo "" >> "$BCFTOOLS_REPORT"
done

ok "bcftools concordance report: $BCFTOOLS_REPORT"

# ---------------------------------------------------------------------------
# Section 3: Full Python Comparison (compare_positions.py)
# ---------------------------------------------------------------------------

log "Running Python position comparison (compare_positions.py)..."

declare -a PYTHON_RESULTS=()

for mode in "${MODES[@]}"; do
    vcf="$PHASE2_DIR/HiSeq.10000.${mode}.vcf"
    log "  Mode: $mode"

    python3 "$COMPARE_SCRIPT" \
        --nirvana2vcf "$vcf" \
        --reference "$ORIGINAL_VCF" \
        > "$PHASE3_RESULTS/${mode}_position_concordance.txt" \
        2> "$PHASE3_RESULTS/${mode}_position_concordance.log"

    # Check for PASS in output
    if grep -q "RESULT: PASS" "$PHASE3_RESULTS/${mode}_position_concordance.txt"; then
        record_pass "Python concordance: $mode"
        PYTHON_RESULTS+=("PASS")
    else
        record_fail "Python concordance: $mode"
        PYTHON_RESULTS+=("FAIL")
    fi
done

ok "Python concordance reports written to $PHASE3_RESULTS/"

# ---------------------------------------------------------------------------
# Section 4: Summary Report
# ---------------------------------------------------------------------------

log "Generating summary report..."

SUMMARY="$PHASE3_RESULTS/summary.txt"
{
    echo "============================================="
    echo "  Phase 3: Position Concordance Summary"
    echo "============================================="
    echo ""
    echo "Original VCF: $ORIGINAL_VCF"
    echo "Original site tuples: $ORIGINAL_COUNT"
    echo ""
    printf "%-15s %-15s %-15s\n" "Mode" "bcftools+comm" "Python"
    printf "%-15s %-15s %-15s\n" "----" "-------------" "------"
    for i in "${!MODES[@]}"; do
        printf "%-15s %-15s %-15s\n" "${MODES[$i]}" "${BCFTOOLS_RESULTS[$i]}" "${PYTHON_RESULTS[$i]}"
    done
    echo ""

    # Overall verdict
    ALL_PASS=true
    for i in "${!MODES[@]}"; do
        if [ "${BCFTOOLS_RESULTS[$i]}" != "PASS" ] || [ "${PYTHON_RESULTS[$i]}" != "PASS" ]; then
            ALL_PASS=false
        fi
    done

    if [ "$ALL_PASS" = true ]; then
        echo "Overall: PASS (100% concordance across all modes)"
    else
        echo "Overall: FAIL (see individual mode results above)"
    fi
    echo "============================================="
} > "$SUMMARY"

ok "Summary report: $SUMMARY"

# ---------------------------------------------------------------------------
# Section 5: Update strategy.md
# ---------------------------------------------------------------------------

if grep -q '| 3\. Position Concordance | Not started |' "$STRATEGY_MD" 2>/dev/null; then
    log "Updating Phase 3 status in strategy.md..."
    sed -i '' 's/| 3\. Position Concordance | Not started |/| 3. Position Concordance | Done |/' "$STRATEGY_MD"
    ok "strategy.md updated"
else
    ok "strategy.md already up to date"
fi

# ---------------------------------------------------------------------------
# Section 6: Console Summary
# ---------------------------------------------------------------------------

echo ""
echo "============================================="
echo "  Phase 3 Complete"
echo "============================================="
echo ""
echo "  Original VCF: $(basename "$ORIGINAL_VCF")"
echo "  Site tuples:   $ORIGINAL_COUNT"
echo ""
echo "  Position concordance results:"
printf "    %-15s %-15s %-15s\n" "Mode" "bcftools+comm" "Python"
printf "    %-15s %-15s %-15s\n" "----" "-------------" "------"
for i in "${!MODES[@]}"; do
    printf "    %-15s %-15s %-15s\n" "${MODES[$i]}" "${BCFTOOLS_RESULTS[$i]}" "${PYTHON_RESULTS[$i]}"
done
echo ""
echo "  Reports:"
echo "    $BCFTOOLS_REPORT"
for mode in "${MODES[@]}"; do
    echo "    $PHASE3_RESULTS/${mode}_position_concordance.txt"
done
echo "    $SUMMARY"
echo ""
if [ "$ERRORS" -eq 0 ]; then
    echo -e "  ${GREEN}Result: PASS${NC} ($WARNINGS warnings, $ERRORS errors)"
else
    echo -e "  ${RED}Result: FAIL${NC} ($WARNINGS warnings, $ERRORS errors)"
fi
echo "============================================="
