#!/usr/bin/env bash
# =============================================================================
# Phase 4: Round-Trip Annotation Concordance
#
# Verifies that annotation values are faithfully converted from Nirvana JSON
# to VCF. This is a self-consistency check — both sides derive from the same
# Nirvana JSON, so the expectation is 100% concordance.
#
# Prerequisites:
#   - Phase 2 complete (mode VCFs exist in data/output/phase2/)
#   - Nirvana JSON available (HiSeq.10000.json.gz in data/output/)
#   - python3 + orjson available
#
# Usage:
#   bash validation/scripts/run_phase4.sh
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
PHASE4_RESULTS="$RESULTS_DIR/phase4"
NIRVANA_JSON="$OUTPUT_DIR/HiSeq.10000.json.gz"
COMPARE_SCRIPT="$SCRIPT_DIR/compare_annotations.py"

# Modes to test.
# Skip decomposed: identical to default for all-SNP data.
MODES=("default" "raw" "no_samples" "csq_only")

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log()   { echo -e "${BLUE}[Phase4]${NC} $*"; }
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

# Nirvana JSON
if [ ! -f "$NIRVANA_JSON" ]; then
    err "Nirvana JSON not found: $NIRVANA_JSON"
    err "Run Phase 1 first: bash validation/scripts/run_phase1.sh"
    exit 1
fi
ok "Nirvana JSON found: $NIRVANA_JSON"

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

# python3
if ! command -v python3 &>/dev/null; then
    err "python3 not found"
    exit 1
fi
ok "python3 available"

# orjson
if ! python3 -c "import orjson" 2>/dev/null; then
    err "orjson not available. Install with: pip install orjson"
    exit 1
fi
ok "orjson available"

# compare_annotations.py
if [ ! -f "$COMPARE_SCRIPT" ]; then
    err "compare_annotations.py not found: $COMPARE_SCRIPT"
    exit 1
fi
ok "compare_annotations.py found"

# ---------------------------------------------------------------------------
# Section 2: Per-mode Annotation Comparison
# ---------------------------------------------------------------------------

log "Running annotation concordance checks..."

mkdir -p "$PHASE4_RESULTS"

declare -a MODE_RESULTS=()

for mode in "${MODES[@]}"; do
    vcf="$PHASE2_DIR/HiSeq.10000.${mode}.vcf"
    report="$PHASE4_RESULTS/${mode}_annotation_concordance.txt"

    log "  Mode: $mode"

    python3 "$COMPARE_SCRIPT" \
        --json "$NIRVANA_JSON" \
        --vcf "$vcf" \
        --mode "$mode" \
        --output "$report" \
        --verbose \
        2> "$PHASE4_RESULTS/${mode}_annotation_concordance.log"

    # Check for PASS/FAIL in output
    if grep -q "RESULT: PASS" "$report"; then
        record_pass "Annotation concordance: $mode"
        MODE_RESULTS+=("PASS")
    else
        record_fail "Annotation concordance: $mode"
        MODE_RESULTS+=("FAIL")
    fi
done

ok "Annotation concordance reports written to $PHASE4_RESULTS/"

# ---------------------------------------------------------------------------
# Section 3: Summary Report
# ---------------------------------------------------------------------------

log "Generating summary report..."

SUMMARY="$PHASE4_RESULTS/summary.txt"
{
    echo "============================================="
    echo "  Phase 4: Annotation Concordance Summary"
    echo "============================================="
    echo ""
    echo "Nirvana JSON: $NIRVANA_JSON"
    echo ""
    echo "Note: decomposed mode skipped (identical to default for all-SNP data)"
    echo ""
    printf "%-15s %-15s\n" "Mode" "Result"
    printf "%-15s %-15s\n" "----" "------"
    for i in "${!MODES[@]}"; do
        printf "%-15s %-15s\n" "${MODES[$i]}" "${MODE_RESULTS[$i]}"
    done
    echo ""

    # Extract key stats from default mode report
    if [ -f "$PHASE4_RESULTS/default_annotation_concordance.txt" ]; then
        echo "Coverage (from default mode):"
        sed -n '/^COVERAGE SUMMARY/,/^$/p' \
            "$PHASE4_RESULTS/default_annotation_concordance.txt" \
            | grep -v "^-" | grep -v "^COVERAGE"
        echo ""
    fi

    # Overall verdict
    ALL_PASS=true
    for i in "${!MODES[@]}"; do
        if [ "${MODE_RESULTS[$i]}" != "PASS" ]; then
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
# Section 4: Update strategy.md
# ---------------------------------------------------------------------------

if grep -q '| 4\. Annotation Concordance | Not started |' "$STRATEGY_MD" 2>/dev/null; then
    log "Updating Phase 4 status in strategy.md..."
    sed -i '' 's/| 4\. Annotation Concordance | Not started |/| 4. Annotation Concordance | Done |/' "$STRATEGY_MD"
    ok "strategy.md updated"
else
    ok "strategy.md already up to date"
fi

# ---------------------------------------------------------------------------
# Section 5: Console Summary
# ---------------------------------------------------------------------------

echo ""
echo "============================================="
echo "  Phase 4 Complete"
echo "============================================="
echo ""
echo "  Nirvana JSON: $(basename "$NIRVANA_JSON")"
echo ""
echo "  Annotation concordance results:"
printf "    %-15s %-15s\n" "Mode" "Result"
printf "    %-15s %-15s\n" "----" "------"
for i in "${!MODES[@]}"; do
    printf "    %-15s %-15s\n" "${MODES[$i]}" "${MODE_RESULTS[$i]}"
done
echo ""
echo "  Note: decomposed mode skipped (identical to default for all-SNP data)"
echo ""
echo "  Reports:"
for mode in "${MODES[@]}"; do
    echo "    $PHASE4_RESULTS/${mode}_annotation_concordance.txt"
done
echo "    $SUMMARY"
echo ""
if [ "$ERRORS" -eq 0 ]; then
    echo -e "  ${GREEN}Result: PASS${NC} ($WARNINGS warnings, $ERRORS errors)"
else
    echo -e "  ${RED}Result: FAIL${NC} ($WARNINGS warnings, $ERRORS errors)"
fi
echo "============================================="
