#!/usr/bin/env bash
# =============================================================================
# Phase 1: Generate Nirvana JSON from Public Data
#
# Idempotent script — safe to re-run after failures, resumes from where it
# left off. Each step checks for its output artifact before running.
#
# Usage:
#   bash validation/scripts/run_phase1.sh
#
# Environment overrides:
#   GENOME_ASSEMBLY=GRCh38  (default: GRCh37 — matches HiSeq.10000 test VCF)
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Section 0: Configuration
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
DATA_DIR="$PROJECT_ROOT/validation/data"
STRATEGY_MD="$PROJECT_ROOT/validation/docs/strategy.md"

GENOME_ASSEMBLY="${GENOME_ASSEMBLY:-GRCh37}"

NIRVANA_DIR="$DATA_DIR/Nirvana"
NIRVANA_DATA="$DATA_DIR/nirvana_data"
OUTPUT_DIR="$DATA_DIR/output"

NIRVANA_DLL="$NIRVANA_DIR/bin/Release/net6.0/Nirvana.dll"
DOWNLOADER_DLL="$NIRVANA_DIR/bin/Release/net6.0/Downloader.dll"

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log()   { echo -e "${BLUE}[Phase1]${NC} $*"; }
ok()    { echo -e "${GREEN}[  OK  ]${NC} $*"; }
warn()  { echo -e "${YELLOW}[ WARN ]${NC} $*"; }
err()   { echo -e "${RED}[ERROR ]${NC} $*" >&2; }

# ---------------------------------------------------------------------------
# Section 1: Prerequisites
# ---------------------------------------------------------------------------

log "Checking prerequisites..."

if ! command -v brew &>/dev/null; then
    err "Homebrew (brew) is required but not found. Install from https://brew.sh"
    exit 1
fi
ok "brew found"

for tool in bcftools bgzip tabix; do
    if ! command -v "$tool" &>/dev/null; then
        case "$tool" in
            bcftools) log "Installing bcftools via brew..."; brew install bcftools ;;
            bgzip|tabix) log "Installing htslib via brew..."; brew install htslib ;;
        esac
    fi
done
ok "bcftools, bgzip, tabix available"

if ! command -v dotnet &>/dev/null; then
    err "dotnet is required but not found. Install .NET SDK from https://dotnet.microsoft.com/download"
    exit 1
fi
# Verify .NET 6.0 runtime is available
if ! dotnet --list-runtimes 2>/dev/null | grep -q 'Microsoft.NETCore.App 6\.'; then
    err ".NET 6.0 runtime is required. Install with: dotnet-install.sh -Runtime dotnet -Channel 6.0"
    exit 1
fi
ok "dotnet with .NET 6.0 runtime available"

if ! command -v git &>/dev/null; then
    err "git is required but not found."
    exit 1
fi
ok "git available"

# ---------------------------------------------------------------------------
# Section 2: Directory Setup
# ---------------------------------------------------------------------------

log "Setting up directories..."
mkdir -p "$DATA_DIR" "$NIRVANA_DATA" "$OUTPUT_DIR"
ok "Directories ready"

# ---------------------------------------------------------------------------
# Section 3: Clone & Build Nirvana
# ---------------------------------------------------------------------------

if [ ! -d "$NIRVANA_DIR/.git" ]; then
    log "Cloning Nirvana..."
    git clone https://github.com/Illumina/Nirvana.git "$NIRVANA_DIR"
else
    ok "Nirvana already cloned"
fi

# Checkout latest release tag for stability
cd "$NIRVANA_DIR"
LATEST_TAG="$(git tag --sort=-v:refname | head -1)"
if [ -n "$LATEST_TAG" ]; then
    CURRENT="$(git describe --tags --exact-match 2>/dev/null || echo '')"
    if [ "$CURRENT" != "$LATEST_TAG" ]; then
        log "Checking out latest release tag: $LATEST_TAG"
        git checkout "$LATEST_TAG" --quiet
    else
        ok "Already on latest tag: $LATEST_TAG"
    fi
else
    warn "No release tags found — using default branch"
fi

# Patch Nirvana's banner display to truncate the git SHA from the version
# string. Without this, the Downloader crashes because author + version
# exceeds the hardcoded 75-char line width in DisplayBanner().
BANNER_SRC="$NIRVANA_DIR/CommandLine/Utilities/CommandLineUtilities.cs"
if grep -q 'InformationalVersion.Length' "$BANNER_SRC" 2>/dev/null; then
    log "Patching Nirvana banner to handle long version strings..."
    sed -i '' \
        -e '/const int lineLength = 75;/i\
            // Truncate version at + to drop git SHA\
            string version = InformationalVersion;\
            int plusIdx = version.IndexOf('"'"'+'"'"');\
            if (plusIdx > 0) version = version.Substring(0, plusIdx);\
' \
        -e 's/lineLength - author\.Length - InformationalVersion\.Length/lineLength - author.Length - version.Length/' \
        -e 's/author, filler2, InformationalVersion/author, filler2, version/' \
        "$BANNER_SRC"
    # Force rebuild after patch
    rm -f "$NIRVANA_DLL"
    ok "Banner patched"
fi

if [ ! -f "$NIRVANA_DLL" ]; then
    log "Building Nirvana (Release)..."
    dotnet build -c Release
else
    ok "Nirvana already built"
fi

if [ ! -f "$NIRVANA_DLL" ]; then
    err "Nirvana build failed — expected binary at $NIRVANA_DLL"
    exit 1
fi
ok "Nirvana binary ready: $NIRVANA_DLL"

# Rebuild the native BlockCompression library for Apple Silicon if needed.
# Nirvana ships an x86_64-only dylib that crashes on arm64.
BLOCK_LIB="$NIRVANA_DIR/bin/Release/net6.0/libBlockCompression.dylib"
if [ "$(uname -m)" = "arm64" ] && file "$BLOCK_LIB" 2>/dev/null | grep -q x86_64; then
    log "Rebuilding libBlockCompression.dylib for arm64..."
    bash "$SCRIPT_DIR/build_block_compression.sh"
fi

cd "$PROJECT_ROOT"

# ---------------------------------------------------------------------------
# Section 4: Download Annotation Data (~30 GB)
# ---------------------------------------------------------------------------

# Reference file as existence check for completed download
REF_FILE="$NIRVANA_DATA/References/Homo_sapiens.${GENOME_ASSEMBLY}.Nirvana.dat"

if [ -f "$REF_FILE" ]; then
    ok "Nirvana annotation data already downloaded for $GENOME_ASSEMBLY"
else
    # Disk space check — require 40 GB free minimum
    AVAILABLE_GB=$(df -g "$DATA_DIR" | awk 'NR==2 {print $4}')
    if [ "$AVAILABLE_GB" -lt 40 ]; then
        err "Insufficient disk space: ${AVAILABLE_GB} GB available, 40 GB required"
        err "Nirvana annotation data is ~30 GB. Free up space and re-run."
        exit 1
    fi

    log "Downloading Nirvana annotation data for $GENOME_ASSEMBLY (~30 GB)..."
    log "This may take 30–60 minutes depending on your connection."
    dotnet "$DOWNLOADER_DLL" --ga "$GENOME_ASSEMBLY" --out "$NIRVANA_DATA"
    ok "Annotation data downloaded"
fi

# ---------------------------------------------------------------------------
# Section 5: Download Test VCF & Run Nirvana
# ---------------------------------------------------------------------------

INPUT_VCF="$DATA_DIR/HiSeq.10000.vcf.gz"
OUTPUT_JSON="$OUTPUT_DIR/HiSeq.10000.json.gz"

if [ ! -f "$INPUT_VCF" ]; then
    log "Downloading HiSeq.10000.vcf.gz test file..."
    curl -L -o "$INPUT_VCF" \
        "https://illumina.github.io/NirvanaDocumentation/files/HiSeq.10000.vcf.gz"
    ok "Test VCF downloaded"
else
    ok "Test VCF already present"
fi

if [ ! -f "$OUTPUT_JSON" ]; then
    # Nirvana's -c flag takes a path PREFIX — it appends .transcripts.ndb, etc.
    # Downloader places files at Cache/<assembly>/Both.*.ndb
    CACHE_PREFIX="$NIRVANA_DATA/Cache/$GENOME_ASSEMBLY/Both"

    log "Running Nirvana on HiSeq.10000.vcf.gz..."
    dotnet "$NIRVANA_DLL" \
        -c "$CACHE_PREFIX" \
        --sd "$NIRVANA_DATA/SupplementaryAnnotation/$GENOME_ASSEMBLY" \
        -r "$REF_FILE" \
        -i "$INPUT_VCF" \
        -o "$OUTPUT_DIR/HiSeq.10000"
    ok "Nirvana annotation complete"
else
    ok "Nirvana JSON output already exists"
fi

if [ ! -f "$OUTPUT_JSON" ]; then
    err "Nirvana output not found at $OUTPUT_JSON"
    exit 1
fi

# ---------------------------------------------------------------------------
# Section 6: Smoke Test with nirvana2vcf
# ---------------------------------------------------------------------------

SMOKE_VCF="$OUTPUT_DIR/HiSeq.10000.smoke.vcf"

log "Running nirvana2vcf smoke test..."

# Ensure nirvana2vcf is available (install in dev mode if needed)
if ! command -v nirvana2vcf &>/dev/null; then
    log "Installing nirvana2vcf in dev mode..."
    pip install -e "$PROJECT_ROOT" --quiet
fi

nirvana2vcf -i "$OUTPUT_JSON" -o "$SMOKE_VCF"

VARIANT_COUNT=$(grep -cv '^#' "$SMOKE_VCF" || true)
ok "Smoke test complete: $VARIANT_COUNT variants in output VCF"

# ---------------------------------------------------------------------------
# Section 7: Update Status & Summary
# ---------------------------------------------------------------------------

# Update Phase 1 status in strategy.md
if grep -q '| 1\. Generate Nirvana JSON | Not started |' "$STRATEGY_MD" 2>/dev/null; then
    log "Updating Phase 1 status in strategy.md..."
    sed -i '' 's/| 1\. Generate Nirvana JSON | Not started |/| 1. Generate Nirvana JSON | Done |/' "$STRATEGY_MD"
    ok "strategy.md updated"
else
    ok "strategy.md already up to date"
fi

# Print summary
NIRVANA_VERSION="$(dotnet "$NIRVANA_DLL" --version 2>&1 | head -1 || echo 'unknown')"
INPUT_SIZE="$(du -h "$INPUT_VCF" | cut -f1)"
OUTPUT_SIZE="$(du -h "$OUTPUT_JSON" | cut -f1)"
SMOKE_SIZE="$(du -h "$SMOKE_VCF" | cut -f1)"

echo ""
echo "============================================="
echo "  Phase 1 Complete"
echo "============================================="
echo "  Assembly:         $GENOME_ASSEMBLY"
echo "  Nirvana version:  $NIRVANA_VERSION"
echo "  Nirvana tag:      ${LATEST_TAG:-default branch}"
echo ""
echo "  Input VCF:        $INPUT_VCF ($INPUT_SIZE)"
echo "  Nirvana JSON:     $OUTPUT_JSON ($OUTPUT_SIZE)"
echo "  Smoke test VCF:   $SMOKE_VCF ($SMOKE_SIZE)"
echo "  Variant count:    $VARIANT_COUNT"
echo "============================================="
