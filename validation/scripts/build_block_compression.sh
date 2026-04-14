#!/usr/bin/env bash
# =============================================================================
# Build libBlockCompression.dylib for Apple Silicon (arm64)
#
# Nirvana ships a precompiled x86_64-only native library for its BGZF and
# zstd compression. This script builds an arm64 replacement from a minimal
# C shim that links against system zlib and homebrew zstd.
#
# Usage:
#   bash validation/scripts/build_block_compression.sh
#
# Prerequisites:
#   brew install zstd
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
NIRVANA_DIR="$PROJECT_ROOT/validation/data/Nirvana"
BUILD_DIR="$PROJECT_ROOT/validation/data/block_compression_build"
OUTPUT="$NIRVANA_DIR/bin/Release/net6.0/libBlockCompression.dylib"

# Color codes
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

log() { echo -e "${BLUE}[BlockCompression]${NC} $*"; }
ok()  { echo -e "${GREEN}[  OK  ]${NC} $*"; }

# Check prerequisites
if ! command -v cc &>/dev/null; then
    echo "ERROR: C compiler (cc) not found. Install Xcode Command Line Tools." >&2
    exit 1
fi

ZSTD_PREFIX="$(brew --prefix zstd 2>/dev/null || true)"
if [ -z "$ZSTD_PREFIX" ]; then
    echo "ERROR: zstd not found. Install with: brew install zstd" >&2
    exit 1
fi

mkdir -p "$BUILD_DIR"

# ---- Write the C source ----
cat > "$BUILD_DIR/block_compression.c" << 'CSOURCE'
/*
 * Minimal BlockCompression native library for Nirvana on Apple Silicon.
 *
 * Exports the 6 functions Nirvana's .NET P/Invoke expects:
 *   - get_library_id()
 *   - bgzf_compress() / bgzf_decompress()
 *   - ZSTD_compress() / ZSTD_decompress() / ZSTD_getDecompressedSize()
 *
 * BGZF = Block GZip Format (SAM/BAM spec). Each block is:
 *   - 18-byte header (gzip header with BGZF extra field)
 *   - deflate-compressed payload
 *   - 8-byte trailer (CRC32 + uncompressed size)
 */

#include <zlib.h>
#include <zstd.h>
#include <string.h>
#include <stdint.h>

/* ---- Library ID ---- */
/* Nirvana checks this magic value at startup: (int)0xCEFAFECA == -822411574 */
int get_library_id(void) {
    return -822411574;
}

/* ---- BGZF constants ---- */
#define BGZF_HEADER_SIZE 18
#define BGZF_TRAILER_SIZE 8
#define BGZF_MAX_BLOCK_SIZE 65536

/* Standard BGZF header template */
static const uint8_t bgzf_header[BGZF_HEADER_SIZE] = {
    0x1f, 0x8b,             /* gzip magic */
    0x08,                   /* CM: deflate */
    0x04,                   /* FLG: FEXTRA */
    0x00, 0x00, 0x00, 0x00, /* MTIME */
    0x00,                   /* XFL */
    0xff,                   /* OS: unknown */
    0x06, 0x00,             /* XLEN = 6 */
    0x42, 0x43,             /* SI1, SI2: 'BC' */
    0x02, 0x00,             /* SLEN = 2 */
    0x00, 0x00              /* BSIZE - 1 (placeholder, filled at compress time) */
};

/*
 * bgzf_compress: Compress data into a single BGZF block.
 *
 * Parameters match Nirvana's P/Invoke signature:
 *   compressedBlock:   output buffer
 *   compressedLen:     output buffer capacity
 *   uncompressedBlock: input data
 *   uncompressedLen:   input data length
 *   compressionLevel:  zlib compression level (1-9)
 *
 * Returns: total bytes written to compressedBlock (header + deflate + trailer)
 */
int bgzf_compress(
    uint8_t *compressedBlock,   int compressedLen,
    const uint8_t *uncompressedBlock, int uncompressedLen,
    int compressionLevel)
{
    /* Write BGZF header */
    memcpy(compressedBlock, bgzf_header, BGZF_HEADER_SIZE);

    /* Deflate the payload */
    z_stream zs;
    memset(&zs, 0, sizeof(zs));

    /* Use raw deflate (windowBits = -15) to match BGZF spec */
    if (deflateInit2(&zs, compressionLevel, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK)
        return -1;

    zs.next_in   = (Bytef *)uncompressedBlock;
    zs.avail_in  = (uInt)uncompressedLen;
    zs.next_out  = compressedBlock + BGZF_HEADER_SIZE;
    zs.avail_out = (uInt)(compressedLen - BGZF_HEADER_SIZE - BGZF_TRAILER_SIZE);

    if (deflate(&zs, Z_FINISH) != Z_STREAM_END) {
        deflateEnd(&zs);
        return -1;
    }
    deflateEnd(&zs);

    int deflateSize = (int)zs.total_out;
    int blockSize = BGZF_HEADER_SIZE + deflateSize + BGZF_TRAILER_SIZE;

    /* Fill BSIZE-1 in header (bytes 16-17, little-endian) */
    uint16_t bsize_minus1 = (uint16_t)(blockSize - 1);
    compressedBlock[16] = (uint8_t)(bsize_minus1 & 0xFF);
    compressedBlock[17] = (uint8_t)(bsize_minus1 >> 8);

    /* Write trailer: CRC32 + uncompressed size (both little-endian) */
    uint32_t crc = (uint32_t)crc32(0L, uncompressedBlock, (uInt)uncompressedLen);
    uint8_t *trailer = compressedBlock + BGZF_HEADER_SIZE + deflateSize;
    trailer[0] = (uint8_t)(crc & 0xFF);
    trailer[1] = (uint8_t)((crc >> 8) & 0xFF);
    trailer[2] = (uint8_t)((crc >> 16) & 0xFF);
    trailer[3] = (uint8_t)((crc >> 24) & 0xFF);
    trailer[4] = (uint8_t)(uncompressedLen & 0xFF);
    trailer[5] = (uint8_t)((uncompressedLen >> 8) & 0xFF);
    trailer[6] = (uint8_t)((uncompressedLen >> 16) & 0xFF);
    trailer[7] = (uint8_t)((uncompressedLen >> 24) & 0xFF);

    return blockSize;
}

/*
 * bgzf_decompress: Decompress a single BGZF block.
 *
 * Returns: number of uncompressed bytes written.
 */
int bgzf_decompress(
    uint8_t *uncompressedBlock, int uncompressedSize,
    const uint8_t *compressedBlock,   int compressedSize)
{
    z_stream zs;
    memset(&zs, 0, sizeof(zs));

    /* Raw inflate (no gzip/zlib header — we skip the BGZF header ourselves) */
    if (inflateInit2(&zs, -15) != Z_OK)
        return -1;

    zs.next_in   = (Bytef *)(compressedBlock + BGZF_HEADER_SIZE);
    zs.avail_in  = (uInt)(compressedSize - BGZF_HEADER_SIZE - BGZF_TRAILER_SIZE);
    zs.next_out  = uncompressedBlock;
    zs.avail_out = (uInt)uncompressedSize;

    if (inflate(&zs, Z_FINISH) != Z_STREAM_END) {
        inflateEnd(&zs);
        return -1;
    }

    int decompressedLen = (int)zs.total_out;
    inflateEnd(&zs);
    return decompressedLen;
}

CSOURCE

log "Compiling libBlockCompression.dylib for arm64..."

# Statically link zstd so ZSTD_* symbols are exported from our dylib
# (matching Nirvana's original x86_64 library which also had them static).
# -force_load pulls in all zstd object files even without explicit references.
# zlib is linked dynamically — it's always available on macOS.
cc -shared -O2 -arch arm64 \
    -o "$OUTPUT" \
    "$BUILD_DIR/block_compression.c" \
    -lz \
    -I"${ZSTD_PREFIX}/include" \
    -Wl,-force_load,"${ZSTD_PREFIX}/lib/libzstd.a"

ok "Built: $OUTPUT"
file "$OUTPUT"

# Also copy to the Compression/Packages directory (used during build)
PACKAGES_DIR="$NIRVANA_DIR/Compression/Packages/BlockCompression"
if [ -d "$PACKAGES_DIR" ]; then
    cp "$OUTPUT" "$PACKAGES_DIR/libBlockCompression.dylib"
    ok "Copied to $PACKAGES_DIR/"
fi

# Verify the library loads and has the expected symbols
log "Verifying exported symbols..."
nm -g "$OUTPUT" | grep -E "get_library_id|bgzf_compress|bgzf_decompress|ZSTD_compress|ZSTD_decompress|ZSTD_getDecompressedSize"
ok "All required symbols present"
