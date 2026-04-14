# Nirvana BlockCompression Apple Silicon Fix

## Problem

After building Nirvana from source and patching the banner crash (see `nirvana_banner.md`), running Nirvana or the Downloader on Apple Silicon (arm64) fails immediately with:

```
ERROR: Unable to find the block GZip compression library (BlockCompression)
```

This prevents Nirvana from reading any annotation data files, which are BGZF-compressed.

## Root Cause

Nirvana uses a native shared library (`libBlockCompression.dylib`) for BGZF and Zstandard compression via .NET P/Invoke (`DllImport`). The library ships precompiled in the Nirvana source tree at `Compression/Packages/BlockCompression/` and is copied to the build output during `dotnet build`.

The shipped binary is **x86_64 only**:

```
$ file libBlockCompression.dylib
Mach-O 64-bit dynamically linked shared library x86_64
```

On Apple Silicon (arm64), the .NET runtime cannot load this x86_64 library. `LibraryUtilities.CheckLibrary()` wraps the native `get_library_id()` call in a try/catch — any failure (load error or wrong return value) throws `MissingCompressionLibraryException("BlockCompression")`, which Nirvana reports as the generic "Unable to find" message.

No C source code is included in the Nirvana repository for this library — only precompiled binaries for x86_64 Linux (.so) and macOS (.dylib).

## What the Library Provides

Nirvana's .NET code imports 6 functions from `BlockCompression`:

| Function | Source | Purpose |
|----------|--------|---------|
| `get_library_id()` | Custom | Returns magic `-822411574` (`0xCEFAFECA`) for load validation |
| `bgzf_compress()` | Custom (wraps zlib) | BGZF block compression (SAM/BAM spec) |
| `bgzf_decompress()` | Custom (wraps zlib) | BGZF block decompression |
| `ZSTD_compress()` | libzstd | Zstandard compression |
| `ZSTD_decompress()` | libzstd | Zstandard decompression |
| `ZSTD_getDecompressedSize()` | libzstd | Get decompressed size from frame header |

The original x86_64 library statically links zlib-ng, zstd, and libdeflate into a single binary. The BGZF functions implement the standard Block GZip Format: an 18-byte gzip header with `BC` extra field, raw deflate payload, and 8-byte trailer (CRC32 + uncompressed size).

## Fix Applied

`build_block_compression.sh` compiles a minimal C replacement that:

1. **Implements `get_library_id()`** — returns the expected magic value `-822411574` (`0xCEFAFECA`). Note: the Nirvana source comments this as "cafeface" but the actual signed int32 value is `0xCEFAFECA` — using the wrong hex literal (`0xCAFEFACE`) causes the validation check to fail silently.

2. **Implements `bgzf_compress()` and `bgzf_decompress()`** — standard BGZF block format using system zlib's `deflateInit2`/`inflateInit2` with raw deflate (`windowBits = -15`). Matches the P/Invoke signatures exactly:
   ```c
   int bgzf_compress(uint8_t *dst, int dstLen,
                     const uint8_t *src, int srcLen, int level);
   int bgzf_decompress(uint8_t *dst, int dstSize,
                       const uint8_t *src, int srcSize);
   ```

3. **Statically links libzstd** via `-Wl,-force_load` so that `ZSTD_compress`, `ZSTD_decompress`, and `ZSTD_getDecompressedSize` are exported directly from the dylib (matching the original library's symbol layout).

### Build command

```bash
cc -shared -O2 -arch arm64 \
    -o libBlockCompression.dylib \
    block_compression.c \
    -lz \
    -I"$(brew --prefix zstd)/include" \
    -Wl,-force_load,"$(brew --prefix zstd)/lib/libzstd.a"
```

### Dependencies

- **zlib**: Linked dynamically (`-lz`). Always present on macOS.
- **zstd**: Linked statically from Homebrew (`brew install zstd`). The static archive is force-loaded so all ZSTD symbols appear as exports.
- **Xcode Command Line Tools**: Provides `cc` (Apple Clang).

## Integration with run_phase1.sh

`run_phase1.sh` detects the architecture mismatch automatically after building Nirvana:

```bash
BLOCK_LIB="$NIRVANA_DIR/bin/Release/net6.0/libBlockCompression.dylib"
if [ "$(uname -m)" = "arm64" ] && file "$BLOCK_LIB" | grep -q x86_64; then
    bash "$SCRIPT_DIR/build_block_compression.sh"
fi
```

This is idempotent — once the arm64 dylib is in place, the `grep -q x86_64` check fails and the build is skipped on subsequent runs.

## Verification

After building, the library exports all required symbols:

```
$ nm -g libBlockCompression.dylib | grep -E 'get_library|bgzf|ZSTD_(com|decom|getDe)'
00000000000004e8 T _get_library_id
000000000000067c T _bgzf_compress
00000000000007b4 T _bgzf_decompress
000000000000a970 T _ZSTD_compress
0000000000052c8c T _ZSTD_decompress
0000000000052130 T _ZSTD_getDecompressedSize
```

Nirvana loads and annotates correctly:

```
$ dotnet Nirvana.dll -c ... --sd ... -r ... -i HiSeq.10000.vcf.gz -o output
---------------------------------------------------------------------------
Nirvana                                             (c) 2022 Illumina, Inc.
Stromberg, Roy, Platzer, Siddiqui, Ouyang, et al                     3.18.1
---------------------------------------------------------------------------

Initialization                                         Time     Positions/s
---------------------------------------------------------------------------
Cache                                               00:00:00.6
SA Position Scan                                    00:00:00.0      298,541

Reference                                Preload    Annotation   Variants/s
---------------------------------------------------------------------------
chr1                                    00:00:00.1  00:00:00.9       10,531

Time: 00:00:03.3
```

### BGZF header gotcha

The standard BGZF header is 18 bytes. The third byte must be `0x08` (CM: deflate compression method). A common mistake is to omit this byte and place the FLG byte (`0x04`) in the CM position, which shifts the entire header and produces files that gzip/bgzip/Python's `gzip` module reject with "Unknown compression method". The correct byte sequence starts:

```
1f 8b 08 04 ...
^^    ^^ ^^
magic CM FLG(FEXTRA)
```

### Magic value gotcha

The Nirvana source comments the library ID constant as "cafeface":

```csharp
const int expectedLibraryId = -822411574; // cafeface
```

This is misleading. The signed int32 `-822411574` is `0xCEFAFECA` in hex, not `0xCAFEFACE`. Using the literal `0xCAFEFACE` in C produces `-889259314`, which fails the validation. The correct C implementation must use the decimal literal directly:

```c
int get_library_id(void) { return -822411574; }
```

## Affected Versions

Observed on Nirvana v3.18.1 (latest release as of 2026-04-13) on macOS with Apple Silicon (M-series). Affects any arm64 macOS system building Nirvana from source. Not an issue on x86_64 macOS or Linux x86_64 (which has a matching `.so`).
