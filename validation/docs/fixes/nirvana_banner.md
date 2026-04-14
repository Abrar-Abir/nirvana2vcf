# Nirvana Banner Crash Fix

## Problem

When Nirvana is built from source and then run (either `Nirvana.dll` or `Downloader.dll`), the program crashes immediately with:

```
System.InvalidOperationException: Unable to display the program banner,
the author name and version string is too long.
```

This prevents the Downloader from fetching annotation data (Phase 1, Section 4).

## Root Cause

`CommandLine/Utilities/CommandLineUtilities.cs` has a `DisplayBanner()` method that formats author and version strings into a fixed 75-character-wide banner:

```csharp
const int lineLength = 75;
int fillerLength2 = lineLength - author.Length - InformationalVersion.Length;

if (fillerLength2 < 1)
    throw new InvalidOperationException("...");
```

When building from a git checkout, .NET automatically appends the full commit SHA to `InformationalVersion` via the `SourceRevisionId` MSBuild property:

| Component | Value | Length |
|-----------|-------|--------|
| Authors | `Stromberg, Roy, Platzer, Siddiqui, Ouyang, et al` | 50 |
| Version | `3.18.1+05f88047f4c9e79e54ed6b6482549258596bb146` | 47 |
| **Total** | | **97** |

97 > 75, so the banner throws before any useful work begins.

## Why MSBuild flags don't help

The obvious fix — `dotnet build -p:SourceRevisionId=` or `-p:InformationalVersion=3.18.1` — does not work. The .NET SDK's `InitializeSourceControlInformation` target repopulates `SourceRevisionId` from git after command-line properties are evaluated, and the SHA persists in the compiled assembly regardless.

## Fix Applied

`run_phase1.sh` patches the Nirvana source after cloning and before building. The sed command in Section 3 of the script:

1. Inserts a version-truncation block before the `lineLength` constant:
   ```csharp
   string version = InformationalVersion;
   int plusIdx = version.IndexOf('+');
   if (plusIdx > 0) version = version.Substring(0, plusIdx);
   ```
2. Replaces `InformationalVersion.Length` with `version.Length` in the filler calculation.
3. Replaces `InformationalVersion` with `version` in the `Console.WriteLine` call.

This strips the `+<sha>` suffix at display time only — the full version remains in assembly metadata.

### Idempotency

The patch is guarded by `grep -q 'InformationalVersion.Length'`, which matches only the unpatched source (the patch replaces that exact string with `version.Length`). Re-running the script after patching is a no-op.

## Affected Versions

Observed on Nirvana v3.18.1 (latest release tag as of 2026-04-13) built with .NET SDK 8.0.419 on macOS. Likely affects any source build where git history is present, since the SHA append is a .NET SDK default behavior.
