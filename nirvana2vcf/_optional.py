"""Lazy imports for optional pysam-backed features."""

import sys


def require_pysam(feature: str):
    try:
        import pysam
        return pysam
    except ImportError:
        print(
            f"Error: {feature} requires pysam. "
            f"Install with: pip install 'nirvana2vcf[full]'",
            file=sys.stderr,
        )
        sys.exit(1)
