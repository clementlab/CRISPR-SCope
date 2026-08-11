"""I/O helpers shared across CRISPRSCope modules."""

from __future__ import annotations

import gzip
from pathlib import Path


_GZIP_MAGIC = b"\x1f\x8b"


def open_text_maybe_gzip(path, mode="rt"):
    """Open a text file that may be gzip-compressed regardless of suffix."""
    if mode != "rt":
        raise ValueError(f"Unsupported mode for open_text_maybe_gzip: {mode!r}")

    path = Path(path)
    with path.open("rb") as handle:
        magic = handle.read(2)

    if magic == _GZIP_MAGIC:
        return gzip.open(path, mode)
    return path.open(mode)
