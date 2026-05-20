"""SHA-256 of files + argv config hashing.

File sha256 uses mmap for files > 16 MB (fast path). Smaller files use
chunked read. Returns hex digest.

``normalized_config_hash`` takes an argv list and an ignore-list of flag
names; filters out the ignored flags + their values; sorts the result; hashes
it. The result is stable across runs that differ only in path-related or
threading flags.
"""
import hashlib
import mmap
import os
from pathlib import Path
from typing import List, Optional, Set


def sha256_of_file(path: "Path | str", *, chunk_size: int = 1 << 20) -> str:
    """Compute SHA-256 hex digest of a file.

    For files > 16 MB, mmap-based read; otherwise chunked.

    Args:
        path: Path to the file.
        chunk_size: Chunk size for the non-mmap path (default 1 MiB).

    Returns:
        Hex digest string (no ``sha256:`` prefix — callers add the prefix
        when embedding in sidecar JSON).

    Raises:
        FileNotFoundError: if ``path`` doesn't exist.
    """
    path = Path(path)
    size = path.stat().st_size  # propagates FileNotFoundError
    h = hashlib.sha256()
    if size == 0:
        return h.hexdigest()
    if size > (16 << 20):
        with open(path, "rb") as fh:
            with mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                h.update(mm)
    else:
        with open(path, "rb") as fh:
            while True:
                buf = fh.read(chunk_size)
                if not buf:
                    break
                h.update(buf)
    return h.hexdigest()


def normalized_config_hash(
    argv: List[str],
    *,
    ignore: Optional[Set[str]] = None,
) -> str:
    """Hash an argv list after removing flags listed in ``ignore`` (with their values).

    Args:
        argv: Full argv (typically ``sys.argv`` or
            ``prior['invocation']['argv_template']``).
        ignore: Set of flag NAMES to filter out (e.g. ``{"--n-threads",
            "--tmp-dir"}``). For each flag, the following positional value
            (if any) is also dropped.

    Returns:
        ``"sha256:" + SHA-256 hex`` of the canonical-form normalized argv.

    Canonicalization:
        Filtered argv items are joined with U+001F (unit separator) so two
        distinct argvs with different splits can't accidentally collide.
        Order is NOT sorted — order is preserved. Two argvs that are
        otherwise identical but differ in flag order produce different hashes
        (intentional; preserves pipeline semantics).
    """
    ignore = ignore or set()
    filtered: List[str] = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        # ``--foo=bar`` form (one token)
        if "=" in arg and arg.startswith("--"):
            flag = arg.split("=", 1)[0]
            if flag in ignore:
                i += 1
                continue
            filtered.append(arg)
            i += 1
            continue
        # ``--foo bar`` form (two tokens)
        if arg in ignore:
            i += 1
            # Skip following non-flag token if present (the flag's value).
            if i < len(argv) and not argv[i].startswith("-"):
                i += 1
            continue
        filtered.append(arg)
        i += 1

    canonical = "\x1f".join(filtered)
    return "sha256:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()
