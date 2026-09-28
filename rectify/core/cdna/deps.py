"""Runtime check for edlib, which cDNA adapter and SSP anchoring needs.

edlib is a core dependency. This guards environments built some other way (a conda env without
it, a vendored tree). Without edlib every fuzzy anchor search in ``read_info`` and ``walkback``
falls back to EXACT matching: a molecule with one basecall error in its SSP is classified Type 2
(no UMI, no deduplication, adapter left as a 5' soft clip) and the poly(A) anchor search misreads
tail lengths. A default ``pip install`` shipped exactly that until 2026-09-28, silently.
"""
import logging


def edlib_available() -> bool:
    try:
        import edlib  # noqa: F401
    except ImportError:
        return False
    return True


def warn_if_edlib_missing(log: logging.Logger, command: str) -> bool:
    """Log a WARNING that names the consequence when edlib is absent. Returns availability."""
    if edlib_available():
        return True
    log.warning(
        "%s: edlib is NOT installed, so SSP and adapter anchoring fall back to EXACT matching. "
        "Molecules with a basecall error in the SSP become Type 2 (no UMI, no deduplication, "
        "adapter left as a 5' soft clip) and poly(A) tail lengths are misread. "
        "Install it (pip install edlib) and rerun.", command)
    return False
