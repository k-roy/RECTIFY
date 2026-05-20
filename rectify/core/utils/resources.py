"""Resource-policy helpers for pipeline concurrency decisions."""

import os
from typing import Optional


def _detect_total_mem_gb() -> Optional[float]:
    """Best-effort physical memory detection without optional dependencies."""
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
    except (AttributeError, ValueError, OSError):
        return None
    if pages <= 0 or page_size <= 0:
        return None
    return pages * page_size / 1024**3


def detect_machine_class(
    total_mem_gb: Optional[float] = None,
    cpu_count: Optional[int] = None,
) -> str:
    """Return ``m1_laptop`` or ``cluster_node`` for concurrency policy."""
    if total_mem_gb is None:
        total_mem_gb = _detect_total_mem_gb()
    if cpu_count is None:
        cpu_count = os.cpu_count()

    if (
        total_mem_gb is not None
        and cpu_count is not None
        and total_mem_gb < 16
        and cpu_count < 8
    ):
        return "m1_laptop"
    return "cluster_node"


def resolve_aligner_concurrency(
    value: str,
    total_threads: int,
    machine_class: str,
) -> int:
    """Resolve ``--aligner-concurrency`` to a bounded worker budget."""
    if total_threads < 1:
        total_threads = 1

    if value == "auto":
        if machine_class == "m1_laptop":
            return 1
        return max(1, total_threads - 2)

    try:
        requested = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "--aligner-concurrency must be 'auto' or a positive integer"
        ) from exc

    if requested < 1:
        raise ValueError("--aligner-concurrency must be >= 1")
    return min(requested, total_threads)

