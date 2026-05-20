"""Version helpers for RECTIFY.

Thin wrappers that normalise the return type and provide the canonical
``get_rectify_git_sha()`` function expected by stage-sidecar wiring
(Commit B, 2026-05-20).

The underlying implementation lives in ``rectify.utils.bam_provenance``.
This module re-exports it as a non-Optional form so sidecar callers always
get a str (falling back to ``"unknown"``).
"""

from __future__ import annotations


def get_rectify_git_sha() -> str:
    """Return current HEAD sha of the rectify install's git repo, or 'unknown'.

    Used in PROVENANCE sidecars for skip-check git_sha matching.
    Falls back to 'unknown' if rectify was pip-installed without a .git dir
    (e.g., wheel install) — in that case, sidecar's git_sha mismatch will
    always fire, which is the safe behavior.

    Implementation note: delegates to ``rectify.utils.bam_provenance.get_rectify_git_sha``
    which caches the result after the first subprocess call.
    """
    from rectify.utils.bam_provenance import get_rectify_git_sha as _impl
    sha = _impl()
    return sha if sha else "unknown"
