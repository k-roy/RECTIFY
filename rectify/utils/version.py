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
    """Identify the code that produced an output. NEVER returns a bare 'unknown'
    when anything better is obtainable.

    This is THE field that pins which rectify produced a dataset, so it is
    mission-critical that it resolve. Preference order (see
    ``rectify.core.provenance.sidecar.resolve_code_fingerprint``):

      1. ``git:<40-hex>[+dirty]`` — a real git checkout. **This is the goal.**
         ``+dirty`` means the work tree had uncommitted changes, so the sha
         alone does NOT reproduce the run.
      2. ``pkghash:<12-hex>`` — sha256 over every ``.py`` in the installed
         package. Backstop for a deployment that is a plain COPY of the repo.
         It still moves whenever the code moves, so drift is detectable — but
         it is NOT resolvable back to a commit. A WARNING is logged.
      3. ``"unknown"`` — only if even hashing fails. Logged as an ERROR.

    History (2026-08-02): this used to return ``"unknown"`` whenever no
    enclosing ``.git`` existed, and every rectify deployment on Hoffman2 was a
    non-git copy. Result: all 33 provenance sidecars behind the MS2 3'-end
    tables recorded ``rectify_git_sha: "unknown"``, and a ~2-month-stale rectify
    produced production data with nothing flagging it. Silent 'unknown' is the
    exact failure mode this must never repeat.

    ⚠️ Deployments SHOULD be git checkouts. ``pkghash`` proves *that* the code
    differs, not *which* code it is; only a git sha is reproducible.
    """
    from rectify.core.provenance.sidecar import resolve_code_fingerprint
    return resolve_code_fingerprint()
