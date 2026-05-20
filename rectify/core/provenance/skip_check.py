#!/usr/bin/env python3
"""
Skip-check decision tree for stage-level resume.

For a given stage of ``rectify run-all``, decide whether the prior run's
sidecar is valid enough to skip recomputing this stage, or whether
something has changed and we need to RUN.

Implements the §6.5.3 decision tree from
``dev/specs/parallel_refactor_plan.md``.

Correctness invariants
----------------------
- **Transient sub_outputs are NEVER sha256-checked.** Region BAMs in
  ``$L_SCRATCH``/``$TMPDIR`` are job-ephemeral by design; missing them is
  expected on re-run, not a reason to invalidate.
- **Persistent outputs ARE sha256-checked.** If a recorded sha256 mismatches
  the current file's sha256, RUN. Never trust mtime — files can be touched
  without content change, and content can change without mtime update.
- **Inputs are checked by role + sha256.** The path may have moved; the
  hash is the truth. If the caller didn't supply a role this run, we
  attempt to resolve from the prior sidecar's recorded path; if THAT
  doesn't exist either, RUN.
- **argv_template + config_hash are stable across ``$L_SCRATCH``/``$TMPDIR``
  drift.** The current run's argv is tokenized via
  :func:`tokenize_argv_paths` before hashing, then compared against the
  prior sidecar's stored ``config_hash``.
- **rectify_git_sha mismatch → RUN** unless ``accept_prior_provenance``
  is True. Default behavior is conservative: any change in rectify source
  invalidates prior outputs.

One wrong branch in this decision tree silently corrupts resume — either by
skipping a stage that should rerun (stale outputs) or by rerunning a stage
that should skip (wasted compute). This module is Opus-owned per
``dev/specs/parallel_refactor_plan.md`` §11.1.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

from .cluster import detect_current_cluster
from .hashing import normalized_config_hash, sha256_of_file
from .path_resolver import (
    PathResolutionError,
    PortablePath,
    tokenize_argv_paths,
)

if TYPE_CHECKING:
    # sidecar.py is Sonnet-owned (Commit A.5 boilerplate). We avoid a
    # runtime import to break the would-be cycle (sidecar will import
    # path_resolver, which is in the same package). At runtime, we import
    # ``read_stage_sidecar`` lazily inside ``can_skip_stage``.
    from .sidecar import ProvenanceRecord

logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------
# Return type
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class SkipDecision:
    """The result of :func:`can_skip_stage`.

    Two factory classmethods (``SKIP``, ``RUN``) build the two variants;
    callers branch on ``.skip``. When SKIP, ``.prior`` carries the parsed
    prior sidecar so the caller can chain outputs into the next stage's
    inputs without re-hashing.
    """

    skip: bool
    reason: str
    prior: Optional["ProvenanceRecord"] = None

    @classmethod
    def SKIP(cls, *, reason: str, prior: "ProvenanceRecord") -> "SkipDecision":
        return cls(skip=True, reason=reason, prior=prior)

    @classmethod
    def RUN(cls, *, reason: str) -> "SkipDecision":
        return cls(skip=False, reason=reason, prior=None)

    def __repr__(self) -> str:
        tag = "SKIP" if self.skip else "RUN"
        return f"SkipDecision({tag}, reason={self.reason!r})"


# ----------------------------------------------------------------------
# Decision tree
# ----------------------------------------------------------------------


def can_skip_stage(
    *,
    stage: str,
    sample_output: Union[Path, str],
    sample_id: str,
    current_inputs: dict,            # role -> Path (or str path)
    current_argv: list[str],
    rectify_git_sha: str,
    force: bool = False,
    force_stages: Optional[set] = None,
    accept_prior_provenance: bool = False,
) -> SkipDecision:
    """Decide whether to skip ``stage`` based on its prior sidecar.

    Returns :class:`SkipDecision.SKIP` if and only if ALL of these hold:

    1. Prior sidecar exists at
       ``<sample_output>/<sample_id>.<stage>.provenance.json``.
    2. Prior sidecar's ``exit_status == 0``.
    3. All declared NON-TRANSIENT outputs exist on disk and their sha256
       matches what the prior sidecar recorded. (Optional outputs are
       allowed to be absent.)
    4. All declared inputs (matched by role, NOT by path) sha256-match the
       prior sidecar.
    5. Current argv's normalized ``config_hash`` matches the prior.
    6. Prior ``rectify_git_sha == rectify_git_sha`` (current run's), OR
       ``accept_prior_provenance`` is True.

    Otherwise returns :class:`SkipDecision.RUN` with a specific reason that
    can be surfaced to the user in a ``[RESUME]`` log line.

    Args:
        stage: e.g. ``"correct"``, ``"analyze"``, ``"drs_trim"``.
        sample_output: the sample's output directory (where sidecars live).
        sample_id: the sample's id; combined with ``stage`` to form the
            sidecar filename.
        current_inputs: mapping ``role -> Path`` of inputs the caller has
            for this run. Roles missing from this dict will be resolved
            from the prior sidecar's recorded path (with the usual
            PortablePath fallback chain) if needed.
        current_argv: typically ``sys.argv``.
        rectify_git_sha: current rectify git sha (e.g. from
            ``git log -1 --format=%H``).
        force: if True, return RUN unconditionally with reason
            ``"--force-all"``.
        force_stages: set of stage names to force-rerun even when
            ``force`` is False (corresponds to ``--force-stage``).
        accept_prior_provenance: if True, ignore git_sha mismatch
            (corresponds to ``--accept-prior-provenance``).

    Returns:
        :class:`SkipDecision`.
    """
    sample_output = Path(sample_output)
    force_stages = force_stages or set()

    # GATE 0: --force-* takes precedence over everything.
    if force:
        return SkipDecision.RUN(reason="--force-all")
    if stage in force_stages:
        return SkipDecision.RUN(reason=f"--force-stage {stage}")

    sidecar_path = sample_output / f"{sample_id}.{stage}.provenance.json"
    if not sidecar_path.exists():
        return SkipDecision.RUN(reason=f"no prior sidecar at {sidecar_path.name}")

    # Lazy import to defer the circular-import risk. By the time
    # can_skip_stage runs, sidecar.py will have been loaded (it's part of
    # the same package).
    from .sidecar import read_stage_sidecar

    try:
        prior = read_stage_sidecar(sidecar_path)
    except Exception as e:
        return SkipDecision.RUN(
            reason=f"sidecar at {sidecar_path.name} unreadable: {e}"
        )

    if prior.exit_status != 0:
        return SkipDecision.RUN(
            reason=f"prior run exited with status {prior.exit_status}"
        )

    # Cross-cluster warning (non-blocking; surfaces here, not at end).
    current_cluster = detect_current_cluster()
    if prior.cluster and prior.cluster != current_cluster:
        logger.warning(
            "[RESUME] cluster mismatch for stage=%s: sidecar was written "
            "on %s, current host is %s; env vars may resolve differently "
            "(expect cached_absolute fallback warnings).",
            stage, prior.cluster, current_cluster,
        )

    # GATE 1: declared NON-TRANSIENT outputs exist + sha256 match.
    # Sub-outputs (region BAMs in $L_SCRATCH / $TMPDIR with transient=True)
    # are intentionally NOT checked here — they're job-ephemeral and may
    # legitimately not exist on the current host.
    for out in prior.outputs:
        try:
            p = out.path.resolve_now(sample_output)
        except PathResolutionError as e:
            return SkipDecision.RUN(
                reason=f"output {out.role!r} cannot be resolved on current host: {e}"
            )
        if not p.exists():
            # `resolve_now` already tried cached_absolute fallback. If we
            # got here, it returned a missing path — shouldn't happen, but
            # be defensive.
            if getattr(out, "optional", False):
                continue
            return SkipDecision.RUN(
                reason=f"output {out.role!r} missing on disk: {p}"
            )
        try:
            current_sha = "sha256:" + sha256_of_file(p)
        except Exception as e:
            return SkipDecision.RUN(
                reason=f"output {out.role!r} cannot be hashed ({p}): {e}"
            )
        if current_sha != out.sha256:
            return SkipDecision.RUN(
                reason=f"output {out.role!r} sha256 changed since prior run"
            )

    # GATE 2: inputs match by role + sha256.
    ignored_roles = set(prior.skip_check_config.get("ignore_input_roles", []))
    for inp in prior.inputs:
        if inp.role in ignored_roles:
            continue

        current_path = current_inputs.get(inp.role)
        if current_path is None:
            # Caller didn't supply this role for the current run. Fall
            # back to the path the prior sidecar recorded.
            try:
                current_path = inp.path.resolve_now(sample_output)
            except PathResolutionError as e:
                return SkipDecision.RUN(
                    reason=f"input {inp.role!r} not supplied this run and "
                           f"cannot be resolved from prior sidecar: {e}"
                )

        current_path = Path(current_path)
        if not current_path.exists():
            return SkipDecision.RUN(
                reason=f"input {inp.role!r} path does not exist: {current_path}"
            )
        try:
            current_sha = "sha256:" + sha256_of_file(current_path)
        except Exception as e:
            return SkipDecision.RUN(
                reason=f"input {inp.role!r} cannot be hashed ({current_path}): {e}"
            )
        if current_sha != inp.sha256:
            return SkipDecision.RUN(
                reason=f"input {inp.role!r} sha256 changed (path={current_path})"
            )

    # GATE 3: argv config_hash matches.
    ignored_argv = set(prior.skip_check_config.get("ignore_argv", []))
    current_template = tokenize_argv_paths(current_argv, sample_output)
    current_hash = normalized_config_hash(current_template, ignore=ignored_argv)
    prior_hash = prior.invocation.get("config_hash")
    if prior_hash is None:
        return SkipDecision.RUN(
            reason="prior sidecar missing invocation.config_hash"
        )
    if current_hash != prior_hash:
        return SkipDecision.RUN(
            reason="invocation argv differs (use --dry-run-resume to see diff)"
        )

    # GATE 4: rectify git_sha matches (or override).
    if prior.rectify_git_sha != rectify_git_sha and not accept_prior_provenance:
        return SkipDecision.RUN(
            reason=(
                f"rectify git_sha changed: "
                f"{prior.rectify_git_sha[:7]} -> {rectify_git_sha[:7]} "
                f"(use --accept-prior-provenance to override)"
            )
        )

    # All gates passed.
    wall = prior.stats.get("wall_seconds", "?") if prior.stats else "?"
    return SkipDecision.SKIP(
        reason=(
            f"prior valid run completed_at={prior.completed_at}, "
            f"wall_seconds={wall}"
        ),
        prior=prior,
    )
