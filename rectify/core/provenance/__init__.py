"""Stage-level provenance sidecars + resume infrastructure.

This package is the MACHINE-READABLE skip-cache for ``rectify run-all``.
Each stage emits a sidecar JSON to
``<sample_output>/<sample_id>.<stage>.provenance.json`` on successful
completion. Subsequent run-all invocations consult sidecars via
:func:`can_skip_stage` to decide what to recompute.

Distinct from ``rectify.utils.provenance.ProvenanceTracker``, which is a
USER-FACING per-output-directory run log + README writer. The two coexist
intentionally; the legacy accumulator is unchanged. Code that imports from
``rectify.utils.provenance`` for run-log purposes should continue to do so.
Code that wants to read/write stage-level sidecars or implement skip-check
should import from this package (``rectify.core.provenance``).

Public API:
    ProvenanceRecord       — dataclass mirroring the v1.0 sidecar schema
    write_stage_sidecar    — atomic write of a sidecar
    read_stage_sidecar     — load + validate a sidecar
    can_skip_stage         — the skip-check decision tree (see spec §6.5.3)
    SkipDecision           — return type of can_skip_stage
    PortablePath           — path envelope (see spec §6.5.2.5)
    PathResolutionError    — raised by PortablePath.resolve_now
    tokenize_argv_paths    — replace absolute paths in argv with templates
    detect_current_cluster — returns "sherlock"/"hoffman2"/"m1"/"other"
    add_resume_args        — argparse helper for --force-*/--dry-run-resume
    sha256_of_file         — SHA-256 hex digest of a file
    normalized_config_hash — hash an argv list, ignoring selected flags
"""
from .sidecar import ProvenanceRecord, write_stage_sidecar, read_stage_sidecar
from .skip_check import can_skip_stage, SkipDecision
from .path_resolver import PortablePath, PathResolutionError, tokenize_argv_paths
from .cluster import detect_current_cluster
from .cli import add_resume_args
from .hashing import sha256_of_file, normalized_config_hash

__all__ = [
    "ProvenanceRecord",
    "write_stage_sidecar",
    "read_stage_sidecar",
    "can_skip_stage",
    "SkipDecision",
    "PortablePath",
    "PathResolutionError",
    "tokenize_argv_paths",
    "detect_current_cluster",
    "add_resume_args",
    "sha256_of_file",
    "normalized_config_hash",
]
