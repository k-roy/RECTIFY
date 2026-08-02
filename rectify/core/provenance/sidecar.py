"""Stage-level provenance sidecar — dataclass + atomic read/write.

Schema version: 1.0

Each stage writes one sidecar JSON to::

    <sample_output>/<sample_id>.<stage>.provenance.json

on successful completion. The JSON follows the schema defined in
``dev/specs/parallel_refactor_plan.md §6.5.2``.

All path fields store :class:`~rectify.core.provenance.path_resolver.PortablePath`
envelopes, NEVER bare strings. Attempting to read or write bare strings in path
fields raises :class:`ValueError`.
"""

from __future__ import annotations

import json
import os
import socket
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, ClassVar, Dict, List, Optional, Union

import logging

from .cluster import detect_current_cluster
from .hashing import normalized_config_hash, sha256_of_file
from .path_resolver import PortablePath, tokenize_argv_paths

_log = logging.getLogger(__name__)

# Cached so repeated sidecar writes in one process don't re-hash the package.
_CODE_FINGERPRINT: Optional[str] = None


def resolve_code_fingerprint() -> str:
    """Identify the *code* that produced an output, as precisely as available.

    Returns, in order of preference:

    1. ``rectify.__git_sha__`` if a build injected one.
    2. ``git:<sha>[+dirty]`` when the installed package sits in a git work tree.
    3. ``pkghash:<12 hex>`` — sha256 over the sorted contents of every ``.py``
       in the package. This is the important fallback: a deployment that is a
       plain *copy* of the repo (no ``.git``) still gets a fingerprint that
       changes whenever the code changes.
    4. ``"unknown"`` only if even that fails.

    Why this exists (2026-08-02): ``__git_sha__`` was read here but **never set
    anywhere in the codebase**, so every provenance sidecar ever written — on
    every machine — recorded ``rectify_git_sha: "unknown"``. The one field meant
    to pin the code version was inert, and that is precisely how a ~2-month-stale
    rectify on Hoffman2 (a non-git copy) produced the MS2 3'-end tables without
    anything flagging the drift. ``pkghash`` closes the non-git case, which is
    the case that actually bit us.
    """
    global _CODE_FINGERPRINT
    if _CODE_FINGERPRINT is not None:
        return _CODE_FINGERPRINT

    fp = "unknown"
    try:
        import rectify as _rect
        injected = getattr(_rect, "__git_sha__", None)
        if injected:
            _CODE_FINGERPRINT = str(injected)
            return _CODE_FINGERPRINT
        pkg_dir = Path(_rect.__file__).resolve().parent

        # (2) git work tree?
        try:
            import subprocess
            # NOTE: `cwd=`, NOT `git -C`. `-C` requires git >= 1.8.5 and
            # Hoffman2 ships git 1.8.3.1, where it silently fails -- which sent
            # a REAL git checkout down the pkghash fallback. Caught only by
            # running this on the actual cluster.
            sha = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=str(pkg_dir), capture_output=True, text=True, timeout=5,
            )
            if sha.returncode == 0 and sha.stdout.strip():
                fp = "git:" + sha.stdout.strip()
                # `+dirty` must mean "the CODE that ran differs from this
                # commit" -- nothing else. Two scopings make that true on real
                # deployments:
                #   --untracked-files=no : deployments keep data/outputs beside
                #       the checkout (the H2 tree is 37 GB); untracked data is
                #       not a code change.
                #   -- .  (cwd is pkg_dir) : restrict to the package subtree.
                #       On H2, 498 tracked docs/scripts/tests differ from HEAD
                #       while rectify/ itself is byte-identical to it -- so an
                #       unscoped check would cry dirty on a run that IS exactly
                #       reproducible from the commit.
                dirty = subprocess.run(
                    ["git", "status", "--porcelain", "--untracked-files=no",
                     "--", "."],
                    cwd=str(pkg_dir), capture_output=True, text=True, timeout=30,
                )
                if dirty.returncode == 0 and dirty.stdout.strip():
                    fp += "+dirty"
                _CODE_FINGERPRINT = fp
                return fp
        except Exception:
            pass

        # (3) content hash of the installed package
        try:
            import hashlib
            h = hashlib.sha256()
            for f in sorted(pkg_dir.rglob("*.py")):
                h.update(str(f.relative_to(pkg_dir)).encode())
                h.update(f.read_bytes())
            fp = "pkghash:" + h.hexdigest()[:12]
            _log.warning(
                "PROVENANCE: rectify at %s is NOT a git checkout — recording %s. "
                "This fingerprint detects code drift but is NOT resolvable back to a "
                "commit, so the run is not reproducible from provenance alone. "
                "Deployments should be git checkouts.", pkg_dir, fp,
            )
        except Exception:
            pass
    except Exception:
        pass

    if fp == "unknown":
        _log.error(
            "PROVENANCE: could not identify the rectify code version at all; "
            "outputs of this run will NOT be traceable to a code state."
        )
    _CODE_FINGERPRINT = fp
    return fp

# ---------------------------------------------------------------------------
# Entry dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class InputEntry:
    """One input file declared for a stage.

    Attributes:
        role: Logical name for this input (e.g. ``"bam"``, ``"genome_fasta"``).
        sha256: ``"sha256:HEX"`` digest of the file at write time.
        size_bytes: File size in bytes (optional; None for virtual inputs).
        path: :class:`PortablePath` envelope; never a bare string.
    """

    role: str
    sha256: str
    size_bytes: Optional[int]
    path: PortablePath

    def to_json(self) -> dict:
        return {
            "role": self.role,
            "sha256": self.sha256,
            "size_bytes": self.size_bytes,
            "path": self.path.to_json(),
        }

    @classmethod
    def from_json(cls, d: dict) -> "InputEntry":
        _require_dict_path_field(d, "path", "InputEntry")
        return cls(
            role=d["role"],
            sha256=d["sha256"],
            size_bytes=d.get("size_bytes"),
            path=PortablePath.from_json(d["path"]),
        )


@dataclass(frozen=True)
class OutputEntry:
    """One primary output file declared for a stage.

    Attributes:
        role: Logical name for this output (e.g. ``"corrected_bam"``).
        sha256: ``"sha256:HEX"`` digest of the file at write time.
        size_bytes: File size in bytes.
        path: :class:`PortablePath` envelope.
        optional: If True, the output may legitimately be absent without
            invalidating the sidecar (e.g. a per-strand split BAM when the
            sample has no reads on that strand).
    """

    role: str
    sha256: str
    size_bytes: Optional[int]
    path: PortablePath
    optional: bool = False

    def to_json(self) -> dict:
        return {
            "role": self.role,
            "sha256": self.sha256,
            "size_bytes": self.size_bytes,
            "path": self.path.to_json(),
            "optional": self.optional,
        }

    @classmethod
    def from_json(cls, d: dict) -> "OutputEntry":
        _require_dict_path_field(d, "path", "OutputEntry")
        return cls(
            role=d["role"],
            sha256=d["sha256"],
            size_bytes=d.get("size_bytes"),
            path=PortablePath.from_json(d["path"]),
            optional=d.get("optional", False),
        )


@dataclass(frozen=True)
class SubOutputEntry:
    """One sub-output (e.g. a per-region BAM in ``$L_SCRATCH``).

    Sub-outputs are typically transient (job-ephemeral) or optional. When
    ``transient=True``, skip-check NEVER sha256-checks them.

    Attributes:
        role: Logical name (e.g. ``"region_bam"``).
        region_id: Optional region identifier (e.g. ``"chrI:1-10000"``).
        sha256: ``"sha256:HEX"`` digest at write time. May be
            ``"sha256:UNKNOWN"`` for transient outputs recorded without
            checking.
        transient: If True, this file lives in ``$L_SCRATCH`` / ``$TMPDIR``
            and is expected to be absent on the next run.
        path: :class:`PortablePath` envelope.
        sidecar: Optional path to a per-region sub-sidecar JSON.
    """

    role: str
    region_id: Optional[str]
    sha256: str
    transient: bool
    path: PortablePath
    sidecar: Optional[PortablePath] = None

    def to_json(self) -> dict:
        d: dict = {
            "role": self.role,
            "region_id": self.region_id,
            "sha256": self.sha256,
            "transient": self.transient,
            "path": self.path.to_json(),
        }
        if self.sidecar is not None:
            d["sidecar"] = self.sidecar.to_json()
        return d

    @classmethod
    def from_json(cls, d: dict) -> "SubOutputEntry":
        _require_dict_path_field(d, "path", "SubOutputEntry")
        sidecar_raw = d.get("sidecar")
        sidecar = None
        if sidecar_raw is not None:
            _require_dict_path_field({"sidecar": sidecar_raw}, "sidecar", "SubOutputEntry.sidecar")
            sidecar = PortablePath.from_json(sidecar_raw)
        return cls(
            role=d["role"],
            region_id=d.get("region_id"),
            sha256=d["sha256"],
            transient=d.get("transient", False),
            path=PortablePath.from_json(d["path"]),
            sidecar=sidecar,
        )


# ---------------------------------------------------------------------------
# Main ProvenanceRecord dataclass
# ---------------------------------------------------------------------------


@dataclass(frozen=False)
class ProvenanceRecord:
    """Full stage-level provenance sidecar (schema_version = ``"1.0"``).

    This dataclass mirrors the on-disk JSON. Use :meth:`from_components` to
    construct from raw stage code, or :meth:`from_json` to deserialize.
    Use :func:`write_stage_sidecar` / :func:`read_stage_sidecar` for I/O.

    Attributes:
        schema_version: Always ``"1.0"``.
        stage: Stage name, e.g. ``"correct"`` or ``"analyze"``.
        stage_subtype: Optional protocol variant (``"drs"``, ``"cdna"``,
            ``"netseq"``, …) or ``None``.
        sample_id: The sample identifier.
        rectify_git_sha: ``git log -1 --format=%H`` at time of writing.
        rectify_version: The ``rectify.__version__`` string.
        started_at: ISO 8601 timestamp of stage start.
        completed_at: ISO 8601 timestamp of stage completion.
        exit_status: Integer exit code. ``0`` = success.
        host: ``socket.gethostname()`` at write time.
        cluster: One of ``"sherlock"``, ``"hoffman2"``, ``"m1"``, ``"other"``.
        scheduler_job_id: SGE ``$JOB_ID`` or SLURM ``$SLURM_JOB_ID``, or
            ``None`` if not running under a scheduler.
        invocation: Dict with keys ``argv``, ``argv_template``,
            ``config_hash``.
        inputs: List of :class:`InputEntry`.
        outputs: List of :class:`OutputEntry`.
        sub_outputs: List of :class:`SubOutputEntry`.
        stats: Free-form dict of stage-specific metrics (e.g.
            ``{"wall_seconds": 42.1, "reads_corrected": 5000}``).
        skip_check_config: Dict with optional keys ``"ignore_input_roles"``
            (list of roles to skip during skip-check) and ``"ignore_argv"``
            (list of flag names to ignore when comparing config_hash).
    """

    SCHEMA_VERSION: ClassVar[str] = "1.0"

    schema_version: str
    stage: str
    stage_subtype: Optional[str]
    sample_id: str
    rectify_git_sha: str
    rectify_version: str
    started_at: str
    completed_at: str
    exit_status: int
    host: str
    cluster: str
    scheduler_job_id: Optional[str]
    invocation: Dict[str, Any]
    inputs: List[InputEntry]
    outputs: List[OutputEntry]
    sub_outputs: List[SubOutputEntry]
    stats: Dict[str, Any]
    skip_check_config: Dict[str, Any]

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_components(
        cls,
        *,
        stage: str,
        sample_id: str,
        sample_output_dir: Union[Path, str],
        started_at: str,
        completed_at: str,
        exit_status: int,
        # Inputs: role -> path (must exist and be hashable)
        inputs: Dict[str, Union[Path, str]],
        # Outputs: role -> (path, optional) tuples, or just role -> path
        outputs: Dict[str, Union[Path, str, tuple]],
        # Sub-outputs: list of dicts with keys:
        #   role, path, transient, [region_id], [sidecar], [sha256]
        sub_outputs: Optional[List[Dict[str, Any]]] = None,
        stats: Optional[Dict[str, Any]] = None,
        skip_check_config: Optional[Dict[str, Any]] = None,
        stage_subtype: Optional[str] = None,
        argv: Optional[List[str]] = None,
        rectify_git_sha: Optional[str] = None,
        rectify_version: Optional[str] = None,
        scheduler_job_id: Optional[str] = None,
    ) -> "ProvenanceRecord":
        """Convenience constructor for stage code.

        Computes sha256 of each input/output, wraps paths in PortablePath,
        tokenizes argv, and builds the full ProvenanceRecord.

        Args:
            stage: Stage name (e.g. ``"correct"``).
            sample_id: The sample identifier.
            sample_output_dir: The stage's output directory (used for
                sample-relative path tokenization).
            started_at: ISO 8601 timestamp; typically captured at stage start.
            completed_at: ISO 8601 timestamp; typically
                ``datetime.now(timezone.utc).isoformat()``.
            exit_status: 0 for success.
            inputs: ``{role: path}`` mapping. Each path must exist.
            outputs: ``{role: path}`` or ``{role: (path, optional_bool)}``.
                Each path must exist.
            sub_outputs: List of dicts describing sub-outputs. Each dict
                must have ``role``, ``path``, ``transient`` keys. Optional
                keys: ``region_id``, ``sidecar``, ``sha256`` (if not
                provided, sha256 is computed if the file exists and
                ``transient=False``; otherwise ``"sha256:UNKNOWN"``).
            stats: Free-form metrics dict.
            skip_check_config: Dict with optional ``"ignore_input_roles"``
                and ``"ignore_argv"`` keys.
            stage_subtype: Protocol variant (``"drs"``, ``"cdna"``, …).
            argv: Raw argv; defaults to ``sys.argv``.
            rectify_git_sha: Current git sha. If ``None``, attempted from
                ``rectify.__git_sha__`` or ``"unknown"``.
            rectify_version: rectify version string. If ``None``, attempted
                from ``rectify.__version__`` or ``"unknown"``.
            scheduler_job_id: SGE/SLURM job id from env (auto-detected if
                ``None``).

        Returns:
            A :class:`ProvenanceRecord` ready for :func:`write_stage_sidecar`.
        """
        sample_output_dir = Path(sample_output_dir)
        argv = argv if argv is not None else sys.argv[:]

        # ---- resolve git sha + version ----
        if rectify_git_sha is None:
            rectify_git_sha = resolve_code_fingerprint()

        if rectify_version is None:
            try:
                import rectify as _rect
                rectify_version = getattr(_rect, "__version__", "unknown")
            except ImportError:
                rectify_version = "unknown"

        # ---- scheduler job id (auto-detect) ----
        if scheduler_job_id is None:
            scheduler_job_id = (
                os.environ.get("SLURM_JOB_ID")
                or os.environ.get("JOB_ID")  # SGE
            )

        # ---- build InputEntry list ----
        input_entries: List[InputEntry] = []
        for role, p in inputs.items():
            p = Path(p)
            hex_digest = sha256_of_file(p)
            sha = "sha256:" + hex_digest
            size = p.stat().st_size
            pp = PortablePath.from_path(p, sample_output_dir=sample_output_dir)
            input_entries.append(InputEntry(
                role=role,
                sha256=sha,
                size_bytes=size,
                path=pp,
            ))

        # ---- build OutputEntry list ----
        output_entries: List[OutputEntry] = []
        for role, val in outputs.items():
            optional_flag = False
            if isinstance(val, tuple):
                p, optional_flag = val[0], val[1]
            else:
                p = val
            p = Path(p)
            hex_digest = sha256_of_file(p)
            sha = "sha256:" + hex_digest
            size = p.stat().st_size
            pp = PortablePath.from_path(p, sample_output_dir=sample_output_dir)
            output_entries.append(OutputEntry(
                role=role,
                sha256=sha,
                size_bytes=size,
                path=pp,
                optional=optional_flag,
            ))

        # ---- build SubOutputEntry list ----
        sub_output_entries: List[SubOutputEntry] = []
        for sub in (sub_outputs or []):
            p = Path(sub["path"])
            transient = sub.get("transient", False)
            # If sha256 provided, use it; otherwise compute if file exists.
            if "sha256" in sub:
                sha = sub["sha256"]
                if not sha.startswith("sha256:"):
                    sha = "sha256:" + sha
            elif p.exists():
                sha = "sha256:" + sha256_of_file(p)
            else:
                sha = "sha256:UNKNOWN"
            pp = PortablePath.from_path(
                p,
                sample_output_dir=sample_output_dir,
                transient=transient,
            )
            sidecar_pp: Optional[PortablePath] = None
            if sub.get("sidecar"):
                sidecar_pp = PortablePath.from_path(
                    sub["sidecar"],
                    sample_output_dir=sample_output_dir,
                    transient=transient,
                )
            sub_output_entries.append(SubOutputEntry(
                role=sub["role"],
                region_id=sub.get("region_id"),
                sha256=sha,
                transient=transient,
                path=pp,
                sidecar=sidecar_pp,
            ))

        # ---- invocation block ----
        argv_template = tokenize_argv_paths(argv, sample_output_dir)
        ignore_argv_set = set(
            (skip_check_config or {}).get("ignore_argv", [])
        )
        config_hash = normalized_config_hash(argv_template, ignore=ignore_argv_set)
        invocation = {
            "argv": argv,
            "argv_template": argv_template,
            "config_hash": config_hash,
        }

        return cls(
            schema_version=cls.SCHEMA_VERSION,
            stage=stage,
            stage_subtype=stage_subtype,
            sample_id=sample_id,
            rectify_git_sha=rectify_git_sha,
            rectify_version=rectify_version,
            started_at=started_at,
            completed_at=completed_at,
            exit_status=exit_status,
            host=socket.gethostname(),
            cluster=detect_current_cluster(),
            scheduler_job_id=scheduler_job_id,
            invocation=invocation,
            inputs=input_entries,
            outputs=output_entries,
            sub_outputs=sub_output_entries,
            stats=stats or {},
            skip_check_config=skip_check_config or {},
        )

    # ------------------------------------------------------------------
    # JSON serialization / deserialization
    # ------------------------------------------------------------------

    def to_json(self) -> dict:
        """Serialize to a JSON-compatible dict."""
        return {
            "schema_version": self.schema_version,
            "stage": self.stage,
            "stage_subtype": self.stage_subtype,
            "sample_id": self.sample_id,
            "rectify_git_sha": self.rectify_git_sha,
            "rectify_version": self.rectify_version,
            "started_at": self.started_at,
            "completed_at": self.completed_at,
            "exit_status": self.exit_status,
            "host": self.host,
            "cluster": self.cluster,
            "scheduler_job_id": self.scheduler_job_id,
            "invocation": self.invocation,
            "inputs": [e.to_json() for e in self.inputs],
            "outputs": [e.to_json() for e in self.outputs],
            "sub_outputs": [e.to_json() for e in self.sub_outputs],
            "stats": self.stats,
            "skip_check_config": self.skip_check_config,
        }

    @classmethod
    def from_json(cls, d: dict) -> "ProvenanceRecord":
        """Deserialize from a dict loaded from a sidecar file.

        Raises:
            ValueError: On schema violations (bare string in path field,
                missing required key, wrong schema_version, etc.).
        """
        _require_key(d, "schema_version", "ProvenanceRecord")
        if d["schema_version"] != cls.SCHEMA_VERSION:
            raise ValueError(
                f"ProvenanceRecord: unsupported schema_version "
                f"{d['schema_version']!r} (expected {cls.SCHEMA_VERSION!r}). "
                "Schema migration is not yet implemented."
            )
        for key in (
            "stage", "sample_id", "rectify_git_sha", "rectify_version",
            "started_at", "completed_at", "exit_status", "host",
            "cluster", "invocation", "inputs", "outputs", "sub_outputs",
            "stats", "skip_check_config",
        ):
            _require_key(d, key, "ProvenanceRecord")

        inputs = [InputEntry.from_json(e) for e in d["inputs"]]
        outputs = [OutputEntry.from_json(e) for e in d["outputs"]]
        sub_outputs = [SubOutputEntry.from_json(e) for e in d["sub_outputs"]]

        return cls(
            schema_version=d["schema_version"],
            stage=d["stage"],
            stage_subtype=d.get("stage_subtype"),
            sample_id=d["sample_id"],
            rectify_git_sha=d["rectify_git_sha"],
            rectify_version=d["rectify_version"],
            started_at=d["started_at"],
            completed_at=d["completed_at"],
            exit_status=d["exit_status"],
            host=d["host"],
            cluster=d["cluster"],
            scheduler_job_id=d.get("scheduler_job_id"),
            invocation=d["invocation"],
            inputs=inputs,
            outputs=outputs,
            sub_outputs=sub_outputs,
            stats=d["stats"],
            skip_check_config=d["skip_check_config"],
        )


# ---------------------------------------------------------------------------
# I/O functions
# ---------------------------------------------------------------------------


def write_stage_sidecar(record: ProvenanceRecord, sample_output: Union[Path, str]) -> Path:
    """Atomically write a stage sidecar JSON to ``sample_output``.

    Filename: ``<sample_output>/<sample_id>.<stage>.provenance.json``.

    Atomic write: write to ``<filename>.tmp`` then :func:`os.replace` (POSIX
    rename, atomic). A crash between the two steps leaves a ``*.tmp`` file but
    never a corrupt final sidecar.

    Args:
        record: The :class:`ProvenanceRecord` to serialize.
        sample_output: The sample's output directory (must already exist).

    Returns:
        The final sidecar path.

    Raises:
        ValueError: If any path field in ``record`` is a bare string rather
            than a :class:`PortablePath` (schema violation).
    """
    sample_output = Path(sample_output)
    _validate_record_paths(record)

    sidecar_name = f"{record.sample_id}.{record.stage}.provenance.json"
    sidecar_path = sample_output / sidecar_name
    tmp_path = sidecar_path.with_suffix(".json.tmp")

    json_dict = record.to_json()
    with open(tmp_path, "w", encoding="utf-8") as fh:
        json.dump(json_dict, fh, indent=2)

    os.replace(tmp_path, sidecar_path)
    return sidecar_path


def read_stage_sidecar(path: Union[Path, str]) -> ProvenanceRecord:
    """Read and validate a stage sidecar JSON.

    Args:
        path: Path to the sidecar JSON file.

    Returns:
        A parsed :class:`ProvenanceRecord`.

    Raises:
        ValueError: On any schema violation (bad version, missing keys,
            bare strings in path fields, malformed JSON).
        FileNotFoundError: If ``path`` doesn't exist.
    """
    path = Path(path)
    try:
        with open(path, "r", encoding="utf-8") as fh:
            raw = fh.read()
    except FileNotFoundError:
        raise

    try:
        d = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"read_stage_sidecar: invalid JSON in {path}: {exc}"
        ) from exc

    if not isinstance(d, dict):
        raise ValueError(
            f"read_stage_sidecar: expected top-level JSON object, "
            f"got {type(d).__name__} in {path}"
        )

    try:
        return ProvenanceRecord.from_json(d)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            f"read_stage_sidecar: schema violation in {path}: {exc}"
        ) from exc


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _require_key(d: dict, key: str, context: str) -> None:
    if key not in d:
        raise ValueError(
            f"{context}: required key {key!r} missing from dict. "
            f"Got keys: {sorted(d.keys())}"
        )


def _require_dict_path_field(d: dict, field_name: str, context: str) -> None:
    """Assert that ``d[field_name]`` is a dict (PortablePath envelope), not a bare string."""
    val = d.get(field_name)
    if val is None:
        raise ValueError(
            f"{context}: required path field {field_name!r} is missing."
        )
    if isinstance(val, str):
        raise ValueError(
            f"{context}: path field {field_name!r} must be a PortablePath "
            f"dict envelope, not a bare string {val!r}. "
            "Schema violation: all path fields must use PortablePath.to_json()."
        )
    if not isinstance(val, dict):
        raise ValueError(
            f"{context}: path field {field_name!r} must be a dict "
            f"(PortablePath envelope), got {type(val).__name__}={val!r}."
        )


def _validate_record_paths(record: ProvenanceRecord) -> None:
    """Raise if any path field in the record is not a PortablePath.

    This is a defence-in-depth check: normal construction via
    :meth:`ProvenanceRecord.from_components` always produces PortablePath
    objects. But callers who manually build the dataclass could slip in a
    bare string.
    """
    for entry in record.inputs:
        if not isinstance(entry.path, PortablePath):
            raise ValueError(
                f"write_stage_sidecar: InputEntry role={entry.role!r} "
                f"has path type {type(entry.path).__name__!r} (expected PortablePath). "
                "Schema violation: wrap paths via PortablePath.from_path()."
            )
    for entry in record.outputs:
        if not isinstance(entry.path, PortablePath):
            raise ValueError(
                f"write_stage_sidecar: OutputEntry role={entry.role!r} "
                f"has path type {type(entry.path).__name__!r} (expected PortablePath). "
                "Schema violation: wrap paths via PortablePath.from_path()."
            )
    for entry in record.sub_outputs:
        if not isinstance(entry.path, PortablePath):
            raise ValueError(
                f"write_stage_sidecar: SubOutputEntry role={entry.role!r} "
                f"has path type {type(entry.path).__name__!r} (expected PortablePath). "
                "Schema violation: wrap paths via PortablePath.from_path()."
            )
        if entry.sidecar is not None and not isinstance(entry.sidecar, PortablePath):
            raise ValueError(
                f"write_stage_sidecar: SubOutputEntry role={entry.role!r} "
                f"sidecar has type {type(entry.sidecar).__name__!r} (expected PortablePath). "
                "Schema violation: wrap paths via PortablePath.from_path()."
            )
