#!/usr/bin/env python3
"""
Per-BAM provenance sidecar — strict reuse gate for ``rectify run-all``.

When ``rectify run-all`` (or ``rectify align``) writes a BAM it expects callers
to be able to reuse later (per-aligner BAMs, ``{sample}.multialigned.bam``), it
also drops a sidecar JSON: ``{bam}.prov.json``. The sidecar records the
rectify version, the git SHA, the aligner name + version, and the timestamp.

On a later invocation, the reuse gate reads the sidecar and refuses to reuse
the BAM unless the recorded ``(rectify_git_sha, aligner_name, aligner_version)``
triple matches the current run's expected values exactly.

The match policy is intentionally strict — strict equality on all three
fields, no semver tolerance. Override with the ``--trust-existing-bams`` flag
when an explicit human-mediated bypass is needed (e.g. reusing BAMs produced
before this feature shipped).

Public surface — keep it small:

- :func:`compute_run_provenance` — capture once at run start; pass into every
  write_sidecar call for that run.
- :func:`write_sidecar` — drops ``{bam}.prov.json`` next to ``bam``.
- :func:`read_sidecar` — returns dict or ``None``.
- :func:`matches_strict` — ``(ok, reason)`` comparing stored vs current.
- :func:`get_rectify_git_sha` — best-effort HEAD SHA of the rectify install.
- :func:`get_aligner_version` — runs ``<aligner> --version`` and parses.
"""

from __future__ import annotations

import json
import logging
import os
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .. import __version__ as _RECTIFY_PKG_VERSION

logger = logging.getLogger(__name__)

SIDECAR_SUFFIX = ".prov.json"
SIDECAR_SCHEMA_VERSION = 1


_GIT_SHA_CACHE: Optional[str] = None
_ALIGNER_VERSION_CACHE: Dict[str, Optional[str]] = {}


def _rectify_repo_root() -> Optional[Path]:
    """Return the directory containing the rectify package source, if it
    lives inside a git checkout. ``None`` if the install is from a wheel/
    sdist (no .git present in any parent of ``rectify/``)."""
    here = Path(__file__).resolve().parent.parent  # rectify/utils/.. = rectify/
    for ancestor in [here, *here.parents]:
        if (ancestor / ".git").exists():
            return ancestor
    return None


def get_rectify_git_sha() -> Optional[str]:
    """Best-effort: full SHA of ``HEAD`` of the rectify install's git
    checkout. Returns ``None`` if rectify was installed from a wheel/sdist
    (no enclosing .git directory)."""
    global _GIT_SHA_CACHE
    if _GIT_SHA_CACHE is not None:
        return _GIT_SHA_CACHE or None
    root = _rectify_repo_root()
    if root is None:
        _GIT_SHA_CACHE = ""
        return None
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(root), capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0:
            sha = result.stdout.strip()
            _GIT_SHA_CACHE = sha
            return sha
    except Exception as exc:
        logger.debug("get_rectify_git_sha failed: %s", exc)
    _GIT_SHA_CACHE = ""
    return None


# Per-aligner version-probe commands. Each entry is ``(argv, regex)``: we run
# ``argv`` and search the combined stdout+stderr for ``regex`` (case-sensitive,
# multiline). The first capture group becomes the version string.
#
# Many of these tools print their version to stderr, not stdout — we always
# combine both streams before searching.
_ALIGNER_VERSION_PROBES: Dict[str, Tuple[List[str], str]] = {
    "minimap2":  (["minimap2", "--version"],          r"^([0-9][\w.\-]+)"),
    "mapPacBio": (["bbversion.sh"],                   r"BBMap version ([\d.]+)"),
    "gapmm2":    (["gapmm2", "--version"],            r"gapmm2[, ]+version[: ]+([\w.]+)"),
    "uLTRA":     (["uLTRA", "--version"],             r"uLTRA[ ]+([\w.]+)"),
    "deSALT":    (["deSALT", "--version"],            r"^([\w.]+)"),
    "bbmap":     (["bbversion.sh"],                   r"BBMap version ([\d.]+)"),
    "bwa":       (["bwa"],                            r"Version:\s*([\w.\-]+)"),
}


def get_aligner_version(aligner_name: str) -> Optional[str]:
    """Best-effort: version string of ``aligner_name``. Returns ``None`` if
    the aligner is not on PATH or its version probe fails. Cached per name.

    Special case: ``'consensus'`` is the rectify consensus-selection step,
    not a third-party binary. Its version is the rectify package version."""
    if aligner_name in _ALIGNER_VERSION_CACHE:
        return _ALIGNER_VERSION_CACHE[aligner_name]
    if aligner_name == "consensus":
        _ALIGNER_VERSION_CACHE["consensus"] = _RECTIFY_PKG_VERSION
        return _RECTIFY_PKG_VERSION
    probe = _ALIGNER_VERSION_PROBES.get(aligner_name)
    if probe is None:
        _ALIGNER_VERSION_CACHE[aligner_name] = None
        return None
    argv, pattern = probe
    if shutil.which(argv[0]) is None:
        _ALIGNER_VERSION_CACHE[aligner_name] = None
        return None
    try:
        result = subprocess.run(
            argv, capture_output=True, text=True, timeout=15,
        )
    except Exception as exc:
        logger.debug("aligner version probe failed for %s: %s", aligner_name, exc)
        _ALIGNER_VERSION_CACHE[aligner_name] = None
        return None
    combined = (result.stdout or "") + "\n" + (result.stderr or "")
    m = re.search(pattern, combined, flags=re.MULTILINE)
    version = m.group(1) if m else None
    _ALIGNER_VERSION_CACHE[aligner_name] = version
    return version


def compute_run_provenance(
    command: Optional[List[str]] = None,
    extra: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    """Capture once at run start. Returns a dict that should be passed into
    every :func:`write_sidecar` call for the duration of the run.

    The dict has all fields *except* ``aligner_name`` / ``aligner_version``,
    which are added by :func:`write_sidecar` at per-aligner write time.

    Args:
        command: ``sys.argv``-style invocation. Stored for debugging only;
            not part of the match key.
        extra: optional extra fields to merge in. Caller-controlled.
    """
    payload: Dict[str, object] = {
        "schema_version": SIDECAR_SCHEMA_VERSION,
        "rectify_pkg_version": _RECTIFY_PKG_VERSION,
        "rectify_git_sha": get_rectify_git_sha(),
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "command": list(command) if command is not None else None,
        "hostname": _hostname(),
    }
    if extra:
        payload.update(extra)
    return payload


def _hostname() -> Optional[str]:
    try:
        import socket
        return socket.gethostname()
    except Exception:
        return None


def sidecar_path(bam_path: Path) -> Path:
    """Conventional sidecar path for ``bam_path``."""
    return Path(str(bam_path) + SIDECAR_SUFFIX)


def write_sidecar(
    bam_path: Path,
    run_provenance: Dict[str, object],
    aligner_name: str,
    aligner_version: Optional[str] = None,
    extra: Optional[Dict[str, object]] = None,
) -> Path:
    """Write ``{bam_path}.prov.json`` next to ``bam_path``.

    ``run_provenance`` is the dict returned by :func:`compute_run_provenance`
    (captured once at run start; shared across sidecars in the same run).
    ``aligner_name`` is required. ``aligner_version`` is auto-probed via
    :func:`get_aligner_version` when omitted.

    Returns the sidecar path. Best-effort: a failure to write logs a warning
    and re-raises only ``OSError`` (so the caller can decide whether a
    missing sidecar should abort the BAM write).
    """
    if aligner_version is None:
        aligner_version = get_aligner_version(aligner_name)
    sidecar = dict(run_provenance)
    sidecar["aligner_name"] = aligner_name
    sidecar["aligner_version"] = aligner_version
    if extra:
        sidecar.update(extra)
    target = sidecar_path(bam_path)
    # Atomic write: tmp + rename.
    tmp = target.with_suffix(target.suffix + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(sidecar, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, target)
    return target


def read_sidecar(bam_path: Path) -> Optional[Dict[str, object]]:
    """Return parsed sidecar dict, or ``None`` if absent / unreadable."""
    target = sidecar_path(bam_path)
    if not target.exists():
        return None
    try:
        with open(target) as fh:
            return json.load(fh)
    except Exception as exc:
        logger.warning("Failed to parse sidecar %s: %s", target, exc)
        return None


_MATCH_KEY_FIELDS: Tuple[str, ...] = (
    "rectify_git_sha",
    "aligner_name",
    "aligner_version",
)


def matches_strict(
    stored: Optional[Dict[str, object]],
    current: Dict[str, object],
) -> Tuple[bool, str]:
    """Return ``(ok, reason)``.

    Match policy v1: strict equality on every field in
    :data:`_MATCH_KEY_FIELDS`. Missing sidecar (``stored is None``) is a
    mismatch with a distinct reason so callers can render a helpful
    "missing sidecar" message.
    """
    if stored is None:
        return False, "no sidecar present"
    differences: List[str] = []
    for field in _MATCH_KEY_FIELDS:
        s_val = stored.get(field)
        c_val = current.get(field)
        if s_val != c_val:
            differences.append(f"{field}: stored={s_val!r} current={c_val!r}")
    if differences:
        return False, "; ".join(differences)
    return True, "match"


def expected_provenance_for_aligner(
    run_provenance: Dict[str, object],
    aligner_name: str,
    aligner_version: Optional[str] = None,
) -> Dict[str, object]:
    """Build the expected-provenance dict that :func:`matches_strict` will
    compare against, using the same fields :func:`write_sidecar` writes."""
    if aligner_version is None:
        aligner_version = get_aligner_version(aligner_name)
    expected = dict(run_provenance)
    expected["aligner_name"] = aligner_name
    expected["aligner_version"] = aligner_version
    return expected


__all__ = [
    "SIDECAR_SUFFIX",
    "SIDECAR_SCHEMA_VERSION",
    "compute_run_provenance",
    "expected_provenance_for_aligner",
    "get_aligner_version",
    "get_rectify_git_sha",
    "matches_strict",
    "read_sidecar",
    "sidecar_path",
    "write_sidecar",
]
