"""Universal per-command provenance: EVERY rectify invocation records its version.

Why this exists
---------------
An audit on 2026-08-02 found that only 5 of 26 command modules wrote a provenance
sidecar, and that ``rectify_git_sha`` resolved to ``"unknown"`` in all 33 sidecars
behind the MS2 3'-end tables. Net effect: for most rectify operations there was no
record at all of which code produced the output, and for the rest the record was
inert. A ~2-month-stale rectify on Hoffman2 produced production data with nothing
flagging it.

Rather than hand-wire 21 more modules (21 chances to get it subtly wrong, and no
guarantee for the 22nd command someone adds next month), this hooks the single
CLI dispatch point. Every subcommand therefore leaves a record, whether or not
its module knows anything about provenance.

Relationship to the rich per-stage sidecars
-------------------------------------------
Stages that write their OWN sidecar (correct, analyze, drs_trim, restore_polya,
trim_cdna_polya, align, cdna_correct, run/stages) keep doing so — those carry
input/output sha256s, stats and skip-check config, which this cannot. This is the
FLOOR, not a replacement: it guarantees the version question is always answerable.
Its filename is distinct (``<command>.command.provenance.json``) so it never
collides with or overwrites a stage sidecar.

Commands that produce no data artifact
--------------------------------------
``test``, ``install-aligners`` and ``verify`` have no output location and produce
nothing to attribute, so they are skipped by ``SKIP_COMMANDS``. Everything else
writes, falling back to the current working directory when no output argument can
be resolved — a record in the wrong-ish place beats no record.
"""
from __future__ import annotations

import json
import logging
import os
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

_log = logging.getLogger(__name__)

#: Commands that transform nothing and have nowhere meaningful to write.
SKIP_COMMANDS = {"test", "install-aligners", "verify"}

#: argparse destinations that plausibly name an output location, best first.
_OUTPUT_ARG_NAMES = (
    "output_dir", "outdir", "out_dir",
    "output", "out",
    "output_path", "output_prefix", "output_bam", "output_tsv",
    "o",
)


def resolve_output_dir(args: Any) -> Optional[Path]:
    """Best-effort output directory for *args*, or None.

    Takes the first ``_OUTPUT_ARG_NAMES`` attribute that is set; if it looks
    like a file (has a suffix) its parent is used. Never raises.
    """
    for name in _OUTPUT_ARG_NAMES:
        val = getattr(args, name, None)
        if not val or not isinstance(val, (str, Path)):
            continue
        try:
            p = Path(val)
        except Exception:
            continue
        # A path with a suffix is a file; its directory is the output location.
        return p.parent if p.suffix else p
    return None


def write_command_sidecar(
    command: str,
    args: Any,
    started_at: str,
    wall_seconds: float,
    exit_status: int,
    argv: Optional[list] = None,
) -> Optional[Path]:
    """Write ``<outdir>/<command>.command.provenance.json``. Never raises.

    Returns the path written, or None if skipped/failed.
    """
    if command in SKIP_COMMANDS:
        return None
    try:
        from rectify.utils.version import get_rectify_git_sha

        out_dir = resolve_output_dir(args) or Path.cwd()
        try:
            out_dir.mkdir(parents=True, exist_ok=True)
        except Exception:
            out_dir = Path.cwd()

        try:
            import rectify as _rect
            version = getattr(_rect, "__version__", "unknown")
        except Exception:
            version = "unknown"

        rec = {
            "schema_version": 1,
            "record_type": "command_invocation",
            "stage": command,
            "rectify_git_sha": get_rectify_git_sha(),
            "rectify_version": version,
            "started_at": started_at,
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "wall_seconds": round(wall_seconds, 3),
            "exit_status": exit_status,
            "host": socket.gethostname(),
            "scheduler_job_id": (
                os.environ.get("SLURM_JOB_ID")
                or os.environ.get("JOB_ID")
                or os.environ.get("PBS_JOBID")
            ),
            "sge_task_id": os.environ.get("SGE_TASK_ID"),
            "cwd": str(Path.cwd()),
            "argv": list(argv if argv is not None else sys.argv),
            # Full resolved argparse namespace — the knobs that shaped the run.
            "args": {
                k: (str(v) if isinstance(v, Path) else v)
                for k, v in sorted(vars(args).items())
                if isinstance(v, (str, int, float, bool, type(None), Path))
            },
        }
        path = out_dir / f"{command}.command.provenance.json"
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(rec, indent=2, default=str))
        tmp.replace(path)          # atomic: a reader never sees a partial file
        return path
    except Exception as exc:            # never let provenance break a command
        _log.debug("command sidecar not written for %s: %s", command, exc)
        return None
