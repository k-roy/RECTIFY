# Sonnet Briefing — Commit A.5 (Provenance + resume infrastructure)

**You are a Sonnet 4.6 subagent.** Goal: land Commit A.5 on M1 working tree. Opus owns the high-judgment pieces (`path_resolver.py`, `skip_check.py`, the migration decision); **you own the boilerplate** (`sidecar.py`, `hashing.py`, `cluster.py`, `cli.py`, `__init__.py`, and all the tests).

**Hard prerequisite:** Commit A must be landed + merged on `drs-validation-rebuild` before you start. Verify `git log --oneline -3` shows Commit A's commit at HEAD before doing anything else.

---

## 1. Pre-recorded Opus decisions (do NOT re-litigate)

### 1.1 Migration decision: COEXIST + cleanup

`rectify/utils/provenance.py` already exists with a `ProvenanceTracker` accumulator class. Opus inspected it on 2026-05-19. Decision recorded in this briefing:

- **Keep `rectify/utils/provenance.py` UNTOUCHED.** It's a user-facing per-directory run log + README writer. Lab processed-output directories already depend on its output format. Do not modify, do not move, do not rename.
- **Create new package `rectify/core/provenance/`** for stage-level sidecars + skip-check (this commit's scope).
- **Delete `ProvenanceTracker.register_staged()` + its 2 callsites in `rectify/core/commands/run/single_sample.py:507, 529`.** That method stores a scratch→Oak path mapping in an in-memory `_staged_map` dict that is NEVER persisted to disk and NEVER read back. It's dead code. Removing the callsites + the method definition is safe.

Naming clarification: yes, there will be both a `rectify.utils.provenance` module (legacy accumulator) and a `rectify.core.provenance` package (new stage-level). Import paths are distinct. Document this in `rectify/core/provenance/__init__.py` with a paragraph explaining the split.

### 1.2 Schema decisions (already in spec §6.5.2)

- Stage-level sidecar lives at `<sample_output>/<sample_id>.<stage>.provenance.json`.
- Schema version `1.0`. Field policy in spec §6.5.2.
- All path fields are `PortablePath` envelopes (§6.5.2.5), not bare strings. Sidecar I/O code REFUSES to read or write bare-string paths in path fields — schema violation → raise.
- SHA-256, not MD5. (Existing `ProvenanceTracker` uses MD5; the new package does not.)
- Sub-outputs flagged `transient: true` are NOT sha256-checked at skip-check time.
- Persistent outputs (`transient: false` or unset) MUST be sample-relative (`kind="sample_relative"`) or env-relative against env vars OTHER than `L_SCRATCH`/`TMPDIR`. `PortablePath` constructor REFUSES to build a non-transient `env_relative` with `env_var ∈ {L_SCRATCH, TMPDIR}`. Sidecar writer respects this; trying to violate it raises `ValueError`.

### 1.3 File ownership

| File | Owner | Approx LOC |
|---|---|---|
| `rectify/core/provenance/__init__.py` | Sonnet | ~40 |
| `rectify/core/provenance/cluster.py` | Sonnet | ~50 |
| `rectify/core/provenance/hashing.py` | Sonnet | ~80 |
| `rectify/core/provenance/cli.py` | Sonnet | ~60 |
| `rectify/core/provenance/sidecar.py` | Sonnet | ~250 |
| `rectify/core/provenance/path_resolver.py` | **Opus (pre-written)** | ~250 |
| `rectify/core/provenance/skip_check.py` | **Opus (pre-written)** | ~250 |
| `rectify/core/commands/run/single_sample.py` (delete 2 lines + check imports) | Sonnet | ~10 |
| `rectify/utils/provenance.py` (delete `register_staged` method) | Sonnet | ~15 |
| `tests/test_provenance_sidecar.py` | Sonnet | ~150 |
| `tests/test_provenance_skip_check.py` | Sonnet | ~100 |
| `tests/test_portable_path.py` | Sonnet | ~200 |
| `tests/test_cluster_detect.py` | Sonnet | ~40 |
| `tests/test_hashing.py` | Sonnet | ~60 |

Opus writes `path_resolver.py` and `skip_check.py` FIRST (before you start), commits them via the same Commit A.5 commit, and includes them in the diff under review. You write everything else.

---

## 2. Critical rules

1. **No commits, no pushes.** Working tree only. Opus reviews + commits at the end.
2. **Verify Commit A landed.** First action: `cd /Users/kevinroy/work/rectify && git log --oneline -3`. If Commit A's commit message ("...region-parallel BAM write infrastructure..." or similar) is not in HEAD or HEAD~1, STOP and report.
3. **`path_resolver.py` and `skip_check.py` are Opus-pre-written.** Read them but DO NOT modify. Your modules import from them and depend on their APIs being as written.
4. **Verify, don't recite.** Per CLAUDE.md.
5. **No `git add -A`.** Surgical staging only (and you don't commit anyway).
6. **`pytest -m "not slow"` MUST pass.** New tests must also pass.
7. **No wiring through to `correct_command.py` or any stage yet.** This commit ADDS infrastructure. Commit B does the first wiring.

---

## 3. Module specs (Sonnet-owned files)

### 3.1 `rectify/core/provenance/__init__.py`

Re-exports the public API + a docstring explaining the split with `rectify.utils.provenance`. Example:

```python
"""Stage-level provenance sidecars + resume infrastructure.

This package is the MACHINE-READABLE skip-cache for `rectify run-all`.
Each stage emits a sidecar JSON to <sample_output>/<sample_id>.<stage>.provenance.json
on successful completion. Subsequent run-all invocations consult sidecars
via `can_skip_stage()` to decide what to recompute.

Distinct from `rectify.utils.provenance.ProvenanceTracker`, which is a
USER-FACING per-output-directory run log + README writer. The two coexist;
the old accumulator is unchanged.

Public API:
    ProvenanceRecord       — dataclass mirroring the v1.0 sidecar schema
    write_stage_sidecar    — atomic write of a sidecar
    read_stage_sidecar     — load + validate a sidecar
    can_skip_stage         — the skip-check decision tree (see spec §6.5.3)
    SkipDecision           — return type of can_skip_stage
    PortablePath           — path envelope (see spec §6.5.2.5)
    PathResolutionError    — raised by PortablePath.resolve_now
    detect_current_cluster — returns "sherlock"/"hoffman2"/"m1"/"other"
    add_resume_args        — argparse helper for --force-all etc.
"""
from .sidecar import ProvenanceRecord, write_stage_sidecar, read_stage_sidecar
from .skip_check import can_skip_stage, SkipDecision
from .path_resolver import PortablePath, PathResolutionError, tokenize_argv_paths
from .cluster import detect_current_cluster
from .cli import add_resume_args
from .hashing import sha256_of_file, normalized_config_hash

__all__ = [
    "ProvenanceRecord", "write_stage_sidecar", "read_stage_sidecar",
    "can_skip_stage", "SkipDecision",
    "PortablePath", "PathResolutionError", "tokenize_argv_paths",
    "detect_current_cluster",
    "add_resume_args",
    "sha256_of_file", "normalized_config_hash",
]
```

### 3.2 `rectify/core/provenance/cluster.py`

```python
"""Cluster detection. Used by sidecar.write to record `cluster`; used by
skip_check to warn on cross-cluster mismatch."""

import os
import platform
import socket
from typing import Literal

Cluster = Literal["sherlock", "hoffman2", "m1", "other"]


def detect_current_cluster() -> Cluster:
    """Heuristic detection of which cluster we're running on.

    Detection order matters: check more-specific markers first.
    """
    # Sherlock: $SHERLOCK env var, or hostname matches sh<N>-<N>[<NN>]
    if os.environ.get("SHERLOCK"):
        return "sherlock"
    hostname = socket.gethostname().lower()
    if hostname.startswith("sh") and "stanford" in hostname:
        return "sherlock"
    # Some Sherlock compute nodes have hostnames like 'sh03-07n10' without
    # the .stanford suffix in $HOSTNAME — pattern-match.
    import re
    if re.match(r"^sh\d+-\d+", hostname):
        return "sherlock"

    # Hoffman2: hostname has .hoffman2.idre.ucla.edu suffix, or matches
    # n<NNNN> / login<N> compute/login nodes.
    if "hoffman2" in hostname or hostname.endswith(".idre.ucla.edu"):
        return "hoffman2"
    if re.match(r"^n\d{3,5}$", hostname) or re.match(r"^login\d+$", hostname):
        # These could be H2 OR another cluster; require an additional marker.
        # H2 has /u/project — check for it.
        if os.path.isdir("/u/project"):
            return "hoffman2"

    # M1 laptop: Darwin + arm64
    uname = platform.uname()
    if uname.system == "Darwin" and uname.machine == "arm64":
        return "m1"

    return "other"
```

### 3.3 `rectify/core/provenance/hashing.py`

```python
"""SHA-256 of files + argv config hashing.

File sha256 uses mmap for files > 16 MB (fast path). Smaller files use
chunked read. Returns hex digest.

normalized_config_hash takes an argv list and an ignore-list of flag names;
filters out the ignored flags + their values; sorts the result; hashes it.
The result is stable across runs that differ only in path-related or
threading flags.
"""
import hashlib
import mmap
import os
from pathlib import Path
from typing import Iterable, List, Optional, Set


def sha256_of_file(path: Path | str, *, chunk_size: int = 1 << 20) -> str:
    """Compute SHA-256 hex digest of a file.

    For files > 16 MB, mmap-based read; otherwise chunked.
    Raises FileNotFoundError if path doesn't exist.
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


def normalized_config_hash(argv: List[str], *, ignore: Optional[Set[str]] = None) -> str:
    """Hash an argv list after removing flags listed in `ignore` (with their values).

    Args:
        argv: full argv (typically `sys.argv` or `prior['invocation']['argv_template']`).
        ignore: set of flag NAMES to filter out (e.g. {"--n-threads", "--tmp-dir"}).
                For each flag, the following positional value (if any) is also dropped.

    Returns SHA-256 hex digest of the canonical-form normalized argv.

    Canonicalization: filtered argv items are joined with U+001F (unit separator)
    so two distinct argvs with different splits can't accidentally collide.
    """
    ignore = ignore or set()
    filtered: List[str] = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        # `--foo=bar` form (one token)
        if "=" in arg and arg.startswith("--"):
            flag = arg.split("=", 1)[0]
            if flag in ignore:
                i += 1
                continue
            filtered.append(arg)
            i += 1
            continue
        # `--foo bar` form (two tokens)
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
```

### 3.4 `rectify/core/provenance/cli.py`

```python
"""Shared argparse helpers for resume flags. Used by every stage parser
that wants to honor --force-*/--accept-prior-provenance/--dry-run-resume."""

import argparse


def add_resume_args(parser: argparse.ArgumentParser) -> None:
    """Add the resume/sidecar flags to `parser`.

    Flags added:
      --force-all                    : ignore all sidecars; rerun everything
      --force-stage NAME[,NAME...]   : ignore sidecars for these stages only
      --accept-prior-provenance      : treat git_sha mismatch as non-blocking
      --dry-run-resume               : print SKIP/RUN decisions and exit
    """
    g = parser.add_argument_group("Resume / sidecar")
    g.add_argument(
        "--force-all",
        action="store_true",
        default=False,
        help="Ignore all sidecars; rerun every stage.",
    )
    g.add_argument(
        "--force-stage",
        type=str,
        default=None,
        metavar="NAME[,NAME...]",
        help="Comma-separated list of stage names to force-rerun (e.g. analyze,correct). "
             "Stages downstream of a forced stage are also forced because their input "
             "sha256 will change.",
    )
    g.add_argument(
        "--accept-prior-provenance",
        action="store_true",
        default=False,
        help="Treat rectify git_sha mismatch between prior and current run as non-blocking. "
             "Default behavior is to rerun any stage whose prior run was at a different "
             "git_sha (safer; conservative). Use this for cosmetic commits that don't "
             "affect output bytes (docs, comments).",
    )
    g.add_argument(
        "--dry-run-resume",
        action="store_true",
        default=False,
        help="Print the SKIP/RUN decision for each stage with full diff (which input "
             "or argv changed) and exit. No work done.",
    )
```

### 3.5 `rectify/core/provenance/sidecar.py` (~250 LOC, the big one)

Defines `ProvenanceRecord` dataclass + atomic `write_stage_sidecar` + `read_stage_sidecar`.

Schema follows spec §6.5.2 EXACTLY. Field-by-field. The dataclass should be JSON-roundtrip clean.

Key APIs:

```python
@dataclass(frozen=False)
class ProvenanceRecord:
    schema_version: str           # "1.0"
    stage: str                    # "correct", "analyze", "drs_trim", ...
    stage_subtype: Optional[str]  # "drs", "cdna", "netseq", or None
    sample_id: str
    rectify_git_sha: str
    rectify_version: str
    started_at: str               # ISO 8601
    completed_at: str
    exit_status: int
    host: str                     # socket.gethostname()
    cluster: str                  # "sherlock"/"hoffman2"/"m1"/"other"
    scheduler_job_id: Optional[str]  # $JOB_ID for SGE, $SLURM_JOB_ID for SLURM
    invocation: Dict              # see spec §6.5.2
    inputs: List[InputEntry]
    outputs: List[OutputEntry]
    sub_outputs: List[SubOutputEntry]
    stats: Dict
    skip_check_config: Dict       # {"ignore_input_roles": [...], "ignore_argv": [...]}

@dataclass(frozen=True)
class InputEntry:
    role: str
    sha256: str          # "sha256:HEX..." OR "abc123..." (sha256_of_file returns bare hex; sidecar adds "sha256:" prefix on write)
    size_bytes: Optional[int]
    path: PortablePath   # NOT a bare string

@dataclass(frozen=True)
class OutputEntry:
    role: str
    sha256: str
    size_bytes: Optional[int]
    path: PortablePath
    optional: bool = False

@dataclass(frozen=True)
class SubOutputEntry:
    role: str
    region_id: Optional[str]
    sha256: str
    transient: bool
    path: PortablePath
    sidecar: Optional[PortablePath]  # path to per-region sidecar if any
```

Functions:

- `write_stage_sidecar(record: ProvenanceRecord, sample_output: Path) -> Path`:
  - Computes filename: `<sample_output>/<sample_id>.<stage>.provenance.json`.
  - Atomic write: `.tmp` + rename, mirroring `ProvenanceTracker.save()` at `rectify/utils/provenance.py:209-214`.
  - Returns the final path.
  - JSON serialization: PortablePath has a `to_json()` method (see Opus-pre-written `path_resolver.py`); use it. Refuses bare strings in `path` fields.

- `read_stage_sidecar(path: Path) -> ProvenanceRecord`:
  - Reads the JSON, validates `schema_version == "1.0"` (raise on mismatch — schema bumps will be handled in a future commit), parses all PortablePath fields via `PortablePath.from_json()`.
  - Returns the dataclass.
  - Raises `ValueError` with a precise message on schema violation (bare string in path field, missing required key, etc.).

- `ProvenanceRecord.from_components(...)` (classmethod): convenience constructor for stage code to use at write time. Accepts the "inputs as dict of role->Path" form that stage code naturally has, and converts to `InputEntry` list (computing sha256 + wrapping in PortablePath).

Be careful: the classmethod constructor needs to know how to:
1. Compute sha256 of each input/output file (call `sha256_of_file`).
2. Wrap each path in PortablePath (call `PortablePath.from_path(p, sample_output_dir=...)`).
3. Tokenize argv via `tokenize_argv_paths` (from path_resolver).
4. Compute `config_hash` via `normalized_config_hash`.
5. Detect cluster via `detect_current_cluster()`.

### 3.6 `rectify/utils/provenance.py` — delete `register_staged`

Two changes:
- Delete the `register_staged` method (lines 186-199).
- Delete the `_staged_map` attribute references if any remain (do a final grep).
- Add a comment at the top of `ProvenanceTracker` referencing the new package: `"""... For machine-readable stage-level sidecars + resume, see `rectify.core.provenance`. ProvenanceTracker is the legacy per-directory run log + README writer."""`.

### 3.7 `rectify/core/commands/run/single_sample.py` — delete callsites

Lines 507 and 529: `tracker.register_staged(...)`. Just delete these lines. Verify the `tracker` variable is still used elsewhere in the function (it should be — single_sample.py uses ProvenanceTracker for the existing per-directory run log). If the deletion makes `tracker` unused in a local scope, remove the unused variable too.

### 3.8 Tests

Five new test files. Each PASSES on M1 before reporting done.

**`tests/test_cluster_detect.py`:**
- Detect on M1: should return "m1" (this test only runs locally on Kevin's M1; mark with `@pytest.mark.skipif(platform.system() != "Darwin")`).
- Mock os.environ["SHERLOCK"] = "1" → returns "sherlock".
- Mock socket.gethostname to return "sh03-07n10" → returns "sherlock".
- Mock os.path.isdir to return True for "/u/project" AND hostname "n1821" → returns "hoffman2".
- Default → "other".

**`tests/test_hashing.py`:**
- sha256_of_file on a tiny known file: assert hash matches `hashlib.sha256(content).hexdigest()`.
- sha256_of_file on a 20MB file (mmap path): create a 20MB file with known content, hash matches.
- sha256_of_file on missing file: raises FileNotFoundError.
- normalized_config_hash: same argv → same hash. Adding an ignored flag → same hash. Adding a non-ignored flag → different hash. Reordering `--foo --bar` vs `--bar --foo` → DIFFERENT hash (order matters; the canonicalization is deliberate).

**`tests/test_provenance_sidecar.py`:**
- Round-trip: build a ProvenanceRecord, write_stage_sidecar to a tmpdir, read it back, assert deep-equal.
- Atomic write: simulate a crash during write (write to tmp file, kill before rename); verify final sidecar file isn't corrupt OR doesn't exist.
- Schema violation: try to write a record with a bare-string path field → ValueError.
- read on invalid JSON: ValueError with helpful message.
- read on schema_version != "1.0": ValueError.
- Optional outputs: an output marked `optional=True` is allowed to be absent without invalidating the sidecar.

**`tests/test_provenance_skip_check.py`:**
- Tests `can_skip_stage` for each of the 4 RUN reasons + the SKIP case. Each is a row in a parametrized table.
- Use a tmpdir as sample_output, write a known sidecar there, then test various mutations.
- Don't test cross-cluster resolution here (that's in test_portable_path).

**`tests/test_portable_path.py`:**
- Per spec §6.5.7, six scenarios:
  - (a) sample-relative tokenization with symlinks (use `tmp_path.symlink_to` or `os.symlink`).
  - (b) env_relative tokenization ordering ($L_SCRATCH before $SCRATCH before $HOME). Set all three env vars to different tmpdirs; pass a path that's under multiple. Verify the LONGEST-prefix wins (= $L_SCRATCH, not $HOME).
  - (c) cached_absolute fallback when env var unset on read side. Write sidecar with env_relative path; unset the env var; resolve → falls back to cached_absolute + emits a logging.warning.
  - (d) cross-cluster resolution (write with $OAK=/tmp/fake_oak; unset $OAK; resolve falls back to cached_absolute and warns).
  - (e) ValueError on attempt to construct non-transient `env_relative` PortablePath with env_var in {L_SCRATCH, TMPDIR}.
  - (f) argv_template tokenization round-trip preserves config_hash across runs with different $L_SCRATCH values (this exercises `tokenize_argv_paths`).

---

## 4. Acceptance criteria

- [ ] Commit A's commit message is in `git log --oneline -3` HEAD or HEAD~1.
- [ ] All Sonnet-owned files (3.1-3.4, 3.5) exist with the spec'd APIs.
- [ ] `register_staged` is deleted from `rectify/utils/provenance.py`; callsites in `single_sample.py:507, 529` are also deleted.
- [ ] `pytest -m "not slow"` passes (no regressions to ~934 baseline + your new tests).
- [ ] All new tests pass.
- [ ] No `correct_command.py` / `analyze_command.py` / stage code has been wired through to the new infra. (Commit B does the first wiring.)
- [ ] `rectify.core.provenance.detect_current_cluster()` correctly returns "m1" on Kevin's M1 — verify by running it once before reporting done.
- [ ] No `git add` / `git commit` runs.
- [ ] `git status -s` and `git diff --stat` are captured in your handoff.

---

## 5. Stop-and-ask triggers

- Commit A is not in git log — Commit A hasn't been merged yet. Stop.
- `path_resolver.py` or `skip_check.py` doesn't exist yet — Opus hasn't written them yet. Stop and report.
- Existing tests regress — investigate root cause, do not just skip the test.
- The existing `ProvenanceTracker` has callers besides `single_sample.py:507, 529` that pass `register_staged` — grep + report; don't silently break them.
- Any import-order issue with the new `__init__.py` re-exports (circular imports, etc.) — report.

---

## 6. Reporting back

When done (all acceptance boxes checked):
1. `git status -s` output.
2. `git diff --stat` summary.
3. `pytest -m "not slow"` summary line.
4. The five new test files' pass/fail count.
5. The result of `python -c "from rectify.core.provenance import detect_current_cluster; print(detect_current_cluster())"` (should print "m1").
6. Confirmation that `grep -rn register_staged rectify/` returns NO matches (you deleted it).
7. Confirmation that `rectify/utils/provenance.py` is otherwise unchanged (only `register_staged` removed; `ProvenanceTracker` class still works).
8. Any TODOs you left for Opus before Commit B.

---

## 7. Time budget

- Read briefing, verify Commit A landed, read Opus-pre-written `path_resolver.py` and `skip_check.py`: ~30 min.
- Boilerplate modules (cluster, hashing, cli, sidecar): ~90-120 min.
- Tests: ~90 min.
- Delete `register_staged` + callsites: ~10 min.
- Pytest + iterate: ~30 min.
- **Total: ~4-5 hours.**

**End of briefing.**
