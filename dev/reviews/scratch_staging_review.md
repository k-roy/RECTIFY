# Review Handoff: Scratch Staging for `rectify run-all`

**Date:** 2026-05-19
**Author:** Kevin R. Roy
**Reviewer scope:** correctness, safety, edge-case coverage, and portability of the scratch-staging implementation

**Status (2026-05-19, post-review):** review findings in
`scratch_staging_review_findings_20260519.md`. BLOCKER (sample_id[:12]
collision) + macOS `$TMPDIR` always-on + rsync filter leaks + Step-4
NFS-sort + `_run_single_sample` try/finally gap are all addressed.

---

## What was implemented

`rectify run-all` previously wrote all intermediate BAM I/O (multi-aligner output, unsorted corrected BAM, pysam sort temp files) directly to the NFS output directory. On HPC clusters with concurrent array tasks this causes severe I/O contention. The fix routes these intermediate files through a fast scratch filesystem and rsyncs only the durable outputs back to NFS.

### New CLI flag

```
--scratch-dir DIR   Base directory for intermediate BAM I/O. Auto-detected from
                    $SCRATCH, $SLURM_TMPDIR, or $TMPDIR when running inside an
                    HPC batch job. No staging overhead on local workstations when
                    no HPC scratch is found. Durable outputs (corrected_reads.tsv,
                    bedgraphs, report) always stay in --output-dir.
```

`--use-scratch` is now a suppressed deprecated alias (was the old flag, had no effect).

---

## Files changed

| File | What changed |
|---|---|
| `rectify/slurm.py` | Added `resolve_scratch_dir()` (always returns a Path, `tempfile` fallback) and `sync_to_oak()` (rsync with BAM filtering options). `make_job_scratch_dir()` was pre-existing and unchanged. |
| `rectify/core/commands/run_command.py` | Added `--scratch-dir` argument; deprecated `--use-scratch`. |
| `rectify/core/commands/run/single_sample.py` | Both `_run_single_sample` and `_process_one_sample` updated with scratch staging logic (see below). |

---

## Design decisions to review

### 1. Auto-detect vs always-on

The two scratch-detection functions serve different purposes:

- `make_job_scratch_dir(prefix)` → `Optional[Path]`: returns `None` when no HPC scratch env vars (`$SCRATCH`, `$SLURM_TMPDIR`, `$TMPDIR`) are present. **Used for auto-detection.** Zero overhead on local workstations.
- `resolve_scratch_dir(prefix, base_dir=None)` → `Path`: always returns a valid path, falling back to `tempfile.mkdtemp()` if nothing is available. **Used only when `--scratch-dir` is explicitly passed**, so the user gets exactly the directory they asked for.

The dispatch in both `_process_one_sample` and `_run_single_sample` looks like:

```python
_scratch_arg = getattr(args, 'scratch_dir', None)
if _scratch_arg is not None:
    _work = resolve_scratch_dir(..., base_dir=_scratch_arg)  # always a Path
    _using_scratch = True
else:
    _scratch = make_job_scratch_dir(...)  # None on laptop
    if _scratch is not None:
        _work = _scratch
        _using_scratch = True
    else:
        _work = sample_output   # write directly to NFS, no overhead
        _using_scratch = False
```

**Question for reviewer:** _(answered in review — non-issue.)_ When `--scratch-dir` is passed, `resolve_scratch_dir` calls `Path(base_dir).mkdir(parents=True, exist_ok=True)`. There is no tempfile fallback in that branch; the only way to reach `tempfile.mkdtemp()` is when `base_dir is None` AND `get_scratch_dir()` returns None. A bad explicit `--scratch-dir` raises `PermissionError`/`FileNotFoundError` on mkdir, which is the right behavior.

### 2. Scratch directory naming

Scratch subdirs are named `<prefix>_<job_id>_<task_id>`. The IDs are portable across schedulers:

- `SLURM_JOB_ID` → `JOB_ID` (UGE/SGE) → `PBS_JOBID` → `os.getpid()` for interactive
- `SLURM_ARRAY_TASK_ID` → `SGE_TASK_ID` → `PBS_ARRAY_INDEX` → `'0'`

Concurrent array tasks (e.g., 50-sample manifest split into a SLURM array) each get a unique subdirectory and never collide.

### 3. DRS Step 4 ordering — critical invariant

DRS mode runs:
- **Step 0** (trim poly-A): BAM → trimmed FASTQ (scratch)
- **Step 1** (align): trimmed FASTQ → per-aligner BAMs (scratch)
- **Step 2** (correct): per-aligner BAMs → corrected_reads.tsv (scratch)
- **Early TSV sync**: copy corrected_reads.tsv scratch → NFS immediately (so analysis reads from canonical NFS path)
- **Step 3** (analyze): reads corrected_reads.tsv from NFS
- **Step 4** (poly-A restore): reads per-aligner BAMs from scratch — these must still be present

The early TSV sync happens BEFORE the final scratch sync+rmtree. Step 4 needs per-aligner BAMs (still on scratch at that point). The final `sync_to_oak` + `rmtree` happens AFTER Step 4 completes.

**Verify:** In `_run_single_sample`, `per_aligner_bams` is populated during alignment (Step 1) and used in Step 4 (line ~620). The `_aligner_bams_for_restore` filter `if Path(str(v)).exists()` guards against the case where a BAM wasn't written. Confirm that per-aligner BAMs are still on scratch (not already cleaned up) when Step 4 runs.

### 4. `_process_one_sample` vs `_run_single_sample`

These are parallel implementations for multi-sample and single-sample paths respectively. Both implement the same scratch logic. Key differences:

- `_process_one_sample` runs inside a `ThreadPoolExecutor` — scratch must be thread-safe. Per-job naming (`<prefix>_<job_id>_<task_id>`) ensures isolation.
- `_process_one_sample` has a `try/except/finally` around the entire body. The `finally` calls `rmtree` only when `_using_scratch and _work.exists()`.
- `_run_single_sample` handles cleanup via:
  - Manual `_shutil.rmtree(scratch_dir, ...)` on the correction error path (line ~551)
  - Manual `_shutil.rmtree(scratch_dir, ...)` on the skip-alignment error path (line ~459)
  - Guarded `sync_to_oak` + `rmtree` in the normal exit path
  - Junction aggregation and DRS Step 4 are individually wrapped in `try/except` so exceptions don't escape without cleanup

**Gap:** `_run_single_sample` does NOT have a `try/finally` wrapper around the full body. If the final `sync_to_oak` itself raises (rsync non-zero), scratch is not cleaned. This is arguably acceptable (user needs to investigate the sync failure and may want the data), but a reviewer should flag if they disagree.

### 5. Sync filter: `exclude_aligner_bams`

`sync_to_oak` supports three sync modes:

| `exclude_bam` | `exclude_aligner_bams` | Effect |
|---|---|---|
| `False` | `False` | sync everything |
| `True` | — | skip all `.bam`/`.bai` |
| `False` | `True` | keep `*.rectified.bam` + `.bai`; skip per-aligner BAMs |

Default for `run-all` (when `--keep-aligner-bams` is not set): `exclude_aligner_bams=True`. Per-aligner BAMs (minimap2, mapPacBio, gapmm2) are discarded after consensus selection to save NFS space; only the rectified BAM is synced back.

When `--bam-dir DIR` is set, per-aligner BAMs went directly to `DIR` (not scratch), so the `exclude_aligner_bams` flag is forced to `False` in the scratch→NFS sync — there are no per-aligner BAMs on scratch to exclude.

### 6. `rsync -rlL` vs `shutil` fallback

`sync_to_oak` uses `rsync -rlL` (recursive, copy symlinks as files, dereference symlinks). The flags do NOT include `-u`, so files are not skipped based on age — rsync overwrites when size/mtime differs. Falls back to `shutil.copy2` with manual filtering if `rsync` is not on PATH. The rsync command uses a trailing slash on `scratch_dir/` to copy contents (not the directory itself), matching rsync convention.

**Verify:** `rsync -rlL` without `--delete` means files already on NFS from a previous partial run are NOT removed. This is intentional — a re-run should add/overwrite, not delete. Confirm this matches expected resumption behavior.

---

## Known gaps / not-yet-done

1. **No end-to-end test on cluster.** The test suite (`pytest -m "not slow"`) passes (1002/1004, 2 pre-existing failures in `test_corrected_consensus_tiebreaker.py` unrelated to this change). But no test exercises the full `run-all` pipeline with actual BAMs through scratch. A smoke test on H2 or Sherlock with a small real sample (e.g., `by4742_drs_wt_rep1` first 10k reads) is the right next step.

2. **`resolve_scratch_dir` tempfile fallback on `--scratch-dir`.** If the user passes `--scratch-dir /nonexistent/path` and the path's parent doesn't exist, `resolve_scratch_dir` falls back to `tempfile.mkdtemp()` rather than failing loudly. The user gets scratch staging in `/tmp` instead of an error. This may be surprising.

3. **Multi-sample manifest scratch isolation.** The `_using_scratch` flag per sample worker is correct, but if `--scratch-dir` points to a shared directory (e.g., `$SCRATCH` itself rather than a per-job subdir), two workers for different samples could collide on the prefix-based naming. The current naming is `rectify_<sample_id[:12]>_<job_id>_<task_id>`, which is safe for array jobs (different `task_id`) but not for intra-process threads under a single job ID. The thread safety relies on the sample_id being unique within a manifest.

4. **Pre-existing test failures.** Two tests in `test_corrected_consensus_tiebreaker.py` fail with `AttributeError: 'float' object has no attribute 'split'` in `corrected_consensus.py:830`. These are pre-existing and unrelated to this PR.

---

## How to verify the implementation

### Unit check (M1, no cluster needed)

```python
import os, tempfile
# Simulate no HPC scratch
for v in ('SCRATCH', 'SLURM_TMPDIR', 'TMPDIR', 'SLURM_JOB_ID'):
    os.environ.pop(v, None)
from rectify.slurm import make_job_scratch_dir
assert make_job_scratch_dir('test') is None  # should be None on laptop

# Simulate $SCRATCH available
os.environ['SCRATCH'] = tempfile.mkdtemp()
result = make_job_scratch_dir('test')
assert result is not None and result.exists()
```

### Cluster smoke test (H2 or Sherlock)

```bash
# H2 — interactive node
qrsh -l h_rt=1:00:00,h_data=8G -pe shared 4
source ~/.bashrc
conda activate rectify

# Use a small DRS sample (~10k reads is enough)
rectify run-all /u/project/guillom/shared/raw/by4742_polya_drs_2025/wt_rep1.bam \
    --drs --Scer \
    -o /u/scratch/$USER/rectify_scratch_test/ \
    --threads 4

# Verify:
# 1. Output dir has corrected_reads.tsv, report.html, etc.
ls /u/scratch/$USER/rectify_scratch_test/
# 2. Scratch was cleaned up (should be empty / not exist)
ls $SCRATCH/rectify_single_*/   # should fail (dir was removed)
```

---

## Diff summary

```
rectify/slurm.py                              +76 lines (resolve_scratch_dir, sync_to_oak)
rectify/core/commands/run_command.py          +13 lines (--scratch-dir arg, deprecate --use-scratch)
rectify/core/commands/run/single_sample.py   ~100 lines net change across both functions
```

Total: ~190 lines added/changed. No new dependencies. Existing `make_job_scratch_dir` and `get_scratch_dir` in `slurm.py` are unchanged.
