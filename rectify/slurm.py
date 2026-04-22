#!/usr/bin/env python3
"""
SLURM integration utilities for RECTIFY.

Provides SLURM-aware CPU detection, thread limit management, and scratch
filesystem staging to optimize performance on HPC clusters.

## Thread limits (CRITICAL)

Setting thread limits prevents account bans on HPC clusters. Libraries like
numpy, sklearn, and pydeseq2 auto-spawn threads that can exceed SLURM
allocation if not constrained. Always call set_thread_limits() before
importing numpy or sklearn.

## Scratch staging

On HPC systems, $SCRATCH typically provides high-bandwidth local storage.
Staging large BAMs to $SCRATCH before correction significantly reduces wall
time and avoids I/O contention across concurrent array tasks.

Recommended pattern in SLURM scripts::

    from rectify.slurm import make_job_scratch_dir, copy_outputs_to_oak
    scratch = make_job_scratch_dir()
    if scratch:
        # stage BAM, run on scratch, copy back
        ...
    else:
        # fall back to direct Oak I/O
        ...

Or use the built-in staging via run_command._run_single_sample() when
--use-scratch is passed on the CLI.

Author: Kevin R. Roy
Date: 2026-03
"""

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def get_available_cpus(default: int = 1) -> int:
    """
    Get number of available CPUs, portable across SLURM, UGE/SGE, PBS, and interactive.

    Priority order:
    1. SLURM_CPUS_PER_TASK  (SLURM)
    2. NSLOTS               (UGE/SGE — set by scheduler to requested -pe slots)
    3. PBS_NUM_PPN          (PBS/Torque — processors per node)
    4. LOKY_MAX_CPU_COUNT   (user override, any scheduler)
    5. os.cpu_count() / 2   (conservative default for shared systems)
    6. default parameter

    Args:
        default: Fallback value if no CPU count can be determined

    Returns:
        Number of CPUs to use
    """
    for var in ('SLURM_CPUS_PER_TASK', 'NSLOTS', 'PBS_NUM_PPN', 'LOKY_MAX_CPU_COUNT'):
        val = os.environ.get(var)
        if val:
            try:
                return int(val)
            except ValueError:
                pass  # Malformed value, try next

    # Fall back to system CPU count (halved for safety on shared systems)
    cpu_count = os.cpu_count()
    if cpu_count:
        return max(1, cpu_count // 2)

    return default


def get_job_id() -> str:
    """
    Return the current job ID, portable across SLURM, UGE/SGE, PBS, and local runs.

    Priority: SLURM_JOB_ID → JOB_ID (UGE/SGE) → PBS_JOBID → os.getpid()
    """
    return (
        os.environ.get('SLURM_JOB_ID')
        or os.environ.get('JOB_ID')       # UGE/SGE
        or os.environ.get('PBS_JOBID')    # PBS/Torque
        or str(os.getpid())               # Interactive fallback
    )


def get_task_id() -> str:
    """
    Return the current array task ID, portable across schedulers.

    Priority: SLURM_ARRAY_TASK_ID → SGE_TASK_ID (UGE/SGE) → PBS_ARRAY_INDEX → '0'
    """
    return (
        os.environ.get('SLURM_ARRAY_TASK_ID')
        or os.environ.get('SGE_TASK_ID')       # UGE/SGE
        or os.environ.get('PBS_ARRAY_INDEX')   # PBS/Torque
        or '0'
    )


def set_thread_limits(n_threads: Optional[int] = None) -> int:
    """
    Set thread limits for common parallel libraries.

    CRITICAL: Must be called BEFORE importing numpy, sklearn, etc.

    Sets:
    - OMP_NUM_THREADS (OpenMP)
    - OPENBLAS_NUM_THREADS (OpenBLAS)
    - MKL_NUM_THREADS (Intel MKL)
    - LOKY_MAX_CPU_COUNT (joblib/sklearn)
    - NUMEXPR_MAX_THREADS (numexpr)

    Args:
        n_threads: Number of threads. If None, uses get_available_cpus()

    Returns:
        Number of threads set
    """
    if n_threads is None:
        n_threads = get_available_cpus()

    thread_str = str(n_threads)

    # Set all thread limit environment variables
    os.environ['OMP_NUM_THREADS'] = thread_str
    os.environ['OPENBLAS_NUM_THREADS'] = thread_str
    os.environ['MKL_NUM_THREADS'] = thread_str
    os.environ['LOKY_MAX_CPU_COUNT'] = thread_str
    os.environ['NUMEXPR_MAX_THREADS'] = thread_str

    return n_threads


def get_scratch_dir() -> Optional[Path]:
    """
    Return the best available high-bandwidth scratch directory for this job.

    Priority order (first that exists wins):
      1. $SCRATCH       — Preferred: persistent per-user scratch (Hoffman2, TACC,
                          and similar HPC systems). BAM files survive for post-job inspection and job
                          resumption. Auto-purged after ~90 days.
      2. $SLURM_TMPDIR  — Node-local tmpdir on some SLURM clusters.
      3. $TMPDIR        — POSIX generic; auto-cleaned at job end (not stable).

    Returns None if no scratch filesystem is detected.

    Notes:
        - $SCRATCH is preferred over $TMPDIR because rectify writes large intermediate
          BAMs that users may want to inspect after the job completes, and because
          job resumption (skipping re-alignment) relies on files surviving past the job.
        - $TMPDIR is appropriate only for purely transient files created and consumed
          within a single tool invocation.
        - Always clean up job-specific subdirectories at job end.
    """
    for var in ('SCRATCH', 'SLURM_TMPDIR', 'TMPDIR'):
        val = os.environ.get(var)
        if val:
            p = Path(val)
            if p.exists():
                return p
    if get_job_id() != str(os.getpid()):  # running under a real scheduler
        logger.warning(
            "get_scratch_dir(): running in a batch job but no scratch directory "
            "found ($SCRATCH, $SLURM_TMPDIR, $TMPDIR unset or non-existent). "
            "All I/O will go directly to the output filesystem, which may cause "
            "severe NFS contention under concurrent array tasks."
        )
    return None


def make_job_scratch_dir(prefix: str = 'rectify') -> Optional[Path]:
    """
    Create a unique per-job scratch directory and return its path.

    The directory name encodes the job ID and array task ID (portable across
    SLURM, UGE/SGE, and PBS) so concurrent array tasks never collide.
    Returns None if no scratch filesystem is available.

    Example::

        scratch = make_job_scratch_dir('rectify_correct')
        if scratch:
            work_bam = scratch / 'sample.rectified.bam'
            shutil.copy(oak_bam, work_bam)
            # ... run correction on scratch ...
            sync_to_oak(scratch, oak_output_dir)
            shutil.rmtree(scratch)

    Args:
        prefix: Directory name prefix (default: 'rectify')

    Returns:
        Path to created scratch dir, or None if scratch unavailable.
    """
    scratch_base = get_scratch_dir()
    if scratch_base is None:
        return None

    job_id = get_job_id()
    task_id = get_task_id()
    scratch_dir = scratch_base / f'{prefix}_{job_id}_{task_id}'
    scratch_dir.mkdir(parents=True, exist_ok=True)
    return scratch_dir


def sync_to_oak(
    scratch_dir: Path,
    oak_dir: Path,
    exclude_bam: bool = False,
    exclude_aligner_bams: bool = False,
) -> None:
    """
    Copy all outputs from a scratch directory back to Oak using rsync.

    Uses rsync -rlL (recursive, copy symlinks as files, dereference
    symlinks) to copy only new/changed files without propagating symlinks.

    Args:
        scratch_dir: Source directory on scratch filesystem.
        oak_dir: Destination directory on Oak (created if needed).
        exclude_bam: If True, skip all *.bam and *.bai files (use when the
            BAM already lives on Oak and was only staged temporarily).
        exclude_aligner_bams: If True, skip per-aligner BAMs (*.minimap2.bam,
            *.mapPacBio.bam, *.gapmm2.bam, etc.) but keep the rectified BAM.
            Takes effect only when exclude_bam is False.

    Notes:
        - Always call this before rmtree(scratch_dir).
        - Falls back to shutil.copytree if rsync is unavailable.
    """
    oak_dir.mkdir(parents=True, exist_ok=True)
    cmd = ['rsync', '-rlL', f'{scratch_dir}/', str(oak_dir) + '/']
    if exclude_bam:
        cmd += ['--exclude=*.bam', '--exclude=*.bai']
    elif exclude_aligner_bams:
        # Keep only *.rectified.bam; discard per-aligner BAMs.
        # rsync rule matching: first matching rule wins.
        cmd += [
            '--include=*.rectified.bam',
            '--include=*.rectified.bam.bai',
            '--exclude=*.bam',
            '--exclude=*.bai',
        ]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error("sync_to_oak: rsync failed (exit %d): %s", e.returncode, e)
        raise RuntimeError(f"sync_to_oak failed — rsync exited with status {e.returncode}") from e
    except FileNotFoundError:
        # rsync binary not installed — fall back to shutil with a warning
        logger.warning(
            "sync_to_oak: rsync not found on PATH; falling back to shutil.copy2. "
            "Install rsync for reliable Oak syncing."
        )
        def _should_skip(p: Path) -> bool:
            if exclude_bam and p.suffix in ('.bam', '.bai'):
                return True
            if exclude_aligner_bams:
                n = p.name
                is_bam_like = n.endswith('.bam') or n.endswith('.bai')
                is_rectified = n.endswith('.rectified.bam') or n.endswith('.rectified.bam.bai')
                if is_bam_like and not is_rectified:
                    return True
            return False

        def _copy_tree(src_dir: Path, dst_dir: Path) -> None:
            dst_dir.mkdir(parents=True, exist_ok=True)
            for src in src_dir.iterdir():
                if _should_skip(src):
                    continue
                dst = dst_dir / src.name
                if src.is_dir():
                    _copy_tree(src, dst)
                else:
                    shutil.copy2(src, dst)

        _copy_tree(scratch_dir, oak_dir)


def is_slurm_job() -> bool:
    """Check if currently running inside a SLURM job."""
    return 'SLURM_JOB_ID' in os.environ


def is_hpc_job() -> bool:
    """Check if currently running inside a recognized HPC batch job (any scheduler)."""
    return bool(
        os.environ.get('SLURM_JOB_ID')   # SLURM
        or os.environ.get('JOB_ID')       # UGE/SGE
        or os.environ.get('PBS_JOBID')    # PBS/Torque
    )


def get_slurm_info() -> dict:
    """
    Get SLURM job information if available.

    Returns:
        Dict with job_id, array_task_id, cpus, mem, etc.
        Empty dict if not in SLURM job.
    """
    if not is_slurm_job():
        return {}

    return {
        'job_id': os.environ.get('SLURM_JOB_ID'),
        'array_task_id': os.environ.get('SLURM_ARRAY_TASK_ID'),
        'cpus': os.environ.get('SLURM_CPUS_PER_TASK'),
        'mem': os.environ.get('SLURM_MEM_PER_NODE'),
        'partition': os.environ.get('SLURM_JOB_PARTITION'),
        'nodelist': os.environ.get('SLURM_NODELIST'),
    }
