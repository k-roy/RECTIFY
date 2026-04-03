#!/usr/bin/env python3
"""
SLURM integration utilities for RECTIFY.

Provides SLURM-aware CPU detection, thread limit management, and scratch
filesystem staging to optimize performance on HPC clusters.

## Thread limits (CRITICAL)

Setting thread limits prevents account bans on Sherlock. Libraries like
numpy, sklearn, and pydeseq2 auto-spawn threads that can exceed SLURM
allocation if not constrained. Always call set_thread_limits() before
importing numpy or sklearn.

## Scratch staging

On Sherlock, $SCRATCH has ~75 GB/s bandwidth vs Oak's shared NFS. Staging
large BAMs to $SCRATCH before correction significantly reduces wall time and
avoids I/O contention across concurrent array tasks.

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

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional


def get_available_cpus(default: int = 1) -> int:
    """
    Get number of available CPUs, respecting SLURM allocation.

    Priority order:
    1. SLURM_CPUS_PER_TASK (if running in SLURM job)
    2. LOKY_MAX_CPU_COUNT (if set by user)
    3. os.cpu_count() / 2 (conservative default for shared systems)
    4. default parameter

    Args:
        default: Fallback value if no CPU count can be determined

    Returns:
        Number of CPUs to use
    """
    # Check SLURM environment first
    slurm_cpus = os.environ.get('SLURM_CPUS_PER_TASK')
    if slurm_cpus:
        try:
            return int(slurm_cpus)
        except ValueError:
            pass  # Malformed value, fall through

    # Check LOKY_MAX_CPU_COUNT (user override)
    loky_cpus = os.environ.get('LOKY_MAX_CPU_COUNT')
    if loky_cpus:
        try:
            return int(loky_cpus)
        except ValueError:
            pass  # Malformed value, fall through

    # Fall back to system CPU count (halved for safety on shared systems)
    cpu_count = os.cpu_count()
    if cpu_count:
        return max(1, cpu_count // 2)

    return default


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
      1. $SCRATCH       — Sherlock per-user scratch (~75 GB/s aggregate)
      2. $SLURM_TMPDIR  — node-local tmpdir on some clusters
      3. $TMPDIR        — generic temp dir fallback

    Returns None if no scratch filesystem is detected.

    Notes:
        - On Sherlock, $SCRATCH is auto-purged after 90 days.
        - Always clean up job-specific subdirectories at job end.
        - Do NOT use $L_SCRATCH across nodes (node-local only).
    """
    for var in ('SCRATCH', 'SLURM_TMPDIR', 'TMPDIR'):
        val = os.environ.get(var)
        if val:
            p = Path(val)
            if p.exists():
                return p
    return None


def make_job_scratch_dir(prefix: str = 'rectify') -> Optional[Path]:
    """
    Create a unique per-job scratch directory and return its path.

    The directory name encodes the SLURM job ID and array task ID so
    concurrent array tasks never collide. Returns None if no scratch
    filesystem is available (falls back to Oak I/O).

    Example::

        scratch = make_job_scratch_dir('rectify_correct')
        if scratch:
            work_bam = scratch / 'sample.consensus.bam'
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

    job_id = os.environ.get('SLURM_JOB_ID', 'local')
    task_id = os.environ.get('SLURM_ARRAY_TASK_ID', '0')
    scratch_dir = scratch_base / f'{prefix}_{job_id}_{task_id}'
    scratch_dir.mkdir(parents=True, exist_ok=True)
    return scratch_dir


def sync_to_oak(scratch_dir: Path, oak_dir: Path, exclude_bam: bool = False) -> None:
    """
    Copy all outputs from a scratch directory back to Oak using rsync.

    Uses rsync -a (archive mode: preserves permissions, timestamps,
    symlinks) to copy only new/changed files.

    Args:
        scratch_dir: Source directory on scratch filesystem.
        oak_dir: Destination directory on Oak (created if needed).
        exclude_bam: If True, skip *.bam and *.bai files (use when the
            BAM already lives on Oak and was only staged temporarily).

    Notes:
        - Always call this before rmtree(scratch_dir).
        - Falls back to shutil.copytree if rsync is unavailable.
    """
    oak_dir.mkdir(parents=True, exist_ok=True)
    cmd = ['rsync', '-a', f'{scratch_dir}/', str(oak_dir) + '/']
    if exclude_bam:
        cmd += ['--exclude=*.bam', '--exclude=*.bai']
    try:
        subprocess.run(cmd, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        # rsync not available — fall back to shutil
        for src in scratch_dir.iterdir():
            if exclude_bam and src.suffix in ('.bam', '.bai'):
                continue
            dst = oak_dir / src.name
            if src.is_dir():
                shutil.copytree(src, dst, dirs_exist_ok=True)
            else:
                shutil.copy2(src, dst)


def is_slurm_job() -> bool:
    """Check if currently running inside a SLURM job."""
    return 'SLURM_JOB_ID' in os.environ


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
