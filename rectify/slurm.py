#!/usr/bin/env python3
"""
SLURM integration utilities for RECTIFY.

Provides SLURM-aware CPU detection and thread limit management to prevent
oversubscription of cluster resources.

CRITICAL: Setting thread limits prevents account bans on Sherlock.
Libraries like numpy, sklearn, and pydeseq2 auto-spawn threads that can
exceed SLURM allocation if not constrained.

Author: Kevin R. Roy
Date: 2026-03
"""

import os
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
