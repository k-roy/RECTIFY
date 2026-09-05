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

    Auto-detection only fires inside a recognized HPC batch job (SLURM, SGE/UGE,
    PBS). On local workstations — including macOS, where ``$TMPDIR`` is always
    set by the OS — this returns ``None`` so callers can write directly to the
    output filesystem with no rsync staging overhead. Users who want explicit
    scratch staging on a workstation can pass ``--scratch-dir DIR``.

    Priority order inside an HPC job (first that exists wins):
      1. $SCRATCH       — Preferred: persistent per-user scratch (Hoffman2, TACC,
                          and similar HPC systems). BAM files survive for post-job
                          inspection and job resumption. Auto-purged after ~90 days.
      2. $SLURM_TMPDIR  — Node-local tmpdir on some SLURM clusters.
      3. $TMPDIR        — POSIX generic; auto-cleaned at job end (not stable).

    Returns None if no scratch filesystem is detected (or if not in an HPC job).

    Notes:
        - $SCRATCH is preferred over $TMPDIR because rectify writes large intermediate
          BAMs that users may want to inspect after the job completes, and because
          job resumption (skipping re-alignment) relies on files surviving past the job.
        - $TMPDIR is appropriate only for purely transient files created and consumed
          within a single tool invocation.
        - Always clean up job-specific subdirectories at job end.
    """
    if not is_hpc_job():
        return None
    for var in ('SCRATCH', 'SLURM_TMPDIR', 'TMPDIR'):
        val = os.environ.get(var)
        if val:
            p = Path(val)
            if p.exists():
                return p
    logger.warning(
        "get_scratch_dir(): running in a batch job but no scratch directory "
        "found ($SCRATCH, $SLURM_TMPDIR, $TMPDIR unset or non-existent). "
        "All I/O will go directly to the output filesystem, which may cause "
        "severe NFS contention under concurrent array tasks."
    )
    return None


#: Minimum free space (GiB) a scratch filesystem must have before rectify will stage onto it.
#: An alignment + consensus run writes tens of GB of intermediate BAMs, and a full scratch
#: fails LATE and confusingly: on Hoffman2 (`/u/scratch` at 97 % in 2026-09) a near-full quota
#: produces an empty stderr, a low maxvmem and every array task dying inside a minute
#: (memory `reference-h2-scratch-quota-silent-kill`). Overridable with $RECTIFY_MIN_SCRATCH_GIB.
MIN_SCRATCH_FREE_GIB = 20.0


def scratch_free_gib(path: Path) -> Optional[float]:
    """Free space on *path*'s filesystem in GiB, or None if it cannot be determined.

    Uses ``shutil.disk_usage``, which reports the FILESYSTEM, not a per-user quota — a scratch
    with room but no quota left still looks free here. The check is therefore a guard against
    the common case (a genuinely full scratch), not a guarantee.
    """
    import shutil as _shutil
    try:
        return _shutil.disk_usage(str(path)).free / (1024 ** 3)
    except OSError:
        return None


def _scratch_has_room(path: Path, min_free_gib: Optional[float] = None) -> bool:
    """True when *path* has at least *min_free_gib* free (or free space is unknown)."""
    if min_free_gib is None:
        try:
            min_free_gib = float(os.environ.get('RECTIFY_MIN_SCRATCH_GIB', MIN_SCRATCH_FREE_GIB))
        except ValueError:
            min_free_gib = MIN_SCRATCH_FREE_GIB
    free = scratch_free_gib(path)
    if free is None:
        return True
    if free < min_free_gib:
        logger.warning(
            "scratch %s has only %.1f GiB free (< %.1f GiB); NOT staging there. "
            "Intermediates will be written to the output directory instead. "
            "Pass --scratch-dir DIR to choose another filesystem, or raise the floor with "
            "$RECTIFY_MIN_SCRATCH_GIB.", path, free, min_free_gib,
        )
        return False
    return True


def make_job_scratch_dir(prefix: str = 'rectify') -> Optional[Path]:
    """
    Create a unique per-job scratch directory and return its path.

    The directory name encodes the job ID and array task ID (portable across
    SLURM, UGE/SGE, and PBS) so concurrent array tasks never collide.
    Returns None if no scratch filesystem is available.

    Example::

        scratch = make_job_scratch_dir('rectify_correct')
        if scratch:
            work_bam = scratch / 'sample.multialigned.bam'
            shutil.copy(oak_bam, work_bam)
            # ... run correction on scratch ...
            sync_to_oak(scratch, oak_output_dir)
            shutil.rmtree(scratch)

    Args:
        prefix: Directory name prefix (default: 'rectify')

    Returns:
        Path to created scratch dir, or None if scratch is unavailable **or too full**
        (below ``MIN_SCRATCH_FREE_GIB``) — the caller then writes to the output filesystem.
    """
    scratch_base = get_scratch_dir()
    if scratch_base is None:
        return None

    job_id = get_job_id()
    task_id = get_task_id()
    if not _scratch_has_room(scratch_base):
        return None
    scratch_dir = scratch_base / f'{prefix}_{job_id}_{task_id}'
    scratch_dir.mkdir(parents=True, exist_ok=True)
    return scratch_dir


def resolve_scratch_dir(prefix: str = 'rectify', base_dir: Optional[Path] = None) -> Path:
    """
    Return a writable scratch subdirectory for intermediate BAM I/O.

    Unlike ``make_job_scratch_dir``, this function **always** returns a valid
    ``Path`` — callers should always clean up with ``shutil.rmtree``.

    Resolution order when *base_dir* is ``None``:

    1. ``$SCRATCH``           — persistent per-user scratch (Hoffman2 / Sherlock / TACC)
    2. ``$SLURM_TMPDIR``      — node-local tmpdir on some SLURM clusters
    3. ``$TMPDIR``            — POSIX generic; auto-cleaned at job end
    4. ``tempfile.mkdtemp()`` — always-available fallback (uses ``/tmp``)

    Args:
        prefix:   Directory name prefix used in the subdirectory name.
        base_dir: Explicit base path (from ``--scratch-dir``). When set, a
                  per-job subdirectory is created under *base_dir* instead of
                  the auto-detected path.

    Returns:
        Path to a freshly-created, per-job scratch subdirectory.
    """
    import tempfile

    if base_dir is not None:
        scratch_base: Optional[Path] = Path(base_dir)
        scratch_base.mkdir(parents=True, exist_ok=True)
    else:
        scratch_base = get_scratch_dir()

    if scratch_base is None:
        tmp = tempfile.mkdtemp(prefix=f'{prefix}_')
        return Path(tmp)

    job_id = get_job_id()
    task_id = get_task_id()
    d = scratch_base / f'{prefix}_{job_id}_{task_id}'
    d.mkdir(parents=True, exist_ok=True)
    return d


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
            *.mapPacBio.bam, *.gapmm2.bam, etc.) but keep the rectified BAM
            and the Step-4 *.corrected_polya.bam. Takes effect only when
            exclude_bam is False.

    Notes:
        - Always call this before rmtree(scratch_dir).
        - The ``drs_trim/`` and ``consensus_ckpt/`` subdirectories are always
          excluded — they hold large regenerable intermediates (trim FASTQ,
          consensus checkpoint state) that should never reach the durable
          output filesystem.
        - Falls back to a shutil-based copy if rsync is unavailable.
    """
    oak_dir.mkdir(parents=True, exist_ok=True)
    cmd = ['rsync', '-rlL', f'{scratch_dir}/', str(oak_dir) + '/']
    # Always-skip transient subdirs and aligner-internal scratch.
    # - drs_trim/, consensus_ckpt/: regenerable intermediate state
    # - *.uLTRA_ultra_tmp/: uLTRA prep_splicing index — rebuilt per run
    # - *.paf: aligner PAF dumps — debugging aid, regenerable
    cmd += [
        '--exclude=drs_trim/',
        '--exclude=consensus_ckpt/',
        '--exclude=*_ultra_tmp/',
        '--exclude=*.paf',
    ]
    if exclude_bam:
        cmd += ['--exclude=*.bam', '--exclude=*.bai']
    elif exclude_aligner_bams:
        # Keep the durable BAMs and discard per-aligner BAMs.  rsync filter
        # rules: first match wins, so list every durable include before the
        # broad *.bam exclude.  Durable set:
        #   *.rectified.bam        — FINAL corrected output (was corrected_consensus.bam)
        #   *.multialigned.bam     — pre-correction merged alignment artifact the
        #                            run reuse-gate keys on (was *.rectified.bam)
        #   corrected_consensus.bam — back-compat symlink to *.rectified.bam
        #   *.corrected_polya.bam  — DRS Step-4 output
        cmd += [
            '--include=*.rectified.bam',
            '--include=*.rectified.bam.bai',
            '--include=*.multialigned.bam',
            '--include=*.multialigned.bam.bai',
            '--include=corrected_consensus.bam',
            '--include=corrected_consensus.bam.bai',
            '--include=*.corrected_polya.bam',
            '--include=*.corrected_polya.bam.bai',
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
        _TRANSIENT_DIRS = {'drs_trim', 'consensus_ckpt'}

        def _should_skip(p: Path) -> bool:
            if p.is_dir() and (p.name in _TRANSIENT_DIRS or p.name.endswith('_ultra_tmp')):
                return True
            if p.suffix == '.paf':
                return True
            if exclude_bam and p.suffix in ('.bam', '.bai'):
                return True
            if exclude_aligner_bams:
                n = p.name
                is_bam_like = n.endswith('.bam') or n.endswith('.bai')
                is_durable_bam = (
                    n.endswith('.rectified.bam')
                    or n.endswith('.rectified.bam.bai')
                    or n.endswith('.multialigned.bam')
                    or n.endswith('.multialigned.bam.bai')
                    or n == 'corrected_consensus.bam'
                    or n == 'corrected_consensus.bam.bai'
                    or n.endswith('.corrected_polya.bam')
                    or n.endswith('.corrected_polya.bam.bai')
                )
                if is_bam_like and not is_durable_bam:
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
