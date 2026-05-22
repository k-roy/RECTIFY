#!/usr/bin/env python3
"""
Region-parallel BAM writer for RECTIFY.

Provides ``write_corrected_bam_parallel()``, which applies the same per-read
mutations as the single-threaded ``write_corrected_bam()`` in ``bam_writer.py``
but divides the work across N parallel processes — one per coord-ordered region.

Per-read mutations are a faithful mirror of ``bam_writer.py:write_corrected_bam``
lines 244-368.  Do NOT refactor them — correctness requires that the mutation
sequence is identical to the single-threaded path.

Worker design:
- Each worker opens its own pysam handle (pysam handles do not survive fork).
- Workers write an unsorted region BAM, sort it with pysam.sort, then touch an
  .ok sentinel to signal durable completion.
- The main thread merges all sorted region BAMs with pysam.merge (heap-merge,
  no final sort needed since every input is coord-sorted).
- Unmapped/secondary/supplementary reads from the last region's worker are
  handled in a separate serial pass after region workers complete.

Author: Kevin R. Roy
"""

import logging
import json
import os
import pickle
import random
import shutil
import string
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from .bam_writer import (
    apply_corrected_edits_to_read,
    _load_corrections_from_tsv,
)
from .regions import RegionPlan, plan_regions

logger = logging.getLogger(__name__)


def _fsync_parent_dir(path: Path) -> None:
    try:
        dir_fd = os.open(str(path.parent), os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(dir_fd)
    finally:
        os.close(dir_fd)


def _atomic_write_text(path: Path, text: str) -> None:
    tmp_path = path.with_name(path.name + '.tmp')
    with open(tmp_path, 'w') as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp_path, path)
    _fsync_parent_dir(path)


def _plan_complete(plan: RegionPlan) -> bool:
    return plan.ok_sentinel.exists() and plan.region_bam.exists()


def _atomic_pickle(path: Path, value) -> None:
    tmp_path = path.with_name(path.name + '.tmp')
    with open(tmp_path, 'wb') as fh:
        pickle.dump(value, fh, protocol=pickle.HIGHEST_PROTOCOL)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp_path, path)
    _fsync_parent_dir(path)


def _chunk_plans(plans: List[RegionPlan], n_batches: int) -> List[List[RegionPlan]]:
    batches: List[List[RegionPlan]] = [[] for _ in range(max(1, n_batches))]
    for idx, plan in enumerate(plans):
        batches[idx % len(batches)].append(plan)
    return [batch for batch in batches if batch]


# ---------------------------------------------------------------------------
# Per-read mutation helper (mirrors bam_writer.write_corrected_bam body)
# ---------------------------------------------------------------------------

def _apply_corrections_to_read(
    read: pysam.AlignedSegment,
    correction: Optional[Dict],
    genome: Optional[Dict[str, str]],
) -> bool:
    """Backward-compatible wrapper around the shared corrected-read helper."""
    return apply_corrected_edits_to_read(read, correction, genome)


# ---------------------------------------------------------------------------
# Per-region worker (module-level so multiprocessing can pickle it)
# ---------------------------------------------------------------------------

def _process_region_for_bam_write(
    plan: RegionPlan,
    input_bam_path: str,
    corrections: Dict[str, Dict],
    genome: Optional[Dict[str, str]],
) -> Dict:
    """Process a single region: apply corrections and emit a sorted region BAM.

    Module-level (not a method) so it is picklable for multiprocessing.Pool.
    Workers open their own pysam handle — pysam handles do not survive fork.

    Returns a dict of per-region stats:
        region_id, n_reads_in, n_reads_out, n_reads_skipped_dedup, wall_seconds
    """
    t0 = time.monotonic()
    region_id = plan.region_id

    # Resume via sentinel: if a prior run completed this region durably, skip.
    if _plan_complete(plan):
        logger.info("[region %s] skip: sentinel exists", region_id)
        # Return stats with zeros (actual counts are lost, but sentinel implies success)
        return {
            "region_id": region_id,
            "n_reads_in": 0,
            "n_reads_out": 0,
            "n_reads_skipped_dedup": 0,
            "wall_seconds": 0.0,
            "resumed": True,
        }
    if plan.ok_sentinel.exists() and not plan.region_bam.exists():
        logger.warning(
            "[region %s] ignoring stale sentinel without region BAM: %s",
            region_id,
            plan.ok_sentinel,
        )

    unsorted_path = plan.tmp_dir / f"{region_id}.{os.getpid()}.unsorted.bam.tmp"
    sorted_tmp_path = plan.region_bam.with_name(
        f"{plan.region_bam.name}.{os.getpid()}.tmp"
    )

    n_reads_in = 0
    n_reads_out = 0
    n_reads_skipped_dedup = 0

    try:
        bam_in = pysam.AlignmentFile(input_bam_path, 'rb')
        try:
            bam_out = pysam.AlignmentFile(
                str(unsorted_path), 'wb', header=bam_in.header
            )
            try:
                for read in bam_in.fetch(plan.chrom, plan.start, plan.end):
                    n_reads_in += 1

                    # Region-boundary dedup: only process reads whose alignment
                    # START (reference_start) falls within [plan.start, plan.end).
                    # Reads that START in an adjacent region are that region's
                    # responsibility; they will appear in both fetches but we only
                    # process each read once.
                    if not read.is_unmapped:
                        if (read.reference_start < plan.start
                                or read.reference_start >= plan.end):
                            n_reads_skipped_dedup += 1
                            continue

                    correction = corrections.get(read.query_name) if not (
                        read.is_unmapped or read.is_secondary or read.is_supplementary
                    ) else None

                    _apply_corrections_to_read(read, correction, genome)
                    bam_out.write(read)
                    n_reads_out += 1
            finally:
                bam_out.close()
        finally:
            bam_in.close()

        # Sort to a temp BAM first. Only the final os.replace publishes this
        # region for resume/merge.
        if sorted_tmp_path.exists():
            sorted_tmp_path.unlink()
        pysam.sort('-@', '1', '-o', str(sorted_tmp_path), str(unsorted_path))
        unsorted_path.unlink()

        # fsync the region BAM to ensure durable write before touching the sentinel.
        with open(sorted_tmp_path, 'rb') as fh:
            os.fsync(fh.fileno())
        os.replace(sorted_tmp_path, plan.region_bam)
        _fsync_parent_dir(plan.region_bam)

        # Write the sentinel AFTER the BAM replacement. Existence implies the
        # region BAM was durably committed.
        _atomic_write_text(plan.ok_sentinel, 'ok\n')
    except Exception:
        for path in (unsorted_path, sorted_tmp_path):
            try:
                if path.exists():
                    path.unlink()
            except OSError:
                pass
        raise

    wall_seconds = time.monotonic() - t0
    logger.info(
        "[region %s] done: %d in, %d out, %d dedup-skipped, %.1fs",
        region_id, n_reads_in, n_reads_out, n_reads_skipped_dedup, wall_seconds,
    )
    return {
        "region_id": region_id,
        "n_reads_in": n_reads_in,
        "n_reads_out": n_reads_out,
        "n_reads_skipped_dedup": n_reads_skipped_dedup,
        "wall_seconds": wall_seconds,
        "resumed": False,
    }


# ---------------------------------------------------------------------------
# Unmapped-reads tail pass (serial, runs after region workers complete)
# ---------------------------------------------------------------------------

def _write_unmapped_reads(
    input_bam_path: str,
    unmapped_bam_path: Path,
    corrections: Dict[str, Dict],
    genome: Optional[Dict[str, str]],
) -> Dict:
    """Write all unmapped reads from the input BAM to a separate output BAM.

    The legacy write_corrected_bam writes unmapped reads unchanged (lines
    247-250 of bam_writer.py). We mirror that here: unmapped reads are
    written without applying per-read mutations.

    Returns stats: n_reads_in, n_reads_out.
    """
    n_reads_in = 0
    n_reads_out = 0

    with pysam.AlignmentFile(input_bam_path, 'rb') as bam_in, \
         pysam.AlignmentFile(str(unmapped_bam_path), 'wb', header=bam_in.header) as bam_out:

        for read in bam_in.fetch(until_eof=True):
            if not read.is_unmapped:
                continue
            n_reads_in += 1
            # Mirror legacy: unmapped reads are written unchanged.
            bam_out.write(read)
            n_reads_out += 1

    logger.info(
        "[unmapped] %d unmapped reads written to %s",
        n_reads_out, unmapped_bam_path,
    )
    return {"n_reads_in": n_reads_in, "n_reads_out": n_reads_out}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def write_corrected_bam_parallel(
    input_bam_path: str,
    corrected_tsv_path: str,
    output_bam_path: str,
    n_threads: int = 4,
    genome: Optional[Dict[str, str]] = None,
    tmp_dir: Optional[str] = None,
    allow_resume: bool = True,
    **kwargs,
) -> Dict:
    """Write a corrected BAM using N parallel region workers.

    Produces output byte-equivalent to ``write_corrected_bam()`` (modulo
    sort tie order, which is non-deterministic in pysam for reads with the
    same coordinate). The two functions share the same per-read mutation
    sequence; they differ only in how reads are partitioned and merged.

    Args:
        input_bam_path:      Path to the original input BAM (must be indexed).
        corrected_tsv_path:  Path to the corrected_reads.tsv from rectify correct.
        output_bam_path:     Destination BAM path (overwritten if exists).
        n_threads:           Number of parallel region workers.
        genome:              Optional pre-loaded genome dict for CIGAR surgery.
        tmp_dir:             Scratch directory for intermediate region BAMs.
                             Defaults to a new tempdir under $TMPDIR.
        allow_resume:        If True, skip regions whose .ok sentinel exists.
        **kwargs:            Recognised: ``keep_tmp`` (bool, default False).

    Returns:
        Dict with aggregated stats:
            n_regions, n_reads_in_total, n_reads_out_total, wall_seconds_total,
            stats_per_region (list of per-region stat dicts).
    """
    t_total = time.monotonic()
    keep_tmp = kwargs.get("keep_tmp", False)

    # Resolve or create tmp_dir.
    if tmp_dir is None:
        rand_suffix = ''.join(random.choices(string.ascii_lowercase, k=8))
        tmp_dir_path = Path(tempfile.gettempdir()) / "rectify_regions" / f"{os.getpid()}_{rand_suffix}"
    else:
        tmp_dir_path = Path(tmp_dir)
    tmp_dir_path.mkdir(parents=True, exist_ok=True)

    # Pre-flight disk space check: need at least 1.5× input BAM size free.
    input_size = os.path.getsize(input_bam_path)
    free = shutil.disk_usage(tmp_dir_path).free
    if free < 1.5 * input_size:
        raise RuntimeError(
            f"Insufficient disk space in {tmp_dir_path}: "
            f"{free / 1e9:.1f} GB free, need at least {1.5 * input_size / 1e9:.1f} GB "
            f"(1.5× input size of {input_size / 1e9:.1f} GB)."
        )

    # Load corrections table once in the parent process.
    corrections = _load_corrections_from_tsv(corrected_tsv_path)
    logger.info(
        "write_corrected_bam_parallel: loaded %d corrected positions from %s",
        len(corrections), corrected_tsv_path,
    )

    # Plan regions.
    plans = plan_regions(input_bam_path, tmp_dir_path)
    logger.info(
        "write_corrected_bam_parallel: %d regions, %d threads, tmp=%s",
        len(plans), n_threads, tmp_dir_path,
    )

    # If resume is disabled, remove any existing sentinels so workers re-run.
    if not allow_resume:
        for plan in plans:
            if plan.ok_sentinel.exists():
                plan.ok_sentinel.unlink()

    # Dispatch region workers.
    worker_args = [
        (plan, input_bam_path, corrections, genome)
        for plan in plans
    ]

    stats_per_region: List[Dict] = []

    enable_subprocess_workers = os.environ.get('RECTIFY_ENABLE_PARALLEL_BAM_WRITER') == '1'
    if n_threads == 1 or len(plans) <= 1 or not enable_subprocess_workers:
        # Sequential path: avoids unsafe multiprocessing/subprocess behavior by
        # default. Set RECTIFY_ENABLE_PARALLEL_BAM_WRITER=1 to opt into external
        # worker processes after validating the local runtime.
        if n_threads > 1 and len(plans) > 1 and not enable_subprocess_workers:
            logger.warning(
                "write_corrected_bam_parallel: n_threads=%d requested, but "
                "parallel BAM writer workers are disabled by default because "
                "they can abort in some pysam/OpenMP runtimes. Running region "
                "workers sequentially; set RECTIFY_ENABLE_PARALLEL_BAM_WRITER=1 "
                "to opt in.",
                n_threads,
            )
        for args in worker_args:
            stats_per_region.append(_process_region_for_bam_write(*args))
    else:
        effective_workers = min(n_threads, len(plans))
        logger.info(
            "write_corrected_bam_parallel: launching %d clean worker subprocesses",
            effective_workers,
        )
        stats_per_region.extend(_run_region_batches_in_subprocesses(
            plans,
            input_bam_path,
            corrections,
            genome,
            tmp_dir_path,
            effective_workers,
        ))

    # Verify all regions completed successfully.
    for plan in plans:
        if not plan.ok_sentinel.exists():
            raise RuntimeError(
                f"Region {plan.region_id} ({plan.chrom}:{plan.start}-{plan.end}) "
                f"did not produce a sentinel file — worker may have failed silently."
            )

    # Unmapped-reads tail pass (serial).
    unmapped_bam_path = tmp_dir_path / "unmapped.bam"
    unmapped_stats = _write_unmapped_reads(
        input_bam_path, unmapped_bam_path, corrections, genome
    )

    # Build merge input list: all region BAMs + unmapped BAM (if non-empty).
    merge_inputs = [str(p.region_bam) for p in plans]
    if unmapped_stats["n_reads_out"] > 0:
        merge_inputs.append(str(unmapped_bam_path))

    # Merge sorted region BAMs.  pysam.merge heap-merges coord-sorted inputs
    # so the output is fully sorted without a final sort pass.
    logger.info(
        "write_corrected_bam_parallel: merging %d BAM shards -> %s",
        len(merge_inputs), output_bam_path,
    )
    pysam.merge(
        '-@', str(n_threads),
        '-p',   # preserve @PG chain
        '-c',   # consolidate colliding @PG IDs
        '-f',   # overwrite output if exists
        str(output_bam_path),
        *merge_inputs,
    )

    # Index the merged output.
    pysam.index(str(output_bam_path))
    logger.info("write_corrected_bam_parallel: indexed %s", output_bam_path)

    # Clean up tmp_dir unless keep_tmp is set.
    if not keep_tmp:
        shutil.rmtree(tmp_dir_path, ignore_errors=True)
    else:
        logger.info("write_corrected_bam_parallel: keeping tmp dir: %s", tmp_dir_path)

    wall_total = time.monotonic() - t_total

    # Aggregate stats.
    n_reads_in_total = sum(s.get("n_reads_in", 0) for s in stats_per_region)
    n_reads_out_total = sum(s.get("n_reads_out", 0) for s in stats_per_region)
    n_reads_out_total += unmapped_stats["n_reads_out"]

    result = {
        "n_regions": len(plans),
        "n_reads_in_total": n_reads_in_total,
        "n_reads_out_total": n_reads_out_total,
        "wall_seconds_total": wall_total,
        "stats_per_region": stats_per_region,
    }
    logger.info(
        "write_corrected_bam_parallel: done in %.1fs, %d reads out",
        wall_total, n_reads_out_total,
    )
    return result


def _process_region_for_bam_write_star(args: Tuple) -> Dict:
    """Unpacking wrapper for pool.imap (which passes a single iterable arg)."""
    return _process_region_for_bam_write(*args)


def _run_region_batches_in_subprocesses(
    plans: List[RegionPlan],
    input_bam_path: str,
    corrections: Dict[str, Dict],
    genome: Optional[Dict[str, str]],
    tmp_dir_path: Path,
    n_workers: int,
) -> List[Dict]:
    """Run region workers as clean Python subprocesses, not forked Pool children."""
    worker_dir = tmp_dir_path / "worker_state"
    worker_dir.mkdir(parents=True, exist_ok=True)
    corrections_pkl = worker_dir / "corrections.pkl"
    genome_pkl = worker_dir / "genome.pkl"
    _atomic_pickle(corrections_pkl, corrections)
    _atomic_pickle(genome_pkl, genome)

    batches = _chunk_plans(plans, n_workers)
    worker_env = os.environ.copy()
    for key in (
        'OMP_NUM_THREADS',
        'OPENBLAS_NUM_THREADS',
        'MKL_NUM_THREADS',
        'VECLIB_MAXIMUM_THREADS',
        'NUMEXPR_NUM_THREADS',
        'LOKY_MAX_CPU_COUNT',
    ):
        worker_env[key] = '1'
    worker_env.setdefault('KMP_INIT_AT_FORK', 'FALSE')
    worker_env['RECTIFY_IMPORT_VISUALIZE'] = '0'

    processes = []
    try:
        for batch_idx, batch in enumerate(batches):
            plans_pkl = worker_dir / f"batch_{batch_idx:03d}.plans.pkl"
            stats_json = worker_dir / f"batch_{batch_idx:03d}.stats.json"
            _atomic_pickle(plans_pkl, batch)
            cmd = [
                sys.executable,
                '-m',
                'rectify.core.bam.bam_writer_parallel',
                '--worker-batch',
                '--input-bam',
                input_bam_path,
                '--plans-pkl',
                str(plans_pkl),
                '--corrections-pkl',
                str(corrections_pkl),
                '--genome-pkl',
                str(genome_pkl),
                '--stats-json',
                str(stats_json),
            ]
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=worker_env,
            )
            processes.append((batch_idx, proc, stats_json))

        stats: List[Dict] = []
        for batch_idx, proc, stats_json in processes:
            stdout, stderr = proc.communicate()
            if proc.returncode != 0:
                details = stderr.strip() or stdout.strip() or f"exit code {proc.returncode}"
                raise RuntimeError(f"BAM writer worker batch {batch_idx} failed: {details}")
            with open(stats_json, 'r') as fh:
                stats.extend(json.load(fh))
    except Exception:
        for _batch_idx, proc, _stats_json in processes:
            if proc.poll() is None:
                proc.terminate()
        raise

    return sorted(stats, key=lambda item: item.get('region_id', ''))


def _worker_batch_main(argv: Optional[List[str]] = None) -> int:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--worker-batch', action='store_true')
    parser.add_argument('--input-bam', required=True)
    parser.add_argument('--plans-pkl', required=True)
    parser.add_argument('--corrections-pkl', required=True)
    parser.add_argument('--genome-pkl', required=True)
    parser.add_argument('--stats-json', required=True)
    args = parser.parse_args(argv)

    with open(args.plans_pkl, 'rb') as fh:
        plans = pickle.load(fh)
    with open(args.corrections_pkl, 'rb') as fh:
        corrections = pickle.load(fh)
    with open(args.genome_pkl, 'rb') as fh:
        genome = pickle.load(fh)

    stats = [
        _process_region_for_bam_write(plan, args.input_bam, corrections, genome)
        for plan in plans
    ]
    _atomic_write_text(Path(args.stats_json), json.dumps(stats, sort_keys=True) + '\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(_worker_batch_main())
