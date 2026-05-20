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
import multiprocessing as mp
import os
import random
import shutil
import string
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from .bam_writer import (
    _decode_eq_seq_inplace,
    _load_corrections_from_tsv,
)
from .read_edits import (
    clip_read_to_corrected_3prime,
    extend_read_5prime_for_junction_rescue,
    extend_read_3prime_for_overcall_rescue,
    extend_read_3prime_for_softclip_rescue,
    softclip_intronic_tail_5prime,
    reroute_intronic_tail_5prime_via_junction,
    realign_exon_blocks,
    reanchor_5prime_for_rescue,
    _hardclip_trailing_a_run,
)
from .regions import RegionPlan, plan_regions

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-read mutation helper (mirrors bam_writer.write_corrected_bam body)
# ---------------------------------------------------------------------------

def _apply_corrections_to_read(
    read: pysam.AlignedSegment,
    correction: Optional[Dict],
    genome: Optional[Dict[str, str]],
) -> bool:
    """Apply per-read corrections in the same sequence as write_corrected_bam.

    Returns True if the read was modified, False otherwise.

    This is a module-level function (not a method) so it can be used by both
    the region worker and the unmapped tail pass.
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        # These are written unchanged — mirror bam_writer.py:247-250.
        return False

    # Decode '='-compressed SEQ (see _decode_eq_seq_inplace docstring).
    if genome is not None:
        _decode_eq_seq_inplace(read, genome)

    if correction is None:
        return False

    modified = False

    # 5'-edge reanchor pre-pass: apply before realign_exon_blocks.
    if genome is not None and correction.get('reanchor_clip_len', 0) > 0:
        modified |= reanchor_5prime_for_rescue(read, genome, anchor_min_run=10)

    # Homopolymer CIGAR surgery: re-align exon blocks.
    if genome is not None:
        modified |= realign_exon_blocks(read, genome)

    # 5' junction rescue: extend soft-clip to exon 1 (Cat3).
    if correction['five_prime_rescued'] and correction['five_prime_position'] is not None:
        modified |= extend_read_5prime_for_junction_rescue(
            read,
            correction['five_prime_position'],
            correction['five_prime_soft_clip'],
            correction['strand'],
            exon_cigar_str=correction.get('five_prime_exon_cigar', ''),
            upstream_trim=correction.get('five_prime_upstream_trim', 0),
        )

    # 5' junction rescue: reroute intronic M ops to exon 1 (Cases 1/2/2b/4).
    _icp = correction.get('five_prime_intron_clip_pos', -1)
    _exon_cig = correction.get('five_prime_exon_cigar', '')
    if (_icp >= 0 and _exon_cig and correction.get('five_prime_rescued')
            and correction['five_prime_position'] is not None):
        modified |= reroute_intronic_tail_5prime_via_junction(
            read,
            clip_boundary=_icp,
            five_prime_position=correction['five_prime_position'],
            exon_cigar_str=_exon_cig,
            strand=correction['strand'],
        )
    elif _icp >= 0 and not _exon_cig and correction.get('five_prime_rescued'):
        modified |= softclip_intronic_tail_5prime(
            read,
            clip_boundary=_icp,
            strand=correction['strand'],
        )

    # Cat2 soft-clip rescue: extend 3' alignment outward into homopolymer.
    if correction.get('sc_rescued_seq'):
        modified |= extend_read_3prime_for_softclip_rescue(
            read,
            correction['strand'],
            correction['sc_homopolymer_extension'],
            correction['sc_rescued_seq'],
            correction['sc_original_softclip_len'],
            hard_clip=True,
        )

    # Over-call rescue: convert trailing soft-clip to aligned bases.
    if correction.get('oc_terminal_base'):
        modified |= extend_read_3prime_for_overcall_rescue(
            read,
            correction['strand'],
            correction['oc_homopolymer_extension'],
            correction['oc_overcall_count'],
            correction['oc_terminal_base'],
        )

    # 3' hard-clip to corrected position.
    modified |= clip_read_to_corrected_3prime(
        read, correction['corrected_3prime'], correction['strand']
    )

    # Hard-clip trailing genomic A-run at the 3' end.
    modified |= _hardclip_trailing_a_run(read, correction['strand'])

    # Tag the corrected 3' end position for IGV visibility.
    read.set_tag('cp', correction['corrected_3prime'])

    return modified


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
    if plan.ok_sentinel.exists():
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

    unsorted_path = plan.tmp_dir / f"{region_id}.unsorted.bam"

    n_reads_in = 0
    n_reads_out = 0
    n_reads_skipped_dedup = 0

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

    # Sort the unsorted region BAM to the final region BAM path.
    pysam.sort('-@', '1', '-o', str(plan.region_bam), str(unsorted_path))
    unsorted_path.unlink()

    # fsync the region BAM to ensure durable write before touching the sentinel.
    with open(plan.region_bam, 'rb') as fh:
        os.fsync(fh.fileno())

    # Touch the sentinel AFTER fsync — existence implies region_bam is durable.
    plan.ok_sentinel.touch()

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

    if n_threads == 1 or len(plans) <= 1:
        # Sequential path: avoids ProcessPool spawn cost for small inputs.
        for args in worker_args:
            stats_per_region.append(_process_region_for_bam_write(*args))
    else:
        effective_workers = min(n_threads, len(plans))
        with mp.Pool(processes=effective_workers) as pool:
            for result in pool.imap(_process_region_for_bam_write_star, worker_args):
                stats_per_region.append(result)

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
