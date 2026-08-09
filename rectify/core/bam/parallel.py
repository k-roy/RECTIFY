#!/usr/bin/env python3
"""
Parallel / streaming BAM processing dispatch for RECTIFY.

Wraps the per-read correction core (``bam_processor.correct_read_3prime``)
with region-based parallel workers, streaming TSV writers, and a streaming +
parallel combination for very large BAMs.

Workers use a spawn-style multiprocessing context by default so pysam/htslib
state opened during parent-side BAM scans is never inherited by child
processes. Each worker loads its own genome/model state once at process start
and opens its own BAM handle per region.

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import gzip
import json
import logging
import multiprocessing as mp
import os

import pysam

from ..correct.indel_corrector import VariantAwareHomopolymerRescue
from ..polya.polya_model import PolyAModel, load_model as load_polya_model
from ...slurm import get_available_cpus
from ...utils.genome import load_genome
from .bam_processor import correct_read_3prime, _load_netseq
from .output import CORRECTION_TSV_HEADER, correction_result_to_tsv_row, write_output_tsv
from .processing_stats import ProcessingStats
from ..position_index import write_position_index

logger = logging.getLogger(__name__)


_REGION_WORKER_STATE: Dict = {}


def _warmup_hp_ed_numba() -> None:
    """Pre-fork JIT warmup for the 3'SS-rescue DP kernel.

    ⚠️ **Near-inert under the DEFAULT start method.** ``_get_bam_worker_context``
    defaults to ``spawn``, and spawned workers do NOT inherit the parent's
    compiled code -- they re-import and reload from numba's on-disk cache.  So
    under the default it is ``@njit(cache=True)`` that actually saves the
    workers, not this call.  It is kept because it DOES help when
    ``RECTIFY_BAM_MP_START_METHOD=fork``, and because it surfaces any JIT
    failure in the parent (where it is visible) rather than in N workers.

    Imported lazily: ``splice_aware_5prime`` is not a top-level dependency of
    this module and a module-level import would risk a cycle.  Best-effort by
    design -- a warmup failure must never break a correction run.  No-op when
    numba is absent.
    """
    try:
        from ..splice.splice_aware_5prime import _warmup_hp_ed_numba as _w
        _w()
    except Exception:                                   # pragma: no cover
        pass


def _get_bam_worker_context() -> mp.context.BaseContext:
    """Return the multiprocessing context used for pysam region workers.

    Linux defaults to ``fork``, which can inherit htslib C-level state from
    parent-side BAM scans and corrupt large indexed fetch workloads. ``spawn``
    gives each worker a fresh interpreter and fresh htslib state. The env var
    is intentionally narrow so HPC debugging can compare methods without a
    code patch.
    """
    method = os.environ.get('RECTIFY_BAM_MP_START_METHOD', 'spawn').strip().lower()
    valid = set(mp.get_all_start_methods())
    if method not in valid:
        logger.warning(
            "Invalid RECTIFY_BAM_MP_START_METHOD=%r; falling back to 'spawn' "
            "(available: %s)",
            method,
            ', '.join(sorted(valid)),
        )
        method = 'spawn'
    return mp.get_context(method)


def _init_region_worker_state(
    genome_path: str,
    polya_model_path: Optional[str],
    shared_kwargs: Dict,
) -> None:
    """Initialise per-process state for region workers.

    Holds only the fields shared across aligners in one ``run-all`` invocation:
    ``genome``, ``polya_model``, ``annotated_junctions``, ``pool_chrom_index``,
    ``gene_interval_trees``, ``netseq_dir``, and config flags. Per-aligner
    fields (``bam_path``, ``variant_aware_rescue``) are NOT stored here — they
    travel with each region task tuple so the same pool can be shared across
    multiple per-aligner correction calls (Commit C-II.0 worker state model).

    This avoids sending large objects such as the genome with every task while
    still keeping all pysam handles worker-local.
    """
    global _REGION_WORKER_STATE

    genome = load_genome(genome_path)
    polya_model = load_polya_model(Path(polya_model_path) if polya_model_path else None)
    _REGION_WORKER_STATE = dict(shared_kwargs)
    _REGION_WORKER_STATE['genome'] = genome
    _REGION_WORKER_STATE['polya_model'] = polya_model


def _process_region_worker_from_state(
    task: Tuple,
) -> Union[List[Dict], str]:
    """Pool entrypoint using process-local state installed by initializer.

    The task tuple is ``(region, bam_path, variant_aware_rescue)``. The first
    field is the genome interval; the latter two are the per-aligner sliver
    that must vary across calls submitted to a shared pool.
    """
    if not _REGION_WORKER_STATE:
        raise RuntimeError("Region worker state was not initialised")
    region, bam_path, variant_aware_rescue = task
    return _process_region_worker(
        region=region,
        bam_path=bam_path,
        variant_aware_rescue=variant_aware_rescue,
        **_REGION_WORKER_STATE,
    )


def _write_results_chunk(fh, results: List[Dict]):
    """Write a chunk of results to file handle."""
    for result in results:
        fh.write('\t'.join(correction_result_to_tsv_row(result)) + '\n')


def _region_tsv_path(checkpoint_dir: Path, region_idx: int) -> Path:
    return checkpoint_dir / f'region_{region_idx:04d}.tsv'


def _region_stats_path(checkpoint_dir: Path, region_idx: int) -> Path:
    return checkpoint_dir / f'region_{region_idx:04d}.stats.json'


def _region_done_path(checkpoint_dir: Path, region_idx: int) -> Path:
    return checkpoint_dir / f'region_{region_idx:04d}.done'


def _is_region_checkpoint_complete(checkpoint_dir: Path, region_idx: int) -> bool:
    return (
        _region_done_path(checkpoint_dir, region_idx).exists()
        and _region_tsv_path(checkpoint_dir, region_idx).exists()
    )


def _fsync_parent_dir(path: Path) -> None:
    """Best-effort parent-directory fsync after atomic checkpoint renames."""
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


def _stats_for_results(results: List[Dict]) -> ProcessingStats:
    stats = ProcessingStats()
    for result in results:
        stats.update_from_result(result)
    return stats


def _stats_from_dict(values: Dict) -> ProcessingStats:
    stats = ProcessingStats()
    for key in stats.to_dict():
        setattr(stats, key, int(values.get(key, 0) or 0))
    return stats


def _copy_prescan_counts(source: ProcessingStats, target: ProcessingStats) -> None:
    target.total_reads_in_bam = source.total_reads_in_bam
    target.reads_unmapped = source.reads_unmapped
    target.reads_secondary = source.reads_secondary
    target.reads_supplementary = source.reads_supplementary
    target.spikein_reads_filtered = source.spikein_reads_filtered


def _write_potential_variants(
    rescue: VariantAwareHomopolymerRescue,
    output_variants_path: Optional[str],
) -> None:
    if not output_variants_path:
        return
    variants = rescue.get_potential_variants()
    if not variants:
        return
    logger.info(f"  Writing {len(variants)} potential variants to {output_variants_path}")
    with open(output_variants_path, 'w') as f:
        header = [
            'chrom', 'position', 'ref_base', 'homopolymer_base',
            'total_reads', 'mismatch_fraction',
            'dominant_mismatch_base', 'dominant_mismatch_fraction',
        ]
        f.write('\t'.join(header) + '\n')
        for v in variants:
            row = [
                v['chrom'],
                str(v['position']),
                v['ref_base'],
                v['homopolymer_base'],
                str(v['total_reads']),
                f"{v['mismatch_fraction']:.3f}",
                v.get('dominant_mismatch_base', ''),
                f"{v.get('dominant_mismatch_fraction', 0.0):.3f}",
            ]
            f.write('\t'.join(row) + '\n')


def _regions_from_intervals(
    chrom_lengths: Dict[str, int],
    intervals_by_chrom: Dict[str, List[Tuple[int, int]]],
    *,
    min_gap_size: int,
    max_region_size: int = 100000,
) -> List[Tuple[str, int, int]]:
    regions: List[Tuple[str, int, int]] = []
    for chrom, chrom_len in chrom_lengths.items():
        if chrom_len <= max_region_size:
            regions.append((chrom, 0, chrom_len))
            continue

        intervals = sorted(intervals_by_chrom.get(chrom, []))
        if not intervals:
            continue

        merged = [list(intervals[0])]
        for start, end in intervals[1:]:
            if start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])

        gaps: List[Tuple[int, int]] = []
        prev_end = 0
        for seg_start, seg_end in merged:
            if seg_start - prev_end >= min_gap_size:
                gaps.append((prev_end, seg_start))
            prev_end = seg_end
        if chrom_len - prev_end >= min_gap_size:
            gaps.append((prev_end, chrom_len))

        if not gaps:
            n_sub = max(1, (chrom_len + max_region_size - 1) // max_region_size)
            sub_size = chrom_len // n_sub
            for i in range(n_sub):
                sub_start = i * sub_size
                sub_end = chrom_len if i == n_sub - 1 else (i + 1) * sub_size
                regions.append((chrom, sub_start, sub_end))
            continue

        current_start = 0
        for gap_start, gap_end in gaps:
            if gap_start > current_start:
                regions.append((chrom, current_start, gap_start))
            current_start = gap_end
        if current_start < chrom_len:
            regions.append((chrom, current_start, chrom_len))
    return regions


def _prescan_bam_once(
    bam_path: str,
    genome: Dict[str, str],
    *,
    min_gap_size: int,
    variant_aware: bool = False,
    output_variants_path: Optional[str] = None,
) -> Tuple[ProcessingStats, List[Tuple[str, int, int]], Optional[VariantAwareHomopolymerRescue]]:
    """Single parent-side BAM scan for stats, region planning, and optional variant scan."""
    stats = ProcessingStats()
    intervals_by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    rescue = None
    if variant_aware:
        rescue = VariantAwareHomopolymerRescue(
            min_variant_fraction=0.8,
            min_reads_for_variant_call=5,
            min_homopolymer_len=4,
            max_rescue_bases=3,
        )

    logger.info("Pre-scanning BAM once for read stats and region planning...")
    n_variant_scanned = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        chrom_lengths = dict(zip(bam.references, bam.lengths))
        for chrom in bam.references:
            intervals_by_chrom[chrom] = []

        for read in bam.fetch(until_eof=True):
            stats.total_reads_in_bam += 1
            if read.is_unmapped:
                stats.reads_unmapped += 1
                continue
            if read.is_secondary:
                stats.reads_secondary += 1
                continue
            if read.is_supplementary:
                stats.reads_supplementary += 1
                continue

            chrom = bam.get_reference_name(read.reference_id)
            if chrom and read.reference_start is not None and read.reference_end is not None:
                intervals_by_chrom.setdefault(chrom, []).append((read.reference_start, read.reference_end))

            if rescue is not None:
                strand = '-' if read.is_reverse else '+'
                rescue.scan_read(read, strand, genome, end='3prime')
                n_variant_scanned += 1
                if n_variant_scanned % 500000 == 0:
                    logger.info(f"  Variant-aware scan: {n_variant_scanned:,} reads...")

    if rescue is not None:
        logger.info(f"  Variant-aware scan total: {n_variant_scanned:,} reads")
        rescue.finalize_scan()
        variant_stats = rescue.get_statistics()
        logger.info(f"  Positions with mismatches: {variant_stats['total_positions_scanned']:,}")
        logger.info(f"  With sufficient coverage: {variant_stats['positions_with_sufficient_coverage']:,}")
        logger.info(f"  Potential variants: {variant_stats['potential_variants_detected']:,}")
        _write_potential_variants(rescue, output_variants_path)

    regions = _regions_from_intervals(
        chrom_lengths,
        intervals_by_chrom,
        min_gap_size=min_gap_size,
    )
    logger.info(f"  Total reads: {stats.total_reads_in_bam:,}")
    logger.info(f"  Unmapped: {stats.reads_unmapped:,}")
    logger.info(f"  Secondary: {stats.reads_secondary:,}")
    logger.info(f"  Supplementary: {stats.reads_supplementary:,}")
    return stats, regions, rescue


def _write_region_results_atomic(checkpoint_dir: Path, region_idx: int, results: List[Dict]) -> None:
    """Atomically commit one checkpointed region's rows and stats."""
    tsv_path = _region_tsv_path(checkpoint_dir, region_idx)
    tsv_tmp = tsv_path.with_name(tsv_path.name + '.tmp')
    with open(tsv_tmp, 'w') as fh:
        _write_results_chunk(fh, results)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tsv_tmp, tsv_path)
    _fsync_parent_dir(tsv_path)

    stats_path = _region_stats_path(checkpoint_dir, region_idx)
    stats_text = json.dumps(_stats_for_results(results).to_dict(), sort_keys=True) + '\n'
    _atomic_write_text(stats_path, stats_text)

    # The done sentinel is written last. A crash before this point leaves files
    # that are safe to overwrite on the next run, but not trusted as complete.
    _atomic_write_text(_region_done_path(checkpoint_dir, region_idx), 'done\n')


def _parse_tsv_result(fields: List[str], seen_read_ids: set) -> Dict:
    row = dict(zip(CORRECTION_TSV_HEADER, fields))

    def _int(name: str, default: int = 0) -> int:
        value = row.get(name, '')
        return int(value) if value not in ('', None) else default

    def _float(name: str, default: float = 0.0) -> float:
        value = row.get(name, '')
        return float(value) if value not in ('', None) else default

    corrections = row.get('correction_applied', '')
    correction_applied = [] if corrections in ('', 'none') else corrections.split(',')
    qc_flags = [flag for flag in row.get('qc_flags', '').split(',') if flag]
    read_id = row.get('read_id', '')
    is_primary = read_id not in seen_read_ids
    seen_read_ids.add(read_id)

    return {
        'read_id': read_id,
        'chrom': row.get('chrom', ''),
        'strand': row.get('strand', ''),
        'original_3prime': _int('original_3prime'),
        'corrected_3prime': _int('corrected_3prime'),
        'ambiguity_min': _int('ambiguity_min'),
        'ambiguity_max': _int('ambiguity_max'),
        'ambiguity_range': _int('ambiguity_range'),
        'polya_length': _int('polya_length'),
        'correction_applied': correction_applied,
        'confidence': row.get('confidence', 'high') or 'high',
        'qc_flags': qc_flags,
        'fraction': _float('fraction', 1.0),
        'five_prime_rescued': row.get('five_prime_rescued', '0') == '1',
        'is_primary_result': is_primary,
    }


def _output_tmp_path(output_path: str) -> Path:
    path = Path(output_path)
    return path.with_name(path.name + '.tmp')


def _rebuild_output_from_region_files(
    checkpoint_dir: Path,
    output_path: str,
    region_count: int,
    prescan_stats: ProcessingStats,
) -> Tuple[ProcessingStats, Dict[Tuple[str, int, str], float]]:
    """Rebuild final TSV and position counts from completed region checkpoints."""
    final_stats = ProcessingStats()
    _copy_prescan_counts(prescan_stats, final_stats)
    pos_counts: Dict[Tuple[str, int, str], float] = {}
    seen_read_ids: set = set()
    tmp_output = _output_tmp_path(output_path)

    opener = gzip.open if output_path.endswith('.gz') else open
    try:
        with opener(tmp_output, 'wt') as out_fh:
            out_fh.write('\t'.join(CORRECTION_TSV_HEADER) + '\n')
            for region_idx in range(region_count):
                if not _is_region_checkpoint_complete(checkpoint_dir, region_idx):
                    raise RuntimeError(f"Checkpoint region {region_idx} is not complete")

                stats_path = _region_stats_path(checkpoint_dir, region_idx)
                has_stats_sidecar = stats_path.exists()
                if has_stats_sidecar:
                    with open(stats_path, 'r') as stats_fh:
                        final_stats.merge(_stats_from_dict(json.load(stats_fh)))

                with open(_region_tsv_path(checkpoint_dir, region_idx), 'r') as region_fh:
                    for line in region_fh:
                        out_fh.write(line)
                        fields = line.rstrip('\n').split('\t')
                        if len(fields) != len(CORRECTION_TSV_HEADER):
                            continue
                        result = _parse_tsv_result(fields, seen_read_ids)
                        key = (result['chrom'], result['corrected_3prime'], result['strand'])
                        pos_counts[key] = pos_counts.get(key, 0.0) + float(result.get('fraction', 1.0))
                        if not has_stats_sidecar:
                            final_stats.update_from_result(result)

        os.replace(tmp_output, output_path)
        _fsync_parent_dir(Path(output_path))
    finally:
        try:
            if tmp_output.exists():
                tmp_output.unlink()
        except OSError:
            pass

    return final_stats, pos_counts


def _rebuild_pos_counts_from_partial(output_path: str, pos_counts: dict) -> int:
    """Re-read an existing partial output TSV to rebuild pos_counts. Returns read count."""
    n = 0
    try:
        with open(output_path, 'r') as f:
            header = f.readline().strip().split('\t')
            try:
                ci = header.index('chrom')
                pi = header.index('corrected_3prime')
                si = header.index('strand')
                fi = header.index('fraction')
            except ValueError:
                return 0
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= max(ci, pi, si, fi):
                    continue
                try:
                    key = (parts[ci], int(parts[pi]), parts[si])
                    frac = float(parts[fi]) if parts[fi].strip() else 1.0
                    pos_counts[key] = pos_counts.get(key, 0.0) + frac
                    n += 1
                except (ValueError, IndexError):
                    pass
    except OSError:
        pass
    return n


def _process_region_worker(
    region: Tuple[str, int, int],
    bam_path: str,
    genome: Dict[str, str],
    apply_atract: bool,
    apply_ag_mispriming: bool,
    ag_threshold: float = 0.65,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_dir: Optional[str] = None,
    variant_aware_rescue: Optional[VariantAwareHomopolymerRescue] = None,
    annotated_junctions: Optional[set] = None,
    pool_chrom_index: Optional[Dict] = None,
    apply_3ss_rescue: bool = True,
    gene_interval_trees: Optional[Dict] = None,
    polya_model: Optional[PolyAModel] = None,
    tmp_dir: Optional[str] = None,
    max_reads_for_variant_rescue: int = 500,
    dt_primed_cDNA: bool = False,
    ont_cDNA: bool = False,
    use_dorado_polya: bool = False,
    min_mapq: int = 0,
    min_aligned_length: int = 0,
    exclusion_detector: Optional['ExclusionRegionDetector'] = None,
) -> Union[List[Dict], str]:
    """
    Worker function to process a single region.

    Opens its own BAM handle to avoid multiprocessing issues.

    Args:
        region: Tuple of (chrom, start, end)
        bam_path: Path to BAM file
        genome: Pre-loaded genome dict
        apply_*: Correction module flags
        netseq_dir: Optional NET-seq BigWig directory
        variant_aware_rescue: Optional variant-aware rescue object (from first pass)
        tmp_dir: If set, pickle results to a temp file in this directory and
            return the file path string instead of the results list.  This
            avoids the OS pipe-buffer-full deadlock that occurs when a very
            large region (e.g. rDNA) produces hundreds of MB of results — a
            pipe write blocks when the buffer is full, causing a deadlock
            between the worker (waiting to finish the write) and the main
            process (waiting for the complete pickle before it starts reading).
            A file write goes to the OS page cache and returns immediately
            regardless of result size.

    Returns:
        List of correction result dicts for reads in region, or a file path
        string if tmp_dir is set.
    """
    import pickle as _pickle
    import tempfile as _tempfile

    chrom, start, end = region
    results = []

    # Open BAM for this region
    bam = pysam.AlignmentFile(bam_path, 'rb')

    # Load NET-seq if needed (per-worker to avoid file handle issues)
    netseq_loader = None
    if netseq_dir:
        netseq_loader = _load_netseq(netseq_dir)

    _rescue_reads_count = 0
    try:
        for read in bam.fetch(chrom, start, end):
            # Skip unmapped/secondary/supplementary
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            # Skip low-MAPQ reads (e.g. multi-mapping reads at MAPQ=0)
            if min_mapq > 0 and read.mapping_quality < min_mapq:
                continue
            # Skip reads with insufficient aligned length (e.g. reads mapping
            # only to a low-complexity T-tract with almost no unique sequence)
            if min_aligned_length > 0 and (read.reference_length or 0) < min_aligned_length:
                continue
            # Boundary deduplication: when a chromosome is split by coordinate
            # (not by a coverage gap), pysam.fetch returns reads that *overlap*
            # [start, end), so a read spanning a sub-region boundary would be
            # processed by both workers.  Only process reads that START in this
            # region — reads starting before 'start' belong to the previous worker.
            if read.reference_start < start:
                continue

            # Cap variant-aware rescue for high-depth regions (e.g. rDNA).
            # After max_reads_for_variant_rescue reads the variant dictionary is
            # already well-populated (min_reads_for_variant_call=5), so further
            # scan updates waste CPU.  Remaining reads fall back to plain rescue.
            if variant_aware_rescue is not None:
                if _rescue_reads_count < max_reads_for_variant_rescue:
                    _local_rescue = variant_aware_rescue
                    _rescue_reads_count += 1
                else:
                    _local_rescue = None
            else:
                _local_rescue = None

            # Apply corrections
            read_results = correct_read_3prime(
                read,
                genome,
                apply_atract=apply_atract,
                apply_ag_mispriming=apply_ag_mispriming,
                ag_threshold=ag_threshold,
                apply_polya_trim=apply_polya_trim,
                apply_indel_correction=apply_indel_correction,
                netseq_loader=netseq_loader,
                variant_aware_rescue=_local_rescue,
                annotated_junctions=annotated_junctions,
                pool_chrom_index=pool_chrom_index,
                apply_3ss_rescue=apply_3ss_rescue,
                gene_interval_trees=gene_interval_trees,
                polya_model=polya_model,
                dt_primed_cDNA=dt_primed_cDNA,
                ont_cDNA=ont_cDNA,
                use_dorado_polya=use_dorado_polya,
                exclusion_detector=exclusion_detector,
            )
            results.extend(read_results)

    finally:
        bam.close()
        if netseq_loader:
            netseq_loader.close()

    if tmp_dir is not None:
        # Write to a temp file; return just the path string through the pipe.
        # This avoids the pipe-buffer-full deadlock for large regions.
        fd, tmp_path = _tempfile.mkstemp(suffix='.pkl', dir=tmp_dir)
        with os.fdopen(fd, 'wb') as _fh:
            _pickle.dump(results, _fh, protocol=_pickle.HIGHEST_PROTOCOL)
        return tmp_path

    return results


def process_bam_file_parallel(
    bam_path: str,
    genome_path: str,
    n_threads: int = 0,
    apply_atract: bool = True,
    apply_ag_mispriming: bool = False,
    ag_threshold: float = 0.65,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_dir: Optional[str] = None,
    output_path: Optional[str] = None,
    max_reads: Optional[int] = None,
    min_gap_size: int = 10000,
    show_progress: bool = True,
    return_stats: bool = False,
    variant_aware: bool = False,
    variant_output_path: Optional[str] = None,
    annotated_junctions: Optional[set] = None,
    pool_chrom_index: Optional[Dict] = None,
    apply_3ss_rescue: bool = True,
    gene_interval_trees: Optional[Dict] = None,
    polya_model_path: Optional[str] = None,
    variant_scan_cache: Optional[str] = None,
    max_reads_for_variant_rescue: int = 500,
    dt_primed_cDNA: bool = False,
    ont_cDNA: bool = False,
    use_dorado_polya: bool = False,
    min_mapq: int = 0,
    min_aligned_length: int = 0,
    reuse_pool_container: Optional[list] = None,
    exclusion_detector: Optional['ExclusionRegionDetector'] = None,
) -> Union[List[Dict], Tuple[List[Dict], ProcessingStats]]:
    """
    Process BAM file with parallel region-based processing.

    Splits chromosomes into independent regions at coverage gaps,
    then processes regions in parallel.

    Args:
        bam_path: Path to BAM file
        genome_path: Path to genome FASTA
        n_threads: Number of threads (0 = auto-detect from SLURM/system)
        apply_atract: Apply A-tract ambiguity detection
        apply_ag_mispriming: Apply AG mispriming screening
        apply_polya_trim: Apply poly(A) tail trimming
        apply_indel_correction: Apply indel artifact correction
        netseq_dir: Optional directory with NET-seq BigWig files
        output_path: Optional output TSV path
        max_reads: Optional maximum number of reads (for testing)
        min_gap_size: Minimum gap size for region splitting
        show_progress: Show progress information
        return_stats: Return ProcessingStats alongside results
        variant_aware: Enable variant-aware homopolymer rescue (two-pass)
        variant_output_path: Optional path to write potential variants TSV

    Returns:
        List of correction result dicts, or tuple of (results, stats) if return_stats=True
    """
    # Auto-detect threads if not specified
    if n_threads <= 0:
        n_threads = get_available_cpus()

    logger.info(f"Using {n_threads} thread(s) for parallel processing")

    # Load genome for sequential correction and variant-aware prescan. In the
    # parallel path, spawned workers load clean per-process copies.
    logger.info(f"Loading genome from {genome_path}...")
    genome = load_genome(genome_path)

    # Load poly(A) model for sequential correction. In the parallel path,
    # spawned workers load clean per-process copies.
    polya_model = load_polya_model(Path(polya_model_path) if polya_model_path else None)
    if polya_model is not None:
        logger.info(f"Loaded poly(A) model from {polya_model_path}")

    # Run variant-aware scan if enabled (first pass).
    # When variant_scan_cache points to a pre-built pickle (from `rectify prescan`),
    # load it directly and skip the scan — used for chunked correction pipelines.
    import pickle as _pickle
    variant_aware_rescue = None
    if variant_aware:
        if variant_scan_cache and Path(variant_scan_cache).exists():
            logger.info(f"Loading pre-computed variant scan from cache: {variant_scan_cache}")
            with open(variant_scan_cache, 'rb') as _f:
                variant_aware_rescue = _pickle.load(_f)

    stats, regions, prescan_rescue = _prescan_bam_once(
        bam_path,
        genome,
        min_gap_size=min_gap_size,
        variant_aware=variant_aware and variant_aware_rescue is None,
        output_variants_path=variant_output_path,
    )
    if prescan_rescue is not None:
        variant_aware_rescue = prescan_rescue
    logger.info(f"  Found {len(regions)} processing regions")

    # Single-threaded fallback
    if n_threads == 1:
        logger.info("Processing regions sequentially...")
        all_results = []
        for i, region in enumerate(regions):
            chrom, start, end = region
            if show_progress:
                logger.info(f"  Region {i+1}/{len(regions)}: {chrom}:{start:,}-{end:,}")

            results = _process_region_worker(
                region, bam_path, genome,
                apply_atract, apply_ag_mispriming,
                ag_threshold,
                apply_polya_trim, apply_indel_correction,
                netseq_dir, variant_aware_rescue,
                annotated_junctions,
                pool_chrom_index,
                apply_3ss_rescue,
                gene_interval_trees,
                polya_model,
                max_reads_for_variant_rescue=max_reads_for_variant_rescue,
                dt_primed_cDNA=dt_primed_cDNA,
                ont_cDNA=ont_cDNA,
                use_dorado_polya=use_dorado_polya,
                min_mapq=min_mapq,
                min_aligned_length=min_aligned_length,
            )
            all_results.extend(results)

            if max_reads and len(all_results) >= max_reads:
                logger.info(f"  Reached max_reads limit ({max_reads})")
                break

        logger.info(f"Completed processing {len(all_results):,} reads")

        # Update stats from results
        for result in all_results:
            stats.update_from_result(result)

        if output_path:
            write_output_tsv(all_results, output_path)
            write_position_index(all_results, output_path)

        if return_stats:
            return all_results, stats
        return all_results

    # Multi-process processing. Use spawn-style workers so no htslib state from
    # parent-side BAM scans is inherited by the workers.
    logger.info(f"Processing {len(regions)} regions across {n_threads} workers...")

    # Commit C-II.0 worker state split: shared across aligners vs per-task.
    # The shared bundle is installed once via the pool initializer; per-task
    # fields (bam_path, variant_aware_rescue) are wrapped into each region
    # tuple so the same pool can later be shared across per-aligner correction
    # calls without re-init overhead.
    shared_kwargs = dict(
        apply_atract=apply_atract,
        apply_ag_mispriming=apply_ag_mispriming,
        ag_threshold=ag_threshold,
        apply_polya_trim=apply_polya_trim,
        apply_indel_correction=apply_indel_correction,
        netseq_dir=netseq_dir,
        annotated_junctions=annotated_junctions,
        pool_chrom_index=pool_chrom_index,
        apply_3ss_rescue=apply_3ss_rescue,
        gene_interval_trees=gene_interval_trees,
        max_reads_for_variant_rescue=max_reads_for_variant_rescue,
        dt_primed_cDNA=dt_primed_cDNA,
        ont_cDNA=ont_cDNA,
        use_dorado_polya=use_dorado_polya,
        min_mapq=min_mapq,
        min_aligned_length=min_aligned_length,
        exclusion_detector=exclusion_detector,
    )
    region_tasks = [(region, bam_path, variant_aware_rescue) for region in regions]

    # The parent no longer needs these large objects once worker initializers
    # will load their own copies. Releasing them also reduces RSS before the
    # spawned workers start.
    try:
        del genome
        del polya_model
    except NameError:
        pass
    import gc as _gc
    _gc.collect()

    all_results = []

    # Pool lifecycle — three cases based on reuse_pool_container:
    #   None        → create + destroy within this call (standard single-aligner path)
    #   empty list  → first aligner: create pool, store in container for caller to reuse
    #   non-empty   → subsequent aligner: reuse pool already stored by first aligner
    import contextlib as _contextlib
    if reuse_pool_container is not None and reuse_pool_container:
        _pool_cm = _contextlib.nullcontext(reuse_pool_container[0])
        logger.info("Reusing shared pool (%d workers).", n_threads)
    else:
        ctx = _get_bam_worker_context()
        logger.info(
            "Using multiprocessing start method '%s' for BAM region workers",
            ctx.get_start_method(),
        )
        # Compile the 3'SS-rescue DP kernel BEFORE forking so workers inherit
        # compiled code; otherwise each worker pays the ~1 s JIT independently.
        # Mirrors junction_refiner._warmup_numba_dp. No-op without numba.
        _warmup_hp_ed_numba()
        _new_pool = ctx.Pool(
            n_threads,
            initializer=_init_region_worker_state,
            initargs=(str(genome_path), polya_model_path, shared_kwargs),
        )
        if reuse_pool_container is not None:
            # First-of-shared: store for caller, skip lifecycle management here.
            reuse_pool_container.append(_new_pool)
            _pool_cm = _contextlib.nullcontext(_new_pool)
            logger.info(
                "Created shared pool (%d workers); stored for reuse by subsequent aligners.",
                n_threads,
            )
        else:
            _pool_cm = _new_pool  # Pool is its own context manager; terminates on __exit__

    with _pool_cm as pool:
        try:
            # Try to use tqdm for progress if available
            from tqdm import tqdm
            results_iter = tqdm(
                pool.imap(_process_region_worker_from_state, region_tasks),
                total=len(region_tasks),
                desc="Processing regions"
            )
        except ImportError:
            results_iter = pool.imap(_process_region_worker_from_state, region_tasks)

        for region_results in results_iter:
            all_results.extend(region_results)

            if max_reads and len(all_results) >= max_reads:
                logger.info(f"Reached max_reads limit ({max_reads})")
                break

    logger.info(f"Completed processing {len(all_results):,} reads")

    # Update stats from results
    for result in all_results:
        stats.update_from_result(result)

    # Write output if requested
    if output_path:
        write_output_tsv(all_results, output_path)
        write_position_index(all_results, output_path)

    if return_stats:
        return all_results, stats
    return all_results


# =============================================================================
# Streaming Output Mode (for memory efficiency on large BAMs)
# =============================================================================


def process_bam_streaming(
    bam_path: str,
    genome_path: str,
    output_path: str,
    chunk_size: int = 10000,
    apply_atract: bool = True,
    apply_ag_mispriming: bool = False,
    ag_threshold: float = 0.65,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_dir: Optional[str] = None,
    show_progress: bool = True,
    annotated_junctions: Optional[set] = None,
    pool_chrom_index: Optional[Dict] = None,
    apply_3ss_rescue: bool = True,
    gene_interval_trees: Optional[Dict] = None,
    polya_model_path: Optional[str] = None,
    dt_primed_cDNA: bool = False,
    ont_cDNA: bool = False,
    use_dorado_polya: bool = False,
    min_mapq: int = 0,
    min_aligned_length: int = 0,
) -> ProcessingStats:
    """
    Process BAM file with streaming output to minimize memory usage.

    Writes results directly to output file instead of accumulating in memory.
    Recommended for BAM files > 10GB.

    Args:
        bam_path: Input BAM path
        genome_path: Genome FASTA path
        output_path: Output TSV path (supports .gz compression)
        chunk_size: Number of reads to process before writing
        apply_*: Correction module flags
        netseq_dir: Optional NET-seq BigWig directory
        show_progress: Show progress information

    Returns:
        ProcessingStats object with comprehensive statistics
    """
    # Load genome
    logger.info(f"Loading genome from {genome_path}...")
    genome = load_genome(genome_path)

    # Load NET-seq if provided
    netseq_loader = None
    if netseq_dir:
        logger.info(f"Loading NET-seq data from {netseq_dir}...")
        netseq_loader = _load_netseq(netseq_dir)

    # Load poly(A) model if provided
    polya_model = load_polya_model(Path(polya_model_path) if polya_model_path else None)
    if polya_model is not None:
        logger.info(f"Loaded poly(A) model from {polya_model_path}")

    # Initialize comprehensive stats
    stats = ProcessingStats()

    # Position count accumulator for the compact index
    from collections import defaultdict as _defaultdict
    _pos_counts = _defaultdict(float)

    # Open output file (support gzip)
    if output_path.endswith('.gz'):
        out_fh = gzip.open(output_path, 'wt')
    else:
        out_fh = open(output_path, 'w')

    _failed = False
    try:
        # Write header
        out_fh.write('\t'.join(CORRECTION_TSV_HEADER) + '\n')

        # Open BAM
        bam = pysam.AlignmentFile(bam_path, 'rb')
        chunk = []
        _progress_interval = 100000  # Log every N reads
        _next_progress = _progress_interval
        import time as _time
        _t_start = _time.monotonic()
        _t_last = _t_start

        try:
            for read in bam:
                stats.total_reads_in_bam += 1

                # Track filtered reads
                if read.is_unmapped:
                    stats.reads_unmapped += 1
                    continue
                if read.is_secondary:
                    stats.reads_secondary += 1
                    continue
                if read.is_supplementary:
                    stats.reads_supplementary += 1
                    continue
                if min_mapq > 0 and read.mapping_quality < min_mapq:
                    continue
                if min_aligned_length > 0 and (read.reference_length or 0) < min_aligned_length:
                    continue

                read_results = correct_read_3prime(
                    read, genome,
                    apply_atract=apply_atract,
                    apply_ag_mispriming=apply_ag_mispriming,
                    ag_threshold=ag_threshold,
                    apply_polya_trim=apply_polya_trim,
                    apply_indel_correction=apply_indel_correction,
                    netseq_loader=netseq_loader,
                    annotated_junctions=annotated_junctions,
                    pool_chrom_index=pool_chrom_index,
                    apply_3ss_rescue=apply_3ss_rescue,
                    gene_interval_trees=gene_interval_trees,
                    polya_model=polya_model,
                    dt_primed_cDNA=dt_primed_cDNA,
                    ont_cDNA=ont_cDNA,
                    use_dorado_polya=use_dorado_polya,
                )
                chunk.extend(read_results)

                # Update comprehensive stats and position index accumulator
                for result in read_results:
                    stats.update_from_result(result)
                    _pos_counts[(result['chrom'], result['corrected_3prime'], result['strand'])] += float(result.get('fraction', 1.0))

                # Write chunk when full
                if len(chunk) >= chunk_size:
                    _write_results_chunk(out_fh, chunk)
                    chunk = []

                # Progress log: fires every _progress_interval BAM reads (not rows)
                if show_progress and stats.total_reads_in_bam >= _next_progress:
                    _t_now = _time.monotonic()
                    _elapsed = _t_now - _t_start
                    _rate_per_min = (stats.total_reads_in_bam / _elapsed * 60) if _elapsed > 0 else 0
                    logger.info(
                        f"  Processed {stats.total_reads_in_bam:,} reads  "
                        f"({_rate_per_min / 1000:.0f}k reads/min, {_elapsed / 60:.1f} min elapsed)"
                    )
                    _next_progress += _progress_interval

            # Write remaining
            if chunk:
                _write_results_chunk(out_fh, chunk)

        except Exception:
            _failed = True
            bam.close()
            raise
        else:
            bam.close()

    finally:
        out_fh.close()
        if netseq_loader:
            netseq_loader.close()
        # Remove partial output on failure so callers don't see a truncated file
        if _failed:
            _partial = Path(output_path)
            if _partial.exists():
                _partial.unlink()
                logger.warning(f"Removed partial output file after error: {output_path}")

    logger.info(f"Completed processing {stats.reads_processed:,} reads")
    logger.info(f"  Output written to {output_path}")

    # Write compact position index
    _base = str(output_path)
    if _base.endswith('.tsv.gz'):
        _index_path = _base[:-len('.tsv.gz')] + '_index.bed.gz'
    elif _base.endswith('.tsv'):
        _index_path = _base[:-len('.tsv')] + '_index.bed.gz'
    else:
        _index_path = _base + '_index.bed.gz'

    with gzip.open(_index_path, 'wt') as _idx_fh:
        _idx_fh.write('chrom\tcorrected_3prime\tstrand\tcount\n')
        for (chrom, pos, strand), count in sorted(_pos_counts.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
            _idx_fh.write(f'{chrom}\t{pos}\t{strand}\t{count:.6f}\n')
    logger.info(f"Position index written to {_index_path}")

    return stats


def process_bam_streaming_parallel(
    bam_path: str,
    genome_path: str,
    output_path: str,
    n_threads: int = 0,
    apply_atract: bool = True,
    apply_ag_mispriming: bool = False,
    ag_threshold: float = 0.65,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_dir: Optional[str] = None,
    show_progress: bool = True,
    variant_aware: bool = False,
    variant_output_path: Optional[str] = None,
    annotated_junctions: Optional[set] = None,
    pool_chrom_index: Optional[Dict] = None,
    apply_3ss_rescue: bool = True,
    gene_interval_trees: Optional[Dict] = None,
    polya_model_path: Optional[str] = None,
    min_gap_size: int = 10000,
    checkpoint_dir: Optional[str] = None,
    keep_checkpoints: bool = False,
    variant_scan_cache: Optional[str] = None,
    max_reads_for_variant_rescue: int = 500,
    dt_primed_cDNA: bool = False,
    ont_cDNA: bool = False,
    use_dorado_polya: bool = False,
    min_mapq: int = 0,
    min_aligned_length: int = 0,
) -> ProcessingStats:
    """
    Process BAM file with parallel region workers and streaming output.

    Combines the memory efficiency of process_bam_streaming (results written
    immediately, not accumulated in RAM) with the throughput of
    process_bam_file_parallel (multiple worker processes).

    Regions are computed from genomic coverage gaps (same as
    process_bam_file_parallel). Each worker processes one region and returns
    its results; the main thread writes each batch to disk as it arrives.

    Args:
        bam_path: Input BAM path
        genome_path: Genome FASTA path
        output_path: Output TSV path (supports .gz compression)
        n_threads: Worker processes (0 = auto-detect from SLURM/system)
        apply_*: Correction module flags
        netseq_dir: Optional NET-seq BigWig directory
        show_progress: Log reads/min progress
        variant_aware: Enable variant-aware homopolymer rescue (two-pass)
        variant_output_path: Optional path to write potential variants TSV
        checkpoint_dir: If set, write resumable per-region checkpoints here
                        (region_*.tsv/.stats.json/.done + rescue_scan.pkl).
        keep_checkpoints: If False (default), delete the dead per-region
                          checkpoint files (+ this chunk's rescue_scan.pkl) after
                          the final output TSV is rebuilt successfully, leaving
                          ~0 checkpoint files. If True, retain them for debugging.
                          Only affects the success path; on failure the committed
                          region checkpoints are always preserved for resume.

    Returns:
        ProcessingStats object
    """
    import time as _time

    if n_threads <= 0:
        n_threads = get_available_cpus()

    logger.info(f"Using {n_threads} worker(s) for parallel streaming processing")

    # Load genome for variant-aware prescan. Spawned workers load clean
    # per-process copies for region correction.
    logger.info(f"Loading genome from {genome_path}...")
    genome = load_genome(genome_path)

    # Load poly(A) model for parent-side setup/logging. Spawned workers load
    # clean per-process copies for region correction.
    polya_model = load_polya_model(Path(polya_model_path) if polya_model_path else None)
    if polya_model is not None:
        logger.info(f"Loaded poly(A) model from {polya_model_path}")

    # Checkpoint directory setup
    import pickle as _pickle
    _chk_dir = Path(checkpoint_dir) if checkpoint_dir else None
    if _chk_dir:
        _chk_dir.mkdir(parents=True, exist_ok=True)
    _rescue_pkl = (_chk_dir / 'rescue_scan.pkl') if _chk_dir else None

    # Variant-aware first pass (single-threaded pre-scan).
    # Priority: (1) variant_scan_cache from `rectify prescan` (chunked pipelines),
    #           (2) checkpoint pkl from a prior streaming run (resume),
    #           (3) run the scan fresh.
    variant_aware_rescue = None
    if variant_aware:
        if variant_scan_cache and Path(variant_scan_cache).exists():
            logger.info(f"Loading pre-computed variant scan from cache: {variant_scan_cache}")
            with open(variant_scan_cache, 'rb') as _f:
                variant_aware_rescue = _pickle.load(_f)
        elif _rescue_pkl and _rescue_pkl.exists():
            logger.info(f"Loading pre-computed variant scan from checkpoint: {_rescue_pkl}")
            with open(_rescue_pkl, 'rb') as _f:
                variant_aware_rescue = _pickle.load(_f)

    stats, regions, prescan_rescue = _prescan_bam_once(
        bam_path,
        genome,
        min_gap_size=min_gap_size,
        variant_aware=variant_aware and variant_aware_rescue is None,
        output_variants_path=variant_output_path,
    )
    if prescan_rescue is not None:
        variant_aware_rescue = prescan_rescue
        if _rescue_pkl:
            with open(_rescue_pkl, 'wb') as _f:
                _pickle.dump(variant_aware_rescue, _f, protocol=_pickle.HIGHEST_PROTOCOL)
            logger.info(f"Scan checkpoint saved: {_rescue_pkl}")
    logger.info(f"  {len(regions)} regions across {n_threads} workers")
    _pos_counts: dict = {}
    _t_start = _time.monotonic()
    _next_progress = 100000

    # Checkpoint: find already-completed regions
    _done_region_idxs: set = set()
    if _chk_dir:
        for _sf in _chk_dir.glob('region_*.done'):
            try:
                _idx = int(_sf.stem.split('_')[1])
            except (ValueError, IndexError):
                continue
            if _is_region_checkpoint_complete(_chk_dir, _idx):
                _done_region_idxs.add(_idx)
            else:
                logger.warning(
                    "Ignoring incomplete checkpoint marker for region %04d "
                    "(done sentinel without committed TSV)",
                    _idx,
                )
        if _done_region_idxs:
            logger.info(
                f"Checkpoint: {len(_done_region_idxs)}/{len(regions)} regions already done, resuming..."
            )

    _regions_to_run = [r for i, r in enumerate(regions) if i not in _done_region_idxs]
    _orig_idxs = [i for i, _r in enumerate(regions) if i not in _done_region_idxs]

    # Use a shared temp directory so workers can offload large results to disk
    # instead of sending them through the multiprocessing pipe.  This prevents
    # the pipe-buffer-full deadlock that occurs on high-coverage regions (e.g.
    # rDNA): a pipe write blocks when the OS buffer is full, and since the main
    # process is also blocked waiting for the complete pickle, neither side can
    # make progress.  Workers write their results to a temp .pkl file and return
    # only the path string; the main process loads and deletes each file.
    import tempfile as _tempfile
    _tmp_dir = _tempfile.mkdtemp(prefix='rectify_region_')

    # Commit C-II.0 worker state split (streaming path mirrors the in-memory
    # process_bam_file_parallel layout above): per-aligner fields go into
    # task tuples; everything else stays at init level.
    shared_kwargs = dict(
        apply_atract=apply_atract,
        apply_ag_mispriming=apply_ag_mispriming,
        ag_threshold=ag_threshold,
        apply_polya_trim=apply_polya_trim,
        apply_indel_correction=apply_indel_correction,
        netseq_dir=netseq_dir,
        annotated_junctions=annotated_junctions,
        pool_chrom_index=pool_chrom_index,
        apply_3ss_rescue=apply_3ss_rescue,
        gene_interval_trees=gene_interval_trees,
        tmp_dir=_tmp_dir,
        max_reads_for_variant_rescue=max_reads_for_variant_rescue,
        dt_primed_cDNA=dt_primed_cDNA,
        ont_cDNA=ont_cDNA,
        use_dorado_polya=use_dorado_polya,
        min_mapq=min_mapq,
        min_aligned_length=min_aligned_length,
    )

    # Do not keep parent-loaded genome/model objects alive across worker
    # creation. Spawned workers load clean copies in their initializer.
    try:
        del genome
        del polya_model
    except NameError:
        pass
    import gc as _gc
    _gc.collect()

    # Build the TSV header (identical to process_bam_streaming)
    _header = CORRECTION_TSV_HEADER

    _failed = False
    out_fh = None
    try:
        if not _chk_dir:
            if output_path.endswith('.gz'):
                out_fh = gzip.open(output_path, 'wt')
            else:
                out_fh = open(output_path, 'w')
            out_fh.write('\t'.join(_header) + '\n')

        _map_fn = 'imap' if _chk_dir else 'imap_unordered'
        ctx = _get_bam_worker_context()
        logger.info(
            "Using multiprocessing start method '%s' for BAM region workers",
            ctx.get_start_method(),
        )
        if _regions_to_run:
            _region_tasks = [
                (region, bam_path, variant_aware_rescue) for region in _regions_to_run
            ]
            # Pre-fork JIT warmup — see the sibling call above.
            _warmup_hp_ed_numba()
            with ctx.Pool(
                n_threads,
                initializer=_init_region_worker_state,
                initargs=(str(genome_path), polya_model_path, shared_kwargs),
            ) as pool:
                _iter = getattr(pool, _map_fn)(_process_region_worker_from_state, _region_tasks)
                for _batch_num, _worker_ret in enumerate(_iter):
                    # Workers return a temp-file path (str) when tmp_dir is set.
                    # Load the pickle and delete the file immediately to free space.
                    if isinstance(_worker_ret, str):
                        with open(_worker_ret, 'rb') as _pkl_fh:
                            region_results = _pickle.load(_pkl_fh)
                        os.unlink(_worker_ret)
                    else:
                        region_results = _worker_ret

                    if _chk_dir:
                        orig_idx = _orig_idxs[_batch_num]
                        _write_region_results_atomic(_chk_dir, orig_idx, region_results)
                    else:
                        _write_results_chunk(out_fh, region_results)

                    for result in region_results:
                        stats.update_from_result(result)
                        _pos_counts.setdefault(
                            (result['chrom'], result['corrected_3prime'], result['strand']), 0.0
                        )
                        _pos_counts[(result['chrom'], result['corrected_3prime'], result['strand'])] += float(result.get('fraction', 1.0))

                    if show_progress and stats.reads_processed >= _next_progress:
                        _elapsed = _time.monotonic() - _t_start
                        _rate = (stats.reads_processed / _elapsed * 60) if _elapsed > 0 else 0
                        logger.info(
                            f"  Processed {stats.reads_processed:,} reads  "
                            f"({_rate / 1000:.0f}k reads/min, {_elapsed / 60:.1f} min elapsed)"
                        )
                        _next_progress += 100000
        elif _chk_dir:
            logger.info("Checkpoint: all regions already complete; rebuilding final output")

    except Exception:
        _failed = True
        raise
    finally:
        if out_fh is not None:
            out_fh.close()
        # Clean up temp dir (any leftover .pkl files from failed workers)
        import shutil as _shutil
        try:
            _shutil.rmtree(_tmp_dir, ignore_errors=True)
        except Exception:
            pass
        if _failed and not _chk_dir:
            _partial = Path(output_path)
            if _partial.exists():
                _partial.unlink()
                logger.warning(f"Removed partial output after error: {output_path}")
        elif _failed and _chk_dir:
            logger.warning(
                f"Processing failed — committed region checkpoints preserved in {_chk_dir}"
            )

    if _chk_dir:
        stats, _pos_counts = _rebuild_output_from_region_files(
            _chk_dir,
            output_path,
            len(regions),
            stats,
        )

        # ── Success-path checkpoint cleanup ─────────────────────────────────
        # The final output TSV has now been rebuilt + atomically committed from
        # the per-region checkpoints above, so region_NNNN.{tsv,stats.json,done}
        # and this chunk's variant-scan pickle (rescue_scan.pkl) are DEAD resume
        # state. On a large cluster panel these per-region files dominate the
        # .checkpoints inode count (the driver of the Oak 15M-inode quota blow),
        # so remove them by default. keep_checkpoints=True (--keep-checkpoints)
        # retains them for debugging.
        #
        # Runs ONLY here, after a successful rebuild. The except/finally block
        # above re-raises on failure and explicitly preserves committed region
        # checkpoints ("committed region checkpoints preserved in ..."), so a
        # requeued task still resumes from completed regions.
        if not keep_checkpoints:
            _region_cleanup: List[Path] = []
            for _ridx in range(len(regions)):
                _region_cleanup.append(_region_tsv_path(_chk_dir, _ridx))
                _region_cleanup.append(_region_stats_path(_chk_dir, _ridx))
                _region_cleanup.append(_region_done_path(_chk_dir, _ridx))
            if _rescue_pkl is not None:
                _region_cleanup.append(_rescue_pkl)
            for _p in _region_cleanup:
                try:
                    _p.unlink()
                except OSError:
                    pass

    logger.info(f"Completed processing {stats.reads_processed:,} reads")
    logger.info(f"  Output written to {output_path}")

    # Write compact position index
    _base = str(output_path)
    if _base.endswith('.tsv.gz'):
        _index_path = _base[:-len('.tsv.gz')] + '_index.bed.gz'
    elif _base.endswith('.tsv'):
        _index_path = _base[:-len('.tsv')] + '_index.bed.gz'
    else:
        _index_path = _base + '_index.bed.gz'

    with gzip.open(_index_path, 'wt') as _idx_fh:
        _idx_fh.write('chrom\tcorrected_3prime\tstrand\tcount\n')
        for (chrom, pos, strand), count in sorted(_pos_counts.items(), key=lambda x: (x[0][0], x[0][1])):
            _idx_fh.write(f'{chrom}\t{pos}\t{strand}\t{count:.6f}\n')
    logger.info(f"Position index written to {_index_path}")

    return stats
