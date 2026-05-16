#!/usr/bin/env python3
"""
Parallel / streaming BAM processing dispatch for RECTIFY.

Wraps the per-read correction core (``bam_processor.correct_read_3prime``)
with region-based parallel workers, streaming TSV writers, and a streaming +
parallel combination for very large BAMs.

Workers fork from the parent process so the loaded genome, poly(A) model, and
walkback scratch arrays in ``rectify.core.correct.walkback`` are shared via
copy-on-write. The per-read function itself is imported from
``bam_processor`` rather than re-defined here so that all workers reach the
same code object and walkback state.

Author: Kevin R. Roy
"""

from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import gzip
import logging
import os

import pysam

from ..correct.indel_corrector import VariantAwareHomopolymerRescue
from ..polya.polya_model import PolyAModel, load_model as load_polya_model
from ...slurm import get_available_cpus
from ...utils.genome import load_genome
from .bam_processor import correct_read_3prime, _load_netseq
from .output import write_output_tsv
from .processing_stats import ProcessingStats
from .regions import get_processing_regions
from .variant_scan import run_variant_aware_scan
from ..position_index import write_position_index

logger = logging.getLogger(__name__)


def _write_results_chunk(fh, results: List[Dict]):
    """Write a chunk of results to file handle."""
    for result in results:
        _pt = result.get('pt_tag')
        _ps = result.get('polya_score')
        row = [
            result['read_id'],
            result['chrom'],
            result['strand'],
            str(result['original_3prime']),
            str(result['corrected_3prime']),
            str(result.get('five_prime_position', '')),  # 5' end (TSS)
            '1' if result.get('five_prime_rescued') else '0',  # 5' rescue flag
            result.get('five_prime_exon_cigar') or '',  # exon CIGAR for Cat3
            str(result.get('alignment_start', '')),  # Read body start
            str(result.get('alignment_end', '')),  # Read body end (exclusive)
            str(result['ambiguity_min']),
            str(result['ambiguity_max']),
            str(result['ambiguity_range']),
            str(result.get('polya_length', 0)),  # poly(A) length, default 0 if not computed
            str(result.get('aligned_a_length', 0)),  # Aligned A's
            str(result.get('soft_clip_a_length', 0)),  # Soft-clipped A's
            result.get('junctions_str', ''),  # Junctions as semicolon-separated string
            str(result.get('n_junctions', 0)),  # Number of junctions
            str(result.get('five_prime_soft_clip_length', 0)),  # 5' soft clip
            str(result.get('three_prime_soft_clip_length', 0)),  # 3' soft clip
            str(result.get('mapq', 0)),  # Mapping quality
            ','.join(result['correction_applied']) if result['correction_applied'] else 'none',
            result['confidence'],
            ','.join(result['qc_flags']),
            f"{result.get('fraction', 1.0):.6f}",
            result.get('gene_id') or '',  # Per-read gene attribution (empty if not computed)
            str(_pt) if _pt is not None else '',
            f'{_ps:.4f}' if _ps is not None else '',
            result.get('polya_source', 'none'),
            str(result.get('sc_homopolymer_extension', 0)),  # Cat2 CIGAR surgery
            result.get('sc_rescued_seq', ''),
            str(result.get('sc_original_softclip_len', 0)),
            str(result.get('five_prime_intron_clip_pos', -1)),  # Case 4 BAM clip
        ]
        fh.write('\t'.join(row) + '\n')


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
    gene_interval_trees: Optional[Dict] = None,
    polya_model: Optional[PolyAModel] = None,
    tmp_dir: Optional[str] = None,
    max_reads_for_variant_rescue: int = 500,
    dt_primed_cDNA: bool = False,
    min_mapq: int = 0,
    min_aligned_length: int = 0,
) -> Union[List[Dict], str]:
    """
    Worker function to process a single region.

    Opens its own BAM handle to avoid multiprocessing issues.

    Args:
        region: Tuple of (chrom, start, end)
        bam_path: Path to BAM file
        genome: Pre-loaded genome dict (shared via fork)
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
                gene_interval_trees=gene_interval_trees,
                polya_model=polya_model,
                dt_primed_cDNA=dt_primed_cDNA,
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
    gene_interval_trees: Optional[Dict] = None,
    polya_model_path: Optional[str] = None,
    variant_scan_cache: Optional[str] = None,
    max_reads_for_variant_rescue: int = 500,
    dt_primed_cDNA: bool = False,
    min_mapq: int = 0,
    min_aligned_length: int = 0,
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

    # Load genome (shared across workers via fork)
    logger.info(f"Loading genome from {genome_path}...")
    genome = load_genome(genome_path)

    # Load poly(A) model if provided (once, shared across all workers via fork)
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
        else:
            variant_aware_rescue = run_variant_aware_scan(
                bam_path=bam_path,
                genome=genome,
                min_variant_fraction=0.8,
                min_reads_for_variant_call=5,
                output_variants_path=variant_output_path,
            )

    # Count total reads and filtered reads for stats
    stats = ProcessingStats()
    logger.info("Counting reads in BAM file...")
    bam = pysam.AlignmentFile(bam_path, 'rb')
    for read in bam:
        stats.total_reads_in_bam += 1
        if read.is_unmapped:
            stats.reads_unmapped += 1
        elif read.is_secondary:
            stats.reads_secondary += 1
        elif read.is_supplementary:
            stats.reads_supplementary += 1
    bam.close()
    logger.info(f"  Total reads: {stats.total_reads_in_bam:,}")
    logger.info(f"  Unmapped: {stats.reads_unmapped:,}")
    logger.info(f"  Secondary: {stats.reads_secondary:,}")
    logger.info(f"  Supplementary: {stats.reads_supplementary:,}")

    # Get processing regions
    logger.info("Identifying processing regions...")
    regions = get_processing_regions(bam_path, min_gap_size=min_gap_size)
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
                gene_interval_trees,
                polya_model,
                max_reads_for_variant_rescue=max_reads_for_variant_rescue,
                dt_primed_cDNA=dt_primed_cDNA,
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

    # Multi-threaded processing
    logger.info(f"Processing {len(regions)} regions across {n_threads} workers...")

    # Create partial function with fixed arguments
    worker_func = partial(
        _process_region_worker,
        bam_path=bam_path,
        genome=genome,
        apply_atract=apply_atract,
        apply_ag_mispriming=apply_ag_mispriming,
        ag_threshold=ag_threshold,
        apply_polya_trim=apply_polya_trim,
        apply_indel_correction=apply_indel_correction,
        netseq_dir=netseq_dir,
        variant_aware_rescue=variant_aware_rescue,
        annotated_junctions=annotated_junctions,
        pool_chrom_index=pool_chrom_index,
        gene_interval_trees=gene_interval_trees,
        polya_model=polya_model,
        max_reads_for_variant_rescue=max_reads_for_variant_rescue,
        dt_primed_cDNA=dt_primed_cDNA,
        min_mapq=min_mapq,
        min_aligned_length=min_aligned_length,
    )

    all_results = []

    # Process with pool
    with Pool(n_threads) as pool:
        try:
            # Try to use tqdm for progress if available
            from tqdm import tqdm
            results_iter = tqdm(
                pool.imap(worker_func, regions),
                total=len(regions),
                desc="Processing regions"
            )
        except ImportError:
            results_iter = pool.imap(worker_func, regions)

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
    gene_interval_trees: Optional[Dict] = None,
    polya_model_path: Optional[str] = None,
    dt_primed_cDNA: bool = False,
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
        header = [
            'read_id', 'chrom', 'strand',
            'original_3prime', 'corrected_3prime',
            'five_prime_position',  # TSS end of the read
            'five_prime_rescued',   # 1 if 5' end was corrected by junction rescue (v2.7.9)
            'five_prime_exon_cigar',  # SAM CIGAR for exon segment of Cat3 rescue (v2.8.0)
            'alignment_start', 'alignment_end',  # Full read body interval
            'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
            'polya_length',  # Total observed poly(A) tail length
            'aligned_a_length', 'soft_clip_a_length',  # Breakdown of poly(A)
            'junctions', 'n_junctions',  # Splice junctions
            'five_prime_soft_clip_length', 'three_prime_soft_clip_length',  # Soft clips
            'mapq',  # Mapping quality
            'correction_applied', 'confidence', 'qc_flags', 'fraction',
            'gene_id',  # Per-read gene attribution (optional)
            'pt_tag',      # dorado pt:i signal-level poly(A) length (blank if absent) (v2.9.0)
            'polya_score', # poly(A) model confidence 0-1 (blank if no model) (v2.9.0)
            'polya_source',  # 'pt_tag' | 'model' | 'none' (v2.9.0)
            'sc_homopolymer_extension',  # Cat2: under-called homopolymer bases (v2.9.1)
            'sc_rescued_seq',            # Cat2: non-poly-A bases matched to ref (v2.9.1)
            'sc_original_softclip_len',  # Cat2: original 3' soft-clip length (v2.9.1)
        ]
        out_fh.write('\t'.join(header) + '\n')

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
                    gene_interval_trees=gene_interval_trees,
                    polya_model=polya_model,
                    dt_primed_cDNA=dt_primed_cDNA,
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
    gene_interval_trees: Optional[Dict] = None,
    polya_model_path: Optional[str] = None,
    min_gap_size: int = 10000,
    checkpoint_dir: Optional[str] = None,
    variant_scan_cache: Optional[str] = None,
    max_reads_for_variant_rescue: int = 500,
    dt_primed_cDNA: bool = False,
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

    Returns:
        ProcessingStats object
    """
    import time as _time

    if n_threads <= 0:
        n_threads = get_available_cpus()

    logger.info(f"Using {n_threads} worker(s) for parallel streaming processing")

    # Load genome once — shared across workers via fork
    logger.info(f"Loading genome from {genome_path}...")
    genome = load_genome(genome_path)

    # Load poly(A) model once (shared via fork)
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
        else:
            variant_aware_rescue = run_variant_aware_scan(
                bam_path=bam_path,
                genome=genome,
                min_variant_fraction=0.8,
                min_reads_for_variant_call=5,
                output_variants_path=variant_output_path,
            )
            if _rescue_pkl:
                with open(_rescue_pkl, 'wb') as _f:
                    _pickle.dump(variant_aware_rescue, _f, protocol=_pickle.HIGHEST_PROTOCOL)
                logger.info(f"Scan checkpoint saved: {_rescue_pkl}")

    # Count total reads and filtered reads for stats (pre-scan, single-threaded)
    stats = ProcessingStats()
    logger.info("Counting reads in BAM file...")
    _bam_prescan = pysam.AlignmentFile(bam_path, 'rb')
    for _r in _bam_prescan:
        stats.total_reads_in_bam += 1
        if _r.is_unmapped:
            stats.reads_unmapped += 1
        elif _r.is_secondary:
            stats.reads_secondary += 1
        elif _r.is_supplementary:
            stats.reads_supplementary += 1
    _bam_prescan.close()
    logger.info(f"  Total reads: {stats.total_reads_in_bam:,}")

    # Compute genomic regions
    logger.info("Identifying processing regions...")
    regions = get_processing_regions(bam_path, min_gap_size=min_gap_size)
    logger.info(f"  {len(regions)} regions across {n_threads} workers")
    _pos_counts: dict = {}
    _t_start = _time.monotonic()
    _next_progress = 100000

    # Checkpoint: find already-completed regions
    _done_region_idxs: set = set()
    if _chk_dir:
        for _sf in _chk_dir.glob('region_*.done'):
            try:
                _done_region_idxs.add(int(_sf.stem.split('_')[1]))
            except (ValueError, IndexError):
                pass
        if _done_region_idxs:
            logger.info(
                f"Checkpoint: {len(_done_region_idxs)}/{len(regions)} regions already done, resuming..."
            )
            _n_partial = _rebuild_pos_counts_from_partial(output_path, _pos_counts)
            logger.info(f"  Rebuilt pos_counts from {_n_partial:,} rows in partial output")

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

    worker_func = partial(
        _process_region_worker,
        bam_path=bam_path,
        genome=genome,
        apply_atract=apply_atract,
        apply_ag_mispriming=apply_ag_mispriming,
        ag_threshold=ag_threshold,
        apply_polya_trim=apply_polya_trim,
        apply_indel_correction=apply_indel_correction,
        netseq_dir=netseq_dir,
        variant_aware_rescue=variant_aware_rescue,
        annotated_junctions=annotated_junctions,
        pool_chrom_index=pool_chrom_index,
        gene_interval_trees=gene_interval_trees,
        polya_model=polya_model,
        tmp_dir=_tmp_dir,
        max_reads_for_variant_rescue=max_reads_for_variant_rescue,
        dt_primed_cDNA=dt_primed_cDNA,
        min_mapq=min_mapq,
        min_aligned_length=min_aligned_length,
    )

    # Build the TSV header (identical to process_bam_streaming)
    _header = [
        'read_id', 'chrom', 'strand',
        'original_3prime', 'corrected_3prime',
        'five_prime_position', 'five_prime_rescued', 'five_prime_exon_cigar',
        'alignment_start', 'alignment_end',
        'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
        'polya_length', 'aligned_a_length', 'soft_clip_a_length',
        'junctions', 'n_junctions',
        'five_prime_soft_clip_length', 'three_prime_soft_clip_length',
        'mapq', 'correction_applied', 'confidence', 'qc_flags', 'fraction',
        'gene_id', 'pt_tag', 'polya_score', 'polya_source',
        'sc_homopolymer_extension', 'sc_rescued_seq', 'sc_original_softclip_len',
        'five_prime_intron_clip_pos',
    ]

    _resuming = bool(_done_region_idxs) and Path(output_path).exists()
    if output_path.endswith('.gz'):
        out_fh = gzip.open(output_path, 'at' if _resuming else 'wt')
    else:
        out_fh = open(output_path, 'a' if _resuming else 'w')

    _failed = False
    try:
        if not _resuming:
            out_fh.write('\t'.join(_header) + '\n')

        _map_fn = 'imap' if _chk_dir else 'imap_unordered'
        with Pool(n_threads) as pool:
            _iter = getattr(pool, _map_fn)(worker_func, _regions_to_run)
            for _batch_num, _worker_ret in enumerate(_iter):
                # Workers return a temp-file path (str) when tmp_dir is set.
                # Load the pickle and delete the file immediately to free space.
                if isinstance(_worker_ret, str):
                    with open(_worker_ret, 'rb') as _pkl_fh:
                        region_results = _pickle.load(_pkl_fh)
                    os.unlink(_worker_ret)
                else:
                    region_results = _worker_ret
                _write_results_chunk(out_fh, region_results)

                for result in region_results:
                    stats.update_from_result(result)
                    _pos_counts.setdefault(
                        (result['chrom'], result['corrected_3prime'], result['strand']), 0.0
                    )
                    _pos_counts[(result['chrom'], result['corrected_3prime'], result['strand'])] += float(result.get('fraction', 1.0))

                if _chk_dir:
                    orig_idx = _orig_idxs[_batch_num]
                    (_chk_dir / f'region_{orig_idx:04d}.done').touch()

                if show_progress and stats.reads_processed >= _next_progress:
                    _elapsed = _time.monotonic() - _t_start
                    _rate = (stats.reads_processed / _elapsed * 60) if _elapsed > 0 else 0
                    logger.info(
                        f"  Processed {stats.reads_processed:,} reads  "
                        f"({_rate / 1000:.0f}k reads/min, {_elapsed / 60:.1f} min elapsed)"
                    )
                    _next_progress += 100000

    except Exception:
        _failed = True
        raise
    finally:
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
                f"Processing failed — partial output preserved for checkpoint resume: {output_path}"
            )

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
