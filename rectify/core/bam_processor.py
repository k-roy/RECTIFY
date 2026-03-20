#!/usr/bin/env python3
"""
BAM processing pipeline for RECTIFY.

This module orchestrates the full correction workflow:
1. A-tract ambiguity detection (universal)
2. AG mispriming screening (for oligo-dT)
3. Poly(A) tail trimming (when poly(A) sequenced)
4. Indel artifact correction (when poly(A) sequenced)
5. NET-seq refinement (optional)

Supports parallel processing via:
- Per-chromosome parallelization for small genomes (e.g., yeast)
- Region-based parallelization with gap splitting for large chromosomes

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import Dict, List, Optional, Tuple, Union
from multiprocessing import Pool
from functools import partial
import gzip
import logging
import pysam
from pathlib import Path

from . import atract_detector
from . import ag_mispriming
from . import polya_trimmer
from . import indel_corrector
from . import netseq_refiner
from .processing_stats import ProcessingStats, write_stats_tsv, generate_stats_report
from ..utils.genome import load_genome, standardize_chrom_name
from ..slurm import get_available_cpus

logger = logging.getLogger(__name__)


def get_read_3prime_position(read: pysam.AlignedSegment) -> Tuple[int, str]:
    """
    Get 3' end position and strand from read.

    Args:
        read: pysam AlignedSegment

    Returns:
        Tuple of (position, strand)
    """
    # Determine strand
    strand = '-' if read.is_reverse else '+'

    # Get 3' end position
    if strand == '+':
        position = read.reference_end - 1  # 0-based inclusive
    else:
        position = read.reference_start  # 0-based

    return position, strand


def correct_read_3prime(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    apply_atract: bool = True,
    apply_ag_mispriming: bool = False,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_loader: Optional[netseq_refiner.NetseqLoader] = None
) -> Dict:
    """
    Apply all corrections to a single read.

    Args:
        read: pysam AlignedSegment
        genome: Genome dict from load_genome()
        apply_atract: Apply A-tract ambiguity detection
        apply_ag_mispriming: Apply AG mispriming screening
        apply_polya_trim: Apply poly(A) tail trimming
        apply_indel_correction: Apply indel artifact correction
        netseq_loader: Optional NET-seq loader for refinement

    Returns:
        Dict with correction results
    """
    # Get original position and strand
    original_position, strand = get_read_3prime_position(read)
    chrom = read.reference_name

    # Standardize chromosome name
    chrom_std = standardize_chrom_name(chrom)

    # Initialize result (use standardized chromosome name for output)
    result = {
        'read_id': read.query_name,
        'chrom': chrom_std,  # Use Roman numeral format (chrI, chrII, etc.)
        'strand': strand,
        'original_3prime': original_position,
        'corrected_3prime': original_position,
        'ambiguity_min': original_position,
        'ambiguity_max': original_position,
        'ambiguity_range': 0,
        'correction_applied': [],
        'qc_flags': [],
        'confidence': 'high',
    }

    current_position = original_position

    # Module 1: A-tract ambiguity (always applied by default)
    if apply_atract:
        atract_result = atract_detector.calculate_atract_ambiguity(
            genome, chrom_std, current_position, strand, downstream_bp=10
        )

        result['ambiguity_min'] = atract_result['ambiguity_min']
        result['ambiguity_max'] = atract_result['ambiguity_max']
        result['ambiguity_range'] = atract_result['ambiguity_range']
        result['downstream_a_count'] = atract_result.get('downstream_a_count')

        if atract_result['has_ambiguity']:
            result['correction_applied'].append('atract_ambiguity')

    # Module 2A: AG mispriming (for oligo-dT technologies)
    if apply_ag_mispriming:
        ag_result = ag_mispriming.screen_ag_mispriming(
            genome, chrom_std, current_position, strand
        )

        if ag_result['is_likely_misprimed']:
            qc_flag = ag_mispriming.get_ag_qc_flag(ag_result)
            result['qc_flags'].append(qc_flag)
            result['ag_content'] = ag_result['ag_content']
            result['ag_confidence'] = ag_result['confidence']

            # Reduce confidence for misprimed reads
            if result['confidence'] == 'high':
                result['confidence'] = 'medium' if ag_result['confidence'] == 'medium' else 'low'

    # Module 2B: Poly(A) tail detection (when poly(A) is sequenced)
    # NOTE: This does NOT correct positions - it only measures poly(A) tail length.
    # Position correction for A-tract ambiguity is already handled above.
    polya_shift = 0
    if apply_polya_trim:
        # Pass atract_result so polya_trimmer can include aligned A's in the count
        atract_for_polya = {
            'tract_length': result.get('ambiguity_range', 0),  # Approximate aligned A's
        }
        polya_result = polya_trimmer.trim_polya_from_read(
            read, strand, atract_result=atract_for_polya
        )

        # Always record poly(A) length (even if 0) for completeness
        result['polya_length'] = polya_result['polya_length']
        result['aligned_a_length'] = polya_result['aligned_a_length']
        result['soft_clip_a_length'] = polya_result['soft_clip_a_length']

        # NOTE: We no longer add 'polya_trim' to correction_applied because
        # soft-clips don't affect genomic position. The only position correction
        # comes from A-tract ambiguity detection (atract_ambiguity).

    # Module 2C: Indel artifact correction (when poly(A) is sequenced)
    indel_shift = 0
    if apply_indel_correction:
        indel_result = indel_corrector.correct_indels_from_read(read, strand, genome)

        if indel_result['has_artifacts']:
            result['correction_applied'].append('indel_correction')
            result['n_indel_artifacts'] = len(indel_result['artifacts'])
            indel_shift = indel_result['correction_bp']

            # Update current position
            current_position = indel_result['corrected_3prime']

    # Update corrected position (after poly(A) and indel corrections)
    result['corrected_3prime'] = current_position

    # Update ambiguity window based on corrections
    if polya_shift != 0 or indel_shift != 0:
        if strand == '+':
            result['ambiguity_min'] = current_position - result['ambiguity_range']
            result['ambiguity_max'] = current_position
        else:
            result['ambiguity_min'] = current_position
            result['ambiguity_max'] = current_position + result['ambiguity_range']

    # Module 3: NET-seq refinement (optional)
    if netseq_loader is not None and result['ambiguity_range'] > 0:
        netseq_result = netseq_refiner.refine_with_netseq(
            netseq_loader,
            chrom_std,
            result['ambiguity_min'],
            result['ambiguity_max'],
            strand,
            current_position,
            proportional_split=False,  # Return single dict, not list
        )

        result['corrected_3prime'] = netseq_result['refined_position']
        result['netseq_confidence'] = netseq_result['confidence']
        result['netseq_method'] = netseq_result['method']
        result['netseq_peak_signal'] = netseq_result['peak_signal']

        # Update overall confidence based on NET-seq
        if netseq_result['confidence'] == 'low':
            result['confidence'] = 'low'
        elif netseq_result['confidence'] == 'medium' and result['confidence'] == 'high':
            result['confidence'] = 'medium'

        result['correction_applied'].append('netseq_refinement')

    # Set QC flags if none added
    if not result['qc_flags']:
        result['qc_flags'].append('PASS')

    return result


def process_bam_file(
    bam_path: str,
    genome_path: str,
    apply_atract: bool = True,
    apply_ag_mispriming: bool = False,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_dir: Optional[str] = None,
    output_path: Optional[str] = None,
    max_reads: Optional[int] = None
) -> List[Dict]:
    """
    Process BAM file and apply all corrections.

    Args:
        bam_path: Path to BAM file
        genome_path: Path to genome FASTA
        apply_atract: Apply A-tract ambiguity detection
        apply_ag_mispriming: Apply AG mispriming screening
        apply_polya_trim: Apply poly(A) tail trimming
        apply_indel_correction: Apply indel artifact correction
        netseq_dir: Optional directory with NET-seq BigWig files
        output_path: Optional output TSV path
        max_reads: Optional maximum number of reads to process

    Returns:
        List of correction result dicts
    """
    # Load genome
    print(f"Loading genome from {genome_path}...")
    genome = load_genome(genome_path)

    # Load NET-seq data if provided
    netseq_loader = None
    if netseq_dir:
        print(f"Loading NET-seq data from {netseq_dir}...")
        netseq_loader = netseq_refiner.NetseqLoader()
        netseq_loader.load_directory(netseq_dir, pattern="*.bw")
        print(f"  Loaded {len(netseq_loader.bigwigs)} BigWig file(s)")

    # Open BAM file
    print(f"Processing BAM file: {bam_path}")
    bam = pysam.AlignmentFile(bam_path, 'rb')

    results = []
    n_processed = 0

    for read in bam:
        # Skip unmapped reads
        if read.is_unmapped:
            continue

        # Skip secondary/supplementary alignments
        if read.is_secondary or read.is_supplementary:
            continue

        # Apply corrections
        result = correct_read_3prime(
            read,
            genome,
            apply_atract=apply_atract,
            apply_ag_mispriming=apply_ag_mispriming,
            apply_polya_trim=apply_polya_trim,
            apply_indel_correction=apply_indel_correction,
            netseq_loader=netseq_loader
        )

        results.append(result)
        n_processed += 1

        # Progress reporting
        if n_processed % 10000 == 0:
            print(f"  Processed {n_processed:,} reads...")

        # Limit for testing
        if max_reads and n_processed >= max_reads:
            print(f"  Reached max_reads limit ({max_reads})")
            break

    bam.close()

    print(f"Completed processing {n_processed:,} reads")

    # Write output if requested
    if output_path:
        write_output_tsv(results, output_path)

    return results


def write_output_tsv(results: List[Dict], output_path: str):
    """
    Write correction results to TSV file.

    Args:
        results: List of correction result dicts
        output_path: Path to output TSV file
    """
    print(f"Writing output to {output_path}...")

    with open(output_path, 'w') as f:
        # Write header
        header = [
            'read_id', 'chrom', 'strand',
            'original_3prime', 'corrected_3prime',
            'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
            'polya_length',  # Total observed poly(A) tail length
            'correction_applied', 'confidence', 'qc_flags'
        ]
        f.write('\t'.join(header) + '\n')

        # Write results
        for result in results:
            row = [
                result['read_id'],
                result['chrom'],
                result['strand'],
                str(result['original_3prime']),
                str(result['corrected_3prime']),
                str(result['ambiguity_min']),
                str(result['ambiguity_max']),
                str(result['ambiguity_range']),
                str(result.get('polya_length', 0)),  # poly(A) length, default 0 if not computed
                ','.join(result['correction_applied']) if result['correction_applied'] else 'none',
                result['confidence'],
                ','.join(result['qc_flags']),
            ]
            f.write('\t'.join(row) + '\n')

    print(f"  Wrote {len(results):,} corrected positions")


def generate_summary_report(results: List[Dict]) -> str:
    """
    Generate summary report from correction results.

    Uses single-pass accumulation for efficiency with large datasets.

    Args:
        results: List of correction result dicts

    Returns:
        Formatted report string
    """
    import numpy as np

    n_total = len(results)
    if n_total == 0:
        return "No reads processed."

    # Single-pass accumulation of all metrics
    n_with_ambiguity = 0
    n_corrected = 0
    ambiguity_ranges = []
    corrected_shifts = []
    by_confidence = {'high': 0, 'medium': 0, 'low': 0}

    for r in results:
        # Ambiguity check
        ambig_range = r['ambiguity_range']
        ambiguity_ranges.append(ambig_range)
        if ambig_range > 0:
            n_with_ambiguity += 1

        # Correction check
        shift = abs(r['corrected_3prime'] - r['original_3prime'])
        if shift > 0:
            n_corrected += 1
            corrected_shifts.append(shift)

        # Confidence
        conf = r['confidence']
        if conf in by_confidence:
            by_confidence[conf] += 1

    # Build report
    report = []
    report.append("=" * 60)
    report.append("RECTIFY Correction Summary")
    report.append("=" * 60)
    report.append("")
    report.append("Overall:")
    report.append(f"  Total reads:            {n_total:,}")
    report.append(f"  With ambiguity:         {n_with_ambiguity:,} ({100*n_with_ambiguity/n_total:.1f}%)")
    report.append(f"  Position corrected:     {n_corrected:,} ({100*n_corrected/n_total:.1f}%)")
    report.append("")

    report.append("Ambiguity Ranges:")
    report.append(f"  Mean:                   {np.mean(ambiguity_ranges):.2f} bp")
    report.append(f"  Median:                 {np.median(ambiguity_ranges):.1f} bp")
    report.append(f"  Maximum:                {max(ambiguity_ranges)} bp")
    report.append("")

    if n_corrected > 0:
        report.append("Position Corrections:")
        report.append(f"  Mean shift:             {np.mean(corrected_shifts):.2f} bp")
        report.append(f"  Median shift:           {np.median(corrected_shifts):.1f} bp")
        report.append(f"  Maximum shift:          {max(corrected_shifts)} bp")
        report.append("")

    report.append("Confidence Levels:")
    for conf in ['high', 'medium', 'low']:
        count = by_confidence.get(conf, 0)
        pct = 100.0 * count / n_total
        report.append(f"  {conf:10s}          {count:7,} ({pct:5.1f}%)")

    report.append("")
    report.append("=" * 60)

    return "\n".join(report)


# =============================================================================
# Parallel Processing Functions
# =============================================================================


def find_coverage_gaps(
    bam_path: str,
    chrom: str,
    min_gap_size: int = 10000
) -> List[Tuple[int, int]]:
    """
    Find gaps in read coverage that can be used to split chromosome.

    Large chromosomes can be split at coverage deserts where no reads exist,
    allowing independent parallel processing of regions.

    Args:
        bam_path: Path to BAM file
        chrom: Chromosome name
        min_gap_size: Minimum gap size in bp to consider a split point

    Returns:
        List of (start, end) tuples representing gap regions
    """
    bam = pysam.AlignmentFile(bam_path, 'rb')

    # Get chromosome length
    chrom_length = None
    for idx, ref in enumerate(bam.references):
        if ref == chrom:
            chrom_length = bam.lengths[idx]
            break

    if chrom_length is None:
        bam.close()
        return []

    # Use index to efficiently find gaps
    # Strategy: sample positions and find where no reads map
    gaps = []
    last_covered = 0

    # Iterate through all reads to find covered regions
    covered_positions = set()
    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        # Add start and end positions
        covered_positions.add(read.reference_start)
        covered_positions.add(read.reference_end)

    bam.close()

    if not covered_positions:
        return [(0, chrom_length)]

    # Sort covered positions
    sorted_positions = sorted(covered_positions)

    # Find gaps larger than threshold
    prev_pos = 0
    for pos in sorted_positions:
        if pos - prev_pos >= min_gap_size:
            gaps.append((prev_pos, pos))
        prev_pos = pos

    # Check final gap
    if chrom_length - prev_pos >= min_gap_size:
        gaps.append((prev_pos, chrom_length))

    return gaps


def get_processing_regions(
    bam_path: str,
    min_gap_size: int = 10000,
    max_region_size: int = 500000
) -> List[Tuple[str, int, int]]:
    """
    Get list of regions for parallel processing.

    Splits large chromosomes at coverage gaps to create balanced work units.

    Args:
        bam_path: Path to BAM file
        min_gap_size: Minimum gap to use as split point
        max_region_size: Target maximum region size

    Returns:
        List of (chrom, start, end) tuples representing processing regions
    """
    bam = pysam.AlignmentFile(bam_path, 'rb')
    chromosomes = list(bam.references)
    chrom_lengths = dict(zip(bam.references, bam.lengths))
    bam.close()

    regions = []

    for chrom in chromosomes:
        chrom_len = chrom_lengths[chrom]

        # For small chromosomes, process as single region
        if chrom_len <= max_region_size:
            regions.append((chrom, 0, chrom_len))
            continue

        # For large chromosomes, split at gaps
        gaps = find_coverage_gaps(bam_path, chrom, min_gap_size)

        if not gaps:
            # No gaps found, process as single region
            regions.append((chrom, 0, chrom_len))
            continue

        # Create regions between gaps
        current_start = 0
        for gap_start, gap_end in gaps:
            if gap_start > current_start:
                regions.append((chrom, current_start, gap_start))
            current_start = gap_end

        # Final region after last gap
        if current_start < chrom_len:
            regions.append((chrom, current_start, chrom_len))

    return regions


def _process_region_worker(
    region: Tuple[str, int, int],
    bam_path: str,
    genome: Dict[str, str],
    apply_atract: bool,
    apply_ag_mispriming: bool,
    apply_polya_trim: bool,
    apply_indel_correction: bool,
    netseq_dir: Optional[str] = None
) -> List[Dict]:
    """
    Worker function to process a single region.

    Opens its own BAM handle to avoid multiprocessing issues.

    Args:
        region: Tuple of (chrom, start, end)
        bam_path: Path to BAM file
        genome: Pre-loaded genome dict (shared via fork)
        apply_*: Correction module flags
        netseq_dir: Optional NET-seq BigWig directory

    Returns:
        List of correction result dicts for reads in region
    """
    chrom, start, end = region
    results = []

    # Open BAM for this region
    bam = pysam.AlignmentFile(bam_path, 'rb')

    # Load NET-seq if needed (per-worker to avoid file handle issues)
    netseq_loader = None
    if netseq_dir:
        netseq_loader = netseq_refiner.NetseqLoader()
        netseq_loader.load_directory(netseq_dir, pattern="*.bw")

    try:
        for read in bam.fetch(chrom, start, end):
            # Skip unmapped/secondary/supplementary
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            # Apply corrections
            result = correct_read_3prime(
                read,
                genome,
                apply_atract=apply_atract,
                apply_ag_mispriming=apply_ag_mispriming,
                apply_polya_trim=apply_polya_trim,
                apply_indel_correction=apply_indel_correction,
                netseq_loader=netseq_loader
            )
            results.append(result)

    finally:
        bam.close()
        if netseq_loader:
            netseq_loader.close()

    return results


def process_bam_file_parallel(
    bam_path: str,
    genome_path: str,
    n_threads: int = 0,
    apply_atract: bool = True,
    apply_ag_mispriming: bool = False,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_dir: Optional[str] = None,
    output_path: Optional[str] = None,
    max_reads: Optional[int] = None,
    min_gap_size: int = 10000,
    show_progress: bool = True,
    return_stats: bool = False
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
                apply_polya_trim, apply_indel_correction,
                netseq_dir
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
        apply_polya_trim=apply_polya_trim,
        apply_indel_correction=apply_indel_correction,
        netseq_dir=netseq_dir
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
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_dir: Optional[str] = None,
    show_progress: bool = True
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
        netseq_loader = netseq_refiner.NetseqLoader()
        netseq_loader.load_directory(netseq_dir, pattern="*.bw")

    # Initialize comprehensive stats
    stats = ProcessingStats()

    # Open output file (support gzip)
    if output_path.endswith('.gz'):
        out_fh = gzip.open(output_path, 'wt')
    else:
        out_fh = open(output_path, 'w')

    try:
        # Write header
        header = [
            'read_id', 'chrom', 'strand',
            'original_3prime', 'corrected_3prime',
            'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
            'polya_length',  # Total observed poly(A) tail length
            'correction_applied', 'confidence', 'qc_flags'
        ]
        out_fh.write('\t'.join(header) + '\n')

        # Open BAM
        bam = pysam.AlignmentFile(bam_path, 'rb')
        chunk = []

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

            result = correct_read_3prime(
                read, genome,
                apply_atract=apply_atract,
                apply_ag_mispriming=apply_ag_mispriming,
                apply_polya_trim=apply_polya_trim,
                apply_indel_correction=apply_indel_correction,
                netseq_loader=netseq_loader
            )
            chunk.append(result)

            # Update comprehensive stats
            stats.update_from_result(result)

            # Write chunk when full
            if len(chunk) >= chunk_size:
                _write_results_chunk(out_fh, chunk)
                chunk = []

                if show_progress and stats.reads_processed % 100000 == 0:
                    logger.info(f"  Processed {stats.reads_processed:,} reads...")

        # Write remaining
        if chunk:
            _write_results_chunk(out_fh, chunk)

        bam.close()

    finally:
        out_fh.close()
        if netseq_loader:
            netseq_loader.close()

    logger.info(f"Completed processing {stats.reads_processed:,} reads")
    logger.info(f"  Output written to {output_path}")

    return stats


def _write_results_chunk(fh, results: List[Dict]):
    """Write a chunk of results to file handle."""
    for result in results:
        row = [
            result['read_id'],
            result['chrom'],
            result['strand'],
            str(result['original_3prime']),
            str(result['corrected_3prime']),
            str(result['ambiguity_min']),
            str(result['ambiguity_max']),
            str(result['ambiguity_range']),
            str(result.get('polya_length', 0)),  # poly(A) length, default 0 if not computed
            ','.join(result['correction_applied']) if result['correction_applied'] else 'none',
            result['confidence'],
            ','.join(result['qc_flags']),
        ]
        fh.write('\t'.join(row) + '\n')


def generate_summary_from_stats(stats) -> str:
    """
    Generate summary report from ProcessingStats or legacy dict.

    Args:
        stats: ProcessingStats object or legacy Dict from process_bam_streaming

    Returns:
        Formatted report string
    """
    # Handle both ProcessingStats and legacy dict format
    if isinstance(stats, ProcessingStats):
        return generate_stats_report(stats)

    # Legacy dict format (backwards compatibility)
    n_total = stats.get('total_reads', 0)
    if n_total == 0:
        return "No reads processed."

    report = []
    report.append("=" * 60)
    report.append("RECTIFY Correction Summary (Streaming Mode)")
    report.append("=" * 60)
    report.append("")
    report.append("Overall:")
    report.append(f"  Total reads:            {n_total:,}")
    report.append(f"  With ambiguity:         {stats['with_ambiguity']:,} ({100*stats['with_ambiguity']/n_total:.1f}%)")
    report.append(f"  Position corrected:     {stats['position_corrected']:,} ({100*stats['position_corrected']/n_total:.1f}%)")
    report.append("")
    report.append("Confidence Levels:")
    for conf in ['high', 'medium', 'low']:
        count = stats['by_confidence'].get(conf, 0)
        pct = 100.0 * count / n_total
        report.append(f"  {conf:10s}          {count:7,} ({pct:5.1f}%)")
    report.append("")
    report.append("=" * 60)

    return "\n".join(report)
