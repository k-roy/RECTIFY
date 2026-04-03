"""
NET-seq signal deconvolution for RECTIFY.

Applies NNLS deconvolution to NET-seq signal to recover true CPA positions
from oligo(A)-spread signal.

The deconvolution uses the 0A PSF (Point-Spread-Function) empirically derived
from sites with no downstream genomic A's, where the "true" position is known.

This module operates on position-level aggregated signal, not per-read.
For read-level output, use proportional assignment from deconvolved peaks.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

from ..config import NETSEQ_DECONV_REGULARIZATION, CHROM_TO_GENOME
from .netseq_refiner import (
    DEFAULT_PSF_0A,
    build_convolution_matrix,
    deconvolve_signal,
    identify_deconvolved_peaks,
)


@dataclass
class DeconvolutionResult:
    """Result of deconvolving a genomic region."""
    chrom: str
    strand: str
    region_start: int
    region_end: int
    n_positions: int
    observed_signal: np.ndarray
    deconvolved_signal: np.ndarray
    peaks: List[Dict]  # Each dict has 'position', 'signal', 'fraction'
    fit_quality: float  # R-squared between observed and predicted


def count_downstream_a(
    genome: Dict[str, str],
    chrom: str,
    position: int,
    strand: str,
    window_size: int = 10,
) -> int:
    """
    Count A's (or T's for minus strand) downstream of position.

    For plus strand: Count A's to the RIGHT of position
    For minus strand: Count T's to the LEFT of position

    Args:
        genome: Dict mapping chrom -> sequence (may use NCBI or standard names)
        chrom: Chromosome name (standard format: chrI, chrII, etc.)
        position: 0-based position
        strand: '+' or '-'
        window_size: Window to check (default 10 bp)

    Returns:
        Number of A's (plus) or T's (minus) in downstream window
    """
    # Map chromosome name to genome format if needed
    genome_chrom = CHROM_TO_GENOME.get(chrom, chrom)

    # Try both standard and NCBI format
    if genome_chrom not in genome:
        if chrom not in genome:
            return 0
        genome_chrom = chrom

    seq = genome[genome_chrom].upper()
    seq_len = len(seq)

    if strand == '+':
        # Plus strand: downstream is to the right (higher positions)
        start = position + 1
        end = min(start + window_size, seq_len)
        if start >= seq_len:
            return 0
        window = seq[start:end]
        return window.count('A')
    else:
        # Minus strand: downstream is to the left (lower positions)
        end = position
        start = max(end - window_size, 0)
        if end <= 0:
            return 0
        window = seq[start:end]
        return window.count('T')


def identify_atract_regions(
    position_counts: Dict[Tuple[str, str, int], float],
    genome: Dict[str, str],
    min_downstream_a: int = 3,
    window_size: int = 10,
    min_region_signal: float = 10.0,
) -> List[Tuple[str, str, int, int, List[int]]]:
    """
    Identify genomic regions with A-tract ambiguity requiring deconvolution.

    Regions with >= min_downstream_a adenosines within window are
    candidates for signal spreading and need deconvolution.

    Args:
        position_counts: Position -> count mapping from aggregated reads
        genome: Reference genome dict
        min_downstream_a: Minimum A's to flag region
        window_size: Window size to check for A's
        min_region_signal: Minimum total signal in region

    Returns:
        List of (chrom, strand, start, end, positions) tuples for regions
        needing deconvolution
    """
    # Group positions by (chrom, strand)
    by_chrom_strand: Dict[Tuple[str, str], List[Tuple[int, float]]] = {}
    for (chrom, strand, pos), count in position_counts.items():
        key = (chrom, strand)
        if key not in by_chrom_strand:
            by_chrom_strand[key] = []
        by_chrom_strand[key].append((pos, count))

    regions = []

    for (chrom, strand), positions_counts in by_chrom_strand.items():
        # Sort by position
        positions_counts.sort(key=lambda x: x[0])

        # Find positions in A-tracts
        atract_positions = []
        for pos, count in positions_counts:
            n_downstream_a = count_downstream_a(genome, chrom, pos, strand, window_size)
            if n_downstream_a >= min_downstream_a:
                atract_positions.append((pos, count, n_downstream_a))

        if not atract_positions:
            continue

        # Cluster consecutive positions into regions
        current_region = [atract_positions[0]]
        for i in range(1, len(atract_positions)):
            pos, count, n_a = atract_positions[i]
            prev_pos = current_region[-1][0]

            # Merge if within 20bp of previous position
            if pos - prev_pos <= 20:
                current_region.append((pos, count, n_a))
            else:
                # Save current region if substantial
                total_signal = sum(c for _, c, _ in current_region)
                if total_signal >= min_region_signal:
                    positions_list = [p for p, _, _ in current_region]
                    start = min(positions_list)
                    end = max(positions_list) + 1
                    regions.append((chrom, strand, start, end, positions_list))
                current_region = [(pos, count, n_a)]

        # Don't forget last region
        if current_region:
            total_signal = sum(c for _, c, _ in current_region)
            if total_signal >= min_region_signal:
                positions_list = [p for p, _, _ in current_region]
                start = min(positions_list)
                end = max(positions_list) + 1
                regions.append((chrom, strand, start, end, positions_list))

    return regions


def deconvolve_region(
    position_counts: Dict[Tuple[str, str, int], float],
    chrom: str,
    strand: str,
    region_start: int,
    region_end: int,
    psf: np.ndarray = None,
    regularization: float = None,
) -> DeconvolutionResult:
    """
    Deconvolve signal within a single A-tract region.

    Args:
        position_counts: Full position count dict
        chrom: Chromosome
        strand: Strand
        region_start: Start of A-tract region (inclusive)
        region_end: End of A-tract region (exclusive)
        psf: Point-spread-function (default: DEFAULT_PSF_0A)
        regularization: NNLS regularization (default: from config)

    Returns:
        DeconvolutionResult with observed, deconvolved, and peaks
    """
    if psf is None:
        psf = DEFAULT_PSF_0A
    if regularization is None:
        regularization = NETSEQ_DECONV_REGULARIZATION

    # Build observed signal vector
    n_positions = region_end - region_start
    if n_positions <= 0:
        raise ValueError(
            f"deconvolve_region(): inverted coordinates — "
            f"region_start={region_start} >= region_end={region_end}"
        )
    observed = np.zeros(n_positions)
    positions = list(range(region_start, region_end))

    for i, pos in enumerate(positions):
        key = (chrom, strand, pos)
        if key in position_counts:
            observed[i] = position_counts[key]

    # Build convolution matrix
    # tract_length = n_positions - 1 (first position is "first non-A")
    A = build_convolution_matrix(n_positions - 1, psf, include_first_non_a=True)

    # Deconvolve
    deconvolved = deconvolve_signal(observed, A, regularization)

    # Calculate fit quality
    predicted = A @ deconvolved
    ss_res = np.sum((observed - predicted) ** 2)
    ss_tot = np.sum((observed - observed.mean()) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    # Identify peaks
    peak_dicts = identify_deconvolved_peaks(deconvolved, positions, threshold_frac=0.1)

    # Add fraction of total signal
    total_deconv = deconvolved.sum()
    for peak in peak_dicts:
        peak['fraction'] = peak['signal'] / total_deconv if total_deconv > 0 else 0.0

    return DeconvolutionResult(
        chrom=chrom,
        strand=strand,
        region_start=region_start,
        region_end=region_end,
        n_positions=n_positions,
        observed_signal=observed,
        deconvolved_signal=deconvolved,
        peaks=peak_dicts,
        fit_quality=r_squared,
    )


def deconvolve_all_regions(
    position_counts: Dict[Tuple[str, str, int], float],
    genome: Dict[str, str],
    min_downstream_a: int = 3,
    min_region_signal: float = 10.0,
    show_progress: bool = True,
) -> Tuple[Dict[Tuple[str, str, int], float], List[DeconvolutionResult]]:
    """
    Deconvolve all A-tract regions and return corrected position counts.

    For positions NOT in A-tracts, counts are passed through unchanged.
    For positions IN A-tracts, signal is redistributed to deconvolved peaks.

    Args:
        position_counts: Raw position counts
        genome: Reference genome (dict of chrom -> sequence)
        min_downstream_a: Minimum downstream A's to trigger deconvolution
        min_region_signal: Minimum total signal in region to deconvolve
        show_progress: Print progress messages

    Returns:
        Tuple of:
        - Deconvolved position counts (new dict)
        - List of DeconvolutionResult for each region
    """
    if show_progress:
        print("  Identifying A-tract regions...")

    # Identify regions needing deconvolution
    regions = identify_atract_regions(
        position_counts, genome,
        min_downstream_a=min_downstream_a,
        min_region_signal=min_region_signal,
    )

    if show_progress:
        print(f"  Found {len(regions)} regions with A-tract ambiguity")

    if not regions:
        # No regions to deconvolve - return copy of original
        return dict(position_counts), []

    # Track which positions are in deconvolved regions
    deconvolved_positions = set()
    for chrom, strand, start, end, positions in regions:
        for pos in positions:
            deconvolved_positions.add((chrom, strand, pos))

    # Start with non-deconvolved positions
    deconv_counts = {}
    for key, count in position_counts.items():
        if key not in deconvolved_positions:
            deconv_counts[key] = count

    # Deconvolve each region
    results = []
    for chrom, strand, start, end, positions in regions:
        result = deconvolve_region(
            position_counts, chrom, strand, start, end
        )
        results.append(result)

        # Add deconvolved peaks to output
        for peak in result.peaks:
            key = (chrom, strand, peak['position'])
            deconv_counts[key] = deconv_counts.get(key, 0) + peak['signal']

    if show_progress:
        total_peaks = sum(len(r.peaks) for r in results)
        avg_fit = np.mean([r.fit_quality for r in results]) if results else 0.0
        print(f"  Deconvolution complete: {total_peaks} peaks, mean R²={avg_fit:.3f}")

    return deconv_counts, results


def proportional_read_assignment(
    records: List,
    deconv_results: List[DeconvolutionResult],
    position_field: str = 'three_prime_raw',
) -> List:
    """
    Assign reads proportionally to deconvolved peak positions.

    For reads in A-tract regions, updates three_prime_corrected to the
    assigned peak position based on deconvolution fractions.

    Uses weighted random assignment so that if a peak has 60% of signal,
    ~60% of reads at that position will be assigned there.

    Args:
        records: List of UnifiedReadRecord (will be modified in place)
        deconv_results: Deconvolution results for regions
        position_field: Field to use for matching ('three_prime_raw')

    Returns:
        Same records list (modified in place)
    """
    # Build lookup for deconvolved regions
    # Key: (chrom, strand, position) -> list of (peak_position, fraction)
    region_lookup: Dict[Tuple[str, str, int], List[Tuple[int, float]]] = {}

    for result in deconv_results:
        if not result.peaks:
            continue

        # For all positions in this region, map to peaks
        for pos in range(result.region_start, result.region_end):
            key = (result.chrom, result.strand, pos)
            peaks_with_fractions = [
                (p['position'], p['fraction'])
                for p in result.peaks
            ]
            region_lookup[key] = peaks_with_fractions

    if not region_lookup:
        return records

    # Assign each record
    rng = np.random.default_rng(42)  # Fixed seed for reproducibility

    for record in records:
        raw_pos = getattr(record, position_field)
        key = (record.chrom, record.strand, raw_pos)

        if key not in region_lookup:
            # Not in deconvolved region - keep as-is
            continue

        peaks = region_lookup[key]
        if len(peaks) == 1:
            # Single peak - assign directly
            record.three_prime_corrected = peaks[0][0]
        else:
            # Multiple peaks - weighted random assignment
            positions = [p[0] for p in peaks]
            fractions = [p[1] for p in peaks]
            # Normalize fractions
            total = sum(fractions)
            probs = [f / total for f in fractions] if total > 0 else [1/len(peaks)] * len(peaks)
            assigned_pos = rng.choice(positions, p=probs)
            record.three_prime_corrected = assigned_pos

    return records
