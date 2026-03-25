"""
Coverage extraction and visualization for RECTIFY.

Provides functions for:
- Extracting coverage from BAM files
- Drawing coverage tracks as filled area plots
- Strand-specific coverage visualization

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def extract_coverage_from_bam(
    bam_path: Union[str, Path],
    chrom: str,
    start: int,
    end: int,
    strand: Optional[str] = None,
    min_mapq: int = 20,
    use_3prime_end: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract coverage from BAM file for a genomic region.

    Args:
        bam_path: Path to indexed BAM file
        chrom: Chromosome name
        start: Start position (0-based)
        end: End position (exclusive)
        strand: Optional strand filter ('+', '-', or None for both)
        min_mapq: Minimum mapping quality
        use_3prime_end: If True, count 3' ends instead of full coverage

    Returns:
        Tuple of (positions, depths) as numpy arrays

    Example:
        positions, depths = extract_coverage_from_bam(
            "sample.bam", "chrIV", 530000, 535000
        )
    """
    import pysam

    # Initialize arrays for the region
    region_size = end - start
    positions = np.arange(start, end)
    depths = np.zeros(region_size, dtype=np.int32)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if use_3prime_end:
            # Count 3' end positions
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < min_mapq:
                    continue

                read_strand = '-' if read.is_reverse else '+'

                # Filter by strand if specified
                if strand is not None and read_strand != strand:
                    continue

                # Get 3' end position
                if read_strand == '+':
                    pos_3prime = read.reference_end - 1  # 0-based
                else:
                    pos_3prime = read.reference_start

                # Add to coverage if in range
                if start <= pos_3prime < end:
                    depths[pos_3prime - start] += 1
        else:
            # Use pysam's count_coverage for full read coverage
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < min_mapq:
                    continue

                read_strand = '-' if read.is_reverse else '+'

                # Filter by strand if specified
                if strand is not None and read_strand != strand:
                    continue

                # Add coverage for this read
                read_start = max(read.reference_start, start)
                read_end = min(read.reference_end, end)

                for pos in range(read_start, read_end):
                    depths[pos - start] += 1

    return positions, depths


def extract_coverage_from_array(
    positions: np.ndarray,
    counts: np.ndarray,
    region_start: int,
    region_end: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract coverage for a region from pre-computed arrays.

    Useful when coverage is already loaded into memory.

    Args:
        positions: Array of genomic positions with signal
        counts: Array of counts at each position
        region_start: Start of region to extract
        region_end: End of region to extract

    Returns:
        Tuple of (positions, depths) for the region
    """
    # Create output arrays
    region_positions = np.arange(region_start, region_end)
    region_depths = np.zeros(region_end - region_start, dtype=np.float64)

    # Find positions in range
    mask = (positions >= region_start) & (positions < region_end)
    for pos, count in zip(positions[mask], counts[mask]):
        region_depths[pos - region_start] = count

    return region_positions, region_depths


def draw_coverage_track(
    ax,
    positions: np.ndarray,
    depths: np.ndarray,
    region_start: int,
    region_end: int,
    color: str = '#0072B2',
    fill: bool = True,
    fill_alpha: float = 0.5,
    line_alpha: float = 1.0,
    linewidth: float = 1.0,
    label: Optional[str] = None,
    normalize: bool = False,
    smooth_window: Optional[int] = None,
    log_scale: bool = False,
):
    """
    Draw coverage as filled area plot.

    Args:
        ax: Matplotlib axes
        positions: Array of genomic positions
        depths: Array of coverage depths
        region_start: Start of display region
        region_end: End of display region
        color: Fill and line color
        fill: Whether to fill area under curve
        fill_alpha: Transparency of fill
        line_alpha: Transparency of line
        linewidth: Width of outline
        label: Legend label
        normalize: If True, normalize to max=1
        smooth_window: Optional smoothing window size
        log_scale: If True, use log10 scale for y-axis

    Example:
        positions, depths = extract_coverage_from_bam(bam_path, chrom, start, end)
        draw_coverage_track(ax, positions, depths, start, end, label="NET-seq")
    """
    ax.set_xlim(region_start, region_end)

    # Apply smoothing if requested
    if smooth_window is not None and smooth_window > 1:
        kernel = np.ones(smooth_window) / smooth_window
        depths = np.convolve(depths, kernel, mode='same')

    # Normalize if requested
    if normalize and depths.max() > 0:
        depths = depths / depths.max()

    # Apply log scale
    if log_scale:
        depths = np.log10(depths + 1)

    # Draw fill
    if fill:
        ax.fill_between(
            positions,
            0,
            depths,
            color=color,
            alpha=fill_alpha,
            label=label,
        )

    # Draw line
    ax.plot(
        positions,
        depths,
        color=color,
        linewidth=linewidth,
        alpha=line_alpha,
        label=label if not fill else None,
    )

    # Set y-axis
    if log_scale:
        ax.set_ylabel('Coverage (log10)')
    elif normalize:
        ax.set_ylabel('Normalized')
    else:
        ax.set_ylabel('Coverage')


def draw_strand_coverage(
    ax,
    positions_plus: np.ndarray,
    depths_plus: np.ndarray,
    positions_minus: np.ndarray,
    depths_minus: np.ndarray,
    region_start: int,
    region_end: int,
    color_plus: str = '#0072B2',
    color_minus: str = '#D55E00',
    fill_alpha: float = 0.5,
    labels: Tuple[str, str] = ('+ strand', '- strand'),
    mirror: bool = True,
):
    """
    Draw strand-specific coverage with plus above axis and minus below.

    Args:
        ax: Matplotlib axes
        positions_plus: Positions for plus strand
        depths_plus: Depths for plus strand
        positions_minus: Positions for minus strand
        depths_minus: Depths for minus strand
        region_start: Start of display region
        region_end: End of display region
        color_plus: Color for plus strand
        color_minus: Color for minus strand
        fill_alpha: Fill transparency
        labels: Labels for legend
        mirror: If True, show minus strand below axis (negated)

    Example:
        # Extract strand-specific coverage
        pos_plus, dep_plus = extract_coverage_from_bam(bam, chrom, start, end, strand='+')
        pos_minus, dep_minus = extract_coverage_from_bam(bam, chrom, start, end, strand='-')

        # Draw mirrored
        draw_strand_coverage(ax, pos_plus, dep_plus, pos_minus, dep_minus, start, end)
    """
    ax.set_xlim(region_start, region_end)

    # Plus strand (above axis)
    ax.fill_between(
        positions_plus,
        0,
        depths_plus,
        color=color_plus,
        alpha=fill_alpha,
        label=labels[0],
    )
    ax.plot(positions_plus, depths_plus, color=color_plus, linewidth=0.5)

    # Minus strand (below axis if mirror=True)
    if mirror:
        ax.fill_between(
            positions_minus,
            0,
            -depths_minus,
            color=color_minus,
            alpha=fill_alpha,
            label=labels[1],
        )
        ax.plot(positions_minus, -depths_minus, color=color_minus, linewidth=0.5)
    else:
        ax.fill_between(
            positions_minus,
            0,
            depths_minus,
            color=color_minus,
            alpha=fill_alpha,
            label=labels[1],
        )
        ax.plot(positions_minus, depths_minus, color=color_minus, linewidth=0.5)

    # Add center line
    ax.axhline(y=0, color='black', linewidth=0.5, alpha=0.5)

    ax.set_ylabel('Coverage')
    ax.legend(loc='upper right', fontsize=6)


def compare_coverage_tracks(
    ax,
    coverage_dict: Dict[str, Tuple[np.ndarray, np.ndarray]],
    region_start: int,
    region_end: int,
    colors: Optional[Dict[str, str]] = None,
    normalize: bool = True,
    alpha: float = 0.7,
    show_legend: bool = True,
):
    """
    Compare multiple coverage tracks on the same axes.

    Args:
        ax: Matplotlib axes
        coverage_dict: Dict mapping label -> (positions, depths)
        region_start: Start of display region
        region_end: End of display region
        colors: Optional dict mapping label -> color
        normalize: If True, normalize each track to max=1
        alpha: Line/fill transparency
        show_legend: Whether to show legend

    Example:
        coverage = {
            'WT': (wt_positions, wt_depths),
            'mutant': (mut_positions, mut_depths),
        }
        compare_coverage_tracks(ax, coverage, start, end)
    """
    from .config import WONG_COLORS

    default_colors = list(WONG_COLORS.values())

    ax.set_xlim(region_start, region_end)

    for i, (label, (positions, depths)) in enumerate(coverage_dict.items()):
        color = colors.get(label) if colors else default_colors[i % len(default_colors)]

        if normalize and depths.max() > 0:
            depths = depths / depths.max()

        ax.plot(
            positions,
            depths,
            color=color,
            linewidth=1.5,
            alpha=alpha,
            label=label,
        )

    if show_legend:
        ax.legend(loc='upper right', fontsize=8)

    ax.set_ylabel('Normalized' if normalize else 'Coverage')
