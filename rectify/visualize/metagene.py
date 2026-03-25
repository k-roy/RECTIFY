"""
Metagene signal aggregation for RECTIFY visualization.

Provides efficient signal extraction and aggregation for metagene analysis:
- PositionIndex: O(1) position lookups using dict-of-dicts structure
- MetagenePipeline: Complete workflow with per-locus normalization

Key features:
- Per-locus normalization: Each locus contributes equally regardless of expression
- Percentile capping: Reduces influence of outlier loci
- Strand awareness: Correct handling of + and - strand coordinates
- Variable-length scaling: Feature bodies scaled to fixed width using np.interp

Performance:
- Index build: ~2s for 10M positions
- Window extraction: O(window_size) dict lookups per locus
- Full metagene (10k loci): <1 second

Author: Kevin R. Roy
"""

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd


@dataclass
class MetageneConfig:
    """Configuration for metagene analysis.

    Attributes:
        window_upstream: Number of bp upstream of feature start
        window_downstream: Number of bp downstream of feature end
        scaled_body_width: Scale variable-length feature bodies to this width
        normalize_per_locus: If True, each locus contributes equally
        cap_percentile: Cap values at this percentile (reduces outlier influence)
        min_locus_signal: Minimum total signal for a locus to be included
    """
    window_upstream: int = 100
    window_downstream: int = 100
    scaled_body_width: int = 20
    normalize_per_locus: bool = True
    cap_percentile: float = 90.0
    min_locus_signal: float = 0.0


class PositionIndex:
    """
    Efficient position lookup using dict-of-dicts structure.

    Stores signal counts at genomic positions for fast O(1) lookups.
    Structure: index[(chrom, strand)][position] = count

    Advantages over pkl caching or BAM scanning:
    - O(1) lookup per position
    - Built on-demand from any DataFrame (~2s for 10M positions)
    - No cache invalidation issues
    - Memory efficient (~200MB for 10M unique positions)

    Example:
        # Build from 3' end positions DataFrame
        ends_df = pd.read_csv("3prime_ends.tsv", sep='\\t')
        index = PositionIndex(ends_df, position_col='polished_3prime')

        # Query single position
        count = index.count_at('chrIV', '+', 530000)

        # Extract window as array
        signal = index.extract_window_array('chrIV', '+', 530000, 530200)
    """

    def __init__(
        self,
        df: pd.DataFrame,
        position_col: str = 'position',
        chrom_col: str = 'chrom',
        strand_col: str = 'strand',
        count_col: Optional[str] = None,
    ):
        """Build index from DataFrame.

        Args:
            df: DataFrame with position data
            position_col: Column name for genomic position
            chrom_col: Column name for chromosome
            strand_col: Column name for strand (+/-)
            count_col: Optional column with pre-computed counts.
                       If None, counts occurrences of each position.
        """
        self._index: Dict[Tuple[str, str], Dict[int, int]] = {}
        self._n_positions = 0
        self._total_counts = 0
        self._build_index(df, position_col, chrom_col, strand_col, count_col)

    def _build_index(
        self,
        df: pd.DataFrame,
        position_col: str,
        chrom_col: str,
        strand_col: str,
        count_col: Optional[str],
    ):
        """Build index using Counter for efficiency."""
        if count_col is not None:
            # Use pre-computed counts
            for _, row in df.iterrows():
                key = (row[chrom_col], row[strand_col])
                pos = int(row[position_col])
                count = int(row[count_col])

                if key not in self._index:
                    self._index[key] = {}
                self._index[key][pos] = self._index[key].get(pos, 0) + count
                self._total_counts += count
        else:
            # Count occurrences using Counter (much faster than iterrows)
            keys = list(zip(
                df[chrom_col],
                df[strand_col],
                df[position_col].astype(int)
            ))
            counts = Counter(keys)

            for (chrom, strand, pos), count in counts.items():
                key = (chrom, strand)
                if key not in self._index:
                    self._index[key] = {}
                self._index[key][pos] = count
                self._total_counts += count

        self._n_positions = sum(len(d) for d in self._index.values())

    @property
    def n_positions(self) -> int:
        """Number of unique positions in index."""
        return self._n_positions

    @property
    def total_counts(self) -> int:
        """Total signal counts across all positions."""
        return self._total_counts

    @property
    def chromosomes(self) -> List[str]:
        """List of chromosomes in index."""
        return sorted(set(chrom for chrom, strand in self._index.keys()))

    def count_at(self, chrom: str, strand: str, position: int) -> int:
        """Get count at single position.

        Args:
            chrom: Chromosome name
            strand: Strand ('+' or '-')
            position: Genomic position

        Returns:
            Count at position (0 if not in index)
        """
        return self._index.get((chrom, strand), {}).get(position, 0)

    def sum_in_window(self, chrom: str, strand: str, start: int, end: int) -> int:
        """Sum counts in [start, end] window (inclusive).

        Args:
            chrom: Chromosome name
            strand: Strand ('+' or '-')
            start: Start position (inclusive)
            end: End position (inclusive)

        Returns:
            Sum of counts in window
        """
        pos_counts = self._index.get((chrom, strand), {})
        return sum(pos_counts.get(pos, 0) for pos in range(start, end + 1))

    def extract_window_array(
        self,
        chrom: str,
        strand: str,
        start: int,
        end: int,
    ) -> np.ndarray:
        """Extract signal as numpy array for a window.

        Args:
            chrom: Chromosome name
            strand: Strand ('+' or '-')
            start: Start position (inclusive)
            end: End position (inclusive)

        Returns:
            1D numpy array of counts for each position in [start, end]
        """
        pos_counts = self._index.get((chrom, strand), {})
        return np.array([pos_counts.get(pos, 0) for pos in range(start, end + 1)])

    def get_positions_in_window(
        self,
        chrom: str,
        strand: str,
        start: int,
        end: int,
    ) -> Dict[int, int]:
        """Get dict of positions with non-zero counts in window.

        Args:
            chrom: Chromosome name
            strand: Strand ('+' or '-')
            start: Start position (inclusive)
            end: End position (inclusive)

        Returns:
            Dict mapping position -> count for positions with signal
        """
        pos_counts = self._index.get((chrom, strand), {})
        return {
            pos: count
            for pos, count in pos_counts.items()
            if start <= pos <= end
        }


class MetagenePipeline:
    """
    Complete metagene workflow with per-locus normalization.

    Aggregates signal across multiple genomic loci (e.g., TRTs, TSS, TES)
    to create averaged metagene profiles.

    Key features:
    - Per-locus normalization: Each locus contributes equally regardless of expression
    - Percentile capping: Reduces influence of outlier loci
    - Strand awareness: Correct handling of + and - strand coordinates
    - Variable-length scaling: Feature bodies scaled to fixed width using np.interp

    Example:
        # Load data
        ends_df = pd.read_csv("3prime_ends.tsv", sep='\\t')
        loci_df = pd.read_csv("trts.tsv", sep='\\t')

        # Build index
        index = PositionIndex(ends_df, position_col='position')

        # Compute profile
        config = MetageneConfig(window_upstream=100, window_downstream=100)
        pipeline = MetagenePipeline(config)
        result = pipeline.compute_profile(loci_df, index)

        # Plot
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        pipeline.plot_profile(ax, result)
    """

    def __init__(self, config: Optional[MetageneConfig] = None):
        """Initialize pipeline with configuration.

        Args:
            config: MetageneConfig with analysis parameters.
                    Uses defaults if not provided.
        """
        self.config = config or MetageneConfig()

    def compute_profile(
        self,
        loci: pd.DataFrame,
        position_index: PositionIndex,
        chrom_col: str = 'chrom',
        strand_col: str = 'strand',
        locus_start_col: str = 'start',
        locus_end_col: str = 'end',
    ) -> Dict:
        """
        Compute metagene profile across all loci.

        Args:
            loci: DataFrame with locus coordinates
            position_index: PositionIndex with signal data
            chrom_col: Column name for chromosome
            strand_col: Column name for strand
            locus_start_col: Column name for locus start position
            locus_end_col: Column name for locus end position

        Returns:
            Dict with keys:
                - 'profile': Aggregated metagene profile (1D array)
                - 'profile_matrix': Per-locus profiles (2D array, n_loci x profile_length)
                - 'sem': Standard error of mean at each position
                - 'n_loci': Number of loci included
                - 'x_positions': X-axis positions relative to feature
                - 'body_width': Scaled body width
        """
        cfg = self.config
        upstream_profiles = []
        body_profiles = []
        downstream_profiles = []
        included_loci = []

        for idx, locus in loci.iterrows():
            chrom = locus[chrom_col]
            strand = locus[strand_col]
            locus_start = int(locus[locus_start_col])
            locus_end = int(locus[locus_end_col])

            # Extract regions (strand-aware)
            if strand == '+':
                up = position_index.extract_window_array(
                    chrom, strand,
                    locus_start - cfg.window_upstream,
                    locus_start - 1
                )
                body_raw = position_index.extract_window_array(
                    chrom, strand,
                    locus_start,
                    locus_end
                )
                down = position_index.extract_window_array(
                    chrom, strand,
                    locus_end + 1,
                    locus_end + cfg.window_downstream
                )
            else:
                # Minus strand: coordinates are flipped
                # "upstream" in transcription sense = higher genomic coords
                up = position_index.extract_window_array(
                    chrom, strand,
                    locus_end + 1,
                    locus_end + cfg.window_upstream
                )[::-1]  # Reverse to match transcription direction

                body_raw = position_index.extract_window_array(
                    chrom, strand,
                    locus_start,
                    locus_end
                )[::-1]

                down = position_index.extract_window_array(
                    chrom, strand,
                    locus_start - cfg.window_downstream,
                    locus_start - 1
                )[::-1]

            # Scale body to fixed width
            if len(body_raw) > 1:
                body_scaled = np.interp(
                    np.linspace(0, len(body_raw) - 1, cfg.scaled_body_width),
                    np.arange(len(body_raw)),
                    body_raw
                )
            elif len(body_raw) == 1:
                body_scaled = np.full(cfg.scaled_body_width, body_raw[0] / cfg.scaled_body_width)
            else:
                body_scaled = np.zeros(cfg.scaled_body_width)

            # Calculate total signal for this locus
            total_signal = up.sum() + body_scaled.sum() + down.sum()

            # Skip loci with insufficient signal
            if total_signal < cfg.min_locus_signal:
                continue

            # Per-locus normalization
            if cfg.normalize_per_locus and total_signal > 0:
                up = up / total_signal
                body_scaled = body_scaled / total_signal
                down = down / total_signal

            upstream_profiles.append(up)
            body_profiles.append(body_scaled)
            downstream_profiles.append(down)
            included_loci.append(idx)

        if not upstream_profiles:
            # Return empty result
            profile_length = cfg.window_upstream + cfg.scaled_body_width + cfg.window_downstream
            return {
                'profile': np.zeros(profile_length),
                'profile_matrix': np.zeros((0, profile_length)),
                'sem': np.zeros(profile_length),
                'n_loci': 0,
                'x_positions': np.arange(profile_length) - cfg.window_upstream,
                'body_width': cfg.scaled_body_width,
            }

        # Stack arrays
        up_arr = np.vstack(upstream_profiles)
        body_arr = np.vstack(body_profiles)
        down_arr = np.vstack(downstream_profiles)

        # Apply percentile capping
        if cfg.cap_percentile < 100:
            all_vals = np.concatenate([
                up_arr.flatten(),
                body_arr.flatten(),
                down_arr.flatten()
            ])
            nonzero = all_vals[all_vals > 0]
            if len(nonzero) > 0:
                cap = np.percentile(nonzero, cfg.cap_percentile)
                up_arr = np.clip(up_arr, 0, cap)
                body_arr = np.clip(body_arr, 0, cap)
                down_arr = np.clip(down_arr, 0, cap)

        # Concatenate to full profile matrix
        profile_matrix = np.hstack([up_arr, body_arr, down_arr])

        # Calculate mean and SEM
        profile_mean = np.mean(profile_matrix, axis=0)
        profile_sem = np.std(profile_matrix, axis=0) / np.sqrt(len(profile_matrix))

        # Create x-axis positions
        x_upstream = np.arange(-cfg.window_upstream, 0)
        x_body = np.linspace(0, cfg.scaled_body_width, cfg.scaled_body_width, endpoint=False)
        x_downstream = np.arange(cfg.scaled_body_width, cfg.scaled_body_width + cfg.window_downstream)
        x_positions = np.concatenate([x_upstream, x_body, x_downstream])

        return {
            'profile': profile_mean,
            'profile_matrix': profile_matrix,
            'sem': profile_sem,
            'n_loci': len(included_loci),
            'x_positions': x_positions,
            'body_width': cfg.scaled_body_width,
            'included_loci': included_loci,
        }

    def compute_multiple_profiles(
        self,
        loci: pd.DataFrame,
        condition_indices: Dict[str, PositionIndex],
        **kwargs,
    ) -> Dict[str, Dict]:
        """Compute profiles for multiple conditions.

        Args:
            loci: DataFrame with locus coordinates
            condition_indices: Dict mapping condition name -> PositionIndex
            **kwargs: Additional arguments passed to compute_profile

        Returns:
            Dict mapping condition name -> profile result dict
        """
        results = {}
        for condition, index in condition_indices.items():
            results[condition] = self.compute_profile(loci, index, **kwargs)
        return results

    def plot_profile(
        self,
        ax,
        result: Dict,
        color: str = '#0072B2',
        label: Optional[str] = None,
        show_sem: bool = True,
        sem_alpha: float = 0.2,
        linewidth: float = 1.5,
        show_feature_bounds: bool = True,
        feature_color: str = 'red',
    ):
        """Plot a single metagene profile.

        Args:
            ax: Matplotlib axes
            result: Result dict from compute_profile
            color: Line color
            label: Legend label
            show_sem: Whether to show SEM shading
            sem_alpha: Alpha for SEM shading
            linewidth: Line width
            show_feature_bounds: Whether to show vertical lines at feature start/end
            feature_color: Color for feature boundary lines
        """
        x = result['x_positions']
        y = result['profile']
        sem = result['sem']
        body_width = result['body_width']

        # Plot profile line
        ax.plot(x, y, color=color, linewidth=linewidth, label=label)

        # Add SEM shading
        if show_sem:
            ax.fill_between(
                x,
                y - sem,
                y + sem,
                color=color,
                alpha=sem_alpha,
            )

        # Add feature boundary lines
        if show_feature_bounds:
            ax.axvline(0, color=feature_color, linestyle='--', alpha=0.5, linewidth=1)
            ax.axvline(body_width, color=feature_color, linestyle='--', alpha=0.5, linewidth=1)

        # Labels
        ax.set_xlabel('Position relative to feature')
        ax.set_ylabel('Signal (normalized)' if self.config.normalize_per_locus else 'Signal')

    def plot_multiple_profiles(
        self,
        ax,
        results: Dict[str, Dict],
        colors: Optional[Dict[str, str]] = None,
        show_sem: bool = True,
        show_legend: bool = True,
        show_feature_bounds: bool = True,
    ):
        """Plot multiple metagene profiles on same axes.

        Args:
            ax: Matplotlib axes
            results: Dict mapping condition name -> result dict
            colors: Optional dict mapping condition name -> color
            show_sem: Whether to show SEM shading
            show_legend: Whether to show legend
            show_feature_bounds: Whether to show feature boundary lines

        Example:
            results = pipeline.compute_multiple_profiles(loci, {
                'WT': wt_index,
                'mutant': mut_index,
            })
            pipeline.plot_multiple_profiles(ax, results, colors={
                'WT': '#0072B2',
                'mutant': '#D55E00',
            })
        """
        from .config import WONG_COLORS

        # Default colors from Wong palette
        default_colors = list(WONG_COLORS.values())

        for i, (condition, result) in enumerate(results.items()):
            color = colors.get(condition) if colors else default_colors[i % len(default_colors)]

            # Only show feature bounds for first condition
            self.plot_profile(
                ax,
                result,
                color=color,
                label=f"{condition} (n={result['n_loci']})",
                show_sem=show_sem,
                show_feature_bounds=(i == 0) if show_feature_bounds else False,
            )

        if show_legend:
            ax.legend(loc='upper right', fontsize=8, framealpha=0.8)


def build_index_from_bam(
    bam_path: Union[str, Path],
    chrom_map: Optional[Dict[str, str]] = None,
    min_mapq: int = 20,
    use_3prime: bool = True,
) -> PositionIndex:
    """Build PositionIndex directly from BAM file.

    Convenience function that extracts 3' (or 5') end positions from a BAM file
    and builds a PositionIndex.

    Args:
        bam_path: Path to BAM file
        chrom_map: Optional dict to map BAM chromosome names (e.g., RefSeq -> chrI)
        min_mapq: Minimum mapping quality
        use_3prime: If True, extract 3' ends; if False, extract 5' ends

    Returns:
        PositionIndex built from BAM data
    """
    import pysam

    positions = []

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue

            chrom = read.reference_name
            if chrom_map:
                chrom = chrom_map.get(chrom, chrom)

            strand = '-' if read.is_reverse else '+'

            if use_3prime:
                # 3' end
                if strand == '+':
                    pos = read.reference_end - 1  # 0-based
                else:
                    pos = read.reference_start
            else:
                # 5' end
                if strand == '+':
                    pos = read.reference_start
                else:
                    pos = read.reference_end - 1

            positions.append({
                'chrom': chrom,
                'strand': strand,
                'position': pos,
            })

    df = pd.DataFrame(positions)
    return PositionIndex(df)
