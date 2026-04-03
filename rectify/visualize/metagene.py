"""
Metagene signal aggregation for RECTIFY visualization.

Provides efficient signal extraction and aggregation for metagene analysis:
- PositionIndex: O(1) position lookups using dict-of-dicts structure
- MetagenePipeline: Complete workflow with per-locus normalization
- LociFilter: Consistent loci filtering across multiple datasets

Key features:
- Per-locus normalization: Each locus contributes equally regardless of expression
- Percentile capping: Reduces influence of outlier loci
- Trimmed mean: Excludes top/bottom loci by signal (better than capping)
- Strand awareness: Correct handling of + and - strand coordinates
- Variable-length scaling: Feature bodies scaled to fixed width using np.interp
- Consistent loci filtering: Use one dataset to determine loci, apply to all

Performance:
- Index build: ~2s for 10M positions
- Window extraction: O(window_size) dict lookups per locus
- Full metagene (10k loci): <1 second

Author: Kevin R. Roy
"""

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd


class StrandOrientationError(Exception):
    """Raised when plus and minus strand peaks don't align in metagene analysis.

    Indicates a strand orientation bug: minus strand signals were not reversed
    to 5'→3' orientation, or center coordinates are wrong for one strand.

    Common causes:
    - Missing signal[::-1] reversal for minus strand
    - Using trt_start instead of trt_end - 1 as center for minus strand loci
    - Off-by-one in coordinate conversion (0-based vs 1-based)
    """
    pass


def verify_strand_balance(
    profiles_plus: np.ndarray,
    profiles_minus: np.ndarray,
    x: np.ndarray,
    tolerance: int = 2,
    raise_on_fail: bool = True,
    label: str = "",
) -> Dict:
    """
    Verify that plus and minus strand metagene profiles have aligned peaks.

    A correctly oriented metagene must show the same peak position for both
    strand subsets. Misaligned peaks indicate missing signal reversal
    (signal[::-1]) or wrong center coordinate convention.

    The canonical TRT metagene has its peak at ~-2 bp (CPA position) for both
    strands when signals are correctly oriented. If the minus strand peak is at
    a very different position (e.g., +133 bp), strand handling is broken.

    Args:
        profiles_plus: 2D array (n_plus x profile_length) for plus strand loci.
                       Each row is one locus's signal profile.
        profiles_minus: 2D array (n_minus x profile_length) for minus strand loci.
        x: 1D array of x-axis positions (relative to center, e.g. np.arange(-150, 151))
        tolerance: Max allowed peak position difference in bp (default: 2)
        raise_on_fail: If True, raise StrandOrientationError when peaks diverge
        label: Optional context string included in error messages

    Returns:
        Dict with keys:
            - 'plus_peak': Peak x-position for plus strand loci
            - 'minus_peak': Peak x-position for minus strand loci
            - 'difference': Absolute difference between peaks (bp)
            - 'passed': True if difference <= tolerance
            - 'n_plus': Number of plus strand loci
            - 'n_minus': Number of minus strand loci

    Raises:
        StrandOrientationError: If |plus_peak - minus_peak| > tolerance
                                and raise_on_fail=True

    Example:
        profiles_plus = extract_profiles(plus_strand_loci)
        profiles_minus = extract_profiles(minus_strand_loci)
        x = np.arange(-150, 151)

        result = verify_strand_balance(profiles_plus, profiles_minus, x)
        print(f"Plus peak: {result['plus_peak']:+d} bp")
        print(f"Minus peak: {result['minus_peak']:+d} bp")
    """
    n_plus = len(profiles_plus)
    n_minus = len(profiles_minus)

    if n_plus == 0 or n_minus == 0:
        raise ValueError(
            f"Cannot verify strand balance: need profiles for both strands "
            f"(n_plus={n_plus}, n_minus={n_minus})"
        )

    plus_mean = np.mean(profiles_plus, axis=0)
    minus_mean = np.mean(profiles_minus, axis=0)

    plus_peak = x[np.argmax(plus_mean)]
    minus_peak = x[np.argmax(minus_mean)]
    difference = abs(int(plus_peak) - int(minus_peak))
    passed = difference <= tolerance

    ctx = f" [{label}]" if label else ""
    print(
        f"STRAND VERIFICATION{ctx}: "
        f"Plus peak={plus_peak:+d} bp (n={n_plus}), "
        f"Minus peak={minus_peak:+d} bp (n={n_minus}), "
        f"Difference={difference} bp {'✓' if passed else '✗ FAIL'}"
    )

    if not passed and raise_on_fail:
        raise StrandOrientationError(
            f"Strand orientation bug detected{ctx}! "
            f"Plus strand peak at {plus_peak:+d} bp, "
            f"minus strand peak at {minus_peak:+d} bp "
            f"(difference={difference} bp, tolerance={tolerance} bp). "
            f"Check: (1) signal[::-1] applied for minus strand, "
            f"(2) window coordinates flipped for minus strand, "
            f"(3) center = trt_end - 1 (not trt_start) for minus strand loci."
        )

    return {
        'plus_peak': int(plus_peak),
        'minus_peak': int(minus_peak),
        'difference': difference,
        'passed': passed,
        'n_plus': n_plus,
        'n_minus': n_minus,
    }


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
        trimmed_proportion: Proportion to cut from each end for trimmed mean (0.1 = 10%)
        use_trimmed_mean: If True, use trimmed mean instead of percentile capping
    """
    window_upstream: int = 100
    window_downstream: int = 100
    scaled_body_width: int = 20
    normalize_per_locus: bool = True
    cap_percentile: float = 90.0
    min_locus_signal: float = 0.0
    trimmed_proportion: float = 0.1
    use_trimmed_mean: bool = False


@dataclass
class LociFilter:
    """
    Filter for consistent loci selection across multiple datasets.

    When comparing metagene profiles from different data types (e.g., RECTIFY,
    NET-seq, PAR-CLIP), outlier loci may differ between datasets. Using the
    same loci across all datasets ensures valid comparisons.

    The filter is computed from a "priority" dataset and then applied to all
    other datasets, ensuring they use the exact same set of loci.

    Attributes:
        kept_indices: Indices of loci to keep (relative to original DataFrame)
        excluded_indices: Indices of loci excluded
        priority_dataset: Name of the dataset used to determine filtering
        method: Filtering method used ('trimmed_mean', 'top_percentile', etc.)
        parameters: Dict of method parameters

    Example:
        # Compute filter from RECTIFY WT data
        filter = pipeline.compute_loci_filter(
            loci_df, rectify_wt_index,
            method='trimmed_mean',
            proportion=0.1  # Keep middle 80%
        )

        # Apply same filter to all datasets
        results = {
            'RECTIFY_WT': pipeline.compute_profile(loci_df, rectify_wt_index, loci_filter=filter),
            'RECTIFY_dst1': pipeline.compute_profile(loci_df, rectify_dst1_index, loci_filter=filter),
            'NET-seq_WT': pipeline.compute_profile(loci_df, netseq_wt_index, loci_filter=filter),
        }
    """
    kept_indices: List[int]
    excluded_indices: List[int]
    priority_dataset: str = "unknown"
    method: str = "trimmed_mean"
    parameters: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_kept(self) -> int:
        """Number of loci kept."""
        return len(self.kept_indices)

    @property
    def n_excluded(self) -> int:
        """Number of loci excluded."""
        return len(self.excluded_indices)

    @property
    def fraction_kept(self) -> float:
        """Fraction of loci kept."""
        total = self.n_kept + self.n_excluded
        return self.n_kept / total if total > 0 else 0.0

    def apply_to_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply filter to DataFrame, keeping only selected loci."""
        return df.iloc[self.kept_indices].copy()

    def __repr__(self) -> str:
        return (
            f"LociFilter(n_kept={self.n_kept}, n_excluded={self.n_excluded}, "
            f"method='{self.method}', priority='{self.priority_dataset}')"
        )


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

    def compute_loci_filter(
        self,
        loci: pd.DataFrame,
        position_index: PositionIndex,
        method: str = 'trimmed_mean',
        proportion: float = 0.1,
        priority_name: str = 'priority',
        chrom_col: str = 'chrom',
        strand_col: str = 'strand',
        locus_start_col: str = 'start',
        locus_end_col: str = 'end',
    ) -> LociFilter:
        """
        Compute a loci filter based on signal from a priority dataset.

        This determines which loci to include/exclude based on signal levels
        in the priority dataset. The resulting filter can be applied to other
        datasets to ensure consistent loci selection.

        Methods:
        - 'trimmed_mean': Exclude top and bottom proportion of loci by total signal
        - 'top_percentile': Exclude top proportion of loci (highest signal)
        - 'bottom_percentile': Exclude bottom proportion of loci (lowest signal)

        Args:
            loci: DataFrame with locus coordinates
            position_index: PositionIndex for the priority dataset
            method: Filtering method ('trimmed_mean', 'top_percentile', 'bottom_percentile')
            proportion: Proportion to exclude (0.1 = 10% from each end for trimmed_mean)
            priority_name: Name for the priority dataset (for tracking)
            chrom_col: Column name for chromosome
            strand_col: Column name for strand
            locus_start_col: Column name for locus start position
            locus_end_col: Column name for locus end position

        Returns:
            LociFilter object with kept/excluded indices

        Example:
            # Create filter from RECTIFY WT, excluding top/bottom 10%
            filt = pipeline.compute_loci_filter(
                loci_df, rectify_wt_index,
                method='trimmed_mean',
                proportion=0.1,
                priority_name='RECTIFY_WT'
            )
            print(f"Keeping {filt.n_kept} of {filt.n_kept + filt.n_excluded} loci")
        """
        cfg = self.config

        # Calculate total signal for each locus
        locus_signals = []
        locus_indices = []

        for idx, locus in loci.iterrows():
            chrom = locus[chrom_col]
            strand = locus[strand_col]
            locus_start = int(locus[locus_start_col])
            locus_end = int(locus[locus_end_col])

            # Sum signal in the window around this locus
            if strand == '+':
                window_start = locus_start - cfg.window_upstream
                window_end = locus_end + cfg.window_downstream
            else:
                window_start = locus_start - cfg.window_downstream
                window_end = locus_end + cfg.window_upstream

            signal = position_index.extract_window_array(
                chrom, strand, window_start, window_end
            )
            total = signal.sum()

            locus_signals.append(total)
            locus_indices.append(idx)

        # Convert to arrays
        signals = np.array(locus_signals)
        indices = np.array(locus_indices)
        n_loci = len(signals)

        # Determine which loci to keep based on method
        if method == 'trimmed_mean':
            # Exclude top and bottom proportion
            n_cut = int(n_loci * proportion)
            sorted_order = np.argsort(signals)
            # Keep middle portion
            if n_cut > 0:
                kept_positions = sorted_order[n_cut:-n_cut]
                excluded_positions = np.concatenate([
                    sorted_order[:n_cut],
                    sorted_order[-n_cut:]
                ])
            else:
                kept_positions = sorted_order
                excluded_positions = np.array([], dtype=int)

        elif method == 'top_percentile':
            # Exclude top proportion (highest signal)
            n_cut = int(n_loci * proportion)
            sorted_order = np.argsort(signals)
            if n_cut > 0:
                kept_positions = sorted_order[:-n_cut]
                excluded_positions = sorted_order[-n_cut:]
            else:
                kept_positions = sorted_order
                excluded_positions = np.array([], dtype=int)

        elif method == 'bottom_percentile':
            # Exclude bottom proportion (lowest signal)
            n_cut = int(n_loci * proportion)
            sorted_order = np.argsort(signals)
            if n_cut > 0:
                kept_positions = sorted_order[n_cut:]
                excluded_positions = sorted_order[:n_cut]
            else:
                kept_positions = sorted_order
                excluded_positions = np.array([], dtype=int)

        else:
            raise ValueError(f"Unknown method: {method}. Use 'trimmed_mean', 'top_percentile', or 'bottom_percentile'")

        # Convert positions back to original DataFrame indices
        kept_indices = [indices[i] for i in kept_positions]
        excluded_indices = [indices[i] for i in excluded_positions]

        return LociFilter(
            kept_indices=kept_indices,
            excluded_indices=excluded_indices,
            priority_dataset=priority_name,
            method=method,
            parameters={'proportion': proportion}
        )

    def compute_profile_filtered(
        self,
        loci: pd.DataFrame,
        position_index: PositionIndex,
        loci_filter: LociFilter,
        **kwargs,
    ) -> Dict:
        """
        Compute metagene profile using only loci specified by a filter.

        This is the primary method for computing profiles with consistent
        loci filtering. The filter is typically computed from a priority
        dataset using compute_loci_filter().

        Args:
            loci: DataFrame with locus coordinates
            position_index: PositionIndex with signal data
            loci_filter: LociFilter specifying which loci to include
            **kwargs: Additional arguments passed to compute_profile

        Returns:
            Dict with profile results (same format as compute_profile)

        Example:
            # Compute filter from RECTIFY WT
            filt = pipeline.compute_loci_filter(loci_df, rectify_wt_index)

            # Apply to all datasets
            results = {
                name: pipeline.compute_profile_filtered(loci_df, idx, filt)
                for name, idx in datasets.items()
            }
        """
        # Filter loci DataFrame
        filtered_loci = loci_filter.apply_to_dataframe(loci)

        # Compute profile on filtered loci
        result = self.compute_profile(filtered_loci, position_index, **kwargs)

        # Add filter info to result
        result['loci_filter'] = loci_filter
        result['filter_method'] = loci_filter.method
        result['filter_priority'] = loci_filter.priority_dataset

        return result

    def compute_multiple_profiles_filtered(
        self,
        loci: pd.DataFrame,
        condition_indices: Dict[str, PositionIndex],
        priority_condition: str,
        filter_method: str = 'trimmed_mean',
        filter_proportion: float = 0.1,
        **kwargs,
    ) -> Tuple[Dict[str, Dict], LociFilter]:
        """
        Compute profiles for multiple conditions with consistent loci filtering.

        First computes a loci filter from the priority condition, then applies
        it to all conditions. This ensures all profiles use the exact same loci.

        Args:
            loci: DataFrame with locus coordinates
            condition_indices: Dict mapping condition name -> PositionIndex
            priority_condition: Name of condition to use for determining loci filter
            filter_method: Filtering method ('trimmed_mean', 'top_percentile', 'bottom_percentile')
            filter_proportion: Proportion to exclude
            **kwargs: Additional arguments passed to compute_profile

        Returns:
            Tuple of (results dict, LociFilter used)

        Example:
            results, filt = pipeline.compute_multiple_profiles_filtered(
                loci_df,
                {'WT': wt_index, 'mutant': mut_index},
                priority_condition='WT',
                filter_method='trimmed_mean',
                filter_proportion=0.1
            )
            print(f"Using {filt.n_kept} loci based on WT signal")
        """
        if priority_condition not in condition_indices:
            raise ValueError(
                f"Priority condition '{priority_condition}' not in condition_indices. "
                f"Available: {list(condition_indices.keys())}"
            )

        # Compute filter from priority condition
        loci_filter = self.compute_loci_filter(
            loci,
            condition_indices[priority_condition],
            method=filter_method,
            proportion=filter_proportion,
            priority_name=priority_condition,
            **{k: v for k, v in kwargs.items()
               if k in ['chrom_col', 'strand_col', 'locus_start_col', 'locus_end_col']}
        )

        # Compute profiles for all conditions using the same filter
        results = {}
        for condition, index in condition_indices.items():
            results[condition] = self.compute_profile_filtered(
                loci, index, loci_filter, **kwargs
            )

        return results, loci_filter

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

    def compute_center_profile(
        self,
        loci: List[Dict],
        position_index: "PositionIndex",
        window: Tuple[int, int],
        normalize: bool = True,
        total_reads: Optional[int] = None,
        verify_strands: bool = True,
        strand_tolerance: int = 2,
        cap_percentile: Optional[int] = None,
    ) -> Dict:
        """
        Compute metagene profile centered on a single position per locus.

        Preferred method for TRT/CPA metagene analysis where each locus has a
        single center position (e.g., first T of T-tract). Handles strand-aware
        coordinate transformation internally — callers do NOT need to flip
        coordinates or reverse signal arrays.

        Strand handling (automatic, always applied):
        - Plus strand:  genomic_start = center + window[0]
                        genomic_end   = center + window[1]
        - Minus strand: genomic_start = center - window[1]  (flipped)
                        genomic_end   = center - window[0]
                        signal array reversed [::-1] to 5'→3' orientation

        If verify_strands=True (default), calls verify_strand_balance() after
        extracting all profiles and raises StrandOrientationError if plus and
        minus strand peaks diverge by more than strand_tolerance bp. This is
        the primary safety net against strand bugs.

        Args:
            loci: List of dicts, each with keys:
                      'chrom' (str), 'strand' ('+'/'-'), 'center' (int, 0-based)
            position_index: PositionIndex with signal data
            window: (window_start, window_end) relative to center (both inclusive),
                    e.g. (-150, 150) extracts 301 bp centered on center position.
                    window_start should be negative (upstream); window_end positive.
            normalize: If True and total_reads provided, apply RPM normalization
                       (multiply by 1e6 / total_reads)
            total_reads: Total mapped reads for RPM normalization (required if normalize=True)
            verify_strands: If True (default), run verify_strand_balance() and
                            raise StrandOrientationError if peaks diverge
            strand_tolerance: Max allowed bp difference between + and - strand
                              peaks (default: 2)
            cap_percentile: If set (e.g. 50), apply window-sum capping at this
                            percentile to each strand's profiles before running
                            strand verification. Reduces influence of outlier loci
                            (e.g. extremely high-expression genes) on peak detection.
                            Does NOT affect the returned profile_matrix or plots —
                            only affects which peak is found for strand verification.
                            Recommended: 50 (matches apply_window_sum_capping default).

        Returns:
            Dict with keys:
                - 'profile': Mean metagene profile across all loci (1D array)
                - 'profile_matrix': Per-locus profiles (2D array, n_loci x profile_length)
                - 'sem': Standard error of mean (1D array)
                - 'n_loci': Number of loci included
                - 'x': X-axis positions (1D array, e.g. np.arange(-150, 151))
                - 'window': The window tuple used
                - 'profiles_plus': Per-locus profiles for + strand loci (2D array)
                - 'profiles_minus': Per-locus profiles for - strand loci (2D array)
                - 'strand_verification': Dict from verify_strand_balance() (if run)
                - 'n_plus': Number of plus strand loci
                - 'n_minus': Number of minus strand loci

        Raises:
            StrandOrientationError: If verify_strands=True and peaks diverge

        Example:
            # Load TRT loci cache and RECTIFY position index
            loci = load_trt_loci_cache("cached_trt_signals_v3_no_poliii.pkl")
            index, total_reads = load_rectify_position_index(tsv_files)

            # Compute metagene — strand verification runs automatically
            pipeline = MetagenePipeline()
            result = pipeline.compute_center_profile(
                loci, index,
                window=(-150, 150),
                total_reads=total_reads,
            )
            print(f"Peak at: {result['x'][np.argmax(result['profile'])]:+d} bp")
        """
        window_start, window_end = window
        profile_length = window_end - window_start + 1
        x = np.arange(window_start, window_end + 1)

        profiles_all = []
        profiles_plus = []
        profiles_minus = []

        for locus in loci:
            chrom = locus['chrom']
            strand = locus['strand']
            center = int(locus['center'])

            # Strand-aware coordinate transformation
            if strand == '+':
                gstart = center + window_start
                gend = center + window_end
            else:
                # Minus strand: flip window in genomic coordinates
                gstart = center - window_end
                gend = center - window_start

            signal = position_index.extract_window_array(chrom, strand, gstart, gend)

            # Reverse minus strand signal to 5'→3' orientation
            if strand == '-':
                signal = signal[::-1]

            # RPM normalization
            if normalize and total_reads is not None and total_reads > 0:
                signal = signal * (1e6 / total_reads)

            profiles_all.append(signal)
            if strand == '+':
                profiles_plus.append(signal)
            else:
                profiles_minus.append(signal)

        if not profiles_all:
            return {
                'profile': np.zeros(profile_length),
                'profile_matrix': np.zeros((0, profile_length)),
                'sem': np.zeros(profile_length),
                'n_loci': 0,
                'x': x,
                'window': window,
                'profiles_plus': np.zeros((0, profile_length)),
                'profiles_minus': np.zeros((0, profile_length)),
                'strand_verification': None,
                'n_plus': 0,
                'n_minus': 0,
            }

        arr_all = np.array(profiles_all)
        arr_plus = np.array(profiles_plus) if profiles_plus else np.zeros((0, profile_length))
        arr_minus = np.array(profiles_minus) if profiles_minus else np.zeros((0, profile_length))

        # Strand balance verification — the primary safety net
        strand_verification = None
        if verify_strands and len(arr_plus) > 0 and len(arr_minus) > 0:
            # Optionally cap outlier loci before peak detection to prevent a few
            # high-expression genes from dominating the aggregate profile shape.
            if cap_percentile is not None:
                from .figure_utils import apply_window_sum_capping
                plus_for_verify, _ = apply_window_sum_capping(arr_plus, percentile=cap_percentile)
                minus_for_verify, _ = apply_window_sum_capping(arr_minus, percentile=cap_percentile)
            else:
                plus_for_verify = arr_plus
                minus_for_verify = arr_minus
            strand_verification = verify_strand_balance(
                plus_for_verify,
                minus_for_verify,
                x,
                tolerance=strand_tolerance,
                raise_on_fail=True,
            )

        profile_mean = np.mean(arr_all, axis=0)
        profile_sem = np.std(arr_all, axis=0) / np.sqrt(len(arr_all))

        return {
            'profile': profile_mean,
            'profile_matrix': arr_all,
            'sem': profile_sem,
            'n_loci': len(profiles_all),
            'x': x,
            'window': window,
            'profiles_plus': arr_plus,
            'profiles_minus': arr_minus,
            'strand_verification': strand_verification,
            'n_plus': len(profiles_plus),
            'n_minus': len(profiles_minus),
        }

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
