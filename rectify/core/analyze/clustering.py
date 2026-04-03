#!/usr/bin/env python3
"""
CPA Cluster Formation Module

Forms clusters from individual CPA (cleavage and polyadenylation) sites.
Supports both fixed-distance and adaptive valley-based clustering.

Author: Kevin R. Roy
Date: 2026-03-17
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool, cpu_count
import numpy as np
import pandas as pd

# Optional: IntervalTree for O(log n) cluster lookups
try:
    from intervaltree import IntervalTree
    HAS_INTERVALTREE = True
except ImportError:
    HAS_INTERVALTREE = False

# Default clustering parameters
DEFAULT_CLUSTER_DISTANCE = 25  # bp for fixed-distance clustering
DEFAULT_MAX_CLUSTER_RADIUS = 10  # bp for adaptive clustering
DEFAULT_MIN_PEAK_SEPARATION = 5  # bp between distinct peaks
DEFAULT_MIN_READS = 5  # minimum reads per cluster
DEFAULT_MIN_SAMPLES = 2  # cluster must appear in N samples


def _cluster_chrom_strand_group(args):
    """Worker function for parallel clustering of a single chrom/strand group."""
    chrom, strand, positions, counts, cluster_distance, min_reads, start_id = args

    if len(positions) == 0:
        return []

    # Find cluster boundaries using vectorized gap detection
    gaps = np.diff(positions)
    split_points = np.where(gaps > cluster_distance)[0] + 1

    cluster_starts = np.concatenate([[0], split_points])
    cluster_ends = np.concatenate([split_points, [len(positions)]])

    clusters = []
    local_id = 0

    for start_idx, end_idx in zip(cluster_starts, cluster_ends):
        cluster_positions = positions[start_idx:end_idx]
        cluster_counts = counts[start_idx:end_idx]

        total_reads = cluster_counts.sum()
        if total_reads < min_reads:
            continue

        cluster_start = cluster_positions.min()
        cluster_end = cluster_positions.max()
        modal_position = int(np.median(cluster_positions))

        clusters.append({
            'cluster_id': f'cluster_{start_id + local_id:06d}',
            'chrom': chrom,
            'strand': strand,
            'start': int(cluster_start),
            'end': int(cluster_end),
            'modal_position': modal_position,
            'n_positions': len(cluster_positions),
            'n_reads': int(total_reads),
        })
        local_id += 1

    return clusters


def cluster_cpa_sites(
    positions_df: pd.DataFrame,
    cluster_distance: int = DEFAULT_CLUSTER_DISTANCE,
    min_reads: int = DEFAULT_MIN_READS,
    chrom_col: str = 'chrom',
    strand_col: str = 'strand',
    position_col: str = 'corrected_position',
    count_col: Optional[str] = None,
    n_workers: int = 1,
) -> pd.DataFrame:
    """
    Cluster CPA sites using fixed-distance merging.

    Groups nearby 3' end positions (within cluster_distance bp) into clusters.
    Uses vectorized operations for efficiency. Optionally parallelizes across
    chromosome/strand combinations.

    Args:
        positions_df: DataFrame with CPA positions
        cluster_distance: Maximum distance (bp) to merge positions
        min_reads: Minimum total reads per cluster
        chrom_col: Column name for chromosome
        strand_col: Column name for strand
        position_col: Column name for position
        count_col: Column name for read counts (optional, defaults to 1 per row)
        n_workers: Number of parallel workers (default: 1 = no parallelization)

    Returns:
        DataFrame with cluster definitions:
            - cluster_id: Unique cluster identifier
            - chrom, strand: Genomic location
            - start, end: Cluster boundaries
            - modal_position: Median position (best estimate of true CPA)
            - n_positions: Number of unique positions in cluster
            - n_reads: Total read count
    """
    if positions_df.empty:
        return pd.DataFrame(columns=[
            'cluster_id', 'chrom', 'strand', 'start', 'end',
            'modal_position', 'n_positions', 'n_reads'
        ])

    # Add count column if not present
    if count_col is None or count_col not in positions_df.columns:
        positions_df = positions_df.copy()
        positions_df['_count'] = 1
        count_col = '_count'

    # Prepare worker arguments
    worker_args = []
    start_id = 0

    for (chrom, strand), group in positions_df.groupby([chrom_col, strand_col]):
        group = group.sort_values(position_col)
        positions = group[position_col].values
        counts = group[count_col].values

        # Estimate cluster count for ID spacing
        if len(positions) > 0:
            gaps = np.diff(positions)
            estimated_clusters = np.sum(gaps > cluster_distance) + 1
        else:
            estimated_clusters = 0

        worker_args.append((
            chrom, strand, positions, counts,
            cluster_distance, min_reads, start_id
        ))
        start_id += estimated_clusters + 100  # Buffer for ID spacing

    # Process with or without parallelization
    if n_workers > 1 and len(worker_args) > 1:
        with Pool(processes=min(n_workers, len(worker_args))) as pool:
            results = pool.map(_cluster_chrom_strand_group, worker_args)
    else:
        results = [_cluster_chrom_strand_group(args) for args in worker_args]

    # Flatten results and reassign sequential IDs
    all_clusters = []
    cluster_id = 0
    for cluster_list in results:
        for cluster in cluster_list:
            cluster['cluster_id'] = f'cluster_{cluster_id:06d}'
            all_clusters.append(cluster)
            cluster_id += 1

    return pd.DataFrame(all_clusters)


def cluster_cpa_sites_adaptive(
    positions_df: pd.DataFrame,
    max_cluster_radius: int = DEFAULT_MAX_CLUSTER_RADIUS,
    min_peak_separation: int = DEFAULT_MIN_PEAK_SEPARATION,
    min_reads: int = DEFAULT_MIN_READS,
    chrom_col: str = 'chrom',
    strand_col: str = 'strand',
    position_col: str = 'corrected_position',
    count_col: Optional[str] = None,
) -> pd.DataFrame:
    """
    Cluster CPA sites using adaptive valley-based clustering.

    Identifies peaks (local maxima by count) and extends clusters to
    valleys between peaks, capped at max_cluster_radius.

    Args:
        positions_df: DataFrame with CPA positions
        max_cluster_radius: Maximum distance (bp) from peak to boundary
        min_peak_separation: Minimum distance (bp) between distinct peaks
        min_reads: Minimum total reads per cluster
        chrom_col: Column name for chromosome
        strand_col: Column name for strand
        position_col: Column name for position
        count_col: Column name for read counts

    Returns:
        DataFrame with cluster definitions including cluster_width
    """
    if positions_df.empty:
        return pd.DataFrame(columns=[
            'cluster_id', 'chrom', 'strand', 'start', 'end',
            'modal_position', 'n_positions', 'n_reads', 'cluster_width'
        ])

    # Add count column if not present
    if count_col is None or count_col not in positions_df.columns:
        positions_df = positions_df.copy()
        positions_df['_count'] = 1
        count_col = '_count'

    clusters = []
    cluster_id = 0

    # Process each chrom/strand combination
    for (chrom, strand), group in positions_df.groupby([chrom_col, strand_col]):
        # Aggregate counts by position
        pos_counts = group.groupby(position_col)[count_col].sum().sort_index()

        if len(pos_counts) == 0:
            continue

        positions = pos_counts.index.values
        counts = pos_counts.values

        # Identify peaks (positions with highest counts, spaced apart)
        peaks = _identify_peaks(positions, counts, min_peak_separation)

        if not peaks:
            continue

        # Find valleys between peaks
        valleys = _find_valleys_between_peaks(positions, counts, peaks)

        # Get adaptive cluster boundaries
        boundaries = _get_adaptive_cluster_boundaries(
            peaks, valleys, max_cluster_radius
        )

        # Create clusters
        for peak_pos, (left_bound, right_bound) in zip(peaks, boundaries):
            # Get positions and counts within boundaries
            mask = (positions >= left_bound) & (positions <= right_bound)
            cluster_positions = positions[mask]
            cluster_counts = counts[mask]

            total_reads = cluster_counts.sum()
            if total_reads < min_reads:
                continue

            clusters.append({
                'cluster_id': f'cluster_{cluster_id:06d}',
                'chrom': chrom,
                'strand': strand,
                'start': int(left_bound),
                'end': int(right_bound),
                'modal_position': int(peak_pos),
                'n_positions': len(cluster_positions),
                'n_reads': int(total_reads),
                'cluster_width': int(right_bound - left_bound + 1),
            })
            cluster_id += 1

    return pd.DataFrame(clusters)


def _identify_peaks(
    positions: np.ndarray,
    counts: np.ndarray,
    min_separation: int,
) -> List[int]:
    """
    Identify peak positions by count, ensuring minimum separation.

    Iteratively selects highest-count positions that are at least
    min_separation bp from all previously selected peaks.
    """
    if len(positions) == 0:
        return []

    # Sort indices by count (descending)
    sorted_indices = np.argsort(counts)[::-1]

    peaks = []
    for idx in sorted_indices:
        pos = positions[idx]

        # Check distance to existing peaks
        if all(abs(pos - p) >= min_separation for p in peaks):
            peaks.append(int(pos))

    # Sort peaks by position
    peaks.sort()
    return peaks


def _find_valleys_between_peaks(
    positions: np.ndarray,
    counts: np.ndarray,
    peaks: List[int],
) -> List[int]:
    """
    Find valleys (local minima) between adjacent peaks.
    """
    valleys = []

    for i in range(len(peaks) - 1):
        left_peak = peaks[i]
        right_peak = peaks[i + 1]

        # Get positions between peaks
        mask = (positions > left_peak) & (positions < right_peak)
        between_positions = positions[mask]
        between_counts = counts[mask]

        if len(between_positions) == 0:
            # No positions between peaks - use midpoint
            valleys.append((left_peak + right_peak) // 2)
        else:
            # Find minimum count position
            min_idx = np.argmin(between_counts)
            valleys.append(int(between_positions[min_idx]))

    return valleys


def _get_adaptive_cluster_boundaries(
    peaks: List[int],
    valleys: List[int],
    max_radius: int,
) -> List[Tuple[int, int]]:
    """
    Calculate adaptive cluster boundaries for each peak.

    Boundaries extend to valleys or max_radius, whichever is closer.
    """
    boundaries = []

    for i, peak in enumerate(peaks):
        # Left boundary
        if i == 0:
            left_bound = peak - max_radius
        else:
            valley = valleys[i - 1]
            midpoint = (valley + peak) // 2
            left_bound = max(midpoint, peak - max_radius)

        # Right boundary
        if i == len(peaks) - 1:
            right_bound = peak + max_radius
        else:
            valley = valleys[i]
            midpoint = (peak + valley) // 2
            right_bound = min(midpoint, peak + max_radius)

        boundaries.append((left_bound, right_bound))

    return boundaries


def build_cluster_count_matrix(
    positions_df: pd.DataFrame,
    clusters_df: pd.DataFrame,
    sample_col: str = 'sample',
    chrom_col: str = 'chrom',
    strand_col: str = 'strand',
    position_col: str = 'corrected_position',
    count_col: Optional[str] = None,
    fraction_col: Optional[str] = None,
) -> pd.DataFrame:
    """
    Build cluster × sample count matrix with O(log n) interval lookup.

    Uses IntervalTree when available for efficient position-to-cluster mapping.
    Supports fractional counts when positions have a fraction column.

    Args:
        positions_df: DataFrame with CPA positions and sample information
        clusters_df: DataFrame with cluster definitions
        sample_col: Column name for sample identifier
        chrom_col: Column name for chromosome
        strand_col: Column name for strand
        position_col: Column name for position
        count_col: Column name for read counts (optional)
        fraction_col: Column name for fractional weights (optional, for proportional assignment)

    Returns:
        DataFrame with cluster_id as index, samples as columns, counts as values
    """
    if positions_df.empty or clusters_df.empty:
        return pd.DataFrame()

    # Add count column if not present
    if count_col is None or count_col not in positions_df.columns:
        positions_df = positions_df.copy()
        positions_df['_count'] = 1
        count_col = '_count'

    # Build position-to-cluster mapping using IntervalTree if available
    if HAS_INTERVALTREE:
        # O(log n) lookup using IntervalTree
        interval_trees = {}

        for _, cluster in clusters_df.iterrows():
            key = (cluster['chrom'], cluster['strand'])
            if key not in interval_trees:
                interval_trees[key] = IntervalTree()
            # IntervalTree uses half-open intervals [start, end)
            interval_trees[key][cluster['start']:cluster['end']+1] = cluster['cluster_id']

        def find_cluster(chrom, strand, pos):
            key = (chrom, strand)
            if key not in interval_trees:
                return None
            intervals = interval_trees[key][pos]
            if intervals:
                return intervals.pop().data
            return None
    else:
        # Fallback: sorted list with linear search
        position_to_cluster = {}

        for _, cluster in clusters_df.iterrows():
            key = (cluster['chrom'], cluster['strand'])
            if key not in position_to_cluster:
                position_to_cluster[key] = []
            position_to_cluster[key].append(
                (cluster['start'], cluster['end'], cluster['cluster_id'])
            )

        for key in position_to_cluster:
            position_to_cluster[key].sort(key=lambda x: x[0])

        def find_cluster(chrom, strand, pos):
            key = (chrom, strand)
            if key not in position_to_cluster:
                return None
            for start, end, cluster_id in position_to_cluster[key]:
                if start <= pos <= end:
                    return cluster_id
            return None

    # Build count matrix using vectorized groupby where possible
    samples = positions_df[sample_col].unique()

    # Assign cluster IDs to all positions using list comprehension (faster than apply)
    # Extract arrays once to avoid repeated DataFrame indexing
    chroms = positions_df[chrom_col].values
    strands = positions_df[strand_col].values
    positions = positions_df[position_col].values

    # Vectorized cluster assignment using list comprehension
    cluster_ids = [
        find_cluster(c, s, p)
        for c, s, p in zip(chroms, strands, positions)
    ]

    positions_df = positions_df.copy()
    positions_df['_cluster_id'] = cluster_ids

    # Filter to positions with assigned clusters
    mask = positions_df['_cluster_id'].notna()
    assigned = positions_df.loc[mask].copy()

    if assigned.empty:
        # No positions assigned to clusters
        count_matrix = pd.DataFrame(
            0, index=clusters_df['cluster_id'], columns=sorted(samples)
        )
        count_matrix.index.name = 'cluster_id'
        return count_matrix

    # Calculate effective counts (with fraction weighting if applicable)
    if fraction_col and fraction_col in assigned.columns:
        assigned['_effective_count'] = assigned[count_col] * assigned[fraction_col]
    else:
        assigned['_effective_count'] = assigned[count_col]

    # Aggregate using pandas groupby (much faster than iterrows)
    count_matrix = assigned.pivot_table(
        index='_cluster_id',
        columns=sample_col,
        values='_effective_count',
        aggfunc='sum',
        fill_value=0
    )
    count_matrix.index.name = 'cluster_id'

    # Ensure all samples are present
    for sample in samples:
        if sample not in count_matrix.columns:
            count_matrix[sample] = 0

    # Ensure all clusters are present
    for cluster_id in clusters_df['cluster_id']:
        if cluster_id not in count_matrix.index:
            count_matrix.loc[cluster_id] = 0

    return count_matrix[sorted(count_matrix.columns)]


def annotate_clusters_with_genes(
    clusters_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    upstream_window: int = 500,
    downstream_window: int = 100,
) -> pd.DataFrame:
    """
    Annotate clusters with nearest gene information.

    Args:
        clusters_df: DataFrame with cluster definitions
        annotation_df: DataFrame with gene annotations (must have:
            chrom, strand, start, end, gene_id, gene_name columns)
        upstream_window: Maximum distance upstream to search for genes
        downstream_window: Maximum distance downstream to search

    Returns:
        clusters_df with additional columns:
            gene_id, gene_name, distance_to_gene_3prime
    """
    if clusters_df.empty:
        return clusters_df

    # Add annotation columns
    clusters_df = clusters_df.copy()
    clusters_df['gene_id'] = None
    clusters_df['gene_name'] = None
    clusters_df['distance_to_gene_3prime'] = None

    # Build gene lookup by chrom/strand
    gene_lookup = defaultdict(list)
    for _, gene in annotation_df.iterrows():
        key = (gene['chrom'], gene['strand'])
        # Gene 3' end position
        if gene['strand'] == '+':
            gene_3prime = gene['end']
        else:
            gene_3prime = gene['start']
        gene_lookup[key].append({
            'gene_id': gene['gene_id'],
            'gene_name': gene.get('gene_name', gene['gene_id']),
            'gene_3prime': gene_3prime,
        })

    # Sort by 3' position for efficient lookup
    for key in gene_lookup:
        gene_lookup[key].sort(key=lambda x: x['gene_3prime'])

    # Annotate each cluster
    for idx, cluster in clusters_df.iterrows():
        key = (cluster['chrom'], cluster['strand'])
        if key not in gene_lookup:
            continue

        cluster_pos = cluster['modal_position']
        strand = cluster['strand']

        best_gene = None
        best_distance = float('inf')

        for gene in gene_lookup[key]:
            gene_3prime = gene['gene_3prime']

            # Calculate strand-aware distance
            if strand == '+':
                # Cluster should be downstream of gene 3' end
                distance = cluster_pos - gene_3prime
            else:
                # Cluster should be upstream (lower coordinate) of gene 3' end
                distance = gene_3prime - cluster_pos

            # Check if within window
            if -upstream_window <= distance <= downstream_window:
                if abs(distance) < abs(best_distance):
                    best_distance = distance
                    best_gene = gene

        if best_gene:
            clusters_df.at[idx, 'gene_id'] = best_gene['gene_id']
            clusters_df.at[idx, 'gene_name'] = best_gene['gene_name']
            clusters_df.at[idx, 'distance_to_gene_3prime'] = int(best_distance)

    return clusters_df
