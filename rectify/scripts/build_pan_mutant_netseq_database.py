#!/usr/bin/env python3
"""
Build Pan-Mutant NET-seq Database for A-Tract Refinement

This script aggregates NET-seq signal from all GSE159603 (Couvillion 2022) mutant
samples to create a comprehensive database of NET-seq signal distributions at
A-tract CPA sites. This database is used to refine nanopore 3' end positions
that fall within ambiguous A-tracts.

Key insight: While individual NET-seq samples may have sparse coverage at any
given CPA site, aggregating across 109 mutant samples provides robust signal
distributions that can resolve A-tract ambiguity.

Usage:
    python build_pan_mutant_netseq_database.py \
        --cpa-clusters /path/to/cpa_clusters.tsv \
        --bam-dir /path/to/aligned_hexamer_trimmed \
        --output /path/to/pan_mutant_netseq_database.h5 \
        --min-acount 4 \
        --threads 16

Output:
    HDF5 file containing:
    - /signals/{cluster_id}: NET-seq signal array for each A-tract CPA site
    - /metadata: Cluster metadata (chrom, strand, tract_start, tract_end, etc.)
    - /aggregated_psf/{acount}: Empirical PSF by A-count category

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-24
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings

import numpy as np
import pandas as pd
import pysam

# Suppress pysam truncation warnings
warnings.filterwarnings('ignore', category=UserWarning, module='pysam')


# =============================================================================
# Constants
# =============================================================================

# Chromosome name mapping (NCBI -> canonical)
NCBI_TO_CHROM = {
    'ref|NC_001133|': 'chrI',
    'ref|NC_001134|': 'chrII',
    'ref|NC_001135|': 'chrIII',
    'ref|NC_001136|': 'chrIV',
    'ref|NC_001137|': 'chrV',
    'ref|NC_001138|': 'chrVI',
    'ref|NC_001139|': 'chrVII',
    'ref|NC_001140|': 'chrVIII',
    'ref|NC_001141|': 'chrIX',
    'ref|NC_001142|': 'chrX',
    'ref|NC_001143|': 'chrXI',
    'ref|NC_001144|': 'chrXII',
    'ref|NC_001145|': 'chrXIII',
    'ref|NC_001146|': 'chrXIV',
    'ref|NC_001147|': 'chrXV',
    'ref|NC_001148|': 'chrXVI',
    'ref|NC_001224|': 'chrMito',
}

CHROM_TO_NCBI = {v: k for k, v in NCBI_TO_CHROM.items()}

# Window extension around A-tract for signal extraction
WINDOW_EXTENSION = 10  # bp on each side

# Minimum reads to include a site in the database
MIN_READS_THRESHOLD = 10


# =============================================================================
# Helper Functions
# =============================================================================

def get_bam_chrom_map(bam_path: str) -> Dict[str, str]:
    """
    Build chromosome name mapping from BAM file to our canonical names.

    Handles various chromosome naming conventions:
    - Short names: chrI, chrII, etc.
    - NCBI format: ref|NC_001133|
    - Long headers with NC_* embedded
    """
    chrom_map = {}

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for ref_name in bam.references:
            # Check for NCBI format
            if ref_name in NCBI_TO_CHROM:
                chrom_map[NCBI_TO_CHROM[ref_name]] = ref_name
            # Check for canonical format
            elif ref_name in CHROM_TO_NCBI:
                chrom_map[ref_name] = ref_name
            # Check for NC_* substring (long headers)
            else:
                for ncbi, canonical in NCBI_TO_CHROM.items():
                    nc_id = ncbi.replace('ref|', '').replace('|', '')
                    if nc_id in ref_name:
                        chrom_map[canonical] = ref_name
                        break

    return chrom_map


def extract_atract_from_sequence(
    sequence: str,
    position: int,
    strand: str,
    window: int = 20
) -> Tuple[int, int, int]:
    """
    Extract A-tract boundaries from genomic sequence.

    For minus strand CPA sites, we look UPSTREAM (higher coords) for T's
    (which are A's on the transcript).

    Returns:
        Tuple of (tract_start, tract_end, tract_length)
    """
    if strand == '-':
        # Minus strand: look upstream (higher coords) for T's
        tract_end = position
        tract_start = position
        for i in range(window):
            upstream_pos = position + i + 1
            if upstream_pos >= len(sequence):
                break
            if sequence[upstream_pos].upper() == 'T':
                tract_start = upstream_pos
            else:
                break
        tract_length = tract_start - tract_end
    else:
        # Plus strand: look upstream (lower coords) for A's
        tract_start = position
        tract_end = position
        for i in range(window):
            upstream_pos = position - i - 1
            if upstream_pos < 0:
                break
            if sequence[upstream_pos].upper() == 'A':
                tract_start = upstream_pos
            else:
                break
        tract_length = tract_end - tract_start

    return tract_start, tract_end, tract_length


def extract_netseq_signal_at_site(
    bam: pysam.AlignmentFile,
    chrom: str,
    position: int,
    strand: str,
    window: int = WINDOW_EXTENSION,
) -> np.ndarray:
    """
    Extract NET-seq 3' end signal at a CPA site.

    For NET-seq, we count the 3' end positions of reads (accounting for strand).
    Minus strand reads have 3' end at alignment start.
    Plus strand reads have 3' end at alignment end.

    Args:
        bam: Open pysam AlignmentFile
        chrom: Chromosome name (in BAM format)
        position: CPA position
        strand: Gene strand
        window: Window size on each side

    Returns:
        Signal array of length 2*window + 1, centered on position
    """
    start = position - window
    end = position + window + 1

    signal = np.zeros(2 * window + 1)

    try:
        for read in bam.fetch(chrom, max(0, start - 100), end + 100):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            # Determine read strand
            read_strand = '-' if read.is_reverse else '+'

            # We want reads on the same strand as the gene
            # (NET-seq is strand-specific)
            if read_strand != strand:
                continue

            # Get 3' end position
            if read_strand == '-':
                three_prime = read.reference_start  # 0-based
            else:
                three_prime = read.reference_end - 1  # 0-based

            # Check if within window
            if start <= three_prime <= end - 1:
                idx = three_prime - start
                signal[idx] += 1

    except ValueError:
        # Region out of bounds
        pass

    return signal


def process_cluster(
    cluster_row: pd.Series,
    bam_paths: List[str],
    chrom_maps: Dict[str, Dict[str, str]],
    reference_fasta: Optional[str] = None,
) -> Optional[Dict]:
    """
    Process a single CPA cluster to aggregate NET-seq signal across all BAMs.

    Args:
        cluster_row: Row from CPA clusters DataFrame
        bam_paths: List of BAM file paths
        chrom_maps: Precomputed chromosome maps for each BAM
        reference_fasta: Path to reference genome (for A-tract detection)

    Returns:
        Dict with cluster_id, signal array, and metadata
    """
    cluster_id = cluster_row['cluster_id']
    chrom = cluster_row['chrom']
    strand = cluster_row['strand']
    position = int(cluster_row['modal_position'])

    # Convert chromosome name to canonical format if needed
    # Handle NCBI format (ref|NC_001133|) -> chrI
    canonical_chrom = NCBI_TO_CHROM.get(chrom, chrom)

    # Skip mitochondrial chromosomes (not in most NET-seq BAMs)
    if canonical_chrom in ['chrMito', 'chrmt', 'chrM']:
        return None

    # Aggregate signal across all BAMs
    total_signal = None
    n_contributing_bams = 0

    for bam_path in bam_paths:
        bam_name = Path(bam_path).stem

        # Get BAM-specific chromosome name
        chrom_map = chrom_maps.get(bam_path, {})
        bam_chrom = chrom_map.get(canonical_chrom)

        if bam_chrom is None:
            continue

        try:
            with pysam.AlignmentFile(bam_path, 'rb') as bam:
                signal = extract_netseq_signal_at_site(
                    bam, bam_chrom, position, strand
                )

                if signal.sum() > 0:
                    n_contributing_bams += 1
                    if total_signal is None:
                        total_signal = signal.copy()
                    else:
                        total_signal += signal

        except Exception as e:
            # Skip problematic BAMs
            continue

    if total_signal is None or total_signal.sum() < MIN_READS_THRESHOLD:
        return None

    return {
        'cluster_id': cluster_id,
        'chrom': canonical_chrom,
        'strand': strand,
        'position': position,
        'signal': total_signal,
        'total_reads': int(total_signal.sum()),
        'n_contributing_bams': n_contributing_bams,
    }


def compute_empirical_psf(
    signals: List[np.ndarray],
    center_idx: int,
) -> np.ndarray:
    """
    Compute empirical PSF from aggregated signals.

    The PSF is computed by normalizing each signal to its peak
    and then averaging across sites.

    Args:
        signals: List of signal arrays (all same length)
        center_idx: Index of the center position

    Returns:
        Normalized PSF array
    """
    if not signals:
        return np.array([])

    # Stack signals
    signal_matrix = np.vstack(signals)

    # Normalize each row by its sum
    row_sums = signal_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    normalized = signal_matrix / row_sums

    # Average across sites
    psf = normalized.mean(axis=0)

    return psf


# =============================================================================
# Main Processing
# =============================================================================

def build_database(
    cpa_clusters_path: str,
    bam_dir: str,
    output_path: str,
    min_acount: int = 4,
    max_clusters: Optional[int] = None,
    threads: int = 8,
    verbose: bool = True,
):
    """
    Build pan-mutant NET-seq database for A-tract refinement.

    Args:
        cpa_clusters_path: Path to CPA clusters TSV
        bam_dir: Directory containing hexamer-trimmed NET-seq BAMs
        output_path: Output HDF5 file path
        min_acount: Minimum downstream A-count to include
        max_clusters: Maximum clusters to process (for testing)
        threads: Number of parallel threads
        verbose: Print progress
    """
    import h5py

    if verbose:
        print(f"Loading CPA clusters from {cpa_clusters_path}")

    # Load CPA clusters
    clusters_df = pd.read_csv(cpa_clusters_path, sep='\t')

    if verbose:
        print(f"  Total clusters: {len(clusters_df):,}")

    # Filter out mitochondrial chromosomes first (not in most NET-seq BAMs)
    mito_chroms = ['chrmt', 'chrMito', 'chrM', 'ref|NC_001224|']
    clusters_df = clusters_df[~clusters_df['chrom'].isin(mito_chroms)].copy()
    if verbose:
        print(f"  Nuclear clusters: {len(clusters_df):,}")

    # Filter to clusters with A-tracts if downstream_base column exists
    # Otherwise process all clusters
    if 'downstream_base' in clusters_df.columns:
        # For A-tract refinement, we want sites WITH downstream A's
        # (these are the ambiguous sites that need refinement)
        atract_df = clusters_df[clusters_df['downstream_base'] == 'A'].copy()
        if verbose:
            print(f"  Clusters with downstream A (A-tracts): {len(atract_df):,}")
        # Also include non-A sites for building complete PSF
        non_atract_df = clusters_df[clusters_df['downstream_base'] != 'A'].copy()
        if verbose:
            print(f"  Clusters without downstream A (0A sites): {len(non_atract_df):,}")
    else:
        # No downstream_base column - process all clusters
        atract_df = clusters_df.copy()
        if verbose:
            print(f"  Processing all clusters (no downstream_base column)")

    if max_clusters:
        atract_df = atract_df.head(max_clusters)
        if verbose:
            print(f"  Limited to {max_clusters} clusters for testing")

    # Find all BAM files
    bam_dir_path = Path(bam_dir)
    bam_paths = sorted(bam_dir_path.glob('*.bam'))

    if verbose:
        print(f"\nFound {len(bam_paths)} BAM files in {bam_dir}")

    if not bam_paths:
        print("ERROR: No BAM files found!")
        return

    # Pre-compute chromosome maps for all BAMs
    if verbose:
        print("Building chromosome maps for all BAMs...")

    chrom_maps = {}
    for bam_path in bam_paths:
        chrom_maps[str(bam_path)] = get_bam_chrom_map(str(bam_path))

    # Process clusters
    if verbose:
        print(f"\nProcessing {len(atract_df):,} A-tract clusters...")

    results = []
    bam_path_strs = [str(p) for p in bam_paths]

    # Process in batches to manage memory
    batch_size = 1000
    n_batches = (len(atract_df) + batch_size - 1) // batch_size

    for batch_idx in range(n_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, len(atract_df))
        batch_df = atract_df.iloc[start_idx:end_idx]

        if verbose and batch_idx % 10 == 0:
            print(f"  Processing batch {batch_idx + 1}/{n_batches} "
                  f"(clusters {start_idx:,}-{end_idx:,})")

        # Process each cluster in the batch
        for _, row in batch_df.iterrows():
            result = process_cluster(
                row, bam_path_strs, chrom_maps
            )
            if result:
                results.append(result)

    if verbose:
        print(f"\n  Successfully processed {len(results):,} clusters with signal")

    # Save to HDF5
    if verbose:
        print(f"\nSaving database to {output_path}")

    with h5py.File(output_path, 'w') as h5f:
        # Create groups
        signals_grp = h5f.create_group('signals')
        metadata_grp = h5f.create_group('metadata')

        # Store signals
        cluster_ids = []
        chroms = []
        strands = []
        positions = []
        total_reads_list = []
        n_bams_list = []

        for result in results:
            cluster_id = str(result['cluster_id'])
            signals_grp.create_dataset(cluster_id, data=result['signal'])

            cluster_ids.append(result['cluster_id'])
            chroms.append(result['chrom'])
            strands.append(result['strand'])
            positions.append(result['position'])
            total_reads_list.append(result['total_reads'])
            n_bams_list.append(result['n_contributing_bams'])

        # Store metadata
        metadata_grp.create_dataset('cluster_id', data=np.array(cluster_ids))
        metadata_grp.create_dataset('position', data=np.array(positions))
        metadata_grp.create_dataset('total_reads', data=np.array(total_reads_list))
        metadata_grp.create_dataset('n_contributing_bams', data=np.array(n_bams_list))

        # Store chromosome and strand as string datasets
        dt = h5py.string_dtype(encoding='utf-8')
        metadata_grp.create_dataset('chrom', data=np.array(chroms, dtype=object), dtype=dt)
        metadata_grp.create_dataset('strand', data=np.array(strands, dtype=object), dtype=dt)

        # Compute and store aggregated PSF
        psf_grp = h5f.create_group('aggregated_psf')
        all_signals = [r['signal'] for r in results]

        if all_signals:
            center_idx = len(all_signals[0]) // 2
            empirical_psf = compute_empirical_psf(all_signals, center_idx)
            psf_grp.create_dataset('all_sites', data=empirical_psf)

        # Store metadata
        h5f.attrs['n_clusters'] = len(results)
        h5f.attrs['n_bam_files'] = len(bam_paths)
        h5f.attrs['window_extension'] = WINDOW_EXTENSION
        h5f.attrs['min_reads_threshold'] = MIN_READS_THRESHOLD
        h5f.attrs['creation_date'] = pd.Timestamp.now().isoformat()

    if verbose:
        print(f"\nDatabase created successfully!")
        print(f"  Clusters with signal: {len(results):,}")
        print(f"  Total reads: {sum(total_reads_list):,}")
        print(f"  Mean reads per cluster: {np.mean(total_reads_list):.1f}")
        print(f"  Median contributing BAMs: {np.median(n_bams_list):.0f}")


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Build pan-mutant NET-seq database for A-tract refinement',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        '--cpa-clusters', '-c',
        required=True,
        help='Path to CPA clusters TSV file'
    )
    parser.add_argument(
        '--bam-dir', '-b',
        required=True,
        help='Directory containing hexamer-trimmed NET-seq BAM files'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output HDF5 file path'
    )
    parser.add_argument(
        '--min-acount', '-m',
        type=int,
        default=4,
        help='Minimum downstream A-count to include (default: 4)'
    )
    parser.add_argument(
        '--max-clusters',
        type=int,
        default=None,
        help='Maximum clusters to process (for testing)'
    )
    parser.add_argument(
        '--threads', '-t',
        type=int,
        default=8,
        help='Number of parallel threads (default: 8)'
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress progress output'
    )

    args = parser.parse_args()

    build_database(
        cpa_clusters_path=args.cpa_clusters,
        bam_dir=args.bam_dir,
        output_path=args.output,
        min_acount=args.min_acount,
        max_clusters=args.max_clusters,
        threads=args.threads,
        verbose=not args.quiet,
    )


if __name__ == '__main__':
    main()
