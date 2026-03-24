#!/usr/bin/env python3
"""
Create Compact A-Tract NET-seq Reference for Package Distribution

This script creates a compact reference file containing only NET-seq signal
at A-tract CPA sites for bundling with pip/conda distributions.

The output is a compressed TSV file with:
- CPA positions that have downstream A-tracts (A-count >= 1)
- NET-seq signal distribution within a ±10bp window
- Much smaller than full genome-wide NET-seq data

Usage:
    python create_atract_netseq_reference.py \
        --cpa-clusters /path/to/cpa_clusters.tsv \
        --netseq-pan /path/to/saccharomyces_cerevisiae_netseq_pan.tsv.gz \
        --reference-fasta /path/to/genome.fasta \
        --output /path/to/atract_netseq_reference.tsv.gz

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-24
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import gzip

import numpy as np
import pandas as pd


# Chromosome mappings
NCBI_TO_CHROM = {
    'ref|NC_001133|': 'chrI', 'ref|NC_001134|': 'chrII', 'ref|NC_001135|': 'chrIII',
    'ref|NC_001136|': 'chrIV', 'ref|NC_001137|': 'chrV', 'ref|NC_001138|': 'chrVI',
    'ref|NC_001139|': 'chrVII', 'ref|NC_001140|': 'chrVIII', 'ref|NC_001141|': 'chrIX',
    'ref|NC_001142|': 'chrX', 'ref|NC_001143|': 'chrXI', 'ref|NC_001144|': 'chrXII',
    'ref|NC_001145|': 'chrXIII', 'ref|NC_001146|': 'chrXIV', 'ref|NC_001147|': 'chrXV',
    'ref|NC_001148|': 'chrXVI', 'ref|NC_001224|': 'chrMito',
}

WINDOW = 10  # Window size for signal extraction


def load_reference_genome(fasta_path: str) -> Dict[str, str]:
    """Load reference genome from FASTA file."""
    from Bio import SeqIO

    sequences = {}
    with gzip.open(fasta_path, 'rt') if fasta_path.endswith('.gz') else open(fasta_path) as f:
        for record in SeqIO.parse(f, 'fasta'):
            chrom = record.id
            # Convert to canonical name if needed
            if chrom in NCBI_TO_CHROM:
                chrom = NCBI_TO_CHROM[chrom]
            sequences[chrom] = str(record.seq).upper()

    return sequences


def count_downstream_as(
    sequences: Dict[str, str],
    chrom: str,
    position: int,
    strand: str,
    window: int = 10,
) -> int:
    """Count consecutive A's downstream of position."""
    seq = sequences.get(chrom, '')
    if not seq:
        return 0

    count = 0
    if strand == '+':
        # Plus strand: look downstream (higher coords)
        for i in range(1, window + 1):
            pos = position + i
            if pos >= len(seq):
                break
            if seq[pos] == 'A':
                count += 1
            else:
                break
    else:
        # Minus strand: look downstream (lower coords), check for T
        for i in range(1, window + 1):
            pos = position - i
            if pos < 0:
                break
            if seq[pos] == 'T':
                count += 1
            else:
                break

    return count


def create_atract_reference(
    cpa_clusters_path: str,
    netseq_pan_path: str,
    reference_fasta_path: str,
    output_path: str,
    min_acount: int = 1,
    min_signal: float = 10.0,
    verbose: bool = True,
):
    """
    Create compact A-tract NET-seq reference file.

    Args:
        cpa_clusters_path: Path to CPA clusters TSV
        netseq_pan_path: Path to pan-mutant NET-seq TSV
        reference_fasta_path: Path to reference genome FASTA
        output_path: Output path for compressed TSV
        min_acount: Minimum A-count to include (default: 1)
        min_signal: Minimum total signal to include (default: 10)
        verbose: Print progress
    """
    if verbose:
        print("Loading reference genome...")
    sequences = load_reference_genome(reference_fasta_path)

    if verbose:
        print(f"Loading CPA clusters from {cpa_clusters_path}...")
    clusters_df = pd.read_csv(cpa_clusters_path, sep='\t')

    # Convert chromosome names
    clusters_df['chrom_canonical'] = clusters_df['chrom'].map(
        lambda x: NCBI_TO_CHROM.get(x, x)
    )

    # Filter out mitochondrial
    mito_chroms = ['chrMito', 'chrmt', 'chrM']
    clusters_df = clusters_df[~clusters_df['chrom_canonical'].isin(mito_chroms)]

    if verbose:
        print(f"  Nuclear clusters: {len(clusters_df):,}")

    # Count downstream A's for each cluster
    if verbose:
        print("Counting downstream A's...")

    acounts = []
    for _, row in clusters_df.iterrows():
        acount = count_downstream_as(
            sequences,
            row['chrom_canonical'],
            int(row['modal_position']),
            row['strand'],
        )
        acounts.append(acount)

    clusters_df['downstream_acount'] = acounts

    # Filter to A-tract sites
    atract_df = clusters_df[clusters_df['downstream_acount'] >= min_acount].copy()

    if verbose:
        print(f"  A-tract clusters (acount >= {min_acount}): {len(atract_df):,}")

    # Load NET-seq pan signal
    if verbose:
        print(f"Loading NET-seq pan signal from {netseq_pan_path}...")

    netseq_df = pd.read_csv(netseq_pan_path, sep='\t')

    # Build position -> signal lookup
    if verbose:
        print("Building signal lookup...")

    signal_lookup = {}
    for _, row in netseq_df.iterrows():
        key = (row['chrom'], row['position'], row['strand'])
        signal_lookup[key] = row['signal']

    # Extract signal for each A-tract site
    if verbose:
        print("Extracting signal at A-tract sites...")

    results = []
    for _, row in atract_df.iterrows():
        chrom = row['chrom_canonical']
        pos = int(row['modal_position'])
        strand = row['strand']
        acount = row['downstream_acount']

        # Extract signal in window
        signals = []
        for offset in range(-WINDOW, WINDOW + 1):
            key = (chrom, pos + offset, strand)
            sig = signal_lookup.get(key, 0.0)
            signals.append(sig)

        total_signal = sum(signals)

        if total_signal >= min_signal:
            # Format signals as comma-separated string
            signal_str = ','.join(f'{s:.1f}' for s in signals)

            results.append({
                'chrom': chrom,
                'position': pos,
                'strand': strand,
                'downstream_acount': acount,
                'total_signal': total_signal,
                'signal_window': signal_str,
            })

    if verbose:
        print(f"  Sites with signal >= {min_signal}: {len(results):,}")

    # Create output dataframe
    output_df = pd.DataFrame(results)

    # Sort by total signal (descending)
    output_df = output_df.sort_values('total_signal', ascending=False)

    # Save to compressed TSV
    if verbose:
        print(f"Saving to {output_path}...")

    output_df.to_csv(output_path, sep='\t', index=False, compression='gzip')

    # Report file size
    import os
    size_mb = os.path.getsize(output_path) / 1e6

    if verbose:
        print(f"\nA-tract NET-seq reference created successfully!")
        print(f"  Sites: {len(output_df):,}")
        print(f"  File size: {size_mb:.2f} MB")
        print(f"  Output: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Create compact A-tract NET-seq reference for package distribution',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        '--cpa-clusters', '-c',
        required=True,
        help='Path to CPA clusters TSV file'
    )
    parser.add_argument(
        '--netseq-pan', '-n',
        required=True,
        help='Path to pan-mutant NET-seq TSV file'
    )
    parser.add_argument(
        '--reference-fasta', '-r',
        required=True,
        help='Path to reference genome FASTA'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output path for compressed TSV'
    )
    parser.add_argument(
        '--min-acount',
        type=int,
        default=1,
        help='Minimum downstream A-count to include (default: 1)'
    )
    parser.add_argument(
        '--min-signal',
        type=float,
        default=10.0,
        help='Minimum total signal to include (default: 10)'
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress progress output'
    )

    args = parser.parse_args()

    create_atract_reference(
        cpa_clusters_path=args.cpa_clusters,
        netseq_pan_path=args.netseq_pan,
        reference_fasta_path=args.reference_fasta,
        output_path=args.output,
        min_acount=args.min_acount,
        min_signal=args.min_signal,
        verbose=not args.quiet,
    )


if __name__ == '__main__':
    main()
