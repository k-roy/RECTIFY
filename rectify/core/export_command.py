#!/usr/bin/env python3
"""
Export corrected 3' end positions to bedGraph/bigWig format.

Creates two types of output:
1. Per-replicate: One file per replicate (e.g., wt_by4742_rep1.plus.bw)
2. Per-condition summed: All replicates summed (e.g., wt_by4742.plus.bw)

Author: Kevin R. Roy
Date: 2026-03-19
"""

import argparse
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd

try:
    import pyBigWig
    HAS_PYBIGWIG = True
except ImportError:
    HAS_PYBIGWIG = False


def load_corrected_3prime(input_path: Path, position_col: str = 'position') -> pd.DataFrame:
    """Load corrected 3' end file."""
    print(f"Loading {input_path}...")
    df = pd.read_csv(input_path, sep='\t')
    print(f"  Loaded {len(df):,} reads")
    print(f"  Columns: {list(df.columns)}")

    # Check required columns
    required = ['chrom', position_col, 'strand']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    return df


def aggregate_positions(
    df: pd.DataFrame,
    position_col: str = 'position',
    group_cols: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Aggregate read counts per position.

    Args:
        df: DataFrame with corrected 3' ends
        position_col: Column with position
        group_cols: Additional columns to group by (e.g., ['replicate'])

    Returns:
        DataFrame with chrom, position, strand, count (and group_cols if specified)
    """
    if group_cols is None:
        group_cols = []

    group_by = ['chrom', position_col, 'strand'] + group_cols
    counts = df.groupby(group_by).size().reset_index(name='count')
    counts.rename(columns={position_col: 'position'}, inplace=True)

    return counts


def write_bedgraph(
    counts_df: pd.DataFrame,
    output_path: Path,
    strand: str,
    chrom_sizes: Optional[Dict[str, int]] = None
) -> int:
    """Write strand-specific bedGraph file."""
    strand_df = counts_df[counts_df['strand'] == strand].copy()
    strand_df = strand_df.sort_values(['chrom', 'position'])

    with open(output_path, 'w') as f:
        for _, row in strand_df.iterrows():
            chrom = row['chrom']
            pos = int(row['position'])
            count = row['count']
            # bedGraph: chrom, start (0-based), end, value
            f.write(f"{chrom}\t{pos}\t{pos+1}\t{count}\n")

    print(f"  Wrote {len(strand_df):,} positions to {output_path}")
    return len(strand_df)


def write_bigwig(
    counts_df: pd.DataFrame,
    output_path: Path,
    strand: str,
    chrom_sizes: Dict[str, int]
) -> int:
    """Write strand-specific bigWig file."""
    if not HAS_PYBIGWIG:
        print(f"  Skipping bigWig (pyBigWig not available)")
        return 0

    strand_df = counts_df[counts_df['strand'] == strand].copy()
    strand_df = strand_df.sort_values(['chrom', 'position'])

    # Filter to known chromosomes
    strand_df = strand_df[strand_df['chrom'].isin(chrom_sizes)]

    if len(strand_df) == 0:
        print(f"  No data for strand {strand}")
        return 0

    bw = pyBigWig.open(str(output_path), "w")

    # Add header with chromosome sizes
    header = [(chrom, size) for chrom, size in chrom_sizes.items()]
    bw.addHeader(header)

    # Add entries by chromosome
    for chrom in strand_df['chrom'].unique():
        chrom_df = strand_df[strand_df['chrom'] == chrom]
        if len(chrom_df) == 0:
            continue

        starts = chrom_df['position'].values.astype(np.int64)
        ends = (chrom_df['position'].values + 1).astype(np.int64)
        values = chrom_df['count'].values.astype(np.float64)

        # Ensure positions are within chromosome bounds
        valid = (starts >= 0) & (ends <= chrom_sizes[chrom])
        if not valid.all():
            starts = starts[valid]
            ends = ends[valid]
            values = values[valid]

        if len(starts) > 0:
            bw.addEntries([chrom] * len(starts), starts, ends=ends, values=values)

    bw.close()
    print(f"  Wrote {len(strand_df):,} positions to {output_path}")
    return len(strand_df)


def get_chrom_sizes_from_genome(genome_path: Path) -> Dict[str, int]:
    """Extract chromosome sizes from genome FASTA."""
    from Bio import SeqIO

    chrom_sizes = {}
    for record in SeqIO.parse(genome_path, "fasta"):
        chrom_sizes[record.id] = len(record.seq)

    return chrom_sizes


def process_per_replicate(
    df: pd.DataFrame,
    output_dir: Path,
    output_format: str,
    position_col: str,
    chrom_sizes: Dict[str, int]
) -> None:
    """Generate per-replicate output files."""
    print("\n=== Per-Replicate Output ===")

    if 'replicate' not in df.columns:
        print("  Warning: No 'replicate' column found. Skipping per-replicate output.")
        return

    # Get unique replicates
    replicates = df['replicate'].unique()
    print(f"Found {len(replicates)} replicates")

    for rep in sorted(replicates):
        rep_df = df[df['replicate'] == rep]
        counts = aggregate_positions(rep_df, position_col=position_col)

        # Clean replicate name for filename
        rep_clean = str(rep).replace('/', '_').replace(' ', '_')

        for strand, suffix in [('+', 'plus'), ('-', 'minus')]:
            if output_format == 'bigwig':
                out_path = output_dir / f"{rep_clean}.{suffix}.bw"
                write_bigwig(counts, out_path, strand, chrom_sizes)
            else:
                out_path = output_dir / f"{rep_clean}.{suffix}.bedgraph"
                write_bedgraph(counts, out_path, strand, chrom_sizes)


def process_per_condition(
    df: pd.DataFrame,
    output_dir: Path,
    output_format: str,
    position_col: str,
    chrom_sizes: Dict[str, int]
) -> None:
    """Generate per-condition summed output files."""
    print("\n=== Per-Condition Summed Output ===")

    if 'condition' not in df.columns:
        print("  Warning: No 'condition' column found. Using all data as single condition.")
        conditions = ['all']
        df = df.copy()
        df['condition'] = 'all'
    else:
        conditions = df['condition'].unique()

    print(f"Found {len(conditions)} conditions")

    for cond in sorted(conditions):
        cond_df = df[df['condition'] == cond]
        # Aggregate across all replicates for this condition
        counts = aggregate_positions(cond_df, position_col=position_col)

        # Clean condition name for filename
        cond_clean = str(cond).replace('/', '_').replace(' ', '_')

        print(f"\n{cond}: {len(cond_df):,} reads -> {len(counts):,} positions")

        for strand, suffix in [('+', 'plus'), ('-', 'minus')]:
            if output_format == 'bigwig':
                out_path = output_dir / f"{cond_clean}.{suffix}.bw"
                write_bigwig(counts, out_path, strand, chrom_sizes)
            else:
                out_path = output_dir / f"{cond_clean}.{suffix}.bedgraph"
                write_bedgraph(counts, out_path, strand, chrom_sizes)


def run(args: argparse.Namespace) -> int:
    """Run the export command."""
    print("RECTIFY Export: Generate bedGraph/bigWig from corrected 3' ends")
    print("=" * 60)

    # Validate inputs
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}")
        return 1

    # Get chromosome sizes
    if args.genome:
        print(f"Loading chromosome sizes from {args.genome}...")
        chrom_sizes = get_chrom_sizes_from_genome(args.genome)
    elif args.chrom_sizes:
        print(f"Loading chromosome sizes from {args.chrom_sizes}...")
        chrom_sizes = {}
        with open(args.chrom_sizes) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chrom_sizes[parts[0]] = int(parts[1])
    else:
        # Use bundled S. cerevisiae sizes
        print("Using bundled S. cerevisiae chromosome sizes (NCBI format)")
        chrom_sizes = {
            'ref|NC_001133|': 230218, 'ref|NC_001134|': 813184, 'ref|NC_001135|': 316620,
            'ref|NC_001136|': 1531933, 'ref|NC_001137|': 576874, 'ref|NC_001138|': 270161,
            'ref|NC_001139|': 1090940, 'ref|NC_001140|': 562643, 'ref|NC_001141|': 439888,
            'ref|NC_001142|': 745751, 'ref|NC_001143|': 666816, 'ref|NC_001144|': 1078177,
            'ref|NC_001145|': 924431, 'ref|NC_001146|': 784333, 'ref|NC_001147|': 1091291,
            'ref|NC_001148|': 948066, 'ref|NC_001224|': 85779,
        }

    print(f"  {len(chrom_sizes)} chromosomes")

    # Default: both types
    if not args.per_replicate and not args.per_condition:
        args.per_replicate = True
        args.per_condition = True

    # Create output directories
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    df = load_corrected_3prime(args.input, position_col=args.position_col)

    # Check format
    if args.format == 'bigwig' and not HAS_PYBIGWIG:
        print("Warning: pyBigWig not available. Falling back to bedGraph format.")
        args.format = 'bedgraph'

    # Process
    if args.per_replicate:
        rep_dir = args.output_dir / "per_replicate"
        rep_dir.mkdir(exist_ok=True)
        process_per_replicate(df, rep_dir, args.format, args.position_col, chrom_sizes)

    if args.per_condition:
        cond_dir = args.output_dir / "per_condition"
        cond_dir.mkdir(exist_ok=True)
        process_per_condition(df, cond_dir, args.format, args.position_col, chrom_sizes)

    print("\nDone!")
    return 0
