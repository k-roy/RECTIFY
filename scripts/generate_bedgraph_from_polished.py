#!/usr/bin/env python3
"""
Generate bedgraph files from RECTIFY polished 3' end TSV files.

Creates strand-specific bedgraph files with all replicates combined per sample.
Output: one plus and one minus bedgraph per sample (RPM-normalized).

This script handles both old and new RECTIFY output column formats:
  - Old format: 'position' column
  - New format: 'corrected_3prime' column

Usage:
    python generate_bedgraph_from_polished.py <input_dir> <output_dir>
    python generate_bedgraph_from_polished.py  # uses defaults

Author: Kevin R. Roy
"""

import pandas as pd
from pathlib import Path
from collections import defaultdict
import sys
import argparse

# Chromosome order for sorted output (yeast)
CHROM_ORDER = [
    'chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII',
    'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrmt'
]


def detect_position_column(filepath):
    """Detect which column contains the 3' end position."""
    header = pd.read_csv(filepath, sep='\t', nrows=0).columns.tolist()

    if 'corrected_3prime' in header:
        return 'corrected_3prime'
    elif 'position' in header:
        return 'position'
    else:
        raise ValueError(f"Cannot find position column in {filepath}. "
                         f"Expected 'corrected_3prime' or 'position'. "
                         f"Found columns: {header}")


def polished_to_bedgraph(input_tsv, output_prefix, normalize_rpm=True):
    """
    Convert polished 3' end TSV to strand-specific bedgraph files.

    Args:
        input_tsv: Path to polished 3' end TSV file
        output_prefix: Output path prefix (will add _plus.bedgraph, _minus.bedgraph)
        normalize_rpm: If True, output RPM-normalized values; else raw counts
    """
    print(f"Processing: {input_tsv.name}")

    # Detect column format
    pos_col = detect_position_column(input_tsv)
    print(f"  Position column: {pos_col}")

    # Read only needed columns for efficiency
    usecols = ['chrom', 'strand', pos_col]

    # Read TSV in chunks to handle large files
    chunk_size = 2_000_000
    plus_counts = defaultdict(lambda: defaultdict(int))
    minus_counts = defaultdict(lambda: defaultdict(int))
    total_reads = 0

    for chunk in pd.read_csv(input_tsv, sep='\t', usecols=usecols, chunksize=chunk_size):
        total_reads += len(chunk)

        # Vectorized groupby for speed
        plus_chunk = chunk[chunk['strand'] == '+']
        minus_chunk = chunk[chunk['strand'] == '-']

        # Count per chrom/position
        for (chrom, pos), count in plus_chunk.groupby(['chrom', pos_col]).size().items():
            plus_counts[chrom][pos] += count

        for (chrom, pos), count in minus_chunk.groupby(['chrom', pos_col]).size().items():
            minus_counts[chrom][pos] += count

        print(f"  Processed {total_reads:,} reads...", end='\r')

    print(f"\n  Total reads: {total_reads:,}")

    # Calculate RPM factor
    rpm_factor = 1e6 / total_reads if normalize_rpm and total_reads > 0 else 1.0

    # Write bedgraph files
    for strand_name, counts_dict in [('plus', plus_counts), ('minus', minus_counts)]:
        output_path = Path(f"{output_prefix}_{strand_name}.bedgraph")

        with open(output_path, 'w') as f:
            # Bedgraph header
            sample_name = input_tsv.stem.replace('_polished_3prime', '')
            f.write(f'track type=bedGraph name="{sample_name}_{strand_name}" '
                    f'description="Nanopore 3\' ends ({strand_name} strand, RPM)"\n')

            # Sort by chromosome order, then by position
            for chrom in CHROM_ORDER:
                if chrom not in counts_dict:
                    continue

                positions = sorted(counts_dict[chrom].keys())
                for pos in positions:
                    count = counts_dict[chrom][pos]
                    value = count * rpm_factor if normalize_rpm else count

                    # Bedgraph is 0-based, half-open
                    # Position is 1-based, so start = pos - 1, end = pos
                    start = int(pos) - 1
                    end = int(pos)

                    f.write(f"{chrom}\t{start}\t{end}\t{value:.4f}\n")

        print(f"  Wrote: {output_path.name}")

    return total_reads


def main():
    parser = argparse.ArgumentParser(
        description="Generate bedgraph files from RECTIFY polished 3' end TSV files."
    )
    parser.add_argument(
        'input_dir', nargs='?',
        default=None,
        help="Directory containing polished 3' end TSV files"
    )
    parser.add_argument(
        'output_dir', nargs='?',
        default=None,
        help="Output directory for bedgraph files"
    )
    parser.add_argument(
        '--pattern', default='*_polished_3prime.tsv',
        help="Glob pattern for input files (default: *_polished_3prime.tsv)"
    )
    args = parser.parse_args()

    # Set defaults if not provided
    if args.input_dir is None:
        input_dir = Path("/path/to/polished_3prime")  # set with --input-dir or configure here
    else:
        input_dir = Path(args.input_dir)

    if args.output_dir is None:
        output_dir = Path("/path/to/bedgraph/nanopore_polished")  # set with --output-dir or configure here
    else:
        output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all polished TSV files
    polished_files = sorted(input_dir.glob(args.pattern))

    if not polished_files:
        print(f"No files matching '{args.pattern}' found in {input_dir}")
        sys.exit(1)

    print(f"Found {len(polished_files)} polished 3' end files:")
    for f in polished_files:
        print(f"  - {f.name}")
    print()

    for tsv_file in polished_files:
        sample_name = tsv_file.stem.replace('_polished_3prime', '')
        output_prefix = output_dir / sample_name

        polished_to_bedgraph(tsv_file, output_prefix, normalize_rpm=True)
        print()

    print("Done! Bedgraph files saved to:")
    print(f"  {output_dir}")


if __name__ == "__main__":
    main()
