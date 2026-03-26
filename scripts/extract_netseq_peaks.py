#!/usr/bin/env python3
"""
Extract NET-seq signal data for bundling with RECTIFY.

This script processes BigWig files and extracts sparse signal data
that can be bundled with the package for A-tract ambiguity resolution.

Author: Kevin R. Roy
Date: 2026-03-18
"""

import argparse
import gzip
from pathlib import Path
import numpy as np

try:
    import pyBigWig
except ImportError:
    print("Error: pyBigWig required. Install with: pip install pyBigWig")
    exit(1)


# Chromosome name mappings for S. cerevisiae
# RefSeq accessions -> standard names
YEAST_CHROM_MAP = {
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
    'ref|NC_001224|': 'chrM',
    # Also support SGD format
    'NC_001133': 'chrI',
    'NC_001134': 'chrII',
    'NC_001135': 'chrIII',
    'NC_001136': 'chrIV',
    'NC_001137': 'chrV',
    'NC_001138': 'chrVI',
    'NC_001139': 'chrVII',
    'NC_001140': 'chrVIII',
    'NC_001141': 'chrIX',
    'NC_001142': 'chrX',
    'NC_001143': 'chrXI',
    'NC_001144': 'chrXII',
    'NC_001145': 'chrXIII',
    'NC_001146': 'chrXIV',
    'NC_001147': 'chrXV',
    'NC_001148': 'chrXVI',
    'NC_001224': 'chrM',
}


def extract_sparse_signal(
    plus_bw_path: str,
    minus_bw_path: str,
    output_path: str,
    min_signal: float = 1.0,
):
    """
    Extract sparse NET-seq signal from BigWig files.

    Args:
        plus_bw_path: Path to plus strand BigWig
        minus_bw_path: Path to minus strand BigWig
        output_path: Output TSV path (will be gzipped if ends in .gz)
        min_signal: Minimum signal threshold to include
    """
    print(f"Loading BigWig files...")
    plus_bw = pyBigWig.open(plus_bw_path)
    minus_bw = pyBigWig.open(minus_bw_path)

    # Get chromosome info
    chroms = plus_bw.chroms()

    print(f"Found {len(chroms)} chromosomes")
    print(f"Extracting positions with signal >= {min_signal}...")

    positions = []

    for chrom, length in chroms.items():
        # Map chromosome name to standard format
        std_chrom = YEAST_CHROM_MAP.get(chrom, chrom)

        # Plus strand
        try:
            values = plus_bw.values(chrom, 0, length)
            for i, val in enumerate(values):
                if val is not None and val >= min_signal:
                    positions.append((std_chrom, i, '+', val))
        except Exception as e:
            print(f"  Warning: {chrom} plus strand: {e}")

        # Minus strand
        try:
            values = minus_bw.values(chrom, 0, length)
            for i, val in enumerate(values):
                if val is not None and val >= min_signal:
                    positions.append((std_chrom, i, '-', val))
        except Exception as e:
            print(f"  Warning: {chrom} minus strand: {e}")

    plus_bw.close()
    minus_bw.close()

    print(f"Found {len(positions):,} positions with signal >= {min_signal}")

    # Write output
    output_path = Path(output_path)
    open_func = gzip.open if str(output_path).endswith('.gz') else open

    print(f"Writing to {output_path}...")
    with open_func(output_path, 'wt') as f:
        f.write("chrom\tposition\tstrand\tsignal\n")
        for chrom, pos, strand, signal in positions:
            f.write(f"{chrom}\t{pos}\t{strand}\t{signal:.2f}\n")

    # Report size
    size_bytes = output_path.stat().st_size
    if size_bytes > 1024 * 1024:
        size_str = f"{size_bytes / (1024*1024):.1f} MB"
    elif size_bytes > 1024:
        size_str = f"{size_bytes / 1024:.1f} KB"
    else:
        size_str = f"{size_bytes} bytes"

    print(f"Output size: {size_str}")
    print("Done!")


def main():
    parser = argparse.ArgumentParser(
        description="Extract NET-seq signal for RECTIFY bundling"
    )
    parser.add_argument(
        "--plus", required=True,
        help="Plus strand BigWig file"
    )
    parser.add_argument(
        "--minus", required=True,
        help="Minus strand BigWig file"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV file (use .gz extension for compression)"
    )
    parser.add_argument(
        "--min-signal", type=float, default=1.0,
        help="Minimum signal threshold (default: 1.0)"
    )

    args = parser.parse_args()

    extract_sparse_signal(
        args.plus,
        args.minus,
        args.output,
        args.min_signal,
    )


if __name__ == "__main__":
    main()
