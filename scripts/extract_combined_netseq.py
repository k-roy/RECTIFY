#!/usr/bin/env python3
"""
Extract combined NET-seq signal from multiple WT replicates.

Combines signal from all WT replicates (2011 + 2022) to get consensus CPA usage.

Author: Kevin R. Roy
Date: 2026-03-18
"""

import gzip
from pathlib import Path
import numpy as np

try:
    import pyBigWig
except ImportError:
    print("Error: pyBigWig required. Install with: pip install pyBigWig")
    exit(1)


# Chromosome name mappings for S. cerevisiae
YEAST_CHROM_MAP = {
    'ref|NC_001133|': 'chrI', 'ref|NC_001134|': 'chrII', 'ref|NC_001135|': 'chrIII',
    'ref|NC_001136|': 'chrIV', 'ref|NC_001137|': 'chrV', 'ref|NC_001138|': 'chrVI',
    'ref|NC_001139|': 'chrVII', 'ref|NC_001140|': 'chrVIII', 'ref|NC_001141|': 'chrIX',
    'ref|NC_001142|': 'chrX', 'ref|NC_001143|': 'chrXI', 'ref|NC_001144|': 'chrXII',
    'ref|NC_001145|': 'chrXIII', 'ref|NC_001146|': 'chrXIV', 'ref|NC_001147|': 'chrXV',
    'ref|NC_001148|': 'chrXVI', 'ref|NC_001224|': 'chrM',
    'NC_001133': 'chrI', 'NC_001134': 'chrII', 'NC_001135': 'chrIII',
    'NC_001136': 'chrIV', 'NC_001137': 'chrV', 'NC_001138': 'chrVI',
    'NC_001139': 'chrVII', 'NC_001140': 'chrVIII', 'NC_001141': 'chrIX',
    'NC_001142': 'chrX', 'NC_001143': 'chrXI', 'NC_001144': 'chrXII',
    'NC_001145': 'chrXIII', 'NC_001146': 'chrXIV', 'NC_001147': 'chrXV',
    'NC_001148': 'chrXVI', 'NC_001224': 'chrM',
}


def combine_bigwigs(bigwig_files, output_path, min_signal=5.0):
    """
    Combine multiple BigWig files (plus/minus pairs) into sparse signal TSV.

    Signal is summed across all replicates.
    """
    print(f"Combining {len(bigwig_files)} BigWig file pairs...")

    # Collect signal from all files
    combined_signal = {}  # (chrom, strand, pos) -> total_signal

    for plus_path, minus_path in bigwig_files:
        print(f"  Processing {Path(plus_path).stem}...")

        plus_bw = pyBigWig.open(plus_path)
        minus_bw = pyBigWig.open(minus_path)

        chroms = plus_bw.chroms()

        for chrom, length in chroms.items():
            std_chrom = YEAST_CHROM_MAP.get(chrom, chrom)

            # Plus strand
            try:
                values = plus_bw.values(chrom, 0, length)
                for i, val in enumerate(values):
                    if val is not None and val > 0:
                        key = (std_chrom, '+', i)
                        combined_signal[key] = combined_signal.get(key, 0) + val
            except Exception:
                pass

            # Minus strand
            try:
                values = minus_bw.values(chrom, 0, length)
                for i, val in enumerate(values):
                    if val is not None and val > 0:
                        key = (std_chrom, '-', i)
                        combined_signal[key] = combined_signal.get(key, 0) + val
            except Exception:
                pass

        plus_bw.close()
        minus_bw.close()

    # Filter by threshold
    print(f"\nTotal positions with any signal: {len(combined_signal):,}")
    filtered = {k: v for k, v in combined_signal.items() if v >= min_signal}
    print(f"Positions with combined signal >= {min_signal}: {len(filtered):,}")

    # Write output
    output_path = Path(output_path)
    open_func = gzip.open if str(output_path).endswith('.gz') else open

    print(f"\nWriting to {output_path}...")
    with open_func(output_path, 'wt') as f:
        f.write("chrom\tposition\tstrand\tsignal\n")
        for (chrom, strand, pos), signal in sorted(filtered.items()):
            f.write(f"{chrom}\t{pos}\t{strand}\t{signal:.2f}\n")

    size_bytes = output_path.stat().st_size
    size_str = f"{size_bytes / (1024*1024):.1f} MB" if size_bytes > 1024*1024 else f"{size_bytes / 1024:.1f} KB"
    print(f"Output size: {size_str}")
    print("Done!")


if __name__ == "__main__":
    # All WT NET-seq datasets
    base_dir = "/oak/stanford/groups/larsms/Users/kevinroy/projects/roadblocks/processed_data/bigwig/reprocessed"

    bigwig_pairs = [
        (f"{base_dir}/netseq_2011_WT.plus.bw", f"{base_dir}/netseq_2011_WT.minus.bw"),
        (f"{base_dir}/netseq_2022_wt1.plus.bw", f"{base_dir}/netseq_2022_wt1.minus.bw"),
        (f"{base_dir}/netseq_2022_wt2.plus.bw", f"{base_dir}/netseq_2022_wt2.minus.bw"),
    ]

    output = "/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/saccharomyces_cerevisiae_netseq_wt.tsv.gz"

    combine_bigwigs(bigwig_pairs, output, min_signal=5.0)
