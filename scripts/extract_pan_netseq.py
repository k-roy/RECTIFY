#!/usr/bin/env python3
"""
Extract pan-mutant NET-seq consensus signal from all datasets.

Combines signal from ALL NET-seq replicates (WT + mutants) to get maximum
coverage of CPA sites. WT signal is tracked separately to enable preferential
weighting in apportionment.

Strategy:
- Include ALL NET-seq datasets to maximize site coverage
- Store WT signal separately for preferential weighting
- For sites with both WT and mutant signal, use WT proportions
- For sites only in mutants, include them but flag as mutant-only

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


def process_bigwig_pair(plus_path, minus_path, signal_dict, is_wt=False):
    """
    Process a BigWig plus/minus pair and accumulate signal.

    Args:
        plus_path: Path to plus strand BigWig
        minus_path: Path to minus strand BigWig
        signal_dict: Dict to accumulate signal into
        is_wt: If True, also track WT-specific signal

    Returns:
        wt_signal dict if is_wt=True, else None
    """
    print(f"  Processing {Path(plus_path).stem}... (WT={is_wt})")

    wt_signal = {} if is_wt else None

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
                    signal_dict[key] = signal_dict.get(key, 0) + val
                    if is_wt and wt_signal is not None:
                        wt_signal[key] = wt_signal.get(key, 0) + val
        except Exception:
            pass

        # Minus strand
        try:
            values = minus_bw.values(chrom, 0, length)
            for i, val in enumerate(values):
                if val is not None and val > 0:
                    key = (std_chrom, '-', i)
                    signal_dict[key] = signal_dict.get(key, 0) + val
                    if is_wt and wt_signal is not None:
                        wt_signal[key] = wt_signal.get(key, 0) + val
        except Exception:
            pass

    plus_bw.close()
    minus_bw.close()

    return wt_signal


def combine_pan_netseq(wt_pairs, mutant_pairs, output_path, min_signal=5.0):
    """
    Combine all NET-seq files with WT preference.

    For each position, we track:
    - total_signal: sum across all datasets
    - wt_signal: sum across WT datasets only

    Apportionment strategy:
    - If WT signal exists, use WT-derived ratios
    - If only mutant signal exists, use mutant ratios (but flag)
    """
    print(f"Combining {len(wt_pairs)} WT + {len(mutant_pairs)} mutant BigWig pairs...")

    # Collect all signal
    total_signal = {}  # (chrom, strand, pos) -> total_signal
    wt_signal = {}     # (chrom, strand, pos) -> wt_signal

    # Process WT datasets first
    print("\nProcessing WT datasets:")
    for plus_path, minus_path in wt_pairs:
        wt_contrib = process_bigwig_pair(plus_path, minus_path, total_signal, is_wt=True)
        if wt_contrib:
            for key, val in wt_contrib.items():
                wt_signal[key] = wt_signal.get(key, 0) + val

    # Process mutant datasets
    print("\nProcessing mutant datasets:")
    for plus_path, minus_path in mutant_pairs:
        process_bigwig_pair(plus_path, minus_path, total_signal, is_wt=False)

    # Filter by threshold
    print(f"\nTotal positions with any signal: {len(total_signal):,}")
    filtered = {k: v for k, v in total_signal.items() if v >= min_signal}
    print(f"Positions with combined signal >= {min_signal}: {len(filtered):,}")

    # Count WT vs mutant-only
    wt_positions = sum(1 for k in filtered if k in wt_signal and wt_signal[k] >= 1.0)
    mutant_only = len(filtered) - wt_positions
    print(f"  - Positions with WT signal: {wt_positions:,}")
    print(f"  - Positions with mutant-only signal: {mutant_only:,}")

    # Write output with WT signal column for preferential weighting
    output_path = Path(output_path)
    open_func = gzip.open if str(output_path).endswith('.gz') else open

    print(f"\nWriting to {output_path}...")
    with open_func(output_path, 'wt') as f:
        # Header includes wt_signal for preferential weighting
        f.write("chrom\tposition\tstrand\tsignal\twt_signal\n")
        for (chrom, strand, pos), signal in sorted(filtered.items()):
            wt_sig = wt_signal.get((chrom, strand, pos), 0.0)
            f.write(f"{chrom}\t{pos}\t{strand}\t{signal:.2f}\t{wt_sig:.2f}\n")

    size_bytes = output_path.stat().st_size
    size_str = f"{size_bytes / (1024*1024):.1f} MB" if size_bytes > 1024*1024 else f"{size_bytes / 1024:.1f} KB"
    print(f"Output size: {size_str}")
    print("Done!")


def find_bigwig_pairs(bigwig_dir):
    """
    Find all BigWig plus/minus pairs in a directory and classify as WT or mutant.

    Returns:
        wt_pairs: List of (plus_path, minus_path) for WT samples
        mutant_pairs: List of (plus_path, minus_path) for mutant samples
    """
    bigwig_dir = Path(bigwig_dir)

    wt_pairs = []
    mutant_pairs = []

    # Find all plus strand BigWig files
    plus_files = sorted(bigwig_dir.glob("*.plus.bw"))

    for plus_path in plus_files:
        # Get corresponding minus strand file
        minus_path = plus_path.parent / plus_path.name.replace(".plus.bw", ".minus.bw")
        if not minus_path.exists():
            print(f"  Warning: Missing minus file for {plus_path.name}")
            continue

        # Parse sample name: netseq_2022_<strain>_<replicate>.plus.bw
        # or: netseq_2011_<strain>.plus.bw
        stem = plus_path.stem.replace(".plus", "")
        parts = stem.split("_")

        # Handle different naming conventions
        if len(parts) >= 3:
            strain = parts[2].lower()  # e.g., "wt", "dst1d", "cbc1"
        else:
            strain = parts[-1].lower()

        # Classify as WT or mutant
        if strain == "wt":
            wt_pairs.append((str(plus_path), str(minus_path)))
        else:
            mutant_pairs.append((str(plus_path), str(minus_path)))

    return wt_pairs, mutant_pairs


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extract pan-mutant NET-seq consensus")
    parser.add_argument("--bigwig-dir", type=str, required=True,
                        help="Directory containing BigWig files")
    parser.add_argument("--output", "-o", type=str, required=True,
                        help="Output TSV path (can be .gz)")
    parser.add_argument("--min-signal", type=float, default=5.0,
                        help="Minimum combined signal threshold (default: 5.0)")
    parser.add_argument("--additional-wt-dir", type=str,
                        help="Additional directory with WT-only BigWig files")

    args = parser.parse_args()

    # Find all BigWig pairs
    wt_pairs, mutant_pairs = find_bigwig_pairs(args.bigwig_dir)

    # Optionally add additional WT directory (e.g., 2011 data)
    if args.additional_wt_dir:
        extra_wt, _ = find_bigwig_pairs(args.additional_wt_dir)
        wt_pairs.extend(extra_wt)

    print(f"Found {len(wt_pairs)} WT BigWig pairs")
    print(f"Found {len(mutant_pairs)} mutant BigWig pairs")

    if len(wt_pairs) == 0 and len(mutant_pairs) == 0:
        print("Error: No BigWig pairs found")
        exit(1)

    combine_pan_netseq(wt_pairs, mutant_pairs, args.output, min_signal=args.min_signal)
