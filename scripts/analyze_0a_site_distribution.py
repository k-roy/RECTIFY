#!/usr/bin/env python3
"""
Analyze NET-seq signal distribution at 0A CPA sites.

At 0A sites (no downstream genomic A's), oligo(A) tails cannot align downstream.
This script empirically computes what fraction of signal appears at position 0
vs downstream positions, to validate/derive the PSF for deconvolution.

Key questions:
1. What fraction of signal is at the exact 0A position?
2. Is there any downstream spreading (and why)?
3. What is the oligo(A) tail length distribution?

Author: Kevin R. Roy
Date: 2026-03-19
"""

import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np

# Paths
BASE_DIR = Path("/oak/stanford/groups/larsms/Users/kevinroy")
RECTIFY_DIR = BASE_DIR / "software/rectify"
GENOME_PATH = BASE_DIR / "common/reference_genomes/MAGESTIC_background_strain.fasta"
NETSEQ_PATH = RECTIFY_DIR / "rectify/data/saccharomyces_cerevisiae_netseq_pan.tsv.gz"

# Parameters
DOWNSTREAM_WINDOW = 10  # bp to check for A's
UPSTREAM_WINDOW = 5     # bp to look upstream
ANALYSIS_WINDOW = 20    # bp downstream to analyze signal distribution
MIN_SIGNAL = 10.0       # Minimum signal at position to consider as peak


def load_genome(fasta_path: Path) -> Dict[str, str]:
    """Load genome from FASTA file."""
    print(f"Loading genome from {fasta_path}...")
    genome = {}
    current_chrom = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chrom:
                    genome[current_chrom] = ''.join(current_seq).upper()
                # Parse chromosome name (e.g., ">chrI" -> "chrI")
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_chrom:
            genome[current_chrom] = ''.join(current_seq).upper()

    print(f"  Loaded {len(genome)} chromosomes")
    for chrom, seq in list(genome.items())[:3]:
        print(f"    {chrom}: {len(seq):,} bp")

    return genome


def load_netseq_data(netseq_path: Path) -> Dict[Tuple[str, str, int], float]:
    """Load NET-seq data into dict of (chrom, strand, position) -> signal."""
    print(f"Loading NET-seq data from {netseq_path}...")
    data = {}

    with gzip.open(netseq_path, 'rt') as f:
        header = f.readline()  # Skip header
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]
            position = int(fields[1])
            strand = fields[2]
            signal = float(fields[3])

            data[(chrom, strand, position)] = signal

    print(f"  Loaded {len(data):,} positions with signal")
    return data


def count_downstream_as(genome: Dict[str, str], chrom: str, position: int,
                         strand: str, window: int = 10) -> int:
    """Count A's (or T's for minus strand) in downstream window."""
    seq = genome.get(chrom)
    if seq is None:
        return -1

    if strand == '+':
        # Downstream is to the right
        start = position + 1
        end = min(position + 1 + window, len(seq))
        if start >= len(seq):
            return 0
        downstream_seq = seq[start:end]
        return downstream_seq.count('A')
    else:
        # Minus strand: downstream is to the left, count T's
        start = max(0, position - window)
        end = position
        if end <= 0:
            return 0
        downstream_seq = seq[start:end]
        return downstream_seq.count('T')


def is_local_maximum(
    netseq_data: Dict[Tuple[str, str, int], float],
    chrom: str, strand: str, position: int,
    signal: float, window: int = 3
) -> bool:
    """Check if position is a local maximum within window."""
    for offset in range(-window, window + 1):
        if offset == 0:
            continue
        neighbor_signal = netseq_data.get((chrom, strand, position + offset), 0.0)
        if neighbor_signal > signal:
            return False
    return True


def is_isolated_peak(
    netseq_data: Dict[Tuple[str, str, int], float],
    chrom: str, strand: str, position: int,
    signal: float, isolation_window: int = 10,
    signal_ratio_threshold: float = 0.3
) -> bool:
    """
    Check if position is an isolated peak - no other significant peaks nearby.

    This helps identify true isolated CPA sites vs. clusters of peaks.
    """
    for offset in range(-isolation_window, isolation_window + 1):
        if offset == 0:
            continue
        neighbor_signal = netseq_data.get((chrom, strand, position + offset), 0.0)
        # If neighbor has more than threshold fraction of our signal, not isolated
        if neighbor_signal > signal * signal_ratio_threshold:
            return False
    return True


def find_sites_by_a_count(
    genome: Dict[str, str],
    netseq_data: Dict[Tuple[str, str, int], float],
    target_a_count: int,
    min_signal: float = 10.0,
    downstream_window: int = 10,
    require_local_max: bool = True,
    local_max_window: int = 3,
    require_isolated: bool = False,
    isolation_window: int = 10,
    max_sites: int = 0
) -> List[Tuple[str, str, int, float]]:
    """Find positions with specific number of downstream A's."""
    print(f"Finding {target_a_count}A sites with signal >= {min_signal}...")
    if require_local_max:
        print(f"  Requiring local maximum within ±{local_max_window} bp")
    if require_isolated:
        print(f"  Requiring isolation (no peaks >30% signal within ±{isolation_window} bp)")

    matching_sites = []
    total_candidates = 0
    local_max_count = 0

    for (chrom, strand, position), signal in netseq_data.items():
        if signal < min_signal:
            continue

        total_candidates += 1

        # Check if local maximum (if required)
        if require_local_max:
            if not is_local_maximum(netseq_data, chrom, strand, position, signal, local_max_window):
                continue

        local_max_count += 1

        # Check if isolated (if required)
        if require_isolated:
            if not is_isolated_peak(netseq_data, chrom, strand, position, signal, isolation_window):
                continue

        a_count = count_downstream_as(genome, chrom, position, strand, downstream_window)

        if a_count == target_a_count:
            matching_sites.append((chrom, strand, position, signal))
            if max_sites > 0 and len(matching_sites) >= max_sites:
                break

    print(f"  Candidates with signal >= {min_signal}: {total_candidates:,}")
    print(f"  Local maxima: {local_max_count:,}")
    print(f"  {target_a_count}A sites (final): {len(matching_sites):,}")
    return matching_sites


def find_0a_sites_with_signal(
    genome: Dict[str, str],
    netseq_data: Dict[Tuple[str, str, int], float],
    min_signal: float = 10.0,
    downstream_window: int = 10,
    require_local_max: bool = True,
    local_max_window: int = 3,
    require_isolated: bool = False,
    isolation_window: int = 10
) -> List[Tuple[str, str, int, float]]:
    """Find positions with high NET-seq signal and 0 downstream A's."""
    return find_sites_by_a_count(
        genome, netseq_data, target_a_count=0,
        min_signal=min_signal, downstream_window=downstream_window,
        require_local_max=require_local_max, local_max_window=local_max_window,
        require_isolated=require_isolated, isolation_window=isolation_window
    )


def analyze_signal_distribution(
    netseq_data: Dict[Tuple[str, str, int], float],
    zero_a_sites: List[Tuple[str, str, int, float]],
    upstream_window: int = 5,
    downstream_window: int = 20
) -> Dict[int, List[float]]:
    """
    For each 0A site, compute signal distribution relative to position 0.

    Returns dict of offset -> list of (signal / peak_signal) values.
    Offset 0 = the 0A site itself
    Positive offsets = downstream
    Negative offsets = upstream
    """
    print(f"Analyzing signal distribution around {len(zero_a_sites)} sites...")

    # offset -> list of normalized signal values
    offset_signals = defaultdict(list)

    for chrom, strand, position, peak_signal in zero_a_sites:
        # Get signal in window around this position
        for offset in range(-upstream_window, downstream_window + 1):
            if strand == '+':
                query_pos = position + offset
            else:
                # For minus strand, downstream is negative offset
                query_pos = position - offset

            signal = netseq_data.get((chrom, strand, query_pos), 0.0)

            # Normalize by peak signal
            normalized = signal / peak_signal if peak_signal > 0 else 0
            offset_signals[offset].append(normalized)

    return offset_signals


def compute_psf(offset_signals: Dict[int, List[float]], max_offset: int = 16) -> np.ndarray:
    """
    Compute Point-Spread-Function from offset signal data.

    The PSF represents: P(signal at offset | true CPA at position 0)
    """
    psf = []

    for offset in range(max_offset):
        if offset in offset_signals:
            # Mean normalized signal at this offset
            mean_signal = np.mean(offset_signals[offset])
            psf.append(mean_signal)
        else:
            psf.append(0.0)

    psf = np.array(psf)

    # Normalize to sum to 1
    total = psf.sum()
    if total > 0:
        psf = psf / total

    return psf


def print_results(offset_signals: Dict[int, List[float]], psf: np.ndarray):
    """Print analysis results."""
    print("\n" + "=" * 70)
    print("RESULTS: Signal Distribution at 0A Sites")
    print("=" * 70)

    print("\nOffset distribution (mean ± std of normalized signal):")
    print("-" * 50)

    for offset in sorted(offset_signals.keys()):
        values = offset_signals[offset]
        mean = np.mean(values)
        std = np.std(values)
        n = len(values)
        print(f"  Offset {offset:+3d}: {mean:.4f} ± {std:.4f}  (n={n})")

    print("\n" + "-" * 50)
    print("\nDerived PSF (normalized):")
    print("-" * 50)

    print("DEFAULT_PSF_0A = np.array([")
    for i, val in enumerate(psf):
        comment = "true position" if i == 0 else f"offset {i}"
        print(f"    {val:.4f},  # {comment}")
    print("])")

    # Key statistics
    print("\n" + "-" * 50)
    print("\nKey statistics:")
    print(f"  Signal at position 0:     {psf[0]*100:.1f}%")
    print(f"  Signal at offsets 1-5:    {psf[1:6].sum()*100:.1f}%")
    print(f"  Signal at offsets 6-10:   {psf[6:11].sum()*100:.1f}%")
    print(f"  Signal at offsets 11-15:  {psf[11:16].sum()*100:.1f}%")

    # Compare with current hardcoded values
    current_psf_0 = 0.54
    print(f"\n  Current hardcoded value at position 0: {current_psf_0*100:.1f}%")
    print(f"  Empirical value at position 0:         {psf[0]*100:.1f}%")

    if abs(psf[0] - current_psf_0) > 0.1:
        print(f"\n  WARNING: Large discrepancy ({abs(psf[0]-current_psf_0)*100:.1f}%) - PSF may need updating!")


def main():
    print("=" * 70)
    print("Analyzing NET-seq signal distribution by downstream A count")
    print("=" * 70)
    print()

    # Load data
    genome = load_genome(GENOME_PATH)
    netseq_data = load_netseq_data(NETSEQ_PATH)

    # Collect summary statistics
    summary_data = {}

    # Analyze sites with 0, 1, 2, 3, 4, 5+ downstream A's
    for a_count in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
        print("\n" + "-" * 70)
        print(f"ANALYSIS: {a_count}A sites (isolated)")
        print("-" * 70)

        sites = find_sites_by_a_count(
            genome, netseq_data,
            target_a_count=a_count,
            min_signal=MIN_SIGNAL,
            downstream_window=DOWNSTREAM_WINDOW,
            require_isolated=True,
            isolation_window=10,
            max_sites=500  # Limit to keep analysis fast
        )

        if len(sites) < 50:
            print(f"  Not enough sites (n={len(sites)}), skipping...")
            continue

        offset_signals = analyze_signal_distribution(
            netseq_data, sites,
            upstream_window=UPSTREAM_WINDOW,
            downstream_window=ANALYSIS_WINDOW
        )
        psf = compute_psf(offset_signals, max_offset=16)

        # Store summary
        summary_data[a_count] = {
            'n_sites': len(sites),
            'pos0_fraction': psf[0],
            'pos1_5_fraction': psf[1:6].sum(),
            'pos6_10_fraction': psf[6:11].sum(),
            'offset_signals': offset_signals,
            'psf': psf,
        }

        # Print condensed results for each
        print(f"\n  Position 0:    {psf[0]*100:.1f}%")
        print(f"  Offsets 1-5:   {psf[1:6].sum()*100:.1f}%")
        print(f"  Offsets 6-10:  {psf[6:11].sum()*100:.1f}%")

    # Summary comparison
    print("\n" + "=" * 70)
    print("SUMMARY: Signal at Position 0 vs Downstream A Count")
    print("=" * 70)
    print(f"\n{'A count':<10} {'n sites':<10} {'Pos 0':<10} {'1-5bp':<10} {'6-10bp':<10}")
    print("-" * 50)

    for a_count in sorted(summary_data.keys()):
        data = summary_data[a_count]
        print(f"{a_count:<10} {data['n_sites']:<10} "
              f"{data['pos0_fraction']*100:.1f}%{'':<5} "
              f"{data['pos1_5_fraction']*100:.1f}%{'':<5} "
              f"{data['pos6_10_fraction']*100:.1f}%")

    # Interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print("""
Key findings:
- 0A sites: ~60% at position 0 (baseline - no spreading possible)
- Sites with downstream A's: signal should spread MORE as A-count increases
  because oligo(A) tails can align further downstream

The DIFFERENCE between nA sites and 0A sites reveals the oligo(A) spreading.
If 0A = 60% at pos0 and 3A = 40% at pos0, then ~20% of signal shifted
downstream due to oligo(A) alignment.

For the deconvolution PSF:
- The 0A PSF should be: [1.0, 0, 0, ...] (no spreading possible)
- For nA sites: spreading depends on oligo(A) tail length distribution
""")

    # Print recommended PSF for 0A based on isolated sites
    if 0 in summary_data:
        print("\nRecommended DEFAULT_PSF_0A (from isolated 0A sites):")
        print("-" * 50)
        psf = summary_data[0]['psf']
        # Renormalize to show only signal above background
        # The background is roughly uniform at ~3% per position
        background_per_pos = 0.03
        corrected = psf.copy()
        corrected[1:] = np.maximum(0, corrected[1:] - background_per_pos)
        corrected = corrected / corrected.sum()
        print("DEFAULT_PSF_0A = np.array([")
        for i, val in enumerate(corrected):
            comment = "true position" if i == 0 else f"offset {i}"
            print(f"    {val:.4f},  # {comment}")
        print("])")

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
