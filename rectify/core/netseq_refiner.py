#!/usr/bin/env python3
"""
NET-seq refinement for RECTIFY.

This module implements NET-seq-based refinement of ambiguous CPA positions.
NET-seq captures nascent RNA 3' ends without poly(A) tails, providing ground
truth for CPA site localization.

Problem: A-tract ambiguity creates positional uncertainty (0-4 bp). NET-seq
data can resolve this by identifying the precise position of transcription
termination.

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import Dict, List, Optional, Tuple
from pathlib import Path
from collections import OrderedDict
import numpy as np

from ..config import (
    NETSEQ_PEAK_WINDOW,
    NETSEQ_PEAK_THRESHOLD,
    NETSEQ_SIGNAL_HIGH,
    NETSEQ_SIGNAL_MEDIUM,
)

# Try to import pyBigWig (optional dependency)
try:
    import pyBigWig
    PYBIGWIG_AVAILABLE = True
except ImportError:
    pyBigWig = None
    PYBIGWIG_AVAILABLE = False


class NetseqLoader:
    """
    Loader for NET-seq BigWig files with caching.

    This class handles loading and caching of NET-seq data from BigWig files.
    Multiple files can be loaded and their signals combined.
    """

    # Maximum number of cached signal arrays
    MAX_CACHE_SIZE = 10000

    def __init__(self, max_cache_size: int = None):
        """
        Initialize NET-seq loader with LRU cache.

        Args:
            max_cache_size: Maximum number of cached signal arrays (default: 10000)
        """
        self.bigwigs = {}  # {name: pyBigWig object}
        self._max_cache_size = max_cache_size or self.MAX_CACHE_SIZE
        # Use OrderedDict for LRU cache behavior
        self._cache = OrderedDict()  # {(chrom, start, end, strand): signal_array}

    def load_bigwig(self, filepath: str, name: str = None):
        """
        Load a BigWig file.

        Args:
            filepath: Path to BigWig file
            name: Optional name for this file (defaults to filename)
        """
        if not PYBIGWIG_AVAILABLE:
            raise ImportError("pyBigWig is required for NET-seq refinement. "
                            "Install with: pip install pyBigWig")

        if name is None:
            name = Path(filepath).stem

        bw = pyBigWig.open(str(filepath))
        self.bigwigs[name] = bw

    def load_directory(self, directory: str, pattern: str = "*.bw"):
        """
        Load all BigWig files from directory.

        Args:
            directory: Directory containing BigWig files
            pattern: Glob pattern for files (default: *.bw)
        """
        dir_path = Path(directory)
        for filepath in dir_path.glob(pattern):
            self.load_bigwig(str(filepath))

    def get_signal(
        self,
        chrom: str,
        start: int,
        end: int,
        strand: str = '+'
    ) -> np.ndarray:
        """
        Get NET-seq signal across genomic region.

        Args:
            chrom: Chromosome name
            start: Start position (0-based)
            end: End position (0-based exclusive)
            strand: Strand ('+' or '-')

        Returns:
            Array of signal values (one per base)
        """
        # Check cache (move to end for LRU behavior)
        cache_key = (chrom, start, end, strand)
        if cache_key in self._cache:
            # Move to end (most recently used)
            self._cache.move_to_end(cache_key)
            return self._cache[cache_key]

        # Combine signal from all loaded BigWigs
        length = end - start
        combined_signal = np.zeros(length)

        for name, bw in self.bigwigs.items():
            try:
                values = bw.values(chrom, start, end)
                # Replace None with 0
                values = [v if v is not None else 0.0 for v in values]
                combined_signal += np.array(values)
            except (KeyError, RuntimeError, ValueError):
                # Chromosome not found or invalid region in this BigWig, skip
                continue

        # Cache result with LRU eviction
        self._cache[cache_key] = combined_signal
        if len(self._cache) > self._max_cache_size:
            # Remove oldest (first) item
            self._cache.popitem(last=False)

        return combined_signal

    def clear_cache(self):
        """Clear signal cache."""
        self._cache.clear()

    def close(self):
        """Close all BigWig files."""
        for bw in self.bigwigs.values():
            bw.close()
        self.bigwigs.clear()
        self._cache.clear()


def find_peaks_in_window(
    signal: np.ndarray,
    positions: List[int],
    threshold: float = NETSEQ_PEAK_THRESHOLD
) -> List[Dict]:
    """
    Find peaks in NET-seq signal within window.

    Args:
        signal: Signal array (one value per base)
        positions: Genomic positions corresponding to signal array
        threshold: Minimum signal threshold (fraction of max)

    Returns:
        List of peak dicts with:
            - position: Genomic position of peak
            - signal: Signal value at peak
            - prominence: Peak prominence (height above surroundings)
    """
    if len(signal) == 0 or len(positions) == 0:
        return []

    if len(signal) != len(positions):
        raise ValueError("Signal and positions must have same length")

    # Find local maxima
    peaks = []
    max_signal = np.max(signal)

    if max_signal == 0:
        return peaks  # No signal

    # Threshold in absolute terms
    abs_threshold = max_signal * threshold

    for i in range(len(signal)):
        if signal[i] < abs_threshold:
            continue

        # Check if this is a local maximum
        is_peak = True

        # Check left neighbor
        if i > 0 and signal[i] <= signal[i-1]:
            is_peak = False

        # Check right neighbor
        if i < len(signal) - 1 and signal[i] <= signal[i+1]:
            is_peak = False

        if is_peak:
            # Calculate prominence (height above minimum of surroundings)
            left_min = np.min(signal[:i]) if i > 0 else 0
            right_min = np.min(signal[i+1:]) if i < len(signal) - 1 else 0
            surrounding_min = min(left_min, right_min)
            prominence = signal[i] - surrounding_min

            peaks.append({
                'position': positions[i],
                'signal': signal[i],
                'prominence': prominence,
            })

    return peaks


def select_best_peak(
    peaks: List[Dict],
    ambiguity_min: int,
    ambiguity_max: int
) -> Optional[Dict]:
    """
    Select best peak within ambiguity window.

    Prioritizes peaks by:
    1. Signal strength
    2. Prominence
    3. Leftmost position (for ties)

    Args:
        peaks: List of peak dicts
        ambiguity_min: Minimum position in ambiguity window
        ambiguity_max: Maximum position in ambiguity window

    Returns:
        Best peak dict, or None if no peaks in window
    """
    if not peaks:
        return None

    # Filter peaks within ambiguity window
    window_peaks = [
        p for p in peaks
        if ambiguity_min <= p['position'] <= ambiguity_max
    ]

    if not window_peaks:
        return None

    # Sort by signal (descending), then by prominence (descending)
    window_peaks.sort(
        key=lambda p: (p['signal'], p['prominence']),
        reverse=True
    )

    return window_peaks[0]


def assign_confidence(
    peak: Optional[Dict],
    all_peaks: List[Dict],
    ambiguity_range: int
) -> Tuple[str, str]:
    """
    Assign confidence level and method for refinement.

    Args:
        peak: Selected peak (or None)
        all_peaks: All peaks found in window
        ambiguity_range: Size of ambiguity window (bp)

    Returns:
        Tuple of (confidence, method)
        - confidence: 'high', 'medium', 'low'
        - method: 'netseq_peak', 'leftmost', 'modal', 'no_data'
    """
    if peak is None:
        return ('low', 'no_data')

    # High confidence: Single clear peak with strong signal
    if len(all_peaks) == 1 and peak['signal'] >= NETSEQ_SIGNAL_HIGH:
        return ('high', 'netseq_peak')

    # Medium confidence: Multiple peaks or moderate signal
    if len(all_peaks) > 1 or peak['signal'] >= NETSEQ_SIGNAL_MEDIUM:
        return ('medium', 'netseq_peak')

    # Low confidence: Weak signal
    return ('low', 'netseq_peak')


def refine_with_netseq(
    netseq_loader: NetseqLoader,
    chrom: str,
    ambiguity_min: int,
    ambiguity_max: int,
    strand: str,
    original_position: int
) -> Dict:
    """
    Refine CPA position using NET-seq data.

    This function:
    1. Extracts NET-seq signal across ambiguity window
    2. Finds peaks in the signal
    3. Selects the best peak
    4. Assigns confidence score

    Args:
        netseq_loader: NetseqLoader with BigWig files loaded
        chrom: Chromosome name
        ambiguity_min: Minimum position in ambiguity window
        ambiguity_max: Maximum position in ambiguity window
        strand: Gene strand ('+' or '-')
        original_position: Original CPA position before refinement

    Returns:
        Dict with:
            - refined_position: Best estimate of true CPA
            - confidence: 'high', 'medium', 'low'
            - method: Refinement method used
            - peak_signal: Signal at refined position (or 0)
            - n_peaks: Number of peaks found
    """
    # Expand window slightly for peak detection
    window_start = ambiguity_min - NETSEQ_PEAK_WINDOW
    window_end = ambiguity_max + NETSEQ_PEAK_WINDOW

    # Get NET-seq signal
    signal = netseq_loader.get_signal(chrom, window_start, window_end, strand)

    # Create position array
    positions = np.arange(window_start, window_end)

    # Find peaks
    all_peaks = find_peaks_in_window(signal, positions.tolist())

    # Select best peak within ambiguity window
    best_peak = select_best_peak(all_peaks, ambiguity_min, ambiguity_max)

    # Assign confidence
    confidence, method = assign_confidence(
        best_peak,
        all_peaks,
        ambiguity_max - ambiguity_min
    )

    # Determine refined position
    if best_peak:
        refined_position = best_peak['position']
        peak_signal = best_peak['signal']
    else:
        # No peak found - use leftmost position (most conservative)
        refined_position = ambiguity_min
        peak_signal = 0.0
        method = 'leftmost'

    return {
        'refined_position': refined_position,
        'confidence': confidence,
        'method': method,
        'peak_signal': peak_signal,
        'n_peaks': len(all_peaks),
        'shift_from_original': refined_position - original_position,
    }


def refine_batch(
    netseq_loader: NetseqLoader,
    positions: List[Dict]
) -> List[Dict]:
    """
    Refine batch of CPA positions using NET-seq.

    Args:
        netseq_loader: NetseqLoader with BigWig files loaded
        positions: List of position dicts with keys:
            - chrom: Chromosome
            - ambiguity_min: Min position
            - ambiguity_max: Max position
            - strand: Strand
            - original_position: Original CPA position

    Returns:
        List of refinement result dicts (same order as input)
    """
    results = []

    for pos_data in positions:
        result = refine_with_netseq(
            netseq_loader,
            pos_data['chrom'],
            pos_data['ambiguity_min'],
            pos_data['ambiguity_max'],
            pos_data['strand'],
            pos_data['original_position']
        )
        results.append(result)

    return results


def calculate_netseq_statistics(refinement_results: List[Dict]) -> Dict:
    """
    Calculate summary statistics for NET-seq refinement.

    Args:
        refinement_results: List of refinement result dicts

    Returns:
        Dict with summary statistics
    """
    if not refinement_results:
        return {
            'total': 0,
            'refined': 0,
            'by_confidence': {},
            'by_method': {},
        }

    # Count by confidence
    by_confidence = {}
    for result in refinement_results:
        conf = result['confidence']
        by_confidence[conf] = by_confidence.get(conf, 0) + 1

    # Count by method
    by_method = {}
    for result in refinement_results:
        method = result['method']
        by_method[method] = by_method.get(method, 0) + 1

    # Count positions that were refined (changed)
    refined = sum(1 for r in refinement_results if r['shift_from_original'] != 0)

    # Calculate mean shift
    shifts = [abs(r['shift_from_original']) for r in refinement_results]
    mean_shift = np.mean(shifts) if shifts else 0.0

    return {
        'total': len(refinement_results),
        'refined': refined,
        'refinement_rate': refined / len(refinement_results),
        'mean_shift': mean_shift,
        'by_confidence': by_confidence,
        'by_method': by_method,
    }


def format_netseq_report(stats: Dict) -> str:
    """
    Format NET-seq refinement statistics as human-readable report.

    Args:
        stats: Statistics dict from calculate_netseq_statistics()

    Returns:
        Formatted report string
    """
    report = []
    report.append("=" * 60)
    report.append("NET-seq Refinement Summary")
    report.append("=" * 60)
    report.append("")

    report.append("Overall:")
    report.append(f"  Total positions:        {stats['total']:,}")
    report.append(f"  Refined:                {stats['refined']:,} ({stats['refinement_rate']:.1%})")
    report.append(f"  Mean shift:             {stats['mean_shift']:.1f} bp")
    report.append("")

    report.append("Confidence Levels:")
    for conf in ['high', 'medium', 'low']:
        count = stats['by_confidence'].get(conf, 0)
        pct = 100.0 * count / stats['total'] if stats['total'] > 0 else 0.0
        report.append(f"  {conf:10s}          {count:7,} ({pct:5.1f}%)")
    report.append("")

    report.append("Refinement Methods:")
    for method, count in stats['by_method'].items():
        pct = 100.0 * count / stats['total'] if stats['total'] > 0 else 0.0
        report.append(f"  {method:20s} {count:7,} ({pct:5.1f}%)")

    report.append("")
    report.append("=" * 60)

    return "\n".join(report)
