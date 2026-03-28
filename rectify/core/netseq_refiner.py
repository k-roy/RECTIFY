#!/usr/bin/env python3
"""
NET-seq refinement for RECTIFY.

This module implements NET-seq-based refinement of ambiguous CPA positions.
NET-seq captures nascent RNA 3' ends without poly(A) tails, providing ground
truth for CPA site localization.

Problem: A-tract ambiguity creates positional uncertainty (0-4 bp). NET-seq
data can resolve this by identifying the precise position of transcription
termination.

Key insight: NET-seq captures oligo-adenylated intermediates (~6.6 bp mean tail).
When genomic A's exist downstream of the CPA site, these short tails align to
them, spreading the signal DOWNSTREAM. NNLS deconvolution recovers the true
peak positions by removing this spreading artifact.

Algorithm:
1. Load 0A PSF (Point-Spread-Function) from training data
2. Build convolution matrix for A-tract region
3. Deconvolve observed NET-seq signal using regularized NNLS
4. Assign reads proportionally to deconvolved peak intensities

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
from collections import OrderedDict
import threading
import numpy as np

from ..config import (
    NETSEQ_PEAK_WINDOW,
    NETSEQ_PEAK_THRESHOLD,
    NETSEQ_SIGNAL_HIGH,
    NETSEQ_SIGNAL_MEDIUM,
    NETSEQ_PSF_POSITION_0_FRACTION,
    NETSEQ_PSF_DOWNSTREAM_FRACTION,
    NETSEQ_DECONV_REGULARIZATION,
)

# Try to import pyBigWig (optional dependency)
try:
    import pyBigWig
    PYBIGWIG_AVAILABLE = True
except ImportError:
    pyBigWig = None
    PYBIGWIG_AVAILABLE = False

# Try to import scipy for NNLS (optional but recommended)
try:
    from scipy.optimize import nnls
    SCIPY_AVAILABLE = True
except ImportError:
    nnls = None
    SCIPY_AVAILABLE = False


# Default 0A PSF (Point-Spread-Function)
# Empirically derived from GSE25107 (2011) and GSE159603 (2022) NET-seq data:
#   - 54% at true position (NETSEQ_PSF_POSITION_0_FRACTION)
#   - 46% spreading downstream (NETSEQ_PSF_DOWNSTREAM_FRACTION)
#   - Downstream decay follows oligo-A tail distribution (~Poisson, mean 5.5 bp)
#
# See config.py for empirical constants and their derivation
DEFAULT_PSF_0A = np.array([
    NETSEQ_PSF_POSITION_0_FRACTION,  # offset 0: true position (~54%)
    0.005,  # offset 1
    0.015,  # offset 2
    0.032,  # offset 3
    0.052,  # offset 4
    0.067,  # offset 5
    0.073,  # offset 6 (mode of downstream distribution)
    0.067,  # offset 7
    0.052,  # offset 8
    0.032,  # offset 9
    0.020,  # offset 10
    0.015,  # offset 11+
    0.010,
    0.008,
    0.005,
    0.003,
])
# Normalize to sum to 1 (preserving relative downstream weights)
DEFAULT_PSF_0A = DEFAULT_PSF_0A / DEFAULT_PSF_0A.sum()


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
        Initialize NET-seq loader with thread-safe LRU cache.

        Args:
            max_cache_size: Maximum number of cached signal arrays (default: 10000)
        """
        self.bigwigs = {}  # {name: pyBigWig object}
        self._max_cache_size = max_cache_size or self.MAX_CACHE_SIZE
        # Use OrderedDict for LRU cache behavior, with lock for thread safety
        self._cache = OrderedDict()  # {(chrom, start, end, strand): signal_array}
        self._cache_lock = threading.Lock()

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
        with self._cache_lock:
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
            except (KeyError, RuntimeError, ValueError) as e:
                # Chromosome not found or invalid region in this BigWig
                import logging
                logging.getLogger(__name__).debug(
                    f"No NET-seq signal for {chrom}:{start}-{end} in {name}: {e}"
                )
                continue

        # Cache result with LRU eviction (thread-safe)
        with self._cache_lock:
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


def load_psf_model(filepath: str) -> np.ndarray:
    """
    Load PSF (Point-Spread-Function) model from file.

    The PSF describes how signal from a true CPA position spreads
    downstream due to oligo-adenylation alignment to genomic A's.

    Args:
        filepath: Path to TSV file with columns 'acount', 'offset', 'probability'

    Returns:
        1D array indexed by offset (0 = true position)
    """
    import pandas as pd
    psf_df = pd.read_csv(filepath, sep='\t')
    psf_0a = psf_df[psf_df['acount'] == 0].copy()

    max_offset = int(psf_0a['offset'].max())
    psf = np.zeros(max_offset + 1)
    for _, row in psf_0a.iterrows():
        psf[int(row['offset'])] = row['probability']

    return psf


def build_convolution_matrix(
    tract_length: int,
    psf: np.ndarray,
    include_first_non_a: bool = True
) -> np.ndarray:
    """
    Build convolution matrix for NNLS deconvolution.

    A[i, j] = P(observe signal at position j | true peak at position i)

    CRITICAL: Signal spreading beyond tract boundaries is LOST, not
    accumulated. Row sums < 1.0 at boundaries correctly represent
    the fact that only ~54% of boundary peak signal is observable.

    Args:
        tract_length: Length of the A-tract
        psf: 0A point-spread-function (probability by offset)
        include_first_non_a: Include first non-A position in matrix

    Returns:
        Convolution matrix of shape (n_positions, n_positions)
        Note: Row sums may be < 1.0 for positions near boundaries
    """
    n = tract_length + (1 if include_first_non_a else 0)
    A = np.zeros((n, n))

    for i in range(n):  # True peak position
        for offset in range(len(psf)):
            j = i + offset  # Observed position
            if j < n:
                A[i, j] = psf[offset]
            # Out of bounds: Signal is LOST (correct physical model)

    return A


def deconvolve_signal(
    observed: np.ndarray,
    A: np.ndarray,
    regularization: float = None
) -> np.ndarray:
    """
    Deconvolve observed signal to recover true peak positions.

    Solves: observed = A @ true_peaks
    Using non-negative least squares with L2 regularization.

    Args:
        observed: Observed signal vector (n_positions,)
        A: Convolution matrix (n_positions, n_positions)
        regularization: L2 regularization strength

    Returns:
        Deconvolved signal (true peak intensities)

    Raises:
        ImportError: If scipy is not available
    """
    if not SCIPY_AVAILABLE:
        raise ImportError(
            "scipy is required for NNLS deconvolution. "
            "Install with: pip install scipy"
        )

    # Use default regularization from config if not specified
    if regularization is None:
        regularization = NETSEQ_DECONV_REGULARIZATION

    n = len(observed)

    # Add regularization: [A; sqrt(reg)*I] @ x = [observed; 0]
    A_reg = np.vstack([A, np.sqrt(regularization) * np.eye(n)])
    b_reg = np.concatenate([observed, np.zeros(n)])

    # Solve with non-negative least squares
    deconvolved, _ = nnls(A_reg, b_reg)

    return deconvolved


def identify_deconvolved_peaks(
    deconvolved: np.ndarray,
    positions: List[int],
    threshold_frac: float = 0.1
) -> List[Dict]:
    """
    Identify true peak positions from deconvolved signal.

    Args:
        deconvolved: Deconvolved signal vector
        positions: Genomic positions corresponding to each index
        threshold_frac: Minimum fraction of max signal to call a peak

    Returns:
        List of peak dicts with 'position' and 'signal'
    """
    if len(deconvolved) == 0 or deconvolved.max() == 0:
        return []

    threshold = threshold_frac * deconvolved.max()
    peaks = []

    for i, val in enumerate(deconvolved):
        if val >= threshold and i < len(positions):
            peaks.append({
                'position': positions[i],
                'signal': float(val),
            })

    return peaks


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
    original_position: int,
    use_deconvolution: bool = True,
    psf: Optional[np.ndarray] = None,
    proportional_split: bool = True,
) -> Union[Dict, List[Dict]]:
    """
    Refine CPA position using NET-seq data.

    This function:
    1. Extracts NET-seq signal across ambiguity window
    2. Optionally deconvolves signal using NNLS to remove spreading artifact
    3. Finds peaks in the (deconvolved) signal
    4. Returns proportional assignments OR winner-take-all

    Args:
        netseq_loader: NetseqLoader with BigWig files loaded
        chrom: Chromosome name
        ambiguity_min: Minimum position in ambiguity window
        ambiguity_max: Maximum position in ambiguity window
        strand: Gene strand ('+' or '-')
        original_position: Original CPA position before refinement
        use_deconvolution: Whether to use NNLS deconvolution (default: True)
        psf: Point-spread-function for deconvolution (default: built-in)
        proportional_split: Return proportional assignments (default: True)

    Returns:
        If proportional_split=True:
            List of dicts, each with:
                - assigned_position: Position for this fraction
                - fraction: Fraction of read assigned here (0-1)
                - confidence: 'high', 'medium', 'low', or 'split'
                - method: Refinement method used
                - peak_signal: NET-seq signal at position

        If proportional_split=False:
            Single dict with:
                - refined_position: Best estimate of true CPA
                - confidence: 'high', 'medium', 'low'
                - method: Refinement method used
                - peak_signal: Signal at refined position
                - n_peaks: Number of peaks found
                - shift_from_original: Position change from input
    """
    # Expand window slightly for peak detection
    window_start = ambiguity_min - NETSEQ_PEAK_WINDOW
    window_end = ambiguity_max + NETSEQ_PEAK_WINDOW

    # Get NET-seq signal
    signal = netseq_loader.get_signal(chrom, window_start, window_end, strand)

    # Create position array
    positions = np.arange(window_start, window_end)

    # Choose deconvolution or simple peak finding
    if use_deconvolution and SCIPY_AVAILABLE:
        all_peaks = _refine_with_deconvolution(
            signal, positions, ambiguity_min, ambiguity_max, strand, psf
        )
    else:
        # Fall back to simple peak finding
        all_peaks = find_peaks_in_window(signal, positions.tolist())

    # Filter peaks within ambiguity window
    window_peaks = [
        p for p in all_peaks
        if ambiguity_min <= p['position'] <= ambiguity_max
    ]

    # Handle no peaks case
    if not window_peaks:
        if proportional_split:
            return [{
                'assigned_position': ambiguity_min,
                'fraction': 1.0,
                'confidence': 'low',
                'method': 'leftmost',
                'peak_signal': 0.0,
                'n_peaks': 0,
            }]
        else:
            return {
                'refined_position': ambiguity_min,
                'confidence': 'low',
                'method': 'leftmost',
                'peak_signal': 0.0,
                'n_peaks': 0,
                'shift_from_original': ambiguity_min - original_position,
            }

    # Proportional splitting
    if proportional_split:
        return _proportional_assignment(
            window_peaks, original_position, use_deconvolution
        )
    else:
        # Winner-take-all (legacy behavior)
        best_peak = select_best_peak(all_peaks, ambiguity_min, ambiguity_max)
        confidence, method = assign_confidence(
            best_peak, all_peaks, ambiguity_max - ambiguity_min
        )

        if best_peak:
            refined_position = best_peak['position']
            peak_signal = best_peak['signal']
        else:
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


def _refine_with_deconvolution(
    signal: np.ndarray,
    positions: np.ndarray,
    ambiguity_min: int,
    ambiguity_max: int,
    strand: str,
    psf: Optional[np.ndarray] = None,
) -> List[Dict]:
    """
    Deconvolve NET-seq signal and identify true peaks.

    Args:
        signal: Raw NET-seq signal array
        positions: Genomic positions for each signal index
        ambiguity_min: Start of ambiguity window
        ambiguity_max: End of ambiguity window
        strand: Strand ('+' or '-')
        psf: Point-spread-function (default: built-in)

    Returns:
        List of peak dicts with 'position', 'signal', 'prominence'
    """
    if psf is None:
        psf = DEFAULT_PSF_0A

    # Build convolution matrix for the tract
    tract_length = ambiguity_max - ambiguity_min + 1
    if tract_length < 2:
        # For single-base or empty tracts, deconvolution is underdetermined;
        # fall back to simple peak finding on the raw signal
        return find_peaks_in_window(signal, positions.tolist())
    A = build_convolution_matrix(tract_length, psf, include_first_non_a=True)

    # Extract signal within tract region (plus first non-A)
    tract_start_idx = np.searchsorted(positions, ambiguity_min)
    # +1 to include the first non-A position in the convolution matrix
    tract_end_idx = np.searchsorted(positions, ambiguity_max + 1, side='right')

    if tract_end_idx - tract_start_idx < 2:
        return []

    tract_signal = signal[tract_start_idx:tract_end_idx]
    tract_positions = positions[tract_start_idx:tract_end_idx]

    # Handle size mismatch (warn since it may indicate a boundary issue)
    if len(tract_signal) != A.shape[0]:
        import logging
        logging.getLogger(__name__).warning(
            f"Signal/matrix size mismatch at {ambiguity_min}-{ambiguity_max}: "
            f"signal={len(tract_signal)}, matrix={A.shape[0]}. Truncating."
        )
        n = min(len(tract_signal), A.shape[0])
        A = A[:n, :n]
        tract_signal = tract_signal[:n]
        tract_positions = tract_positions[:n]

    # Skip if no signal
    if tract_signal.sum() < 1.0:
        return []

    # Deconvolve using config regularization
    try:
        deconvolved = deconvolve_signal(tract_signal, A, regularization=NETSEQ_DECONV_REGULARIZATION)
    except Exception:
        # Fall back to observed signal
        deconvolved = tract_signal

    # Identify peaks
    peaks = identify_deconvolved_peaks(
        deconvolved, tract_positions.tolist(), threshold_frac=0.1
    )

    # Add prominence (relative intensity)
    total = sum(p['signal'] for p in peaks) if peaks else 1.0
    for p in peaks:
        p['prominence'] = p['signal'] / total if total > 0 else 0.0

    return peaks


def _proportional_assignment(
    peaks: List[Dict],
    original_position: int,
    use_deconvolution: bool,
) -> List[Dict]:
    """
    Assign reads proportionally to multiple peaks.

    Args:
        peaks: List of peak dicts with 'position' and 'signal'
        original_position: Original position before refinement
        use_deconvolution: Whether deconvolution was used

    Returns:
        List of assignment dicts with position, fraction, confidence
    """
    if not peaks:
        return []

    total_signal = sum(p['signal'] for p in peaks)
    if total_signal == 0:
        total_signal = 1.0

    # Determine confidence based on peak dominance
    max_frac = max(p['signal'] / total_signal for p in peaks)
    if max_frac > 0.9:
        confidence = 'high'
    elif max_frac > 0.7:
        confidence = 'medium'
    else:
        confidence = 'split'

    method = 'deconvolved' if use_deconvolution else 'netseq_peak'

    results = []
    for peak in peaks:
        fraction = peak['signal'] / total_signal
        results.append({
            'assigned_position': peak['position'],
            'fraction': fraction,
            'confidence': confidence,
            'method': method,
            'peak_signal': peak['signal'],
            'n_peaks': len(peaks),
            'shift_from_original': peak['position'] - original_position,
        })

    # Sort by fraction (highest first)
    results.sort(key=lambda x: x['fraction'], reverse=True)

    return results


def refine_batch(
    netseq_loader: NetseqLoader,
    positions: List[Dict],
    use_deconvolution: bool = True,
    psf: Optional[np.ndarray] = None,
    proportional_split: bool = True,
) -> List[Union[Dict, List[Dict]]]:
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
            - count: (optional) Read count for this position
        use_deconvolution: Whether to use NNLS deconvolution
        psf: Point-spread-function for deconvolution
        proportional_split: Return proportional assignments

    Returns:
        List of refinement results (same order as input)
        If proportional_split=True, each result is a List[Dict]
        If proportional_split=False, each result is a Dict
    """
    results = []

    for pos_data in positions:
        result = refine_with_netseq(
            netseq_loader,
            pos_data['chrom'],
            pos_data['ambiguity_min'],
            pos_data['ambiguity_max'],
            pos_data['strand'],
            pos_data['original_position'],
            use_deconvolution=use_deconvolution,
            psf=psf,
            proportional_split=proportional_split,
        )

        # If proportional split and we have a read count, weight the fractions
        if proportional_split and isinstance(result, list):
            read_count = pos_data.get('count', 1.0)
            for r in result:
                r['weighted_count'] = read_count * r['fraction']

        results.append(result)

    return results


def aggregate_proportional_results(
    results: List[List[Dict]],
) -> Dict[Tuple[str, str, int], float]:
    """
    Aggregate proportional assignments into position -> count mapping.

    Args:
        results: List of proportional assignment results from refine_batch

    Returns:
        Dict mapping (chrom, strand, position) -> aggregated weighted count
    """
    from collections import defaultdict
    aggregated = defaultdict(float)

    for result_list in results:
        if not isinstance(result_list, list):
            result_list = [result_list]

        for r in result_list:
            pos = r.get('assigned_position', r.get('refined_position'))
            chrom = r.get('chrom', 'unknown')
            strand = r.get('strand', '+')
            count = r.get('weighted_count', r.get('fraction', 1.0))

            aggregated[(chrom, strand, pos)] += count

    return dict(aggregated)


def calculate_netseq_statistics(
    refinement_results: List[Union[Dict, List[Dict]]],
    flatten: bool = True,
) -> Dict:
    """
    Calculate summary statistics for NET-seq refinement.

    Args:
        refinement_results: List of refinement result dicts or lists (proportional)
        flatten: If True, flatten proportional results before counting

    Returns:
        Dict with summary statistics
    """
    if not refinement_results:
        return {
            'total_input': 0,
            'total_output': 0,
            'refined': 0,
            'refinement_rate': 0.0,
            'mean_shift': 0.0,
            'by_confidence': {},
            'by_method': {},
            'n_proportional': 0,
            'n_multi_peak': 0,
            'total_weighted_count': 0.0,
        }

    # Flatten if needed
    flat_results = []
    n_proportional = 0
    n_multi_peak = 0

    for result in refinement_results:
        if isinstance(result, list):
            n_proportional += 1
            if len(result) > 1:
                n_multi_peak += 1
            if flatten:
                flat_results.extend(result)
            else:
                # Take the first (highest fraction) result
                flat_results.append(result[0] if result else {})
        else:
            flat_results.append(result)

    # Count by confidence
    by_confidence = {}
    for result in flat_results:
        conf = result.get('confidence', 'unknown')
        by_confidence[conf] = by_confidence.get(conf, 0) + 1

    # Count by method
    by_method = {}
    for result in flat_results:
        method = result.get('method', 'unknown')
        by_method[method] = by_method.get(method, 0) + 1

    # Count positions that were refined (changed)
    refined = sum(
        1 for r in flat_results
        if r.get('shift_from_original', 0) != 0
    )

    # Calculate mean shift
    shifts = [abs(r.get('shift_from_original', 0)) for r in flat_results]
    mean_shift = np.mean(shifts) if shifts else 0.0

    # Calculate total weighted counts for proportional assignments
    total_weighted = sum(
        r.get('weighted_count', r.get('fraction', 1.0))
        for r in flat_results
    )

    return {
        'total_input': len(refinement_results),
        'total_output': len(flat_results),
        'refined': refined,
        'refinement_rate': refined / len(flat_results) if flat_results else 0.0,
        'mean_shift': mean_shift,
        'by_confidence': by_confidence,
        'by_method': by_method,
        'n_proportional': n_proportional,
        'n_multi_peak': n_multi_peak,
        'total_weighted_count': total_weighted,
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

    total_in = stats.get('total_input', stats.get('total', 0))
    total_out = stats.get('total_output', stats.get('total', 0))

    report.append("Overall:")
    report.append(f"  Input positions:        {total_in:,}")
    if total_out != total_in:
        report.append(f"  Output assignments:     {total_out:,}")
    report.append(f"  Refined:                {stats['refined']:,} ({stats['refinement_rate']:.1%})")
    report.append(f"  Mean shift:             {stats['mean_shift']:.1f} bp")

    # Proportional splitting stats
    n_prop = stats.get('n_proportional', 0)
    n_multi = stats.get('n_multi_peak', 0)
    if n_prop > 0:
        report.append("")
        report.append("Proportional Assignment:")
        report.append(f"  Proportionally split:   {n_prop:,} ({100.0 * n_prop / total_in:.1f}%)")
        report.append(f"  Multi-peak regions:     {n_multi:,} ({100.0 * n_multi / total_in:.1f}%)")
        if 'total_weighted_count' in stats:
            report.append(f"  Total weighted count:   {stats['total_weighted_count']:,.1f}")

    report.append("")
    report.append("Confidence Levels:")
    for conf in ['high', 'medium', 'split', 'low']:
        count = stats['by_confidence'].get(conf, 0)
        pct = 100.0 * count / total_out if total_out > 0 else 0.0
        if count > 0:
            report.append(f"  {conf:10s}          {count:7,} ({pct:5.1f}%)")
    report.append("")

    report.append("Refinement Methods:")
    for method, count in stats['by_method'].items():
        pct = 100.0 * count / total_out if total_out > 0 else 0.0
        report.append(f"  {method:20s} {count:7,} ({pct:5.1f}%)")

    report.append("")
    report.append("=" * 60)

    return "\n".join(report)
