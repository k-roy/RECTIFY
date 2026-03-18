#!/usr/bin/env python3
"""
A-Tract Deconvolution Module

Deconvolves oligo-adenylation spreading artifacts to recover true CPA positions
within A-tracts using the 0A point-spread-function (PSF).

Key Insight:
At 0A sites (no A-tract), ~54% of signal stays at the true CPA position while
~46% spreads downstream due to poly-A tail alignment. Within A-tracts, this
spreading is constrained by tract boundaries, allowing us to deconvolve the
observed signal to recover true peak positions.

CRITICAL IMPLEMENTATION NOTE:
Signal spreading beyond tract boundaries is LOST, not accumulated at the boundary.
Row sums < 1.0 for boundary positions correctly represent signal loss.
For example, a peak at the last A position has row sum ~0.54 because 46% of
signal spreads beyond the tract and is lost to observation.

DO NOT NORMALIZE ROWS - this is the key fix for correct boundary peak recovery.

Author: Kevin R. Roy
Date: 2026-03-17
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd

# Try to import scipy for NNLS
try:
    from scipy.optimize import nnls
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

# Default PSF based on 0A site analysis
# 54% at position 0, 46% spreading downstream with exponential decay
DEFAULT_PSF = np.array([
    0.54,   # Position 0 (true CPA)
    0.005, 0.015, 0.032, 0.052, 0.067, 0.073, 0.067,  # 1-7 bp downstream
    0.052, 0.038, 0.028, 0.020, 0.015, 0.011, 0.008,  # 8-14 bp downstream
    0.006, 0.005, 0.004, 0.003, 0.002, 0.002,         # 15-20 bp downstream
])

# Normalize to sum to 1
DEFAULT_PSF = DEFAULT_PSF / DEFAULT_PSF.sum()

# Cache for convolution matrices
_MATRIX_CACHE: Dict[int, np.ndarray] = {}


def load_psf(filepath: str) -> np.ndarray:
    """
    Load point-spread-function from file.

    Expected format: TSV with columns 'acount', 'offset', 'probability'
    Extracts the 0A PSF (acount=0).

    Args:
        filepath: Path to PSF file

    Returns:
        PSF array indexed by offset
    """
    df = pd.read_csv(filepath, sep='\t')
    psf_0a = df[df['acount'] == 0].copy()

    max_offset = int(psf_0a['offset'].max())
    psf = np.zeros(max_offset + 1)

    for _, row in psf_0a.iterrows():
        psf[int(row['offset'])] = row['probability']

    return psf


def build_convolution_matrix(
    tract_length: int,
    psf: np.ndarray = None,
    include_first_non_a: bool = True,
) -> np.ndarray:
    """
    Build convolution matrix for A-tract deconvolution.

    A[i, j] = P(observe signal at position j | true peak at position i)

    CRITICAL: Signal spreading beyond tract boundaries is LOST, not accumulated.
    Row sums < 1.0 for boundary positions correctly represent signal loss.

    Args:
        tract_length: Number of A's in the tract
        psf: Point-spread-function (default: use DEFAULT_PSF)
        include_first_non_a: Include the first non-A position after the tract

    Returns:
        Convolution matrix of shape (n_positions, n_positions)
        Note: Row sums may be < 1.0 for positions near boundaries
    """
    if psf is None:
        psf = DEFAULT_PSF

    n = tract_length + (1 if include_first_non_a else 0)
    A = np.zeros((n, n))

    for i in range(n):  # True peak position
        for offset in range(len(psf)):
            j = i + offset  # Observed position

            if j < n:
                # Within bounds - use PSF probability
                A[i, j] = psf[offset]
            # Out of bounds: Signal is LOST, not accumulated at boundary
            # This is the correct physical model - spreading beyond the
            # observable region reduces the observed signal

    # DO NOT NORMALIZE ROWS
    # Row sums < 1.0 at boundaries correctly represent signal loss
    # NNLS will properly account for this when recovering true peaks

    return A


def get_cached_convolution_matrix(
    tract_length: int,
    psf: np.ndarray = None,
) -> np.ndarray:
    """
    Get convolution matrix from cache or compute it.

    Caches matrices for commonly used tract lengths to improve performance.

    Args:
        tract_length: Number of A's in the tract
        psf: Point-spread-function (default: use DEFAULT_PSF)

    Returns:
        Copy of cached convolution matrix
    """
    global _MATRIX_CACHE

    # Use hash of PSF to distinguish different PSFs
    psf_key = hash(psf.tobytes()) if psf is not None else 0
    cache_key = (tract_length, psf_key)

    if cache_key not in _MATRIX_CACHE:
        _MATRIX_CACHE[cache_key] = build_convolution_matrix(tract_length, psf)

    return _MATRIX_CACHE[cache_key].copy()


def clear_matrix_cache():
    """Clear the convolution matrix cache."""
    global _MATRIX_CACHE
    _MATRIX_CACHE.clear()


def deconvolve_signal(
    observed: np.ndarray,
    convolution_matrix: np.ndarray = None,
    tract_length: int = None,
    psf: np.ndarray = None,
    regularization: float = 0.01,
) -> np.ndarray:
    """
    Deconvolve observed signal to recover true peak positions.

    Solves: observed = A @ true_peaks
    Using non-negative least squares with L2 regularization.

    Args:
        observed: Observed signal vector (n_positions,)
        convolution_matrix: Pre-computed convolution matrix (optional)
        tract_length: A-tract length (used if matrix not provided)
        psf: Point-spread-function (used if matrix not provided)
        regularization: L2 regularization strength (default: 0.01)

    Returns:
        Deconvolved signal array (true peak intensities)
    """
    if not SCIPY_AVAILABLE:
        raise ImportError(
            "scipy is required for deconvolution. "
            "Install with: pip install scipy"
        )

    # Build matrix if not provided
    if convolution_matrix is None:
        if tract_length is None:
            raise ValueError("Either convolution_matrix or tract_length must be provided")
        convolution_matrix = build_convolution_matrix(tract_length, psf)

    A = convolution_matrix
    n = len(observed)

    # Ensure dimensions match
    if A.shape[0] != n:
        min_n = min(A.shape[0], n)
        A = A[:min_n, :min_n]
        observed = observed[:min_n]
        n = min_n

    # Add L2 regularization: [A; sqrt(reg)*I] @ x = [observed; 0]
    A_reg = np.vstack([A, np.sqrt(regularization) * np.eye(n)])
    b_reg = np.concatenate([observed, np.zeros(n)])

    # Solve with non-negative least squares
    deconvolved, residual = nnls(A_reg, b_reg)

    return deconvolved


def identify_peaks(
    deconvolved: np.ndarray,
    threshold_fraction: float = 0.1,
    min_intensity: float = 0.0,
) -> List[Tuple[int, float]]:
    """
    Identify true peak positions from deconvolved signal.

    Args:
        deconvolved: Deconvolved signal vector
        threshold_fraction: Minimum fraction of max signal to call a peak
        min_intensity: Minimum absolute intensity to call a peak

    Returns:
        List of (position, intensity) tuples, sorted by intensity descending
    """
    if len(deconvolved) == 0 or deconvolved.max() == 0:
        return []

    threshold = max(threshold_fraction * deconvolved.max(), min_intensity)
    peaks = []

    for i, val in enumerate(deconvolved):
        if val >= threshold:
            peaks.append((i, val))

    # Sort by intensity descending
    peaks.sort(key=lambda x: -x[1])

    return peaks


def compute_fit_quality(
    observed: np.ndarray,
    deconvolved: np.ndarray,
    convolution_matrix: np.ndarray,
) -> Dict[str, float]:
    """
    Compute quality metrics for the deconvolution fit.

    Args:
        observed: Observed signal
        deconvolved: Deconvolved signal
        convolution_matrix: Convolution matrix used

    Returns:
        Dict with quality metrics:
        - r_squared: Coefficient of determination
        - rmse: Root mean squared error
        - signal_recovery: Fraction of observed signal recovered
    """
    predicted = convolution_matrix @ deconvolved

    # R-squared
    residual_ss = np.sum((observed - predicted) ** 2)
    total_ss = np.sum((observed - observed.mean()) ** 2)
    r_squared = 1 - residual_ss / total_ss if total_ss > 0 else 0
    r_squared = max(0, r_squared)

    # RMSE
    rmse = np.sqrt(np.mean((observed - predicted) ** 2))

    # Signal recovery
    signal_recovery = predicted.sum() / observed.sum() if observed.sum() > 0 else 0

    return {
        'r_squared': r_squared,
        'rmse': rmse,
        'signal_recovery': signal_recovery,
    }


def deconvolve_with_confidence(
    observed: np.ndarray,
    tract_length: int,
    psf: np.ndarray = None,
    regularization: float = 0.01,
    threshold_fraction: float = 0.1,
) -> Dict:
    """
    Full deconvolution pipeline with confidence scoring.

    Args:
        observed: Observed signal vector
        tract_length: Length of the A-tract
        psf: Point-spread-function (default: use DEFAULT_PSF)
        regularization: L2 regularization strength
        threshold_fraction: Minimum fraction of max for peak calling

    Returns:
        Dict with:
        - deconvolved: Deconvolved signal array
        - peaks: List of (position, intensity) tuples
        - confidence: Dict of confidence scores per peak
        - fit_quality: Dict with R², RMSE, etc.
        - convolution_matrix: The matrix used
    """
    # Build matrix
    A = build_convolution_matrix(tract_length, psf)

    # Deconvolve
    deconvolved = deconvolve_signal(observed, A, regularization=regularization)

    # Identify peaks
    peaks = identify_peaks(deconvolved, threshold_fraction)

    # Compute fit quality
    fit_quality = compute_fit_quality(observed, deconvolved, A)

    # Compute per-peak confidence
    total_deconv = deconvolved.sum()
    confidence = {}

    for pos, intensity in peaks:
        # Fraction of total signal
        c_intensity = intensity / total_deconv if total_deconv > 0 else 0

        # Position confidence (boundaries tend to have true peaks)
        dist_to_boundary = min(pos, tract_length - pos)
        c_position = 1.0 if pos == tract_length else max(0.3, 1.0 - dist_to_boundary * 0.1)

        # Combined score
        combined = (c_intensity * c_position * max(0.1, fit_quality['r_squared'])) ** 0.33

        # Flag
        if combined >= 0.3:
            flag = 'strong'
        elif combined >= 0.1:
            flag = 'weak'
        else:
            flag = 'marginal'

        confidence[pos] = {
            'intensity_fraction': c_intensity,
            'position_score': c_position,
            'r_squared': fit_quality['r_squared'],
            'combined': combined,
            'flag': flag,
        }

    return {
        'deconvolved': deconvolved,
        'peaks': peaks,
        'confidence': confidence,
        'fit_quality': fit_quality,
        'convolution_matrix': A,
    }


def verify_boundary_recovery():
    """
    Unit test: Verify that boundary peaks are correctly recovered.

    This validates the key fix: no boundary accumulation ensures that
    peaks at tract boundaries are not underestimated.

    Returns:
        True if all tests pass, raises AssertionError otherwise
    """
    print("Testing boundary recovery...")

    # Use default PSF
    psf = DEFAULT_PSF

    # Test: 5-position tract (4A + first non-A)
    tract_length = 4
    A = build_convolution_matrix(tract_length, psf)
    row_sums = A.sum(axis=1)

    # Verify row sums decrease toward boundary
    for i in range(len(row_sums) - 1):
        assert row_sums[i] >= row_sums[i+1], \
            f"Row sums should decrease: {row_sums[i]:.3f} >= {row_sums[i+1]:.3f}"

    # Verify boundary has lowest row sum (approximately psf[0])
    assert row_sums[-1] < row_sums[0], "Boundary should have lowest row sum"
    assert abs(row_sums[-1] - psf[0]) < 0.1, \
        f"Boundary row sum should be ~{psf[0]:.3f}, got {row_sums[-1]:.3f}"

    # Test deconvolution of boundary peak
    true_signal = np.array([0, 0, 0, 0, 1.0])
    observed = A @ true_signal
    deconvolved = deconvolve_signal(observed, A, regularization=0.01)

    # Boundary peak should be recovered approximately
    recovery_error = abs(deconvolved[4] - 1.0)
    assert recovery_error < 0.5, \
        f"Boundary peak recovery error too high: {recovery_error:.3f}"

    print("  ✓ All boundary recovery tests passed")
    return True
