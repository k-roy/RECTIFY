#!/usr/bin/env python3
"""
A-Tract Refiner Module

This module provides the user-facing API for refining nanopore 3' end positions
within A-tracts using the bundled NET-seq reference data.

Usage:
    from rectify.core.analyze.atract_refiner import ATractRefiner

    # Initialize (automatically loads bundled reference)
    refiner = ATractRefiner()

    # Refine a single position
    result = refiner.refine_position('chrI', 6281, '+')
    # Returns: [{'position': 6291, 'fraction': 0.41, ...}, ...]

    # Batch refinement
    results = refiner.refine_batch([
        {'chrom': 'chrI', 'position': 6281, 'strand': '+'},
        {'chrom': 'chrII', 'position': 150000, 'strand': '-'},
    ])

The module uses pre-computed NET-seq signal distributions from 109 GSE159603
mutant samples to resolve positional ambiguity in nanopore data.

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-24
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import importlib.resources
import gzip

import numpy as np
import pandas as pd

from ...config import (
    NETSEQ_PSF_POSITION_0_FRACTION,
    NETSEQ_DECONV_REGULARIZATION,
)


# Try to import scipy for NNLS
try:
    from scipy.optimize import nnls
    SCIPY_AVAILABLE = True
except ImportError:
    nnls = None
    SCIPY_AVAILABLE = False


# Default bundled reference filename
BUNDLED_REFERENCE = 'saccharomyces_cerevisiae_atract_netseq.tsv.gz'

# Window size used in reference (must match create_atract_netseq_reference.py)
WINDOW_SIZE = 10


class ATractRefiner:
    """
    Refines nanopore 3' end positions using bundled NET-seq A-tract reference.

    This class provides the primary user-facing API for A-tract refinement.
    It automatically loads the bundled reference data on initialization.
    """

    def __init__(
        self,
        reference_path: Optional[str] = None,
        use_deconvolution: bool = True,
        regularization: float = None,
    ):
        """
        Initialize the A-tract refiner.

        Args:
            reference_path: Path to A-tract reference TSV (default: bundled)
            use_deconvolution: Whether to apply NNLS deconvolution
            regularization: L2 regularization for NNLS (default: from config)
        """
        self.use_deconvolution = use_deconvolution
        self.regularization = regularization or NETSEQ_DECONV_REGULARIZATION

        # Load reference data
        if reference_path is None:
            reference_path = self._get_bundled_reference_path()

        self._load_reference(reference_path)

    def _get_bundled_reference_path(self) -> str:
        """Get path to bundled reference file."""
        try:
            # Python 3.9+
            with importlib.resources.files('rectify.data').joinpath(BUNDLED_REFERENCE) as path:
                return str(path)
        except AttributeError:
            # Python 3.8 fallback
            import pkg_resources
            return pkg_resources.resource_filename('rectify.data', BUNDLED_REFERENCE)

    def _load_reference(self, filepath: str):
        """Load A-tract reference from file."""
        self.reference_df = pd.read_csv(filepath, sep='\t')

        # Build position index for fast lookup
        # Key: (chrom, strand, position) -> row index
        self._position_index = {}
        for idx, row in self.reference_df.iterrows():
            key = (row['chrom'], row['strand'], row['position'])
            self._position_index[key] = idx

        # Parse signal windows into arrays
        self._signal_cache = {}
        for idx, row in self.reference_df.iterrows():
            signals = [float(x) for x in row['signal_window'].split(',')]
            self._signal_cache[idx] = np.array(signals)

        # Build PSF from config
        self._build_psf()

    def _build_psf(self):
        """Build point-spread-function for deconvolution."""
        # Use same PSF as other modules
        self.psf = np.array([
            NETSEQ_PSF_POSITION_0_FRACTION,  # Position 0 (~54%)
            0.005, 0.015, 0.032, 0.052, 0.067, 0.073, 0.067,  # 1-7
            0.052, 0.038, 0.028, 0.020, 0.015, 0.011, 0.008,  # 8-14
            0.006, 0.005, 0.004, 0.003, 0.002, 0.002,         # 15-20
        ])
        self.psf = self.psf / self.psf.sum()

    def find_nearest_site(
        self,
        chrom: str,
        position: int,
        strand: str,
        max_distance: int = 50,
    ) -> Optional[int]:
        """
        Find nearest A-tract site in reference.

        Args:
            chrom: Chromosome name
            position: Query position
            strand: Strand ('+' or '-')
            max_distance: Maximum search distance

        Returns:
            Row index of nearest site, or None if not found
        """
        # Exact match
        key = (chrom, strand, position)
        if key in self._position_index:
            return self._position_index[key]

        # Search nearby
        for offset in range(1, max_distance + 1):
            for pos in [position - offset, position + offset]:
                key = (chrom, strand, pos)
                if key in self._position_index:
                    return self._position_index[key]

        return None

    def refine_position(
        self,
        chrom: str,
        position: int,
        strand: str,
        proportional: bool = True,
        threshold_frac: float = 0.1,
    ) -> Union[Dict, List[Dict]]:
        """
        Refine a nanopore 3' end position using NET-seq reference.

        Args:
            chrom: Chromosome name
            position: Current estimated 3' end position
            strand: Gene strand ('+' or '-')
            proportional: Return proportional assignments (True) or winner-take-all
            threshold_frac: Minimum fraction of max signal to call a peak

        Returns:
            If proportional=True:
                List of dicts with position, fraction, signal, shift

            If proportional=False:
                Single dict with position, confidence, shift
        """
        # Find nearest reference site
        site_idx = self.find_nearest_site(chrom, position, strand)

        if site_idx is None:
            # No reference data - return original position
            if proportional:
                return [{
                    'position': position,
                    'fraction': 1.0,
                    'signal': 0.0,
                    'method': 'no_reference',
                    'shift': 0,
                }]
            else:
                return {
                    'position': position,
                    'confidence': 'low',
                    'signal': 0.0,
                    'method': 'no_reference',
                    'shift': 0,
                }

        # Get signal and reference position
        row = self.reference_df.iloc[site_idx]
        ref_position = row['position']
        signal = self._signal_cache[site_idx]

        # Apply deconvolution if enabled
        if self.use_deconvolution and SCIPY_AVAILABLE:
            refined_signal = self._deconvolve_signal(signal)
            method = 'atract_deconv'
        else:
            refined_signal = signal
            method = 'atract_raw'

        # Find peaks
        positions = np.arange(
            ref_position - WINDOW_SIZE,
            ref_position + WINDOW_SIZE + 1
        )

        # Identify significant peaks
        max_signal = refined_signal.max()
        threshold = threshold_frac * max_signal

        peaks = []
        for i, (pos, sig) in enumerate(zip(positions, refined_signal)):
            if sig >= threshold:
                peaks.append({
                    'position': int(pos),
                    'signal': float(sig),
                })

        if not peaks:
            # No significant peaks - use maximum
            max_idx = np.argmax(refined_signal)
            peaks = [{
                'position': int(positions[max_idx]),
                'signal': float(refined_signal[max_idx]),
            }]

        # Compute fractions
        total_signal = sum(p['signal'] for p in peaks)
        for peak in peaks:
            peak['fraction'] = peak['signal'] / total_signal if total_signal > 0 else 0

        if proportional:
            results = []
            for peak in sorted(peaks, key=lambda p: -p['fraction']):
                results.append({
                    'position': peak['position'],
                    'fraction': peak['fraction'],
                    'signal': peak['signal'],
                    'method': method,
                    'shift': peak['position'] - position,
                    'downstream_acount': row['downstream_acount'],
                })
            return results
        else:
            best_peak = max(peaks, key=lambda p: p['signal'])
            best_frac = best_peak['fraction']

            if best_frac > 0.8:
                confidence = 'high'
            elif best_frac > 0.5:
                confidence = 'medium'
            else:
                confidence = 'low'

            return {
                'position': best_peak['position'],
                'confidence': confidence,
                'signal': best_peak['signal'],
                'method': method,
                'shift': best_peak['position'] - position,
                'downstream_acount': row['downstream_acount'],
                'n_peaks': len(peaks),
            }

    def _deconvolve_signal(self, signal: np.ndarray) -> np.ndarray:
        """Apply NNLS deconvolution to remove spreading artifact."""
        n = len(signal)

        # Build convolution matrix
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(len(self.psf)):
                if i + j < n:
                    A[i, i + j] = self.psf[j]

        # Add regularization
        A_reg = np.vstack([A, np.sqrt(self.regularization) * np.eye(n)])
        b_reg = np.concatenate([signal, np.zeros(n)])

        # Solve NNLS
        try:
            deconvolved, _ = nnls(A_reg, b_reg)
            return deconvolved
        except Exception:
            return signal

    def refine_batch(
        self,
        positions: List[Dict],
        proportional: bool = True,
    ) -> List[Union[Dict, List[Dict]]]:
        """
        Refine a batch of positions.

        Args:
            positions: List of dicts with chrom, position, strand
            proportional: Return proportional assignments

        Returns:
            List of refinement results (same order as input)
        """
        results = []
        for pos_data in positions:
            result = self.refine_position(
                chrom=pos_data['chrom'],
                position=pos_data['position'],
                strand=pos_data['strand'],
                proportional=proportional,
            )
            results.append(result)
        return results

    def get_statistics(self) -> Dict:
        """Get reference statistics."""
        return {
            'n_sites': len(self.reference_df),
            'total_signal': float(self.reference_df['total_signal'].sum()),
            'mean_signal_per_site': float(self.reference_df['total_signal'].mean()),
            'acount_distribution': self.reference_df['downstream_acount'].value_counts().to_dict(),
            'use_deconvolution': self.use_deconvolution,
            'regularization': self.regularization,
        }


# Convenience function for quick refinement
def refine_atract_position(
    chrom: str,
    position: int,
    strand: str,
    proportional: bool = True,
) -> Union[Dict, List[Dict]]:
    """
    Quick refinement of a single position (creates temporary refiner).

    For batch refinement, create an ATractRefiner instance instead.

    Args:
        chrom: Chromosome name
        position: Current position
        strand: Strand ('+' or '-')
        proportional: Return proportional assignments

    Returns:
        Refinement result(s)
    """
    refiner = ATractRefiner()
    return refiner.refine_position(chrom, position, strand, proportional)
