#!/usr/bin/env python3
"""
Pan-Mutant NET-seq Refiner for A-Tract Ambiguity Resolution

This module uses pre-computed pan-mutant NET-seq distributions to refine
nanopore 3' end positions that fall within ambiguous A-tracts.

Key Concept:
Nanopore oligo-dT priming captures poly-A tails, but when the CPA site falls
within a genomic A-tract, the read's alignment becomes ambiguous - it could
map to any position within the tract. NET-seq data (which captures nascent
RNA 3' ends without poly-A selection bias) provides ground truth for the
true distribution of termination sites within A-tracts.

Algorithm:
1. Load pan-mutant NET-seq database (HDF5 from build_pan_mutant_netseq_database.py)
2. For each nanopore read falling within an A-tract:
   a. Look up the corresponding NET-seq signal distribution
   b. Apply NNLS deconvolution if needed
   c. Assign the read proportionally to peaks OR winner-take-all

Usage:
    from rectify.core.analyze.pan_mutant_refiner import PanMutantRefiner

    refiner = PanMutantRefiner('/path/to/pan_mutant_netseq_database.h5')
    refined = refiner.refine_position(
        chrom='chrXII',
        position=450000,
        strand='-',
        tract_start=450000,
        tract_end=450008,
    )

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-24
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np

# Try to import h5py
try:
    import h5py
    H5PY_AVAILABLE = True
except ImportError:
    h5py = None
    H5PY_AVAILABLE = False

# Try to import scipy for NNLS
try:
    from scipy.optimize import nnls
    SCIPY_AVAILABLE = True
except ImportError:
    nnls = None
    SCIPY_AVAILABLE = False

from ...config import (
    NETSEQ_PSF_POSITION_0_FRACTION,
    NETSEQ_DECONV_REGULARIZATION,
)


class PanMutantRefiner:
    """
    Refines nanopore 3' end positions using pan-mutant NET-seq data.

    This class loads a pre-computed HDF5 database of pan-mutant NET-seq
    signal distributions at A-tract CPA sites and uses them to resolve
    positional ambiguity in nanopore data.
    """

    def __init__(
        self,
        database_path: str,
        use_deconvolution: bool = True,
        regularization: float = None,
        cache_size: int = 10000,
    ):
        """
        Initialize the pan-mutant refiner.

        Args:
            database_path: Path to HDF5 database from build_pan_mutant_netseq_database.py
            use_deconvolution: Whether to apply NNLS deconvolution (recommended)
            regularization: L2 regularization for NNLS (default: from config)
            cache_size: Maximum cached signal arrays
        """
        if not H5PY_AVAILABLE:
            raise ImportError("h5py is required for pan-mutant refinement. "
                             "Install with: pip install h5py")

        self.database_path = Path(database_path)
        self.use_deconvolution = use_deconvolution
        self.regularization = regularization or NETSEQ_DECONV_REGULARIZATION
        self.cache_size = cache_size

        # Load database
        self._load_database()

        # Signal cache (LRU-style)
        self._signal_cache = {}
        self._cache_order = []

    def _load_database(self):
        """Load metadata from HDF5 database."""
        self.h5f = h5py.File(self.database_path, 'r')

        # Load metadata
        metadata = self.h5f['metadata']
        self.cluster_ids = metadata['cluster_id'][:]
        self.chroms = [c.decode() if isinstance(c, bytes) else c
                       for c in metadata['chrom'][:]]
        self.strands = [s.decode() if isinstance(s, bytes) else s
                        for s in metadata['strand'][:]]
        self.positions = metadata['position'][:]
        self.total_reads = metadata['total_reads'][:]

        # Build position index for fast lookup
        # Key: (chrom, strand, position) -> cluster_id
        self._position_index = {}
        for i, (chrom, strand, pos, cid) in enumerate(
            zip(self.chroms, self.strands, self.positions, self.cluster_ids)
        ):
            self._position_index[(chrom, strand, pos)] = str(cid)

        # Load aggregated PSF if available
        if 'aggregated_psf' in self.h5f and 'all_sites' in self.h5f['aggregated_psf']:
            self.aggregated_psf = self.h5f['aggregated_psf']['all_sites'][:]
        else:
            self.aggregated_psf = None

        # Store database attributes
        self.n_clusters = self.h5f.attrs.get('n_clusters', len(self.cluster_ids))
        self.window_extension = self.h5f.attrs.get('window_extension', 10)

    def _get_signal(self, cluster_id: str) -> Optional[np.ndarray]:
        """
        Get signal array for a cluster (with caching).

        Args:
            cluster_id: Cluster ID as string

        Returns:
            Signal array or None if not found
        """
        # Check cache
        if cluster_id in self._signal_cache:
            # Move to end of cache order (most recently used)
            self._cache_order.remove(cluster_id)
            self._cache_order.append(cluster_id)
            return self._signal_cache[cluster_id]

        # Load from HDF5
        if cluster_id in self.h5f['signals']:
            signal = self.h5f['signals'][cluster_id][:]

            # Cache with eviction
            self._signal_cache[cluster_id] = signal
            self._cache_order.append(cluster_id)

            if len(self._cache_order) > self.cache_size:
                # Remove oldest entry
                oldest = self._cache_order.pop(0)
                del self._signal_cache[oldest]

            return signal

        return None

    def find_cluster(
        self,
        chrom: str,
        position: int,
        strand: str,
        search_window: int = 5,
    ) -> Optional[str]:
        """
        Find the cluster ID for a given position.

        Searches for exact match first, then within a window.

        Args:
            chrom: Chromosome name
            position: Genomic position
            strand: Strand ('+' or '-')
            search_window: Search window for nearby clusters

        Returns:
            Cluster ID as string, or None if not found
        """
        # Exact match
        key = (chrom, strand, position)
        if key in self._position_index:
            return self._position_index[key]

        # Search nearby
        for offset in range(1, search_window + 1):
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
        tract_start: Optional[int] = None,
        tract_end: Optional[int] = None,
        proportional: bool = True,
        threshold_frac: float = 0.1,
    ) -> Union[Dict, List[Dict]]:
        """
        Refine a nanopore 3' end position using pan-mutant NET-seq data.

        Args:
            chrom: Chromosome name
            position: Current estimated 3' end position
            strand: Gene strand
            tract_start: Start of A-tract (if known)
            tract_end: End of A-tract (if known)
            proportional: Return proportional assignments (True) or winner-take-all (False)
            threshold_frac: Minimum fraction of max signal to call a peak

        Returns:
            If proportional=True:
                List of dicts with:
                - position: Refined position
                - fraction: Probability of this position
                - signal: NET-seq signal at position
                - method: 'pan_mutant_deconv' or 'pan_mutant_raw'

            If proportional=False:
                Single dict with:
                - position: Best estimate position
                - confidence: 'high', 'medium', 'low'
                - signal: NET-seq signal at position
                - shift: Position change from input
        """
        # Find corresponding cluster
        cluster_id = self.find_cluster(chrom, position, strand)

        if cluster_id is None:
            # No data available - return original position with low confidence
            if proportional:
                return [{
                    'position': position,
                    'fraction': 1.0,
                    'signal': 0.0,
                    'method': 'no_data',
                }]
            else:
                return {
                    'position': position,
                    'confidence': 'low',
                    'signal': 0.0,
                    'shift': 0,
                    'method': 'no_data',
                }

        # Get signal
        signal = self._get_signal(cluster_id)

        if signal is None or signal.sum() == 0:
            if proportional:
                return [{
                    'position': position,
                    'fraction': 1.0,
                    'signal': 0.0,
                    'method': 'no_signal',
                }]
            else:
                return {
                    'position': position,
                    'confidence': 'low',
                    'signal': 0.0,
                    'shift': 0,
                    'method': 'no_signal',
                }

        # Apply deconvolution if enabled
        if self.use_deconvolution and SCIPY_AVAILABLE and self.aggregated_psf is not None:
            refined_signal = self._deconvolve_signal(signal)
            method = 'pan_mutant_deconv'
        else:
            refined_signal = signal
            method = 'pan_mutant_raw'

        # Find peaks
        center_idx = len(refined_signal) // 2
        center_pos = self.positions[np.where(self.cluster_ids == int(cluster_id))[0][0]]

        # Compute position for each index
        positions = np.arange(
            center_pos - self.window_extension,
            center_pos + self.window_extension + 1
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
                    'index': i,
                })

        if not peaks:
            # No significant peaks - use maximum
            max_idx = np.argmax(refined_signal)
            peaks = [{
                'position': int(positions[max_idx]),
                'signal': float(refined_signal[max_idx]),
                'index': max_idx,
            }]

        # Compute fractions
        total_signal = sum(p['signal'] for p in peaks)
        for peak in peaks:
            peak['fraction'] = peak['signal'] / total_signal if total_signal > 0 else 0

        if proportional:
            # Return all peaks with their fractions
            results = []
            for peak in sorted(peaks, key=lambda p: -p['fraction']):
                results.append({
                    'position': peak['position'],
                    'fraction': peak['fraction'],
                    'signal': peak['signal'],
                    'method': method,
                    'shift': peak['position'] - position,
                })
            return results
        else:
            # Return best peak (winner-take-all)
            best_peak = max(peaks, key=lambda p: p['signal'])

            # Determine confidence
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
                'shift': best_peak['position'] - position,
                'method': method,
                'n_peaks': len(peaks),
            }

    def _deconvolve_signal(self, signal: np.ndarray) -> np.ndarray:
        """
        Apply NNLS deconvolution to remove spreading artifact.

        Args:
            signal: Raw signal array

        Returns:
            Deconvolved signal array
        """
        if self.aggregated_psf is None:
            return signal

        n = len(signal)
        psf = self.aggregated_psf

        # Build convolution matrix
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(len(psf)):
                if i + j < n:
                    A[i, i + j] = psf[j]

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
            positions: List of dicts with keys:
                - chrom: Chromosome
                - position: Current position
                - strand: Strand
                - tract_start: Optional A-tract start
                - tract_end: Optional A-tract end
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
                tract_start=pos_data.get('tract_start'),
                tract_end=pos_data.get('tract_end'),
                proportional=proportional,
            )
            results.append(result)

        return results

    def get_statistics(self) -> Dict:
        """Get database statistics."""
        return {
            'n_clusters': self.n_clusters,
            'window_extension': self.window_extension,
            'total_reads': int(self.total_reads.sum()),
            'mean_reads_per_cluster': float(self.total_reads.mean()),
            'has_aggregated_psf': self.aggregated_psf is not None,
            'use_deconvolution': self.use_deconvolution,
            'regularization': self.regularization,
        }

    def close(self):
        """Close HDF5 file."""
        try:
            if hasattr(self, 'h5f') and self.h5f is not None:
                self.h5f.close()
                self.h5f = None
        except Exception:
            pass  # Ignore errors during cleanup

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def __del__(self):
        self.close()
