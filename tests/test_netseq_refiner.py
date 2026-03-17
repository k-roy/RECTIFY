#!/usr/bin/env python3
"""
Tests for NET-seq refinement module.
"""

import pytest
import numpy as np
from unittest.mock import Mock, patch, MagicMock
from rectify.core import netseq_refiner
from rectify import config


class TestFindPeaksInWindow:
    """Test peak finding in signal."""

    def test_single_peak(self):
        """Test finding single clear peak."""
        signal = np.array([0.0, 0.1, 0.5, 1.0, 0.5, 0.1, 0.0])
        positions = list(range(1000, 1007))

        peaks = netseq_refiner.find_peaks_in_window(signal, positions, threshold=0.5)

        assert len(peaks) == 1
        assert peaks[0]['position'] == 1003  # Peak at index 3
        assert peaks[0]['signal'] == 1.0

    def test_multiple_peaks(self):
        """Test finding multiple peaks."""
        signal = np.array([0.0, 1.0, 0.0, 0.8, 0.0])
        positions = list(range(1000, 1005))

        peaks = netseq_refiner.find_peaks_in_window(signal, positions, threshold=0.5)

        assert len(peaks) == 2
        assert peaks[0]['position'] == 1001  # First peak
        assert peaks[1]['position'] == 1003  # Second peak

    def test_no_peaks_above_threshold(self):
        """Test signal below threshold."""
        signal = np.array([0.1, 0.15, 0.1, 0.12, 0.1])
        positions = list(range(1000, 1005))

        # Threshold 0.9: 0.15 * 0.9 = 0.135, peak at 0.15 is > 0.135
        # But with threshold > 1.0, nothing passes
        peaks = netseq_refiner.find_peaks_in_window(signal, positions, threshold=1.1)

        assert len(peaks) == 0  # No peaks above 110% of max

    def test_flat_signal(self):
        """Test flat signal (no peaks)."""
        signal = np.array([0.5, 0.5, 0.5, 0.5])
        positions = list(range(1000, 1004))

        peaks = netseq_refiner.find_peaks_in_window(signal, positions, threshold=0.5)

        assert len(peaks) == 0  # No local maxima

    def test_empty_signal(self):
        """Test empty signal."""
        signal = np.array([])
        positions = []

        peaks = netseq_refiner.find_peaks_in_window(signal, positions)

        assert len(peaks) == 0

    def test_zero_signal(self):
        """Test all-zero signal."""
        signal = np.array([0.0, 0.0, 0.0, 0.0])
        positions = list(range(1000, 1004))

        peaks = netseq_refiner.find_peaks_in_window(signal, positions)

        assert len(peaks) == 0  # No signal


class TestSelectBestPeak:
    """Test peak selection."""

    def test_select_highest_signal(self):
        """Test selection of peak with highest signal."""
        peaks = [
            {'position': 1001, 'signal': 0.5, 'prominence': 0.5},
            {'position': 1003, 'signal': 1.0, 'prominence': 0.8},
            {'position': 1005, 'signal': 0.3, 'prominence': 0.3},
        ]

        best = netseq_refiner.select_best_peak(peaks, 1000, 1010)

        assert best['position'] == 1003  # Highest signal

    def test_select_within_window(self):
        """Test selection only considers peaks within window."""
        peaks = [
            {'position': 999, 'signal': 2.0, 'prominence': 2.0},   # Outside window
            {'position': 1003, 'signal': 1.0, 'prominence': 1.0},  # Inside window
            {'position': 1011, 'signal': 2.0, 'prominence': 2.0},  # Outside window
        ]

        best = netseq_refiner.select_best_peak(peaks, 1000, 1010)

        assert best['position'] == 1003  # Only peak in window

    def test_no_peaks_in_window(self):
        """Test when no peaks fall within window."""
        peaks = [
            {'position': 999, 'signal': 1.0, 'prominence': 1.0},
            {'position': 1011, 'signal': 1.0, 'prominence': 1.0},
        ]

        best = netseq_refiner.select_best_peak(peaks, 1000, 1010)

        assert best is None

    def test_empty_peaks(self):
        """Test with no peaks."""
        best = netseq_refiner.select_best_peak([], 1000, 1010)

        assert best is None


class TestAssignConfidence:
    """Test confidence assignment."""

    def test_high_confidence_single_peak(self):
        """Test high confidence for single strong peak."""
        peak = {'position': 1003, 'signal': 1.5, 'prominence': 1.5}
        all_peaks = [peak]

        confidence, method = netseq_refiner.assign_confidence(
            peak, all_peaks, ambiguity_range=5
        )

        assert confidence == 'high'
        assert method == 'netseq_peak'

    def test_medium_confidence_multiple_peaks(self):
        """Test medium confidence for multiple peaks."""
        peak = {'position': 1003, 'signal': 1.5, 'prominence': 1.5}
        all_peaks = [
            peak,
            {'position': 1005, 'signal': 1.2, 'prominence': 1.0},
        ]

        confidence, method = netseq_refiner.assign_confidence(
            peak, all_peaks, ambiguity_range=5
        )

        assert confidence == 'medium'
        assert method == 'netseq_peak'

    def test_low_confidence_weak_signal(self):
        """Test low confidence for weak signal."""
        peak = {'position': 1003, 'signal': 0.3, 'prominence': 0.2}
        all_peaks = [peak]

        confidence, method = netseq_refiner.assign_confidence(
            peak, all_peaks, ambiguity_range=5
        )

        assert confidence == 'low'
        assert method == 'netseq_peak'

    def test_no_peak(self):
        """Test confidence when no peak found."""
        confidence, method = netseq_refiner.assign_confidence(
            None, [], ambiguity_range=5
        )

        assert confidence == 'low'
        assert method == 'no_data'


class TestNetseqLoader:
    """Test NET-seq BigWig loader."""

    @patch('rectify.core.netseq_refiner.PYBIGWIG_AVAILABLE', True)
    @patch('rectify.core.netseq_refiner.pyBigWig')
    def test_load_bigwig(self, mock_pybigwig):
        """Test loading a BigWig file."""
        # Mock pyBigWig.open()
        mock_bw = MagicMock()
        mock_pybigwig.open.return_value = mock_bw

        loader = netseq_refiner.NetseqLoader()
        loader.load_bigwig('/path/to/file.bw', name='test')

        assert 'test' in loader.bigwigs
        mock_pybigwig.open.assert_called_once()

    def test_get_signal(self):
        """Test getting signal from BigWig."""
        # Mock BigWig object
        mock_bw = MagicMock()
        mock_bw.values.return_value = [0.0, 1.0, 2.0, 3.0, 4.0]

        loader = netseq_refiner.NetseqLoader()
        loader.bigwigs['test'] = mock_bw

        signal = loader.get_signal('chrI', 1000, 1005, '+')

        assert len(signal) == 5
        assert signal[4] == 4.0
        mock_bw.values.assert_called_once_with('chrI', 1000, 1005)

    def test_get_signal_with_nones(self):
        """Test handling None values in BigWig."""
        mock_bw = MagicMock()
        mock_bw.values.return_value = [0.0, None, 2.0, None, 4.0]

        loader = netseq_refiner.NetseqLoader()
        loader.bigwigs['test'] = mock_bw

        signal = loader.get_signal('chrI', 1000, 1005, '+')

        # None values should be replaced with 0
        assert signal[1] == 0.0
        assert signal[3] == 0.0

    def test_signal_caching(self):
        """Test signal caching."""
        mock_bw = MagicMock()
        mock_bw.values.return_value = [1.0, 2.0, 3.0]

        loader = netseq_refiner.NetseqLoader()
        loader.bigwigs['test'] = mock_bw

        # First call - should hit BigWig
        signal1 = loader.get_signal('chrI', 1000, 1003, '+')

        # Second call - should use cache
        signal2 = loader.get_signal('chrI', 1000, 1003, '+')

        # Should only call values() once
        assert mock_bw.values.call_count == 1
        np.testing.assert_array_equal(signal1, signal2)

    def test_combine_multiple_bigwigs(self):
        """Test combining signal from multiple BigWigs."""
        mock_bw1 = MagicMock()
        mock_bw1.values.return_value = [1.0, 2.0, 3.0]

        mock_bw2 = MagicMock()
        mock_bw2.values.return_value = [0.5, 1.0, 1.5]

        loader = netseq_refiner.NetseqLoader()
        loader.bigwigs['bw1'] = mock_bw1
        loader.bigwigs['bw2'] = mock_bw2

        signal = loader.get_signal('chrI', 1000, 1003, '+')

        # Should combine both signals
        expected = np.array([1.5, 3.0, 4.5])
        np.testing.assert_array_equal(signal, expected)


class TestRefineWithNetseq:
    """Test full refinement workflow."""

    def test_refine_with_clear_peak(self):
        """Test refinement with clear NET-seq peak."""
        # Mock loader
        loader = Mock()
        # Window: 1000-1005, expanded to 995-1010 = 15 positions
        # Signal with clear peak at index 8 (position 1003)
        signal = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 2.0, 1.0, 0.8, 0.5, 0.3, 0.2, 0.1, 0.0])
        loader.get_signal.return_value = signal

        result = netseq_refiner.refine_with_netseq(
            loader,
            chrom='chrI',
            ambiguity_min=1000,
            ambiguity_max=1005,
            strand='+',
            original_position=1005
        )

        assert result['refined_position'] == 1002  # Peak position (995 + 7 = 1002)
        assert result['confidence'] in ['high', 'medium']
        assert result['method'] == 'netseq_peak'
        assert result['n_peaks'] >= 1

    def test_refine_no_signal(self):
        """Test refinement with no NET-seq signal."""
        loader = Mock()
        # Window: 1000-1005, expanded to 995-1010 = 15 positions
        loader.get_signal.return_value = np.zeros(15)  # All zeros

        result = netseq_refiner.refine_with_netseq(
            loader,
            chrom='chrI',
            ambiguity_min=1000,
            ambiguity_max=1005,
            strand='+',
            original_position=1005
        )

        # Should fall back to leftmost position
        assert result['refined_position'] == 1000
        assert result['confidence'] == 'low'
        assert result['method'] == 'leftmost'


class TestBatchProcessing:
    """Test batch refinement."""

    def test_refine_batch(self):
        """Test batch refinement of multiple positions."""
        loader = Mock()
        # Return appropriate signal size for each call (15 positions)
        loader.get_signal.return_value = np.array([0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0], dtype=float)

        positions = [
            {
                'chrom': 'chrI',
                'ambiguity_min': 1000,
                'ambiguity_max': 1005,
                'strand': '+',
                'original_position': 1005,
            },
            {
                'chrom': 'chrII',
                'ambiguity_min': 2000,
                'ambiguity_max': 2005,
                'strand': '+',
                'original_position': 2005,
            },
        ]

        results = netseq_refiner.refine_batch(loader, positions)

        assert len(results) == 2
        assert 'refined_position' in results[0]
        assert 'confidence' in results[0]


class TestStatistics:
    """Test NET-seq refinement statistics."""

    def test_calculate_statistics(self):
        """Test statistics calculation."""
        results = [
            {
                'refined_position': 1003,
                'confidence': 'high',
                'method': 'netseq_peak',
                'shift_from_original': -2,
            },
            {
                'refined_position': 2005,
                'confidence': 'medium',
                'method': 'netseq_peak',
                'shift_from_original': 0,
            },
            {
                'refined_position': 3001,
                'confidence': 'low',
                'method': 'leftmost',
                'shift_from_original': -3,
            },
        ]

        stats = netseq_refiner.calculate_netseq_statistics(results)

        assert stats['total'] == 3
        assert stats['refined'] == 2  # Two with non-zero shift
        assert stats['refinement_rate'] == pytest.approx(2/3)
        assert stats['by_confidence']['high'] == 1
        assert stats['by_method']['netseq_peak'] == 2

    def test_statistics_empty(self):
        """Test statistics with empty input."""
        stats = netseq_refiner.calculate_netseq_statistics([])

        assert stats['total'] == 0

    def test_format_report(self):
        """Test report formatting."""
        stats = {
            'total': 100,
            'refined': 85,
            'refinement_rate': 0.85,
            'mean_shift': 1.5,
            'by_confidence': {'high': 50, 'medium': 30, 'low': 20},
            'by_method': {'netseq_peak': 85, 'leftmost': 15},
        }

        report = netseq_refiner.format_netseq_report(stats)

        assert 'NET-seq Refinement' in report
        assert '100' in report
        assert '85.0%' in report
        assert '1.5 bp' in report


class TestEdgeCases:
    """Test edge cases."""

    def test_mismatched_signal_positions(self):
        """Test error when signal and positions length mismatch."""
        signal = np.array([1.0, 2.0, 3.0])
        positions = [1000, 1001]  # Wrong length

        with pytest.raises(ValueError):
            netseq_refiner.find_peaks_in_window(signal, positions)

    def test_ambiguity_window_outside_signal(self):
        """Test when ambiguity window is very large."""
        loader = Mock()
        # Window: 1000-2000, expanded to 995-2005 = 1010 positions
        # Provide signal of correct size
        loader.get_signal.return_value = np.zeros(1010)

        # This should work even with large window
        result = netseq_refiner.refine_with_netseq(
            loader,
            chrom='chrI',
            ambiguity_min=1000,
            ambiguity_max=2000,  # Very large window
            strand='+',
            original_position=1500
        )

        assert 'refined_position' in result
        # With no signal, should use leftmost
        assert result['refined_position'] == 1000
