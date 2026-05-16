#!/usr/bin/env python3
"""
Tests for AG mispriming detection module.
"""

import pytest
from rectify.core.polya import ag_mispriming
from rectify import config


# Mock genome for testing
MOCK_GENOME = {
    'ref|NC_001133|': (  # chrI
        'N' * 1000 +
        # Position 1000-1049: Low AG content (~30%)
        'TTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTT' +
        # Position 1050-1099: Medium AG content (~60%)
        'AGAGCTCTAG' * 5 +
        # Position 1100-1149: High AG content (~80%)
        'AGAGAGAGAG' * 5 +
        # Position 1150-1199: Very high AG content (~90%)
        'AAAGGGAAAG' * 5 +
        'N' * 10000
    )
}


class TestAGContentCalculation:
    """Test AG content calculation."""

    def test_calculate_ag_content_zero(self):
        """Test sequence with no A/G bases."""
        assert ag_mispriming.calculate_ag_content('TTCCTTCC') == 0.0

    def test_calculate_ag_content_half(self):
        """Test sequence with 50% AG."""
        assert ag_mispriming.calculate_ag_content('AGTC' * 10) == 0.5

    def test_calculate_ag_content_full(self):
        """Test sequence with 100% AG."""
        assert ag_mispriming.calculate_ag_content('AGAGAGAG') == 1.0

    def test_calculate_ag_content_mixed(self):
        """Test mixed sequence."""
        # AGATC = A, G, A = 3 A's/G's out of 5 = 0.6
        assert ag_mispriming.calculate_ag_content('AGATC') == pytest.approx(0.6)

    def test_calculate_ag_content_case_insensitive(self):
        """Test case insensitivity."""
        assert ag_mispriming.calculate_ag_content('aGaGaG') == 1.0

    def test_calculate_ag_content_empty(self):
        """Test empty sequence."""
        assert ag_mispriming.calculate_ag_content('') == 0.0


class TestAGMisprimingScreening:
    """Test AG mispriming screening."""

    def test_low_ag_content(self):
        """Test position with low AG content downstream (not misprimed)."""
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrI', 1000, '+', window=50
        )

        assert result['window_size'] == 50
        assert result['ag_content'] < 0.5
        assert not result['is_likely_misprimed']
        assert result['ag_score'] <= config.AG_RICHNESS_SCORE_THRESHOLD
        assert not result['insufficient_data']

    def test_medium_ag_content(self):
        """Position 1050: ~60% AG content, weighted score above threshold."""
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrI', 1050, '+', window=50
        )

        assert result['window_size'] == 50
        assert 0.5 < result['ag_content'] < 0.7
        assert result['ag_score'] > config.AG_RICHNESS_SCORE_THRESHOLD
        assert result['is_likely_misprimed']

    def test_high_ag_content(self):
        """Test position with high AG content (likely misprimed)."""
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrI', 1100, '+', window=50
        )

        assert result['window_size'] == 50
        assert result['ag_content'] >= 0.75
        assert result['ag_score'] > config.AG_RICHNESS_SCORE_THRESHOLD
        assert result['is_likely_misprimed']

    def test_very_high_ag_content(self):
        """Test position with very high AG content."""
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrI', 1150, '+', window=50
        )

        assert result['window_size'] == 50
        assert result['ag_content'] >= 0.85
        assert result['ag_score'] > config.AG_RICHNESS_SCORE_THRESHOLD
        assert result['is_likely_misprimed']

    def test_custom_threshold(self):
        """Test with custom AG-richness threshold."""
        # Position with ~60% AG content
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrI', 1050, '+', window=50, threshold=0.55
        )

        # Should be flagged with lower threshold
        assert result['is_likely_misprimed']

    def test_minus_strand(self):
        """Test AG screening on minus strand."""
        # Minus strand: downstream in gene coords = LEFT in genomic coords
        # Position 1149 looking LEFT (positions 1099-1148) has genomic 'AGAGAGAG...'
        # When reverse complemented: A→T, G→C, so this becomes 'TCTCTCTC...' (0% AG in RNA)
        # For high AG in RNA orientation, we need high T/C in genomic sequence
        # Position 1049 looking LEFT (positions 999-1048) has 'TTCCTTCC...' (~60% TC)
        # Reverse complement: T→A, C→G, so ~60% AG in RNA orientation
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrI', 1049, '-', window=50
        )

        assert result['window_size'] == 50
        # Should have ~60% AG content (from reverse complementing TC-rich region)
        assert result['ag_content'] >= 0.55
        # May or may not be flagged depending on exact threshold

    def test_base_composition(self):
        """Test base composition reporting."""
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrI', 1100, '+', window=50
        )

        comp = result['base_composition']
        assert 'A' in comp
        assert 'G' in comp
        # Position 1100-1149 is 'AGAGAGAGAG' * 5, so only A's and G's
        assert comp['A'] == 25
        assert comp['G'] == 25
        # Check total matches window size
        assert sum(comp.values()) == 50


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_near_chromosome_end(self):
        """Test position near chromosome end (insufficient window)."""
        # Position near end where we can't get full 50bp window
        chrom_len = len(MOCK_GENOME['ref|NC_001133|'])
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrI', chrom_len - 10, '+', window=50
        )

        assert result['insufficient_data']
        assert result['window_size'] < config.AG_RICHNESS_MIN_WINDOW
        assert not result['is_likely_misprimed']  # Don't flag insufficient data
        assert result['ag_content'] is None

    def test_invalid_chromosome(self):
        """Test with invalid chromosome."""
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrXXX', 1000, '+', window=50
        )

        assert result['insufficient_data']
        assert not result['is_likely_misprimed']

    def test_minimum_window_requirement(self):
        """Test minimum window size requirement."""
        # Test with position that gives exactly minimum window
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrI', 1000, '+',
            window=50, min_window=20
        )

        assert result['window_size'] >= 20
        assert not result['insufficient_data']


class TestBatchProcessing:
    """Test batch AG screening."""

    def test_batch_screening(self):
        """Test batch processing of positions."""
        positions = [
            {'chrom': 'chrI', 'position': 1000, 'strand': '+'},  # Low AG
            {'chrom': 'chrI', 'position': 1100, 'strand': '+'},  # High AG
            {'chrom': 'chrI', 'position': 1150, 'strand': '+'},  # Very high AG
        ]

        results = ag_mispriming.screen_ag_mispriming_batch(
            MOCK_GENOME, positions, window=50
        )

        assert len(results) == 3
        assert not results[0]['is_likely_misprimed']
        assert results[1]['is_likely_misprimed']
        assert results[2]['is_likely_misprimed']

    def test_empty_batch(self):
        """Test batch processing with empty input."""
        results = ag_mispriming.screen_ag_mispriming_batch(
            MOCK_GENOME, [], window=50
        )

        assert len(results) == 0


class TestStatistics:
    """Test statistical summaries."""

    def test_calculate_statistics(self):
        """Test statistics calculation."""
        ag_results = [
            {'ag_score': 1.0, 'ag_content': 0.3, 'is_likely_misprimed': False, 'insufficient_data': False},
            {'ag_score': 18.5, 'ag_content': 0.7, 'is_likely_misprimed': True, 'insufficient_data': False},
            {'ag_score': 23.0, 'ag_content': 0.8, 'is_likely_misprimed': True, 'insufficient_data': False},
        ]

        stats = ag_mispriming.calculate_ag_statistics(ag_results)

        assert stats['total'] == 3
        assert stats['valid'] == 3
        assert stats['likely_misprimed'] == 2
        assert stats['mispriming_rate'] == pytest.approx(2/3)
        assert stats['mean_ag_content'] == pytest.approx(0.6)

    def test_statistics_with_insufficient_data(self):
        """Test statistics with some insufficient data."""
        ag_results = [
            {'ag_score': 1.0, 'ag_content': 0.3, 'is_likely_misprimed': False, 'insufficient_data': False},
            {'ag_score': None, 'ag_content': None, 'is_likely_misprimed': False, 'insufficient_data': True},
        ]

        stats = ag_mispriming.calculate_ag_statistics(ag_results)

        assert stats['total'] == 2
        assert stats['valid'] == 1
        assert stats['insufficient_data'] == 1

    def test_statistics_empty(self):
        """Test statistics with empty input."""
        stats = ag_mispriming.calculate_ag_statistics([])

        assert stats['total'] == 0
        assert stats['mispriming_rate'] == 0.0

    def test_format_report(self):
        """Test report formatting."""
        stats = {
            'total': 100,
            'valid': 95,
            'likely_misprimed': 12,
            'mispriming_rate': 0.126,
            'mean_ag_content': 0.55,
            'median_ag_content': 0.52,
            'max_ag_content': 0.85,
            'insufficient_data': 5,
        }

        report = ag_mispriming.format_ag_report(stats)

        assert 'AG Mispriming' in report
        assert '100' in report
        assert '12.6%' in report
        assert '55.0%' in report  # Mean AG content


class TestQCFlags:
    """Test QC flag assignment.

    The flag describes the *genomic* sequence in the 19 bp window downstream
    (in mRNA orientation) of a read's called 3' end — not the read itself.
    AG_RICH means that window looks like an oligo-dT(V) priming substrate,
    so the called 3' end is plausibly an internal-priming artefact.

    Single-band classification: AG_RICH / PASS / INSUFFICIENT_DATA. The
    historical AG_RICH_HIGH / AG_RICH_MEDIUM distinction was removed when
    the two thresholds were collapsed onto the single Youden-optimal cutoff
    at score=17.0.
    """

    def test_qc_flag_pass(self):
        """PASS — downstream genomic window not AG-enriched."""
        result = {'ag_content': 0.4, 'is_likely_misprimed': False, 'insufficient_data': False}
        assert ag_mispriming.get_ag_qc_flag(result) == 'PASS'

    def test_qc_flag_ag_rich(self):
        """AG_RICH — downstream genomic window scores above threshold."""
        result = {'ag_content': 0.8, 'is_likely_misprimed': True, 'insufficient_data': False}
        assert ag_mispriming.get_ag_qc_flag(result) == 'AG_RICH'

    def test_qc_flag_insufficient_data(self):
        """INSUFFICIENT_DATA dominates when the downstream window is shorter than min_window."""
        result = {'ag_content': None, 'is_likely_misprimed': False, 'insufficient_data': True}
        assert ag_mispriming.get_ag_qc_flag(result) == 'INSUFFICIENT_DATA'


class TestScoreClassification:
    """Test ag_score-driven flagging across the MOCK_GENOME positions.

    Each position represents a putative called 3' end whose 19 bp downstream
    genomic context drives the AG-richness score.
    """

    def test_score_above_threshold_flags_high_ag_context(self):
        """High-AG downstream genomic context (pos 1100) → flagged."""
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrI', 1100, '+', window=50
        )
        assert result['ag_content'] >= 0.75
        assert result['ag_score'] > config.AG_RICHNESS_SCORE_THRESHOLD
        assert result['is_likely_misprimed']

    def test_score_above_threshold_flags_medium_ag_context(self):
        """Medium-AG downstream genomic context (pos 1050) → flagged (single-cutoff contract)."""
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrI', 1050, '+', window=50
        )
        assert 0.5 < result['ag_content'] < 0.7
        assert result['ag_score'] > config.AG_RICHNESS_SCORE_THRESHOLD
        assert result['is_likely_misprimed']

    def test_score_below_threshold_not_flagged(self):
        """Low-AG downstream genomic context (pos 1000) → not flagged."""
        result = ag_mispriming.screen_ag_mispriming(
            MOCK_GENOME, 'chrI', 1000, '+', window=50
        )
        assert result['ag_content'] < 0.65
        assert result['ag_score'] <= config.AG_RICHNESS_SCORE_THRESHOLD
        assert not result['is_likely_misprimed']
