#!/usr/bin/env python3
"""
Tests for A-tract ambiguity detection module.
"""

import pytest
from rectify.core import atract_detector
from rectify import config


# Mock genome for testing
# Using 'chrI' as key (atract_detector handles name conversion internally)
MOCK_GENOME = {
    'chrI': (
        'N' * 1000 +  # Padding
        'ATCG' * 10 +  # Position 1000-1039: No A-tract (10 repeats of ATCG = 40bp)
        'AAAAAAAAAAAA' +  # Position 1040-1051: 12 A's
        'TTTTTTTTTTTT' +  # Position 1052-1063: 12 T's (for - strand)
        'ATGATGATGATG' +  # Position 1064-1075: Mixed
        'AAAAAAA' + 'TCGAT' +  # Position 1076-1087: 7 A's then mixed
        'N' * 10000  # More padding
    )
}


class TestCalculateAtractAmbiguity:
    """Test calculate_atract_ambiguity function."""

    def test_low_acount_plus_strand(self):
        """Test position with few downstream A's (+ strand)."""
        # Position 1010: next 10bp in 'ATCG'*10 section
        # Position 1010 = 'C', next 10bp = 'CGATCGATCG' has 2 A's
        result = atract_detector.calculate_atract_ambiguity(
            MOCK_GENOME, 'chrI', 1010, '+', downstream_bp=10
        )

        assert result['downstream_a_count'] == 2
        assert result['expected_shift'] == 0.3  # From config for 2A
        # ambiguity_min = int(1010 - 0.3) = 1009
        assert result['ambiguity_min'] == 1009
        assert result['ambiguity_max'] == 1010
        assert result['ambiguity_range'] == 1
        assert result['has_ambiguity']

    def test_high_acount_plus_strand(self):
        """Test position with 10 downstream A's (+ strand)."""
        # Position 1040 has 10 A's in next 10bp
        result = atract_detector.calculate_atract_ambiguity(
            MOCK_GENOME, 'chrI', 1040, '+', downstream_bp=10
        )

        assert result['downstream_a_count'] == 10
        assert result['expected_shift'] == 3.8  # From config
        assert result['ambiguity_range'] == 4  # 1040 - 1036
        assert result['ambiguity_min'] == 1036  # int(1040 - 3.8) = int(1036.2) = 1036
        assert result['ambiguity_max'] == 1040
        assert result['has_ambiguity']

    def test_seven_acount_plus_strand(self):
        """Test position with 7 downstream A's (+ strand)."""
        # Position 1076 has 7 A's in next 10bp
        result = atract_detector.calculate_atract_ambiguity(
            MOCK_GENOME, 'chrI', 1076, '+', downstream_bp=10
        )

        assert result['downstream_a_count'] == 7
        assert result['expected_shift'] == 2.6  # From config
        assert result['ambiguity_range'] == 3  # 1076 - 1073
        assert result['ambiguity_min'] == 1073  # int(1076 - 2.6)
        assert result['ambiguity_max'] == 1076

    def test_minus_strand_ttract(self):
        """Test position with T-tract on minus strand (RNA A's)."""
        # Position 1063 on minus strand: look LEFT for T's (downstream in gene coords)
        # Positions 1053-1062 contain T's
        result = atract_detector.calculate_atract_ambiguity(
            MOCK_GENOME, 'chrI', 1063, '-', downstream_bp=10
        )

        # Downstream in gene orientation (LEFT in genomic coords) should be reverse-complemented
        # So T's become A's in RNA orientation
        assert result['downstream_a_count'] == 10
        assert result['expected_shift'] == 3.8

        # Minus strand: ambiguity extends RIGHTWARD (higher genomic coords)
        assert result['ambiguity_min'] == 1063
        assert result['ambiguity_max'] == 1063 + 4  # int(1063 + 3.8) + 1
        assert result['ambiguity_range'] > 0

    def test_tract_length_calculation(self):
        """Test contiguous A-tract length calculation."""
        # Position 1076: 'AAAAAAA' starts here (positions 1076-1082 = 7 A's), then 'TCGAT'
        result = atract_detector.calculate_atract_ambiguity(
            MOCK_GENOME, 'chrI', 1076, '+', downstream_bp=10
        )

        assert result['tract_length'] == 7  # Contiguous A's starting at position 1076

    def test_saturated_acount(self):
        """Test A-counts > 10 use saturated shift value."""
        # Position 1040: 12 A's (positions 1040-1051), then T's at 1052+
        # With downstream_bp=12, we count all 12 A's
        result = atract_detector.calculate_atract_ambiguity(
            MOCK_GENOME, 'chrI', 1040, '+', downstream_bp=12
        )

        assert result['downstream_a_count'] == 12
        # Should use saturated value (>10 uses DEFAULT_MAX_SHIFT)
        assert result['expected_shift'] == config.DEFAULT_MAX_SHIFT


class TestAmbiguityCategories:
    """Test ambiguity categorization functions."""

    def test_get_ambiguity_category(self):
        """Test ambiguity range categorization."""
        assert atract_detector.get_ambiguity_category(0) == 'none'
        assert atract_detector.get_ambiguity_category(1) == 'low'
        assert atract_detector.get_ambiguity_category(2) == 'low'
        assert atract_detector.get_ambiguity_category(3) == 'medium'
        assert atract_detector.get_ambiguity_category(4) == 'medium'
        assert atract_detector.get_ambiguity_category(5) == 'high'
        assert atract_detector.get_ambiguity_category(10) == 'high'


class TestBatchProcessing:
    """Test batch ambiguity calculation."""

    def test_calculate_batch(self):
        """Test batch processing of positions."""
        positions = [
            {'chrom': 'chrI', 'position': 1010, 'strand': '+'},  # 2A (in ATCG section)
            {'chrom': 'chrI', 'position': 1040, 'strand': '+'},  # 10A (start of A-tract)
            {'chrom': 'chrI', 'position': 1076, 'strand': '+'},  # 7A (start of 7-A tract)
        ]

        results = atract_detector.calculate_atract_ambiguity_batch(
            MOCK_GENOME, positions, downstream_bp=10
        )

        assert len(results) == 3
        assert results[0]['downstream_a_count'] == 2  # Position 1010 has 2 A's in next 10bp
        assert results[1]['downstream_a_count'] == 10  # Position 1040 has 10 A's in next 10bp
        assert results[2]['downstream_a_count'] == 7  # Position 1076 has 7 A's in next 10bp

    def test_empty_batch(self):
        """Test batch processing with empty input."""
        results = atract_detector.calculate_atract_ambiguity_batch(
            MOCK_GENOME, [], downstream_bp=10
        )

        assert len(results) == 0


class TestSummaryStatistics:
    """Test summary statistics functions."""

    def test_summarize_distribution(self):
        """Test ambiguity distribution summary."""
        ambiguities = [
            {'ambiguity_range': 0, 'downstream_a_count': 0},
            {'ambiguity_range': 3, 'downstream_a_count': 7},
            {'ambiguity_range': 4, 'downstream_a_count': 10},
        ]

        summary = atract_detector.summarize_ambiguity_distribution(ambiguities)

        assert summary['total'] == 3
        assert summary['with_ambiguity'] == 2
        assert summary['mean_range'] == pytest.approx(7/3)
        assert summary['max_range'] == 4
        assert summary['by_category']['none'] == 1
        assert summary['by_category']['medium'] == 2

    def test_summarize_empty(self):
        """Test summary with empty input."""
        summary = atract_detector.summarize_ambiguity_distribution([])

        assert summary['total'] == 0
        assert summary['mean_range'] == 0.0

    def test_format_report(self):
        """Test report formatting."""
        summary = {
            'total': 100,
            'with_ambiguity': 85,
            'mean_range': 2.5,
            'median_range': 2.0,
            'max_range': 5,
            'by_category': {'none': 15, 'low': 30, 'medium': 40, 'high': 15},
            'by_acount': {0: 15, 4: 20, 7: 30, 10: 35},
        }

        report = atract_detector.format_ambiguity_report(summary)

        assert 'A-tract Ambiguity Summary' in report
        assert '100' in report
        assert '85.0%' in report
        assert '2.50 bp' in report  # Format uses 2 decimal places


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_near_chromosome_end(self):
        """Test position near chromosome end (insufficient downstream sequence)."""
        # Position near end of mock genome where we can't get full 10bp window
        # The end of the genome is N's, so downstream_a_count will be 0
        result = atract_detector.calculate_atract_ambiguity(
            MOCK_GENOME, 'chrI', len(MOCK_GENOME['chrI']) - 5, '+', downstream_bp=10
        )

        # Should handle gracefully - partial window of N's gives 0 A's
        assert result['downstream_a_count'] == 0
        assert result['expected_shift'] == 0.0
        assert result['ambiguity_range'] == 0

    def test_invalid_chromosome(self):
        """Test with invalid chromosome name."""
        result = atract_detector.calculate_atract_ambiguity(
            MOCK_GENOME, 'chrXXX', 1000, '+', downstream_bp=10
        )

        # Should return no ambiguity for invalid chrom
        assert result['downstream_a_count'] is None
        assert result['ambiguity_range'] == 0


class TestStrandLogic:
    """Test strand-specific ambiguity calculation."""

    def test_plus_strand_shifts_left(self):
        """Test + strand ambiguity extends leftward (upstream)."""
        # + strand with A-tract: observed position is shifted RIGHT
        # True CPA is LEFTWARD
        result = atract_detector.calculate_atract_ambiguity(
            MOCK_GENOME, 'chrI', 1076, '+', downstream_bp=10
        )

        assert result['ambiguity_min'] < 1076
        assert result['ambiguity_max'] == 1076

    def test_minus_strand_shifts_right(self):
        """Test - strand ambiguity extends rightward (downstream in genomic coords)."""
        # - strand with T-tract: observed position is shifted LEFT
        # True CPA is RIGHTWARD
        result = atract_detector.calculate_atract_ambiguity(
            MOCK_GENOME, 'chrI', 1063, '-', downstream_bp=10
        )

        assert result['ambiguity_min'] == 1063
        assert result['ambiguity_max'] > 1063
