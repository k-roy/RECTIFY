#!/usr/bin/env python3
"""
Tests for indel artifact correction module.
"""

import pytest
from unittest.mock import Mock, patch
from rectify.core import indel_corrector
from rectify import config


class TestIsAtractDeletion:
    """Test A-tract deletion classification."""

    def test_small_deletion_with_a_rich_flanks(self):
        """Test small deletion (2bp) with A-rich flanks (is artifact)."""
        deletion = {
            'pos': 10,
            'length': 2,
            'ref_seq': 'TT',  # Deleted bases are NOT A's
            'read_pos': 10,
        }
        read_seq = 'AAAAAAAAAAAAAAAAAA'  # All A's

        result = indel_corrector.is_atract_deletion(
            deletion, '', read_seq, '+'
        )

        assert result  # Should be classified as artifact

    def test_deletion_of_as(self):
        """Test deletion of A's (not an artifact, expected in A-tract)."""
        deletion = {
            'pos': 10,
            'length': 2,
            'ref_seq': 'AA',  # Deleted bases ARE A's
            'read_pos': 10,
        }
        read_seq = 'AAAAAAAAAAAAAAAAAA'

        result = indel_corrector.is_atract_deletion(
            deletion, '', read_seq, '+'
        )

        assert not result  # Expected deletion, not artifact

    def test_large_deletion(self):
        """Test large deletion (>3bp) should not be artifact."""
        deletion = {
            'pos': 10,
            'length': 5,  # Too large
            'ref_seq': 'TTCCC',
            'read_pos': 10,
        }
        read_seq = 'AAAAAAAAAAAAAAAAAA'

        result = indel_corrector.is_atract_deletion(
            deletion, '', read_seq, '+'
        )

        assert not result  # Too large to be artifact

    def test_low_a_richness_flanks(self):
        """Test deletion with non-A-rich flanks (not an artifact)."""
        deletion = {
            'pos': 10,
            'length': 2,
            'ref_seq': 'TT',
            'read_pos': 10,
        }
        read_seq = 'TCGATCGATCGATCGATC'  # Mixed bases, not A-rich

        result = indel_corrector.is_atract_deletion(
            deletion, '', read_seq, '+'
        )

        assert not result  # No A-rich flanks

    def test_minus_strand_deletion(self):
        """Test deletion classification on minus strand."""
        deletion = {
            'pos': 10,
            'length': 2,
            'ref_seq': 'AA',  # On minus strand, looking for T's
            'read_pos': 10,
        }
        read_seq = 'TTTTTTTTTTTTTTTTTT'  # All T's (poly-A in RNA orientation)

        # On minus strand, deleted bases should NOT be T's for artifact
        # These are A's, so this IS an artifact
        result = indel_corrector.is_atract_deletion(
            deletion, '', read_seq, '-'
        )

        assert result


class TestDetectIndelArtifacts:
    """Test indel artifact detection."""

    def create_mock_read(self, ref_start, ref_end, cigar_tuples, sequence, strand='+'):
        """Create mock pysam read."""
        read = Mock()
        read.reference_start = ref_start
        read.reference_end = ref_end
        read.reference_name = 'chrI'
        read.cigartuples = cigar_tuples
        read.query_sequence = sequence
        read.query_length = len(sequence)
        read.is_reverse = (strand == '-')
        return read

    def test_detect_deletion_near_3prime_plus(self):
        """Test detection of deletion near 3' end (+ strand)."""
        # Mock read: 80M2D18M (deletion at position 80)
        # 3' end is at right, so deletion at 80 is 18bp from end
        read = self.create_mock_read(
            ref_start=1000,
            ref_end=1100,  # 80 + 2(D) + 18 = 100
            cigar_tuples=[(0, 80), (2, 2), (0, 18)],  # 80M, 2D, 18M
            sequence='A' * 98,  # 80 + 18 = 98 bases
            strand='+'
        )

        with patch('rectify.core.indel_corrector.extract_deletions') as mock_del:
            mock_del.return_value = [
                {'read_pos': 80, 'ref_pos': 1080, 'length': 2, 'ref_seq': 'TT'}
            ]

            artifacts = indel_corrector.detect_indel_artifacts(read, '+')

            assert len(artifacts) == 1
            assert artifacts[0]['type'] == 'deletion'
            assert artifacts[0]['length'] == 2

    def test_detect_deletion_near_3prime_minus(self):
        """Test detection of deletion near 3' end (- strand)."""
        # For minus strand, 3' end is at LEFT (start of read)
        read = self.create_mock_read(
            ref_start=1000,
            ref_end=1100,
            cigar_tuples=[(0, 18), (2, 2), (0, 80)],  # 18M, 2D, 80M
            sequence='A' * 98,
            strand='-'
        )

        with patch('rectify.core.indel_corrector.extract_deletions') as mock_del:
            mock_del.return_value = [
                {'read_pos': 18, 'ref_pos': 1018, 'length': 2, 'ref_seq': 'AA'}
            ]

            artifacts = indel_corrector.detect_indel_artifacts(read, '-')

            assert len(artifacts) == 1

    def test_no_deletions(self):
        """Test read with no deletions."""
        read = self.create_mock_read(
            ref_start=1000,
            ref_end=1100,
            cigar_tuples=[(0, 100)],  # 100M
            sequence='A' * 100,
            strand='+'
        )

        with patch('rectify.core.indel_corrector.extract_deletions') as mock_del:
            mock_del.return_value = []

            artifacts = indel_corrector.detect_indel_artifacts(read, '+')

            assert len(artifacts) == 0


class TestCorrectPositionForIndels:
    """Test position correction for indel artifacts."""

    def test_correct_for_deletion_plus_strand(self):
        """Test correction for deletion artifact (+ strand)."""
        # Deletion artifact: aligner removed 2bp that should be there
        # For + strand, move position LEFTWARD (upstream)
        artifacts = [
            {
                'type': 'deletion',
                'length': 2,
                'is_artifact': True,
            }
        ]

        result = indel_corrector.correct_position_for_indels(
            original_pos=1099,
            indel_artifacts=artifacts,
            strand='+'
        )

        assert result['corrected_position'] == 1097  # 1099 - 2
        assert result['correction_bp'] == -2
        assert result['n_deletions'] == 1

    def test_correct_for_deletion_minus_strand(self):
        """Test correction for deletion artifact (- strand)."""
        # For minus strand, move position RIGHTWARD (downstream in genomic coords)
        artifacts = [
            {
                'type': 'deletion',
                'length': 2,
                'is_artifact': True,
            }
        ]

        result = indel_corrector.correct_position_for_indels(
            original_pos=1000,
            indel_artifacts=artifacts,
            strand='-'
        )

        assert result['corrected_position'] == 1002  # 1000 + 2
        assert result['correction_bp'] == 2
        assert result['n_deletions'] == 1

    def test_correct_for_insertion(self):
        """Test correction for insertion artifact."""
        artifacts = [
            {
                'type': 'insertion',
                'length': 1,
                'is_artifact': True,
            }
        ]

        result = indel_corrector.correct_position_for_indels(
            original_pos=1099,
            indel_artifacts=artifacts,
            strand='+'
        )

        # Insertion: aligner added bases, move RIGHTWARD
        assert result['corrected_position'] == 1100  # 1099 + 1
        assert result['correction_bp'] == 1
        assert result['n_insertions'] == 1

    def test_multiple_artifacts(self):
        """Test correction for multiple artifacts."""
        artifacts = [
            {'type': 'deletion', 'length': 2, 'is_artifact': True},
            {'type': 'deletion', 'length': 1, 'is_artifact': True},
        ]

        result = indel_corrector.correct_position_for_indels(
            original_pos=1099,
            indel_artifacts=artifacts,
            strand='+'
        )

        assert result['corrected_position'] == 1096  # 1099 - 2 - 1
        assert result['correction_bp'] == -3
        assert result['n_deletions'] == 2

    def test_non_artifact_ignored(self):
        """Test that non-artifacts are ignored."""
        artifacts = [
            {'type': 'deletion', 'length': 2, 'is_artifact': False},  # Not artifact
            {'type': 'deletion', 'length': 1, 'is_artifact': True},   # Artifact
        ]

        result = indel_corrector.correct_position_for_indels(
            original_pos=1099,
            indel_artifacts=artifacts,
            strand='+'
        )

        # Only 1bp correction (from second deletion)
        assert result['corrected_position'] == 1098
        assert result['n_deletions'] == 1


class TestCorrectIndelsFromRead:
    """Test full indel correction workflow."""

    def create_mock_read(self, ref_start, ref_end, cigar_tuples, sequence, strand='+'):
        """Create mock pysam read."""
        read = Mock()
        read.reference_start = ref_start
        read.reference_end = ref_end
        read.reference_name = 'chrI'
        read.cigartuples = cigar_tuples
        read.query_sequence = sequence
        read.query_length = len(sequence)
        read.is_reverse = (strand == '-')
        return read

    def test_correct_read_with_artifact(self):
        """Test full correction for read with artifact."""
        read = self.create_mock_read(
            ref_start=1000,
            ref_end=1100,
            cigar_tuples=[(0, 80), (2, 2), (0, 18)],
            sequence='A' * 98,
            strand='+'
        )

        # Mock detect_indel_artifacts to return artifact
        with patch('rectify.core.indel_corrector.detect_indel_artifacts') as mock_detect:
            mock_detect.return_value = [
                {
                    'type': 'deletion',
                    'length': 2,
                    'is_artifact': True,
                    'pos': 80,
                    'ref_pos': 1080,
                }
            ]

            result = indel_corrector.correct_indels_from_read(read, '+')

            assert result['original_3prime'] == 1099
            assert result['has_artifacts']
            assert result['corrected_3prime'] == 1097  # Moved upstream by 2bp
            assert result['correction_bp'] == -2

    def test_correct_read_no_artifact(self):
        """Test correction for read with no artifacts."""
        read = self.create_mock_read(
            ref_start=1000,
            ref_end=1100,
            cigar_tuples=[(0, 100)],
            sequence='A' * 100,
            strand='+'
        )

        with patch('rectify.core.indel_corrector.detect_indel_artifacts') as mock_detect:
            mock_detect.return_value = []

            result = indel_corrector.correct_indels_from_read(read, '+')

            assert result['original_3prime'] == 1099
            assert not result['has_artifacts']
            assert result['corrected_3prime'] == 1099  # No change
            assert result['correction_bp'] == 0


class TestStatistics:
    """Test indel correction statistics."""

    def test_calculate_statistics(self):
        """Test statistics calculation."""
        results = [
            {'has_artifacts': True, 'correction_bp': -2},
            {'has_artifacts': True, 'correction_bp': -3},
            {'has_artifacts': False, 'correction_bp': 0},
        ]

        stats = indel_corrector.calculate_indel_statistics(results)

        assert stats['total'] == 3
        assert stats['with_artifacts'] == 2
        assert stats['artifact_rate'] == pytest.approx(2/3)
        assert stats['mean_correction'] == 2.5  # abs(-2) + abs(-3) / 2
        assert stats['median_correction'] == 2.5

    def test_statistics_empty(self):
        """Test statistics with empty input."""
        stats = indel_corrector.calculate_indel_statistics([])

        assert stats['total'] == 0
        assert stats['artifact_rate'] == 0.0

    def test_format_report(self):
        """Test report formatting."""
        stats = {
            'total': 100,
            'with_artifacts': 15,
            'artifact_rate': 0.15,
            'mean_correction': 2.3,
            'median_correction': 2.0,
            'max_correction': 3,
        }

        report = indel_corrector.format_indel_report(stats)

        assert 'Indel Artifact Correction' in report
        assert '100' in report
        assert '15.0%' in report
        assert '2.3 bp' in report


class TestEdgeCases:
    """Test edge cases."""

    def test_empty_deletion_info(self):
        """Test deletion with missing ref_seq."""
        deletion = {
            'pos': 10,
            'length': 2,
            'ref_seq': '',  # Empty
            'read_pos': 10,
        }
        read_seq = 'AAAAAAAAAAAAAAAAAA'

        result = indel_corrector.is_atract_deletion(
            deletion, '', read_seq, '+'
        )

        assert not result  # Can't classify without ref_seq

    def test_short_read_flanks(self):
        """Test deletion in very short read."""
        deletion = {
            'pos': 2,
            'length': 1,
            'ref_seq': 'T',
            'read_pos': 2,
        }
        read_seq = 'AAAA'  # Too short for proper flanks

        result = indel_corrector.is_atract_deletion(
            deletion, '', read_seq, '+'
        )

        assert not result  # Read too short
