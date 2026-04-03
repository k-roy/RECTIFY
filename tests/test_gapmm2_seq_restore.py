#!/usr/bin/env python3
"""
Tests for Bug 3 fix: gapmm2 PAF→BAM conversion drops read sequences (SEQ=*).

gapmm2 outputs PAF format which does not carry read sequences.  _paf_to_bam()
therefore leaves query_sequence=None on every gapmm2 BAM record.  When gapmm2
wins consensus selection the output BAM would contain SEQ=* records that break
all downstream steps (indel correction, poly-A trimming, etc.).

The fix in _process_and_write_batch() calls _restore_sequence_from_aligner_reads()
before writing, copying SEQ+QUAL from another aligner's record for the same read.

Author: Kevin R. Roy
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from unittest.mock import MagicMock, patch

from rectify.core.consensus import _restore_sequence_from_aligner_reads


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_mock_read(query_name="read1", query_sequence=None, query_qualities=None):
    """Return a minimal mock pysam.AlignedSegment."""
    read = MagicMock()
    read.query_name = query_name
    read.query_sequence = query_sequence
    read.query_qualities = query_qualities
    return read


# ---------------------------------------------------------------------------
# Tests for _restore_sequence_from_aligner_reads
# ---------------------------------------------------------------------------

class TestRestoreSequenceFromAlignerReads:
    """Tests for the sequence-restore helper."""

    def test_gapmm2_seq_none_restored_from_minimap2(self):
        """When gapmm2 wins with SEQ=None, sequence is copied from minimap2."""
        seq = "ACGTACGTACGT"
        quals = [30] * len(seq)

        gapmm2_read = make_mock_read("read1", query_sequence=None, query_qualities=None)
        minimap2_read = make_mock_read("read1", query_sequence=seq, query_qualities=quals)

        aligner_reads = {
            "gapmm2": gapmm2_read,
            "minimap2": minimap2_read,
        }

        _restore_sequence_from_aligner_reads(gapmm2_read, aligner_reads)

        assert gapmm2_read.query_sequence == seq
        assert gapmm2_read.query_qualities == quals

    def test_gapmm2_seq_none_restored_from_mapPacBio(self):
        """Sequence is restored from mapPacBio when minimap2 is absent."""
        seq = "TTTTAAAACCCC"
        quals = [25] * len(seq)

        gapmm2_read = make_mock_read("read2", query_sequence=None, query_qualities=None)
        mpb_read = make_mock_read("read2", query_sequence=seq, query_qualities=quals)

        aligner_reads = {
            "gapmm2": gapmm2_read,
            "mapPacBio": mpb_read,
        }

        _restore_sequence_from_aligner_reads(gapmm2_read, aligner_reads)

        assert gapmm2_read.query_sequence == seq
        assert gapmm2_read.query_qualities == quals

    def test_no_aligner_has_sequence_does_not_crash(self):
        """When all aligners have SEQ=None, function logs a warning and returns."""
        gapmm2_read = make_mock_read("read3", query_sequence=None)
        other_read = make_mock_read("read3", query_sequence=None)

        aligner_reads = {
            "gapmm2": gapmm2_read,
            "minimap2": other_read,
        }

        import logging
        with patch.object(logging.getLogger("rectify.core.consensus"), "warning") as mock_warn:
            _restore_sequence_from_aligner_reads(gapmm2_read, aligner_reads)
            mock_warn.assert_called_once()
            assert "read3" in mock_warn.call_args[0][0]

        # best_read is left unchanged (still None) — caller writes SEQ=* as fallback
        assert gapmm2_read.query_sequence is None

    def test_read_with_existing_sequence_is_not_overwritten(self):
        """A read that already has SEQ set is not modified (non-gapmm2 winner)."""
        original_seq = "GGGGCCCC"
        original_quals = [40] * len(original_seq)
        donor_seq = "AAAAAAAAAAAA"

        winner_read = make_mock_read("read4", query_sequence=original_seq,
                                     query_qualities=original_quals)
        donor_read = make_mock_read("read4", query_sequence=donor_seq)

        aligner_reads = {
            "minimap2": winner_read,
            "gapmm2": donor_read,
        }

        # _restore_sequence_from_aligner_reads should only be called when SEQ is None;
        # verify it correctly skips the winner itself and finds the first non-None donor
        # (in this test the winner already has SEQ, so this call is a no-op guard check)
        # The function returns on the first non-None donor — which may be the winner itself
        # if it is iterated first.  This test checks the return value is deterministic.
        _restore_sequence_from_aligner_reads(winner_read, aligner_reads)

        # Winner's sequence stays as-is (it was already non-None; first non-None found
        # is winner_read itself, so sequence is set to its own value — idempotent)
        assert winner_read.query_sequence == original_seq

    def test_single_aligner_no_sequence_does_not_crash(self):
        """Edge case: only one aligner present and it has no sequence."""
        gapmm2_read = make_mock_read("read5", query_sequence=None)
        aligner_reads = {"gapmm2": gapmm2_read}

        import logging
        with patch.object(logging.getLogger("rectify.core.consensus"), "warning") as mock_warn:
            _restore_sequence_from_aligner_reads(gapmm2_read, aligner_reads)
            mock_warn.assert_called_once()

        assert gapmm2_read.query_sequence is None

    def test_quality_scores_transferred_alongside_sequence(self):
        """Both sequence and quality scores are copied atomically."""
        seq = "ACGT"
        quals = [10, 20, 30, 40]

        gapmm2_read = make_mock_read("read6", query_sequence=None, query_qualities=None)
        minimap2_read = make_mock_read("read6", query_sequence=seq, query_qualities=quals)

        _restore_sequence_from_aligner_reads(
            gapmm2_read,
            {"gapmm2": gapmm2_read, "minimap2": minimap2_read}
        )

        assert gapmm2_read.query_sequence == seq
        assert gapmm2_read.query_qualities == quals
