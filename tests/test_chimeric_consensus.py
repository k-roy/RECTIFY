#!/usr/bin/env python3
"""
Tests for chimeric consensus CIGAR construction.

Covers the reference-continuity bug (Bug 7) where segments from different
aligners produced phantom large insertions (e.g., 947I) or implausible N
operations (e.g., 49426N) when their reference coordinates did not meet
cleanly at segment boundaries.

Author: Kevin R. Roy
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from unittest.mock import MagicMock

from rectify.core.chimeric_consensus import (
    CigarEvent,
    ChimericSegment,
    build_chimeric_cigar,
    _validate_chimeric_cigar,
    _merge_cigar_ops,
    cigar_to_events,
    find_sync_points,
    build_query_ref_map,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_read(cigar_tuples, ref_start, query_seq=None, is_reverse=False):
    """Return a minimal mock pysam.AlignedSegment."""
    read = MagicMock()
    read.cigartuples = cigar_tuples
    read.reference_start = ref_start
    read.query_length = sum(
        length for op, length in cigar_tuples
        if op in (0, 1, 4, 7, 8)  # query-consuming ops
    )
    read.query_sequence = query_seq or ('A' * read.query_length)
    read.query_qualities = None
    read.is_reverse = is_reverse
    return read


def seg(q_start, q_end, winner):
    s = ChimericSegment(q_start=q_start, q_end=q_end, position='interior')
    s.winning_aligner = winner
    return s


# ---------------------------------------------------------------------------
# _validate_chimeric_cigar
# ---------------------------------------------------------------------------

class TestValidateChimericCigar:
    def test_normal_cigar_passes(self):
        # 50M 200N 50M — normal splice junction
        cigar = [(0, 50), (3, 200), (0, 50)]
        assert _validate_chimeric_cigar(cigar, read_length=100)

    def test_large_insertion_rejected(self):
        # 947I in a 1869-bp read: max_insertion = 1869//4 = 467
        cigar = [(0, 23), (3, 24716), (1, 947), (0, 213)]
        assert not _validate_chimeric_cigar(cigar, read_length=1869)

    def test_insertion_at_threshold_accepted(self):
        # Exactly at limit: read=400, max=100; insertion=100 → OK
        cigar = [(0, 50), (1, 100), (0, 50)]
        assert _validate_chimeric_cigar(cigar, read_length=400)

    def test_insertion_one_over_threshold_rejected(self):
        cigar = [(0, 50), (1, 101), (0, 50)]
        assert not _validate_chimeric_cigar(cigar, read_length=400)

    def test_huge_intron_rejected(self):
        # 49426N exceeds default max_intron=10000
        cigar = [(0, 100), (3, 49426), (0, 100)]
        assert not _validate_chimeric_cigar(cigar, read_length=200)

    def test_intron_at_threshold_accepted(self):
        cigar = [(0, 100), (3, 10000), (0, 100)]
        assert _validate_chimeric_cigar(cigar, read_length=200)

    def test_custom_max_intron(self):
        # Yeast: max intron ~2000 bp
        cigar = [(0, 100), (3, 2001), (0, 100)]
        assert not _validate_chimeric_cigar(cigar, read_length=200, max_intron=2000)
        assert _validate_chimeric_cigar(cigar, read_length=200, max_intron=5000)


# ---------------------------------------------------------------------------
# build_chimeric_cigar — reference continuity
# ---------------------------------------------------------------------------

class TestBuildChimericCigar:
    """
    build_chimeric_cigar receives FULL-READ events for each aligner (not per-segment
    events). extract_events_for_query_range clips them to the segment's query range.
    Tests must construct events from a full-read CIGAR covering all query positions.
    """

    def test_continuous_single_segment_passthrough(self):
        """Single segment, no stitching — CIGAR should pass through unchanged."""
        read_a = make_read([(0, 100)], ref_start=1000)
        events_a = cigar_to_events([(0, 100)], 1000)

        segments = [seg(0, 100, 'a')]
        aligner_events = {'a': events_a}
        aligner_reads = {'a': read_a}

        ref_start, cigar = build_chimeric_cigar(segments, aligner_events, aligner_reads)
        assert ref_start == 1000
        assert cigar == [(0, 100)]

    def test_two_continuous_segments_no_gap(self):
        """
        Two segments from different aligners meeting at a sync point.

        aligner_a: 100M at ref 1000  → q[0..100] → r[1000..1100]
        aligner_b: 100M at ref 1000  → q[0..100] → r[1000..1100]

        Segment A (q=[0,50]) winner=a  → 50M  r[1000..1050]
        Segment B (q=[50,100]) winner=b → 50M  r[1050..1100]

        No gap; should merge to 100M.
        """
        events_a = cigar_to_events([(0, 100)], 1000)
        events_b = cigar_to_events([(0, 100)], 1000)

        segments = [seg(0, 50, 'a'), seg(50, 100, 'b')]
        aligner_events = {'a': events_a, 'b': events_b}
        aligner_reads = {
            'a': make_read([(0, 100)], 1000),
            'b': make_read([(0, 100)], 1000),
        }

        ref_start, cigar = build_chimeric_cigar(segments, aligner_events, aligner_reads)
        assert ref_start == 1000
        assert cigar == [(0, 100)]

    def test_reference_gap_bridged_with_N(self):
        """
        Segment A (winner=a) ends at ref 1050.
        Segment B (winner=b) starts at ref 1350 — 300-bp reference gap.

        aligner_a: 100M at ref 1000  → q[0..100]  r[1000..1100]
        aligner_b: 100M at ref 1300  → q[0..100]  r[1300..1400]

        When extracting segment B (q=[50,100]) from aligner_b:
          → first event starts at r=1350 (q=50 → r=1300+50=1350)
          → cur_ref after A = 1050
          → gap = 1350 - 1050 = 300  → N bridge inserted

        Result: 50M  300N  50M
        """
        events_a = cigar_to_events([(0, 100)], 1000)   # q[0..100] → r[1000..1100]
        events_b = cigar_to_events([(0, 100)], 1300)   # q[0..100] → r[1300..1400]

        segments = [seg(0, 50, 'a'), seg(50, 100, 'b')]
        aligner_events = {'a': events_a, 'b': events_b}
        aligner_reads = {
            'a': make_read([(0, 100)], 1000),
            'b': make_read([(0, 100)], 1300),
        }

        ref_start, cigar = build_chimeric_cigar(segments, aligner_events, aligner_reads)
        assert ref_start == 1000
        assert cigar == [(0, 50), (3, 300), (0, 50)]

    def test_reference_regression_returns_sentinel(self):
        """
        Segment B's aligner starts at ref 800, before where A ended (ref 1050).
        Reference regression → return (None, []) sentinel.

        aligner_a: 100M at ref 1000  → ends at r=1100; segment A [0..50] ends at r=1050
        aligner_b: 100M at ref 750   → segment B [50..100] starts at r=800 < 1050
        """
        events_a = cigar_to_events([(0, 100)], 1000)
        events_b = cigar_to_events([(0, 100)], 750)

        segments = [seg(0, 50, 'a'), seg(50, 100, 'b')]
        aligner_events = {'a': events_a, 'b': events_b}
        aligner_reads = {
            'a': make_read([(0, 100)], 1000),
            'b': make_read([(0, 100)], 750),
        }

        ref_start, cigar = build_chimeric_cigar(segments, aligner_events, aligner_reads)
        assert ref_start is None
        assert cigar == []

    def test_phantom_insertion_prevented(self):
        """
        Reproduces the class of bug seen in ysh1_rep2.

        aligner_a: 20I + 80M at ref 1000 → q[0..100], ends at r=1080
          Segment A (q=[0..50]) → 20I + 30M, ends at r=1030

        aligner_b: 100M at ref 1200 → q[0..100]
          Segment B (q=[50..100]) → 50M starting at r=1250

        Without fix: phantom ~220I from reference gap (1250-1030=220).
        With fix: 30M bridge inserted as 220N, result is 20I 30M 220N 50M.
        """
        # aligner_a: 20I then 80M at ref 1000
        events_a = cigar_to_events([(1, 20), (0, 80)], 1000)
        # aligner_b: 100M at ref 1200
        events_b = cigar_to_events([(0, 100)], 1200)

        segments = [seg(0, 50, 'a'), seg(50, 100, 'b')]
        aligner_events = {'a': events_a, 'b': events_b}
        aligner_reads = {
            'a': make_read([(1, 20), (0, 80)], 1000),
            'b': make_read([(0, 100)], 1200),
        }

        ref_start, cigar = build_chimeric_cigar(segments, aligner_events, aligner_reads)

        # Must not contain any I larger than the segment's own 20I
        for op, length in cigar:
            if op == 1:  # I
                assert length <= 20, (
                    f"Phantom insertion detected: {length}I in chimeric CIGAR {cigar}"
                )

        # Gap between r=1030 (end of A) and r=1250 (start of B's events) = 220N
        assert (3, 220) in cigar, f"Expected 220N bridge in {cigar}"

    def test_real_world_ysh1_like_cigar_rejected_by_validation(self):
        """
        A CIGAR resembling the ysh1_rep2 pathological case:
        100S 23M 24716N 20I 10M 1I 22M 164I 34107N 15M 194I ...

        This should fail _validate_chimeric_cigar (164I and 194I exceed
        read_len // 4 for a 1869-bp read, and 34107N exceeds max_intron=10000).
        """
        pathological = [
            (4, 100),   # 100S
            (0, 23),    # 23M
            (3, 24716), # 24716N — exceeds max_intron
            (1, 20),    # 20I
            (0, 10),    # 10M
            (1, 164),   # 164I — exceeds max_insertion
            (0, 15),    # 15M
        ]
        read_length = 100 + 23 + 20 + 10 + 164 + 15  # = 332

        assert not _validate_chimeric_cigar(
            pathological, read_length=read_length, max_intron=10_000
        )


# ---------------------------------------------------------------------------
# Integration: select_best_chimeric falls back on bad assembly
# ---------------------------------------------------------------------------

class TestSelectBestChimericFallback:
    """
    Test that select_best_chimeric falls back correctly when build_chimeric_cigar
    returns invalid results. _fallback_simple_selection uses a relative import
    (rectify.core.chimeric_consensus) not available in the test environment, so
    we mock it out.
    """

    def _make_chimeric_result(self, aligner_name, read):
        from rectify.core.chimeric_consensus import ChimericResult
        return ChimericResult(
            read_id=read.query_name,
            is_chimeric=False,
            segment_winners=[('whole', aligner_name, 0, read.query_length)],
            chimeric_cigar=list(read.cigartuples),
            chimeric_ref_start=read.reference_start,
            confidence='medium',
            n_segments=1,
            n_aligners_used=1,
            five_prime_aligner=aligner_name,
            interior_aligners=[aligner_name],
            three_prime_aligner=aligner_name,
        )

    def test_fallback_on_reference_regression(self):
        """
        When build_chimeric_cigar returns (None, []), select_best_chimeric
        must fall back to single-aligner rather than crashing.
        """
        from rectify.core.chimeric_consensus import select_best_chimeric
        import unittest.mock as mock

        genome = {'chrI': 'A' * 250_000}
        read_a = make_read([(0, 50)], 1000)
        read_a.reference_name = 'chrI'
        read_a.is_reverse = False
        read_a.query_name = 'test_read'
        read_a.query_sequence = 'A' * 50

        read_b = make_read([(0, 50)], 1050)
        read_b.reference_name = 'chrI'
        read_b.is_reverse = False
        read_b.query_name = 'test_read'
        read_b.query_sequence = 'A' * 50

        aligner_reads = {'a': read_a, 'b': read_b}
        fallback_result = self._make_chimeric_result('a', read_a)

        with mock.patch(
            'rectify.core.chimeric_consensus.build_chimeric_cigar',
            return_value=(None, []),
        ), mock.patch(
            'rectify.core.chimeric_consensus._fallback_simple_selection',
            return_value=fallback_result,
        ) as mock_fallback:
            result = select_best_chimeric(aligner_reads, genome)

        mock_fallback.assert_called_once()
        assert result is fallback_result

    def test_fallback_on_invalid_cigar(self):
        """
        When build_chimeric_cigar returns a CIGAR that fails validation
        (947I in a 50-bp read context), select_best_chimeric must call
        _fallback_simple_selection rather than returning the bad CIGAR.
        """
        from rectify.core.chimeric_consensus import select_best_chimeric
        import unittest.mock as mock

        genome = {'chrI': 'A' * 250_000}
        read_a = make_read([(0, 50)], 1000)
        read_a.reference_name = 'chrI'
        read_a.is_reverse = False
        read_a.query_name = 'test_read'
        read_a.query_sequence = 'A' * 50

        read_b = make_read([(0, 50)], 1050)
        read_b.reference_name = 'chrI'
        read_b.is_reverse = False
        read_b.query_name = 'test_read'
        read_b.query_sequence = 'A' * 50

        aligner_reads = {'a': read_a, 'b': read_b}
        fallback_result = self._make_chimeric_result('a', read_a)

        # 947I far exceeds max_insertion for a 50-bp read (max = max(50//4=12, 100)=100)
        bad_cigar = [(0, 23), (3, 24716), (1, 947), (0, 213)]

        with mock.patch(
            'rectify.core.chimeric_consensus.build_chimeric_cigar',
            return_value=(1000, bad_cigar),
        ), mock.patch(
            'rectify.core.chimeric_consensus._fallback_simple_selection',
            return_value=fallback_result,
        ) as mock_fallback:
            result = select_best_chimeric(aligner_reads, genome)

        mock_fallback.assert_called_once()
        assert result.chimeric_cigar != bad_cigar
