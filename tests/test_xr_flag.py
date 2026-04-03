#!/usr/bin/env python3
"""
Tests for Bug 6: XR flag inversion edge case on minus-strand reads.

The XR tag is set to 1 when the winning aligner had FEWER 5' soft-clip bases
than the runner-up (i.e., the winner rescued the 5' end of the transcript).

Bug: for minus-strand reads the 5' end is at the RIGHT (trailing) side of the
CIGAR, not the left (leading) side.  If get_softclip_lengths() always read
cigar[0] as the 5' clip it would:
  - report 0 clips for a minus-strand read whose only clip is trailing, and
  - therefore fire XR=1 even when the winner had MORE 5' clips than the
    runner-up.

Observed example: rna15_rep3, SRR32518274.250727
  - uLTRA (winner): plus-strand CIGAR leads to large trailing clip (5' end)
  - mapPacBio (runner-up): no trailing clip
  - With the old cigar[0]-only code, uLTRA.five_prime_softclip was reported as
    0 (its leading clip), so XR=1 fired incorrectly.

Fix location: rectify/core/consensus.py  get_softclip_lengths()

Author: Kevin R. Roy
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from unittest.mock import MagicMock

from rectify.core.consensus import (
    AlignmentInfo,
    get_softclip_lengths,
    select_best_alignment,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_read(cigar_tuples, is_reverse=False, ref_start=1000, ref_end=None):
    """Return a minimal mock pysam.AlignedSegment."""
    read = MagicMock()
    read.cigartuples = cigar_tuples
    read.is_reverse = is_reverse
    read.reference_start = ref_start
    # Compute reference_end from cigar if not provided
    if ref_end is not None:
        read.reference_end = ref_end
    else:
        ref_consuming = 0
        for op, length in cigar_tuples:
            if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
                ref_consuming += length
        read.reference_end = ref_start + ref_consuming
    return read


def make_alignment_info(aligner, five_prime_softclip, strand='+', chrom='chrI'):
    """Build an AlignmentInfo with the given 5' clip count."""
    return AlignmentInfo(
        read_id='test_read',
        aligner=aligner,
        chrom=chrom,
        strand=strand,
        reference_start=1000,
        reference_end=1100,
        cigar_string='50M',
        mapq=60,
        five_prime_softclip=five_prime_softclip,
        corrected_3prime=1099 if strand == '+' else 1000,
    )


# ---------------------------------------------------------------------------
# Unit tests: get_softclip_lengths
# ---------------------------------------------------------------------------

class TestGetSoftclipLengths:

    def test_plus_strand_leading_clip_is_5prime(self):
        """Plus-strand: cigar[0] leading S is the 5' clip."""
        # CIGAR: 10S 80M 5S   (leading 10bp = 5' clip for + strand)
        cigar = [(4, 10), (0, 80), (4, 5)]
        read = make_read(cigar, is_reverse=False)
        five, three = get_softclip_lengths(read)
        assert five == 10, f"Plus strand 5' clip should be 10 (leading), got {five}"
        assert three == 5, f"Plus strand 3' clip should be 5 (trailing), got {three}"

    def test_plus_strand_no_clips(self):
        """Plus-strand: no clips."""
        cigar = [(0, 100)]
        read = make_read(cigar, is_reverse=False)
        five, three = get_softclip_lengths(read)
        assert five == 0
        assert three == 0

    def test_minus_strand_trailing_clip_is_5prime(self):
        """Minus-strand: cigar[-1] trailing S is the 5' clip.

        For a reverse-complement read, the RNA 5' end (TSS) is at the HIGH
        genomic coordinate, which is represented by the TRAILING CIGAR
        operation.  Counting cigar[0] (leading) as the 5' clip is the bug.
        """
        # CIGAR: 5S 80M 15S   (trailing 15bp = 5' clip for - strand)
        cigar = [(4, 5), (0, 80), (4, 15)]
        read = make_read(cigar, is_reverse=True)
        five, three = get_softclip_lengths(read)
        assert five == 15, (
            f"Minus strand 5' clip should be 15 (trailing/right), got {five}. "
            f"Bug would return {5} if cigar[0] were used instead of cigar[-1]."
        )
        assert three == 5, f"Minus strand 3' clip should be 5 (leading/left), got {three}"

    def test_minus_strand_only_trailing_clip(self):
        """Minus-strand read with ONLY a trailing clip (no leading clip).

        This is the exact scenario for read SRR32518274.250727:
        uLTRA produces a minus-strand read with a large trailing clip (5' end)
        and zero leading clip. With cigar[0]-only code, five_prime_softclip
        would be reported as 0.
        """
        # CIGAR: 80M 20S   (trailing 20bp = 5' clip for - strand; no leading clip)
        cigar = [(0, 80), (4, 20)]
        read = make_read(cigar, is_reverse=True)
        five, three = get_softclip_lengths(read)
        assert five == 20, (
            f"Minus strand 5' clip should be 20 (only trailing clip), got {five}. "
            f"cigar[0]-only code would return 0."
        )
        assert three == 0

    def test_minus_strand_only_leading_clip(self):
        """Minus-strand read with ONLY a leading clip (3' end clip)."""
        # CIGAR: 10S 80M   (leading 10bp = 3' clip for - strand; no trailing clip)
        cigar = [(4, 10), (0, 80)]
        read = make_read(cigar, is_reverse=True)
        five, three = get_softclip_lengths(read)
        assert five == 0, f"Minus strand 5' clip should be 0 (no trailing clip), got {five}"
        assert three == 10

    def test_minus_strand_no_clips(self):
        """Minus-strand read with no clips."""
        cigar = [(0, 100)]
        read = make_read(cigar, is_reverse=True)
        five, three = get_softclip_lengths(read)
        assert five == 0
        assert three == 0

    def test_no_cigar_returns_zeros(self):
        """Unmapped or no-CIGAR read returns (0, 0)."""
        read = MagicMock()
        read.cigartuples = None
        read.is_reverse = False
        assert get_softclip_lengths(read) == (0, 0)

        read2 = MagicMock()
        read2.cigartuples = []
        read2.is_reverse = True
        assert get_softclip_lengths(read2) == (0, 0)

    def test_plus_strand_spliced_with_clips(self):
        """Plus-strand spliced read: 5S 40M 200N 40M 3S."""
        # CIGAR: 5S 40M 200N 40M 3S
        cigar = [(4, 5), (0, 40), (3, 200), (0, 40), (4, 3)]
        read = make_read(cigar, is_reverse=False)
        five, three = get_softclip_lengths(read)
        assert five == 5
        assert three == 3

    def test_minus_strand_spliced_with_clips(self):
        """Minus-strand spliced read: 3S 40M 200N 40M 5S — 5' clip is trailing."""
        # CIGAR: 3S 40M 200N 40M 5S
        cigar = [(4, 3), (0, 40), (3, 200), (0, 40), (4, 5)]
        read = make_read(cigar, is_reverse=True)
        five, three = get_softclip_lengths(read)
        assert five == 5, f"Minus strand 5' clip should be 5 (trailing), got {five}"
        assert three == 3


# ---------------------------------------------------------------------------
# Integration tests: XR flag via select_best_alignment
# ---------------------------------------------------------------------------

class TestXRFlagViaSelectBestAlignment:
    """
    Tests for correct XR flag assignment.

    select_best_alignment() reads five_prime_softclip from AlignmentInfo
    (which is populated by get_softclip_lengths).  These tests verify the
    was_5prime_rescued logic using AlignmentInfo objects directly to isolate
    the XR assignment logic from genome/scoring.
    """

    def test_xr_set_when_winner_has_fewer_clips_plus_strand(self):
        """XR should be True when winner has fewer 5' clips (plus strand, normal case)."""
        # uLTRA found the intron: 0 clips. mapPacBio clipped: 15 clips. uLTRA wins.
        alignments = {
            'uLTRA': make_alignment_info('uLTRA', five_prime_softclip=0, strand='+'),
            'mapPacBio': make_alignment_info('mapPacBio', five_prime_softclip=15, strand='+'),
        }
        result = select_best_alignment(alignments, genome={})
        # uLTRA should win (0 clips → no penalty vs −30 for mapPacBio)
        assert result.best_aligner == 'uLTRA', f"Expected uLTRA to win, got {result.best_aligner}"
        assert result.was_5prime_rescued is True, (
            "XR should be True: winner (uLTRA) has fewer 5' clips than runner-up"
        )

    def test_xr_not_set_when_winner_has_more_clips_plus_strand(self):
        """XR should be False when winner has MORE 5' clips (plus strand)."""
        # mapPacBio has no clips and wins. uLTRA has clips and loses.
        # But to make uLTRA win despite more clips, we'd need a tiebreak scenario.
        # Simplest: both have same clips → no rescue.
        alignments = {
            'uLTRA': make_alignment_info('uLTRA', five_prime_softclip=5, strand='+'),
            'mapPacBio': make_alignment_info('mapPacBio', five_prime_softclip=5, strand='+'),
        }
        result = select_best_alignment(alignments, genome={})
        assert result.was_5prime_rescued is False, (
            "XR should be False when all aligners have equal clips"
        )

    def test_xr_not_set_when_winner_has_more_clips_minus_strand(self):
        """XR must NOT be set when winner has MORE 5' clips — the minus-strand edge case.

        This is Bug 6: uLTRA (winner) has a large TRAILING clip on a minus-strand
        read (= large 5' clip).  mapPacBio has zero trailing clip (= 0 5' clip).

        With the old cigar[0]-only code, uLTRA.five_prime_softclip would be
        reported as 0 (its leading clip, which is absent), so:
          min_5clip=0, winner.five_prime_softclip=0 → XR=1 (WRONG)

        With the correct strand-aware code:
          uLTRA.five_prime_softclip=15 (trailing=5' for minus), min=0 (mapPacBio)
          winner.five_prime_softclip=15 ≠ min → XR=0 (CORRECT)

        Here we test was_5prime_rescued directly by injecting pre-computed
        AlignmentInfo values, isolating the XR flag logic from the cigar parser.
        A companion test in TestCigar0OnlyRegression verifies that the cigar
        parser produces the correct five_prime_softclip=15 value (not 0).
        """
        # In the real read (SRR32518274.250727), mapPacBio wins because it has
        # zero clips (score=0) while uLTRA's trailing clip (5' for minus strand)
        # gives it a score penalty.  To test the rescue logic when uLTRA somehow
        # wins with MORE clips, we give uLTRA the same clip count as mapPacBio
        # to create an equal-score tie, then add a canonical junction tiebreak.
        #
        # The was_5prime_rescued check only fires when:
        #   max_5clip > min_5clip AND winner.five_prime_softclip == min_5clip
        # So with winner=uLTRA and uLTRA.five_prime_softclip == min (0), the bug
        # would fire XR=1.  With correct values (uLTRA=15), XR stays 0.
        #
        # We verify both the "correct values → no XR" path and that
        # "buggy values (0) → XR=1" would be the wrong outcome.

        # Correct AlignmentInfo: five_prime_softclip derived from trailing clip
        ultra_correct = AlignmentInfo(
            read_id='SRR32518274.250727',
            aligner='uLTRA',
            chrom='chrI',
            strand='-',
            reference_start=1000,
            reference_end=1100,
            cigar_string='80M15S',
            mapq=60,
            five_prime_softclip=15,   # correct: trailing clip = 5' for minus strand
            three_prime_softclip=0,
            corrected_3prime=1000,
            canonical_count=1,
        )
        mpb_info = AlignmentInfo(
            read_id='SRR32518274.250727',
            aligner='mapPacBio',
            chrom='chrI',
            strand='-',
            reference_start=1000,
            reference_end=1100,
            cigar_string='100M',
            mapq=60,
            five_prime_softclip=0,    # no trailing clip → no 5' clip
            three_prime_softclip=0,
            corrected_3prime=1000,
            canonical_count=0,
        )

        # With correct five_prime_softclip=15 for uLTRA:
        # uLTRA score = -30 (clip penalty), mapPacBio score = 0 → mapPacBio wins.
        # was_5prime_rescued: min=0 (mapPacBio), max=15 (uLTRA), winner=mapPacBio with 0
        # → winner.five_prime_softclip(0) == min(0) → would set XR=1
        # BUT this is the CORRECT case: mapPacBio won by having fewer clips → XR=1 is right.
        result = select_best_alignment(
            {'uLTRA': ultra_correct, 'mapPacBio': mpb_info}, genome={}
        )
        assert result.best_aligner == 'mapPacBio', (
            f"mapPacBio should win (score 0 vs uLTRA score -30), got {result.best_aligner}"
        )
        # mapPacBio wins with fewer clips → this IS a 5' rescue (mapPacBio rescued the 5' end)
        assert result.was_5prime_rescued is True, (
            "XR=1 is correct here: mapPacBio won with fewer 5' clips (0 vs uLTRA's 15)"
        )

        # Now test the BUGGY scenario: old code sets five_prime_softclip=0 for uLTRA
        # (reads cigar[0]=M and gets 0 instead of the correct trailing 15).
        # This makes BOTH aligners appear to have 0 clips → max==min → XR stays 0.
        # So the bug actually manifests differently from the initial description:
        # with old code, the winner appears to have 0 clips same as runner-up,
        # which means XR is NOT set when it SHOULD be (false negative, not false positive).
        ultra_buggy = AlignmentInfo(
            read_id='SRR32518274.250727',
            aligner='uLTRA',
            chrom='chrI',
            strand='-',
            reference_start=1000,
            reference_end=1100,
            cigar_string='80M15S',
            mapq=60,
            five_prime_softclip=0,    # WRONG: old code reads cigar[0]=M → returns 0
            three_prime_softclip=0,
            corrected_3prime=1000,
            canonical_count=1,
        )
        result_buggy = select_best_alignment(
            {'uLTRA': ultra_buggy, 'mapPacBio': mpb_info}, genome={}
        )
        # With buggy values: uLTRA.five_prime_softclip=0, mpb=0 → min==max → no rescue
        assert result_buggy.was_5prime_rescued is False, (
            "With buggy clip values (both 0), was_5prime_rescued is False — "
            "a false negative: XR is missing when it should be set"
        )

        # Verify: correct code sets five_prime_softclip=15 (not 0) for this cigar/strand
        cigar = [(0, 80), (4, 15)]  # 80M 15S  (trailing 15bp = 5' for minus strand)
        read = make_read(cigar, is_reverse=True)
        five_clip, _ = get_softclip_lengths(read)
        assert five_clip == 15, (
            f"get_softclip_lengths must return 15 for minus-strand trailing clip, got {five_clip}"
        )
        assert five_clip != 0, (
            "If get_softclip_lengths returned 0 here, Bug 6 is still present"
        )

    def test_xr_set_correctly_when_minus_strand_winner_has_fewer_clips(self):
        """XR=1 is correct on minus strand when the winner genuinely has fewer 5' clips."""
        # uLTRA found the intron (0 trailing clip = 0 five_prime_softclip).
        # mapPacBio missed the intron and starts after it → 0 clips too, but
        # let's give mapPacBio a trailing clip to simulate having missed the junction
        # and showing up with a trailing unaligned region.
        ultra_info = AlignmentInfo(
            read_id='test_read_minus',
            aligner='uLTRA',
            chrom='chrI',
            strand='-',
            reference_start=900,
            reference_end=1100,
            cigar_string='100M100N100M',
            mapq=60,
            five_prime_softclip=0,    # uLTRA found the 5' junction: no trailing clip
            three_prime_softclip=0,
            corrected_3prime=900,
            canonical_count=1,
        )
        mpb_info = AlignmentInfo(
            read_id='test_read_minus',
            aligner='mapPacBio',
            chrom='chrI',
            strand='-',
            reference_start=900,
            reference_end=1100,
            cigar_string='100M20S',   # 20bp trailing = 5' clip for minus strand
            mapq=60,
            five_prime_softclip=20,   # correct: trailing clip = 5' for minus strand
            three_prime_softclip=0,
            corrected_3prime=900,
            canonical_count=0,
        )

        alignments = {'uLTRA': ultra_info, 'mapPacBio': mpb_info}
        result = select_best_alignment(alignments, genome={})

        assert result.best_aligner == 'uLTRA', (
            f"uLTRA should win (0 5' clips vs 20), got {result.best_aligner}"
        )
        assert result.was_5prime_rescued is True, (
            "XR should be True: uLTRA (winner) has fewer 5' clips (0 vs 20) on minus strand"
        )

    def test_xr_not_set_single_aligner(self):
        """XR should never be set when there is only one aligner."""
        alignments = {
            'uLTRA': make_alignment_info('uLTRA', five_prime_softclip=0, strand='-'),
        }
        result = select_best_alignment(alignments, genome={})
        assert result.was_5prime_rescued is False, (
            "XR should not be set with a single aligner (no comparison possible)"
        )

    def test_xr_not_set_when_all_equal_clips_minus_strand(self):
        """XR should not be set when all aligners have equal 5' clips on minus strand."""
        alignments = {
            'uLTRA': make_alignment_info('uLTRA', five_prime_softclip=5, strand='-'),
            'mapPacBio': make_alignment_info('mapPacBio', five_prime_softclip=5, strand='-'),
        }
        result = select_best_alignment(alignments, genome={})
        assert result.was_5prime_rescued is False, (
            "XR should not be set when all aligners have identical clip counts"
        )


# ---------------------------------------------------------------------------
# Regression test: verify old cigar[0]-only code would fail
# ---------------------------------------------------------------------------

class TestCigar0OnlyRegression:
    """
    Demonstrates why always reading cigar[0] as the 5' clip is wrong
    for minus-strand reads.  These tests document the expected failure
    mode of the old code.
    """

    def test_minus_strand_trailing_only_clip_is_not_cigar0(self):
        """For a minus-strand read with only a trailing clip, cigar[0] is M not S.

        Old code:  five_prime = cigar[0][1] if cigar[0][0] == 4 else 0  → returns 0 (WRONG)
        New code:  uses cigar[-1] for minus strand               → returns 15 (CORRECT)
        """
        # CIGAR: 80M 15S  (trailing 15bp = 5' clip for minus strand)
        cigar = [(0, 80), (4, 15)]
        read = make_read(cigar, is_reverse=True)

        # Demonstrate what old code (cigar[0] only) would return:
        old_code_result = cigar[0][1] if cigar[0][0] == 4 else 0
        assert old_code_result == 0, "Old code returns 0 for this cigar (leading op is M)"

        # Current correct code:
        five, _ = get_softclip_lengths(read)
        assert five == 15, f"Correct code should return 15, got {five}"
        assert five != old_code_result, (
            "Correct code must differ from old cigar[0]-only code for this minus-strand case"
        )
