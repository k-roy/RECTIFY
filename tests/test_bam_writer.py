"""
Unit tests for bam_writer CIGAR surgery functions.

Focuses on correctness of clip_read_to_corrected_3prime and
softclip_read_to_corrected_3prime, especially edge cases involving
deletions (D) or intron skips (N) at the alignment terminus.

Run with:
    pytest tests/test_bam_writer.py -v

Author: Kevin R. Roy
"""

import pytest
import pysam
from typing import List, Tuple

from rectify.core.bam_writer import (
    clip_read_to_corrected_3prime,
    softclip_read_to_corrected_3prime,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_read(
    chrom: str,
    start: int,
    cigar: List[Tuple[int, int]],
    strand: str,
    seq: str,
) -> pysam.AlignedSegment:
    hdr = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6'},
        'SQ': [{'SN': chrom, 'LN': 3_000_000}],
    })
    read = pysam.AlignedSegment(hdr)
    read.query_name = 'test_read'
    read.reference_name = chrom
    read.reference_start = start
    read.cigartuples = cigar
    read.is_reverse = (strand == '-')
    read.is_unmapped = False
    read.is_secondary = False
    read.is_supplementary = False
    read.mapping_quality = 60
    read.query_sequence = seq
    return read


def _cigar_str(read: pysam.AlignedSegment) -> str:
    """Return CIGAR string from read."""
    return read.cigarstring or ''


def _has_terminal_D_or_N(read: pysam.AlignedSegment, strand: str) -> bool:
    """Return True if the alignment's 3' terminal op is D or N (invalid)."""
    cigar = list(read.cigartuples or [])
    if not cigar:
        return False
    # Strip terminal H/S (they're valid clip ops at the boundary)
    while cigar and cigar[-1][0] in (4, 5):
        cigar.pop()
    while cigar and cigar[0][0] in (4, 5):
        cigar.pop(0)
    if strand == '+':
        return bool(cigar) and cigar[-1][0] in (2, 3)
    else:
        return bool(cigar) and cigar[0][0] in (2, 3)


# ---------------------------------------------------------------------------
# Tests: clip_read_to_corrected_3prime — plus strand
# ---------------------------------------------------------------------------

class TestClipReadPlusStrand:
    """Tests for plus-strand hard-clipping."""

    def test_simple_clip(self):
        """Basic clip of trailing matches."""
        # CIGAR: 10= (ref 1000-1009), clip to 1005
        read = _make_read('chrI', 1000, [(8, 10)], '+', 'A' * 10)
        clipped = clip_read_to_corrected_3prime(read, 1005, '+')
        assert clipped
        assert read.reference_end == 1006  # last base at 1005, end = 1006
        assert 'H' in _cigar_str(read)
        assert not _has_terminal_D_or_N(read, '+')

    def test_no_clip_needed(self):
        """corrected_3prime at current boundary — no change."""
        read = _make_read('chrI', 1000, [(8, 10)], '+', 'A' * 10)
        clipped = clip_read_to_corrected_3prime(read, 1009, '+')
        assert not clipped

    def test_deletion_before_polya_trailing_match(self):
        """
        CIGAR: 3=4D6= at ref 1000.
        Deletion spans [1003, 1006]; last match ends at 1012.
        corrected_3prime = 1002 (last base of 3=).
        Walk-back removes 6= (query) and 4D (ref), leaving [3=].
        Result must NOT be [3=4D6H] — the 4D must be stripped before H.
        """
        # ref layout: 1000-1002 match (3=), 1003-1006 deleted (4D), 1007-1012 match (6=)
        seq = 'A' * 3 + 'G' * 6   # 9 query bases
        read = _make_read('chrI', 1000,
                          [(8, 3), (2, 4), (8, 6)],  # 3=, 4D, 6=
                          '+', seq)
        assert read.reference_end == 1013

        clipped = clip_read_to_corrected_3prime(read, 1002, '+')
        assert clipped
        assert not _has_terminal_D_or_N(read, '+'), (
            f'Terminal D/N found after clipping: {_cigar_str(read)}'
        )
        # The 3' alignment boundary should be at most 1002
        assert read.reference_end - 1 <= 1002

    def test_corrected_3prime_within_deletion_span(self):
        """
        corrected_3prime falls inside a deletion span (no query base there).
        CIGAR: 3=4D6= (same as above), corrected_3prime = 1005 (inside 4D).
        Walk removes 6= only (6 ref consumed = n_ref_clip), leaving [3=, 4D].
        Fix must strip the terminal 4D before appending H.
        """
        seq = 'A' * 3 + 'G' * 6
        read = _make_read('chrI', 1000,
                          [(8, 3), (2, 4), (8, 6)],
                          '+', seq)
        # corrected_3prime = 1006 → n_ref_clip = 1012 - 1006 = 6
        clipped = clip_read_to_corrected_3prime(read, 1006, '+')
        assert clipped
        assert not _has_terminal_D_or_N(read, '+'), (
            f'Terminal D after clip within deletion span: {_cigar_str(read)}'
        )

    def test_intron_skip_not_left_terminal(self):
        """
        CIGAR: 10=300N10= — spliced read.  If corrected to the first exon,
        the N must not be left as a terminal op.
        """
        seq = 'A' * 20
        read = _make_read('chrI', 1000,
                          [(8, 10), (3, 300), (8, 10)],  # 10=, 300N, 10=
                          '+', seq)
        # clip to end of first exon (ref 1009)
        clipped = clip_read_to_corrected_3prime(read, 1009, '+')
        assert clipped
        assert not _has_terminal_D_or_N(read, '+'), (
            f'Terminal N after clip: {_cigar_str(read)}'
        )

    def test_regression_299e1402_pattern(self):
        """
        Regression for read 299e1402 (chrII, plus strand).
        CIGAR pattern: ...5X3=4D6= where 6= is the poly-A region.
        corrected_3prime set to end of 3= block.
        Result must NOT end in [4D6H].
        """
        # Simplified version of the actual read geometry
        seq = 'T' * 5 + 'G' * 3 + 'A' * 6  # 5X + 3= + 6= = 14 query bases
        read = _make_read('chrII', 168400,
                          [(8, 5), (8, 3), (2, 4), (8, 6)],  # 5X, 3=, 4D, 6=
                          '+', seq)
        # corrected_3prime = end of 3= block = 168400+5+3-1 = 168407
        clipped = clip_read_to_corrected_3prime(read, 168407, '+')
        assert clipped
        cigar_str = _cigar_str(read)
        assert '4D' not in cigar_str or not cigar_str.endswith('H'), (
            f'Terminal 4D...H pattern still present: {cigar_str}'
        )
        assert not _has_terminal_D_or_N(read, '+'), (
            f'Terminal D/N after clip: {cigar_str}'
        )


# ---------------------------------------------------------------------------
# Tests: clip_read_to_corrected_3prime — minus strand
# ---------------------------------------------------------------------------

class TestClipReadMinusStrand:
    """Tests for minus-strand hard-clipping (clips from the left)."""

    def test_simple_clip_minus(self):
        """Basic clip of leading matches for minus strand."""
        # CIGAR: 10= (ref 1000-1009), clip 3' end (left) to 1005.
        # For minus strand, reference_start == corrected_3prime after clipping.
        read = _make_read('chrI', 1000, [(8, 10)], '-', 'A' * 10)
        clipped = clip_read_to_corrected_3prime(read, 1005, '-')
        assert clipped
        assert read.reference_start == 1005
        assert not _has_terminal_D_or_N(read, '-')

    def test_leading_deletion_not_left_terminal_minus(self):
        """
        CIGAR: 4D6= (ref 1000-...) — deletion at the very start.
        corrected_3prime = 1005 (inside the deletion after some matches are clipped).
        Leading D must not be left after clipping.
        """
        seq = 'G' * 6
        read = _make_read('chrI', 1000,
                          [(2, 4), (8, 6)],  # 4D, 6=
                          '-', seq)
        # corrected_3prime = 1006 → n_ref_clip = 1006 - 1000 = 6
        # Loop: process 4D (ref=4, query=0, need=6 → 4<=6, pop), n_ref=4
        #        process 6= (ref=6, query=6, need=2 → partial: trim to 4=), n_ref=6
        # cigar = [4=] (since D was fully popped and 6= trimmed to 4=)
        clipped = clip_read_to_corrected_3prime(read, 1006, '-')
        assert clipped
        assert not _has_terminal_D_or_N(read, '-'), (
            f'Terminal D/N on minus strand: {_cigar_str(read)}'
        )

    def test_deletion_within_deletion_span_minus(self):
        """
        CIGAR: 6=4D6= (minus strand), corrected_3prime within the 4D span.
        After loop removes left 6= (query), leading 4D should be stripped.
        """
        seq = 'A' * 6 + 'G' * 6
        read = _make_read('chrI', 1000,
                          [(8, 6), (2, 4), (8, 6)],  # 6=, 4D, 6=
                          '-', seq)
        # corrected_3prime = 1006 (last base of 4D span, inside deletion)
        # n_ref_clip = 1006 - 1000 = 6
        # Loop removes 6= (ref=6), exits. cigar = [4D, 6=]. Leading D!
        clipped = clip_read_to_corrected_3prime(read, 1006, '-')
        assert clipped
        assert not _has_terminal_D_or_N(read, '-'), (
            f'Leading D/N on minus strand after clip: {_cigar_str(read)}'
        )


# ---------------------------------------------------------------------------
# Tests: softclip_read_to_corrected_3prime
# ---------------------------------------------------------------------------

class TestSoftclipRead:
    """Tests for soft-clip variant (same D/N stripping required)."""

    def test_softclip_trailing_deletion_stripped_plus(self):
        """
        Same geometry as terminal-D test but using softclip variant.
        Result must not be [3=4D6S].
        """
        seq = 'A' * 3 + 'G' * 6
        read = _make_read('chrI', 1000,
                          [(8, 3), (2, 4), (8, 6)],
                          '+', seq)
        clipped = softclip_read_to_corrected_3prime(read, 1006, '+')
        assert clipped
        assert not _has_terminal_D_or_N(read, '+'), (
            f'Terminal D/N in softclip result: {_cigar_str(read)}'
        )

    def test_softclip_leading_deletion_stripped_minus(self):
        """Minus strand softclip: leading D must not survive."""
        seq = 'G' * 6 + 'A' * 6
        read = _make_read('chrI', 1000,
                          [(8, 6), (2, 4), (8, 6)],
                          '-', seq)
        clipped = softclip_read_to_corrected_3prime(read, 1006, '-')
        assert clipped
        assert not _has_terminal_D_or_N(read, '-'), (
            f'Leading D/N in softclip minus result: {_cigar_str(read)}'
        )


# ---------------------------------------------------------------------------
# Tests: trailing/leading I stripping (v3.2.7)
# ---------------------------------------------------------------------------

class TestTrailingInsertionStripping:
    """
    An alignment cannot end (hard-clip path) or have a trailing/leading I
    immediately adjacent to a soft-clip (soft-clip path) — the insertion bases
    must be absorbed into the clip op.

    Pattern: CIGAR ends in ...MIM and we clip to remove the last M.
    After the walk the rightmost op is I; it must be absorbed.
    """

    def test_clip_trailing_insertion_plus(self):
        """Plus strand hard-clip: trailing I must be absorbed into H.

        CIGAR 10M2I5M, start=1000, ref_end=1015, current_end=1014.
        Clip to corrected_3prime=1009 (remove last 5 ref bases of the 5M block).
        After the walk the 5M is gone; the 2I is now the rightmost op.
        Expected result: 10M7H (I absorbed), NOT 10M2I5H.
        """
        # CIGAR: 10M 2I 5M  →  clip last 5 ref bases  →  expect 10M7H (not 10M2I5H)
        seq = 'A' * 10 + 'T' * 2 + 'C' * 5
        read = _make_read('chrI', 1000,
                          [(0, 10), (1, 2), (0, 5)],
                          '+', seq)
        clipped = clip_read_to_corrected_3prime(read, 1009, '+')
        assert clipped
        cigar = list(read.cigartuples)
        # Last op must be H (hard-clip), not I or M+I
        assert cigar[-1][0] == 5, f'Expected H at end, got {_cigar_str(read)}'
        # I must not appear adjacent to H
        assert not any(op == 1 for op, _ in cigar[:-1][-1:]), (
            f'Trailing I before H in: {_cigar_str(read)}'
        )
        # Sequence length must match new CIGAR query span
        query_span = sum(l for op, l in cigar if op in (0, 1, 4, 7, 8))
        assert len(read.query_sequence) == query_span, (
            f'seq len {len(read.query_sequence)} != query span {query_span}'
        )

    def test_clip_leading_insertion_minus(self):
        """Minus strand hard-clip: leading I must be absorbed into H."""
        # CIGAR: 5M 2I 10M  →  clip first 5 ref bases from left  →  expect H for 7
        seq = 'C' * 5 + 'T' * 2 + 'A' * 10
        read = _make_read('chrI', 1000,
                          [(0, 5), (1, 2), (0, 10)],
                          '-', seq)
        clipped = clip_read_to_corrected_3prime(read, 1005, '-')
        assert clipped
        cigar = list(read.cigartuples)
        assert cigar[0][0] == 5, f'Expected leading H, got {_cigar_str(read)}'
        assert not any(op == 1 for op, _ in cigar[1:2]), (
            f'Leading I after H in: {_cigar_str(read)}'
        )
        query_span = sum(l for op, l in cigar if op in (0, 1, 4, 7, 8))
        assert len(read.query_sequence) == query_span

    def test_softclip_trailing_insertion_plus(self):
        """Plus strand soft-clip: trailing I must be absorbed into S."""
        seq = 'A' * 10 + 'T' * 2 + 'C' * 5
        read = _make_read('chrI', 1000,
                          [(0, 10), (1, 2), (0, 5)],
                          '+', seq)
        clipped = softclip_read_to_corrected_3prime(read, 1009, '+')
        assert clipped
        cigar = list(read.cigartuples)
        assert cigar[-1][0] == 4, f'Expected trailing S, got {_cigar_str(read)}'
        # No I adjacent to the S
        if len(cigar) >= 2:
            assert cigar[-2][0] != 1, f'I immediately before S in: {_cigar_str(read)}'

    def test_softclip_leading_insertion_minus(self):
        """Minus strand soft-clip: leading I must be absorbed into S."""
        seq = 'C' * 5 + 'T' * 2 + 'A' * 10
        read = _make_read('chrI', 1000,
                          [(0, 5), (1, 2), (0, 10)],
                          '-', seq)
        clipped = softclip_read_to_corrected_3prime(read, 1005, '-')
        assert clipped
        cigar = list(read.cigartuples)
        assert cigar[0][0] == 4, f'Expected leading S, got {_cigar_str(read)}'
        if len(cigar) >= 2:
            assert cigar[1][0] != 1, f'I immediately after leading S in: {_cigar_str(read)}'
