#!/usr/bin/env python3
"""Unit tests for ``rectify.core.correct.walkback``.

The walk-back algorithm under test compares the read base to the reference
base at each aligned position and stops at the first non-stop-base match.
These tests use synthetic ``pysam.AlignedSegment`` objects so no real BAM
or FASTA on disk is required.
"""
from __future__ import annotations

import pysam
import pytest

from rectify.core.correct.walkback import (
    APPLIED_NONE,
    APPLIED_WALKBACK,
    THREE_PRIME_SIDE_LEFT,
    THREE_PRIME_SIDE_RIGHT,
    walkback_3prime,
    walkback_drs,
)


def _make_read(
    *,
    chrom: str = "chrI",
    start: int = 1000,
    seq: str = "",
    cigar=((0, 0),),
    is_reverse: bool = False,
) -> pysam.AlignedSegment:
    """Build a minimal ``pysam.AlignedSegment`` for unit testing."""
    hdr = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": chrom, "LN": 100_000}],
        }
    )
    r = pysam.AlignedSegment(hdr)
    r.query_name = "test"
    r.reference_name = chrom
    r.reference_start = start
    r.cigartuples = list(cigar)
    r.is_reverse = is_reverse
    r.is_unmapped = False
    r.is_secondary = False
    r.is_supplementary = False
    r.mapping_quality = 60
    r.query_sequence = seq
    return r


# ---------------------------------------------------------------------------
# Terminal-position gating (no walk-back needed)
# ---------------------------------------------------------------------------
class TestTerminalGate:
    def test_terminal_nonstop_match_no_walkback(self):
        """Terminal non-A read base matching genome → correctly anchored."""
        # 10 bp aligned, read ends in C, reference at that pos is also C.
        ref = "X" * 1000 + "ACGTACGTAC" + "X" * 100
        read = _make_read(start=1000, seq="ACGTACGTAC", cigar=((0, 10),))
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT)
        assert applied == APPLIED_NONE
        assert orig == corr == 1009

    def test_terminal_stop_base_matching_genome_walkback_fires(self):
        """Terminal A on genomic A → walkback fires to find true CPA.

        Genome:  ...CGTACGT AAA
        Read:    ...CGTACGT AAA
        The last 3 A's could be poly-A tail aligned over the genomic A-run.
        Walkback scans inward and anchors at T (pos 1006), the last non-A
        read-genome agreement.
        """
        ref = "X" * 1000 + "CGTACGT" + "AAA" + "X" * 100
        read = _make_read(start=1000, seq="CGTACGT" + "AAA", cigar=((0, 10),))
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT)
        assert applied == APPLIED_WALKBACK
        assert orig == 1009
        assert corr == 1006  # T immediately before the A-stretch


# ---------------------------------------------------------------------------
# Walk-back actually fires when poly-A over-extends into genomic A-stretch
# ---------------------------------------------------------------------------
class TestWalkbackFires:
    def test_polya_extended_into_genomic_a_stretch_right(self):
        """Core RECTIFY case: poly-A tail aligning into a genomic A-tract.

        Genome:  ...ACGTC AAAAA...
        Read:    ...ACGTC AAAAA   (last 5 A's are basecalled tail over genomic A)
        The terminal A matches the genomic A, but walkback MUST still fire —
        non-genomically encoded poly-A routinely lands on genomic A-runs.
        The scan walks back to C (the last non-A read-genome agreement).
        """
        ref = "X" * 1000 + "ACGTC" + "AAAAA" + "X" * 100
        read = _make_read(start=1000, seq="ACGTC" + "AAAAA", cigar=((0, 10),))
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT)
        assert applied == APPLIED_WALKBACK
        assert orig == 1009
        assert corr == 1004  # C immediately before the A-stretch

    def test_polya_extends_past_genomic_tract_mismatching_terminal(self):
        """Unambiguous over-extension: read has more A's than the genomic A-run.

        Genome:  ...ACGTC AAA T  X...   (3 A's then T at position 1008)
        Read:    ...ACGTC AAAA A         (4 A's then A at position 1008 — mismatch!)
        The last A in the read is basecalled tail aligned over a genomic T;
        walkback finds the canonical anchor at the C (position 1004).
        """
        ref = "X" * 1000 + "ACGTC" + "AAA" + "TT" + "X" * 100
        read = _make_read(start=1000, seq="ACGTC" + "AAA" + "AA", cigar=((0, 10),))
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT)
        assert applied == APPLIED_WALKBACK
        assert corr == 1004  # the C before the A-stretch
        assert orig == 1009

    def test_polya_v_primer_tip_artifact(self):
        """QuantSeq V-primer artifact: read ends in ...AAAAAAAG (terminal G)
        over a genomic A-run. Terminal G mismatches genomic A → walkback fires,
        scan continues through the A's to the canonical anchor.
        """
        ref = "X" * 1000 + "ACGTC" + "AAAAAA" + "X" * 100  # 5 ref bases + 6 A's
        read = _make_read(
            start=1000,
            seq="ACGTC" + "AAAAA" + "G",  # 5 ref + 5 A's + 1 G (artifact)
            cigar=((0, 11),),
        )
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT)
        assert applied == APPLIED_WALKBACK
        assert corr == 1004  # the C before the A-stretch

    def test_polya_extension_left_side_v_primer_artifact(self):
        """Mirror of the V-primer-tip artifact on the LEFT side of the alignment.

        For a minus-strand gene under a left-3'-end protocol, the read's
        terminal base (at ``reference_start``) is the V-primer tip (a G or
        other non-A) sitting over a genomic A. That terminal mismatch fires
        Case 3 of the gate and the walkback proceeds.

        ref pos 1000+:  A A A A A C G T A C
        read         :  G A A A A A G T A C   (terminal G is the artifact)
        Scan from left:
          (0,1000) G/A mismatch → skip
          (1,1001) A/A stop_base → skip
          (2,1002) A/A stop_base → skip
          (3,1003) A/A stop_base → skip
          (4,1004) A/A stop_base → skip
          (5,1005) A/C mismatch → skip
          (6,1006) G/G MATCH, non-stop → anchor at 1006.
        """
        ref = "X" * 1000 + "AAAAAC" + "GTAC" + "X" * 100
        read = _make_read(start=1000, seq="GAAAAA" + "GTAC", cigar=((0, 10),))
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_LEFT)
        assert applied == APPLIED_WALKBACK
        assert orig == 1000
        assert corr == 1006


# ---------------------------------------------------------------------------
# Bad input handling
# ---------------------------------------------------------------------------
class TestInputValidation:
    def test_bad_three_prime_side_raises(self):
        ref = "A" * 100
        read = _make_read(start=10, seq="ACGT", cigar=((0, 4),))
        with pytest.raises(ValueError, match="three_prime_side"):
            walkback_3prime(read, ref, "middle")

    def test_bad_stop_base_raises(self):
        ref = "A" * 100
        read = _make_read(start=10, seq="ACGT", cigar=((0, 4),))
        with pytest.raises(ValueError, match="stop_base"):
            walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT, stop_base="X")

    def test_unmapped_read_returns_unchanged(self):
        """Read with no aligned pairs → no correction, no crash."""
        ref = "A" * 100
        # CIGAR all soft-clip → no aligned (qp, rp) pairs.
        read = _make_read(start=10, seq="AAAA", cigar=((4, 4),))
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT)
        assert applied == APPLIED_NONE
        assert orig == corr


# ---------------------------------------------------------------------------
# DRS wrapper: BAM strand passes through to gene strand
# ---------------------------------------------------------------------------
class TestDrsWrapper:
    def test_drs_plus_strand_uses_right_side(self):
        ref = "X" * 1000 + "ACGTC" + "AAA" + "TT" + "X" * 100
        read = _make_read(
            start=1000, seq="ACGTC" + "AAAAA", cigar=((0, 10),), is_reverse=False
        )
        orig, corr, applied, strand = walkback_drs(read, ref)
        assert strand == "+"
        assert applied == APPLIED_WALKBACK
        assert corr == 1004

    def test_drs_minus_strand_uses_left_side(self):
        """Same V-primer-tip pattern as the left-side core test, exercised
        through the DRS wrapper (``is_reverse=True`` → minus-strand gene)."""
        ref = "X" * 1000 + "AAAAAC" + "GTAC" + "X" * 100
        read = _make_read(
            start=1000, seq="GAAAAA" + "GTAC", cigar=((0, 10),), is_reverse=True
        )
        orig, corr, applied, strand = walkback_drs(read, ref)
        assert strand == "-"
        assert applied == APPLIED_WALKBACK
        assert corr == 1006
