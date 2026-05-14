#!/usr/bin/env python3
"""Unit tests for the QuantSeq REV walkback wrapper.

QuantSeq REV chemistry is *antisense*: the cDNA reads are reverse-complementary
to the mRNA. After pysam-driven reverse-complementation:

    is_reverse = False  ->  gene is '-' strand, 3' end at the LEFT
    is_reverse = True   ->  gene is '+' strand, 3' end at the RIGHT

These tests mirror the structure of ``tests/test_walkback_readvsref.py``'s
``TestDrsWrapper`` but with the strand inversion that QuantSeq REV adds.
"""
from __future__ import annotations

import pysam
import pytest

from rectify.core.correct.protocols.quantseq_rev import walkback_quantseq_rev
from rectify.core.correct.walkback import APPLIED_NONE, APPLIED_WALKBACK


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


class TestQuantSeqRevStrandInversion:
    """Verify the BAM-strand -> gene-strand inversion."""

    def test_is_reverse_true_yields_plus_strand_and_right_side(self):
        """is_reverse=True -> gene is '+', 3' end at RIGHT side.

        Same V-primer-tip pattern as the core right-side test but routed
        through the QuantSeq REV wrapper.
        """
        # genome + 5 ref bases + 6 A's; read carries 5 ref + 5 A's + V-primer G
        ref = "X" * 1000 + "ACGTC" + "AAAAAA" + "X" * 100
        read = _make_read(
            start=1000,
            seq="ACGTC" + "AAAAA" + "G",
            cigar=((0, 11),),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_WALKBACK
        assert corr == 1004  # the C before the A-stretch
        assert orig == 1010  # rightmost ref position (reference_end - 1)

    def test_is_reverse_false_yields_minus_strand_and_left_side(self):
        """is_reverse=False -> gene is '-', 3' end at LEFT side.

        Mirror of the core left-side V-primer-tip test routed through the
        QuantSeq REV wrapper.
        """
        ref = "X" * 1000 + "AAAAAC" + "GTAC" + "X" * 100
        read = _make_read(
            start=1000,
            seq="GAAAAA" + "GTAC",
            cigar=((0, 10),),
            is_reverse=False,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "-"
        assert applied == APPLIED_WALKBACK
        assert orig == 1000
        assert corr == 1006

    def test_inversion_is_opposite_of_drs(self):
        """Sanity check: a read with is_reverse=True should map to gene '+'
        in QuantSeq REV but to gene '-' in DRS. Use the DRS wrapper for
        the contrast.
        """
        from rectify.core.correct.walkback import walkback_drs

        ref = "X" * 1000 + "ACGTC" + "AAA" + "TT" + "X" * 100
        rd_qs = _make_read(
            start=1000, seq="ACGTC" + "AAAAA", cigar=((0, 10),), is_reverse=True
        )
        rd_drs = _make_read(
            start=1000, seq="ACGTC" + "AAAAA", cigar=((0, 10),), is_reverse=True
        )
        _, _, _, qs_strand = walkback_quantseq_rev(rd_qs, ref)
        _, _, _, drs_strand = walkback_drs(rd_drs, ref)
        assert qs_strand == "+"
        assert drs_strand == "-"
        assert qs_strand != drs_strand


class TestQuantSeqRevWalkbackFires:
    """Walkback fires on internal-priming over a genomic A-run."""

    def test_over_extension_into_genomic_a_tract_right_side(self):
        """is_reverse=True, gene '+', poly-A over-extension at the RIGHT.

        Genome:  ...ACGTC AAA TT X...  (3 genomic A's then T at pos 1008)
        Read:    ...ACGTC AAA AA       (5 A's; last A over a genomic T -> mismatch)
        """
        ref = "X" * 1000 + "ACGTC" + "AAA" + "TT" + "X" * 100
        read = _make_read(
            start=1000,
            seq="ACGTC" + "AAA" + "AA",
            cigar=((0, 10),),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_WALKBACK
        assert corr == 1004  # the C before the A-stretch
        assert orig == 1009

    def test_terminal_match_no_walkback(self):
        """Non-A terminal matching genome -> no walkback (Case 1 gate)."""
        ref = "X" * 1000 + "ACGTACGTAC" + "X" * 100
        read = _make_read(
            start=1000,
            seq="ACGTACGTAC",
            cigar=((0, 10),),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_NONE
        assert orig == corr == 1009

    def test_internal_priming_terminal_a_on_genomic_a_fires_walkback(self):
        """Core RECTIFY case for QuantSeq REV: poly-A aligning to genomic A-tract.

        is_reverse=True → gene '+', 3' end at RIGHT.
        Genome:  ...ACGTC AAAAA X...
        Read:    ...ACGTC AAAAA      (all A's on genomic A; walkback must fire)
        The scan anchors at C (pos 1004), the last non-A read-genome match.
        """
        ref = "X" * 1000 + "ACGTC" + "AAAAA" + "X" * 100
        read = _make_read(
            start=1000,
            seq="ACGTC" + "AAAAA",
            cigar=((0, 10),),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_WALKBACK
        assert orig == 1009
        assert corr == 1004  # C before the A-stretch

    def test_internal_priming_left_side_terminal_a_on_genomic_a(self):
        """Core RECTIFY case on LEFT side: gene '-', poly-A aligning to genomic A.

        is_reverse=False → gene '-', 3' end at LEFT.
        Genome:  ...AAAAA CGTAC X...
        Read:    ...AAAAA CGTAC       (all leftmost A's on genomic A; walkback fires)
        Scan from left anchors at C (pos 1005), the first non-A inward agreement.
        """
        ref = "X" * 1000 + "AAAAA" + "CGTAC" + "X" * 100
        read = _make_read(
            start=1000,
            seq="AAAAA" + "CGTAC",
            cigar=((0, 10),),
            is_reverse=False,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "-"
        assert applied == APPLIED_WALKBACK
        assert orig == 1000
        assert corr == 1005  # C immediately after the A-stretch


class TestQuantSeqRevBoundary:
    """Edge cases the wrapper should pass through to the core."""

    def test_unmapped_no_aligned_pairs(self):
        """Read with all soft-clip CIGAR -> wrapper returns APPLIED_NONE."""
        ref = "A" * 100
        read = _make_read(start=10, seq="AAAA", cigar=((4, 4),), is_reverse=False)
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert applied == APPLIED_NONE
        assert orig == corr
        assert gene_strand == "-"
