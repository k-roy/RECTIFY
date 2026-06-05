#!/usr/bin/env python3
"""Unit tests for the NET-seq walkback wrapper (rectify.core.correct.protocols.netseq).

NET-seq (yeast single-end + human mNET-seq read 2) is antisense; its poly-A
walkback geometry is IDENTICAL to QuantSeq REV. These tests mirror
``tests/test_quantseq_rev_walkback.py`` and assert ``walkback_netseq`` reproduces
the QuantSeq REV walkback exactly (same corrected CPA, same gene strand) for both
``is_reverse`` states — locking the protocol contract so future edits can't
silently diverge the validated geometry.
"""
from __future__ import annotations

import pysam

from rectify.core.correct.protocols.netseq import walkback_netseq
from rectify.core.correct.protocols.quantseq_rev import walkback_quantseq_rev
from rectify.core.correct.walkback import APPLIED_WALKBACK


def _make_read(*, chrom="chrI", start=1000, seq="", cigar=((0, 0),), is_reverse=False):
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": 100_000}]}
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


def test_netseq_plus_strand_right_side():
    """is_reverse=True -> gene '+', poly-A at RIGHT; CPA walks back to the C before the A-tract."""
    ref = "X" * 1000 + "ACGTC" + "AAAAAA" + "X" * 100
    read = _make_read(start=1000, seq="ACGTC" + "AAAAA" + "G", cigar=((0, 11),), is_reverse=True)
    orig, corr, applied, gene_strand = walkback_netseq(read, ref)
    assert gene_strand == "+"
    assert applied == APPLIED_WALKBACK
    assert corr == 1004   # the C immediately 5' of the genomic A-stretch = true CPA
    assert orig == 1010   # rightmost aligned ref position (reference_end - 1)


def test_netseq_matches_quantseq_rev_both_strands():
    """NET-seq walkback must reproduce QuantSeq REV exactly (identical antisense geometry)."""
    ref = "X" * 1000 + "ACGTC" + "AAAAAA" + "X" * 100
    for is_rev in (True, False):
        r_net = _make_read(start=1000, seq="ACGTC" + "AAAAA" + "G", cigar=((0, 11),), is_reverse=is_rev)
        r_qsr = _make_read(start=1000, seq="ACGTC" + "AAAAA" + "G", cigar=((0, 11),), is_reverse=is_rev)
        assert walkback_netseq(r_net, ref) == walkback_quantseq_rev(r_qsr, ref)
