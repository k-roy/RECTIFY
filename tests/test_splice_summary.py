#!/usr/bin/env python3
"""Unit tests for ``rectify.core.analyze.splice_summary``."""
from __future__ import annotations

import pysam
import pytest

from rectify.core.analyze.splice_summary import (
    CLASS_ALTERNATIVE,
    CLASS_ANNOTATED,
    CLASS_NO_SPAN,
    CLASS_NOVEL,
    CLASS_UNSPLICED,
    GeneIntronSet,
    classify_read,
    extract_n_ops,
    snap_n_op_to_annotated,
)


def _read(*, start: int, cigar, is_reverse: bool = False, chrom: str = "chrI") -> pysam.AlignedSegment:
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": 3_000_000}]}
    )
    r = pysam.AlignedSegment(hdr)
    r.query_name = "t"
    r.reference_name = chrom
    r.reference_start = start
    r.cigartuples = list(cigar)
    r.is_reverse = is_reverse
    r.is_unmapped = False
    r.is_secondary = False
    r.is_supplementary = False
    r.mapping_quality = 60
    # Build a query sequence of the right length (1 per M/I/S, 0 per D/N).
    qlen = sum(L for op, L in cigar if op in (0, 1, 4, 7, 8))
    r.query_sequence = "A" * qlen
    return r


def _gene(introns):
    return GeneIntronSet(
        gene_id="YDR381W",
        chrom="chrI",
        strand="+",
        introns=frozenset(introns),
    )


# ---------------------------------------------------------------------------
# extract_n_ops
# ---------------------------------------------------------------------------
class TestExtractNOps:
    def test_no_n_op(self):
        r = _read(start=100, cigar=[(0, 50)])
        assert extract_n_ops(r) == []

    def test_single_n_op(self):
        # 100 + 30M -> 130; N=200 -> intron [130, 330); 20M -> 350
        r = _read(start=100, cigar=[(0, 30), (3, 200), (0, 20)])
        assert extract_n_ops(r) == [(130, 330)]

    def test_two_n_ops(self):
        r = _read(start=100, cigar=[(0, 30), (3, 200), (0, 50), (3, 100), (0, 20)])
        assert extract_n_ops(r) == [(130, 330), (380, 480)]

    def test_deletion_advances_ref_pos(self):
        # 100 + 10M -> 110; 5D -> 115; N=100 -> [115, 215)
        r = _read(start=100, cigar=[(0, 10), (2, 5), (3, 100), (0, 20)])
        assert extract_n_ops(r) == [(115, 215)]


# ---------------------------------------------------------------------------
# snap_n_op_to_annotated
# ---------------------------------------------------------------------------
class TestSnap:
    def test_exact_match_returns_same(self):
        assert snap_n_op_to_annotated((100, 200), frozenset([(100, 200)]), 5) == (100, 200)

    def test_within_tolerance_both_ends(self):
        # N-op (102, 198), annotated (100, 200) — within 5 bp on both sides
        assert snap_n_op_to_annotated((102, 198), frozenset([(100, 200)]), 5) == (100, 200)

    def test_one_end_out_of_tolerance(self):
        # Donor 7 bp off — exceeds tolerance — no snap
        assert snap_n_op_to_annotated((107, 200), frozenset([(100, 200)]), 5) == (107, 200)

    def test_tolerance_zero_no_snap(self):
        assert snap_n_op_to_annotated((101, 200), frozenset([(100, 200)]), 0) == (101, 200)

    def test_picks_closest_when_two_candidates(self):
        annotated = frozenset([(100, 200), (120, 220)])
        # (115, 215) is 5+5=10 from (120,220) and 15+15=30 from (100,200)
        assert snap_n_op_to_annotated((115, 215), annotated, 10) == (120, 220)


# ---------------------------------------------------------------------------
# classify_read
# ---------------------------------------------------------------------------
class TestClassifyRead:
    # Use intron [1000, 1200) for a hypothetical gene.
    GENE = _gene([(1000, 1200)])

    def test_unspliced_intron_retained(self):
        # Read fully covers the intron, no N-cigar.
        r = _read(start=900, cigar=[(0, 400)])  # 900..1300
        assert classify_read(r, self.GENE) == CLASS_UNSPLICED

    def test_annotated_spliced_exact_match(self):
        # exon1: 900..1000, intron N: 1000..1200, exon2: 1200..1300
        r = _read(start=900, cigar=[(0, 100), (3, 200), (0, 100)])
        assert classify_read(r, self.GENE) == CLASS_ANNOTATED

    def test_alternative_5ss_one_end_matches(self):
        # Donor 1 bp downstream of annotated, acceptor at annotated
        # exon1 ends at 1001 (not 1000), intron 1001..1200, exon2 starts 1200
        r = _read(start=900, cigar=[(0, 101), (3, 199), (0, 100)])
        # N-op = (1001, 1200): donor 1001 not annotated, acceptor 1200 IS annotated
        assert classify_read(r, self.GENE) == CLASS_ALTERNATIVE

    def test_alternative_3ss_one_end_matches(self):
        # Donor at annotated 1000, acceptor at 1205 (5 bp downstream of annotated 1200)
        r = _read(start=900, cigar=[(0, 100), (3, 205), (0, 100)])
        # N-op = (1000, 1205): donor 1000 IS annotated, acceptor 1205 is not
        assert classify_read(r, self.GENE) == CLASS_ALTERNATIVE

    def test_novel_neither_end_matches(self):
        # Both ends are 10 bp off
        r = _read(start=900, cigar=[(0, 110), (3, 180), (0, 100)])
        # N-op = (1010, 1190): neither end annotated
        assert classify_read(r, self.GENE) == CLASS_NOVEL

    def test_no_intron_span(self):
        # Read entirely upstream of intron
        r = _read(start=500, cigar=[(0, 100)])
        assert classify_read(r, self.GENE) == CLASS_NO_SPAN

    def test_snapping_rescues_basecaller_noise(self):
        # N-op 2 bp off — without snapping it's novel, with snapping it's annotated.
        r = _read(start=900, cigar=[(0, 102), (3, 196), (0, 100)])
        # N-op = (1002, 1198) — within 5 bp tolerance of (1000, 1200)
        assert classify_read(r, self.GENE, snap_tolerance=0) == CLASS_NOVEL
        assert classify_read(r, self.GENE, snap_tolerance=5) == CLASS_ANNOTATED

    def test_two_intron_gene_one_alt_one_annotated(self):
        gene = _gene([(1000, 1100), (1200, 1300)])
        # First N matches annotated; second N has acceptor 1305 (alt 3'SS)
        # exon1: 900..1000, N1: 1000..1100, exon2: 1100..1200, N2: 1200..1305, exon3: 1305..1400
        r = _read(start=900, cigar=[(0, 100), (3, 100), (0, 100), (3, 105), (0, 95)])
        assert classify_read(r, gene) == CLASS_ALTERNATIVE

    def test_two_intron_gene_both_annotated(self):
        gene = _gene([(1000, 1100), (1200, 1300)])
        r = _read(start=900, cigar=[(0, 100), (3, 100), (0, 100), (3, 100), (0, 100)])
        assert classify_read(r, gene) == CLASS_ANNOTATED
