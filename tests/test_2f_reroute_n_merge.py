"""ISSUE-007: the reroute's N-merge must not fuse a PARTIALLY TRIMMED junction.

`reroute_intronic_tail_5prime_via_junction` trims the CIGAR back to
``clip_boundary`` and appends/prepends the new N, then merges consecutive N ops.
The merge was written for the legitimate Case-5 abutting case — a FULLY INTACT
pre-existing N ending exactly at ``clip_boundary``. When ``clip_boundary`` falls
INSIDE an aligner-called junction the trim shortens it and the merge then fuses
the remnant with the new N, silently deleting every junction in between: read
eca6079d had a 1,241 bp N cut to 755 bp and merged into a 97,213 bp AG-AG
"intron" that replaced six real junctions.

The guard was keyed on the N being CUT (``op == 3 and trim < length``), not on
the terminal op being N — the latter would break Case 5. `grep -rn
n_boundary_adjust` found only the two code comments, so nothing else covers that
case; the abutting test below is the first. ISSUE-026 invariant B (2026-09-05)
widened the guard to ``op == 3``: an N wholly inside the trimmed run is refused
too (a 5' rescue may never remove an aligner-called N-op). Case 5's abutting N
still merges — it ends exactly at ``clip_boundary`` and the trim loop breaks
before ever reaching it.

Run with:
    pytest tests/test_2f_reroute_n_merge.py -v
"""

from typing import List, Tuple

import pysam
import pytest

from rectify.core.bam.read_edits import (
    reroute_intronic_tail_5prime_via_junction,
    softclip_intronic_tail_5prime,
)


def _make_read(start: int, cigar: List[Tuple[int, int]], strand: str,
               name: str = 'r') -> pysam.AlignedSegment:
    hdr = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6'},
        'SQ': [{'SN': 'chrN', 'LN': 1_000_000}],
    })
    read = pysam.AlignedSegment(hdr)
    read.query_name = name
    read.reference_name = 'chrN'
    read.reference_start = start
    read.cigartuples = cigar
    read.is_reverse = (strand == '-')
    read.is_unmapped = False
    read.is_secondary = False
    read.is_supplementary = False
    read.mapping_quality = 60
    read.query_sequence = 'A' * sum(l for op, l in cigar if op in (0, 1, 4, 7, 8))
    return read


def _n_ops(read):
    pos = read.reference_start
    out = []
    for op, length in read.cigartuples:
        if op == 3:
            out.append((pos, pos + length))
        if op in (0, 2, 3, 7, 8):
            pos += length
    return out


class TestMinusStrand:
    """Read at 100: 50M (100-150) 200N (150-350) 30M (350-380).
    The 5' end is the HIGH-coordinate end."""

    def _read(self):
        return _make_read(100, [(0, 50), (3, 200), (0, 30)], '-')

    def test_abutting_intact_n_still_merges(self):
        """Case 5 `n_boundary_adjust`: clip_boundary == the N's exon-2 edge, so
        the trim never touches the N and the two N ops are one boundary-adjusted
        junction. This must keep working."""
        read = self._read()
        assert reroute_intronic_tail_5prime_via_junction(
            read, clip_boundary=350, five_prime_position=400,
            exon_cigar_str='30M', strand='-') is True
        assert _n_ops(read) == [(150, 400)]
        assert read.cigarstring == '50M250N30M'

    def test_cut_n_is_refused(self):
        """clip_boundary 300 falls inside the 150-350 junction."""
        read = self._read()
        before_cigar = read.cigartuples
        before_start = read.reference_start
        assert reroute_intronic_tail_5prime_via_junction(
            read, clip_boundary=300, five_prime_position=400,
            exon_cigar_str='30M', strand='-') is False
        assert read.cigartuples == before_cigar
        assert read.reference_start == before_start

    def test_refusal_leaves_no_fused_giant_intron(self):
        """Without the guard the read would carry a single 150-400 N in place of
        its real 150-350 junction."""
        read = self._read()
        reroute_intronic_tail_5prime_via_junction(
            read, clip_boundary=300, five_prime_position=400,
            exon_cigar_str='30M', strand='-')
        assert _n_ops(read) == [(150, 350)]

    def test_eca6079d_shape_is_refused(self):
        """The reported carrier: a 1,241 bp N cut and fused into 97,213 bp."""
        read = _make_read(0, [(0, 100), (3, 1241), (0, 200)], '-')
        assert reroute_intronic_tail_5prime_via_junction(
            read, clip_boundary=100 + 755, five_prime_position=97_500,
            exon_cigar_str='200M', strand='-') is False
        assert _n_ops(read) == [(100, 1341)]


class TestPlusStrand:
    """Read at 100: 30M (100-130) 200N (130-330) 50M (330-380).
    The 5' end is the LOW-coordinate end."""

    def _read(self):
        return _make_read(100, [(0, 30), (3, 200), (0, 50)], '+')

    def test_cut_n_is_refused(self):
        read = self._read()
        before_cigar = read.cigartuples
        before_start = read.reference_start
        assert reroute_intronic_tail_5prime_via_junction(
            read, clip_boundary=300, five_prime_position=49,
            exon_cigar_str='30M', strand='+') is False
        assert read.cigartuples == before_cigar
        assert read.reference_start == before_start
        assert _n_ops(read) == [(130, 330)]

    def test_intact_n_fully_inside_the_trimmed_region_is_refused(self):
        """ISSUE-026 invariant B (2026-09-05): an N the trim would consume
        ENTIRELY (trim == length) is refused too — a 5' rescue may never remove
        an aligner-called N-op. (Until B this was the `_get_intronic_query_bases`
        pass-through for a "spurious micro-junction inside the target intron";
        T0 chrX 771560c4 showed the same geometry deleting a REAL 4,073-nt
        annotated junction and its 123-nt terminal exon — an exon-skip collapse.)"""
        read = self._read()
        before_cigar, before_start = read.cigartuples, read.reference_start
        assert reroute_intronic_tail_5prime_via_junction(
            read, clip_boundary=330, five_prime_position=49,
            exon_cigar_str='30M', strand='+') is False
        assert read.cigartuples == before_cigar and read.reference_start == before_start
        assert _n_ops(read) == [(130, 330)]


class TestInvariantBNeverRemovesAnAlignerN:
    """ISSUE-026 invariant B: neither writer helper may delete an aligner-called
    N-op, whether the boundary CUTS it (ISSUE-007) or the re-placed run CONTAINS
    it whole. The minus-strand shape is T0 chrX 771560c4: read 50M 200N 30M at
    100 (junction 150-350, terminal exon 350-380), candidate intron 120-500 —
    clip_boundary 120 sits before the read's own junction."""

    def test_minus_reroute_refuses_a_contained_junction(self):
        read = _make_read(100, [(0, 50), (3, 200), (0, 30)], '-')
        before = read.cigartuples
        assert reroute_intronic_tail_5prime_via_junction(
            read, clip_boundary=120, five_prime_position=500,
            exon_cigar_str='60M', strand='-') is False
        assert read.cigartuples == before
        assert _n_ops(read) == [(150, 350)]

    def test_minus_softclip_refuses_a_contained_junction(self):
        read = _make_read(100, [(0, 50), (3, 200), (0, 30)], '-')
        before = read.cigartuples
        assert softclip_intronic_tail_5prime(read, clip_boundary=120, strand='-') is False
        assert read.cigartuples == before
        assert _n_ops(read) == [(150, 350)]

    def test_plus_softclip_refuses_a_contained_junction(self):
        read = _make_read(100, [(0, 30), (3, 200), (0, 50)], '+')
        before = read.cigartuples
        assert softclip_intronic_tail_5prime(read, clip_boundary=360, strand='+') is False
        assert read.cigartuples == before
        assert _n_ops(read) == [(130, 330)]


class TestNoQueryCorruption:
    """A refusal must not consume query bases: the read is returned untouched."""

    def test_refusal_preserves_query_length_accounting(self):
        read = _make_read(100, [(0, 50), (3, 200), (0, 30)], '-')
        q_before = sum(l for op, l in read.cigartuples
                       if op in (0, 1, 4, 7, 8))
        reroute_intronic_tail_5prime_via_junction(
            read, clip_boundary=300, five_prime_position=400,
            exon_cigar_str='30M', strand='-')
        q_after = sum(l for op, l in read.cigartuples if op in (0, 1, 4, 7, 8))
        assert q_after == q_before == len(read.query_sequence)


class TestSoftclipIntronicTailAlsoRefusesACutJunction:
    """ISSUE-007's second call site, found via Sumner hold-out read 34625d8e.

    When the icp gate routes away from `extend` and the reroute refuses, the
    writer falls back to `softclip_intronic_tail_5prime` — which trims to the
    same clip_boundary with the same `min(length, excess)` and so could shorten
    an aligner-called junction on its own. On 34625d8e that turned an annotated
    CT-AC 91378775-91382812 into 91378775-91380924; the canonical-destination
    guard could not see it, because both edges are still CT-AC.
    """

    def test_minus_strand_cut_is_refused(self):
        read = _make_read(100, [(0, 50), (3, 200), (0, 30)], '-')
        before = read.cigartuples
        assert softclip_intronic_tail_5prime(read, clip_boundary=300, strand='-') is False
        assert read.cigartuples == before
        assert _n_ops(read) == [(150, 350)]

    def test_plus_strand_cut_is_refused(self):
        read = _make_read(100, [(0, 30), (3, 200), (0, 50)], '+')
        before = read.cigartuples
        assert softclip_intronic_tail_5prime(read, clip_boundary=300, strand='+') is False
        assert read.cigartuples == before
        assert _n_ops(read) == [(130, 330)]

    def test_a_boundary_outside_any_junction_still_soft_clips(self):
        """The guard must not disable the helper: a clip_boundary that falls in
        an M op trims normally."""
        read = _make_read(100, [(0, 50), (3, 200), (0, 30)], '-')
        assert softclip_intronic_tail_5prime(read, clip_boundary=360, strand='-') is True
        assert _n_ops(read) == [(150, 350)]      # junction intact
        assert read.cigartuples[-1][0] == 4      # trailing soft clip added
