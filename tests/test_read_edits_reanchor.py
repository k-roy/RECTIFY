"""
Unit tests for _apply_reanchor_from_clip_len vs reanchor_5prime_for_rescue.

For every synthetic read, both functions must produce byte-identical
(cigartuples, reference_start) output. Covers all edge cases identified
in the plan: leading H, leading/trailing S already present, anchor at op
boundary, anchor mid-op, D/N before anchor — on both strands.

Run with:
    pytest tests/test_read_edits_reanchor.py -v
"""

import copy
import pytest
import pysam
from typing import Dict, List, Tuple, Optional

from rectify.core.bam.read_edits import (
    reanchor_5prime_for_rescue,
    _apply_reanchor_from_clip_len,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_CHROM = 'chrI'
_CHROM_LEN = 2_000_000


def _make_read(
    start: int,
    cigar: List[Tuple[int, int]],
    seq: str,
    strand: str = '+',
) -> pysam.AlignedSegment:
    hdr = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6'},
        'SQ': [{'SN': _CHROM, 'LN': _CHROM_LEN}],
    })
    r = pysam.AlignedSegment(hdr)
    r.query_name = 'test_read'
    r.reference_name = _CHROM
    r.reference_start = start
    r.cigartuples = cigar
    r.is_reverse = (strand == '-')
    r.is_unmapped = False
    r.is_secondary = False
    r.is_supplementary = False
    r.mapping_quality = 60
    r.query_sequence = seq
    return r


def _copy_read(read: pysam.AlignedSegment) -> pysam.AlignedSegment:
    """Return a fresh read with the same state (cigar + start + seq + strand)."""
    r = _make_read(
        read.reference_start,
        list(read.cigartuples),
        read.query_sequence,
        '-' if read.is_reverse else '+',
    )
    return r


def _build_genome(seq: str, start: int = 0) -> Dict[str, str]:
    """Build a minimal genome dict with a single chromosome."""
    chrom_seq = 'N' * start + seq + 'N' * (max(0, _CHROM_LEN - start - len(seq)))
    return {_CHROM: chrom_seq[:_CHROM_LEN]}


def _assert_reanchor_equivalent(
    read_orig: pysam.AlignedSegment,
    genome: Dict[str, str],
    clip_len: int,
    label: str,
) -> None:
    """Run both implementations and assert identical output."""
    r_full = _copy_read(read_orig)
    r_fast = _copy_read(read_orig)

    changed_full = reanchor_5prime_for_rescue(r_full, genome, anchor_min_run=10)
    changed_fast = _apply_reanchor_from_clip_len(r_fast, clip_len)

    assert changed_full is True, f'{label}: full reanchor did not fire'
    assert changed_fast is True, f'{label}: fast reanchor did not fire'

    assert list(r_fast.cigartuples) == list(r_full.cigartuples), (
        f'{label}: cigar mismatch\n'
        f'  full: {r_full.cigarstring}\n'
        f'  fast: {r_fast.cigarstring}'
    )
    assert r_fast.reference_start == r_full.reference_start, (
        f'{label}: reference_start mismatch: full={r_full.reference_start} fast={r_fast.reference_start}'
    )


# ---------------------------------------------------------------------------
# Genome for + strand tests: starts at position 1000.
# First 5 bases are mismatches (X cluster), then 20 matching bases (anchor).
#
# Read seq:  XXXXXAAAAAAAAAAAAAAAAAAA  (5 mismatches + 20 A's)
# Ref  seq:  BBBBBAAAAAAAAAAAAAAAAAAA  (5 B's + 20 A's at ref pos 1000)
# Cigar: 25M — the first 5 M bases are mismatches; bases 5-24 match.
# Expected reanchor: leading 5S, then 20M, reference_start shifts +5.
# ---------------------------------------------------------------------------

_REF_START = 1000
_MISMATCH_LEN = 5
_ANCHOR_LEN = 20
_ANCHOR_BASE = 'A'
_MISMATCH_BASE_READ = 'C'   # mismatch in read
_MISMATCH_BASE_REF = 'T'    # mismatch in genome


def _plus_genome_basic() -> Dict[str, str]:
    """Genome: _MISMATCH_LEN T's then _ANCHOR_LEN A's at _REF_START."""
    ref_region = _MISMATCH_BASE_REF * _MISMATCH_LEN + _ANCHOR_BASE * _ANCHOR_LEN
    return _build_genome(ref_region, start=_REF_START)


def _plus_seq_basic() -> str:
    return _MISMATCH_BASE_READ * _MISMATCH_LEN + _ANCHOR_BASE * _ANCHOR_LEN


# ---------------------------------------------------------------------------
# Test cases — + strand
# ---------------------------------------------------------------------------

class TestPlusStrand:

    def test_basic_mismatch_cluster(self):
        """Basic + strand: 5-base mismatch cluster in 25M, anchor at offset 5."""
        genome = _plus_genome_basic()
        seq = _plus_seq_basic()
        read = _make_read(_REF_START, [(0, 25)], seq, '+')  # 25M
        _assert_reanchor_equivalent(read, genome, clip_len=_MISMATCH_LEN, label='plus_basic')

    def test_anchor_at_op_boundary(self):
        """+ strand: mismatch cluster spans entire first M op, anchor starts at next op."""
        # CIGAR: 5M (all mismatches) + 20M (all matches)
        # Anchor starts at op boundary: clip_len = 5, run_start_off = 0 in second op.
        genome = _plus_genome_basic()
        seq = _plus_seq_basic()
        read = _make_read(_REF_START, [(0, 5), (0, 20)], seq, '+')  # 5M 20M
        _assert_reanchor_equivalent(read, genome, clip_len=_MISMATCH_LEN, label='plus_op_boundary')

    def test_with_leading_h(self):
        """+ strand: hard-clip at start (H op skipped by both implementations)."""
        # CIGAR: 3H 5X 20=  — H at front, then mismatch cluster, then match run.
        # Hard clip not in query_sequence; query = 25 bases (5 mismatch + 20 match).
        genome = _plus_genome_basic()
        seq = _plus_seq_basic()  # 25 bases
        read = _make_read(_REF_START, [(5, 3), (8, 5), (7, 20)], seq, '+')  # 3H 5X 20=
        _assert_reanchor_equivalent(read, genome, clip_len=_MISMATCH_LEN, label='plus_leading_h')

    def test_with_leading_s(self):
        """+ strand: pre-existing S at start; mismatch cluster follows."""
        # CIGAR: 3S 5M(mismatch) 20M(match)
        # Seq: 3 arbitrary + 5 mismatch + 20 matching = 28 bases.
        # Clip for S + mismatches: run_start_q = 8.
        genome = _build_genome(
            'GGG' + _MISMATCH_BASE_REF * 5 + _ANCHOR_BASE * 20,
            start=_REF_START - 3,
        )
        # The 3S bases don't consume reference; reference_start = 1000 as before.
        seq = 'GGG' + _MISMATCH_BASE_READ * 5 + _ANCHOR_BASE * 20
        # CIGAR: 3S at front (query-only), 25M for the rest
        read = _make_read(_REF_START, [(4, 3), (0, 25)], seq, '+')  # 3S 25M
        _assert_reanchor_equivalent(read, genome, clip_len=3 + _MISMATCH_LEN, label='plus_leading_s')

    def test_with_deletion_before_anchor(self):
        """+ strand: D op between leading mismatches and anchor; r_acc must include D."""
        # CIGAR: 3X 2D 20=
        # Read: 3 mismatch bases + 20 matching bases (D consumes ref only).
        # Genome at ref: 3 mismatches, 2 deleted-from-read (any non-A), 20 A's.
        ref_seq = _MISMATCH_BASE_REF * 3 + 'GG' + _ANCHOR_BASE * 20
        genome = _build_genome(ref_seq, start=_REF_START)
        seq = _MISMATCH_BASE_READ * 3 + _ANCHOR_BASE * 20  # 23 bases
        read = _make_read(_REF_START, [(8, 3), (2, 2), (7, 20)], seq, '+')  # 3X 2D 20=
        _assert_reanchor_equivalent(read, genome, clip_len=3, label='plus_del_before_anchor')

    def test_anchor_mid_m_op(self):
        """+ strand: mismatch cluster mid-M op, anchor starts at offset k > 0."""
        # Single 30M op: first 7 bases are mismatch, next 23 are match.
        # Anchor starts at offset 7 within the M op.
        ref_seq = _MISMATCH_BASE_REF * 7 + _ANCHOR_BASE * 23
        genome = _build_genome(ref_seq, start=_REF_START)
        seq = _MISMATCH_BASE_READ * 7 + _ANCHOR_BASE * 23
        read = _make_read(_REF_START, [(0, 30)], seq, '+')  # 30M
        _assert_reanchor_equivalent(read, genome, clip_len=7, label='plus_mid_m_op')


# ---------------------------------------------------------------------------
# Test cases — - strand
# ---------------------------------------------------------------------------
# For - strand reads, the 5' RNA end is at the BAM-right (high coordinate).
# Mismatch cluster at the BAM-right end → trailing S in new CIGAR.
# reference_start is unchanged.
#
# Layout: read maps [_REF_START, _REF_START + 25).
#   Low-coord end (BAM left): 20 A's that match the genome (body).
#   High-coord end (BAM right): 5 mismatched bases (= RNA 5' cluster).
# Genome: 20 A's then 5 T's.
# Read seq (BAM order, low→high coord): 20 A's + 5 C's.
# CIGAR: 25M.
# Expected reanchor: 20M 5S, reference_start unchanged.
# ---------------------------------------------------------------------------

_MINUS_REF_START = 2000


def _minus_genome_basic() -> Dict[str, str]:
    """Genome: _ANCHOR_LEN A's then _MISMATCH_LEN T's at _MINUS_REF_START."""
    ref_region = _ANCHOR_BASE * _ANCHOR_LEN + _MISMATCH_BASE_REF * _MISMATCH_LEN
    return _build_genome(ref_region, start=_MINUS_REF_START)


def _minus_seq_basic() -> str:
    return _ANCHOR_BASE * _ANCHOR_LEN + _MISMATCH_BASE_READ * _MISMATCH_LEN


class TestMinusStrand:

    def test_basic_mismatch_cluster(self):
        """Basic - strand: 5-base mismatch cluster at BAM-right end of 25M."""
        genome = _minus_genome_basic()
        seq = _minus_seq_basic()
        read = _make_read(_MINUS_REF_START, [(0, 25)], seq, '-')
        _assert_reanchor_equivalent(read, genome, clip_len=_MISMATCH_LEN, label='minus_basic')

    def test_anchor_at_op_boundary(self):
        """- strand: cluster spans last M op entirely; anchor at op boundary."""
        genome = _minus_genome_basic()
        seq = _minus_seq_basic()
        read = _make_read(_MINUS_REF_START, [(0, 20), (0, 5)], seq, '-')  # 20M 5M
        _assert_reanchor_equivalent(read, genome, clip_len=_MISMATCH_LEN, label='minus_op_boundary')

    def test_with_trailing_s(self):
        """- strand: pre-existing trailing S; mismatch cluster before it."""
        # CIGAR: 20M(match) 5M(mismatch) 3S
        # Seq: 20 match + 5 mismatch + 3 arbitrary = 28 bases.
        # Genome: 20 A's + 5 T's (mismatches) at _MINUS_REF_START.
        # The 3S doesn't consume reference; reference_start = _MINUS_REF_START.
        genome = _minus_genome_basic()
        seq = _ANCHOR_BASE * _ANCHOR_LEN + _MISMATCH_BASE_READ * _MISMATCH_LEN + 'GGG'
        read = _make_read(_MINUS_REF_START, [(0, 25), (4, 3)], seq, '-')  # 25M 3S
        _assert_reanchor_equivalent(read, genome, clip_len=_MISMATCH_LEN + 3, label='minus_trailing_s')

    def test_anchor_mid_m_op(self):
        """- strand: cluster mid-M op at BAM-right; anchor at offset from end."""
        # 30M: last 7 bases (BAM-right) are mismatch, first 23 are anchor (match).
        ref_seq = _ANCHOR_BASE * 23 + _MISMATCH_BASE_REF * 7
        genome = _build_genome(ref_seq, start=_MINUS_REF_START)
        seq = _ANCHOR_BASE * 23 + _MISMATCH_BASE_READ * 7
        read = _make_read(_MINUS_REF_START, [(0, 30)], seq, '-')
        _assert_reanchor_equivalent(read, genome, clip_len=7, label='minus_mid_m_op')

    def test_with_deletion_before_cluster(self):
        """- strand: D op between anchor and cluster (D consumes ref, not query)."""
        # CIGAR: 20= 2D 5X  (BAM order: match, deletion, mismatch at right end)
        # Genome: 20 A's + 2 G's (deleted from read) + 5 T's.
        # Read (BAM seq): 20 A's + 5 C's = 25 bases (D not in seq).
        ref_seq = _ANCHOR_BASE * 20 + 'GG' + _MISMATCH_BASE_REF * 5
        genome = _build_genome(ref_seq, start=_MINUS_REF_START)
        seq = _ANCHOR_BASE * 20 + _MISMATCH_BASE_READ * 5  # 25 bases
        read = _make_read(_MINUS_REF_START, [(7, 20), (2, 2), (8, 5)], seq, '-')  # 20= 2D 5X
        _assert_reanchor_equivalent(read, genome, clip_len=5, label='minus_del_before_cluster')


# ---------------------------------------------------------------------------
# No-op guard: _apply_reanchor_from_clip_len with degenerate inputs.
# ---------------------------------------------------------------------------

class TestDegenerate:

    def test_clip_len_zero(self):
        genome = _plus_genome_basic()
        seq = _plus_seq_basic()
        read = _make_read(_REF_START, [(0, 25)], seq, '+')
        assert _apply_reanchor_from_clip_len(read, 0) is False

    def test_clip_len_gte_seq_len(self):
        genome = _plus_genome_basic()
        seq = _plus_seq_basic()
        read = _make_read(_REF_START, [(0, 25)], seq, '+')
        assert _apply_reanchor_from_clip_len(read, 25) is False
        assert _apply_reanchor_from_clip_len(read, 30) is False

    def test_unmapped_read(self):
        hdr = pysam.AlignmentHeader.from_dict({
            'HD': {'VN': '1.6'},
            'SQ': [{'SN': _CHROM, 'LN': _CHROM_LEN}],
        })
        r = pysam.AlignedSegment(hdr)
        r.is_unmapped = True
        r.query_sequence = 'ACGT'
        assert _apply_reanchor_from_clip_len(r, 2) is False
