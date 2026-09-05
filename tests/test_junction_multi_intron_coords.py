"""Multi-intron N-op coordinate regression tests (ISSUE-004).

The genomic cursor of every position-tracking CIGAR walker must advance across
an ``N`` op — an intron consumes reference.  Before this guard existed the
walkers used ``_REF_CONSUMING`` (M/D/=/X only), so on a read with more than one
intron every junction after the first was reported short by the summed length of
all preceding introns.  Human DRS exposed it (81/81 multi-intron panel reads
wrong, 85 % of the observed junction pool phantom); yeast never did, because
~95 % of intron-bearing yeast genes have a single intron.

These tests derive the expected coordinates independently from the CIGAR (SAM
spec: M/D/N/=/X consume reference) and from ``pysam.AlignedSegment.get_blocks``,
so they cannot be satisfied by the same off-by-intron arithmetic they guard.

Author: Kevin R. Roy (agent S1)
"""

import sys
from pathlib import Path

import pysam
import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.utils.genome import standardize_chrom_name  # noqa: E402
from rectify.core.splice.junction_refiner import _iter_n_ops, _n_op_ref_start  # noqa: E402
from rectify.core.splice.junction_scoring import (  # noqa: E402
    _N,
    _REF_CONSUMING,
    _REF_CONSUMING_POS,
    _collect_junction_counts_core,
    collect_junctions_from_bam,
)

CHROM = "chr5"


def _pool_chrom():
    """The pool builders key junctions by ``standardize_chrom_name``.

    These tests are about COORDINATES, so they compare against the standardized
    name rather than the raw one — chromosome-name policy is a separate defect
    (ISSUE-001).  Resolved at CALL time, not import time: the contig registry is
    process-global module state that another test module may have populated.
    """
    return standardize_chrom_name(CHROM)
CHROM_LEN = 200_000
REF_START = 1_000

# 4 exons / 3 introns, deliberately unequal intron lengths so an off-by-intron
# walker cannot land on a right answer by accident.
MULTI_CIGAR = [(4, 12), (0, 40), (3, 771), (0, 30), (3, 11_964),
               (0, 25), (3, 1_389), (0, 35), (4, 8)]
SINGLE_CIGAR = [(0, 50), (3, 300), (0, 50)]
# Indels on both sides of an intron: D advances the cursor, I does not.
INDEL_CIGAR = [(0, 20), (2, 3), (0, 10), (3, 500), (0, 15), (1, 4), (0, 20),
               (3, 800), (0, 30)]

_SAM_REF_CONSUMING = {0, 2, 3, 7, 8}   # M D N = X — straight from the SAM spec


def _expected_n_ops(ref_start, cigar):
    """(cigar_idx, intron_start, intron_end) per the SAM spec, computed here."""
    pos = ref_start
    out = []
    for i, (op, length) in enumerate(cigar):
        if op == _N:
            out.append((i, pos, pos + length))
        if op in _SAM_REF_CONSUMING:
            pos += length
    return out


def _query_len(cigar):
    return sum(l for op, l in cigar if op in (0, 1, 4, 7, 8))


def _make_read(cigar, name="r", ref_start=REF_START):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"},
         "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
    )
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.reference_id = 0
    read.reference_start = ref_start
    read.mapping_quality = 60
    read.cigartuples = cigar
    n = _query_len(cigar)
    read.query_sequence = "ACGT" * (n // 4) + "ACGT"[: n % 4]
    read.query_qualities = pysam.qualitystring_to_array("I" * n)
    return read


def _write_bam(path, reads):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"},
         "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
    )
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.index(str(path))


# ---------------------------------------------------------------------------
# The op sets themselves
# ---------------------------------------------------------------------------

def test_ref_consuming_pos_includes_n_and_ref_consuming_does_not():
    """The two sets are NOT interchangeable — do not 'unify' them.

    ``_REF_CONSUMING`` is the EXON-flank set used by
    ``_apply_junction_replacement`` to find and consume the aligned ops abutting
    an N-op; widening it to include N would let those walkers step through an
    intron.  ``_REF_CONSUMING_POS`` is the cursor set.
    """
    assert _N not in _REF_CONSUMING
    assert _N in _REF_CONSUMING_POS
    assert _REF_CONSUMING_POS == _REF_CONSUMING | {_N}
    assert _REF_CONSUMING_POS == frozenset(_SAM_REF_CONSUMING)


# ---------------------------------------------------------------------------
# _iter_n_ops / _n_op_ref_start
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cigar", [MULTI_CIGAR, SINGLE_CIGAR, INDEL_CIGAR])
def test_iter_n_ops_matches_sam_spec_walk(cigar):
    read = _make_read(cigar)
    expected = _expected_n_ops(REF_START, cigar)
    got = [(i, s, e) for i, s, e, _q in _iter_n_ops(read)]
    assert got == expected


def test_iter_n_ops_matches_pysam_blocks():
    """Cross-check against pysam's own block decomposition (independent code)."""
    read = _make_read(MULTI_CIGAR)
    blocks = read.get_blocks()          # aligned blocks, gaps = N (and D)
    introns = [(a_end, b_start) for (_, a_end), (b_start, _) in zip(blocks, blocks[1:])
               if b_start > a_end]
    got = [(s, e) for _i, s, e, _q in _iter_n_ops(read)]
    assert got == introns
    assert len(got) == 3


def test_iter_n_ops_third_intron_is_not_short_by_preceding_introns():
    """The exact ISSUE-004 signature, stated as an inequality."""
    read = _make_read(MULTI_CIGAR)
    got = [(s, e) for _i, s, e, _q in _iter_n_ops(read)]
    buggy = []
    pos = REF_START
    for op, length in MULTI_CIGAR:
        if op == _N:
            buggy.append((pos, pos + length))
        if op in _REF_CONSUMING:        # the old (wrong) cursor set
            pos += length
    assert got[0] == buggy[0]           # first intron was never affected
    assert got[1] != buggy[1]
    assert got[2] != buggy[2]
    assert got[1][0] - buggy[1][0] == 771            # 1st intron length
    assert got[2][0] - buggy[2][0] == 771 + 11_964   # 1st + 2nd


def test_n_op_ref_start_matches_iter_n_ops():
    read = _make_read(MULTI_CIGAR)
    for idx, start, _end, _q in _iter_n_ops(read):
        assert _n_op_ref_start(read, idx) == start


def test_iter_n_ops_q_split_unchanged_by_the_cursor_fix():
    """q_split counts QUERY bases; introns consume none, so it must not move."""
    read = _make_read(MULTI_CIGAR)
    q_splits = [q for _i, _s, _e, q in _iter_n_ops(read)]
    # 12S leading (not counted) + 40M, then +30M, then +25M
    assert q_splits == [40, 70, 95]


def test_iter_n_ops_leading_softclip_still_excluded_from_q_split():
    read = _make_read([(4, 9), (0, 20), (3, 100), (0, 20)])
    assert [q for _i, _s, _e, q in _iter_n_ops(read)] == [20]


# ---------------------------------------------------------------------------
# Pool builders
# ---------------------------------------------------------------------------

def test_pool_builders_report_only_real_coordinates(tmp_path):
    bam = tmp_path / "multi_intron.bam"
    reads = [_make_read(MULTI_CIGAR, name="multi"),
             _make_read(SINGLE_CIGAR, name="single", ref_start=150_000)]
    _write_bam(bam, reads)

    truth = {(_pool_chrom(), s, e) for r in reads
             for _i, s, e in _expected_n_ops(r.reference_start, r.cigartuples)}
    assert len(truth) == 4

    pool = collect_junctions_from_bam(str(bam))
    assert pool == truth

    anchor, raw = _collect_junction_counts_core(str(bam), min_anchor_overhang=0)
    assert set(raw) == truth
    assert set(anchor) <= truth


def test_pool_builder_single_intron_read_unchanged(tmp_path):
    """A single-intron read was always right — prove the fix did not move it."""
    bam = tmp_path / "single_intron.bam"
    _write_bam(bam, [_make_read(SINGLE_CIGAR, name="single")])
    assert collect_junctions_from_bam(str(bam)) == {
        (_pool_chrom(), REF_START + 50, REF_START + 350)
    }
