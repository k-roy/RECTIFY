"""Regression tests for the HP-drift guard in junction_refiner.

The guard controls the ubiquitous-ONT-homopolymer-undercall failure mode: a
purely evidence-driven re-placer can slide a junction boundary *into* a
homopolymer run (re-partitioning identical bases between intron and exon),
fabricating a false non-canonical junction.  The guard (``hp_drift_margin``)
requires a move that lands a boundary INSIDE a homopolymer run to clear an extra
evidence margin, while leaving moves to genuine sequence transitions (real
non-canonical acceptors) untouched — so motif-blind discovery is preserved and
only the homopolymer failure mode is reined in.

Two layers:
  1. ``_hp_run_across`` — the detector (boundary inside a run vs a transition).
  2. ``refine_read_junctions`` end-to-end — the guard flips a real drift decision
     (holds when hp_drift_margin > 0) without changing anything at margin 0
     (byte-identical) and without blocking a move to a transition site.

Author: Kevin R. Roy
"""
from __future__ import annotations

import sys
from pathlib import Path

import pysam
import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.splice.junction_refiner import _hp_run_across, refine_read_junctions
from rectify.core.splice.junction_scoring import _D, _EQ, _N


# ---------------------------------------------------------------------------
# 1. Detector: _hp_run_across
# ---------------------------------------------------------------------------
class TestHpRunAcross:
    # ...CG A A A A A A A A GT...  a clean 8-A run at [3, 11)
    SEQ = "CG" + "A" * 8 + "GT" + "C" * 5   # A-run = indices [2, 10)

    def test_inside_run_returns_length(self):
        # boundaries strictly inside the run: seq[pos-1]==seq[pos]=='A'
        for pos in (4, 6, 9):
            assert _hp_run_across(self.SEQ, pos, 4) == 8, pos

    def test_run_start_is_a_transition(self):
        # pos=2: seq[1]='G', seq[2]='A' -> transition (the run is entirely exonic)
        assert _hp_run_across(self.SEQ, 2, 4) == 0

    def test_run_end_is_a_transition(self):
        # pos=10: seq[9]='A', seq[10]='G' -> transition
        assert _hp_run_across(self.SEQ, 10, 4) == 0

    def test_short_run_below_min_returns_zero(self):
        s = "CG" + "A" * 3 + "GT"          # only a 3-A run
        assert _hp_run_across(s, 3, 4) == 0   # inside it, but 3 < min_run 4
        assert _hp_run_across(s, 3, 3) == 3   # min_run 3 -> now it counts

    def test_out_of_bounds_returns_zero(self):
        assert _hp_run_across(self.SEQ, 0, 4) == 0
        assert _hp_run_across(self.SEQ, len(self.SEQ), 4) == 0

    def test_non_A_homopolymer_also_detected(self):
        s = "AT" + "G" * 6 + "AT"          # a 6-G run at [2, 8)
        assert _hp_run_across(s, 5, 4) == 6
        assert _hp_run_across(s, 2, 4) == 0   # G-run start is a transition


# ---------------------------------------------------------------------------
# 2. End-to-end: the guard flips a drift decision in refine_read_junctions
# ---------------------------------------------------------------------------
# Canonical GT-AG junction with an 8-A run IMMEDIATELY after the acceptor.  A read
# undercalled by 1 A in the run is placed at the true junction (with a 1-base
# deletion in the run); the intron-grown placement (acceptor +1, INTO the run)
# explains the read with NO deletion, so it scores strictly better and an
# unguarded re-placer drifts there.  The +1 acceptor lands inside the A-run, so
# the guard must veto it.
#
#   [lpad 30] [exon1 50] GT [intron body 86] AG | [A x8] G [exon2 tail 30] [rpad]
#            ^30         ^80=D              ^168 ^170=A (exon2 start)
LPAD, EXON1_LEN, RPAD = 30, 50, 20
INTRON = "GT" + "C" * 86 + "AG"            # 90 bp, canonical; acceptor dinuc "AG"
EXON2 = "A" * 8 + "GATC" * 10              # A-run [0,8) of exon2, then a run-free tail
GENOME = "T" * LPAD + "C" * EXON1_LEN + INTRON + EXON2 + "T" * RPAD
D = LPAD + EXON1_LEN                        # 80  donor / intron start
A = D + len(INTRON)                         # 170 acceptor / intron end == exon2 start
CHROM = "chrTEST"


def _make_read(cigar_ops, ref_start):
    header = pysam.AlignmentHeader.from_references([CHROM], [len(GENOME)])
    read = pysam.AlignedSegment(header)
    read.query_name = "hp_drift_read"
    read.reference_id = 0
    read.reference_start = ref_start
    read.mapping_quality = 60
    read.is_unmapped = False
    read.is_reverse = False
    # build a query consistent with the CIGAR (copy ref bases for =, skip for D/N)
    q, r = [], ref_start
    for op, ln in cigar_ops:
        if op == _EQ:
            q.append(GENOME[r:r + ln]); r += ln
        elif op in (_D, _N):
            r += ln
    read.query_sequence = "".join(q)
    read.cigartuples = cigar_ops
    return read


def _acceptors_after_refine(hp_drift_margin):
    """Refine the undercalled read; return the acceptor the refiner chose."""
    # true placement: exon1(50=) intron(90N) 7 run-A's (7=) 1 deleted A (1D) tail(30=)
    cigar = [(_EQ, EXON1_LEN), (_N, len(INTRON)), (_EQ, 7), (_D, 1), (_EQ, 30)]
    read = _make_read(cigar, ref_start=LPAD)
    pool = {CHROM: [(D, A), (D, A + 1)]}   # true + the into-run drift candidate
    reps = refine_read_junctions(
        read, pool, set(), GENOME, "+",
        motif_blind=True, hp_drift_margin=hp_drift_margin,
    )
    # replacements are (cigar_idx, old_ns, old_ne, new_js, new_je)
    new_acc = A
    for _idx, ons, one, njs, nje in reps:
        if ons == D and one == A:
            new_acc = nje
    return new_acc


def test_true_acceptor_is_a_transition_and_drift_is_inside_run():
    """Frame check: the true acceptor is a transition; +1 lands inside the run."""
    assert GENOME[D:D + 2] == "GT"
    assert GENOME[A - 2:A] == "AG"
    assert _hp_run_across(GENOME, A, 4) == 0          # true acceptor: transition
    assert _hp_run_across(GENOME, A + 1, 4) == 8       # drift +1: inside the 8-A run


def test_unguarded_replacer_drifts_into_the_homopolymer():
    """With no guard (margin 0), the better-scoring intron-grown placement wins."""
    assert _acceptors_after_refine(hp_drift_margin=0.0) == A + 1


def test_guard_vetoes_the_drift():
    """With the guard, the into-run drift is held at the true acceptor."""
    assert _acceptors_after_refine(hp_drift_margin=3.0) == A


def test_guard_margin_is_hp_specific():
    """The motif-blind-discovery half of the design, as a contract: the emit gate
    adds the extra margin ONLY where ``_hp_run_across > 0``.  The true canonical
    acceptor and every downstream transition are 0 (no margin → a move there is
    judged on evidence alone, exactly like arm-B), while only the into-run drift
    carries the guard.  Combined with test_guard_vetoes_the_drift (the into-run
    move IS vetoed) and test_unguarded_replacer_drifts (the move happens without
    the guard), this pins that the guard touches the drift and nothing else."""
    assert _hp_run_across(GENOME, A, 4) == 0          # true acceptor: transition, no margin
    assert _hp_run_across(GENOME, A + 12, 4) == 0      # a downstream transition: no margin
    assert _hp_run_across(GENOME, A + 1, 4) == 8       # the drift: guarded (inside the run)


def test_margin_zero_is_byte_identical_to_no_guard():
    """hp_drift_margin=0.0 must not change any decision (the shipped default off).
    The undercalled read drifts to A+1 both with the param unset and =0.0."""
    # param unset uses the 0.0 default; explicit 0.0 must match.
    cigar = [(_EQ, EXON1_LEN), (_N, len(INTRON)), (_EQ, 7), (_D, 1), (_EQ, 30)]
    pool = {CHROM: [(D, A), (D, A + 1)]}

    def acc(**kw):
        read = _make_read(cigar, ref_start=LPAD)
        reps = refine_read_junctions(read, pool, set(), GENOME, "+", motif_blind=True, **kw)
        a = A
        for _i, ons, one, njs, nje in reps:
            if ons == D and one == A:
                a = nje
        return a

    assert acc() == acc(hp_drift_margin=0.0) == A + 1
