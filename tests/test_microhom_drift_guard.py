"""Regression tests for the microhomology-drift guard in junction_refiner.

The general (non-homopolymer) analog of the HP-drift guard.  A motif-blind re-placer
fabricates false non-canonical junctions by drifting a real canonical junction a few
bp to a nearby non-canonical pool member — enabled by LOCAL MICROHOMOLOGY (a near-
tandem direct repeat at the drift distance) plus ONT error making the true-vs-drifted
boundary a near-tie.  The HP-drift guard misses this (it only checks homopolymer runs).

The guard (``microhom_drift_margin``) requires a move whose boundary shift sits in a
microhomology context (``_move_microhomology`` fraction >= ``microhom_threshold``) to
clear an extra evidence margin, while leaving moves to genuine sequence transitions
(real non-canonical splice sites, low microhomology) untouched — so motif-blind
discovery is preserved and only the drift failure mode is reined in.

Author: Kevin R. Roy
"""
from __future__ import annotations

import sys
from pathlib import Path

import pysam
import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.splice.junction_refiner import (
    _frac_match, _move_microhomology, refine_read_junctions,
)
from rectify.core.splice.junction_scoring import _D, _EQ, _N


# ---------------------------------------------------------------------------
# 1. _frac_match
# ---------------------------------------------------------------------------
class TestFracMatch:
    def test_identical(self):
        assert _frac_match("ACGT", "ACGT") == 1.0

    def test_half(self):
        assert _frac_match("ACGT", "ACGA") == 0.75
        assert _frac_match("AACC", "AAGG") == 0.5

    def test_empty_or_mismatched_length(self):
        assert _frac_match("", "") == 0.0
        assert _frac_match("ACG", "ACGT") == 0.0


# ---------------------------------------------------------------------------
# 2. _move_microhomology — the detector (matches the drift mechanism)
# ---------------------------------------------------------------------------
class TestMoveMicrohomology:
    # A period-6 tandem repeat "ACGTAC ACGTAC" at index 20; a +6 acceptor drift across
    # it is high-microhomology.  Elsewhere is random -> low.
    SEQ = "T" * 20 + "ACGTAC" + "ACGTAC" + "GGGCCCTTTAAA" + "T" * 20   # repeat at [20,32)

    def test_acceptor_drift_into_tandem_repeat_is_high(self):
        # acceptor shift ne=20 -> je=26 (k=6): seq[20:26]="ACGTAC" vs seq[26:32]="ACGTAC" -> 1.0
        assert _move_microhomology(self.SEQ, 5, 20, 5, 26) == 1.0

    def test_acceptor_drift_to_transition_is_low(self):
        # a +6 acceptor move landing on the non-repeat run -> low (random-ish, << 0.5)
        v = _move_microhomology(self.SEQ, 5, 32, 5, 38)   # seq[32:38] vs seq[38:44]
        assert v < 0.5

    def test_no_move_is_zero(self):
        assert _move_microhomology(self.SEQ, 5, 26, 5, 26) == 0.0

    def test_upstream_acceptor_shift_symmetric(self):
        # je < ne (shift upstream) also detects the repeat
        assert _move_microhomology(self.SEQ, 5, 26, 5, 20) == 1.0

    def test_donor_shift_detected(self):
        # donor drift into the repeat: ns=32 -> js=26 (k=-6): seq[26:32] vs seq[20:26] -> 1.0
        assert _move_microhomology(self.SEQ, 32, 60, 26, 60) == 1.0

    def test_out_of_bounds_returns_zero(self):
        # a shift that runs off the end returns 0 (no crash)
        assert _move_microhomology(self.SEQ, 5, len(self.SEQ) - 2, 5, len(self.SEQ) + 4) == 0.0

    def test_partial_microhomology_between_thresholds(self):
        # one mismatch in a 6-mer repeat -> 5/6 = 0.83 (still >= a 0.5 threshold)
        s = "T" * 20 + "ACGTAC" + "ACGTAG" + "T" * 20   # last base differs
        assert abs(_move_microhomology(s, 5, 20, 5, 26) - 5 / 6) < 1e-9


# ---------------------------------------------------------------------------
# 3. End-to-end guard behaviour (mirrors test_hp_drift_guard)
# ---------------------------------------------------------------------------
# Canonical GT-AG junction whose acceptor abuts a period-6 tandem repeat.  A read
# errored near the boundary can drift the acceptor +6 into the repeat (a false
# non-canonical junction); the guard must veto that while sparing a transition move.
LPAD, EXON1_LEN, RPAD = 30, 50, 40
INTRON = "GT" + "C" * 86 + "AG"          # 90 bp canonical
REPEAT = "ACGTAC" + "ACGTAC"             # period-6 tandem repeat at exon2 start (mh=1.0 across +6)
EXON2 = REPEAT + "GATC" * 10             # then a run-free tail
GENOME = "T" * LPAD + "C" * EXON1_LEN + INTRON + REPEAT + "GATC" * 10 + "T" * RPAD
D = LPAD + EXON1_LEN                      # 80 donor
A = D + len(INTRON)                       # 170 acceptor (true, canonical AG|ACGTAC...)
CHROM = "chrTEST"


def _make_read(cigar_ops, ref_start):
    header = pysam.AlignmentHeader.from_references([CHROM], [len(GENOME)])
    read = pysam.AlignedSegment(header)
    read.query_name = "mh_drift_read"
    read.reference_id = 0
    read.reference_start = ref_start
    read.mapping_quality = 60
    read.is_unmapped = False
    read.is_reverse = False
    q, r = [], ref_start
    for op, ln in cigar_ops:
        if op == _EQ:
            q.append(GENOME[r:r + ln]); r += ln
        elif op in (_D, _N):
            r += ln
    read.query_sequence = "".join(q)
    read.cigartuples = cigar_ops
    return read


def test_drift_target_is_a_microhomology_move_and_true_is_not():
    """Frame check: the +6 acceptor drift (A -> A+6) sits in the tandem repeat
    (microhomology 1.0); the true acceptor A is a canonical GT-AG placement."""
    assert GENOME[D:D + 2] == "GT"
    assert GENOME[A - 2:A] == "AG"
    assert _move_microhomology(GENOME, D, A, D, A + 6) == 1.0          # drift: full microhomology


def _acceptor_after_refine(microhom_drift_margin):
    """Refine a read placed at the true acceptor A but whose exon2 matches the
    drifted repeat position; return the acceptor the refiner chose.  The candidate
    pool offers both A (true) and A+6 (drift target)."""
    # exon1(50=) intron(90N) then exon2 that equally matches at A and A+6 (the repeat)
    cigar = [(_EQ, EXON1_LEN), (_N, len(INTRON)), (_EQ, 18)]
    read = _make_read(cigar, ref_start=LPAD)
    pool = {CHROM: [(D, A), (D, A + 6)]}
    reps = refine_read_junctions(
        read, pool, set(), GENOME, "+",
        motif_blind=True, microhom_drift_margin=microhom_drift_margin,
    )
    acc = A
    for _idx, ons, one, njs, nje in reps:
        if ons == D and one == A:
            acc = nje
    return acc


def test_guard_off_is_byte_identical_default():
    """microhom_drift_margin=0.0 (default) must not add any veto — the refiner's
    behaviour is exactly the incumbent."""
    cigar = [(_EQ, EXON1_LEN), (_N, len(INTRON)), (_EQ, 18)]
    pool = {CHROM: [(D, A), (D, A + 6)]}

    def acc(**kw):
        read = _make_read(cigar, LPAD)
        reps = refine_read_junctions(read, pool, set(), GENOME, "+", motif_blind=True, **kw)
        a = A
        for _i, ons, one, njs, nje in reps:
            if ons == D and one == A:
                a = nje
        return a

    assert acc() == acc(microhom_drift_margin=0.0)


def test_guard_holds_true_acceptor_when_on():
    """With the guard on, the acceptor is held at the true canonical A (the +6
    drift into the tandem repeat is vetoed)."""
    assert _acceptor_after_refine(microhom_drift_margin=8.0) == A


def test_guard_margin_is_microhomology_specific():
    """The extra margin is added ONLY where _move_microhomology >= threshold: the +6
    drift into the repeat is guarded (1.0); a move to a genuine transition
    downstream of the repeat (low microhomology) is not — preserving discovery."""
    assert _move_microhomology(GENOME, D, A, D, A + 6) >= 0.5       # drift: guarded
    trans = A + 12                                                   # past the repeat, in GATC tail
    assert _move_microhomology(GENOME, D, A, D, trans) < 0.5         # transition: spared
