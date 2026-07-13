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
    _frac_match, _hp_run_across, _move_microhomology, _effective_veto_margin,
    refine_read_junctions,
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

    def test_ambiguity_bases_do_not_match(self):
        # audit A5: N==N (and any non-ACGT, incl. lowercase) must NOT count as agreement,
        # else a genome ambiguity run scores phantom microhomology 1.0 and falsely vetoes.
        assert _frac_match("NNNN", "NNNN") == 0.0
        assert _frac_match("nnnn", "nnnn") == 0.0
        assert _frac_match("XXXX", "XXXX") == 0.0
        assert _frac_match("acgt", "acgt") == 0.0          # lowercase are not real bases here
        assert _frac_match("NANA", "NANA") == 0.5          # only the two real A's agree


# ---------------------------------------------------------------------------
# 1b. _hp_run_across — audit A5: an N-run is NOT a homopolymer
# ---------------------------------------------------------------------------
class TestHpRunAcrossAmbiguity:
    def test_n_run_is_not_a_homopolymer(self):
        assert _hp_run_across("NNNNN", 2, 4) == 0            # A5: was 5 (false HP run)
        assert _hp_run_across("XXXXXX", 3, 4) == 0

    def test_real_homopolymer_still_detected(self):
        assert _hp_run_across("AAAAA", 2, 4) == 5            # real A-run unaffected
        assert _hp_run_across("ACGTC", 2, 4) == 0            # transition, no run


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

    def test_both_boundary_transition_not_masked_by_other_repeat(self):
        # audit A8: a both-boundary move whose DONOR abuts a (CAG)n repeat (frac 1.0) but
        # whose ACCEPTOR is a genuine transition (frac 0.0) must be SPARED (min over moved
        # boundaries = 0.0), not vetoed by max() masking the transition with the repeat.
        seq = "CAG" + "CAG" + "GTAAGTACTAAC" + "GCC" + "TTA" + "CTGCTG"   # no N, pure ACGT
        # donor shift ns=3->js=6 (CAG|CAG, frac 1.0); acceptor shift ne=18->je=21 (GCC|TTA, 0.0)
        assert _move_microhomology(seq, 3, 18, 6, 21) == 0.0            # min → spared (A8 fixed)
        # sanity: the donor leg alone IS a full repeat, the acceptor leg alone IS a transition
        assert _move_microhomology(seq, 3, 18, 6, 18) == 1.0           # donor-only move: repeat
        assert _move_microhomology(seq, 3, 18, 3, 21) == 0.0           # acceptor-only move: transition


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


# ---------------------------------------------------------------------------
# 4. _effective_veto_margin — the read-evidence near-tie cap on the read-blind
#    drift margins (bounds the discovery-loss the audit flagged; MICROHOM_AUDIT_SYNTHESIS)
# ---------------------------------------------------------------------------
class TestEffectiveVetoMargin:
    def test_cap_disabled_is_identity(self):
        # cap <= 0.0 → eff_margin returned unchanged (byte-identical to pre-cap veto)
        assert _effective_veto_margin(hold_margin=0.0, eff_margin=8.0, drift_near_tie_cap=0.0) == 8.0
        assert _effective_veto_margin(hold_margin=1.0, eff_margin=9.0, drift_near_tie_cap=-1.0) == 9.0

    def test_no_drift_added_cap_is_noop(self):
        # eff_margin == hold_margin (no hp/microhom drift flagged) → cap irrelevant,
        # hold_margin (read-agnostic) is never capped even with a tiny cap
        assert _effective_veto_margin(hold_margin=2.0, eff_margin=2.0, drift_near_tie_cap=1.0) == 2.0

    def test_cap_bounds_the_drift_margin(self):
        # hold < cap < eff → veto margin is capped at `cap`
        assert _effective_veto_margin(hold_margin=0.0, eff_margin=8.0, drift_near_tie_cap=2.0) == 2.0
        assert _effective_veto_margin(hold_margin=1.0, eff_margin=8.0, drift_near_tie_cap=3.0) == 3.0

    def test_cap_above_eff_is_inactive(self):
        # cap >= eff_margin → the cap never binds, eff_margin returned unchanged
        assert _effective_veto_margin(hold_margin=0.0, eff_margin=3.0, drift_near_tie_cap=8.0) == 3.0

    def test_hold_margin_is_never_capped(self):
        # THE interaction the audit flags: a cap BELOW hold_margin must NOT drop the
        # veto below hold_margin — hold (read-agnostic blunt prior) is a floor.
        assert _effective_veto_margin(hold_margin=4.0, eff_margin=8.0, drift_near_tie_cap=2.0) == 4.0
        # formula check across the regime: max(hold, min(eff, cap))
        assert _effective_veto_margin(hold_margin=4.0, eff_margin=8.0, drift_near_tie_cap=6.0) == 6.0


def test_near_tie_cap_byte_identical_when_disabled():
    """drift_near_tie_cap=0.0 (default) must not change the refiner's output vs omitting
    it — the cap is inert at default (independent of any drift margin setting)."""
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

    # with the microhom guard on, cap=0.0 == omitting cap == same veto
    assert acc(microhom_drift_margin=8.0) == acc(microhom_drift_margin=8.0, drift_near_tie_cap=0.0)
