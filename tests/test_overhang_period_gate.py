"""Tandem-period self-match gate (the 2026-08-11 reporter-construct fix).

The I_eff estimators are blind to periodicity with unit > ~3 nt; an MS2
stem-loop array passed the gate at 76-81 bits and exploded resolver cost
>500x on reporter constructs. The period term caps W_max at (period - 1),
generalizing the poly(A)/(AG)n refusals.
"""

import random

from rectify.core.splice.overhang_informativeness import (
    COUNTERS,
    assess_overhang,
    max_search_window_bp,
    min_self_match_period,
    reset_counters,
)

MS2 = "ACATGAGGATCACCCATGT"          # canonical stem-loop unit, 19 nt
RANDOMISH = "ACGTGATCCATGCTTACGCTGACTATCGGACTTCAGATCCGTACTGACGATCCATGCATC"
ACT1_60 = "ATGGATTCTGAGGTTGCTGCTTTGGTTATTGACAACGGTTCTGGTATGTGTAAAGCCGGT"


def test_ms2_array_detected_and_capped():
    a = assess_overhang(MS2 * 3, max_window=5000)
    assert a.period == 19
    assert a.w_max_bp == 18          # was 5000 before the fix
    assert not a.refused             # small boundary search stays allowed


def test_error_tolerant_at_10pct():
    rng = random.Random(644)
    seq = list(MS2 * 3)
    for i in rng.sample(range(len(seq)), int(0.1 * len(seq))):
        seq[i] = rng.choice("ACGT")
    a = assess_overhang("".join(seq), max_window=5000)
    assert a.period == 19
    assert a.w_max_bp <= 18


def test_generic_kmer_repeat_capped():
    a = assess_overhang("GATCCTAG" * 8, max_window=5000)
    assert a.period == 8
    assert a.w_max_bp == 7


def test_homopolymer_and_dinuc_still_refused_via_period():
    # the period term generalizes the old refusals: period 1 and 2 -> W_max 0
    assert min_self_match_period("A" * 60) == 1
    assert min_self_match_period("AG" * 30) == 2
    assert assess_overhang("A" * 60).refused
    assert assess_overhang("AG" * 30).refused


def test_aperiodic_sequences_unaffected():
    for seq in (RANDOMISH, ACT1_60):
        a = assess_overhang(seq, max_window=5000)
        assert a.period is None
        assert a.w_max_bp == 5000
        assert not a.refused


def test_short_sequences_safe():
    # below min_overlap no shift can qualify; nothing raises
    for seq in ("", "ACGT", "ACGTACGTACGT"[:11]):
        assert min_self_match_period(seq) is None
    # a 12-mer of period 4 with only 8 overlap at shift 4: still needs
    # n - shift >= min_overlap, so not reported
    assert min_self_match_period("ACGTACGTACGT") is None


def test_max_search_window_bp_capped():
    assert max_search_window_bp(MS2 * 3, max_window=5000) == 18
    assert max_search_window_bp(RANDOMISH, max_window=5000) == 5000


def test_counter_increments():
    reset_counters()
    assess_overhang(MS2 * 3, max_window=5000)
    assess_overhang(RANDOMISH, max_window=5000)
    assert COUNTERS['period_bounded'] == 1
    assert COUNTERS['assessed'] == 2


def test_period_cap_never_widens():
    # the period term may only shrink the window
    for seq in (MS2 * 3, RANDOMISH, "GATCCTAG" * 8):
        with_cap = max_search_window_bp(seq, max_window=5000)
        assert with_cap <= 5000
