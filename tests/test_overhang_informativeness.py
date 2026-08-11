"""Unit tests for the shared informativeness gate (planning/641 §2).

Covers the SPEC's own worked examples (T1 semantics: the poly(A) overhang is
refused before any candidate is touched), the alpha knob's monotonicity (T7),
the ambiguity canonicalisation (planning/621), and the exact-cutoff contract
of the bounded HP edit distance (planning/596 T2 argument).
"""

import math
import random

import pytest

from rectify.core.splice.overhang_informativeness import (
    COUNTERS,
    OverhangAssessment,
    ambiguity_window,
    assess_overhang,
    canonicalize_junction,
    effective_information_bits,
    hp_edit_distance_bounded,
    max_search_window_bp,
    reset_counters,
    same_junction,
)

POLYA16 = 'AAAAAAAAAAAATAAA'   # the SPEC's canonical low-information overhang


@pytest.fixture(autouse=True)
def _fresh_counters():
    reset_counters()
    yield
    reset_counters()


def _rand(n, seed):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


# ---------------------------------------------------------------------------
# I_eff / W_max — the SPEC §2 table
# ---------------------------------------------------------------------------

class TestInformationContent:
    def test_spec_polya_example_is_about_5_bits(self):
        i = effective_information_bits(POLYA16)
        assert 3.0 < i < 6.0  # SPEC: "I_eff ~ 5 bits"

    def test_polya_refused_at_default_alpha(self):
        assert max_search_window_bp(POLYA16, alpha=0.01, max_window=2000) == 0

    def test_pure_homopolymer_refused(self):
        assert max_search_window_bp('A' * 30, alpha=0.01) == 0
        assert max_search_window_bp('T' * 100, alpha=0.01) == 0

    def test_dinucleotide_repeat_refused(self):
        # Composition entropy alone would score (AG)8 at 16 bits; the
        # conditional term is what kills it.
        assert max_search_window_bp('AG' * 8, alpha=0.01) == 0

    def test_random_16mer_clamps_to_max_window(self):
        assert max_search_window_bp(_rand(16, 7), alpha=0.01, max_window=2000) == 2000

    def test_random_30mer_clamps_to_max_window(self):
        assert max_search_window_bp(_rand(30, 11), alpha=0.01, max_window=5000) == 5000

    def test_short_overhangs_refused(self):
        # A <=4 nt overhang cannot distinguish its placement anywhere.
        for s in ('', 'A', 'AC', 'ACG', 'ACGT'):
            assert max_search_window_bp(s, alpha=0.01) == 0

    def test_non_acgt_contributes_nothing(self):
        assert effective_information_bits('N' * 30) == 0.0
        # N-padding must not raise the information content
        base = effective_information_bits(POLYA16)
        assert effective_information_bits(POLYA16 + 'NNNN') <= base + 1e-9

    def test_case_insensitive(self):
        assert effective_information_bits(POLYA16.lower()) == \
            effective_information_bits(POLYA16)

    def test_mixed_12mer_gets_small_window(self):
        # Informative but weakly so: a bounded, sub-max window.
        w = max_search_window_bp('ACGGTTACAGTC', alpha=0.01, max_window=5000)
        assert 0 < w < 5000


class TestAlphaKnob:
    """641 T7: alpha is a published knob with a monotone effect."""

    @pytest.mark.parametrize('seq', [
        POLYA16, 'ACGGTTACAGTC', _rand(16, 3), _rand(10, 5), 'AAG' * 10,
    ])
    def test_window_monotone_in_alpha(self, seq):
        w1 = max_search_window_bp(seq, alpha=0.001, max_window=10**9)
        w2 = max_search_window_bp(seq, alpha=0.01, max_window=10**9)
        w3 = max_search_window_bp(seq, alpha=0.1, max_window=10**9)
        assert w1 <= w2 <= w3
        from rectify.core.splice.overhang_informativeness import (
            min_self_match_period,
        )
        period_capped = (min_self_match_period(seq) is not None
                         and w3 <= (min_self_match_period(seq) or 1))
        if 0 < w1 and w3 < 10**9 and not period_capped:
            # In the unclamped regime the ratio is exactly the alpha ratio.
            # A tandem-periodic sequence (e.g. (AAG)n) is instead capped at
            # period-1 INDEPENDENT of alpha — the 2026-08-11 period gate —
            # so only weak monotonicity applies there.
            assert w3 >= 10 * w2 >= 100 * w1 or w1 == 0


class TestAssessOverhang:
    def test_refusal_counted(self):
        a = assess_overhang(POLYA16, alpha=0.01, max_window=2000)
        assert isinstance(a, OverhangAssessment)
        assert a.refused and a.w_max_bp == 0
        assert COUNTERS['assessed'] == 1
        assert COUNTERS['refused'] == 1

    def test_window_bounded_counted(self):
        a = assess_overhang('ACGGTTACAGTC', alpha=0.01, max_window=5000)
        assert not a.refused and 0 < a.w_max_bp < 5000
        assert COUNTERS['window_bounded'] == 1
        assert COUNTERS['refused'] == 0

    def test_complex_overhang_not_bounded(self):
        a = assess_overhang(_rand(30, 13), alpha=0.01, max_window=5000)
        assert not a.refused and a.w_max_bp == 5000
        assert COUNTERS['window_bounded'] == 0


# ---------------------------------------------------------------------------
# Junction ambiguity canonicalisation (planning/621)
# ---------------------------------------------------------------------------

class TestAmbiguity:
    def test_no_ambiguity(self):
        #        0123456789
        seq = 'TTGTACCCAGCA'  # intron [2, 10): GT...AG; seq[1]=T vs seq[9]=G etc.
        l, r = ambiguity_window(seq, 2, 10)
        assert (l, r) == (0, 0)
        assert canonicalize_junction(seq, 2, 10) == (2, 10)

    def test_right_slide(self):
        # G[s] == G[e] => +1 is a no-op.
        #      s=2      e=8
        seq = 'TTGTAAAAGTCC'
        # seq[2]='G', seq[8]='G' -> r_amb >= 1; seq[3]='T', seq[9]='T' -> >= 2
        l, r = ambiguity_window(seq, 2, 8)
        assert r >= 2 and l == 0
        assert canonicalize_junction(seq, 2, 8) == (2, 8)
        # a +1-slid version canonicalises back to the leftmost form
        assert canonicalize_junction(seq, 3, 9) == (2, 8)

    def test_same_junction_within_class(self):
        seq = 'TTGTAAAAGTCC'
        assert same_junction(seq, (2, 8), (3, 9))
        assert not same_junction(seq, (2, 8), (2, 9))   # different length
        assert not same_junction(seq, (2, 8), (4, 11))  # different length

    def test_length_invariant_under_canonicalisation(self):
        seq = _rand(500, 21)
        for s in range(50, 400, 37):
            e = s + 60
            cs, ce = canonicalize_junction(seq, s, e)
            assert ce - cs == e - s


# ---------------------------------------------------------------------------
# Bounded HP edit distance — the planning/596 exactness contract
# ---------------------------------------------------------------------------

def _hp_ed_reference(s1, s2):
    """Unpruned reference implementation (mirrors _hp_edit_distance)."""
    return hp_edit_distance_bounded(s1, s2, cutoff=-1.0)


class TestBoundedEditDistance:
    def test_basics(self):
        assert hp_edit_distance_bounded('', '') == 0.0
        assert hp_edit_distance_bounded('ACGT', '') == 4.0
        assert hp_edit_distance_bounded('', 'ACG') == 3.0
        assert hp_edit_distance_bounded('ACGT', 'ACGT') == 0.0
        assert hp_edit_distance_bounded('ACGT', 'ACCT') == 1.0

    def test_hp_indel_half_cost(self):
        # deleting an A inside an AA run costs 0.5
        assert hp_edit_distance_bounded('AAAG', 'AAG') == 0.5
        assert hp_edit_distance_bounded('AAG', 'AAAG') == 0.5
        # a non-HP indel costs 1.0
        assert hp_edit_distance_bounded('ACG', 'AG') == 1.0

    def test_matches_splice_aware_implementation(self):
        from rectify.core.splice.splice_aware_5prime import _hp_edit_distance
        rng = random.Random(99)
        for _ in range(300):
            n = rng.randint(0, 40)
            m = rng.randint(0, 40)
            s1 = ''.join(rng.choice('AACGT') for _ in range(n))  # A-biased: HP runs
            s2 = ''.join(rng.choice('AACGT') for _ in range(m))
            assert hp_edit_distance_bounded(s1, s2) == _hp_edit_distance(s1, s2)

    def test_cutoff_contract(self):
        """true <= cutoff => exact; true > cutoff => result > cutoff."""
        rng = random.Random(17)
        for _ in range(300):
            n = rng.randint(1, 30)
            m = rng.randint(1, 30)
            s1 = ''.join(rng.choice('AACGGT') for _ in range(n))
            s2 = ''.join(rng.choice('AACGGT') for _ in range(m))
            true = _hp_ed_reference(s1, s2)
            for cutoff in (0.0, true - 0.5, true, true + 0.5, 10.0):
                if cutoff < 0:
                    continue
                got = hp_edit_distance_bounded(s1, s2, cutoff=cutoff)
                if true <= cutoff:
                    assert got == true, (s1, s2, cutoff)
                else:
                    assert got > cutoff, (s1, s2, cutoff)

    def test_equal_to_cutoff_never_pruned(self):
        # ED exactly at the cutoff must come back exact (tiebreakers depend on it)
        assert hp_edit_distance_bounded('ACGT', 'ACCT', cutoff=1.0) == 1.0
