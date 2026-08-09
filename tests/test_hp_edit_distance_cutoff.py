"""Equivalence tests for the exact early-exit cutoff in ``_hp_edit_distance``.

The cutoff is an OPTIMISATION: it must never change which candidate a caller
selects.  ``_rescue_3ss_truncation_body`` compares each candidate's edit
distance against its running best with

    if _ed < _best_local_ed or (_ed == _best_local_ed and _cur < _best):

so an ED *tie* is a live candidate resolved by a 4-level tiebreaker.  The
contract the cutoff must therefore satisfy is:

  * if the true distance is ``<= cutoff``  -> return it EXACTLY
  * if the true distance is ``>  cutoff``  -> return anything ``> cutoff``

Both halves are asserted below.  A cutoff that pruned at equality would satisfy
the second clause but break the first, silently changing junction selection --
which is precisely what these tests exist to catch.
"""

import random

import pytest

from rectify.core.splice import splice_aware_5prime as _mod
from rectify.core.splice.splice_aware_5prime import (
    _HP_ED_MAX_LEN,
    _HP_ED_NUMBA,
    _HP_ED_NUMBA_MIN_CELLS,
    _hp_edit_distance,
)


def _py(s1: str, s2: str, cutoff: float = -1.0) -> float:
    """Force the pure-Python branch, whatever numba is doing.

    ``_hp_edit_distance`` reads the module-level ``_HP_ED_NUMBA`` at call time,
    so flipping it here selects the reference implementation.  This is what lets
    the numba kernel be compared against the loop it must reproduce exactly.
    """
    saved = _mod._HP_ED_NUMBA
    _mod._HP_ED_NUMBA = False
    try:
        return _hp_edit_distance(s1, s2, cutoff)
    finally:
        _mod._HP_ED_NUMBA = saved


def _check(s1: str, s2: str) -> None:
    """Assert the cutoff contract across a sweep of cutoffs for one pair."""
    truth = _hp_edit_distance(s1, s2)

    # No cutoff (default) must be bit-identical to the historical behaviour.
    assert _hp_edit_distance(s1, s2, -1.0) == truth

    # A cutoff at or above the true distance must reproduce it EXACTLY.
    for c in (truth, truth + 0.5, truth + 1.0, truth + 1000.0):
        got = _hp_edit_distance(s1, s2, c)
        assert got == truth, (
            f"cutoff={c} changed an admissible result: {got} != {truth} "
            f"for {s1!r} vs {s2!r}"
        )

    # A cutoff below the true distance may return anything, so long as it is
    # still > cutoff (so the caller's `<` and `==` both correctly fail).
    c = truth - 0.5
    if c >= 0.0:
        got = _hp_edit_distance(s1, s2, c)
        assert got > c, f"pruned return {got} must exceed cutoff {c}"


def test_randomised_pairs_exact():
    """2000 randomised pairs — the bulk equivalence sweep."""
    rng = random.Random(20260807)
    for _ in range(2000):
        n = rng.randint(1, 60)
        m = rng.randint(1, 60)
        s1 = "".join(rng.choice("ACGT") for _ in range(n))
        s2 = "".join(rng.choice("ACGT") for _ in range(m))
        _check(s1, s2)


def test_homopolymer_rich_pairs():
    """HP runs drive the 0.5-cost path — bias the alphabet to exercise it."""
    rng = random.Random(11)
    for _ in range(500):
        s1 = "".join(rng.choice("AAAACG") for _ in range(rng.randint(1, 50)))
        s2 = "".join(rng.choice("AAAACG") for _ in range(rng.randint(1, 50)))
        _check(s1, s2)


@pytest.mark.parametrize(
    "s1,s2",
    [
        ("", ""),
        ("", "ACGT"),
        ("ACGT", ""),
        ("A", "A"),
        ("A", "C"),
        ("ACGT", "ACGT"),                 # identical
        ("AAAA", "AAAA"),                 # HP run, identical
        ("AAAAC", "CAAAA"),               # HP run at each end
        ("AAAAAAAAAA", "A"),              # long HP vs single base
        ("ACGTACGTAC", "TGCATGCATG"),     # maximally dissimilar, equal length
    ],
)
def test_edge_cases(s1, s2):
    _check(s1, s2)


def test_truncation_boundary():
    """Sequences at/over _HP_ED_MAX_LEN are truncated — cutoff must respect that."""
    for n in (_HP_ED_MAX_LEN - 1, _HP_ED_MAX_LEN, _HP_ED_MAX_LEN + 37):
        s1 = "ACGT" * (n // 4 + 1)
        s2 = "ACGA" * (n // 4 + 1)
        _check(s1[:n], s2[:n])


def test_cutoff_actually_prunes():
    """Guard against a no-op: an impossible cutoff must trigger the early exit.

    Uses two long dissimilar sequences whose true distance is large, with a
    cutoff of 0.0 — the row minimum must exceed it well before the final row.
    """
    s1 = "A" * 150
    s2 = "C" * 150
    truth = _hp_edit_distance(s1, s2)
    assert truth > 0.0
    got = _hp_edit_distance(s1, s2, 0.0)
    assert got > 0.0
    assert got != truth or truth <= 1.0  # pruned value need not equal the truth


# ===========================================================================
# Numba-kernel equivalence — the JIT path must reproduce the Python loop EXACTLY
# ===========================================================================

pytestmark_numba = pytest.mark.skipif(
    not _HP_ED_NUMBA, reason="numba not installed in this environment"
)


def _check_numba(s1: str, s2: str) -> None:
    """Numba result must be bit-identical to the Python loop, cutoff or not."""
    ref = _py(s1, s2)
    got = _hp_edit_distance(s1, s2)          # takes the numba path when eligible
    assert got == ref, f"numba {got} != python {ref} for len({len(s1)},{len(s2)})"

    # And with pruning active, across cutoffs either side of the true value.
    for c in (ref, ref + 1.0, max(0.0, ref - 0.5)):
        r_ref = _py(s1, s2, c)
        r_got = _hp_edit_distance(s1, s2, c)
        if c >= ref:
            # Admissible: both must return the exact true distance.
            assert r_got == ref and r_ref == ref, (
                f"cutoff={c}: numba={r_got} python={r_ref} truth={ref}"
            )
        else:
            # Pruned: both need only exceed the cutoff.
            assert r_got > c and r_ref > c


@pytestmark_numba
def test_numba_matches_python_randomised():
    """1000 pairs large enough to cross _HP_ED_NUMBA_MIN_CELLS."""
    rng = random.Random(20260807)
    side = int(_HP_ED_NUMBA_MIN_CELLS ** 0.5) + 5      # guarantees n*m >= threshold
    for _ in range(1000):
        n = rng.randint(side, side + 80)
        m = rng.randint(side, side + 80)
        s1 = "".join(rng.choice("ACGT") for _ in range(n))
        s2 = "".join(rng.choice("ACGT") for _ in range(m))
        _check_numba(s1, s2)


@pytestmark_numba
def test_numba_matches_python_homopolymer():
    """HP-biased alphabet — exercises the 0.5-cost rule in the JIT kernel."""
    rng = random.Random(99)
    for _ in range(400):
        n = rng.randint(25, 90)
        m = rng.randint(25, 90)
        s1 = "".join(rng.choice("AAAACG") for _ in range(n))
        s2 = "".join(rng.choice("AAAACG") for _ in range(m))
        _check_numba(s1, s2)


@pytestmark_numba
def test_numba_matches_python_at_truncation_boundary():
    for n in (_HP_ED_MAX_LEN - 1, _HP_ED_MAX_LEN, _HP_ED_MAX_LEN + 53):
        s1 = ("ACGTT" * (n // 5 + 2))[:n]
        s2 = ("ACGAT" * (n // 5 + 2))[:n]
        _check_numba(s1, s2)


@pytestmark_numba
def test_numba_path_is_actually_taken():
    """Anti-no-op: prove the dispatcher really routes large DPs to the kernel."""
    calls = {"n": 0}
    real = _mod._hp_ed_dp_numba

    def _spy(a1, a2, cutoff):
        calls["n"] += 1
        return real(a1, a2, cutoff)

    _mod._hp_ed_dp_numba = _spy
    try:
        big = "ACGT" * 30          # 120x120 = 14400 cells, well over the threshold
        _hp_edit_distance(big, big)
        assert calls["n"] == 1, "large DP did NOT reach the numba kernel"
        # And small inputs must NOT pay the marshalling cost.
        before = calls["n"]
        _hp_edit_distance("ACGT", "ACGA")     # 16 cells, under the threshold
        assert calls["n"] == before, "small DP wrongly routed to numba"
    finally:
        _mod._hp_ed_dp_numba = real


@pytestmark_numba
def test_warmup_is_safe_and_idempotent():
    _mod._warmup_hp_ed_numba()
    _mod._warmup_hp_ed_numba()
