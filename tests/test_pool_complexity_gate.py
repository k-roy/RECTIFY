"""Pool-admission complexity gate (planning/648 §1 form 2, planning/649).

The gate answers a different question from the search-radius gate that shares
the module: given a junction's OBSERVED span D and its exonic flank
complexity, is the placement better than chance (``D * 2**-I_eff <= alpha``)?
The D-dependence is the property under test — a low-complexity flank is only
disqualifying for a LONG-range junction.
"""
import pytest

from rectify.core.splice.overhang_informativeness import (
    DEFAULT_POOL_FLANK_BP,
    filter_pool_by_flank_complexity,
    flank_information_bits,
    junction_chance_expectation,
)

# A chromosome with two intron-flank contexts at known offsets:
#   exon1 complex  [0:100), intron [100:200), exon2 complex [200:300)
#   poly(A)/T flank block at [300:400) so a junction ending there is low-info.
COMPLEX_L = "GCATTACGGATCCTAGCATCGGATTCAGCTACGATCGATTGCACGTATCAG" * 2   # 102 nt
COMPLEX_R = "TTGACCGATCAGCATGGCTAACGTTCAGGCATCAGTTACGGCATTGACCAT" * 2   # 102 nt
LOWCOMP = "T" * 60 + "A" * 60
CHROM = COMPLEX_L + ("N" * 100) + COMPLEX_R + LOWCOMP
GENOME = {"chrTest": CHROM}

I0 = len(COMPLEX_L)                 # intron start (exon-1 side flank = complex)
I1 = I0 + 100                       # intron end   (exon-2 side flank = complex)
LOW0 = len(COMPLEX_L) + 100 + len(COMPLEX_R) + 20   # inside the low-complexity block


def test_complex_flanks_score_far_above_lowcomplexity():
    complex_bits = flank_information_bits(CHROM, I0, I1)
    low_bits = flank_information_bits(CHROM, LOW0, LOW0 + 20)
    assert complex_bits > 20.0
    assert low_bits < 6.0
    assert complex_bits > 3 * low_bits


def test_worse_of_the_two_flanks_is_used():
    """One complex side must not rescue a low-complexity partner side."""
    # intron ENDS at the low-complexity block: exon-1 flank complex, exon-2 poor
    mixed = flank_information_bits(CHROM, I0, LOW0)
    poor_only = flank_information_bits(CHROM, LOW0 - 100, LOW0)
    assert mixed == pytest.approx(poor_only)


def test_D_dependence_short_intron_survives_a_flank_that_kills_a_long_one():
    """THE point of the E_chance form: the same flank is fine for a short span."""
    short_e = junction_chance_expectation(CHROM, LOW0 - 40, LOW0 - 20)   # D=20
    long_e = junction_chance_expectation(CHROM, LOW0 - 8000, LOW0)       # D=8000
    assert short_e < long_e
    assert long_e / max(short_e, 1e-12) == pytest.approx(8000 / 20, rel=1e-6)


def test_filter_refuses_long_lowcomplexity_and_keeps_complex():
    good = ("chrTest", I0, I1)                      # complex flanks, D=100
    bad = ("chrTest", LOW0 - 8000, LOW0)            # low-complexity flank, D=8000
    kept, n_refused = filter_pool_by_flank_complexity(
        {good, bad}, set(), GENOME, alpha=0.01,
    )
    assert n_refused == 1
    assert good in kept and bad not in kept


def test_annotated_junctions_are_never_refused():
    bad = ("chrTest", LOW0 - 8000, LOW0)
    kept, n_refused = filter_pool_by_flank_complexity(
        {bad}, {bad}, GENOME, alpha=0.01,
    )
    assert n_refused == 0 and bad in kept


def test_missing_chromosome_is_retained_not_dropped():
    """Never drop a junction for lack of evidence."""
    j = ("chrAbsent", 10, 9000)
    kept, n_refused = filter_pool_by_flank_complexity({j}, set(), GENOME, alpha=0.01)
    assert n_refused == 0 and j in kept


def test_alpha_is_monotone():
    """Larger alpha admits a superset — the knob must be well-behaved."""
    pool = {("chrTest", LOW0 - d, LOW0) for d in (50, 200, 1000, 8000)}
    pool |= {("chrTest", I0, I1)}
    sizes = []
    prev = None
    for alpha in (0.0001, 0.01, 1.0, 100.0):
        kept, _ = filter_pool_by_flank_complexity(pool, set(), GENOME, alpha=alpha)
        if prev is not None:
            assert prev <= kept, "a larger alpha must admit a superset"
        prev = kept
        sizes.append(len(kept))
    assert sizes == sorted(sizes)
    assert sizes[0] < sizes[-1], "the sweep must actually move the pool"


def test_input_sets_are_not_mutated():
    pool = {("chrTest", LOW0 - 8000, LOW0), ("chrTest", I0, I1)}
    annot = set()
    before = set(pool)
    filter_pool_by_flank_complexity(pool, annot, GENOME, alpha=0.01)
    assert pool == before and annot == set()


def test_flank_bp_default_matches_465b_prototype():
    assert DEFAULT_POOL_FLANK_BP == 15
