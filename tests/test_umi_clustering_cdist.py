"""651: the batched-cdist pair-set path must be byte-identical to the per-pair loop.

The pre-651 behaviour IS the fallback loop (unchanged code), so A/B-ing the batched path
against the fallback is exactly "new == old". Star splitting with precomputed neighbour
sets is A/B-ed against its own distance-verifying form directly.
"""
from __future__ import annotations

import random

import pytest

from rectify.core.umi import clustering
from rectify.core.umi.clustering import (
    _levenshtein_neighbour_sets,
    _split_into_stars,
    umi_components_directional,
)
from rapidfuzz.distance import Levenshtein

FIXED_T = [0, 1, 6, 7, 12, 13, 18, 19, 24, 25, 26]
UMI_LEN = 27


def _make_structured_umi(rng: random.Random) -> str:
    return "".join("T" if j in FIXED_T else rng.choice("ACG") for j in range(UMI_LEN))


def _corrupt(u: str, rng: random.Random, sub_p: float, indel_p: float,
             fixed_len: bool) -> str:
    out = []
    for ch in u:
        r = rng.random()
        if r < indel_p / 2:
            continue  # deletion
        if r < indel_p:
            out.append(rng.choice("ACGT"))  # insertion before ch
        if rng.random() < sub_p:
            out.append(rng.choice([c for c in "ACGT" if c != ch]))
        else:
            out.append(ch)
    s = "".join(out)
    if fixed_len:
        # the cdna extractor slices exactly UMI_LEN, so indels shift content, not length
        s = (s + "AACG")[:UMI_LEN] if len(s) < UMI_LEN else s[:UMI_LEN]
    return s


def _make_bucket(seed: int, n_mol: int, dup: float, fixed_len: bool) -> list:
    rng = random.Random(seed)
    reads = []
    for _ in range(n_mol):
        true = _make_structured_umi(rng)
        k = 1 + (1 if rng.random() < (dup - 1) else 0)
        for _ in range(k):
            reads.append(_corrupt(true, rng, sub_p=0.013, indel_p=0.004,
                                  fixed_len=fixed_len))
    rng.shuffle(reads)
    return reads


@pytest.mark.parametrize("seed", [1, 2, 3])
@pytest.mark.parametrize("max_edit", [2, 3])
@pytest.mark.parametrize("n_mol,dup", [(60, 1.5), (400, 1.13), (1500, 1.13)])
@pytest.mark.parametrize("fixed_len", [True, False])
def test_cdist_path_identical_to_perpair_loop(monkeypatch, seed, max_edit, n_mol,
                                              dup, fixed_len):
    umis = _make_bucket(seed, n_mol, dup, fixed_len)
    got_new = umi_components_directional(umis, max_edit)
    monkeypatch.setattr(clustering, "_levenshtein_neighbour_sets",
                        lambda unique, k: None)
    got_old = umi_components_directional(umis, max_edit)
    assert got_new == got_old


def test_cdist_path_identical_at_ed1_mixed_length(monkeypatch):
    # mixed lengths disable fast_hamming even at ed=1 -> exercises cdist vs loop there too
    umis = _make_bucket(7, 300, 1.3, fixed_len=False)
    umis.append(umis[0] + "AC")
    got_new = umi_components_directional(umis, 1)
    monkeypatch.setattr(clustering, "_levenshtein_neighbour_sets",
                        lambda unique, k: None)
    got_old = umi_components_directional(umis, 1)
    assert got_new == got_old


@pytest.mark.parametrize("seed", [11, 12])
@pytest.mark.parametrize("max_edit", [1, 2, 3])
def test_star_split_with_neighbours_identical(seed, max_edit):
    rng = random.Random(seed)
    uniq = list({_corrupt(_make_structured_umi(rng), rng, 0.05, 0.01, True)
                 for _ in range(300)})
    counts = [rng.randint(1, 5) for _ in uniq]
    neighbours = _levenshtein_neighbour_sets(uniq, max_edit)
    assert neighbours is not None, "compiled rapidfuzz + numpy expected in the test env"
    component = list(range(len(uniq)))
    with_sets = _split_into_stars(component, uniq, counts, max_edit,
                                  neighbours=neighbours)
    without = _split_into_stars(component, uniq, counts, max_edit)
    assert with_sets == without


def test_neighbour_sets_match_bruteforce():
    rng = random.Random(99)
    uniq = list({_corrupt(_make_structured_umi(rng), rng, 0.08, 0.02, False)
                 for _ in range(200)})
    for max_edit in (1, 2, 3):
        got = _levenshtein_neighbour_sets(uniq, max_edit)
        assert got is not None
        expect = {}
        for i in range(len(uniq)):
            for j in range(i + 1, len(uniq)):
                if Levenshtein.distance(uniq[i], uniq[j],
                                        score_cutoff=max_edit) <= max_edit:
                    expect.setdefault(i, set()).add(j)
                    expect.setdefault(j, set()).add(i)
        assert got == expect
