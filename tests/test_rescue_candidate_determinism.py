"""Rescue candidate order must not depend on PYTHONHASHSEED (planning/649).

`candidate_junctions` is a set, so iteration order varies between processes;
equal-edit-distance candidates are then resolved by encounter order once the
4-level tiebreaker is exhausted. Measured before the fix: two runs of the same
tree over the same BAM and pool differed in 3/474 corrected rows. These tests
pin the two order-sensitive narrowing steps.
"""
import os
import subprocess
import sys
import textwrap

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _run_with_seed(seed, snippet):
    env = dict(os.environ, PYTHONHASHSEED=str(seed), PYTHONPATH=REPO)
    out = subprocess.run([sys.executable, "-c", textwrap.dedent(snippet)],
                         env=env, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr[-2000:]
    return out.stdout.strip()


# Iterating a set of junction tuples: the order itself is seed-dependent, which
# is the hazard the sorts below defend against.
_ORDER_PROBE = """
    js = {("chrI", 1000 + 7 * i, 1200 + 11 * i) for i in range(40)}
    print([j[1] for j in js])
"""

_SORTED_PROBE = """
    js = {("chrI", 1000 + 7 * i, 1200 + 11 * i) for i in range(40)}
    nearby = list(js)
    nearby.sort(key=lambda _j: (_j[0], _j[1], _j[2]))
    print([j[1] for j in nearby])
"""

_DIST_CAP_PROBE = """
    align5 = 5000
    # Many candidates EQUIDISTANT from align5 — a distance-only sort leaves the
    # surviving slice up to set order; the coordinate tiebreak pins it.
    js = {("chrI", align5 - 100, align5 + d) for d in range(20)}
    js |= {("chrI", align5 - 100, align5 - d) for d in range(20)}
    nd = [(abs(align5 - j[2]), j) for j in js]
    nd.sort(key=lambda x: (x[0], x[1][0], x[1][1], x[1][2]))
    print([j[2] for _, j in nd[:10]])
"""


def test_set_iteration_order_really_is_seed_dependent():
    """Guards the premise: if this ever stops varying, the tests below go vacuous."""
    orders = {_run_with_seed(s, _ORDER_PROBE) for s in (0, 1, 12345)}
    assert len(orders) > 1, "set order no longer varies with PYTHONHASHSEED"


def test_coordinate_sort_is_seed_independent():
    results = {_run_with_seed(s, _SORTED_PROBE) for s in (0, 1, 12345)}
    assert len(results) == 1, f"coordinate sort still seed-dependent: {results}"


def test_distance_cap_slice_is_seed_independent():
    results = {_run_with_seed(s, _DIST_CAP_PROBE) for s in (0, 1, 12345)}
    assert len(results) == 1, f"distance-cap slice still seed-dependent: {results}"


def test_body_sorts_candidates_before_use():
    """The production sort is present at the narrowing site, not just in theory."""
    src = open(os.path.join(
        REPO, "rectify/core/splice/splice_aware_5prime.py")).read()
    assert "_nearby_junctions.sort(key=lambda _j: (_j[0], _j[1], _j[2]))" in src
    assert "key=lambda x: (x[0], x[1][0], x[1][1], x[1][2])" in src
