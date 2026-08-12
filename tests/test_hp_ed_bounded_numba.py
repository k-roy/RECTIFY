"""Numba twin of hp_edit_distance_bounded — bit-identity + guard behavior.

The kernel is opt-in (RECTIFY_HP_ED_NUMBA, same knob and spawn-RSS rationale
as splice_aware_5prime's kernel). These tests exercise identity via a
subprocess with the env set, so the default test run stays numba-free.
"""

import os
import subprocess
import sys

import pytest

from rectify.core.splice.overhang_informativeness import hp_edit_distance_bounded


def test_default_off_no_kernel():
    if os.environ.get('RECTIFY_HP_ED_NUMBA', '0').strip() in ('', '0', 'false',
                                                              'False', 'no', 'off'):
        from rectify.core.splice import overhang_informativeness as oi
        assert oi._hp_ed_bounded_numba is None


def test_python_path_reference_values():
    # anchors the reference semantics the kernel must reproduce
    assert hp_edit_distance_bounded("", "ACGT") == 4.0
    assert hp_edit_distance_bounded("ACGT", "ACGT") == 0.0
    assert hp_edit_distance_bounded("AAAA", "AA") == 1.0   # HP 0.5-cost deletions
    v = hp_edit_distance_bounded("ACGTACGT", "TTTTTTTT", cutoff=1.0)
    assert v == 2.0  # cutoff + 1.0 on prune


@pytest.mark.skipif(
    subprocess.run([sys.executable, "-c", "import numba"],
                   capture_output=True).returncode != 0,
    reason="numba not installed",
)
def test_kernel_bit_identity_subprocess():
    code = r"""
import random
from rectify.core.splice import overhang_informativeness as oi
assert oi._hp_ed_bounded_numba is not None
def py_ref(s1, s2, cutoff):
    saved = oi._hp_ed_bounded_numba
    oi._hp_ed_bounded_numba = None
    try:
        return oi.hp_edit_distance_bounded(s1, s2, cutoff)
    finally:
        oi._hp_ed_bounded_numba = saved
rng = random.Random(1234)
bases = "ACGTacgtN"
for _ in range(800):
    s1 = "".join(rng.choice(bases) for _ in range(rng.randint(0, 120)))
    s2 = "".join(rng.choice(bases) for _ in range(rng.randint(0, 260)))
    if rng.random() < 0.4:
        s1 += "A" * rng.randint(0, 30)
    cutoff = rng.choice([-1.0, 0.0, 3.0, 25.0, 60.0])
    a = oi.hp_edit_distance_bounded(s1, s2, cutoff)
    b = py_ref(s1, s2, cutoff)
    assert a == b, (s1, s2, cutoff, a, b)
print("OK")
"""
    env = dict(os.environ, RECTIFY_HP_ED_NUMBA="1")
    r = subprocess.run([sys.executable, "-c", code], env=env,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "OK" in r.stdout
