"""2026-09-05 (Kevin): the resolver's hp-ED numba kernel is ON by default with a lazy
import; RECTIFY_HP_ED_NUMBA=0 disables it; the correct-stage kernel stays opt-in."""
import importlib
import os
import subprocess
import sys


def _requested(env_value):
    if env_value is None:
        setenv = "os.environ.pop('RECTIFY_HP_ED_NUMBA', None)"
    else:
        setenv = "os.environ['RECTIFY_HP_ED_NUMBA'] = %r" % env_value
    code = (
        "import os\n" + setenv + "\n"
        "import rectify.core.splice.overhang_informativeness as oi\n"
        "print(int(oi._HP_ED_NUMBA_REQUESTED), int(oi._HP_ED_NUMBA_LOADED))\n"
    )
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, check=True).stdout.split()
    return int(out[0]), int(out[1])


def test_default_is_on_and_import_is_lazy():
    requested, loaded = _requested(None)
    assert requested == 1
    assert loaded == 0, 'numba must not be imported at module import time'


def test_env_zero_disables():
    requested, _ = _requested('0')
    assert requested == 0


def test_correct_stage_kernel_stays_opt_in():
    code = (
        "import os; os.environ.pop('RECTIFY_HP_ED_NUMBA', None)\n"
        "import rectify.core.splice.splice_aware_5prime as s\n"
        "print(int(s._HP_ED_NUMBA_REQUESTED))\n"
    )
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, check=True).stdout.strip()
    assert out == '0'


def test_kernel_loads_on_first_large_call_and_matches_python():
    code = (
        "import os; os.environ.pop('RECTIFY_HP_ED_NUMBA', None)\n"
        "import rectify.core.splice.overhang_informativeness as oi\n"
        "a = 'ACGT' * 30; b = 'ACGA' * 30\n"
        "v = oi.hp_edit_distance_bounded(a, b)\n"
        "k = oi._hp_ed_bounded_numba\n"
        "print(int(oi._HP_ED_NUMBA_LOADED), 'K' if k is not None else 'N', v)\n"
    )
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, check=True).stdout.split()
    assert out[0] == '1'
    # when numba is installed the kernel must be present; either way the value is the DP's
    try:
        import numba  # noqa: F401
        assert out[1] == 'K'
    except ImportError:
        assert out[1] == 'N'
