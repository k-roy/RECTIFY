"""Tests for rectify.core.provenance.cluster.detect_current_cluster."""

import platform
import sys
from unittest import mock

import pytest

from rectify.core.provenance.cluster import detect_current_cluster


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _patch_env(**env_overrides):
    """Context manager: replace os.environ with overrides (full replacement)."""
    return mock.patch.dict("os.environ", env_overrides, clear=True)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    platform.system() != "Darwin",
    reason="M1 detection test only runs on macOS (Kevin's laptop)",
)
def test_detect_m1_on_darwin():
    """On Kevin's M1 MacBook, detect_current_cluster() returns 'm1'."""
    # No env manipulation needed; just confirm the live result.
    result = detect_current_cluster()
    assert result == "m1", (
        f"Expected 'm1' on macOS/arm64 (or Rosetta), got {result!r}. "
        "If this machine is Darwin but not Apple Silicon, update the test."
    )


def test_detect_sherlock_via_env_var():
    """$SHERLOCK=1 in environment → 'sherlock'."""
    with _patch_env(SHERLOCK="1"):
        assert detect_current_cluster() == "sherlock"


def test_detect_sherlock_via_hostname_pattern():
    """Hostname matching sh03-07n10 → 'sherlock'."""
    with mock.patch("socket.gethostname", return_value="sh03-07n10"):
        with _patch_env():
            assert detect_current_cluster() == "sherlock"


def test_detect_sherlock_via_hostname_with_stanford():
    """Hostname containing .stanford. → 'sherlock'."""
    with mock.patch("socket.gethostname", return_value="sh02-11n01.stanford.edu"):
        with _patch_env():
            assert detect_current_cluster() == "sherlock"


def test_detect_hoffman2_via_hostname():
    """Hostname containing 'hoffman2' → 'hoffman2'."""
    with mock.patch("socket.gethostname", return_value="login1.hoffman2.idre.ucla.edu"):
        with _patch_env():
            assert detect_current_cluster() == "hoffman2"


def test_detect_hoffman2_via_node_hostname_and_u_project():
    """n<NNNN> hostname + /u/project exists → 'hoffman2'."""
    with mock.patch("socket.gethostname", return_value="n1821"):
        with mock.patch("os.path.isdir", side_effect=lambda p: p == "/u/project"):
            with _patch_env():
                assert detect_current_cluster() == "hoffman2"


def test_detect_hoffman2_login_node_with_u_project():
    """login<N> hostname + /u/project exists → 'hoffman2'."""
    with mock.patch("socket.gethostname", return_value="login3"):
        with mock.patch("os.path.isdir", side_effect=lambda p: p == "/u/project"):
            with _patch_env():
                assert detect_current_cluster() == "hoffman2"


def test_detect_other_when_no_markers():
    """Unknown hostname + no env vars + no /u/project → 'other' (on non-Darwin)."""
    if platform.system() == "Darwin":
        pytest.skip("Darwin always resolves to m1; can't test 'other' here.")
    with mock.patch("socket.gethostname", return_value="somerandommachine"):
        with mock.patch("os.path.isdir", return_value=False):
            with _patch_env():
                assert detect_current_cluster() == "other"
