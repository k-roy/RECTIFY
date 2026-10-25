"""abPOA abort isolation: a cluster that kills abPOA must not kill the library.

Context. abPOA sizes its DP allocation from the partial-order graph and, for some inputs, asks for a
single 2^38-byte block and **aborts the process** when ``posix_memalign`` fails. Inside
``correct-cdna``'s process pool that killed a worker and failed the whole library; two libraries of a
35-library campaign were lost (Chanfreau ``[[1190]]`` tasks 8, 28). No cap fixes it — a cluster-size
cap is the wrong axis (16,384 x 2 kb is fine, 4 x 100 kb is not) and no length cap is viable because
a 286,560 nt consensus succeeded in production (``[[1195b]]``). So the fix catches the abort.

These tests do not need abPOA installed: the point of the design is the process boundary, so the
tests exercise the boundary with a stubbed child.
"""
from __future__ import annotations

import sys

import pytest

from rectify.core.cdna import consensus as C


@pytest.fixture(autouse=True)
def _reset_counter():
    C.poa_aborted_count(reset=True)
    yield
    C.poa_aborted_count(reset=True)


def _seqs(n, length):
    return ["ACGT" * (length // 4) for _ in range(n)]


def test_short_cluster_stays_in_process(monkeypatch):
    """Behaviour preservation: below the gate nothing changes, and no child is spawned."""
    monkeypatch.setattr(C, "HAS_POA", True)
    monkeypatch.setattr(C, "_poa_in_process", lambda s: "INPROC")
    monkeypatch.setattr(C, "_poa_isolated", lambda s: pytest.fail(
        "a short cluster must NOT be isolated; that would pay a process per cluster"))
    assert C.poa_consensus_from_strings(_seqs(4, 1000)) == "INPROC"
    assert C.poa_aborted_count() == 0


def test_long_cluster_is_isolated(monkeypatch):
    """At or above the gate the call goes through the child, carrying the same sequences."""
    monkeypatch.setattr(C, "HAS_POA", True)
    monkeypatch.setattr(C, "_poa_in_process", lambda s: pytest.fail(
        "a long cluster must NOT run in-process; an abort would kill the worker"))
    seen = {}
    monkeypatch.setattr(C, "_poa_isolated", lambda s: (seen.update(n=len(s)), "ISOLATED")[1])
    assert C.poa_consensus_from_strings(_seqs(4, C.POA_ISOLATE_MIN_LEN)) == "ISOLATED"
    assert seen["n"] == 4


def test_child_death_falls_back_and_is_counted(monkeypatch):
    """A dying child yields None (so callers fall back) and increments the counter."""
    monkeypatch.setattr(C, "HAS_POA", True)

    class _Dead:
        returncode = -6              # SIGABRT, which is how abPOA actually dies
        stdout = b""
        stderr = b"[SIMDMalloc] posix_memalign fail!\nSize: 274877906944, Error: ENOMEM\n"

    monkeypatch.setattr("subprocess.run", lambda *a, **k: _Dead())
    assert C._poa_isolated(_seqs(4, 100_000)) is None
    assert C.poa_aborted_count() == 1, "an abort that is not counted is an invisible fallback"


def test_child_success_returns_the_consensus(monkeypatch):
    monkeypatch.setattr(C, "HAS_POA", True)

    class _Ok:
        returncode = 0
        stdout = b"ACGTACGT"
        stderr = b""

    monkeypatch.setattr("subprocess.run", lambda *a, **k: _Ok())
    assert C._poa_isolated(_seqs(4, 100_000)) == "ACGTACGT"
    assert C.poa_aborted_count() == 0


def test_child_empty_payload_is_no_consensus(monkeypatch):
    """Exit 0 with no payload means abPOA declined — NOT an abort, so it must not be counted."""
    monkeypatch.setattr(C, "HAS_POA", True)

    class _Empty:
        returncode = 0
        stdout = b""
        stderr = b""

    monkeypatch.setattr("subprocess.run", lambda *a, **k: _Empty())
    assert C._poa_isolated(_seqs(4, 100_000)) is None
    assert C.poa_aborted_count() == 0


def test_spawn_failure_is_also_counted(monkeypatch):
    monkeypatch.setattr(C, "HAS_POA", True)

    def _boom(*a, **k):
        raise OSError("no fork for you")

    monkeypatch.setattr("subprocess.run", _boom)
    assert C._poa_isolated(_seqs(4, 100_000)) is None
    assert C.poa_aborted_count() == 1


def test_counter_reset_is_per_region(monkeypatch):
    """io.write_stage1_fastq reads-and-resets, so a region reports only its own aborts."""
    monkeypatch.setattr(C, "HAS_POA", True)
    C._POA_ABORTED[0] = 3
    assert C.poa_aborted_count(reset=True) == 3
    assert C.poa_aborted_count() == 0


@pytest.mark.skipif(not C.HAS_POA, reason="pyabpoa not installed in this interpreter")
def test_real_abort_is_survived():
    """🔴 THE ACCEPTANCE TEST. A real 4 x 100,000 nt cluster aborts abPOA in production.

    This asserts the parent SURVIVES it, which is the whole defect: before the fix the same input
    killed the pool worker and the library with it. It must not be skipped on a machine with abPOA,
    because every other test here stubs the boundary it is meant to prove.
    """
    import random
    random.seed(1195)
    base = "".join(random.choice("ACGT") for _ in range(100_000))

    def mut(s):
        b = list(s)
        for _ in range(500):
            b[random.randrange(len(b))] = random.choice("ACGT")
        return "".join(b)

    out = C.poa_consensus_from_strings([mut(base) for _ in range(4)])
    # Either abPOA coped (then we have a consensus) or it aborted (then we fell back and counted
    # it). What is NOT acceptable is this process dying, which is what the test proves by reaching
    # the assert at all.
    assert out is None or isinstance(out, str)
    if out is None:
        assert C.poa_aborted_count() >= 1, (
            "abPOA produced nothing yet no abort was counted — the fallback would be silent")
