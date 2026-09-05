#!/usr/bin/env python3
"""A BAM-region worker that dies must FAIL the run, not hang it (BUGS_TO_FIX D8, 2026-08-22).

`multiprocessing.Pool` replaces a worker killed by SIGABRT/SIGSEGV but never re-queues the task
it was running, so ``for r in pool.imap(...)`` blocks forever. Four DRS samples sat idle for
9.7 h each this way (glibc ``free(): invalid next size`` in a native extension). The guard
``_iter_pool_results`` polls the worker roster and raises ``RegionWorkerDied``.

The test really aborts a worker (``os.abort()`` -> SIGABRT, same signal as glibc's fatal error)
under the production ``spawn`` context and bounds the whole thing with a thread deadline so a
regression shows up as a failure, not a hung test session.
"""
from __future__ import annotations

import multiprocessing as mp
import os
import threading
import time

import pytest

from rectify.core.bam.parallel import RegionWorkerDied, _iter_pool_results


def _task(x):
    if x == 3:
        os.abort()          # SIGABRT, exactly what glibc raises on heap corruption
    time.sleep(0.05)
    return x * 10


def _run(items, poll_s):
    ctx = mp.get_context("spawn")
    out, err = [], []
    with ctx.Pool(2) as pool:
        try:
            for r in _iter_pool_results(pool, pool.imap(_task, items), poll_s=poll_s):
                out.append(r)
        except RegionWorkerDied as e:
            err.append(e)
    return out, err


def _with_deadline(fn, seconds):
    box = {}
    t = threading.Thread(target=lambda: box.setdefault("r", fn()), daemon=True)
    t.start(); t.join(seconds)
    assert not t.is_alive(), f"still running after {seconds}s — the dead-worker guard did not fire (hang regression)"
    return box["r"]


def test_dead_worker_raises_instead_of_hanging():
    out, err = _with_deadline(lambda: _run(list(range(8)), poll_s=0.5), seconds=60)
    assert len(err) == 1, f"expected RegionWorkerDied, got results={out} err={err}"
    msg = str(err[0])
    assert "worker(s) died" in msg and ("signal 6" in msg or "replaced" in msg), msg
    # Results before the crash may or may not have been yielded (imap is ordered); what matters
    # is that item 3 never silently disappears into an infinite wait.
    assert 30 not in out


def test_healthy_pool_is_unaffected():
    out, err = _with_deadline(lambda: _run([0, 1, 2, 4, 5], poll_s=0.5), seconds=60)
    assert err == [] and out == [0, 10, 20, 40, 50]
