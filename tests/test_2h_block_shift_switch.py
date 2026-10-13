"""ISSUE-047 whole-block shifts: RECTIFY_2H_BLOCK_SHIFT=on (default) | report | off (queue r1 D3).

Master applies whole-block shifts whenever Module 2H runs (cd200a9); the filing asked for a read review before any
default-on. 'report' plans each shift and counts it (whole_block_shifts_reported_not_written) without writing it,
so the class can be sized and carded; 'off' skips the planning. The default stays 'on' until Kevin rules.
"""
from collections import Counter

from rectify.core.splice.block_shift import AtomicReadPlan
from rectify.core.splice.junction_refiner import refine_read_junctions

from tests.test_ultracode_block_shift import block_input


def _run(monkeypatch, mode):
    if mode is None:
        monkeypatch.delenv('RECTIFY_2H_BLOCK_SHIFT', raising=False)
    else:
        monkeypatch.setenv('RECTIFY_2H_BLOCK_SHIFT', mode)
    read, expected, g, pool = block_input()
    counters = Counter()
    changes = refine_read_junctions(read, pool, set(), g, '+', counters=counters)
    return changes, counters


def test_default_writes_the_shift(monkeypatch):
    changes, counters = _run(monkeypatch, None)
    assert len(changes) == 1 and isinstance(changes[0], AtomicReadPlan)
    assert counters['whole_block_shifts'] == 1


def test_report_counts_without_writing(monkeypatch):
    changes, counters = _run(monkeypatch, 'report')
    assert not any(isinstance(c, AtomicReadPlan) for c in changes)
    assert counters['whole_block_shifts_reported_not_written'] == 1 and counters['whole_block_shifts'] == 0


def test_off_neither_writes_nor_counts(monkeypatch):
    changes, counters = _run(monkeypatch, 'off')
    assert not any(isinstance(c, AtomicReadPlan) for c in changes)
    assert not counters
