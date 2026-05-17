"""Tests for protocol-aware bundled penalty table resolution.

Pinned to the S. cerevisiae bundle:
  - The flat (DRS) keys must remain resolvable for back-compat callers.
  - The protocols map must route --Scer --dT-primed-cDNA → cdna table and
    --Scer --short-read --dT-primed-cDNA → qsrev table, falling back to DRS
    when a protocol-specific file is not yet bundled.
"""
from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import pytest

from rectify import data as d


ORG = 'saccharomyces_cerevisiae'


# ---------------------------------------------------------------------------
# Back-compat: protocol=None must return the DRS table
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('getter', [
    d.get_bundled_junction_penalty_table,
    d.get_bundled_str_penalty_table,
    d.get_bundled_junction_overhang_table,
])
def test_back_compat_no_protocol_returns_drs_table(getter):
    """Callers that don't pass `protocol` get the same files as before."""
    p = getter(ORG)
    assert p is not None, f'{getter.__name__} returned None'
    assert p.exists(), f'{p} does not exist'
    # DRS tables have no protocol suffix
    assert '_cdna' not in p.name
    assert '_qsrev' not in p.name


# ---------------------------------------------------------------------------
# Explicit protocol kwarg
# ---------------------------------------------------------------------------

def test_drs_protocol_matches_back_compat():
    """protocol='drs' resolves to the same files as protocol=None."""
    for getter in (d.get_bundled_junction_penalty_table,
                   d.get_bundled_str_penalty_table,
                   d.get_bundled_junction_overhang_table):
        assert getter(ORG, protocol='drs') == getter(ORG)


def test_cdna_protocol_resolves():
    """protocol='cdna' returns a path (either cdna-specific or DRS fallback)."""
    for getter in (d.get_bundled_junction_penalty_table,
                   d.get_bundled_str_penalty_table,
                   d.get_bundled_junction_overhang_table):
        p = getter(ORG, protocol='cdna')
        assert p is not None, f'{getter.__name__}(protocol="cdna") returned None'
        assert p.exists()


def test_qsrev_protocol_resolves():
    """protocol='qsrev' returns a path. Overhang has no qsrev entry per handoff
    — must fall back to the DRS table."""
    for getter in (d.get_bundled_junction_penalty_table,
                   d.get_bundled_str_penalty_table,
                   d.get_bundled_junction_overhang_table):
        p = getter(ORG, protocol='qsrev')
        assert p is not None
        assert p.exists()

    # Verify the QSrev overhang specifically falls back to the DRS file
    qsrev_overhang = d.get_bundled_junction_overhang_table(ORG, protocol='qsrev')
    drs_overhang = d.get_bundled_junction_overhang_table(ORG, protocol='drs')
    assert qsrev_overhang == drs_overhang, (
        'QSrev overhang should fall back to DRS table (handoff §6 — short-read '
        'junctions are too sparse to calibrate)'
    )


# ---------------------------------------------------------------------------
# Protocol-from-args helper
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('short_read,dt_cdna,expected', [
    (False, False, 'drs'),
    (False, True,  'cdna'),
    (True,  True,  'qsrev'),
    (True,  False, 'drs'),   # short_read alone is unspecified → DRS fallback
])
def test_protocol_for_args(short_read, dt_cdna, expected):
    args = Namespace(short_read=short_read, dT_primed_cDNA=dt_cdna)
    assert d._protocol_for_args(args) == expected


def test_protocol_for_args_missing_attrs_defaults_to_drs():
    """Subcommands without short_read / dT_primed_cDNA attrs default to DRS."""
    args = Namespace()
    assert d._protocol_for_args(args) == 'drs'


# ---------------------------------------------------------------------------
# resolve_reference_paths integration
# ---------------------------------------------------------------------------

def test_resolve_reference_paths_picks_cdna_protocol():
    """--Scer --dT-primed-cDNA should pick the cdna-tagged table when it
    exists; otherwise gracefully fall back to the DRS table."""
    args = Namespace(
        organism='saccharomyces_cerevisiae',
        genome=None, annotation=None,
        short_read=False, dT_primed_cDNA=True,
        junction_penalty_table=None,
        str_penalty_table=None,
        junction_overhang_table=None,
    )
    d.resolve_reference_paths(args, require_genome=False, verbose=False)

    for attr in ('junction_penalty_table', 'str_penalty_table',
                 'junction_overhang_table'):
        val = getattr(args, attr)
        assert val is not None, f'{attr} not filled by resolver'
        assert Path(val).exists(), f'{attr} -> {val} does not exist'


def test_resolve_reference_paths_back_compat_no_dt_cdna():
    """Without --dT-primed-cDNA, must pick DRS tables (back-compat behavior)."""
    args = Namespace(
        organism='saccharomyces_cerevisiae',
        genome=None, annotation=None,
        short_read=False, dT_primed_cDNA=False,
        junction_penalty_table=None,
        str_penalty_table=None,
        junction_overhang_table=None,
    )
    d.resolve_reference_paths(args, require_genome=False, verbose=False)

    assert args.junction_penalty_table is not None
    assert 'penalty_scores.tsv' in args.junction_penalty_table
    assert '_cdna' not in args.junction_penalty_table
    assert '_qsrev' not in args.junction_penalty_table
