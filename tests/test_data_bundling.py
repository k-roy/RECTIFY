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


# ---------------------------------------------------------------------------
# Phase B: nested-shape resolver + per-UMI-bin fallback
# ---------------------------------------------------------------------------

def test_cdna_pooled_when_bin_not_specified():
    """Without a ``bin`` kwarg, cdna resolves to the pooled penalty_scores_cdna.tsv."""
    p = d.get_bundled_junction_penalty_table(ORG, protocol='cdna')
    assert p is not None
    assert p.exists()
    assert p.name == 'penalty_scores_cdna.tsv'


@pytest.mark.parametrize('umi_bin', ['umi1', 'umi2', 'umi3plus'])
def test_cdna_per_bin_resolves_to_specific_file(umi_bin):
    """Per-bin TSVs are bundled (Phase D complete); resolver must return the
    bin-specific file, not fall back to the pooled table."""
    bin_specific = (d.BUNDLED_DATA_DIR /
                    f'genomes/{ORG}/penalty_tables/penalty_scores_cdna_{umi_bin}.tsv')
    assert bin_specific.exists(), f'{bin_specific} missing — re-run Phase D bundling'

    p = d.get_bundled_junction_penalty_table(ORG, protocol='cdna', bin=umi_bin)
    assert p is not None
    assert p.name == f'penalty_scores_cdna_{umi_bin}.tsv'


def test_qsrev_flat_protocols_still_resolves():
    """Flat-shape detection: qsrev keeps its ``key`` directly on the protocol
    map (no ``pooled`` nesting). Resolver must still return the qsrev table."""
    p = d.get_bundled_junction_penalty_table(ORG, protocol='qsrev')
    assert p is not None
    assert p.exists()
    assert p.name == 'penalty_scores_qsrev.tsv'


def test_drs_no_protocol_unchanged():
    """Org-level flat key (no protocol) returns the DRS penalty_scores.tsv."""
    p = d.get_bundled_junction_penalty_table(ORG)
    assert p is not None
    assert p.exists()
    assert p.name == 'penalty_scores.tsv'


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


# ---------------------------------------------------------------------------
# Phase C smoke test: ONT-cDNA per-UMI-bin PenaltyTableSet
# ---------------------------------------------------------------------------

def test_ont_cdna_penalty_table_set_routes_by_xc_tag():
    """Smoke test: bundled per-bin tables load and route correctly by XC tag.

    Exercises the full chain used by ``rectify correct --ONT-cDNA``:
      1. get_bundled_junction_penalty_table resolves all four bin paths on disk
      2. HpPenaltyTable.from_tsv loads each table (with shared STR fallback path)
      3. PenaltyTableSet.for_read routes XC=1→umi1, XC=2→umi2, XC≥3→umi3plus,
         absent→pooled
    """
    import pysam
    from rectify.core.splice.junction_refiner import PenaltyTableSet
    from rectify.core.splice.hp_penalty import HpPenaltyTable

    _str_p = d.get_bundled_str_penalty_table(ORG, protocol='cdna')
    _str_path = str(_str_p) if _str_p is not None else None

    bins = {}
    for bin_name in ('umi1', 'umi2', 'umi3plus', 'pooled'):
        p = d.get_bundled_junction_penalty_table(ORG, protocol='cdna', bin=bin_name)
        assert p is not None, f"Bundled per-bin table missing for bin='{bin_name}'"
        assert p.exists(), f"Bundled per-bin file not on disk: {p}"
        bins[bin_name] = HpPenaltyTable.from_tsv(str(p), str_path=_str_path)

    pts = PenaltyTableSet(bins=bins, umi_tag='XC')

    def _read(xc):
        r = pysam.AlignedSegment()
        r.query_name = 'smoke'
        r.query_sequence = 'ACGT'
        r.flag = 0
        if xc is not None:
            r.set_tag('XC', xc)
        return r

    assert pts.for_read(_read(1))    is bins['umi1'],     "XC=1 → umi1"
    assert pts.for_read(_read(2))    is bins['umi2'],     "XC=2 → umi2"
    assert pts.for_read(_read(3))    is bins['umi3plus'], "XC=3 → umi3plus"
    assert pts.for_read(_read(10))   is bins['umi3plus'], "XC=10 → umi3plus"
    assert pts.for_read(_read(None)) is bins['pooled'],   "missing XC → pooled"
