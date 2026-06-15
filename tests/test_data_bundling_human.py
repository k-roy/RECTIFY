"""Tests for the homo_sapiens bundled penalty-table arm.

Deliberately NOT a clone of test_data_bundling.py (yeast). The assertions here
pin the places the human arm *differs* from yeast, so a registry edit that
silently breaks human resolution fails loudly:

  - Human ``cdna`` is FLAT-shaped (``{'junction_penalty_table': ...}``), whereas
    yeast ``cdna`` is NESTED (``pooled`` / ``umi1`` / ``umi2`` / ``umi3plus``).
    There are no human UMI-bin tables → a ``bin=`` request must fall back to the
    pooled cdna table, the OPPOSITE of yeast's per-bin resolution.
  - Human ``cdna`` / ``qsrev`` carry only ``junction_penalty_table``; their STR
    and overhang tables fall back to the human DRS tables.
  - ``qsrev`` junction table exists for human (not a fallback).
  - Aliases human / hg38 / grch38 resolve to the same files as homo_sapiens.

Genome + annotation are intentionally NOT bundled for human (too large); only the
penalty tables ship. These tests therefore never touch get_bundled_genome_path.
"""
from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import pytest

from rectify import data as d


ORG = 'homo_sapiens'


# ---------------------------------------------------------------------------
# DRS-family tables are all bundled (HP, STR, overhang)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('getter,expected_name', [
    (d.get_bundled_junction_penalty_table, 'penalty_scores.tsv'),
    (d.get_bundled_str_penalty_table,      'str_penalty_scores.tsv'),
    (d.get_bundled_junction_overhang_table, 'junction_overhang_table.tsv'),
])
def test_drs_tables_bundled(getter, expected_name):
    """All three DRS-family tables ship for human and resolve with no protocol."""
    p = getter(ORG)
    assert p is not None, f'{getter.__name__} returned None'
    assert p.exists(), f'{p} does not exist'
    assert p.name == expected_name
    assert 'homo_sapiens' in str(p)


def test_drs_protocol_matches_back_compat():
    """protocol='drs' resolves to the same files as protocol=None."""
    for getter in (d.get_bundled_junction_penalty_table,
                   d.get_bundled_str_penalty_table,
                   d.get_bundled_junction_overhang_table):
        assert getter(ORG, protocol='drs') == getter(ORG)


# ---------------------------------------------------------------------------
# cDNA: FLAT shape, junction-only, no UMI bins (differs from yeast)
# ---------------------------------------------------------------------------

def test_cdna_junction_table_is_human_specific():
    """protocol='cdna' returns the human cdna table, not a DRS fallback."""
    p = d.get_bundled_junction_penalty_table(ORG, protocol='cdna')
    assert p is not None and p.exists()
    assert p.name == 'penalty_scores_cdna.tsv'


@pytest.mark.parametrize('umi_bin', ['umi1', 'umi2', 'umi3plus'])
def test_cdna_bin_falls_back_to_pooled_cdna(umi_bin):
    """Human has NO per-UMI-bin cdna tables (unlike yeast). A bin= request must
    fall back to the pooled cdna table, never to a bin-specific file."""
    # Guard: confirm the human arm genuinely lacks the per-bin file, so this test
    # exercises the fallback rather than silently passing if one is later added.
    bin_specific = (d.BUNDLED_DATA_DIR /
                    f'genomes/{ORG}/penalty_tables/penalty_scores_cdna_{umi_bin}.tsv')
    assert not bin_specific.exists(), (
        f'{bin_specific} now exists — human grew per-UMI-bin cdna tables; update '
        'this test and the registry to route bins to the specific files.'
    )

    p = d.get_bundled_junction_penalty_table(ORG, protocol='cdna', bin=umi_bin)
    assert p is not None and p.exists()
    assert p.name == 'penalty_scores_cdna.tsv', (
        f"bin={umi_bin} should fall back to pooled cdna table, got {p.name}"
    )


def test_cdna_str_and_overhang_fall_back_to_drs():
    """Human cdna defines only junction_penalty_table; STR + overhang fall back
    to the human DRS tables."""
    assert (d.get_bundled_str_penalty_table(ORG, protocol='cdna')
            == d.get_bundled_str_penalty_table(ORG, protocol='drs'))
    assert (d.get_bundled_junction_overhang_table(ORG, protocol='cdna')
            == d.get_bundled_junction_overhang_table(ORG, protocol='drs'))


# ---------------------------------------------------------------------------
# QuantSeq REV: junction table exists; STR + overhang fall back to DRS
# ---------------------------------------------------------------------------

def test_qsrev_junction_table_is_human_specific():
    p = d.get_bundled_junction_penalty_table(ORG, protocol='qsrev')
    assert p is not None and p.exists()
    assert p.name == 'penalty_scores_qsrev.tsv'


def test_qsrev_str_and_overhang_fall_back_to_drs():
    assert (d.get_bundled_str_penalty_table(ORG, protocol='qsrev')
            == d.get_bundled_str_penalty_table(ORG, protocol='drs'))
    assert (d.get_bundled_junction_overhang_table(ORG, protocol='qsrev')
            == d.get_bundled_junction_overhang_table(ORG, protocol='drs'))


# ---------------------------------------------------------------------------
# Organism aliases resolve to identical paths
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('alias', ['human', 'hg38', 'grch38'])
def test_aliases_resolve_to_homo_sapiens(alias):
    for getter in (d.get_bundled_junction_penalty_table,
                   d.get_bundled_str_penalty_table,
                   d.get_bundled_junction_overhang_table):
        assert getter(alias) == getter(ORG), (
            f'alias {alias!r} did not resolve like {ORG!r} for {getter.__name__}'
        )


# ---------------------------------------------------------------------------
# resolve_reference_paths integration (the real consumer path)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Loadability: tables must load through the real consumer, not just resolve.
# The human tables were written by a newer profiler than the loader was last
# exercised against — a renamed column / stray NaN passes every path-existence
# check above but dies at `rectify correct` runtime. This closes the
# "confirm hp_penalty.py loads the human tables" open item (status 2026-05-27).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('protocol', ['drs', 'cdna', 'qsrev'])
def test_human_tables_load_via_consumer(protocol):
    from rectify.core.splice.hp_penalty import HpPenaltyTable

    jp = d.get_bundled_junction_penalty_table(ORG, protocol=protocol)
    str_p = d.get_bundled_str_penalty_table(ORG, protocol=protocol)
    assert jp is not None and jp.exists()

    table = HpPenaltyTable.from_tsv(
        str(jp), str_path=(str(str_p) if str_p is not None else None)
    )
    assert table is not None
    # Exercise a real lookup — the consumer's actual contract, not just construction.
    cost = table.del_cost(4, ref_base='A')
    assert isinstance(cost, float)
    assert cost >= 0.0


def test_resolve_reference_paths_drs_default():
    """--organism homo_sapiens with no protocol flag fills DRS tables."""
    args = Namespace(
        organism='homo_sapiens',
        genome=None, annotation=None,
        short_read=False, dT_primed_cDNA=False,
        junction_penalty_table=None,
        str_penalty_table=None,
        junction_overhang_table=None,
    )
    d.resolve_reference_paths(args, require_genome=False, verbose=False)
    assert args.junction_penalty_table is not None
    assert Path(args.junction_penalty_table).exists()
    assert args.junction_penalty_table.endswith('penalty_scores.tsv')
    assert '_cdna' not in args.junction_penalty_table
    assert '_qsrev' not in args.junction_penalty_table


def test_resolve_reference_paths_cdna_protocol():
    """--organism homo_sapiens --dT-primed-cDNA picks the human cdna junction
    table, with STR + overhang falling back to the human DRS tables."""
    args = Namespace(
        organism='homo_sapiens',
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
    assert args.junction_penalty_table.endswith('penalty_scores_cdna.tsv')
    # STR + overhang have no human cdna variant → DRS files
    assert args.str_penalty_table.endswith('str_penalty_scores.tsv')
    assert args.junction_overhang_table.endswith('junction_overhang_table.tsv')
