"""Bundled penalty-table selection for `rectify correct` (ISSUE-005).

`correct` used to load NO empirical penalty table at all unless the user passed
one — the only auto-selection was the ONT-cDNA per-UMI-bin branch — so every DRS
run scored junctions with flat edit costs and no STR context, silently.  The
tables are per organism and per protocol; the one thing that must never happen
is quietly scoring one organism's reads with another organism's calibration.

Author: Kevin R. Roy (agent S1)
"""

import logging
import sys
from pathlib import Path

import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.commands.correct_command import (  # noqa: E402
    penalty_table_protocol,
    select_penalty_tables,
)


# ---------------------------------------------------------------------------
# protocol mapping
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("config,expected", [
    ({}, 'drs'),
    ({'ont_cDNA': True}, 'cdna'),
    ({'dt_primed_cDNA': True}, 'qsrev'),
    ({'is_netseq': True, 'dt_primed_cDNA': True}, 'qsrev'),
    # ONT PCR-cDNA wins over the antisense flag it also sets
    ({'ont_cDNA': True, 'dt_primed_cDNA': True}, 'cdna'),
])
def test_penalty_table_protocol(config, expected):
    assert penalty_table_protocol(config) == expected


# ---------------------------------------------------------------------------
# selection
# ---------------------------------------------------------------------------

def test_human_drs_selects_the_bundled_human_table():
    jpt, spt = select_penalty_tables({'organism': 'homo_sapiens'})
    assert jpt is not None and 'homo_sapiens' in jpt
    assert Path(jpt).exists()
    assert spt is not None and 'homo_sapiens' in spt


def test_yeast_drs_selects_the_bundled_yeast_table():
    jpt, _spt = select_penalty_tables({'organism': 'saccharomyces_cerevisiae'})
    assert jpt is not None and 'saccharomyces_cerevisiae' in jpt


def test_organism_aliases_resolve():
    a, _ = select_penalty_tables({'organism': 'human'})
    b, _ = select_penalty_tables({'organism': 'homo_sapiens'})
    assert a == b


def test_protocol_changes_the_selected_table():
    drs, _ = select_penalty_tables({'organism': 'homo_sapiens'})
    cdna, _ = select_penalty_tables({'organism': 'homo_sapiens', 'ont_cDNA': True})
    qsrev, _ = select_penalty_tables({'organism': 'homo_sapiens',
                                      'dt_primed_cDNA': True})
    assert drs != cdna != qsrev
    assert 'cdna' in Path(cdna).name
    assert 'qsrev' in Path(qsrev).name


def test_unbundled_organism_never_borrows_another_organisms_table(caplog):
    """The whole point: no silent cross-organism calibration."""
    with caplog.at_level(logging.WARNING):
        jpt, spt = select_penalty_tables({'organism': 'danio_rerio'})
    assert jpt is None and spt is None
    assert 'danio_rerio' in caplog.text, caplog.text
    assert 'saccharomyces' not in caplog.text


def test_unknown_organism_and_no_genome_warns_and_selects_nothing(caplog):
    with caplog.at_level(logging.WARNING):
        jpt, spt = select_penalty_tables({})
    assert (jpt, spt) == (None, None)
    assert 'organism unknown' in caplog.text


def test_explicit_user_table_wins():
    jpt, spt = select_penalty_tables({
        'organism': 'homo_sapiens',
        'junction_penalty_table': '/tmp/mine.tsv',
        'str_penalty_table': '/tmp/mine_str.tsv',
    })
    assert (jpt, spt) == ('/tmp/mine.tsv', '/tmp/mine_str.tsv')


def test_selection_is_logged_at_info(caplog):
    """A silent change of cost scale is the defect; the log line is the fix."""
    with caplog.at_level(logging.INFO):
        select_penalty_tables({'organism': 'homo_sapiens'})
    assert 'penalty table' in caplog.text.lower()
    assert 'homo_sapiens' in caplog.text
