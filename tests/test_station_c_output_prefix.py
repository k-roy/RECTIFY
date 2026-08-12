"""`pool-gate -o` is a PREFIX, and a dotted prefix must not be truncated.

Two sessions (668-drs-arm and cdna-trim-fix) each lost a pass to `-o` being a
prefix rather than a directory, which prompted a look at the writer. The help
string was already correct; the WRITER was not.

`write_pool_gate_outputs` used `Path.with_suffix('.pool_gate.tsv')`, and
`with_suffix` REPLACES an existing suffix rather than appending. So a versioned
prefix silently collapsed:

    -o run/sample.v1  ->  run/sample.pool_gate.tsv
    -o run/sample.v2  ->  run/sample.pool_gate.tsv   # same file, overwritten

Two runs, one output, no warning — a silent-overwrite bug, not a cosmetic one.
"""

import os
from pathlib import Path

import pytest

from rectify.core.consensus.station_c import write_pool_gate_outputs

_ROWS: list = []
_SUMMARY: dict = {'n_junctions': 0}


def test_dotted_prefix_is_not_truncated(tmp_path):
    tsv, js = write_pool_gate_outputs(_ROWS, _SUMMARY, tmp_path / 'sample.v2')
    assert tsv.name == 'sample.v2.pool_gate.tsv'
    assert js.name == 'sample.v2.pool_gate.json'


def test_versioned_prefixes_do_not_collide(tmp_path):
    """The actual failure mode: v1 and v2 writing to one file."""
    tsv1, _ = write_pool_gate_outputs(_ROWS, _SUMMARY, tmp_path / 'sample.v1')
    tsv2, _ = write_pool_gate_outputs(_ROWS, _SUMMARY, tmp_path / 'sample.v2')
    assert tsv1 != tsv2
    assert tsv1.exists() and tsv2.exists()


def test_plain_prefix_still_works(tmp_path):
    tsv, js = write_pool_gate_outputs(_ROWS, _SUMMARY, tmp_path / 'wt_rep1')
    assert tsv.name == 'wt_rep1.pool_gate.tsv'
    assert js.name == 'wt_rep1.pool_gate.json'
    assert tsv.parent == tmp_path


def test_nested_prefix_creates_its_parent(tmp_path):
    tsv, _ = write_pool_gate_outputs(_ROWS, _SUMMARY, tmp_path / 'a' / 'b' / 's')
    assert tsv == tmp_path / 'a' / 'b' / 's.pool_gate.tsv'
    assert tsv.exists()


@pytest.mark.parametrize('as_str', [False, True])
def test_directory_prefix_is_refused_with_an_actionable_message(tmp_path, as_str):
    """Silently writing `<dir>.pool_gate.tsv` BESIDE the directory is what
    confused two sessions. Refuse, and say what to pass instead."""
    d = tmp_path / 'outdir'
    d.mkdir()
    target = str(d) + os.sep if as_str else d
    with pytest.raises(ValueError) as exc:
        write_pool_gate_outputs(_ROWS, _SUMMARY, target)
    msg = str(exc.value)
    assert 'PREFIX' in msg and 'not a directory' in msg
    assert 'pool_gate.tsv' in msg, "message should show the resulting filename"
