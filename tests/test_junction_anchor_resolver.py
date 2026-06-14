"""Tests for the junction-anchor gate CLI resolver + per-organism default + the
split-script generation that bakes the resolved value into the chunk-merge job.

The anchor gate itself (the consensus winner-selection branch) is covered by
``test_junction_anchor_gate.py``. THIS file covers the *wiring* the advisor
flagged as the real risk: yeast→0 is byte-identical by construction, so a green
yeast suite says nothing about whether the human flip actually fires. These
tests assert the human path resolves to 10 AND that 10 reaches the generated
chunk-merge script through the production ``_generate_scripts`` path.
"""
import argparse
from pathlib import Path

import pytest

from rectify.data import (
    resolve_min_junction_anchor_bp,
    default_min_junction_anchor_bp,
    add_junction_anchor_args,
)
from rectify.core.commands.split_command import (
    _chunk_merge_body,
    create_split_parser,
    _generate_scripts,
)


def _ns(**kw):
    base = dict(organism=None, min_junction_anchor_bp=None)
    base.update(kw)
    return argparse.Namespace(**base)


# ── per-organism default ────────────────────────────────────────────────────
@pytest.mark.parametrize("organism,expected", [
    ('homo_sapiens', 10),
    ('human', 10),
    ('hg38', 10),
    ('grch38', 10),
    ('saccharomyces_cerevisiae', 0),
    ('yeast', 0),
    ('saccer', 0),
    (None, 0),
    ('mus_musculus', 0),      # unknown organism → conservative gate-off
])
def test_default_per_organism(organism, expected):
    assert default_min_junction_anchor_bp(organism) == expected


# ── resolver: explicit wins, else per-organism ─────────────────────────────
def test_resolver_human_default_on():
    assert resolve_min_junction_anchor_bp(_ns(organism='human')) == 10


def test_resolver_yeast_default_off():
    assert resolve_min_junction_anchor_bp(_ns(organism='saccharomyces_cerevisiae')) == 0


def test_resolver_custom_genome_default_off():
    # human full-genome run that forgot --organism → gate stays off (footgun the
    # split-time guard warns about), NOT silently 10.
    assert resolve_min_junction_anchor_bp(_ns(organism=None)) == 0


def test_resolver_explicit_overrides_organism():
    assert resolve_min_junction_anchor_bp(_ns(organism='human', min_junction_anchor_bp=5)) == 5
    # explicit 0 must defeat the human default (operator deliberately disabling)
    assert resolve_min_junction_anchor_bp(_ns(organism='human', min_junction_anchor_bp=0)) == 0


def test_resolver_bare_namespace_is_safe():
    assert resolve_min_junction_anchor_bp(argparse.Namespace()) == 0


# ── arg declaration ─────────────────────────────────────────────────────────
def test_add_arg_default_is_none_sentinel():
    p = argparse.ArgumentParser()
    add_junction_anchor_args(p)
    assert p.parse_args([]).min_junction_anchor_bp is None
    assert p.parse_args(['--min-junction-anchor-bp', '12']).min_junction_anchor_bp == 12


# ── chunk-merge body bakes the value ────────────────────────────────────────
def test_chunk_merge_body_bakes_value():
    b = _chunk_merge_body(2, ['mapPacBio', 'minimap2'], 's', Path('/tmp/o'),
                          '/g.fa', '/a.gff', 'python', '/src',
                          min_junction_anchor_bp=10)
    assert 'min_junction_anchor_bp=10,' in b


# ── integration: real generation path flips human ON ───────────────────────
def _gen(tmp_path, argv_extra):
    sub = argparse.ArgumentParser(prog='rectify').add_subparsers(dest='command')
    sp = create_split_parser(sub)
    argv = ['in.fastq', '-o', str(tmp_path), '--generate-slurm'] + argv_extra
    args = sp.parse_args(argv)
    args.output_dir = tmp_path
    _generate_scripts(args, n_chunks=2, sample_prefix='smpl')
    return (tmp_path / 'run_array_chunk_merge.sh').read_text()


def test_generate_human_bakes_10(tmp_path):
    script = _gen(tmp_path, ['--organism', 'human'])
    assert 'min_junction_anchor_bp=10,' in script


def test_generate_yeast_bakes_0(tmp_path):
    script = _gen(tmp_path, ['--Scer'])
    assert 'min_junction_anchor_bp=0,' in script


def test_generate_explicit_overrides(tmp_path):
    script = _gen(tmp_path, ['--organism', 'human', '--min-junction-anchor-bp', '7'])
    assert 'min_junction_anchor_bp=7,' in script


# ── co-activation guard: mapPacBio + gate-off warns ────────────────────────
def test_guard_warns_mappacbio_gate_off(tmp_path, capsys):
    _gen(tmp_path, [])                       # organism None → gate 0, mapPacBio in by default
    err = capsys.readouterr().err
    assert 'mapPacBio is in the aligner panel but the junction-anchor' in err


def test_guard_silent_when_human(tmp_path, capsys):
    _gen(tmp_path, ['--organism', 'human'])  # gate 10 → no warning
    err = capsys.readouterr().err
    assert 'spurious-intron gaming vector' not in err


def test_guard_silent_when_mappacbio_skipped(tmp_path, capsys):
    _gen(tmp_path, ['--skip-map-pacbio'])    # gate 0 but no mapPacBio → no warning
    err = capsys.readouterr().err
    assert 'spurious-intron gaming vector' not in err
