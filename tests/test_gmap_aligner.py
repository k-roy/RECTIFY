"""Smoke tests for the GMAP junction-aligner integration.

GMAP is wired as an opt-in splice-aware aligner alongside uLTRA/deSALT. These
tests cover the wiring (CLI choices + namespace forwarding) and the two
fail-fast guards in run_gmap (missing binary, missing/incomplete db) without
needing a real gmap install or a built database.
"""
import argparse

import pytest

from rectify.cli import create_parser
from rectify.core.align import multi_aligner


def test_gmap_in_run_all_junction_choices():
    p = create_parser()
    ns = p.parse_args([
        'run-all', 'x.bam', '--drs', '--Scer',
        '--base-aligners', 'minimap2', 'mapPacBio',
        '--junction-aligners', 'deSALT', 'uLTRA', 'gmap',
        '--gmap-db', '/tmp/db', '-o', '/tmp/out',
    ])
    assert ns.junction_aligners == ['deSALT', 'uLTRA', 'gmap']
    assert ns.gmap_db == '/tmp/db'
    assert ns.gmap_path == 'gmap'


def test_gmap_in_align_junction_choices():
    p = create_parser()
    ns = p.parse_args([
        'align', 'x.fastq', '--genome', 'g.fa',
        '--junction-aligners', 'gmap', '-o', '/tmp/out',
    ])
    assert 'gmap' in ns.junction_aligners


def test_run_gmap_missing_binary(monkeypatch, tmp_path):
    monkeypatch.setattr(multi_aligner.shutil, 'which', lambda _x: None)
    with pytest.raises(FileNotFoundError, match='gmap not found'):
        multi_aligner.run_gmap(
            reads_path=str(tmp_path / 'r.fastq'),
            genome_path=str(tmp_path / 'g.fa'),
            output_bam=str(tmp_path / 'out.bam'),
        )


def test_run_gmap_missing_db(monkeypatch, tmp_path):
    # Binary "found" so we reach the db check; db dir does not exist.
    monkeypatch.setattr(multi_aligner.shutil, 'which', lambda _x: '/usr/bin/gmap')
    with pytest.raises(FileNotFoundError, match='gmap_build'):
        multi_aligner.run_gmap(
            reads_path=str(tmp_path / 'r.fastq'),
            genome_path=str(tmp_path / 'g.fa'),
            output_bam=str(tmp_path / 'out.bam'),
            gmap_db=str(tmp_path / 'nonexistent_db'),
        )


def test_run_gmap_incomplete_db_rejected(monkeypatch, tmp_path):
    # An existing but empty db dir (no *.genomecomp) must be rejected, not used.
    monkeypatch.setattr(multi_aligner.shutil, 'which', lambda _x: '/usr/bin/gmap')
    empty_db = tmp_path / 'db'
    empty_db.mkdir()
    with pytest.raises(FileNotFoundError, match='missing or incomplete'):
        multi_aligner.run_gmap(
            reads_path=str(tmp_path / 'r.fastq'),
            genome_path=str(tmp_path / 'g.fa'),
            output_bam=str(tmp_path / 'out.bam'),
            gmap_db=str(empty_db),
        )
