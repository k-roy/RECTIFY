#!/usr/bin/env python3
"""Paired-end short-read COMPASS panel wiring for `rectify run-all`.

`rectify align --short-read --read2 --aligners all` has always expanded to the
7-member COMPASS panel; `run-all` had no --read2 at all, so paired short-read
runs had to drop down to `align` + `aggregate` by hand. These tests pin the
run-all integration.

The subtlety worth guarding: run-all does NOT hand align the literal "all"
sentinel — `run/stages.py::_run_alignment` passes an EXPLICIT aligner list, so
align_command's own "all"-expansion never fires from the run-all path. Selecting
the COMPASS panel therefore has to happen in _run_alignment too; a naive
"plumb --read2 through" would silently run bbmap+bwa *paired* instead. The panel
tests below fail on exactly that mistake.

No aligner binaries are needed: run_align is intercepted, so nothing shells out.
"""

import argparse
from pathlib import Path

import pytest

from rectify.core.commands import align_command
from rectify.core.commands.align_command import COMPASS_PE_ALIGNERS
from rectify.core.commands.run_command import (
    _validate_paired_end_args,
    create_run_parser,
)
from rectify.core.commands.run.stages import _run_alignment


EXPECTED_COMPASS_PANEL = [
    'bbmap',
    'STAR_default',
    'STAR_noncanonical',
    'HISAT2_default',
    'HISAT2_noncanonical',
    'magicblast',
    'gsnap',
]


def _run_all_parser():
    sub = argparse.ArgumentParser().add_subparsers(dest='command')
    create_run_parser(sub)
    return sub.choices['run-all']


# ── The panel constant (single source of truth shared with `align`) ──────────

def test_compass_panel_constant_is_the_seven_member_set():
    """6 tools; STAR + HISAT2 each contribute a canonical and non-canonical mode."""
    assert COMPASS_PE_ALIGNERS == EXPECTED_COMPASS_PANEL
    assert len(COMPASS_PE_ALIGNERS) == 7


def test_align_all_expansion_uses_the_shared_constant():
    """`align`'s --short-read --read2 "all" expansion must stay tied to the constant.

    Guards the align/run-all pair against drifting apart — the whole point of
    the integration is that both select the identical panel.
    """
    import inspect
    src = inspect.getsource(align_command.run_align)
    assert 'COMPASS_PE_ALIGNERS' in src


# ── run-all parser ───────────────────────────────────────────────────────────

class TestRunAllPairedEndParser:
    def test_accepts_read2_and_read_length(self):
        args = _run_all_parser().parse_args([
            'R1.fastq.gz', '--short-read', '--read2', 'R2.fastq.gz',
            '--read-length', '100', '-o', 'out/',
        ])
        assert args.read2 == Path('R2.fastq.gz')
        assert args.read_length == 100
        assert args.short_read is True

    def test_positional_r1_may_trail_read2(self):
        """The documented invocation puts R1 AFTER --read2:

            rectify run-all --short-read --read2 R2.fastq.gz R1.fastq.gz -o results/

        --read2 is single-valued so it can't swallow R1, but this interleaving of
        an optional and a nargs='?' positional is exactly where argparse surprises
        people — pin it.
        """
        args = _run_all_parser().parse_args([
            '--short-read', '--read2', 'R2.fastq.gz', 'R1.fastq.gz', '-o', 'results/',
        ])
        assert args.input == Path('R1.fastq.gz')
        assert args.read2 == Path('R2.fastq.gz')
        assert _validate_paired_end_args(args) is None

    def test_dash_2_is_an_alias_for_read2(self):
        args = _run_all_parser().parse_args([
            'R1.fastq.gz', '--short-read', '-2', 'R2.fastq.gz', '-o', 'out/',
        ])
        assert args.read2 == Path('R2.fastq.gz')

    def test_defaults_are_backwards_compatible(self):
        """Existing invocations keep working: read2 absent, read_length 150."""
        args = _run_all_parser().parse_args(['s.fastq.gz', '-o', 'out/'])
        assert args.read2 is None
        assert args.read_length == 150

    def test_base_aligners_accepts_compass_members(self):
        args = _run_all_parser().parse_args([
            'R1.fastq.gz', '--short-read', '-2', 'R2.fastq.gz', '-o', 'out/',
            '--base-aligners', 'STAR_default', 'gsnap',
        ])
        assert args.base_aligners == ['STAR_default', 'gsnap']

    def test_short_read_help_describes_both_panels(self):
        """The help used to claim bbmap+bwa unconditionally, which was wrong."""
        action = next(
            a for a in _run_all_parser()._actions if a.dest == 'short_read'
        )
        assert 'COMPASS' in action.help
        assert '--read2' in action.help
        assert 'bbmap + bwa' in action.help


# ── Panel selection (the part a naive integration gets wrong) ─────────────────

class TestRunAllPanelSelection:
    """_run_alignment must pick the COMPASS panel itself, not defer to align."""

    @staticmethod
    def _captured_align_args(tmp_path, monkeypatch, **kwargs):
        captured = {}

        class _Stop(Exception):
            pass

        def _fake_run_align(args):
            captured['args'] = args
            raise _Stop()

        # _run_alignment imports run_align from the module at call time.
        monkeypatch.setattr(align_command, 'run_align', _fake_run_align)

        with pytest.raises(_Stop):
            _run_alignment(
                input_path=Path('R1.fastq.gz'),
                sample_id='s',
                sample_output_dir=tmp_path,
                genome_path=Path('genome.fa'),
                annotation_path=None,
                threads=4,
                **kwargs,
            )
        return captured['args']

    def test_paired_short_read_selects_the_full_compass_panel(self, tmp_path, monkeypatch):
        args = self._captured_align_args(
            tmp_path, monkeypatch,
            short_read=True, read2=Path('R2.fastq.gz'),
        )
        assert list(args.aligners) == EXPECTED_COMPASS_PANEL

    def test_paired_short_read_plumbs_read2_and_read_length(self, tmp_path, monkeypatch):
        """STAR needs read_length for --sjdbOverhang; R2 must reach the wrappers."""
        args = self._captured_align_args(
            tmp_path, monkeypatch,
            short_read=True, read2=Path('R2.fastq.gz'), read_length=100,
        )
        assert args.read2 == Path('R2.fastq.gz')
        assert args.read_length == 100

    def test_single_end_truseq_selects_the_compass_se_panel(self, tmp_path, monkeypatch):
        """Plain --short-read (TruSeq-style RNA-seq) = splice-aware COMPASS
        SE subset (Kevin 2026-08-17: QuantSeq vs TruSeq panels split)."""
        args = self._captured_align_args(tmp_path, monkeypatch, short_read=True)
        assert list(args.aligners) == [
            'bbmap', 'STAR_default', 'STAR_noncanonical',
            'HISAT2_default', 'HISAT2_noncanonical',
        ]
        assert args.read2 is None

    def test_single_end_dt_primed_selects_the_3prime_panel(self, tmp_path, monkeypatch):
        """--short-read + dT-primed (QuantSeq REV) keeps the 3'-end panel."""
        args = self._captured_align_args(
            tmp_path, monkeypatch, short_read=True, dt_primed_cdna=True)
        assert list(args.aligners) == ['bbmap', 'bwa']
        assert args.read2 is None

    def test_long_read_default_is_unchanged(self, tmp_path, monkeypatch):
        """The DRS/PCR-cDNA production base panel: minimap2 alone —
        mapPacBio/gapmm2 de-paneled 2026-08-17 (opt-in via --base-aligners)."""
        args = self._captured_align_args(tmp_path, monkeypatch)
        assert list(args.aligners) == ['minimap2']
        assert args.read2 is None
        assert args.read_length == 150

    def test_explicit_base_aligners_still_override_the_panel(self, tmp_path, monkeypatch):
        args = self._captured_align_args(
            tmp_path, monkeypatch,
            short_read=True, read2=Path('R2.fastq.gz'), base_aligners=['bbmap'],
        )
        assert list(args.aligners) == ['bbmap']


# ── Validation of nonsensical combinations ───────────────────────────────────

class TestPairedEndValidation:
    @staticmethod
    def _args(**kw):
        base = dict(
            read2=None, short_read=False, manifest=None,
            input=None, chunked_alignment=False,
        )
        base.update(kw)
        return argparse.Namespace(**base)

    def test_no_read2_is_always_valid(self):
        assert _validate_paired_end_args(self._args()) is None
        assert _validate_paired_end_args(
            self._args(input=Path('s.bam'), manifest=Path('m.tsv'))
        ) is None

    def test_valid_paired_end_combo_passes(self):
        err = _validate_paired_end_args(self._args(
            read2=Path('R2.fastq.gz'), short_read=True, input=Path('R1.fastq.gz'),
        ))
        assert err is None

    def test_read2_without_short_read_is_rejected(self):
        err = _validate_paired_end_args(self._args(
            read2=Path('R2.fastq.gz'), input=Path('R1.fastq.gz'),
        ))
        assert err is not None
        assert '--short-read' in err

    def test_read2_with_bam_input_is_rejected(self):
        err = _validate_paired_end_args(self._args(
            read2=Path('R2.fastq.gz'), short_read=True, input=Path('s.bam'),
        ))
        assert err is not None
        assert 'BAM' in err

    def test_read2_with_manifest_is_rejected(self):
        """The manifest schema (sample_id, path, condition) carries no mate column."""
        err = _validate_paired_end_args(self._args(
            read2=Path('R2.fastq.gz'), short_read=True, manifest=Path('m.tsv'),
        ))
        assert err is not None
        assert '--manifest' in err

    def test_read2_with_chunked_alignment_fails_loudly_with_the_alternative(self):
        """Must NOT silently drop R2 — and must name the supported chunked route."""
        err = _validate_paired_end_args(self._args(
            read2=Path('R2.fastq.gz'), short_read=True,
            input=Path('R1.fastq.gz'), chunked_alignment=True,
        ))
        assert err is not None
        assert '--chunked-alignment' in err
        assert 'rectify split' in err
        assert '--generate-slurm' in err
