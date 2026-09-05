"""`run-all`'s Station B/C wiring for ISSUE-009 and ISSUE-013.

Both defects were the same shape: the fix landed on the `rectify triage` /
`rectify pool-gate` CLIs and `run-all` kept calling the old thing, so the
one-command path a lab member actually uses could not reach it. That is exactly
the failure mode `test_run_all_protocol_flags.py` exists for, and these are its
Station B/C siblings.

- **ISSUE-013.** `_run_station_bc` derived Station C's length pre-gate from
  `multi_aligner.derive_max_intron` — 2x the LONGEST annotated intron, clamped
  to 500,000. GENCODE v48 chr5's longest intron is 772,519 bp, so the bound
  saturated the clamp and no human junction could ever trip the pre-gate. The
  `pool-gate` CLI moved to `derive_pool_gate_max_intron` (2x p99.5 = 310,100 on
  the same annotation, 5,000 unchanged on yeast); `run-all` must follow.
- **ISSUE-009.** Station B's junction leg proposes with Module 2H, the engine
  that produced the corrected BAM, so it re-derives its own fixed point.
  `run-all` knows where this run's pre-correction BAM is; the user should not
  have to assemble paths.
"""
import argparse
import inspect
from pathlib import Path

import pytest

from rectify.core.commands.run import stages


def _run_all_parser():
    """The Station B/C flags are attached in `rectify.cli` AFTER
    `create_run_parser` (cli.py: `add_station_bc_args(subparsers.choices[...])`),
    so a locally-built run-all parser does not carry them — go through the real
    top-level parser, as test_station_c.py does."""
    from rectify.cli import create_parser
    top = create_parser()
    for action in top._actions:
        choices = getattr(action, 'choices', None)
        if isinstance(choices, dict) and 'run-all' in choices:
            return choices['run-all']
    raise AssertionError('run-all subparser not found')


BASE = ['in.bam', '-o', 'out']


# ── ISSUE-009: the flag exists, is opt-in, and reaches the triage namespace ──

def test_run_all_exposes_triage_original_bams():
    opts = [o for a in _run_all_parser()._actions for o in a.option_strings]
    assert '--triage-original-bams' in opts
    assert '--triage-correction-regression-ratio' in opts


def test_triage_original_bams_is_opt_in_and_has_three_modes():
    p = _run_all_parser()
    # absent -> None: the leg is inert and the run is unchanged
    assert p.parse_args(BASE).triage_original_bams is None
    # present with no paths -> [] : auto-discover (falsy, so callers MUST test
    # `is not None` rather than truthiness)
    assert p.parse_args(BASE + ['--triage-original-bams']).triage_original_bams == []
    # present with paths -> explicit override
    assert p.parse_args(
        BASE + ['--triage-original-bams', 'a.bam', 'b.bam']
    ).triage_original_bams == ['a.bam', 'b.bam']


def test_correction_regression_ratio_defaults_to_one():
    p = _run_all_parser()
    assert p.parse_args(BASE).triage_correction_regression_ratio == 1.0
    assert p.parse_args(
        BASE + ['--triage-correction-regression-ratio', '0.5']
    ).triage_correction_regression_ratio == 0.5


def test_stages_threads_original_bams_into_the_triage_namespace():
    """The flag must be ASSIGNED into the namespace handed to run_triage, not
    merely read — the invisible failure mode protocol_flags guards against."""
    src = inspect.getsource(stages)
    assert 'original_bams=_orig_bams' in src, (
        'stages.py reads the flag but never passes it to run_triage')
    assert "getattr(args, 'triage_original_bams', None) is not None" in src, (
        'the auto-discover mode is an EMPTY LIST, so a truthiness test would '
        'silently disable it')
    assert 'max_correction_regression_ratio=' in src


# ── ISSUE-009: what auto-discovery resolves to ──────────────────────────────

def _touch(p: Path, size: int = 64):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(b'\0' * size)
    return p


def test_auto_discovery_prefers_the_multialigned_bam(tmp_path):
    """<sample>.multialigned.bam is the merged pre-correction artifact that
    `correct` actually consumed — the exact "before" of this run's "after", and
    the shape the fix was measured on. Per-aligner BAMs are one step further
    back, so they are the fallback."""
    _touch(tmp_path / 's1.multialigned.bam')
    _touch(tmp_path / 's1.minimap2.bam')
    args = argparse.Namespace(triage_original_bams=[])
    got = stages._pre_correction_bams(tmp_path, 's1', args)
    assert got == [str(tmp_path / 's1.multialigned.bam')]


def test_auto_discovery_falls_back_to_the_per_aligner_bams(tmp_path):
    _touch(tmp_path / 's1.minimap2.bam')
    _touch(tmp_path / 's1.gapmm2.bam')
    args = argparse.Namespace(triage_original_bams=[])
    got = stages._pre_correction_bams(tmp_path, 's1', args)
    assert len(got) == 2
    assert all(Path(p).name.startswith('s1.') for p in got)
    assert not any('multialigned' in p for p in got)


def test_explicit_paths_override_discovery(tmp_path):
    _touch(tmp_path / 's1.multialigned.bam')
    args = argparse.Namespace(triage_original_bams=['/elsewhere/orig.bam'])
    assert stages._pre_correction_bams(tmp_path, 's1', args) == [
        '/elsewhere/orig.bam']


def test_nothing_on_disk_yields_an_empty_list(tmp_path):
    """Which leaves the leg inert rather than crashing a fail-soft post-stage."""
    args = argparse.Namespace(triage_original_bams=[])
    assert stages._pre_correction_bams(tmp_path, 's1', args) == []


# ── ISSUE-013: Station C's pre-gate in run-all ──────────────────────────────

def test_stages_uses_the_quantile_pre_gate_not_the_aligner_bound():
    src = inspect.getsource(stages._run_station_bc)
    assert 'derive_pool_gate_max_intron' in src, (
        'run-all must derive Station C\'s pre-gate the way the pool-gate CLI '
        'does; derive_max_intron saturates its 500,000 clamp on human and the '
        'pre-gate then never fires')
    assert 'derive_max_intron' not in src.replace(
        'derive_pool_gate_max_intron', ''), (
        'the aligner bound must not be used for the Station C pre-gate')
    assert 'length pre-gate' in src, 'the resolved bound must be logged'


@pytest.mark.parametrize('annotation,expected', [('yeast', 5000)])
def test_the_wired_bound_is_unchanged_on_yeast(annotation, expected):
    """The regression that matters most: whatever run-all now derives, yeast
    must still get exactly 5,000."""
    import rectify
    from rectify.core.consensus.station_c import derive_pool_gate_max_intron
    gff = (Path(rectify.__file__).parent / 'data' / 'genomes'
           / 'saccharomyces_cerevisiae'
           / 'saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz')
    if not gff.exists():
        pytest.skip('bundled yeast annotation absent')
    assert derive_pool_gate_max_intron(gff) == expected
