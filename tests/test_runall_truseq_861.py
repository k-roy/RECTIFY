"""861 — regressions found by the H2 TruSeq/COMPASS smoke (planning/861).

Two defects on master fd2e2d2:
  1. `run-all --chunked-alignment` died with
     `TypeError: generate_alignment_scripts() got an unexpected keyword argument
     'write_per_aligner_corrected_bams'` — on EVERY datatype (the call at
     run/chunked_batch.py:221 is unconditional).
  2. A COMPASS arm launched against a MISSING prebuilt index did not degrade to
     DROPPED-ALIGNER: `hisat2`'s Perl wrapper orphaned a grandchild that kept the
     captured stdout/stderr pipes open, so `subprocess.run(capture_output=True)`
     blocked for ALIGNER_TIMEOUT (6 h) at 0 % CPU.  The pre-flight check turns
     that into a clean FileNotFoundError before the binary is launched.
"""
import inspect
from pathlib import Path

import pytest

from rectify.core.commands import split_command
from rectify.core.align import multi_aligner as ma


# ── 1. the chunked-alignment TypeError ──────────────────────────────────────

def test_generate_alignment_scripts_accepts_write_per_aligner_corrected_bams():
    sig = inspect.signature(split_command.generate_alignment_scripts)
    assert 'write_per_aligner_corrected_bams' in sig.parameters
    assert sig.parameters['write_per_aligner_corrected_bams'].default is False


def test_chunked_batch_call_kwargs_are_all_accepted():
    """Every kwarg run/chunked_batch.py passes must exist on the callee."""
    src = Path(inspect.getsourcefile(
        __import__('rectify.core.commands.run.chunked_batch', fromlist=['x'])
    )).read_text()
    start = src.index('generate_alignment_scripts(')
    depth, i = 0, start + len('generate_alignment_scripts')
    while True:
        if src[i] == '(':
            depth += 1
        elif src[i] == ')':
            depth -= 1
            if depth == 0:
                break
        i += 1
    call = src[start:i]
    accepted = set(inspect.signature(split_command.generate_alignment_scripts).parameters)
    import re
    passed = set(re.findall(r'^\s*(\w+)\s*=', call, flags=re.M))
    unknown = {k for k in passed if k not in accepted}
    assert not unknown, f"chunked_batch passes kwargs the callee rejects: {unknown}"


def test_generate_alignment_scripts_forwards_the_flag(monkeypatch, tmp_path):
    seen = {}

    def fake(args, n_chunks, sample_prefix, **kw):
        seen.update(kw)

    monkeypatch.setattr(split_command, '_generate_scripts', fake)
    split_command.generate_alignment_scripts(
        n_chunks=2, sample_prefix='s', output_dir=tmp_path,
        genome=None, annotation=None, python_path='python', rectify_src='.',
        write_per_aligner_corrected_bams=True,
    )
    assert seen.get('write_per_aligner_corrected_bams') is True


# ── 2. the missing-index pre-flight ─────────────────────────────────────────

def test_star_missing_index_raises_before_launch(tmp_path):
    with pytest.raises(FileNotFoundError, match='COMPASS index missing'):
        ma._require_compass_index('STAR_default', 'star', tmp_path / 'nope')


def test_star_present_index_passes(tmp_path):
    d = tmp_path / 'STAR_annotated_75_bp_SJDB_index'
    d.mkdir()
    (d / 'SAindex').write_text('x')
    (d / 'genomeParameters.txt').write_text('x')
    ma._require_compass_index('STAR_default', 'star', d)  # must not raise


def test_star_half_built_index_still_raises(tmp_path):
    d = tmp_path / 'STAR_annotated_75_bp_SJDB_index'
    d.mkdir()
    (d / 'SAindex').write_text('x')          # genomeParameters.txt absent
    with pytest.raises(FileNotFoundError):
        ma._require_compass_index('STAR_default', 'star', d)


@pytest.mark.parametrize('suffix', ['.1.ht2', '.1.ht2l'])
def test_hisat2_index_accepts_both_ht2_flavours(tmp_path, suffix):
    p = tmp_path / 'gv'
    (tmp_path / f'gv{suffix}').write_text('x')
    ma._require_compass_index('HISAT2_default', 'hisat2', p)


def test_hisat2_missing_index_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match='HISAT2 index prefix'):
        ma._require_compass_index('HISAT2_default', 'hisat2', tmp_path / 'gv')


@pytest.mark.parametrize('ext', ['.nin', '.nsq', '.nal'])
def test_blast_db_variants(tmp_path, ext):
    (tmp_path / f'gv{ext}').write_text('x')
    ma._require_compass_index('magicblast', 'blast', tmp_path / 'gv')


def test_blast_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match='BLAST database'):
        ma._require_compass_index('magicblast', 'blast', tmp_path / 'gv')


def test_gsnap_empty_dir_raises(tmp_path):
    d = tmp_path / 'GSNAP' / 'gv'
    d.mkdir(parents=True)
    with pytest.raises(FileNotFoundError, match='GSNAP genome dir'):
        ma._require_compass_index('gsnap', 'gsnap', d)


def test_gsnap_populated_dir_passes(tmp_path):
    d = tmp_path / 'GSNAP' / 'gv'
    d.mkdir(parents=True)
    (d / 'gv.ref153positions').write_text('x')
    ma._require_compass_index('gsnap', 'gsnap', d)


def test_every_compass_wrapper_preflights_its_index():
    """A wrapper that skips the pre-flight can still deadlock the whole run."""
    for fn in (ma.run_star, ma.run_hisat2, ma.run_magicblast, ma.run_gsnap):
        src = inspect.getsource(fn)
        assert '_require_compass_index(' in src, f"{fn.__name__} has no index pre-flight"
        assert src.index('_require_compass_index(') < src.index('subprocess.run('), \
            f"{fn.__name__} pre-flights AFTER launching the binary"
