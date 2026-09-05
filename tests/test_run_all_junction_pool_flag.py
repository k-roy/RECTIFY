"""Tests for `--junction-pool-max-intron-len` on `rectify run-all` (dev/BUGS_TO_FIX.md C3).

Why this file exists: the flag existed only on `rectify consensus`
(consensus_command.py) — `run-all` and `align` had no way to reach it at all, so
the ONLY way to bound the non-annotated-junction candidate pool used by the
align-stage consensus selection's 5' soft-clip rescue was to hand-run `align`
then `consensus` as separate commands.

Where this flag actually lands (traced by reading the code, not assumed): the
`--junction-pool-max-intron-len` mechanism (`pool_max_intron_len` in
`select_best_alignment`) is exclusively the ALIGN-stage multi-aligner consensus
selection that produces `<sample>.multialigned.bam` -- `run-all`'s stages.py calls
`align_command.run_align`, which calls `run_consensus_selection` itself
(align_command.py, "Run consensus selection" block) when >1 aligner succeeds and
`--no-consensus` was not passed. The LATER "merge -> consensus" step of the
correct-first pipeline (`write_corrected_consensus_bam`) has no candidate-pool
concept at all and cannot be what this flag means. So the plumbing is a two-hop
chain and both hops are tested here:

  run_command.py (parser) -> single_sample.py -> stages.py::_run_alignment
      (Namespace attr `junction_pool_max_intron_len`, mirroring how
      resolver_acceptor_classes/resolver_atac are threaded)
    -> align_command.py::run_align reads it via getattr(...) and passes it as
       the `pool_max_intron_len` kwarg to run_consensus_selection.

`align_command.py` also gained the same flag on its own parser (one `add_argument`
block) so `rectify align` standalone can set it too -- this is NOT required for
run-all (run_align reads via getattr(..., 0) regardless of whether align's own
parser defines the flag), it's the brief-sanctioned "same one-line plumbing"
extra.

No aligner binaries or real BAMs are needed for the second hop: the
`--trust-existing-bams` per-aligner-BAM checkpoint (already in align_command.py)
is used to short-circuit `_run_one_aligner` with two pre-seeded dummy files, the
same technique as test_run_all_compass_pe.py's run_align-interception style, one
level deeper (through run_align for real instead of intercepting it).
"""
import argparse
from pathlib import Path
from unittest import mock

import pytest

from rectify.core.commands import align_command
from rectify.core.commands.run.stages import _run_alignment


def _run_all_parser():
    from rectify.core.commands.run_command import create_run_parser
    sub = argparse.ArgumentParser(prog='rectify').add_subparsers(dest='command')
    create_run_parser(sub)
    return sub.choices['run-all']


def _align_parser():
    from rectify.core.commands.align_command import create_align_parser
    sub = argparse.ArgumentParser(prog='rectify').add_subparsers(dest='command')
    create_align_parser(sub)
    return sub.choices['align']


def _consensus_parser():
    from rectify.core.commands.consensus_command import create_consensus_parser
    sub = argparse.ArgumentParser(prog='rectify').add_subparsers(dest='command')
    create_consensus_parser(sub)
    return sub.choices['consensus']


# ── the flag exists on run-all, parses, and matches consensus's default ─────

def test_run_all_exposes_the_flag():
    opts = [o for a in _run_all_parser()._actions for o in a.option_strings]
    assert '--junction-pool-max-intron-len' in opts, (
        "run-all must offer --junction-pool-max-intron-len (was consensus-only)"
    )


def test_run_all_default_matches_consensus_default():
    """Compare parser to parser -- a hardcoded literal would silently drift
    the moment either side's default changes."""
    run_default = _run_all_parser().get_default('junction_pool_max_intron_len')
    consensus_default = _consensus_parser().get_default('junction_pool_max_intron_len')
    assert run_default == consensus_default
    assert run_default == 0, "consensus_command.py's own default is 0 (off)"


def test_run_all_parses_the_value():
    p = _run_all_parser()
    args = p.parse_args(
        ['in.bam', '-o', 'out', '--junction-pool-max-intron-len', '3000']
    )
    assert args.junction_pool_max_intron_len == 3000


def test_align_also_exposes_the_flag_with_the_same_default():
    """Elective (brief-sanctioned) extra: align_command.run_align is where the
    value actually reaches run_consensus_selection, so `rectify align`
    standalone gets the same knob for free as one more add_argument block."""
    opts = [o for a in _align_parser()._actions for o in a.option_strings]
    assert '--junction-pool-max-intron-len' in opts
    align_default = _align_parser().get_default('junction_pool_max_intron_len')
    consensus_default = _consensus_parser().get_default('junction_pool_max_intron_len')
    assert align_default == consensus_default


# ── hop 1: stages.py._run_alignment -> the args Namespace handed to run_align ─

class _Stop(Exception):
    """Deliberate short-circuit once run_align has been called and captured --
    mirrors TestRunAllPanelSelection._captured_align_args in
    test_run_all_compass_pe.py (same interception technique, new file so that
    gate file stays untouched)."""


def _captured_align_args(tmp_path, monkeypatch, **kwargs):
    captured = {}

    def _fake_run_align(args):
        captured['args'] = args
        raise _Stop()

    # _run_alignment imports run_align from the module at call time
    # ("from ..align_command import run_align" inside the function), so the
    # patch target is the ORIGIN module, not rectify.core.commands.run.stages.
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


def test_align_stage_passes_the_value_through(tmp_path, monkeypatch):
    args = _captured_align_args(
        tmp_path, monkeypatch, junction_pool_max_intron_len=1234,
    )
    assert args.junction_pool_max_intron_len == 1234


def test_align_stage_default_is_off(tmp_path, monkeypatch):
    """Omitting the kwarg (the pre-C3 call shape) must not crash and must
    default to 0 -- _run_alignment's own default, matching consensus's."""
    args = _captured_align_args(tmp_path, monkeypatch)
    assert args.junction_pool_max_intron_len == 0


# ── hop 2: align_command.run_align -> run_consensus_selection's kwarg ───────

def test_run_align_passes_the_value_to_run_consensus_selection(tmp_path):
    """The load-bearing fact: everything threaded onto args upstream is inert
    unless run_align actually forwards it into the pool_max_intron_len kwarg.

    Drives the REAL (unmocked) run_align with two pre-seeded per-aligner BAMs
    under --trust-existing-bams so _run_one_aligner's existing checkpoint
    short-circuits before shelling out to any real aligner binary or
    `samtools sort` -- only check_aligner_available and run_consensus_selection
    itself are mocked.
    """
    reads = tmp_path / 'reads.fastq'
    reads.touch()
    genome_path = tmp_path / 'genome.fa'
    genome_path.write_text('>chr1\n' + 'ACGT' * 50 + '\n')

    prefix = 'sampleX'
    for aligner in ('minimap2', 'gapmm2'):
        (tmp_path / f'{prefix}.{aligner}.bam').write_bytes(b'X' * 2001)

    args = argparse.Namespace(
        reads=reads, genome=genome_path, output_dir=tmp_path, annotation=None,
        threads=1, aligners=['minimap2', 'gapmm2'], short_read=False,
        read2=None, read_length=150, junction_aligners=[], max_intron=5000,
        resolver_acceptor_classes='canonical', resolver_atac=False,
        junction_pool_max_intron_len=4321, no_consensus=False,
        chimeric_consensus=False, junc_bonus=9, junc_bed=None,
        parallel_aligners=False, minimap2_path='minimap2',
        mapPacBio_path='mapPacBio.sh', gapmm2_path='gapmm2',
        ultra_path='uLTRA', desalt_path='deSALT', gmap_path='gmap',
        gmap_db=None, mapPacBio_chunks=1, mapPacBio_chunk_idx=None,
        prefix=prefix, keep_sam=False, sort=True, index=True, verbose=False,
        checkpoint_dir=None, trust_existing_bams=True,
    )

    captured = {}

    def _fake_run_consensus_selection(**kwargs):
        captured['kwargs'] = kwargs
        # Deliberately incomplete: run_align reads stats['consensus_high'] /
        # stats['5prime_rescued'] right after this call for logging, so an
        # empty dict makes it fail with a KeyError that run_align's own
        # try/except turns into rc=1. That's fine -- the kwarg is already
        # captured by the time that happens.
        return {}

    with mock.patch(
        'rectify.core.align.multi_aligner.check_aligner_available',
        return_value=True,
    ), mock.patch(
        'rectify.core.consensus.consensus.run_consensus_selection',
        side_effect=_fake_run_consensus_selection,
    ):
        align_command.run_align(args)

    assert 'kwargs' in captured, (
        "run_consensus_selection was never called -- the >1-successful-aligner "
        "consensus branch was not reached"
    )
    assert captured['kwargs']['pool_max_intron_len'] == 4321
