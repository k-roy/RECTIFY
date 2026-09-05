"""Tests for `--junction-pool-max-intron-len` / `--junction-pool-min-anchor-bp`
on `rectify run-all` (dev/BUGS_TO_FIX.md C3, plus a coordinator follow-up that
threaded the sibling `--junction-pool-min-anchor-bp` flag through the identical
path).

Why this file exists: both flags existed only on `rectify consensus`
(consensus_command.py) — `run-all` and `align` had no way to reach either at
all, so the ONLY way to bound the non-annotated-junction candidate pool used by
the align-stage consensus selection's 5' soft-clip rescue was to hand-run
`align` then `consensus` as separate commands.

Where these flags actually land (traced by reading the code, not assumed): the
mechanism they configure (`pool_max_intron_len` / `pool_min_anchor_bp` in
`select_best_alignment`) is exclusively the ALIGN-stage multi-aligner consensus
selection that produces `<sample>.multialigned.bam` -- `run-all`'s stages.py calls
`align_command.run_align`, which calls `run_consensus_selection` itself
(align_command.py, "Run consensus selection" block) when >1 aligner succeeds and
`--no-consensus` was not passed. The LATER "merge -> consensus" step of the
correct-first pipeline (`write_corrected_consensus_bam`) has no candidate-pool
concept at all and cannot be what either flag means. So the plumbing is a
two-hop chain and both hops are tested here, for both flags:

  run_command.py (parser) -> single_sample.py -> stages.py::_run_alignment
      (Namespace attrs `junction_pool_max_intron_len` /
      `junction_pool_min_anchor_bp`, mirroring how
      resolver_acceptor_classes/resolver_atac are threaded)
    -> align_command.py::run_align reads them via getattr(...) and passes them
       as the `pool_max_intron_len` / `pool_min_anchor_bp` kwargs to
       run_consensus_selection.

`align_command.py` also gained both flags on its own parser (one `add_argument`
block each) so `rectify align` standalone can set them too -- this is NOT
required for run-all (run_align reads via getattr(..., 0) regardless of
whether align's own parser defines the flag), it's the brief-sanctioned "same
one-line plumbing" extra.

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


# ── the sibling flag, --junction-pool-min-anchor-bp, threaded the same way ──
# (coordinator follow-up: identical path, identical test shapes as above)

def test_run_all_exposes_the_min_anchor_flag():
    opts = [o for a in _run_all_parser()._actions for o in a.option_strings]
    assert '--junction-pool-min-anchor-bp' in opts, (
        "run-all must offer --junction-pool-min-anchor-bp (was consensus-only)"
    )


def test_run_all_min_anchor_default_matches_consensus_default():
    run_default = _run_all_parser().get_default('junction_pool_min_anchor_bp')
    consensus_default = _consensus_parser().get_default('junction_pool_min_anchor_bp')
    assert run_default == consensus_default
    assert run_default == 0, "consensus_command.py's own default is 0 (off)"


def test_run_all_parses_the_min_anchor_value():
    p = _run_all_parser()
    args = p.parse_args(
        ['in.bam', '-o', 'out', '--junction-pool-min-anchor-bp', '8']
    )
    assert args.junction_pool_min_anchor_bp == 8


def test_align_also_exposes_the_min_anchor_flag_with_the_same_default():
    opts = [o for a in _align_parser()._actions for o in a.option_strings]
    assert '--junction-pool-min-anchor-bp' in opts
    align_default = _align_parser().get_default('junction_pool_min_anchor_bp')
    consensus_default = _consensus_parser().get_default('junction_pool_min_anchor_bp')
    assert align_default == consensus_default


# ── hop 0: single_sample.py -> both flags reach BOTH _run_alignment call ────
# sites. Closes a gap the coordinator named after 854db78: the parser and
# hop-1/hop-2 tests above exercise stages.py and align_command.py directly,
# so none of them would notice if single_sample.py itself dropped a flag
# before ever calling _run_alignment. single_sample.py has TWO call sites
# (_process_one_sample for manifest/multi-sample mode, _run_single_sample
# for the single-sample CLI entry — see LOG.md's C3 multi_sample-coverage
# checkpoint) and each must pass BOTH flags, or the flag parses on run-all
# and is silently dropped on the floor -- the exact defect class this whole
# file exists to catch, just one hop earlier than hop 1 tests for.
#
# House pattern (test_run_all_protocol_flags.py::test_correct_stage_receives_ont_cdna)
# checks a single call site with `'X' in src` / `'X=' in src`. That is not
# enough here: with two call sites, deleting the kwarg line from ONE of them
# would still leave the substring present once, and a bare `in` check would
# not notice. Count occurrences instead so a single deletion at EITHER site
# fails the assertion.

def test_single_sample_both_call_sites_pass_both_flags():
    """Both _run_alignment call sites in single_sample.py must pass BOTH
    junction-pool flags through as getattr(args, '<flag>', 0) kwargs.

    Without this, a flag parses on run-all's own parser (proven above) and
    reaches stages.py._run_alignment's default value ONLY by accident of
    single_sample.py's own default matching -- an attacker-free but very
    real failure mode: someone edits one call site (e.g. adding a new
    kwarg) and forgets the other, and no test catches it because the
    parser tests never call single_sample.py, and hop 1/hop 2 tests call
    _run_alignment / run_align directly, bypassing single_sample.py
    entirely.
    """
    import inspect
    from rectify.core.commands.run import single_sample

    src = inspect.getsource(single_sample)

    for flag_name in ('junction_pool_max_intron_len', 'junction_pool_min_anchor_bp'):
        # Assigned (not merely mentioned) as a kwarg at both call sites.
        assign_count = src.count(f'{flag_name}=getattr(')
        assert assign_count == 2, (
            f"{flag_name}=getattr(...) must appear at BOTH _run_alignment "
            f"call sites in single_sample.py (_process_one_sample and "
            f"_run_single_sample); found {assign_count} occurrence(s) -- "
            f"a call site is dropping this flag."
        )
        # And it must read the MATCHING args attribute at both sites --
        # guards a copy-paste that keeps the right kwarg name but reads the
        # wrong (or misspelled) source attribute.
        read_count = src.count(f"args, '{flag_name}', 0")
        assert read_count == 2, (
            f"getattr(args, '{flag_name}', 0) must appear at BOTH call "
            f"sites; found {read_count} occurrence(s)."
        )


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


def test_align_stage_passes_the_min_anchor_value_through(tmp_path, monkeypatch):
    args = _captured_align_args(
        tmp_path, monkeypatch, junction_pool_min_anchor_bp=8,
    )
    assert args.junction_pool_min_anchor_bp == 8


def test_align_stage_min_anchor_default_is_off(tmp_path, monkeypatch):
    args = _captured_align_args(tmp_path, monkeypatch)
    assert args.junction_pool_min_anchor_bp == 0


# ── hop 2: align_command.run_align -> run_consensus_selection's kwargs ──────

def _build_align_args(tmp_path, prefix='sampleX', **overrides):
    """The Namespace run_align needs to reach its "Run consensus selection"
    block with two pre-seeded per-aligner BAMs under --trust-existing-bams
    (see test_run_align_passes_the_value_to_run_consensus_selection's
    docstring for why this avoids real aligner binaries / samtools).
    `overrides` lets each test set just the field(s) it cares about."""
    reads = tmp_path / 'reads.fastq'
    reads.touch()
    genome_path = tmp_path / 'genome.fa'
    genome_path.write_text('>chr1\n' + 'ACGT' * 50 + '\n')

    for aligner in ('minimap2', 'gapmm2'):
        bam_path = tmp_path / f'{prefix}.{aligner}.bam'
        if not bam_path.exists():
            bam_path.write_bytes(b'X' * 2001)

    base = dict(
        reads=reads, genome=genome_path, output_dir=tmp_path, annotation=None,
        threads=1, aligners=['minimap2', 'gapmm2'], short_read=False,
        read2=None, read_length=150, junction_aligners=[], max_intron=5000,
        resolver_acceptor_classes='canonical', resolver_atac=False,
        junction_pool_max_intron_len=0, junction_pool_min_anchor_bp=0,
        no_consensus=False, chimeric_consensus=False, junc_bonus=9,
        junc_bed=None, parallel_aligners=False, minimap2_path='minimap2',
        mapPacBio_path='mapPacBio.sh', gapmm2_path='gapmm2',
        ultra_path='uLTRA', desalt_path='deSALT', gmap_path='gmap',
        gmap_db=None, mapPacBio_chunks=1, mapPacBio_chunk_idx=None,
        prefix=prefix, keep_sam=False, sort=True, index=True, verbose=False,
        checkpoint_dir=None, trust_existing_bams=True,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _run_align_and_capture_consensus_kwargs(args):
    """Drives the REAL (unmocked) run_align; only check_aligner_available and
    run_consensus_selection itself are mocked (both at their ORIGIN module,
    per the function-local-import trap noted in _captured_align_args above).
    Returns the kwargs run_consensus_selection was called with."""
    captured = {}

    def _fake_run_consensus_selection(**kwargs):
        captured['kwargs'] = kwargs
        # Deliberately incomplete: run_align reads stats['consensus_high'] /
        # stats['5prime_rescued'] right after this call for logging, so an
        # empty dict makes it fail with a KeyError that run_align's own
        # try/except turns into rc=1. That's fine -- the kwargs are already
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
    return captured['kwargs']


def test_run_align_passes_the_value_to_run_consensus_selection(tmp_path):
    """The load-bearing fact: everything threaded onto args upstream is inert
    unless run_align actually forwards it into the pool_max_intron_len kwarg.
    """
    args = _build_align_args(tmp_path, junction_pool_max_intron_len=4321)
    kwargs = _run_align_and_capture_consensus_kwargs(args)
    assert kwargs['pool_max_intron_len'] == 4321


def test_run_align_passes_the_min_anchor_value_to_run_consensus_selection(tmp_path):
    """Same load-bearing fact for the sibling flag's pool_min_anchor_bp kwarg."""
    args = _build_align_args(tmp_path, junction_pool_min_anchor_bp=8)
    kwargs = _run_align_and_capture_consensus_kwargs(args)
    assert kwargs['pool_min_anchor_bp'] == 8


def test_both_flags_reach_run_consensus_selection_independently(tmp_path):
    """Targets the risk a SECOND flag on the same path introduces: a
    copy-paste error wiring one flag's CLI value into the other's kwarg (or
    into both). Distinct, asymmetric values on each so a swap would be
    caught by either assertion."""
    args = _build_align_args(
        tmp_path,
        junction_pool_max_intron_len=3000,
        junction_pool_min_anchor_bp=8,
    )
    kwargs = _run_align_and_capture_consensus_kwargs(args)
    assert kwargs['pool_max_intron_len'] == 3000
    assert kwargs['pool_min_anchor_bp'] == 8
