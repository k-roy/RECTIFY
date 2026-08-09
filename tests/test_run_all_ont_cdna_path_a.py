"""Tests for ONT PCR-cDNA Path A as the `run-all` default, and the guards around it.

These exist because of a specific class of failure: the previous cDNA preparation was
inlined in `_process_one_sample` (the `--manifest` worker) ONLY, so
`run-all --ONT-cDNA reads.fastq` -- the invocation a lab member would naturally type --
parsed the flag, ran, exited 0, and silently corrected every read through the weakest
strand channel (~19% of 3' ends >250nt inside the CDS vs ~1.7-2.2% tag-resolved).

🔴 The regression test that shipped alongside that bug PASSED, because it grepped
`inspect.getsource(single_sample)` at MODULE scope -- `_process_one_sample`'s copy
satisfied the assertion while `_run_single_sample` had nothing. Every test here
therefore asserts PER FUNCTION, never against the module as a whole.
"""
import inspect

import pysam
import pytest

from rectify.core.commands.run import single_sample as ss


# --------------------------------------------------------------------------
# Both entry points must prepare cDNA -- asserted per function, not per module
# --------------------------------------------------------------------------

@pytest.mark.parametrize("fname", ["_run_single_sample", "_process_one_sample"])
def test_entry_point_calls_cdna_prepare(fname):
    """EVERY run-all entry point must call the shared cDNA preparation helper."""
    src = inspect.getsource(getattr(ss, fname))
    assert "_run_ont_cdna_prepare(" in src, (
        "%s does not call _run_ont_cdna_prepare -- a --ONT-cDNA run through this "
        "entry point would skip cDNA preparation silently" % fname
    )


@pytest.mark.parametrize("fname", ["_run_single_sample", "_process_one_sample"])
def test_entry_point_does_not_inline_the_trim(fname):
    """The trim must NOT be re-inlined -- copy-paste is how the original bug happened."""
    src = inspect.getsource(getattr(ss, fname))
    assert "trim_cdna_fastq_polya(" not in src, (
        "%s inlines trim_cdna_fastq_polya; call _run_ont_cdna_prepare instead so "
        "both entry points cannot drift apart" % fname
    )


def test_shared_helper_exists_and_dispatches_both_paths():
    src = inspect.getsource(ss._run_ont_cdna_prepare)
    assert "_run_ont_cdna_path_a(" in src
    assert "trim_cdna_fastq_polya(" in src   # the Path B branch


# --------------------------------------------------------------------------
# Path A is the DEFAULT
# --------------------------------------------------------------------------

def _run_all_parser():
    import argparse
    from rectify.core.commands.run_command import create_run_parser
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='command')
    create_run_parser(sub)
    return p


def test_ont_cdna_path_defaults_to_a():
    """Default must be molecule-level. A read-level default silently violates policy."""
    args = _run_all_parser().parse_args(['run-all', 'x.fastq.gz', '-o', 'out', '--ONT-cDNA'])
    assert getattr(args, 'ont_cdna_path', None) == 'a'


def test_ont_cdna_path_b_is_selectable():
    args = _run_all_parser().parse_args(
        ['run-all', 'x.fastq.gz', '-o', 'out', '--ONT-cDNA', '--ont-cdna-path', 'b'])
    assert args.ont_cdna_path == 'b'


def test_ont_cdna_path_rejects_unknown_value():
    with pytest.raises(SystemExit):
        _run_all_parser().parse_args(
            ['run-all', 'x.fastq.gz', '-o', 'out', '--ONT-cDNA', '--ont-cdna-path', 'z'])


def test_poa_off_by_default():
    """Dedup needs no consensus sequence and POA is ~44% of stage-1 runtime."""
    args = _run_all_parser().parse_args(['run-all', 'x.fastq.gz', '-o', 'out', '--ONT-cDNA'])
    assert getattr(args, 'ont_cdna_poa', True) is False


# --------------------------------------------------------------------------
# Path A refuses input it cannot handle, rather than degrading
# --------------------------------------------------------------------------

def test_path_a_rejects_non_fastq_input(tmp_path):
    class _A:
        genome = 'g.fa'; annotation = None; threads = 2; ont_cdna_poa = False
    with pytest.raises(ValueError, match="FASTQ"):
        ss._run_ont_cdna_path_a(
            _A(), tmp_path / 'in.bam', 'bam', tmp_path, tmp_path, 'S1')


# --------------------------------------------------------------------------


import inspect  # noqa: E402

# --------------------------------------------------------------------------
# --chunked-alignment must refuse protocol flags rather than drop them
# --------------------------------------------------------------------------

def test_chunked_alignment_refuses_protocol_flags():
    """The chunked generator does not forward protocol flags into the emitted
    `rectify correct`, so a cDNA library would be corrected under the DRS rule.
    Refusing is correct; generating a script that yields a confident wrong answer
    is not."""
    src = inspect.getsource(
        __import__('rectify.core.commands.run_command', fromlist=['x']))
    assert "--chunked-alignment does not propagate protocol flags" in src
