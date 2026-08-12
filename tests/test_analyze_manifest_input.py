"""`rectify analyze --manifest` must not require the unused positional input.

Reported by the 668-drs-arm session: a correct manifest-mode invocation died
with `error: the following arguments are required: input`, four minutes into a
cohort run. The positional is never read in manifest mode — argparse simply
made it mandatory.

The workaround in the tree is `run/multi_sample.py:43`, which passes the
manifest path as BOTH `input=` and `manifest=`. That must keep working: it is
what existing callers do.
"""

import argparse

import pytest

from rectify.core.commands.analyze_command import create_analyze_parser


def _parse(argv):
    root = argparse.ArgumentParser()
    sub = root.add_subparsers(dest='command')
    create_analyze_parser(sub)
    return root.parse_args(argv)


def test_manifest_mode_needs_no_positional(tmp_path):
    """The reported failure: this invocation used to be rejected by argparse."""
    args = _parse(['analyze', '--manifest', 'samples.tsv', '-o', str(tmp_path)])
    assert args.manifest == 'samples.tsv'
    assert args.input is None


def test_positional_mode_is_unchanged(tmp_path):
    args = _parse(['analyze', 'corrected.tsv', '-o', str(tmp_path)])
    assert args.input == 'corrected.tsv'
    assert args.manifest is None


def test_both_still_accepted(tmp_path):
    """run/multi_sample.py passes the manifest as both — must not break."""
    args = _parse(['analyze', 'samples.tsv', '--manifest', 'samples.tsv',
                   '-o', str(tmp_path)])
    assert args.input == 'samples.tsv' and args.manifest == 'samples.tsv'


def test_neither_is_rejected_with_a_message_naming_both_options(tmp_path, caplog):
    """Making the positional optional must not let a no-input run start."""
    import logging
    from rectify.core.commands.analyze_command import run_analyze

    args = _parse(['analyze', '-o', str(tmp_path)])
    caplog.set_level(logging.ERROR)
    assert run_analyze(args) == 2
    msg = caplog.text
    assert '--manifest' in msg
    assert 'positional' in msg
