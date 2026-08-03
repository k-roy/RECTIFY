"""Tests for `rectify run-all`'s protocol flags — in particular `--ONT-cDNA`.

Why this file exists: until 2026-08-03 `--ONT-cDNA` was reachable ONLY from
`rectify correct`. `run-all` and `split` offered just `--drs` / `--dT-primed-cDNA`,
so the one-command path a lab member would actually use could not select the ONT
PCR-cDNA protocol at all — and `--dT-primed-cDNA` (QuantSeq REV) applies a *fixed*
antisense strand rule, which for a both-orientations PCR-cDNA library puts the RNA
3' end at the wrong terminus for roughly half the reads. That is exactly the defect
`planning/541` fixed inside `correct`; leaving `run-all` unable to ask for it made
the fix unreachable in practice.

These tests lock down (a) the flag exists and reaches `args.ONT_cDNA`, (b) it is
threaded into the correct-stage namespace, and (c) two protocols are rejected
loudly rather than one silently winning.
"""
import argparse

import pytest

from rectify.core.commands.run_command import (
    _validate_protocol_flags,
    _PROTOCOL_FLAGS,
)


def _ns(**kw):
    base = dict(drs=False, dT_primed_cDNA=False, ONT_cDNA=False, short_read=False)
    base.update(kw)
    return argparse.Namespace(**base)


# ── the flag exists on run-all and lands on args.ONT_cDNA ───────────────────
def _run_all_parser():
    # create_run_parser() wires the subcommand in and returns None, so recover
    # the parser from the subparsers action's choices.
    from rectify.core.commands.run_command import create_run_parser
    sub = argparse.ArgumentParser(prog='rectify').add_subparsers(dest='command')
    create_run_parser(sub)
    return sub.choices['run-all']


def test_run_all_exposes_ont_cdna():
    p = _run_all_parser()
    opts = [o for a in p._actions for o in a.option_strings]
    assert '--ONT-cDNA' in opts, "run-all must offer --ONT-cDNA (was correct-only)"


def test_run_all_ont_cdna_defaults_off_and_parses():
    p = _run_all_parser()
    base = ['in.bam', '-o', 'out']
    assert p.parse_args(base).ONT_cDNA is False
    assert p.parse_args(base + ['--ONT-cDNA']).ONT_cDNA is True


def test_run_all_still_offers_the_other_protocols():
    # Guard against a rename silently dropping one.
    opts = [o for a in _run_all_parser()._actions for o in a.option_strings]
    for flag, _dest in _PROTOCOL_FLAGS:
        assert flag in opts, f"run-all lost {flag}"


# ── the flag is threaded into the correct stage ─────────────────────────────
def test_correct_stage_receives_ont_cdna():
    """stages.py must put ONT_cDNA on the namespace it hands to correct.

    Without this the flag parses and is then dropped on the floor — the failure
    mode that is invisible until someone checks a strand column.
    """
    import inspect
    from rectify.core.commands.run import stages
    src = inspect.getsource(stages)
    assert 'ONT_cDNA' in src, "stages.py never mentions ONT_cDNA"
    # It must be *assigned into* the correct_args namespace, not merely read.
    assert 'ONT_cDNA=' in src, "ONT_cDNA is read but never passed to correct_args"


# ── mutually exclusive protocols ───────────────────────────────────────────
def test_no_protocol_is_valid():
    assert _validate_protocol_flags(_ns()) is None


@pytest.mark.parametrize("dest", [d for _f, d in _PROTOCOL_FLAGS])
def test_exactly_one_protocol_is_valid(dest):
    assert _validate_protocol_flags(_ns(**{dest: True})) is None


def test_ont_cdna_plus_dt_primed_is_rejected():
    err = _validate_protocol_flags(_ns(ONT_cDNA=True, dT_primed_cDNA=True))
    assert err is not None
    assert '--ONT-cDNA' in err and '--dT-primed-cDNA' in err
    # the message must say which one PCR-cDNA wants, or it isn't actionable
    assert 'PCR-cDNA' in err


def test_drs_plus_ont_cdna_is_rejected():
    assert _validate_protocol_flags(_ns(drs=True, ONT_cDNA=True)) is not None


def test_bare_namespace_is_safe():
    # Callers that build a partial namespace must not crash the validator.
    assert _validate_protocol_flags(argparse.Namespace()) is None
