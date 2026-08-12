"""The dropped-aligner gate: a partial panel must not pass as exit 0.

Every "aligner unavailable" path in `align_command._run_one_aligner` logs a
warning and returns ``(aligner, None)``, letting the run continue and exit 0.
That is a false green: a requested 3-aligner panel silently becomes a
2-aligner one, and any acceptance gate that checks only the exit code records
a pass. Observed live on SCG — deSALT absent from the env produced a green
"3-aligner" run with 2 arms.

`--require-aligners` makes that fatal. The DROPPED-ALIGNER summary line is
emitted either way, so the condition is greppable even in legacy runs.
"""

import argparse
import logging

from rectify.core.commands.align_command import create_align_parser


def _parse(argv):
    root = argparse.ArgumentParser()
    sub = root.add_subparsers(dest='command')
    create_align_parser(sub)
    return root.parse_args(argv)


def test_flag_parses_and_defaults_off(tmp_path):
    """Default stays permissive — existing runbooks must not change behaviour."""
    base = ['align', 'reads.fastq', '--genome', 'g.fa', '-o', str(tmp_path)]
    assert _parse(base).require_aligners is False
    assert _parse(base + ['--require-aligners']).require_aligners is True


def _drop_gate(aligners, results, require, caplog):
    """Mirror of the gate in align_command (see 'Dropped-aligner gate').

    Kept as a local reimplementation because the real gate sits mid-way through
    a ~500-line function that needs a genome, a FASTQ and real aligner binaries
    to reach. The assertions below pin the CONTRACT (which aligners count as
    dropped, and when it is fatal); `test_gate_source_is_wired` pins that the
    real code path still carries it.
    """
    dropped = [a for a in aligners if not results.get(a)]
    if dropped:
        logging.getLogger('rectify').error(
            "DROPPED-ALIGNER: %d of %d requested aligners produced no BAM: %s",
            len(dropped), len(aligners), ', '.join(dropped),
        )
        if require:
            return 1
    return 0


def test_missing_arm_is_fatal_under_require(caplog):
    caplog.set_level(logging.ERROR)
    rc = _drop_gate(
        ['minimap2', 'uLTRA', 'deSALT'],
        {'minimap2': 'a.bam', 'uLTRA': 'b.bam', 'deSALT': None},
        require=True,
        caplog=caplog,
    )
    assert rc == 1
    assert 'DROPPED-ALIGNER' in caplog.text
    assert 'deSALT' in caplog.text


def test_missing_arm_warns_but_passes_without_require(caplog):
    """The legacy behaviour, but now greppable rather than silent."""
    caplog.set_level(logging.ERROR)
    rc = _drop_gate(
        ['minimap2', 'deSALT'],
        {'minimap2': 'a.bam', 'deSALT': None},
        require=False,
        caplog=caplog,
    )
    assert rc == 0
    assert 'DROPPED-ALIGNER' in caplog.text


def test_full_panel_is_silent(caplog):
    caplog.set_level(logging.ERROR)
    rc = _drop_gate(
        ['minimap2', 'deSALT'],
        {'minimap2': 'a.bam', 'deSALT': 'b.bam'},
        require=True,
        caplog=caplog,
    )
    assert rc == 0
    assert 'DROPPED-ALIGNER' not in caplog.text


def test_gate_source_is_wired():
    """Guard the real call site — the local mirror above could otherwise drift.

    Asserts the gate exists in `align_command`, keys off `require_aligners`,
    and sits BEFORE the resolver post-pass (an absent deSALT must fail before
    the resolver spends an arm's worth of work).
    """
    import inspect
    from rectify.core.commands import align_command

    src = inspect.getsource(align_command.run_align)
    assert 'DROPPED-ALIGNER' in src
    assert "require_aligners" in src
    gate_at = src.index('DROPPED-ALIGNER')
    resolver_at = src.index('if want_resolver:')
    assert gate_at < resolver_at, "gate must precede the resolver post-pass"
