"""planning 873 / 861 G-5, G-10, G-21 + the scratch free-space guard.

Four independent regressions, each of which shipped as a silent wrong answer rather than a crash:

1. ``run-all --short-read`` had NO adapter-trim stage, and the adapter was invisible in the BAM
   because HISAT2 runs ``--no-softclip`` and STAR ``EndToEnd`` on that path (planning 861 §S0c:
   13.0 % of R1 carried the TruSeq adapter, 11.5 % were full-length dimers, 0.01 % of records
   showed any terminal clip).
2. ``--bwa-path`` / ``--bbmap-path`` did not exist on ``run-all``, so a cluster env with those
   executables outside ``PATH`` silently lost the two short-read arms (planning 862 G-10).
3. ``correct`` printed and wrote ``Protocol: DRS (direct RNA-seq)`` for every short-read and
   NET-seq run (planning 861 G-21) — a log reader auditing protocol mode was told the opposite
   of the truth.
4. ``make_job_scratch_dir`` staged onto ``$SCRATCH`` however full it was; on Hoffman2 a near-full
   ``/u/scratch`` kills every array task inside a minute with an EMPTY stderr.

Everything here is a pure-function test: no cutadapt, no aligner, no BAM.
"""
import argparse
import os
from pathlib import Path

import pytest

from rectify.core.commands.run.short_read_trim import (
    TRUSEQ_R1_ADAPTER,
    TRUSEQ_R2_ADAPTER,
    build_short_read_cutadapt_cmd,
    parse_cutadapt_report,
)
from rectify.core.bam.processing_stats import ProcessingStats, generate_stats_report, write_stats_tsv
from rectify import slurm


# ──────────────────────────────────────────────────────────────────────────────
# 1. the trim command builder (G-5)
# ──────────────────────────────────────────────────────────────────────────────
def test_single_end_cmd_has_adapter_and_min_length():
    cmd = build_short_read_cutadapt_cmd(Path('r1.fastq.gz'), Path('out.fastq.gz'))
    assert cmd[0] == 'cutadapt'
    assert '-a' in cmd and cmd[cmd.index('-a') + 1] == TRUSEQ_R1_ADAPTER
    # -m 20, not 0: a 0-length record is rejected by bbmap/magicblast and silently
    # dropped by STAR, and 11.5 % of the 861 library were full-length dimers.
    assert cmd[cmd.index('-m') + 1] == '20'
    assert cmd[-1] == 'r1.fastq.gz'
    assert '-A' not in cmd and '-p' not in cmd


def test_paired_end_cmd_trims_both_mates_in_one_pass():
    cmd = build_short_read_cutadapt_cmd(
        Path('r1.fq.gz'), Path('o1.fq.gz'), read2=Path('r2.fq.gz'), output2=Path('o2.fq.gz'))
    assert cmd[cmd.index('-A') + 1] == TRUSEQ_R2_ADAPTER
    assert cmd[cmd.index('-p') + 1] == 'o2.fq.gz'
    # both inputs are on ONE command line; trimming the mates separately desynchronises
    # the FASTQs as soon as one mate falls below -m.
    assert cmd[-2:] == ['r1.fq.gz', 'r2.fq.gz']


def test_paired_end_needs_both_output_and_read2():
    with pytest.raises(ValueError):
        build_short_read_cutadapt_cmd(Path('r1'), Path('o1'), read2=Path('r2'))
    with pytest.raises(ValueError):
        build_short_read_cutadapt_cmd(Path('r1'), Path('o1'), output2=Path('o2'))


def test_quality_and_nextseq_are_mutually_exclusive():
    with pytest.raises(ValueError):
        build_short_read_cutadapt_cmd(Path('r1'), Path('o'), quality_cutoff=20, nextseq_trim=20)
    assert '-q' in build_short_read_cutadapt_cmd(Path('r1'), Path('o'), quality_cutoff=20)
    assert '--nextseq-trim=20' in build_short_read_cutadapt_cmd(Path('r1'), Path('o'), nextseq_trim=20)


def test_custom_adapter_and_executable_are_honoured():
    cmd = build_short_read_cutadapt_cmd(
        Path('r1'), Path('o'), adapter='ACGTACGT', cutadapt_path='/opt/env/bin/cutadapt', threads=8)
    assert cmd[0] == '/opt/env/bin/cutadapt'
    assert cmd[cmd.index('-a') + 1] == 'ACGTACGT'
    assert cmd[cmd.index('-j') + 1] == '8'


def test_netseq_linker_is_not_the_short_read_default():
    """The NET-seq stage has its own adapter and its own no-``-u`` rule; the two must not merge."""
    from rectify.core.commands.run.netseq_pipeline import CHURCHMAN_3P_LINKER
    assert TRUSEQ_R1_ADAPTER != CHURCHMAN_3P_LINKER
    assert '-u' not in build_short_read_cutadapt_cmd(Path('r1'), Path('o'))


def test_report_parser_reads_the_minimal_report_and_tolerates_junk():
    txt = 'status\tin_reads\tin_bp\tout_reads\tw/adapters\tout_bp\nOK\t200000\t15000000\t176984\t25911\t12000000\n'
    rep = parse_cutadapt_report(txt)
    assert rep['in_reads'] == 200000 and rep['out_reads'] == 176984 and rep['w/adapters'] == 25911
    # a report that could not be read must not raise: the FASTQ is still valid.
    assert parse_cutadapt_report('') == {}
    assert parse_cutadapt_report('only one line') == {}
    assert parse_cutadapt_report('a\tb\n1\n') == {}


# ──────────────────────────────────────────────────────────────────────────────
# 2. run-all exposes the flags and forwards the aligner paths (G-10)
# ──────────────────────────────────────────────────────────────────────────────
def _run_all_parser():
    from rectify.core.commands.run_command import create_run_parser
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='command')
    create_run_parser(sub)
    return p


@pytest.mark.parametrize('flag,dest,value,expected', [
    ('--short-read-adapter', 'short_read_adapter', 'ACGT', 'ACGT'),
    ('--short-read-adapter-r2', 'short_read_adapter_r2', 'TTTT', 'TTTT'),
    ('--short-read-min-length', 'short_read_min_length', '31', 31),
    ('--short-read-nextseq-trim', 'short_read_nextseq_trim', '20', 20),
    ('--bwa-path', 'bwa_path', '/opt/bwa', '/opt/bwa'),
    ('--bbmap-path', 'bbmap_path', '/opt/bbmap.sh', '/opt/bbmap.sh'),
])
def test_run_all_accepts_the_new_flags(flag, dest, value, expected):
    args = _run_all_parser().parse_args(['run-all', 'r.fastq.gz', '--short-read', flag, value])
    assert getattr(args, dest) == expected


def test_trim_defaults_are_on_with_a_skip_switch():
    args = _run_all_parser().parse_args(['run-all', 'r.fastq.gz', '--short-read'])
    assert args.skip_short_read_trim is False           # the stage runs by default
    assert args.short_read_min_length == 20
    assert args.short_read_adapter is None              # None -> the TruSeq stem
    args2 = _run_all_parser().parse_args(
        ['run-all', 'r.fastq.gz', '--short-read', '--skip-short-read-trim'])
    assert args2.skip_short_read_trim is True


def test_run_alignment_forwards_bbmap_and_bwa_paths_into_align_args():
    """The Namespace handed to run_align is hand-built, so a missing key is silently dropped."""
    import inspect
    from rectify.core.commands.run import stages
    sig = inspect.signature(stages._run_alignment)
    assert 'bbmap_path' in sig.parameters and 'bwa_path' in sig.parameters
    src = inspect.getsource(stages._run_alignment)
    assert 'bbmap_path=bbmap_path or' in src and 'bwa_path=bwa_path or' in src
    # and the run-all entry point must actually pass them through
    from rectify.core.commands.run import single_sample
    ss = inspect.getsource(single_sample)
    assert "bbmap_path=getattr(args, 'bbmap_path', None)" in ss
    assert "bwa_path=getattr(args, 'bwa_path', None)" in ss


def test_single_sample_runs_the_trim_before_alignment():
    """The stage must sit ahead of _run_alignment, or the aligner sees untrimmed reads."""
    import inspect
    from rectify.core.commands.run import single_sample
    src = inspect.getsource(single_sample)
    i_trim = src.index('run_short_read_trim(')
    i_align = src.index('_run_alignment(')
    assert i_trim < i_align
    assert "skip_short_read_trim" in src


# ──────────────────────────────────────────────────────────────────────────────
# 3. the protocol label (G-21)
# ──────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize('protocol,needle', [
    ('drs', 'DRS (direct RNA-seq)'),
    ('dt_cdna', 'oligo-dT-primed cDNA'),
    ('ont_cdna', 'ONT PCR-cDNA'),
    ('short_read', 'short-read Illumina (--short-read)'),
    ('netseq', 'NET-seq (--netseq)'),
])
def test_stats_report_names_the_real_protocol(protocol, needle):
    report = generate_stats_report(ProcessingStats(), protocol=protocol)
    assert needle in report


def test_short_read_and_netseq_are_no_longer_labelled_drs(tmp_path):
    for protocol in ('short_read', 'netseq'):
        assert 'DRS (direct RNA-seq)' not in generate_stats_report(ProcessingStats(), protocol=protocol)
        out = tmp_path / f'{protocol}.tsv'
        write_stats_tsv(ProcessingStats(), out, protocol=protocol)
        head = out.read_text()
        assert '# Protocol: DRS' not in head
    assert '# Protocol: DRS' in (lambda p: (write_stats_tsv(ProcessingStats(), p, protocol='drs'),
                                            p.read_text())[1])(tmp_path / 'drs.tsv')


def test_correct_command_resolves_short_read_protocol():
    """--short-read --dT-primed-cDNA is a QuantSeq library and KEEPS the oligo-dT label."""
    import inspect
    from rectify.core.commands import correct_command
    src = inspect.getsource(correct_command)
    assert "else 'short_read' if _is_short_read" in src
    i_dt = src.index("else 'dt_cdna' if is_dt_primed")
    i_sr = src.index("else 'short_read' if _is_short_read")
    assert i_dt < i_sr


# ──────────────────────────────────────────────────────────────────────────────
# 4. the scratch free-space guard
# ──────────────────────────────────────────────────────────────────────────────
def test_scratch_free_gib_reports_a_number_for_a_real_path(tmp_path):
    free = slurm.scratch_free_gib(tmp_path)
    assert free is None or free >= 0.0


def test_scratch_free_gib_returns_none_for_a_bad_path():
    assert slurm.scratch_free_gib(Path('/definitely/not/a/mount/point/xyzzy')) is None


def test_full_scratch_is_refused_and_unknown_space_is_allowed(monkeypatch, tmp_path):
    monkeypatch.setattr(slurm, 'scratch_free_gib', lambda p: 1.0)
    assert slurm._scratch_has_room(tmp_path) is False
    monkeypatch.setattr(slurm, 'scratch_free_gib', lambda p: 500.0)
    assert slurm._scratch_has_room(tmp_path) is True
    # unknown free space must NOT block the run — the check is a guard, not a gate
    monkeypatch.setattr(slurm, 'scratch_free_gib', lambda p: None)
    assert slurm._scratch_has_room(tmp_path) is True


def test_floor_is_env_overridable(monkeypatch, tmp_path):
    monkeypatch.setattr(slurm, 'scratch_free_gib', lambda p: 50.0)
    monkeypatch.setenv('RECTIFY_MIN_SCRATCH_GIB', '100')
    assert slurm._scratch_has_room(tmp_path) is False
    monkeypatch.setenv('RECTIFY_MIN_SCRATCH_GIB', 'not-a-number')   # must not raise
    assert slurm._scratch_has_room(tmp_path) is True


def test_make_job_scratch_dir_falls_back_when_scratch_is_full(monkeypatch, tmp_path):
    monkeypatch.setattr(slurm, 'get_scratch_dir', lambda: tmp_path)
    monkeypatch.setattr(slurm, 'scratch_free_gib', lambda p: 0.5)
    assert slurm.make_job_scratch_dir('rectify_test') is None      # caller writes to output dir
    monkeypatch.setattr(slurm, 'scratch_free_gib', lambda p: 999.0)
    made = slurm.make_job_scratch_dir('rectify_test')
    assert made is not None and made.exists() and made.name.startswith('rectify_test_')
