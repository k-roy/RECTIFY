"""The generated long-read chunk scripts, executed — prescan -> per-arm correct -> chunk merge.

Three defects hid here because nothing ran the scripts (their text was tested, never their effect):

* ISSUE-079 — the chunk merge replayed each corrected TSV onto the RAW chunk BAM, undoing every
  Module 2H placement in ``corrected_consensus.bam``.
* ISSUE-081 — ``correct`` stopped emitting the flat ``corrected_reads.tsv`` (manifest-only default),
  which is the only name the scripts look for: the skip-check never fired and the merge found no
  TSVs.
* ISSUE-082 — a ``$L_SCRATCH`` in a comment inside an UNQUOTED heredoc; under ``set -u`` the merge
  died on any host that does not define it.
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from rectify.core.commands import split_command as sc
from tests.test_issue079_writer_input_handoff import ARMS, build_two_arm_library

pytestmark = pytest.mark.skipif(shutil.which('bash') is None, reason='needs bash')
REPO = Path(__file__).resolve().parents[1]


def _run(script, out, scratch):
    env = {k: v for k, v in os.environ.items() if k not in ('L_SCRATCH', 'SGE_TASK_ID', 'PBS_ARRAY_INDEX')}
    env.update(SCRATCH=str(scratch), SLURM_ARRAY_TASK_ID='0', SLURM_CPUS_PER_TASK='1',
               PYTHONPATH=str(REPO))
    return subprocess.run(['bash', str(out / script)], cwd=str(out), env=env,
                          capture_output=True, text=True)


def test_generated_chunk_route_keeps_2h_placements_and_resumes(tmp_path):
    (tmp_path / 'lib').mkdir()
    bams, fasta, annotation, _, truth = build_two_arm_library(tmp_path / 'lib')
    out, scratch = tmp_path / 'out', tmp_path / 'scratch'
    scratch.mkdir()
    (out / 'merged_bams').mkdir(parents=True)
    for arm, bam in bams.items():
        chunk = out / 'aligner_chunks' / arm / 'chunk_000'
        chunk.mkdir(parents=True)
        for dest in (out / 'merged_bams' / f'S.{arm}.bam', chunk / f'S_chunk_000.{arm}.bam'):
            shutil.copy(bam, dest)
            shutil.copy(str(bam) + '.bai', str(dest) + '.bai')
    sc._generate_scripts(argparse.Namespace(
        output_dir=out, genome=fasta, annotation=annotation, scheduler='slurm',
        python_path=sys.executable, rectify_src=str(REPO), other_aligners=list(ARMS),
        skip_map_pacbio=True, slurm_partition=None, slurm_account=None,
        uge_queue=None, uge_pe='smp', pbs_queue='workq'), n_chunks=1, sample_prefix='S')

    done = _run('run_prescan.sh', out, scratch)
    assert done.returncode == 0, done.stderr[-2000:]
    for arm in ARMS:
        done = _run(f'run_array_correct_{arm}.sh', out, scratch)
        assert done.returncode == 0, done.stderr[-2000:]
        chunk = out / 'aligner_chunks' / arm / 'chunk_000'
        assert (chunk / 'corrected_reads.tsv').exists()              # ISSUE-081
        assert (chunk / 'corrected_reads.writer_input.bam').exists()  # ISSUE-079
        assert 'already corrected' in _run(f'run_array_correct_{arm}.sh', out, scratch).stdout

    done = _run('run_array_chunk_merge.sh', out, scratch)             # ISSUE-082: L_SCRATCH unset
    assert done.returncode == 0, done.stderr[-2000:]
    with pysam.AlignmentFile(str(out / 'chunks' / 'chunk_000' / 'corrected_consensus.bam')) as bam:
        got = {r.query_name: (r.reference_start, r.cigarstring) for r in bam}
    assert got == truth

    # The retained inputs are released once the chunk's products stand; both stages then skip.
    assert not list(out.rglob('*.writer_input.bam'))
    assert 'already merged' in _run('run_array_chunk_merge.sh', out, scratch).stdout
    assert 'already corrected' in _run(f'run_array_correct_{ARMS[0]}.sh', out, scratch).stdout

    # A lost merge product cannot be rebuilt from an evicted pair — and must not pretend to be.
    (out / 'chunks' / 'chunk_000' / 'corrected_consensus.bam.bai').unlink()
    redo = _run('run_array_chunk_merge.sh', out, scratch)
    assert redo.returncode != 0 and 'evicted' in redo.stderr
    # ... the correct task sees the same thing and rebuilds its pair instead of skipping.
    assert 'already corrected' not in _run(f'run_array_correct_{ARMS[0]}.sh', out, scratch).stdout
    assert (out / 'aligner_chunks' / ARMS[0] / 'chunk_000' / 'corrected_reads.writer_input.bam').exists()
