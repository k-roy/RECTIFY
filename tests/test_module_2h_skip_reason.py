"""Module 2H must say when it is not running (ISSUE-003).

`rectify correct` gates junction refinement on `--aligner-bams` (or
`--junction-pool-cache`) plus `--annotation`, and used to skip it in total
silence: the corrected BAM looked like it had been refined. On the tester's
62k-read human slice that is 6 reads changed instead of 30,703.

`module_2h_skip_reason` is the predicate behind the log line; the run() call
site must agree with `_has_junction_context`, which is asserted below.

Author: Kevin R. Roy (agent S1)
"""

import sys
from pathlib import Path

import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.commands.correct_command import module_2h_skip_reason  # noqa: E402

RUNS = {'aligner_bams': ['a.bam'], 'annotation_path': 'genes.gtf'}


def test_the_configuration_that_actually_runs_2h_has_no_reason():
    assert module_2h_skip_reason(RUNS) is None


def test_pool_cache_alone_is_enough():
    assert module_2h_skip_reason(
        {'junction_pool_cache': 'pool.pkl', 'annotation_path': 'genes.gtf'}
    ) is None


def test_a_bare_correct_run_names_both_missing_inputs():
    reason = module_2h_skip_reason({})
    assert reason is not None
    assert '--aligner-bams' in reason
    assert '--annotation' in reason


def test_missing_aligner_bams_is_named():
    reason = module_2h_skip_reason({'annotation_path': 'genes.gtf'})
    assert reason is not None and '--aligner-bams' in reason
    assert '--annotation' not in reason


def test_missing_annotation_is_named():
    reason = module_2h_skip_reason({'aligner_bams': ['a.bam']})
    assert reason is not None and '--annotation' in reason


def test_both_skip_flags_together_disable_the_whole_block():
    reason = module_2h_skip_reason(
        {**RUNS, 'skip_junction_refinement': True, 'apply_3ss_rescue': False})
    assert reason is not None and '--skip-junction-refinement' in reason


def test_skip_refinement_alone_still_builds_the_pool_for_2f():
    """--skip-junction-refinement keeps the pool index that Module 2F needs."""
    assert module_2h_skip_reason({**RUNS, 'skip_junction_refinement': True}) is None


def test_predicate_matches_the_gate_in_run():
    """The log line must not be able to disagree with the branch it explains."""
    src = (RECTIFY_ROOT / 'rectify/core/commands/correct_command.py').read_text()
    gate = src.index('_has_junction_context = (')
    call = src.rindex('module_2h_skip_reason(config)')   # the call site, not the def
    assert gate < call
    assert 'if not _has_junction_context:' in src
