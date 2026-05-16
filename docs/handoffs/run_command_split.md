# Split `rectify/core/commands/run_command.py`

**Current:** 2,949 lines (largest single file in the repo). The
`rectify run-all` orchestrator — drives the whole pipeline,
single-sample or via manifest, with optional batch (SLURM) dispatch
and resume-on-restart.

## Goal

Carve into per-stage runners. The current monolith intermixes
pipeline-stage logic (align → correct → analyze → junction-aggregate)
with sample dispatch (single vs manifest vs chunked-batch) with
helper functions (path resolution, BAM integrity checks). Splitting
along stage boundaries makes each piece independently inspectable.

## DO THIS LAST

Of the seven handoff briefs, run this one **after** the others have
settled. Almost every other split touches a function that
`run_command.py` calls. Doing this first means rewriting your work
mid-refactor.

## Proposed split

| New file | Contents | Approx lines |
|---|---|---|
| `commands/run_command.py` (kept) | `run()`, `create_run_parser()`, top-level dispatch (single vs manifest vs chunked-batch) | ~500 |
| `commands/run/single_sample.py` (new) | `_run_single_sample`, `_process_one_sample`, `_rectified_bam_path` | ~1,000 |
| `commands/run/multi_sample.py` (new) | `_run_multi_sample`, `_run_analysis_manifest` | ~250 |
| `commands/run/chunked_batch.py` (new) | `_generate_chunked_pipeline` (currently 900 lines — by far the biggest function in this file) | ~900 |
| `commands/run/stages.py` (new) | `_run_alignment`, `_run_correction`, `_run_correction_per_aligner`, `_combine_corrected_tsvs`, `_run_analysis`, `_run_junction_aggregation` | ~580 |
| `commands/run/helpers.py` (new) | `_resolve_reference_paths`, `_collect_per_aligner_bams`, `_bam_has_md_tags`, `_validate_bam_integrity` | ~130 |

Note: this creates a `commands/run/` sub-subpackage. Make sure
`commands/run/__init__.py` re-exports the symbols `run_command.py`
needs.

## Migration plan

1. Create `commands/run/__init__.py`.
2. Move helpers first (`commands/run/helpers.py`) — these are
   used by everything else and have no further dependencies.
3. Move `commands/run/stages.py` — per-pipeline-stage runners that
   the sample-level functions call.
4. Move `commands/run/single_sample.py`,
   `commands/run/multi_sample.py`,
   `commands/run/chunked_batch.py`.
5. Update `run_command.py` to import from the new locations.

## Test gate

- `tests/test_run_command_wiring.py` — direct tests of the wiring
- `tests/test_correct_command_drs.py` — end-to-end via subprocess
- `tests/test_validation_reads.py` — full pipeline
- `tests/test_cdna_correct.py` — multi-stage cDNA pipeline (slow)
- `tests/test_cdna_chain_canary.py` — full chain canary (slow)

## Critical pitfalls

- The 900-line `_generate_chunked_pipeline` is a single function that
  emits shell scripts. Don't try to break it up; just move it.
- `_run_single_sample` shares ~30 local variables with
  `_process_one_sample`. Audit which need to be parameters versus
  re-derived inside the moved function.
- `_run_correction_per_aligner` runs the correct step once per
  aligner. The aligner-name list is built in `_run_alignment` and
  passed through; preserve the contract.
- The resume-on-restart logic (skip-if-{output}-exists) is sprinkled
  throughout `_process_one_sample`. After the split, that logic
  should remain in `single_sample.py` — don't push it down into
  `stages.py`, because the resume decision is sample-level.
- `commands/batch_command.py` imports from `run_command` for
  manifest parsing. After the split, update its imports.
