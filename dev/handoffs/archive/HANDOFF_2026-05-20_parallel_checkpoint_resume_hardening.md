# Parallel Checkpoint/Resume Hardening

Date: 2026-05-20

Status: implemented in Codex working tree on 2026-05-20. Kept as a design note
for review and future agents.

## Problem

`rectify.core.bam.parallel.process_bam_streaming_parallel()` writes each
completed region's rows directly into the shared corrected TSV, then touches a
`region_XXXX.done` sentinel.

Current relevant code:

- `rectify/core/bam/parallel.py:_write_results_chunk()`
- `rectify/core/bam/parallel.py:process_bam_streaming_parallel()`
- `rectify/core/bam/parallel.py` around the loop that writes `region_results`
  and then touches `.done`

Failure modes:

- Crash after rows are appended but before `.done` is touched: resume reruns the
  region and duplicates rows.
- Crash after `.done` is touched but before output data is durable: resume skips
  the region and can drop rows.
- Shared append output makes correctness depend on process and filesystem
  ordering rather than explicit per-region artifacts.

## Desired Behavior

Checkpoint state should be derived from atomic, durable, per-region outputs.
Resume should be idempotent:

- A region is either absent and rerunnable, or complete and reusable.
- Re-running after a crash cannot duplicate or drop rows.
- The final corrected TSV and position index are rebuilt deterministically from
  completed region files.

## Proposed Design

Use a checkpoint directory with one region file per region:

- `region_XXXX.tsv.tmp`
- `region_XXXX.tsv`
- `region_XXXX.done`

Worker or parent flow:

1. Serialize region rows to `region_XXXX.tsv.tmp` using the canonical
   `CORRECTION_TSV_HEADER` / `correction_result_to_tsv_row()` schema.
2. Flush and `os.fsync()` the temp file.
3. `os.replace(tmp, final_region_tsv)`.
4. `fsync()` the checkpoint directory if practical on the platform.
5. Touch/write `region_XXXX.done.tmp`, fsync, then `os.replace()` to
   `region_XXXX.done`.

Resume flow:

- Treat a region as complete only if both `region_XXXX.done` and
  `region_XXXX.tsv` exist.
- If `.tmp` files exist on startup, ignore or remove them; never treat them as
  complete.
- Rebuild the final output TSV by concatenating completed region TSVs in
  original region order, writing one canonical header at the top.
- Rebuild `ProcessingStats` and `_pos_counts` from the concatenated rows, or
  store per-region stats as separate atomic JSON/TSV sidecars and sum them.

This can be implemented without changing `correct_read_3prime()`.

## Files To Touch

- `rectify/core/bam/parallel.py`
- Possibly `rectify/core/bam/output.py` if a helper like
  `write_results_rows(fh, results)` is useful.
- Tests under `tests/test_resume_correctness.py` or a new
  `tests/test_parallel_checkpoint_resume.py`.

## Test Plan

Add deterministic crash-injection tests without using real multiprocessing:

1. **Crash before done**
   - Simulate a completed `region_0000.tsv` without `.done`.
   - Resume should rerun/rewrite the region once, not append duplicate rows.

2. **Crash after done missing data**
   - Simulate `.done` present but `region_0000.tsv` absent.
   - Resume should treat the region as incomplete and rerun it.

3. **Tmp file ignored**
   - Simulate `region_0000.tsv.tmp` and no final file.
   - Resume should ignore tmp and rerun.

4. **Deterministic final concat**
   - Create region files out of completion order.
   - Rebuild final TSV and assert rows are emitted in region index order with a
     single header.

5. **Position index parity**
   - Rebuilt final TSV and index should match a clean non-resumed run for the
     same synthetic regions.

## Validation Commands

Targeted:

```bash
pytest tests/test_resume_correctness.py tests/test_parallel_processing.py -q
```

Then:

```bash
pytest tests/test_correct_command_parallel_default.py tests/test_validation_reads.py -q
```

## Notes

The 2026-05-20 audit-fix tranche already unified the TSV schema used by
streaming and parallel output. Any checkpoint-region writer should reuse that
schema; do not introduce a second header list.
