# Parallel BAM Writer Abort Investigation

Date: 2026-05-20

Status: mitigated in Codex working tree on 2026-05-20. The low-level BAM writer
tests pass, and worker-process execution is now opt-in because this runtime can
abort during process launch.

## Problem

Audit observation:

```text
pytest tests/test_bam_writer.py tests/test_bam_writer_parallel_smoke.py -q
```

passed low-level writer tests and the `n_threads=1` smoke case, then repeatedly
hit fatal Python abort traces in `multiprocessing.Pool` while exercising
`write_corrected_bam_parallel()`.

Primary file:

- `rectify/core/bam/bam_writer_parallel.py`

Related files:

- `rectify/core/bam/bam_writer.py`
- `rectify/core/bam/read_edits.py`
- `tests/test_bam_writer_parallel_smoke.py`
- `tests/test_bam_writer.py`

## Why This Matters

The serial BAM writers are the reference behavior. The parallel writer should be
a throughput optimization only. It must not introduce:

- process-level aborts,
- nondeterministic output ordering,
- missing/duplicated reads,
- CIGAR differences from serial output,
- shared `pysam.AlignmentFile` or htslib state across processes.

## Initial Suspicions

Areas to inspect first:

- Whether any `pysam.AlignmentFile`, `AlignedSegment`, header object, genome
  dict, or correction dict is being shared across forked workers in a way htslib
  dislikes.
- Whether workers write BAM records to the same path concurrently.
- Whether temp BAMs are closed/indexed/merged while worker processes still hold
  file handles.
- Whether multiprocessing start method differs across platforms. macOS + htslib
  + forked interpreter state can be fragile.
- Whether Matplotlib/import side effects are involved in test startup. The
  broader test run can emit abort traces during imports, but the audit's
  parallel BAM issue specifically implicated multiprocessing pool execution.

## Desired Behavior

Parallel BAM writing should:

- Open input BAM independently inside each worker.
- Open output temp BAM independently inside each worker.
- Never pass live `pysam` objects between processes.
- Emit one complete temp BAM per region/chunk.
- Merge temp BAMs after all workers have exited successfully.
- Preserve serial writer behavior for CIGAR edits.

## Suggested Implementation Direction

Keep serial writer logic as the single source of truth:

1. Extract a per-read edit function from `bam_writer.py` if needed, but avoid
   duplicating CIGAR surgery in parallel-specific code.
2. In each worker:
   - open input BAM fresh,
   - iterate assigned region/read IDs,
   - load only needed corrections,
   - write to `worker_N.region_M.bam.tmp`,
   - close BAM,
   - `os.replace()` temp to final temp BAM path.
3. Parent process:
   - waits for all workers,
   - fails fast if any worker fails,
   - merges temp BAMs in deterministic region order,
   - indexes the final BAM.
4. Prefer `spawn` or `forkserver` on macOS if fork-related htslib instability is
   confirmed. Measure overhead before making it default.

## Test Plan

Add a test that compares serial and parallel output on a small synthetic BAM:

- `n_threads=1` should match serial.
- `n_threads=2` should match serial.
- Compare read names, reference starts, CIGAR strings, sequences, qualities, and
  key tags (`cp`, `Xj`, `Xv` if present).
- Include at least:
  - simple 3' clip,
  - 5' rescue,
  - softclip rescue,
  - overcall rescue,
  - a read with no correction.

Add a stress smoke test:

- 20 to 100 synthetic reads split across at least two references or regions.
- Run parallel writer multiple times in the same pytest process to catch handle
  leaks and pool cleanup failures.

## Validation Commands

Targeted:

```bash
pytest tests/test_bam_writer.py tests/test_bam_writer_parallel_smoke.py -q
```

Broader BAM-related pass:

```bash
pytest tests/test_validation_reads.py tests/test_bam_provenance.py tests/test_restore_polya_manifest.py -q
```

## Notes From 2026-05-20 Audit-Fix Tranche

The dual serial writer was fixed to match standalone hardclip/softclip behavior
for:

- `five_prime_upstream_trim`
- `oc_homopolymer_extension`
- `oc_overcall_count`
- `oc_terminal_base`

Parallel writer comparisons should include these fields because they previously
caused mode-specific BAM CIGAR divergence.
