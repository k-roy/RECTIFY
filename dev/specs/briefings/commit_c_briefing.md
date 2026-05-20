# Opus Briefing — Commit C (Dynamic aligner-queue concurrency)

**Goal:** add Axis B concurrency for multi-aligner correction without breaking
the correct-first pipeline order or the Commit B manifest/sidecar contract.

This is an Opus-owned commit. Do not hand the whole implementation to a
mechanical subagent: the hard part is deciding the execution boundary, not
writing a pool loop.

---

## 1. Current State To Verify First

Run these before editing:

```bash
cd /Users/kevinroy/work/rectify
git log --oneline -5
git status --short rectify/core/commands/run/stages.py rectify/core/commands/correct_command.py rectify/core/bam/parallel.py rectify/core/bam/bam_writer_parallel.py
ssh sherlock 'echo ok'
ssh hoffman2 'bash -lc "cd /u/home/k/kevinroy/software/rectify && git log --oneline -1"'
ssh sherlock 'bash --norc --noprofile -c "cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify && git log --oneline -1"'
```

Expected baseline at briefing time:

- M1/GitHub/H2/Sherlock are synced at `6c7f846` or newer.
- `rectify/core/commands/run/stages.py` is clean.
- `rectify/core/commands/correct_command.py` is clean.
- The repo may have Kevin WIP elsewhere. Stage only explicit files.

Also read:

- `dev/specs/parallel_refactor_plan.md` §0 and Commit C.
- `dev/specs/profile_results.md` §7 for the protocol-conditional payoff.
- `HANDOFF.md` and `AGENT_FIXES.md` before any cluster run.
- Drive `rectify/CLAUDE.md` for HPC chunking/scratch rules.

---

## 2. Important Boundary Decision

Commit C is **not** just "run five existing `_run_correction()` calls in a
ProcessPoolExecutor."

That simpler implementation would multiply the current per-aligner region
workers: five aligners × `--threads` workers. On a 16-core node with five
aligners, that can accidentally spawn ~80 correction workers plus writer
workers, violating the resource discipline this refactor is supposed to add.

The planned design is a single shared worker budget:

```text
total active correction workers <= aligner_queue_workers
```

Workers pull `(aligner_name, region)` tasks from one queue, so fast aligners do
not leave static worker slices idle and slow aligners do not over-subscribe the
node.

If the exact shared-region queue is too invasive for one commit, land a
smaller preparatory commit instead:

1. Add resource detection and `--aligner-concurrency`.
2. Keep default `auto` resolving to `1` until the region queue exists.
3. Add tests proving the parser/resource policy.
4. Do **not** claim Axis B performance.

---

## 3. Files And Likely Shape

### 3.1 Add `rectify/core/utils/resources.py`

Add a small, dependency-light helper:

```python
def detect_machine_class(total_mem_gb: float | None = None, cpu_count: int | None = None) -> str:
    """Return 'm1_laptop' or 'cluster_node'."""
```

Rules from the spec:

- `m1_laptop` when `total_mem_gb < 16 and cpu_count < 8`.
- Otherwise `cluster_node`.
- Make it injectable for tests so CI does not depend on the runner.

Add:

```python
def resolve_aligner_concurrency(value: str, total_threads: int, machine_class: str) -> int:
```

Suggested policy:

- `"1"` -> `1`.
- integer string `N` -> clamp to `[1, min(N, total_threads)]`.
- `"auto"` on `m1_laptop` -> `1`.
- `"auto"` on `cluster_node` -> `max(1, total_threads - 2)`.

This helper should not import numpy, psutil, pandas, or pysam.

### 3.2 Modify `run-all` parser

Add `--aligner-concurrency {auto,1,N}` to `rectify/core/commands/run_command.py`.

Default: `auto`.

Help text should be explicit:

- `auto` disables Axis B on small laptops.
- On cluster nodes, `auto` reserves two CPUs for merge/main-process overhead.
- `1` preserves sequential-aligner behavior.

### 3.3 Modify `_run_correction_per_aligner`

Current verified boundary:

- `rectify/core/commands/run/stages.py:_run_correction_per_aligner` loops over
  `per_aligner_bams`.
- It calls `_run_correction(...)`, which builds an argparse namespace and calls
  `correct_command.run(correct_args)`.
- Commit B's manifest-only default means callers must locate
  `corrected_reads.manifest.tsv` when `corrected_reads.tsv` is absent.

Minimum safe behavior:

- Always honor `--aligner-concurrency 1` by preserving today's serial loop.
- For any concurrent implementation, pass a reduced per-aligner thread count so
  aggregate workers do not exceed the resolved budget.
- Keep `aligner_bams=list(per_aligner_bams.values())` for every aligner so
  Module 2H still has the cross-aligner candidate junction pool.
- Treat both `corrected_reads.tsv` and `corrected_reads.manifest.tsv` as valid
  per-aligner outputs.

Preferred full behavior:

- Introduce an internal region-task API below `correct_command.run`, so region
  tasks can be queued across aligners while sharing one worker pool.
- Final per-aligner merge/manifest/sidecar steps happen after all region tasks
  for that aligner complete.
- Do not write a stage sidecar before final BAM + manifest durability.

---

## 4. Tests

Add focused unit tests before cluster validation:

- `tests/test_resources.py`
  - machine-class detection with injected CPU/memory values.
  - `resolve_aligner_concurrency("auto", ...)` on laptop vs cluster.
  - numeric clamp and invalid-value errors.
- Parser test for `rectify run-all --aligner-concurrency`.
- `_run_correction_per_aligner` test with monkeypatched `_run_correction`:
  - `--aligner-concurrency 1` preserves serial behavior.
  - manifest-only output is accepted as success.
  - aggregate thread budget is not multiplied by number of aligners.

Do not require real aligners in unit tests.

---

## 5. Cluster Acceptance

Only run when queues are sane. Check `AGENT_FIXES.md` first.

H2 acceptance target from the spec:

- Han wt_R1, five aligners, 16-core node.
- Correct phase wall time <= slowest aligner alone × 1.3.
- No OOM on 64 GB H2 node.
- `--aligner-concurrency 1` output equivalent to the prior sequential-aligner
  path.

M1 acceptance:

- `auto` falls back to sequential aligners.
- Active process count should peak near `N+1`, not `5N+1`.

If full Han wt_R1 remains blocked, run only unit tests and commit as a
preparatory Commit C. Do not mark the performance acceptance box complete.

---

## 6. Non-Goals

- Do not change the Commit B default manifest policy.
- Do not touch `rectify/utils/provenance.py`.
- Do not implement analyze streaming; that is Commit D.
- Do not submit the 6.7M-read heap-corruption smoke from a saturated queue.
- Do not commit Kevin WIP files shown by `git status`.

