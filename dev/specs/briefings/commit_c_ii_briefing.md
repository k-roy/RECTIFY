# Briefing — Commit C-II (Shared per-aligner pool with work-stealing)

**Status:** DRAFT v2 — scratch outside the rectify repo. Supersedes the
2026-05-22 v1 (which was based on a wrong premise — see §0). Move into
`rectify/dev/specs/briefings/commit_c_ii_briefing.md` when concurrent
agent activity is clear and the index is safe.

**Predecessor:** Preparatory Commit C landed at `d762b91` on
2026-05-22 — resolver plumbed, manifest-output regression fixed, loop
stays sequential. Commit D (analyze sidecars) landed at `8c9f80c`. The
reanchor perf fix landed at `55089f7`. HEAD as of writing: `55089f7`.

**Owner:** Opus. Advisor-gated. Split into two commits:

- **C-II.0** — Refactor `bam_parallel.py` worker state model. Pure
  refactor, semantic no-op. Pytest-only acceptance.
- **C-II.1** — Add shared cross-aligner pool with work-stealing.
  Cluster-acceptance gated (Han wt_R1 ≤ slowest×1.3).

The split is also race-safer given 5 active claude sessions on M1 —
each commit assembled and pushed atomically before starting the next.

---

## 0. Correction to v1 — what changed

The v1 briefing claimed Commit C-II would "extract per-region
correction work out of `correct_command.run()`." **That premise was
wrong.** `correct_command.run()` already delegates region-parallel
correction to `bam_parallel.process_bam_file_parallel` (line 885 in
`correct_command.py`). Region decomposition is *already done*.

What is NOT done is sharing the worker pool across per-aligner
correction calls. Today each call to `process_bam_file_parallel`
spawns its own pool (line 769 of `bam_parallel.py`:
`with ctx.Pool(n_threads, ...) as pool:`) with workers initialized
via `_init_region_worker_state(genome_path, polya_model_path,
worker_kwargs)`. The worker init bakes in **per-aligner state**
(`bam_path` + `variant_aware_rescue`), so the pool can't be shared
across aligners without first making those fields per-task.

Confirmed by grep on 2026-05-22:

- **Per-aligner fields** (need to move to per-task): `bam_path`,
  `variant_aware_rescue`.
- **Shared fields** (stay at init): `genome` (loaded from `genome_path`),
  `polya_model` (loaded from `polya_model_path`), `annotated_junctions`
  (cross-aligner junction pool, same set for all aligners),
  `pool_chrom_index`, `gene_interval_trees`, all config flags
  (`apply_atract`, `ag_threshold`, `polya_trim`, `indel_correction`,
  `apply_3ss_rescue`, `dt_primed_cDNA`, `min_mapq`,
  `min_aligned_length`, `max_reads_for_variant_rescue`), `netseq_dir`.

So 2 fields per-aligner, ~14 fields shared. The refactor is much
smaller than v1 implied — closer to 100-150 LOC than 400-800.

## 1. C-II.0 — Worker state refactor

### Files touched

- `rectify/core/bam/parallel.py` — primary
- (Maybe) `rectify/core/commands/correct_command.py` — only if the call
  site changes shape

### Concrete shape

Today (`bam_parallel.py:64-89`):

```python
def _init_region_worker_state(genome_path, polya_model_path, worker_kwargs):
    global _REGION_WORKER_STATE
    _REGION_WORKER_STATE = dict(worker_kwargs)  # includes bam_path + rescue
    _REGION_WORKER_STATE['genome'] = load_genome(genome_path)
    _REGION_WORKER_STATE['polya_model'] = load_polya_model(...)

def _process_region_worker_from_state(region):
    return _process_region_worker(region=region, **_REGION_WORKER_STATE)
```

Tasks today are bare `(chrom, start, end)` tuples; per-aligner config is
baked into the worker state at init.

After C-II.0:

```python
def _init_region_worker_state(genome_path, polya_model_path, shared_kwargs):
    global _REGION_WORKER_STATE
    _REGION_WORKER_STATE = dict(shared_kwargs)  # only shared fields
    _REGION_WORKER_STATE['genome'] = load_genome(genome_path)
    _REGION_WORKER_STATE['polya_model'] = load_polya_model(...)

def _process_region_worker_from_state(task):
    # task = (region, bam_path, variant_aware_rescue)
    region, bam_path, variant_aware_rescue = task
    return _process_region_worker(
        region=region,
        bam_path=bam_path,
        variant_aware_rescue=variant_aware_rescue,
        **_REGION_WORKER_STATE,
    )
```

Call site in `process_bam_file_parallel`:

```python
# Before:
results_iter = pool.imap(_process_region_worker_from_state, regions)

# After:
tasks = [(r, bam_path, variant_aware_rescue) for r in regions]
results_iter = pool.imap(_process_region_worker_from_state, tasks)
```

Same change in `process_bam_streaming_parallel` if it follows the same
init pattern.

### Acceptance — C-II.0

- New test `tests/test_bam_parallel_state.py` exercises
  `process_bam_file_parallel` on the validation BAM with `n_threads=2`,
  captures the sorted output row hash, and asserts a recorded golden
  hash. Pre-refactor and post-refactor must produce identical hashes.
- Full fast suite (`pytest -m "not slow"`) passes (baseline 1239 + new
  tests). Two known unrelated WIP failures in `test_read_edits_reanchor`
  are now resolved by `55089f7`, so target is ~1241+ passing.
- No cluster run required for C-II.0.

### Stop-and-ask triggers — C-II.0

- If grep shows more than 3 fields varying between aligners (i.e., the
  shared/per-aligner split is wrong), surface and re-scope.
- If `process_bam_streaming_parallel` has a different worker init pattern
  that doesn't compose with C-II.0's design, drop streaming from C-II.0
  scope and handle it separately.
- If `n_threads=1` single-threaded fallback (line 683) has its own
  region-processing path that doesn't go through the worker state
  machinery — verify it still works post-refactor (it should, since it
  doesn't call `_init_region_worker_state`).

## 2. C-II.1 — Shared pool across aligners + work-stealing

### Files touched

- `rectify/core/commands/run/stages.py` — `_run_correction_per_aligner`
  becomes a coordinator that builds the full `(aligner, region)` task
  list and dispatches through a single shared pool.
- `rectify/core/bam/parallel.py` — add optional `executor` /
  `shared_pool` parameter to `process_bam_file_parallel` and
  `process_bam_streaming_parallel`. When provided, use the injected
  pool instead of creating one; skip the `initializer` argument since
  the shared pool's init must happen once at construction time.
- `rectify/core/commands/correct_command.py` — accept an optional
  injected pool through `args` or a new kwarg.
- New file (only if needed): `rectify/core/correct/region_dispatcher.py`
  for the shared queue + work-stealing harness. Avoid this if the
  coordination can live cleanly inside `stages.py`.

### Concrete shape

In `_run_correction_per_aligner`:

```python
n_workers = _resolved_ac  # already computed in Prep-C
ctx = _get_bam_worker_context()
shared_init_args = (genome_path, polya_model_path, shared_kwargs)

with ctx.Pool(n_workers,
              initializer=_init_region_worker_state,
              initargs=shared_init_args) as pool:
    futures_by_aligner = {}
    for aligner, bam_path in per_aligner_bams.items():
        # Build per-aligner task list. Tasks carry (region, bam_path,
        # variant_aware_rescue) — the per-aligner sliver.
        tasks = [(r, str(bam_path), per_aligner_rescue[aligner]) for r in plan]
        futures_by_aligner[aligner] = pool.imap_unordered(
            _process_region_worker_from_state, tasks
        )

    # Drain each aligner's results into its own writer. Work-stealing
    # is automatic — Pool drains the union of all submitted tasks in
    # any order.
    for aligner, results_iter in futures_by_aligner.items():
        write_per_aligner_outputs(aligner, results_iter)
        write_per_aligner_sidecar(aligner)
```

### Acceptance — C-II.1

- Existing `--aligner-concurrency 1` path remains byte-identical (no
  semantic change in serial mode).
- New `--aligner-concurrency auto` actually parallelizes — measurable
  via the sequential vs concurrent wall time of the test suite's
  validation pipeline.
- Cluster: Han wt_R1, 5 aligners, 16-core H2 node.
  **Target: correct phase wall ≤ slowest_aligner_alone × 1.3.**
  Submit via `h2-qsub`; poll every 5 min per CLAUDE.md "Cluster job
  monitoring" until past the first compute-heavy stage. If it fails
  the 1.3× gate, investigate work-stealing tuning before reverting.
- Re-run the heap-corruption smoke from
  `AGENT_FIXES.md "[2026-05-20] rectify correct pysam heap
  corruption ... STILL OPEN"`. C-II.0 changes the worker state model;
  C-II.1 changes the pool sharing pattern. Either could affect the
  heap-corruption signature. Don't claim Axis B safe-at-scale until
  the smoke clears.

### Stop-and-ask triggers — C-II.1

- If `process_bam_file_parallel`'s internal pre-pool logic (genome
  load at line 652, prescan at line 671) is per-aligner-stateful in
  ways that prevent injecting a shared pool, surface and re-design.
- If the `variant_aware` prescan adds a per-aligner first-pass step
  that itself uses parallelism, decide whether the shared pool also
  carries that prescan or whether prescans stay per-aligner-serial.
- If the 1.3× gate fails on H2, do NOT widen the gate — diagnose first
  (probably mapPacBio-tail starvation, fixable with finer region
  granularity for slow aligners).

## 3. Heap-corruption coordination

`AGENT_FIXES.md "[2026-05-20] rectify correct pysam heap corruption
on Han 2023 wt_R1 — STILL OPEN"` documents a heap corruption that
appears at ≥4M-read scale. The mitigation
(`--checkpoint-dir + threads=8`) silent-hung; the structural fix is
still open.

C-II.0 and C-II.1 both alter the worker state model. They may help
(workers now own their own bam_path → fewer chances for cross-task
state contamination) or hurt (any new state-sharing bugs). Plan to
re-run the smoke after C-II.0 lands (local pytest only — the bug
shows at cluster scale, so failure there must be tested on H2).

If the C-II.0 byte-equivalence test passes but the cluster smoke
still hangs, that's evidence the bug lives below the worker state
layer — probably in htslib state inheritance via fork (already
mitigated by `spawn` context at `_get_bam_worker_context` line 51,
but worth re-verifying).

## 4. Non-goals

- Do not change the Commit B manifest-only default.
- Do not touch the analyze pipeline (Commit D territory — landed at
  `8c9f80c`).
- Do not touch `restore_polya_command.py` (Commit D landed sidecar
  wiring there).
- Do not implement NETSEQ parallel path (Commit F).
- Do not implement Commit E shared-pool integration here — leave a
  TODO at E's integration point so E can land its non-pool work
  in parallel.
- Do not push without explicit user OK given the concurrent-session
  race history (see [[concurrent-sessions-git-race]]).

## 5. Sequencing relative to other commits

- **After Commit D (DONE at `8c9f80c`)** ✓
- **Before Commit E** — E uses the same shared-pool primitive. Doing
  E first would require a second shared-pool implementation in
  `align_command.py`, then merging the two later.
- **Cluster heap-corruption smoke** — re-run after C-II.0, again after
  C-II.1 cluster validation, before claiming Axis B safe-at-scale.

## 6. Open questions for Kevin / advisor

1. Should `process_bam_file_parallel`'s `variant_aware` prescan
   (line 660-679) remain per-aligner-serial, or also use the shared
   pool? Prescan touches the BAM differently from correction; sharing
   the pool might force two different worker init patterns.
2. Is `_run_correction_per_aligner`'s sort-step (today: pysam.sort on
   the input BAM when no .bai exists, line 388-391 of `stages.py`)
   safe to run concurrently across aligners, or should it stay serial
   before the shared pool spins up? It's per-aligner so concurrent
   should be safe, but worth confirming.
3. M1 policy: should C-II.1 actually run on M1 (`auto → 1` keeps it
   serial per Prep-C's resolve_aligner_concurrency), or stay strictly
   sequential on M1 forever and only exercise the shared queue on
   cluster nodes? Mac BLAS/OpenMP libraries have shown contention
   issues with mp.Pool in the past.
