# Spec — Default-parallel RECTIFY across DRS, ONT cDNA, and NETSEQ/QSrev

**Status:** Design — not yet implemented. Commit Zero (profile + audit) gates everything.
**Owner:** Kevin R. Roy (M1 orchestrator); per-commit execution split between Opus 4.7 (design / judgment) and Sonnet 4.6 subagents (mechanical work). See §11 for the delegation policy.
**Branch:** `drs-validation-rebuild` (do not target `master` — frozen at 0.9.0).
**Created:** 2026-05-19
**Supersedes / depends on:**
- Triggering issue: `project_status_markdowns/ISSUE_parallel_bam_write_and_lazy_tsv.md`
- Triggering observation: Han wt_R1 on H2 (12.49M reads, 16 cores, `pod_smp.q@n1821`, job 13446575) showed `rectify correct` finishing parallel correction in ~5 min then sitting another ~2.5 hr in single-threaded `write_corrected_bam` at 93.8% CPU on one core, 0% I/O wait.
- Related memory: `feedback_small_chunks_owners_preemption.md`, `feedback_m1_memory_discipline.md`, `feedback_provenance_alongside_outputs.md`, `feedback_sherlock_sbatch_set_u.md`.

---

## 0. Goal

Make region-parallel BAM write + per-region TSV emission the **default** path across all three pipeline arms (DRS, ONT cDNA, NETSEQ/QSrev), add a second parallelism axis (concurrent aligners via dynamic work queue), and restructure `analyze` to consume per-region partials so the merge point is pushed all the way down to genuinely-global substeps (cluster calling + DESeq2 count matrix).

**Companion requirement — stage-level resume.** Every stage in `rectify run-all` emits a PROVENANCE.json sidecar on completion that captures inputs (with sha256), outputs (with sha256), config args, rectify git_sha, and exit status. A subsequent `run-all` invocation reads each stage's prior sidecar, validates it against the current inputs + config, and skips the stage if and only if all checks pass. The pipeline is resumable from any interruption point top-to-bottom without recomputing already-valid stages. See §6.5 for schema + skip-check semantics; see Commit A.5 for the resume infrastructure and per-stage commits for the wiring.

**Success metric (PROTOCOL-CONDITIONAL — see `dev/specs/profile_results.md` §7).** Commit Zero profiling on 2026-05-19 found that the parallelism payoff differs sharply by protocol:

- **QSrev / short-read / dT-primed cDNA:** end-to-end Han wt_R1 (12.95M reads) on a 16-core H2 node drops from `align + correct + analyze ≈ 5+ hr` to **under 75 min** on a clean run. Driver: region-parallel `write_corrected_bam` (which today runs single-threaded and accounts for ~90% of wall time at -j 16; `realign_exon_blocks` is 88% of WRITE phase).
- **DRS / ONT-cDNA / long-read:** the dominant cost (`realign_exon_blocks` inside the already-parallel correct-region loop, ~91% of wall at -j 1) is NOT meaningfully reduced by this refactor's region-parallel write infrastructure. Expected speedup from Commits A-F on DRS is ~1.1-1.3× on wall time (driven primarily by Axis B concurrent aligners + small write-phase wins). The architectural payoff for DRS is resume + provenance + analyze partial-streaming, not raw wall time.

No correctness regression on the 934-test suite. Byte-equivalent (modulo deterministic sort order) corrected BAM/TSV outputs compared to the current pipeline for both protocols.

**Secondary success metric (both protocols):** a `run-all` invocation that follows a successful prior run on the same inputs completes in under 60 seconds total (all stages skipped via sidecars). A `run-all` resumed from an interruption after stage K runs only stages K+1..N.

**Out-of-scope follow-up (noted for future work, NOT this PR):** for QSrev, `realign_exon_blocks` runs TWICE — once in correct-region (12.9% of correct phase) and again inside write_corrected_bam (88% of write phase) on the same reads. Caching the post-realign CIGAR in the correct-phase TSV and skipping the second realign in write would collapse QSrev WRITE-phase work by ~50× and would compose with the parallelism refactor. For DRS, no equivalent quick win — the dominant cost IS the first realign pass.

**Non-goals (out of scope for this PR):**
- GPU acceleration of any aligner.
- Streaming pysam writes via Rust/C extensions.
- Replacing `merge_corrected_tsvs` consensus logic — only parallelizing it.
- Cluster scheduler integration changes (sbatch / qsub templates stay as they are).
- A new public release version bump. Land on `drs-validation-rebuild`; release decision happens later.

---

## 1. Architectural overview

### 1.1 Two parallelism axes

| Axis | Scope | Implementation | Bound |
|---|---|---|---|
| **A — regions** | Within a single `rectify correct` invocation; split the input BAM into gap-delimited coord regions | Reuse existing `get_processing_regions()` from `rectify/core/bam/parallel.py:244-458` (extract to `rectify/core/bam/regions.py`) | CPU cores |
| **B — aligners** | Across the 5 per-aligner correct passes inside `_run_correction_per_aligner` (`stages.py:294-415`) | Single shared `multiprocessing.Queue` of `(aligner_id, region_id)` task tuples; N total workers pull dynamically | CPU cores × RAM |

Axes compose: total active workers = `min(total_cores - reserve, sum of region tasks across aligners)`. Static partitioning (`n_workers_per_aligner = total_cores // n_aligners`) was rejected because minimap2 finishes 3-5× faster than deSALT/gapmm2 — static splits leave fast-aligner workers idle. The dynamic queue keeps all N workers saturated.

### 1.2 Merge points (latest possible)

The merge point for each output type is pushed as far downstream as physically necessary:

| Output | New merge point | Reason |
|---|---|---|
| Per-region corrected BAM | `samtools merge` of coord-sorted region BAMs, run once at end of `rectify correct` | Final BAM is a single artifact for IGV / downstream tooling |
| Per-region corrected TSV | **Manifest only by default** (`corrected_reads.manifest.tsv` listing region TSVs); single concatenated TSV produced only when `--emit-merged-tsv` is set (default ON during back-compat window — see §6) | Downstream `analyze` reads manifest; external scripts that want one file still get it |
| Per-aligner TSVs → consensus | `merge_corrected_tsvs` runs once after all aligners complete; may itself be region-parallel (decision deferred to Commit Zero profile) | Consensus calling is per-position; can be parallelized by chrom |
| Bedgraphs (per-condition) | Per-chrom streaming; concat at end | No global view needed |
| Correction stats, histograms | Per-region stats → reduce-sum | Trivially merge-summable |
| Position index (tabix) | Per-chrom slices → final tabix on concatenated bed | Tabix requires sorted input but tolerates concatenated per-chrom shards |
| Cluster calling | Two-pass: (1) global coverage quantiles streamed from manifest in constant RAM; (2) per-chrom cluster calling with pre-computed thresholds, parallelized via `ProcessPoolExecutor` | Adaptive thresholds genuinely require global view; cluster calling itself does not |
| DESeq2 count matrix | Runs once on collapsed cluster table after cluster calling | Inherently global; small input (~10⁴ rows) so single-threaded is fine |

### 1.3 Sort-then-merge (not merge-then-sort)

Each region worker emits a **coord-sorted** region BAM (per-region sort is parallel-trivial since regions are gap-delimited and small). The final `samtools merge` with `--threads N` heap-merges pre-sorted inputs in O(N log K) and produces a fully-sorted output directly. **No final global sort pass is needed.** This is what `samtools merge` was designed for.

This corrects the original issue doc's "emit in `imap` order → merge → cheap final sort" — the final sort on a 12.5 GB BAM is not cheap (5-10 min single-threaded).

### 1.4 Shared infrastructure under `rectify/core/bam/`

New modules:

- `rectify/core/bam/regions.py` — extends existing 148-line file. Adds `RegionPlan` dataclass (`region_id`, `chrom`, `start`, `end`, `expected_read_count`, `tmp_dir`). Factors `get_processing_regions()` out of `parallel.py` so `bam_writer_parallel`, NETSEQ, cDNA, and `variant_scan` can all consume it.
- `rectify/core/bam/bam_writer_parallel.py` — `write_corrected_bam_parallel(input_bam, corrections_table, output_bam, n_threads, genome, tmp_dir, *, allow_resume=True)`. Workers apply per-read mutations (the same set as today's `bam_writer.py:200-368`) and emit coord-sorted region BAMs with `.ok` sentinels for resume support.
- `rectify/core/bam/tsv_partition.py` — `RegionTsvWriter` class + `corrected_reads.manifest.tsv` schema.
- `rectify/core/bam/netseq_bam_parallel.py` — **duplicated** (~150 LOC) from the corrected-bam path. NETSEQ produces per-base 3'-end counts, has no `realign_exon_blocks`, has its own exclusion-stats path; per CLAUDE.md "don't introduce abstractions beyond what the task requires", duplicating is preferred to generalizing across three arms with different output shapes.

New test utility:

- `tests/utils/bam_compare.py` — QNAME-sort both BAMs, compare record-by-record on `(CIGAR, FLAG, POS, MAPQ, all tags)`. Reused across DRS, cDNA, NETSEQ equivalence tests.

---

## 2. Pipeline-stage audit — what gets parallelized

This table is the **definitive list of stages this PR touches.** Commit Zero (§3) measures each on Han wt_R1 to confirm the bottleneck attribution before infrastructure work begins.

| # | Stage | Current code | Today's parallelism | New parallelism | Commit |
|---|---|---|---|---|---|
| 1 | DRS poly-A trim | `drs_trim_command.py:trim_drs_bam_polya` | Single-threaded fetch loop | Region-parallel via `regions.py` | E (if profile confirms ≥10% of wall time) |
| 2 | Prescan: variant scan | `bam/variant_scan.py:run_variant_aware_scan` (line 54-65 fetch loop) | None. `--threads` flag at `prescan_command.py:113-118` is only for genome load | Region-parallel; merge = per-position frequency sum on a `VariantAwareHomopolymerRescue.merge(other)` new method | D |
| 3 | Prescan: junction pool | `splice/junction_scoring.py:200-236` `for bam_path in aligner_bams:` | Sequential across 5 aligner BAMs | `ProcessPoolExecutor(max_workers=n_aligners)`; merge = `set.union` | A (~30 LOC freebie alongside infrastructure) |
| 4 | Multi-aligner `align` | `stages.py:_run_alignment` | Verify in Commit Zero — currently sequential? | If sequential, parallel-across-aligners via shared queue | E (only if confirmed as bottleneck) |
| 5 | Per-aligner correct: region loop | `bam/parallel.py:process_bam_file_parallel` / `process_bam_streaming_parallel` | Already parallel (`mp.Pool(n_threads)`) | Switch TSV emission to per-region writes + manifest (no algorithmic change) | A + B |
| 6 | Per-aligner correct: `write_corrected_bam` | `bam/bam_writer.py:200-368`, called from `correct_command.py:848` | **Single-threaded — the bottleneck** | Region-parallel via new `write_corrected_bam_parallel` | A + B |
| 7 | `restore_polya` BAM rewrite | `correct_command.py:789` | Single-threaded soft-clip rewrite | Region-parallel via `bam_writer_parallel` infra (same shape) | B (near-free if same machinery) |
| 8 | Per-aligner correct across 5 aligners | `stages.py:_run_correction_per_aligner` `:294-415` | Sequential | Dynamic `(aligner_id, region_id)` work queue | C |
| 9 | `merge_corrected_tsvs` | `consensus/corrected_consensus.py:577` | Single-threaded; loads all 5 aligner TSVs into RAM | Per-chrom `ProcessPoolExecutor`; per-chrom consensus tables → final concat | D (if profile confirms it's the next wall) |
| 10 | Analyze: bedgraphs, stats, position index, histograms | `commands/analyze_command.py:95` (`load_corrected_positions`) | Loads single concatenated TSV; mostly single-threaded substeps | Per-chrom streaming from manifest; reduce at end | D |
| 11 | Analyze: cluster calling | `cluster_cpa_sites_adaptive` | Single-threaded global pass | Two-pass: streamed coverage quantiles → per-chrom cluster calling via pool | D |
| 12 | Analyze: DESeq2 | `build_cluster_count_matrix` → `run_deseq2_*` | Single-threaded (R subprocess) | Unchanged — inherently global, runs on small input | — |
| 13 | ONT cDNA stage-1 consensus FASTQ | `cdna_correct_command.py:run`, `cdna/io.py:write_stage1_fastq` line 145 | Single-threaded cluster loop | `ProcessPoolExecutor` over UMI clusters; per-shard FASTQ concat at end | E |
| 14 | ONT cDNA post-align correct | Reuses `correct_command.py` | Inherits §5-§7 wins automatically | (transparent) | — |
| 15 | NETSEQ 3'-end count | `commands/netseq_command.py:run_netseq`, `bam/netseq_bam_processor.py` (640 LOC) | Single-threaded fetch loop | Region-parallel via duplicated `netseq_bam_parallel.py` | F |

**Stages explicitly NOT in scope:**
- Index BAM (`pysam.index`) — fast, already runs once at end.
- Genome load — fast (~1-2 s for yeast); already cached across stages.
- Calibration table loads (penalty tables, junction overhang) — fast, run once.

---

## 3. Commit Zero — Profile + audit (blocking gate)

**No code changes.** Output: `/Users/kevinroy/work/rectify/dev/specs/profile_results.md`.

### 3.1 Procedure

On H2 (preferred over Sherlock for agent autonomy — no 2FA):

```bash
ssh hoffman2 'bash -lc "
  cd /u/home/k/kevinroy/scratch/profile_runs
  module load conda/23.11.0 && conda activate rectify
  # 200k-read slice of Han wt_R1 — already staged under han2023_chunks_v2/chunk_001.bam
  python -m cProfile -o write_bam.prof -m rectify correct \
    /u/project/guillom/shared/processed/han2023/chunks_v2/wt_R1.chunk_001.bam \
    --Scer --output-dir profile_runs/write_bam_test \
    --n-threads 1 --no-parallel  # force single-threaded path for clean attribution
  python -c \"import pstats; p = pstats.Stats('write_bam.prof'); p.sort_stats('cumulative').print_stats(40)\"
"'
```

### 3.2 Stage timings to capture (all on Han wt_R1, 12.49M reads, 16-core H2 node)

Wall time and CPU% per stage. Each measured via `time` + `top -bn1 -d 5 -p $PID` snapshots written to a log file (existing `h2-queue-probe` skill captures the pattern).

1. `drs_trim_bam_polya` — full run.
2. `prescan` Step 1 (variant scan) — full run.
3. `prescan` Step 2 (junction pool) — full run.
4. `_run_alignment` — per-aligner timings.
5. Per-aligner `correct` region phase — already known (~5 min). Re-confirm.
6. Per-aligner `write_corrected_bam` — already known (~2.5 hr). Re-confirm.
7. `merge_corrected_tsvs` — full run.
8. `analyze` substep timings.

### 3.3 cProfile decomposition of `write_corrected_bam`

**STATUS: DONE 2026-05-19.** See `dev/specs/profile_results.md` for the full attribution. Summary:

- QSrev WRITE phase: REALIGN 88.0%, BGZF 1.9%, SET_TAG 0.2% — `realign_exon_blocks` dominates as predicted.
- DRS WRITE phase: only 0.5% of total wall (CORRECT_REGION takes 98.4%). The bottleneck for DRS lives upstream and is already region-parallel via `process_bam_file_parallel`.

**Refuted escape hatches (do NOT implement):**
- `compresslevel=1` Commit Zero-bis: BGZF is 1.9% (QSrev) of WRITE, not the proposed ≥30%. Skip.
- Batched `set_tag`: SET_TAG is 0.2% of WRITE. Skip.

**Confirmed decision:** proceed with Commits A-F as written. Protocol-specific speedup expectations are in §0.

### 3.4 Acceptance criteria for Commit Zero (Phase Zero-A: DONE 2026-05-19)

- [x] `profile_results.md` exists with bucket attribution for QSrev + DRS comparison (Phase Zero-B full-pipeline timing still pending — confirmation only, not a gate).
- [x] py-spy native+raw decomposition of `write_corrected_bam` into the four buckets (§3.3).
- [x] Go/rescope decision recorded: **PROCEED with Commits A-F**, with protocol-conditional success metric.
- [x] This spec updated to reflect the protocol-conditional outcome.

### 3.5 Estimated wall time

- Profile run on 200k slice: ~5 min.
- Full-sample timing run on Han wt_R1: ~3-4 hr (this is the bottleneck the PR is fixing — one full slow pass for baseline).
- Analysis + writeup: ~30 min Opus.
- **Total: ~4 hr cluster wall + ~30 min agent time.**

---

## 4. Per-commit plan

Each commit is independently mergeable and passes the existing `pytest -m "not slow"` suite (~1 min, 934 tests). Slow markers (cDNA smoke + chain canary) pass at the end of each commit that touches their code path.

### Commit A — Shared infrastructure (no behavior change)

**Goal:** introduce the new modules without changing any default code path. Existing tests pass unchanged.

**Files:**
- **Modify** `rectify/core/bam/parallel.py` — extract `get_processing_regions()` to `regions.py`, import it back. No behavior change.
- **Modify** `rectify/core/bam/regions.py` — add `RegionPlan` dataclass + `get_processing_regions()`.
- **New** `rectify/core/bam/bam_writer_parallel.py` — `write_corrected_bam_parallel()` with sort-then-merge, idempotent `.ok` sentinels for preemption resume, and an internal `_process_region_for_bam_write()` worker that mirrors `bam_writer.py:200-368` mutations per-region.
- **New** `rectify/core/bam/tsv_partition.py` — `RegionTsvWriter` + manifest schema.
- **Modify** `rectify/core/splice/junction_scoring.py` — `build_junction_pool` swaps the `for bam_path in aligner_bams:` loop for `ProcessPoolExecutor`; result is the same set union. ~30 LOC change. Test added.
- **New** `tests/utils/bam_compare.py` — QNAME-sort + record-by-record comparison helper.
- **New** `tests/test_bam_writer_parallel_smoke.py` — calls `write_corrected_bam_parallel()` on a 1k-read fixture and confirms output matches single-threaded `write_corrected_bam` via `bam_compare.py`.

**Acceptance:**
- [ ] `pytest -m "not slow"` passes (no regressions).
- [ ] New parallel-write smoke test passes.
- [ ] No existing code path is rerouted to the new functions yet.
- [ ] Diff reviewed by Opus before merge to `drs-validation-rebuild`.

**Estimated LOC:** ~500 new, ~50 modified.
**Delegation:** Sonnet subagent (mechanical work, well-specified). Opus reviews diff.

---

### Commit A.5 — Provenance + resume infrastructure (no stage wiring yet)

**Goal:** introduce the provenance/sidecar/skip-check module and a `tracker` migration decision, without wiring any stage through it. Existing run-all behavior unchanged. This lands before Commit B so that B can be the first commit to actually emit a stage-level sidecar.

**Files:**
- **New** `rectify/core/provenance/__init__.py` — public API: `ProvenanceRecord`, `write_stage_sidecar`, `read_stage_sidecar`, `can_skip_stage`, `SkipDecision`.
- **New** `rectify/core/provenance/sidecar.py` (~250 LOC) — `ProvenanceRecord` dataclass matching §6.5.2 schema; `write_stage_sidecar(stage, sample_output, inputs, outputs, sub_outputs, argv, stats, ignore_argv, ignore_input_roles)` does atomic write via `tmp + rename`; `read_stage_sidecar` validates schema version + returns parsed record.
- **New** `rectify/core/provenance/skip_check.py` (~200 LOC) — `can_skip_stage(...)` implementing §6.5.3; `SkipDecision` is a small dataclass with `.skip: bool`, `.reason: str`, `.prior: Optional[ProvenanceRecord]`.
- **New** `rectify/core/provenance/hashing.py` (~80 LOC) — `sha256_of_file(path)` with mmap fast-path for large files; `normalized_config_hash(argv_template, ignore)` for stable argv hashing across runs.
- **New** `rectify/core/provenance/path_resolver.py` (~200 LOC) — `PortablePath` dataclass (§6.5.2.5), `from_path(p, sample_output_dir)` tokenization, `resolve_now(sample_output_dir)` resolution with fallback-to-cached-absolute, `tokenize_argv_paths(argv, sample_output)` for invocation tokenization, `PathResolutionError` exception. Refuses to construct a non-transient `env_relative` PortablePath with `env_var ∈ {L_SCRATCH, TMPDIR}` (enforces the §6.5.3 invariant).
- **New** `rectify/core/provenance/cluster.py` (~50 LOC) — `detect_current_cluster()` returns `"sherlock" | "hoffman2" | "m1" | "other"` per §6.5.2.5 heuristics. Used by sidecar writer and skip-check.
- **New** `rectify/core/provenance/cli.py` (~60 LOC) — shared argparse helpers `add_resume_args(parser)` that adds `--force-all`, `--force-stage`, `--accept-prior-provenance`, `--dry-run-resume` to any stage parser. Each command's `create_*_parser` calls this helper.
- **Modify** `rectify/core/commands/run/helpers.py` (or wherever `tracker` is defined — verify with `grep -rn "register_staged" rectify/`) — apply the Commit A.5 prereq decision (replace / extend / coexist). Document the decision in the commit message.
- **New** `tests/test_provenance_sidecar.py` (~150 LOC) — round-trip: write a sidecar, read it back, mutate one field, check `can_skip_stage` returns RUN with the right reason. Includes a test for each of the 4 RUN-reasons (missing output, sha256 mismatch on input/output, config diff, git_sha mismatch).
- **New** `tests/test_provenance_skip_check.py` (~100 LOC) — table-driven test of `SkipDecision` outcomes.
- **New** `tests/test_portable_path.py` (~200 LOC) — covers: (a) sample-relative tokenization with symlinks, (b) env_relative tokenization ordering (`$L_SCRATCH` before `$SCRATCH` before `$HOME`), (c) `cached_absolute` fallback when env var unset on read side, (d) cross-cluster resolution (set `$OAK=/tmp/fake_oak`, write sidecar, unset `$OAK`, re-resolve — should fall back to cached_absolute and warn), (e) `ValueError` on attempt to construct non-transient `env_relative` PortablePath with `env_var ∈ {L_SCRATCH, TMPDIR}`, (f) `argv_template` tokenization round-trip preserves `config_hash` across runs with different `$L_SCRATCH` values.
- **Modify** `dev/specs/parallel_refactor_plan.md` (this file) — mark Commit A.5 status DONE on landing.

**Acceptance:**
- [ ] All new unit tests pass — including the 6 path-portability scenarios in `test_portable_path.py`.
- [ ] `pytest -m "not slow"` passes (no regressions).
- [ ] `rectify run-all --help` displays the new resume flags.
- [ ] No existing stage emits or consumes sidecars yet (Commit B and later wire them in).
- [ ] `tracker.register_staged` decision recorded in commit message + this spec.
- [ ] PortablePath refuses to construct a non-transient `env_relative` path with `env_var ∈ {L_SCRATCH, TMPDIR}` (enforces the §6.5.3 invariant that persistent outputs live in `$SAMPLE_OUTPUT`).
- [ ] `detect_current_cluster()` correctly identifies H2, Sherlock, and M1 in unit tests (mock `os.uname` + hostname env).
- [ ] Opus reviews diff — especially the `path_resolver` resolution-order logic and the invariant enforcement.

**Estimated LOC:** ~950 new, ~30 modified.
**Delegation:** Sonnet for the boilerplate sidecar/hashing/cli modules; Opus owns the `tracker` migration decision (requires reading existing code + judgment), the skip-check decision tree implementation (correctness-critical), and the `path_resolver` design (cross-cluster correctness traps).

---

### Commit B — DRS default flip + back-compat decision

**Goal:** wire `correct_command.py` through the new infra; flip the default to parallel write; keep escape hatches.

**Prerequisite:** back-compat grep for `corrected_reads.tsv` consumers (Opus does this before Sonnet starts):
```bash
cd /Users/kevinroy/Library/CloudStorage/GoogleDrive-kevinrjroy@gmail.com/My\ Drive/Work/Chanfreau\ Lab
grep -rn "corrected_reads\.tsv" scripts/ methods/ h2_review/ rectify/dev/ --include="*.py" --include="*.sh" --include="*.md" | head -50
ssh hoffman2 'grep -rn "corrected_reads\.tsv" /u/project/guillom/shared/processed/han2023/ 2>/dev/null | head -20'
```
The grep result determines the default for `--emit-merged-tsv`:
- If ≤5 external consumers and all easily migrated → manifest-only by default; add `--emit-merged-tsv` opt-in.
- If >5 consumers or any are external (collaborator scripts) → keep `--emit-merged-tsv` ON by default; flip in a later PR after consumers migrate.

**Files:**
- **Modify** `rectify/core/commands/correct_command.py`:
  - `:848` replace `write_corrected_bam` call with `write_corrected_bam_parallel` (default path).
  - Add `--legacy-single-threaded` escape hatch routing to old code.
  - Add small-input gate: if input BAM < 100 MB OR estimated reads < 500k, route to sequential path automatically (prevents test-suite regression from ProcessPool spawn cost).
  - `:670`/`:727` — switch TSV emission to `RegionTsvWriter`.
  - `:789` `restore_polya` BAM rewrite — also routes through `bam_writer_parallel` infra (~free with the new machinery).
- **Modify** `rectify/core/commands/analyze_command.py:95` — detect manifest vs single-TSV input; load accordingly.
- **Modify** `rectify/core/bam/loaders.py` — manifest-aware `load_corrected_positions`.
- **New** `tests/test_correct_command_parallel_default.py` — runs DRS smoke through new default path; compares output to legacy via `bam_compare.py` + `tsv_compare`.

**Acceptance:**
- [ ] Han wt_R1 end-to-end smoke (`scripts/smoke/han_wt_R1_drs_e2e.sh` or equivalent) completes in target wall time (< 75 min on 16-core H2 node).
- [ ] Output BAM is QNAME-sort-equivalent to legacy.
- [ ] Output `corrected_reads.tsv` (concatenated from manifest if `--emit-merged-tsv` is on) is row-set-equivalent to legacy (sort columns, then `cmp`).
- [ ] `pytest -m "not slow"` passes.
- [ ] `pytest -m slow` cDNA smoke passes.
- [ ] `--legacy-single-threaded` flag produces identical output to current main branch.
- [ ] Small-input gate verified: a 50k-read test BAM does NOT spawn workers (instrument with a sentinel log message).
- [ ] Per-region PROVENANCE.json sidecars present for region BAMs.
- [ ] **Stage-level sidecars emitted** for `correct` and `restore_polya` stages (per §6.5.2 schema): inputs (BAM + genome + annotation + prescan caches), outputs (final corrected BAM + index + manifest), sub_outputs (per-region BAMs).
- [ ] **Resume verified:** rerun the smoke after a successful run; resume layer skips both stages with `decision=SKIP`; total wall time < 30 seconds.
- [ ] **Resume verified after input change:** touch the input BAM (re-hash); rerun; both stages rerun with `reason=input bam sha256 changed`.
- [ ] **Resume verified after argv change:** rerun with a different `--n-threads` (should SKIP — listed in `ignore_argv`); rerun with a different `--min-coverage` (should RUN — affects output).
- [ ] Opus reviews diff + smoke results + resume traces before merge.

**Estimated LOC:** ~250 new, ~170 modified (resume wiring adds ~50 LOC vs. original estimate).
**Delegation:** mostly Sonnet, but Opus owns:
- The back-compat grep + default decision.
- Interpreting any divergence the equivalence test surfaces (real bug vs. acceptable sort-tie artifact).
- The first stage-level sidecar wiring (sets the pattern for D, E, F).

---

### Commit C — Dynamic aligner-queue concurrency

**Goal:** add the second parallelism axis. Multiple aligners' regions processed concurrently via a single shared work queue.

**Files:**
- **Modify** `rectify/core/commands/run/stages.py:_run_correction_per_aligner` `:294-415`:
  - Replace per-aligner serial loop with a `multiprocessing.Queue` of `(aligner_id, region_id, region_plan)` tuples and N workers pulling dynamically.
  - N = `total_cores - reserve` (reserve = 2 for main + housekeeping).
  - Each worker: pop task → run correct for that aligner's region → emit region BAM + region TSV in that aligner's subdir.
  - When all tasks for an aligner finish, that aligner's `write_corrected_bam_parallel` merge + final TSV manifest writes are kicked off in a follow-up small thread pool (one merge thread per aligner; doesn't compete with the main work queue).
- **New** `rectify/core/utils/resources.py` — `detect_machine_class()` returns one of `{"m1_laptop", "cluster_node"}` based on `total_mem_gb < 16 AND os.cpu_count() < 8`. Used to disable Axis B on M1 (preserves `feedback_m1_memory_discipline.md` constraint).
- **Modify** `_run_correction_per_aligner` to honor `--aligner-concurrency {auto,1,N}`. `auto` picks based on `detect_machine_class()`.

**Acceptance:**
- [ ] On 16-core H2 node: Han wt_R1 5-aligner correct phase runs in wall time ≤ (slowest_aligner_alone × 1.3). The 30% slack accounts for queue overhead.
- [ ] On M1 (autodetected): falls back to sequential aligners, full N region workers each. Verified by checking the active process count peaks at N+1, not 5N+1.
- [ ] No OOM on 64 GB H2 node with 5 aligners concurrent.
- [ ] `--aligner-concurrency 1` produces output identical to Commit B (sequential aligners).
- [ ] Opus reviews diff + resource autodetect logic.

**Estimated LOC:** ~250 new, ~50 modified.
**Delegation:** Opus owns this commit. Boundary-condition heavy (RAM autosizing, queue lifecycle, M1 vs cluster), hard to specify completely up front for a Sonnet agent.

---

### Commit D — Analyze partial-streaming + merge_corrected_tsvs parallelism

**Goal:** push the merge point all the way down. Most analyze substeps run on per-region partials.

**Conditional on Commit Zero results:** if `merge_corrected_tsvs` shows up as <10% of wall time, defer to a later PR and keep this commit smaller (analyze-only). If it's a major fraction, include both.

**Files:**
- **New** `rectify/core/analyze/streaming.py` — `load_corrected_positions_partitioned(manifest_path) → Iterator[ChromShard]` + per-substep helpers (`run_bedgraph_per_shard`, `accumulate_position_index_per_shard`, etc.).
- **Modify** `rectify/core/commands/analyze_command.py`:
  - Detect manifest input → use streaming loader.
  - Bedgraphs: per-chrom workers → final concat.
  - Stats: per-region dicts → reduce-sum.
  - Position index: per-chrom tabix shards → final concat-then-tabix.
  - Histograms: per-shard bins → sum.
  - Cluster calling: two-pass. Pass 1 streams from manifest to compute global coverage quantiles in constant RAM (no positions dataframe loaded). Pass 2 calls clusters per-chrom via `ProcessPoolExecutor` with pre-computed thresholds.
  - DESeq2: unchanged.
- **Modify** `rectify/core/consensus/corrected_consensus.py:merge_corrected_tsvs` — split per-chrom, run consensus calling in `ProcessPoolExecutor`, final concat. Only if Commit Zero flags it.
- **Modify** `rectify/core/commands/cdna_analyze_command.py` — share the same streaming loader.
- **New** `tests/test_analyze_streaming_equivalence.py` — runs analyze through both paths on a small fixture; compares all outputs (bedgraphs, cluster table, stats, position index, histograms).

**Acceptance:**
- [ ] Analyze on Han wt_R1 outputs match legacy outputs (file-by-file diff).
- [ ] Peak RAM during analyze cluster calling stays under 8 GB (versus current ~16-20 GB on Han wt_R1).
- [ ] Total analyze wall time reduced by ≥3× on a 16-core node.
- [ ] `pytest -m "not slow"` passes.
- [ ] **Stage-level sidecars emitted** for `merge_corrected_tsvs` (if in scope) and `analyze`. Analyze sidecar's outputs list: bedgraphs (per-condition), cluster table, DESeq2 results, position index, 3′ histograms.
- [ ] **Resume verified:** rerun analyze immediately after a successful run; SKIP with wall time < 5 seconds.
- [ ] **Resume verified after correct sidecar invalidation:** modify `corrected_reads.manifest.tsv` (re-hash); rerun run-all; both `merge_corrected_tsvs` and `analyze` rerun because their input sha256 changed.
- [ ] Opus reviews diff + cluster-calling two-pass logic + equivalence outputs + resume traces.

**Estimated LOC:** ~400 new, ~170 modified.
**Delegation:** Opus owns the cluster-calling two-pass refactor and equivalence interpretation. Sonnet owns the per-substep streaming helpers, the analyze sidecar wiring (follows Commit B pattern), and the test fixture.

---

### Commit E — DRS poly-A trim + multi-aligner align + ONT cDNA stage-1

**Goal:** catch the remaining stages flagged by Commit Zero as bottlenecks. Each is region-parallel (or cluster-parallel for cDNA stage-1) in shape.

**Conditional on Commit Zero results.** Skip any sub-stage whose measured wall time is < 10% of pipeline total.

**Files:**
- **Modify** `rectify/core/commands/drs_trim_command.py:trim_drs_bam_polya` — region-parallel via `regions.py`. Output is a trimmed BAM + metadata parquet; both shard-then-merge naturally.
- **Modify** `rectify/core/commands/run/stages.py:_run_alignment` — if today's code is sequential across aligners (verify in Commit Zero), wrap in `ProcessPoolExecutor(max_workers=n_aligners)`. Each aligner is its own subprocess (samtools / mapPacBio / etc.), so this is process-level concurrency, not Python-thread concurrency.
- **Modify** `rectify/core/commands/cdna_correct_command.py` + `rectify/core/cdna/io.py:write_stage1_fastq` line 145 — `ProcessPoolExecutor` over UMI clusters; per-shard FASTQ files; concat at end via `cat`. UMI clusters are independent (no cross-cluster state) so this is trivially correct.
- **New** `tests/test_cdna_stage1_parallel_equivalence.py` — verifies stage-1 FASTQ output is identical (sorted by cluster ID) between sequential and parallel paths.

**Acceptance:**
- [ ] DRS poly-A trim wall time reduced ≥4× on 16-core node.
- [ ] Multi-aligner align wall time reduced to ≈max(per_aligner_time).
- [ ] cDNA stage-1 FASTQ wall time reduced ≥4× on 16-core node.
- [ ] All output files identical (modulo cluster order) to legacy.
- [ ] cDNA smoke test passes.
- [ ] **Stage-level sidecars emitted** for each newly-parallelized stage: `drs_trim`, `align` (one per aligner, or one combined sidecar with sub_outputs — Opus decides at commit time), `cdna_stage1_consensus`.
- [ ] **Resume verified:** rerun run-all after a successful run; every stage in this commit SKIPs.
- [ ] Opus reviews diff.

**Estimated LOC:** ~250 new, ~120 modified.
**Delegation:** Sonnet for each sub-piece; Opus reviews and gates on Commit Zero scope. Opus decides align sidecar shape (one-per-aligner vs combined).

---

### Commit F — NETSEQ duplicated parallel path

**Goal:** apply the region-parallel pattern to NETSEQ via duplication (not abstraction).

**Files:**
- **New** `rectify/core/bam/netseq_bam_parallel.py` — ~150 LOC. Mirrors `bam_writer_parallel.py` shape but emits per-base 3'-end count TSVs and handles NETSEQ exclusion stats. Shares only `RegionPlan` + `get_processing_regions()` from `regions.py`.
- **Modify** `rectify/core/commands/netseq_command.py:run_netseq` — route default through `netseq_bam_parallel`. Add `--legacy-single-threaded` escape.
- **New** `tests/test_netseq_parallel_equivalence.py` — runs both paths on a small fixture; compares 3'-end count TSVs + exclusion stats.

**Acceptance:**
- [ ] NETSEQ on Han 2023 cluster_revisit test sample completes in wall time ≤ legacy/4.
- [ ] Output TSVs identical (sorted by chrom, pos, strand).
- [ ] Exclusion stats identical.
- [ ] `pytest -m "not slow"` passes.
- [ ] **Stage-level sidecar emitted** for the NETSEQ stage.
- [ ] **Resume verified:** rerun after a successful run; SKIP with wall time < 5 seconds.

**Estimated LOC:** ~180 new, ~60 modified.
**Delegation:** Sonnet. Opus reviews diff.

---

### Commit F.5 — Run-all resume end-to-end integration tests

**Goal:** verify that the resume infrastructure works correctly across the full pipeline with failure injection. This is the commit that proves the §6.5 design end-to-end before Commit G publishes it.

**Prerequisite:** Commits A.5 + B + D + E + F have all landed (every stage emits a sidecar).

**Files:**
- **New** `tests/integration/test_run_all_resume_drs.py` — uses a small DRS fixture (~50k reads, runs in ~30 s end-to-end). Test matrix:
  1. **Clean run → resume:** run-all to completion; rerun-all; every stage SKIPs; total wall < 10 s.
  2. **Interrupt at each stage K (1..N):** SIGTERM after stage K's sidecar is written but before stage K+1 starts; rerun-all; stages 1..K SKIP, stages K+1..N RUN; final outputs match the clean-run outputs.
  3. **Interrupt mid-stage (no sidecar):** SIGTERM during stage K's work (before sidecar write); rerun-all; stage K reruns from per-region sentinels (fast); stages K+1..N RUN.
  4. **Input mutation:** clean run; touch input BAM (re-hash); rerun-all; every stage reruns with `reason=input ... sha256 changed`.
  5. **Argv mutation (ignored flag):** clean run; rerun with different `--n-threads`; every stage SKIPs.
  6. **Argv mutation (load-bearing flag):** clean run; rerun with different `--min-coverage`; analyze reruns; correct SKIPs (analyze's input — corrected TSV — unchanged).
  7. **Force flag:** clean run; rerun with `--force-stage analyze`; analyze reruns; everything else SKIPs.
  8. **Git sha mismatch:** clean run; manually edit sidecar's `rectify_git_sha`; rerun without `--accept-prior-provenance`; all stages rerun. Rerun with `--accept-prior-provenance`; all SKIP.
  9. **Corrupt sidecar:** truncate a sidecar; rerun; that stage reruns with `reason=invalid sidecar`.
  10. **Missing output file:** delete one declared output; rerun; that stage reruns with `reason=missing output ...`.
  11. **Sample directory moved:** clean run; `mv $SAMPLE_OUTPUT $SAMPLE_OUTPUT_new`; rerun-all with `--output-dir $SAMPLE_OUTPUT_new`; every stage SKIPs (sample-relative paths resolve correctly under the new sample dir).
  12. **Env var drift ($SCRATCH changed):** clean run with `$SCRATCH=/tmp/scratch_A`; rerun with `$SCRATCH=/tmp/scratch_B` (after copying files); every stage SKIPs via `cached_absolute` fallback with the documented drift warning logged.
  13. **Env var unset:** clean run with `$OAK=/tmp/oak`; unset `$OAK`; copy the genome file to its `cached_absolute` location; rerun; stage SKIPs via cached_absolute fallback.
  14. **Env var unset AND cached path missing:** clean run; unset `$OAK`; delete the cached_absolute file; rerun; stage RUNs with `reason=cannot resolve input genome: ...`.
  15. **Transient sub_outputs missing:** clean run with region BAMs in `$L_SCRATCH`; verify `$L_SCRATCH` is empty after; rerun; every stage SKIPs (sub_outputs flagged transient, not sha256-checked).
  16. **Cross-cluster sidecar (informational):** copy a sidecar's `cluster: hoffman2` into a fake-Sherlock test env; rerun; warning logged; resolution falls back to `cached_absolute`; SKIPs if `cached_absolute` path exists on the fake-cluster filesystem, else RUNs cleanly.
- **New** `tests/integration/test_run_all_resume_cdna.py` — minimal version of the matrix for the cDNA arm.
- **New** `tests/integration/test_run_all_resume_netseq.py` — minimal version for NETSEQ.
- **Modify** `dev/specs/parallel_refactor_plan.md` — mark this commit DONE.

**Acceptance:**
- [ ] All 10 DRS resume scenarios pass.
- [ ] cDNA and NETSEQ resume scenarios pass.
- [ ] Resume test suite runs in < 5 min wall (acceptable for `slow` marker; not blocking the fast suite).
- [ ] `--dry-run-resume` output is human-readable: includes per-stage decision + reason + (when RUN) the specific diff that triggered it.
- [ ] Opus reviews integration test design + reads through resume traces.

**Estimated LOC:** ~400 new test code.
**Delegation:** Sonnet drafts the test scenarios from the matrix above; Opus reviews + designs the failure-injection harness (`SIGTERM`-then-resume in a controlled way).

---

### Commit G — Docs + CHANGELOG + deprecation roadmap

**Files:**
- **Modify** `CHANGELOG.md` — describe the parallelism additions, the manifest-as-canonical-artifact change, the `--legacy-single-threaded` flag, the back-compat plan for `corrected_reads.tsv`, the new resume model + sidecar schema + CLI flags.
- **Modify** `docs/PARALLELISM.md` — new doc (or update existing) describing the two axes, the merge points, when each is active, how to override for M1 vs cluster.
- **New** `docs/RESUME.md` — user-facing doc on the resume model: when stages skip vs run, the sidecar schema, the `--force-*` flags, `--dry-run-resume`, what to do when a sidecar is rejected. Include a worked example: "you killed the run after correct finished, here's what `run-all` resumes."
- **Modify** `docs/ALIGNER_RECOMMENDATIONS.md` — if Commit Zero surfaces aligner-specific timing tradeoffs, update.
- **Modify** `dev/specs/parallel_refactor_plan.md` (this file) — mark each commit's status; archive the doc to `dev/specs/archive/` once Commit G lands.

**Acceptance:**
- [ ] Docs build cleanly (mkdocs).
- [ ] CHANGELOG references the original issue (`ISSUE_parallel_bam_write_and_lazy_tsv.md`).
- [ ] Spec doc marked DONE.

**Estimated LOC:** ~200 doc.
**Delegation:** Sonnet drafts; Opus + user finalizes wording.

---

## 5. Validation strategy

### 5.1 Equivalence tests (Commits A, B, C, D, E, F each add their own)

The `tests/utils/bam_compare.py` helper QNAME-sorts both BAMs and compares record-by-record. The bar:

- **CIGAR, FLAG, POS, MAPQ identical.**
- **All tags identical** — including the rectify-specific `cp` tag set in `bam_writer.py`.
- **Read count identical** — no silent drops at region boundaries.

For TSVs: sort by `(chrom, pos, strand, read_id)`, then `cmp -s`. Any divergence is a failure unless explained (e.g., float formatting differences from a re-sum order — vanishingly unlikely with integer counts but possible with `mean_q_score`-style aggregates).

### 5.2 Region-boundary correctness

**Risk:** a read spanning a region boundary might resolve soft-clips or junction rescues differently in the partitioned vs. unified pass.

**Mitigation:** `get_processing_regions()` already splits at coverage gaps (low-coverage regions where no read spans the boundary). Verify this property in Commit A: assert that for every region boundary `pos`, `bam.count(chrom, pos-1, pos+1) == 0`. If violated, the boundary is unsafe and regions must be merged.

**Backup:** if reads do cross boundaries despite the gap split, add overlap zones (read fetched from both adjoining regions, deduped by `read_name` at merge time).

### 5.3 Smoke runs on real data

Required before each default-flip:

| Commit | Smoke sample | Expected wall time |
|---|---|---|
| B (DRS default) | Han wt_R1 (12.49M reads) on H2 16-core | < 75 min E2E |
| C (aligner queue) | Han wt_R1 on H2 16-core, 5 aligners | < 45 min for correct phase |
| D (analyze) | Han wt_R1 corrected output → analyze | < 15 min E2E analyze |
| E (cDNA stage-1) | by4742 cDNA PCB114 chunk | < 30 min |
| F (NETSEQ) | Han 2023 cluster_revisit fixture | < 5 min |

All smoke runs submitted via the `h2-qsub` skill.

### 5.4 Regression suite

After every commit:
```bash
ssh hoffman2 'bash -lc "
  cd /u/home/k/kevinroy/software/rectify
  source activate rectify
  pytest -m \"not slow\" 2>&1 | tail -30
"'
```
Slow markers run at the end of Commits B, E, F (anywhere cDNA or full-pipeline code paths change).

---

## 6. Backwards compatibility

### 6.1 `corrected_reads.tsv` consumers

The grep (run in Commit B prereqs) catalogs external consumers. Until that grep returns, the default is **`--emit-merged-tsv` ON** — `correct` produces both the manifest and a single concatenated TSV. The concatenated TSV is a cheap `cat` of region files; cost is ~5% of total runtime.

**Migration path post-grep:**
- If 0-5 external consumers: write a 1-page migration note, ping consumers, flip default to manifest-only in next minor version.
- If >5 or any external (collaborator scripts): keep emitted-by-default for one full release cycle (≥6 months), add a deprecation warning when read by the legacy path.

### 6.2 `--legacy-single-threaded` flag

Added in Commit B, removed in a future PR (NOT this one). Provides escape hatch for:
- Memory-constrained environments (e.g., debugging on M1 with other heavy work).
- Bisection of correctness regressions.
- CI environments without enough cores for the parallel path to win.

### 6.3 BAM file format

No change. Output BAM is bit-equivalent (modulo sort tie order) to legacy. Existing `corrected.bam` / `softclipped.bam` / etc. naming preserved.

### 6.4 PROVENANCE.json — two layers

Per `feedback_provenance_alongside_outputs.md`, every output gets a sidecar. Two layers coexist:

- **Sub-stage layer (per-region):** each region BAM/TSV gets a small sidecar listing region coords, worker_id, git_sha, invocation timestamp, sha256. Used by the parent stage's resume logic and by `--keep-tmp` debugging. Lives next to the temp file: `<tmp_dir>/region_NNN.bam.provenance.json`.
- **Stage layer (per pipeline stage):** each `run-all` stage (drs_trim, prescan, align, correct, restore_polya, merge_corrected_tsvs, write_corrected_bam, analyze, …) emits one top-level sidecar to the sample output directory. This is the artifact §6.5's skip-check consults. Lives at `<sample_output>/<sample_id>.<stage>.provenance.json`.

The stage-level sidecar references its sub-stage sidecars by path + sha256, so a chain-of-custody walk is possible. The temp dir is preserved if `--keep-tmp` is set; otherwise cleaned after merge — but the **stage-level sidecar persists in `<sample_output>` regardless**, since it's load-bearing for resume.

### 6.5 Stage-level resume via provenance sidecars

#### 6.5.1 Why this exists

`rectify run-all` chains 6-9 stages (DRS arm) or 4-6 stages (cDNA / NETSEQ arms). A run interrupted at stage K — by cluster preemption, SIGTERM at the 24-hr campus.q cap, Ctrl-C, machine reboot, or a fix-and-rerun cycle on a bug in stage K+1 — currently has no clean way to resume. The user either reruns from stage 1 (wastes hours of compute) or hand-edits skip flags into the invocation (fragile, undocumented). This section specifies the resume protocol.

#### 6.5.2 Sidecar schema (v1.0)

Stage-level sidecar at `<sample_output>/<sample_id>.<stage>.provenance.json`:

```json
{
  "schema_version": "1.0",
  "stage": "correct",
  "stage_subtype": "drs",
  "sample_id": "han_wt_R1",
  "rectify_git_sha": "b1565ce",
  "rectify_version": "0.9.0-dev",
  "started_at": "2026-05-19T14:32:11.234Z",
  "completed_at": "2026-05-19T14:54:27.891Z",
  "exit_status": 0,
  "host": "n1821.h2.local",
  "cluster": "hoffman2",
  "scheduler_job_id": "13446575",
  "invocation": {
    "argv": ["rectify", "correct", "input.bam", "--drs", "--Scer", "..."],
    "argv_template": ["rectify", "correct", "$SAMPLE_OUTPUT/wt_R1.aligned.bam", "--drs", "--Scer", "..."],
    "config_hash": "sha256:...",
    "cwd_template": "$SAMPLE_OUTPUT",
    "cwd_cached": "/u/project/guillom/shared/processed/han2023/wt_R1"
  },
  "inputs": [
    {"role": "bam", "sha256": "...", "size_bytes": 12500000000,
     "path": {"template": "wt_R1.aligned.bam",
              "cached_absolute": "/u/project/guillom/shared/processed/han2023/wt_R1/wt_R1.aligned.bam",
              "kind": "sample_relative"}},
    {"role": "genome", "sha256": "...",
     "path": {"template": "$OAK/references/scer/genome.fsa.gz",
              "cached_absolute": "/oak/stanford/groups/larsms/references/scer/genome.fsa.gz",
              "kind": "env_relative",
              "env_var": "OAK"}},
    {"role": "annotation", "sha256": "...",
     "path": {"template": "$OAK/references/scer/genes.gff.gz", "cached_absolute": "...", "kind": "env_relative", "env_var": "OAK"}},
    {"role": "prescan_variant_cache", "sha256": "...",
     "path": {"template": "prescan/rescue_scan.pkl", "cached_absolute": "...", "kind": "sample_relative"}},
    {"role": "prescan_junction_pool", "sha256": "...",
     "path": {"template": "prescan/junction_pool.pkl", "cached_absolute": "...", "kind": "sample_relative"}}
  ],
  "outputs": [
    {"role": "corrected_bam", "sha256": "...", "size_bytes": 12000000000,
     "path": {"template": "wt_R1.corrected.bam", "cached_absolute": "...", "kind": "sample_relative"}},
    {"role": "corrected_bam_index", "sha256": "...",
     "path": {"template": "wt_R1.corrected.bam.bai", "cached_absolute": "...", "kind": "sample_relative"}},
    {"role": "corrected_tsv_manifest", "sha256": "...",
     "path": {"template": "corrected_reads.manifest.tsv", "cached_absolute": "...", "kind": "sample_relative"}},
    {"role": "corrected_tsv_merged", "sha256": "...", "optional": true,
     "path": {"template": "corrected_reads.tsv", "cached_absolute": "...", "kind": "sample_relative"}}
  ],
  "sub_outputs": [
    {"role": "region_bam", "region_id": "r000", "sha256": "...",
     "transient": true,
     "path": {"template": "$L_SCRATCH/rectify_regions_${JOB_ID}/region_000.bam",
              "cached_absolute": "/lscratch/13446575/rectify_regions_13446575/region_000.bam",
              "kind": "env_relative",
              "env_var": "L_SCRATCH"},
     "sidecar": {"template": "$L_SCRATCH/rectify_regions_${JOB_ID}/region_000.bam.provenance.json",
                 "cached_absolute": "...", "kind": "env_relative", "env_var": "L_SCRATCH"}}
  ],
  "stats": {
    "n_reads_input": 12490000,
    "n_reads_corrected": 12490000,
    "n_reads_dropped": 0,
    "wall_seconds": 1334.6,
    "peak_rss_gb": 8.4
  },
  "skip_check_config": {
    "ignore_input_roles": [],
    "ignore_argv": ["--n-threads", "--tmp-dir", "--verbose", "--keep-tmp",
                    "--aligner-concurrency", "--legacy-single-threaded",
                    "--scratch-dir", "--temp-dir", "--output-dir"]
  }
}
```

**Field policy:**
- `inputs` and `outputs` use **role-keyed** entries, not positional. Stage code declares the role; the resume layer matches by role, not by path. This lets users move sample directories without breaking resume (sha256 is the truth; path is a locator hint).
- Every `path` field is a **PortablePath envelope** (see §6.5.2.5), not a bare string. Schema validator rejects bare strings.
- `skip_check_config.ignore_argv` lists flags whose value does NOT affect output bytes. **Path-related flags (`--output-dir`, `--tmp-dir`, `--scratch-dir`, `--temp-dir`) are always ignored** — paths are checked via the PortablePath resolution layer, not via argv string compare.
- `argv` (resolved) is recorded for diagnostics; `argv_template` (with env-var substitutions) is what feeds `config_hash` after `ignore_argv` filtering. Two runs that differ only in scratch path resolution produce the same `config_hash`.
- `cluster` field (`"hoffman2"` / `"sherlock"` / `"m1"` / `"other"`) is recorded so cross-cluster resume can warn loudly even if the user manually copies sidecars.

#### 6.5.2.5 PortablePath envelope and path resolution

Cluster environments make absolute-path storage unsafe for sidecars:

| Variable | Stability | Example |
|---|---|---|
| `$HOME` | per-user, stable | `/u/home/k/kevinroy` (H2), `/home/users/kevinroy` (Sherlock) |
| `$OAK` (Sherlock) | per-group, stable | `/oak/stanford/groups/larsms` |
| `$SCRATCH` (Sherlock) | per-user, stable | `/scratch/users/kevinroy` |
| `$L_SCRATCH` (Sherlock) | **per-job, ephemeral** | `/lscratch/13446575` (gone after job exits) |
| `$TMPDIR` (H2) | **per-job, ephemeral** | `/scratch/job/<jobid>` |
| `$SAMPLE_OUTPUT` (rectify-internal) | per-sample, mobile | wherever `--output-dir` points |

A naive sidecar storing `/lscratch/13446575/...` for a region BAM will never resolve again on the next run (different `$L_SCRATCH`). A naive sidecar storing `/u/project/guillom/...` for the sample dir will break if the user later mounts the data on Sherlock at `/oak/.../han2023/`. The PortablePath envelope solves both.

**Schema:**

```python
@dataclass(frozen=True)
class PortablePath:
    template: str            # tokenized form, e.g. "$OAK/refs/genome.fsa.gz" or "wt_R1.corrected.bam"
    cached_absolute: str     # the literal absolute path at sidecar-write time
    kind: Literal["sample_relative", "env_relative", "absolute"]
    env_var: Optional[str]   # set iff kind == "env_relative", e.g. "OAK"
```

**Tokenization at sidecar-write time** (`PortablePath.from_path(p, sample_output_dir)`):

1. Resolve `p` and `sample_output_dir` via `os.path.realpath` (handles symlinks).
2. If resolved `p` is inside resolved `sample_output_dir` → `kind="sample_relative"`, `template = relpath(p, sample_output_dir)`.
3. Else, walk `KNOWN_ENV_VARS = ("SAMPLE_OUTPUT", "L_SCRATCH", "TMPDIR", "SCRATCH", "OAK", "WORK", "HOME")` in that order (longest-prefix-first to avoid `$HOME` swallowing `$SCRATCH`). For the first env var whose realpath value is a prefix of `realpath(p)`, set `kind="env_relative"`, `env_var=<var>`, `template=f"${VAR}/{relpath(p, env_value)}"`.
4. Else → `kind="absolute"`, `template = str(p)`.
5. In every case, `cached_absolute = str(p)` records the literal resolved path at write time.

**Transient sub-outputs.** Region BAMs in `$L_SCRATCH` / `$TMPDIR` are job-scoped and disappear at job exit. Their sidecar entries have `"transient": true` and are explicitly **excluded from the skip-check's sha256 verification** (§6.5.3 step 1). They are recorded for chain-of-custody and `--keep-tmp` debugging only. A stage's persistent outputs (the merged BAM, manifest, etc.) always live in `$SAMPLE_OUTPUT` (never in scratch) and are non-transient.

**Resolution at sidecar-read time** (`PortablePath.resolve_now(sample_output_dir) -> Path`):

1. By `kind`:
   - `sample_relative` → `(sample_output_dir / template).resolve()`.
   - `env_relative` → expand `template` against current env via `os.path.expandvars`, then `realpath`. If `$VAR` is unset in the current environment, fall through to step 2.
   - `absolute` → `Path(template).resolve()`.
2. **Fallback to `cached_absolute`** if step 1 returned a non-existent path. Log a warning: `"path drift detected: template=$X/foo resolved to /current/y/foo (doesn't exist), falling back to cached_absolute=/old/x/foo"`. If `cached_absolute` exists, use it.
3. If neither exists → raise `PathResolutionError` with both candidates in the message; the skip-check converts this to `RUN(reason="cannot resolve input <role>: tried <template-resolved> and <cached_absolute>, neither exists")`.

**Cross-cluster guard.** When `prior.cluster != current_cluster`, the resume layer logs a loud warning `[RESUME] cluster mismatch: sidecar from hoffman2, current host is sherlock — env vars may resolve differently`. It does NOT auto-fail (cross-cluster resume CAN work if the same data has been replicated to compatible paths), but flags every PathResolutionError with `(cross-cluster: env vars likely don't match)`.

**Detection of current cluster** (`rectify/core/provenance/cluster.py`):
- Sherlock if `$SHERLOCK` set OR hostname matches `sh[0-9]+-[0-9]+.*`.
- Hoffman2 if hostname matches `n[0-9]+` OR `login[0-9]+` AND `$HOSTNAME` ends with `.hoffman2.idre.ucla.edu` (or similar marker).
- M1 if `os.uname().sysname == "Darwin"` AND `platform.machine() == "arm64"`.
- Else `"other"`.

#### 6.5.3 Skip-check decision tree

At each stage entry in `run-all`, before doing any work:

```python
def can_skip_stage(stage, sample_output, current_inputs, current_argv, rectify_git_sha,
                   *, force=False, accept_prior_provenance=False) -> SkipDecision:
    if force:
        return RUN(reason="--force-stage or --force-all")

    sidecar_path = sample_output / f"{sample_id}.{stage}.provenance.json"
    if not sidecar_path.exists():
        return RUN(reason="no prior sidecar")

    prior = load_sidecar(sidecar_path)              # raises on schema mismatch -> RUN
    if prior["exit_status"] != 0:
        return RUN(reason=f"prior run exited with status {prior['exit_status']}")

    # Cross-cluster warning (non-blocking; surfaces here, not at end)
    current_cluster = detect_current_cluster()
    if prior.get("cluster") and prior["cluster"] != current_cluster:
        logger.warning(f"[RESUME] cluster mismatch: sidecar from {prior['cluster']}, "
                       f"current = {current_cluster}; env vars may resolve differently")

    # 1. All declared NON-TRANSIENT outputs exist + sha256 match.
    #    Sub-outputs flagged transient (region BAMs in $L_SCRATCH / $TMPDIR) are SKIPPED here:
    #    they're recorded for chain-of-custody but don't gate skip decisions.
    for out in prior["outputs"]:
        portable = PortablePath.from_json(out["path"])
        try:
            p = portable.resolve_now(sample_output)
        except PathResolutionError as e:
            return RUN(reason=f"cannot resolve output {out['role']}: {e}")
        if out.get("optional") and not p.exists():
            continue
        if not p.exists():
            return RUN(reason=f"missing output {out['role']}: tried {portable.template} -> {p}")
        if sha256_of(p) != out["sha256"]:
            return RUN(reason=f"output {out['role']} sha256 changed since prior run")

    # (Sub-outputs intentionally NOT checked here. They are transient and may not exist.
    #  If a stage genuinely depends on a region BAM persisting, mark its sidecar entry
    #  transient=false and promote it to `outputs`.)

    # 2. Inputs match prior recorded inputs (sha256, not path).
    #    Path resolution tries template first, then cached_absolute, with cross-cluster warning.
    ignored_roles = set(prior["skip_check_config"].get("ignore_input_roles", []))
    for inp in prior["inputs"]:
        if inp["role"] in ignored_roles:
            continue
        current_path = current_inputs.get(inp["role"])
        if current_path is None:
            # Caller didn't supply this role this run. Fall back to prior's recorded path.
            try:
                current_path = PortablePath.from_json(inp["path"]).resolve_now(sample_output)
            except PathResolutionError as e:
                return RUN(reason=f"input {inp['role']} not provided and cannot resolve from sidecar: {e}")
        if not Path(current_path).exists():
            return RUN(reason=f"input {inp['role']} path resolved to {current_path} but file missing")
        if sha256_of(current_path) != inp["sha256"]:
            return RUN(reason=f"input {inp['role']} sha256 changed")

    # 3. Config args match (path-related flags excluded via ignore_argv).
    #    Compare against argv_template (env-var-tokenized), not raw argv.
    ignored_argv = set(prior["skip_check_config"].get("ignore_argv", []))
    current_template = tokenize_argv_paths(current_argv, sample_output)
    current_hash = normalized_config_hash(current_template, ignore=ignored_argv)
    if current_hash != prior["invocation"]["config_hash"]:
        return RUN(reason="config args changed (compare via --dry-run-resume for diff)")

    # 4. Rectify version.
    if prior["rectify_git_sha"] != rectify_git_sha and not accept_prior_provenance:
        return RUN(reason=f"rectify git_sha changed: {prior['rectify_git_sha'][:7]} -> {rectify_git_sha[:7]} "
                          f"(pass --accept-prior-provenance to override)")

    return SKIP(reason=f"prior valid run at {prior['completed_at']}, "
                       f"wall={prior['stats']['wall_seconds']:.1f}s")
```

A `SkipDecision.SKIP` returns the prior sidecar so downstream stages can use its `outputs` as their `inputs` without recomputing sha256.

**Key invariants:**
- Transient sub-outputs (region BAMs in `$L_SCRATCH` / `$TMPDIR`) are **never sha256-checked** at skip time. They are job-ephemeral by design.
- Persistent outputs that the next stage needs **must live in `$SAMPLE_OUTPUT`** (kind=`sample_relative`). Stages MUST NOT declare a `$L_SCRATCH` path as a non-transient output. The sidecar writer enforces this: attempting to write a non-transient `env_relative` output with `env_var ∈ {L_SCRATCH, TMPDIR}` raises `ValueError`.
- `argv` is tokenized to `argv_template` before hashing. Two invocations that differ only in `$L_SCRATCH` resolution (different job IDs, same logical command) produce the same `config_hash`.

#### 6.5.4 CLI flags

Added to `rectify run-all` (and the per-stage commands, for parity):

| Flag | Effect |
|---|---|
| (default) | Check sidecars; skip valid stages; rerun the rest. Print a one-line `[RESUME] stage=correct decision=SKIP reason=...` for each stage. |
| `--force-all` | Ignore all sidecars; rerun the entire pipeline. |
| `--force-stage <name>[,<name>...]` | Ignore sidecars for the listed stages only. Common pattern: `--force-stage analyze` after fixing an analyze bug. Stages downstream of a forced stage are also forced (their input sha256 will change). |
| `--accept-prior-provenance` | Treat rectify git_sha mismatch as non-blocking. Useful for cosmetic commits (docs / comments). Default: rerun on any git_sha change. |
| `--dry-run-resume` | Print the SKIP/RUN decision for each stage with full diff (which input or argv changed) and exit. No work done. |

#### 6.5.5 Interaction with per-region resume

The two resume layers are independent:

- **Stage-level resume (this section):** at stage K entry, decides whether to skip K entirely.
- **Region-level resume (§7.3, Commit A's `.ok` sentinels):** if stage K is running but was interrupted mid-flight, the per-region sentinels let workers skip already-done regions and resume from the unfinished ones.

When a stage runs (not skipped), the stage-level sidecar is **NOT** written until the stage's final merge completes successfully. A mid-flight crash leaves region sentinels but no stage-level sidecar → next `run-all` will rerun the stage, but the region sentinels make that rerun fast (workers skip done regions). After the stage's final merge writes the stage-level sidecar, region sentinels become redundant and are cleaned up unless `--keep-tmp` is set.

#### 6.5.6 Multi-sample manifests

For `rectify run-all --manifest samples.tsv` (multi-sample mode driven by `commands/run/multi_sample.py`), each sample's stages are independently sidecared in that sample's output directory. A `[RESUME] sample=wt_R1 stage=correct decision=SKIP` log line per (sample, stage) gives a clean picture. The manifest-driver itself does not have a sidecar — only the per-sample stages.

#### 6.5.7 Existing `tracker.register_staged` infrastructure

`commands/run/single_sample.py:507,529` calls `tracker.register_staged(...)` — a legacy tracking mechanism. **Commit A.5 prereq:** read `commands/run/helpers.py` (or wherever `tracker` is defined), decide whether to:
- Replace it entirely with the new sidecar layer (cleanest).
- Extend it to write sidecars as a side effect (lowest-risk).
- Keep it for its current purpose and layer sidecars alongside (most defensive).

Default to replace if `tracker` only handles path registration; extend if it handles other concerns (e.g., cleanup ordering). Decision recorded in Commit A.5 commit message.

---

## 7. Cross-cutting concerns

### 7.1 RAM budget

| Machine | Cores | RAM | Mode | Per-worker peak |
|---|---|---|---|---|
| M1 laptop | 8 | 8 GB | Axis A only (max 2 region workers, sequential aligners) | ~700 MB |
| H2 pod_smp.q@n1821 | 16 | 256 GB | Axes A + B; total workers = 14 | ~700 MB |
| Sherlock larsms (AMD Milan) | 64 | 256 GB | Axes A + B; total workers = 60 | ~700 MB |

Auto-detection: `rectify/core/utils/resources.py:detect_machine_class()` returns `"m1_laptop"` if `total_mem_gb < 16 AND os.cpu_count() < 8` — gates Axis B off.

### 7.2 Disk budget

Per-region BAMs total approximately 1× input BAM size. For a 12.5 GB input, scratch needs ~13 GB. Workers write to `$TMPDIR`:
- H2: `$SCRATCH` (typically `/u/scratch/<user>/`) or job-local `$TMPDIR`
- Sherlock: `$L_SCRATCH` (node-local, fastest, per-job ephemeral) or `$SCRATCH` (shared, larger, per-user stable)
- M1: `/tmp` (RAM-backed on macOS — adds memory pressure! Override to `$HOME/tmp/rectify_regions` on M1)

Pre-flight check at `bam_writer_parallel` entry: `shutil.disk_usage(tmp_dir).free` must be ≥ `1.5 × input_bam_size`. Hard error otherwise.

**Path portability and resume.** Per-region BAMs live in scratch dirs whose paths differ across runs (especially `$L_SCRATCH` / `$TMPDIR` which are per-job). The provenance layer (§6.5.2.5) handles this:
- Region BAMs are recorded as **transient sub_outputs** with `kind="env_relative"`, `env_var ∈ {L_SCRATCH, TMPDIR}`. Their sha256 is recorded but the skip-check (§6.5.3 step 1) explicitly does NOT verify them — they're job-ephemeral by design.
- The persistent merged BAM, manifest, and final TSV ALWAYS land in `$SAMPLE_OUTPUT` (which is on `$OAK` for Sherlock, `/u/project/` for H2, or wherever `--output-dir` points). The skip layer verifies these by sha256.
- The PortablePath envelope (§6.5.2.5) ensures sample-directory moves and env-var drift don't break resume — sample-relative paths re-resolve under the new sample dir; env-relative paths fall back to `cached_absolute` with a warning.

### 7.3 Sherlock `larsms,owners` preemption

Per-region `.ok` sentinel files enable resume: on rerun, workers skip regions whose `region_NNN.unsorted.bam.ok` exists. The `samtools merge` step at the end is NOT idempotent — if preempted mid-merge, must restart the merge from scratch (~30-60 s cost on 12.5 GB).

Document `--keep-tmp` as the recommended mode for `larsms,owners` runs.

### 7.4 Determinism

The dynamic work queue (Axis B) processes tasks in non-deterministic order. This affects:
- BAM sort tie order — broken by QNAME in the comparison helper, not in the output (output is coord-sorted).
- TSV row order within a region — deterministic (per-region workers process reads in fetch order).
- Manifest order — deterministic (regions written in coord order).

No user-facing nondeterminism.

### 7.5 Logging

Each worker logs to a per-region file under `$TMPDIR/rectify_regions/region_NNN.log`. Main thread tails errors at completion. Existing top-level `[TIMING]` log lines preserved.

---

## 8. Open questions / decision gates

Resolved before each commit:

| Question | Resolved by | Default if unresolved |
|---|---|---|
| Does BGZF compression dominate `write_corrected_bam`? | Commit Zero cProfile | Assume no; proceed with region parallelism |
| Is `_run_alignment` sequential across aligners today? | Commit Zero timing | Assume yes; include in Commit E |
| Is `merge_corrected_tsvs` a major bottleneck? | Commit Zero timing | Defer to a future PR if <10% wall |
| Does `VariantAwareHomopolymerRescue` have mergeable state? | Read `core/correct/indel_corrector.py` in Commit D | If no, region overlap zones; if yes, simple frequency-sum merge |
| How many external consumers depend on `corrected_reads.tsv`? | Back-compat grep in Commit B prereqs | Keep emitted-by-default if >5 |
| Does region-boundary read crossover occur in practice? | Assertion in Commit A | If yes, add overlap zones |

---

## 9. Risks (in order of impact)

1. **Profile mis-attribution.** If we proceed without Commit Zero and `write_corrected_bam` turns out to be BGZF-bound, the entire architecture is overbuilt. Mitigation: Commit Zero is hard-gated.
2. **Region-boundary correctness divergence.** Reads spanning boundaries could resolve mutations differently in partitioned vs unified passes. Mitigation: §5.2.
3. **Sherlock `larsms` preemption mid-merge.** Acceptable cost (~30-60 s), but document it.
4. **M1 fallback misses a case.** If `detect_machine_class()` doesn't fire (e.g., user runs on a 4-core 16-GB cloud VM), they get the full Axis B treatment and OOM. Mitigation: also gate on `os.cpu_count()` AND `total_mem_gb`, log the detected class loudly at startup, document the override flag.
5. **Worker startup cost regresses small tests.** Mitigated by the small-input gate in Commit B (§3 row 5).
6. **External consumer breakage** from manifest-only default. Mitigated by §6.1 back-compat policy.
7. **Path-portability bugs in sidecars** — e.g., a stage records `$L_SCRATCH/foo` as a non-transient output, or `cached_absolute` falls back silently when env var resolution actually changed semantics (different data at the same env-relative path). Mitigated by: (a) PortablePath construction refuses non-transient `env_relative` paths in `$L_SCRATCH`/`$TMPDIR` (Commit A.5 invariant), (b) sha256 verification is the source of truth — if cached_absolute resolves to different bytes, skip-check returns RUN, (c) cross-cluster mismatch logs a loud warning. The residual risk: a user copies a sidecar across machines AND copies the underlying file but with different bytes — sha256 catches this.
8. **DRS users will see far smaller speedup than QSrev users.** Commit Zero profiling (2026-05-19) showed DRS bottleneck (`realign_exon_blocks` inside the already-parallel correct-region loop, 91% of wall) is not addressed by this refactor's region-parallel write infrastructure. Expected DRS wall-time improvement: 1.1-1.3× (vs ~6× for QSrev). Mitigation: communicate clearly in CHANGELOG + docs/RESUME.md that the parallelism refactor delivers (a) QSrev/cDNA-short ~6× wall-time speedup, (b) DRS/cDNA-long: architectural improvements (resume, sidecars, manifest TSVs, analyze partial-streaming) but minimal raw wall-time win, (c) future algorithmic work (cross-phase realign caching) needed for DRS to see comparable speedup. Don't let users feel misled by the "default-parallel RECTIFY" framing if they're running DRS workloads.

---

## 10. Estimated effort

| Phase | Wall (incl cluster waits) | Agent time |
|---|---|---|
| Commit Zero (profile) | 4 hr | 30 min Opus |
| Commit A | 1 day | ~3 hr Sonnet + 30 min Opus review |
| Commit A.5 (resume infra) | 1 day | ~3 hr Sonnet + 1.5 hr Opus (tracker decision + skip-check correctness) |
| Commit B | 1.5 days | ~4 hr Sonnet + 1.5 hr Opus (incl. first sidecar wiring) + smoke wall |
| Commit C | 1 day | ~2 hr Opus + smoke wall |
| Commit D | 1.5 days | ~3 hr Opus + 2.5 hr Sonnet + smoke wall |
| Commit E | 1 day | ~3.5 hr Sonnet + 30 min Opus + smoke wall |
| Commit F | 0.5 day | ~2 hr Sonnet + 30 min Opus |
| Commit F.5 (resume integration tests) | 0.75 day | ~3 hr Sonnet + 1 hr Opus |
| Commit G | 0.25 day | ~1.5 hr Sonnet + 30 min Opus |
| **Total** | **~9 days wall** | **~30 hr agent time** |

---

## 11. Delegation policy

### 11.1 Opus owns

- **Commit Zero entirely** (profile interpretation drives all downstream decisions).
- **Commit A.5 high-judgment pieces:** `tracker.register_staged` migration decision (replace / extend / coexist) and the skip-check decision tree implementation (correctness-critical, one wrong branch silently corrupts resume).
- **Commit B first stage-level sidecar wiring** — sets the template every later commit follows.
- **Commit C** (concurrency primitives, RAM autosizing, M1 vs cluster boundary).
- **Cluster-calling two-pass refactor in Commit D** (correctness traps).
- **Commit F.5 failure-injection harness design** (SIGTERM-then-resume in a controlled way is tricky to get right; Sonnet implements the scenarios but Opus designs the framework).
- **All correctness-divergence diagnostics** when equivalence or resume tests fail.
- **Back-compat grep + default decision** in Commit B prereqs.
- **Per-commit diff review** before merge to `drs-validation-rebuild`.

### 11.2 Sonnet subagents own (one commit at a time)

- Commit A scaffolding (well-specified, low judgment).
- Commit A.5 boilerplate: `sidecar.py`, `hashing.py`, `cli.py`, unit tests for each.
- Commit B mechanical wiring (after Opus runs the grep, picks the default, and wires the first stage sidecar).
- Commit D per-substep streaming helpers + test fixtures.
- Commit E sub-pieces (each is a self-contained transform of an existing function).
- Commit F (pure copy-of-template work).
- Commit F.5 test scenarios (drafted from the §F.5 matrix; Opus reviews + the harness is Opus-designed).
- Commit G doc drafts (including `docs/RESUME.md`).

### 11.3 Briefing policy

Per-commit Sonnet briefings are written **at the time that commit starts**, not all up front. Reason: Commit Zero's profile reshapes the rest of the plan. Pre-writing Commit C while Commit Zero is in flight risks invalidation.

Each briefing includes:
- Goal (one paragraph from this spec).
- File-and-line references for every touched location.
- Acceptance criteria (copy from §4).
- Verify-don't-rewrite instructions per CLAUDE.md: "prove every claim with a tool call (grep / import / --help) rather than recite from memory."
- The relevant section of this spec attached verbatim.
- A no-write-without-grep rule: if the briefing says "modify `foo.py:42`", the subagent must `Read foo.py` first and confirm the line content matches before editing.

### 11.4 Subagent guardrails

- **No `git add -A` / `git add .`** — surgical staging per CLAUDE.md "Surgical staging" rule.
- **No commits to `master`** — only `drs-validation-rebuild`.
- **No skipping pre-commit hooks** unless explicitly authorized.
- **Cross-cluster commit-status reconciliation** before any subagent that runs on H2 or Sherlock — verify M1 / GitHub / cluster heads match per CLAUDE.md "Cross-cluster commit-status check" rule.
- **No phantom-flag punch lists** — if claiming a file has a bug or a function is wrong, the subagent must grep-verify the claim and quote the offending lines.

---

## 12. Resumption checklist (if this work is interrupted mid-stream)

After context loss / session switch, the orchestrator should:

1. `git log --oneline -10` on `drs-validation-rebuild` to see which commits landed.
2. Read this spec (`dev/specs/parallel_refactor_plan.md`) — find each commit's checkboxes; the unchecked ones are the work remaining.
3. Read `dev/specs/profile_results.md` — confirm Commit Zero decision still holds (if Commit Zero is done).
4. Cross-cluster commit-status check (M1 = GitHub = H2 = Sherlock heads).
5. Re-read this section (§4) for the next unfinished commit; resume from its acceptance-criteria list.

If `profile_results.md` does not exist, **Commit Zero has not run** — do not start any other commit. This is the hard gate.

**For mid-flight `rectify run-all` invocations** (a separate concern from this spec's own implementation resume): once Commit A.5 + B have landed, the user can resume any interrupted run by reinvoking the same `rectify run-all` command. The resume layer reads sidecars and skips valid stages. Use `--dry-run-resume` to preview the decision tree, `--force-stage <name>` to force-rerun a specific stage, or `--force-all` to nuke and restart. The pipeline's resumability is a user-facing feature (§6.5); this spec's resumability (the checklist above) is an internal orchestration tool for the implementation work.

---

**End of spec.**
