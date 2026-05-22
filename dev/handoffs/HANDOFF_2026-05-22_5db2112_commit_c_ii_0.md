# Handoff — Commit C-II.0 landed, C-II.1 ready to start

**Date:** 2026-05-22
**Branch:** `drs-validation-rebuild`
**HEAD:** `5db2112` — `refactor(bam): split worker state into shared init + per-task tuple`
**Origin:** `5db2112` (pushed)
**Clusters:** H2 has `5db2112` files staged via rsync (no `.git`); Sherlock untouched (`1ab71f0` rsync from earlier sessions, sync before next cluster run).

---

## 1. What landed this session

Three commits, in order:

- `d762b91` — **Preparatory Commit C**: accept manifest output + plumb `--aligner-concurrency` in per-aligner driver. Fixes a latent Commit B regression in `_run_correction_per_aligner` (skip-check + success-check now accept both `corrected_reads.manifest.tsv` and `corrected_reads.tsv`). Resolver from `resources.py` is now invoked but the loop stays sequential — Axis B perf explicitly deferred. 8 new tests in `tests/test_run_correction_per_aligner.py`.

- `5db2112` — **Commit C-II.0**: split `bam_parallel.py` worker state into shared init bundle + per-task tuple. The shared bundle (`genome`, `polya_model`, `annotated_junctions`, `pool_chrom_index`, `gene_interval_trees`, config flags, `netseq_dir`) installs once per pool initializer; per-aligner sliver (`bam_path`, `variant_aware_rescue`) travels in each region task tuple. Both `process_bam_file_parallel` and `process_bam_streaming_parallel` updated. Pure semantic no-op proven by byte-equivalence test (`tests/test_bam_parallel_state.py`) with recorded golden hash `b267b1da6b8533a84b3701c7b1974f7666fba80fd602667d62bbfb80f1017c95`.

Co-committed with sibling agents this session:

- `6d1c83e` — docs: README_KR_edits merge + CLAUDE.md path updates + doc consolidation proposal.
- `7ecf534` — fix(analyze): missing `cluster_gene_attribution.py` module (fresh-clone-fixing).
- `8c9f80c` — feat(analyze): sidecar wiring for analyze + restore_polya stages (Commit D).
- `55089f7` — perf(bam): reanchor O(CIGAR-ops) shortcut.

---

## 2. Verified

- Fast suite slice for C-II.0 blast radius: 148 passed + 8 skipped + 0 failures (1:17 wall).
- Tests: `test_bam_parallel_state.py`, `test_correct_command_drs.py`, `test_correct_command_parallel_default.py`, `test_resume_correctness.py`, `test_run_correction_per_aligner.py`, `test_validation_reads.py`, `test_resources.py`.
- Byte-equivalence: validation BAM at `n_threads=2` produces the same JSON-canonical-sorted hash before and after the worker-state split. The `n_threads=1` (no Pool) and `n_threads=2` (Pool) paths converge — no serial/parallel divergence.
- Origin = local = `5db2112`; H2 carries the two C-II.0 files at matching md5.

---

## 3. Open items

### Commit C-II.1 — shared pool across aligners + work-stealing

**Owner:** Opus, advisor-gated.

**Updated design briefing (v2):** `Chanfreau Lab/commit_c_ii_briefing_DRAFT.md` in Drive cwd. The v1 (`dev/specs/briefings/commit_c_ii_briefing_DRAFT.md` if a sibling agent moved it — check both locations) was based on a wrong premise; v2 reflects the actual `bam_parallel.py` shape and proposes a 2-commit C-II.0 + C-II.1 split.

**Discovered constraint not in the briefing:** `process_bam_file_parallel` does the prescan + region planning + variant_aware first-pass INSIDE the function, before the pool spawns. For a true cross-aligner shared pool, those pre-pool steps need to move out so they can run per-aligner serially BEFORE the shared pool starts, building up the full task list. This is a real lift on top of the C-II.0 worker state split.

**Files C-II.1 will touch:**
- `rectify/core/commands/run/stages.py` — `_run_correction_per_aligner` becomes a coordinator.
- `rectify/core/bam/parallel.py` — expose pre-pool prep separately from pool execution.
- `rectify/core/commands/correct_command.py` — accept injected pool.
- Possibly new `rectify/core/correct/region_dispatcher.py`.

**Acceptance for C-II.1:**
- Existing `--aligner-concurrency 1` path remains byte-identical (verify via existing tests).
- `--aligner-concurrency auto` actually parallelizes — measurable via the validation pipeline wall time.
- **Cluster gate:** Han wt_R1 on H2 (16-core node), correct phase wall ≤ slowest_aligner_alone × 1.3. Submit via h2-qsub; poll every 5 min per CLAUDE.md.

### Heap-corruption smoke (AGENT_FIXES "STILL OPEN" entry, 2026-05-20)

Not run this session. Deferred per prior AGENT_FIXES note ("Han wt_R1 6.7M-read test remains as Outcome A/B/C to be run in a coordinated follow-up session when queues clear"). C-II.0 changes the worker state model — could help or hurt the heap-corruption signature. **Don't claim Axis B safe-at-scale until this smoke clears.**

Setup notes for the smoke:

- BAM candidate path on H2: `/u/project/guillom/shared/processed/han_2023_ysh1aa_vs_wt_rectify_v0.9.0/results/wt_R1/wt_R1.rectified.bam` (325 MB; need to confirm this is the 6.7M-read post-consensus BAM, not a smaller earlier-pipeline version).
- Invocation pattern from prior run (note: that run was QSrev, not DRS — DRS smoke needs the DRS Han wt_R1 BAM, different file): in the README at the same path.
- Queue: H2 pod_smp.q. As of session end, no current queue check captured.
- Resource sizing: 16 threads, 16-core node, h_data=8G per slot, h_rt=12h (the chrI-V subset took 61:40 wall, so full BAM may take 2-3 hours plus queue).

### Documentation chase

- AGENT_FIXES.md needs an entry documenting C-II.0 (the worker state split) so future sessions see the design choice without rereading commit messages.
- HANDOFF.md (the main file) has Kevin's WIP from earlier sessions; not touched this session to avoid the concurrent-session race documented in `feedback_concurrent_sessions_git_race.md` (memory). Resume agents can fold this dated handoff into HANDOFF.md proper once the race-active window passes.

---

## 4. Resume command

```bash
# M1 preflight
cd /Users/kevinroy/work/rectify
git log --oneline -5
git status -s | head -15

# Confirm cross-cluster heads
git ls-remote origin drs-validation-rebuild | awk '{print substr($1,1,7)}'
ssh hoffman2 'cd /u/home/k/kevinroy/software/rectify && md5sum rectify/core/bam/parallel.py tests/test_bam_parallel_state.py'
md5 /Users/kevinroy/work/rectify/rectify/core/bam/parallel.py /Users/kevinroy/work/rectify/tests/test_bam_parallel_state.py

# Read the updated C-II briefing before starting C-II.1
cat "/Users/kevinroy/Library/CloudStorage/GoogleDrive-kevinrjroy@gmail.com/My Drive/Work/Chanfreau Lab/commit_c_ii_briefing_DRAFT.md"

# Sanity check pre-existing tests still pass
python -m pytest tests/test_bam_parallel_state.py tests/test_run_correction_per_aligner.py tests/test_resources.py -q
```

When starting C-II.1:

1. Read the briefing's §1 "Concrete shape" and §2 "Files".
2. Call advisor before substantive code.
3. Implement in stages: extract pre-pool prep from `process_bam_file_parallel` → modify `_run_correction_per_aligner` to use shared pool → demux per-aligner results → finalize per-aligner.
4. Verify `--aligner-concurrency 1` path remains byte-identical via the existing `test_bam_parallel_state.py` golden.
5. **Atomic stage+commit** (one Bash invocation chained with `&&`) — concurrent sessions may still be active.

---

## 5. Hard rules for the resolver

Same as prior handoff:

- NEVER `git add -A` or `git add .` — working tree has WIP for multiple upcoming commits.
- NEVER commit on a cluster — M1 is the git source of truth.
- Run `pytest -m "not slow"` slice (the 148-test blast radius) before each commit touching `bam_parallel.py`.
- Branch is `drs-validation-rebuild`. `master` is frozen at 0.9.0.
- Multiple claude sessions still active on M1 — stage and commit must be one atomic Bash call per `feedback_concurrent_sessions_git_race.md`.

---

## 6. Memory updates this session

- `feedback_concurrent_sessions_git_race.md` (new) — stage+commit atomicity rule when sibling sessions are active.
- `project_rectify_commit_c_status.md` (new) — Preparatory C is landed; Axis B perf split into C-II.0 (done) + C-II.1 (pending).
- Both linked from `MEMORY.md`.
