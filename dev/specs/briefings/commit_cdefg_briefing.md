# Follow-up Briefing — Commits C, D, E, F, F.5, G

**Audience:** the next agent picking up the parallelism refactor.
**Status at handoff (2026-05-20 evening):** Commits Zero, A, A.5, B landed; multiple parallel-session commits have shifted the baseline since. Current HEAD: `d790c00 docs(handoffs): archive active session notes` (probe to confirm).

**Important:** before consulting this umbrella briefing for Commit C, read `dev/specs/briefings/commit_c_briefing.md` (which already exists, landed in `4200d34`) — it is the canonical detailed plan for Commit C. This umbrella briefing summarizes C briefly and focuses on D/E/F/F.5/G which don't have per-commit briefings yet.

**This is ONE briefing covering FIVE commits.** Don't try to land them all in one session — that's 3-5 days of work. **Pick one commit at a time** based on the sequencing in §0. Each commit's section below is self-contained and can be extracted as a standalone Sonnet brief when that commit is actually dispatched.

**Master spec to read first:** `dev/specs/parallel_refactor_plan.md`, especially §4 Commit C / D / E / F / F.5 / G subsections and §11 (Opus vs Sonnet ownership split). This briefing FILLS IN current-state deltas; it does NOT replace the master spec.

---

## 0. Critical orientation — state since `ba61d31`

The session that landed Commit B ended at `ba61d31`. Multiple parallel sessions have committed since then. **PROBE GIT STATE FIRST** before doing anything else:

```bash
cd /Users/kevinroy/work/rectify
git log --oneline -20
git status -s | head -30
git ls-remote origin drs-validation-rebuild | awk '{print substr($1,1,7)}'
ssh hoffman2 'bash -lc "cd /u/home/k/kevinroy/software/rectify && git log --oneline -1"'
ssh sherlock 'echo ok && bash --norc --noprofile -c "cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify && git log --oneline -1"'
```

### Known-landed work since `ba61d31` that affects C–G (read in commit order, oldest first)

- **`cfdca4a fix(analyze,scripts): bedgraph 1-bp left-shift in 3 emitters`** — corrected_3prime is 0-based; 3 bedgraph emitters had an off-by-one. Touches `analyze_command.py`-adjacent code; verify your Commit D analyze refactor preserves the fix.
- **`45cbc13 fix(qname): cross-aligner QNAME hardening + auto-sanitize validator`** — touches `consensus.py` (+51), `corrected_consensus.py` (+45), per-aligner wrapper QNAME sanitization, plus RN aux tag plumbing (read_num sidecar architecture — see `dev/specs/briefings/read_num_sidecar_briefing.md` and HANDOFF.md). **Affects D scope:** `merge_corrected_tsvs` consumes per-aligner BAMs whose QNAMEs are reliably normalized AND may carry `RN:i:` aux tags; D's parallelization must respect both keying paths.
- **`6c7f846 fix(restore-polya): accept corrected TSV manifests`** — `restore_polya_command.py` is now manifest-aware (+99 LOC, 53-test fixture). **CLOSES one of two TODOs from Commit B.** Still pending: the `restore_polya` STAGE SIDECAR emission (separate from the manifest-loader piece).
- **`4200d34 docs(parallel): add Commit C briefing`** — `dev/specs/briefings/commit_c_briefing.md`. Read it before doing C; it's the canonical detailed plan.
- **`5235e59 feat(parallel): add aligner concurrency policy`** — `rectify/core/utils/resources.py` (~63 LOC, `resolve_aligner_concurrency()`) + integration in `run_command.py` + 52-test `test_resources.py`. **This is the AUTO-SIZING POLICY piece of Commit C, NOT the dynamic work queue.** `_run_correction_per_aligner` at `stages.py:294` is still sequential across aligners. Commit C's remainder is the queue refactor.
- **`a4e2933 fix(cli): validate aligner concurrency argument`** — small followup validator.
- **`12816fe fix(consensus): accept corrected TSV manifests`** — `corrected_consensus.py:545` adds `_read_corrected_tsv_or_manifest()`; `merge_corrected_tsvs` is now manifest-aware. **CLOSES the second TODO from Commit B.** D's `merge_corrected_tsvs` parallelization still pending.
- **`d73c1f2 docs(dev): triage dirty worktree disposition`** — WIP cleanup pass; check working-tree state with a fresh `git status -s` before staging.
- **`717e99e docs(agent_fixes): set2 cascade scale entry`** — AGENT_FIXES.md updated: the heap-corruption bug ALSO cascades at CHUNK scale via merged-BAM lookups in `--aligner-bams` paths, not just at 6.7M-read full-sample scale. **Affects E (multi-aligner align) and D (merge_corrected_tsvs)**: any code that opens the merged consensus BAM via pysam fetch in workers may trigger the same htslib state corruption.
- **`e0c9dc3 docs(briefing): axis 2 htslib state corruption`** — `dev/specs/briefings/axis2_htslib_state_corruption_briefing.md`. **READ THIS** before touching anything in `--aligner-bams`-consuming workers. The bug is in Axis 2 (htslib internals); the parallelism refactor (Axis 1) is hypothesized to side-step it via per-region pre-partitioning but is unproven at scale.
- Subsequent docs commits (`cb27a53`, `8ca5795`, `a358011`, `d022398`, `2187b8f`, `d790c00`) — docs/handoff updates, no code impact on C-G.

### Working tree expectations

Lots of WIP from parallel sessions. **Before staging any commit:** run `git status -s`, identify which files are YOURS to touch vs which are parallel-session WIP, and stage surgically (per `CLAUDE.md` "Surgical staging" rule). Never `git add -A` or `git add .`.

### Cross-cluster sync hygiene

Per `CLAUDE.md` "Cross-cluster commit-status check at every session start" rule. Probe Sherlock ControlMaster with `ssh sherlock 'echo ok'` before assuming it's live (memory `feedback-sherlock-controlmaster-test`).

### Recommended commit ordering

C → D → E → F → F.5 → G.

(The restore_polya + merge_corrected_tsvs manifest-loader cleanup from Commit B's TODO list is **already done** via `6c7f846` + `12816fe`. The restore_polya STAGE SIDECAR emission is still pending and folds naturally into D.)

Reasoning:
- C is mechanically self-contained and independent of D/E/F/G — picks off the work-queue refactor that the existing `resources.py` helper enables. See `commit_c_briefing.md` for the detailed plan.
- D restructures analyze to consume per-region partials AND parallelizes `merge_corrected_tsvs` AND emits the restore_polya stage sidecar. After D, E/F can wire their stages' sidecars + manifest paths against a stable analyze contract.
- E and F are mechanical Sonnet-style work; their order doesn't matter to each other.
- F.5 is the resume integration test matrix — needs all prior commits done so every stage emits a sidecar.
- G is final docs + CHANGELOG; happens last.

### Hard rules across all commits

1. **NO `git add -A` / `git add .`.** Surgical staging only.
2. **Verify, don't recite.** Per CLAUDE.md. Grep + Read before claiming a function's location/signature.
3. **`pytest -m "not slow"` MUST pass** after every commit. Current baseline is ~1133 + whatever's landed since.
4. **Watchdog mitigation:** prior agents stalled at 600s pytest transcript silence. Use `pytest -v --tb=line`; print ~60s status pings during long commands.
5. **No commits / no pushes by Sonnet subagents.** Opus reviews + commits + pushes.
6. **No changes to `rectify.utils.provenance` (legacy run log) or `rectify.core.provenance` (the A.5 package).** They're done.
7. **The full-scale Han wt_R1 6.7M-read structural smoke remains DEFERRED** until the axis-2 investigation produces a verdict OR cluster queues clear AND the user gives an explicit go. Do not opportunistically dispatch it from this commit chain.

---

## 1. Commit C — Dynamic aligner-queue concurrency (REMAINDER)

**Canonical detailed briefing:** `dev/specs/briefings/commit_c_briefing.md`. **READ THAT FIRST.** This section is the executive summary.

**Status at handoff:** AUTO-SIZING POLICY landed in `5235e59` (`resources.py:resolve_aligner_concurrency`); CLI flag landed in `a4e2933`. The DYNAMIC WORK QUEUE in `_run_correction_per_aligner` at `stages.py:294` is still pending.

**Owner per spec §11.1:** Opus (concurrency primitives, RAM autosizing, M1 vs cluster boundary). The remaining work-queue refactor is correctness-critical; if dispatched to Sonnet, Opus must review the queue lifecycle + error propagation logic carefully.

### 1.1 What's left

Replace the sequential per-aligner correction loop at `rectify/core/commands/run/stages.py:294-415` (`_run_correction_per_aligner`) with a single shared `multiprocessing.Queue` of `(aligner_id, region_id)` task tuples and N workers pulling dynamically. Reason: minimap2 finishes 3-5× faster than deSALT/gapmm2 on the same regions; static aligner partitioning idles the fast workers.

### 1.2 Design

- Total worker count: from `resolve_aligner_concurrency()` (already lands the M1 vs cluster decision). Workers can range from `n_aligners` (1 worker per aligner, mostly idle) down to a dynamic pool sized by `total_cores - reserve`.
- Single `multiprocessing.Queue` holds tasks `(aligner_id, region_id, region_plan)`. Each worker `get`s a task, runs the correction for that (aligner, region) pair, emits the per-region BAM/TSV under the aligner's subdir, then loops for the next task.
- When all tasks for a given aligner complete, kick off that aligner's final `samtools merge` + index + sidecar emission in a small thread (one thread per aligner; non-competing with the work queue).
- Errors: a worker exception poisons the queue; remaining workers drain gracefully; the main thread re-raises with the offending `(aligner_id, region_id)` in the message.

### 1.3 Files to modify

- `rectify/core/commands/run/stages.py` — replace `_run_correction_per_aligner` body. Function signature unchanged so callers in `single_sample.py` and the chunked-batch driver don't need updates.
- `rectify/core/commands/run/stages.py` — possibly extract the queue management into a helper `_run_aligner_correction_queue()` to keep the main function readable.
- `tests/test_aligner_concurrency.py` — NEW (or extend `test_resources.py`). At least 4 tests:
  1. Single aligner → falls back to direct call (no queue overhead).
  2. Multi-aligner sequential mode (`--aligner-concurrency 1`) → same result as pre-Commit-C.
  3. Multi-aligner dynamic mode → all aligners' outputs equivalent to sequential mode.
  4. Worker exception → main thread receives the error with `(aligner_id, region_id)` context.

### 1.4 Acceptance

- [ ] On 16-core H2 node: Han wt_R1 5-aligner correct phase runs in wall time ≤ `(slowest_aligner_alone × 1.3)`. Sentinel: log line `[aligner_queue] finished aligner_id=<X> in <wall>s, total queue wall <total>s`.
- [ ] On M1 (auto-detected): falls back to sequential aligners per `resolve_aligner_concurrency()`. Verify by checking active process count peaks at `N+1`, not `5N+1`.
- [ ] `--aligner-concurrency 1` produces output identical to pre-Commit-C (sequential aligners). Use `tests/utils/bam_compare.py` for cross-aligner equivalence.
- [ ] `pytest -m "not slow"` passes.
- [ ] No OOM on 64 GB H2 node with 5 aligners concurrent at default sizing.

### 1.5 Estimated LOC

~150-200 LOC modified in `stages.py`; ~200 LOC tests. The auto-sizing piece is already done.

---

## 2. Restore_polya + merge_corrected_tsvs manifest cleanup — MOSTLY DONE

**Status:** the manifest-aware loader pieces both landed during parallel sessions:

- **`6c7f846 fix(restore-polya): accept corrected TSV manifests`** — `restore_polya_command.py` now manifest-aware. Tests at `tests/test_restore_polya_manifest.py`.
- **`12816fe fix(consensus): accept corrected TSV manifests`** — `corrected_consensus.py:545` `_read_corrected_tsv_or_manifest` wraps `merge_corrected_tsvs` callers.

**Still pending (folds into Commit D):**
- `<sample_id>.restore_polya.provenance.json` sidecar emission. Add at restore_polya stage entry: `can_skip_stage` + `[RESUME]` log line. At success exit: `write_stage_sidecar(...)` with inputs (corrected TSV/manifest + original poly-A metadata parquet from `rectify trim-polya`) and outputs (restored BAM + index). Roughly 50 LOC in `restore_polya_command.py`.

No standalone commit needed; D absorbs this.

---

## 3. Commit D — Analyze partial-streaming + merge_corrected_tsvs parallelism

**Owner per spec §11.1:** Opus owns the cluster-calling two-pass refactor + equivalence interpretation. Sonnet can own the per-substep streaming helpers + analyze sidecar wiring + test fixtures.

**Conditional on Commit Zero results (recorded in `profile_results.md`):** if `merge_corrected_tsvs` shows up as <10% of wall time at scale, defer that piece to a later PR and ship D as analyze-only.

### 3.1 Scope

Restructure `analyze_command.py:95` (`load_corrected_positions`) and downstream substeps to consume per-region partials from a manifest, only collapsing to a global view where mathematically required.

Per-substep migration (spec §1.2 table):

| Substep | Current | After D |
|---|---|---|
| `generate_bedgraphs` | global df | per-chrom streaming + final concat |
| Stats (`corrected_reads_stats.tsv`) | aggregate over df | per-region reduce-sum |
| Position index (`*_index.bed.gz`) | sorted global → tabix | per-chrom tabix shards + final concat |
| 3′ histograms | global | per-shard summable bins |
| Cluster calling (`cluster_cpa_sites_adaptive`) | global df | TWO-PASS: pass 1 streamed coverage quantiles in constant RAM; pass 2 per-chrom calls via `ProcessPoolExecutor` with pre-computed thresholds |
| DESeq2 count matrix | small input | unchanged |

### 3.2 Files

- `rectify/core/analyze/streaming.py` — NEW, `load_corrected_positions_partitioned(manifest_path) → Iterator[ChromShard]` + per-substep helpers.
- `rectify/core/commands/analyze_command.py` — detect manifest vs single-TSV; dispatch substeps to streaming forms.
- `rectify/core/bam/loaders.py` — already manifest-aware from Commit B; verify the partial-streaming forms match.
- `rectify/core/consensus/corrected_consensus.py:merge_corrected_tsvs` — split per-chrom, run consensus calling in `ProcessPoolExecutor`, final concat. ONLY if Commit Zero flags it.
- `rectify/core/commands/cdna_analyze_command.py` — share the same streaming loader.
- `<sample_id>.analyze.provenance.json` sidecar emission.
- `tests/test_analyze_streaming_equivalence.py` — NEW, ~5 tests comparing analyze outputs byte-for-byte across streaming vs legacy paths on a small fixture.

### 3.3 Cluster-calling two-pass (Opus owns; correctness-critical)

The adaptive cluster caller `cluster_cpa_sites_adaptive` uses global coverage quantiles. Two-pass refactor:
1. **Pass 1:** stream the manifest, compute per-chrom coverage quantiles incrementally. Constant RAM (no positions dataframe materialized).
2. **Pass 2:** for each chrom, call clusters with pre-computed thresholds in a `ProcessPoolExecutor` worker. Concat resulting cluster table at end.

**Correctness invariant:** the streaming quantile computation MUST agree with the legacy global quantile (within float tolerance). Test fixture must verify this on a non-trivial sample (≥10k positions across multiple chroms).

### 3.4 Acceptance

- [ ] Analyze on Han wt_R1 outputs file-by-file diff-equivalent to legacy outputs.
- [ ] Peak RAM during analyze cluster calling < 8 GB (versus current ~16-20 GB on Han wt_R1).
- [ ] Total analyze wall time ≥3× faster on a 16-core node.
- [ ] `<sample_id>.analyze.provenance.json` sidecar emitted; rerun returns SKIP.
- [ ] Cluster-calling two-pass quantile equivalence test passes.
- [ ] `pytest -m "not slow"` + `pytest -m slow` cDNA smoke pass.
- [ ] Opus reviews diff + cluster-calling two-pass logic + equivalence outputs.

### 3.5 LOC

~400 new (`streaming.py`), ~170 modified (`analyze_command.py` + `cdna_analyze_command.py` + `corrected_consensus.py`), ~250 tests. **Total ~800 LOC** — D is the biggest of the remaining commits.

---

## 4. Commit E — DRS poly-A trim + multi-aligner align + ONT cDNA stage-1

**Owner:** Sonnet for each sub-piece; Opus reviews and gates on Commit Zero scope.

**Conditional on Commit Zero results:** skip any sub-stage whose measured wall time is <10% of pipeline total.

### 4.1 DRS poly-A trim

- `rectify/core/commands/drs_trim_command.py:trim_drs_bam_polya` — currently single-threaded fetch loop. Region-parallel via `regions.py:plan_regions` + a worker that mirrors the per-read mutation pattern (lighter than `bam_writer_parallel` — just poly-A trim + metadata emission).
- Output: trimmed BAM + metadata parquet. Both shard-then-merge naturally.
- Emit `<sample_id>.drs_trim.provenance.json` sidecar.

### 4.2 Multi-aligner align

- `rectify/core/commands/run/stages.py:_run_alignment` — Commit Zero TODO: confirm whether 5 aligners run concurrently today (probably sequentially). If sequential, wrap in `ProcessPoolExecutor(max_workers=n_aligners)`. Each aligner is its own subprocess (samtools / mapPacBio / etc.), so this is process-level concurrency, not Python-thread concurrency.
- Sidecar shape decision (Opus): one sidecar per aligner (`<sample_id>.align.<aligner>.provenance.json`) OR one combined sidecar with `sub_outputs` listing per-aligner BAMs? Recommend per-aligner sidecars for clean skip-check (a single aligner failing shouldn't invalidate the others' work).
- **Heap-corruption hazard:** the merged-BAM lookups in `--aligner-bams` paths now also crash at chunk scale per `717e99e`. If `_run_alignment` opens the merged consensus BAM via pysam in workers, that's a known risk. Read the axis-2 briefing (`docs/briefings/axis2_htslib_state_corruption*.md` or similar) BEFORE designing the worker pattern.

### 4.3 ONT cDNA stage-1 consensus FASTQ

- `rectify/core/commands/cdna_correct_command.py` + `rectify/core/cdna/io.py:write_stage1_fastq` — stage-1 produces consensus FASTQ from UMI clusters; currently iterates clusters serially. `ProcessPoolExecutor` over UMI clusters; per-shard FASTQ outputs; concat at end. UMI clusters are independent (no cross-cluster state) so this is trivially correct.
- Emit `<sample_id>.cdna_stage1.provenance.json` sidecar.
- Tests: `tests/test_cdna_stage1_parallel_equivalence.py` — NEW, verifies stage-1 FASTQ output identical (sorted by cluster ID) between sequential and parallel paths.

### 4.4 Acceptance

- [ ] DRS poly-A trim wall ≥4× faster on 16-core node (only required if Commit Zero showed >10% of total wall on this stage).
- [ ] Multi-aligner align wall time ≈ `max(per_aligner_time)`.
- [ ] cDNA stage-1 FASTQ wall ≥4× faster on 16-core node.
- [ ] All three sidecars emitted; resume verified.
- [ ] cDNA smoke test passes.
- [ ] Opus reviews diff.

### 4.5 LOC

~250 new, ~120 modified. Three sub-pieces; can be split across multiple commits if cleaner.

---

## 5. Commit F — NETSEQ duplicated parallel path

**Owner:** Sonnet. Per CLAUDE.md "don't introduce abstractions beyond what the task requires" — duplicate ~150 LOC into `rectify/core/bam/netseq_bam_parallel.py` rather than generalize the abstraction across three arms with different output shapes.

### 5.1 Scope

`rectify/core/commands/netseq_command.py:run_netseq` and `rectify/core/bam/netseq_bam_processor.py` (640 LOC) — region-parallel via the same `regions.py:plan_regions` infrastructure but with NETSEQ-specific output (per-base 3'-end counts, not corrected BAM). No `realign_exon_blocks` involvement; no per-read CIGAR mutations beyond what NETSEQ already does.

### 5.2 Files

- `rectify/core/bam/netseq_bam_parallel.py` — NEW, ~150 LOC. Mirrors `bam_writer_parallel.py` shape but emits per-base count TSVs + handles NETSEQ exclusion stats. Shares only `RegionPlan` + `plan_regions()` from `regions.py`.
- `rectify/core/commands/netseq_command.py:run_netseq` — route default through `netseq_bam_parallel`. Small-input gate (same threshold as Commit B: 100 MB OR 500k reads). `--legacy-single-threaded` escape.
- `<sample_id>.netseq.provenance.json` sidecar.
- `tests/test_netseq_parallel_equivalence.py` — NEW. Both paths on a small fixture; compare 3'-end count TSVs + exclusion stats.

### 5.3 Acceptance

- [ ] NETSEQ on Han 2023 cluster_revisit test sample wall ≤ legacy/4 on 16-core node.
- [ ] Output TSVs identical (sorted by chrom, pos, strand).
- [ ] Exclusion stats identical.
- [ ] Sidecar emitted; resume verified.
- [ ] `pytest -m "not slow"` passes.

### 5.4 LOC

~180 new, ~60 modified.

---

## 6. Commit F.5 — Run-all resume end-to-end integration tests

**Owner:** Sonnet drafts the 16 scenarios from spec §F.5 (briefing draft was in `dev/specs/parallel_refactor_plan.md`); Opus designs the SIGTERM failure-injection harness (correctness-critical).

**Hard prerequisite:** every prior commit (B, restore_polya cleanup, D, E, F) has emitted its stage sidecar. F.5 verifies the resume protocol end-to-end.

### 6.1 Scope

`tests/integration/test_run_all_resume_drs.py` + `test_run_all_resume_cdna.py` + `test_run_all_resume_netseq.py`. Use small DRS/cDNA/NETSEQ fixtures (~50k reads each, run in ~30s end-to-end). 16-scenario matrix per spec §F.5 covering:
1-10: clean rerun, interrupt-at-each-stage, mid-stage interrupt, input mutation, ignored-argv mutation, load-bearing-argv mutation, force flags, git_sha mismatch, corrupt sidecar, missing output.
11-16: sample directory moved, $SCRATCH drift, env var unset, env var unset AND cached missing, transient sub_outputs missing, cross-cluster sidecar.

### 6.2 Failure-injection harness (Opus designs)

The hard part: a `with sigterm_after(seconds=N)` context manager that sends SIGTERM to the rectify subprocess after N seconds, captures partial state, then validates resume-from-partial. Must be robust to:
- pytest worker isolation
- platform differences (macOS vs Linux SIGTERM handling)
- BAM index file partial state (`.bai` half-written)

### 6.3 Acceptance

- [ ] All 16 DRS resume scenarios pass.
- [ ] cDNA + NETSEQ scenarios pass.
- [ ] Test suite runs in <5 min wall (acceptable for `slow` marker).
- [ ] `--dry-run-resume` output human-readable: per-stage decision + reason + RUN diff.

### 6.4 LOC

~400 new test code.

---

## 7. Commit G — Docs + CHANGELOG + deprecation roadmap

**Owner:** Sonnet drafts; Opus + user finalizes wording.

### 7.1 Files

- `CHANGELOG.md` — describe the parallelism additions, the manifest-as-canonical-artifact change, `--legacy-single-threaded`, `--emit-merged-tsv` default flip, new resume model + sidecar schema + CLI flags, `rectify export-merged-tsv` shim.
- `docs/PARALLELISM.md` — NEW. Two parallelism axes, merge points, M1 vs cluster sizing, override flags.
- `docs/RESUME.md` — NEW. User-facing doc on resume model: when stages skip vs run, sidecar schema (`<sample_id>.<stage>.provenance.json`), `--force-*` flags, `--dry-run-resume`. Include a worked example: "you killed the run after correct finished, here's what `run-all` resumes."
- `docs/ALIGNER_RECOMMENDATIONS.md` — if Commit Zero / E surfaces aligner-specific timing tradeoffs, update.
- `dev/specs/parallel_refactor_plan.md` — mark all commit checkboxes DONE; archive to `dev/specs/archive/` once Commit G lands.
- `rectify/AGENT_FIXES.md` — finalize the heap-corruption entry based on whatever the axis-2 investigation produced.

### 7.2 Acceptance

- [ ] mkdocs build clean.
- [ ] CHANGELOG references the original issue (`ISSUE_parallel_bam_write_and_lazy_tsv.md`).
- [ ] Master spec marked DONE.

### 7.3 LOC

~250 doc.

---

## 8. Time budget across C–G

| Commit | Wall (incl smoke) | Agent work | Notes |
|---|---|---|---|
| C remainder (work queue) | ~1.5 days | ~3 hr Opus + ~2 hr Sonnet + smoke wall | |
| restore_polya cleanup | ~0.5 day | ~2 hr Sonnet | Fold into D if cleaner |
| D (analyze + merge_corrected_tsvs) | ~2 days | ~4 hr Opus + ~4 hr Sonnet + smoke wall | Biggest remaining commit |
| E (drs_trim + align + cdna_stage1) | ~1.5 days | ~5 hr Sonnet + 1 hr Opus + smoke wall | Three sub-pieces |
| F (NETSEQ) | ~0.5 day | ~2 hr Sonnet | |
| F.5 (resume integration tests) | ~1 day | ~3 hr Sonnet + ~1.5 hr Opus harness | Needs B/D/E/F sidecars in place |
| G (docs) | ~0.5 day | ~2 hr Sonnet + 30 min Opus finalize | |
| **Total** | **~7-8 days wall** | **~30 hr agent** | One commit at a time |

---

## 9. Stop-and-ask triggers (applies to any commit)

Halt and report to Opus / user immediately if:

- Any `git log --oneline -3` HEAD doesn't match what this briefing assumes (HEAD has moved during your work; rebase/handle).
- The axis-2 heap-corruption investigation produces new findings that contradict the design here (e.g., the bug turns out to be in `bam_writer_parallel.py` itself rather than in worker fetch patterns).
- An equivalence test fails and you can't quickly localize the divergence.
- `pytest -m "not slow"` regresses for a reason not explained by your changes.
- Cluster ControlMaster dies mid-work (Sherlock Duo prompt suddenly appearing).
- Any uncommitted parallel-session WIP looks like it conflicts with your commit's files — STOP, surface it, don't silently overwrite.

---

## 10. Quick-start checklist for a follow-up agent

When you pick this up:

1. ☐ Read this briefing in full (~15 min).
2. ☐ Read `dev/specs/parallel_refactor_plan.md` §4 entries for the commit you're picking (~15 min).
3. ☐ Read `rectify/AGENT_FIXES.md` for the latest heap-corruption status (~5 min).
4. ☐ Read the axis-2 htslib briefing if your commit touches `--aligner-bams` workers (~10 min).
5. ☐ Probe current state: cross-cluster HEADs + `git status -s` + Sherlock ControlMaster.
6. ☐ Pick ONE commit. Don't try to land multiple in one session.
7. ☐ Follow that commit's §1-§7 section above (scope, files, acceptance, LOC).
8. ☐ Run the acceptance criteria before claiming done.
9. ☐ Report back per the standard 7-item structure (git status, diff stat, pytest, per-test outcomes, smoke if applicable, sidecar verification, TODOs).

**End of briefing.**
