# Session Handoff — Commits Zero / A / A.5 / B landed; structural smoke deferred

**Date:** 2026-05-20 (session ran 2026-05-19 evening → 2026-05-20 morning Pacific)
**Branch:** `drs-validation-rebuild`
**M1 HEAD:** `ba61d31 feat(parallel): Commit B - DRS default flip + first stage sidecar wiring`
**Cross-cluster sync:** M1 = GitHub = H2 = Sherlock all at `ba61d31` (verified at session end)
**Prior handoff:** `handoffs/HANDOFF_2026-05-19_manifest_correction_divergence.md`
**Next session expected:** tomorrow evening (Pacific) per Kevin

---

## 1. What was done

**Parallelism refactor — Commits Zero through B (load-bearing 60% of `dev/specs/parallel_refactor_plan.md` scope):**

- **Commit Zero — Profile attribution:** py-spy native+raw decomposition of `write_corrected_bam` on QSrev 200k slice + DRS 20k slice. Result in `dev/specs/profile_results.md` (committed in `16148a4`). Decision recorded: PROCEED with Commits A-F as written; BGZF `compresslevel=1` and batched `set_tag` escape hatches both REFUTED. Discovered the **protocol-conditional** payoff — QSrev WRITE phase 65% (realign 88% within); DRS WRITE phase 0.5% (correct phase 98% dominates).
- **Commit A — Shared parallel BAM-write infrastructure** (`39fbb63`). New: `bam_writer_parallel.py` (region-parallel write with idempotent `.ok` sentinels + sort-then-merge), `tsv_partition.py` (manifest + RegionTsvWriter), `regions.py:RegionPlan + plan_regions()`, `junction_scoring.py:build_junction_pool` ProcessPoolExecutor, `tests/utils/bam_compare.py`, smoke test parametrized over n_threads ∈ {1,2,4}.
- **Commit A.5 — Provenance + resume infrastructure** (`0379dde`). New `rectify.core.provenance` package: `path_resolver.py` (PortablePath envelope with sample-relative / env-relative / absolute + cross-cluster cached_absolute fallback + invariant against non-transient `$L_SCRATCH`/`$TMPDIR`), `skip_check.py` (4-gate decision tree + SkipDecision), `sidecar.py` (atomic write + schema validation), `cluster.py`, `hashing.py`, `cli.py`. Migration decision recorded as `COEXIST + cleanup` — kept `rectify.utils.provenance.ProvenanceTracker` (legacy run log), deleted dead `register_staged` method + callsites.
- **Commit B — DRS default flip + first stage sidecar wiring** (`ba61d31`). `correct_command.py` now defaults to region-parallel write for cluster-sized inputs (small-input gate at 100 MB OR 500k reads), emits `<sample_id>.correct.provenance.json` sidecar at stage exit, accepts `--legacy-single-threaded` escape + resume flags via `add_resume_args`. Manifest-only TSV output by default (Option B); `--emit-merged-tsv` opts back in. New `rectify export-merged-tsv <manifest>` shim subcommand (back-compat). `_load_corrections_from_tsv` auto-detects manifest vs single-TSV — all internal callers work transparently.

**In-session side-quests / coordination:**

- **Sherlock `PenaltyTableSet` resync** (`2a492f9`) — committed TestPenaltyTableSet (7 tests) + 5 per-UMI penalty TSVs that were locally-untracked on M1; cleared the "lost diff" concern from the other agent.
- **NaN-junctions tiebreaker bug** (`2894a61`) — `corrected_consensus.py:826` `_eff_key` was crashing on pandas NaN (truthy float) bypassing `or ''` guard; replaced with `isinstance(_juncs, str)` check. 2 previously-failing tests now pass.
- **`AGENT_FIXES.md` integration** — read the cross-agent coordination log, tempered the Commit B hypothesis entry (the `--checkpoint-dir + threads=8` mitigation ALREADY silent-hung at job 25510786; Commit B's pre-partitioning is meaningfully different but unproven). Saved memory `reference-agent-fixes-md` so future sessions read this file first.
- **CLAUDE.md updates** — added new section "Sherlock: direct SSH to active allocation (queue bypass)" documenting the SSH-from-login-to-compute-node pattern (Kevin's suggestion, lines 215-258 of `/Users/kevinroy/Library/CloudStorage/GoogleDrive-kevinrjroy@gmail.com/My Drive/Work/Chanfreau Lab/CLAUDE.md`).
- **Memory additions** (in `/Users/kevinroy/.claude/projects/.../memory/`): `feedback_sherlock_controlmaster_test.md` (probe ControlMaster at session start; don't assume blocked), `reference_sherlock_rectify_env.md` (Sherlock conda activation path: `source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh; conda activate rectify`), `reference_agent_fixes_md.md`.

**Documentation artifacts (committed in `16148a4`):**

- `dev/specs/parallel_refactor_plan.md` — master spec, 947 lines, with stage-level resume protocol and PortablePath envelope for cross-cluster sidecar correctness.
- `dev/specs/profile_results.md` — Commit Zero attribution + decision.
- `dev/specs/briefings/commit_zero_briefing.md`, `commit_a_briefing.md`, `commit_a5_briefing.md`, `commit_b_briefing.md` — Sonnet briefings for each commit.
- `docs/handoffs/HANDOFF_2026-05-20_corrected_consensus_nan_junctions_bug.md`, `HANDOFF_2026-05-20_sherlock_penaltytableset_resync.md` — handoffs for the two side-quests (both now resolved by `2894a61` and `2a492f9` respectively).

---

## 2. What's verified

- `pytest -m "not slow"` on M1 after Commit B → **1133 passed, 34 skipped, 4 deselected, 1 xfailed, 0 failed** in 41 s. Baseline was 966; +167 from new Commit A/A.5/B tests plus fixture migrations.
- `pytest -m slow` cDNA chain canary → **4 passed in 636 s** (after Commit B; verifies the manifest-only default doesn't break the cDNA pipeline path).
- **All 10 Commit B-specific tests pass** (`test_correct_command_parallel_default.py` 8 + `test_resume_correctness.py` 2). Includes the critical `test_sidecar_not_written_on_crash_before_sidecar` — verifies the §6.5 ordering invariant (sidecar emission MUST follow output durability).
- **8 Commit A tests** (3 `bam_writer_parallel_smoke` × n_threads + 5 `junction_scoring_parallel`) pass on M1.
- **64 Commit A.5 tests** (provenance/sidecar + skip-check + portable_path + cluster_detect + hashing) pass on M1.
- `assert_bams_equivalent` confirms parallel-vs-legacy correctness via `test_legacy_single_threaded_equivalent_output` — byte-identical modulo sort-tie ordering.
- `rectify export-merged-tsv --help` works; round-trip equivalence verified in `test_emit_merged_tsv_produces_both`.
- **Cross-cluster sync at session end:** `git log --oneline -1` on M1, GitHub, H2 (`/u/home/k/kevinroy/software/rectify`), and Sherlock (`/oak/stanford/groups/larsms/Users/kevinroy/software/rectify`) all return `ba61d31`.
- **Live-code check on all 3 clusters:** the NaN-junctions fix (and all Commit A/A.5/B code) is importable + runnable on M1, H2, and Sherlock — each cluster has rectify installed as an editable install pointing at the git working tree.
- **Sherlock chrI-V smoke (2,009,670 reads, 16 threads)** ran clean at 61:40 wall, exit 0 in a prior session at this scale — already validated Commit B's architecture below the heap-corruption threshold (per AGENT_FIXES.md scale table).

**NOT VERIFIED:**

- **Han wt_R1 6.7M-read structural-resolution smoke** — H2 pod_smp.q had 25,600 jobs waiting at submit time; Sherlock larsms saturated with Kevin's own `wt_R1_co` array (25515911) + owners had 7,555 queued. The full-scale outcome (A/B/C from briefing §4) for the AGENT_FIXES.md heap-corruption "STILL OPEN" entry remains untested. Commit B's `write_corrected_bam_parallel` pre-partitioning hypothesis is unverified at 6.7M-read scale.
- **`rectify restore-polya` sidecar emission** — only `correct` emits a sidecar in this commit. `restore_polya_command.py` still uses `pd.read_csv` directly on the corrected TSV (works with `--emit-merged-tsv` or the shim, but not manifest-aware natively). Deferred to Commit D/E.
- **`corrected_consensus.merge_corrected_tsvs` manifest-awareness** — same limitation. Deferred to Commit D/E.

---

## 3. Open items

- **Han wt_R1 6.7M-read structural smoke is the highest-value follow-up.** Tests whether Commit B's per-region pre-partitioning resolves OR mitigates OR fails-the-same-way as the AGENT_FIXES.md heap-corruption bug. **Why deferred:** queue saturation on both clusters at submission time; the prior `--checkpoint-dir + threads=8` mitigation already silent-hung at job 25510786 so the hypothesis is uncertain and needs a real test. Outcome buckets A/B/C are documented in `dev/specs/briefings/commit_b_briefing.md` §4 with explicit instructions on how to update AGENT_FIXES.md for each.
- **`restore_polya_command.py` not manifest-aware** + **`<sample_id>.restore_polya.provenance.json` sidecar not emitted.** **Why deferred:** Sonnet trimmed scope on Commit B to land the bigger `correct` stage cleanly; restore_polya stays on legacy path until Commit D/E. Users running `rectify restore-polya` after a default manifest-only `rectify correct` need `--emit-merged-tsv` or to run the shim first. Documented in Commit B's commit message as a Commit D/E follow-up.
- **`corrected_consensus.merge_corrected_tsvs` uses `pd.read_csv` directly** — same pattern, same workaround, same deferral to D/E.
- **Remaining parallelism commits (C, D, E, F, F.5, G).** Per `dev/specs/parallel_refactor_plan.md`:
  - C — Dynamic aligner-queue concurrency (Opus owns; ~250 LOC; orthogonal to the deferred Han wt_R1 test)
  - D — Analyze partial-streaming + `merge_corrected_tsvs` parallelism (Opus owns the cluster-calling two-pass; would also pick up the manifest-aware restore_polya + merge_corrected_tsvs TODOs above)
  - E — DRS poly-A trim + multi-aligner align + ONT cDNA stage-1 (Sonnet)
  - F — NETSEQ duplicated parallel path (Sonnet)
  - F.5 — Run-all resume end-to-end integration tests (Sonnet drafts; Opus designs failure-injection harness)
  - G — Docs + CHANGELOG + deprecation note
  **Why deferred:** session length + queue blockage on the structural smoke. Natural sequencing for tomorrow: either dispatch C while queues clear, OR clear the Han wt_R1 smoke first to resolve the AGENT_FIXES.md open bug, OR pick off the restore_polya/merge_corrected_tsvs TODOs in a small consolidation commit.
- **`tasks/#5 Phase Zero-B full Han wt_R1 stage timings`** remains pending. Originally intended as the "before measurement" for Commit B's speedup claim. **Why deferred:** same heap-corruption blocker as the structural smoke — the full-sample run is what would crash. When the structural smoke runs, Zero-B's per-stage timings come along for free.
- **Lots of Kevin's pre-existing WIP in the working tree** (figures, validation BAM rebuilds, scripts/calibration, docs/architecture, etc.) — not my scope; Kevin owns the disposition. Listed under §5 with `[uncommitted, Kevin's WIP]`.

---

## 4. Resume command

**Resume:** start by reading `dev/specs/parallel_refactor_plan.md` §0 (success metric is protocol-conditional) + `rectify/AGENT_FIXES.md` (cross-agent coordination, especially the "STILL OPEN" heap-corruption entry which Commit B may or may not resolve).

Then probe state in this order:

```bash
# 1. Cross-cluster head reconciliation (per CLAUDE.md rule). Should all return ba61d31 or newer.
cd /Users/kevinroy/work/rectify && git log --oneline -1
git ls-remote origin drs-validation-rebuild | awk '{print substr($1,1,7)}'
ssh hoffman2 'bash -lc "cd /u/home/k/kevinroy/software/rectify && git log --oneline -1"'
ssh sherlock 'bash --norc --noprofile -c "cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify && git log --oneline -1"'

# 2. Probe Sherlock ControlMaster (per memory feedback-sherlock-controlmaster-test).
ssh sherlock 'echo ok'   # should return within ~2s without Duo prompt

# 3. Check cluster queue depth (last seen: H2 disastrously backlogged, Sherlock larsms saturated by Kevin's own wt_R1_co).
ssh hoffman2 'bash -lc "qstat -g c | grep -E pod_smp; qstat -u kevinroy | head -5"'
ssh sherlock 'bash --norc --noprofile -c "squeue -u kevinroy | head -10; squeue -p larsms -h | wc -l"'
```

**Branch then by outcome:**

- **If Sherlock larsms is free** (Kevin's `wt_R1_co` 25515911 / 25514775 has drained):
  → kick off the **Han wt_R1 6.7M-read structural smoke** per `dev/specs/briefings/commit_b_briefing.md` §4 (run `rectify correct` on the full merged BAM at `-j 8 --tmp-dir $L_SCRATCH/rectify_regions`). Watch for outcome A/B/C; update `rectify/AGENT_FIXES.md` per the briefing's instructions; commit the update.

- **If H2 pod_smp.q has drained but Sherlock is still busy**:
  → run the smoke on H2 instead (same command; H2 was the original crash environment so this is the more authentic test). Same outcome handling.

- **If BOTH clusters are still saturated:**
  → pick from one of these productive alternatives, each independent of the smoke:
    - **Commit C** (dynamic aligner queue) — Opus-owned, ~250 LOC, orthogonal to the smoke. Read `dev/specs/parallel_refactor_plan.md` §4 Commit C; design pending.
    - **Restore_polya + merge_corrected_tsvs cleanup** — pick off the two §3 Open-item TODOs in a small consolidation commit (~50 LOC: switch their `pd.read_csv` direct reads to use the manifest-aware loader from `rectify.core.bam.bam_writer:_load_corrections_from_tsv`).
    - **Commit D briefing draft** — restructures `analyze_command.py` for partial-streaming; bigger lift; would benefit from the smoke result first because it informs the cluster-calling two-pass design.

- **If anything is unclear** (cluster sync drift, ControlMaster dead, unexpected work in `git status`):
  → STOP, re-read this HANDOFF, then ask Kevin before action. Especially: don't unilaterally commit anything in his pre-existing WIP listed in §5.

---

## 5. Files touched

**Committed in this session** (all on `drs-validation-rebuild`, pushed to origin, pulled on H2 + Sherlock):

- `rectify/cli.py` — register `export-merged-tsv` subcommand (`ba61d31`)
- `rectify/core/analyze/loaders.py` — manifest-aware `load_corrected_positions` (`ba61d31`)
- `rectify/core/bam/bam_writer.py` — manifest-aware `_load_corrections_from_tsv` (`ba61d31`)
- `rectify/core/bam/bam_writer_parallel.py` — NEW, region-parallel write infrastructure (`39fbb63`)
- `rectify/core/bam/regions.py` — `RegionPlan` + `plan_regions()` (`39fbb63`)
- `rectify/core/bam/tsv_partition.py` — NEW, manifest + RegionTsvWriter (`39fbb63`)
- `rectify/core/commands/correct_command.py` — dispatch + sidecar emission + resume flags + manifest-only default (`ba61d31`)
- `rectify/core/commands/export_merged_tsv_command.py` — NEW, back-compat shim (`ba61d31`)
- `rectify/core/consensus/corrected_consensus.py` — NaN-junctions fix at `_eff_key` (`2894a61`)
- `rectify/core/provenance/__init__.py`, `cli.py`, `cluster.py`, `hashing.py`, `path_resolver.py`, `sidecar.py`, `skip_check.py` — NEW package (`0379dde`)
- `rectify/core/splice/junction_scoring.py` — `build_junction_pool` ProcessPoolExecutor (`39fbb63`)
- `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna{,_umi1,_umi2,_umi3plus}.tsv`, `penalty_scores_qsrev.tsv` (`2a492f9`)
- `rectify/utils/provenance.py` — deleted dead `register_staged` method (`0379dde`)
- `rectify/utils/version.py` — NEW, `get_rectify_git_sha()` helper (`ba61d31`)
- `rectify/AGENT_FIXES.md` — Commit B hypothesis entry (`ba61d31`)
- `tests/test_bam_writer_parallel_smoke.py`, `test_junction_scoring_parallel.py`, `tests/utils/bam_compare.py`, `tests/utils/__init__.py` (`39fbb63`)
- `tests/test_cluster_detect.py`, `test_hashing.py`, `test_portable_path.py`, `test_provenance_sidecar.py`, `test_provenance_skip_check.py` (`0379dde`)
- `tests/test_correct_command_parallel_default.py`, `test_resume_correctness.py` — NEW (`ba61d31`)
- `tests/test_correct_command_drs.py`, `test_quantseq_rev_integration.py`, `test_validation_reads.py`, `test_junction_refiner.py` — fixture migrations / TestPenaltyTableSet additions (`2a492f9`, `ba61d31`)
- `dev/specs/parallel_refactor_plan.md`, `profile_results.md`, `briefings/commit_{zero,a,a5,b}_briefing.md` (`16148a4`, `ba61d31`)
- `docs/handoffs/HANDOFF_2026-05-20_corrected_consensus_nan_junctions_bug.md`, `HANDOFF_2026-05-20_sherlock_penaltytableset_resync.md` (`16148a4`)
- `handoffs/HANDOFF_2026-05-19_manifest_correction_divergence.md` — `[uncommitted]` archive of the prior HANDOFF before writing this one

**Edited outside the rectify repo:**

- `/Users/kevinroy/Library/CloudStorage/GoogleDrive-kevinrjroy@gmail.com/My Drive/Work/Chanfreau Lab/CLAUDE.md` — added "Sherlock: direct SSH to active allocation (queue bypass)" section (lines 215-258). `[uncommitted; lab-level doc, Kevin's call whether to commit/push elsewhere]`
- `/Users/kevinroy/.claude/projects/.../memory/feedback_sherlock_controlmaster_test.md`, `reference_sherlock_rectify_env.md`, `reference_agent_fixes_md.md` — new memory entries (auto-persisted, no commit needed)
- `/Users/kevinroy/.claude/projects/.../memory/MEMORY.md` — index updated with the three new entries

**Working tree at handoff time** (NOT my scope unless flagged; mostly Kevin's WIP from other workstreams running in parallel):

- `[uncommitted, Kevin's WIP]` `AGENT_FIXES.md` — empty-SEQ filter recovery section (lines 69-103). Documents the workaround for tainted pre-`e8c8070` mapPacBio BAMs. Worth Kevin committing.
- `[uncommitted, Kevin's WIP]` `CLAUDE.md` (the rectify-local one, distinct from the lab Drive one I edited)
- `[uncommitted, Kevin's WIP]` `README.md`, `TODO.md`, `dev/BUGS_TO_FIX.md`, `docs/ALIGNER_RECOMMENDATIONS.md`, `docs/troubleshooting.md`
- `[uncommitted, Kevin's WIP]` `docs/figures/*.svg` + `*.png` (12 figure regenerations)
- `[uncommitted, Kevin's WIP]` `docs/images/*.png` (5 deletions)
- `[uncommitted, Kevin's WIP]` `rectify/core/consensus/corrected_consensus.py` — Kevin added `_cigar_hp_edit_distance` + `_cigar_aligned_bases` helpers since I last committed this file; not my scope to touch
- `[uncommitted, Kevin's WIP]` `rectify/data/validation/aligners/*.bam`, `PROVENANCE.json`, `README.md`, `VALIDATION_READS.md` (validation BAM rebuild wave)
- `[uncommitted, Kevin's WIP]` `scripts/calibration/*`, `scripts/validation_data/*` (calibration + validation utilities)
- `[uncommitted, Kevin's WIP]` many untracked: `_h2_rebaseline/`, `README_KR_edits.md`, `RECTIFY_SHERLOCK_HANDOFF_*.md`, `docs/aligners/`, `docs/architecture/`, `docs/handoffs/cdna_umi_*.md`, `docs/penalty_tables_quickref.md`, `docs/protocols/`, `pA_tail_DRS_*.md`, `handoffs/HANDOFF_2026-05-18_*.md`, etc.

If next session: **do not commit any `[uncommitted, Kevin's WIP]` items without explicit Kevin OK.** They are visible in `git status -s` and may look stale, but they are active workstreams.

---
