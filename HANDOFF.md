# Handoff — drs-validation-rebuild (session 2026-05-24 → 25)

**Date:** 2026-05-25
**Branch:** `drs-validation-rebuild`
**HEAD:** `93d99b8` (this doc was written at `961c844`; the perf thread has since landed — see banner below)

> ✅ **PERF THREAD LANDED (2026-05-25, perf agent) — supersedes the "concurrent
> agent / do NOT commit" warnings in this doc.** The perf work referenced
> throughout is now committed and verified:
> - `ed3df74` pool anchor floor + cross-family concordance relaxation;
>   `93d99b8` docs. **HEAD is now `93d99b8`.**
> - `rectify/core/splice/junction_scoring.py` and
>   `tests/test_junction_anchor_filter.py` are **COMMITTED** (no longer in-flight —
>   the do-NOT-commit warnings below are resolved for these two).
> - Result: Module 2F 3'SS-rescue candidate-collapse **155×–911× across all 5
>   aligners**, 66–75% of reads now skip rescue entirely. No real-junction loss
>   (floor drops only single-family junctions; relaxation recovers ≥2-family;
>   `minimap2`+`gapmm2` count as one family since gapmm2 wraps minimap2).
> - Periodicity/complexity dimension **DEFERRED to the Human DRS agent**
>   (`dev/handoffs/HANDOFF_2026-05-25_human_drs_validation.md` §7). **No splice-motif
>   gating** in rectify (memory: `feedback_no_motifs_unbiased_discovery`).
> - Full write-up: `dev/PERF_AUDIT.md` (resolution + 10k verify numbers).
> - **Still OPEN — this doc's actual focus, untouched by the perf thread:** the
>   upf1D validation-set doubling and master-plan execution below.

> ⚠️ **UNCOMMITTED WORK IN TREE — perf agent wrapped up; you are now the SOLE
> agent.** As of 2026-05-25 the perf/aligner agent finished and left its work
> UNCOMMITTED on top of HEAD `93d99b8`: `git status` shows modified **code**
> (`multi_aligner.py`, `align_command.py`) and **docs**
> (`docs/aligners/*`, `ALIGNER_RECOMMENDATIONS.md`) plus dev artifacts (full list
> in §5). Nobody is editing it now, but it is NOT your task. **Stage ONLY your own
> files with explicit paths; NEVER `git add -A`/`.`. Recommend confirming with
> Kevin whether to commit the leftover perf/aligner work first, so you start from
> a clean tree.**

> 🎯 **Why this handoff exists (Kevin, 2026-05-25):** the run-all manifest test
> hit a perf bug, and chasing it consumed the whole session. **That perf thread
> has since LANDED (top banner) — it is done, not your blocker.** The **original
> session goals — the upf1D validation-set doubling and master-plan execution —
> were NOT completed.** This handoff (for a fresh Opus agent) deliberately
> re-centers on those: the perf detail below is context, your task is §3/§4.

---

## 1. What was done

- **Release strategy decided: stay on v0.9.x; v1.0.0 deferred to beta-cohort
  feedback.** Recorded in memory (`project_rectify_release_strategy`); plan doc
  retitled to "0.9.x hardening" (`6c1d4ac`). Manifest mode flagged experimental.
- **Module 2F terminal peel** shipped earlier in the session (`32c6f2e`,
  `2718c07`), then the perf cascade it surfaced at production scale:
  - peel cost-bound, eliminates the OOM (`8e8dc8c`)
  - 3'SS-rescue candidate narrowing + pool-radius shrink (`bd20f9e`)
  - dual-site interval-tree fetch — intron-length-independent, human-ready
    (`961c844`)
  - pool anchor floor + cross-family concordance relaxation (`ed3df74`, perf
    agent 2026-05-25) — candidate-collapse 155×–911×; see top banner + §2
- **Diagnosed the run-all manifest hang.** py-spy (29,999 samples) → ~87% CPU in
  `_hp_edit_distance` from **per-read scans of the full junction pool**. Three
  root causes: (a) per-read global-collection scan, (b) vestigial 10 kb radius
  (sized for a removed `junction_proximity_bp=5000`), (c) bound on the wrong axis
  (radius on `intron_start` bounds intron *length* → breaks for human). Write-up:
  `AGENT_FIXES.md` [2026-05-24] + the perf-audit doc.
- **NEW-080** minus-strand artifact-junction flank fix (`cc86cc0`); backlog/plan
  reconciliation (`64fd716`, `4c05e07`).

## 2. What's verified

- **Full fast suite green at HEAD `93d99b8`:** `pytest -m "not slow"` →
  **1315 passed, 35 skipped, 4 deselected, 1 xfailed** (~5 min). (The earlier
  focused 218-test run at `961c844` is also still green.)
- OOM eliminated: live `/proc` probe of the full-pool run showed RSS **64 GB → ~14 GB**.
- **VERIFIED — candidate-collapse + no real-junction loss (perf agent, 2026-05-25).**
  On 10k 5-aligner data with committed code, dual-site fetch + anchor floor cut
  rescue candidates/read from ~390 to 0.4–2.5 (**155×–911×**; 66–75% of reads now
  skip rescue entirely). Pool-composition equivalence shows the floor drops ONLY
  single-family junctions and the relaxation recovers ONLY ≥2-family ones → zero
  real-junction loss. The earlier 800-read end-to-end correct verify showed all
  dropped rescues spurious.
- **STILL NOT RUN — the named `dev/verify_dualsite_correct.sbatch` full-pool
  completion (<32 GB) run, and the exact end-to-end "793 s → X" wall-time.** The
  above strongly implies completion (candidate-driven HP-ED work near-eliminated),
  but the specific full-pool wall/RSS job was not executed — optional gold-standard
  confirmation.

## 3. Open items  ← the original goals live here

- **[ORIGINAL GOAL — NOT DONE] upf1D validation-set doubling.** Add **2 more +/-
  strand reads per category** (9 categories) to `rectify/data/validation/`,
  sourced from **upf1D** (NMD-deficient → rich in novel/alternative intron usage;
  PMID 24722551, Kawashima/Chanfreau 2014). Choose reads where a `rectify
  correct` step was load-bearing for the best alignment, reflecting a diversity
  of winning aligners; key by XV tag. cat5 candidate loci hint: 2-intron genes
  YGL076C / YBR111W-A / YGL033W. **Blocked on:** upf1D being re-trimmed/re-mapped
  with CURRENT rectify (existing processed upf1D is stale-trim-era) — which is
  exactly what the run-all manifest test below was for. **Why it stalled:** the
  manifest test hung on the perf bug and the session became perf debugging.
- **[ORIGINAL GOAL — PARTIALLY DONE] Master-plan execution.**
  `dev/specs/v1_0_0_master_plan_20260523.md` has parallel (P1–P7, Sonnet) +
  sequential (S1–S9, Opus) buckets. The plan was created/reconciled and NEW-080 +
  2F landed, but most P/S items (e.g. P1 doc reconciliation, P2 NEW-079 off-by-1)
  are unexecuted. Why deferred: the perf detour.
- **run-all manifest test (wt+upf1 5% subset)** — the vehicle to reprocess upf1D.
  Trimmed+aligned BAMs already staged on Sherlock at
  `$SCRATCH/rectify_wt_by4742_rep1_25846844_0/`; the previously-failing
  inline-correct stage is what the perf fix targets. Re-attempt once the
  dual-site fix is scale-verified.
- **Perf-audit follow-ups** — the primary 2F/HP-ED hotspot is **LANDED** (`ed3df74`,
  see §2). Remaining secondary hotspots are identified + documented in
  `dev/PERF_AUDIT.md` (Part III; the 2026-05-24 findings were merged in and the
  separate FINDINGS file deleted) but NOT fixed: the 28 s non-HP-ED
  softclip-realignment outlier; the per-aligner setup / `--annotation` rebuild
  (Finding #3); `2C` indel correction on human. Periodicity/complexity dimension
  DEFERRED to the Human DRS agent (`dev/handoffs/HANDOFF_2026-05-25_human_drs_validation.md`
  §7). `_real_junctions` pool-floor bypass is a documented open item. **No
  splice-motif gating in rectify** (memory: `feedback_no_motifs_unbiased_discovery`).
- **[CORRECTNESS — perf agent, 2026-05-25] mapPacBio `intronlen`/`maxindel`
  misconfigured.** `run_map_pacbio` (`multi_aligner.py` ~L723) sets
  `intronlen=max_intron` — backwards: `intronlen` is a MINIMUM D→N relabel
  threshold, not a max — and leaves `maxindel` at BBMap's 16 kb default. On human
  chr5 DRS, mapPacBio emitted ~1k introns vs 200k–418k for the other aligners.
  The 2026-05-24 commit `b83e537` ("wire --max-intron into mapPacBio intronlen")
  was a misread and made labeling WORSE; masked for yeast (short introns). Fix
  (TODO): `intronlen=10`–`20`, `maxindel=max(200000, max_intron)`, fix the
  comment. Full entry: `AGENT_FIXES.md` [2026-05-25]. Relevant to the upf1D
  remapping above, though masked for yeast so not a hard blocker for it.
- **2F mapPacBio case-1 screening** (`RECTIFY_PEEL_SCREEN`) — in progress,
  deferred behind perf. Goal: find reads where the peel gets a strictly-better
  alignment (hunch: mapPacBio).
- **Commit-Zero-B profiling** (size perf items C-II.2 / D / F) — pending.

## 4. Resume command

**Resume by re-centering on the ORIGINAL goals. The perf thread has LANDED
(`ed3df74`/`93d99b8`, committed + verified per §2) — it is no longer a blocker.**

1. **Perf gate is effectively clear.** The dual-site + anchor-floor fix is committed
   and the candidate-collapse / no-real-loss verification is done (§2). The only
   remaining named check is the full-pool completion run
   (`dev/verify_dualsite_correct.sbatch`) — run it if you want the exact full-pool
   wall/RSS, but it is NOT required to proceed.
2. **Re-submit the run-all manifest test** (wt+upf1 5% subset) — the inline-correct
   stall it hit is exactly the bug the landed fix targets. **First sync the cluster
   rectify to `ed3df74`+** — H2/Sherlock are non-git rsync copies currently BEHIND
   this commit, so an unsynced run would execute the OLD pre-fix code. Staged
   inputs: Sherlock `$SCRATCH/rectify_wt_upf1_manifest_test_20260523/` (sbatch +
   manifest) and `$SCRATCH/rectify_wt_by4742_rep1_25846844_0/` (trimmed/aligned
   BAMs). On completion, upf1D is freshly processed → do the **validation-set
   doubling**: per category in `rectify/data/validation/VALIDATION_READS.md`, pick
   2 new +/- reads from upf1D where correction was load-bearing.
3. **Then** return to master-plan execution: open
   `dev/specs/v1_0_0_master_plan_20260523.md`, run the P1–P7 parallel bucket and
   the un-done S-bucket items.

## 5. Files touched (this session)

Committed (through HEAD `93d99b8`):
- `rectify/core/splice/splice_aware_5prime.py` — 2F peel + candidate narrowing
  (`8e8dc8c`, `bd20f9e`, `2718c07`, `32c6f2e`, `cc86cc0`)
- `rectify/core/bam/bam_processor.py` — dual-site interval-tree fetch
  (`bd20f9e`, `961c844`)
- `rectify/core/splice/junction_scoring.py` + `tests/test_junction_anchor_filter.py`
  + `tests/test_junction_scoring_parallel.py` — pool anchor floor + cross-family
  concordance relaxation + graded periodicity helpers (`ed3df74`, perf agent)
- `dev/PERF_AUDIT.md` (resolution + 10k verify) +
  `dev/handoffs/HANDOFF_2026-05-25_human_drs_validation.md` §7 (`93d99b8`)
- `AGENT_FIXES.md` — diagnosis entry + perf-audit pointer
- `dev/specs/v1_0_0_master_plan_20260523.md` — 0.9.x retitle (`6c1d4ac`)
- `tests/test_terminal_peel.py` — peel + cost-bound tests (`8e8dc8c`)

Memory (not in repo): `project_rectify_release_strategy.md` + `MEMORY.md` index.

`[uncommitted]`:
- `HANDOFF.md` — this file. Stale prior handoff archived to
  `handoffs/HANDOFF_2026-05-23_04effa5.md`.

**Perf/junction work — COMMITTED at `ed3df74`/`93d99b8`** (no longer in-flight):
`rectify/core/splice/junction_scoring.py`, `tests/test_junction_anchor_filter.py`,
`tests/test_junction_scoring_parallel.py`, `dev/PERF_AUDIT.md`,
`dev/handoffs/HANDOFF_2026-05-25_human_drs_validation.md`.

> ⚠️ **WORKING TREE IS NOT CLEAN — perf agent's leftover uncommitted work (do NOT
> sweep into your commits).** The perf/aligner agent wrapped up but left these
> uncommitted on top of HEAD `93d99b8` (`git status`, 2026-05-25):
> - modified code: `rectify/core/align/multi_aligner.py`,
>   `rectify/core/commands/align_command.py`
> - modified docs: `docs/ALIGNER_RECOMMENDATIONS.md`,
>   `docs/aligners/{minimap2,gapmm2}.md`
> - untracked: `docs/aligners/{deSALT,mapPacBio,minisplice,winnowmap2}.md`,
>   `dev/analyze_ab.py`, `dev/verify_dualsite_correct.sbatch`,
>   `dev/specs/briefings/`,
>   `dev/validation_review/`, `dev/branding/`
> - also modified (shared docs): `AGENT_FIXES.md`, `dev/PERF_AUDIT.md`,
>   `HANDOFF.md` (this file)
>
> These are NOT yours. Stage only the files your validation-doubling / master-plan
> work creates or changes, with explicit `git add <paths>`, atomic stage+commit.
