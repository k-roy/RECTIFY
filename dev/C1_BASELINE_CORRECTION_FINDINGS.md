# FINDINGS — the native-aligner "C1 build" task was already done; I was on the wrong branch

**Written:** 2026-07-02 (session that was asked to execute `dev/HANDOFF_NATIVE_ALIGNER_BUILD.md`).
**Status:** No code shipped. A wrong-baseline reimplementation was built, caught before any commit,
and fully discarded. This doc is the durable record so the correction is not lost.

---

## TL;DR

- **C1 (the "inner −logP scorer" the handoff asks to build) is ALREADY BUILT, LOCKED, and VALIDATED**
  on branch `worktree-agent-a25a2c1e784ad37dc` (checked out at
  `.claude/worktrees/agent-a25a2c1e784ad37dc`). See `dev/C1_DESIGN.md` (locked 2026-06-29) +
  `tests/test_c1_lengthlaw.py` + the committed `align_exon_block_global(penalty_table, lam, ins_lengthlaw)`.
- **I was dropped into the MAIN worktree `/Users/kevinroy/work/rectify` on branch `drs-validation-rebuild`,
  which has C1 REMOVED.** I read that branch's `local_aligner.py` (no C1), concluded C1 was unbuilt, and
  built a reimplementation — on the branch the handoff says NEVER to commit to.
- **The reimplementation was wrong on every axis** (see below) and would have reintroduced a bug real SIRV
  already caught. It was discarded; the main worktree is restored pristine; nothing under
  `.claude/worktrees/` was touched.

## Why my reimplementation was wrong (all three refuted by `dev/C1_DESIGN.md` on the agent branch)

1. **Wrong cost source.** I fed `HpPenaltyTable.del_cost/ins_cost` (the `penalty_score` column) into the DP.
   C1_DESIGN.md §"The cost": *"`penalty_score` is NOT −logP … Do NOT feed `penalty_score` into the DP."*
   The correct signal is **`rate_mean`** (sums to ≈1 over {M,X,D,I} per column = a calibrated emission dist).
2. **Wrong mechanism.** I replaced the per-position gap-**EXTEND** term. C1_DESIGN.md proves
   *"scale-mixing is provably incoherent in a global, full-consumption DP"* (match-count varies across paths).
   The only coherent partial is a **baseline-anchored log-odds DELTA on gap-OPEN**:
   `Δ_D(hp,bc) = λ·ln(P_D(hp,bc)/P_D(1,bc))`, zero at hp=1, λ≈1, gated on `homo_mask`. Already implemented.
3. **Reintroduced a known-harmful bug.** I applied the insertion law unconditionally. C1_DESIGN.md
   §"OVER-CALL RESULT": real SIRV (LRGASP RNA002 + SG-NEx) showed the ins delta **hallucinates indels
   3–7%, growing with run length**; it is **GATED OFF by default** (`ins_lengthlaw=False`) with a
   regression test. My version had no such gate.

## Branch / worktree topology (the trap)

| | path | branch | C1 state |
| --- | --- | --- | --- |
| my cwd (session start) | `/Users/kevinroy/work/rectify` | `drs-validation-rebuild` (9c183f6) | **REMOVED** |
| the handoff's branch | `.claude/worktrees/agent-a25a2c1e784ad37dc` | `worktree-agent-a25a2c1e784ad37dc` (8dccd4e) | **BUILT + validated** |

- The two branches have **diverged** (neither is an ancestor of the other). `local_aligner.py` and
  `hp_penalty.py` differ between them — the agent branch has the C1 code, drs does not.
- The session-start `gitStatus` claimed the current branch was `worktree-agent-a25a2c1e784ad37dc`, but the
  actual cwd HEAD is `drs-validation-rebuild`. That discrepancy is the whole trap.
- The agent worktree has **active, recent uncommitted WIP** (files touched 2026-07-02 ~21:58, minutes before
  the 22:34 handoff): `scripts/benchmark/{addressability_precheck,human_noncanon_sim,real_diff_composition,
  tie_logp_probe}.py`, `dev/DIRECTOR_ALGO_EVAL_FABLE*.md`. In this heavily-concurrent setup (~12 worktrees,
  some locked, one-agent-per-workspace) it may be another agent's live territory — do NOT edit/commit there
  without confirming it is safe.

## Corrected "done / open" map (from code + commits + the C-series design docs, NOT the prose handoff)

**DONE (committed on the agent branch):**
- C1 HP-length-law gap cost (locked 2026-06-29): gap-OPEN log-odds delta from `rate_mean`, deletion-only
  (ins gated off), byte-identical when `penalty_table=None`, Claim-A proven (in-run vs out-of-run placement:
  boundary_sub 0.00→0.55/0.78), clean-run false-indel ≈0. Real-SIRV over-call SAFETY result done.
- A large investigation + benchmark apparatus (`scripts/benchmark/c1_*.py`, `c2_gate.py`…`c6_headroom.py`,
  panel/blindspot scripts) and design docs `dev/C{1..6}_DESIGN.md`, `dev/ALIGNER_MEMBER_DESIGN.md`,
  `dev/NATIVE_ALIGNER_OVERVIEW.md`, `dev/HANDOFF_ALIGNER_BENCHMARK.md`.

**OPEN (the handoff's genuine, and DANGEROUS, next steps):**
- **Outer junction re-aligner** (the "defining purpose"): enumerate candidate junctions in the panel window,
  zero the motif term, neutralize annotation priors, score flanks with the inner scorer, pick best. This
  layer MOVES junctions → discards C1's boundary-invariance → owns a stricter **false-JUNCTION** FDR (the C6
  variant-adjacent trap). NOT covered by C1's false-indel tests.
- **Make-or-break, do-no-harm half (real-data-testable NOW):** snap-or-hold on the 3 cross-platform anchors
  (SLC35A4/TMED9/SQSTM1; `dev/COMPASS_2corroborated_CROSSPLATFORM.md`). Either snap or hold is a PASS if
  defensible; the only failure is silently inheriting the panel's choice.
- **Make-or-break, addressability half (C>B) — UNVALIDATED, no real-data substrate.** The handoff forces a
  named choice: (1) build an INDEPENDENT-error-model sim (SIRV/published ONT profile, NOT RECTIFY's own
  tables — else C wins by construction), or (2) label arm C explicitly "build-to-understand, unvalidated."
- C1_DESIGN.md's own open item: the C1 code was PRE-`aligner_bench_live`; needs syncing to Sherlock for the
  real-SIRV / cluster benchmarks.

## Recommendation / questions for the PI (Kevin)

1. **Where should the native-aligner work happen?** The handoff's branch (`worktree-agent-a25a2c1e784ad37dc`)
   has active recent WIP — is that worktree free for me, or is another agent on it? I was placed in the main
   worktree on the forbidden branch, which itself needs clarifying.
2. **Given C1 is done, the real next step is the outer junction re-aligner + the addressability make-or-break
   — which the handoff itself flags as unvalidated/dangerous and requires a named honest path.** Confirm
   whether to proceed with that (and which addressability path: independent-error sim vs labeled
   build-to-understand), or a different priority (e.g. the snap-or-hold do-no-harm test first, which IS
   real-data-testable now).
