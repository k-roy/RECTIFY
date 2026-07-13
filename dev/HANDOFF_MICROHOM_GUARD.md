# Handoff — Microhomology-drift guard (session 2026-07-11)

**Branch:** `worktree-agent-a25a2c1e784ad37dc` (worktree). Guard is COMMITTED here
(`3a716aa` feat + `17e57ba` tests/geometry-fix). **NEVER commit to `drs-validation-rebuild`.**
**Constraint:** microhom_drift_margin default stays 0.0 (OFF, byte-identical) until the
triple-audit is all-clear AND the COMPASS real-data threshold is confirmed.

---

## Done & verified
- **Built the microhomology-drift guard** — the general (non-HP) analog of the shipped HP-drift
  guard. New pure helpers `_frac_match`, `_move_microhomology` in `junction_refiner.py`; params
  `microhom_drift_margin` (op point 8.0), `microhom_threshold` (0.5) threaded through
  `refine_read_junctions → _run_sequential/_run_parallel (worker kwargs) → refine_bam_junctions`.
  Byte-identical when margin=0 (default).
- **Ground-truth validated** (spike-in fab panel + R3 discovery panel, per-read truth, seed excluded):
  - Fabrication FDR (canonical→drift): **1.31% → 0.00%** at m=8 (m=3 → 0.05%, 96% cut).
  - R3 discovery recall: **0.284 → 0.284 PRESERVED EXACTLY** (b=0, c=0) at both m=3 and m=8.
  - Canonical (plain) recall: **+0.010** (0.931→0.941) — guard also prevents spurious canonical drift.
  - Composition with HP-guard on `mix_fair_out` (HP m=3 + microhom m=8): no regression, slight
    improvement (ACC_A_D0 +0.007, BOT_A_D1 +0.005). Zero discovery cost — the HP-guard bar.
- **Byte-identity at scale:** full `not slow` suite at default = **1653 passed** (== pre-guard
  baseline; the 1 error is the known pre-existing missing-fixture file, present before the guard).
  64 refiner/guard tests + 14 new `tests/test_microhom_drift_guard.py` all pass.
- **COMPASS BBMap survived the outage:** job 33531828 COMPLETED, 4/4 `*.bbmap.junc.tsv` present;
  LR junc done; `recall_analyze.py` present on Sherlock.

## ⛔ AUDIT VERDICT: HOLD — default stays OFF (2026-07-11, workflow `wte43x5rc`)
`dev/MICROHOM_AUDIT_SYNTHESIS.md` (STATUS: COMPLETE). **Do NOT flip the `microhom_drift_margin`
default.** Guard ships as code, default OFF (byte-identical, safe — no action needed). Two blocking
reasons: (1) **gate unsatisfied** — 2 of 3 auditors stalled on API errors (discovery-loss at CP3,
detector-correctness at PLAN); incomplete ≠ all-clear. (2) **read-blind detector fault is
mechanically LIVE** (confirmed by synthesizer's independent code read): `_move_microhomology` is
GENOME-only, so the mh≥0.5 veto trigger and the read-evidence delta are INDEPENDENT — a real cryptic
the read distinguishes (delta>0) can still trip mh≥0.5 and be vetoed whenever delta < margin. byte-
identity auditor CLEARED inertness but explicitly declined to endorse the flip.

### ▶ FIX BUILT + AUDIT-V2 LAUNCHED (2026-07-12)
- **Fix committed `05664bc`** (branch worktree-agent-...): near-tie read-evidence cap
  `_effective_veto_margin(hold, eff, cap) = max(hold, min(eff, cap))`, param `drift_near_tie_cap`
  (default 0.0 = OFF = byte-identical) threaded through all 4 refine fns. Bounds the read-blind
  discovery-loss (move with delta_improve≥cap never drift-vetoed); hold_margin never capped. HONEST
  SCOPE: cap BOUNDS but does NOT close the fault (same score axis, no discriminating signal in (0,cap)).
- **Internal check PASSED** (per advisor — NOT efficacy validation): byte-identity at default
  (613 + 170 + 20 targeted tests green; refiner/validation tests would break if default changed),
  `_effective_veto_margin` arithmetic + hold-uncapped interaction unit-tested (20 in
  test_microhom_drift_guard.py). Default UNCHANGED — no flip.
- **Triple Opus-Max audit LAUNCHED: workflow `we02p54lj` (run wf_8994db9c-70a)**, script
  `dev/microhom_audit_v2.workflow.js`. 3 orthogonal Opus-Max auditors + synthesis, heavy incremental
  checkpointing (a stall costs one leg, not the verdict). Durable → `dev/MICROHOM_AUDIT_V2_*.md` +
  `dev/MICROHOM_AUDIT_V2_SYNTHESIS.md`. **RESUME:** on notification read the V2 synthesis — if a
  leg stalled, relaunch it pointed at its partial `.md` (or `Workflow({scriptPath, resumeFromRunId:
  "wf_8994db9c-70a"})` to reuse cached completed legs). The discovery-loss leg is load-bearing: it
  INDEPENDENTLY builds the real-cryptic-microhomology panel at varied delta_improve and decides
  whether the cap (at margin=3) bounds discovery loss OR the (0,cap) overlap is irreducible (→ needs a
  positional-distinctiveness signal, not aggregate delta). Only if all_clear → flip default at the
  audit's recommended (margin, threshold, cap).

### ✅ AUDIT-V4 DONE → HOLD + DETECTOR FIXES LANDED (2026-07-13)
`dev/MICROHOM_AUDIT_V4_SYNTHESIS.md`. 4 Opus-Max, 2/task; 7/9 agents stalled but the redundant design gave a
robust CONSENSUS. Results: cap helper CONSENSUS CLEAR; **A5 CONSENSUS HOLDING FAULT** (`_frac_match` N==N
phantom microhomology; shared by `_hp_run_across`); A8 (max-over-boundaries masks a genuine transition — A
fault / B collapses, agree on mechanics); discovery-loss RATE still unmeasured (both stalled) but the
delta>0∧mh≥0.5 case confirmed real + (0,cap) overlap IRREDUCIBLE → positional signal needed to CLOSE.
**FIXES IMPLEMENTED (Phase 5, this branch):** A5 — `_frac_match`/`_hp_run_across` count only ACGT (const
`_ACGT`); A8 — `_move_microhomology` combines moved boundaries by **min** not max. Byte-identical off (guard
still default OFF); 24 microhom tests + guard/refiner suites green (70 passed refiner-level).
**COMMITTED `d1fd08d`** (broad suite green: 645 passed, 33 skipped, 1 xfailed; 24 microhom tests). Detector
fixes + full audit trail durable on branch worktree-agent-... Cap = 05664bc. Guard STILL default OFF. NOT-YET-fixed / open: (i) discovery-loss rate quantification (stalls); (ii)
positional-distinctiveness signal to CLOSE (cap only bounds); (iii) min-k sensitivity floor; (iv) COMPASS
real-data (independent prereq).

### ▶ AUDIT-V4 REDUNDANT LAUNCHED (2026-07-13, workflow w1oi74gmz / run wf_8e872a33-eda)
Script `dev/microhom_audit_v4.workflow.js`. Per user: **4 Opus-Max auditors, TWO INDEPENDENT per task**
(discovery-loss A+B, detector-correctness A+B) for redundancy AND independent verdicts → **robust consensus**.
Each builds its OWN harness (own record `dev/MICROHOM_AUDIT_V4_<task>_<A|B>.md`, own scratch
`…/scratchpad/audit_v4/<task>_<variant>/`; do NOT couple the twins), retry-on-stall (2 attempts). Synthesis
computes CONSENSUS per task (A vs B agree?) + folds in V2 byte-identity (cleared). all_clear ONLY if consensus
CLEAR on every lens. (V3 2-agent version was STOPPED and superseded by this.)
**RESUME on notification:** read `dev/MICROHOM_AUDIT_V4_SYNTHESIS.md`. Check `discovery_loss_consensus` +
`detector_consensus`: if A&B AGREE CLEAR on both AND byte-identity cleared → all_clear (then decide the flip
at the recommended (margin,threshold,cap), still ALSO gated on COMPASS real-data). If A&B DISAGREE or either
HOLDs → not clear; reconcile per the stronger evidence. If a V4 leg STILL stalled both attempts, relaunch it
(`Workflow({scriptPath:".../microhom_audit_v4.workflow.js", resumeFromRunId:"wf_8e872a33-eda"})` reuses
completed legs' cached results) or measure inline. Key open Qs the numbers resolve: (a) is any (margin,cap)
defensible or is a positional-distinctiveness signal required to CLOSE (cap only bounds); (b) detector
A5(N==N)/A8(masking) real fault or not.

### ⛔ AUDIT-V2 RESULT (2026-07-12, workflow we02p54lj) → HOLD AGAIN
`dev/MICROHOM_AUDIT_V2_SYNTHESIS.md`. Same failure mode as V1: **byte-identity + structural-honesty leg
COMPLETED & CLEARED both rounds; the 2 load-bearing empirical legs (discovery-loss, detector-correctness)
STALLED on API errors AGAIN** (discovery-loss at CP1/orientation, detector at PLAN). Gate unsatisfied → HOLD.
- **CLEARED (byte-identity + structural honesty):** cap byte-identical off (incl. the guard_on_cap_off config
  that RUNS the replaced veto line: refined counts 815→308 yet output identical → `_effective_veto_margin(
  ...,cap=0)==eff_margin` exactly). Structural honesty HONEST: one-scalar-axis proven, cap only-lowers-
  never-raises (0 violations/200k), clamp identity. Flagged (not faults): P7 ergonomic sharp edge (any tiny
  positive cap ≈ near-max guard-weakening); the `test_restore_cat3_plus_2` pytest -x stop is an orthogonal
  missing-DATA fixture outside the cap commit range (guard suites 20/20 green).
- **★ DECISIVE STRUCTURAL FINDING:** the near-tie cap is a **BOUND, not a CLOSE**. Because delta_improve /
  eff_margin / cap share ONE score axis, a real cryptic and a fab drift with the SAME delta_improve are
  treated IDENTICALLY inside (0,cap). The cap protects delta_improve≥cap moves but CANNOT separate real
  discovery from fabrication in the near-tie band. To CLOSE the read-blind fault needs a POSITIONAL-
  DISTINCTIVENESS signal (distinctive in-window exon2 bases the incumbent soft-clip cannot absorb) — a
  deeper change than the cap.
- **UNMEASURED (the stalled legs):** the empirical discovery-loss RATE vs (margin, cap, microhomology_frac);
  detector A5 (`_frac_match` scores N==N as a match → false-pos veto on genome ambiguity runs) + A8
  (max-over-both-boundaries masking). Named for relaunch.
- **Root cause of stalls = API infrastructure, not the charter.** 2 rounds lost the load-bearing legs.
  Robust alternative: measure the discovery-loss rate INLINE (no workflow-agent stall risk), then adversarial-
  review the completed result with a SHORT Opus-Max agent.

### Remediation path (the audit's prescription — all keep default OFF)
1. **Complete the 2 stalled audits**, resumed from partial `dev/MICROHOM_AUDIT_{discovery-loss,
   detector-correctness}.md` (relaunch pointed at the partials, not restart).
2. **Add a read-evidence safeguard: near-tie veto gate** — only apply the veto when `|delta| < small`
   (a true near-tie), so a strongly-supported cryptic is NEVER vetoed regardless of genomic
   microhomology. Directly neutralizes the read-blind risk.
3. **Enlarge the tuning/validation panel** (broad R3/cryptic-3'SS at varied microhomology_frac).
4. **Revised operating point when the gate is satisfied: `microhom_drift_margin=3.0` NOT 8.0** (m=3
   already 96% fab-suppression, discovery-neutral on tuning panel; m=8 unproven on discovery loss),
   `microhom_threshold=0.5` unchanged. COMPASS real-data confirmation (§4b) is a ship prerequisite
   REGARDLESS.

## Open / in flight
2. **COMPASS real-data recall** (Sherlock bg `bm5g94ea5`, re-run of `bdkwin9ud`): `recall_analyze.py
   --min-split 3 --min-samples 2` — confirms the 0.5 threshold separates real SMA drift (fabrication)
   from real discovery on non-simulated data. Non-circular replacement for the retracted Snaptron
   verdicts (Snaptron is STAR-built/motif-filtered → confirm-only, cannot prove fabrication).
   **DATA ALL LOADED CLEAN** (durable): 4 SR samples SRR6376956/58/60/62 (360–407k junctions each,
   files OK), 13 LR paired samples. First run's stdout (sections 1–4b: truth set / feasibility gate /
   recall / reverse lead-check) was lost to ssh pipe block-buffering (inventory prints to stderr,
   analysis to stdout); re-run writes stdout to `/scratch/users/kevinroy/sma_recall/recall_result.txt`.

**Audit key finding so far (discovery-loss lens, the load-bearing one):** in a PURE microhomology
tie the guard is PROVABLY IRRELEVANT — the scorer soft-clips the drift-distance prefix, so a read
consistent with both incumbent and drift target scores delta=0, and the refiner's tie-break holds
the incumbent REGARDLESS of the guard (guard-OFF already fails to discover). Real discovery-loss
attributable to the guard therefore requires a cryptic site that scores STRICTLY BETTER (real
distinguishing evidence) yet still trips mh≥0.5 — the auditor is constructing exactly that case to
settle whether 0.5/8.0 is safe. Byte-identity auditor already confirmed (code read) inert at default.

## Resume (concrete)
- **Audit:** `ls dev/MICROHOM_AUDIT_SYNTHESIS.md` — when it exists, read it.
  - If **all-clear** (all 3 auditors find no fault) → flip the default: set
    `microhom_drift_margin: float = 8.0` in `refine_read_junctions` (and matching defaults in
    `_run_sequential/_run_parallel/refine_bam_junctions`) + surface it alongside `hp_drift_margin`
    in `correct_command.py`; re-run `pytest -m "not slow"`; ship alongside the HP-guard.
  - If **discovery-loss auditor finds the 0.5 threshold suppresses real discovery** → re-tune
    threshold/margin per its recommendation, re-run fab+R3 panels, re-audit before flipping.
- **COMPASS recall:** check `tasks/bdkwin9ud.output` (Sherlock `recall_analyze.py`). If it errored,
  re-run on Sherlock: `PYTHONPATH=/scratch/users/kevinroy/rectify_guard python
  /scratch/users/kevinroy/sma_recall/recall_analyze.py --min-split 3 --min-samples 2`. Interpret:
  does COMPASS 6-aligner short-read support separate the flagged drift from genuine discovery?

## Files
- `rectify/core/splice/junction_refiner.py` — guard impl (`_frac_match`, `_move_microhomology`, params)
- `tests/test_microhom_drift_guard.py` — 14 regression tests
- `dev/MICROHOMOLOGY_DRIFT_GUARD_DESIGN.md` — design + Phase 1–3 validation
- `dev/MICROHOM_AUDIT_*.md` — the 3 auditor lenses + synthesis (in flight)
- Sherlock `/scratch/users/kevinroy/sma_recall/` — COMPASS recall harness; `rectify_guard` = deployed guard code
