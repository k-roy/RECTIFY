# HANDOFF — del_open_delta arm-F staircase fix + recovery-axis panel

Branch: `worktree-agent-a25a2c1e784ad37dc`. Python: `/Users/kevinroy/miniconda3/bin/python`.

## Done
- Confirmed baseline green (`pytest tests/test_junction_refiner.py tests/test_hp_drift_guard.py` = 46 passed, 17 skip).
- Applied `dev/arm_f_del_open_delta.patch` (arm-F wiring) — tests still green (lam=0 path untouched).
- **Fixed the staircase in `rectify/core/splice/hp_penalty.py::_score_hp_affine_del`** (Option B, run-latch
  three-state DP: H=avail / Hs=spent / D). Discount is ONE-TIME PER HP RUN; a genuine run boundary re-arms it.
- **Proof** (`scratchpad/staircase_proof.py`): headline k=2..5 absorbing in-run bases →
  BUGGY = k·open, FIXED = open+(k-1)·extend (BRAKE OK). Random 4000/4000: FIXED==once-per-run oracle,
  FIXED>=BUGGY always, no-discount FIXED==BUGGY (byte-identity fence).
- **λ = 0.1** chosen (del units, del_hp=0.5): open-cost graded 0.41→0.18 across hp4..12, monotone, never floored.
  (λ=0.2 floors at hp>=6; λ=1.0 = the buggy panels' value floors almost everywhere.) Rate ratios monotone (data).

## FAIR PANEL RESULT (mix_fair_out, λ=0.1) — canonical-drift, del_open's designed domain
Refined counts: arm_Ff (standalone) 1165/5600; arm_Ffg (guard) 169/5600 ≈ arm_E 171.
Recovery (overall R0flank): B=0.839  E_m3p0=0.911  Ff=0.760  Ffg=0.911.
- arm_Ff standalone SIGNIFICANTLY WORSE than arm-B (over-shifts): ACC_A_D2 −0.562 p=0 (225 lost/0 gained),
  BOT_A_D1 −0.405 p=0, ACC_A_D1 −0.110 p=0. Worst in ACC_A_D2 = del_open's OWN designed C1 domain.
- arm_Ffg vs arm_E (McNemar, paired): NO cell significant (all p>=0.25), mixed sign — ACC_A_D0 +3 (p=0.25),
  BOT_A_D1 +2 (p=0.5), ACC_A_D2 −1, ACC_A_D8 −2. Overall 0.911 == 0.911. A WASH.
STRUCTURAL WHY: del_open lowers a gap-open by <=0.32 cost units (open 0.18-0.41 vs del_hp 0.5); guard margin
is 3.0 in the SAME units → del_open is 3-10% of the margin, cannot flip a guard veto. Raise λ → floors;
lower margin → over-shift returns (arm_Ff proves it). No operating point helps. => RECOVERY axis REJECT.

## R3B PANEL RESULT (mix_r3b_out, λ=0.1) — R3 discovery
Refined: arm_E_m3p0 292/4800; arm_Ff (standalone) 694/4800; arm_Ffg (guard) 291/4800 (guard clamps to ≈E).
R3 recovery: B=0.608  E=0.608  Ff=0.604  Ffg=0.608 (HP cell all 0.284; arm_Ff over-shifts a couple → 0.278).
McNemar E_m3p0 vs Ffg on R3: b=0 c=0 p=1.0 on BOTH HP and plain cells — arm_Ffg IDENTICAL to arm_E.
R0/R1 HP-drift cells: E=0.968 = Ffg=0.968 > B=0.954; arm_Ff 0.952 (again slightly worse than B).
Guarded refined 292(E) vs 291(Ffg): del_open perturbs exactly ONE guarded decision of 4800.

## VERDICT — REJECT CONFIRMED DEFINITIVELY (all 3 axes)
- D0/D2 recovery: Ffg ≈ E, no significant cell (all p>=0.25), mixed sign, WORSE at ACC_A_D2 (designed domain).
- Over-shift: Ff standalone over-shifts (−0.562 @ACC_A_D2 p=0); Ffg = E.
- R3 discovery: Ffg == E exactly (McNemar b=c=0). No Pareto win anywhere. Coherent one-time law does NOT belong.
Tests: 146 passed / 17 skip / 1 xfail (junction_refiner + hp_drift_guard + junction_scoring_parallel + splice_junction).

## Working tree (do NOT commit; task-sanctioned to leave applied)
- arm_F patch APPLIED + Option-B run-latch fix in hp_penalty.py + runners (_make_arm_ff/_make_r3b_arms/_score_arms_multi)
  + generated BAMs (arm_Ff*/arm_Ffg*/arm_E_m3p0 in mix_fair_out & mix_r3b_out). To restore pristine:
  `git checkout rectify/core/splice/hp_penalty.py rectify/core/splice/junction_refiner.py rectify/core/splice/junction_scoring.py`
  then remove the untracked runners/BAMs if desired.

## Single least-sure caveat
Verdict is stated at the guard operating point (hp_drift_margin=3.0). del_open (<=0.32) is swamped there by
design; a much SMALLER margin would readmit del_open's influence — but arm_Ff (margin 0) proves that path just
restores over-shift. Did NOT sweep the full margin×λ plane (M1 refines are expensive; mechanism covers it).

## Files
- FIX: `rectify/core/splice/hp_penalty.py` (_score_hp_affine_del run-latch).
- Runners: `scripts/benchmark/noncanon_sim/_make_arm_ff.py`, `_make_r3b_arms.py`, `_score_arms_multi.py`.
- Proof: `scratchpad/staircase_proof.py`. Patch: `dev/arm_f_del_open_delta.patch` (currently APPLIED).
- Do NOT commit. Working tree left with patch applied + fix + runners (note in final report).
