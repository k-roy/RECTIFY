# Coherent `del_open_delta` log-odds law in the re-placer — REJECTED (2026-07-06)

**Question (PI):** the `penalty_score` column we rejected (arm-C) is "not −logP, incoherent to
sum." The *coherent* quantity is `del_open_delta = λ·ln(rate_mean(hp)/rate_mean(1))` — a baseline-
anchored log-odds **gap-OPEN** delta (zero at hp=1), validated for exon-block indel attribution
(C1) but never wired into the junction *search*. Does it behave differently — help, or also hurt?

**Verdict: NO — it also over-shifts. The HP-drift guard is required regardless of the del-cost law.**
Keep the empirical table out of the re-placer search in **either** form (penalty_score OR del_open_delta).

## The two quantities (genuinely different — not an accidental re-wire of arm-C)
| | `penalty_score` (arm-C) | `del_open_delta` (arm-F) |
|---|---|---|
| formula | reciprocal-rate `c/rate_mean` | `λ·ln(rate_mean(hp)/rate_mean(1))`, zero at hp=1 |
| values (AT del, hp1→8) | 0.44→0.03 (falls) | 0→2.60 (rises), λ=1 |
| applied as | additive **per-base** del cost | one-time **gap-OPEN** discount (needs an affine-del DP) |

## Results (arm-F vs the ship arm-E = motif-blind + guard @m3.0)
- **Canonical drift (mix_fair_out), recovery:** ACC_A_D0 arm-E 0.993 vs **arm-F 0.312 (−0.680)**;
  ACC_A_D2 1.000 vs 0.825 (−0.175). arm-F never beats arm-E anywhere.
- **Over-shift (reads moved off the aligner placement):** arm-B 580 · arm-C 1085 · **arm-E 171** ·
  **arm-F 654** (~3.8× the surgical guard; net-harmful moves).
- **Distinct from arm-C:** at ACC_A_D0 arm-C got +0.69 vs arm-B, arm-F only +0.005 → the wiring is
  real, not an accidental arm-C. Yet it still over-shifts.
- **Discovery (mix_r3b_out):** arm-F merely *preserves* R3 (HP 0.290 vs arm-B/E 0.284, n.s.) while
  moving 495 reads vs arm-E's 294 — extra moves, zero benefit.
- **λ sweep:** λ=1 floors to 0.05 for all hp≥4 (delta 2.6–3.3 >> the re-placer's 0.5–1.0 del scale);
  the honest graded run **λ=0.2** over-shifts identically. **No λ separates beneficial from harmful
  moves** — the same "no sweet spot" as the arm-D hold_margin, restated on the del-cost axis.
- **Mechanism (verified):** ACC_A_D2 context `...CCAG|GT AAAAAA`; arm-F slides the acceptor 270→273,
  absorbing G,T,A into the intron and landing **inside the AAAAAA run** — the exact boundary-
  absorption over-shift the guard was built to stop. arm-E's `_hp_run_across` holds at 270.

## The deep reason
The degeneracy lives in the boundary **search** (a run-absorbing shift is a 0-error alignment), not
in the del-cost calibration. C1's coherence fix addresses indel *attribution within fixed boundaries*;
it cannot address the boundary *search*. So no del-cost term — coherent or not — fixes it; only a
prior (the guard) can. This closes the "is any cost term worth it in the search" thread: **no.**

## Provenance / recoverability
- Byte-identical off, verified 3 ways (affine reduces to linear exactly 0/3000; narrow suite 46 passed
  /17 skipped; broad splice/junction/penalty 507 passed). Investigated by a Completeness Agent
  (Opus), ~2.5M tokens, 66 tool calls.
- The exact arm-F wiring (hp_penalty.py +139 affine-del DP, junction_scoring.py +46, junction_refiner.py
  +50) is **reverted from production** and preserved as `dev/arm_f_del_open_delta.patch` (re-appliable);
  `scripts/benchmark/noncanon_sim/_make_arm_f.py` is the runner (needs the patch applied).
- **Caveats the agent flagged:** (1) *why* arm-F fails to fix D0 while arm-C fixes it is inferred (one-
  time vs per-base cheapening), only the harmful moves were mechanistically verified; (2) λ was ported
  from a gap_open=−4 DP into a ~0.5-scale DP — mitigated by the graded λ=0.2 giving the same verdict,
  but a fully scale-matched re-derivation of λ in the re-placer's own units was not done. Neither
  affects the verdict (arm-E dominance + the D2 over-shift stand on the verified mechanism).

## FAITHFUL RE-RUN — proof gap CLOSED, REJECT confirmed definitively (2026-07-07)
The triple adversarial pass (Workflow whf5mx4qj) held the verdict but flagged a real fault: the arm-F affine-DP
had a STAIRCASE (del→match(cost0)→del re-collects the one-time gap-OPEN discount per in-run base), so arm-F was
never a faithful one-time law in HP runs, and the RECOVERY axis was never faithfully tested. A follow-up agent
FIXED it and re-ran:
- STAIRCASE FIX (hp_penalty.py::_score_hp_affine_del, Option-B run-latch 3-state DP H=avail/Hs=spent/D): discount
  is ONE-TIME PER RUN; a genuine run boundary re-arms it. PROVEN (scratchpad/staircase_proof.py): k=2..5 absorb →
  BUGGY=k·open, FIXED=open+(k-1)·extend (brake); random 4000/4000 == once-per-run oracle; byte-identity off.
- SCALE-MATCHED λ=0.1 (del units, del_hp=0.5): open-cost graded 0.41→0.18 across hp4..12, monotone, NEVER floored.
- RESULT (faithful, λ=0.1), INDEPENDENTLY RE-SCORED here on the agent's bams:
  * arm_Ffg (faithful + guard) vs arm_E: R3 mix_r3b_out b=0 c=0 p=1.0 (IDENTICAL); fair mix_fair_out all cells
    non-sig (p>=0.25), mixed sign, overall 0.911==0.911 — a WASH. del_open adds NOTHING on top of the guard.
  * arm_Ff (faithful, STANDALONE no guard) vs arm_B: over-shifts — ACC_A_D2 -0.562 p=0 (225 lost, del_open's OWN
    C1 designed domain!), BOT_A_D1 -0.405, ACC_A_D1 -0.110. A faithful one-time law standalone is WORSE than arm-B.
  STRUCTURAL WHY: del_open lowers a gap-open by <=0.32 cost units; the guard margin is 3.0 in the SAME units, so
  del_open is 3-10% of the margin and cannot flip a guard veto. Raise λ→floors; lower margin→over-shift returns.
=> REJECT CONFIRMED DEFINITIVELY across all 3 axes (over-shift closed analytically by arbiter 2's del-cost-
   invariance theorem; recovery + discovery closed empirically here). The guard ALONE is the answer; the empirical
   table stays OUT of the re-placer search in every form (penalty_score column AND coherent del_open_delta).
Faithful wiring preserved: dev/arm_f_del_open_delta_faithful.patch (staircase-fixed; supersedes the buggy
arm_f_del_open_delta.patch). Detailed agent handoff: dev/HANDOFF_DEL_OPEN_STAIRCASE.md.
