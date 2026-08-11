# INSCOST AUDIT — LENS: IMPACT-MEASUREMENT VALIDITY

Adversarial auditor. Goal: find the mis-measurement that makes "full-run ins_cost is sounder / safe to switch" WRONG.
Lens = are the investigator's IMPACT numbers trustworthy? Re-derive witness; check call-change-rate methodology
(per-junction vs per-read; right panels; direction); could change be net-HARMFUL (dropping REAL non-canonical
recoveries)?; independently re-score guard-interaction (R3 discovery + drift-fix preserved?).

## PLAN
- [ ] Read INSCOST_INVESTIGATION.md (durable record)
- [ ] Read junction_scoring.py _score_junction / _score_hp_anchored / ins_cost / _hp_run_length
- [ ] Re-derive fabrication witness (all-A(12)) independently, flag OFF vs ON
- [ ] Locate the call-change-rate harness / panels (mix_fair_out, mix_r3b_out) — how measured?
- [ ] Check: per-junction or per-read? right denominator? direction attribution?
- [ ] Check net-harm: are the 0->4 moves REAL recoveries lost, not just fabrications removed?
- [ ] Re-score guard-interaction claim independently
- [ ] Verdict

## CHECKPOINTS
- CK1: Read INSCOST_INVESTIGATION.md + junction_scoring.py (_score_junction/_score_hp_anchored) + hp_penalty.py (ins_cost,_hp_run_length). Flag = env RECTIFY_FULL_RUN_INS=="1", default OFF. full_ins_costs precomputed on full rescue window, sliced [k:] and [:k][::-1]. Confirmed code matches report.
- CK2: Measurement harness = _ins_compare.py. CALL-CHANGE is PER-READ (iterates truth.items(), denom=considered=reads w/ has_true_junction in rung set). Direction via _canonical_tier + recovery-vs-truth. Arms built by _make_ins_arms.py (motif_blind + penalty_table, guard=0 or 3.0), flag from env per subprocess.
- CK2a: STRUCTURAL FLAW forming: anti-fabrication headline (fair 35:6, +0.005) is NO-GUARD config. On fair(drift) panel truth==canonical BY CONSTRUCTION, so "became_MORE_canonical" == "toward truth" trivially — cannot distinguish fabrication-removal from canonical-bias. Production uses guard=3.0 where SAME flip is net-NEG (fair 0.899->0.898, all 4 changes recov_LOST). Need to verify by re-score.
- CK3: REPRODUCED investigator numbers EXACTLY via _ins_compare.py.
  FAIR no-guard: 42/5600=0.75%, 35 reads 4->0, 6 reads 0->4, R0flank +0.005. (panel is 100% R0flank => truth==canonical => "MORE_canonical"=="toward truth" tautology.)
  FAIR guard=3.0 (PRODUCTION): 4/5600=0.07%, ALL 4 changes 0->4 recov_LOST, R0flank -0.001. NET NEGATIVE.
  R3B no-guard: 10/4800=0.21%, recov_GAINED=5 recov_LOST=4, R3 +0.001 R0 -0.001.
  R3B guard=3.0 (PRODUCTION): 5/4800=0.10%, recov_LOST=3 recov_GAINED=1, R3 -0.001 R0 -0.001 WT -0.001.
  => The "+0.005 / 35:6 anti-fabrication" headline is NO-GUARD only. Production (guard) config: flip is net-NEG on BOTH panels.
- CK4: READ-LEVEL changed-read dump (changed_dump.py). SMOKING GUN for net-harm direction:
  R3B PRODUCTION (guard=3.0): 5 changed = 3 LOST (R0 true_tier0 canonical dropped 0->4 FABRICATED; WT true0 recovered->lost 0->4 FABRICATED; R3 true_tier8 REAL non-canonical recovery DROPPED 4->8) + 1 GAINED (R1) + 1 neutral. => full-run INTRODUCES fabrications (0->4) and DROPS a real R3 non-canonical discovery; net -2.
  R3B no-guard: LOST by rung {R0:1, WT:2, R3:1}; GAINED {R1:1, R3:2, WT:2}. R3 net +1 (WASH). Full-run BOTH removes and introduces fabrications ~balanced on the panel that actually has non-canonical truth.
  => The clean "4->0 anti-fabrication" direction ONLY appears on fair (100% R0flank canonical-truth). On the discovery panel full-run is not directionally anti-fabrication.
- CK4a: CONFIRMED fair panel = 5600 reads ALL R0flank (rung_rec printed only R0flank; R0/R1/R3/WT n=0). So the 35:6 headline is measured on a panel with ZERO non-canonical truth -> "MORE_canonical"=="toward truth" tautologically, and it is structurally INCAPABLE of detecting dropped-real-non-canonical-recovery. CONFOUNDED.
- CK5: Fabrication WITNESS reproduced EXACTLY (investigator witness.py AND my own in-process toggle): all-A(12) A-free genome OFF=1.7604 ON=8.2584, 4.69x. The run-splitting gaming + its removal are REAL and correctly measured. Witness is NOT the flaw.

## VERDICT (impact-measurement-validity lens)
FAULT FOUND. The witness + call-change RATES are correctly measured, but the DIRECTION/RECOVERY evidence
cited to support "full-run is sounder / anti-fabrication / safe to switch" is MIS-MEASURED in two compounding ways
that inflate "sounder", and REVERSES sign under the deployment configuration.

STRONGEST CHALLENGE (holds):
(1) CONFOUNDED HEADLINE PANEL. The "35:6 toward canonical / +0.005 recovery / anti-fabrication" signal is measured
   ENTIRELY on mix_fair_out, which is 100% R0flank == canonical-truth. On a canonical-only panel
   "became_MORE_canonical" is DEFINITIONALLY "moved toward truth", so the metric cannot distinguish "removed a
   fabrication" from "biased toward canonical", and it is structurally INCAPABLE of detecting the failure mode
   that matters for a DISCOVERY tool: dropping a REAL non-canonical recovery (there is no non-canonical truth on
   this panel). The corroboration is an artifact of panel construction.
(2) OFF-DEPLOYMENT CONFIG. The net-positive numbers are the NO-GUARD arms. Production ships hp_drift_margin=3.0
   (the investigator's own guard section calls the guard the DOMINANT, deployed drift-fix; guard+full-run 0.898
   ~= guard+per-cut 0.899). In the guard-on (production) config the SAME flip is net-NEGATIVE on BOTH panels:
   fair R0flank 0.899->0.898 (ALL 4 changed reads are recov_LOST, 0 gained); r3b R3 0.604->0.604 with changed
   reads 3 LOST : 1 GAINED.
(3) ON THE ONLY NON-CANONICAL-TRUTH PANEL (r3b) FULL-RUN INTRODUCES FABRICATIONS AND DROPS A REAL DISCOVERY.
   Production config, read-level: of 5 changed reads, 3 are recoveries LOST — a canonical R0 truth demoted
   0->4 (fabricated non-canonical), a canonical WT truth demoted 0->4 (fabricated), and a REAL R3 non-canonical
   truth recovery DROPPED (4->8) — vs 1 gained. The clean 4->0 anti-fabrication direction does NOT appear here.
   No-guard r3b: R3 net +1 (wash), full-run both removes AND introduces fabrications ~balanced.
CONSEQUENCE: The recovery/anti-fabrication CORROBORATION for "sounder" is not trustworthy (confounded + sign-
   reverses in production), and "safe to switch the DEFAULT" is contradicted by the deployment-config measurement
   (mildly net-negative + introduces fabrications on discovery). Either the effect is real (then production is
   net-negative, so don't flip the default on recovery grounds) or it is noise at +-1-2 reads/1600 (then the
   +0.005/35:6 headline is ALSO noise and corroborates nothing). Either way the empirical case for "sounder /
   safe to flip default" collapses. NOTE: this lens does NOT refute the PRINCIPLED model argument or the DP-unlock
   (structural) — those ride on other lenses; but the report's stated EMPIRICAL corroboration is inflated.
SECONDARY: call-change is measured PER-READ; the discovery unit is the JUNCTION. Per-read weights by coverage,
   under-representing rare cryptic junctions (the exact discoveries the tool exists for). Panels are rung-balanced
   sim reads so it is not the primary flaw, but the "0.21%/0.75%" rates are per-read, not per-junction.
- CK6: ADVISOR BLOCKER resolved -> claim #2 CORRECTED. Verified production call site correct_command.py:746 does NOT pass hp_drift_margin -> default 0.0 => GUARD IS OPT-IN, not shipped. So no-guard arms = deployment config (guard axis); on no-guard full-run is net-POSITIVE fair +0.005, wash r3b. My "reverses/net-negative in production" claim is FALSIFIED. Retract it.
- CK6a: BUT two NEW verified config mismatches keep the fault alive, independent of the guard:
  * motif_blind default=False (junction_refiner.py:507); correct_command.py does NOT pass motif_blind => PRODUCTION = motif_blind=False. ALL impact arms built motif_blind=True (canonical-motif prior DISABLED). The prior in production already suppresses non-canonical fabrications (tier_beats_alt, junction_refiner.py:655) — same subsumption logic the investigator applied to the guard, but the motif prior IS the production default. => the +0.005/35:6 anti-fabrication is measured with the prior OFF = UPPER BOUND, expected to shrink under motif_blind=False.
  * junction_penalty_table default='' (split_command.py:1245) => flag is a COMPLETE NO-OP in default pipeline; only acts under opt-in --junction-penalty-table (and even then motif_blind=False).
- CK6b: Testing motif_blind=False empirically on fair panel to confirm subsumption (decisive check).
- CK7: SUBSUMPTION strengthener (advisor), from reproduced per-rung recovery, fair R0flank:
    off_ng(per-cut,noguard)=0.777  on_ng(full-run,noguard)=0.781  [full-run marginal +0.004]
    off_g3(per-cut,guard) =0.899  on_g3(full-run,guard) =0.898  [full-run marginal -0.001]
    guard effect: per-cut 0.777->0.899=+0.122 ; full-run 0.781->0.898=+0.117.
  => The guard prior delivers ~+0.12 under BOTH flag states; full-run's marginal benefit (~+0.004) is subsumed to ~0/neg once ANY canonical prior is on. This is the mechanism that predicts the motif_blind=False (production) benefit shrinks too.

## REVISED VERDICT (supersedes the earlier VERDICT block; reconciles per advisor)
RETRACTION: my earlier claim #2 ("production ships guard=3.0, flip is net-negative") is FALSE. Verified the
guard is OPT-IN (default 0.0; correct_command.py:746 does not pass it). So for the guard axis the deployment
config == no-guard, where full-run is net-POSITIVE on canonical truth (+0.005/35:6) and neutral on discovery.
The earlier "empirical case COLLAPSES" wording is OVERSTATED — corrected to "config-scoped and one-sided."

FAULT THAT SURVIVES (impact-measurement validity):
The recovery/anti-fabrication evidence cited to CORROBORATE "sounder" is measured off the SHIPPED configuration
and is one-sided, so it OVER-FRAMES a narrow benefit as general "net-positive anti-fabrication":
 (A) CONFIG MISMATCH (motif prior). All impact arms use motif_blind=True (canonical-motif prior OFF). Production
    ships motif_blind=False (correct_command.py passes no motif_blind; default False). The motif prior already
    suppresses non-canonical fabrications (tier_beats_alt) — the SAME subsumption the investigator documented for
    the guard (CK7: any canonical prior flattens full-run's +0.004 to ~0/neg). So the 35:6/+0.005 is an UPPER
    BOUND measured with the production prior disabled; it is expected to shrink under motif_blind=False. [empirical
    test in flight: mbF_build.log]
 (B) CONFOUNDED / ONE-SIDED PANEL. The 35:6 headline is 100% R0flank (canonical truth) where "MORE_canonical" ==
    "toward truth" tautologically. It correctly measures fabrication-REMOVAL but CANNOT measure the discovery COST.
    The discovery panel (r3b, the only non-canonical truth) shows full-run is a WASH on real R3 recovery (0.589
    no-guard, net +1 read) AND itself INTRODUCES fabrications (read-level: canonical truths demoted 0->4). So the
    net anti-fabrication direction claimed does NOT generalize to the discovery use case RECTIFY exists for.
 (C) BLAST RADIUS. junction_penalty_table default='' => the flag is a COMPLETE no-op in the default pipeline;
    it only acts under opt-in --junction-penalty-table (and then motif_blind=False). The impact — positive or
    negative — is confined to that opt-in path, not the shipped default.
 (D) SECONDARY: call-change is PER-READ; the discovery unit is the JUNCTION. Per-read weights by coverage,
    under-representing rare cryptic junctions.

WHAT SURVIVES: the fabrication WITNESS (gaming removal 1.76->8.26) is REAL, reproduced independently in-process.
The DP-unlock (structural, cut-independence) and the PRINCIPLED model argument are NOT touched by this lens. So
through the impact-validity lens the CORE verdict (full-run = sounder model + DP unlock) SURVIVES. The defensible
fault is narrower: the recovery/anti-fabrication CORROBORATION is config-scoped (motif prior off, guard off, sim
only) and one-sided (canonical-truth panel), and the report over-frames it as general net-positive anti-fabrication
— it does not establish a discovery-use benefit, and by the investigator's own subsumption logic it will shrink
under the shipped motif_blind=False config.

## CK8 — DECISIVE: motif_blind=False (production prior) test FALSIFIES my subsumption fault (A)
Built fair arms with motif_blind=False (production setting) + penalty_table + no guard, flag OFF vs ON:
  CALL-CHANGE 42/5600 = 0.75% ; 35 reads 4->0, 6 reads 0->4 ; R0flank 0.778 -> 0.783 = +0.005.
  => IDENTICAL to motif_blind=True (+0.005, 35:6). The canonical-motif PRIOR does NOT subsume full-run's benefit
  (unlike the GUARD, which does). Full-run's ins-cost mechanism is ORTHOGONAL to the motif tie-break prior, so the
  anti-fabrication benefit TRANSFERS to the production motif_blind=False config. My fault (A) is EMPIRICALLY FALSE.

## FINAL RECONCILED VERDICT (impact-measurement-validity lens)
After adversarial probing, the investigator's impact numbers are TRUSTWORTHY. Every headline reproduced EXACTLY
(witness 1.76->8.26 in-process; fair 42/5600 +0.005 35:6; r3b 10/4800; guard configs). My two strongest
quantitative attacks were EMPIRICALLY FALSIFIED by checking the shipped config:
  - claim #2 (guard ships on -> net-negative): FALSE. Guard is OPT-IN (default 0.0), so no-guard = deployment;
    full-run is net-POSITIVE there.
  - fault A (motif prior subsumes benefit): FALSE (CK8). Benefit transfers to production motif_blind=False.
NET-HARM CHECK (task's core question): r3b R3 real non-canonical discovery 0.589 no-guard, delta +0.001
  (read-level R3 gained 2 / lost 1). Discovery is PRESERVED, NOT net-harmful. The lone "R3 lost" read has
  tier/recovery disagreement (normalization slide) = rebuttable noise at n=1600.
RESIDUAL LEGITIMATE CAVEATS (do not overturn "sounder"; mostly disclosed by investigator):
  (B) The "anti-fabrication 35:6" is canonical-truth-only (R0flank); it demonstrates fabrication-REMOVAL, not a
      discovery gain. On non-canonical truth (r3b) full-run is a WASH and itself introduces a few 0->4 fabrications
      (balanced by removals). So "net-positive anti-fabrication" should read "removes canonical-truth fabrications
      without harming discovery" — the investigator's own hedge ("recovery small/corroborating, not primary";
      flagged the 0->4 counter-signal) is appropriate.
  (C) BLAST RADIUS: junction_penalty_table default='' => the flag is a COMPLETE no-op in the DEFAULT pipeline;
      impact (pos or neg) is confined to the opt-in --junction-penalty-table path.
  (D) call-change is PER-READ not per-junction (discovery unit); panels rung-balanced so minor.
  (E) sim-only; human ONT transfer UNEXERCISED (investigator flagged as revalidation #3).
CONCLUSION: fault_found = NO mis-measurement that inflates "sounder" or makes the switch unsafe on impact grounds.
The impact evidence SURVIVES this lens. The verdict (full-run sounder + DP unlock) survives. Only framing/scope
caveats remain (B/C/D/E), which the report largely discloses.

## CK9 — advisor caught asymmetry: benefit validated prior-ON (fair), cost only ever prior-OFF (r3b mb=True).
Building the DECISIVE config: r3b, motif_blind=False (production), no-guard, flag off vs on -> R3 recovery delta.
This is opt-in-penalty-table production on a NON-CANONICAL-truth panel where the motif prior FIGHTS the truth.

## CK10 — DECISIVE DISCOVERY-COST TEST (r3b, motif_blind=False=production, no-guard, flag OFF vs ON)
CALL-CHANGE 7/4800 = 0.15%. Per-rung recovery delta:
  R3 (cryptic non-canonical discovery) 0.258 -> 0.258 = +0.000  <-- DISCOVERY EXACTLY PRESERVED
  R0 -0.001, R1 +0.001, WT +0.000. Direction balanced: 3 LOST / 3 GAINED / 1 neutral.
NB R3 baseline 0.258 (mb=False) < 0.589 (mb=True) because the canonical prior itself biases away from cryptic
  recovery — but the FLAG delta on top is 0.000. Full-run neither helps nor harms cryptic discovery in production.
=> Advisor decision rule satisfied: R3 delta neutral (0.000) on the production prior-ON non-canonical panel.
  The "not net-harmful to discovery" half of the verdict is AIRTIGHT. No read-level dig needed (aggregate flat).

## SIGN-OFF (impact-measurement-validity lens)
fault_found = FALSE. No mis-measurement inflates "sounder"; no net-harm to discovery.
- All investigator headline numbers reproduced EXACTLY (witness in-process; fair/r3b call-change & recovery).
- Both my quantitative attacks empirically FALSIFIED: guard is opt-in (no-guard=deployment); motif prior does NOT
  subsume the benefit (fair mb=False = +0.005/35:6, identical).
- Decisive production-config discovery test: R3 delta 0.000. Discovery preserved. Change is NOT net-harmful.
- Witness (gaming removal) REAL & independent; DP-unlock & principled argument untouched by this lens.
VERDICT SURVIVES. Residual = framing/scope caveats only (B: 35:6 is fabrication-REMOVAL on canonical truth, not a
  discovery gain — discovery is flat; C: no-op in DEFAULT pipeline, opt-in --junction-penalty-table only; D:
  per-read not per-junction; E: sim-only, human transfer unexercised — all largely disclosed by the investigator).
