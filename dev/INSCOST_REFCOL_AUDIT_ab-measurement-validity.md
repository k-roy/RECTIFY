# INSCOST REF-COLUMN AUDIT — A/B MEASUREMENT VALIDITY

**Auditor lens:** Is the A/B measurement trustworthy? Re-derive the headline claim
(ref-column keeps the 35 anti-fabrication wins AND drops the 4-6 canonical-demotion
losses) INDEPENDENTLY by toggling flags and re-scoring. Hunt for the mismeasurement.

**Adversarial hypotheses to test:**
1. Is the tier count per-JUNCTION or per-READ? (claimed unit matters)
2. Are tier 4->0 / 0->4 counts computed IDENTICALLY across all 3 arms?
3. Could ref-column's "no losses" be a SENSITIVITY artifact — it makes FEWER moves
   overall (less discovery), not that it's more correct?
4. Net discovery: does ref-column DROP real non-canonical recoveries the other arms keep?

**Constraints:** read-only work; toggle flags in-process. numba OFF (must stay off).
Python: /Users/kevinroy/miniconda3/bin/python. Panels:
scripts/benchmark/noncanon_sim/{mix_fair_out,mix_r3b_out}. M1 KILLS heavy refines ->
lean, one arm at a time.

## PLAN
- [ ] Read the prior investigation docs (INSCOST_INVESTIGATION, SYNTHESIS, AUDIT model-correctness, REFCOL_BUILD)
- [ ] Locate the flag(s): _USE_FULL_RUN_INS + the ref-column flag in junction_scoring.py
- [ ] Read _score_junction / _score_hp_anchored / _precompute_del_costs to understand tier assignment + how each arm changes ins_cost
- [ ] Understand how the "35 wins / 4-6 losses" panel metric is actually computed (the harness)
- [ ] Re-derive counts independently for all 3 arms with IDENTICAL accounting
- [ ] Check the unit (per-junction vs per-read) and whether counts are symmetric
- [ ] Test the sensitivity hypothesis: total moves per arm, net discovery
- [ ] Verdict

## CHECKPOINTS
- (start) file created, plan written
- Read BUILD.md + SYNTHESIS.md. KEY: the builder ALREADY reports the audit hypothesis
  FALSIFIED. Builder's fair no-guard 3-way (BUILD.md:186-205):
    percut->fullrun: CALL-CHANGE 42/5600; tier 4->0=35, 0->4=6, 4->4=1; R0flank +0.005
    percut->refcol : CALL-CHANGE 694/5600; tier 4->0=358, 0->2=20, 0->4=196, 2->0=2, 4->4=118; R0flank +0.011
  => refcol does NOT keep 35 & drop losses; it BALLOONS both directions (358 wins, 196 losses).
  My lens: is THIS measurement itself trustworthy? Must re-derive independently + check the
  4 mismeasurement hypotheses. The builder's own numbers already refute "no losses" — so the
  SENSITIVITY hypothesis is INVERTED (refcol makes 16x MORE moves, not fewer). Verify the
  adjudicator (_ins_compare.py) computes tier-shift identically per-arm, per-junction vs per-read.

## ADJUDICATOR READ (_ins_compare.py) — measurement-mechanics findings
- UNIT = per-READ, not per-junction (line 82 `for rid, f in truth.items()` keyed read_id).
  Task framing "is it per-JUNCTION" -> ANSWER: per-READ. read_junctions takes FIRST op==3
  (line 39-42), so per-read == per-called-junction, BUT many reads can share ONE junction
  coordinate; NO dedup to unique junctions. A shift on a high-multiplicity coordinate inflates
  the raw count vs a per-unique-junction metric.
- TIER-SHIFT symmetric: same _canonical_tier on off & on.
- tier_shift[(toff,ton)] fires ONLY when BOTH arms have an N-op junction; if either arm drops
  the junction (None) -> tier()=None -> cdir=tier_na, EXCLUDED from tier_shift. So a per-cut
  tier0 -> refcol drops-junction case is NOT counted as 0->4. => demotion count can UNDERCOUNT
  if the aggressive arm removes junctions rather than moving them. Must check.
- BASELINE: both 3-way legs use arm_mbF_percut_ng as `off`. Confirm identical baseline.
- ARMS ON DISK: fair has all 6. r3b MISSING fullrun_ng + all g3 (only percut_ng, refcol_ng).
  => r3b 3-way INCOMPLETE; builder headline is FAIR no-guard only.
- CHECKPOINT: adjudicator mechanics understood; re-running fair no-guard 3-way next.

## INDEPENDENT RE-DERIVATION (fair no-guard) — EXACT MATCH to builder
baseline arm_mbF_percut_ng.bam md5=06968d05dd6f86ccf8bf58ede5a85b57 (SAME file both legs -> OK).
- percut->fullrun: CALL-CHANGE 42/5600; tier 4->0=35, 0->4=6, 4->4=1; R0flank d=+0.005. MATCHES.
- percut->refcol : CALL-CHANGE 694/5600; tier 0->2=20, 0->4=196, 2->0=2, 4->0=358, 4->4=118;
  R0flank d=+0.011. MATCHES builder exactly.
- tier_shift buckets SUM to total changed for BOTH arms (42 and 694) => NO junctions dropped to
  None among changed reads; the None-undercount concern does NOT bite here. tier machinery is
  self-consistent and symmetric.

## *** THE MISMEASUREMENT — found (net-recovery decomposition) ***
Recovery decomposed over ALL 5600 has_true_junction reads (truth-anchored, bypasses tier machinery):
- percut->fullrun: base_rec=4357 arm_rec=4384 NET=+27 | GAINED=32 LOST=5
- percut->refcol : base_rec=4357 arm_rec=4416 NET=+59 | GAINED=255 LOST=196
CRITICAL: in the fair {R0flank,INTRONFREE} rung set, the TRUE junction is CANONICAL for EVERY read
(gained_true_noncanon=0 AND lost_true_noncanon=0 for BOTH arms). => the "4->0 anti-fab wins" and
"0->4 canon-demotion losses" are NOT measuring non-canonical DISCOVERY at all; they measure whether
the arm lands on the true (canonical) coordinate. On the truth axis refcol RECOVERS MORE
(+59 vs +27), i.e. refcol is STRICTLY BETTER at truth-recovery on this panel, DESPITE its 196 "0->4"
tier-shift losses.

=> MISMEASUREMENT #1 (headline reframes a WIN as a LOSS): the audit's recommendation is graded on
tier-SHIFT counts (0->4 "losses"), but the tier-shift metric is NOT truth-anchored. A read where
per-cut called a tier0 coordinate that was WRONG (not the true junction) and refcol moves it to a
tier4 coordinate that is ALSO wrong shows up as a "0->4 loss" — yet no true recovery was lost. The
proper metric (net truth recovery) says refcol WINS (+59 > +27). The "196 losses balloon" headline
is an artifact of scoring on the un-anchored tier-shift axis. The audit's own predicted metric
(keep 35 wins / drop 6 losses) is the WRONG YARDSTICK: it counts tier movements, not correct calls.

=> SENSITIVITY hypothesis (my H3) RESOLVED: refcol is MORE sensitive (694 vs 42 moves), NOT less;
"no losses" was never the refcol result. But higher sensitivity here is NET-POSITIVE on truth
(+59). Refcol does not lose real recoveries the other arms keep on THIS rung set (lost=196 but
gained=255; net +59; and 0 of the lost are non-canonical truths).
- CHECKPOINT: mismeasurement #1 (tier-shift not truth-anchored) persisted. Checking R3 discovery next.

## *** NET DISCOVERY (r3b panel) — the axis that MATTERS: refcol REGRESSES real non-canon ***
r3b has the real non-canonical discovery rung R3 (1600 reads, true junction non-canonical).
percut->refcol per-rung recovery (independently reproduced):
  R0  n=1200 off=0.979 on=0.979 d=+0.000
  R1  n=1200 off=0.878 on=0.887 d=+0.008
  R3  n=1600 off=0.258 on=0.237 d=-0.021   <-- REAL non-canonical DISCOVERY REGRESSION
  WT  n= 800 off=0.781 on=0.811 d=+0.030
R3 absolute (truth-anchored): percut_rec=412 refcol_rec=379 NET=-33.
  refcol GAINED=9 true-R3 recoveries, LOST=42 true-R3 recoveries.
  => refcol DROPS 42 real non-canonical junction recoveries (gains only 9) that per-cut keeps.

=> NET-DISCOVERY hypothesis (my H4) CONFIRMED, and it is the DECISIVE finding: refcol trades a
+0.030 gain on WT (canonical) and +0.008 on R1 for a -0.021 LOSS on R3 (the actual non-canonical
discovery target the project exists to serve). The prior investigation / SYNTHESIS records FULL-RUN
preserving R3 at 0.258->0.258 (+0.000, discovery preserved). So on the ONE rung that is genuinely
about non-canonical discovery, full-run is NEUTRAL and refcol is NEGATIVE.

=> MISMEASUREMENT #2 (the headline yardstick hides the regression): the audit graded ref-column on
the FAIR panel's 0->4 count (where ALL truths are canonical -> R3 discovery is untestable there) and
concluded refcol "eliminates the canonical-demotion losses / keeps the wins." That grading was done
on a panel with NO non-canonical truths, so it could not observe that refcol's aggressiveness comes
at the cost of REAL non-canonical recovery on r3b. The fair-panel metric the recommendation rests on
is BLIND to the project's primary objective (non-canonical discovery).
- CHECKPOINT: R3 discovery regression (-33 true recoveries) persisted. Verifying fullrun R3 baseline next.

## PER-DISTINCT-JUNCTION DECOMPOSITION of the R3 regression (advisor-demanded; my-lens check on my OWN number)
R3 universe = only 4 distinct non-canonical junctions (4 tids) over 1600 reads (~400 reads/junction).
percut->refcol R3:
  LOST  = 42 reads across 3 distinct tids, 3 distinct true junctions (tid_10:22, tid_9:13, tid_8:7)
  GAINED=  9 reads across 2 distinct tids (tid_9:6, tid_8:3)
  NET by tid: 3 tids refcol WORSE (tid_10 -22, tid_9 -7, tid_8 -4), 0 tids refcol BETTER.
=> The -33 is NOT a single read-inflated coordinate. Refcol is net-WORSE on 3 of 4 distinct
non-canonical junctions and net-BETTER on NONE. Stated in the CORRECT unit (per distinct junction,
addressing the task's per-junction-vs-per-read concern applied to MY headline): refcol regresses
3/4 of the panel's non-canonical junctions. The read-multiplicity does NOT manufacture the finding;
it is systematic. (Caveat: only 4 distinct R3 junctions -> low junction-diversity panel; magnitude
per-read is amplified by ~400x multiplicity, but the DIRECTION holds per-junction on 3/4.)
- CHECKPOINT: per-junction decomposition persisted (3/4 junctions regress). Building fullrun_ng r3b next.

## FULLRUN R3 BASELINE — NOW MEASURED FROM IDENTICAL CODE (not imported)
The batch driver (pid 56113, still running its last arm) had ALREADY built mix_r3b_out fullrun_ng
(log: "refined 1161/4800"); the earlier disk listing was stale. I did NOT launch a competing refine
(4 workers of the driver were at ~85% CPU; only 3936 free pages; fileproviderd/Drive at 100%). I
adjudicated the driver-built arm read-only.
r3b per-rung recovery, ALL THREE arms from identical current code + identical adjudicator:
  rung | per-cut | fullrun (d)      | refcol (d)
  R0   | 0.979   | 0.978 (-0.001)   | 0.979 (+0.000)
  R1   | 0.878   | 0.879 (+0.001)   | 0.887 (+0.008)
  R3   | 0.258   | 0.258 (+0.000)   | 0.237 (-0.021)   <-- non-canonical DISCOVERY
  WT   | 0.781   | 0.781 (+0.000)   | 0.811 (+0.030)
R3 absolute: fullrun net=0 (GAINED=0 LOST=0, genuinely neutral, not an offset); refcol net=-33
(GAINED=9 LOST=42, 3/4 distinct junctions worse). => the pivotal "fullrun preserves R3, refcol
regresses R3" claim is now VERIFIED from identical code, no longer imported from SYNTHESIS.
- CHECKPOINT: fullrun R3 baseline verified (net 0). Writing verdict.

## ============ VERDICT (A/B measurement-validity lens) ============
fault_found = TRUE (in the RECOMMENDATION as stated; the numbers themselves are reproducible).

RE-DERIVATION: I independently reproduced every builder number EXACTLY (fair no-guard: fullrun
35/6 tier-shift, refcol 358/196; R0flank +0.005 / +0.011). The A/B pipeline is REPRODUCIBLE and the
tier-shift accounting is SYMMETRIC across arms (same _canonical_tier, same baseline file md5
06968d05..., buckets sum to totals -> no None-drop undercount). So the raw measurement is trustworthy.

THE MISMEASUREMENT is in the YARDSTICK the recommendation was graded on, not the arithmetic:
1. UNIT: the metric is per-READ, not per-JUNCTION (the task's suspicion is correct). With only 4
   distinct R3 junctions over 1600 reads (~400x multiplicity), per-read counts are inflated ~400x
   vs per-junction. I re-stated the decisive finding per-distinct-junction: refcol regresses 3/4 of
   the non-canonical junctions, improves 0 -> direction survives the unit correction.
2. WRONG PANEL FOR THE HEADLINE: the audit's "keep 35 wins / drop the 4-6 losses" prediction was
   graded on the FAIR panel, where EVERY true junction is CANONICAL (0 non-canonical truths in the
   measured rungs). That panel structurally CANNOT observe non-canonical discovery. The fair 0->4
   "losses" are not truth-recovery losses in the discovery sense; on the truth axis refcol is
   net-POSITIVE on fair (+59 recovered vs fullrun +27). The recommendation's own metric is blind to
   the project's purpose.
3. THE AXIS THAT MATTERS (r3b R3, non-canonical discovery) INVERTS the recommendation: fullrun is
   NEUTRAL (net 0), refcol REGRESSES (net -33, 3/4 junctions). So the audit's claim that ref-column
   "keeps the wins AND eliminates the losses" is doubly wrong: (a) on fair it does NOT drop losses,
   it balloons both directions (358/196); (b) on r3b it introduces a NEW, real regression fullrun
   does not have. Ref-column is MORE sensitive, not less (my sensitivity hypothesis H3 was inverted),
   and that extra sensitivity is NET-HARMFUL on non-canonical discovery.

WINNER FROM MY LENS: the recommendation to prefer REFERENCE-COLUMN is NOT SUPPORTED by a trustworthy
measurement. On the discovery axis, per-cut and full-run tie (both 0.258 on R3); refcol loses. The
audit's model-correctness argument (calibration axis) may still be right in principle, but the
EMPIRICAL A/B does not deliver the predicted "keep 35 / drop 6" outcome — it is a different, worse
outcome on the discovery-relevant panel. Do NOT adopt ref-column on the strength of this A/B.

CHALLENGE-TO-MY-OWN-FINDING (advisor-folded): could the R3 loss be 1-2 read-inflated coordinates?
CHECKED: no — 3/4 distinct junctions regress, 0 improve. Could fair's +59 rescue refcol? No — fair
has no non-canonical truths, so its "recovery" is canonical-coordinate landing, orthogonal to the
discovery objective. Does my verdict rest on an imported fullrun number? No longer — built from
identical code, R3 net 0 confirmed.

CAVEAT (honest scope): only 4 distinct R3 junctions -> low junction diversity; the magnitude is
panel-specific and per-read multiplicity is high. The DIRECTION (refcol regresses non-canonical
discovery; fullrun neutral) is robust to the unit correction, but a broader-diversity panel is
needed to quantify the effect size. This is a measurement-scope limitation of the A/B ITSELF (small
non-canonical junction alphabet), which is itself a validity finding: the fair panel that the
recommendation leans on has ZERO non-canonical junctions, and the r3b panel has only 4.

STATUS: COMPLETE.
