# INSCOST Synthesis — full-run vs per-cut ins_cost switch

Status: IN PROGRESS. System of record for orchestrator recovery.

## Task
Synthesize investigator + 2 audits on switching default from per-cut ins_cost
to full-run ins_cost (`_USE_FULL_RUN_INS`) in `rectify/core/splice/junction_scoring.py`.
Decide: (1) sounder model? (2) re-validation complete / at-risk-of-flipping?
(3) PROCEED vs HOLD + checklist. (4) result most likely to change materially.

## Inputs summary
- INVESTIGATOR: full-run IS sounder (principled error-model + DP-unlock). Did NOT flip default,
  did NOT build concat DP. Witness reproduced: OFF=1.7604 (per-cut split k=6), ON=8.2584 (full-run 12*ic12).
- AUDITOR-1 (impact-validity): fault_found=FALSE, verdict survives. Falsified own attacks:
  guard is OPT-IN (hp_drift_margin default 0.0); motif prior does NOT subsume; ran
  motif_blind=False production config -> discovery preserved (r3b R3 0.258->0.258, +0.000).
- AUDITOR-2 (revalidation-completeness): fault_found=TRUE, "safe to switch" half does NOT survive.
  Claims mainline {table ON, motif_blind=False, guard OFF} measured on ZERO reads; flag x
  _CANONICAL_HP_PRIOR (fixed 0.5 noise floor vs 4.7x inflated scores) untested; plan re-tunes
  only hp_drift_margin, not _CANONICAL_HP_PRIOR.

## KEY CONFLICT to resolve
Auditor-1 says it RAN motif_blind=False no-guard config (discovery preserved).
Auditor-2 says that config measured on ZERO reads and _CANONICAL_HP_PRIOR interaction untested.
=> Crux: is _CANONICAL_HP_PRIOR (junction_scoring.py:310) gated by motif_blind, or does it fire
   regardless? refiner:655 `tier_beats_alt=(current_tier>=4) and not motif_blind` gates a
   DIFFERENT thing (tier-beats-alt), not necessarily the HP prior.
   Need to read the code to break the tie.

## Plan
1. Read junction_scoring.py around :310 (_CANONICAL_HP_PRIOR), :60 (flag), :710/:839 (scoring fns).
2. Read junction_refiner.py :655, :729-741 (motif_blind gate, move-gate).
3. Read correct_command.py :746 (mainline call — confirm motif_blind/guard defaults).
4. Resolve: does motif_blind=False actually leave _CANONICAL_HP_PRIOR active AND independent of
   the flag rescale? Did Auditor-1's run exercise it?
5. Synthesize verdict.

## Checkpoints
- (start) file created, plan written.
- CODE VERIFIED on shipping branch `drs-validation-rebuild` (flag uncommitted, lives on
  investigator's benchmark worktree; motif_blind ALSO absent here — task #8 added it there):
  * Flag/precomputed_ins_costs/full_ins_costs: NOT in this worktree (grep empty).
  * `motif_blind`: absent from entire rectify/ tree on this branch.
  * junction_refiner.py:616 `tier_beats_alt = current_tier >= 4` — NO `and not motif_blind`.
    => On the SHIPPING branch, `_CANONICAL_HP_PRIOR` (0.5) is ALWAYS applied when current
       N-op is non-canonical (tier>=4) and a canonical (tier<4) alternative exists (:647).
  * PER-CUT MECHANISM CONFIRMED (junction_scoring.py:766-769): ins_costs computed as
    `penalty_table.ins_cost(_hp_run_length(query, i), query[i])` where `query` = the SLICED
    segment (rescue[k:] for t1, rescue[:k][::-1] for t2). _hp_run_length is measured WITHIN
    the segment => a 12-A run split at k=6 becomes two 6-A runs, each charged ic(6)=0.1467
    (table min), total 12*ic(6)=1.76, vs full-run 12*ic(12)=8.26. Gaming reproduced exactly.
  * Only `_score_junction`'s k-sweep exploits per-segment run length (Auditor-2's confirmation
    that the other ins_cost consumers don't k-sweep stands). Blast radius correct.

## CONFLICT RESOLUTION (the two auditors)
- Auditor-2 LITERAL claim "mainline motif-aware config measured on ZERO reads" is REBUTTED by
  Auditor-1, who ran motif_blind=False (prior ACTIVE) + table + no-guard on BOTH panels:
  fair 42/5600 (35 reads 4->0, R0flank +0.005), r3b R3 0.258->0.258 (+0.000), 7/4800 balanced.
  On the SHIPPING branch the prior is always-on, so Auditor-1's motif_blind=False runs ARE the
  correct production proxy. Empirical net-harm question for SIM panels = answered (neutral+).
- Auditor-2's SURVIVING valid point (narrower, real): the switch rescales HP-INSERTION scores
  ~4.7x, and `_CANONICAL_HP_PRIOR=0.5` was calibrated to the per-cut scale ("one HP-DELETION
  equivalent"). The flag doesn't touch DEL cost, so the prior's ANCHOR is unchanged — BUT the
  prior is a fixed offset on score gaps that (when HP-insertion-driven) now rescale, so its
  EFFECTIVE strength vs those gaps erodes. A single neutral measurement at legacy 0.5 does NOT
  prove 0.5 is still right on the new scale, nor that best-tuned-full-run >= best-tuned-per-cut.
  The 5-item plan re-tunes only hp_drift_margin, omitting _CANONICAL_HP_PRIOR + best-vs-best.
  => genuine RE-VALIDATION-COMPLETENESS gap, gating the DEFAULT-ON flip (not the model).

## VERDICT (draft, pre-advisor)
(1) full_run_is_sounder = TRUE. UNANIMOUS on model: Auditor-2 EXPLICITLY concedes the model
    half ("charging an over-call at an arbitrary cut k has no biological meaning; calibration-
    independent"). Per-cut is a search/scoring degeneracy exploiting the U-shaped table min at
    hp6. Investigator: per-cut is never MORE correct than full-run.
(2) revalidation_complete = FALSE. Missing: (a) _CANONICAL_HP_PRIOR re-tune to new scale +
    best-tuned-full-run >= best-tuned-per-cut in shipping (prior-always-on) mode; (b) del_open
    arm-F re-run (scale shifted 4.7x); (c) human ONT transfer (direction UNKNOWN); (d) numba-ON
    cluster path flag-ON; (e) guard re-sweep (may drop <3.0).
(3) final_call = PROCEED-TO-SWITCH, GATED: adopt full-run as reference model / land flag, but
    do NOT ship default-ON until the ordered checklist below is green (esp. the _CANONICAL_HP_
    PRIOR re-tune + best-vs-best, and the human transfer). Investigator already did NOT flip.
(4) most_likely_to_change = the HUMAN ONT DRS transfer (investigator item #3).
    READING USED: "the result that can overturn the GO/NO-GO bottom line." It is the ONLY
    external-validity test; every neutral/positive number to date is SIM; the investigator
    itself flags human as direction-UNKNOWN (human HP-length distribution differs => could
    shift the fabrication/recovery balance — direction unknown, NOT assumed worse). The
    disclosed 0->4 canonical-demotion counter-signal is the specific quantity to watch under
    shipping always-on _CANONICAL_HP_PRIOR.
    SEPARABLE NOTE: del_open / arm-F is the result MOST LIKELY TO LITERALLY FLIP (an existing
    verdict made on the OLD score scale; scale shifted 4.69x, 1.76->8.26; investigator says it
    "must be re-run, not assumed to carry") — but it is a separable sub-verdict, not the go/no-go.

## ADVISOR-FOLDED CORRECTIONS
- CROSS-BRANCH GAP (neither auditor caught): flag is UNCOMMITTED on the benchmark worktree,
  which ALSO carries motif_blind (task #8); the shipping/merge target (master) lacks both. All
  credited impact numbers were measured on the benchmark branch. The motif_blind=False ==
  shipping-always-on-prior equivalence holds for THAT ONE line only; it assumes the rest of the
  scoring code is identical across branches. => CHECKLIST ITEM: commit the flag and re-run the
  gate on the ACTUAL merge target; confirm no other scoring-relevant divergence.
- GATE FRAMING FIX: do NOT adopt Auditor-2's "best-tuned-full-run >= best-tuned-per-cut". Per-cut
  is a degeneracy (exploits the U-shaped table hp6 min) — tuning _CANONICAL_HP_PRIOR up to let it
  win a sim metric would overfit the artifact. Correct gate: (i) confirm 0.5 still represents the
  intended noise floor on the new full-run scale (anchor = one HP-DELETION equivalent, and DEL
  cost is UNCHANGED by the flag, so the anchor HOLDS; only the prior's effective strength vs
  HP-insertion-driven gaps erodes), and (ii) confirm NO net discovery regression (esp. human).
- DP-unlock / ~30x speedup is NOT a switch justification: unbuilt, admittedly not byte-identical
  (needs MECH1 boundary fix + MECH3 FP tolerance), numba path unexercised. Contingent downstream
  benefit only. The switch stands on model-correctness alone.

## ORDERED RE-VALIDATION CHECKLIST (gates the DEFAULT-ON flip; investigator did NOT flip)
1. Commit the flag on / rebase to the ACTUAL merge target (master); diff benchmark vs master for
   any other scoring-relevant divergence beyond the flag + motif_blind line. Re-establish the
   baseline on the branch the switch actually lands on.
2. FULL suite flag-ON on the merge target (`pytest -m "not slow"`, ~1603) incl. cat1-cat9
   fixtures — local green gate.
3. numba-ON DP path, flag-ON, on the CLUSTER build (verify reversed-slice list->float64 array
   feeds _score_hp_dp_numba correctly) BEFORE trusting any cluster run.
4. _CANONICAL_HP_PRIOR re-confirm on the full-run scale: verify 0.5 still = intended noise floor
   (NOT "beat per-cut"); check no net discovery regression on fair + r3b, no-guard (prior always
   on = shipping proxy).
5. Guard hp_drift_margin re-sweep (task #16): predicted smallest zero-discovery-cost margin may
   DECREASE below 3.0 since full-run does part of the drift suppression itself.
6. del_open / arm-F verdict RE-RUN under the full-run scale (DEL_OPEN_DELTA_FINDING) — do not
   assume the per-cut-era verdict carries.
7. Yeast DRS real-data HP-drift transfer (task #17, Sherlock).
8. Human ONT DRS transfer (task #18, GENCODE truth) — the decisive external-validity test;
   direction unknown; watch the 0->4 canonical-demotion balance.
Flip default-ON to ship ONLY after 2/3/4 are green and 7/8 show no regression.

## STATUS: COMPLETE.
