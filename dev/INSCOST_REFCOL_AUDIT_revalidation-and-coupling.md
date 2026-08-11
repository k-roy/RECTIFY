# INSCOST REFCOL AUDIT — Downstream Coupling + Re-Validation Lens

**Auditor role:** adversarial. Lens = DOWNSTREAM COUPLING + RE-VALIDATION.
**Question:** Whichever ins_cost variant wins (full-run vs reference-column) changes
the SCORE SCALE. Constants tuned to the per-cut regime must be re-checked:
- `_CANONICAL_HP_PRIOR` (0.5 noise floor, junction_scoring.py:~310)
- `hp_drift_margin` (3.0)
- del_open / arm-F verdict (scale shift could flip it)

**Core suspicions to test:**
1. Does the builder's A/B hold `_CANONICAL_HP_PRIOR` FIXED at 0.5 while ins scale moved
   → biasing the comparison toward whichever variant happens to sit near per-cut scale?
2. Is ref-column scale CLOSER to per-cut (less re-tuning) or FURTHER (more)?
3. Enumerate everything that must be re-validated for the recommended winner.
4. Is the recommendation safe to act on BEFORE that re-validation?
5. Cross-branch: flags uncommitted on THIS worktree (has motif_blind); master has neither.

## Plan
- [ ] Read the three investigation docs (INVESTIGATION, SYNTHESIS, AUDIT_model-correctness)
- [ ] Read INSCOST_REFCOL_BUILD.md (builder durable record — may be null)
- [ ] Read junction_scoring.py around: _score_junction(807), _score_hp_anchored(693),
      _CANONICAL_HP_PRIOR(~310), hp_drift_margin, del_open, _precompute_del_costs
- [ ] Determine how the A/B was run: were coupled constants held fixed?
- [ ] Reason about scale: per-cut vs full-run vs ref-column magnitude of ins_cost sum
- [ ] Enumerate re-validation checklist for the winner
- [ ] Verdict: is the recommendation safe to act on now?

## Checkpoints
- [read] All 4 source docs read: INVESTIGATION, SYNTHESIS, AUDIT_model-correctness, BUILD.
- [KEY] BUILD.md line 186-205: the builder ALREADY RAN the 3-way A/B. The audit hypothesis
  is **FALSIFIED**. percut->fullrun fair no-guard: 35 (4->0 wins) : 6 (0->4 losses), churn
  42/5600=0.75%, R0flank +0.005. percut->refcol fair no-guard: 358 (4->0) : 196 (0->4),
  churn 694/5600=12.39% (16x MORE), R0flank +0.011. Ref-column does NOT eliminate the
  0->4 demotion class — it BALLOONS it (196 vs 6). Recovery higher but via massive churn
  both directions. Ref-column is a far more AGGRESSIVE re-scorer, NOT a clean dominator.
- [CODE] junction_refiner.py:685-689 — `score_cmp = score - canonical_discount` where
  canonical_discount = _CANONICAL_HP_PRIOR = 0.5 (ADDITIVE on the raw _score_junction sum).
  :729-738 — `eff_margin = hold_margin (+ hp_drift_margin=3.0 if into_hp)`; veto if
  `best_score_cmp > incumbent_score - eff_margin`. So BOTH coupled constants (0.5 prior,
  3.0 margin) are FIXED ADDITIVE OFFSETS on the same scale as the ins/del-cost sum. When
  the ins-cost scale moves, these offsets do NOT rescale -> their EFFECTIVE stringency vs
  score gaps changes. This is the coupling my lens targets, and it is REAL.
- [SCALE from witnesses] Wm (genuine 4-base over-call at 8-A ref acceptor):
  per-cut 0.7772 | full-run 2.7528 | refcol 0.7888. W1 (pure-gaming all-A12):
  per-cut 1.7604 | full-run 8.2584 | refcol 6.8748.
- [HARNESS] _make_ins_arms_refcol.py: motif_blind=False (prior ACTIVE — shipping proxy),
  varies ONLY hp_drift_margin (0/3.0) + ins flag. _CANONICAL_HP_PRIOR is a HARDCODED
  module constant (junction_scoring.py:338) imported into refiner:136 — NOT a parameter.
  => CONFIRMED: the A/B held _CANONICAL_HP_PRIOR FIXED at 0.5 across all 3 variants
  while the ins scale moved. This is exactly the bias my lens flagged. The question is
  DIRECTION + magnitude.

## SCALE ANALYSIS (my lens core finding)
The 0.5 prior is applied ADDITIVELY: `score_cmp = score - 0.5` for canonical alternatives
(non-canonical candidates get 0). A non-canonical current placement is DEMOTED to a
canonical alt when (canonical_score - 0.5) <= (noncanon_score). The prior's EFFECTIVE
tipping power depends on the score SCALE of the ins-cost-driven gaps it competes with.

Witness Wm (genuine over-call) magnitudes tell the scale story:
- per-cut  0.7772  <- prior 0.5 is ~64% of this score. Big lever.
- refcol   0.7888  <- prior 0.5 is ~63%. ESSENTIALLY THE SAME SCALE AS PER-CUT (+1.5%).
- full-run 2.7528  <- prior 0.5 is ~18%. Prior is 3.5x WEAKER relative to the score.

CRUCIAL INFERENCE for the bias question:
- FULL-RUN moved the scale 3.5x UP -> the fixed 0.5 prior became RELATIVELY WEAKER ->
  the prior tips FEWER ties toward canonical -> full-run's demotion count should be
  SUPPRESSED by the un-retuned prior (a prior that's too weak lets non-canonical hold,
  it does NOT create canonical demotions). Consistent with full-run's LOW 0->4 count (6).
- REF-COL scale ~= per-cut scale (Wm 0.7888 vs 0.7772). So the 0.5 prior sits at
  ~the SAME relative strength as it was tuned for. Ref-col is the LEAST-mistuned variant
  w.r.t. the prior. => LESS re-tuning of the prior needed for ref-col than for full-run.

## CORRECTION (advisor-folded) — the single-witness scale claim is FRAGILE; do NOT lean on it
- I was about to make "refcol scale == per-cut (Wm 0.7888 vs 0.7772) => prior correctly
  scaled => the 196 losses are INTRINSIC" load-bearing. WRONG. It rests on ONE witness.
  W1 CONTRADICTS it: refcol 6.8748 vs per-cut 1.7604 (~4x, ~= full-run 8.2584). Ref-col's
  scale is BIMODAL BY DESIGN: ~per-cut in the genuine-over-call regime (Wm), ~full-run in
  the gaming/A-free regime (W1). So "refcol is close to per-cut" is only half true.
- The REAL resolution of the tension (scale looks per-cut-ish yet 196 demotions) is the
  **12.39% churn** itself: refcol is BEHAVIORALLY THE FARTHEST of the three from per-cut
  (694/5600 calls change vs full-run's 42/5600). So the task's scale sub-question INVERTS:
  ref-col is NOT "closer to per-cut -> less re-tuning". Behaviorally it demands the MOST
  re-validation of the three. That is the clean lens finding.

## BIAS QUESTION — ANSWER (the task's central question)
Q: Does holding _CANONICAL_HP_PRIOR FIXED at 0.5 while ins scale moved BIAS the comparison?
A: For the recommendation AS STATED (a bare flag-flip to refcol default), NO — because 0.5
   IS the shipping value. The A/B at prior=0.5 is a VALID production-proxy of exactly the
   flip being recommended. Holding it fixed does not bias the test of that flip; it tests it
   faithfully. And under that faithful test refcol LOSES on its own terms (196 demotions vs
   6; 16x churn). The audit hypothesis (refcol keeps ~35 wins AND eliminates the 0->4 losses)
   is EMPIRICALLY FALSIFIED — needs nothing more to reject the recommendation.
   The ONLY way the fixed prior "biases" anything is against a DIFFERENT, un-recommended
   variant: refcol+retuned-prior. That is a separate proposal, untested, and if pursued it
   would itself need best-vs-best + del_open/arm-F re-run + guard re-sweep. So the fixed
   prior does not rescue the recommendation; it just means "refcol-as-recommended fails, and
   the escape hatch is unvalidated."

## MECHANISM of the 196 demotions (reasoned; arms did not persist, do NOT re-run)
Prior fires ONLY when current_tier >= 4 (refiner:655; discount at :685-687). So the 0->4
demotions split into two subsets:
  (A) non-canonical incumbents where per-cut made the canonical fix and refcol did NOT
      -> prior ACTIVE for these; a prior re-tune COULD touch them.
  (B) canonical incumbents refcol demoted to non-canonical on RAW score
      -> prior is OFF (current_tier < 4 => tier_beats_alt False => score_cmp = raw score).
      A prior re-tune canNOT touch subset (B) at all.
=> The demotion class is NOT curable by prior re-tune alone. This is WHY refcol is not safe
   to act on: not "losses provably intrinsic" (I retract that phrasing), but "the only
   escape hatch (prior re-tune) is unvalidated, touches only subset (A), and is a DIFFERENT
   recommendation than the one on the table."

## RE-VALIDATION CHECKLIST — what must be green BEFORE acting on ANY winner
(Every item is made RISKIER for refcol by its 12.39% churn vs full-run's 0.75%: 16x more
calls change => 16x more surface for each check to flip.)
1. _CANONICAL_HP_PRIOR (0.5) re-confirm — SCOPED to the prior-active regime ONLY (current_tier
   >= 4). NOT "beat per-cut" (per-cut is a degeneracy). Confirm 0.5 still = intended noise
   floor on the chosen variant's scale. NOTE: for full-run the prior is 3.5x RELATIVELY WEAKER
   (0.5/2.75 vs 0.5/0.78); for refcol it is bimodal (per-cut-like in over-call, full-run-like
   in gaming). NEITHER re-confirmed on-scale.
2. hp_drift_margin (3.0) re-sweep — 3.0 is an ABSOLUTE additive constant tuned to per-cut
   score deltas. full-run rescales HP scores ~4.7x (W1 1.76->8.26); refcol rescales bimodally.
   Guard effective stringency is UNKNOWN on either new scale. Under the guard the flag-effect
   collapsed to 0.07% (fair) for full-run — the guard MASKS the ins-cost differences, so the
   no-guard A/B (where refcol showed 196 demotions) is the discriminating config, and the
   guard=3.0 config is NOT a safe substitute for re-tuning.
3. del_open / arm-F verdict RE-RUN — reached under per-cut ins scale; both variants shift the
   scale (full-run 4.7x; refcol bimodal). Score-scale shift can FLIP this sub-verdict. Do not
   assume it carries.
4. FULL suite flag-ON (pytest -m "not slow", ~1603) + cat1-cat9 fixtures on the chosen variant.
   Only the narrow 46-suite + a 201-subset were run flag-ON for full-run; refcol flag-ON full
   suite NOT run. (refcol byte-identity OFF confirmed = 46 tests; that only proves default path
   unchanged, NOT flag-ON correctness.)
5. numba-ON DP path flag-ON on the CLUSTER build. numba is OFF locally (must stay off here).
   refcol FORCES pure-Python (BUILD line 44-46: its ic vector is length R+1 column-indexed,
   WRONG for the Q-indexed _score_hp_dp_numba kernel) => refcol has NO numba path yet; the
   ~30x DP-unlock for refcol requires a NEW column-indexed numba kernel, unwritten/unverified.
6. Human ONT DRS transfer (task #18 truth=GENCODE) — the decisive external-validity test;
   direction UNKNOWN. refcol's 196-demotion class is the specific quantity to watch; on human
   HP-length distributions the demotion balance could worsen.
7. CROSS-BRANCH — the ENTIRE A/B ran on THIS benchmark worktree, which carries motif_blind
   (task #8) + both uncommitted flags. master has NEITHER. The motif_blind=False == shipping-
   always-on-prior equivalence must be RE-ESTABLISHED on the actual merge target; diff
   benchmark vs master for any other scoring-relevant divergence before trusting the numbers.

## VERDICT (my lens: downstream coupling + re-validation)
- fault_found = TRUE. The recommendation to build+adopt REFERENCE-COLUMN as the winner does
  NOT survive: the builder's own A/B (prior held at the shipping value 0.5) FALSIFIES the
  audit's load-bearing prediction (keep ~35 wins, eliminate 0->4 losses). Refcol instead
  produces 196 demotions (33x more than full-run's 6) and 12.39% churn (16x full-run's 0.75%).
- The task's bias hypothesis (fixed prior biases the comparison) is REFUTED as a rescue for
  refcol: 0.5 is the shipping value, so the A/B faithfully tests the flip-as-recommended.
- The task's scale sub-question INVERTS: refcol is NOT closer to per-cut (bimodal; behaviorally
  the FARTHEST at 12.39% churn) -> it needs the MOST re-validation, not the least.
- winner_from_my_lens: refcol is OUT. No default flip of EITHER remaining variant is safe until
  items 1-7 are green. full-run carries FAR less re-validation exposure (0.75% churn, numba
  path structurally present, already partially suite-tested) and remains the credible
  candidate; status-quo per-cut is the safe hold. Recommendation is NOT safe to act on now.
- verdict_survives = FALSE (the "adopt reference-column" verdict does not survive).

STATUS: COMPLETE.
