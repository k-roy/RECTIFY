# INSCOST AUDIT — Re-validation Completeness Lens

Adversarial auditor of the claim: "FULL-RUN ins_cost is sounder AND safe to switch (default ON)."
Lens: RE-VALIDATION COMPLETENESS — switching changes results, so EVERY validated finding must be re-checked.
Task: find the flaw that makes "full-run is sounder / safe to switch" WRONG.

## Plan
1. Read dev/INSCOST_INVESTIGATION.md (durable record).
2. Read junction_scoring.py:_score_junction, _score_hp_anchored, ins_cost machinery.
3. Read junction_refiner.py HP-drift guard (hp_drift_margin) — the ONE subtlety.
4. Enumerate EVERYTHING that depends on the score scale / ins_cost:
   - yeast DRS transfer
   - human chr5 / GENCODE transfer
   - hp_drift_margin=3.0 tuning (fit to OLD scoring?)
   - del_open REJECT verdict (arm-F)
   - the 46 refiner/guard tests (pin exact CIGARs/scores?)
   - sim panel calibration (penalty tables + margins)
   - cat1-cat9 fixtures
5. For each: does the investigator's revalidation_plan cover it? Could it FLIP?
6. Toggle flag in-process to gather concrete witnesses where cheap.
7. Verdict.

## Checkpoints
- (start) file created, beginning orientation.
- read INSCOST_INVESTIGATION.md + junction_scoring.py (_score_junction/_score_hp_anchored) + junction_refiner.py move-gate + DEL_OPEN_DELTA_FINDING.md.

## ENUMERATION — what depends on the score scale, and plan coverage

Investigator 5-item plan: (1) full test suite flag-ON; (2) yeast DRS transfer;
(3) human ONT DRS transfer; (4) del_open/arm-F re-run; (5) hp_drift_margin re-tune.

COVERED: full suite (#1), yeast (#2), human (#3), arm-F (#4), hp_drift_margin (#5).

MISSED by the plan (score-scale-coupled, NOT re-checked):
- M1. **`_CANONICAL_HP_PRIOR = 0.5`** (junction_scoring.py:310). Docstring: "0.5 equals
  one Nanopore HP-deletion equivalent — THE EXPECTED NOISE FLOOR for splice-site scoring."
  It is an edit-distance discount compared to raw scores in the move decision
  (refiner:686-687 `score_cmp = score - 0.5`). Full-run inflates HP-insertion scores
  ~4.7x (witness 1.76->8.26); the 0.5 "noise floor" stays fixed -> now ~4.7x too small
  relative to the scores it arbitrates. NOT in the re-tune plan (#5 re-tunes ONLY
  hp_drift_margin).
- M1-CRITICAL. **The A/B panels are ALL `motif_blind`.** refiner:655
  `tier_beats_alt = (current_tier>=4) and not motif_blind`. In motif_blind mode
  `_CANONICAL_HP_PRIOR`, tier, and novel priors are STRUCTURALLY DISABLED. So the
  investigator's ENTIRE validation (mix_fair/mix_r3b, arm-C style = motif_blind +
  penalty_table) NEVER exercised the canonical noise-floor discount. The DEFAULT
  production pipeline is motif-AWARE (motif_blind=False, arm-A), where the fixed 0.5
  discount IS compared against full-run's 4.7x-inflated HP scores. => the flag's
  interaction with the default mode's noise floor was measured on ZERO reads, and none
  of the 5 plan items exercise it (they reuse the motif_blind discovery harness).
- M2. **`hold_margin`** (the BLUNT move-prior on ALL moves). Only hp_drift_margin is
  in #5. The guard-on counter-signal reads (4 reads/panel, 0->4 canonical->non-canonical,
  R0flank 0.899->0.898) are NOT into-HP moves, so the guard misses them; only hold_margin
  (blunt) would catch them. Omitted from re-tune.
- M3. **Two LIVE per-substring ins-cost consumers NOT switched:**
  (a) `corrected_consensus._cigar_hp_edit_distance` — the PRIMARY aligner/winner-selection
      sort key (CLAUDE.md: merge settled on hp_edit_distance). Uses penalty_table ins_cost
      the OLD way.
  (b) `splice_aware_5prime` Cat3 5' rescue exon-vs-intron placement via `_hp_edit_distance`
      (lines 1542/1704/1732/1746/2033/2034).
  If full-run is "THE sounder insertion model," these are now provably inconsistent with it,
  and their prior validation (cat3 fixtures, winner selection) is stale under the new model.
  The flag doesn't change them (so no "changed result"), but the "sounder model" CLAIM is
  applied to 1 decision and silently contradicted at 2 live others. Plan neither switches
  nor re-validates them.
- M4. **Ship-config reframing (from investigator's OWN numbers).** The intended production
  config is arm-E = motif_blind + hp_drift_margin=3.0 (DEL_OPEN doc: "ship arm-E").
  Guard-ON flag OFF vs ON: 0.07% (fair) / 0.10% (r3b), R0flank 0.899->0.898 (-0.001),
  R3 preserved. => in the SHIPPED config full-run's recovery/anti-fabrication benefit is
  ZERO-to-slightly-NEGATIVE. The "net-positive 35:6, +0.005" is a NO-GUARD artifact that
  evaporates in production. Corroborated by DEL_OPEN's deep reason: HP boundary-absorption
  is a SEARCH degeneracy that cost-term changes (arm-F del, and now full-run ins) do NOT
  fix — "only a prior (the guard) can." Full-run is another HP-axis cost-term tweak, so
  the guard subsumes it (exactly what the guard-on data shows). Leaves ONLY the (unbuilt,
  still needs MECH1+MECH3) DP speedup as the real justification.

CONFIRMED NON-ISSUES (empirically checked / not a flip risk):
- 46 refiner/guard tests: investigator RAN them flag-ON = 46 passed (same as OFF). Not a flip.
- `_score_junction_fallback` per-cut split: grep shows NO live caller (dead/tests only). Weak.
- 0.284 anchor: penalty_table=None => flag is structural no-op. Trivially preserved (correct).
- score_bin "floor(score)": refiner docstring line 24 is STALE; code uses raw score_cmp. No live floor.

## STRONGEST CHALLENGE (holds)
The full-run A/B was validated EXCLUSIVELY in motif_blind mode, which structurally
DISABLES the scale-coupled canonical priors (_CANONICAL_HP_PRIOR "noise floor"=0.5,
tier, novel). The DEFAULT production pipeline is motif-AWARE, where 0.5 is actively
compared against full-run's ~4.7x-inflated HP scores. So the switch changes the score
scale by ~4.7x on HP runs while leaving the coupled constants fixed AND un-tested in the
mode where they fire. The plan re-tunes 1 of >=3 scale-coupled constants (hp_drift_margin
only; omits _CANONICAL_HP_PRIOR + hold_margin), tests 0 of them in default mode, and
re-validates 0 of the 2 live per-substring ins consumers (consensus winner-key, Cat3
5' rescue). The A/B is therefore CONFOUNDED (full-run-with-stale-margins vs
per-cut-with-margins-tuned-for-per-cut), and the counter-signal (guard-missed 0->4
canonical demotions) is the exact symptom of an un-re-tuned noise-floor constant.
=> "safe to switch" is UNPROVEN; a valid switch needs best-tuned-full-run >= best-tuned-
per-cut across the FULL margin set in the DEFAULT motif-aware mode, plus consistency
of the 2 un-switched consumers. Following the 5-item plan would NOT establish this.

## ADVISOR DISCRIMINATOR RESOLVED (config intersection) — CONFIRMED 2026-07-08
Advisor asked: does a REAL config have BOTH penalty_table set (flag active) AND
motif_blind=False (_CANONICAL_HP_PRIOR active)? Grep of the production path answers YES,
and it is the MAINLINE `rectify correct` path:
- correct_command.py:746 `refine_bam_junctions(...)` passes
  `penalty_table_path=config.get('junction_penalty_table')` and AUTO-BUNDLES a junction
  penalty table for ont_cDNA (lines 706-734) => **penalty_table ACTIVE => flag BITES.**
- It does NOT pass `motif_blind` => defaults **False => motif-AWARE =>
  _CANONICAL_HP_PRIOR(0.5)+tier+novel priors ALL ACTIVE.**
- It does NOT pass `hold_margin` or `hp_drift_margin` => both default **0.0 => the HP-drift
  GUARD is OFF in the mainline product.**
- `motif_blind` is NOT a CLI flag (grep of commands/ = empty) => users cannot get the
  discovery arm; the shipped `rectify correct` is ALWAYS motif-aware.

=> The MAINLINE production config = {penalty_table ON, motif_blind=False, guard OFF}.
The investigator validated the flag ONLY on motif_blind panels (priors OFF) and reasoned
about safety via guard-ON (=3.0) numbers. NEITHER matches the mainline. So the flag's
default-ON behavior in the shipped `rectify correct` was measured on ZERO reads:
  * flag x _CANONICAL_HP_PRIOR interaction (fixed 0.5 "noise floor" vs 4.7x-inflated HP
    scores) is UNTESTED;
  * mainline has NO guard and NO hold_margin (both 0.0) to absorb the score-scale shift,
    so the 0->4 canonical->non-canonical demotion risk is MAXIMAL there (protection is
    only the exact-tie is_alt breaker) — and the "guard subsumes full-run" reassurance
    (M4) is about a config the mainline does not run.
M1 and M4 are therefore TWO DISTINCT REAL CONFIGS, both unvalidated for net benefit:
config-1 (mainline motif-aware+table+guard-off) untested at the flag×prior intersection;
config-2 (discovery motif_blind+guard@3.0) shows the benefit subsumed by the guard.

## M3 CORRECTED (advisor catch — DROPPED as a per-cut inconsistency)
Read `_cigar_hp_edit_distance` (corrected_consensus.py) + `_hp_edit_distance` usage:
they do NOT k-sweep. `_cigar_hp_edit_distance` charges insertions as
`length * ins_cost(_hp_run_length(GENOME, ref_pos), ...)` (genome-context, per-CIGAR-op);
splice_aware compares full-sequence `_hp_edit_distance(rescue, exon)` vs
`_hp_edit_distance(rescue, intron)` (no per-cut split). => NOT the per-cut gaming artifact;
"same-kind inconsistency" claim is WRONG. DROPPED from the strongest challenge. (Residual,
non-load-bearing: consensus uses a genome-context ins model, a third formulation — a minor
soundness nuance, not a re-validation gap.)

## VERDICT SPLIT
- "sounder MODEL" (over-call = property of the read's one physical run): SURVIVES —
  principled, calibration-independent. Concede.
- "safe to switch / default ON": DOES NOT SURVIVE as established — plan incomplete on
  >=3 coupled constants + default-mode interaction + 2 live consumers; a listed-but-
  unmeasured comparison (re-tuned margins in motif-aware mode) COULD flip net-positive to
  net-negative. Conjunction fails => verdict_survives = False.
