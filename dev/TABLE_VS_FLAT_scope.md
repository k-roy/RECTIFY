# TABLE vs FLAT (SCOPED) — does ANY empirical-table formulation + HP-drift guard beat flat + guard?

## The challenge
PI not convinced dropping the empirical DRS error penalty table (in favor of flat hand-set costs)
is justified for the native junction re-placer's boundary search.

## My angle (SCOPED rejection, not general)
arm-C / del_open_delta rejected table formulations that CHEAPEN IN-RUN INDELS (grease the drift).
Build the STRONGEST table use that AVOIDS that:
- (a) table for SUBSTITUTION / mismatch-context costs ONLY (leave del/ins flat) — mismatch cost
  doesn't move a boundary into a run.
- (b) a table variant that does NOT make in-run deletions cheaper than flat.
- (c) table + HP-drift GUARD together vs flat + GUARD (maybe guard neutralizes the over-shift so
  the table's accuracy benefit shows).

A/B each vs flat, GUARD ON, on mix_fair_out + mix_r3b_out, recovery vs truth.
Question: does ANY table formulation + guard BEAT flat + guard? If yes, "drop the table" is WRONG.

## Constraints
- M1 kills heavy/concurrent refines. REUSE existing aligned.sorted.bam. Single lean refines. Small panels.
- Python /Users/kevinroy/miniconda3/bin/python. numba OFF.

## Plan
1. Read junction_scoring.py (_score_junction, _score_hp_anchored, flat constants, table path).
2. Read HpPenaltyTable / penalty_scores.tsv structure — what columns exist, what's measured.
3. Read paired_arm_test.py + _make_arm_e.py + refine_bam_junctions — how recovery-vs-truth is scored.
4. Locate the sim panels (mix_fair_out, mix_r3b_out) + their existing bams + truth.
5. Design table formulations (a)/(b)/(c) that avoid cheapening in-run indels.
6. Run lean A/B refines GUARD ON: flat+guard vs each table+guard. Recovery vs truth.
7. Verdict: does any table+guard beat flat+guard?

## Checkpoints
(append one line per sub-step)

## CHECKPOINT 1 — orientation done (mechanics confirmed)
- `_score_junction`->`_score_hp_anchored` (junction_scoring.py). FLAT costs: _FLAT_SUB=1.0, _FLAT_DEL_NORMAL=1.0,
  _FLAT_DEL_HP=0.5, _FLAT_INS=1.25. Guard = hp_drift_margin in junction_refiner.py:730 (adds margin to any move
  whose new boundary lands in a >=4 HP run). Native re-aligner: penalty_table=None, motif_blind=True.
- HpPenaltyTable overrides ONLY del/ins (del_cost, ins_cost). NO sub_cost method exists — substitution is ALWAYS
  _FLAT_SUB=1.0 in _score_hp_anchored (line 874/912). So "table for substitution only" needs a NEW code path.
- Bundled Scer penalty_scores.tsv HAS op_type X (mismatch) + M rows, currently UNUSED by the re-placer.
- CRITICAL NUMBERS (Scer, base_class AT):
  * DEL empirical: hp1=0.44 hp2=0.42 hp3=0.29 hp4=0.17 hp5=0.12 hp6=0.08 hp7=0.06 hp8=0.03 ... hp12=0.02
    vs FLAT del_hp=0.5 (hp>=4). => EMPIRICAL DEL IN-RUN IS MUCH CHEAPER THAN FLAT. This IS the drift grease
    (arm-C's rejection mechanism). A table del_cost CANNOT be used without cheapening in-run deletions.
  * X (mismatch) empirical: hp1=1.00 hp2=0.89 hp3=0.71 hp4=0.64 hp5=0.72 hp6=0.80 hp7=0.93 hp8=1.24
    vs FLAT sub=1.0. Mismatch is a DIAGONAL move (consumes q+r) -> does NOT move a boundary into a run.
  * INS empirical: hp1=1.25 hp2=0.41 hp3=0.22 ... (also cheaper in-run; cheapening ins can also grease drift on
    the insertion side, so leave ins FLAT for the strict test).

## ANGLE REFINED — the ONE table use that provably avoids cheapening the drift axis:
(a) TABLE-SUB: use empirical X-stratum for substitution/mismatch cost ONLY; del + ins stay FLAT. Mismatch is a
    diagonal DP move; it cannot pull a boundary into an HP run. This is the strongest scoped table formulation.
(c) TABLE-SUB + GUARD ON vs FLAT + GUARD ON, on mix_fair_out + mix_r3b_out, recovery vs truth.
Also as a NEGATIVE control / completeness: (b) FULL TABLE (del+ins+sub) + guard — expected to over-shift, but
test whether guard neutralizes it enough to beat flat.
