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

## CHECKPOINT 2 (resumed run) — measurement design locked
- STR path confirmed as the ONE guard-blind table axis: hp_penalty.py:420-426, str_del_cost() uses
  str_penalty_scores.tsv (di/tri repeat, >=3 copies). MANY STR del penalties < 1.0 (AC:0.14-0.78,
  AT:0.15-0.81, AG:0.08-0.68) => STR table ALSO cheapens in-run deletions (same drift-grease mechanism
  as arm-C/HP-del), but guard-blind (guard is HP-only, _hp_run_across).
- refine_bam_junctions takes penalty_table_path + str_penalty_table_path + motif_blind + hp_drift_margin.
  So TABLE+guard arm = pass both table paths + motif_blind=True + hp_drift_margin=3.0. FLAT+guard = arm_E.
- PANEL PROBE: mix_fair_out + mix_r3b_out each have only 14 unique true junctions; STR-at-acceptor=0,
  STR-at-donor=1, HP4+-at-acceptor=4/1. => existing panels barely exercise STR path. Full-panel A/B will
  show table effect dominated by HP-del path (which IS guarded). STR path needs a TARGETED synthetic probe.
- Substitution path: DEAD (no sub row in either table; _FLAT_SUB=1.0 always). Confirmed.
- PLAN: (1) full-panel A/B TABLE+guard vs FLAT+guard on both panels, recovery vs truth (reuse bams, 2 lean
  refines/panel, n_workers=4). (2) STR-path direct probe: synthetic STR-boundary junction, score truth-place
  vs drifted-into-STR under FLAT vs TABLE(str) — does cheaper STR-del flip the argmin to the drifted place?

## CHECKPOINT 3 (resumed) — cost-direction + prevalence (advisor-corrected: STR table is BIMODAL)
COST DIRECTION (measured on bundled Scer tables):
- HP-del: UNIFORMLY cheaper than flat. run>=4 min=0.019 max=0.350 (flat del_hp=0.5). base A:0.44->0.02,
  base C:0.85->0.03 as run grows. => pure drift-grease; guard is the real fix on HP axis. SUPPORTS drop.
- STR-del: BIMODAL. 89 entries: 67 (75%) CHEAPER than flat (<1.0, grease); 22 (25%) MORE EXPENSIVE (>=1.0,
  forbid). Most-expensive (rare-del, FOR-table): AAC3=10, AAG5=10, ACC3=10, GGT3=10, GTT3=10, AAT3=8.8,
  CTT4=8.5, ATT3=6.7. Cheapest (grease): AG7=0.08, AGT7=0.12, AC11=0.14. Flat 1.0 cannot express EITHER;
  guard (HP-only) reaches NEITHER. => STR rare-del is the ONE place the table has a specificity capability
  neither flat nor guard has.
PREVALENCE (370 real annotated Scer introns, boundary +/-10bp):
- HP4+ near boundary: 273 (73.8%) [guard covers]
- STR near boundary:  43 (11.6%) [guard-BLIND]
    RARE-del STR (>=1.0, FOR-table forbid-drift): 13 (3.5%)   <- the table's helpful niche
    COMMON-del STR (<1.0, table GREASES drift):   30 (8.1%)   <- table LIABILITY, 2.3x more common
=> The table's helpful direction (rare-del STR) fires at 3.5% of boundaries; its harmful direction (common-del
   STR, greases drift like arm-C) fires 2.3x MORE often (8.1%). A whole-table del cost is NET-NEGATIVE on
   STR axis. Only a STR-del table RESTRICTED to the rare-del/expensive direction could help, at 3.5% ceiling.

## CHECKPOINT 4 (resumed) — two-direction STR drift probe (real scoring primitives, _precompute_del_costs)
Cost a 1-unit boundary DRIFT into an STR run must pay (this is the del the DP charges to move the boundary):
  RARE-del STR  AACx4: FLAT=3.00  TABLE=10.84  -> TABLE FORBIDS drift (WIN)
  RARE-del STR  GGTx3: FLAT=3.00  TABLE=11.69  -> TABLE FORBIDS drift (WIN)
  RARE-del STR  ATTx3: FLAT=3.00  TABLE= 7.51  -> TABLE FORBIDS drift (WIN)
  COMMON-del AC x6:    FLAT=2.00  TABLE= 0.38  -> TABLE GREASES drift (LOSE)
  COMMON-del AT x9:    FLAT=2.00  TABLE= 0.35  -> TABLE GREASES drift (LOSE)
  COMMON-del AG x5:    FLAT=2.00  TABLE= 0.45  -> TABLE GREASES drift (LOSE)
  HP-del contrast (A6): FLAT(del_hp)=0.50 TABLE=0.08 -- guard vetoes move regardless (guarded axis).
=> Confirmed mechanism BOTH directions. Table's ONLY win is rare-del STR (forbid drift), guard-blind, 3.5% prev.
   Table LOSES on common-del STR (grease, 8.1% prev) exactly like arm-C. A whole STR-del table is NET-NEGATIVE.
