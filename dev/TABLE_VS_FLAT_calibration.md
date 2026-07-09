# TABLE vs FLAT — are the hand-set flat costs even RIGHT?

**Angle:** The native re-aligner scores candidate junctions with FLAT hand-set costs
(sub=1.0, del_normal=1.0, del_hp=0.5, ins=1.25) — never fit. Question: on ground-truth
sim (known junctions), does flat MIS-PLACE a measurable fraction that calibrated costs
(derived from the empirical DRS `rate_mean`) would fix? Is del_hp=0.5 / ins=1.25 near an
accuracy optimum, or arbitrary?

## Plan
1. Read `junction_scoring.py` — confirm the flat constants + how sub/del/ins costs enter
   `_score_hp_anchored`. Identify which axis each cost touches (drift vs mismatch).
2. Read the empirical table (`penalty_scores.tsv` / `HpPenaltyTable`) — extract `rate_mean`
   per (motif, hp_len, error_kind). Derive calibrated -log-likelihood costs.
3. Locate ground truth: `scripts/benchmark/noncanon_sim/` panels (mix_fair_out, mix_r3b_out),
   `paired_arm_test.py` recovery-vs-truth. Understand the scoring harness.
4. Measure junction-placement ACCURACY under flat costs on ground truth (recovery / exact
   placement / fabrication rate), guard ON. This is the baseline.
5. Derive calibrated costs from rate_mean; A/B vs flat, guard on. Does calibrated place more
   accurately?
6. Sweep flat costs (del_hp, ins ratio) — is 1.25/0.5 near an optimum or arbitrary?
7. Verdict: is DROPPING the table justified from the accuracy angle?

## Checkpoints
(append one line per sub-step below)

## Checkpoint log
- [orient] Read junction_scoring.py: flat constants _FLAT_SUB=1.0 _FLAT_DEL_NORMAL=1.0 _FLAT_DEL_HP=0.5 _FLAT_INS=1.25. In DP: `sub` charges mismatch (diag), `del_normal/del_hp` charge deletion (left move, ref base skipped), `ins` charges insertion (above move, query base). del_hp applies when ref base is in HP run >= hp_min_run(=4). With penalty_table: del/ins costs looked up per-HP-length; sub UNCHANGED at 1.0 (table has no substitution row — confirm).
- [orient] Empirical table penalty_scores.tsv (132 rows): cols op_type/base_class/hp_length/rate_mean/count_total/penalty_score/low_count. Only D and I ops (NO substitution). penalty_score = the reciprocal-rate heuristic (arm-C's over-shift column). rate_mean = the MEASURED per-base error rate — this is what I calibrate -logLR from.
- [orient] Ground truth: mix_fair_out + mix_r3b_out. read_truth.tsv (5800 rows fair) has true_donor/true_acceptor/motif_rung/context/has_true_junction. paired_arm_test.py scores AMBIGUITY-AWARE recovery (normalize_junction slide within HP/repeat) of refined junction vs truth. Scored cells: motif_rung R3 (cryptic non-canon) + R0flank (canonical flanking-HP).
- [orient] Existing arms: arm_A=motif_blind False (incumbent+table), arm_B=motif_blind flat table-free (SHIPPED re-aligner, NO guard), arm_C=motif_blind+empirical table (NO guard), arm_E=arm_B+hp_drift_margin (guard on, table-free = SHIPPED config). arm_E_m{0p5..4p0}=guard sweep. So table-vs-flat WITH guard needs a motif_blind+table+guard arm (arm_ins_* / arm_mbF_* may be it — check).
