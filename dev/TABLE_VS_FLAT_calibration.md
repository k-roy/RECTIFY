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

## CHECKPOINT — code facts VERIFIED (2026-07-08 resume) + refined plan (advisor)
- VERIFIED flat constants junction_scoring.py:113-116: _FLAT_SUB=1.0 _FLAT_DEL_NORMAL=1.0 _FLAT_DEL_HP=0.5 _FLAT_INS=1.25.
- VERIFIED sub DEAD END: HpPenaltyTable.from_tsv (hp_penalty.py:217-244) reads ONLY op in {D,I} into del/ins tables;
  X (mismatch) rows in penalty_scores.tsv are DROPPED at load. No sub_cost method. Substitution is a DIAGONAL DP move
  (consumes q+r) -> placement-neutral to first order. => "table for mismatch context" cannot even be expressed. 1 line, done.
- VERIFIED STR path (the real guard-blind candidate): _precompute_del_costs (hp_penalty.py:420-426) calls
  penalty_table.str_del_cost(unit,copies) ONLY when penalty_table is not None AND run==1 AND _str_repeat_info fires
  (di/tri STR). Native re-aligner runs penalty_table=None => STR path NEVER reached; flat charges del_normal=1.0 there.
  Guard (_hp_run_across, homopolymer run>=4) does NOT fire at STR boundaries (run==1). So STR is the ONE axis table
  expresses that flat cannot AND guard cannot reach.
- FLAG (cuts AGAINST table): Scer str_del costs for common dinucs are CHEAPER than flat 1.0: AC 3cp=0.78 6cp=0.19
  11cp=0.14; AT 3cp=0.81 -> 0.14; AG 3cp=0.55 -> 0.08. Same in-run cheapening that sank arm-C, now on guard-blind axis.
  => table-STR at least as likely to FABRICATE (over-shift into repeat) as to fix. Must measure FABRICATION not just recovery.

## REFINED PLAN (advisor, cheapest discriminator first)
1. FLAT-RATIO SWEEP (core angle, NO refine): candidate pool is cost-independent (built from BAM + canonical priors).
   Score-level probe: existing panel truth + candidate pools, hold del_normal=sub=1, sweep del_hp x ins, guard margin
   applied, measure EXACT-PLACEMENT accuracy. Flat surface near (0.5,1.25) => hand-set ok => drop justified.
   Peak away => under-justified. This fills `measurement` + answers the angle directly.
2. STR FREQUENCY SCAN (annotation-only, zero compute): fraction of annotated junctions (.junc.bed / real) whose
   boundary sits in a di/tri STR context. Near-zero => niche empty => drop justified.
3. If STR freq nonzero: tiny targeted STR probe, measure FABRICATION + recovery table-STR+guard vs flat+guard.

## MEASUREMENT 1 — STR-boundary FREQUENCY scan (real Scer annotated introns) [str_freq_scan.py]
Zero-refine annotation scan of all 346 resolvable annotated S. cerevisiae introns (378 total, chrmt unresolved).
Boundary = donor pos (intron_start) and acceptor pos (intron_end-1), the bases whose deletion the DP charges.
- STR-boundary (di/tri repeat >=3 copies, _str_repeat_info): donor 1/346, acceptor 0/346, EITHER 1/346 = **0.29%**.
  The ONE case: chrV donor 307746 in an AAG x3 triplet.
- HP-boundary (mono run >=4, the GUARD's domain): EITHER 13/346 = 3.76% (10x larger, AND guard-covered).
VERDICT on STR lead: the table's guard-blind STR niche is essentially EMPTY (0.29%) on real yeast junctions. Even
if the table placed all STR-boundary junctions perfectly and flat mis-placed all of them, the ceiling benefit is
0.29% of junctions. And the str_del costs are CHEAPER than flat (AAG 3cp=2.47 is actually >1 here, but AC/AT dinucs
0.14-0.78 << 1) -> fabrication risk on that tiny niche. STR does NOT rescue the table.
