# TABLE vs FLAT (guard-blind axis) — does dropping the empirical DRS penalty table over-generalize from the HP-drift result?

## Angle
The HP-drift guard fires ONLY on moves INTO a homopolymer run (`_hp_run_across`).
The claim "guard fixes the drift so no cost model matters" is only as general as the guard.
Find junction-placement errors OUTSIDE the guard's domain (non-HP boundary ambiguity,
substitution-heavy, mismatch context near junction, STR/di-tri repeat boundaries,
adjacent-run) where a calibrated/table cost helps and flat costs fail. Measure on ground truth.

## Plan
1. Orient: read junction_scoring.py (_score_junction, _score_hp_anchored, flat constants, hp_drift guard).
2. Read HpPenaltyTable + penalty_scores.tsv — what dimensions does the table actually encode? (sub context? del context? or only HP?)
3. Read paired_arm_test.py + sim panel layout (mix_fair_out, mix_r3b_out) — the truth format + recovery metric.
4. Identify the guard-blind strata in the panels (non-HP boundary ambiguity, STR, substitution-heavy).
5. Build/probe: score truth vs mis-placed candidates under FLAT and under TABLE, guard ON, on those strata.
6. Measure accuracy/recovery/fabrication flat vs table, guard on, guard-blind strata.
7. Verdict: is dropping the table justified? strongest case FOR table (if any)? residual risk?

## Checkpoints
- [init] file created, plan written.

## Checkpoint: orientation done (code read)
- FLAT costs: _FLAT_SUB=1.0, _FLAT_DEL_NORMAL=1.0, _FLAT_DEL_HP=0.5 (run>=4), _FLAT_INS=1.25. junction_scoring.py:113-116.
- TABLE overrides ONLY del_cost + ins_cost (per HP-run-length, base-class AT/CG). SUBSTITUTION cost is FLAT 1.0 in BOTH modes — the DP hardcodes cost_sub=sub=_FLAT_SUB (junction_scoring.py:874,912,967; numba :72). => TABLE CANNOT help pure-substitution/mismatch context. Negative result to nail.
- TABLE ALSO consults STR sub-table: _precompute_del_costs (hp_penalty.py:420-426) — when run==1 (NOT a homopolymer) AND position in a di/tri STR, uses penalty_table.str_del_cost(unit,copies) instead of flat del_normal=1.0. FLAT mode always charges del_normal=1.0 there. => STR-boundary is the axis where TABLE has a real, distinct, calibrated cost the guard cannot reach AND flat cannot express.
- GUARD (_hp_run_across, junction_refiner.py:462): returns 0 (guard does NOT fire) whenever seq[pos-1]!=seq[pos] — i.e. at ANY sequence TRANSITION. STR boundaries (...AT|AT... boundary base != neighbor), adjacent-run boundaries (...AAAA|CCCC...), non-canonical acceptors at a real transition => all GUARD-BLIND.
- Panel strata: mix_fair_out = flanking-HP A-runs (ACC/DON/BOT x D0..D8) + base variants ACC_{C,G,T}_D3 (non-A homopolymers) + INTRONFREE. mix_r3b_out = HP / plain / INTRONFREE (R0/R1/R3/WT). NEITHER panel has a di/tri STR-boundary stratum. => must BUILD STR + adjacent-run cases to test the guard-blind axis where table's STR path lives.
- Refine reuse: _make_arm_e.py / refine_bam_junctions reuse aligned.sorted.bam; _make_ins_arms* show ins-cost axis probes already done (percut/fullrun/refcol).

## Next: call advisor before building STR/adjacent-run guard-blind cases.
