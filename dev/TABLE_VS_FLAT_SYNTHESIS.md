# TABLE vs FLAT — synthesis (2026-07-08): DROP the table from the re-placer is JUSTIFIED

VERDICT: dropping_table_justified=TRUE; flat_costs_defensible=TRUE. Grounded in the 3 angles' MEASURED durable
records (agents killed by 2 API outages but persisted the numbers; durability policy).

## Why (3 orthogonal angles converge)
1. SUBSTITUTION unreachable: the table has ONLY D and I ops (no sub row); DP hardcodes cost_sub=_FLAT_SUB=1.0 in
   BOTH modes. Mismatch is a diagonal move (consumes q+r) → placement-neutral. "Table for mismatch context" is a
   DEAD END — cannot even be expressed.
2. HP-drift axis owned by the TUNED GUARD (hp_drift_margin=3.0, task #16), not del_hp's magnitude. Empirical in-run
   del is MUCH cheaper than flat (hp8=0.03 vs flat 0.5) → pure drift-grease → the table HURTS here (= arm-C).
3. STR-del (the ONE distinct guard-blind table axis) fails on 3 stacked grounds:
   - RANK-EQUIVALENT to flat (measured): inside an STR a full-unit drift gives identical downstream seq → score 0
     in BOTH models (genuine ambiguity, a TIE → refiner needs strictly-better to move, so no fabrication). Where
     truth differs from drift, truth already wins under flat; the table only COMPRESSES the winning margin (e.g.
     2.0→0.34), NEVER flips the argmin winner. So it changes magnitudes, not placements.
   - NET-NEGATIVE prevalence: str_del is BIMODAL (67/89=75% CHEAPER than flat = grease; 22/89=25% forbid). The
     grease direction fires ~2.3x more than forbid. And str_del loads the penalty_score COLUMN = the SAME
     reciprocal-rate heuristic REJECTED as arm-C → cheapens in-repeat del → greases drift into STR (fabrication).
   - NEAR-EMPTY niche: STR boundary on 346 real Scer introns = 1/346 = 0.29% (HP boundary, the guard's domain,
     = 3.76%). The table's only distinct axis barely exists on real junctions.

## RECOMMENDATION: KEEP FLAT + guard; DROP the table from the native re-placer. Already validated on real DRS
(task #17, 17 fix/0 harm) + human ONT vs GENCODE (task #18, 153 fix/0 harm). Flat costs are ORDINALLY justified
(del_hp=0.5 < del_normal=1.0 < ins=1.25; Nanopore HP dropout cheap, insertions rarer).

## HONEST GAP (the PI's "are the costs even RIGHT"): the flat-RATIO optimality sweep was NOT executed. Verdict
"drop the table" does not depend on it, but "flat costs are OPTIMAL (not just ordinally ok)" is UNMEASURED.
FOLLOW-UP (cheap, score-level, no refine): hold del_normal=sub=1.0, sweep del_hp x ins with guard ON, measure
exact-placement accuracy on mix_fair_out + mix_r3b_out. Flat plateau near (0.5,1.25) → hand-set OK → upgrades
"ordinally defensible" to "measured plateau." Peak away → flat under-tuned (still not a table argument).
