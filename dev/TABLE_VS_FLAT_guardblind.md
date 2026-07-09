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

## Checkpoint: code + table values RE-VERIFIED (2026-07-08 resume)
- FLAT constants re-confirmed junction_scoring.py:113-116 (_FLAT_SUB=1.0, _FLAT_DEL_NORMAL=1.0, _FLAT_DEL_HP=0.5, _FLAT_INS=1.25).
- STR del path re-confirmed hp_penalty.py:420-426: run==1 & in di/tri STR -> penalty_table.str_del_cost(unit,copies); FLAT always del_normal=1.0.
- GUARD re-confirmed junction_refiner.py:462-490 + application 730-734: _hp_run_across returns 0 at any transition (b!=seq[pos-1]); STR boundary = transition => GUARD-BLIND. Guard only adds margin for moves INTO a mononuc run>=4.
- BUNDLED STR TABLE EXISTS: rectify/data/genomes/{saccharomyces_cerevisiae,homo_sapiens}/penalty_tables/str_penalty_scores.tsv.
- **KEY**: str_del loads the penalty_score COLUMN (hp_penalty.py:262) = the SAME reciprocal-rate heuristic (c/rate_mean, additive) that was REJECTED as arm-C for over-shifting into runs. Loader drops low_count & count<min_count & penalty<=0.
- STR del penalty_score values: AC repeats 0.14-0.78 (MUCH cheaper than flat 1.0), AAG/AAT 1.7-10 (capped). So at an AC/AT/AG dinuc-repeat boundary, TABLE makes STR-base deletions CHEAP -> same over-shift mechanism as arm-C. This is a STRIKE AGAINST the table on the guard-blind axis, not for it.
- IMPLICATION forming: the guard-blind STR axis is where the table has a DISTINCT cost, but that cost is the rejected reciprocal-rate heuristic which CHEAPENS in-repeat deletions => greases boundary drift into the STR (fabrication), the very failure the guard was built to stop but which the (mononuc-only) guard does NOT reach. So table on the guard-blind axis is plausibly HARMFUL, not helpful.

## Next: MEASURE. Build STR-boundary ground-truth cases, score truth vs drifted candidate under FLAT vs TABLE (str_del), guard ON. Check: does TABLE mis-rank (prefer drifted-into-STR) where FLAT+guard holds truth?

## Checkpoint: PROBE 1 RUN (str_probe.py) — first numbers on ground-truth-style STR geometry
STR table sanity (str_del_cost vs flat 1.0):
  AC/6=0.189, AC/4=0.390, AT/6=0.286, AG/6=0.223 (dinuc STR del MUCH cheaper than flat 1.0)
  AAG/4=5.20, AAT/5=8.75 (trinuc STR del MORE expensive — capped rare)
Geometry: (AC)6 STR at true acceptor je=68. GUARD-BLIND confirmed: _hp_run_across(je)=0, _hp_run_across(je+2)=0. _str_repeat_info(je)=('AC',7) => STR path ACTIVE.

DIRECTION A (FABRICATION, read matches truth):
  je_true=68  FLAT=0.000 TABLE=0.000
  je+2        FLAT=0.000 TABLE=0.000
  je+4        FLAT=0.000 TABLE=0.000
  => INSIDE an STR, a full-unit drift yields IDENTICAL downstream seq => ALL score 0 in BOTH models. Genuine sequence ambiguity; neither model breaks it. IMPORTANT: a TIE does NOT trigger a move (refiner requires strictly-better candidate to move; keeps incumbent). So FLAT does not fabricate here, and TABLE does not help. Fabrication direction: NEUTRAL, no table advantage.

DIRECTION B (RECOVERY, read lost 1 AC unit to contraction; TRUE junction je_true needs a 1-unit STR del):
  je_true (needs del) FLAT=2.000 TABLE=0.340
  je+2   (no del)     FLAT=0.000 TABLE=0.000
  => Under FLAT, WRONG drifted junction wins by 2.0. Under TABLE, drift STILL wins (0.0<0.340) but margin shrinks to 0.34. TABLE does NOT flip the winner => truth still loses. So even the pro-table recovery case does not make truth win. BUT margin compression is real (2.0 -> 0.34).
  => Caveat: in an STR, the contracted read is genuinely ambiguous — je+2-with-no-del IS an equally-valid parse of the observed bases. So "drift wins" is arguably CORRECT parsimony, not an error. Need to test a case where truth is NOT sequence-ambiguous.

## Next: (1) test whether ANY realistic geometry makes table FLIP truth-to-win (partial STR / STR abutting non-repeat so truth is unambiguous). (2) Confirm the tie=>no-move claim in refiner. (3) Test trinuc STR (AAG/AAT) where table del is EXPENSIVE (>flat) — does table there HURT recovery vs flat?

## Checkpoint: PROBE 2 RUN (str_probe2.py) — bounded-STR recovery, dinuc + trinuc, + HP control
Tie-in-move confirmed via code (junction_refiner.py:701-710): motif_blind tuple = (score_cmp, is_alt, ...); incumbent is_alt=0 so a TIE keeps incumbent (no move). Verified.

Bounded-STR recovery (read lost 1 unit; truth=STR-start needs 1-unit del; drift=+1unit):
  dinuc AC/6 distinct tail : FLAT truth=2.00 drift=0.00->drift ; TABLE truth=0.34 drift=0.00->drift ; NO FLIP
  dinuc AC/6 short tail     : same, NO FLIP
  dinuc AT/6                : FLAT truth=2.00 drift=0.00 ; TABLE truth=0.572 drift=0.00 ; NO FLIP
  trinuc AAG/4 (del=5.2)    : FLAT truth=3.00 drift=0.00 ; TABLE truth=2.267 drift=0.00 ; NO FLIP
  trinuc AAT/5 (del=8.75)   : FLAT truth=3.00 drift=0.00 ; TABLE truth=2.267 drift=0.00 ; NO FLIP
HP control (polyA-6, guard DOES fire): guard(je_drift)=6>0 -> margin (not cost) resolves.

KEY MECHANISTIC FINDING: an STR-contracted read is FUNDAMENTALLY AMBIGUOUS. Because the STR is periodic, skipping one ref unit (drift) leaves exactly n-1 units that align PERFECTLY (score 0.0) to the contracted read's remaining STR+tail. The tail does NOT frame-shift into a mismatch (periodicity). So drift wins 0.0 in BOTH models; TABLE only COMPRESSES truth's penalty (2.0->0.34 dinuc; 3.0->2.27 trinuc) but NEVER FLIPS the winner. => Even the advisor's strongest pro-table geometry does NOT make truth win. The ambiguity is intrinsic to STR contraction, unbreakable by any per-base cost.
Note trinuc: TABLE del (5.2/8.75) is MORE than flat (1.0) yet truth still scores LOWER under table (2.27<3.0) — because the DP finds a cheaper mixed del/sub path, not the raw STR del. So even "expensive" trinuc table doesn't make table hurt here; it's dominated by drift=0 anyway.

## Next: the ONE remaining pro-table geometry to test — ANCHORED del: truth unambiguous because a downstream base breaks STR periodicity (drift forced into a real mismatch). Only there can cost model flip the winner. Build it.

## Checkpoint: PROBE 3 RUN (str_probe3.py) — anchored/periodicity-broken STR geometries
Case i (single-base del in (AC)4): je_true/je+1/je+2 ALL identical (FLAT=1.0, TABLE=0.277). Single-base del is placeable anywhere in the periodic STR at equal cost => junction position doesn't change the score => tie among candidates => no directed move. Table only scales magnitude (1.0->0.277), same across candidates.
Case ii (perfect read): je_true=0, je+2=0 (tie by periodicity), je-2 (ref gains unit, forces indel)=FLAT 2.0/TABLE 0.555 -> TRUTH WINS both.
Case iii (acceptor at STR-end, exon2=pure anchor): je_true=0, je+1=0, je-1=FLAT 1.0/TABLE 0.277, je-2=FLAT 2.0/TABLE 0.555 -> TRUTH WINS both.

UNIVERSAL FINDING (all 3 probes): in EVERY STR geometry the candidate RANKING (winner or tie) is IDENTICAL between FLAT and TABLE. The table only SCALES the magnitude of a losing/deleting candidate's penalty; it NEVER REORDERS which candidate is best. When truth-vs-drift genuinely differ (periodicity broken, je+-k), truth ALREADY wins under FLAT. When they don't differ, it's a real sequence tie neither model breaks. => The STR table NEVER rescues a case flat loses, and NEVER fabricates a case flat holds. On the guard-blind STR axis, TABLE and FLAT are RANK-EQUIVALENT for junction placement.

## Next: airtight confirmation — systematic flip-search sweep over many random STR units/copies/tails/error positions; count any case where FLAT-winner != TABLE-winner. Expect ZERO.
