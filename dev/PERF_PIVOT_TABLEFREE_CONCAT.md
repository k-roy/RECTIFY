# PERF PIVOT — concat-DP is byte-identical for the SHIPPING (table-free) re-placer (2026-07-08)

## The realization (verified H1/H2/H3)
The full-run ins_cost re-validation may be addressing a config we DON'T ship. Chain:
- H1: the native re-aligner (shipping = motif-blind + guard) AND the measured-slow workloads (Sumner
  genome-wide sumner_gw_discover.py, human/yeast transfer real_drs_hp_drift.py) all call refine with
  penalty_table_path=None. The empirical table was REJECTED from the re-placer (del_open/arm-C conclusion).
- H2: with penalty_table=None, _score_hp_anchored ins_costs = [ins]*Q (flat 1.25) — per-cut / full-run / refcol
  are ALL IDENTICAL. The full-run flag is a NO-OP in the table-free config. (Confirmed: item-5 margin re-sweep
  recovery curve came back IDENTICAL to flag-off — because arm-B/arm-E use no table.)
- H3 (decisive): the concat-DP audit's OWN control — for penalty_table=None, concat with the MECH1 boundary-column
  exclusion is BYTE-IDENTICAL, 0/20000 diverged (dev/CONCAT_DP_AUDIT_end-to-end-integration.md). Costs
  (ins=1.25, del=1.0/0.5, sub=1.0) are exactly float-representable -> no FP tolerance needed. MECH2 (per-cut ins
  truncation) — the ONLY divergence source — DOES NOT EXIST when ins is flat.

## Implication
The ~30x concat-DP speedup is available NOW, BYTE-IDENTICALLY, for the config we ship (motif-blind + guard,
table-free) and the config that's slow (Sumner genome-wide, human transfer). NO ins_cost change needed.
Full-run ins_cost only matters IF a penalty table is passed to the re-placer — which (a) del_open/arm-C rejected
and (b) is not what's slow.

## Consequence for the re-validation checklist
- Items 4/5/6/7/8 (full-run re-validation) are MOOT for the shipping table-free re-placer (flag is a no-op there).
  They'd only matter for a table-in-re-placer config we don't ship.
- The running full-run transfers (sbatch 33182801) are a NO-OP (real_drs_hp_drift uses table=None) -> will
  reproduce flag-off exactly (0.9884 / 0.7914). Confirming control at best; candidate to cancel.
- CORRECT perf path: implement the concat-DP for the table-free path (MECH1 boundary-col exclusion), byte-identical,
  ~30x. Meticulous build + triple audit (byte-identity 0/20000 is the gate).

## Caveat to verify
The PRODUCTION mainline (correct_command, master) reportedly uses table ON + motif_blind OFF (incumbent arm-A) —
a DIFFERENT code path. If a future deployment puts the table IN the re-placer, full-run would re-matter there. But
the native re-aligner this session built (and everything measured slow) is table-free.

## ADVISOR-CORRECTED (2026-07-08) — three fixes to the framing above
1. H2 was ASSERTED, now COMPUTED: per-cut vs full-run vs refcol _score_junction, penalty_table=None, 4000 random
   crossing-junction cases with injected HP runs -> 4000/4000 EXACT-identical for BOTH variants. Flag is a genuine
   no-op without a table (scalar-score level). The running transfers (33182801) are the END-TO-END confirmation —
   KEEP them (do NOT cancel): flag-on must reproduce flag-off EXACTLY (yeast 0.9884 / human 0.7914); if not, the
   pivot premise is broken. That is the cheap empirical test of the linchpin.
2. "Available NOW byte-identically" OVERSTATES: the 0/20000 None-path result was an auditor's SCALAR-score
   reconstruction; the concat-DP was NEVER built (every builder stalled). Real remaining work = BUILD it + verify
   byte-identical END-TO-END on the None path (same chosen (js,je), same CIGAR, same emitted BAM) with the MECH1
   boundary-column exclusion correctly implemented. Not "confirm identical" (that's the harness) — "build + verify."
3. SCOPE SMALLER: on the None path MECH2 (per-cut ins) + MECH3 (FP tolerance) both VANISH (flat, dyadic ins);
   only MECH1 remains. Far narrower than the general penalty-table concat. Proportionate = focused build + a
   byte-identity harness (20k None-path cases, scalar AND end-to-end BAM diff) + ONE adversarial pass on the
   boundary handling. NOT another 5-agent meticulous-build/triple-audit.

## PREMISE TO STATE EXPLICITLY TO THE PI (what makes this valid)
The whole pivot rests on: THE DEPLOYED NATIVE RE-ALIGNER IS TABLE-FREE (motif-blind + guard, penalty_table=None).
True for everything measured slow (Sumner genome-wide, human transfer) and for the shipping arm-E. IF a future
deployment ever puts the empirical table back INTO the re-placer, full-run stops being moot there and the shelved
investigation (dev/INSCOST_*.md) is the resumption point. Production mainline (correct_command/master) reportedly
runs table-ON + motif-blind-OFF (incumbent arm-A) — a SEPARATE path, not the native re-aligner.

## BUILD RESULT (2026-07-08) — table-free vectorized DP is BYTE-IDENTICAL + 14.3x
Prototype dev/concat_dp_prototype.py (inline build, harness-gated). Design: _score_junction's min_k[t1(k)+t2(k)]
computed via TWO vectorized DP passes instead of ~2L(<=60) per-k calls:
- all_t1(rescue, ref, del_costs): all t1(k)=score(rescue[k:], ref) in ONE reversed-prefix DP (query-suffix family
  via the reversal: t1(k) = D[L-k][R] over reverse(rescue) vs reverse(ref)). VERIFIED 0/3000 vs per-k reference.
- t2(k) = score(reverse(rescue)[L-k:], ref_intron_end_rev) = the SAME query-suffix shape -> reuse all_t1 on
  reverse(rescue). best = min_k(t1[k]+t2[k]), same early-exit, k in [0,L) (MECH1 boundary handled — no k=L col).
VERIFIED: score_junction_fast vs the REAL _score_junction, table=None, 8000/8000 EXACT (mismatches=0).
SPEEDUP: 400 candidate scorings 26.26ms/call -> 1.83ms/call = 14.3x; 2 DP passes/candidate vs ~60.
REMAINING: port behind a flag (gate on penalty_table is None), byte-identity harness in-repo + tests, END-TO-END
BAM diff on a sim panel, 1 adversarial audit. Then cluster deploy + re-run Sumner genome-wide/human at ~14x.

## PORT VERIFIED — 3 independent byte-identity proofs (2026-07-08)
Ported into _score_junction (flag _USE_CONCAT_DP, gated penalty_table is None):
- Standalone prototype: 8000/8000 exact scalar.
- In-repo 20k harness (default vs flag-on, table=None): 20000/20000 EXACT, 0 mismatches.
- END-TO-END BAM diff (refine mix_fair_out flag off vs on): 0/5800 reads differ (same reference_start + CIGAR).
- Suites: 70 passed flag-ON / 46 flag-OFF (refiner+guard+slide+false-junction+anchor-gate).
- penalty-table path flag ON vs OFF: unchanged (fast path strictly gated on penalty_table is None).
14.3x (26.26 -> 1.83 ms/call), 2 DP passes vs ~60. REMAINING: 1 adversarial audit (boundary/edge cases) -> flip
default + cluster deploy + re-run Sumner genome-wide/human at 14x.

## ADVERSARIAL AUDIT — HELD; exact-by-construction (2026-07-08)
Single adversarial auditor (dev/CONCAT_DP_TABLEFREE_AUDIT.md): NO divergence. Attacked the del_costs-reversal
weak spot the random harness misses (HP runs abutting intron_end, del_hp=0.5 entering the score; mid-k winners
exercising t2_vec[L-k]) with 10k targeted cases + all edge classes (L=1, R=0 both flanks, truncated windows,
N/lowercase, all-HP) — all EXACT. Gating CORRECT (flag ON + penalty_table present == legacy; only fires when
RECTIFY_CONCAT_DP==1 AND table is None -> cannot corrupt production table path; per-worker flag mismatch is safe
since both paths identical). FP: all costs dyadic -> every partial sum exactly representable -> 0 ULP GUARANTEED
(exact, not merely observed). Fixed the one flagged hazard: factored _FLAT_INS/_FLAT_SUB/_FLAT_DEL_* into shared
module constants (both fns reference them) so a future retune can't silently break identity. Re-verified post-
refactor 5000/5000 exact + suites green. => 4 byte-identity proofs (8k/20k/5800-e2e/10k-adversarial) + audit HELD.
READY to flip default + cluster-deploy (PI review).

## PIVOT PREMISE CONFIRMED END-TO-END (2026-07-08, sbatch 33182801)
Full-run flag-ON transfers reproduce flag-OFF EXACTLY on real data (table-free config): yeast Bguard 0.9884 ==
0.9884, human Bguard 0.7914 == 0.7914. Confirms H2 (full-run is a genuine no-op without a penalty table) at the
real-data pipeline level, not just scalar. Full-run correctly shelved; concat-DP is the right perf answer for the
table-free re-placer.

## CLUSTER CONFIRMATION (2026-07-09, sbatch 33235941) — byte-identical + faster on REAL data
Real Sumner SMA_7.12 @2%, concat ON (default) vs OFF: every summary count IDENTICAL (raw_annot_canon 95482,
ref_novel_noncanon 6990, revealed_novel_noncanon 1788, same reads md5) — the only "diff" was the cosmetic sample-
name line. Byte-identity HOLDS end-to-end on real human data. Wall-clock 4:29 (off) -> 1:14 (on) = 3.6x on this
tiny 2% job (fixed overhead: 3.2GB genome load + pool build dominates); refine-only carries the full 14x, so a
full-depth run (refine-dominated) approaches 14x. Concat-DP PROMOTED + validated on cluster.
