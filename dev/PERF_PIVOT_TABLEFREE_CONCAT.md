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
