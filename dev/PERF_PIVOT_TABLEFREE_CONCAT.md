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
