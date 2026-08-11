# Flat-cost optimality sweep — the PI's "are 0.5/1.25 the right costs?" (2026-07-08)

Guard-OFF (raw cost effect), mix_fair_out, recovery vs truth (n=400/cell). Shipped hand-set = del_hp=0.5, ins=1.25.

## del_hp axis (ins=1.25 fixed) — PERFECTLY FLAT
  del_hp  ACC_A_D0  ACC_A_D2  BOT_A_D1
  0.25    0.307     0.930     0.757
  0.50    0.307     0.930     0.757   <- shipped
  0.75    0.307     0.930     0.757
  1.00    0.307     0.930     0.757
=> a 4x range in the homopolymer-deletion cost changes placement accuracy by ZERO reads.

## ins axis (del_hp=0.5 fixed) — NEARLY FLAT
  ins     ACC_A_D0  ACC_A_D2  BOT_A_D1
  1.00    0.307     0.920     0.757
  1.25    0.307     0.930     0.757   <- shipped
  2.00    0.312     0.935     0.767
=> across ins 1.0->2.0 (2x), recovery moves 1-3 reads/400 (~1%); ins=2.0 marginally best, within noise.

## VERDICT (closes the PI's Q3)
Placement accuracy is a BROAD FLAT PLATEAU around the shipped (0.5, 1.25). The exact cost MAGNITUDES are
essentially irrelevant to placement — the cost ORDERING (del_hp < del_normal=1.0 < ins; HP-dropout cheap,
insertions rarer) does the work, and the guard handles the drift. So the hand-set 0.5/1.25 are validated as a
MEASURED plateau, not a fragile lucky guess: even the extreme ins=2.0 differs by ~1-2 reads. Combined with the
3-angle adversarial result (the empirical table does NOT beat flat: no sub row, HP-axis=arm-C grease owned by the
guard, STR rank-equivalent+net-negative+near-empty), the complete answer: DROP the table, KEEP flat — both the
choice-to-drop and the specific magnitudes are now measured, not asserted.
CAVEAT (honest): guard-OFF, one canonical-HP panel; establishes cost-magnitude INSENSITIVITY, not global
optimality across all data. Sufficient given the ordering rationale + table-doesn't-beat-flat. (ins=1.5 point not
run; trend between 1.25 and 2.0 is monotone-flat.)
