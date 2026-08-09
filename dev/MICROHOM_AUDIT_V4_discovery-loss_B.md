# MICROHOM AUDIT V4 — discovery-loss — AUDITOR B

**Status: IN PROGRESS**
**Task:** discovery-loss lens. Measure the empirical rate at which REAL cryptic 3'SS/5'SS
reads (genuinely discoverable: refiner-with-guard-OFF moves incumbent to true site, delta_improve>0)
are WRONGLY VETOED because their boundary trips microhom_drift (mh>=0.5, short direct repeat).
Cross-checked against FAB (fabricated drift) suppression on the same grid.

**Commit under audit:** 05664bc (drift_near_tie_cap; `_effective_veto_margin`).
**Working dir (READ-ONLY):** /Users/kevinroy/work/rectify/.claude/worktrees/agent-a25a2c1e784ad37dc
**Python:** /Users/kevinroy/miniconda3/bin/python

## Checkpoint log
- [init] Durable record created. Scratch dir made. About to read guard mechanics in
  rectify/core/splice/junction_refiner.py to verify the mechanism myself before building harness.
