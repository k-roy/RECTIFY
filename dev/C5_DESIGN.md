# C5 — panel-failure tail / FracMinHash fallback: gate result (DRAFT 2026-06-30)

C5's facet (OVERVIEW row C5) is a **FracMinHash containment fallback localizer** for
the panel-failure tail — reads with **no acceptable window at all**. It is EXPLICITLY
GATED behind a *measured depletion trigger*: build nothing until a measurement proves a
recoverable tail EXISTS at realistic error. This is that measurement (decoder-free;
Stage-2 only PROTOTYPES containment for the kill-gate — no production index code).

Provenance: recovered from the C5 fan-out agent's `c5_tail_measure.py` (the agent
stalled before running it; the director fixed a one-line Stage-2 bug — `candidates`→
`cand_meta` — and ran it inline).

## The load-bearing distinction (window-level, per OVERVIEW lines 88-92)
The discovery ceiling is at the **WINDOW** level: a read in the RIGHT window with a
mis-called junction/indel inside it is recoverable by **realignment** (C3/refiner) —
NOT C5. C5 (the localizer) earns its keep ONLY on reads with no acceptable window.
So the tail is split:
- **empty-union** — no panel member placed the read at all → C5.
- **wrong-window** (non-paralog) — placed in the WRONG locus → C5.
- **right-window-wrong-internal (rwwi)** → C3/refiner, NOT C5.
- **paralog wrong-window** → C4 (containment can't disambiguate near-identical copies).
- **zero-evidence** — noised beyond any k-mer signal to the true window → unrecoverable
  by ANY method (never fabricate a placement).

## RESULT — Stage 1 (the dep-commit gate), minimap2-alone upper bound, reps=60, 6210 reads/rate

| hot-read rate | empty | wrongW | rwwi | correct | **C5 tail** | recov | zeroEv | C5tail% |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.00 | 0 | 13 | 311 | 5886 | 0 | 0 | 0 | **0.00%** |
| 0.05 | 0 | 15 | 1913 | 4282 | 0 | 0 | 0 | **0.00%** |
| 0.10 | 2 | 21 | 2006 | 4181 | 2 | 2 | 0 | **0.03%** |
| 0.20 | 24 | 23 | 2021 | 4142 | 24 | 24 | 0 | **0.39%** |
| 0.35 | 1062 | 25 | 994 | 4129 | 1062 | 995 | 67 | 17.10% |
| 0.50 | 1523 | 23 | 544 | 4120 | 1523 | 1268 | 255 | 24.52% |

(`hot_frac=0.3`; the rate is the hot sub-population's per-base error. Result file:
`scripts/benchmark/c5_tail_measure_result.txt`.)

## VERDICT — C5 DEFERRED (no trigger at realistic error)

**At realistic error (rate ≤ 0.20) the C5-addressable tail is ~0–0.4%.** It only
becomes substantial (17–24%) at EXTREME injected hot-read rates (0.35–0.50 per-base
error), which is far beyond the measured RNA004 hot-tail (~12%, ~1-2% of reads;
SESSION-10). And the bulk of below-ceiling reads at every rate is **rwwi** (right
window, wrong internal — C3/refiner territory, 311→2021 reads) plus a small steady
**paralog wrong-window** (C4, ~13–25), NOT C5. So per the pre-committed gate:

> tail ~0 at realistic error → **DEFER** (no depletion trigger; record the fraction;
> do NOT build the FracMinHash index/fallback). This is the correct outcome of a
> *measured* trigger, not a refute — the dep-commit gate is doing its job.

The window-level split is itself a finding: **the panel-failure problem at realistic
error is overwhelmingly a within-window realignment problem (C3/refiner) and a paralog
problem (C4), not a localization problem (C5).** A single-aligner panel (minimap2) is
the upper bound on the tail; the real 5-aligner panel would shrink it further.

## Residue / next (gated on a real trigger)
- **Stage-2 positive control (in flight):** at the elevated-error regime where the
  tail IS large, does the prototype FracMinHash containment localize the recoverable
  tail to the true window better than a genome-scale random-window null (with
  `zeroEv_acc ~ 0` proving it's containment specificity, not the recoverability
  filter)? This is a mechanism control only — it does NOT change the DEFER, because
  Stage-1 shows no realistic-error trigger. (Director re-running after the one-line fix.)
- **The real trigger** would be a MEASURED panel-failure tail on real data (multi-
  aligner, real reads) — e.g. the ~12% unmapped/all-herd fraction the OVERVIEW cites
  as unverified. Gate C5 on that measurement, not on injected extreme error.
- No member/index code built (the gate-refuted/deferred idiom). Smoke unaffected (the
  probe adds no shared-file edits).
