# MICROHOM-DRIFT GUARD + NEAR-TIE CAP (05664bc) — V4 REDUNDANT-AUDIT SYNTHESIS

**Date:** 2026-07-13
**Synthesizer:** subagent (agent-a25a2c1e784ad37dc), READ-ONLY, worktree `worktree-agent-a25a2c1e784ad37dc`.
**Question:** Does the redundant (2×2) V4 audit clear the microhomology-drift guard + near-tie cap so the
default `microhom_drift_margin` can be flipped from 0.0 (inert) to enabled?
**Consensus rule (adversarial):** `all_clear=true` ONLY if byte-identity cleared (V2) AND **both**
discovery-loss auditors reach CLEAR with consistent numbers AND **both** detector auditors reach CLEAR. A HOLD
or a real fault from **either** auditor on a lens blocks that lens. Unresolved A/B disagreement = NOT-clear (a
split verdict cannot license a flip).

---

# ============ BOTTOM LINE ============

## `all_clear = FALSE` → **HOLD** (do not flip the default; ship as code, default OFF)

The verdict is **overdetermined** by two independent, individually-sufficient blockers, neither of which
depends on the one open A/B disagreement (A8):

1. **DETECTOR-CORRECTNESS is NOT clear — A5 is a CONSENSUS holding fault.** Both detector auditors
   independently confirmed, with byte-consistent numbers, that `_frac_match` scores `N==N` (any `x==x`, incl.
   `n`, `X`, lowercase) as a MATCH → phantom microhomology ≥ 0.5 → **false-positive veto of a genuinely
   read-supported move** (`delta_improve = 6.0` edit units). A CLEARED-value detector (N excluded) returns
   < 0.5 and would SPARE the move. Confirmed at unit level AND end-to-end decision-flip by BOTH auditors.

2. **DISCOVERY-LOSS reached NO verdict from EITHER auditor.** Auditor B produced **zero** persisted work
   (empty scratch dir, no `.md`). Auditor A stalled at CP2: `harness.py` + `sweep.py` written, one decisive
   read confirmed, but **no results TSV ever persisted and no verdict written**. The gate requires BOTH
   discovery-loss auditors to reach CLEAR with consistent numbers; neither completed. **Audit-incomplete fails
   the gate by construction** (same principle the V1 synthesis applied) — it is not "all clear," it is
   "not audited."

Byte-identity (V2) is the ONLY lens that is cleared. One cleared lens out of three does not clear the flip.

---

# ============ PER-LENS CONSENSUS ============

## LENS 1 — BYTE-IDENTITY (V2, prerequisite) — **CLEARED** (single completed auditor, robust)
Not re-litigated here (per brief). V2 `dev/MICROHOM_AUDIT_V2_byte-identity.md` established, with independent
sequential + parallel BAM diffs HEAD vs the pre-cap parent `05664bc^`:
- `junction_refiner.py` is git-IDENTICAL between `17e57ba` (prior cleared guard state) and `05664bc^` → the
  cap commit is a PURE additive delta on an already-inert file; only removed code line is the veto compare,
  replaced by the `_effective_veto_margin(...)` variant.
- Sequential (n_workers=1, RAW order, SHA256 over full SAM text): all 9 cells MATCH across
  mix_r3b / mix_fair / fab_sweep × {default, motifblind, **guard_on_cap_off**}. The `guard_on_cap_off`
  config genuinely ENTERS the veto block (refined counts drop 815→308, 585→188, 152→92) yet stays
  byte-identical → proves `_effective_veto_margin(hold, eff, cap=0.0) == eff_margin` exactly on the live path.
- Parallel (n_workers=3, position-sorted): all 9 cells MATCH + HEAD-self determinism.
- Structural honesty of the cap docstring: HONEST. The cap shares ONE score axis with `delta_improve`
  (P1), only ever LOWERS the veto margin (P2, 0 violations in 200k+grid), floors at `hold` (P3), equals
  `clamp(cap, hold, eff)` on the active branch (P4), is identity when disabled (P5), monotone in cap (P6).
  Re-admission is on `delta_improve` ALONE → real cryptics and fab drift with the same delta are treated
  IDENTICALLY inside `(0, cap)` → the cap BOUNDS but does not close, and adds **no discriminating signal**.
- **Ergonomic flag (P7, not a fault):** discontinuity at cap=0 — `cap=0` is OFF (full veto) but any tiny
  `cap=0+ε` gives near-MAX guard-weakening (`max(hold, ε)`). "Small cap" ≠ "small effect." Document before
  exposing the knob to tuners.

**Verdict: CLEARED on inertness + honest-framing.** Does NOT endorse the flip — the operating-point efficacy
is out of this lens.

---

## LENS 2 — DISCOVERY-LOSS — **NOT CLEAR (audit-incomplete; NO verdict from either auditor)**

| Auditor | Persisted state | Verdict reached? |
|---|---|---|
| A (`MICROHOM_AUDIT_V4_discovery-loss_A.md`) | Truncated at **CP2**. `harness.py` + `sweep.py` written to scratch. ONE decisive read confirmed (k=6, mh_frac=0.8: `delta_improve=6.0`, `mh_measured=0.667≥0.5`, discoverable guard-OFF → a real, discoverable cryptic that trips mh≥0.5). **No results TSV persisted; no verdict.** | **NO** |
| B (`…discovery-loss_B.md`) | **File does not exist.** Scratch dir `…/discovery-loss_B/` is **empty** — zero harness, zero results, zero `.md`. Total stall, nothing persisted. | **NO** |

**Agreement level: N/A — neither auditor completed.** A single confirmed read (A's CP2) demonstrates the
delta>0-yet-mh≥0.5 case is *real and constructible* (the exact case the V1 synthesis flagged as never-built),
but A never ran the sweep that would QUANTIFY the discovery-loss rate vs mh_frac × (margin, cap), and B
produced nothing. **Per the adversarial gate, this lens cannot be CLEAR.**

**Structural finding (carried from V2 CP2 + the detector records, labeled as adjacent evidence — the V4
discovery-loss auditors did NOT produce it):** the `(0, cap)` overlap is **IRREDUCIBLE**. V2's re-admission
harness proved the cap re-admits a flagged move precisely when `clamp(cap, hold, eff) ≤ delta_improve`, i.e.
on `delta_improve` alone; a REAL cryptic and a FAB drift with the SAME `delta_improve` are re-admitted or
vetoed IDENTICALLY. The cap therefore only **BOUNDS** the read-blind loss (protects strong-evidence moves);
inside the near-tie band it adds no separating signal. **→ A positional-distinctiveness / read-evidence signal
is REQUIRED to CLOSE the discovery-loss fault. The cap alone cannot.** (This is a design conclusion, not a
V4 discovery-loss clearance.)

---

## LENS 3 — DETECTOR-CORRECTNESS — **NOT CLEAR (A5 = consensus HOLDING fault; A8 disagreement unresolved)**

Both auditors completed and agree on the two decisive facts. `_effective_veto_margin` is CLEARED by both.

### A5 (`N==N` phantom microhomology) — **CONSENSUS HOLDING FAULT** (both auditors, consistent numbers)
`_frac_match` uses raw `x==y` with NO `.upper()` and NO non-ACGT exclusion, so `N==N` (and `n`, `X`,
lowercase-self) counts as a MATCH. `_move_microhomology` inherits it; `genome_seq` reaches the refiner with N
intact (`load_genome` uppercases but does NOT strip N — `genome.py:152/155`); no upstream pre-filter blocks an
N-region boundary in motif-blind mode. `_hp_run_across` shares the flaw (`HP('NNNNN',2,4)=5`).

| Evidence | Auditor A | Auditor B |
|---|---|---|
| Unit false-positive | `FM('NNNN','NNNN')=1.0`, `('nnnn','nnnn')=1.0`, `('XXXX','XXXX')=1.0` | `FM("NNNNNN","NNNNNN")=1.0`; tip `ACNTTT` vs `ACNGGG` = **0.5 trips** (correct 0.333 spares) |
| End-to-end decision FLIP (honest read, `delta_improve=6.0`) | Constructed masked-acceptor move: `MH=1.0` vetoes; unmasked identical geometry `MH=0.0` spared | `e2e2_results.json` / `b_e2e_independent.json`: **V2_N_phantom** actual mh=0.5 TRIPS (correct 0.333 spare), **V3_full_N** 1.0 TRIPS (correct 0.0 spare); both `FAULT_phantom_veto=true`; guard-OFF je=176 (moved) → guard-ON(m=8) je=170 (held) |
| Cap rescue | `EVM(0.5,1.5,{0,.25,1,5})={1.5,.5,1,1.5}` — nonzero veto every cap; cap on same axis as delta, cannot fix a wrong SCORE | cap=2 restores je=176 in V2/V3 (delta 6 > cap) but DEFAULT cap=0 leaves phantom veto UNBOUNDED across `(hold, hold+margin)` band |

The two records are numerically consistent (identical constructs, identical mh=0.5/1.0 vs correct 0.333/0.0
flips, identical guard-OFF→ON je 176→170). This is a **robust, byte-consistent consensus fault.**

**Note on B's record form:** B's `.md` ends after the RESULTS sections without a formal "VERDICT:" line, but
its RESULTS header states verbatim **"A5 (N==N phantom microhomology) — HOLDING-CLASS DETECTOR FAULT,"** and
the persisted JSON (`e2e2_results.json`, `b_e2e_independent.json`, `sparse_n.json`) backs it. Counted as a
COMPLETED holding finding.

**Realism / severity (honest scoping, does NOT clear the fault):** the bundled `--Scer` yeast S288C R64-5-1
has **N=0** (both auditors: total 12,157,105 bp, 0 N-runs) → A5 cannot fire on the default genome. User GRCh38
has ~150 Mb N but at assembly gaps where reads don't align, so a candidate junction boundary essentially never
lands in an N-run. B's `sparse_n.json` sharpens the trigger threshold: with random real bases, min N to trip
0.5 is 1 (k=4/6), 2 (k=8), 3 (k=10) — i.e. even a FEW N in the flank suffice. Net: A5 is a genuine SCORING
defect with a **narrow but non-empty** realistic trigger. Per the task's explicit decision rule ("if A5/A8 are
real, they are HOLDING faults"), narrow scope is severity context, NOT grounds to clear.

### A8 (`max` over both boundaries) — **A/B DISAGREE ON LABEL (not on fact) → treated as NOT-clear; moot (A5 already blocks)**
Both agree on the MECHANICS (verified against source `junction_refiner.py:499-532`): `_move_microhomology`
returns `max()` over the acceptor leg and the donor leg; a move shifting BOTH boundaries where the donor leg
sits in a repeat (frac 1.0) and the acceptor leg is a genuine transition (frac 0.0) yields combined MH=1.0 →
VETO, whereas acceptor-only (donor held) → MH=0.0 → spared.

- **Auditor A: A8 is a DISTINCT HOLDING fault.** Pure-ACGT, k=3 (within default `max_boundary_shift=50`), NO
  N required → fires on the DEFAULT `--Scer` genome. A donor-side `(CAG)n` microsatellite (ubiquitous at real
  splice boundaries) masks a genuine acceptor cryptic. **Violates the detector's own docstring** ("A move to a
  genuine sequence transition … is untouched, preserving motif-blind discovery") — the genuine acceptor
  transition IS touched. Reproducer persisted (`a8_reproducer.py`, asserts pass). Also: single-boundary
  acceptor into a 2 bp coincidental `AT|AT` (freq ~1/16) → MH=1.0 → veto, so even without both-boundary
  geometry short coincidental repeats over-trip 0.5.
- **Auditor B: A8 is NOT a separate fault** — over-broad SCOPE by design that collapses into the
  ALREADY-CLEARED read-blind structural point. `max` can only INFLATE (→ over-veto, never accept-a-fake); the
  return value 1.0 is a CORRECT report of the donor leg's real microhomology; the docstring states "returns
  the max over the moved boundaries." A real acceptor cryptic co-moving with a donor repeat is (on one score
  axis) indistinguishable from drift — the same structural limitation V2 cleared as honest.

**Reconciliation:** A and B disagree on no *fact* — only on whether the `max`-masking is a distinct
detector-correctness defect (A: the docstring's per-boundary spare promise is broken) or a re-description of
the cleared one-axis read-blindness (B: max is documented and numerically truthful). Under the adversarial
rule, **unresolved disagreement is NOT-clear and cannot license a flip.** It is **moot for the verdict**: A5
alone already blocks the detector lens. A8's operational value is that it sharpens the FIX (evaluate boundaries
independently rather than `max`), regardless of the label.

### `_effective_veto_margin` — **CONSENSUS CLEAR**
- A: 12/12 regime probes PASS (cap≤0 / eff==hold no-ops byte-identical; active = `max(hold,min(eff,cap))`;
  veto never below hold; monotone in cap; range `[hold,eff]`).
- B: 10-row truth table ALL PASS; invariant `hold ≤ got ≤ eff` holds every regime (hold>cap→hold, cap≤0→eff,
  eff==hold→eff, cap≥eff→eff, degenerate equals).
Both auditors: **`_effective_veto_margin` is CORRECT.** The near-tie cap does NOT rescue detector correctness —
it lives on the same score axis and cannot fix a wrong `_frac_match` value (both auditors, matching logic).

---

# ============ RECONCILED FINDINGS (the answers the brief asks for) ============

- **Detector consensus:** A5 = CONSENSUS HOLDING fault (both auditors, byte-consistent numbers, unit +
  end-to-end flip). `_effective_veto_margin` = CONSENSUS CLEAR (12/12 and 10/10). A8 = A/B disagree on label
  (not fact) → NOT-clear, but moot (A5 already blocks). **Detector lens: NOT CLEAR.**
- **Discovery-loss consensus:** NONE POSSIBLE — B persisted nothing, A truncated at CP2 with no results/verdict.
  **Discovery-loss lens: NOT CLEAR (audit-incomplete).** Adjacent evidence (V2 CP2) shows the `(0,cap)` overlap
  is irreducible on `delta_improve` alone → **a positional-distinctiveness signal is REQUIRED to close; the cap
  only BOUNDS.**
- **`(margin, threshold, cap)`:** **do NOT emit a shippable operating point.** No cap value rescues A5 (wrong
  frac SCORE, same axis). The blocking gate before any operating point:
  1. **Fix `_frac_match` to exclude non-ACGT** (`x==y and x in 'ACGT'`) — kills the `N==N`/`n`/`X` false
     positive in `_move_microhomology` AND the shared `_hp_run_across` N-run. Recommended INDEPENDENTLY by both
     auditors.
  2. **Evaluate the two boundaries INDEPENDENTLY (guard only the moved boundary), not `max` over both** — so an
     unrelated donor repeat cannot mask a genuine acceptor transition (addresses A8). Consider a minimum k
     (≥3–4) to stop 2 bp coincidental direct repeats tripping 0.5.
  3. **Then** run the never-completed discovery-loss quantification (A's `sweep.py` design: discovery-loss rate
     + fab FDR across mh_frac × (margin, cap)) on a real cryptic/R3 panel.
  4. Add the read-evidence / positional-distinctiveness gate the `(0,cap)` overlap requires.
- **Contingent-only prior operating point (NOT a clearance):** IF/WHEN the gate is later satisfied, the V1
  synthesis tentatively favored `microhom_threshold=0.5`, `microhom_drift_margin=3.0` (the data-supported
  discovery-neutral point; treat 8.0 as unproven) PLUS an absolute-near-tie gate. Marked explicitly
  contingent — it is NOT licensed by this V4 audit.
- **COMPASS real-data confirmation remains an INDEPENDENT ship prerequisite regardless** of any synthetic-audit
  outcome — the non-circular replacement for the retracted Snaptron verdicts. Required before the default is
  ever flipped, whether or not the synthetic audits pass.

---

# ============ THE ONE FAULT TO FIX FIRST ============

**`_frac_match` counts `N==N` (and `X`, lowercase-self) as a match → phantom microhomology ≥ 0.5 → false-positive
veto of a read-supported move (`delta_improve=6.0`).** This is the CONSENSUS detector fault (both auditors, unit
+ end-to-end decision-flip, backed by persisted JSON). Fix: `_frac_match` must treat any non-ACGT as a NON-match
(`sum(1 for x,y in zip(a,b) if x==y and x in 'ACGT')/len(a)`), which also protects the shared `_hp_run_across`.
Then evaluate boundaries independently (A8 fix), then quantify discovery-loss, then COMPASS. No cap value
substitutes for these — the cap bounds, it does not correct the score.

STATUS: COMPLETE.
