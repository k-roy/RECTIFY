# MICROHOM-DRIFT GUARD + NEAR-TIE CAP — AUDIT V2 SYNTHESIS

**Date:** 2026-07-12
**Synthesizer:** subagent (agent-a25a2c1e784ad37dc), Opus 4.8. READ-ONLY (no edits/commits/default flips).
**Question:** Should we flip the default `microhom_drift_margin` from 0.0 (inert) to enabled — shipping the
microhomology-drift guard, now with the Phase-4 near-tie read-evidence cap (`drift_near_tie_cap`), as an
active default?
**Under audit:** commit `05664bc` "feat(refiner): near-tie read-evidence cap on the read-blind drift guards"
(HEAD of `worktree-agent-a25a2c1e784ad37dc`). Cap commit range = junction_refiner.py + tests + design.md ONLY
(git-confirmed).
**Policy:** `all_clear=true` (→ flip permitted) ONLY if all THREE auditors COMPLETED (reached a verdict) AND
none found a holding fault. A stalled/incomplete leg is NOT all-clear.
**Continuity:** supersedes the prior HOLD in `dev/MICROHOM_AUDIT_SYNTHESIS.md` (the pre-cap V1 round, which
also HELD, on the same gate-unsatisfied grounds + the read-blind discovery-loss concern).

---

## THE FIX UNDER AUDIT (Phase 4)

`_effective_veto_margin(hold_margin, eff_margin, drift_near_tie_cap)` (junction_refiner.py:535-558):
```
veto_margin = max(hold_margin, min(eff_margin, drift_near_tie_cap))   when cap active (>0 AND eff>hold)
            = eff_margin                                              when cap<=0 OR eff==hold (byte-identical)
```
Threaded through `refine_read_junctions` / `_run_sequential` / `_run_parallel` / `refine_bam_junctions`.
- **CLAIM (code + docstring):** a move with `delta_improve >= cap` is NEVER drift-vetoed → the cap *bounds*
  the read-blind discovery loss.
- **CONCESSION (docstring, 549-554):** `delta_improve`, `eff_margin`, and the cap share the SAME score axis,
  so the cap "does not add a discriminating signal inside the `(0, cap)` near-tie band, where real cryptics
  and fabricated drift still overlap." The cap **bounds but does not close** the fault.

---

## COMPLETION STATUS OF THE THREE V2 AUDITORS (the gate)

The three durable `dev/MICROHOM_AUDIT_V2_*.md` records are the source of truth (returned findings were
`[null,null,null]` — all three agents dropped their stream; only the on-disk `.md` persists). Read in full:

| Leg | File (V2) | State | Verdict? |
| --- | --- | --- | --- |
| **discovery-loss** (LOAD-BEARING) | `MICROHOM_AUDIT_V2_discovery-loss.md` (80 ln, mtime 14:51) | **STALLED at CHECKPOINT 1** (orientation only) | **NONE** |
| **detector-correctness** | `MICROHOM_AUDIT_V2_detector-correctness.md` (63 ln, mtime 14:50) | **STALLED at PLAN** (0 checkpoints) | **NONE** |
| **byte-identity + structural-honesty** | `MICROHOM_AUDIT_V2_byte-identity.md` (175 ln, mtime 15:13) | **substantively COMPLETE** (CP0–CP3 done) | CLEAR / HONEST |

**→ 2 of 3 legs never reached a verdict.** By policy this fails the gate by construction — regardless of any
further compute, `all_clear` CANNOT be true. This is the DECISIVE, blocking reason for the call below.

### discovery-loss (LOAD-BEARING) — STALLED at CHECKPOINT 1
Did orientation only: verified guard + `_effective_veto_margin` + scorer + fab-harness code, restated the
irreducibility question, and named the decisive panel. Then STALLED. **PLAN steps 3–5 — the entire empirical
core — never ran:** it never built the read-cryptic panel (real cryptic query + incumbent N-op at the wrong
canonical site) and never measured `DISCOVERY-LOSS RATE = f(margin ∈ {3,8}) × (cap ∈ {off,1,2,3}) ×
microhomology_frac`. So the **empirical magnitude** of discovery loss under the cap is **UNMEASURED**. No verdict.

### detector-correctness — STALLED at PLAN (zero checkpoints)
Transcribed `_move_microhomology` + `_effective_veto_margin` and laid out a 12-case adversarial matrix (A1–A12
+ the `_effective_veto_margin` regime grid). Appended ZERO checkpoints. The detector's correctness is therefore
UNVERIFIED by adversarial test. The two named risks worth a concrete relaunch brief:
- **A5** — `_frac_match` compares by raw char equality, so on ambiguity/`N` runs `N==N` scores as a MATCH →
  a genome N-run near a boundary could spuriously trip `mh >= 0.5` (false-positive veto). UNTESTED.
- **A8** — `_move_microhomology` returns `max` over BOTH moved boundaries; a benign transition shift on one
  boundary could be MASKED by an unrelated repeat on the other boundary that also moved → spurious veto.
  UNTESTED. (This is the headline detector concern.)
No verdict.

### byte-identity + structural-honesty — substantively COMPLETE (CLEAR / HONEST)
Reached verdicts on both duties; the only loose end is a trailing "count below" on a `-x` pytest re-run (see
"not a fault" below). Findings (independently reproduced by that auditor, and the two load-bearing code facts
re-confirmed by me from source):
- **A. BYTE-IDENTITY = CLEAR.** HEAD (`05664bc`, cap present, default) vs the pre-cap parent (`05664bc^`,
  git-proven identical to the earlier-cleared guard state `17e57ba`) is byte-identical:
  - SEQUENTIAL (n_workers=1, RAW record order, SHA256 over full SAM text): MATCH on mix_r3b / mix_fair /
    fab_sweep × {default, motifblind, **guard_on_cap_off**}. The `guard_on_cap_off` config (motif-blind +
    `microhom_drift_margin=8` + `cap=0`) drives `eff_margin>0` so `_effective_veto_margin` ACTUALLY RUNS and
    the veto fires (refined counts drop, e.g. 815→308) — yet output stays identical, proving
    `_effective_veto_margin(...,cap=0.0)==eff_margin` exactly. STRONGER than the V1 audit, which never
    executed the replaced veto line.
  - PARALLEL (n_workers=3, POSITION-SORTED digest): all 9 cells MATCH; raw-order parallel divergence is
    guard-independent scheduling noise (HEAD-vs-HEAD self shows the same), correctly not flagged.
- **B. STRUCTURAL HONESTY = HONEST (neither over- nor under-claimed).** Verified on both the pure function and
  end-to-end:
  - **P1 one scalar axis:** `_effective_veto_margin` signature carries NO read/seq input; the caller uses the
    result ONLY as `best_score_cmp > incumbent_score - veto_margin` (⟺ `delta_improve < veto_margin`). The cap
    adds no new term, no read feature, no extra comparison.
  - **P2 cap ONLY lowers, never raises** (`veto_margin <= eff_margin`, 0 violations in 200k+grid) → relative
    to no-cap the cap can only ADD moves, never remove a real discovery.
  - **P4 clamp identity / P5 disabled-branch identity / P6 monotone-in-cap** all confirmed.
  - **P7 ergonomic sharp edge (flag, not a fault):** `f(cap=0)=eff` (full veto) but `f(cap=0+ε)=max(hold,ε)`
    (near-max protection). "Small cap" ≠ "small effect" — any tiny positive cap is a near-maximal
    guard-weakening. Worth documenting for whoever tunes it.
  - **End-to-end re-admission (harness_readmit2.py) — the decisive mechanism check:** on a constructed case
    (period-6 repeat so `+6` is drift-FLAGGED `mh=0.667`, built from the `+6` frame so it scores strictly
    better, `delta_improve≈6`, with a 1× acceptor boundary error so the N-op is actually scored), guard-OFF
    MOVES (discovers +6), guard-ON-no-cap VETOES (the read-blind fault, live), and the cap re-admits the
    flagged move **iff `clamp(cap,hold,eff) <= delta_improve`**. Adversarial correction to the task's phrasing:
    it is LOWERING the cap toward `hold` that re-admits; RAISING it toward `eff` restores full suppression
    (task said "raising" — directionally inverted; reported as a correction, not a code fault).
- **Not a fault:** `pytest -m "not slow"` stopped under `-x` at `test_restore_polya_from_parquet.py::
  test_restore_cat3_plus_2` — a missing-DATA fixture assertion (`validation_reads_polya_trim_metadata.tsv`
  absent in this worktree). Git-confirmed OUTSIDE the cap commit range; references junction_refiner/microhom/
  cap 0 times. Orthogonal to the cap. The guard-relevant suites are green:
  `test_microhom_drift_guard.py` 20/20 (14 pre-cap + 6 new cap), `test_junction_refiner.py` +
  `test_hp_drift_guard.py` 46 passed / 17 skipped.

---

## THE DISCOVERY-LOSS QUESTION — SPLIT INTO PROVEN vs UNMEASURED

The task's conditional ("if discovery-loss shows the cap is INSUFFICIENT — the `(0,cap)` overlap is
irreducible — the call is HOLD + a positional-distinctiveness signal is needed") must be answered carefully,
because the discovery-loss leg that was *supposed* to answer it STALLED. The question has two halves, and they
landed on opposite sides:

### PROVEN (structural) — the cap CANNOT CLOSE the read-blind fault; a positional signal IS required to close it
This is **established**, but by the **COMPLETED byte-identity/structural-honesty leg** + direct source, NOT by
the stalled discovery-loss leg:
- `_move_microhomology` is **read-blind** — signature `(seq, ns, ne, js, je)`, genome + coordinates only, never
  sees the read/query (junction_refiner.py:499-532, re-confirmed by me from source).
- The cap is a pure **aggregate-delta** operation on the **one shared score axis** (P1), and re-admits a move
  purely on `delta_improve` (end-to-end re-admission finding). Therefore a REAL cryptic and a FABRICATED drift
  that happen to share the same `delta_improve` are treated **IDENTICALLY** — both re-admitted iff
  `veto_margin <= delta_improve`, both vetoed otherwise.
- Consequence: **no aggregate-delta cap can discriminate real from fabricated inside the `(0, cap)` band.** The
  overlap is irreducible *with an aggregate-delta signal*. The code's own docstring concedes exactly this. So
  "close the fault → need a POSITIONAL-DISTINCTIVENESS signal the incumbent soft-clip cannot absorb" is a
  **defensible affirmative design conclusion.** (The soft-clip mechanism — scorer free-prefixes the first
  `drift` rescue bases, which is why a pure-microhomology tie scores `delta=0` and is guard-irrelevant — is
  precisely why distinguishing evidence must be POSITIONAL / outside that prefix. Both V2 orientation notes
  established this.)

### UNMEASURED (empirical) — the discovery-loss RATE/MAGNITUDE under the cap is NOT quantified
The cap **bounds** the loss (moves with `delta_improve >= cap` always survive — P2 + re-admission), but **how
often real cryptics actually fall in `(0, cap)`** on a realistic panel — the discovery-loss *rate* as a
function of `(margin, cap, microhomology_frac)` — was the entire job of discovery-loss PLAN steps 3–5, which
**never ran**. So there is **no measured number** bounding the residual discovery loss, and therefore **no
empirically-validated `(margin, threshold, cap)` operating point.**

**Net (state exactly this, do not blur it):** the *structural* half of discovery-loss's question is answered —
the cap bounds but cannot close, and closing it needs a positional signal. The *empirical magnitude* half is
UNANSWERED because the load-bearing leg stalled. Do NOT write "discovery-loss showed the cap is insufficient"
(that leg reached no verdict); do write "the completed structural analysis proves the cap cannot close the
read-blind fault, and the empirical loss rate that would justify a specific operating point was never measured."

---

# ===== VERDICT =====

## SHIP-ENABLED (flip the default) vs HOLD  →  **HOLD**

**DECISIVE (blocking) reason — the audit gate is UNSATISFIED.** Policy: flip only if ALL THREE auditors
COMPLETED with no holding fault. **Two of three never reached a verdict** (discovery-loss STALLED at CP1 /
orientation; detector-correctness STALLED at PLAN / zero checkpoints). An incomplete audit is not "all clear" —
it "failed to run." This alone forbids the flip and fixes `all_clear=false`, independent of any further compute.

**Reinforcing (not required for the call):**
1. The **load-bearing** leg is the one that stalled. It exists specifically to measure the discovery-loss rate
   the flip's safety depends on; that number does not exist.
2. **Structurally, the cap cannot close the read-blind fault** (proven above): it bounds discovery loss but adds
   no discriminating signal inside `(0, cap)`; a POSITIONAL-DISTINCTIVENESS signal (one the incumbent soft-clip
   cannot absorb) is required to actually separate real cryptics from fabricated drift there. Shipping the guard
   as an active default while the fault is only bounded — not closed — with an unmeasured residual rate is
   premature.
3. The operating point in play (`microhom_threshold=0.5`, `microhom_drift_margin=8.0`, and the prior synthesis's
   fallback `m=3`) is tuned on **5 fab / 2 real** cases — exactly the under-powering the stalled discovery-loss
   leg was chartered to attack.

**Ship the guard + cap as CODE with default OFF (already done — byte-identically inert, CLEARED). Do NOT flip
the default.**

## RELAUNCH / FIX PLAN (before any re-audit or flip)

1. **PROCEDURAL (primary) — relaunch the two stalled legs from their partial `.md`, no threshold change until
   then.** They persisted their orientation/plan, so resume (not restart):
   - **discovery-loss** ← `dev/MICROHOM_AUDIT_V2_discovery-loss.md` @ CHECKPOINT 1: execute PLAN steps 3–5 —
     build the read-cryptic panel (real cryptic query, incumbent N-op at the wrong canonical site, distinguishing
     evidence placed OUTSIDE the soft-clipped drift-prefix), then measure DISCOVERY-LOSS RATE across
     `margin ∈ {3,8} × cap ∈ {off,1,2,3} × microhomology_frac`, and re-measure fab suppression on the same grid.
     Deliver a concrete `(margin, cap)` with fab≈0 AND a *bounded, measured* discovery loss — or state it is
     irreducible with an aggregate-delta signal.
   - **detector-correctness** ← `dev/MICROHOM_AUDIT_V2_detector-correctness.md` @ PLAN: run the A1–A12 matrix,
     prioritizing **A5** (`N==N` false-positive on ambiguity runs) and **A8** (max-over-both-boundaries masking).
2. **DESIGN (secondary, structural) — to CLOSE (not merely bound) the read-blind fault, add a
   POSITIONAL-DISTINCTIVENESS signal.** The cap is on the wrong axis (aggregate delta). Candidate directions:
   (a) gate the veto on the candidate exon2 carrying distinctive in-window bases the incumbent soft-clip cannot
   absorb (breaks the soft-clip tie-masking); (b) an absolute near-tie gate (`|delta| < small`) so a
   strongly-supported cryptic is never vetoed regardless of genomic microhomology. Keep the cap as the bound;
   the positional signal is what actually separates real from fabricated inside `(0, cap)`.
3. **REAL-DATA (independent gate, required regardless) — the COMPASS real-data threshold check.** 0.5/8.0 (and
   m=3) were tuned on synthetic spike-ins with 2 real cases. A synthetic discovery-loss test, even completed,
   cannot substitute for confirming on non-simulated data that the threshold separates REAL drift (fabrication)
   from REAL discovery. This is a ship prerequisite, not optional, and the non-circular replacement for the
   retracted STAR-built/motif-filtered Snaptron verdicts.

## RECOMMENDED OPERATING POINT  →  **NONE is defensible yet — the blocking gap is the unmeasured discovery-loss rate.**

No validated `(margin, threshold, cap)` can be recommended: the empirical discovery-loss rate that would
justify one was never measured (discovery-loss stalled at orientation). If any point is carried forward as a
*provisional* starting hypothesis for the relaunch, it is `threshold=0.5, margin=3.0` (the prior synthesis's
discovery-neutral-on-the-tuning-panel point) — but flagged as **TUNING-PANEL-ONLY (5 fab / 2 real,
under-powered)**, which is exactly the under-powering the stalled leg exists to attack, so it is NOT a
recommendation. Any real operating point additionally requires: (i) the two stalled legs completed with
verdicts, (ii) a POSITIONAL-DISTINCTIVENESS signal to close (not just bound) the read-blind fault, and (iii)
COMPASS real-data confirmation.

STATUS: COMPLETE.
