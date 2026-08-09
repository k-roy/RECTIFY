# MICROHOM-DRIFT GUARD — AUDIT SYNTHESIS

**Date:** 2026-07-11
**Synthesizer:** subagent (agent-a25a2c1e784ad37dc)
**Question:** Should we flip the default `microhom_drift_margin` from 0.0 (inert) to enabled, shipping the
microhomology-drift guard alongside the HP-drift guard?
**Policy:** Flip the default ONLY if all three auditors are all-clear (no holding fault).

---

## PLAN

1. Read the three auditor findings + the byte-identity durable record (`dev/MICROHOM_AUDIT_byte-identity.md`).
2. Establish ground truth on the operating point: threshold, margin, and the size/composition of the tuning set.
3. Evaluate each auditor's strongest challenge; decide whether it holds.
4. Weigh the DISCOVERY-LOSS concern (0.5 threshold set on 2 real cases suppressing real non-canonical junctions).
5. SHIP vs HOLD verdict + specific design/threshold change if HOLD.
6. Decide whether the COMPASS real-data threshold check is still needed.
7. Persist verdict here.

## CHECKPOINTS
(appended as each sub-step completes)

### CP1 — inputs read (2026-07-11)
- Only ONE auditor finding was passed in the task prompt: the BYTE-IDENTITY lens (verdict: CLEARED, no fault).
- Found on disk THREE auditor durable records:
  - `dev/MICROHOM_AUDIT_byte-identity.md` (206 lines) — COMPLETE, verdict CLEARED.
  - `dev/MICROHOM_AUDIT_discovery-loss.md` (99 lines) — **STALLED at CHECKPOINT 3** (no verdict written).
  - `dev/MICROHOM_AUDIT_detector-correctness.md` (43 lines) — **STALLED at PLAN** (zero checkpoints, no verdict).
- Design/result docs: `dev/MICROHOMOLOGY_DRIFT_GUARD_DESIGN.md`, `dev/HANDOFF_MICROHOM_GUARD.md`.
- Operating point under consideration: **microhom_threshold=0.5, microhom_drift_margin=8.0**.
- Tuning set: Phase-1 threshold-separation = **5 fab cases vs 2 real R3 cases** (under-powered, per design doc CAVEAT).
- Phase 2-3: fab FDR 1.31%->0.00% at m=8; R3 discovery recall 0.284 PRESERVED EXACTLY; canonical +0.010.

### CP2 — guard code verified (byte-identity leg confirmed independently)
- `junction_refiner.py:788-805`: guard block gated on `microhom_drift_margin > 0.0`; at default it is skipped,
  `_move_microhomology` never called, `eff_margin` stays `hold_margin` (0.0), and the veto at line 799 requires
  `eff_margin > 0.0` — short-circuits before `incumbent_score`/`best_score_cmp` are compared. BYTE-IDENTITY HOLDS.
- Veto rule when ON: `if _move_microhomology(...) >= microhom_threshold: eff_margin += microhom_drift_margin`,
  then `if best_score_cmp > incumbent_score - eff_margin: moves = False`. (score_cmp: lower = better, `<` in tuple cmp.)
- KEY: the discovery-loss auditor's mechanistic finding (scorer soft-clips the drift-distance prefix => pure
  microhomology ties score delta=0 => refiner tie-break holds incumbent REGARDLESS of guard) is the crux of
  whether 0.5 is safe. That auditor STALLED before constructing the "cryptic scores strictly better yet trips
  mh>=0.5" case that would actually demonstrate (or refute) discovery loss.

### CP3 — detector is read-blind (the discovery-loss path is mechanically LIVE, not defused)
- `_move_microhomology(seq=genome_seq, ns, ne, js, je)` (junction_refiner.py:499-532) operates PURELY on the
  genome and the two candidate coordinate pairs. It NEVER sees the read/query evidence.
- Consequence: the veto trigger `mh >= 0.5` is decided by GENOMIC context alone; the score delta that must clear
  the margin is decided by READ evidence. These are INDEPENDENT. So a cryptic acceptor that the READ genuinely
  distinguishes (scores strictly better, delta>0, discovered with guard OFF) can STILL sit in a genomic
  partial-repeat that trips mh>=0.5 and get vetoed whenever its evidence delta < margin (8.0 edit-dist units over
  the 30bp scored window — a very high bar; a real cryptic differing by only 2-3 in-window bases has delta << 8).
- The HANDOFF's "guard is PROVABLY IRRELEVANT" claim is TRUE ONLY for the pure-tie case (delta=0, incumbent held
  regardless of guard). It does NOT cover the delta>0-yet-mh>=0.5 case, which the discovery-loss auditor NAMED at
  its checkpoint 3 as "the case that actually matters" and then STALLED before building. That case is
  mechanistically live => the flagged discovery-loss fault is OPEN, not closed.

---

# ===== VERDICT =====

## (1) Per-auditor: strongest challenge + does it hold?

### A. BYTE-IDENTITY (COMPLETE — CLEARED, but does NOT endorse the flip)
- **Strongest challenge:** the parallel path (n_workers>1) produced sha256 DIFFER between HEAD (guard present,
  default off) and the pre-guard commit on the same input — an apparent byte-identity violation.
- **Does it hold? NO.** The DIFFER is guard-INDEPENDENT scheduling-order noise: a HEAD-vs-HEAD self-consistency
  run also gives raw-order DIFFER but position-sorted MATCH. Sequential (n_workers=1) HEAD-vs-PRE is byte-identical
  WITH raw order preserved across all 6 panel/config cells; parallel is byte-identical after coordinate-sort
  normalization. Confounder (concat-DP) isolated by comparing against the DIRECT PARENT b6f07f7 (only
  junction_refiner.py changed b6f07f7..HEAD). Structural: microhom feat commit has ZERO deleted lines (pure +74);
  at default the guard block is skipped and the veto short-circuits before any incumbent/score read. pytest
  -m "not slow" = 1653 passed; test_microhom_drift_guard.py 14/14.
- **VERDICT: CLEARED on the inertness/byte-identity leg.** BUT the auditor EXPLICITLY declined to endorse enabling
  the guard: "this audit does NOT endorse enabling the guard — the 0.5 / 8.0 operating point (set on 5 fab / 2
  real, under-powered) is an efficacy concern outside this lens." So even the one cleared auditor does not clear
  the FLIP — only the "off = inert" property.

### B. DISCOVERY-LOSS (INCOMPLETE — STALLED at checkpoint 3, NO VERDICT)
- **Strongest challenge (the one the task flags as most-likely-genuine):** the 0.5 threshold, set on only 2 real
  R3 cases, could veto a REAL non-canonical (cryptic 3'SS) junction that incidentally sits next to a short direct
  repeat (mh >= 0.5) — suppressing genuine discovery.
- **Does it hold? UNRESOLVED, and leaning OPEN.** The auditor established two useful mechanistic facts:
  (i) in a PURE microhomology tie the scorer soft-clips the drift-distance prefix so delta=0 and the refiner's
  tie-break holds the incumbent REGARDLESS of the guard (guard is irrelevant there — not a discovery loss caused
  by the guard); (ii) genuine guard-caused loss therefore requires a cryptic site that scores STRICTLY BETTER
  (delta>0) yet still trips mh>=0.5. The auditor NAMED this decisive case at checkpoint 3 and then STALLED before
  building it. Independent code inspection here (CP3) confirms this path is MECHANICALLY LIVE: `_move_microhomology`
  is READ-BLIND (genome-only), so mh>=0.5 (genomic) and delta>0 (read evidence) are independent conditions that
  can co-occur — a real cryptic differing by 2-3 in-window bases (delta << 8) in a partial-repeat context WOULD
  be vetoed. **The challenge is not defeated; it is unadjudicated with the mechanism pointing toward it holding.**

### C. DETECTOR-CORRECTNESS (INCOMPLETE — STALLED at PLAN, NO CHECKPOINTS, NO VERDICT)
- **Intended strongest challenge:** find an input where `_move_microhomology` scores WRONG — false-neg (misses
  real drift) or false-pos (flags a real discovery). Planned coverage: donor-side drift, genome-edge k-mers,
  k=0/1/huge, non-ACGT/N bases, overlapping/nested repeats, both-boundary shifts, max-over-both masking, and
  acceptor-vs-donor geometry.
- **Does it hold? UNRESOLVED — the audit never ran.** Zero checkpoints were written. The detector's correctness
  (esp. the max-over-both-boundaries masking case and non-ACGT handling) is therefore UNVERIFIED by adversarial test.

## (2) SHIP or HOLD?  ==> **HOLD.**

The decision is fixed by POLICY alone, independent of any further compute:
- **PRIMARY (blocking) reason — the audit gate is unsatisfied.** Policy: "flip the default only if ALL CLEAR (no
  holding fault)." Two of three auditors never reached a verdict (discovery-loss stalled at CP3; detector-
  correctness stalled at PLAN). Absence of a completed audit is NOT "all clear" — it is "audit incomplete," which
  fails the gate by construction. This alone forbids the flip.
- **SECONDARY (reinforcing) — even the completed auditor declined to endorse the flip.** Byte-identity cleared
  only the off=inert property and explicitly named the 0.5/8.0 operating point as an out-of-lens efficacy concern.
- **SECONDARY (reinforcing) — the flagged discovery-loss fault is mechanically live.** The read-blind detector
  makes the delta>0-yet-mh>=0.5 veto path real; no evidence exists that it never fires on genuine cryptic sites,
  and the tuning set (5 fab / 2 real) is far too small to bound the discovery-loss rate.

**Ship the guard as CODE with default OFF (already done, byte-identical) — do NOT flip the default.**

## (3) If HOLD — the specific design/threshold change to pursue before re-auditing

Ranked, cheapest-first:

1. **Complete the two stalled audits FIRST (no threshold change until then).** The gate is procedural: finish
   discovery-loss (build the delta>0-yet-mh>=0.5 case; sweep discovery-loss rate vs microhomology_frac and vs
   margin) and detector-correctness (run the planned edge-case matrix). Relaunch each pointed at its partial
   `dev/MICROHOM_AUDIT_{discovery-loss,detector-correctness}.md` (resume, not restart).

2. **Prefer the smaller margin m=3 over m=8 at the SAME threshold 0.5.** Existing tuning data already show m=3
   gives 96% fab suppression (1.31% -> 0.05%, 1 residual read) with R3 discovery PRESERVED at 0.284 exactly and
   canonical improved — i.e. m=3 is discovery-neutral ON THE TUNING PANEL. The m=3 -> m=8 jump buys only the last
   0.05% of fab suppression while ~2.7x-ing the evidence bar a real cryptic must clear to survive a mh>=0.5 flag.
   In edit-distance units over a 30bp window, m=8 vetoes any real cryptic that differs by fewer than ~8 in-window
   bases if it trips the (genomic) threshold — a plausibly common case. **m=3 is the safer operating point the
   current data already support; do not adopt m=8 without discovery-loss data quantifying its extra cost.**

3. **Add a discovery-loss safeguard that consults READ evidence, not just genomic context.** The root risk is the
   detector's read-blindness. Options (design-time): (a) only apply the veto when the evidence delta is SMALL in
   absolute terms (a true near-tie), not merely smaller than a fixed margin — i.e. gate on `|delta| < small` so a
   strongly-supported cryptic is never vetoed regardless of genomic microhomology; (b) exempt moves whose
   candidate exon2 carries distinctive in-window bases the incumbent soft-clip cannot absorb (breaks the
   soft-clip-tie masking); (c) raise microhom_threshold above 0.5 toward the fab/real separation midpoint once
   a larger real panel pins the real-cryptic microhomology distribution (Phase-1 measured real=0.33 on n=2 only).

4. **Enlarge the tuning/validation set before locking ANY threshold.** 5 fab / 2 real cannot bound either the fab
   FDR or the discovery-loss rate. Need a broad R3/cryptic-3'SS panel at varied microhomology_frac (the discovery-
   loss auditor's step 5) PLUS the real-data COMPASS confirmation (item 4 below).

## (4) Is the COMPASS real-data threshold check still needed regardless?  ==> **YES, required regardless.**

The 0.5 threshold and 8.0 margin were tuned on synthetic spike-in data with only 2 real cases. A synthetic
discovery-loss test (even completed) cannot substitute for confirming, on non-simulated data, that the threshold
separates REAL SMA drift (fabrication) from REAL discovery. COMPASS is a separate, independent gate — needed
whether or not the synthetic audits pass, and needed before the default is ever flipped. It is also the
non-circular replacement for the retracted Snaptron verdicts (Snaptron is STAR-built/motif-filtered => confirm-
only). Keep it in flight; it is a ship prerequisite, not optional.

---

## RECOMMENDED OPERATING POINT (if/when the gate is later satisfied)
- microhom_threshold: **0.5** (unchanged, pending a larger real panel to pin real-cryptic microhomology dist).
- microhom_drift_margin: **3.0**, NOT 8.0 — the data-supported discovery-neutral point; treat 8.0 as unproven.
- PLUS an absolute-near-tie gate on the veto (`|delta| < small`) to neutralize the read-blind-detector risk.
- All contingent on: two stalled audits completed + COMPASS real-data confirmation.

STATUS: COMPLETE.
