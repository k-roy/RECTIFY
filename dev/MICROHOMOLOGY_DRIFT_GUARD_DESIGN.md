# Microhomology-drift guard — design + Phase-1 validation (2026-07-11)

**Problem (spike-in + Sumner + COMPASS-recall-pending):** the motif-blind re-placer fabricates non-canonical
junctions via NON-HOMOPOLYMER drift (~1.2% FDR, spike-in ground truth) — it shifts a real canonical junction a
few-to-many bp to a nearby non-canonical pool member. The HP-drift guard catches ZERO of this (it only checks
homopolymer runs). Mechanism (explained + spike-in-confirmed): motif-blindness REMOVED the specificity prior; at
loci with LOCAL MICROHOMOLOGY (a near-tandem-repeat at the drift distance) + ONT error, the true vs drifted boundary
is a near-tie the error tips wrong. Pool-dependent (needs the false site in the pool = paralog/high-expression loci).

## The fix: generalize the HP-drift guard to a MICROHOMOLOGY-drift guard
A homopolymer is the MAXIMAL microhomology (every base identical). The non-HP drift is the SAME phenomenon at
PARTIAL microhomology. So generalize the guard's detector from "boundary inside a homopolymer run" to "the MOVE
sits in a microhomology context."

### Detector (MOVE-aware — matches what the guard sees: shift ns/ne -> js/je)
`move_microhomology_frac(genome, ne, je)` for an acceptor shift (donor symmetric): k = je - ne; the drift is
enabled iff the k-mer absorbed into the intron (genome[ne:je]) near-matches the new-exon k-mer (genome[je:je+k]).
Return the FRACTION of the k-mer that matches at the drift distance. This is EXACTLY the spike-in's installed
repeat genome[A:A+k] ~ genome[A+k:A+2k]. Check both donor and acceptor shifts; take the max.

### Guard rule (mirrors the HP-drift guard)
On a move, if `move_microhomology_frac >= mh_threshold` (~0.5), require the move to beat the incumbent by
`microhom_drift_margin` (a new param, default 0 = OFF = byte-identical). Real non-canonical discovery has LOW
microhomology -> no margin -> discovery preserved. ADD alongside `hp_drift_margin` (keep the validated HP-guard
intact; the new guard adds non-HP coverage). Byte-identical when the new param is 0.

## ★ PHASE 1 VALIDATION (spike-in ground truth, /tmp/mh_char2.py)
Move-aware microhomology_frac on the DRIFT move (true->drift) vs the REAL move (decoy->true non-canonical):
| set | microhomology_frac | 
| FAB-drift (installed, non-HP, k=6..28) | 0.60-0.71 (mean 0.66), consistent across all k |
| R3 genuine non-canonical (real splice) | 0.33 (no repeat) |
=> a threshold ~0.5 SEPARATES fabricated drift (>=0.60) from real discovery (0.33), across all drift distances.
The detector catches exactly the non-HP drift the HP-guard misses, without flagging real discovery. DESIGN VALIDATED.
CAVEAT: small n (FAB 5, R3 2); must re-validate the threshold on REAL SMA drift (SNRPN etc.) via COMPASS recall
(in flight) + a broader R3/cryptic-3'SS panel before locking the threshold.

## BUILD PLAN (rigor = HP-guard + concat-DP discipline)
1. Implement `move_microhomology_frac` + wire a `microhom_drift_margin` / `mh_threshold` into junction_refiner
   (mirror the hp_drift_margin plumbing through refine_read_junctions -> _run_sequential/_run_parallel ->
   refine_bam_junctions). Byte-identical when margin=0 (default).
2. Regression tests (mirror test_hp_drift_guard): detector unit tests (repeat vs random; both shift dirs; edges),
   guard vetoes the fab drift, spares the real transition, byte-identical at 0.
3. VALIDATE: re-run the spike-in fab panel -> does arm-Bguard-micro cut the ~1.2% non-HP drift FDR? AND the R3/
   discovery panel -> is recall preserved (like the HP-guard's zero-discovery-cost)? Tune mh_threshold + margin.
4. REAL-DATA: test the detector's threshold on the COMPASS-labeled SMA drift (real fabrication vs real discovery).
5. Adversarial audit (byte-identity + does it kill any real discovery) before flipping default.

## ★★★ PHASE 2-3 VALIDATION (2026-07-11) — DECISIVE, CLEAN WIN
Implemented (junction_refiner.py, committed): _move_microhomology detector + microhom_drift_margin/microhom_threshold
threaded through all 4 refine fns, byte-identical when margin=0 (refiner+guard suites 50 passed).
GROUND-TRUTH VALIDATION (spike-in fab panel + R3 discovery panel; per-read truth, seed excluded):
| arm | fabrication FDR (canon->drift) | R3 discovery (HP cell) | canonical (plain) |
| arm-B no-guard | 1.31% (24 reads; == spike-in ~1.2%) | 0.284 | 0.931 |
| microhom m=3   | 0.05% (1 read; 96% cut) | 0.284 PRESERVED (b=0,c=0) | 0.941 (+0.010, improved) |
| microhom m=8   | 0.00% (ELIMINATED) | 0.284 PRESERVED | 0.941 (improved) |
=> the general non-HP drift the spike-in/Sumner/Snaptron exposed is ELIMINATED, discovery preserved EXACTLY,
canonical slightly IMPROVED (guard also prevents spurious canonical drift). Zero discovery cost — the HP-guard bar.
Operating point: microhom_drift_margin=8.0, microhom_threshold=0.5 (full elimination, zero discovery cost).
Measurement note: the initial 5.7%->4.4% was a MEASUREMENT ARTIFACT (counted seed reads legitimately at drift as
fabrication); per-read-truth (canonical-origin reads moved to drift) is the correct FDR -> 1.31%->0%.
REMAINING: regression tests (mirror test_hp_drift_guard); adversarial audit (byte-identity + no-discovery-loss);
COMPASS real-data threshold confirmation (in flight); then flip default / ship alongside the HP-guard.

## ⛔ PHASE 3 AUDIT → HOLD (2026-07-11, workflow wte43x5rc; dev/MICROHOM_AUDIT_SYNTHESIS.md)
The triple adversarial audit returned **HOLD — do NOT flip the default.** byte-identity CLEARED (guard inert
at default, 1653 passed) but the audit found a **genuine design fault the synthetic validation above MISSED:**
`_move_microhomology` is **READ-BLIND** (genome-only). The mh≥0.5 veto trigger (genomic) and the score
delta_improve that must clear the margin (read evidence) are INDEPENDENT → a real cryptic the read distinguishes
(delta_improve>0, discovered with guard OFF) can still trip mh≥0.5 and be vetoed whenever delta_improve < margin.
At m=8 that suppresses any real cryptic differing by <~8 in-window bases. The m=8 "zero discovery cost" above is
an ARTIFACT of the under-powered panel (5 fab / 2 real): the panel doesn't populate the (0, margin) delta band
where the fault lives. Two auditors (discovery-loss, detector-correctness) STALLED on API errors → gate unsatisfied.

## ★ PHASE 4 FIX (2026-07-12) — near-tie read-evidence cap (`drift_near_tie_cap`)
Added `_effective_veto_margin(hold_margin, eff_margin, drift_near_tie_cap)` +
`drift_near_tie_cap` param threaded through all 4 refine fns (default 0.0 = disabled = byte-identical).
Rule: `veto_margin = max(hold_margin, min(eff_margin, drift_near_tie_cap))` when the cap is active —
caps the read-BLIND drift margins (hp_drift + microhom_drift) so a move with `delta_improve >= cap` is
NEVER drift-vetoed, while NEVER capping `hold_margin` (read-agnostic blunt prior; it stays a floor).
**HONEST SCOPE (advisor, load-bearing):** delta_improve, eff_margin and the cap share the SAME score
axis, so the cap **BOUNDS** the read-blind discovery-loss (protects moves with strong read evidence) but
adds **NO discriminating signal inside the (0, cap) near-tie band** — real cryptics and fabricated drift
still overlap there. Its value is STRUCTURAL: decouple hold_margin from the drift margins + make the
discovery-loss ceiling explicit and tunable. If the (0, cap) overlap proves fatal to discovery, the only
real remedy is a DIFFERENT signal (per-base positional distinctiveness the incumbent soft-clip cannot
absorb) — NOT built speculatively; the independent audit panel decides whether it's needed.
Provisional operating point (to be validated/tuned by the independent audit, NOT self-tuned): margin=3.0
(m=8→m=3 per audit; m=3 already 96% fab-suppression), microhom_threshold=0.5, drift_near_tie_cap=TBD.
Internal check ONLY (per advisor — not efficacy validation): byte-identical at default, suite green,
_effective_veto_margin arithmetic + hold-uncapped interaction unit-tested (test_microhom_drift_guard.py).
The load-bearing validation (real-cryptic-microhomology panel at varied delta_improve → discovery-loss vs
(margin, cap)) is DELEGATED to the independent Opus-Max audit so it isn't a self-tuned rubber stamp.

## ⚖ AUDIT-V4 (2026-07-13, 4 Opus-Max, 2/task, robust consensus) → HOLD + DETECTOR FIXES
`dev/MICROHOM_AUDIT_V4_SYNTHESIS.md`. Despite 7/9 agents hitting API stalls, the redundant (2-per-task)
design produced a robust CONSENSUS. Findings:
- **`_effective_veto_margin` (the cap): CONSENSUS CLEAR** (A 12/12 probes; B 10-row truth table). The cap is
  arithmetically correct; it just can't rescue a wrong detector `frac`.
- **A5 (CONSENSUS HOLDING FAULT, both detector auditors):** `_frac_match` scored `N==N` (+ non-ACGT,
  lowercase) as a MATCH → phantom microhomology ≥0.5 → falsely vetoed a read-supported move
  (delta_improve=6.0). Shared by `_hp_run_across` (NNNNN→run 5). Pipeline realization narrow today (bundled
  yeast N=0; human N at gaps) but intrinsic.
- **A8 (A: distinct fault / B: collapses into read-blind point — agree on mechanics):** `_move_microhomology`
  `max()` over both boundaries let a donor-side (CAG)n repeat MASK a genuine acceptor transition → over-veto,
  fires on the default genome. Moot for the verdict (A5 blocks) but sharpens the fix.
- **discovery-loss RATE: still UNMEASURED** (both auditors stalled) — but A confirmed the delta>0 ∧ mh≥0.5
  case is REAL and constructible; the (0,cap) overlap is IRREDUCIBLE → a positional-distinctiveness signal is
  REQUIRED to CLOSE (cap only BOUNDS).

## ★ PHASE 5 DETECTOR FIXES (2026-07-13) — A5 + A8 (byte-identical off; guard still default OFF)
1. **A5:** `_frac_match` now counts a match only when `x == y and x in _ACGT` (module const `frozenset("ACGT")`);
   `_hp_run_across` returns 0 when the boundary base ∉ _ACGT. Kills phantom microhomology / phantom HP on N/gap
   runs. Tests: `TestFracMatch.test_ambiguity_bases_do_not_match`, `TestHpRunAcrossAmbiguity`.
2. **A8:** `_move_microhomology` combines the MOVED boundaries by **min** (was max): a move is drift-suspect
   only if EVERY shifted boundary sits in microhomology; a genuine transition on any moved boundary spares the
   whole move (an edge-truncated moved boundary contributes 0.0 → spares). Single-boundary moves unchanged
   (min == that frac), so all prior detector tests + the validated single-boundary fab suppression are
   preserved. Test: `TestMoveMicrohomology.test_both_boundary_transition_not_masked_by_other_repeat` (A's exact
   pure-ACGT reproducer). Deferred: min-k sensitivity floor (needs the discovery-loss panel to tune).
STILL OPEN before ANY enable: (i) quantify discovery-loss rate (stalls); (ii) add the positional-
distinctiveness signal to CLOSE (cap only bounds); (iii) COMPASS real-data confirmation (independent).
