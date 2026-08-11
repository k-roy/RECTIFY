# MICROHOM AUDIT V5 — SYNTHESIS (8-auditor redundant vet)

**Synthesizer role:** reconcile FOUR tasks × TWO independent Opus-Max auditors each into a robust
consensus on RECTIFY's **positional-distinctiveness CLOSE** for the motif-blind junction re-placer's
microhomology-drift guard.

**Under audit:** `_positional_signal(genome, q, q_split, ne, new_je, W=28)` +
`_semiglobal_ed(query, ref)` in `rectify/core/splice/junction_refiner.py` (lines 580-629), wired as
`drift_positional_gate` (lines 915-930, default 0.0 = OFF). A drift-flagged would-be-veto is SPARED when
`_positional_signal >= gate`. CLAIM under audit (dev/DISCOVERY_LOSS_PANEL_RESULT.md CP4/CP5): the
read-blind microhomology fault is **CLOSED** — WIRED m3/cap2/gate1 gives ~0.4% discovery loss / ~4.3%
fab-residual, ed_signal separates the delta overlap band at 98-99% balanced accuracy.

**Evidence provenance (source of truth = the 8 durable .md on disk):**
- signal-correctness_A — **STALLED at CP0** (16 lines, no findings, no verdict). Contributes nothing.
- signal-correctness_B — **STALLED at CP5** (138 lines, no final verdict line) but with complete
  unit-level probe results (CP2-CP4) + decisive reachability. Direction is unambiguous.
- independent-remeasure_A — COMPLETE, verdict HOLD (construction-dependent).
- independent-remeasure_B — COMPLETE, verdict HOLD (construction-dependent).
- byte-identity-architecture_A — COMPLETE, verdict CLEAR (safety) + fault_found=TRUE (architecture).
- byte-identity-architecture_B — COMPLETE, verdict CLEAR (safety), fault_found=FALSE (calls it debt).
- strategic-realdata_A — COMPLETE, verdict HOLD (enabling claim unproven).
- strategic-realdata_B — COMPLETE, verdict HOLD (enabling claim unproven).

---

## THE UNIFYING FACT — one mechanism, four vantage points: the W=28 horizon

Every auditor hit the **same mechanism** from a different angle. `_positional_signal` can only see the
base that distinguishes a real cryptic from a fabrication if that base sits **within W=28 bp** of the
acceptor. Outside W, both reference windows `genome[ne:ne+34]` and `genome[new_je:new_je+34]` become
identical (or the discriminating base is truncated away) → `_positional_signal` returns **0** → `0 < gate`
→ the read is NOT spared → the veto fires → real discovery is LOST.

| Auditor | Angle | Same fact restated |
| --- | --- | --- |
| signal-correctness B (B7) | sign axis | "blind past W — conservative, NOT a wrong sign" (sig=0 ≠ wrong sign) |
| independent-remeasure A | span sweep | knee at tandem-array span ≈ 28; 55-100% discovery loss past it |
| independent-remeasure B | length sweep | knee at shared-microhom length L ≈ 24-28; 14-62% loss, large denom |
| strategic A | distinguishing-base **density** | signal ≈ ±d; at d=1 (paralog PSV) ~2% loss at shipped gate1 |
| strategic B | distinguishing-base **position** | in-window SNP → esig=+2 (robust); SNP at je+40 (>W) → esig=0 → vetoed |

**This is why `close_is_correct=TRUE` and `close_is_complete=FALSE` are CONSISTENT, not contradictory.**
sig=0 on long microhomology is not a wrong *sign* (correctness passes — 0 is a legitimate "no in-window
evidence" value), but it *is* the discovery-loss *mechanism* (completeness fails — the read carried
decisive full-length evidence that the fixed window discarded).

**This also dissolves a false numeric conflict.** strategic's "~2% loss at d=1" and remeasure's
"14-100% loss at span≥28" are **two points on ONE curve**, not conflicting measurements:
- strategic sweeps *density within the window* (d = # distinguishing bases inside W). Low d, in-window.
- remeasure sweeps *distance to the window edge* (how far the divergence is pushed out). Divergence
  reaching/exceeding W.
Both are the horizon. As you push the discriminating information away from the junction (longer
microhomology) OR thin it out (fewer PSVs), the signal degrades along the same axis.

---

## PER-TASK CONSENSUS

### Task 1 — signal-correctness → CLOSE IS CORRECT (on the sign axis), single-auditor caveat

**A/B agreement:** cannot be a two-auditor consensus — A died at CP0. **B's evidence affirmatively
supports correctness:**
- `_semiglobal_ed` primitive: **4000/4000** random cross-checks agree with an independent brute-force
  min-over-ref-prefixes Levenshtein; hard-anchor at ref[0] and free-suffix both confirmed. No fault.
- `_positional_signal` sign: TRUE-POSITIVE (read==moved) → +18 spare; TRUE-NEGATIVE (read==incumbent) →
  −16 no-spare; ambiguous → 0. **Signs correct in the normal case.**
- **Minus strand: geometry CORRECT** — verified with a real pysam is_reverse read. Forward-genome frame
  (pysam q, CIGAR, ne all forward) makes the ed comparison strand-consistent; docstring "acceptor"
  wording is biologically loose on minus strand (ne is the donor there) but produces **no sign error**.

**Two latent unit-level wrong-sign bugs — both UNREACHABLE on real refs (B's CP4 is decisive):**
- B5 (genome-end truncation of the incumbent window → spurious spare) and B6 (rescue off the moved end →
  false veto) are deterministic wrong signs at the unit level, but require a junction within ~84 bp of a
  contig end. On the shipped Scer reference **min(len − intron_end) = 3929 bp**; **0/370 junctions within
  84 bp**; smallest nuclear contig 230 kb. → LATENT / theoretical only; worth a cheap defensive clamp
  (return None when the ref window truncates below rescue length), non-blocking.
- The both-boundary donor blind spot (B's probe d) is BY DESIGN (acceptor-only; docstring argues a donor
  term would hurt genuine cryptics) — a real hole only if both-boundary moves reach the gate, which B was
  mid-way proving does not happen at small δ (the two windows collide → tie → move not selected).

**Reconciled finding:** the signal is **implemented correctly** — no sign fault, no edge bug reachable on
chromosome-scale references, minus-strand correct. The W=28 blindness is NOT a sign fault (it is the
completeness fault, Task 2). **CAVEAT: this rests on ONE partial auditor** (A dead, B stalled before
writing its verdict line). Correct on the merits; recommend a re-run for two-auditor parity, and do not
launder absence-of-verdict into a "clean consensus pass."

### Task 2 — independent-remeasure → NOT COMPLETE: 0.4% is construction-dependent

**A/B agreement: STRONG, independent, converging.** Two auditors built two DIFFERENT constructions
(A: repeat-block + controllable-complexity downstream context + focused microsatellite panel; B:
period-p tandem microsatellite spine swept across W). **Both first reproduced the claim to validate
their harness, then broke it at the same knee:**

| | reproduce (easy regime) | break regime | knee |
| --- | --- | --- | --- |
| A | 0.0%/4.3% + 98% band-sep ✓ | disc-loss 7.2% agg → **55-100%** at span≥28 | span ≈ 28, sharp |
| B | 0.7%/0.0% (also re-ran the ORIGINAL panel → 98% ✓) | disc-loss **14-62%** at L~24-28, large denom | L ≈ 24-28 |

- **Bias confirmed decisively:** the original `EXON2 = U(k≤10) + U'(≤10) + random TAIL` caps contiguous
  microhomology at k≤10 and appends a maximally-distinctive random tail *immediately* after — the
  discriminating divergence always sits deep INSIDE W=28. The regime that breaks the signal (contiguous
  microhomology ≥ 28 bp) is **structurally excluded** from the panel that "proves the close."
- **The break is real discovery, not a vanishing denominator:** B's genuine-discovery check (L=26, 10%
  err, n=200) shows 70% of cryptics discoverable guard-off, 18% LOST by the shipped close. A's hand case
  (period-6 copies-6, span-30): full-read ed 0 (cryptic) vs 6 (incumbent) — **decisive evidence in the
  read** — yet `_positional_signal = 0` because the distinctive tail is at read-offset 30, past W=28.
- **Generalizes across unit periods incl. homopolymers** (B's period-generalization: p1/HP DLclose
  18-20% & fab-resid 75-88% once the HP run is long — the close FAILS on the very HP case that motivated
  the guard). Error-rate-independent (the signal is a *difference* of two eds vs the same corrupted
  rescue → common-mode error cancels; the blinding is pure window-truncation geometry).
- **Fails BOTH ways at intermediate span** (A's CP3): span-16 fab-residual 50% — inside the partially
  ambiguous window, error-driven fabrications score esig>0 and are wrongly SPARED. Not a safe monotone
  guard in the microhomology regime.

**Reconciled finding:** "0.4% discovery loss" is a **property of short-repeat microhomology only** — the
close holds where the divergence is inside W=28 (exactly what the original panel exclusively samples) and
degrades to 14-100% discovery loss where contiguous microhomology reaches the window length (biologically
ubiquitous: (CAG)n, poly-pyrimidine tracts, (AT)n, long HP). The "fault CLOSED" claim is **construction-
dependent, not robust.** Both auditors: HOLD on the *claim language*; code is correct + default-off so
non-blocking, but CP4/CP5's headline must be downscoped.

### Task 3 — byte-identity-architecture → CODE SAFETY CLEAR; architecture is design debt

**A/B agreement on the SAFETY property: complete and independent.** Both proved byte-identity three ways
(static diff / trip-wire / real-BAM SHA256), both green pytest, both confirm unwired:
- Real-BAM diff parent-vs-HEAD **identical SHA256** at default (block bypassed) AND when the veto block is
  forced live (margin8) — A: `e81d03a8…` / `033c0976…`; B: `3d93cd22…` (raw) / `e976505c…` (sorted,
  incl. parallel n_workers=2 fork path).
- Trip-wire: `_positional_signal` / `_semiglobal_ed` called **0×** at default (even with drift+cap on),
  called **1×/2×** only under a positive control (gate=1.0) — proving the 0-counts are real nulls.
- Unwired: `correct_command.py:746` passes NONE of the 4 drift kwargs (B adds: NOR `motif_blind`, default
  False — the *entire* motif-blind arm is dark code); no CLI flag, no config key, zero non-test caller.
- pytest -m "not slow": **1670 passed**; the 1 ERROR (`test_restore_cat3_plus_2`) is a pre-existing
  missing-fixture-file that references the guard 0× and **fails identically on the pre-close parent** →
  CLOSE-independent, orthogonal.

**The A(true)/B(false) split is a LABELING difference, NOT a substantive disagreement** — I explicitly do
NOT trip the "unresolved A/B disagreement = NOT-clear" rule on it. Both agree on every fact:
- byte-identical, green, unwired, inert at default;
- `_semiglobal_ed` is a flat-cost (HP-blind) near-clone of the scorer's HP-aware `_score_hp_anchored`
  (same "left ref fixed, free right suffix" contract);
- t1(k=0) vs `_semiglobal_ed` **disagree 34% of the time** (102/300) on HP-adjacent acceptors (structural,
  not literal, redundancy);
- **the discrimination BELONGS IN THE SCORER**: expose a k=0 hard-anchored distinctiveness term from
  `_score_junction` (reusing the existing HP-aware primitive, keeping the load-bearing free-k score
  alongside) so `delta_improve` separates on its own axis → likely makes BOTH `drift_positional_gate` AND
  the second O(m·n) `_semiglobal_ed` alignment deletable.
The only difference: A calls this an architecture **fault** (band-aid outside a scorer whose blind spot it
leaves intact, redundant with a primitive the scorer already contains); B calls it design **debt**. Same
recommendation. **The safety property is CLEAR from both; the architecture finding feeds worth_enabling,
it does not block the code as-is.**

**Reconciled finding:** the code is **safe to carry (default-off, byte-identical, unwired, inert, no test
regression) — there is nothing to block; you cannot HOLD an inert change.** The architecture point (the
discrimination should live in the scorer, HP-aware, not as a 4th/5th stacked move-gate parameter + a
cruder second alignment) is real and feeds the enablement decision.

### Task 4 — strategic-realdata → HOLD: enabling claim NOT validated on anything real

**A/B agreement: STRONG, independent, converging on four independently-sufficient strategic gaps.** Both
verified the ed_signal is NOT technically fragile (B: adversarial tight-paralog null — a single in-window
SNP holds esig≈+2 despite 6% error because the difference cancels common-mode error) — so the HOLD is
**strategic (untested-where-it-counts), not a math defect.**

1. **The guard's TARGET metric was never measured on real data.** COMPASS §4b — the recovery/drift/
   inconclusive classification that IS the real-data fabrication rate — OOM'd (exit 137), then TIMED OUT
   (30 min; algorithmic O(142k × per-chrom-list)), then was declared "context, not a gate." The single
   number the close exists to reduce has NO real-data value.
2. **The close was never RUN on COMPASS.** The headline 32× recall win (raw-mm2 0.54% → arm-B 17.46% vs
   independent motif-agnostic BBMap short-read truth) is **arm-E (`hp_drift_margin=3.0`) = the HP guard,
   NOT the microhom margin + cap + positional gate** (RESULT.md L176: "microhom guard NOT yet applied").
   `refined.bam` exercises neither the precision nor the recall side of the feature under audit.
3. **The close guards a feature that cannot be turned on** (dark code — see Task 3).
4. **All close evidence is one synthetic panel whose signal is manufactured by a property real target
   loci lack.** The panel's TAIL is RANDOM (~0.75/base divergence starting at offset k) → strong esig.
   Real paralog / rDNA / SMA microsatellite loci (the fabrication-enriched set) are near-identical tandems
   where distinctiveness is SPARSE and often OUT-OF-WINDOW (>28 bp) → esig weak-to-zero at exactly the
   hardest loci. A's d-sweep quantifies it: the panel's own at-risk cryptics have **ndisc median = 7, only
   9.6% at d≤2** — it structurally undersamples the high-mh-AND-low-d paralog regime (SMN1/SMN2 ~1 PSV /
   1-2 kb) the guard TARGETS. And real paralog origin is **label-unknowable** to single-read ONT (6%
   error, no phasing) — untestable by any ground-truth panel (B's CP5).

**Foot-gun flag (A):** the panel presents **gate=2** as a "balanced knob" (6.2%/2.1% — looks like a mild
trade). A's end-to-end sweep: at d=1 paralog loci gate=2 costs **~45% discovery loss**. Do not expose
gate=2 as a user knob without a loud warning.

**Reconciled finding:** the acceptor close is a legitimate, well-engineered, byte-identical-off piece of
CODE at a legitimate stopping point — but its EFFICACY claim ("the fault CLOSED") is synthetic-only,
over-favored by construction, and the real-data gate (COMPASS with the guard ON) was never run. Both
auditors: HOLD on enablement / the "CLOSED" claim. For a discovery tool, adding 4 params + a foot-gun
operating point to suppress an UNMEASURED real fabrication rate is not net-positive until COMPASS shows
the real fab-removed AND confirms recall is held. "Honest stop (acceptor close done, default-OFF, gate on
COMPASS)" beats "keep building" — both strategic auditors strongly endorse.

---

## SYNTHESIZED VERDICT

### (a) close_is_correct — **YES** (signal implemented correctly; single-auditor caveat)
No sign fault, no edge bug reachable on chromosome-scale references, minus-strand geometry correct,
primitive verified 4000/4000. The W=28 blindness is NOT a sign fault. The two latent wrong-sign bugs
(genome-end truncation) are unreachable on real refs (min 3929 bp vs 84 needed). **Caveat: rests on one
partial auditor (A dead at CP0, B stalled before its verdict line) — the affirmative evidence is solid but
a re-run for two-auditor parity is warranted. A cheap defensive clamp closes the latent B5/B6 bugs.**

### (b) close_is_complete — **NO** (0.4% is construction-dependent, does NOT robustly close on independent panels)
Two independent auditors, two different constructions, same knee at span/L ≈ 28 bp. Both reproduced
0.4%/4.3% FIRST (harness validated), then showed it is a property of the U+U'+random-TAIL construction
that caps microhomology at k≤10 inside the window. On tandem microsatellites / long HP the close degrades
to **14-100% discovery loss**, generalizing across unit periods including homopolymers, and admits
fabrication at intermediate span. Strong consensus: the fault does NOT robustly close.

### (c) worth_enabling — **GATED-ON-COMPASS** (NOT flat NO)
Split the two dispositions cleanly:
- **Carrying the code default-off: SAFE and at a legitimate stopping point.** Byte-identical, unwired,
  inert, no test regression. Nothing to block; you cannot HOLD an inert change. Mark the acceptor close
  DONE + default-OFF.
- **Enabling the guard + presenting the fault as "CLOSED": BLOCKED.** Both strategic auditors HOLD; the
  efficacy claim is synthetic-only and over-favored by construction; the real-data gate was never run;
  both byte-identity auditors say the discrimination belongs in the scorer, not as a 4th/5th stacked gate
  parameter + a redundant second alignment. Not net-positive for a discovery tool on the current evidence.

**Do NOT read the HOLD as "revert the code."** No auditor said that. The concern blocks *enablement* and
the *"CLOSED" claim*, not the code carried default-off.

### (d) real_data_gate — what COMPASS must show before enabling
A **NEW** COMPASS run with the close **ENABLED** — `microhom_drift_margin=3.0 / drift_near_tie_cap=2.0 /
drift_positional_gate=1.0` — on **both arms**, NOT arm-E (`hp_drift_margin=3.0`, which is what `refined.bam`
actually ran). Specifically:
1. Fix §4b's O(142k × per-chrom-list) scan with a **position-binned index** so the recovery/drift/
   inconclusive classification COMPLETES (feasible; the fix is already specified).
2. Report the **real fabrication rate the guard removes** (guard-ON vs guard-OFF) — the number the close
   exists to reduce, which has no real-data value today.
3. Confirm the guard does **NOT drop** motif-agnostic-BBMap-confirmed non-canonical recall — directly
   rebutting the recall doc's *assumed* line-176 preservation claim (if the guard vetoes real low-d
   cryptics BBMap can see, guard-on recall could drop; that is exactly what §4b-with-guard-ON would
   measure and never did).
4. Run over the SMA/paralog real drift set: **SNRPN ~6 bp, UBA1 ~26 bp, PCBP2 ~28 bp** (the spike-in
   doc's non-HP drift set — precisely the near-W regime where Task 2 shows the signal degrades).

---

## FAULTS TO FIX (priority order)

1. **[DECISIVE] The W=28 blind spot.** `_positional_signal` returns 0 whenever the discriminating
   divergence is pushed to/past offset 28 (contiguous microhomology ≥ ~28 bp: (CAG)n, poly-pyrimidine
   tracts, (AT)n, long HP). Fix: adapt W to the detected repeat span (extend the window past a detected
   tandem run), OR stop treating sig=0 as "no evidence to spare" (fall through to margin/cap when both
   windows are identical rather than vetoing). Note this only shifts the knee; past the *scorer's own*
   ~W horizon nothing is discoverable — so pair it with (2).
2. **[ARCHITECTURE] Move the discrimination INTO `_score_junction`.** Expose a k=0 hard-anchored
   distinctiveness term (reusing the existing HP-aware `_score_hp_anchored`, keeping the load-bearing
   free-k score alongside) so `delta_improve` separates on its own axis. This likely **deletes both**
   `drift_positional_gate` AND the second O(m·n) `_semiglobal_ed` alignment. If the external metric is
   kept instead: make it **HP-aware** (call `_score_hp_anchored`, not the flat-cost `_semiglobal_ed`
   which disagrees 34% on HP-adjacent acceptors — the DRS default) and reconcile the window mismatch
   (scorer rescue=30 vs signal W=28).
3. **[LATENT] Defensive clamp for B5/B6** genome-end truncation: return None when the ref window truncates
   below rescue length. Unreachable on Scer today, cheap insurance for fragmented assemblies.
4. **[CLAIM] Downscope "the fault CLOSED"** in dev/DISCOVERY_LOSS_PANEL_RESULT.md CP4/CP5 to "closes for
   contiguous microhomology whose discriminating divergence lies within ~W=28 bp; degrades to 14-100%
   residual discovery loss for longer tandem microhomology / long HP." Add a **long-microsatellite / long-HP
   arm** to dev/discovery_loss_panel.py (tandem array with divergence swept past W) so the frontier lives
   in the repo's own validation, not just this audit.
5. **[FOOT-GUN] Do not expose gate=2** as a user knob without a loud warning — on the panel's numbers it
   looks like a mild 6.2%/2.1% trade but costs ~45% discovery at d=1 paralog loci.
6. **[ORTHOGONAL] Fix the pre-existing missing-fixture test** (`test_restore_cat3_plus_2`, missing
   `scripts/validation_data/rebuild_2026_05/trimmed/validation_reads_polya_trim_metadata.tsv`) so
   `pytest -m "not slow"` is a clean 0-error green. CLOSE-independent (fails identically on the parent).

---

## FINAL CALL

**The positional CLOSE is CORRECT (signal implemented correctly, no reachable sign/edge fault) but NOT
COMPLETE (the 0.4% discovery-loss result is construction-dependent — it holds only for short-repeat
microhomology inside the W=28 window and degrades to 14-100% loss on tandem microsatellites / long HP, the
exact regime the guard targets). The code is SAFE to carry default-OFF (byte-identical, unwired, inert, no
test regression) — nothing to block. But the "fault CLOSED" claim is OVERCLAIMED and enabling the guard is
NOT justified: all efficacy evidence is one synthetic panel over-favored by construction, and the
real-data gate (a close-ENABLED COMPASS run with §4b completed) was never run. Disposition: mark the
acceptor close DONE + default-OFF; downscope "CLOSED"; do the scorer-level fix (which likely deletes the
gate + the redundant alignment) rather than stacking a 5th parameter; enable ONLY after the specified
close-enabled COMPASS measurement over the SMA/paralog loci.**

Consistent with prior audit rounds V1/V2/V4 (all HOLD). The single new synthesis contribution: the four
tasks are one mechanism (the W=28 horizon), which is why correct-yet-incomplete is coherent, not
contradictory.
