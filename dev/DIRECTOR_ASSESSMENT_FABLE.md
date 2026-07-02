# Director assessment (Fable) — RECTIFY native-aligner + correct-step program

*Fresh-Director strategic assessment. Render date 2026-07-01. Do NOT commit until
reviewed. Basis: code-grounded for walkback / scoring / select (read in full);
documentation-grounded for the green refiner/Cat3/Cat6 modules (CLAUDE.md +
HANDOFF); sim harness via a read-only Explore subagent + advisor at both forks.*

---

## Verdict: on the right track on PROCESS, off on SEQUENCING

The gate-first discipline ("prove against truth before building; fitness is the
truth set, never an internal score") is sound and has earned its keep — 4 plausible
directions refuted/deferred, 1 confirmed, a real consensus bug (the Discovery
tiebreaker) surfaced AND fixed. I confirmed the C3-refute and the Discovery fix by
reading the code, not just trusting the banked verdict. Keep the discipline.

**But the program sequenced its measurements wrong.** The ONE measurement that
discriminates the entire "is the native aligner justified?" question — the
**novel-junction blind-spot surface** (how often the REAL 5-aligner panel fails to
produce the true site at all) — was never made. Every C-gate to date tested
arbitration ("given a truth member, does the arbiter pick it?"), which is a
different question. The headline motivation (isoform flattening) is therefore
still **unmeasured**. That is the gap to close first, and it can pivot the program:
empty surface -> pivot to correct-step + transfer; large surface -> build signal.

Steelman "we are on the WRONG track": every verdict is sim-grounded on a
placeholder-magnitude injector (overdisp_v ~0.27 vs real 0.70), and real-data
transfer (Deliverable B) is BLOCKED — so it is fair to say "0 facets proven on
real data, 1 sim-proven shipped on principle." This is the strongest bear case.
It does NOT invalidate the process; it means the transfer loop is a first-class
dependency, and the decisive next measurement must be one that DOESN'T depend on
injector calibration. Fortunately one does (below).

---

## Angle 1 — trim / correct / consensus

**Highest-leverage correct-step target is NOT dirty-pA in isolation.** dirty-pA is
real and addressable but bounded, and my own walkback read predicts it re-refutes
on the naive stratum:
- The core scan anchors on a genome match; intra-tail A->non-A subs never trigger
  a stop (they aren't genome matches). So the UNIFORM-sub stratum shows ceiling and
  misleads exactly as C2-clean did.
- The ONLY mechanism that bites: `_max_non_stop_in_tail = max(1, k//4) = 1`
  (walkback.py:632). 2+ adjacent non-A in the k=4 window (a force-aligned adapter
  stub `...TC`, or a clustered intra-tail burst) defeats the tail-context guard ->
  the walk short-stops before the true CPA.
- Addressability is aligner-specific and already settled by the code comments:
  **mapPacBio force-aligns** past pA (the cat1_plus_1 `rb=A,T,A,A` pattern) so the
  stub reaches the scan; **minimap2 soft-clips** it (out of scan). So leverage is
  real but concentrated in force-aligning members.

So dirty-pA is worth ONE targeted gate — on the **clustered-≥2-non-A / adapter-stub
stratum only**, addressability-first — not a broad build. The green modules
(junction refiner, Cat3 5' rescue, hp_edit_distance arbiter) are documented
resolved; the one live correct-step enhancement (precomputed up_amb/down_amb for
soft-clip flexing) is bounded and lower-leverage than the blind-spot measurement.

## Angle 2 — de-novo aligner features (the real leverage)

The instrument (graded ladder + variant-injection + dirty-pA + cryptic-pA) is the
RIGHT design, but the ladder DOESN'T EXIST yet — only a single JUNCTION_DISCOVERY
2x2, no deviation gradient, no exon-size axis, no error stratification. The scorer
already measures native panel recovery (no truth-member injection), so the missing
piece is the driver stratum, not the scoring machinery.

**Missing beyond the ladder** (genuine gaps, in priority): dirty-pA + cryptic-pA
(extend/truncate) strata; alt-TSS / 5'-terminal-exon variation (ties to Cat3);
**isoform-level truth** (today per-junction + per-indel only — a read with all
junctions right but a wrong TSS/TES scores fine, so isoform-discovery FDR is
invisible); quantification bias (alt-TSS/TES); chimera / trans-splicing /
RT-template-switching (all absent — flag as scope, not necessarily build).

## Cross-cutting — what re-opens

- **C3-as-accuracy: STAYS refuted** (confirmed by reading select/scoring).
- **Flat-Q: STAYS null.** C1-del ship stands (with a real-data spot-check owed).
- **C2 and C6 do NOT re-open as standalone builds** — they become STRATA INSIDE
  the blind-spot measurement program (C2-dirty = clustered-stub CPA rung; C6 =
  variant-adjacent FP rung). The measurement tells you what to build; don't
  pre-commit a facet.
- The earlier "gate complete / opportunities exhausted" framing WAS premature (the
  PI's correction is right): novel-junction discovery was tested only via injected
  members, never as real-panel-failure sizing.

## Prioritization (leverage x risk x cost; simulation central)

1. **Angle-2 clean-rung blind-spot gate (FIRST — discriminates the whole program,
   injector-FREE).** Rungs 1-2 (unannotated-canonical, then GC-AG/AT-AC/1-off/2-off)
   on **error-free** reads. The overview's own claim is minimap2 snaps *even on
   error-free reads* -> this measures pure motif/annotation-snapping with ZERO
   dependence on the placeholder injector, and transfers to real data immediately.
   Addressability-clean (clean canonical evidence = strictly recoverable).
   Metric: real 5-panel native recovery of the true site (ambiguity-aware) over
   (deviation x exon-size). Guards: addressability (true site strictly wins on the
   motif-blind primary score) vs zero-evidence; FDR (no over-call on clean annotated
   reads). M1-light. THEN add error-overlay stratification (carries the calibration
   caveat).
2. **Dirty-pA clustered/stub gate (SECOND, cheap correct-step).** Addressability
   first (which members force-align the stub), stratum = clustered-≥2-non-A / stub,
   NOT uniform-sub. Truth = constructed CPA.
3. **Transfer-loop unblock (ENABLING dependency, parallel).** Unblock Deliverable B
   / lock the SIRV injector magnitude — required before ANY error-overlaid or
   magnitude-sensitive conclusion is trustworthy. Every new stratum carries a
   real-data spot-check plan.

Respect the budget lessons: <=3 heavy agents, no nested panel fan-out, director-run
inline gates (the c3/c4/c5 method never failed), write-to-files reliability.
