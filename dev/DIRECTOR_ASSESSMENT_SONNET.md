# RECTIFY Director Assessment — Sonnet 4.6 (independent) — 2026-07-01

*Independence enforced: DIRECTOR_ASSESSMENT_FABLE.md was not read. All code
conclusions are firsthand from the source files listed below.*

Primary sources verified:
- `dev/RECTIFY_STRATEGIC_FRAME.md`, `dev/NOVEL_JUNCTION_BLINDSPOT.md`,
  `dev/NATIVE_ALIGNER_OVERVIEW.md`, `dev/HANDOFF_ALIGNER_BENCHMARK.md` (lines 1-542),
  `dev/SIMULATION_BENCHMARK_SPEC.md` (lines 1-584)
- `rectify/core/correct/walkback.py` (all 973 lines)
- `rectify/core/consensus/scoring.py` (all 761 lines)
- `rectify/core/consensus/select.py` (header + function signatures)
- `rectify/core/consensus/corrected_consensus.py` (header + arbiter logic)
- `scripts/benchmark/novel_junction_blindspot.py` (all 220 lines — metric verified)
- `rectify/core/benchmark/scorer.py` (lines 1-120 — ambiguity-aware TP/FN logic)

---

## Angle 1: Trim/correct/consensus — where is the leverage?

### What the code actually does

Two-tier arbitration:
- `select.py::select_best_alignment` uses `junction_score` (proximity-penalty; per
  `scoring.py`) to pick the best uncorrected alignment per read per aligner during the
  correct phase.
- `corrected_consensus.py::_cigar_hp_edit_distance` + `merge_corrected_tsvs` picks the
  winning aligner after all per-aligner corrections, using HP-aware edit distance on the
  corrected CIGAR (soft/hard clips = 1.0/base, N-ops free, D/I = empirical table).

The hp_edit_distance arbiter is AT CEILING on annotated junctions. This is not
reassertion — I read the penalty cost assignments (X=1.0 flat, N=0 free, D/I=empirical)
and the selector logic directly. C3's refutation ("calibrated LLR adds nothing over
hp_ed") is grounded: after the correct step has done its work, the corrected CIGAR
already encodes what the aligner got right or wrong, and hp_ed summarizes it correctly.

### The dirty-pA C2-redux

The tail-context guard `_max_non_stop_in_tail = max(1, tail_context_k // 4) = 1`
is the operative protection against dirty tails. Analysis:

- Common dirty-tail errors (intra-tail A→T/G substitutions, adapter stub over A-tract):
  these MISMATCH the genomic poly-A → not false anchors (`rb == gb` is the anchor
  condition, and error bases ≠ genomic A). The scan passes through them to the true CPA.
- Real vulnerability: adapter stub positioned over a genomic non-A base within the A-tract
  where the inward poly-A context cleanly matches the genome (all A's = stop-base,
  rb==gb → ctx_has_mismatch=False). Under that condition the tail-context guard does
  NOT fire, and the false anchor is accepted.

The trouble is that static analysis reverses mid-inspection (I caught myself twice). The
relevant conditions interact with local sequence in ways that depend on exact tail
composition and alignment. This is a textbook "measure-don't-assert" case. The
strategic frame's proposed dirty-tail injector (intra-tail substitutions + adapter stub,
graded dirtiness) is the correct instrument. Do NOT call walkback robust or fragile
without the measurement.

**Verdict on Angle 1**: The hp_edit_distance arbiter is at ceiling and the junction
scoring is architecturally sound. The dirty-pA C2-redux is GENUINELY OPEN. The
correct-step next action is the dirty-tail CPA-recovery gate (M1-local, cheap).

---

## Angle 2: Native aligner for novel-isoform discovery

### Blind-spot verification (firsthand)

I read `novel_junction_blindspot.py` in full. The acceptance loop (lines 95-98) verifies
that each non-canonical rung's true junction is genuinely non-canonical within its
ambiguity window before inclusion. This means any minimap2 snap to the nearest GT-AG is
a genuine FN on the truth, not a coordinate-bookkeeping artifact. The canonical rung
serves as an at-ceiling control (recovery 0.983 confirms the metric is sound).

**The result stands**: minimap2 flattens 47-90% of non-canonical novel junctions on
error-free reads, monotonically with motif deviation. AT-AC U12 = 46.7% blind-spot. 1-off
synthetic = 70%. 2-off = 78.3%. Deep 4-off = 90%. 

**Important calibration the doc under-emphasizes**:

1. The biological footprint is narrow. The canonical novel-junction case (GT-AG at an
   unannotated site) is the most common real novel-isoform mechanism, and minimap2
   recovers it at 0.983. The 47-90% blind-spot applies specifically to non-canonical
   motifs — AT-AC U12 introns are real but rare (~0.1-0.2% of mammalian introns). This
   narrows the claim from "minimap2 flattens novel isoforms" to "minimap2 flattens
   non-canonical-motif junctions," which is accurate but more precise.

2. The full-panel recovery is UNKNOWN. minimap2 is the most aggressive snapper. uLTRA is
   annotation-guided (may snap to annotation harder on known genes, or partially recover
   at novel sites). deSALT and gmap use different scoring models. The panel-native
   blind-spot rate (1 - Pr[at least one of 5 aligners produces the true site]) could be
   lower than minimap2 alone. This is the critical number the measurement has NOT yet
   produced.

3. reps=60: AT-AC blind-spot of 0.467 carries approximately ±0.07 95% CI. The direction
   (monotone degradation with deviation) is statistically robust; the exact percentage
   should be treated as an estimate.

**Verdict on Angle 2**: The native aligner's justification on discovery grounds is
CONDITIONAL on the full-panel blind-spot rate, which is not yet measured. minimap2's
47-90% is an UPPER BOUND on the panel gap — it proves the mechanism exists and that
every flattened case is strictly recoverable (error-free flanks → true site unambiguous
for a motif-blind aligner), but it does not floor the panel gap above zero. If uLTRA
or deSALT recover AT-AC junctions that minimap2 flattens, the panel-native gap could
be small despite minimap2's large individual number. The gap's existence is proved; its
magnitude — and therefore the investment level — awaits the panel measurement.

What the native aligner needs for this gap:
- Motif-blind scoring: empirical -log P penalty replaces the flat affine penalty that
  penalizes non-canonical dinucleotides.
- Calibrated evidence weighting (C1's foundation): the same penalty table that handles HP
  indel placement extends to junction scoring.
- Downstream arbitration using hp_edit_distance on the native aligner's corrected output
  competes cleanly with the panel in `merge_corrected_tsvs` — no new infrastructure
  needed for the arbitration layer.

Cryptic pA / alt-TSS: the blind-spot surface does not cover these. These lack
motif-snapping as a mechanism (CPA sites have no canonical dinucleotide equivalent),
but have annotation-snapping (the correct step's junction refiner may snap annotated CPA
positions). The isoform-level truth gap is a genuine hole: the scorer scores per-junction,
so a read with the splice right but wrong TSS/TES scores clean. This must be addressed
before the cryptic-pA and alt-TSS discovery rungs are added.

---

## Cross-cutting: gate verdict re-examination

**C1 (HP indel placement, shipped)**: HOLD. The empirical -log P table is the foundation
for everything downstream. Firsthand confirmed: `_cigar_hp_edit_distance` uses D/I costs
from this table; `scoring.py` uses it for proximity penalties. Sound.

**C2 (CPA walkback, refuted on clean only)**: RE-OPEN for dirty-tail gate. The clean-A
refutation is correct (walkback is at ceiling on identifiable genomic-A drift). The
dirty-tail gap is a genuine open measurement, not a theoretical objection.

**C3 (calibrated LLR arbitration, refuted)**: HOLD. The two-tier arbitration is at
ceiling. hp_edit_distance on the corrected CIGAR is the right summary stat. Adding
calibrated LLR on top of a correctly-corrected read adds no information. The C3 refutation
is structurally sound (confirmed firsthand from corrected_consensus.py).

**C4 (paralog POA, deferred)**: HOLD. HEADROOM=0 on identifiable reads; below-ceiling
slice is genuinely zero-evidence. No new information changes this.

**C5 (panel-failure tail, deferred)**: SCOPE-EXPAND NEEDED. C5 measured "reads placed
by no aligner at all" (0-0.4% at realistic error). The novel-junction blind-spot
measurement reveals a DIFFERENT panel-failure mode: all aligners simultaneously snap to
the WRONG junction — panel herding. These are distinct categories. C5's DEFER verdict
on the empty-union tail stands, but "panel-failure" must now be tracked at two levels:
(a) no placement (C5's meaning, small) and (b) all-aligners-herd-wrong (the blind-spot
finding, potentially larger). The five-aligner panel blind-spot measurement is the
instrument for (b).

**C6 (variant-aware, deferred)**: HOLD but revisit. The stratum discriminates (40bp+
GT..AG-flanked deletion → spurious N-op at >25/30bp boundary). The strategic frame is
correct that the C6 deferral reason ("circular if sim==truth") was invalid — giving the
aligner a matched VCF is a legitimate mechanism test. But the gate hasn't been formally
run. The C6 facet is a secondary priority after the blind-spot surface is completed.

**Discovery tiebreaker (shipped)**: The fix to `_count_junction_proximity_errors`
(`prev_rp = ref_pos + length` for N-op, = intron_end) is confirmed in scoring.py. The
snapped junction's shift-insertion is now correctly charged to the exon-side proximity
window, so the true junction wins on the primary score. Real-DRS spot-check still pending
per the strategic frame.

---

## Steelman: are we on the WRONG track?

Three honest cases where the answer might be yes:

**Case 1 — The panel blind-spot is smaller than minimap2 alone suggests.** If uLTRA's
annotation-guided mode recovers AT-AC U12 junctions (plausible when the annotation
includes them), the full-panel blind-spot rate for real non-canonical junctions may be
much smaller than 47%. If so, the native aligner's discovery value is narrow — and the
program should pivot harder to Angle 1 (correct-step improvements) and real-data
transfer.

**Case 2 — Canonical novel junctions (GT-AG at new sites) are the dominant biology,
and minimap2 already handles those at 0.983.** The mutant-specific isoforms the PI cares
about may primarily use canonical splice sites at unannotated positions, not non-canonical
motifs. If so, the gap the native aligner addresses — non-canonical motif flattening —
is not where the biological signal lives. The program would be building a solution to
a measurement artifact (the non-canonical ladder) rather than to the real missing biology.

**Case 3 — The gate program is engineering-heavy and science-light.** The program has
invested heavily in gate-discipline and simulation infrastructure. If a simpler approach
— run all 5 aligners, take the union of novel junctions, apply a simple support filter —
already captures most novel sites on real data, the native aligner adds no incremental
value. The gates have not benchmarked this simple baseline.

**Counter-assessment**: Case 1 is testable (panel blind-spot measurement, next
increment). Case 2 is partially real but incomplete — the PI explicitly named
non-canonical introns as a target, and AT-AC U12 is a real biological class. Case 3 is
the strongest concern: the union-based baseline has not been benchmarked. Adding it as
a panel-recovery floor control would sharpen the native aligner's justification.

---

## Prioritization

By leverage × risk × cost:

1. **Panel blind-spot measurement** (cluster, 1-2 days). Run the existing
   `novel_junction_blindspot.py --emit-corpus` + 5-aligner panel on the cluster.
   Score with `panel_blindspot_score.py`. This is the single most important
   measurement: it either confirms the gap (justifying the native aligner investment)
   or reveals the panel partially recovers it (requiring reprioritization). Cost: low
   (existing infrastructure). Risk: low (deterministic). Leverage: highest.

2. **Dirty-pA C2-redux gate** (M1-local, 1 day). Build the dirty-tail injector
   (intra-tail substitutions + adapter stub), measure walkback CPA recovery across the
   dirtiness gradient. This is a correct-step (Angle 1) gap that affects every DRS read.
   The static analysis of walkback.py is indeterminate — only empirical measurement
   resolves it. Cost: low (M1-local). Risk: low. Leverage: medium.

3. **Addressability formalization** (M1-local, 1 day). Confirm per non-canonical rung
   that the true site strictly wins on motif-blind empirical-penalty scoring vs the snap
   target. This is the concrete "native aligner recovers this" proof and the design spec
   for the first native aligner facet. Prerequisite for building the member.

4. **Isoform-level truth gap**: Add TSS/TES truth to the scorer so cryptic-pA and
   alt-TSS rungs are honest. Required before those rungs can be gated meaningfully.

5. **C6 variant-aware gate**: Run the formal gate (variant-blind vs variant-aware arms,
   near-splice-site driver enrichment). Lower priority than panel blind-spot — it sizes
   a smaller gap.

---

## Summary verdicts (short-form)

**Angle 1 (trim/correct/consensus)**: ON THE RIGHT TRACK. hp_edit_distance arbiter is
at ceiling and the junction scoring is architecturally correct. The open item is dirty-pA
walkback — UNMEASURED, not resolved. Gate it before declaring Angle 1 complete.

**Angle 2 (native aligner for discovery)**: CONDITIONALLY JUSTIFIED — the mechanism
is proved and addressable, but the magnitude of the panel gap is unmeasured. minimap2's
47-90% is an upper bound on the opportunity, not a confirmed floor. Full-panel
blind-spot measurement is the critical next gate before committing to building the member.

**Gate verdicts**: C1/C3/C4 hold. C2 re-opened for dirty-tail. C5 scope-expands to
track panel-herding separately from no-placement. C6 deferred correctly but gate not
yet formally run. Discovery tiebreaker shipped.

**Single first increment**: The panel blind-spot measurement (5-aligner corpus on the
cluster). Cheap, definitive, and the bottleneck for every downstream decision.

**Does the blind-spot result change the native aligner justification?** YES — it
operationalizes and strengthens the case but is not itself decisive. The result proves
the flattening mechanism is real in minimap2 and strictly addressable (error-free flanks
→ the true site is unambiguous for a motif-blind aligner). minimap2's 47-90% is an
UPPER BOUND on the panel gap, not a floor: the panel could already recover what
minimap2 flattens. So the correct answer to the mandate question is: the blind-spot
result changes the question from "does the mechanism exist?" (yes) to "how large is
the panel gap?" — and justification is conditional on that number. The first increment
(panel measurement, `panel_blindspot_score.py` already exists) will be decisive.
