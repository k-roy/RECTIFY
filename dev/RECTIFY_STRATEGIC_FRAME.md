# RECTIFY — strategic frame for a fresh-Director assessment (2026-07-01)

*Purpose: situate the native-aligner gate program inside the ENTIRE RECTIFY
package, incorporate the PI's redirect (novel-isoform discovery + dirty-pA + variant
simulation are under-tested and higher-value than paralogs), and hand a fresh Fable
Director agent a neutral mandate to assess — from multiple angles — whether we are
on the right track. The Director may conclude we are NOT; steelman that explicitly.*

Read alongside (do not duplicate): `dev/NATIVE_ALIGNER_OVERVIEW.md` (plain-language),
`dev/SIMULATION_BENCHMARK_SPEC.md` (benchmark spec), `dev/C1_DESIGN.md`..`dev/C6_DESIGN.md`,
`dev/C3_DESIGN.md`, `dev/DISCOVERY_TIEBREAK.md`, `dev/HANDOFF_ALIGNER_BENCHMARK.md`
(gate history + verdicts), `rectify/CLAUDE.md` (pipeline index).

---

## 1. The whole RECTIFY package (the frame)

RECTIFY is a per-read ONT (DRS / cDNA) transcript-correction pipeline:

```
trim → align (5-aligner PANEL: minimap2, gapmm2, uLTRA, deSALT, gmap)
     → correct (per aligner): indel/HP · junction refine (2H) · 5' Cat3 rescue
                              · 3' CPA walkback (poly-A) → per-aligner corrected TSV/BAM
     → merge/consensus (pick winning aligner per read: hp_edit_distance arbiter)
     → analyze (aggregate 5'/3', clustering, DESeq2)
```

The **native de-novo aligner** is a proposed 6th, *orthogonal* member — NOT a 6th
correlated placer, but a calibrated-likelihood realignment/arbitration layer that
runs downstream of the panel to recover what the panel's shared flat-affine,
quality-blind, motif-snapping bias flattens away.

**Two distinct improvement surfaces (the PI's two assessment angles):**
- **(1) Can the trim / correct / consensus steps be improved?** (walkback, junction
  refine, 5' rescue, the hp_edit_distance arbiter — the SHIPPING code.)
- **(2) What should the de-novo aligner exhibit** to maximize not just orthogonal
  diversity but the *discovery* of novel isoforms — alternate TSS, unannotated +
  non-canonical introns, cryptic pA sites that extend OR truncate a transcript —
  that may appear ONLY in specific mutants?

---

## 2. Honest state of the gate program — including where we UNDER-tested

Solid + committed (do not re-litigate): C1 (HP indel placement) confirmed+shipped;
the **Discovery tiebreaker bug fixed+shipped** (a real consensus bug the gate
surfaced: a snapped junction hid its shift-insertion in an intron-side blind spot of
`_count_junction_proximity_errors` and won via the tiebreaker); WS-1 confirmed the
RNA004 divergence is a **hot-read-tail** chemistry property (validates the injector's
3-layer design). Flat-Q null (per-base quality adds nothing over the error-type table).

**Where the earlier "gate complete / opportunities exhausted" conclusion was
PREMATURE (the PI's correction, now accepted):**
- **Novel-junction discovery was tested only narrowly.** The C3 junction probe
  *injected* a truth-site member and asked "does the arbiter pick it over minimap2's
  snap?" (yes) — it did NOT measure the thing the native aligner exists for: across a
  GRADED ladder of novelty, **how often does the REAL 5-aligner panel fail to produce
  the true site at all**, and how does that scale with read quality? "Isoform
  flattening" — the headline motivation — was never systematically sized. C5's "no
  panel-failure tail" used a crude uniform-error injection, not a novelty gradient, so
  it does not speak to this.
- **C2 (3'/CPA) was refuted only on CLEAN genomic-A drift** (walkback at ceiling). It
  was NEVER pressure-tested against **dirty tails**: empirically, yeast Dorado leaves
  an **adapter stub** on the pA tail and intra-tail errors, so a true A(20) reads as
  e.g. `AAAAATAAAAAGAAAAAAAAATC` (internal T/G = tail errors; trailing TC = adapter
  stub). Whether walkback finds the true CPA through that is untested → a CORRECT-step
  (angle 1) question, re-opened.
- **C6 (variant-aware) was deferred too fast** as "read-alone can't; needs a VCF which
  on sim == truth = circular." That objection was wrong: injecting a real CATALOGUE
  variant AND giving the aligner the catalog is a legitimate *mechanism* test — its
  real-world validity is simply that matched-DNA VCFs routinely exist. Simulation with
  constructed truth is *cleaner* than raw variant-transcriptome data.

The chicken-or-egg that makes simulation CENTRAL: novel isoforms (and mutant-specific
ones especially) are exactly what the incumbent aligners bias against reporting, so
**real data carries no ground truth for the junctions/ends they flatten.** Only
constructed-truth simulation can measure the blind spot. This is not a fallback; it is
the correct instrument.

---

## 3. The plan under the two angles

### Angle 1 — trim/correct/consensus improvements
- **Dirty-pA walkback robustness (C2-redux, empirical).** Build a dirty-tail injector
  from the empirical Dorado failure modes: (i) intra-tail A→non-A substitutions at an
  empirical rate, (ii) a trailing **adapter-stub** of empirical length/composition,
  (iii) graded by read quality (the WS-1 hot-tail model). Truth = the constructed CPA.
  Measure walkback's CPA recovery across the dirtiness gradient; where it breaks =
  a correct-step fix (e.g. guard tuning), NOT an aligner facet.
- **Discovery tiebreaker fix** — shipped; needs a real-DRS recall spot-check before
  production reliance.
- Open the floor: any other correct-step gap the Director surfaces (junction-refine,
  5' Cat3 rescue interactions, the arbiter itself).

### Angle 2 — de-novo aligner features (simulation-central)
- **Novel-junction discovery ladder** (graded, error-overlaid):
  - Rung 1: multi-intron gene, move one intron's 5′/3′ SS to a nearby *unannotated but
    canonical* GT-AG (annotation-snapping, motif held).
  - Rung 2: non-canonical sites at graded deviation from GT-AG (GC-AG → AT-AC → 1-off
    → 2-off), parameterized by consensus distance.
  - Rung 3: rungs 1–2 on short internal exons AND **5′ terminal exons near a TSS**
    (short-anchor drop/soft-clip; overlaps existing Cat3 rescue → attribute correctly).
  - Rung 4: two independent shifted introns in one transcript (compounding).
  - Overlay: the WS-1-validated error injector, graded clean-bulk → low-Q, so recovery
    is stratified by read quality.
  - Metric: the **blind-spot surface** — real 5-aligner panel native recovery of the
    true site (ambiguity-aware) over (deviation × exon-size × error). Guards:
    ADDRESSABILITY (does the true site strictly win on a motif-blind empirical-penalty
    score → recoverable) vs ZERO-EVIDENCE (deep deviation + high error → read can't
    distinguish → honest ceiling); and FDR (no over-call on clean/annotated reads).
- **Cryptic pA / alternate 3′ ends** (extend OR truncate): simulate cryptic CPA sites
  that lengthen/shorten a transcript vs the annotated end; test discovery + FDR (the
  discovery analog of the dirty-pA correct-step test).
- **Alternate TSS / 5′ ends**: 5′-terminal-exon variants (ties to Cat3 rescue).
- **Variant-aware realignment**: sample coding variants from a catalogue (gnomAD/
  ClinVar) at realistic SNV:indel ratio, enrich the near-splice-site drivers, keep a
  far-from-splice control; inject into constructed reads. Arms: variant-BLIND
  (fabricates spurious intron) vs variant-AWARE (recovers true placement given the
  catalog); specificity control (must not suppress real near-variant novel junctions).

The unifying discovery target (PI): isoforms that appear ONLY in specific mutants
(alt TSS, unannotated+non-canonical introns, cryptic extending/truncating pA) — which
is exactly why constructed-truth simulation must be central, and why "diversity" is
necessary but not sufficient: the member must be *evidence-weighing*, not motif/
annotation-snapping.

---

## 4. Empirical error models (the realism anchors — already partly built)
- **Hot-read tail** (WS-1 validated): RNA004 = clean bulk ~1% + a ~1-2% tail at ~12%;
  the injector's per-read over-dispersion + burst-HMM reproduce this. Magnitude knob
  targets the tail, not the decent bulk.
- **HP/indel profile**: del-dominant empirical penalty table (C1).
- **Dirty-pA-tail model** (to build): adapter-stub + intra-tail error from the yeast
  Dorado observation.
- **Variant spectrum** (to build): catalogue-sampled coding variation.
All strata score against CONSTRUCTED truth; every claim ambiguity-aware; fitness is
the truth set, never an internal score (the standing discipline).

---

## 5. Mandate for the fresh Fable Director (the assessment)
Explore the ACTUAL code (trim, correct/walkback, junction refine, 5' Cat3 rescue, the
hp_edit_distance consensus, the aligner panel) and this plan from MULTIPLE angles, and
assess — neutrally, steelmanning "we are on the WRONG track" — :
1. **Angle 1**: where can trim/correct/consensus most improve? Is dirty-pA walkback
   the right first target, or is there a higher-leverage correct-step gap?
2. **Angle 2**: what features should the de-novo aligner exhibit to maximize
   NOVEL-ISOFORM discovery (alt TSS, unannotated+non-canonical introns, cryptic
   extending/truncating pA) — especially mutant-specific isoforms — beyond mere
   orthogonal diversity? Is the graded-ladder + variant-injection + dirty-pA
   simulation the right instrument, and is anything missing (e.g. isoform-level vs
   junction-level truth, quantification bias, chimera/trans-splicing, RT artifacts)?
3. **Cross-cutting**: are the gate verdicts we banked (C1 ship, C2/C3 refute, C4/C5/C6
   defer, Discovery fix) still correct under this wider frame, or did any refute/defer
   test too narrow a slice (as C2-clean-only and the novel-junction case already did)?
4. Prioritization: sequence Angle-1 vs Angle-2 work by leverage × risk × cost, with
   simulation central and the budget/concurrency lessons (≤3 heavy agents, no nested
   panel fan-out, write-to-files reliability) respected.
Deliverable: a strategic assessment written to `dev/DIRECTOR_ASSESSMENT_FABLE.md`
(are we on the right track? refinements? risks? the recommended first increment),
plus a short report. Do NOT commit; report back for integration.
