# Novel-junction blind-spot surface — the discriminating measurement (2026-07-01)

The gate the Fable Director (`dev/DIRECTOR_ASSESSMENT_FABLE.md`) identified as the ONE
that discriminates the entire "is the native aligner justified?" question — and that
the program had never made. Prior C-gates tested ARBITRATION ("given a truth member,
does the consensus pick it?"). This tests the thing the native aligner exists for:
**how often the panel fails to produce the true novel junction at all** — the
isoform-flattening rate — as a function of splice-motif deviation from GT-AG.

Injector-FREE by design (error-free reads): minimap2 snaps non-canonical junctions
even on error-free reads (the overview's own claim), so this measures PURE
motif-snapping with ZERO dependence on the placeholder-magnitude error injector (the
program's standing bear case), and it transfers to real data immediately.

Script: `scripts/benchmark/novel_junction_blindspot.py` (result:
`novel_junction_blindspot_result.txt`). Generalizes `gen_junction_discovery_stratum`
into a graded deviation ladder; scored by the shipped ambiguity-aware scorer
(`score_bam`); the canonical rung is the at-ceiling control; addressability is clean by
construction (error-free exact-match exon flanks ⇒ the true site is strictly
recoverable by a motif-blind aligner, so non-recovery is snapping bias, not evidence
loss). No shared-file edits; smoke GREEN.

## RESULT — minimap2, error-free, reps=60

| rung | motif | dev(GT-AG) | real motif | native recovery | **blind-spot** | snap |
| --- | --- | --- | --- | --- | --- | --- |
| control | GT-AG | 0 | Y | 0.983 | 0.017 | 0 |
| U12 | AT-AC | 2 | **Y** | 0.533 | **0.467** | 11 |
| 1-off | GA-AG | 1 | n | 0.300 | **0.700** | 35 |
| 2-off | CT-AC | 2 | n | 0.217 | **0.783** | 16 |
| deep | CA-TC | 4 | n | 0.100 | **0.900** | 26 |

(GC-AG dropped out: the scorer correctly treats GC..AG as canonical, so it is not a
non-canonical discovery case — i.e. NOT a blind spot. That is the right answer; GC-AG
is in the recognized set.)

## VERDICT — BUILD SIGNAL (the native aligner IS justified on discovery grounds)

**minimap2 natively flattens 47–90% of non-canonical novel junctions, monotonically
with motif deviation — including ~47% of REAL AT-AC U12 minor-spliceosome junctions —
on ERROR-FREE reads.** The canonical control recovers at ceiling (0.983), so the
metric is sound. Because the reads are error-free with clean flanks, every flattened
call is *strictly recoverable* by an evidence-weighing, motif-blind member → this is
precisely the addressable regime the native de-novo aligner targets. This is the
first measurement in the program that directly demonstrates the headline motivation
(isoform flattening), and it is **injector-independent and real-data-transferable** —
it does NOT ride the placeholder injector magnitude that blocks the other verdicts.

This **reverses the earlier "opportunities exhausted" framing** decisively: the panel
does flatten novel isoforms, at a large and graded rate, and the gap is addressable.

## Honest caveats + the next increments
1. **minimap2-only.** This is the canonical snapper and represents the panel's shared
   flat-affine bias, but the true panel-native recovery is "≥1 of the 5 aligners gets
   it." uLTRA/gmap are annotation-guided (may snap to annotation harder OR recover
   annotated novels); deSALT is another engine. **NEXT (cluster): run the full
   5-aligner panel on this ladder** and report panel-native recovery — the honest
   flattening number. (minimap2 alone flattening 47–90% already sets the direction.)
2. **Snap-target opportunism.** The random core means the canonical motif minimap2
   snaps TO arises opportunistically in the flanking/core sequence (as in real
   genomes). The RECOVERY number (does minimap2 call the TRUE site) is robust to this;
   the `snap` count is construction-dependent and reported as context only.
3. **Add the remaining ladder axes** (per the strategic frame): exon-size (short
   internal + 5′-terminal-exon-near-TSS — ties to Cat3 rescue), multi-intron
   (compounding), then the error overlay (which re-introduces the injector-calibration
   caveat — do it AFTER the clean surface + panel confirmation).
4. **Addressability formalization.** Confirm, per rung, that the true site strictly
   wins on a motif-blind empirical-penalty score vs the snap target (the C1/hp_penalty
   −logP scale) — the concrete "a native member recovers this" proof, and the design
   spec for that member.
5. **isoform-level truth gap** (Director): the scorer scores per-junction; a read with
   the junction right but a wrong TSS/TES still scores clean, so isoform-discovery FDR
   is invisible. Needed before the cryptic-pA / alt-TSS discovery rungs.

## Cross-model correction (Sonnet outward Director, 2026-07-01) — the 47% is a FLOOR
Independent (verified-Sonnet) review adds a biology point that STRENGTHENS the case: AT-AC is only ~25%
of U12-type minor-spliceosome introns; the MAJORITY of U12 introns are GT-AG in sequence and are therefore
INVISIBLE to this ladder (they read as canonical). So the measured 47% AT-AC blind-spot UNDER-states the real
minor-spliceosome flattening, not overstates it. Also: minimap2 has TWO mechanistically distinct snapping
biases — `--splice-flank` (motif snapping) and `--junc-bonus`/`--junc-bed` (annotation snapping) — which earlier
descriptions conflated; the ladder (no annotation) isolates the motif-snapping mechanism. And a literature
check (ESPRESSO/Science Adv 2023 retains a canonical-motif gate; 2passtools/Genome Biol 2021 re-aligns with an
annotation guide; IsoQuant) found NO published long-read tool that removes the motif prior at the scoring level
or uses HP-law empirical deletion costs for junction scoring → the calibrated-DP-no-motif-prior member concept
is genuinely novel, not a reinvention. Decision threshold (both Sonnet + Opus reviewers converge): native
member JUSTIFIED if panel-native AT-AC recovery < ~60%; if deSALT/mapPacBio independently recover most AT-AC
(> ~80%), the gain is arbitration/union not a new placer → shift to correct-step + C4/C6.

## PANEL RESULT (job 32422876, COMPLETED 2026-07-01) — mapPacBio BREAKS the herd (plot-twist, verified)
Full-panel native recovery per rung (only minimap2 + mapPacBio produced BAMs — deSALT + uLTRA FAILED on the
synthetic per-contig corpus: "All requested aligners failed"; likely index-build/annotation setup on tiny
contigs, to be fixed + re-run):
| rung | minimap2 | mapPacBio | PANEL(union) | panel blind-spot |
| --- | --- | --- | --- | --- |
| GT-AG canon | 0.983 | 1.000 | 1.000 | 0.000 |
| AT-AC (U12) | 0.533 | 1.000 | 1.000 | 0.000 |
| GA-AG 1off | 0.300 | 0.983 | 0.983 | 0.017 |
| CT-AC 2off | 0.217 | 1.000 | 1.000 | 0.000 |
| CA-TC deep | 0.100 | 1.000 | 1.000 | 0.000 |
VERIFIED (not a scorer artifact): mapPacBio CIGAR on a CA-TC_deep read = `199=200N201=` — a genuine N-op
intron at the true coordinates (all 60 deep-rung reads carry the N-op). mapPacBio (BBMap, splice-aware) does
NOT snap non-canonical junctions to GT-AG; it places the intron on read evidence regardless of motif. The
per-contig construction is "easy," but minimap2 fails on the SAME easy setup (0.10-0.53), so mapPacBio's win
is attributable to its SCORING model (no GT-AG prior), not the construction.
IMPLICATION (pre-committed rule: panel-native recovery >~80% => gain is arbitration/union, not a new placer):
the panel ALREADY contains a non-snapping member (mapPacBio) that recovers the non-canonical novels minimap2
flattens. This WEAKENS the "native placer justified for non-canonical intron discovery" case on THIS
measurement, and shifts the question to: (a) does the CONSENSUS/arbiter actually PICK mapPacBio's true
placement over minimap2's snap end-to-end (the corrected-consensus + Discovery-fix should — TEST it); (b) does
the `_CANONICAL_HP_PRIOR` CORRECTOR then RE-SNAP mapPacBio's recovered non-canonical junction back to canonical
(the double-prior ablation is now the RELEVANT next test, run on mapPacBio-recovered reads); (c) does this hold
under READ ERROR + a realistic multi-decoy GENOME + the OTHER discovery targets (alt-TSS, cryptic pA,
variant-adjacent, SMN). CAVEATS: error-free, per-contig, 2-of-4 aligners, non-canonical-INTRON axis only.
Interpretation (build vs pivot) is advisor-gated; do not declare the native aligner unjustified on this alone.

## CORRECTION (adversarial Sonnet reviewer) — retract "47% is a FLOOR"
The earlier "Cross-model correction" claim that the 47% AT-AC blind-spot UNDER-states real U12 flattening is a
LOGICAL ERROR (conflates spliceosome-TYPE identification with COORDINATE placement). For RECTIFY's mission
(recovering junction COORDINATES), GT-AG-motif U12 introns are placed CORRECTLY (they sit at canonical
coordinates — captured by the 0.983 control rung); only the AT-AC minority is a coordinate problem, and 47% IS
its scope, not an understatement. The majority-U12-are-GT-AG argument makes the real coordinate-flattening rate
LOWER than 47%, not higher. => Treat 47% as the AT-AC-coordinate blindspot scope, NOT a floor on U12 flattening.

## PI CAVEAT (2026-07-01, load-bearing) — mapPacBio has MAJOR weaknesses on real human data
mapPacBio's ~100% non-canonical recovery is on ERROR-FREE, single-contig SYNTHETIC reads. Kevin (PI): mapPacBio
has major weaknesses on REAL HUMAN data — so its synthetic success does NOT transfer, and the native-aligner
case stands even where mapPacBio recovers the synthetic non-canonical sites we hand it. => the ERROR OVERLAY
(and, later, real-genome/real-read) test is DECISIVE, not confirmatory: if mapPacBio's advantage collapses under
RNA004-bulk error / real complexity, the calibrated-−logP native member is justified. Do NOT conclude "panel
already covers it" from the clean synthetic panel result; weight the error-overlay + real-data result far more.

## RECONCILIATION (Sumner human eval, 2026-07-02) — mapPacBio's synthetic "win" is its real-data PATHOLOGY
Real-human evidence (dev/SUMNER_HUMAN_ALIGNER_EVAL.md): mapPacBio recovers the synthetic non-canonical junction
for the SAME reason it fails on real human ONT-DRS — it has NO splice-model gate and emits any scored gap as an
intron. On clean synthetic reads that lands on the true site (looks like recovery); on real human data the same
behavior gives a 97.7% SPURIOUS-novel rate (SMA chr5) and consensus that SELECTS its artifacts. => the "panel
already covers it via mapPacBio" PIVOT is REFUTED by real data (mapPacBio's coverage is illusory — it emits true
AND false indiscriminately). The panel has NO member that discovers novel non-canonical junctions at acceptable
PRECISION (workhorses mm2/uLTRA/deSALT are annotation/GT-AG-biased → flatten novel; de-novo mapPacBio/gmap are
spurious-dominated → both effectively rejected on human). That precision-recall gap = the native member's target.
=> the native aligner IS JUSTIFIED. CRITICAL HARNESS GAP: the blindspot ladder measures RECOVERY only, which is
GAMED by indiscriminate emitters (mapPacBio "recovers" by emitting everything). NEXT BUILD (highest priority):
add an FDR/PRECISION axis — per-aligner spurious-non-canonical N-op rate + FP-junction rate on reads with NO
true novel junction (the over-call control, the fp_canonical_snap track applied to panel recovery). Only the
recovery+FDR pair reconciles synthetic vs real and correctly scores mapPacBio.
