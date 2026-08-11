# ADVERSARIAL REVIEW — DIRECTOR_ALGO_EVAL_SYNTHESIS.md
*Reviewer: Sonnet (independent, adversarial). EARLY DRAFT — do not act on until final.*
*Date: 2026-07-01. Report-back only; no commits.*

---

## Executive summary

Three attacks LAND with varying severity; the synthesis's top-3 ordering largely
survives with one demotion. Details follow.

---

## ATTACK 1 — The ablation is NOT co-equal with the blindspot surface (LANDS)

**The claim:** The `_CANONICAL_HP_PRIOR` / `is_novel` double-prior ablation is ranked
as "a co-equal first experiment" with the blindspot surface, and the two "share one
harness."

**The steelman against:**

The corrector (`junction_refiner.py`) is a POST-ALIGNMENT step. It operates on reads
the aligner has already placed. If minimap2 snaps a CA-TC junction to GT-AG during
alignment (the 90% blindspot case), the corrector receives a canonical N-op — it has
no non-canonical junction to evaluate, no `tier_beats_alt=True` branch fires, and the
`canonical_discount` of 0.5 is never applied. The ablation is structurally invisible
for the large-blindspot population.

More precisely: the ablation's leverage is bounded by the aligner's native recovery
rate. For the deep-deviation rung (90% blindspot), only 10% of reads arrive at the
corrector with the correct non-canonical junction. If the corrector re-snaps half of
those, the ablation recovers 5% additional reads. That is SECONDARY leverage, not
co-equal.

The tier-gate (`tier < 4 else 0.0`) makes this worse, not better: at deep deviation
the CANDIDATE canonical junctions would have tier < 4, so the discount DOES fire on
the surviving 10% — but only on that 10%. The synthesis uses the tier-gate to argue
the designers "already disabled it where deviation is large," but this misreads the
code: the tier in line 647 is the tier of the CANDIDATE junction, not the tier of the
TRUE non-canonical junction. The discount fires on canonical candidates regardless of
how deep the true non-canonical deviation is.

**What the synthesis gets right:** "shared harness" is technically accurate — the
same ladder reads go through the aligner AND the corrector. You CAN measure both in
one pass. But "co-equal" implies comparable leverage.

**Verdict:** LANDS. The ablation is a SECONDARY measurement, bounded to the
non-snapped fraction (~10–53% depending on rung). Running it in the same harness is
efficient; elevating it to "co-equal first experiment" is an overclaim. Correct
framing: "measure aligner recall first; within the recovered (non-snapped) reads,
measure corrector re-snapping as a secondary experiment."

---

## ATTACK 2 — Same-model non-independence (LANDS, acknowledged)

**The claim:** The synthesis opens with a corrected provenance note: "three
independent same-model (Opus 4.8) runs." It then presents 11 "CONSENSUS FINDINGS"
and explicitly states items 1, 4, 5, 10, 11 are "largely two faithful restatements
of the program's own priors (low surprise — agreement is expected, not
triangulated)."

**The steelman against:**

The synthesis labels these 11 items "CONSENSUS FINDINGS" despite acknowledging that 5
of them are doc echoes and that all three runs read the same docs on the same model.
The "consensus" framing still imparts the rhetorical weight of independent
confirmation to items that are autocorrelated by construction. A reader scanning the
§1 header ("CONSENSUS FINDINGS — where the independent assessments agree") gets a
false impression of triangulation for the full set.

The genuinely new, non-doc-echoed findings are: Report B's code trace
(`_CANONICAL_HP_PRIOR` path, §2) and Report A's NM≥1/NM==0 operationalization (§3).
Everything else had already been asserted in `NATIVE_ALIGNER_OVERVIEW.md` and the
SPEC before the evaluations ran.

**Verdict:** LANDS. The synthesis acknowledges the problem but does not fully quarantine
it. The fix is to split §1 into: "DOC-ECHOED PRIORS (low-information; listed for
completeness)" and "INDEPENDENTLY COMPUTED FINDINGS (high-information)." Items 2, 3,
6, 7, 8, 9 may have independent value; items 1, 4, 5, 10, 11 do not. Using "CONSENSUS"
as the section header for all 11 conflates these.

---

## ATTACK 3 — "47% is a floor" overclaims; <60% threshold underdetermined (LANDS, partial)

**The claim (NOVEL_JUNCTION_BLINDSPOT.md, adopted by synthesis):** "AT-AC is only
~25% of U12-type minor-spliceosome introns; the MAJORITY of U12 introns are GT-AG in
sequence and are therefore INVISIBLE to this ladder. So the measured 47% AT-AC
blind-spot UNDER-states the real minor-spliceosome flattening, not overstates it."

**The steelman against:**

This confuses spliceosome-type identification with junction coordinate placement. The
majority of U12 introns that are GT-AG in sequence are CORRECTLY PLACED by minimap2
at the right coordinates. For RECTIFY's mission (novel isoform discovery by junction
coordinate), placing a GT-AG U12 junction correctly IS a success — the coordinate is
right, even though the spliceosome type is unlabeled. The ladder's control rung
(GT-AG: 98.3% recovery) captures this: the panel handles GT-AG U12 introns
identically to GT-AG major-spliceosome introns. The "invisible to this ladder" framing
misleads: they're invisible because they're NOT a problem for coordinate recovery.

The 47% AT-AC blindspot is exactly the problem for the minority of U12 introns that
ARE AT-AC. The synthesis's "floor" argument would hold only if RECTIFY's mission
included spliceosome-type annotation, which it does not.

For the <60% threshold: the document attributes it to "both Sonnet + Opus reviewers
converge" without derivation. 60% is a round number. The entire build/no-build
decision for the motif-blind realigner (Candidate #2) turns on whether the
five-aligner panel clears this threshold. If deSALT adds 25 pp to minimap2's 53%
AT-AC recovery (reaching ~78%), the gate closes — but 60% vs 70% would give opposite
answers with the same data. No FDR model, no power calculation, no sensitivity
analysis anchors it.

**Verdict:** PARTIAL LANDS. The "floor" claim is biologically confused and should be
removed or clarified. The <60% threshold is the pivotal go/no-go number and needs
either a principled derivation or explicit acknowledgment that it is a stake-in-the-
ground judgment call. As written it reads like a derived result when it is a prior.

---

## ATTACK 4 — Reject-list: pangenome/VG too hastily dismissed? (FAILS)

**The claim:** WFA/strobemer/genome-wide pangenome/r-index are "dependency-violating"
and rejected.

**The steelman against:**

For SMN1/SMN2 specifically, a local variation graph over the chr5 SMN locus is the
natural tool — it's what HapSMA and ctyper effectively do. The synthesis simultaneously
recommends "build the local 2-copy graph only if truth demands it," which IS a
pangenome approach for the locus. The "genome-wide" qualifier does real work here.
Local-locus VG is NOT dependency-violating if it's scoped to 500 kb and never invokes
the genome-wide vg toolkit.

BUT: the synthesis's routing is correct in practice. It says test winnowmap2 + abPOA
first; build local graph only if they fail. The rejection of "genome-wide pangenome" is
a scope/dependency argument that applies to the global case. The local case is
allowed under Candidate #5.

**Verdict:** FAILS. The reject-list applies to genome-wide pangenomes; the local
locus graph for SMN is explicitly preserved as Candidate #5. No overclaim here.

---

## ATTACK 5 — "Measurement → member" as indefinite deferral (FAILS)

**The claim:** The build-order principle ("measurement → member, never member →
search for a win") risks infinite deferral.

**The steelman against:**

C1 was confirmed and shipped; C2/C3 were rejected quickly; the blindspot surface now
returns a quantitative build signal (47-90% flattening). The program HAS shipped a
win (C1) and HAS closed gates (C2, C3). The synthesis explicitly sets a threshold
(panel AT-AC recovery < 60%) at which the gate passes. That is a falsifiable,
executable decision criterion, not a deferral strategy.

**Verdict:** FAILS. The threshold is stated; there is a next concrete step (run the
full five-aligner panel on the ladder); C1 shipped. The measurement-first principle is
operationalized here, not used to defer.

---

## VERDICT: What survives, what needs correction

**TOP-3 ORDERING — MOSTLY SURVIVES:**

1. Blindspot surface measurement SURVIVES as #1. The data is real and compelling.
   The ablation should be DEMOTED from "co-equal" to "secondary measurement within
   the same harness run" — it runs in the same pass but has secondary leverage.

2. Motif-blind empirical realignment SURVIVES as #2, properly gated on step 1.
   Survives intact.

3. Variant-aware track (cis-variant / SMN2) SURVIVES as #3. No attack lands here.

**THREE CORRECTIONS NEEDED:**

A. In §7 (top-3 description): demote the ablation from "co-equal first" to
   "secondary within the same harness." Recommended text: "measure native panel
   recall; within the non-snapped recovered reads, measure corrector re-snapping
   (ON/OFF/evidence-scaled) as a secondary probe."

B. In §1: relabel the section header to distinguish doc-echoed priors (items 1,4,5,
   10,11) from independently computed findings. "CONSENSUS FINDINGS" is misleading
   for items that are autocorrelated by construction.

C. In §3 (single-source additions) and NOVEL_JUNCTION_BLINDSPOT.md: remove or
   qualify the "47% is a floor" claim. The floor argument conflates spliceosome-type
   identification with junction coordinate placement. The 47% IS the problem scope,
   not an understatement of it. Also: either derive the <60% threshold from a model
   or flag it explicitly as a judgment stake-in-the-ground, not a convergent result.

---

---

## POST-ADVISOR ADJUSTMENTS (2026-07-01)

**Attack 1 tightening:** "Structurally invisible / disjoint populations" was too
strong. The refiner does search candidates within a window; if the aligner snap lands
inside that window, disabling the prior could recover the true non-canonical site.
Populations are not provably disjoint; leverage is bounded and secondary (the correct
claim). Softened above. Second prong added: even after ablating `_CANONICAL_HP_PRIOR`
and `is_novel`, the `tier` tiebreaker remains at priority 2/3 in the sort tuple
(lines 660, 663) — it still biases canonical. So the ablation is INCOMPLETE as well
as low-leverage. That makes the "co-equal" ranking doubly wrong.

**Attack 3a elevated:** The "47% is a floor" conflation is a substantive logical error
that propagated through the same-model chain uncritically. Elevated to co-lead with
Attack 1.

**Attack 2 downgraded:** The synthesis pre-concedes the same-model problem explicitly.
Listing it as a hard LAND is padding. Downgraded to "pre-conceded; residual is
editorial."

**Final ranking of attacks:**
- LEAD (both LAND): Attack 1 (ablation not co-equal) + Attack 3a (floor conflation)
- SECONDARY (partial): Attack 3b (<60% threshold underdetermined)
- EDITORIAL (pre-conceded): Attack 2 (same-model echo)
- FAILS: Attacks 4, 5

*Adversarial reviewer: Sonnet 4.6, independent, REPORT-BACK ONLY — no code changes.*
