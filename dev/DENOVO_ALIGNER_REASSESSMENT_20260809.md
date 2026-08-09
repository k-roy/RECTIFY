# REASSESSMENT — the de novo aligner "too gameable" verdict + the consensus-triage fold-in (2026-08-09)

**Question (PI):** revisit the standalone de novo aligner (the orthogonal panel complement). It was
previously deemed too prone to gaming our scoring system — reassess that conclusion. And: should its
ideas instead be folded into consensus selection as a **triage** layer (high-confidence alignments
bypass; distressed alignments — junction-proximal mismatch overhangs, 5'/3' soft clips — get targeted
Rectify re-alignment)?

**Sources:** `dev/ALIGNER_MEMBER_DESIGN.md`, `dev/HANDOFF_NATIVE_ALIGNER_BUILD.md`,
worktree `agent-a25a2c1e784ad37dc` docs (`HANDOFF_NATIVE_ALIGNER_VETTING.md`,
`ALIGNER_BENCH_STATE_AUDIT_20260721.md` — the definitive state audit, `HANDOFF_GUARD_REEVAL.md`),
`dev/SNAP_OR_HOLD_ANCHOR_COORD_FINDING.md`, and the shipped consensus code on the current branch.

---

## 1. The "gaming" verdict — CONFIRMED, but narrower than remembered

Four distinct gaming episodes are on the record, all real, all reproduced:

| # | What was gamed | Evidence | How it was closed |
|---|---|---|---|
| 1 | **Internal consensus score as fitness** | a pure re-weighting flipped GMAP annotated win/loss 0.09→1.07 with no aligner change | external ground-truth benchmark gate; "fitness = truth set, NEVER the internal score" |
| 2 | **−logP empirical table in the junction-boundary SEARCH** (arm-C) | over-shifts junctions into HP runs (ACC_A_D1: 73 reads moved 2 bp into the A-run; refined 442 reads vs arm-B's 164); del_open_delta variant REJECTED by an impossibility argument — a shift into a run is emission-identical to truth except the gap's label (free `N` vs penalized `D`), so the wrong placement is the global min **regardless of the table** | table confined to RANKING fixed candidates (consensus `hp_edit_distance`; C1-del exon-interior), never boundary search |
| 3 | **Per-cut `ins_cost` run-splitting** | the per-k truncation split homopolymers to exploit the non-monotonic DRS ins table (A×12: 1.76 vs 8.26) | full-run ins_cost, then mooted by the table-free concat-DP (14.3×, byte-identical) |
| 4 | **The junction-coordinate METRIC itself** (compensating-indel exploit) | the re-placer slid junction coordinates with matched I/D pairs while the read's exon sequence never moved; the "32× recall / 27% fabrication" story was ~85–95% this artifact (corrected: **~4% vs minimap2 0.5% ≈ 8×**); an 8-reviewer guard program was built against the artifact | one structural invariant — *a junction move must not raise the read's indel burden* (`e40ca00`, always-on) — reverted ~94% of moves, removed ~91% of apparent fabrication at zero real-discovery cost |

**Reassessed conclusion.** The verdict is *earned and should stand* — for the specific formulation
tried: **likelihood-guided junction-boundary search.** That avenue is closed by more than empirics
(the episode-2 degeneracy argument is essentially an impossibility theorem: as long as N-ops are free
relative to indels, any pure-likelihood boundary search is structurally exploitable). Do not
relitigate arm-C / del_open_delta.

But "too gameable" was **never a verdict against the de novo aligner's ideas.** Every episode was a
gaming of a *search* or a *metric*, and every fix was a **structural constraint, not a better score**
(truth-set fitness; motif-blind flat costs; the compensating-indel invariant — while score-based
hold-margins provably had *no sweet spot*). What survived that discipline is validated and real:

- **motif-blind re-placement + compensating-indel invariant** = the discovery engine:
  corrected real recall **~4% vs minimap2 0.5% (≈8×)** on Sumner chr5 DRS with independent
  short-read corroboration; yeast BY4742 DRS transfer 17 fix / 0 harm; human chr5 transfer
  +956 junctions moved onto annotations / 0 harm; Q1 Illumina adjudication 153 FIX / 0 HARM.
- **C1-del** (HP-length-law deletion cost) confirmed for exon-interior indel attribution
  (real-SIRV over-call 0.0000); C1-ins refuted and gated off.
- The 2026-07-21 state audit's own conclusion anticipates the PI's present question verbatim:
  *"panel-membership is the wrong frame — it manifested as a refiner."*

## 2. The strongest single argument FOR the triage fold-in: C3's refute

The C3 (calibrated-LLR arbiter) keystone was refuted with **headroom = 0**: given a truth-containing
candidate set, the *already-shipped* `hp_edit_distance` arbiter picks truth **100%** — even on
strata where the members disagree 100%. Translation: **the panel's failures are
missing-candidate failures, not arbitration failures.** A better judge buys nothing; what's missing
is a **generator of the missing candidate**. The PI's triage-and-re-align proposal is precisely
that generator, feeding the already-validated judge. This is the correct division of labor and the
core reason to endorse the fold-in.

## 3. Assessment of the proposed architecture — endorse, with three cautions

**Proposed flow (endorsed):** after per-aligner correction + `hp_edit_distance` ranking:

1. **High-confidence bypass** — clean junction-proximal overhangs, no unexplained 5'/3' soft clips,
   annotated junctions, low hp_ed → label `high_confidence`, skip. (Yeast WT: ~98.8% of junction
   reads are already correctly placed — the do-no-harm data say these are at ceiling; touching them
   is pure downside + compute.)
2. **Triage → targeted re-alignment** — junction-proximal mismatch/indel clusters, soft-clipped
   ends, short-anchor junctions → run the validated machinery: motif-blind refiner (internal
   junction placement), `align_clip_to_exon`/Cat3-rescue path (5' overhang placement), the standing
   `up_amb`/`down_amb` ambiguity-window flexing at annotated junctions (the open "future
   enhancement" in CLAUDE.md), C1-del table for exon-interior indel attribution only.
3. **Re-entry, not auto-accept** — the re-aligned candidate re-enters the hp_ed arbiter as one more
   member. The re-aligner *proposes*; the validated ranking *disposes*. This preserves the exact
   asymmetry that already works (ranking is safe; searching is the gameable act) and is the
   architectural anti-gaming control.

**Grounding:** every triage signal already exists in shipped code —
`scoring.py::_count_junction_proximity_errors` (:575), `get_softclip_lengths` + the 5'/3' clip
machinery (:401–:507), `corrected_consensus.py` hp_edit_distance primary sort. The triage layer is
mostly wiring, not new science. Triage also dissolves the standing perf blocker: the refiner at
~0.24–0.34 s/junction-read made full-depth genome-wide human infeasible (~83 hr/arm/sample);
concat-DP already gives 14.3×, and triaging to the distressed few-% multiplies that into
production-feasible territory.

**Caution 1 — the invisible failure class (the big one).** A flattened non-canonical junction looks
CLEAN: minimap2 snaps it to the nearby canonical motif with tidy overhangs and no soft clips (the
blindspot rises to 0.90 for deep non-canonical). Read-level distress triage will *bypass* exactly
the reads the discovery engine exists for. So keep **two distinct gates**:
- *accuracy triage* (read-level distress → re-align) as proposed; and
- *discovery triage* (pool-level): loci where a minority of reads shows junction-proximal distress,
  recurrent low-support alternate placements, or cross-sample recurrence → re-examine the whole
  pool with the motif-blind pass. Do not let the high-confidence bypass silently cap discovery.

**Caution 2 — keep the anti-gaming constraints on the re-align phase.** Compensating-indel
invariant always-on; flat costs (no −logP) in any boundary search; motif/annotation as soft
tie-breakers, never hard filters, in discovery mode; abstain by default. The residual *real* risk is
microhomology-rich drift (~1.2% pool-dependent FDR; the drift guards were SHELVED 2026-07-17 because
`_positional_signal` could not separate recovery from drift on real data — a coin flip). It is
bounded by triage scope + candidate re-entry + recurrence, not by resurrecting the guards.

**Caution 3 — the 3' end.** C2 (CPA change-point decoder) was refuted: the guarded walkback is at
ceiling and the decoder over-calls A-ending bodies. "3' soft-clip triage" should route to the
existing walkback/A-tract machinery, not a new 3' decoder.

## 4. Practical state + the honest gates

- **Everything validated is stranded on `worktree-agent-a25a2c1e784ad37dc`** (503 commits ahead of
  master; last commit 2026-07-17; never merged): `motif_blind`, the compensating-indel invariant
  (`e40ca00`), concat-DP (`e1ed90c`, default ON), the C1-del `penalty_table` knob, the shelved guard
  machinery (dormant, default 0.0 — candidate for removal in the port). First engineering step is a
  **surgical port of the keeper set** to a fresh branch off master.
- **Two open validations from the state audit remain the honest gates** before the triage layer
  ships as default behavior:
  1. *(cheap, M1-inline)* the `scoring.py::_count_junction_proximity_errors` post-N `prev_rp`
     surgical fix — in-probe validated to kill the canonical-snap residual with zero truth-recall
     cost, never applied. It also makes the triage signal itself more accurate at exactly the
     junction-proximal positions the triage keys on.
  2. *(cluster)* reproduce the ~8× corrected recall on a **second independent DRS corpus** with
     orthogonal short-read truth (data staged on Sherlock). If it doesn't reproduce, the discovery
     half of the triage layer stays research-only.
- **Do not relitigate:** arm-C / del_open_delta (−logP in boundary search), the C2 decoder, the C3
  LLR arbiter, Flat-Q, blunt hold-margins.

## Bottom line

The "too prone to gaming" conclusion stands **for likelihood-guided boundary search** — that door is
closed by an impossibility argument, not just bad runs. It does **not** condemn the de novo
aligner's ideas, most of which survived under structural (non-score) constraints and carry real-data
validation. The proposed fold-in — triage in consensus selection, high-confidence bypass, targeted
re-alignment, candidate re-entry under `hp_edit_distance` — is the right architecture: it is where
the program's own 2026-07-21 audit already landed ("it manifested as a refiner"), it converts the
validated-but-stranded worktree work into a shippable stage, it fixes the perf wall, and C3's
headroom-zero refute is direct evidence that a candidate-generator-plus-existing-judge is the
correct complement to the panel. The one design addition it needs is the pool-level discovery gate,
so the clean-looking flattened-junction class isn't bypassed out of existence.
