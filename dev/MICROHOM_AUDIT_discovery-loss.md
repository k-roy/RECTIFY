# MICROHOMOLOGY-DRIFT GUARD — ADVERSARIAL AUDIT: DOES IT KILL REAL DISCOVERY?

Auditor: adversarial (discovery-loss lens). Branch: worktree-agent-a25a2c1e784ad37dc
Started: (see git log). Model: Opus 4.8. READ-ONLY audit (no edits/commits).

## THE LENS
The guard adds `microhom_drift_margin` (op point 8.0) extra evidence for any junction
move whose `_move_microhomology >= microhom_threshold` (0.5). Threshold 0.5 was set on
ONLY 2 real R3 cases — under-powered. RISK: real non-canonical (cryptic 3'SS) sometimes
sits next to a short direct repeat -> the true move has incidental microhomology >= 0.5
-> guard demands +8 evidence -> a REAL discovery is vetoed. Load-bearing question:
IS 0.5 TOO LOW (kills real discovery), and is there a threshold/margin that keeps fab~0
WITHOUT discovery loss, or is the tradeoff irreducible?

## PLAN
1. Read `_move_microhomology` / `_frac_match` / guard application in junction_refiner.py.
   Understand EXACTLY what the guard gates on (which moves, what evidence, what margin).
2. Read design doc + existing tests + benchmark harness (build_panel.py, gen_reads.py).
3. Reproduce the shipped validation (fab FDR, R3 discovery) to establish baseline.
4. BUILD adversarial case: REAL non-canonical junction whose move has microhom >= 0.5
   (cryptic 3'SS adjacent to a short direct repeat, non-HP). Refine w/ guard (m=8) vs
   off. Does it SUPPRESS the real discovery?
5. SWEEP: broader R3/cryptic panel at varied microhomology_frac. Quantify discovery-loss
   rate vs microhomology_frac. Find the microhom_frac at which real discovery starts dying.
6. Verdict: is 0.5 too low? Is there a margin/threshold combo with fab~0 AND no discovery
   loss, or an irreducible tradeoff?

## CHECKPOINTS (appended as I go)

---
### CHECKPOINT 1 — orientation complete (guard mechanism fully mapped)

GUARD MECHANISM (junction_refiner.py:765-798): a MOVE (incumbent (ns,ne) -> best (new_js,new_je))
is vetoed (kept at incumbent) iff BOTH:
  (a) `_move_microhomology(genome, ns,ne, new_js,new_je) >= microhom_threshold` (0.5), AND
  (b) `best_score_cmp > incumbent_score - microhom_drift_margin` (8.0)  [i.e. move doesn't beat
      incumbent by >= 8.0 evidence units].
=> DISCOVERY DIES iff: the READ is from a REAL cryptic non-canonical transcript, its incumbent
   N-op is at the WRONG (canonical/nearby) site, the refiner would MOVE it to the true cryptic
   site, that move is flagged mh>=0.5, AND the evidence delta (incumbent_score - cryptic_score)
   is < 8.0.

SCORER (junction_scoring.py:983 _score_junction): scores only the FIRST 30 bp of exon2
(`_MAX_RESCUE=30`), intron_end-anchored bilateral DP (t1+t2). Units ~edit-distance
(`_CANONICAL_HP_PRIOR`=0.5). => distinctive exon2 bases BEYOND 30bp do NOT contribute evidence.
Margin 8.0 = 16x the canonical prior.

LOAD-BEARING EMPIRICAL CLAIM (mh_characterize.py): "real cryptic = 0.33 microhomology" measured
on mix_r3b_out — a panel whose exon2 is `random.Random`. 0.33 ~ E[match] of a random 6-mer vs a
random neighbor (~1/4/base + noise). REAL cryptic 3'SS are NOT random: alternative/cryptic sites
frequently arise BECAUSE of short local direct repeats. This is the n=2 under-powering to attack.

ADVISOR STEER (critical): `_make_read` in the existing test builds a CANONICAL read (fills query
from GENOME at incumbent coords) -> wrong polarity. I MUST inject the CRYPTIC exon2 as the query
while placing the incumbent N-op at the canonical/wrong site, so the read genuinely WANTS to move.
Discriminating axis = repeat-length (local detector) vs anchor-length/distinctiveness (global
evidence). Sweep the MARGIN too, not just threshold.

NEXT: build adversarial harness (own read builder, cryptic exon2 injected), first check =
k=6 boundary repeat + distinctive 60bp cryptic exon2, margin=8.0, does the move survive?

---
### CHECKPOINT 2 — first case built; geometry subtlety found

First adversarial construction (k=6 perfect tandem repeat CTGCAG at the T<->C boundary,
60bp distinctive tail): mh(C->T)=1.0 BUT delta(inc-true)=0.0 and the read fails to discover T
even with the guard OFF (m=0). ROOT CAUSE: I built exon2 = R + R + tail with the true acceptor
T and canonical C BOTH inside a PERFECT tandem repeat, so placing the read at C vs T is
genuinely SEQUENCE-EQUIVALENT (both see R + tail). That is a TRUE tie, not a discoverable
cryptic site -> the guard holding here is arguably CORRECT (advisor point 3: extended/perfect
repeat => hold is defensible).

NEED: a construction where the boundary microhomology is >=0.5 (so the guard flags the move)
BUT the cryptic exon2 is genuinely DISTINGUISHABLE from the canonical placement (real evidence
delta > 0). The discriminating axis is SHORT boundary repeat (local mh>=0.5) + distinctive exon2
that DIFFERS at the placement. Reconsidering geometry: the perfect-tandem trap must be broken
by making the repeat PARTIAL (mh in [0.5,1.0)) and the tail immediately distinctive.

---
### CHECKPOINT 3 — scorer soft-clip discovery: WHY delta=0 (key mechanistic finding)

Instrumented _score_junction at C (incumbent) vs T (cryptic) for the partial-repeat cases:
BOTH score 0.0 for ALL mh values. ROOT CAUSE (mechanistic, important): the scorer's t1(k) does a
FREE-PREFIX soft-clip — rescue[k:] aligned to genome[acc:acc+buf]. The read's cryptic exon2 is
genome[T:], and genome[T:] shifted right by `drift` overlaps genome[C:] wherever a repeat exists.
So at C the scorer soft-clips the first `drift` bases (CTGCAG) and rescue[6:]=CAGCTATTATCC...
matches genome[C:]=CAGCTATTATCC... PERFECTLY -> score 0. The scorer genuinely cannot distinguish
C from T because the read IS consistent with both at the drift distance (that IS the microhomology).

=> With scores TIED (delta=0), the refiner's is_alt tie-break keeps the read at the INCUMBENT C
   REGARDLESS of the guard (guard-OFF m=0 already gives C(LOST)). So in the pure soft-clip-tie case
   the guard is IRRELEVANT: the motif-blind refiner already fails to discover T.

IMPLICATION for the audit: discovery-LOSS ATTRIBUTABLE TO THE GUARD requires a case where T scores
STRICTLY BETTER than C (delta > 0) with guard OFF -> refiner discovers T -> then turning the guard
ON vetoes it. Need to BREAK the soft-clip tie: the cryptic exon2 must contain a base within the
first `drift` bases (the soft-clipped region) OR a length/indel signal that makes T strictly win.
Reconsidering: the distinguishing evidence must live in a region the C-placement soft-clip canNOT
absorb. NEXT: construct T-strictly-wins case, confirm guard-OFF discovers, then test guard veto.
