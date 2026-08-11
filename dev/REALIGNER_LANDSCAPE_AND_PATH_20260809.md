# The Rectify Re-aligner — landscape assessment + the durable path (2026-08-09)

**Question (PI):** the mapPacBio re-implementation became the "resolver" — it re-arbitrates the
minimap2 arm's soft clips and junctions, and pressure now exists to extend it to terminal CLEAN
ends. The original "better mapPacBio" intent has evolved into the Re-aligner program. Assess the
landscape; propose a durable path.

**Companion docs:** `dev/DENOVO_ALIGNER_REASSESSMENT_20260809.md` (why triage-and-re-align is the
right frame; the score-gaming history), `dev/ALIGNER_BENCH_STATE_AUDIT_20260721.md` (the native
re-aligner program's state of record), Chanfreau `planning/641/643/644` (the resolver SPEC,
implementation, and the T3–T8 acceptance campaign — CLOSED 2026-08-09 with all six verdicts).

---

## 1. The landscape: three lineages that converged on ONE design pattern

| lineage | station in the pipeline | what it touches | validation state |
|---|---|---|---|
| **Resolver** (`overhang_resolver`, 641/643/644) | pre-consensus, panel-arm level (fills mapPacBio's slot on the minimap2 arm) | terminal SOFT CLIPS only — re-places them across junctions from a splice-site index, under the information bound (`W_max = α·2^I_eff`; refusal first-class) | **644 campaign CLOSED (2026-08-09):** recovers 51% of mapPacBio's beyond-mm2 gold (58% of the alternative-intron gap) at **17× precision** (313 vs 5,456 junk) and **~2% of compute** (~52× arm-vs-arm; ~14× at production scale); T6 synthetic: information diagonal exact, poly(A)/repeat refusals 0-FP, decoy FP 2.8% (homologous copies — mappability mitigation named); T8: the shared gate withdraws 251 phantom rescues |
| **Native re-aligner / refiner** (Module 2H arc, worktree branch — now merged into `feat/realigner-triage`) | post-consensus, read level | INTERNAL junction placements on the consensus read — motif-blind re-placement inside the panel window; compensating-indel invariant always-on; concat-DP | corrected real recall ~4% vs minimap2 0.5% (≈8×) on Sumner chr5 DRS w/ SR corroboration; yeast DRS 17 fix/0 harm; human chr5 +956/0 harm; Illumina adjudication 153 FIX/0 HARM; −logP-in-search REFUTED (impossibility argument) |
| **Triage layer** (built 2026-08-09, `feat/realigner-triage`) | consensus orchestration | NOTHING directly — classifies (read-evidence only), routes distressed reads to legs, arbitrates re-entry by hp_ed | classifier + junction leg + strict re-entry gate: 11 tests green incl. end-to-end bundled smoke; policy thresholds pending corpus tuning |

**The convergence is real and it is the finding.** All three lineages independently arrived at the
same architecture the score-gaming history forced: **a constrained proposer + an evidence arbiter +
refusal as a first-class outcome.** The resolver's information bound, the refiner's
compensating-indel invariant and motif-blind flat costs, and the triage layer's
bypass-plus-re-entry are the same idea at three altitudes. That shared pattern — not any single
component — is the durable asset, and it is why "a better mapPacBio" stopped being a general-purpose
aligner: the general-purpose form (`maxindel=200000`, ~12,500 chance acceptors per window) IS the
unbounded search the pattern exists to forbid.

**What mapPacBio's residual value actually is** (measured, 644 T3/T4): (a) **overriding minimap2's
confident UNCLIPPED placements** (the SRC1 class — 112 reads); (b) **non-canonical dinucleotide
classes** (AT-AC etc., outside the resolver's canonical-by-construction index); (c) **whole-read
relocations**; (d) the 14 Gould-S6 novels no minimap2-based arm reaches. Everything else it did is
now done better, 17× cleaner and ~50× cheaper, by the resolver.

## 2. The clean-ends question — recommend AGAINST extending the resolver's contract

The pressure to let the resolver reconsider terminal **clean** ends is the SRC1 class knocking on
the wrong door. Three reasons to refuse:

1. **The resolver's entire validated safety case rests on its input contract.** It touches only
   soft clips — positions where the aligner has ALREADY declared uncertainty. Its α-budget FP
   model, its refusal semantics, and the 644 verdicts are all conditional on that contract.
   Re-arbitrating confident placements converts "resolve declared uncertainty" into "search against
   the aligner's confidence" — a different, strictly harder FDR problem, and precisely the shape of
   unconstrained search that arm-C/del_open taught us is structurally gameable.
2. **A confident-but-wrong placement is invisible read-by-read** — clean ends, tidy overhangs.
   Read-level machinery (resolver, triage classifier) cannot see it BY CONSTRUCTION. It is only
   visible at the **locus level**: 112 reads at SRC1 where another algorithm family disagrees;
   recurrence across samples. The right owner is the pool-level discovery gate (Station C below),
   which triggers on locus evidence and only then re-examines reads.
3. **The contested-GCR1 lesson cuts the other way too:** when the resolver's candidate DID compete
   at a gold locus, it won on sequence with 62-read support vs the catalogue's 4 — flagged as a
   possible real alternative 3′SS. Locus-level adjudication, not blanket clean-end search, is how
   that class gets decided.

## 3. The durable architecture: ONE Re-aligner, three stations, one contract

```
                     ┌──────────────────────────────────────────────────────────┐
  Station A          │ RESOLVER (pre-consensus, per-arm)                        │
  arm level          │   input: declared uncertainty (terminal soft clips)      │
                     │   tool: splice-site index + information bound            │
                     └──────────────────────────────────────────────────────────┘
                     ┌──────────────────────────────────────────────────────────┐
  Station B          │ TRIAGE + REFINER (post-consensus, per-read)              │
  read level         │   input: read-evidence distress (junction-proximal       │
                     │   errors, unexplained clips, unannotated junctions)      │
                     │   legs: motif-blind junction re-placement · Cat3/clip    │
                     │   rescue (shared informativeness gate) · 3' walkback     │
                     │   arbiter: strict hp_ed re-entry; bypass the clean ~98%  │
                     └──────────────────────────────────────────────────────────┘
                     ┌──────────────────────────────────────────────────────────┐
  Station C          │ DISCOVERY GATE (pool/locus level) — NOT YET BUILT        │
  locus level        │   triggers: minority distress, recurrent alternate       │
                     │   placements, cross-family confident disagreement        │
                     │   (the SRC1 detector), cross-sample recurrence           │
                     │   actions: targeted motif-blind re-examination of the    │
                     │   locus pool; optionally a BOUNDED full-realignment      │
                     │   probe (mapPacBio's residual role) on that locus only   │
                     └──────────────────────────────────────────────────────────┘

  Shared contract (all stations): constrained proposer · evidence arbiter ·
  refusal first-class · ambiguity-canonicalised scoring · compensating-indel
  invariant · NO likelihood-guided boundary search · no internal-score fitness.
```

**mapPacBio's endgame:** demote from always-on panel member to **Station-C probe** — invoked only
on pool-flagged loci and unresolved triage reads. That keeps its measured residual (confident-
disagreement discovery, non-canonical classes, relocations) while collapsing its ~80 core-h/panel
cost to the flagged slice. The pre-registered `[[617]]` leave-one-out on CORRECTED output remains
the formal adjudicator for removing it from the default panel — do not skip it.

## 3b. Clarification (PI question, 2026-08-09): probe ≠ arbiter, and the trigger-power problem

**Is the proposal "don't run mapPacBio on all reads — only on Station-C-flagged loci, as an
arbiter"?** Half yes, with one correction and one honest weakness:

- **Probe, not arbiter.** At a flagged locus mapPacBio GENERATES an alternative full-realignment
  candidate; it never judges. The arbiter stays the evidence layer (hp_ed re-entry,
  ambiguity-canonicalised scoring, recurrence/corroboration). An instrument whose beyond-minimap2
  junction set is 98.1% junk cannot adjudicate anything — but its rare genuine disagreements
  (SRC1: 112 reads) are exactly the candidates the arbiter needs to be SHOWN. Same asymmetry as
  every station: proposal is cheap to forgive, judgment is not.
- **The trigger-power problem (the weakness the question exposes).** The SRC1 class was DISCOVERED
  via cross-family disagreement — i.e. using mapPacBio's panel-wide output. In the strong form
  (mapPacBio never runs broadly), Station C loses its best trigger for the confident-but-wrong
  class and must rely on mapPacBio-independent signals: minority distress, recurrent low-support
  alternate placements, cross-sample recurrence, and a per-locus error-density anomaly (a spliced
  read forced THROUGH a real intron should carry an error-dense segment). Plausible, unmeasured.
- **The middle path — SCOUT MODE.** Run mapPacBio on a ~5–10% subsample (or only reads over loci
  with any weak distress) purely as a trigger source. Fabrication is per-read random and does not
  recur; genuine disagreement recurs at a locus — so locus-level trigger power survives
  subsampling at ~5% of the compute (the same logic as the Sumner 5% genome-wide first pass).
  Full-depth probe then runs only on flagged loci.
  **⚠ MEASURED CORRECTION (2026-08-10, `dev/STATIONC_MAPPACBIO_HARVEST_20260810.md`):** the
  recurrence statistics above hold only for RANDOM fabrication — 33% of mapPacBio's junk
  junctions recur (≥2 reads; systematic repeat/homology misalignment), so recurrence is a first
  gate (3× cut), never a sufficient one. Admission must be TWO-TRACK: canonical-class +
  support≥2 ≈ 48% precision (harvests half the canonical residual as a ~25-candidate
  shortlist); the non-canonical track (where the unique Gould-class discovery lives, 0/14
  canonical) stays ~0.3% under any support threshold and requires orthogonal evidence
  (short-read corroboration, cross-sample recurrence, co-localized mm2-side distress).
- **The decision stays empirical:** the pre-registered `[[617]]` leave-one-out on CORRECTED output
  decides whether mapPacBio's arm still contributes to the final product with resolver + triage in
  place. Nothing → out of the consensus panel, retained as scout + probe. Something → it stays,
  and Station C consumes its disagreements panel-wide.

## 4. The path (sequenced, falsifiable)

1. **Land the resolver** (641/643 owns it): rebase `feat/overhang-resolver-641` onto the updated
   line; before human deployment add the T6 decoy mitigation (uniqueness/mappability check).
   Follow-ups already named in 644: support-stratified re-census; adjudicate the GCR1 contested
   3′SS; optional AT-AC index arrays if U12 coverage is wanted.
2. **Land `feat/realigner-triage`** (this branch: native re-aligner merge + triage MVP) once the
   full suite is green on it. Then wire the triage clip legs: 5′ clips → the Cat3/rescue path
   (which now carries the shared informativeness gate from T8), terminal clips → the resolver
   machinery. One refusal discipline across both.
3. **Tune the triage policy on real corpora** (the classifier thresholds are placeholders):
   measure bypass fraction + leg yield on the BY4742 DRS and upf1Δ gold windows; the 644 gold
   catalogue is the ready-made truth set for triage-leg recall/precision.
4. **Build Station C** (the genuinely new work): the SRC1-class detector (cross-family confident
   disagreement census), minority-distress and recurrence triggers, locus-scoped re-examination,
   and the bounded mapPacBio probe hook. The open validation gate from the state audit lives here
   too: reproduce the ~8× discovery recall on a second independent DRS corpus.
5. **Run the `[[617]]` leave-one-out** on corrected output → decide mapPacBio's default-panel
   status on evidence.
6. **Retire the shelved drift-guard machinery** (dormant params) in a cleanup pass once the triage
   layer is the specificity control — one system, not two.

### 4b. Clip-leg wiring — the one-sided-anchor division of labor (641 proposal, 2026-08-11; ACCEPTED as the design)

641's 644f forensics (H2 14283115) closed the mechanism question on the resolver's 14-junction
Gould residual and handed the acceptance decision to this layer, with a measured case:

- All 14 are non-canonical junctions with exactly ONE canonical side (5 canonical-donor,
  9 canonical-acceptor); index_visible 0/14.
- On the 39 supporting reads, **minimap2 itself performs the canonical-preference flatten** the
  noncanon control measured in the rearb pass: 32/39 reads assert a canonical-in-class variant
  with ONE boundary shifted 1–9 bp (25 off by exactly 1 bp; chrXI:93365 = 11/11 reads uniformly
  off-by-1). 7/39 express the junction as a 14–22 bp terminal clip.
- Per-read the snap moves junction-proximal ED by only ~1–2 — under arb_margin (2.0), so no
  per-read contest can displace the canonical incumbent. **Recurrence across reads is the only
  decisive signal** — the "18 margin-beats need pool-level evidence" argument, now with a
  concrete 32-read, 13-locus instance.

**Division of labor (CORRECTED per Kevin 2026-08-11 — the governing spec is
`~/work/UCLA/Chanfreau_Lab/planning/644h_realigner_integration.md` §3; 641-impl has handed the
whole sequence to this agent):**
1. **Enumeration (resolver library, post-landing work):** per flagged junction test BOTH
   hypotheses — (A) donor held, acceptor scanned ±w; (B) acceptor held, donor scanned ±w —
   grammar-free on the scanned side, w ≈ 16 bp, scored with the existing bounded DP. The
   original one-sided proposal is superseded: **the incumbent does not reveal which side was
   snapped** (both its ends look canonical — that is what the snap did). Both-free excluded
   initially; A∪B contains every measured truth.
2. **Trigger (Kevin's criterion): side-split flank error density.** A wrong boundary
   concentrates mismatch/indel density in the first bases of ONE junction-adjacent flank; the
   pooled ed_cur misses 1-bp snaps (~1–2 ED, under the clean-skip threshold). Flag liberally —
   acceptance is downstream.
3. **Adjudication (triage layer):** per-read deltas are sub-margin by construction, so no
   per-read rule can take these without reopening the FP door. The resolver emits per-read
   evidence records (read, incumbent, hypothesis, best_alt, ed_delta, flank_density_L/R); the
   triage layer aggregates per locus and **accepts on pool recurrence** (chrXI:93365 = 11/11
   reads voting the same off-by-1) — same discipline as the Cat3/rescue plan and §3b's
   two-track admission. The 7 clipped reads route through the terminal-clip leg with the same
   two-sided candidate class. Phase-2 note: the 644h overhang-quality likelihood applies
   unchanged — candidates are admitted on (recurrence × short-side overhang quality), never on
   margin alone.

Also riding with the resolver branch: `124c84a` RECTIFY_SKIP_REGIONS (yeast-rdna shorthand;
resolver + rescue bypass for rDNA reads; 2.6× on chrXII, byte-identical non-rDNA output, zero
gold cost — the 10 rDNA Gould entries are 9,139-bp repeat-unit artifacts hit by no arm).

## 5. What NOT to do (the closed doors, so they stay closed)

- No likelihood/−logP-guided boundary search, in any station (impossibility argument; arm-C,
  del_open_delta, per-cut ins run-splitting all refuted).
- No blanket clean-end re-arbitration in the resolver (SRC1 belongs to Station C).
- No internal consensus score as a fitness function — truth sets and orthogonal corroboration only.
- No un-gated searches: every proposer keeps an explicit refusal path and a bounded window.
