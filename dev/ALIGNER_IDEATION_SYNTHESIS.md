# RECTIFY native-aligner — first-principles ideation synthesis (2026-06-18)

Source: 6-lens creativity panel + high-effort synthesis judge (workflow `w03wt9tmh`, run
`wf_bb59e428-529`; full result `tasks/w03wt9tmh.output`). Charter:
`dev/workflows/alignment_firstprinciples_ideation.js`. This is the **design-doc seed** — every concept
is a HYPOTHESIS PENDING ABLATION on the simulation benchmark, never a finding.

## PI localization-hypothesis verdict: **ADOPT_WITH_FALLBACK**
Reuse the 5 aligners' placement cluster as the localization window; native RECTIFY does **local
realignment** in that window (the entire novelty budget → refinement, where the 87-90% indel/residue
reality lives). Dependency-light, tractable now (`local_aligner.py` is a live 791-line affine-gap DP).
Operating downstream of placement diversity = an arbitration/refinement layer, NOT a 6th correlated
placer → friendly to consensus-collapse.

**Sharpened discovery ceiling (the key insight):** the cap is at the **WINDOW level, not the junction
level.** A read in the correct window but with a mis-called/panel-disagreed junction *inside* it IS
discoverable by native realignment — so GMAP-class novel junctions in touched windows are reachable.
Only reads with **NO acceptable window** (all-5-misplace, GMAP-unmapped ~12%, empty/wrong union) are
out of reach by construction → that residual is the ONLY place an independent localizer earns its dep.
Build that fallback (FracMinHash containment lead candidate) but **gate it behind a measured depletion
trigger**; size the tail on the benchmark before paying the index cost.

## Meta-lever (4 of 6 lenses converged independently)
Convert RECTIFY's **hard/flat/quality-blind costs → CALIBRATED LIKELIHOODS, and emit POSTERIORS**, so
consensus arbitrates on a validated scale. This is simultaneously (a) the orthogonality source — the
panel shares a flat-affine / quality-blind error family — and (b) the structural defense against the
scoring-artifact family that already bit the project (the 0.09→1.07 integer-score flip). It is NOT a
new stage. Verified grounding: `local_aligner.py:73-76,604` main recurrence is flat
(_MATCH=2/_MISMATCH=-4/_GAP_OPEN=-4/_GAP_EXTEND=-1, binary homo_mask homo_mismatch=-2); `scoring.py`
consumes no phred; `HpPenaltyTable` (hp_penalty.py:143,179,261,270) exists with per-length del/ins
costs but the **gap recurrence does not consult it** — real internal headroom.

## Top 5 concepts for the design doc (each cleared orthogonality + dep-light + a position-exact ablation)
1. **Calibrated HP-length-law emission cost** wired into the live affine-gap recurrence (the 5-panel
   convergent concept). Replace the flat gap cost inside HP runs with −log P(obs_run|true_run,base)
   from the already-in-tree empirical table. Targets the dominant HP/poly-A/STR micro-indel slice.
   Dep light (KB table, numba path exists). Ablation: position-exact indel concordance at HP
   boundaries, must beat BOTH flat affine AND current homo_mismatch=-2; held-out region split.
   Novelty: uncertain (pair-HMM lineage; honest delta = applying RECTIFY's own table in the
   arbitration DP, not just the splice refiner).
2. **Channel-aware 3'-end / poly-A CPA boundary decoder** — JOINT localization+refinement (the
   inherited panel end IS the drifted quantity, can't inherit-then-refine). 2-state templated-vs-tail
   change-point under the A-run length law. RECTIFY's core mission; fires on ~every poly-A read.
   Ablation: |est−true CPA| stratified by genomic-A abutment vs heuristic walkback. Dep light
   (sequence + optional Q, no signal).
3. **Posterior-emitting refinement + likelihood-ratio consensus arbitration** — refiner emits
   posterior + runner-up; consensus compares aligner paths by calibrated LR instead of integer-max.
   Attacks the no-hill-climb-into-artifact constraint at its root. Ablation: replicate the 0.09→1.07
   scenario, show LR does not flip. Novel for the arbitration framing (SPRT/BCJR applied to consensus).
4. **POA-pooled per-locus consensus refinement (project-back)** — strongest-orthogonality reframe:
   change the unit from per-read to per-locus pool (cluster → POA consensus → align once → project
   back), denoising stochastic HP/indel error by majority BEFORE the genome decision. New orthogonal
   error family = mis-clustering/paralog collapse (no panel member has it). Dep light (abPOA/spoa KB
   C lib). Ablation: per-read HP-boundary accuracy independent vs pooled-then-projected at ≥5-read
   loci; refuted if mis-clustering errors (test on paralogs) exceed indel errors removed.
5. **FracMinHash containment fallback localizer** — the discovery-ceiling escape hatch; the only
   mechanism for reads with NO acceptable panel window. Containment is robust where exact-seed chains
   vanish at 5-10% ONT error. Orthogonal by construction (fires only where no panel placement exists).
   **Gated behind a measured depletion trigger**; ablation must also SIZE the tail before committing
   the dep.

## Credible reframes (challenge the 2-stage split)
- **Per-locus POOLED refinement** (concept 4 generalized): per-read noise is baked in before either
  stage; pool-then-place dissolves the localization/refinement debate for deep loci. Two-stage-needing
  case shrinks to the panel-dropped tail + singletons.
- **SLAM loop-closure / cross-read landmark joint estimation**: junctions/CPA as noisy observations of
  a shared latent landmark across reads; solve a sparse factor graph, snap ambiguous placements to
  consensus. Injects population/coverage structure no per-read aligner can access.
- **Particle-filter pipeline with explicit DEPLETION trigger**: panel's 5 placements seed a weighted
  hypothesis cloud; when all particles have low likelihood, SIGNAL depletion → trigger the fallback.
  Converts "when is novel localization needed" from assumption to measured trigger. (Cloud without the
  trigger = SCE-with-weights → cut.)
- **Joint EM over transcriptome + alignments**: novel junctions EMERGE as new isoform models when a
  read cluster resists all current models. Research-grade (convergence/identifiability on paralogs).

## Rejected (carry the reasons into the design doc so they aren't re-proposed)
- Pangenome/personalized-reference/variation-graph: dep violation + paradigm-rename.
- Strobemer fallback seeder / MI seed-selection: standard_paradigm_renamed (seed-chain-extend reskin).
- WFA-banded engine as a standalone member: not orthogonal (shares minimap2's affine optimum) — it is
  enabling INFRASTRUCTURE (O(Ns) whole-read realignment, build AVX2/scalar to dodge the AVX-512 SIGILL
  trap), NOT a concept.
- r-index pan-isoform localization: not orthogonal. Learned k-mer reweight as fallback: wrong slice.
  GC-stem structure-aware mismatch: untractable.

## Confabulation flags (from the synthesis self-audit — do NOT treat as findings)
- VERIFIED at file:line: the flat-recurrence / no-phred / HpPenaltyTable-exists facts above.
- NOT verified (taken on faith, not load-bearing): the ~92% uncorrected-QuantSeq-drift figure, the
  ~127 GMAP recurrent novels, the 12%-fewer-reads GMAP figure — these are exactly what the benchmark /
  Deliverable B must MEASURE; do not cite as established.
- Novelty downgrades: Trellis-MAP = pair-HMM (→uncertain); SLAM = overlaps FLAIR/StringTie collapse +
  CPA pileup (only the factor-graph-feedback framing is thin-novel); per-read self-calibration is thin.
- Inaccuracy flagged but concept survives: "HpPenaltyTable walled in the splice refiner" is wrong (it's
  imported widely) — but the gap RECURRENCE is flat and doesn't consult it, so the concept holds.

## Benchmark requirements (THIS IS the simulation-benchmark spec — Deliverable A)
**Framing metric (most important):** the discriminator is **EXACT INDEL POSITION CONCORDANCE WITH
TRUTH, not edit distance** — at every contested position edit distance is tied by construction; only
"which tied placement matches truth" separates the concepts. Build the whole contain/measure design
around position-exact truth.

CONTAIN: (2) HP runs A/C/G/T × len 1-12, del-dominant length-dependent miscall, known true run length;
(3) genomic-A-tract abutting known true CPA (the drift confound); (4) NIC/NNC novel junctions not in
GTF, BOTH canonical and non-canonical (separate discovery from annotation-snapping; measure
anti-minisplice FDR); (5) STR/tandem-repeat positions where indel placements are edit-distance-TIED;
(6) paralog loci (SMN1/SMN2-style) with known true locus (window-selection + POA mis-clustering);
(7) CONSTRUCTED panel-failure reads with known origin (all-5-misplace / GMAP-unmapped analogues — size
the discovery-ceiling tail AND the fallback recovery); (8) coverage strata singleton→deep + read-Q
deciles with BOTH calibrated and deliberately-miscalibrated phred (isolate soft-decision gain from
recalibration gain).

MEASURE: exact-indel-position accuracy at HP boundaries; corrected-run-length accuracy; |est−true CPA|
stratified by genomic-A; novel-junction RECALL and FDR stratified annotated/de-novo × canonical/
non-canonical; tied-indel placement agreement on STR; window-selection accuracy on paralogs; fallback
recovery vs random-window null; consensus accuracy integer-score vs LR (replicate 0.09→1.07, show LR
holds).

TWO META-REQUIREMENTS every ablation inherits: (1) **held-out train/test split** for any calibrated
table (else the win is overfitting); (2) **fitness = THIS truth set, NEVER the current internal
score** (the internal score was provably artifact-prone). Use the ambiguity-aware match
(`normalize_junction` / `_canonical_within_window`) so a correct call one bp into the donor/acceptor
repeat is not charged FP.

## Open questions (Q2 answered; others below)
- Q1 tail size (panel-failure residual): **benchmark-measured** (requirement 7), not a pre-build user call.
- Q2 phred plumbing: **ANSWERED** — Q is propagated to output (`chimeric_consensus.py:952`,
  `consensus.py:532`) but NOT consumed in any scoring decision; soft-decision concept must wire it in,
  ONT-Q-miscalibration caveat stands.
- Q3 calibration set / chemistry match: the HP table is currently IVT-calibrated; simulation-first
  supplies truth, held-out split mandatory. Confirm production chemistry matches the table's source.
- Q4 POA bin purity: **benchmark-measured** (requirement 6).
- Q5 **USER CALL:** is de-novo NON-canonical junction discovery a PRIMARY or secondary deliverable? The
  anti-minisplice prior carries real FDR risk and is anti-correlated with annotation-reward — its
  priority determines the first design-doc cut.
