# Improving Splice-Junction Modeling & Refinement in RECTIFY

**Author role:** Algorithmic Discovery — splice-junction proposals for RECTIFY's correct-first
long-read ensemble (ONT DRS/cDNA, *S. cerevisiae*).
**Date:** 2026-06-19.
**Grounded in:** `02_synthesis/{COMPARISON,DEEP_DIVE}.md`; dossiers deSALT / uLTRA / gmap / minimap2;
`03_adversarial/{redteam_denovo,redteam_winrates_selection,redteam_annotation_ecosystem}.md`;
RECTIFY source `rectify/core/splice/{junction_refiner,junction_scoring,hp_penalty,
calibrate_junction_overhang}.py` and `corrected_consensus.py`; CLAUDE.md v3.1.x / v3.3.0 notes.

> Verified against `origin/drs-validation-rebuild` @ 366c885 (2026-06-19).

**Convention.** Each proposal is tagged **ESTABLISHED** (the mechanism exists in a published
aligner / standard method, we are porting it) or **NOVEL** (not done this way in the surveyed
tools). Every proposal carries: *mechanism · failure mode addressed · feasibility · expected
impact · validation · risk*.

---

## 0. Two framing facts that reshape every proposal

Before the proposals, two findings from the adversarial pass change what is worth building:

**(F1) Winner-selection runs HP-edit-distance (Path A) in production; the popularity vote is the fallback.**
Both production call sites of `merge_corrected_tsvs` — `single_sample.py:238-250` and
`split_command.py:1085-1094` — pass per-aligner **raw** BAMs + genome, and `merge_corrected_tsvs`
activates the lazy raw-BAM HP path without materializing corrected BAMs:
`use_hp_ed = bool(per_aligner_corrected_bams or per_aligner_raw_bams)` (`corrected_consensus.py:1262`;
lazy path `:1264-1364`). So **Path A (HP-edit-distance on corrected 3'-ends) is the default production
metric.** The **legacy 5-level sort** (`_five_rescued, _chimera_ok, _conf_rank, _n_agree, _span,
_n_junc`) is the **fallback only** — used when no BAMs/genome are supplied, or if BAM staging yields no
rows (`use_hp_ed=False`, `corrected_consensus.py:1377`). In that fallback the dominant non-flag
discriminator `_n_agree` is a **popularity vote**: the aligner sitting in the majority cluster wins
ties, so deSALT's "homogeneous junctions" can win by **herd bias** (it makes other aligners agree with
it) rather than demonstrated accuracy, and uLTRA's GTF-snapping is **annotation-circular** ("snapped-to-GTF
= correct" is assumed, not tested). These herd-bias and annotation-circularity critiques apply to the
fallback path and to any selection that leans on cross-aligner agreement — they do not describe the
default production metric.

**Consequence:** A better junction *scorer* (Proposals 1–2) is **already wired into the default
selection metric on the build** (Path A scores HP-aware corrected positions); the remaining
high-leverage work is to make that score *calibrated and principled* (Proposal 2) and to charge
N-ops their true cost (Proposal 4 / J1a), rather than to "turn a dead path on." Several proposals
below therefore sharpen — rather than replace — the live HP-ED selection.

**(F2) "8-nt exon" is unsupported; the real micro-exon floor is ~6 nt** (`redteam_annotation`
U9: minimap2 issue #253 is a 6-nt exon; the uLTRA paper bins ≤10 nt at ~60% accuracy). Micro-exon
proposals (Proposal 3) must be validated against **≤10 nt** reality, not a headline number, and on
yeast — where micro-exons are *rare*, capping the achievable gain.

---

## Proposal 1 — Cross-read junction consensus as a first-class step

**(ESTABLISHED mechanism — deSALT Pass-2 exon pooling / GMAP-absent; NOVEL as a decoupled
post-alignment RECTIFY stage operating across all 5 aligners.)**

### What RECTIFY does today
Module 2H (`junction_refiner.refine_read_junctions`) is **per-read**: for each N-op it scores
candidate junctions from a pool and re-snaps that read. The pool (`build_junction_pool`) already
unions annotated + per-aligner-BAM junctions — so RECTIFY *has* the raw material for cross-read
consensus but never computes a **consensus site** or a **cross-read support count**. Each read is
refined in isolation against a static candidate list with no notion of "how many reads/aligners
independently support junction *J*."

### Mechanism
Insert a **junction-consensus build** step between alignment and per-read refinement, generalizing
deSALT's cross-read exon pooling to the *whole ensemble* (deSALT only pools its own reads; no panel
aligner pools across aligners):

1. **Pool pass (one stream over all 5 aligner BAMs + annotation).** For every N-op in every read of
   every aligner, emit `(chrom, donor, acceptor, strand, aligner_id)`. Aggregate into a
   `JunctionConsensus` table keyed by `(chrom, strand, donor, acceptor)` with fields:
   `n_reads`, `n_aligners` (distinct aligners ≥1 read), `canonical_tier`, `is_annotated`,
   `mean_flank_overhang`, `splice_strength` (Proposal 2).
2. **Consensus collapse.** Cluster nearby junctions sharing a strand within a small window
   (`±max_boundary_shift = 50 bp`, already a 2H constant). Within a cluster, elect a **consensus
   site** = the member maximizing a calibrated support score (Proposal 2's probability), *not* raw
   read count (read count alone is the herd-bias trap from F1). Store the cluster → consensus map.
3. **Re-snap pass.** Per-read refinement (existing 2H surgery) now snaps each N-op to its cluster's
   **consensus** site, with the per-read sequence score (`_score_junction`) as the *veto*: a read is
   only moved if its own sequence supports the consensus at least as well as its current placement
   (prevents forcing a genuinely-different isoform onto the majority).

### Data structures & complexity
- `JunctionConsensus`: `dict[(chrom,strand) -> IntervalTree(donor) ]`; each node carries the support
  record. Build is **O(R · j)** (R reads, j N-ops/read) — one linear pass, same cost as the
  existing pool build.
- Consensus collapse: per (chrom,strand), sort cluster members O(k log k), k = junctions in window;
  total **O(J log J)**, J = distinct junctions (~thousands in yeast). Negligible.
- Re-snap: unchanged from current 2H — **O(R · candidates)** but `candidates` shrinks to the
  cluster consensus + the read's current site (≈2), so it is *faster* than today's radius scan.
- Memory: O(J) ≈ a few MB for yeast; trivially fits the chunk-array workers (fork-shared, like the
  existing `_WORKER_POOL_STATE`).

### Failure mode addressed
Per-read 2H can snap two reads of the *same* isoform to two different nearby sites (the inhomogeneity
deSALT fixes internally and minimap2/gapmm2/mapPacBio do not). Cross-read consensus enforces the
homogeneity directly — and does it **across aligners**, so a junction one aligner placed precisely
rescues the same junction in reads where every aligner was noisy. It also **breaks F1's herd bias**:
the consensus site is chosen by a calibrated score, not by which aligner family is largest.

### Feasibility
Medium. All primitives exist (`build_junction_pool`, `_candidates_near`, IntervalTree usage in
`calibrate_junction_overhang.py`, the fork-shared worker pool). New code ≈ one module
(`junction_consensus.py`) + a wiring change in the correction driver. The prescan path already builds
a `prebuilt_junction_pool` — extend it to a `prebuilt_junction_consensus`.

### Expected impact
**High.** This is the single mechanism the win-rate analysis attributes to the best aligner (deSALT),
generalized to all five and de-biased. Directly improves junction homogeneity → tighter corrected-CPA
clusters at clean exon edges → the actual quantity NET-seq validates.

### Validation
- **cat9** (Module 2H junction-refine reads, 4 reads in `validation_reads.bam`): consensus must not
  regress the 4 known corrections (00a1c9b3, 00a1e01e, 0b3b593b, d3357db5); ideally tightens them.
- **cat7** (chimeric/`XU=1`): consensus + Proposal 2 overhang/strength must keep cat7 junctions.
- **Homogeneity metric** (new): per annotated intron, stdev of donor/acceptor across reads should
  drop vs per-read 2H. Measure on `wt_by4742_rep1`.
- **NET-seq concordance** of corrected 3' ends (orthogonal gold standard) before/after.

### Risk
Medium. Consensus can **wrongly collapse two genuinely distinct nearby junctions** (alt-5'SS/3'SS
isoforms < 50 bp apart — real in yeast, e.g. some RP genes). Mitigations: (a) the per-read sequence
veto in step 3; (b) cluster only junctions sharing one boundary exactly (alt-donor and alt-acceptor
isoforms differ at one end only, so anchoring on the shared end avoids merging them); (c) keep
`max_boundary_shift` tight. The redteam's herd-bias warning means we must **measure against truth, not
just internal agreement** — internal homogeneity is exactly what can be gamed.

---

## Proposal 2 — Unified calibrated junction scorer (probabilistic, replaces heuristic tiers)

**(ESTABLISHED components — MaxEnt/GeneSplicer splice models, empirical HP penalties, GMAP's
canonical-as-reward; NOVEL as a single calibrated log-odds combining sequence support + splice
strength + cross-read count + annotation prior, replacing 2H's hand-ordered tuple.)**

### What RECTIFY does today
2H selection is a **lexicographic tuple** (`redteam` confirms it): `(score_cmp, tier/is_alt,
is_alt/tier, is_novel, abs_delta)` with a hard `_CANONICAL_HP_PRIOR = 0.5` discount and an
order that flips on `current_tier >= 4`. This is brittle (the v3.1.1 `range(L)` and `int(score)`
binning bugs were direct symptoms) and **not a probability** — it cannot be compared across reads or
fed into selection as a calibrated quantity. The HP penalty table (`hp_penalty.HpPenaltyTable`,
`penalty_scores.tsv`) is already empirical and calibrated; the *combination* with splice strength and
support is heuristic.

### Mechanism
Define, for a candidate junction *J* supporting read *r*, a **log-odds score**

```
  S(r, J) =  w_seq · logP_seq(r | J)          # HP-aware sequence support (existing _score_hp_anchored, as a likelihood)
           + w_don · logP_donor(J)            # MaxEnt 5'SS strength (MaxEntScan MES:5' or GeneSplicer)
           + w_acc · logP_acceptor(J)         # MaxEnt 3'SS strength (MES:3')
           + w_xr  · log(1 + n_aligners(J))   # cross-read/aligner support (Proposal 1), saturating
           + w_ann · 1[J annotated]           # annotation prior — TIE-BREAKER weight, small, never gating
```

- `logP_seq` reuses the existing HP-anchored DP but emitted as a **per-base log-likelihood** under
  the empirical del/ins/sub costs (the costs are already −log-probabilities in spirit; calibrate the
  constant so they are true log-odds).
- `logP_donor/acceptor` = **MaxEntScan** (Yeo & Burge 2004) 9-mer donor / 23-mer acceptor maximum
  entropy scores, or **GeneSplicer**. These are the field-standard sequence splice-strength models
  GMAP's modern builds adopted (MaxEnt, per gmap.md §3). **Train yeast-specific PWM/MaxEnt tables
  from the SGD annotation** (compact, well-annotated — yeast is *ideal* for this; see Proposal 5).
- Weights `w_*` fit by **logistic regression** on a labeled set: positives = annotated junctions
  with ≥4-aligner concordance (the same concordant-read definition `calibrate_junction_overhang.py`
  already uses to calibrate overhang thresholds); negatives = junctions that lost to a concordant
  alternative. This is a **direct generalization of the existing overhang calibration** — same data,
  richer model.

The canonical/annotation signals enter as **additive priors with small fitted weights**, never as
gates — preserving the v3.1.x "no candidate guards, canonical/annotation as tie-breaker only" policy
(which is marked PERMANENT in CLAUDE.md) while making the trade-off *calibrated* instead of a magic
0.5.

### Failure mode addressed
1. The hand-tuned tuple's brittleness (F1: the metric is not even a comparable scalar; the binning
   bugs). 2. The herd-bias of raw read counts (saturating `log(1+n_aligners)` weight, fit so it can
   never override clear sequence evidence). 3. minimap2's `--splice-flank=no` deliberately discards
   the flanking-base splice prior and hands it to RECTIFY — but RECTIFY currently has **no positive
   splice-strength model**, only canonical/non-canonical tiers. MaxEnt fills exactly that gap.

### Feasibility
Medium. MaxEntScan is a small, public, well-understood model (5'/23-mer score matrices); a yeast PWM
is hundreds of lines to build from the GFF. The logistic fit is one scikit-learn call on the existing
concordance data. The scorer plugs into `_score_junction`'s return as an added term. Main cost is the
**calibration harness + held-out evaluation**, not the model.

### Expected impact
**High.** HP-ED selection (Path A) is the live production metric (see F1), so a *calibrated* score
`S(r,J)` upgrades the metric that already drives selection — no separate wire-in step is needed. As a
pure 2H refinement upgrade it is
**medium** (2H already gets most cat9 reads right). The big win is replacing the lexicographic HP-ED
tuple / `_n_agree` fallback with a calibrated `S(r,J)`, which is the redteam's #1 recommended fix.

### Validation
- All 146 `test_junction_refiner.py` + `test_consensus_selection.py` tests must pass with default and
  empirical penalty table (the current bar).
- cat9 (4 reads) refined identically or better; cat7 chimeric junctions retained.
- **Calibration curves**: predicted P(correct) vs observed concordance on held-out junctions
  (reliability diagram) — the scorer must be *calibrated*, not just discriminative.
- **Ablation** (redteam experiment 5): re-rank with only `S(r,J)` and random tie-break; quantify how
  many wins flip vs the popularity tuple. A large flip count proves the old metric was a tiebreaker
  artifact.

### Risk
Medium. MaxEnt is trained on *metazoan* splice sites in its standard form — **must retrain for yeast**
(yeast 5'SS is more constrained: strong GUAUGU branch context, near-invariant). Risk of over-fitting
the logistic weights to one dataset (`wt_by4742_rep1`); mitigate with cross-validation and a second
DRS sample (redteam experiment 7). Do **not** let `w_ann` grow large — annotation-circularity
(redteam U8/claim 8) is the failure to guard against on yeast where most reads are FSM.

---

## Proposal 3 — De-novo micro-exon recovery (≤10 nt) without annotation

**(ESTABLISHED — GMAP's splice-probability-gated microexon search + uLTRA's all-exon injection;
NOVEL for RECTIFY's de-novo aligners, which currently inherit whatever the aligner emitted.)**

### What RECTIFY does today
Nothing dedicated. uLTRA recovers micro-exons by injecting *annotated* exons into chaining
(annotation-circular — useless for novel micro-exons); deSALT's k=8 local-hash rescues *some*; the
others drop exons < ~13 nt. RECTIFY has no post-hoc micro-exon insertion.

### Mechanism (GMAP-style, annotation-free)
When two consensus exons flank an N-op and the read has **unexplained internal mismatch/indel mass
near the junction** (an `_has_boundary_error` hit, already computed), probe for a hidden micro-exon:
slide a candidate ≤10 nt exon block inside the intron and accept it **only if** both flanking
splice sites clear a MaxEnt **splice-probability gate** (Proposal 2's `logP_donor/acceptor` >
threshold) *and* inserting it lowers the read's sequence score by > the cost of the two new junctions.
This is exactly GMAP's `--microexon-spliceprob` gate (0.90–0.95), reproduced from RECTIFY's own
splice model. Cross-read support (Proposal 1) gates it further: only insert a micro-exon seen
concordantly in ≥N reads — preventing per-read noise from manufacturing spurious 6-mers.

### Failure mode addressed
Reads truncated/mis-joined across a true micro-exon (the uLTRA niche) — recovered **de novo**, so it
works on *novel* micro-exons that uLTRA's annotation-injection cannot, avoiding the circularity
(redteam claim 8).

### Feasibility
Medium-low. Needs Proposal 2's splice model first. The search is bounded (≤10 nt × small intron
window) so cheap, but CIGAR surgery to *insert* an exon (M-N-M-N-M) is more invasive than 2H's
boundary-shift surgery.

### Expected impact
**Low on yeast** (per F2 + uLTRA dossier: yeast has few introns, fewer micro-exons; uLTRA wins only
2% partly for this reason). Higher value as the metazoan-porting story. On yeast, expect a handful of
reads. Rank accordingly.

### Validation
No current validation read exercises micro-exon insertion — would need a **new cat10** built from a
yeast gene with a documented short internal exon (rare; may need a synthetic/SIRV control). NET-seq
does not help here (it is a 3'-end assay). This is the weakest-validated proposal on yeast.

### Risk
Medium-high: spurious micro-exon insertion creates false isoforms; the splice-prob + cross-read
double gate is essential. Low blast radius if gated conservatively (few reads affected).

---

## Proposal 4 — Better novel-junction handling (avoid minimap2 noise AND uLTRA circularity)

**(NOVEL synthesis — calibrated novelty acceptance via Proposal 2 + Proposal 1, replacing both the
"free N-op" loophole and annotation-snapping.)**

### What RECTIFY does today
Two opposite failure modes coexist. minimap2 emits noisy novel junctions (±1–few bp jitter; and the
redteam notes **N-ops cost 0 in HP-edit-distance**, so an aligner can invent a free false intron to
*lower* its score). uLTRA hard-snaps to the GTF (circular: rewards canonicality over truth). 2H's
v3.1.x policy correctly refuses to *gate* on canonical/annotation, but production selection (F1) has
**only the binary overhang flag** (`_chimera_ok`) as a defense against false introns.

### Mechanism
1. **Close the free-N loophole.** In the HP-edit-distance / unified scorer, charge an N-op a **small
   calibrated open cost** = −logP(novel intron of this length at this site) from Proposal 2, instead
   of 0. A *supported* canonical novel junction pays ~nothing (high splice prob); an *invented* one
   to skip mismatches pays its true improbability. This directly fixes the redteam's "N=0 free-pass"
   (winrates note (a)).
2. **Graded novelty acceptance**, not binary snapping: a novel junction is accepted on its **own**
   MaxEnt strength + cross-read support (Proposals 1–2), with annotation only a small prior. A
   strongly-supported novel GT-AG site beats a weakly-supported annotated one (preserves the
   PERMANENT no-gate policy); a 1-bp-jittered novel site loses to the consensus cluster (Proposal 1
   collapses the jitter).

### Failure mode addressed
minimap2 jitter (collapsed by consensus), the free-false-intron exploit (charged its real cost), and
uLTRA circularity (novelty judged by sequence+support, not GTF membership).

### Feasibility
Medium — it is the natural consequence of Proposals 1+2, not separate machinery. The N-op cost change
is a few lines in `_cigar_hp_edit_distance` / the scorer.

### Expected impact
**Medium-high** for selection correctness (it is the redteam's most-emphasized correctness gap), even
if it moves few *validation* reads. Mostly protects against regressions and bad-aligner gaming.

### Validation
- **False-intron stress test** (redteam experiment 8): inject reads with a known-false long intron;
  confirm the N-cost + overhang filter prevent deSALT/gapmm2 from winning them.
- **Annotation-circularity test** (redteam experiment 6): restrict to reads whose true junction is
  novel (orthogonal truth); uLTRA's snapping advantage must vanish.
- cat7/cat9 unchanged.

### Risk
Low-medium. Over-charging N-ops could suppress *real* novel introns; the cost must be calibrated
(Proposal 2's logP), not a fixed penalty — a fixed penalty would reintroduce a gate.

---

## Proposal 5 — Learned sequence splice-site model (SpliceAI-style CNN): worth it for yeast?

**(ESTABLISHED — SpliceAI / Pangolin CNNs; assessed, mostly NOT recommended for yeast.)**

### Assessment
**Not worth a deep-learning CNN for yeast; worth a lightweight PWM/MaxEnt (Proposal 2) instead.**

- **Yeast is the wrong regime for SpliceAI.** SpliceAI's value is *long-range* context (10 kb
  flanks) capturing exon-definition signals in large metazoan introns. Yeast introns are < 1 kb,
  mostly single-intron genes, with a **highly constrained 5'SS (GUAUGU) and branchpoint (UACUAAC)**
  — almost fully captured by a short PWM/MaxEnt window. The compact, well-annotated genome means a
  **position weight matrix / MaxEnt model trained on SGD junctions is near-saturating**; a CNN would
  add cost and over-fitting risk for marginal gain.
- **Training data is thin.** ~280 yeast introns is plenty for a PWM/MaxEnt but **small for a CNN**;
  transfer-learning a metazoan SpliceAI onto yeast inherits metazoan biases (the same retrain problem
  as Proposal 2, worse).
- **Determinism / reproducibility** matter for a validation panel; a PWM is deterministic and
  inspectable, a CNN is neither without pinning.

**Where a learned model *would* pay off: the metazoan port.** If RECTIFY targets human/mouse DRS,
SpliceAI/Pangolin scores as the `logP_donor/acceptor` term in Proposal 2 (drop-in: they output
per-base donor/acceptor probabilities) would beat MaxEnt — exon-definition over long introns is
exactly their strength. Architect Proposal 2's splice term behind an **interface** (`splice_strength(chrom,
pos, strand) -> (logP_donor, logP_acceptor)`) so the yeast PWM and a metazoan CNN are
interchangeable backends. That keeps the door open at zero present cost.

### Feasibility / impact / risk
Feasibility low for yeast (high effort, low payoff). Expected impact **low (yeast) / high (metazoan)**.
Risk: over-engineering; metazoan-bias transfer. **Recommendation: build the PWM/MaxEnt backend now
(Proposal 2), define the interface, defer the CNN to the metazoan milestone.**

### Validation
If ever built: standard SpliceAI top-k accuracy on held-out SGD junctions vs the PWM/MaxEnt baseline —
the CNN must *beat* the PWM on yeast to justify itself (it likely won't).

---

## Ranking by impact / effort

| Rank | Proposal | Impact | Effort | Est. tag | One-line rationale |
|---|---|---|---|---|---|
| **1** | **P2 — Unified calibrated scorer (sharpen the live HP-ED metric)** | High | Med | NOVEL combo | HP-ED already runs in production on the build; P2 makes it *calibrated* (quality, not lexicographic tuple) and adds the splice-strength model RECTIFY lacks; everything else builds on it. |
| **2** | **P1 — Cross-read junction consensus (first-class)** | High | Med | ESTABLISHED/NOVEL | Generalizes deSALT's winning mechanism to all 5 aligners, de-biased; tightens homogeneity → better CPA. |
| **3** | **P4 — Calibrated novel-junction handling** | Med-High | Low | NOVEL synthesis | Mostly free given P1+P2; closes the N=0 free-intron loophole + kills annotation-circularity. |
| **4** | **P3 — De-novo micro-exon recovery** | Low (yeast) | Med | ESTABLISHED | GMAP/uLTRA capability, annotation-free; but yeast has few micro-exons — defer, or save for metazoan. |
| **5** | **P5 — Learned CNN splice model** | Low (yeast) | High | ESTABLISHED | Not worth it for compact well-annotated yeast; build PWM/MaxEnt now, keep a CNN-swappable interface for the metazoan port. |

**Cross-cutting prerequisite (do first, ~free):** the **provenance + Path A-vs-B-fallback ablation**
from `redteam_winrates_selection.md` (experiments 1–2). Path A (HP-ED) is the production default (both
call sites pass raw BAMs + genome — `corrected_consensus.py:1262`, `single_sample.py:238-250`,
`split_command.py:1085-1094`) and Path B is the fallback, so this is an ablation of **two live paths**.
Confirm the cited 78.9/18.2/2/0.8/0.1 win rates are post-v3.3.0-fix and **commit `aligner_summary.tsv`**
(the win-rate provenance caution stands: single dataset, un-committed artifact). Without this we cannot
tell whether any proposal *improves* selection or merely changes the metric. This single audit gates the
value claims of P1–P4.

---

## Honest caveats carried from the adversarial pass

- **Win rates are one dataset, possibly pre-fix.** Under the fallback popularity sort, deSALT's
  "homogeneity → wins" could be herd bias (`_n_agree` popularity), not accuracy; P1 must be validated
  against **orthogonal truth (NET-seq, novel-junction sets)**, not internal agreement — internal
  homogeneity is the exact quantity that can be gamed.
- **uLTRA's micro-exon / annotation strength is circular** on yeast; P3/P4 deliberately route around
  it with de-novo splice-strength gating.
- **N-ops are currently free** in HP-edit-distance (cost 0, `corrected_consensus.py:142-143`) — any
  scorer change must charge introns their calibrated improbability or it inherits the same exploit. On
  yeast the overhang filter is the main defense; the junction-anchor gate defaults OFF at 0 bp for
  yeast (10 bp for human).
- **MaxEnt/SpliceAI are metazoan-trained** — yeast needs a retrained PWM/MaxEnt; the CNN is deferred.
- **All penalty/splice tables are R10.4.1 + S. cerevisiae-specific** and must not transfer to HiFi or
  other organisms without recalibration. These tables are **bundled** under
  `rectify/data/genomes/*/penalty_tables/` for S. cerevisiae and H. sapiens with protocol variants,
  and `--Scer` auto-resolves them — so the open work is "validate/version the bundled tables" (J2),
  not "regenerate absent tables."
