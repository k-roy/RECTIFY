# Learned / ML Approaches for Long-Read RNA-seq Alignment & 3'-End Determination in RECTIFY

**Author role:** Algorithmic-discovery researcher proposing *learned* methods to advance
RECTIFY's alignment + 3'-end/CPA pipeline.
**Inputs read:** `02_synthesis/{COMPARISON.md, DEEP_DIVE.md}`, `01_investigation/ont_drs_ecosystem.md`,
all three `03_adversarial/` reports. Data inventory verified against the repo
(`rectify/data/`: NET-seq WT/pan/atract tables, 36-read validation set, junction_refiner.py,
corrected_consensus.py selection code).

**Framing discipline.** Every proposal below names a **concrete baseline it must beat**, the
**exact training data RECTIFY already has**, **inference cost** under the chunked-SLURM mandate,
and an **overfit/reproducibility risk** rating. Each is tagged **ESTABLISHED** (the method class
is published and works on this data type) or **NOVEL/SPECULATIVE** (no precedent on this exact
problem). I am explicit about where ML is *overkill* for a compact, well-annotated 12 Mb yeast
genome versus where it genuinely buys accuracy.

**One caveat that governs everything (from `redteam_winrates_selection.md`).** The headline win
rates (deSALT 78.9 / mapPacBio 18.2 / …) are a **single-dataset popularity vote** produced by a
legacy 5-level sort, not a measured accuracy comparison. There is **no committed ground-truth
concordance table**. *Before any ML model is trained or any heuristic is declared "beaten," a
held-out truth set and a committed `aligner_summary.tsv` must exist* (Experiments 1–3 of that
red-team). **ML cannot fix an unmeasured objective.** This is the precondition for every idea
here, and it is itself the highest-value, lowest-effort work — but it is *not* ML.

---

## 0. The Honest Map: Where ML Helps vs Where It Is Overkill

| Problem | ML verdict | Why |
|---|---|---|
| **Winner/consensus selection** (replace heuristic sort) | **ML genuinely helps** | The current sort is a hand-tuned 5-key tiebreaker with a documented herd-bias flaw; features are cheap; a GBDT over them is the textbook fix. **#1 ranked.** |
| **Poly-A / CPA boundary from raw signal** | **ML helps IF signal is retained** | The information that disambiguates the CPA (dwell time) is *destroyed by basecalling*; only a signal/dwell model recovers it. But it requires squiggle access RECTIFY does not currently keep. **High impact, high plumbing cost.** |
| **Splice-site scoring (SpliceAI/Pangolin-style)** | **Overkill for yeast, real for metazoan** | Yeast has ~300 introns, nearly all annotated, strongly canonical GT-AG. A learned splice CNN is solving a problem the GTF + a 2-line canonical prior already solves. Banked as a **metazoan-only** future module. |
| **Seq2seq alignment post-correction (transformer over the window)** | **Mostly overkill / high overfit risk** | RECTIFY's deterministic HP-aware DP modules (2C/2E/2G/2H) already encode the error model interpretably. A transformer would be a less-auditable reimplementation trained on 36 reads. **Not worth it now.** |
| **DNA/RNA foundation-model embeddings** | **Mostly hype here** | Genomic LMs help where labels are scarce and motifs are subtle. Yeast 3'-end determination is label-rich (NET-seq) and motif-simple. Honest assessment: **not worth it**, with one narrow speculative exception (chimera detection). |

---

## 1. Learned Per-Read Consensus / Winner SELECTOR  — **ESTABLISHED · RANK #1 (impact/effort)**

**Directly answers the adversarial selection finding.**

### The model
A **gradient-boosted decision tree** (LightGBM/XGBoost) — or a tiny 2-layer MLP — that, for each
`(read_id, aligner)` candidate row, outputs P(this aligner's corrected call is the true CPA /
true junction). Replace the legacy 5-level sort
(`corrected_consensus.py::merge_corrected_tsvs` Path B) with `argmax` over this score. GBDT, not
deep net: tabular features, tens-of-thousands of training rows, full feature interpretability
(SHAP), and CPU inference in microseconds/row — mandatory under the SLURM no-GPU norm.

### Features (all already computed per row)
- `hp_edit_distance` (`_cigar_hp_edit_distance`) — already the intended Path-A key.
- splice strength: canonical-tier of each N-op, annotation membership, junction overhang
  (`junction_overhang_table`), 2H refiner score.
- cross-read support: `_n_agree` (kept as a *feature*, demoted from being the decider),
  cross-read junction count at this locus.
- aligner identity (one-hot) — lets the model *learn* deSALT-junction / mapPacBio-CPA
  specialization instead of it being asserted in CLAUDE.md.
- poly-A: `pt:i` dorado tail length, `polya_score`, distance from soft-clip boundary to genomic
  homopolymer, walk-back magnitude (2E), soft-clip-rescue magnitude (2G).
- self-confidence flag `confidence` — kept as a feature, but the model learns *how much* to trust
  it rather than using it as a raw tiebreaker (fixes red-team confounder #7).

### Training data RECTIFY has
- **Labels = NET-seq.** `saccharomyces_cerevisiae_netseq_wt.tsv.gz` is a genome-wide, per-base,
  orthogonal 3'-end signal. **Label a candidate "correct" if its corrected 3' end sits at/near a
  NET-seq peak; "wrong" otherwise.** This is exactly the field's gold-standard validation signal
  (per `ont_drs §5.1`), and it makes the selector trainable at genome scale on the user's own
  real DRS reads — *not* the 36 validation reads.
- **Junction labels = annotation + multi-aligner concordance** (≥4 aligners agreeing, the same
  concordance signal already used to calibrate the overhang table in `calibrate_junction_overhang.py`).
- **Replicates** provide a reproducibility label: a call reproduced across `wt_rep1/rep2` is more
  trustworthy than a singleton.

### Baseline to beat
The legacy 5-level sort (`_five_rescued, _chimera_ok, _conf_rank, _n_agree, _span, _n_junc`).
**Metric:** fraction of reads whose selected CPA lands within ±N bp of a NET-seq peak on a
**held-out chromosome** (chromosome-level split — see overfit note). Must beat the heuristic by a
margin that survives replicate cross-validation.

### Inference cost & feasibility
Trivial. A LightGBM model is a few hundred KB, scores millions of rows/sec on one CPU core, drops
into `merge_corrected_tsvs` with zero new dependencies of consequence (LightGBM is pip-installable
and pure-CPU). **Feasibility: HIGH.** This is the single most actionable ML idea.

### Validation
- Held-out **chromosome** split (never a random row split — rows from the same locus leak).
- Compare against the heuristic on NET-seq concordance *and* on replicate reproducibility.
- The selector must be A/B-ablatable (red-team Experiment 2): ship it behind a flag so the legacy
  sort remains the falsifiable baseline.

### Overfit / reproducibility risk: **LOW–MEDIUM**
- Risk: GBDT can memorize locus-specific quirks. Mitigation: chromosome-split CV, shallow trees
  (max_depth ≤ 5), monotonic constraints (e.g., higher hp_edit_distance must never increase the
  score) for auditability.
- Reproducibility: pin seed + library version; commit the trained model **and** the training
  script + label-generation code. A model that can't be regenerated from committed code is
  banned.

### Expected impact: **HIGH.** This is the one place the adversarial review *proves* the current
mechanism is weak (herd-bias popularity vote, self-reported confidence). A learned selector trained
on an orthogonal truth signal (NET-seq) replaces a popularity vote with a calibrated accuracy
estimate — exactly the gap the red-team identified.

---

## 2. Signal-Level / Basecaller-Integrated CPA & Poly-A Model  — **ESTABLISHED (components) · NOVEL (fusion) · RANK #2**

### The problem this uniquely solves
Per `ont_drs §1.3`: the basecaller **destroys** the information that fixes the CPA. A poly-A tract
is a flat low-variance current segment whose *dwell time* encodes length; segmentation cannot
recover the base count from the sequence alone. Every alignment-only correction (2E/2C) is fighting
with one hand tied — the truth lives in the squiggle. This is the **one place** where a
signal-level model is not overkill but *necessary* for bp-accuracy.

### The model(s)
- **Established component A — dorado `--estimate-poly-a`** (`pt:i` tag). Already a learned/HMM
  estimator. RECTIFY already ingests `pt:i` into parquet. **Low-hanging fruit: actually *use*
  `pt:i` as a prior in 2E walk-back**, not just for trim/restore bookkeeping. The walk-back target
  length should be Bayesian-shrunk toward the dorado estimate.
- **Established component B — Remora-style dwell/modified-base model.** Remora (ONT's signal model
  framework) trains small NNs on raw-signal chunks anchored to a reference position. A Remora-style
  head can predict, per candidate CPA column, P(this reference base is the true 3' terminus | local
  dwell features). This is the **NOVEL fusion**: feed signal-derived dwell features *as extra
  columns* into the Module-2 correction (or into the §1 selector) so the CPA is chosen jointly from
  alignment geometry **and** raw dwell.

### Training data
- **NET-seq peaks** as CPA labels (as in §1).
- **dorado `pt:i`** as a weak/auxiliary label and feature.
- Signal: requires **POD5/squiggle retention**, which RECTIFY does **not** currently keep (it
  starts from basecalled BAM/FASTQ). This is the gating cost.

### Baseline to beat
Modules 2E + 2C (alignment-only HP walk-back). Metric: ±bp NET-seq concordance at the CPA, on
held-out reads where 2E/2C are known to be ambiguous (long A-tracts terminating in a genomic A-run
— exactly the `find_polya_boundary` false-stop cases in CLAUDE.md v3.0.3).

### Inference cost & feasibility
- **`pt:i`-as-prior: HIGH feasibility, near-zero cost** — do this first, it's not even a new model.
- **Full Remora dwell fusion: MEDIUM-LOW feasibility.** Requires (a) POD5 retention through the
  pipeline (storage + plumbing on `$SCRATCH`), (b) signal-to-reference anchoring (move-table from
  dorado), (c) per-read signal I/O that does not blow the chunked-SLURM I/O budget. GPU helpful for
  training, not strictly needed for the tiny inference head.

### Overfit / reproducibility risk: **MEDIUM**
Remora models are chemistry-specific (R10.4.1 / RNA004). A model trained on one pore generation
must be re-trained for the next (same caveat as the empirical penalty tables, which are already
"R10.4.1-specific"). Reproducibility is fine *if* the squiggle inputs are archived; otherwise the
training set is irreproducible.

### Expected impact: **HIGH ceiling, gated by plumbing.** The `pt:i`-prior sub-step is high-value
and immediate. The full dwell-fusion is the only path to sub-bp CPA in long A-tracts but costs a
real signal-retention engineering investment. **Recommend: ship the `pt:i` prior now; scope the
dwell fusion as a funded R&D track, not a quick win.**

---

## 3. Learned Splice-Site Model (SpliceAI / Pangolin-style)  — **ESTABLISHED · OVERKILL for yeast, RANK #1 for metazoan**

### The model
A 1-D dilated-convolution sequence model (SpliceAI: ~10 kb receptive field; Pangolin: tissue/usage
extension) that outputs per-base P(donor)/P(acceptor) from genomic sequence alone. Used as a
**junction scorer/tie-breaker inside Module 2H** (`junction_refiner.py`), replacing/augmenting the
hard-coded canonical GT-AG `_CANONICAL_HP_PRIOR = 0.5`.

### Honest feasibility & value: **yeast = OVERKILL; metazoan = the headline future module**
- **Yeast.** ~300 introns, almost all annotated, overwhelmingly canonical GT-AG, short, few
  alternative splice sites. Module 2H already uses (a) the annotation and (b) a canonical prior as
  tie-breakers. A learned splice CNN would, on yeast, mostly reproduce "GT-AG + annotated" — which
  two lines of code already encode. **Training data is also thin**: ~300 positive junctions is far
  below SpliceAI's training regime (GENCODE, ~10⁵ junctions). **Verdict: do not build a yeast
  splice CNN.** It is solving a solved problem and would overfit ~300 sites.
- **Metazoan (the stated future target).** Here it inverts. Human/mouse have abundant non-canonical
  junctions, micro-exons, and dense alternative splicing where the annotation is incomplete and the
  GT-AG prior is too weak. **SpliceAI/Pangolin are pre-trained and downloadable** — no training
  needed; use as a frozen scorer. This is the natural answer to the COMPARISON.md observation that
  GMAP's "splice-model-first" sandwich DP is the closest analogue to Module 2H: a learned splice
  model *is* the modern splice-model-first scorer, and it slots into 2H's existing
  "score-first, canonical-as-tiebreaker" interface cleanly.

### Training data
- **Yeast: none needed and none sufficient** (use the GTF).
- **Metazoan: none needed** — SpliceAI/Pangolin ship pre-trained weights. RECTIFY's own multi-aligner
  concordance could *fine-tune* a usage head later, but that is speculative.

### Baseline to beat
Module 2H's current `_CANONICAL_HP_PRIOR = 0.5` + annotation tie-break. On yeast it will not beat
it. On a held-out metazoan dataset with novel/non-canonical junctions, the SpliceAI score should
beat the flat canonical prior on junction-placement accuracy.

### Inference cost
SpliceAI scores can be **precomputed once per genome** (like the existing `annotation.junc.bed`
cache) into a per-base track, then looked up at O(1) during 2H. **No per-read NN inference** —
fully compatible with chunked SLURM. Feasibility: **HIGH for metazoan (precompute-and-cache).**

### Overfit / reproducibility risk: **LOW** (frozen public model, genome-level precompute).

### Expected impact: **NIL for yeast (skip), HIGH for metazoan (the priority module when RECTIFY
moves to human/mouse).** Mark explicitly as a *deferred* module gated on the metazoan roadmap.

---

## 4. Learned Error-Model / Alignment Post-Correction (seq2seq / transformer over the window)  — **NOVEL/SPECULATIVE · MOSTLY NOT WORTH IT**

### The idea
A small encoder-decoder or windowed transformer that takes the alignment window (read bases +
reference + CIGAR) around the 3' terminal exon and emits a corrected CIGAR / corrected 3' position
— a learned replacement for the deterministic 2C/2E/2G/2H surgery.

### Why it is mostly overkill here — be blunt
- RECTIFY's correction modules are **interpretable, validated DP** that already encode the ONT HP
  error model via empirical penalty tables. A transformer would be a **less-auditable
  reimplementation** of a working, debuggable system.
- **Training data is fatally thin for a sequence model.** The labeled set with bp-exact corrected
  CIGARs is the **36 validation reads** (`validation_reads.bam` + `corrected_reads.tsv`). A
  transformer trained on 36 examples is a memorizer, not a model. NET-seq gives CPA *position*
  labels but not corrected-CIGAR labels, so it can't supervise CIGAR surgery directly.
- **Reproducibility/overfit risk: HIGH.** Non-determinism, chemistry-specificity, and the
  CLAUDE.md history of CIGAR-surgery edge cases (terminal D/N stripping, in-deletion-span CPAs)
  show this domain is full of rare structural cases a small learned model will silently get wrong —
  and unlike the DP modules, you can't read off *why*.

### The one defensible narrow version
A **learned HP-indel placement scorer** (not a full seq2seq): a tiny model that, given an
ambiguous homopolymer run with several near-equal DP placements, picks the placement — i.e. learn
the penalty *values* that 2C/2H currently read from the empirical table. This is really "learn the
penalty table from data" and is **strictly an extension of the existing empirical-table approach**,
trainable on multi-aligner concordance + NET-seq, low-risk, interpretable. **That** is worth a
small experiment; the full transformer is not.

### Baseline / impact
Baseline: the empirical penalty tables + DP. Expected impact of the full transformer: **low and
high-risk → NOT recommended.** Expected impact of the narrow "learn the penalty table" version:
**modest, low-risk → optional.**

---

## 5. Foundation-Model / DNA-RNA LM Embeddings (Nucleotide Transformer, DNABERT-2, Evo, etc.)  — **SPECULATIVE · MOSTLY HYPE HERE**

### Honest assessment
Genomic LMs earn their keep where **labels are scarce and the signal is subtle/distributed**
(regulatory-grammar prediction, variant-effect, cross-species transfer). RECTIFY's core problem is
the opposite: **labels are abundant** (genome-wide NET-seq), the **genome is tiny and fully known**
(no need to "learn" yeast sequence statistics a 12 Mb FASTA already contains), and the decisive
features (homopolymer length, dwell time, canonical dinucleotide) are **explicit and
low-dimensional**. Embedding a 6-mer window into a 1024-d transformer space to then predict "is
this a CPA" is using a sledgehammer where a GBDT over 15 hand-features (§1) wins on accuracy,
speed, interpretability, *and* reproducibility.

- **Inference cost** is also disqualifying for the chunked-SLURM mandate: per-read LM embedding is
  GPU-hungry and orders of magnitude slower than the §1 GBDT — antithetical to the
  "$SCRATCH 75 GB/s, no NFS contention, finish in minutes" design philosophy.
- **Overfit/repro risk: HIGH** (giant frozen models, version drift, opaque failure modes).

### The single narrow place it *might* help (speculative)
**Chimera artifact detection.** `ont_drs §4.4` cites a 2024 **genomic-language-model** study that
flags DRS chimeras (two molecules basecalled as one) — a real 3'/junction confounder that
RECTIFY's `build_chimeric_cigar` + overhang filter handle heuristically. A pretrained genomic LM
that scores "does this read look like one coherent transcript or a stitch of two?" is a *plausible*
auxiliary feature for the chimera flag — because here the signal (a discontinuity in sequence
"grammar") is exactly the distributed/subtle kind LMs are good at, and the heuristic overhang
filter is admittedly crude. **Even so: a feature, not a pipeline, and only after §1 is in place.**

### Verdict: **NOT worth it** as a core method. One speculative auxiliary use (chimera scoring),
ranked last.

---

## 6. Ranking (Impact / Effort) and What NOT To Build

| Rank | Proposal | Tag | Impact | Effort | Verdict |
|---|---|---|---|---|---|
| **1** | **Learned winner SELECTOR (GBDT) trained on NET-seq labels** (§1) | ESTABLISHED | **High** | **Low** | **BUILD FIRST.** Directly fixes the proven selection-mechanism flaw; cheap; CPU-only; interpretable. |
| **2a** | **Use dorado `pt:i` as a Bayesian prior in 2E walk-back** (§2) | ESTABLISHED | Med-High | **Very low** | **BUILD NOW.** No new model; data already ingested. |
| **2b** | **Signal/dwell (Remora-style) CPA fusion** (§2) | NOVEL fusion | **High ceiling** | **High (POD5 plumbing)** | **Scope as funded R&D**, not a quick win. The only path to sub-bp CPA in long A-tracts. |
| **3** | **SpliceAI/Pangolin junction scorer in Module 2H** (§3) | ESTABLISHED | **Metazoan: High** / Yeast: ~Nil | Low (precompute) | **DEFER to metazoan.** Skip for yeast. |
| **4** | **"Learn the penalty table" HP-placement scorer** (§4 narrow) | NOVEL, low-risk | Modest | Low-Med | Optional small experiment. |
| — | **Seq2seq/transformer alignment post-correction** (§4 full) | SPECULATIVE | Low | High | **DO NOT BUILD.** 36-read training set; high overfit; replaces auditable DP with a black box. |
| — | **DNA/RNA foundation-model embeddings as core method** (§5) | SPECULATIVE | Low | High | **DO NOT BUILD.** Label-rich, motif-simple problem; GPU cost violates SLURM design; one speculative chimera-feature aside. |

### The precondition that outranks all ML (from `redteam_winrates_selection.md`)
**Build the ground-truth concordance harness first** (NET-seq / replicate-based per-aligner
accuracy, committed `aligner_summary.tsv`, Path A vs Path B ablation). It is the objective every
model in §1–§4 trains and validates against. **No ML proposal here is creditable until a held-out
truth metric exists** — and that harness is itself low-effort, high-value, and not ML. Said plainly:
the most valuable next step is *measurement*, and the highest-impact *learned* method (§1 selector)
is the natural first consumer of that measurement.

---

## 7. Cross-Cutting Risks (apply to every learned proposal)

1. **Single-dataset, single-chemistry training.** All RECTIFY truth signals come from one
   `wt_by4742_rep1` R10.4.1 yeast run. Any model risks learning that dataset, not the biology.
   *Mitigation:* chromosome-level held-out splits, replicate cross-validation, and a hard rule that
   models retrain per chemistry (RNA002 → RNA004 → DNA pores differ — same caveat as the penalty
   tables).
2. **Label leakage from the heuristic.** If labels are derived from multi-aligner concordance and
   the model also sees aligner identity, it can learn the herd, not the truth (red-team confounder
   #1). *Mitigation:* prefer **NET-seq** (aligner-independent) labels; treat concordance labels as
   weak supervision only.
3. **Reproducibility is a release gate.** Commit training script + label-generation code + pinned
   model artifact. A model regenerable only from un-archived squiggle/intermediate data is banned —
   this matches RECTIFY's existing discipline (committed penalty tables, committed validation BAMs).
4. **Interpretability ≥ marginal accuracy.** Given the CLAUDE.md history of subtle CIGAR/strand
   correctness bugs, prefer GBDT + SHAP / monotonic constraints over deep nets wherever the accuracy
   gap is small. An auditable model that is 1% worse beats an opaque one that is 1% better in a
   pipeline whose entire value proposition is bp-level correctness.
```
