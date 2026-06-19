# Improving Seeding & Indexing for RECTIFY — Algorithmic Discovery Proposals

**Goal.** Improve the state of the art in *seeding & indexing* for long-read RNA-seq alignment
in RECTIFY's context: ONT DRS / dT-primed cDNA (R10.4.1) + PacBio HiFi, where the deliverable
is a **base-accurate 3' end (CPA)** and **sequence-supported splice junctions**, not raw mapping
score. The pain points are (1) **short terminal exons** that minimap2/gapmm2 miss and deSALT/uLTRA
recover, (2) **homopolymer indel noise** in the A-run immediately 5' of the CPA, and (3) **junction
wobble** that correct-first / Module 2H must clean up downstream.

**Stack reality (feasibility anchor).** RECTIFY shells out to the `minimap2` *binary*
(`multi_aligner.py::run_minimap2`, L246) with an `extra_args` passthrough, *not* mappy. deSALT is a
vendored binary; mapPacBio is JVM; uLTRA uses namfinder (strobemers) + minimap2; gapmm2 wraps
mappy/minimap2. So the cheapest interventions are **CLI-flag / index-build changes routed through the
existing wrappers**; the most expensive are new Python seeding engines. There is no in-process seed
API to monkey-patch — everything is subprocess + BAM/PAF parsing.

---

## ⚠️ Ground-truth caveats that constrain every proposal below

Two findings from `03_adversarial/` bound what we can claim and how we must validate:

1. **The win-rate causal story is not load-bearing.** Production winner selection
   (`merge_corrected_tsvs`, both call sites) passes **no per-aligner BAMs**, so it runs the **legacy
   5-level popularity sort** (`_five_rescued, _chimera_ok, _conf_rank, _n_agree, _span, _n_junc`) —
   *not* HP-edit-distance and *not* 3'-end accuracy. The "deSALT 78.9% because cross-read junctions"
   chain is an UNTESTED hypothesis on a single yeast DRS sample with unconfirmed pre/post-`index_col`-fix
   provenance. **Therefore: do NOT justify a seeding change by "it will raise aligner X's win rate."**
   Justify it by a *ground-truth concordance* metric (below) measured independently of selection.

2. **There is no committed accuracy oracle yet.** The single most valuable enabling work for *any*
   seeding proposal is a **held-out 3'-end / junction truth set** (NET-seq CPA peaks — already bundled
   in `data/`; Quant-seq / 3'-seq; annotated intron coords for spliced reads). Every proposal's
   validation section assumes this oracle exists; building it (Proposal 0) gates the rest.

**Validation metrics used throughout** (all *selection-independent*, per-aligner, on the truth set):
- **`seed_exon_recall`** — fraction of true terminal/internal exons that received ≥1 seed (the direct
  short-exon-miss metric; measure by counting exons with zero anchors in the chain).
- **`cpa_mae`** — median absolute error of the *corrected* 3' end vs NET-seq/3'-seq CPA.
- **`junction_concordance`** — fraction of N-ops whose boundaries match an orthogonal junction truth
  within ±0 / ±2 bp.
- **`unmapped_rate`, `wall_time`, `peak_rss`** — cost guards.

---

## Proposal 0 (ENABLER) — Commit a selection-independent accuracy oracle + per-aligner seed/exon-recall harness

**Mechanism.** A small harness that, for a fixed read set, runs each aligner, and computes
`seed_exon_recall`, `cpa_mae`, `junction_concordance` against (a) bundled NET-seq tables, (b) the GFF
intron coordinates, on reads whose true structure is known. Crucially it logs, per read, **which
annotated exons received an anchor** (parse the BAM CIGAR + the chain, or instrument the aligner's
seed dump where available — minimap2 `--print-seeds`, namfinder NAM output).

**Benefit / failure mode.** Not a seeding change itself — it is the measuring instrument. Without it,
every proposal below is unfalsifiable (per caveat 1). Directly addresses "win rates are not accuracy."

**Feasibility.** HIGH. Pure Python + pysam + existing wrappers. minimap2 binary already supports
`--print-seeds`; deSALT/mapPacBio expose no seed dump, so for those infer exon coverage from final
CIGAR (sufficient for `seed_exon_recall` at exon granularity).

**Risk/cost.** Low. Main cost is curating the truth set. Risk: NET-seq CPA peaks are an *assay*, not
absolute truth — treat as concordance, not accuracy.

**Established vs novel.** ESTABLISHED methodology (standard aligner benchmarking); the contribution is
committing it to the repo so numbers stop being CLAUDE.md assertions.

---

## Proposal 1 (HIGH IMPACT / LOW EFFORT) — Syncmer / open-syncmer seeding for minimap2 via `--print-seeds`-validated k/w + denser sketch

**Mechanism.** minimap2's minimizers select the *window-minimum* k-mer, so seed spacing is
data-dependent and a window can place its only seed *outside* a short terminal exon → the exon gets no
anchor and is soft-clipped. **Syncmers** (Edgar 2021) and **open-syncmers** select a k-mer as a seed
iff its minimal-hash *s*-mer (s<k) sits at a fixed offset within the k-mer — a *context-free*,
position-independent decision. This gives (a) more uniform seed spacing (no long seedless stretches),
(b) higher conservation under substitution/indel error than minimizers at equal density, and (c)
*window-guarantee* coverage of short exons. Published results (Shaw & Yu 2021, Dutta et al. 2022) show
open-syncmers improve seed *conservation* by ~10–20% at matched density on ~10–15% error reads — exactly
the ONT regime.

**The pragmatic catch + the real proposal.** minimap2 does **not** expose a syncmer flag. So the
*shippable* version is two-tier:
- **1a (this release, trivial):** systematically tune `-k` / `-w` for terminal-exon recall using the
  Proposal-0 harness. RECTIFY already forces `-k14` (the documented DRS recipe); sweep `-k {11,12,13,14}`
  × `-w {3,5}` and pick the (k,w) that maximizes `seed_exon_recall` on short exons without inflating
  `unmapped_rate`/`wall_time`. Lowering w to 3 (denser minimizers) is the closest in-tool approximation
  to syncmer uniformity. **Zero code beyond `extra_args`.**
- **1b (novel, opt-in aligner):** add a Python `run_syncmer_seeder` that computes open-syncmer anchors
  (a ~50-line numpy/`mmh3` kernel), chains them with the existing collinear-chaining logic, and emits a
  PAF that flows through the *same* `_paf_to_bam` path gapmm2 already uses (sequence injection,
  cs→CIGAR). This is a *seeding-only* aligner whose niche is **short-exon recall**, complementing deSALT.

**Benefit / failure mode.** Directly attacks **short terminal-exon miss** (the minimap2/gapmm2 gap).
Expected: +5–15% `seed_exon_recall` on exons ≤20 bp from 1a alone; 1b targets the residual ≤10 bp tail.

**Feasibility.** 1a: HIGH (flag sweep). 1b: MEDIUM — reuse `_paf_to_bam`; the syncmer kernel and a
banded base-DP (edlib, already a dep) are the new code (~1–2 wk).

**Risk/cost.** 1a: negligible (could marginally raise runtime at w=3). 1b: a from-scratch aligner has a
long correctness tail (CIGAR edge cases, minus-strand RC) — the same bugs the dossiers catalogue for
gapmm2. Gate 1b behind validation showing it beats 1a's recall.

**Established vs novel.** Syncmers = ESTABLISHED. Tuning (1a) = ESTABLISHED. A standalone syncmer
*aligner inside RECTIFY* (1b) = NOVEL engineering, established primitives.

**Validation.** Proposal-0 harness; ship 1a if `seed_exon_recall` ↑ and `cpa_mae`/`unmapped_rate` not
worse. Promote 1b only if it recovers exons 1a still misses *and* its `junction_concordance` ≥ deSALT.

---

## Proposal 2 (HIGH IMPACT / MEDIUM EFFORT) — Homopolymer-collapsed (HPC) seeding *for the 3'/poly-A region only*

**Mechanism.** ONT's dominant residual error (even R10.4.1) is **homopolymer-length miscall** —
indels inside A-runs immediately 5' of the CPA. In HPC space, every maximal homopolymer run collapses to
one symbol (`AAAAA → A`), so seeds become *invariant* to HP-length error: a read with `AAAAAA` and the
genome with `AAAAAAA` produce the **same** HPC seed. Assemblers (HiCanu, mdBG/minimap2 `map-pb`'s
`MM_I_HPC`) use this to defeat HP noise. **The benefit for RECTIFY is seed *anchoring*, not the final
alignment**: HPC seeds let the chainer place an anchor *across* the noisy A-run so the terminal exon is
seeded even when its A-tail length is wrong, then the base-DP runs in *normal* space to keep CPA precision.

**Why minimap2's global `map-pb` HPC is the wrong tool.** The dossiers correctly note HPC *hurts* ONT
splice mapping when applied genome-wide (collapses real sequence signal, breaks junction k-mers). So the
proposal is **scoped HPC**: build a *second, HPC-collapsed* minimizer/syncmer index used **only to
recover anchors in the 3' poly-A-proximal window** of reads that minimap2 soft-clipped. Concretely:
a Python pass that (i) finds reads soft-clipped at the 3' end near a genomic A-run (RECTIFY already
detects these for Module 2G), (ii) re-seeds the clipped tail + flank in HPC space against an HPC index
of the local genome, (iii) emits a corrected anchor that feeds the existing soft-clip-rescue / walk-back
modules with a *seed-supported* 3' boundary instead of a heuristic walk.

**Benefit / failure mode.** Attacks **HP-indel-driven 3'-end drift** and **A-run terminal-exon
soft-clipping** (the minimap2/deSALT failure mode). Gives Module 2E/2G a *seed-anchored* CPA rather than
a pure reference walk-back, reducing over/under-correction in long A-runs. Expected `cpa_mae` ↓ on reads
ending in/near genomic A-runs (the hardest CPA class).

**Feasibility.** MEDIUM. The HPC transform + a k-mer hash over the *yeast* genome is tiny (Python/numpy;
yeast is 12 Mb). Hardest part is the run-length↔genomic-coordinate back-mapping (collapsing loses
position; keep a run-length table to invert). Integrates at the Module 2G boundary, not the aligner.

**Risk/cost.** HPC anchors are ambiguous in A-rich regions (many genomic A-runs) — must constrain the
HPC search to the *local* alignment neighborhood (±intron/±gene window), never genome-wide, or false
anchors will pull CPAs into the wrong A-run. This is the central correctness risk; the local-window
constraint is mandatory.

**Established vs novel.** HPC seeding = ESTABLISHED (Canu, mdBG, minimap2). **Scoped HPC re-seeding of
only the poly-A-proximal soft-clip to feed a downstream 3'-end corrector = NOVEL** (no published tool
does HPC seeding *selectively* for CPA recovery; this is RECTIFY-specific).

**Validation.** On reads whose true CPA is known and which end in/near genomic A-runs: compare `cpa_mae`
of (current 2E/2G walk-back) vs (HPC-anchored 2E/2G). Must show ↓ `cpa_mae` with no increase in
gross-mislocalization (CPA moved >50 bp).

---

## Proposal 3 (MEDIUM IMPACT / MEDIUM EFFORT) — Strobemer / fuzzy-seed aligner as a short-exon specialist (productionize what uLTRA already proves)

**Mechanism.** **Strobemers** (Sahlin 2021) link two or more shorter *strobes* sampled from a window
by a hash-minimization rule, producing a single seed that spans a longer reference span while tolerating
mutations/indels *between* strobes — i.e. a **fuzzy, gapped seed**. NAMs (non-overlapping approximate
matches from strobemers, via namfinder) are uLTRA's *current default seeder* and are a published reason
uLTRA recovers micro-exons. Because strobes are short, at least one strobe can land inside a short
terminal exon even when no full k-mer matches; because the seed is gapped, an indel in the exon does not
destroy it. Strobealign (Sahlin & Marçais 2024) shows strobemer seeding matches/beats minimizer
sensitivity on noisy data at comparable speed.

**Proposal.** Expose namfinder (already vendored with uLTRA) as a **standalone seeder** decoupled from
uLTRA's annotation index: NAMs → existing collinear chaining → `_paf_to_bam`. This gives a **de-novo,
annotation-free** strobemer aligner — capturing uLTRA's micro-exon strength *without* its GTF
dependency, poly-A stripping, and operational brittleness (the dossiers' three biggest uLTRA complaints).

**Benefit / failure mode.** **Short / micro-exon recall** without annotation — fills the gap between
deSALT (de-novo but minimizer-grade short-exon recall) and uLTRA (great short-exon recall but
annotation-bound, no CPA). Independent failure mode → genuine ensemble value.

**Feasibility.** MEDIUM. namfinder binary already ships; wrapper + chain + PAF→BAM reuses gapmm2's path.
Risk concentrated in chaining strobemer NAMs correctly (overlap handling) and in CPA handling (namfinder
has no poly-A model — must NOT strip poly-A, unlike uLTRA, so the tail reaches Module 2E/2G).

**Risk/cost.** MEDIUM. Strobemer parameter space (strobe count/length/window) needs tuning; wrong
params regress to worse-than-minimizer. Another from-scratch chainer = correctness tail.

**Established vs novel.** Strobemers/NAMs = ESTABLISHED (and already in-tree via namfinder). A
**de-novo, poly-A-preserving strobemer aligner in RECTIFY's panel** = NOVEL composition of established
parts.

**Validation.** Proposal-0 harness restricted to multi-exon reads with short internal/terminal exons:
must match uLTRA's `seed_exon_recall` while keeping a non-stripped 3' end (`cpa_mae` defined, unlike
uLTRA's which is undefined post-strip).

---

## Proposal 4 (MEDIUM IMPACT / MEDIUM EFFORT) — Unified de-novo + annotation-prior index (RdBG breadth × GTF/CPA priors)

**Mechanism.** The four index structures trade off differently:
- **minimap2 minimizer hash** — fast, de-novo, but samples (can miss exons).
- **deSALT RdBG (deBGA)** — graph of collapsed unipaths; *one lookup returns every genomic copy* →
  breadth-of-candidates (paralogs, repeats, short exons) — the structural root of its sensitivity.
- **GMAP sampled-oligomer hash + IIT** — fixed-stride sampling (deterministic), IIT carries optional
  splice/SNP priors.
- **uLTRA annotation-segment index** — exact annotated exon/junction coordinates, but caps at the GTF.

**The unified idea.** No single index simultaneously gives (a) de-novo breadth and (b) annotation/CPA
priors *as soft anchors*. Build a composite index: an **RdBG-style / de-Bruijn de-novo seed layer**
(breadth) **augmented with an auxiliary interval index (IIT-like)** carrying *annotated exon boundaries,
annotated junctions, and empirical CPA/NET-seq sites* as **priors that seed candidate anchors but never
gate** (exactly RECTIFY's Module 2H "tie-breaker not gate" philosophy, and GMAP's IIT/canonical-as-reward
design). At seed time: query the de-novo layer for breadth, AND inject annotated-exon and CPA-site
anchors into the chain (à la uLTRA's "inject seedless annotated exons") — but de-novo reads still chain
freely, so novel structure is never lost.

**Why this best supports cross-read junction consensus.** deSALT's win is *cross-read exon pooling* —
projecting all reads' match blocks to consensus boundaries. A graph index that exposes all genomic copies
of a unipath is what makes that pooling cheap. Adding annotation+CPA priors to that same graph gives the
pool a *seeded prior* at known junctions/CPAs without forcing reads onto them. This is the index structure
most aligned with RECTIFY's downstream consensus + Module 2H.

**Benefit / failure mode.** Combines deSALT-style sensitivity (short/repeat exons) with uLTRA-style
junction sharpness *and* a 3'-end-aware seed prior (Proposal 5), in one index — attacking *all three*
failure modes at the seed stage rather than post-hoc.

**Feasibility.** LOW–MEDIUM (this is the ambitious one). Building a new graph index in Python is large;
the *pragmatic* route is to **not rebuild deSALT** but to (i) keep deSALT as the breadth layer, and (ii)
add the annotation/CPA prior as **injected `--junc-bed`-style anchors + a CPA anchor file** consumed by a
thin chaining post-pass. I.e. realize the "unified index" as a *layered ensemble + shared prior file*,
not a monolithic new data structure.

**Risk/cost.** Highest of the set. A true unified graph index is a multi-month effort with deSALT's
fragility risks (NFS temp files, dedup, SIGSEGV on yeast GTF). The layered approximation is far cheaper
but is "ensemble glue," not a new index — manage expectations accordingly.

**Established vs novel.** RdBG, IIT, annotation injection are each ESTABLISHED. **A single index unifying
de-novo graph breadth with annotation + empirical-CPA soft priors = NOVEL** (no surveyed tool does it; it
is the natural index for RECTIFY's correct-first + Module-2H design).

**Validation.** Compare the layered prototype against the current ensemble on `seed_exon_recall`,
`junction_concordance`, AND `cpa_mae` jointly. Only worth the full graph build if the layered version
shows the priors materially help all three.

---

## Proposal 5 (HIGH IMPACT / LOW–MEDIUM EFFORT) — 3'-END-AWARE seeding: anchor candidate alignments near annotated/empirical CPA sites

**Mechanism.** Today no aligner in the panel uses *3'-end priors* at seed time — CPA accuracy is purely a
*post-hoc correction* problem. But RECTIFY ships exactly the priors needed: the GFF, and three NET-seq
tables (`saccharomyces_cerevisiae_netseq_{wt,pan}.tsv.gz`, `..._atract_netseq.tsv.gz`). Proposal: build a
**CPA-anchor index** — a sorted interval set of annotated 3' ends + empirical NET-seq CPA peaks — and use
it to **bias seeding/extension toward terminating the terminal exon at a known CPA**, as a *soft anchor*
(reward, never gate, mirroring `--junc-bonus`/Module 2H). Two delivery routes:
- **5a (cheap, today):** generate a **3'-end "anchor BED"** analogous to the existing `--junc-bed`
  junction BED, and pass annotated/empirical CPA coordinates into the soft-clip-rescue / walk-back
  modules as *candidate stop sites*: when a soft-clipped or HP-ambiguous 3' end is within ±W of a CPA-BED
  site, prefer that site (tie-break), seeding the corrected CPA on real biology instead of a pure
  reference walk. This is a *seeding-of-the-corrector* change, fully inside RECTIFY's Python.
- **5b (novel):** a true seed-time CPA anchor — inject a synthetic anchor at the CPA site into the chain
  so the terminal exon is *extended to* the CPA during DP (GMAP-style joint boundary search, but for the
  unspliced 3' end instead of a junction).

**Benefit / failure mode.** Attacks **CPA imprecision** — the core deliverable — at its source. Where
the genome offers several near-equal 3' columns (HP runs), a CPA prior breaks the tie toward the
biologically observed site. Expected `cpa_mae` ↓, especially on the A-run-terminal class that defeats
every aligner.

**Feasibility.** 5a: HIGH — reuses the junc-BED machinery and Module 2E/2G hooks; the CPA-anchor file is
trivially built from bundled NET-seq + GFF. 5b: MEDIUM — needs anchor injection into a chainer (pairs
naturally with Proposal 1b/3/4's custom chain).

**Risk/cost.** **Circularity is the dominant risk** (the adversarial reports' uLTRA "annotation-snapping
= assumed correct" warning applies in full). A CPA prior that *gates* would manufacture concordance with
the very annotation we validate against. Mitigations: (i) prior is a **bounded tie-break only** (±W,
small), never overrides a strongly sequence-supported end; (ii) **validate on novel/un-annotated CPAs**
(reads whose true 3' end is NOT in the prior set) to prove it doesn't just snap to known sites; (iii) keep
NET-seq as concordance, not truth.

**Established vs novel.** Using annotated/assay priors as soft alignment hints = ESTABLISHED (minimap2
`--junc-bed`, GMAP IIT). **3'-end/CPA-specific seed/anchor priors from NET-seq = NOVEL** — no aligner does
3'-end-aware seeding; this is the most RECTIFY-distinctive idea here.

**Validation.** Two-arm: (a) annotated CPAs — `cpa_mae` ↓ expected; (b) **held-out novel CPAs** — `cpa_mae`
must NOT degrade and gross-mislocalization must not rise (the anti-circularity test). Ship only if (b) is
clean.

---

## Proposal 6 (LOW–MEDIUM IMPACT / LOW EFFORT) — Error-model-aware / spaced seeds tuned to R10.4.1; + A/B-test deSALT's missing `-x ont2d`

**Mechanism.** Two distinct, cheap levers:
- **6a — deSALT preset A/B (the dossiers' #1 flagged item).** deSALT runs with **no `-x`** → the `null`
  (~13% error) model on ONT DRS data, never compared against `-x ont2d`. The seed parameters
  (`-l,-s,-a`) and gap/Z-drop are preset-dependent. A one-line A/B (`-x ont2d` vs default) on the
  Proposal-0 harness either confirms `null` is fine or unlocks free sensitivity. **Pure config.**
- **6b — spaced / error-profiled seeds.** Spaced seeds (fixed don't-care positions in the seed mask)
  are provably more sensitive than contiguous seeds at fixed weight under a substitution model; tuning
  the mask to R10.4.1's *error spectrum* (substitution-dominant at certain contexts, HP-deletion-dominant
  in runs) can raise seed hit-rate on noisy reads. minimap2 has no spaced-seed flag, so 6b lands in a
  custom seeder (Proposal 1b/3). Lower priority than 6a.

**Benefit / failure mode.** 6a: potential free sensitivity / better short-exon recall from deSALT — the
panel's strongest aligner — at zero engineering cost. 6b: marginal seed-sensitivity gains on the noisiest
reads.

**Feasibility.** 6a: HIGH (config + harness run). 6b: MEDIUM (needs custom seeder; defer).

**Risk/cost.** 6a: trivial — worst case, the A/B shows no change and we document `null` is intentional
(closing the dossiers' open question). 6b: spaced-seed design is finicky; learned masks risk overfitting
one chemistry (the penalty tables are already "R10.4.1-specific" — same overfitting hazard).

**Established vs novel.** Spaced seeds = ESTABLISHED (PatternHunter, 2002). Read-type presets =
ESTABLISHED. R10.4.1-*learned* seed mask = NOVEL but speculative and overfit-prone.

**Validation.** 6a: harness diff of `seed_exon_recall`/`cpa_mae`/`junction_concordance` with/without
`-x ont2d`; commit the result either way to close the flag. 6b: only after a custom seeder exists.

---

## Ranked summary

Ranking is by **expected impact ÷ effort**, conditioned on Proposal 0 existing (it gates measurability).

| Rank | Proposal | Mechanism (1-line) | Primary failure mode fixed | Status | Impact | Effort | Risk |
|---|---|---|---|---|---|---|---|
| — | **0. Accuracy oracle + seed-recall harness** | Selection-independent truth set + per-aligner exon/CPA/junction metrics | *(enabler — unblocks all)* | ESTABLISHED method | — | Low | Low |
| **1** | **6a. deSALT `-x ont2d` A/B** | Use the ONT error preset instead of `null` | Sensitivity / short-exon recall on the best aligner | ESTABLISHED | Med–High | **Trivial** | Trivial |
| **2** | **1a. minimap2 k/w sweep (syncmer-approx)** | Denser, more uniform minimizers via tuned `-k`/`-w` | Short terminal-exon seed miss | ESTABLISHED | Med–High | Low | Low |
| **3** | **5a. CPA-anchor BED (3'-end-aware corrector seeding)** | NET-seq/GFF CPA sites as bounded tie-break in 2E/2G | CPA imprecision in HP/A-runs | NOVEL (priors estab.) | High | Low–Med | Med (circularity) |
| **4** | **2. Scoped HPC re-seeding of poly-A-proximal soft-clips** | HPC-invariant anchors across A-runs feeding 2E/2G | HP-indel 3'-end drift; A-run exon clip | NOVEL (HPC estab.) | High | Med | Med (local-window bound) |
| **5** | **3. De-novo strobemer aligner (namfinder, poly-A-preserving)** | Fuzzy gapped seeds; uLTRA recall w/o annotation/strip | Micro/short-exon recall, annotation-free | NOVEL composition | Med | Med | Med |
| **6** | **1b. Standalone open-syncmer seeder** | Context-free seeds → uniform coverage of short exons | Residual ≤10 bp exon miss | NOVEL eng. | Med | Med | Med |
| **7** | **5b. Seed-time CPA anchor injection** | Synthetic CPA anchor extends terminal exon in DP | CPA imprecision (seed-level) | NOVEL | Med–High | Med | Med (circularity) |
| **8** | **4. Unified de-novo + annotation/CPA-prior index** | Graph breadth × IIT-style soft priors in one index | All three, at seed stage | NOVEL | High (long-term) | High | High |
| **9** | **6b. R10.4.1-learned spaced seeds** | Error-profiled seed mask | Noisy-read seed sensitivity | Speculative | Low–Med | Med | Med (overfit) |

**Sequencing recommendation.** Do **0 → 1 (6a) → 2 (1a) → 3 (5a)** first: an enabler plus three
low-effort, high-leverage wins that need no new aligner. Treat 2 (HPC) and 3/6 (strobemer) as the next
tier of genuine seeding novelty, and 4 (unified index) as a research bet justified only if the priors in
5a/2 demonstrably help — and only via the *layered* approximation before any monolithic graph build.

**Cross-cutting risk (applies to 5 and 4).** Any annotation/CPA prior must be a **bounded reward, never
a gate**, and every such proposal must include the **novel-site arm** in validation, or we will
re-manufacture the annotation-circularity the adversarial reports already flagged.
