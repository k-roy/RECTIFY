# RECTIFY native aligner + its gating benchmark — overview

*A living document, plain-language first, headed for the README in some form.
Revisit and polish as the work evolves. Working state:
`dev/HANDOFF_ALIGNER_BENCHMARK.md`; technical spec:
`dev/SIMULATION_BENCHMARK_SPEC.md`; design seed:
`dev/ALIGNER_IDEATION_SYNTHESIS.md`; design refinements:
`dev/ALIGNER_INVESTIGATION_SYNTHESIS.md`.*

**Last updated:** 2026-06-29

---

## The one-sentence version

We are building a **new, native RECTIFY alignment algorithm** that complements
the existing panel of aligners with a deliberately *orthogonal* approach — one
that recovers **novel transcript isoforms** the current aligners flatten away and
that is robust to the **specific error modes of Nanopore direct-RNA sequencing
(DRS)** — and, *before* we build it, we are building a **ground-truth benchmark**
that can prove the new approach actually earns its place.

---

## Background: what RECTIFY does, and why aligners matter

Nanopore (ONT) sequencing reads RNA directly (direct-RNA, "DRS") or as cDNA. The
first thing any analysis must do is **align** each read to the genome — and for
RNA, find its **splice junctions** (where introns were removed). Getting the
junctions right is what lets you discover transcript **isoforms** (alternative
splice forms of a gene); getting the **3′ end / cleavage-and-polyadenylation
(CPA)** site right is core to RECTIFY's mission.

RECTIFY doesn't trust one aligner. It runs a **panel** (minimap2, gapmm2, uLTRA,
deSALT, gmap, …), corrects characteristic ONT errors, and forms a **consensus**,
picking the best-supported alignment per read. Different aligners fail in
different ways, so combining beats any one alone. (This per-read multi-aligner
arbitration idea is the PI's published **COMPASS** method, Roy et al. 2023 NAR,
now living inside RECTIFY.)

But the panel has **shared blind spots** — and that is the opening for a new
algorithm.

---

## The gap we are filling — the main point

Three facts define the opportunity:

1. **The panel "herds."** Several panel members share the same underlying
   engine family — a **flat, affine-gap, quality-blind** cost model. So they tend
   to make the *same* mistakes and agree with each other for the *wrong* reason.
   Counting their agreement therefore overcounts confidence (the "herd trap"). A
   new member only adds real information if its errors are **independent** of that
   shared family.

2. **Novel isoforms get flattened.** Existing aligners are biased toward the
   *expected*: minimap2 **snaps** a real but non-canonical splice junction onto
   the nearest canonical `GT-AG` motif (or onto an annotated junction) instead of
   reporting the true, novel site — and it does this *even on error-free reads*.
   So genuine discoveries get silently turned into known isoforms.

3. **DRS has its own error modes.** Direct-RNA reads mis-length **homopolymers**
   (`AAAA…`), carry errors in **bursts** on globally **"hot"** reads (not
   uniformly), and have hard-to-place **poly-A / CPA** ends. Generic affine-gap
   aligners handle these poorly.

### How the new algorithm is *orthogonal* (the design's central lever)

The unifying move (four of six independent ideation lenses converged on it):
**replace hard / flat / quality-blind costs with calibrated likelihoods, and emit
posteriors.** Instead of an integer "best score," the native aligner scores
placements on an **empirical −log P penalty scale** (RECTIFY already ships such
HP/error tables) and produces **probabilities**. This is simultaneously:

- **the orthogonality source** — its error axis differs from the panel's shared
  flat-affine family, so its mistakes are independent and it genuinely *de-herds*
  the consensus; and
- **a structural defense against scoring artifacts** — see the 0.09→1.07 story
  below.

### What the new algorithm *is*, mechanically

Not a 6th correlated placer. It is primarily a **native realignment / arbitration
layer that runs *downstream* of the panel**: reuse the 5 aligners' placement
cluster as a **localization window**, then do **local realignment** inside that
window on a calibrated scale. A key insight bounds the ambition honestly: **the
discovery ceiling is at the *window* level, not the junction level** — a read in
the right window but with a mis-called junction *inside* it is still recoverable
by realignment; only reads with **no acceptable window at all** (all aligners
misplace; ~12% are unmapped) are out of reach, and those are the *only* place an
independent localizer earns its keep.

---

## Why a benchmark has to come first — "the GATE"

To claim the new aligner *improves* anything you must measure it against a
**known-correct answer**. With real biological data you don't know the right
answer — that's the whole problem you're solving. So we build a benchmark with
**ground truth** first, and enforce a hard rule:

> **No native-aligner ("member") code is written until a validated ground-truth
> benchmark proves, against truth, that the new approach beats the incumbent on
> the case it targets.**

Two disciplines make the gate trustworthy:

- **Fitness is the truth set, NEVER the internal score.** The cautionary tale:
  re-weighting the *internal* consensus score once flipped an aligner's apparent
  quality from **0.09 to 1.07 with no change to any alignment**. An internal
  score can be gamed; truth cannot. (Tuning the error model against *simulated*
  errors and trusting a green number is the same trap one level up —
  "hill-climbing into the simulator's error model.")
- **A stratum only counts if the incumbent is BELOW CEILING on it *and* the gap
  is ADDRESSABLE by the proposed member.** A benchmark where the incumbent already
  scores 100% cannot separate the concepts (the "vertical-slice finding": an
  isolated, cleanly-flanked homopolymer is non-discriminating — *both* minimap2
  and the candidate DP score 1.000, because any in-run indel placement is
  ambiguity-equivalent). The real signal lives at the **hard boundary cases**
  (indel-vs-substitution at run edges, run-bleeds-into-flank, adjacent runs,
  background noise) — and we guard against the opposite error too (the "paralog
  zero-evidence trap": a case the incumbent fails but *no* method could recover
  doesn't count either).

### The framing metric

**Exact indel-position concordance with truth — not edit distance.** At every
contested position edit distance is *tied by construction*, so it can't separate
the concepts; only *which tied placement matches truth* can. Scoring is
**ambiguity-aware**: a call one base into a donor/acceptor repeat, or an indel
anywhere inside a homopolymer run, is *not* charged as an error — only genuine
mistakes count.

### Two kinds of ground truth (we use both, deliberately)

- **Simulation** — reads where we *know* the truth because we built the
  transcript, then injected errors. Absolute per-base truth — but only as
  realistic as the injected error model. Validates **placement mechanics**.
- **Real spike-ins (SIRV / Sequins / ERCC)** — synthetic RNA of *known sequence*
  spiked into a real run → **real ONT errors on known sequences**. The gold
  standard for whether our *simulated* errors are realistic, and the only truth
  valid at homopolymers and native CPA.

A simulation win is **necessary, not sufficient** — it must *transfer* to real
data. This is the separate, complementary **Deliverable B**: real-data
corroboration, whose orthogonal junction truth set is the **COMPASS short-read
detector** (an independent, non-ONT view of which junctions are real).

### Two-tier benchmark

- **Tier 1 — controlled micro-benchmark:** hand-built mini-loci with known truth
  per failure mode → *discriminates the concepts* (where position-exact
  concordance is scored). Light enough to run on a laptop.
- **Tier 2 — realistic transcriptome simulation:** whole-transcriptome reads with
  per-read origin → global novel-junction recall/FDR and **sizing the
  panel-failure tail**. Heavy → runs on the cluster.

---

## The targeted capabilities (the "C-facets") — the new aligner's to-do list

Each is a benchmark stratum **and** a candidate member capability; each had to
clear three bars to qualify — **orthogonal**, **dependency-light**, and **has a
position-exact ablation**:

| Facet | Status | Capability | Mechanism | Incumbent weakness it targets |
|---|---|---|---|---|
| **C1** | ✅ built + Claim-A-proven (del-only; ins gated off after real SIRV caught it) | Homopolymer / STR indel correction | calibrated **HP-length-law** emission cost wired into the gap recurrence (−log P(obs_run\|true_run)) | flat affine misplaces indels out of the run / "repairs" a mismatch with a spurious indel |
| **C2** | ❌ refuted *as placement* (gate, 2026-06-29): shipped guarded walkback already at ceiling on the identifiable genomic-A drift; a soft CPA posterior is deferred behind C3 | 3′ poly-A **CPA** placement | 2-state templated-vs-tail **change-point** under the A-run length law (joint localize+refine) | 3′ ends drift; genomic-A tracts confound |
| **C3** | ⬜ **the active next facet** (the keystone) | Calibrated arbitration | refiner emits **posterior + runner-up**; consensus compares paths by **likelihood ratio (LLR)**, not integer-max | hard, quality-blind scores → the 0.09→1.07 artifact |
| **C4** | ⬜ future | Paralog / multi-copy loci (e.g. SMN1/SMN2) | **POA-pooled** per-locus consensus (cluster → consensus → align once → project back) | per-read placement ambiguous at the lone distinguishing base; mis-clustering |
| **C5** | ⬜ future (gated on a measured tail) | The **panel-failure tail** | **FracMinHash** containment fallback localizer — the only mechanism for reads with no acceptable window; **gated** behind a measured depletion trigger | reads no incumbent places are invisible to a single-aligner run |
| **C6** | ⬜ future (stratum built + discriminating) | Variant-aware junctions | variant/haplotype-aware emission | a deletion near a splice site gets re-expressed as a spurious intron |
| **Discovery** | ✅ guardrail prototype done (WS-3; soft prior) | Novel-junction recovery | evidence-weighing instead of motif/annotation snapping; novel-call support must not be read-quality-tail-enriched | the headline isoform-flattening bias |

A spin-off research idea threads through these: a read's **"hotness"** (how
error-prone it is, *estimable from its own well-aligned exon regions*) as a new
signal to **down-weight unreliable reads when discovering novel junctions** — but
only as a **soft** down-weight (a read can be clean in its exons yet bursty at the
junction), and only if the benchmark proves the FDR lift.

### Read-quality structure & the "novel-feature support" principle (pan-dataset view)

A guiding, dataset-agnostic principle for both the aligner *and* whole-pipeline QC,
distilled from looking at real ONT error structure:

- **Reads are not one population.** Plotting each read's error rate often shows
  *structure* — frequently a **bimodal** distribution: a large peak of decent-to-good
  reads plus a smaller, well-separated **error-rich tail**. (A clean way to see it:
  fit a two-component Gaussian mixture to per-read error rate and ask whether the
  minor peak separates.) That minor peak is itself worth labeling — is it
  **contamination** (reads that aren't really from the target — they map better to
  another reference) or **genuinely error-rich reads** (degraded molecules / poor
  pores)? Tracking that split **across runs and chemistries** is a window into
  flow-cell/pore behavior (e.g. RNA002 vs RNA004), and is exactly the lens needed to
  read the RNA004 "hot tail" finding (is it a separable minor peak, and which kind?).

- **For building the aligner, we can largely *disregard* the low-Q tail.** The job
  is robust non-canonical splice sites and 5′/3′ ends recovered from the **bulk of
  decent-to-good reads.** The error-rich tail neither defines nor should drive those
  calls.

- **The principle that makes this actionable — novel features must be supported by
  decent reads in proportion to the overall quality spectrum.** If a dataset is, say,
  90% decent / 10% low-Q and that error is spread roughly uniformly across reads,
  then a *real* novel isoform (novel junction or novel 3′ end) should **sample from
  that same spectrum** — i.e. be supported by ~90% decent reads. If instead a putative
  novel 3′ end or splice site is **90% low-Q and only 10% decent** — and especially
  if the errors are **enriched in exactly the read segments that dictate that novel
  call** — be **immediately suspicious**: it is far more likely a low-quality artifact
  than real biology. This is a per-read-reliability-weighted FDR control on discovery,
  and it operationalizes the "hotness" idea above into a concrete, testable check:
  *novel-feature support, stratified by read quality, must not be tail-enriched.*

This is the high-level, pan-dataset assessment lens for both the de-novo aligner and
the whole pipeline — a benchmarkable metric (novel-call support vs the read-quality
spectrum) and a guardrail against discovering isoforms that exist only in the noise.
(Caveat to honor: distinguishing a *genuine* low-prevalence isoform that happens to
sit in hard-to-sequence sequence from a noise artifact is itself subtle — the test
is a *prior*, not a hard gate, consistent with the soft-down-weight discipline.)

### Design discipline: what was *rejected*, and why (so it isn't re-proposed)

The orthogonality bar is a real gatekeeper. A WFA-banded engine as a standalone
member was rejected — it "shares minimap2's affine optimum," so it is enabling
*infrastructure*, not an orthogonal concept. Pangenome/variation-graph,
strobemer reseeding, and r-index localization were rejected as
dependency-violating or paradigm-renames. Several figures motivating the work
(e.g. uncorrected-3′-drift %, the count of recurrent GMAP-only novels) are
explicitly **unverified** — they are exactly what the benchmark and Deliverable B
must *measure*, not assert.

---

## Where we are now (2026-06-29)

**The gate (benchmark) is built, green, and red-teamed** — it generates controlled
strata (homopolymers, junctions, paralogs, variants), simulates reads, runs the
incumbent aligner, and scores against truth (ambiguity-aware; separate annotated/
novel and canonical/non-canonical tracks).

**The first aligner facet, C1, is built and proven** — an empirical in-run gap
discount that measurably improves homopolymer indel *placement* against truth
(hand-verified at the CIGAR level), without hallucinating on clean reads. Notably,
**real SIRV data caught a bug in C1** (an insertion-discount that hallucinated
indels, 3–7% on two independent real datasets); we gated it off, and the
deletion-only law is what ships — the spike-in grounding working exactly as
intended. A measurement bug (introns inflating the per-read error span) was also
fixed so error rates are now per-exonic-base.

**Real-data grounding is live, with one open puzzle.** Four independent real
spike-in sources are measured. Two RNA002-era SIRV sources agree the clean,
mod-free error clustering is *low*; the current-chemistry **RNA004** source looked
divergently *burstier* — but the signature points to a **hot-read subpopulation /
contamination artifact**, not a chemistry effect, and a competitive-re-alignment
job (WS-1) is settling it now. The injector's magnitude knobs remain
**placeholders** pending that clean number.

**Claim B was correctly ruled out on real SIRV.** The per-length error-*shape*
cannot be honestly validated on SIRV (too few long homopolymers; and a
"truth-confounded" trap where the only discriminating reads have unknowable
per-read truth). The genuine validation path — a held-out injection simulator — is
specified and deferred (multi-night). This is the adversarial process catching a
would-be fake win, not a setback.

**Three workstreams are in flight** (see Priorities): WS-1 (settle RNA004 + a fair
four-way + the biological-excess contrast), C2 (the next facet), and WS-3 (the
novel-feature discovery guardrail).

---

## Priorities (2026-06-29) — the scoping rule for directing this work

Two tracks run in parallel. **The facet track is the deliverable; the realism
track exists only to keep the gate trustworthy.** The one-line rule (advisor):
***aim, not breadth*** — facet-building is the goal; realism work must not crowd it
out, and the RNA004 puzzle needs one competitive re-alignment, not a swarm.

**Track A — aligner facets (the product):**
1. **C1** ✅ done (homopolymer indel placement).
2. **C2** ❌ **refuted as a placement facet** (gate, 2026-06-29): the shipped guarded
   walkback is already at ceiling on the identifiable genomic-A drift, so there's no
   placement headroom to build into — the gate killed a redundant facet. Tellingly,
   the CPA *uncertainty* (`ambiguity_range`) already ships **with no consumer**, so a
   soft CPA posterior only earns its keep once C3 exists → this independently
   promotes C3.
3. **C3 — now THE active next facet, the keystone** (calibrated LLR arbitration):
   C1 (and any future soft CPA/discovery signal) only improves the *pipeline* if the
   consensus can *use* calibrated, probabilistic output instead of integer-max. C3
   converts per-read facet wins into consensus wins, consumes the already-shipped
   `ambiguity_range`, and is the structural fix for the 0.09→1.07 artifact.
4. **C4 / C5 / C6** — gated, later (each needs its enabling measurement: C5 a
   measured panel-failure tail; C4/C6 strata already exist and discriminate).
- **Discovery guardrail** (novel-feature support, the read-quality-tail FDR prior):
  prototype done (WS-3); folds into C3 and the discovery-FDR machinery.

*Pattern worth noting:* the gate has now **refuted or deferred three plausible
directions** (C1's insertion-discount, the SIRV length-shape "Claim B," and
C2-as-placement) and **confirmed one** (C1 deletion placement). That hit rate is the
prove-don't-assert discipline working as intended — most plausible ideas don't
survive contact with truth, and catching that *before* shipping code is the point.

**Track B — error realism / gate calibration (in service of the gate):**
1. **Settle RNA004** (WS-1, in flight): competitive alignment to strip
   contamination → real-vs-artifact; plus the fair four-way and the same-sample
   SIRV-vs-endogenous *biological-excess* contrast (the systematic-error-vs-biology
   axis).
2. **Calibrate the injector** to the cleaned real targets once WS-1 lands (the
   magnitude is still a placeholder).
3. **Genuine length-shape validation (Claim B):** the held-out injection simulator
   on natural long-homopolymer templates at SIRV-measured rates. Multi-night.

---

## The end game

1. A **trustworthy, ground-truth benchmark** that can decide whether any aligner
   change is a real improvement on the targeted capabilities.
2. A **realistic, calibrated read simulator** (anchored to SIRV absolute truth,
   cross-checked against a second simulator and real reads) so simulation results
   transfer to real data — with real-data corroboration (Deliverable B / COMPASS)
   as the independent check.
3. A **native RECTIFY aligner** — an orthogonal, calibrated-likelihood
   realignment/arbitration layer (plus a gated FracMinHash fallback for the tail)
   — that **complements the panel's blind spots**, measurably improving recovery
   of **novel isoforms** and handling of **DRS error modes**, shipped only where
   the benchmark proves the gain.
4. Generalization from the **yeast** proving ground to **human** loci (SMN1/SMN2
   and other paralog / novel-isoform-rich regions), with organism differences (RNA
   modifications, longer poly-A) accounted for — one engine, per-organism parameter
   sets.

The throughline: **every claimed improvement is backed by measurement against
truth, not by intuition.**

---

## Mini-glossary

- **Panel / consensus** — the set of existing aligners RECTIFY runs, and the
  per-read arbitration (COMPASS-style) that picks the best.
- **Member / native aligner** — the new orthogonal algorithm being added.
- **Orthogonal** — making *independent* errors vs the panel's shared flat-affine,
  quality-blind family (so it adds real information / de-herds the consensus).
- **Junction classes** — ANNOTATED (in the reference), NIC (novel combination of
  known sites), NNC (novel site). Recall and FDR are tracked separately for
  canonical vs non-canonical.
- **CPA** — cleavage-and-polyadenylation site (the transcript 3′ end).
- **Tier 1 / Tier 2** — controlled discrimination vs realistic external validity.
- **SIRV / Sequins / ERCC** — synthetic spike-in RNAs of known sequence
  (absolute truth).
- **Deliverable A / B** — the simulation benchmark (this GATE) / real-data
  corroboration (COMPASS orthogonal junction truth).

---

*Provenance: written by the benchmark-builder agent on branch
`worktree-agent-a25a2c1e784ad37dc`; render date 2026-06-29. Sources in this repo:
`dev/SIMULATION_BENCHMARK_SPEC.md`, `dev/ALIGNER_IDEATION_SYNTHESIS.md`,
`dev/ALIGNER_INVESTIGATION_SYNTHESIS.md`, `dev/HANDOFF_ALIGNER_BENCHMARK.md`, and
the RECTIFY architecture docs. Figures flagged "unverified" above are pending
measurement by the benchmark / Deliverable B.*
