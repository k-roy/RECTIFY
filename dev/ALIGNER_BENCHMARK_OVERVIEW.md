# RECTIFY de-novo aligner + its gating benchmark — plain-language overview

*A living document. Revisit and polish as the work evolves. Pairs with the
working state in `dev/HANDOFF_ALIGNER_BENCHMARK.md` and the technical spec in
`dev/SIMULATION_BENCHMARK_SPEC.md`.*

**Last updated:** 2026-06-28

---

## The one-sentence version

We are building a **new, from-scratch alignment algorithm** that joins RECTIFY's
existing panel of aligners with a deliberately *different* (orthogonal) approach —
one that recovers **novel transcript isoforms** the current aligners miss and that
is robust to the **specific error modes of Nanopore direct-RNA sequencing** — and,
*before* we build it, we are building a **ground-truth benchmark** that can prove
the new aligner actually earns its place.

---

## Background: what RECTIFY does, and why aligners matter here

Nanopore (ONT) sequencing reads RNA molecules directly (direct-RNA sequencing,
"DRS") or as cDNA copies. The first thing any analysis must do is **align** each
read to the genome: figure out where it came from, and — crucially for RNA —
where its **splice junctions** are (the boundaries where introns were removed).
Getting the junctions right is what lets you discover transcript **isoforms**
(alternative ways a gene is spliced).

RECTIFY doesn't trust a single aligner. It runs a **panel** of them (minimap2,
gapmm2, uLTRA, deSALT, gmap, …), then **corrects** characteristic ONT errors and
forms a **consensus** — picking, per read, the best-supported alignment. The
panel is RECTIFY's strength: different aligners fail in different ways, so
combining them beats any one alone.

But the panel has **shared blind spots**. The incumbents were largely designed
for DNA or for short, well-behaved introns, and they lean on built-in assumptions
that hurt exactly the cases we care about most.

---

## The gap we are filling — the *main point* of this line of work

Two blind spots motivate a **new aligner member** with an orthogonal algorithm:

1. **Novel-isoform / novel-junction recovery.** Existing aligners are biased
   toward the *expected*. The clearest measured example: minimap2 **snaps** a
   real but non-canonical splice junction onto the nearest canonical `GT-AG`
   motif (or onto an annotated junction), rather than reporting the true, novel
   site. It does this *even on error-free reads*. So genuinely new isoforms — the
   discoveries that matter biologically — get silently flattened into known ones.
   A new member that **weighs evidence** (the read's own bases) instead of
   **snapping to motifs/annotation** can recover these. (In benchmark terms this
   is the JUNCTION_DISCOVERY stratum; we proved the bias is *member-addressable* —
   the true site always carries strictly more support, so an evidence-weighing
   algorithm can win it back.)

2. **ONT DRS-specific error modes.** Direct-RNA reads have characteristic
   mistakes the incumbents handle poorly: long **homopolymer** runs (`AAAA…`) get
   mis-lengthened; errors come in **bursts** and on globally **"hot"** (error-rich)
   reads rather than uniformly; and the **3′ poly-A / cleavage site (CPA)** is
   hard to place. A member built around **run-/length-aware penalties** and
   **calibrated, error-aware scoring** can correct these where a generic affine-gap
   aligner cannot.

The new aligner is **not** meant to replace the panel — it is meant to
**complement** it: contribute the alignments the others systematically get wrong,
so the consensus improves on novel isoforms and DRS error modes.

---

## Why a benchmark has to come first — "the GATE"

To claim the new aligner *improves* anything, you must measure it against a
**known-correct answer**. But with real biological data you don't know the right
answer — that's the whole problem you're trying to solve. So we build a benchmark
with **ground truth** first, and we enforce a hard rule:

> **No new-aligner ("member") code is written until the benchmark can prove,
> against truth, that the new approach beats the incumbent on the cases it
> targets.** Optimize against a broken ruler and you ship fake improvements.

That benchmark is "Deliverable A, the GATE." Its job is to be *trustworthy* and
*discriminating*: able to tell a real improvement from noise, on each targeted
capability, in an **ambiguity-aware** way (e.g. an indel placed anywhere inside a
homopolymer run, or a junction shifted one base into a repeat, is *not* charged as
an error — only genuine mistakes count).

### Two kinds of ground truth (we use both, for different reasons)

- **Simulation** — generate reads where we *know* the truth because we built the
  transcript, then add errors. Gives absolute, per-base truth — but only as
  realistic as the error model we inject.
- **Real spike-ins (SIRV / Sequins / ERCC)** — synthetic RNA of *known sequence*
  physically spiked into a real sequencing run. Gives **real ONT errors on known
  sequences** — the gold standard for checking whether our *simulated* errors are
  realistic.

The current frontier is making simulated errors **realistic**, then **anchoring
their magnitude to the real spike-ins** so simulation conclusions transfer to real
data.

---

## The targeted capabilities (the "C-facets") — the new aligner's to-do list

Each is a benchmark stratum + a future member capability. The benchmark proves
the incumbent is below ceiling *and* that the failure is addressable, before we
build the fix:

| Facet | The capability | Incumbent weakness the benchmark shows |
|---|---|---|
| **C1** | Homopolymer / short-tandem-repeat indel correction (run-/length-aware cost) | flat affine gap misplaces indels out of the run, "repairs" a mismatch with a spurious indel |
| **C2** | 3′ cleavage-site (CPA) placement on native poly-A | 3′ ends mis-set; genomic-A runs confound |
| **C3** | Calibrated, error-aware posterior (a "soft" abstain band, not a hard gate) | aligners commit hard even on unreliable reads |
| **C4** | Paralog / multi-copy locus resolution (e.g. SMN1 vs SMN2) by pooling weak evidence | per-read placement ambiguous at the lone distinguishing base |
| **C5** | The "panel-failure tail" — reads *no* incumbent places | invisible to a single-aligner run |
| **C6** | Variant-aware junction calling (don't fabricate junctions near real variants) | a deletion near a splice site gets re-expressed as a spurious intron |
| **Discovery** | Novel-junction recovery without motif/annotation snapping | the headline bias above |

A spin-off research idea threads through these: a read's **"hotness"** (how
error-prone it is, *estimable from its own well-aligned regions*) could be a new
signal to **down-weight unreliable reads when discovering novel junctions**,
cutting false discoveries — but only if the benchmark proves the lift, and only as
a *soft* down-weight (a hot read can be clean in its exons yet bursty at the
junction).

---

## Where we are now (2026-06-28)

The benchmark machinery is **built and green**: it generates controlled cases
(homopolymers, junctions, paralogs, variants), simulates realistic reads, runs the
incumbent aligner, and scores against truth — ambiguity-aware, with separate
tracks for annotated vs novel and canonical vs non-canonical junctions.

The active work is **error realism**. We measured that real ONT errors are *not*
uniform — they **cluster in bursts** and concentrate on **hot reads** — and built
a 3-layer error injector to reproduce this. Its magnitude knobs are still
**placeholders**, calibrated to a contaminated upper-bound estimate
(read-vs-genome conflates RNA modifications and alignment artifacts). **SIRV
spike-ins** are how we replace placeholders with the real, clean number — and
settle an open question (does error clustering really sit near 3× or 5× a uniform
baseline?).

**This session's contributions:**
- Built and validated the **exon-GTF loader** — the format SIRV/Sequins/LRGASP
  use — which was the one piece blocking all real-spike-in integration.
- **Submitted the SIRV measurement job** on the cluster: align real ONT DRS reads
  to the SIRV reference, keep the spike-in reads (known sequence = absolute
  truth), and measure their clean error structure → the real calibration target.
- Built the **LRGASP/NanoSim truth join** — lets us bring in a *second*,
  independent simulator. With our pbsim3 + LRGASP's NanoSim + real SIRV we get a
  three-way triangulation: if both simulators miss the same real behavior, that's
  a shared blind spot to engineer around.

---

## The end game

1. A **trustworthy, ground-truth benchmark** that can decide whether *any*
   aligner change is a real improvement on the targeted capabilities.
2. A **realistic, calibrated read simulator** (anchored to SIRV absolute truth,
   cross-checked against a second simulator and real reads) so simulation results
   transfer to real data.
3. A **new de-novo aligner member** with an orthogonal algorithm — evidence-
   weighing junction discovery + run-/length-aware, error-calibrated scoring —
   that **complements the panel's blind spots**, measurably improving recovery of
   **novel isoforms** and handling of **DRS error modes**, shipped only where the
   benchmark proves the gain.
4. Generalization from the **yeast** proving ground to **human** loci (SMN1/SMN2
   and other paralog/novel-isoform-rich regions), with organism differences (RNA
   modifications, longer poly-A tails) properly accounted for — one engine,
   per-organism parameter sets.

The throughline: **every claimed improvement is backed by measurement against
truth, not by intuition.**

---

*Provenance: written by the benchmark-builder agent on branch
`worktree-agent-a25a2c1e784ad37dc`; render date 2026-06-28. Sources: this
codebase's `dev/SIMULATION_BENCHMARK_SPEC.md`, `dev/HANDOFF_ALIGNER_BENCHMARK.md`,
and the RECTIFY architecture docs.*
