# GMAP — Source-Level Investigation (cDNA / Iso-Seq Spliced Aligner)

**Tool:** GMAP (Genomic Mapping and Alignment Program)
**Primary paper:** Wu TD, Watanabe CK. *GMAP: a genomic mapping and alignment program for
mRNA and EST sequences.* **Bioinformatics 21(9):1859–1875 (2005).** DOI 10.1093/bioinformatics/bti310.
**Methods chapter:** Wu TD, Reeder J, Lawrence M, Becker G, Brauer MJ. *GMAP and GSNAP for
Genomic Sequence Alignment: Enhancements to Speed, Accuracy, and Functionality.* **Methods in
Molecular Biology 1418:283–334 (2016).** PMID 27008021.
**Source / docs:** http://research-pub.gene.com/gmap/ ; unofficial archive
https://github.com/juliangehring/GMAP-GSNAP (`src/gmap.c`, `src/stage1.c`, `src/stage2.c`,
`src/stage3.c`, `src/intron.c/.h`, `src/dynprog*.c`, `util/gmap_build.pl.in`).

**Role in RECTIFY:** **NOT in the current 5-aligner panel** (minimap2, mapPacBio, gapmm2, uLTRA,
deSALT). GMAP is investigated as (a) the **historically dominant annotation-FREE spliced
aligner** — the de-facto PacBio Iso-Seq mapper before minimap2 — and (b) an **algorithmic
contrast**: its splice placement is driven by an explicit *sequence-based splice-site model*
plus *sandwich DP*, not by minimizer-chain heuristics (minimap2) or a de-Bruijn index (deSALT)
or annotation (uLTRA). It is a candidate aligner for the panel.

Claim tagging, consistent with the sibling dossiers:
- **[FACT-PAPER]** — stated in Wu & Watanabe 2005 or the 2016 MiMB chapter.
- **[FACT-SRC]** — verifiable in GMAP source (`gmap.c` defaults), the `gmap` man page, the
  GMAP README, or RECTIFY source.
- **[INFERENCE]** — my analytical synthesis; not asserted by the authors.

> **Defaults caveat.** Two default sources disagree because they are different vintages: the
> Trusty/Xenial **man page** gives `--chimera-margin=40` and `--microexon-spliceprob=0.90`;
> the archived **`src/gmap.c`** gives `chimera_margin = 30` and `microexon_spliceprob = 0.95`,
> and `maxintronlen = 200000` (man page `-K` default text says 1000000). Numbers below are
> tagged with the source they came from; re-verify against the installed binary's `--help`.

---

## 1. Index & Oligomer Seeding

### What is built

**[FACT-SRC/README]** `gmap_build` preprocesses the genome into a set of binary files. The core
of the index is an **oligomer (k-mer) hash / offset structure** plus a positions array — a
*sampled* genomic k-mer lookup table, not a suffix tree.

- **Default k-mer size `-k = 15`** (allowed ≤ 16; historically 12–15). **[FACT-SRC, man page
  `--kmer`; README]**
- **Default base size `--basesize = 12`** — GMAP uses a two-level scheme: the *offset* table is
  keyed on the 12-mer "base" prefix, and the remaining 3 positions (`15 − 12`) are resolved by a
  compressed "gamma pointer" layer. **[FACT-SRC/README]**
- **Default sampling interval `--sampling = 3`** (allowed 1–3) — the genome's oligomers are
  **sampled every 3 bp** rather than at every position. **[FACT-SRC/README]** This is the
  "minimal sampling" idea of the 2005 paper: store ~1/3 of positions, accept a slightly larger
  search to recover full sensitivity at lookup time.
- **Index file types [FACT-README]:** `<genome>.ref<...>gammaptrs` (~64 MB for hg, k=15),
  `<genome>.ref<...>offsetscomp` (~350 MB), positions file, plus modern variants
  `<genome>.sarray` (enhanced suffix array) and `ref12153bitpackptrs/bitpackcomp`
  (SIMD-bitpacked compressed hash). **IIT** ("Interval Index Tree") files hold splice-site /
  intron / SNP auxiliary data for SNP- and splice-tolerant modes.

### Memory footprint

**[FACT-README]** RAM to **build** the index scales with k: k=12 → 64 MB, k=13 → 256 MB,
k=14 → 1 GB, **k=15 → 4 GB**. The on-disk compressed index for a mammalian genome is a few
hundred MB (gammaptrs 64 MB + offsetscomp ~350 MB). **[INFERENCE]** For RECTIFY's
*S. cerevisiae* R64 (~12 Mb) the index is trivially small (a few MB) and builds in seconds —
the large-genome memory cost is a non-issue at yeast scale, exactly as for deSALT.

### The sampled-k-mer tradeoff

**[FACT-PAPER, "Minimal sampling strategy"]** Two design levers set sensitivity/speed:
**(i) long oligomers (k=15)** for high *specificity* — a 15-mer hit is almost always a true
local match, so few spurious anchors enter chaining; **(ii) sparse sampling (every 3 bp)** for
low *memory and lookup cost*. The paper frames this as minimizing the number of oligomer
lookups needed to map a cDNA reliably, using an **adaptive sampling scheme** (sample more
densely only when the sparse pass is ambiguous). **[INFERENCE]** This is the same
specificity-vs-density tradeoff minimap2 exposes via `(w,k)`-minimizers — but GMAP samples on a
*fixed stride* of the genome rather than picking a window-minimum, so its seeds are positionally
deterministic, not hash-selected.

**[FACT-PAPER, 2016 chapter]** The 2016 release re-engineered this layer: **compressed genomic
hash tables with SIMD (SSE/AVX) fast access**, **linear-genome SIMD representation**,
**enhanced suffix arrays (ESAs)**, and **>4 Gbp large-genome** support — speed/memory upgrades
that did not change the core three-stage algorithm.

---

## 2. Three-Stage Algorithm (map → refine/splice-DP → align)

GMAP's source organizes alignment into **Stage 1 / Stage 2 / Stage 3** (`stage1.c`,
`stage2.c`, `stage3.c`; the `gmap.c` comments and timers `Stage2_compute`, `Stage3_compute`
confirm the naming). **[FACT-SRC]** Conceptually: *find candidate genomic regions → connect
seeds into an exon/intron skeleton with splice-aware DP → produce base-level alignment + CIGAR.*

### Stage 1 — approximate genomic mapping (oligomer hashing + diagonalization)

**[FACT-PAPER, "Genomic mapping"]** Each sampled query 15-mer is looked up in the genomic hash,
yielding candidate genomic positions. GMAP **diagonalizes** these hits — groups them by their
(genomic − query) offset (the *diagonal*), so all seeds belonging to one ungapped exon segment
fall on (nearly) the same diagonal. Diagonals with enough support define **candidate genomic
regions** ("chromosomal segments") where the cDNA plausibly originates. **[FACT-PAPER, 2016
chapter: "GMAP uses diagonalization to find exon regions."]**

**Ranking of regions/paths [FACT-PAPER + FACT-SRC].** Multiple candidate regions can survive
(paralogs, gene families, chimeras). GMAP scores each by seed coverage and produces up to
**`-n / --npaths = 5`** reported paths (`maxpaths_report = 5` in `gmap.c`; `-n 0` ⇒ report 1 if
non-chimeric, 2 if a chimera is detected). **[FACT-SRC]** The best path is the primary
alignment; lower-ranked ones are alternates.

### Stage 2 — refinement: oligomer chaining + sandwich DP (the distinctive stage)

**[FACT-PAPER, "Oligomer chaining" + "Sandwich DP"]** Within a candidate region, Stage 2 turns
the scattered same-diagonal seeds into an ordered **chain** of high-scoring segment pairs (HSPs
= ungapped exon blocks), allowing diagonal *jumps* between them. A jump where the genomic offset
grows much faster than the query offset is a **candidate intron**; a small balanced jump is an
indel.

The defining mechanism is **sandwich DP** (the paper's term; the literature also calls it
"attack from both sides"). When two anchored exon segments flank an unresolved gap that should
contain an intron, GMAP runs **two dynamic programs simultaneously from both anchored ends
toward the middle** — one extending rightward from the upstream exon's 3′ edge (the donor side),
one extending leftward from the downstream exon's 5′ edge (the acceptor side) — and solves for
the single intron position that **maximizes (alignment score + splice-site score)**. **[FACT-PAPER]**

**[INFERENCE]** This is the structural contrast with every other panel aligner: minimap2/deSALT
place the intron as a *long gap op* inside one left-to-right banded DP (cost driven by the
two-piece affine `q2`/`e2` gap and a GT-AG dinucleotide penalty); GMAP instead *brackets the
gap between two anchored exon segments and searches the boundary jointly from both sides*. The
sandwich formulation is what lets GMAP slide the exon/intron boundary to the position that is
simultaneously the best base alignment **and** the strongest splice signal, even when basecall
errors sit right at the junction.

`-K / --intronlength` bounds the intron a sandwich gap may span (man page default text 1000000;
`gmap.c` `maxintronlen = 200000`, used in Stage 2 and `Stage3_mergeable`). `--min-intronlength`
(default 9, man page) sets the threshold below which a genomic gap is called a **deletion**
rather than an intron — i.e. GMAP's own N-vs-D decision boundary. **[FACT-SRC/man page]**

### Stage 3 — base-level alignment, CIGAR, microexons, chimeras

**[FACT-PAPER, "Sandwich DP / nucleotide-level alignment"; 2016 chapter]** Stage 3 runs
**nucleotide-level DP** (`dynprog*.c`) to resolve the remaining mismatches, indels, and exact
intron boundaries within the Stage-2 skeleton, producing the final base-perfect alignment and
the SAM CIGAR (with `N` for introns). The 2016 chapter summarizes it as "dynamic programming at
the nucleotide level to resolve mismatches, indels and intron boundaries."

- **SNP / indel tolerance [FACT-PAPER].** A headline 2005 result: GMAP "generates accurate gene
  structures even in the presence of substantial polymorphisms and sequence errors." The DP
  absorbs mismatches and small indels without breaking the exon chain, and a SNP-tolerant mode
  (IIT-driven) can align against known variants. This robustness is *why* GMAP was usable on
  noisy/divergent ESTs and later on PacBio reads.
- **Microexon detection [FACT-PAPER + FACT-SRC].** Very short internal exons (down to a few bp)
  are missed by seed-based mapping. GMAP has a dedicated microexon search that inserts a tiny
  exon **only if a flanking splice-site probability exceeds `--microexon-spliceprob`**
  (`gmap.c` default **0.95**; man page 0.90). **[FACT-SRC]** This statistical-significance gate
  is the same competency uLTRA markets (exons < 30 nt) but GMAP achieves it *de novo* from the
  splice model, not from annotation.
- **Chimera / read-through detection [FACT-SRC].** If a margin of query sequence
  (`--chimera-margin`; `gmap.c` **30**, man page 40) cannot be aligned to the same locus, GMAP
  searches for the remainder elsewhere and reports a **chimeric (two-part) alignment** — used
  for fusion/read-through transcripts. `-n 0` forces 2 paths when a chimera is found.
  **[FACT-SRC]** GMAP/GSNAP also support **circular-chromosome** alignment for organelle/plasmid
  references. **[FACT-README]**

---

## 3. Sequence-Based Splice-Site Model (the algorithmic signature)

This is GMAP's defining property and the main reason it merits a contrast dossier.

**[FACT-PAPER]** GMAP places introns using a **sequence-based splice-site model** evaluated on
the *genome*, **without requiring annotation and without (originally) a learned probabilistic
model**. The 2005 abstract is explicit: it generates accurate gene structures "**without using
probabilistic splice site models**." Instead, the sandwich DP adds a **reward** when the chosen
intron boundary yields a recognized splice dinucleotide pair, and selects the
exon/intron boundary that maximizes *combined* (base-alignment + splice) score.

### Canonical vs semi-canonical vs non-canonical

**[FACT-SRC, `intron.h`; man page `--canonical-mode`]** GMAP recognizes three splice-signal
classes, in priority order:
1. **Canonical `GT…AG`** (the dominant U2 intron).
2. **Semi-canonical `GC…AG`** and **`AT…AC`** (U12/minor and GC-AG variants).
3. **Non-canonical** (anything else) — allowed, but not rewarded.

`--canonical-mode` controls the reward magnitude: **0 = low, 1 = high (default), 2 =
conditional** (low reward for high-identity sequence, high reward otherwise). **[FACT-SRC/man
page]** The "conditional" mode is the SNP-aware refinement: when the read aligns nearly
perfectly, GMAP trusts the base alignment and *down-weights* the canonical bonus (so it won't
bend a clean alignment to manufacture a GT-AG site); when the read is noisy, it *up-weights* the
splice prior to pull the ambiguous boundary onto a real signal. **[INFERENCE]** This is a more
adaptive policy than minimap2's fixed non-canonical penalty `-C` or deSALT's fixed `-R 9`.

### Splice strength / probability and `--cross-species`

**[FACT-SRC/man page]** GMAP scores the *strength* of donor/acceptor signals (position-specific
sequence content around the dinucleotide), and the microexon gate above is expressed as a
**splice-site probability** threshold (0.90–0.95). For divergent targets, `--cross-species`
switches to "a more sensitive search for canonical splicing," widening the boundary search so a
splice site shifted by polymorphism is still found. **[FACT-PAPER, 2016 chapter]** Later
releases integrated **MaxEnt (maximum-entropy) scoring of novel splice sites** — first in GSNAP,
then ported to GMAP — adding a genuinely probabilistic donor/acceptor model on top of the
original reward scheme. So the "no probabilistic model" claim is true of the 2005 core but
*superseded* by the MaxEnt-enabled modern builds. **[FACT-PAPER, 2016 chapter]**

**[INFERENCE]** Net contrast for RECTIFY: GMAP's junction placement is **splice-model-first**
(boundary chosen to jointly maximize base score + sequence splice strength, sandwich-bracketed
from both sides), where minimap2 is **gap-cost-first** (intron emitted when the two-piece affine
break-even is crossed, GT-AG only as a per-position bonus) and uLTRA is **annotation-first**
(exon coordinates supplied by the GTF). GMAP sits closest in *spirit* to RECTIFY's own Module 2H
junction refiner, which is also sequence-first and canonical-aware — but GMAP does it *inside*
the aligner, per read, at alignment time.

---

## 4. De-novo (annotation-free) vs annotation-guided; determinism; threading

### Annotation-free operation

**[FACT-PAPER]** GMAP is **fundamentally annotation-free**: it discovers exon/intron structure
purely from genome sequence + the splice model. Optional known-splicesite / known-intron IIT
files can *assist* (and GSNAP's `-s`/`-N` flags govern known-vs-novel splicing), but the core
GMAP cDNA workflow needs only the genome. **[INFERENCE]** Contrast:
- **uLTRA** is annotation-*required* (builds candidate-exon collinear chains from the GTF;
  excels at small exons *because* it knows where they are). GMAP gets small-exon (microexon)
  competence *without* annotation via its splice-probability gate — a different route to the
  same capability.
- **minimap2** is heuristic/annotation-optional (`--junc-bed` is only an additive bonus). GMAP's
  sandwich DP is a stronger *de novo* splice solver but far slower (Section 5).
- **deSALT** is also de novo (cross-read exon pooling). GMAP refines each read *independently*
  against the genome+model, with **no cross-read consensus** — a key difference from deSALT's
  batch exon-inference, and **[INFERENCE]** a reason deSALT can produce more *homogeneous*
  junctions across reads while GMAP optimizes each read in isolation.

### Determinism

**[INFERENCE]** With fixed sampling stride, a deterministic diagonalization, and DP tie-breaks,
GMAP's output is **deterministic per read** (no random minimizer hashing, no
unordered-thread-dependent merge). Thread count does not change alignments — only throughput.
This is a mild advantage over hash-seeded mappers for reproducibility. (Not asserted in the
paper; inferred from the fixed-stride/DP design.)

### Multithreading model

**[FACT-SRC/README]** GMAP parallelizes **per read**: `-t/--nthreads` worker threads each align
whole reads independently; the README notes total threads = 2 + `-t` (one reader, one writer).
**Each individual alignment is single-threaded** (the Stage 2/3 DP is not internally
parallelized). **[INFERENCE]** So GMAP scales with read count across cores but cannot speed up a
single hard alignment — and its per-read DP cost is high, so wall-time on large long-read sets is
dominated by the slow per-read path even when fully threaded. This is the practical root of the
speed gap with minimap2 (Section 5/6).

---

## 5. Strengths

- **[FACT-PAPER] Historically very accurate splice sites, de novo.** The sandwich-DP +
  sequence-splice-model placement produces precise canonical/semi-canonical junctions without
  annotation — the property that made GMAP the reference Iso-Seq mapper for a decade.
- **[FACT-PAPER] Strong SNP/indel and divergence tolerance.** Accurate gene structures "even in
  the presence of substantial polymorphisms and sequence errors"; SNP-tolerant IIT mode and
  `--cross-species` extend this to variant/cross-organism targets. **[INFERENCE]** Directly
  relevant to noisy long reads — GMAP was usable on PacBio CLR/CCS precisely because its DP
  absorbs the error.
- **[FACT-SRC] De novo microexon detection** gated by splice probability (0.90–0.95) — recovers
  few-bp internal exons without a GTF (uLTRA-class capability, annotation-free).
- **[FACT-SRC] Chimera / read-through / circular detection** built in.
- **[FACT-PAPER] No annotation dependency** and **deterministic** output **[INFERENCE]** — a
  clean, reproducible baseline.
- **[FACT-SRC] Adaptive canonical reward** (`--canonical-mode 2`, conditional) — won't bend a
  clean alignment to fake a canonical site, but rescues noisy junctions onto real signals.

## 6. Weaknesses (speed / scaling / 3′-end / homopolymer)

- **[FACT-PAPER/benchmark] SPEED is the dominant weakness.** GMAP runs **nucleotide-level
  sandwich/Stage-3 DP per read**, single-threaded per alignment. Independent benchmarks put
  **minimap2 ≥ ~30× faster** than GMAP-class long-read mappers (Li 2018), and the GraphMap2
  evaluation found **GMAP the slowest tool tested, ~3–7× slower than GraphMap2** on long
  RNA-seq. **[INFERENCE]** Against minimap2's SIMD KSW2 + minimizer chaining, GMAP is roughly
  **one to two orders of magnitude slower** on whole long-read datasets — the single reason
  minimap2 (2018) displaced it as the Iso-Seq standard.
- **[INFERENCE] Scaling on large long-read sets.** Per-read DP cost × millions of reads makes
  GMAP wall-time-prohibitive at modern Iso-Seq/ONT depths even when fully threaded; throughput is
  read-parallel only, with no within-alignment SIMD acceleration of the DP in the classic path.
- **[FACT-README] Build memory** 4 GB for k=15 (mammalian); irrelevant at yeast scale but a
  cost on large genomes (cf. deSALT's ~35 GB human index — different magnitude, same "heavy on
  big genomes" theme).
- **[INFERENCE] Noisy ONT homopolymers / poly-A 3′ ends.** GMAP has **no poly-A awareness** — it
  aligns an in-read A-run as ordinary sequence. On ONT DRS, systematic homopolymer-deletion bias
  lets the banded DP slip the 3′ exon end into a genomic A-run, leaving the reported CPA
  imprecise — *exactly the artifact RECTIFY corrects* (Modules 2B/2E/2G/`find_polya_boundary`).
  GMAP's splice model helps internal junctions but does nothing for the unspliced 3′ terminal
  exon / CPA, where RECTIFY's value lies. (Same limitation noted for deSALT/minimap2; GMAP is no
  better here.) The sandwich DP was designed and validated on **mRNA/EST and PacBio**, not on
  ONT's specific homopolymer error profile.
- **[INFERENCE] Per-read isolation.** No cross-read exon consensus (unlike deSALT) → junctions
  can be placed a few bp apart between reads from the same isoform, the very inhomogeneity
  RECTIFY's correct-first selection penalizes.
- **[FACT-SRC] Default drift across versions.** `chimera_margin`, `microexon_spliceprob`,
  `maxintronlen` differ between man page and source — a reproducibility footgun if defaults are
  assumed rather than pinned.

---

## 7. Iso-Seq Ecosystem Role (and why minimap2 displaced it)

**[FACT — cDNA_Cupcake wiki / Iso-Seq history]** Before minimap2, **GMAP was the default genome
aligner for PacBio Iso-Seq**. It sat **between clustering and collapse** in the ToFU / Cupcake
pipeline: `pbtranscript`/ICE clustering → **GMAP genome alignment** →
`collapse_isoforms_by_sam.py` (Cupcake) → SQANTI. The recommended invocation (also recorded in
RECTIFY's `pacbio_ecosystem.md` §5.1):

```
gmap -d <db> -D <dbdir> -f samse -n 0 --no-chimeras -t <N> input.fasta > aln.sam
```

Key Iso-Seq flags: `-f samse` (single-end SAM), **`-n 0`** (report a single best path; 2 if
chimeric), `--no-chimeras` (≡ `--chimera-margin 0`, disable fusion search for clean isoform
collapse), `-z` (cDNA direction), `-K` (max intron for the target genome). **[FACT-SRC]**

**[FACT — cDNA_Cupcake "Best practice" wiki]** The current Cupcake guidance is explicit:
**"minimap2 is recommended; GMAP is a second recommendation, which sometimes has better
alignments but is much slower."** GMAP, deSALT, STAR(long), and BLAT are all listed as
acceptable, with minimap2 the default. **[FACT — IsoSeq3/SQANTI3 docs]** The modern stack
replaced GMAP entirely with **`pbmm2 align --preset ISOSEQ`** (minimap2 `splice:hq`); collapse
moved into `isoseq collapse`; QC into SQANTI3/pigeon. GMAP also appears in **Cogent** (de novo
coding-genome reconstruction → "fake genome" → GMAP-align poorly-mapped reads) and **TAMA**
workflows as an alternate aligner.

**Why minimap2 won [INFERENCE]:** identical-or-better Iso-Seq junction accuracy on *clean HiFi*
reads (where GMAP's error-tolerance advantage evaporates) at **~30×+ the speed**, with native
SAM/BAM, presets, and active maintenance. GMAP's accuracy edge persisted mainly on *noisy /
divergent / cross-species* inputs — a shrinking niche as HiFi accuracy rose. GMAP remains
offered (Cupcake's "sometimes better") but is no longer the default anywhere in the modern
PacBio stack.

---

## 8. Relevance to RECTIFY

**[INFERENCE throughout this section.]**

1. **As an algorithmic contrast (primary value).** GMAP is the cleanest example of
   **splice-model-first** junction placement among long-read aligners: sandwich DP brackets the
   intron from both anchored exon ends and chooses the boundary that jointly maximizes base score
   + sequence splice strength, with an adaptive canonical reward (`--canonical-mode 2`) that
   *down-weights* the splice prior on clean reads. This is conceptually the **same policy
   RECTIFY's Module 2H** (`junction_refiner.py`) applies *post hoc* — sequence-first scoring,
   canonical as a tie-breaker / prior (RECTIFY's `_CANONICAL_HP_PRIOR = 0.5`), never a hard gate.
   GMAP demonstrates that this philosophy can live *inside* the aligner. RECTIFY's design choice
   — let aligners place junctions, then refine sequence-first downstream — is validated by GMAP's
   approach but decoupled from GMAP's speed cost.

2. **As a candidate aligner — likely a poor fit for the production panel.** GMAP's ~1–2 orders
   of magnitude speed disadvantage is disqualifying for RECTIFY's chunked-alignment HPC workflow
   (the whole `rectify split` mandate exists to make alignment I/O/compute tractable; adding a
   per-read-DP aligner that is 30×+ slower than minimap2 would dominate wall time). Its lack of
   poly-A / ONT-homopolymer awareness means it would *not* fix the 3′-end artifacts RECTIFY
   targets — it would be just another upstream aligner whose junctions RECTIFY's correct-first
   logic must clean.

3. **Where it might still earn a slot:** GMAP's **de novo microexon detection** and **strong
   divergence/SNP tolerance** could add value as an *occasional* aligner for difficult loci
   (small internal exons, polymorphic strains) — competing with uLTRA's small-exon niche but
   **annotation-free**, which matters because RECTIFY's deSALT path already disables annotation
   (`-G` SIGSEGV) and its yeast workflow is GTF-light. If added, expect a **low win rate**
   (like uLTRA's ~2% / gapmm2's ~0.8%) — a specialist tie-breaker, not a workhorse.

4. **Determinism / reproducibility.** GMAP's fixed-stride, hash-free-selection seeding gives
   deterministic per-read output independent of thread count — a minor plus for a validation
   panel, contrasting with minimizer/thread-order effects elsewhere.

5. **Net assessment.** GMAP belongs in RECTIFY's investigation set as the **historical and
   algorithmic reference point** for annotation-free splice alignment, and as the lineage that
   minimap2 displaced. It is **not recommended for the production 5-aligner panel** (speed +
   no 3′-end benefit), but its sandwich-DP splice-model-first design is the strongest external
   precedent for RECTIFY's own sequence-first junction-refinement philosophy.

---

## 9. Source / Paper References

**Papers / chapters**
- Wu TD, Watanabe CK. *GMAP: a genomic mapping and alignment program for mRNA and EST sequences.*
  Bioinformatics 21(9):1859–1875 (2005). DOI 10.1093/bioinformatics/bti310. Sections used:
  *Minimal sampling strategy* (k=15, sampling 3 bp, adaptive sampling), *Genomic mapping*
  (oligomer hashing + diagonalization, Stage 1), *Oligomer chaining* and *Sandwich DP* (Stage 2,
  splice-site rewarded boundary search), *nucleotide-level alignment* (Stage 3, SNP/indel),
  *microexon identification with statistical significance*, *chimera detection*.
- Wu TD, Reeder J, Lawrence M, Becker G, Brauer MJ. *GMAP and GSNAP for Genomic Sequence
  Alignment: Enhancements to Speed, Accuracy, and Functionality.* Methods Mol Biol 1418:283–334
  (2016). PMID 27008021. Used: "diagonalization to find exon regions, oligomer chaining of short
  k-mers, nucleotide-level DP for mismatches/indels/intron boundaries"; **MaxEnt novel-splice
  scoring** (GSNAP→GMAP); SIMD compressed genomic hash tables; linear-genome SIMD; ESAs;
  >4 Gbp genomes.

**GMAP source / docs**
- Genentech GMAP/GSNAP: http://research-pub.gene.com/gmap/ (README — index k/basesize/sampling
  defaults, build-memory table, gammaptrs/offsetscomp/IIT/sarray/bitpack file types, threading).
- Archived source: https://github.com/juliangehring/GMAP-GSNAP — `src/gmap.c`
  (`maxpaths_report=5`, `chimera_margin=30`, `canonical_mode=1`, `microexon_spliceprob=0.95`,
  `maxintronlen=200000`), `src/stage1.c`/`stage2.c`/`stage3.c` (stage structure),
  `src/intron.c/.h` (canonical/semi-canonical GT-AG/GC-AG/AT-AC), `src/dynprog*.c`
  (sandwich/base DP), `util/gmap_build.pl.in`.
- `gmap(1)` man page (Ubuntu Trusty/Xenial, manpages.debian.org): `--kmer`, `--basesize`,
  `--sampling`, `-K/--intronlength` (1000000), `--min-intronlength` (9), `--chimera-margin`
  (40), `-n/--npaths`, `-z/--direction`, `--cross-species`, `--microexon-spliceprob` (0.90),
  `--canonical-mode` (0/1/2), `-f/--format` (samse/sampe), `-t/--nthreads`.

**Benchmarks / ecosystem**
- Li H. *Minimap2.* Bioinformatics 34(18):3094 (2018) — "minimap2 is ≥30× faster than
  mainstream long-read mappers (BLASR, BWA-MEM, NGMLR, GMAP) at higher accuracy."
- Marić J et al. *GraphMap2 — splice-aware RNA-seq mapper for long reads.* bioRxiv 720458 (2019)
  — GMAP the slowest tool, ~3–7× slower than GraphMap2.
- Magdoll/cDNA_Cupcake wiki — *Best practice for aligning Iso-Seq* ("minimap2 recommended; GMAP
  second, sometimes better but much slower"), ToFU/Cogent GMAP usage.
- PacBio IsoSeq3 / SQANTI3 docs — GMAP→pbmm2/minimap2 displacement.
- RECTIFY `dev/aligner_investigation/01_investigation/pacbio_ecosystem.md` §5.1 (GMAP Iso-Seq
  invocation), `minimap2.md`, `deSALT.md`, `uLTRA.md` (panel context, win rates).

---

## 10. Confidence & Caveats

- **High confidence (FACT):** k=15 / basesize=12 / sampling-every-3-bp defaults; build-memory
  table (4 GB @ k=15); three-stage structure; sandwich DP as the splice-placement mechanism;
  canonical/semi-canonical GT-AG/GC-AG/AT-AC recognition and `--canonical-mode 0/1/2`;
  microexon-spliceprob gate; chimera margin; npaths=5; "no probabilistic model" (2005) → MaxEnt
  (2016); annotation-free operation; per-read threading; minimap2 ≥30× faster.
- **Version-sensitive (re-verify with installed `--help`):** `chimera_margin` (30 src vs 40 man),
  `microexon_spliceprob` (0.95 src vs 0.90 man), `maxintronlen` (200000 src vs 1000000 man-page
  `-K` text).
- **From search-snippet extraction (re-verify against the paper PDF):** exact sandwich-DP
  scoring constants and splice-strength formula — the Oxford/PMC full text and Springer chapter
  were 403-blocked; the mechanism is FACT, the numeric splice-score weights were not recovered.
- **All "Relevance to RECTIFY" (§8), the speed-magnitude framing, and the deSALT/minimap2/uLTRA
  contrasts are INFERENCE** — analytical synthesis grounded in the cited GMAP facts plus
  RECTIFY's `CLAUDE.md` and sibling dossiers, not direct quotes from the GMAP authors.
