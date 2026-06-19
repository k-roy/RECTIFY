# deSALT — Source-Level Technical Investigation

> **Build note:** code-level claims here were verified against `master`; see `../CORRECTIONS_vs_DRS_BUILD.md` for re-verification vs `origin/drs-validation-rebuild` (flag sets CONFIRMED — deSALT no `-x`/`-G`; penalty tables present & auto-resolved on the build).

**Tool:** deSALT (De Bruijn graph-based Spliced Aligner for Long Transcriptome reads)
**Version in RECTIFY:** v1.5.6 (vendored Linux/x86_64 binary, 773 KB)
**Repo:** https://github.com/ydLiu-HIT/deSALT (depends on a bundled fork of deBGA)
**Primary paper:** Liu B, Liu Y, Li J, Guo H, Zang T, Wang Y. *deSALT: fast and accurate long transcriptomic read alignment with de Bruijn graph-based index.* **Genome Biology 20, 274 (2019).** PMC6913027; bioRxiv 612176.
**Index dependency:** Liu B, Guo H, Brudno M, Wang Y. *deBGA: read alignment with de Bruijn graph-based seed and extension.* **Bioinformatics 32(21):3224–3232 (2016).**

**Role in RECTIFY:** Tier-2 aligner in the 5-aligner correct-first ensemble. Invoked by
`run_desalt()` in `rectify/core/align/multi_aligner.py` (~line 1506). In RECTIFY's correct-first
win rates, **deSALT wins 78.9% of reads — the single best aligner** (vs mapPacBio 18.2%, uLTRA 2%,
gapmm2 0.8%, minimap2 0.1%). This report explains the algorithmic reasons.

A label is attached to every non-trivial claim:
- **[FACT-PAPER]** stated in the Genome Biology 2019 / deBGA 2016 papers.
- **[FACT-SRC]** verifiable in the deSALT/deBGA source, README, or RECTIFY source.
- **[INFERENCE]** my reasoning connecting the mechanism to RECTIFY's 78.9% win rate. Not asserted by the authors.

---

## 1. RdBG Index Construction

### What is built

deSALT does not build its own index format — `deSALT index ref.fa <index_route>` **wraps deBGA**
to produce a **Reduced de Bruijn Graph index (RdBG-Index)** [FACT-SRC, README: "deSALT requires
deBGA for generation of RdBG-index"]. deSALT ships a patched deBGA whose notable change is
`START_POS_REF = 0` replacing the upstream `START_POS_REF = 2048` [FACT-SRC, README].

**Input constraint** [FACT-SRC, README]: reference FASTA lines must be ≤ 500 bp wide;
`changelinewidth.py` reformats if needed. (Mirrored by RECTIFY's expectation that the index
directory contains a non-empty `ref.seq` — see `_is_built_desalt_index()` in `multi_aligner.py`.)

### How the RdBG is constructed (from deBGA 2016) [FACT-PAPER]

1. **De Bruijn graph of the reference.** deBGA builds a dBG over the genome with a user-defined
   k-mer length. Vertices and edges are derived from all **(k+1)-mers** of the reference.
   Default `-k = 22`, restricted to the 20–28 bp range [FACT-SRC, README `-k/--index-kmer`].
2. **Unipath (unitig) extraction.** deBGA inspects in-/out-degrees of vertices to find the start
   vertices of **non-branching paths (unipaths)**, traverses them, and **collapses each unipath's
   vertices into a single unitig sequence**. This is the "reduction" in *Reduced* dBG: long
   non-branching k-mer chains become one node, so repeated genomic substrings are stored once.
3. **Multi-copy bookkeeping.** For each unitig that occurs multiple times in the genome, deBGA
   records the **start positions of all copies**.
4. **Hash-table index over unipaths.** The RdBG-Index is a **hash-table data structure** indexing
   the *unipaths* (not the raw sequence). It exposes three lookups [FACT-PAPER]:
   (i) k-mer → unipath coordinate + the k-mer's **offset within that unipath**;
   (ii) unipath coordinate → **all genomic positions** of that unipath;
   (iii) unipath coordinate → its unitig sequence.

### Why the RdBG matters for accuracy [INFERENCE]

The hash returns *every* genomic copy of a matched unipath in one lookup. For a spliced read
crossing a repeat or a paralog, this surfaces *all* candidate loci cheaply, instead of relying on
minimizer sampling that can miss the true exon. A unipath that is unique in the genome yields a
single unambiguous anchor, which is exactly the high-confidence seed you want at exon ends and
near the 3′ CPA. This breadth-of-candidates property is the root of deSALT's high sensitivity.

### Index files & memory [FACT-SRC / FACT-PAPER]

- RECTIFY treats `ref.seq` (deBGA's main reference output) as the existence/non-empty sentinel for
  a built index [FACT-SRC, `_is_built_desalt_index()`]. The full deBGA index is a *directory* of
  several files (reference sequence, unipath sequences/offsets, and the k-mer hash); the README
  does not enumerate filenames but distributes prebuilt indexes by memory footprint.
- **Prebuilt index memory** [FACT-SRC, README]: Human GRCh38 ≈ **35 GB**, Mouse GRCm38 ≈ 31 GB,
  Fruit fly DM6 ≈ 3.5 GB. deBGA build of hg19 (k=22): ~4.4 h, ~69.55 GB build memory; alignment
  peak ~40.32 GB for human [FACT-PAPER, deBGA].
- **For yeast (RECTIFY's S. cerevisiae R64, ~12 Mb)** the index is tiny (a few hundred MB) and
  loads/builds in seconds — so the large-genome memory cost (a Weakness in general) is **not a
  practical concern in RECTIFY's yeast workflow** [INFERENCE].

---

## 2. Two-Pass Alignment Algorithm

deSALT's central contribution is a **tailored two-pass strategy** built on the RdBG index
[FACT-PAPER]. The two passes are: (1) build coarse **alignment skeletons** for every read; pool
their exon evidence genome-wide; (2) **re-align each read** against a refined, exon-aware reference.

### Pass 1 — Alignment skeleton generation [FACT-PAPER]

1. **Match Blocks (MBs).** Using the RdBG hash, deSALT finds **match blocks (MBs)** — maximal exact
   matches between the read and genome, recovered as unipath hits. Seeding parameters:
   seeding l-mer `-l = 15`, seed step `-s = 5` bp, local-hash k-mer `-a = 8`, max hits per seed
   `-n = 50` [FACT-SRC, README].
2. **Sparse Dynamic Programming (SDP) chaining.** MBs are composed into an optimized
   **alignment skeleton** by an **SDP** approach: DP-style alignment is applied only *between*
   non-overlapping seeds, estimating edits cheaply rather than doing full base-level DP across the
   whole read [FACT-PAPER]. Chain controls: min chain/skeleton score `-c = 30`, max read gap
   `-g = 2000`, max intron `-I = 200000`, min-frag merge distance `-i = 20`, strand-diff `-d = 10`
   [FACT-SRC, README].
3. The skeleton is a **coarse, gapped placement** of the read across the genome — it pinpoints the
   approximate exon structure (which genomic segments the read spans, with large gaps = candidate
   introns) but is not yet base-accurate.

### Pass 2 — Refined alignment [FACT-PAPER]

1. **Exon inference.** deSALT **integrates ALL skeletons** by projecting their constituent MBs onto
   the reference genome; the projections are analyzed to **infer exon regions**. This is *de novo*,
   **GTF-free** exon discovery (a GTF can optionally seed it via `-G`, but RECTIFY does not use it).
2. **Spliced reference / refined re-alignment.** Each read is re-aligned using the inferred exons.
   deSALT searches for **additional local matches on exons between or near** the skeleton's exons,
   using a **local hash query** with the very small `-a = 8` k-mer to catch **small/short exons**
   and exons with weak seed coverage that Pass 1 missed [FACT-PAPER]. `-e = 5` sets the exon
   extension range [FACT-SRC, README].
3. The result is a **full-length, base-level spliced alignment** that recovers exons and splice
   junctions along the entire read [FACT-PAPER].

### Why two passes beat single-pass long-read mappers [INFERENCE]

A single-pass minimizer mapper (minimap2) commits to one chain per read from sparse seeds; noisy
ONT bases at exon boundaries and short terminal exons routinely cause it to truncate or mis-place
the junction. deSALT's first pass is deliberately *coarse and tolerant*; the second pass then does
the precise work **after** it already knows where exons are — both for this read and (crucially)
from every other read (Section 3). That ordering is why minimap2 wins only **0.1%** of reads in
RECTIFY while deSALT wins **78.9%**.

---

## 3. Cross-Read Junction Refinement (the key accuracy mechanism)

This is the property that most directly explains RECTIFY's win rate.

**[FACT-PAPER]** deSALT's exon-inference step **"integrates all the generated alignment skeletons"**
— it does **not** treat each read in isolation. The MBs from *every read in the batch* are projected
onto the genome, and their pooled projections define the exon regions used to refine *each individual*
read in Pass 2. Batch size `-B = 655350` reads per loop [FACT-SRC, README] sets the pool size.

### What cross-read pooling buys you

- **Consensus exon boundaries.** A single noisy read may have a frayed exon end. When dozens of
  reads from the same locus contribute MBs, their **projections agree** on the true exon/intron
  boundary, and that consensus boundary is applied back to the noisy read. The junction a read is
  snapped to is supported by the *collective* evidence, not just that read's error-prone bases.
- **Small-exon rescue.** A micro-exon may be missed in any one read's seeds but is collectively
  visible across the pool; Pass 2's local-hash query (k=8) then finds it in each read.
- **Homogeneous alignments** [FACT-PAPER]: because all reads at a locus are refined against the
  *same* inferred exon set, reads spanning the same junction are aligned **consistently** (same
  CIGAR junction placement) rather than each landing a few bp apart. The paper explicitly markets
  "sensitive, accurate, and **homogeneous** alignments."

### Why this drives the 78.9% win rate in RECTIFY [INFERENCE]

RECTIFY's correct-first ensemble chooses a winning aligner per read **based on the corrected 3′-end
position**, and its junction-refinement Module 2H rewards reads whose splice junctions are
sequence-supported and consistent. deSALT's cross-read homogeneity is *precisely* the trait that
correction logic rewards:

1. **Sharp, consensus junctions** mean the corrected 3′ end (CPA) sits at a clean exon boundary,
   not in an intron — directly improving the 3′-end metric RECTIFY scores on.
2. **Consistent junction placement across reads** means deSALT's corrected positions cluster
   tightly, so its reads repeatedly beat aligners whose junctions are placed idiosyncratically per
   read (minimap2/gapmm2). This compounding effect — best at exactly the property being scored —
   is why one aligner takes ~4/5 of all reads rather than a more even split.

deSALT effectively performs a mini cross-read junction-correction *internally* before RECTIFY's own
correction even runs, giving it a structural head start over the other four aligners.

---

## 4. Base-Level DP & Splice Scoring

### Scoring scheme [FACT-SRC, README defaults]

| Parameter | Flag | Default |
|---|---|---|
| Match score | `-m / --match-score` | **1** |
| Mismatch score | `-M / --mis-score` | **2** |
| Gap open penalty (two-tier) | `-O / --open-pen` | **2, 32** |
| Gap extension penalty (two-tier) | `-E / --ext-pen` | **1, 0** |
| Z-drop | `-z / --zdrop` | **400** |
| Band width | `-w / --band-width` | **500** |
| Non-canonical splice penalty | `-R / --noncan` | **9** |

### Engine [FACT-SRC / INFERENCE]

- **Dual affine gap penalties** (`-O 2,32` / `-E 1,0`) are the classic **two-cost gap model used by
  minimap2/ksw2**: a cheap gap (open 2, extend 1) for ordinary indels and a near-flat long gap
  (open 32, extend 0) so that **introns** — very long reference gaps — are penalized essentially
  by a constant open cost, not per base. This lets a single banded DP emit a long `N`/intron op
  without the score collapsing [INFERENCE; the parameter shape is the FACT].
- **SIMD-banded extension with Z-drop** (`-w 500`, `-z 400`) is the standard ksw2-style banded
  Smith-Waterman; Z-drop terminates extension when the score drops too far, the mechanism that
  decides where an exon ends and an intron begins [INFERENCE on engine identity; defaults are FACT].
- **GT-AG handling via the non-canonical penalty** (`-R 9`) [FACT-SRC]: canonical **GT-AG** donor/
  acceptor dinucleotides incur **no** penalty; non-canonical junctions are penalized by 9. The DP
  is thus biased toward placing introns at canonical splice signals — a **travel-cost / spliced-SDP**
  formulation where the intron-skip op carries the splice-signal-dependent cost. `-T/--trans-strand`
  restricts splice-site search to the transcript strand.

### Read-type presets (`-x`) tune the error model [FACT-SRC, README]

`ccs` (PacBio CCS, ~1% err), `clr` (PacBio CLR, ~15%), `ont1d` (ONT 1D, >20%), `ont2d` (ONT 2D,
>12%), `null` (default, ~13%). For DRS nanopore data the appropriate preset is an ONT profile;
**RECTIFY invokes `deSALT aln` without an explicit `-x`, so it runs the `null` default (~13% error
model)** [FACT-SRC, `run_desalt()` builds `['aln','-t',…,'-f',…,'-o',…]` plus optional `extra_args`].

---

## 5. Why deSALT Wins (accuracy mechanisms, condensed)

1. **RdBG breadth of candidates** — every genomic copy of a matched unipath is returned in one
   hash lookup, so true exons (including repeated/paralogous and short ones) are rarely missed
   [FACT-PAPER → INFERENCE].
2. **Coarse-then-precise two passes** — commits to base-level junctions only *after* exon
   structure is known, avoiding the truncation/mis-junction failures of single-pass mappers on
   noisy ends [FACT-PAPER → INFERENCE].
3. **Cross-read exon pooling → homogeneous junctions** — the decisive factor: junctions are placed
   by *collective* evidence, producing consistent, sharp boundaries that put corrected 3′ ends at
   clean exon edges — exactly what RECTIFY's correct-first selection scores [FACT-PAPER → INFERENCE].
4. **Small-exon local-hash rescue (k=8)** — recovers micro-exons missed by sparse seeding, so reads
   are full-length and junction-complete [FACT-PAPER].
5. **GT-AG-biased spliced DP (`-R 9`)** — anchors introns at canonical signals, improving junction
   precision and, downstream, 3′-end accuracy [FACT-SRC → INFERENCE].

---

## 6. Strengths

- **Best-in-class sensitivity & junction homogeneity** for noisy long RNA reads; recovers
  full-length exon/junction structure [FACT-PAPER].
- **GTF-free** de novo exon/junction inference — no annotation dependency [FACT-PAPER]. RECTIFY
  relies on this: it deliberately does **not** pass `-G` (yeast GTF triggers a SIGSEGV) and finds
  deSALT's de novo splice detection sufficient [FACT-SRC, `run_desalt()` docstring].
- **Tunable splice scoring** (`-R`, two-tier gap costs) and per-platform presets [FACT-SRC].
- **Tiny index/fast on small genomes** (yeast) — the large-genome memory cost does not apply
  [INFERENCE].

## 7. Weaknesses

- **Memory of the RdBG index on large genomes** — ~35 GB for human GRCh38 [FACT-SRC]; deBGA build
  ~70 GB / ~4.4 h for hg19 [FACT-PAPER]. Heavy for mammalian targets (irrelevant for yeast).
- **Slower than minimap2** — two passes + base-level SIMD DP cost more wall time than a single-pass
  minimizer mapper [INFERENCE; consistent with paper's accuracy-over-speed positioning].
- **Parameter sensitivity** — many knobs (`-k`, `-l`, `-a`, `-s`, `-I`, `-R`, `-x`); wrong preset or
  intron cap (`-I 200000`) can degrade junctions. RECTIFY pins defaults and avoids `-x` [FACT-SRC].
- **Engineering fragility observed in RECTIFY** [FACT-SRC, `run_desalt()` and `_dedup_desalt_bam()`]:
  - emits **each primary alignment N times** → RECTIFY deduplicates post-hoc;
  - **misparses gzipped FASTQ** (reads gz header as text; exit 0 with garbage) and stops at the
    first empty-sequence record → RECTIFY decompresses and cleans the FASTQ first;
  - **`-G` yeast GTF → SIGSEGV** → annotation disabled;
  - **`-f` temp file must be on local (non-NFS) disk** or "double free / corruption" crash → RECTIFY
    writes the temp file under `$TMPDIR`;
  - **conda allocator incompatibility** ("double free or corruption" under parallel launch) →
    RECTIFY strips `LD_LIBRARY_PATH` so it uses system glibc.
- **DRS poly-A / homopolymer tails** [INFERENCE]: deSALT has **no poly-A awareness** — it aligns the
  homopolymer tail like ordinary sequence. On nanopore A/T-homopolymer runs the banded DP can extend
  or slip the 3′ end into a genomic A-run, leaving the reported CPA imprecise. This is *precisely the
  artifact RECTIFY's correction stack exists to fix* (Modules 2B/2E/2G), so deSALT's strong junction
  placement combined with RECTIFY's downstream 3′-end correction is complementary — deSALT nails the
  exon structure, RECTIFY nails the tail.

---

## 8. Source / Paper References

**Papers**
- Liu B. et al. *deSALT: fast and accurate long transcriptomic read alignment with de Bruijn
  graph-based index.* Genome Biology 20, 274 (2019). PMC6913027 / bioRxiv 612176.
  Methods used: *alignment skeleton generation* (match blocks + SDP), *exon inference by integrating
  all skeletons*, *refined second-pass alignment* (local hash query, small-exon detection).
- Liu B., Guo H., Brudno M., Wang Y. *deBGA: read alignment with de Bruijn graph-based seed and
  extension.* Bioinformatics 32(21):3224–3232 (2016). RdBG construction, unipath/unitig extraction,
  RdBG-Index hash table, the three index lookups, memory figures.

**deSALT source / README** (https://github.com/ydLiu-HIT/deSALT)
- Commands: `deSALT index ref.fa <index_route>` (wraps deBGA), `deSALT aln [opts] -f <tmp> <index> <reads>`.
- Default parameters (Section 4) verbatim from README options table.
- deBGA fork note: `START_POS_REF = 0` (was 2048); prebuilt index memory figures.

**deBGA source** (https://github.com/HongzheGuo/deBGA, https://github.com/hitbc/deBGA)
- `deBGA index [options] reference.fasta <index_route>`, `-k` 21–28 bp.

**RECTIFY source** (`/home/user/RECTIFY/rectify/core/align/multi_aligner.py`)
- `run_desalt()` (~L1506): builds `deSALT aln -t -f -o`, no `-x`, no `-G`; FASTQ cleaning;
  local-disk temp file; `LD_LIBRARY_PATH` stripped; SAM→sorted BAM→dedup→calmd.
- `_dedup_desalt_bam()` (~L1480): removes deSALT's duplicate primary alignments.
- `_is_built_desalt_index()` (~L1553): requires non-empty `ref.seq`.
- `rectify/core/commands/install_aligners_command.py`: `DESALT_VERSION = '1.5.6'`, vendored at
  `rectify/data/bin/linux_x86_64/deSALT`, source build from the v1.5.6 GitHub tag.

---
*FACT-PAPER / FACT-SRC claims are traceable to the cited papers, the deSALT/deBGA READMEs, or
RECTIFY source. INFERENCE claims (especially the 78.9%-win-rate causal chain in Sections 2–5 and the
DRS poly-A behavior in Section 7) are my reasoning, not assertions by the deSALT authors.*
