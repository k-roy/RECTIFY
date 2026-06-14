# Round-2 Alignment — Designer 1 slice: Isoform Reference Library + Round-2 Alignment

**Author:** Designer 1 (2-designer + adversary panel)
**Date:** 2026-06-14
**Status:** DESIGN ONLY — for a later builder. No production code here.
**Scope:** Items 1–4 of the Round-2 feature. Coordinate lift-back (transcript→genome
realignment) and per-read Round-1-vs-Round-2 scoring are **Designer 2's slice**; this doc
defines the interface (block-map format) Designer 2 consumes but does not implement it.

---

## 0. Framing and the one-paragraph thesis

Round 1 (existing) is **junction DISCOVERY**: a consensus of splice-aware genome aligners
(minimap2, uLTRA, deSALT, mapPacBio, +winnowmap2/gmap being added) finds where introns are,
motif-agnostically. Round 2 (new) is **alignment-quality MAXIMIZATION over already-known
junctions**: build a library of full-length cDNA sequences (introns removed) padded with
genomic flank, and align the subset of reads that Round 1 handled poorly with *non-spliced*
aligners that produce a single contiguous semi-global alignment. Because those aligners no
longer have to *discover* introns, they can place noisy reads better. The recommended way to
keep the library finite and unbiased is to **instantiate only the cDNA paths that reads
actually traverse** (read-evidence-bounded), not to enumerate splice-graph paths
combinatorially.

The load-bearing, adversary-facing assumption is in §1.5: that we can get a read-evidence
isoform set that **includes Round-1 novel/non-canonical junctions** without the WT-annotation
gating that off-the-shelf isoform tools apply by default. Flagged loudly there.

---

## 1. Isoform reference library construction

### 1.1 The arithmetic that kills naive path enumeration (rejected alternative)

Per-gene splice-graph path enumeration is the obvious idea and the wrong one. The cost is not
in the median gene (most human genes have <4 annotated transcripts and the median is small),
it is in the **tail**: a splice graph over N exons with skip/alternative edges has up to
~2^(N−2) source→sink paths. TTN (titin) has ~363 exons; even a tiny fraction of skip edges
makes its path set astronomically large. A "cap of K isoforms/gene" is arbitrary and silently
**discards biology at exactly the genes with the most alternative splicing** — the genes we
most care about in mutants. Path-length bounds and "only paths using ≥1 observed junction"
help but do not bound the blow-up: a single observed exon-skip combined with annotated
alternative segments still multiplies.

**Verdict: REJECT naive enumeration.** It re-introduces an arbitrary bias (the cap) precisely
where RECTIFY's unbiased-discovery philosophy matters most, and it manufactures isoforms no
read supports (which then become spurious multi-mapping targets — see §4.2).

### 1.2 Recommended: read-evidence-bounded isoform instantiation

Instantiate **one cDNA sequence per distinct read-supported exon-chain (intron chain)**,
where the chain may use any junction from `annotation ∪ Round-1-observed`. This bounds the
library at **O(distinct intron chains observed) ≤ O(reads)** — naturally finite, scales with
the data not with combinatorics, and lets novelty enter *only where a read provides evidence*
(no WT-annotation cap, no phantom isoforms).

This is exactly the "collapse" operation that mature long-read isoform-reconstruction tools
already implement (FLAIR, StringTie2 long/mix mode, IsoQuant, Bambu, ESPRESSO). **Build-vs-reuse
recommendation: REUSE the collapse step, do not reimplement it.** RECTIFY already lives next to
FLAIR (the Sumner-lab side of this workspace uses FLAIR→OARFISH), so FLAIR is the path of least
operational friction. The library builder consumes a collapsed transcript model set
(GTF + per-transcript exon structure) and turns it into the padded-cDNA FASTA + block map (§2).

> **Why a graph still appears in the design:** we DO build a per-gene splice graph as the
> internal representation (§1.3) — but as a *deduplication/union/validation* structure over
> observed chains, NOT as an object we enumerate paths from. The graph answers "is this
> observed chain consistent / a sub-path / novel vs annotation?" It is fed; it is not unrolled.

### 1.3 Data structure: per-locus splice graph (representation, not enumerator)

Per gene locus (or per Round-1 connected component of overlapping alignments, to stay
annotation-light):

```
SpliceGraph(locus_id):
  nodes:  exonic segments — boundaries are the SORTED UNION of all exon/junction
          breakpoints from {annotation, Round-1 observed junctions, read TSS/TES}.
          Each node = (genome_start, genome_end, strand).
  edges:  two kinds:
          - ADJACENCY edge  (segment i .. segment i+1 contiguous on genome) = "intron retained"
          - JUNCTION edge   (donor_segment_end -> acceptor_segment_start, with N-op span) =
                            "spliced"; carries source ∈ {annotation, round1}, support_read_count,
                            donor/acceptor dinucleotides (recorded, NOT gated — §1.5), novel flag.
  observed_chains: list of (ordered edge sequence, supporting_read_ids, support_count)
                   — THIS is what becomes a library entry, one per distinct chain above
                   a support floor (default support>=1; see §1.5/§4 for the floor knob).
```

The graph is JSON/parquet-serializable; nodes and edges are integer-indexed; observed_chains
references edge indices. No path enumeration routine exists in this structure — the only
"paths" are the ones already present in `observed_chains`.

### 1.4 Enumeration algorithm + bounds (it's collapse, not enumeration)

```
For each locus L:
  1. Collect all reads whose Round-1 alignment overlaps L.
  2. For each read, extract its intron chain = ordered list of (donor, acceptor) junctions
     from its Round-1 alignment (after the existing consensus + anchor-gate).
  3. Canonicalize the chain (sort by coordinate, snap junctions within +/- W bp (default W=3)
     to the highest-support junction at that locus to absorb wobble — this is the ONLY
     "correction" and it is support-driven, not motif-driven).
  4. Group reads by canonical chain -> each group with support >= floor is one cDNA isoform.
  5. Also emit each ANNOTATED transcript chain for L as an isoform (so well-covered WT
     isoforms are present even if a particular run is shallow) — tagged source=annotation.
  6. De-duplicate identical chains across the annotation/observed union.
```

Bounds, concretely:
- **Per-locus isoform count** = (# distinct observed chains ≥ floor) + (# annotated transcripts).
  Both are data-bounded; no 2^N term anywhere.
- **Library size** ≈ Σ_loci that count. For a human ONT run this is on the order of
  (annotated transcripts ≈ 2.5×10^5 in GENCODE) + (distinct novel chains actually seen,
  typically far smaller because deep coverage collapses to few chains/gene). Order-of-magnitude:
  a few × 10^5 cDNAs — see §4.1 for index cost.
- **Single-exon reads contribute no new isoform** beyond the trivial gene body (they have no
  intron chain); they are also exactly the reads the Round-2 gate excludes (§4.3).

### 1.5 ⚠ LOAD-BEARING ASSUMPTION — how novelty enters WITHOUT WT bias

This is the assumption the adversary should attack hardest.

The whole point is to represent **mutant/pathological exon combinations not in WT annotation**.
Read-evidence instantiation gives us that *only if the collapse step does not silently
re-impose WT/canonical gating*. Off-the-shelf isoform tools violate this by default:

- **FLAIR**: the default pipeline runs `flair correct`, which **snaps read junctions to the
  reference annotation/genome splice sites** and `flair collapse --check_splice` enforces
  GT-AG-ish coverage around junctions and ≥3-read support — exactly the WT bias RECTIFY exists
  to avoid. ([FLAIR modules docs](https://flair.readthedocs.io/en/latest/modules.html),
  [FLAIR issue #371 — non-canonical junctions dropped](https://github.com/BrooksLabUCSC/flair/issues/371))
- **StringTie2 / IsoQuant / Bambu**: all carry reference-guided correction and min-support
  defaults that can drop low-support non-canonical chains.

**Therefore the design does NOT take a tool's default output.** Two acceptable implementations,
both of which must be VERIFIED (by the builder, with a spiked non-canonical junction) before use:

- **(Preferred) Skip the tool's correction step entirely.** Feed RECTIFY's *own*
  Round-1 consensus intron chains directly into the collapse-only step:
  - FLAIR: bypass `flair correct`; run `flair collapse` with `--support 1`,
    **no** `--check_splice`, and `--no_gtf_end_adjustment` so TSS/TES are not snapped to the GTF.
    Verify it does not silently drop non-canonical introns (issue #371 indicates it can —
    hence the spike test is mandatory, not optional).
  - This keeps collapse (the part we want — grouping reads into chains) and discards correction
    (the WT-snapping we do not want).
- **(Fallback) Implement the §1.4 collapse ourselves.** It is ~150 lines: group Round-1 chains,
  snap within W bp to highest-support, emit one cDNA/chain. This is the safest route for the
  unbiased-discovery guarantee because no external default can leak in. Recommend the builder
  start here if the FLAIR spike test shows any non-canonical attrition.

**Junction wobble snapping (W) is the one place bias could sneak in.** It snaps to the
**highest-support junction within W bp**, never to a motif and never to the annotation
preferentially — a novel high-support junction wins over a low-support annotated one. W default
3 bp; expose as a flag; W=0 disables snapping for paranoid runs.

**Tiering (recommended combination), in priority order into the library:**
1. Annotated transcripts (GENCODE) — base layer, always present.
2. Read-supported chains that match annotation — fold into (1) by dedup.
3. Read-supported chains with ≥1 novel junction (the mutant payload) — support floor default 1,
   tagged `source=round1_novel`, dinucleotides recorded but not gated.

---

## 2. Padding strategy + transcript→genome block map (Designer-2 interface)

### 2.1 Why pad, and how much

Each isoform's terminal exons are extended with **genomic flank** so that reads whose TSS
starts upstream of the annotated 5' exon, or whose 3' end / poly-A region extends past the
annotated TES, still have reference to align into rather than soft-clipping off the end.

- **Default pad: 1–2 kb each side.** NOT 10 kb by default. Rationale below.
- **Adaptive upper bound up to ~10 kb**, but **hard-clipped at the nearest intergenic
  boundary** (nearest neighboring annotated exon on the same strand, minus a 50 bp guard). A
  naive 10 kb pad routinely swallows a neighboring gene's exon in the compact human genome,
  creating spurious "this read aligns great to gene A's pad" hits that are really neighbor-gene
  transcripts. Clipping at intergenic boundaries is mandatory.
- **3' pad connects to the existing poly-A walkback.** RECTIFY's walkback needs genomic
  sequence *past* the CPA to walk through the A-tract; the 3' pad supplies exactly that, so the
  pad length on the 3' side should be ≥ the walkback's max look-ahead window (reuse that
  constant, do not invent a new one).

### 2.2 Behavior of reads whose TSS/CPA fall inside the pad

- A read whose 5' end maps into the 5' pad: fine — it aligns contiguously through pad→exon1.
  At lift-back (Designer 2), the pad portion is plain genomic (no N-op), so the genome
  coordinates are a straight 1:1 block. The TSS is then just the read's leftmost genomic
  coordinate; Round 2 has not biased it.
- A read whose 3' end / poly-A falls in the 3' pad: the aligned 3' end lands in genomic pad;
  poly-A walkback then runs on genomic sequence as usual. The pad must be genomic (real
  sequence), never N's, or walkback breaks.
- A read extending *beyond* the pad: it soft-clips past the pad end. That read is a candidate
  for pad extension on the next build, or simply keeps its Round-1 alignment (Round 2 is
  allowed to lose — see §3/§4).

### 2.3 Block-map format (the interface Designer 2 lifts coordinates from)

One record per library cDNA. Emitted as a sidecar (parquet preferred for scale; JSONL
acceptable). **This is the contract; get orientation right or lift-back is a bug farm.**

```jsonc
{
  "cdna_id": "ENSG00000155657.TTN.chain00042.r1novel",   // unique; encodes gene + chain hash + source
  "gene_id": "ENSG00000155657",
  "strand": "+",                  // GENE/RNA strand. See orientation rule below.
  "cdna_length": 104231,          // total length of the padded cDNA sequence in the FASTA
  "source": "round1_novel",       // annotation | round1_match | round1_novel
  "blocks": [
    // ordered 5'->3' ALONG THE cDNA (transcript) sequence.
    // each block is a contiguous run that maps 1:1 (no indels) to a genome interval.
    {"t_start": 0,      "g_start": 100000, "length": 1500, "is_pad": true,  "kind": "pad5"},
    {"t_start": 1500,   "g_start": 101500, "length": 320,  "is_pad": false, "kind": "exon"},
    {"t_start": 1820,   "g_start": 140200, "length": 410,  "is_pad": false, "kind": "exon"},
    // ... interior exons ...
    {"t_start": 103900, "g_start": 250000, "length": 331,  "is_pad": true,  "kind": "pad3"}
  ],
  "junctions": [
    // EXPLICIT N-op insertion sites for Designer 2. One per exon-exon boundary.
    // Between block[i] and block[i+1], the transcript coordinate is CONTIGUOUS
    // (block[i].t_start+length == block[i+1].t_start) while the genome coordinate JUMPS.
    // THAT JUMP **IS** the intron Designer 2 emits as an N-op of length = g_gap.
    {"after_block_index": 1, "t_pos": 1820, "g_donor": 101820, "g_acceptor": 140200,
     "intron_len": 38380, "source": "round1", "novel": true,
     "donor_dinuc": "GT", "acceptor_dinuc": "AG"}   // recorded, NOT gated
  ]
}
```

**Orientation rule (the minus-strand trap):**
- The cDNA FASTA sequence is always written **5'→3' in RNA/gene orientation**. For a
  minus-strand gene this means the cDNA is the reverse-complement of the plus-strand genome.
- `t_start` increases 5'→3' along that RNA-oriented cDNA.
- `g_start` is always a **plus-strand genome coordinate** (the canonical reference frame). For a
  minus-strand gene, consecutive blocks have **decreasing** `g_start`, and the per-block 1:1 map
  is `cdna[t_start + k]  ==  revcomp(genome[g_start + length - 1 - k])`. Designer 2 must read
  `strand` and apply revcomp; do not assume g_start increases with t_start.
- Recommendation: include a redundant `g_end` per block so Designer 2 never has to recompute
  the genome interval under a strand. (Cheap, removes a whole class of off-by-one bugs.)

**The N-op site is the transcript-contiguous / genome-jump boundary.** Stated explicitly in the
`junctions` array so Designer 2 does not have to re-derive it from the blocks. `intron_len` is
the genome gap; for minus strand it is still a positive genomic distance (g_donor/g_acceptor
given in plus-strand coords with donor>acceptor flagged by strand).

### 2.4 Provenance sidecar (lab convention)

Every library build emits `PROVENANCE.json`: GENCODE version + release date, Round-1 junction
source (run id, commit, BAM paths), collapse method (FLAIR collapse-only vs in-house) + tool
versions, pad parameters (default len, max len, intergenic-clip guard), W snap value, support
floor, build script filename, render date. Per the lab's provenance-alongside-every-output rule.

---

## 3. Round-2 aligners + invocation

Assessed, not just tabulated. The reference is already-spliced cDNA, so the job is
**contiguous, semi-global mapping of noisy ONT reads to a transcript**, letting the read's
5'/3' ends soft-clip into/through the pad. NO splice mode in Round 2 — that is Round 1's job.

| Aligner | Suited? | Invocation (Round-2, cDNA library) | Installed in RECTIFY? |
|---|---|---|---|
| **minimap2** | **YES — primary workhorse.** `map-ont` is *built* for noisy-contiguous-to-transcript; no splice. | `minimap2 -ax map-ont --for-only -N 100 -p 0.8 --secondary=yes -t <n> lib.mmi reads.fq`. `--for-only` (forward-transcript-only, library is RNA-oriented) and high `-N` because paralogous isoforms of the same gene are near-identical. **No `-G`/`--max-intron`** (no splice). | ✅ |
| **mapPacBio (BBMap)** | **YES — secondary workhorse.** Long-read mode, contiguous. Drop the mammalian-intron flags used in Round 1. | `mapPacBio.sh ref=lib.fa in=reads.fq ... maxindel=100 intronlen=0` (no intron-finding; large deletions not expected within a cDNA). `pt:i` QNAME handling per existing RECTIFY convention. | ✅ (compute-node only) |
| **GMAP→cDNA** | **MAYBE — candidate.** GMAP can do cDNA-vs-cDNA, but its value-add (spliced alignment) is moot here; run **`-n 1`, no `--nofails`** (the `--nofails` artifact is documented in ALIGNER_RECOMMENDATIONS — do not repeat it). Evaluate as a tie-breaker only. | `gmap -d lib_db -n 1 -f samse ...` | ✅ (added to panel 2026-06-14) |
| **winnowmap2** | **MAYBE.** Strongest genome aligner for chr5 repeats; in non-splice contiguous mode it behaves like minimap2. Worth trying as a minimap2 alternative if paralog repeats cause minimap2 mis-mapping. | `winnowmap -ax map-ont ...` (no splice) | ✅ (added to panel 2026-06-14) |
| **BLAT** | **NO (default) — relegate.** BLAT is a high-identity (>95%) mRNA/EST tool; raw ONT at ~87–95% identity is outside its design point and it is slow on long noisy reads. Only meaningful on **consensus-corrected** reads (cDNA pipeline UMI consensus), where identity is high. Keep OUT of the default Round-2 panel; offer as opt-in for corrected-consensus inputs. | (opt-in only) `blat lib.2bit reads.fa out.psl -minIdentity=90` | ✅ not standard in RECTIFY |
| **Magic-BLAST** | **NO (default) — relegate.** Designed as an RNA-seq aligner incl. long reads, but no evidence it beats minimap2 on noisy ONT-to-transcript, and BLAST-family is slow at scale. Park as a research comparison, not production. ([Magic-BLAST paper](https://www.researchgate.net/publication/334691704_Magic-BLAST_an_accurate_RNA-seq_aligner_for_long_and_short_reads)) | (research only) | ❌ |

**Net Round-2 panel recommendation:** **minimap2 `map-ont` (primary) + mapPacBio (secondary)**,
with winnowmap2/GMAP as optional tie-breakers. BLAT/Magic-BLAST relegated to corrected-consensus
or research use. This mirrors the genome panel's structure (a small consensus) but with
splice DISABLED everywhere.

**Common flags across the panel:** output SAM with full CIGAR + soft-clips preserved
(`-Y`/equivalent so 5'/3' clips into pad are visible to Designer 2), keep secondaries for
paralog disambiguation, no intron/splice mode anywhere.

**⚠ Orientation flag is PER-PROTOCOL, not blanket (do not force `--for-only` everywhere):**
- **DRS:** read is RNA-sense → `--for-only` is correct (map forward-transcript only).
- **PCR-cDNA:** reads come off BOTH strands of the cDNA. `--for-only` would silently drop the
  antisense-oriented half. Use `--for-only` ONLY if the upstream `ont_cdna` pipeline has already
  oriented reads (it has XS/strand-tag machinery — **verify whether Round-2 cDNA input is
  pre-oriented**). If not pre-oriented, map BOTH orientations (omit `--for-only`). This is a
  per-protocol branch in the invocation builder, not a constant.

---

## 4. Index / scale / cost

### 4.1 Library size and indexing

- **Library size:** annotated GENCODE transcripts (~2.0–2.5×10^5 across coding+lncRNA;
  [GENCODE stats](https://pmc.ncbi.nlm.nih.gov/articles/PMC11024653/)) + observed novel chains
  (data-bounded, typically ≪ annotated because deep coverage collapses to few chains/gene) +
  padding. Padded total sequence ≈ transcriptome length + 2×pad×N_cdnas. With 1–2 kb pads and
  ~3×10^5 cDNAs this is on the order of a few Gb of sequence — comparable to indexing the genome
  once, done once per run/condition, not per read.
- **minimap2 index:** `minimap2 -x map-ont -d lib.mmi lib.fa` once; tiny build, reused across
  all reads. mapPacBio builds its index at run start from `ref=lib.fa`. Index files are a build
  artifact with their own PROVENANCE.
- **Rebuild cadence:** the library is a function of (annotation, Round-1 junctions). Rebuild only
  when Round-1 output changes (new sample/condition), not per chunk. Cache keyed on
  (GENCODE version hash, Round-1 junction-set hash).

### 4.2 Multi-mapping across paralogous / overlapping isoforms

This is the real correctness hazard, and the read-evidence library **mitigates it by
construction**: a collapsed, non-redundant library has far fewer near-duplicate isoforms than a
combinatorially enumerated one, so fewer ambiguous targets. Remaining handling:
- Keep secondaries (`-N 100`) so all near-equal cDNA hits are visible; do **not** trust a single
  primary (the cDNA pipeline's duplicate-primary hazard, documented in minimap2.md, applies).
- Disambiguate by **score margin**: if the best two cDNA hits are within ε of each other, this is
  a Designer-2 scoring decision — Round 2 should hand both candidates to Designer-2 scoring
  rather than silently picking one. Designer 1's job is only to *surface* the ambiguity (emit all
  within-margin hits), not resolve it.
- Same-gene isoforms differing only by an exon-skip will both align well to a read that doesn't
  span the distinguishing region — that is correct and benign; lift-back maps them to the same
  genome span anyway.

### 4.3 Read-subsetting gate (keep Round 2 from dominating runtime)

Round 2 runs on a **subset**, never all reads. Gate = **(A) spans ≥1 junction AND (B) Round-1
alignment is weak.**

- **(A) Multi-exon only.** A single-exon read gains *nothing* from a pre-spliced reference —
  there are no introns to avoid discovering, so Round 2 ≡ Round 1 for it. Exclude all reads whose
  Round-1 alignment has zero N-ops. (Sharp, cheap, eliminates a large fraction of reads.)
- **(B) Weak Round-1, using signals RECTIFY ALREADY computes** (do not invent new metrics):
  - large terminal/internal **soft-clip** length (read body not placed),
  - high **HP-aware edit distance (HP-ED)** for the read's best Round-1 consensus alignment,
  - **failed the junction-local perfect-match anchor gate** (the read's spliced alignment lost
    the anchor test — the exact case Round 2 might rescue),
  - low **MAPQ** / **multi-aligner disagreement** in the Round-1 consensus (aligners disagreed on
    the chain).
  A read passing **any** of these (above per-signal thresholds, exposed as flags) AND condition
  (A) enters Round 2.

- **Default thresholds** are knobs; start conservative (e.g. soft-clip ≥ 20 bp OR anchor-gate
  fail OR consensus disagreement) and tune on a smoke set. The gate is an OR of weakness signals
  so we never *miss* a rescuable read; precision (not over-running) is secondary because Round 2
  is allowed to lose (Designer-2 scoring keeps the Round-1 alignment if Round-2 isn't better).

- **Cost envelope:** with (A)+(B), Round 2 typically touches a small minority of reads. minimap2
  `map-ont` at ~9 s / 10k reads (per the chr5 benchmark in ALIGNER_RECOMMENDATIONS) means even a
  generous subset is cheap relative to the 5-aligner genome round. Chunk Round-2 alignment per
  the lab's mandatory chunking rule when run on a cluster.

---

## 5. Assumptions flagged for the adversary

1. **(§1.5, biggest)** That a collapse step (FLAIR collapse-only, or in-house) can ingest
   Round-1 novel/non-canonical chains **without** re-imposing WT/GT-AG/min-support gating. FLAIR
   issue #371 shows the default path drops non-canonical junctions — the spike test is mandatory.
   If neither route preserves non-canonical chains cleanly, the in-house collapse (§1.4) is the
   fallback and the feature still stands.
2. **(§1.2)** That distinct read-supported intron chains are few enough per gene to keep the
   library ~10^5-scale. True for typical depth; a pathological hyper-noisy sample could inflate
   distinct chains via junction wobble — W-snapping (§1.4 step 3) is the mitigation; verify it
   collapses wobble rather than fragmenting it.
3. **(§2.1)** That 1–2 kb default pad (intergenic-clipped) is enough for TSS/CPA variation in
   practice, and that the 3' pad ≥ walkback look-ahead is the right coupling. If real TSS spread
   exceeds the pad often, adaptive extension is the lever, not a flat 10 kb (which swallows
   neighbors).
4. **(§3)** That minimap2 `map-ont`-no-splice + mapPacBio actually align a *non-trivial* subset
   of reads BETTER than the Round-1 genome consensus — i.e. that the whole feature pays off. This
   must be demonstrated on a smoke set before scaling; if Round 2 never wins, the gate (§4.3) and
   library are wasted work. This is the feature's empirical go/no-go.
5. **(§4.2)** That surfacing within-margin multi-hits to Designer 2 (rather than resolving them)
   is the right division of labor. If Designer 2 expects a single resolved alignment, the
   interface must change.
7. **(§3) Orientation flag `--for-only` is DRS-correct but cDNA-hazardous.** PCR-cDNA reads come
   off both cDNA strands; forcing forward-only drops the antisense half unless `ont_cdna` already
   oriented them. Assumes the invocation builder branches per protocol and the cDNA pre-orientation
   state is verified. If wrong, Round 2 silently loses half of every cDNA sample.
8. **Strand orientation (§2.3)** is correct as specified (RNA-oriented cDNA, plus-strand g_start,
   revcomp mapping for minus genes). A single off-by-one or un-revcomp'd block silently corrupts
   every minus-strand lift-back — Designer 2 should assert the 1:1 identity per block on a test
   cDNA at build time.

---

## 6. Builder starting checklist (concrete first steps)

1. Read `docs/ALIGNER_RECOMMENDATIONS.md` + `methods/aligner_environments.md` for current panel
   + env activation. minimap2/mapPacBio/winnowmap2/gmap all present.
2. Spike test FIRST: take one read with a known non-canonical Round-1 junction; run FLAIR
   collapse-only (`--support 1`, no `--check_splice`, `--no_gtf_end_adjustment`) on RECTIFY
   Round-1 chains; confirm the non-canonical chain survives into the GTF. If it does NOT,
   implement the §1.4 in-house collapse instead.
3. Build the SpliceGraph (§1.3) for one small locus; emit observed_chains.
4. Emit padded cDNA FASTA + block map (§2.3) for that locus; assert per-block 1:1 identity
   (forward and revcomp) against the genome.
5. Index with minimap2 `-x map-ont -d`; align a Round-2-gated read subset; preserve soft-clips.
6. Hand block map + per-read candidate hits (with within-margin multis) to Designer 2.
7. PROVENANCE.json + chunked cluster submission throughout.
```

Sources referenced inline: GENCODE stats (PMC11024653), minimap2 transcriptome mapping (lh3/minimap2 issues #340/#790, janbio blog), FLAIR docs + issue #371, Magic-BLAST paper, long-read isoform-tool assessment (Nat Commun 2024 s41467-024-48117-3).
