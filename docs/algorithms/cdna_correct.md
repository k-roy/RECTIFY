# `rectify correct-cdna` — UMI consensus for ONT PCR-cDNA

`rectify correct-cdna` collapses PCR copies of the same starting RNA molecule
into a single per-molecule consensus read, for libraries that carry a
structured UMI on the cDNA adapter (currently Oxford Nanopore SQK-PCB114.24
is the documented configuration; the algorithm generalizes to any kit with a
fixed-length UMI adjacent to a known anchor motif).

<figure markdown>
  ![PCB114.24 read architecture: 5' adapter — SSP motif — 27-nt structured UMI ((TT-VVVV)×4 + TTT) — cDNA insert — VNP poly(T) — 3' adapter. The same molecule can be sequenced sense or antisense, with the UMI preserved in either orientation.](../figures/cdna_umi_architecture.png){ width="720" }
  <figcaption>Oxford Nanopore SQK-PCB114.24 read architecture. The 27-nt structured UMI sits immediately downstream of the SSP motif and supplies ~43 M molecular barcodes. Sense and antisense reads of the same starting molecule share the same UMI.</figcaption>
</figure>

## Goal

For SQK-PCB114.24 PCR-cDNA libraries, every starting RNA molecule produces:
- N PCR copies (typically 2–10) of the same dsDNA cDNA
- Each PCR copy can be sequenced as either strand of the dsDNA → "sense" or
  "antisense" reads relative to the original mRNA
- Within each orientation, basecalled-3' ends are deterministic (set by RT
  priming or template-switching position); basecalled-5' ends fall off
  stochastically due to nanopore translocation truncation

`rectify correct-cdna` is **Stage 1** of the cDNA pipeline: it collapses reads
sharing `(locus, orientation, UMI)` into one consensus per UMI cluster and
emits a FASTQ (`stage1_consensus.fastq.gz`) for re-alignment, carrying the
cluster metadata as FASTQ comment tags. The downstream stages — re-alignment
(`rectify align`) and gene assignment, isoform clustering, and Type-1↔Type-2
reconciliation (`rectify cdna-analyze`) — operate on post-align coordinates.

> **Note:** earlier internal designs performed a cross-orientation
> sense+antisense pair-overlap merge inside `correct-cdna`. In the shipped
> command this is deferred: orientation pairing and isoform reconciliation now
> live in `rectify cdna-analyze`, which has post-alignment positions precise
> enough for the ±5 bp grouping. The cross-orientation merge sections below
> describe that design and are retained for background; the live command stops
> at the per-cluster Stage-1 FASTQ.

---

## Empirical foundation (validated 2026-05-10)

Established by analysis of `wt_rep1.bam` chrI, 40,788 primary reads:

1. **UMI architecture confirmed.** PCB114.24 has a 27-nt structured UMI
   `(TT-VVVV)×4 + TTT` immediately downstream of the SSP motif
   `TTTCTGTTGGTGCTGATATTGCT`. ~16 informational V positions × 3 bases ≈ 43 M
   diversity. UMI-extractable from 81.7% of reads (the rest have basecaller
   errors in the SSP region or were truncated before reaching it).

2. **Sense and antisense reads of the same starting molecule share the same
   UMI.** Out of 1,769 multi-read clusters at same `(anchor, orientation)`,
   65.5% (1,158) showed both FLAG=0 and FLAG=16 reads with matching UMIs.
   This is a strict-matching lower bound; with edit-distance tolerance the
   true rate is higher.

3. **PCR siblings have IDENTICAL UMIs (Hamming=0).** In size-2 same-anchor
   clusters: 6.5% at d=0, 0% at d=1–5, the rest at d=6+ matching the random-
   pair distribution. Basecaller is remarkably accurate on the structured
   UMI; there's essentially no smear of d=1–3 substitution variants.

4. **Indels are common in homopolymer T-runs (57% of pairs show Lev<Ham).**
   Therefore Levenshtein, not Hamming, is mandatory for accurate UMI
   clustering.

5. **Random-pair Levenshtein distance distribution:**
   - d=0: 0/30000 (0 %)
   - d=3: 0.003 %
   - d=4: 0.04 %
   - d=5: 0.21 %
   - d=6: 0.93 %

6. **Recommended UMI clustering threshold: Levenshtein ≤ 3.** Captures all
   true PCR sibs (d=0 + indel-shifted) with ~0.003 % false-positive rate.
   Lev ≤ 4 is acceptable; Lev ≥ 5 starts admitting random collisions.

---

## CLI

```
rectify correct-cdna INPUT.bam -o OUTDIR \
    [--umi-edit-distance 3]      # max Levenshtein distance within a UMI cluster (default 3)
    [--anchor-window-bp 25]      # alignment-end anchor window for same-locus grouping (default 25)
    [--per-cluster-cap 200]      # drop clusters above N reads (PCR-jackpot / polyA pileups)
    [--gff genes.gff3]           # used for rDNA masking (rRNA_gene loci excluded by default)
    [--reference REF.fa]         # enables walk-back anchor canonicalization + pA length
    [--umi-clustering directional|components]   # default: directional (umi_tools 2× rule)
    [--strand-aware-consensus]   # split reads by is_reverse, build per-strand sub-consensus, merge
    [--no-poa]                   # force pileup-only consensus even if abPOA is available
    [--no-mask-rdna]             # disable rDNA masking
    [--region chrI]              # restrict to one BAM region (testing)
    [--workers N]                # parallel worker processes over BAM regions (default 1)
    [-v]
```

Downstream: run `rectify align` on `stage1_consensus.fastq.gz`, then
`rectify cdna-analyze` on the resulting BAM.

---

## Inputs

- **Aligned BAM** (PCR-cDNA reads, sorted + indexed). Auto-indexes if `.bai` is absent.
- `--reference` FASTA (optional) — enables walk-back anchor canonicalization and
  poly(A) tail-length measurement during UMI extraction.
- `--gff` (optional) — used only for rDNA masking (rRNA_gene loci excluded by
  default to avoid the O(n²) UMI-clustering blow-up on the chrXII tandem repeats).

## Outputs

```
OUTDIR/
└── stage1_consensus.fastq.gz   # one consensus record per UMI cluster
```

Each FASTQ record carries TAB-separated SAM-style comment tags so
`rectify align -y` (`minimap2 -y`) propagates them into the post-align BAM:

| Tag | Meaning | Tag | Meaning |
|-----|---------|-----|---------|
| `XU` | canonical UMI | `XO` | orientation (fwd/rev) |
| `XC` | cluster size | `XR` | input read IDs |
| `XM` | consensus method | `XF` | full-length tier |
| `XA` | pre-align poly(A) length | `XT` | read type (1/2) |
| `XY` | read subtype | `XQ` | 5' pre-trim length |
| `XK` | 3' pre-trim length | `XB` | strand split (n_top/n_bot) |

Gene assignment, isoform clustering, and Type-1↔Type-2 pairing run downstream
in `rectify cdna-analyze` (which emits `clusters.tsv`, a per-molecule
`corrected_reads.tsv`, `isoforms.tsv`, `t1t2_pairs.tsv`, and `consensus_tagged.bam`).

---

## Stage 1 — Per-orientation cluster + consensus

<figure markdown>
  ![UMI consensus pipeline: reads sharing (locus, orientation, canonical UMI) form a cluster; within the cluster, abPOA computes a per-orientation consensus sequence; the consensus is re-aligned to the reference and emitted as a single high-quality read per starting strand.](../figures/cdna_poa_consensus.png){ width="720" }
  <figcaption>Per-orientation UMI consensus pipeline. Reads sharing <code>(anchor, orientation, UMI)</code> within Levenshtein ≤ 3 are clustered; abPOA computes a partial-order-alignment consensus; the consensus is emitted to the Stage-1 FASTQ (<code>stage1_consensus.fastq.gz</code>). Re-alignment to the reference is a separate downstream <code>rectify align</code> step — Stage 1 itself writes no BAM.</figcaption>
</figure>

### 1a. Per-read UMI extraction

For each primary alignment in the input BAM:
1. Get `seq = read.query_sequence` (BAM SEQ, in reference orientation).
2. Search for SSP motif (forward) in seq:
   - if found at position `p`: `umi_basecalled = seq[p+len(SSP) : p+len(SSP)+27]`
   - mark `seq_has_ssp_at_5prime_of_bam_seq = True`
3. Else search for SSP_RC motif:
   - if found at position `p`: `umi_basecalled = revcomp(seq[p-27:p])`
   - mark `seq_has_ssp_at_5prime_of_bam_seq = False`
4. Else: skip (this read does not contribute; ~18 % of reads).

The `umi_basecalled` is in the original basecalled orientation regardless of
which BAM orientation the read was stored in.

### 1b. Anchor + orient labels

`orient` indicator for clustering: based on whether SSP appeared as forward or
RC in BAM SEQ. This separates reads that came from one strand of the dsDNA from
reads that came from the other strand. (Independent of FLAG.)

`anchor`: basecalled-3' genomic position. The basecalled-3' is the end opposite
the SSP/UMI side, i.e., the polyT/VNP side. In genomic coordinates:
- If SSP was found at forward (5' of BAM seq): basecalled-3' is at the OTHER
  end of the read sequence. Whether that maps to genomic high or low coord
  depends on FLAG: `aln_end - 1` if FLAG=0, `pos - 1` if FLAG=16.
- If SSP was found as RC (3' of BAM seq): basecalled-3' is at BAM SEQ START,
  which is `pos - 1` (FLAG=0) or `aln_end - 1` (FLAG=16).

(Implementation will table-ify this; the rule is "basecalled-3' is the end
of the read OPPOSITE the SSP, mapped through the FLAG.")

### 1c. Cluster initialization

Group reads by `(chrom, anchor // anchor_window * anchor_window, orient)`. This
gives "candidate molecule clusters" — reads that *could* be from the same
starting molecule, gated by position and orientation but not yet by UMI.

Drop clusters with > `--per-cluster-cap` reads (these are polyA-pileup hot
spots where many distinct mRNAs share the same anchor; UMI clustering inside
these gives correct results but is expensive and the per-molecule meaning is
preserved by the UMI groups).

### 1d. UMI subclustering within each anchor cluster

Within each anchor cluster, build a Levenshtein graph: two UMIs share an edge
if their Levenshtein distance ≤ `--umi-edit-distance` (default 3). Find
connected components. Each component is a "molecule cluster" — the reads
assigned to it are the PCR siblings of one starting strand.

Use a `umi_tools`-style **directional** algorithm to handle UMI-error chains
correctly: a UMI A is considered the same molecule as UMI B if there's a path
A→B through intermediate UMIs each at d≤threshold. Implemented with custom
connected-components clustering (`umi_components` / `umi_components_directional`,
the latter a umi_tools-style directional method) over rapidfuzz Levenshtein
edit-distances — no networkx dependency.

For each molecule cluster:
- If 1 read: the read passes through unchanged (singleton).
- If ≥2 reads: build a POA consensus.

### 1e. Per-cluster consensus

Use **abPOA** (`pip install pyabpoa`) to build a multiple-sequence-alignment
consensus across the reads in the cluster. Operate on the basecalled-3'
portion of each read (the part downstream of SSP+UMI+GGG, i.e., the actual
cDNA insert). Strip SSP/UMI/GGG/barcode/adapter regions before consensus.

Per-base quality for the consensus: use depth of agreement (e.g., 30 + 5×
depth for fully-agreeing positions, capped; lower for tie-positions).

### 1f. Emit consensus FASTQ

Write each per-cluster consensus to `stage1_consensus.fastq.gz`, with the
cluster metadata as FASTQ comment tags. Re-alignment is a separate step:
`rectify align` runs the multi-aligner panel on this FASTQ and propagates the
comment tags into the BAM via `minimap2 -y`.

---

## Cross-orientation pair-overlap merge (deferred to `cdna-analyze`)

> The pairing and merge described in this section is **not performed by the
> shipped `correct-cdna`**; orientation pairing and Type-1↔Type-2 reconciliation
> now run in `rectify cdna-analyze` on post-align coordinates. The design is
> retained here for background.

<figure markdown>
  ![Isoform clustering and T1↔T2 reconciliation schematic: Stage 1 emits one consensus per (locus, orientation, UMI); Stage 2 pairs sense+antisense consensuses sharing a UMI within a span filter and computes a position-overlap merge that recovers the full molecule from its two complementary basecalls.](../figures/cdna_isoform_clustering.png){ width="720" }
  <figcaption>Stage 2 pairs Stage 1 consensuses that share a UMI across orientation and whose anchors fall within <code>--max-cross-orient-span</code> bp of each other. The pair-overlap merge recovers the full molecule, with each orientation contributing the half it sequenced most accurately.</figcaption>
</figure>

For each pair of Stage 1 consensuses (A: orient=fwd, B: orient=rev) that share
the same UMI canonical (Lev ≤ threshold across orient) **AND whose anchors are
within `--max-cross-orient-span` bp of each other** (default 3000 for yeast;
tune for organism mRNA size): they are putatively from the same starting RNA
molecule.

**The span filter is mandatory.** Empirical chrI analysis (2026-05-10) found
1086 cross-orient UMI matches at Lev≤3 across the chromosome, of which only
51 spanned ≤ 3 kb. The other 1035 were random UMI collisions (consistent with
the ~330 expected from the random-pair Lev≤3 baseline of 0.003 % times the
~11 M cross-pair candidate space). Without the span filter, Stage 2 would
mistakenly merge thousands of unrelated molecules.

The two consensuses cover **opposite ends** of the molecule:
- The `orient=fwd` consensus has reliable basecalled-3' at the polyT/VNP side
  → its mapping anchor is at one end of the original mRNA span.
- The `orient=rev` consensus has reliable basecalled-3' at the SSP/cap side
  → its mapping anchor is at the OTHER end of the mRNA span.
- Their alignments should overlap in the middle of the mRNA span.

Merge by alignment-aware pair-overlap consensus:
1. Re-align both consensuses to the reference (already done at end of Stage 1).
2. Identify the genomic overlap window (where both alignments cover the same
   reference positions).
3. In the overlap: take per-position majority base, breaking ties by quality;
   in the non-overlap regions: keep each consensus's own base.
4. Emit one combined record per molecule, with:
   - `XU:Z:<canonical_UMI>` tag
   - `XS:i:<n_reads_sense>` and `XA:i:<n_reads_antisense>` tags
   - Quality scores reflecting consensus depth (sum over both orientations
     in overlap region, single orientation otherwise)

If only one orientation of a molecule was sequenced (singleton in Stage 2
pairing), the Stage 1 consensus passes through unchanged.

---

## Algorithm complexity

Asymptotic cost of Stage 1 (`correct-cdna`):
- Stage 1a (UMI extraction): O(N) reads (motif find per read).
- Stage 1c (anchor cluster): O(N).
- Stage 1d (UMI subcluster): per anchor cluster O(K²) Levenshtein +
  connected/directional components, with K = cluster size (mostly 1–10; rare
  poly(A) hot-spots are capped by `--per-cluster-cap`). O(N) overall at K = O(1)
  with rapidfuzz Levenshtein.
- Stage 1e (POA consensus): per-cluster abPOA, only for clusters with ≥ 2 reads.

Re-alignment (`rectify align`) and orientation pairing (`cdna-analyze`) are
separate downstream steps with their own cost. (Indicative wall-clock figures
depend heavily on sample size, UMI-extractable fraction, and `--workers`; run a
small `--region` first to gauge throughput rather than relying on a single
quoted number.)

---

## Dependencies

- `pysam` (already in the rectify env)
- `rapidfuzz` — Levenshtein distance for UMI clustering
- `pyabpoa` — abPOA partial-order-alignment consensus (optional; `--no-poa`
  falls back to a pileup-only consensus when it is unavailable)

Re-alignment of the Stage-1 FASTQ is handled by `rectify align` (the standard
multi-aligner panel), not by `correct-cdna` itself.

---

## Design notes

- **Cross-orientation merge uses position overlap, not sequence overlap.**
  Sequence overlap would be more forgiving in mid-mRNA regions where the two
  consensuses don't align cleanly, but position overlap (using alignment
  coordinates) is simpler and was empirically sufficient on the validation set.
- **Reads without an extractable UMI are dropped, with a counter logged.**
  Position-only dedup would admit inter-molecule collisions and is not the goal
  of this subcommand. (The ~18% drop rate cited above is from the 2026-05-10
  chrI validation and will vary by sample/basecaller.)
- **Orientation-ambiguous reads** (SSP motif not found cleanly) are dropped in
  Stage 1a and reported in stats.
- **Exon-aware clustering is intentionally not performed.** Reads that start
  at the same exon boundary but follow different splice forms should be
  treated as distinct starting molecules — which is what the
  `(locus, orientation, UMI)` cluster key already produces. The ±25 bp
  anchor window absorbs minor priming jitter without merging splice variants.
- **Sample multiplexing is out of scope.** The input BAM is assumed to be
  pre-demultiplexed per sample.
