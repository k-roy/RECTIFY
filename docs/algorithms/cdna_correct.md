# `rectify correct-cdna` — Design Spec

**Status:** draft v1, 2026-05-10
**Author:** Kevin Roy + Claude
**Scope:** new subcommand for the rectify pipeline that merges PCR-cDNA reads
originating from the same starting RNA molecule into a single per-molecule
consensus read.

---

## Goal

For SQK-PCB114.24 PCR-cDNA libraries, every starting RNA molecule produces:
- N PCR copies (typically 2–10) of the same dsDNA cDNA
- Each PCR copy can be sequenced as either strand of the dsDNA → "sense" or
  "antisense" reads relative to the original mRNA
- Within each orientation, basecalled-3' ends are deterministic (set by RT
  priming or template-switching position); basecalled-5' ends fall off
  stochastically due to nanopore translocation truncation

`rectify correct-cdna` collapses this constellation of reads back into one
high-quality, full-length consensus read per starting RNA molecule, by:

1. Per-orientation clustering of reads sharing `(locus, orientation, UMI)`
   → one consensus per starting strand
2. Cross-orientation pair-overlap merge of matching sense+antisense
   consensuses → one final consensus per starting RNA molecule

The result is a BAM with one record per starting molecule, where each record
is a high-quality merged read covering the full mRNA span (or as much of it
as the input reads collectively spanned).

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
rectify correct-cdna INPUT.bam --genome REF.fa -o OUTDIR \
    [--umi-edit-distance 3]              # default 3 (Lev), validated empirically
    [--anchor-window-bp 25]              # default ±25 bp on basecalled-3'
    [--ssp-motif TTTCTGTTGGTGCTGATATTGCT]  # PCB114-specific, override for other kits
    [--umi-length 27]                    # PCB114-specific
    [--min-cluster-size 1]               # 1 = also emit singletons; 2 = require ≥2 reads
    [--cross-orient-merge / --no-cross-orient-merge]  # Stage 2 toggle (default on)
    [--threads 8]
    [--per-cluster-cap 200]              # safety: drop polyA-pileup loci with >N reads
                                          # (these are multi-molecule, not single-mol PCR sibs)
```

---

## Inputs

- **Aligned BAM** (any short-read aligner; minimap2 / bwa / bbmap for cDNA)
- **Reference genome** FASTA (used for re-aligning the merged consensus)
- (no annotation file required — purely sequence + position based)

## Outputs

```
OUTDIR/
├── cdna_consensus.bam        # one record per starting molecule (after Stage 2 merge)
├── cdna_consensus.bam.bai
├── cdna_clusters.tsv         # cluster manifest:
│   #   cluster_id  chrom  anchor  orient  umi_canonical  n_reads  consensus_read_id
├── cdna_per_orient_consensus.bam  # intermediate: Stage 1 output (one per (locus,orient,UMI) cluster)
├── cdna_stats.tsv            # per-sample QC (see below)
└── cdna_correct.log
```

`cdna_stats.tsv` fields:
- input_reads_total, input_reads_with_umi
- distinct_clusters, mean_cluster_size, median, p95
- pct_reads_in_size1_cluster (= dedup ratio inverse)
- pct_clusters_paired_cross_orient (Stage 2 success rate)
- mean_consensus_length, mean_consensus_qual

---

## Stage 1 — Per-orientation cluster + consensus

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
A→B through intermediate UMIs each at d≤threshold. Implementation via the
networkx `connected_components` algorithm or rapidfuzz-based pairwise
edit-distance.

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

### 1f. Re-align consensus

Re-align each per-cluster consensus to the reference (minimap2 `-ax splice
-G 3000` for spliced cDNA) and emit as Stage 1 output BAM.

---

## Stage 2 — Cross-orientation pair-overlap merge

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

For a yeast-scale BAM (~1.5 M reads) with ~50% UMI-extractable:
- Stage 1a (UMI extraction): O(N) reads, ~1 μs each (motif find) → ~1 s
- Stage 1c (anchor cluster): O(N), ~1 s
- Stage 1d (UMI subcluster): per anchor cluster O(K²) Lev distance + connected
  components, with K = avg cluster size (mostly 1–10, rare polyA hot-spots
  capped). Total O(N) at K=O(1), ~10 s with rapidfuzz Levenshtein.
- Stage 1e (POA consensus): per cluster ~5 ms × ~30 % of N clusters → ~30 s.
- Stage 1f (re-alignment): minimap2 streaming, ~30 s.
- Stage 2: O(M log M) where M = Stage 1 output count, ~10 s.

Total expected: **~90 s for a full sample of 1.5 M reads** on M1 with 8 cores.

---

## Dependencies

- `pysam` (already in rectify env)
- `rapidfuzz` (Levenshtein) — `pip install rapidfuzz`
- `pyabpoa` (POA consensus) — `pip install pyabpoa`
- `mappy` (in-process minimap2) — `pip install mappy`
- (no networkx — implement connected components in pure Python; simpler)

---

## Open design questions

1. **Should the cross-orient merge use sequence overlap or position overlap?**
   Position overlap is simpler (use alignment coordinates); sequence overlap
   handles cases where the consensuses don't align cleanly in the middle. Plan:
   start with position overlap; revisit if mid-mRNA accuracy is poor.

2. **What about reads that don't have an extractable UMI?** ~18 % of reads.
   Options: (a) drop them; (b) fall through to position-only dedup with a
   "no-UMI" cluster label. **Default: drop, but log the count.** Position-only
   dedup admits inter-molecule collisions and isn't the design's purpose.

3. **What about reads that are orientation-ambiguous (SSP not found cleanly)?**
   These are dropped in 1a. Flagged in stats.

4. **Do we cluster by exon-aware position?** Reads that start at the same exon
   but at slightly different splice forms would have similar but not identical
   anchors. The 25 bp window handles small differences; for splice form
   diversity at the same gene, those should cluster as separate molecules
   (different starting RNAs), which is what we want.

5. **How do we handle barcodes (sample multiplexing)?** Out of scope for this
   subcommand — assume the input BAM is already demultiplexed per sample.

---

## Validation checklist (before merging to rectify master)

- [ ] On `wt_rep1.chrI.bam`: emits cdna_consensus.bam with N reads, where N is
      between (size-1 cluster count) and (input read count) — gives the
      empirical dedup ratio.
- [ ] On the same data: cdna_stats.tsv pct_clusters_paired_cross_orient is in
      the 30–70 % range (mostly singletons due to undersampling, but a real
      fraction of cross-orient pairs as observed in the empirical analysis).
- [ ] At least 2 spot-checked clusters where Stage 2 paired sense+antisense
      consensuses give a longer / cleaner alignment than either input.
- [ ] Test suite: unit tests for UMI extraction, anchor computation,
      Levenshtein clustering, POA consensus.
- [ ] Cross-cluster integration test: synthetic "5 PCR sibs of 1 molecule, 3
      sense + 2 antisense" → must collapse to 1 record.

---

## Implementation roadmap

- **v0 (this session):** Standalone Python module under `~/work/cdna_dev/src/`,
  testable on `wt_rep1.chrI.bam`. Stages 1a + 1b + 1c + 1d only (UMI
  extraction + clustering, no consensus). Validates the cluster output by
  printing per-cluster stats.
- **v1:** Add Stage 1e (abPOA consensus) and 1f (minimap2 re-align). Output
  per-orientation consensus BAM.
- **v2:** Add Stage 2 (cross-orient merge). Output final cdna_consensus.bam.
- **v3:** Integrate into rectify as `rectify correct-cdna` subcommand;
  follow rectify's CLI / logging / provenance conventions.
- **v4:** Production run on the full 6 PCB114.24 samples on H2; deposit to
  `/u/project/guillom/shared/processed/by4742_cdna_pcb114_2026_consensus/`.
