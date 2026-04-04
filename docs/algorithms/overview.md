# Algorithm Overview

RECTIFY processes long-read RNA-seq data in three stages: **alignment**, **correction**, and **analysis**.

---

## Pipeline diagram

```
Input: FASTQ or BAM
        │
        ▼
┌───────────────────────────────────────────────────────────────┐
│  Stage 1: Alignment  (skipped if BAM input)                   │
│                                                               │
│   minimap2 ──┐                                                │
│   mapPacBio ─┼─→ Consensus selection by junction scoring ──→ BAM
│   gapmm2  ───┘                                                │
│   [+ uLTRA, deSALT if enabled]                               │
└───────────────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────────────┐
│  Stage 2: Correction (per read)                               │
│                                                               │
│   5' end                                                      │
│   ├─ Junction rescue (soft-clip → upstream exon matching)     │
│   └─ 3'SS truncation rescue                                   │
│                                                               │
│   3' end                                                      │
│   ├─ Poly(A) trimming (--polya-sequenced)                     │
│   ├─ False junction removal                                   │
│   ├─ Indel artifact correction (walk-back algorithm)          │
│   ├─ A-tract ambiguity detection                              │
│   ├─ AG-mispriming screening                                  │
│   └─ NET-seq refinement (NNLS deconvolution, optional)        │
│                                                               │
│   Output: corrected_3ends.tsv + corrected_3ends_index.bed.gz  │
└───────────────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────────────┐
│  Stage 3: Analysis  (multi-sample only)                       │
│                                                               │
│   ├─ Adaptive CPA clustering (valley-based)                   │
│   ├─ Count matrix building (streaming, O(clusters × samples)) │
│   ├─ DESeq2 (gene-level + cluster-level)                      │
│   ├─ APA shift analysis (Jensen-Shannon divergence)           │
│   ├─ GO enrichment (hypergeometric + BH FDR)                  │
│   ├─ Motif discovery (STREME)                                 │
│   └─ HTML report                                              │
└───────────────────────────────────────────────────────────────┘
```

---

## Algorithm pages

| Page | Description |
|------|-------------|
| [3' End Indel Correction](3prime_indel_correction.md) | Walk-back algorithm for poly(A)/A-tract artifacts |
| [Soft-Clip Rescue](softclip_rescue.md) | 5' junction rescue and 3' homopolymer rescue |
| [Multi-Aligner Consensus](multi_aligner_consensus.md) | Scoring and selection across up to 5 aligners |
| [False Junction Handling](false_junction_handling.md) | Removing spurious junctions from poly(A) artifacts |
| [NET-seq Refinement](netseq_refinement.md) | NNLS deconvolution of oligo(A)-spreading |
| [Adaptive Clustering](adaptive_clustering.md) | Valley-based CPA site grouping |
| [APA Detection](apa_detection.md) | Alternative polyadenylation isoform quantification |

---

## Key design decisions

### Why not just trim poly(A)?

Simple poly(A) trimming fails when poly(A) tails land on genomic A-tracts: the boundary between tail A's and genomic A's is ambiguous, and aligners introduce indels to maximize alignment score, shifting the apparent 3' end downstream of the true CPA site. RECTIFY's walk-back algorithm handles this by comparing read and genome sequences at each position.

### Why multi-aligner consensus?

No single aligner handles all read types optimally. minimap2 is fast but aggressive about placing junctions; mapPacBio forces mismatches at junction boundaries; gapmm2 handles large indels better. Running all three and selecting per-read avoids systematic biases.

### Why cluster-level DESeq2?

Gene-level expression counts can stay constant while CPA site usage shifts proximal/distal. Cluster-level analysis catches these isoform-specific changes that a simple sum would miss.
