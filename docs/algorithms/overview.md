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
│   N-op refinement (Module 2H)                                 │
│   ├─ Junction re-scoring (HP-aware split-alignment)           │
│   └─ [--junction-penalty-table] empirical HP costs (optional) │
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

## Empirical HP penalty calibration

Module 2H (N-op junction refinement) scores candidates using per-HP-length deletion and insertion
costs. The defaults are heuristic constants derived from known Nanopore error behaviour. For
production runs, these can be replaced with empirically derived values using
`empirical_cigar_error_profiler.py`.

### How the profiler works

1. **Multi-aligner agreement** — for each read present in all aligner BAMs, find genomic intervals
   that every aligner classifies as exonic (non-N-op). Within these agreed exonic regions, collect
   per-position CIGAR ops from all aligners.

2. **Consensus op calling** — a position is emitted only if all sequence-bearing aligners (those
   with `query_sequence`) agree on the same op type (M/D/X/I). Aligners without a stored query
   sequence (e.g. gapmm2, which uses explicit `=`/`X` CIGAR ops) contribute to D/I/X counts but
   abstain from M calls, which require per-base sequence comparison.

3. **HP context lookup** — for each agreed op, the reference homopolymer run length at that
   position is recorded using the same `_hp_run_length()` used by Module 2H.

4. **Rate and penalty computation** — error rates per (op, HP length) are computed across all
   observed positions; the penalty table normalises so that `sub(HP=1) = 1.0` (reference), then
   `penalty(op, hp) = rate(sub, HP=1) / rate(op, hp)`. Higher error rate → lower penalty
   (the operation is more expected and should be discounted in scoring).

### Key empirical finding (wt_by4742_rep1, 100K reads)

| HP length | del rate | empirical penalty | heuristic penalty |
|-----------|----------|-------------------|-------------------|
| 1         | ~2.4%    | 1.00              | 1.00              |
| 4         | ~5.8%    | 0.41              | 0.50              |
| 8         | ~19.7%   | 0.12              | 0.25              |

The current heuristic step-function (`del_hp=0.5` for HP≥4) is systematically ~2× too strict:
long homopolymer deletions are even cheaper than assumed. The empirical table corrects this
per-HP-length with sub-1% CV across dataset chunks, indicating stable estimates from ~50K reads.

### Integration with Module 2H

`HpPenaltyTable` (in `junction_refiner.py`) wraps the `penalty_scores.tsv` output. It is
threaded through `_score_hp_anchored` and `_hp_edit_distance` as an optional argument; when
`None`, existing heuristic behaviour is unchanged. Load via `--junction-penalty-table PATH` in
`rectify correct`.

---

## Key design decisions

### Why not just trim poly(A)?

Simple poly(A) trimming fails when poly(A) tails land on genomic A-tracts: the boundary between tail A's and genomic A's is ambiguous, and aligners introduce indels to maximize alignment score, shifting the apparent 3' end downstream of the true CPA site. RECTIFY's walk-back algorithm handles this by comparing read and genome sequences at each position.

### Why multi-aligner consensus?

No single aligner handles all read types optimally. minimap2 is fast but aggressive about placing junctions; mapPacBio forces mismatches at junction boundaries; gapmm2 handles large indels better. Running all three and selecting per-read avoids systematic biases.

### Why cluster-level DESeq2?

Gene-level expression counts can stay constant while CPA site usage shifts proximal/distal. Cluster-level analysis catches these isoform-specific changes that a simple sum would miss.
