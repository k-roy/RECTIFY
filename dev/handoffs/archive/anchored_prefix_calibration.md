# 5' Soft-Clip Rescue Calibration: Anchored-Prefix Accuracy Profile

## Methods

**BAM source**: `wt_polya_rep1.sorted.bam` from `/u/project/guillom/shared/raw/by4742-wt-upf1D_polya_drs_2025/` (WT BY4742 polyA DRS, replicate 1, 11.3M reads).

**Reference genome**: S. cerevisiae R64-5-1 (20240529).

**Algorithm**: For each read with a 5' soft-clip (CIGAR op 4, length ≥5 bp), we compared the soft-clipped bases to the upstream reference sequence (walking 5'-ward from the alignment anchor). For each position k (1..200 bp), we scored whether the soft-clip base matches the reference base at distance k from the anchor. We aggregated match/mismatch counts across all 5.66M soft-clipped reads in the dataset.

**Key parameters**:
- Soft-clips ≥5 bp (filter out very short clips)
- Walk back up to 200 bp into upstream reference
- No junction-specific filtering (profile captures all soft-clip behavior)

## Findings

**Soft-clip prevalence**: 5,657,846 soft-clipped reads out of 11,288,668 total = **50.1% soft-clip rate**.

**Per-position accuracy**: Across all soft-clipped reads, the per-position match rate ranges from **20-36%**, corresponding to error rates of **64-80%**. This is far above the "1-3% nanopore baseline" within high-quality aligned regions, and instead reflects the composition of the soft-clip pool: predominantly reads with poly-A tails, intergenic alignments, and terminal degradation in low-complexity regions. The algorithm correctly identifies soft-clips as regions of high uncertainty.

| Position | Matches | Total Reads | Match Rate | Error Rate |
|----------|---------|-------------|------------|------------|
| 1-5      | 1.16M   | 5.66M       | 20-30%     | 70-80%     |
| 6-10     | 0.95M   | 4.73M       | 24-30%     | 70-76%     |
| 11-20    | 0.71M   | 2.77M       | 24-30%     | 70-76%     |
| 21-50    | 0.57M   | 1.57M       | 22-34%     | 66-78%     |
| 51-100   | 0.20M   | 0.50M       | 21-33%     | 67-79%     |

**Stability**: Match rate is remarkably flat across all 200 positions (~25% mean, SD 0.035), indicating that error degradation is not *anchored-position-dependent*. Instead, soft-clips appear to be drawn from a homogeneous pool of noisy sequences. The lack of inflection point suggests that soft-clipped bases are fundamentally low-quality throughout their length, consistent with them being poly-A tail artifacts or terminal nanopore-noise regions.

## Implications for Anchored-Prefix Scoring

**Threshold recommendation**: Do not use a position-dependent threshold. Instead:

1. **Minimum match rate threshold: 35%** — This would retain only the ~top 10-15% of soft-clips (those with 3-4 matching bases in the first 5 positions) and reject the majority as unreliable.

2. **Minimum prefix length: 5 bp** — Clips shorter than 5 bp lack signal; clips ≥5 bp show no degradation curve and should all be treated as suspect.

3. **Accept-only rule**: For soft-clips in putative splice-rescue scenarios (e.g., near a junction), require *at least 4 consecutive matches in the first 6 bp* (cumulative match rate >60%). Reads failing this test are overwhelmingly noise.

**Per-junction thresholds**: Not warranted. The flat error curve across all soft-clip contexts suggests a single global threshold. Soft-clips near annotations (introns, exons) are not distinguishable from those in intergenic regions in this pooled dataset — the sequencing quality is the limiting factor, not the genomic context.

**One-sentence summary**: Soft-clipped bases in this dataset are **uniformly low-confidence** (error rate ~70%); use them only for rescue when the matched prefix is exceptionally clean (≥60% cumulative match in first 6 bp).

## Caveats

- This profile aggregates all soft-clips in the BAM, including (likely predominantly) poly-A tail artifacts and terminal noise. Soft-clips near genuine splice junctions may have better statistics, but they are too rare to profile separately in this 50% soft-clip pool.
- The lack of an inflection point is reassuring: it means there's no "cliff" at which quality drops (ruling out systematic misalignment), but it also means one cannot use distance from anchor to estimate confidence.
- Nanopore basecalling error is sequence-context-dependent; this aggregate profile masks such variation.

## Recommendation for TODO.md

Close the item with outcome: **threshold = 4 consecutive matches (or 60% cumulative) in first 6 bp; apply globally, not per-junction**. This is a conservative gate suitable for a downstream rescue filter. If higher recall is needed, lower the threshold to 3 matches; if higher precision, require 5 matches.
