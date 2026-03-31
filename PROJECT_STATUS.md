# RECTIFY — Project Status

**Version:** 2.6.0
**Last updated:** 2026-03-31

---

## What rectify does

RECTIFY is a framework for correcting systematic 3' end positioning errors in
nanopore direct RNA-seq data. Errors arise from poly(A) tail homopolymer
ambiguity, AG mispriming, and sequencing indels near cleavage/polyadenylation
sites (CPA). The corrected positions enable accurate CPA site mapping and
differential CPA usage analysis between conditions.

---

## Current state: production-ready

### Core correction pipeline — stable

| Module | Status | Notes |
|---|---|---|
| A-tract ambiguity correction | Stable | NET-seq-guided repositioning within homopolymer runs |
| Poly(A) trimming | Stable | Default on in `run-all` |
| Indel artifact correction | Stable | |
| AG mispriming detection | Stable | |
| Variant-aware rescue | Stable | |
| NET-seq refinement | Stable | Bundled S. cerevisiae data (Churchman 2011, Harlen 2016) |
| Spike-in filter | Stable | Fixed 2026-03-29 (false positive bug) |

### Multi-aligner consensus — stable

Three aligners run in parallel (minimap2, mapPacBio/BBMap, gapmm2) and vote
on the best alignment per read. Splice junction rescue via `--junc-bed` is
wired into minimap2 alignment.

Known gap: 5' soft-clip rescue is length-based, not sequence-based. Can cause
intron-spanning reads to lose to shorter-mapping reads. See TODO.md.

### Streaming memory pipeline — new (2026-03-31)

**Background:** The previous multi-sample analysis path combined all per-sample
corrected TSVs into one file before analyzing. For 21 samples / 150M reads
this required 64+ GB RAM and OOM-killed repeatedly on Sherlock.

**What was built:**

- **Tier 1** (`analyze_command.py`): Drop unused columns on load; lower
  chunked-loading threshold 5 GB → 500 MB.

- **Tier 2** (`analyze_command.py`, `run_command.py`): `rectify analyze --manifest`
  two-pass streaming mode. Pass 1 aggregates positions for clustering; Pass 2
  accumulates per-cluster counts in a dict. Peak RAM = O(clusters × samples)
  regardless of read depth. Validated: 21 samples / 150M reads on 16 GB node.

- **Tier 3** (`bam_processor.py`): `rectify correct` now writes
  `corrected_3ends_index.bed.gz` alongside every `corrected_3ends.tsv`. This
  is a pre-aggregated (chrom, pos, strand, count) file ~300× smaller than the
  full per-read TSV. Manifest mode uses it automatically when present, making
  both passes near-instant.

See `CLAUDE.md` for usage patterns and index back-fill recipe.

### Analysis outputs — stable

| Output | Status |
|---|---|
| CPA clusters (`cpa_clusters.tsv`) | Stable |
| Count matrix (`cluster_counts.tsv`) | Stable |
| Bedgraphs | Stable (single-file mode); skipped in manifest mode — see TODO |
| PCA plot | Stable |
| Sample heatmap | Stable |
| DESeq2 gene + cluster results | Stable (requires ≥2 conditions) |
| GO enrichment | Stable (requires `--go-annotations`) |
| De novo motif discovery | Stable (requires `--genome --run-motif`) |
| Cluster shift analysis | Stable |
| HTML report | Stable |

---

## Active experiment: roadblocks project (Sherlock)

21-sample nanopore run (PRJNA1229592 + in-house BY4742):
- **Per-sample correction:** All 21 samples complete (`minimap2_run/` on Oak)
- **Combined analysis:** Job 20063057 queued (2 tasks; see below)
- **Position indices:** Being generated for remaining samples (background)

### Analysis split (two strain backgrounds)

| Task | Samples | Reference | Output dir |
|---|---|---|---|
| 1 | rna15_rep1-3, wt_rep1-3, ysh1_rep1-3 | WT | `combined_cpa/` |
| 2 | wt_by4742_rep1-3, dst1d_rep1-3, ski7d_rep1-3, 4nqo_rep1-3 | WT_BY4742 | `combined_by4742/` |

Script: `minimap2_run/run_combined_analysis.sh` (uses `--manifest`, 16 GB RAM)

---

## Branches

| Branch | Purpose | State |
|---|---|---|
| `master` | Main development | 2 commits ahead of origin; other in-progress work unstaged |
| `beta` | Parallel development track | Has streaming pipeline + CLAUDE.md (cherry-picked 2026-03-31) |

---

## Key file locations (Sherlock)

| Resource | Path |
|---|---|
| rectify source | `/oak/.../software/rectify/` |
| roadblocks run dir | `/oak/.../projects/roadblocks/rectify_alpha_output/minimap2_run/` |
| manifest | `minimap2_run/manifest.tsv` |
| combined analysis script | `minimap2_run/run_combined_analysis.sh` |
| sherlock SLURM profile | `rectify/slurm_profiles/sherlock_larsms.yaml` |
