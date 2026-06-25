# Adaptive Clustering

After per-read correction, RECTIFY groups nearby CPA positions into **clusters** using a valley-based adaptive algorithm. Clusters represent distinct CPA sites used by a gene.

**Implementation:** `rectify/core/analyze/clustering.py`

<figure markdown>
  ![Adaptive clustering schematic: histogram of corrected 3' ends along a gene; a fixed-distance cluster naively merges two real peaks, while the valley-based adaptive cluster splits them at the local minimum between peaks.](../figures/adaptive_clustering.png){ width="600" }
  <figcaption>Fixed-distance clustering (top) merges nearby CPA peaks that should be reported as distinct isoforms. The adaptive valley-based algorithm (bottom) splits peaks at local minima, preserving the multi-site structure that gene-level counts would miss.</figcaption>
</figure>

---

## The problem

A single gene typically has multiple CPA sites within a few hundred bp of each other. Some sites are major (high read count), others are minor. Grouping these into biologically meaningful clusters requires:

1. Not merging two distinct CPA sites used by different isoforms
2. Not splitting a single site into multiple clusters due to position scatter

A fixed-distance merge (e.g., "merge any two positions within 25 bp") handles case 2 but can incorrectly merge case 1.

---

## Valley-based algorithm

The adaptive algorithm treats CPA usage as a 1D density function and finds natural boundaries:

**Step 1: Find peaks**

Scan the position histogram for local maxima (positions with more reads than both neighbors).

**Step 2: Find valleys**

Between each pair of adjacent peaks, find the local minimum (the "valley"). The valley position is the cluster boundary.

**Step 3: Assign boundaries**

Each cluster extends from valley to valley. The cluster's representative position (`modal_position`) is the median position of the cluster's reads (an unweighted median — despite the historical "modal" name, it is not the most-used position).

```
Position:  ...100  101  102  103  104  105  106  107  108  109  110...
Read count:   0    2    15   8    2    1    3    19   7    2    0

Peak at 102 (count=15), Peak at 107 (count=19)
Valley at 105 (count=1) → cluster boundary here

Cluster 1: positions 100-105, peak=102
Cluster 2: positions 105-110, peak=107
```

**Implementation:**

```python
def cluster_cpa_sites_adaptive(
    positions_df,             # DataFrame of CPA positions (chrom, strand, corrected_position, …)
    max_cluster_radius=10,    # Max distance (bp) from peak to cluster boundary
    min_peak_separation=5,    # Min distance (bp) between distinct peaks
    min_reads=5,              # Minimum total reads to keep a cluster
    chrom_col='chrom',
    strand_col='strand',
    position_col='corrected_position',
    count_col=None,
):
    """
    Valley-based CPA clustering.

    Returns: DataFrame with columns cluster_id, chrom, strand, start, end,
    modal_position, cluster_com, n_positions, n_reads, cluster_width
    """
```

---

## Count matrix building

After clustering, a count matrix (clusters × samples) is built using a streaming two-pass approach:

**Pass 1:** Aggregate positions across all samples → find peaks → cluster

**Pass 2:** For each sample, stream positions in 100k-row chunks (the manifest
loop in `analyze/manifest.py`), look up cluster membership via `IntervalTree`,
and accumulate per-`(cluster, sample)` counts in a dict.

`build_cluster_count_matrix()` assembles the final matrix from the streamed
positions and cluster definitions:

```python
def build_cluster_count_matrix(
    positions_df,                 # CPA positions with a `sample` column
    clusters_df,                  # cluster definitions from cluster_cpa_sites_adaptive
    sample_col='sample',
    position_col='corrected_position',
    count_col=None,
    fraction_col=None,            # optional fractional (proportional) weights
):
    """Cluster × sample count matrix via O(log n) IntervalTree lookup."""
```

The IntervalTree lookup is O(log n) per read, making this efficient for large datasets.

---

## Cluster parameters

| Flag (`rectify analyze`) | Default | Effect |
|--------------------------|---------|--------|
| `--cluster-distance` | 25 bp | Fixed-distance merge for the simple clustering mode |
| `--min-reads` | 5 | Discard clusters with fewer total reads |
| `--max-cluster-radius` | 10 bp | Adaptive: max distance from peak to cluster boundary |
| `--min-peak-sep` | 5 bp | Adaptive: minimum separation between distinct CPA peaks |
| `--min-cluster-samples` | 2 | Discard clusters seen in fewer than N samples |

---

## Output

`cpa_clusters.tsv` (adaptive mode adds `cluster_width`):

```tsv
cluster_id      chrom  strand  start   end     modal_position  cluster_com  n_positions  n_reads
cluster_000001  chrI   +       34500   34530   34521           34520        14           1205
cluster_000002  chrI   +       34580   34620   34601           34600        9            342
```

`modal_position` is the **median** position of the cluster's reads (an unweighted
median; the "modal" name is historical and it is *not* the most-used position).
`cluster_com` is the read-weighted center of mass (signal centroid). Clusters are
attributed to genes by the gene-attribution module during analysis.

---

## See also

- [APA Detection](apa_detection.md) — how clusters become APA isoforms
- [Multi-Sample Analysis](../user_guide/multi_sample.md) — the two-pass streaming pipeline
- `rectify analyze` — the command that runs clustering
