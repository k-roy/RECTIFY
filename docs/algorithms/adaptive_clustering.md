# Adaptive Clustering

After per-read correction, RECTIFY groups nearby CPA positions into **clusters** using a valley-based adaptive algorithm. Clusters represent distinct CPA sites used by a gene.

**Implementation:** `rectify/core/analyze/clustering.py`

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

Each cluster extends from valley to valley. The cluster peak is the most-used position within the cluster.

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
    positions,            # Dict: chrom → [{'pos': int, 'strand': str, 'count': float}]
    min_cluster_size=10,  # Minimum reads to keep a cluster
    prominence_threshold=0.3,  # Min peak prominence relative to neighbors
):
    """
    Valley-based CPA clustering.

    Returns: DataFrame with cluster_id, chrom, strand, start, end, peak_pos, count
    """
```

---

## Count matrix building

After clustering, a count matrix (clusters × samples) is built using a streaming two-pass approach:

**Pass 1:** Aggregate positions across all samples → find peaks → cluster

**Pass 2:** For each sample, stream reads in 100k-row chunks, look up cluster membership via `IntervalTree`, accumulate in a dict

```python
def build_cluster_count_matrix(manifest_path, clusters, chunk_size=100_000):
    """
    Streaming count matrix. Peak RAM: O(clusters × samples).

    Returns: count_matrix (DataFrame), sample_ids, cluster_ids
    """
```

The IntervalTree lookup is O(log n) per read, making this efficient for large datasets.

---

## Cluster parameters

| Parameter | Default | Effect |
|-----------|---------|--------|
| `--merge-distance` | 25 bp | Fixed-distance merge for simple clustering mode |
| `--min-reads` | 5 | Discard clusters with fewer total reads |
| `--min-samples` | 1 | Discard clusters seen in fewer samples |
| `prominence_threshold` | 0.3 | Discard peaks shorter than 30% of the tallest neighbor |

---

## Output

`cpa_clusters.tsv`:

```tsv
cluster_id   chrom  strand  start   end     peak_pos  total_count
CLU_chrI_1   chrI   +       34500   34530   34521     1205
CLU_chrI_2   chrI   +       34580   34620   34601     342
```

Each cluster is mapped to a gene via `cluster_gene_mapping.tsv` using the gene attribution module.

---

## See also

- [APA Detection](apa_detection.md) — how clusters become APA isoforms
- [Multi-Sample Analysis](../user_guide/multi_sample.md) — the two-pass streaming pipeline
- `rectify analyze` — the command that runs clustering
