# APA Detection

Alternative polyadenylation (APA) detection identifies genes that use multiple CPA sites and quantifies how usage shifts between conditions.

**Implementation:** `rectify/core/analyze/apa_detection.py`, `rectify/core/analyze/shift_analysis.py`

---

## What is APA?

Many genes in yeast and other organisms have multiple cleavage and polyadenylation (CPA) sites, producing mRNA isoforms that differ in their 3' UTR length. Genes using CPA sites at different distances from the stop codon are said to undergo **alternative polyadenylation**.

The two archetypal isoforms are:

| Isoform | CPA site | 3' UTR | Regulatory consequence |
|---------|---------|--------|----------------------|
| **Proximal** | Upstream (near stop codon) | Short | Shorter 3' UTR; escapes miRNA/RBP regulation |
| **Distal** | Downstream (far from stop) | Long | Longer 3' UTR; more regulatory potential |

---

## APA isoform detection

RECTIFY groups CPA clusters into APA isoforms per gene:

```python
@dataclass
class APAIsoform:
    isoform_id: str
    gene_id: str
    cluster_ids: List[str]              # Which clusters constitute this isoform
    usage_per_sample: Dict[str, float]  # sample_id → fraction of gene's reads

@dataclass
class GeneAPAProfile:
    gene_id: str
    isoforms: List[APAIsoform]
    mean_read_count: int
    has_proximal_distal: bool           # Does this gene have clear prox/dist sites?
```

```python
def detect_apa_isoforms(gene_clusters, count_matrix):
    """
    Group clusters into isoforms by proximity and usage correlation.

    Returns: Dict[gene_id, GeneAPAProfile]
    """
```

---

## Proximal / distal site identification

For each gene with ≥ 2 clusters:

```python
def identify_proximal_distal_tes(gene_id, clusters, threshold_bp=200):
    """
    Classify clusters as proximal or distal based on distance from stop codon.

    threshold_bp: clusters within this distance of the stop codon are "proximal"

    Returns: (proximal_cluster_id, distal_cluster_id) or (None, None)
    """
```

For plus-strand genes:
- **Proximal** = cluster closest to stop codon (smallest genomic coordinate past stop)
- **Distal** = cluster furthest from stop codon (largest coordinate)

---

## APA shift analysis

`shift_analysis.py` quantifies whether the distribution of reads across CPA sites changes between conditions.

### Jensen-Shannon divergence

For each gene with ≥ 2 clusters, the shift is measured by Jensen-Shannon divergence between the condition and reference usage distributions:

```
JSD(P || Q) = ½ KL(P || M) + ½ KL(Q || M)
where M = ½(P + Q)
```

JSD is symmetric and bounded in [0, 1]:
- **JSD = 0**: identical usage distribution
- **JSD = 1**: completely different distributions

```python
def analyze_cluster_shifts(count_matrix, metadata, reference_condition):
    """
    Compute per-gene APA shifts.

    Returns: DataFrame with gene_id, js_divergence, direction, proximal_delta, distal_delta, padj
    """
```

### Shift direction

For genes with proximal and distal sites:

| `direction` | `proximal_delta` | Meaning |
|-------------|-----------------|---------|
| `proximal_shift` | positive | More reads at proximal site |
| `distal_shift` | negative | More reads at distal site |

### Statistical significance

P-values from a permutation test (1000 permutations of condition labels) are adjusted using Benjamini-Hochberg FDR. The default significance threshold is `padj < 0.05`.

---

## Dual-resolution DESeq2

RECTIFY complements APA analysis with DESeq2 at two resolutions:

| Level | Count unit | Detects |
|-------|-----------|---------|
| Gene | Sum of all cluster counts | Overall expression change |
| Cluster | Count per CPA cluster | Isoform-specific changes |

A gene where total expression is unchanged but CPA usage shifts will show:
- Non-significant in gene-level DESeq2
- Significant in cluster-level DESeq2 (one cluster up, another down)

---

## Output

| File | Content |
|------|---------|
| `tables/shift_results.tsv` | Gene-level shift analysis (JSD, direction, padj) |
| `tables/deseq2_clusters_*.tsv` | Cluster-level differential expression |
| `cpa_clusters.tsv` | Cluster definitions with peak positions |
| `cluster_gene_mapping.tsv` | Cluster → gene attribution |

---

## See also

- [Adaptive Clustering](adaptive_clustering.md) — how CPA positions become clusters
- [Output Formats](../user_guide/output_formats.md) — column definitions for `shift_results.tsv`
