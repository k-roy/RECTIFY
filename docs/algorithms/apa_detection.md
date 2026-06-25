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

RECTIFY groups reads into APA isoforms per gene by the triplet
`(gene, junction signature, 3' cluster)`:

```python
@dataclass
class APAIsoform:
    isoform_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    tes_position: int                   # Transcript end site (3' end) position
    tes_cluster_id: Optional[str]       # CPA cluster this TES belongs to
    junction_signature: Tuple[Tuple[int, int], ...]  # Sorted (donor, acceptor) pairs
    supporting_reads: List[str]         # Read IDs supporting this isoform
    is_canonical_tes: bool              # TES matches the annotated gene end
```

```python
def detect_apa_isoforms(records, gene_attributions,
                        cluster_assignments=None, min_reads_per_isoform=3,
                        annotated_tes=None, tes_tolerance=50):
    """
    Detect APA isoforms by grouping reads by (gene, junction signature, TES
    cluster).

    Returns: List[APAIsoform]
    """
```

---

## Proximal / distal site identification

For genes with multiple transcript end sites (TES):

```python
def identify_proximal_distal_tes(profiles, min_tes_difference=100):
    """
    Classify each gene's TES usage as proximal (shorter 3' UTR) vs distal
    (longer 3' UTR) based on genomic position and strand.

    profiles: Dict[gene_id, GeneAPAProfile]
    min_tes_difference: minimum bp between TES to classify prox/dist

    Returns: DataFrame of per-gene proximal/distal TES assignments
    """
```

- **Proximal** = the TES giving the shorter 3' UTR (closer to the stop codon)
- **Distal** = the TES giving the longer 3' UTR (farther from the stop codon)

(Strand is accounted for: "closer to the stop codon" means smaller genomic
coordinate on the plus strand and larger on the minus strand.)

---

## APA shift analysis

`shift_analysis.py` quantifies whether the distribution of reads across CPA sites changes between conditions.

### Jensen-Shannon divergence

For each gene with ≥ 2 clusters, the shift magnitude is measured by the
Jensen-Shannon divergence between the two conditions' cluster-usage
distributions. RECTIFY computes it via `scipy.spatial.distance.jensenshannon`
(base 2) and squares the returned distance, so the reported
`distribution_divergence` is the JSD itself, bounded in [0, 1]:

- **0**: identical usage distribution
- **1**: completely disjoint distributions

```python
def analyze_cluster_shifts(count_matrix, clusters_df, condition_a, condition_b,
                           sample_metadata, cluster_gene_attributions=None,
                           min_total_counts=20):
    """
    Compare cluster-usage distributions between two conditions, per gene.

    Returns a DataFrame with: gene_id, gene_name, n_clusters,
    major_cluster_a, major_cluster_b, shift_bp, shift_direction,
    distribution_divergence.
    """
```

### Shift direction

`shift_bp` is the change in the major cluster's position between conditions,
and `shift_direction` summarizes it:

| `shift_direction` | Meaning |
|-------------------|---------|
| `downstream` | major CPA site moved downstream (3' UTR lengthening) |
| `upstream` | major CPA site moved upstream (3' UTR shortening) |
| `same` | major cluster unchanged |
| `unknown` | a major cluster could not be determined in one condition |

A boolean `same_major` column is also emitted (true when the major cluster id is
identical across the two conditions).

This function reports shift magnitude and direction; it does not itself compute
a per-gene p-value. (Cluster- and gene-level differential-usage significance —
`padj`, `log2FoldChange` — comes from the DESeq2 step below.)

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
| `tables/shift_analysis_{condition}.tsv` | Gene-level shift analysis (`shift_bp`, `shift_direction`, `distribution_divergence`) |
| `tables/deseq2_clusters_{condition}.tsv` | Cluster-level differential usage (`padj`, `log2FoldChange`, `baseMean`) |
| `tables/deseq2_genes_{condition}.tsv` | Gene-level differential usage |
| `cpa_clusters.tsv` | Cluster definitions (`modal_position`, `cluster_com`, …) |

---

## See also

- [Adaptive Clustering](adaptive_clustering.md) — how CPA positions become clusters
- [Output Formats](../user_guide/output_formats.md) — column definitions for `shift_analysis_{condition}.tsv`
