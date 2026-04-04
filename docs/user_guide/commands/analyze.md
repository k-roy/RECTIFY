# rectify analyze

Downstream analysis: CPA clustering, DESeq2, APA detection, GO enrichment, and motif discovery.

Usually called automatically by `rectify run --manifest`. Call directly when you have already corrected all samples and want to re-run analysis.

---

## Usage

```bash
rectify analyze /dev/null \
    --manifest manifest.tsv \
    --annotation genes.gff.gz \
    --reference wt \
    --run-deseq2 \
    --run-motif \
    -o results/combined/
```

!!! note "The positional input"
    The positional argument is a legacy placeholder. Pass `/dev/null` (or any path) when using `--manifest` mode.

---

## Examples

```bash
# Full analysis — bundled yeast
rectify analyze /dev/null \
    --manifest manifest.tsv \
    --Scer \
    --reference wt \
    --run-deseq2 \
    --run-motif \
    -o results/combined/

# Analysis without motif discovery (faster, no MEME required)
rectify analyze /dev/null \
    --manifest manifest.tsv \
    --Scer \
    --reference wt \
    --run-deseq2 \
    -o results/combined/

# Custom genome + GO annotations
rectify analyze /dev/null \
    --manifest manifest.tsv \
    --genome genome.fa.gz \
    --annotation genes.gff.gz \
    --go-annotations go.tsv.gz \
    --reference wt \
    --run-deseq2 \
    --run-motif \
    -o results/combined/
```

---

## Arguments

### Input

| Argument | Description |
|----------|-------------|
| `input` | Positional placeholder (use `/dev/null` with `--manifest`) |
| `--manifest, -m` | Sample manifest TSV (sample_id, path, condition) |

### Reference data

| Argument | Description |
|----------|-------------|
| `--annotation` | Gene annotation GFF/GTF |
| `--Scer` | Bundled *S. cerevisiae* data |
| `--genome` | Reference genome FASTA (needed for motif sequence extraction) |
| `--go-annotations` | GO annotation TSV for enrichment |

### Required

| Argument | Description |
|----------|-------------|
| `-o, --output-dir` | Output directory |
| `--reference` | Reference condition name (case-insensitive) |

### Analysis modules

| Argument | Default | Description |
|----------|---------|-------------|
| `--run-deseq2` | off | Run DESeq2 differential expression |
| `--run-motif` | off | Run STREME motif discovery (requires MEME Suite) |
| `--run-go` | off | Run GO enrichment |
| `--run-shift` | on | APA shift analysis (proximal/distal usage changes) |

### Clustering parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--merge-distance` | 25 | Max bp distance between sites to merge into one cluster |
| `--min-reads` | 5 | Minimum total reads per cluster |
| `--min-samples` | 1 | Minimum samples contributing reads to a cluster |

### Performance

| Argument | Default | Description |
|----------|---------|-------------|
| `--threads` | 4 | Threads for parallel operations |

---

## Output files

| File | Description |
|------|-------------|
| `cpa_clusters.tsv` | CPA site clusters with read counts per sample |
| `cluster_gene_mapping.tsv` | Cluster → gene attribution |
| `tables/deseq2_genes_*.tsv` | Gene-level DE results |
| `tables/deseq2_clusters_*.tsv` | Cluster-level DE results |
| `tables/shift_results.tsv` | APA shift analysis |
| `tables/go_enrichment_*.tsv` | GO enrichment |
| `tables/motif_summary_*.tsv` | Enriched motifs |
| `motifs/*/` | STREME motif files and logos |
| `plots/pca_samples.png` | PCA of samples by cluster usage |
| `plots/heatmap_top_clusters.png` | Heatmap of top DE clusters |
| `analysis_report.html` | Comprehensive HTML summary |

---

## Notes

- DESeq2 requires at least 2 replicates per condition
- Motif discovery (`--run-motif`) requires the MEME Suite (`streme` on PATH); install via `conda install -c bioconda meme`
- Bedgraph generation is skipped in manifest mode; use `rectify export` for per-sample bigWigs
- The streaming pipeline never loads more than one sample's data into RAM — validated on 21 samples / 150M reads with 16 GB peak RAM
