# rectify analyze

Downstream analysis: CPA clustering, DESeq2, APA detection, GO enrichment, and motif discovery.

Usually called automatically by `rectify run-all --manifest`. Call directly when you have already corrected all samples and want to re-run analysis.

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
| `--run-motif` | off | Run STREME motif discovery (requires MEME Suite + `--genome`) |
| `--motif-upstream` | 100 | Window upstream of CPA for motif discovery (bp) |
| `--motif-downstream` | 50 | Window downstream of CPA for motif discovery (bp) |

GO enrichment runs automatically when `--go-annotations` is provided; the bundled yeast file is auto-selected with `--Scer`.

### Clustering parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--cluster-distance` | 25 | Max bp distance between sites to merge into one cluster |
| `--min-reads` | 5 | Minimum total reads per cluster |
| `--min-cluster-samples` | 2 | Minimum samples a cluster must appear in (filters singleton-sample noise) |
| `--max-cluster-radius` | 10 | Adaptive clustering: max bp from peak to cluster boundary |
| `--min-peak-sep` | 5 | Adaptive clustering: min bp between distinct CPA peaks |

### Performance

| Argument | Default | Description |
|----------|---------|-------------|
| `--threads` | 4 | Threads for parallel operations |

---

## Output files

| File | Description |
|------|-------------|
| `cpa_clusters.tsv` | CPA site clusters with read counts per sample |
| `tables/deseq2_genes_{condition}.tsv` | Gene-level DE results |
| `tables/deseq2_clusters_{condition}.tsv` | Cluster-level DE results |
| `tables/shift_analysis_{condition}.tsv` | APA shift analysis |
| `tables/go_enrichment_up_{condition}.tsv` / `go_enrichment_down_{condition}.tsv` | GO enrichment (up / down) |
| `tables/motif_summary_{condition}.tsv` | Enriched motifs |
| `motifs/*/` | STREME motif files and logos |
| `plots/pca_samples.png` | PCA of samples by cluster usage |
| `plots/sample_heatmap.png` | Sample-clustering heatmap |
| `report.html` | Comprehensive HTML summary |

---

## Notes

- DESeq2 requires at least 2 replicates per condition
- Motif discovery (`--run-motif`) requires the MEME Suite (`streme` on PATH); install via `conda install -c bioconda meme`
- Bedgraph generation is skipped in manifest mode; use `rectify export` for per-sample bigWigs
- The streaming pipeline never loads more than one sample's data into RAM — validated on 21 samples / 150M reads with 16 GB peak RAM
