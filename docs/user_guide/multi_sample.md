# Multi-Sample Analysis

RECTIFY's multi-sample mode enables differential expression, APA detection, GO enrichment, and motif discovery across conditions. It uses a streaming pipeline that keeps peak RAM to a few MB regardless of dataset size.

---

## The manifest file

Provide a tab-separated manifest with three required columns:

```tsv
sample_id	path	condition
wt_rep1	/data/wt_rep1.fastq.gz	wt
wt_rep2	/data/wt_rep2.fastq.gz	wt
ko_rep1	/data/ko_rep1.fastq.gz	ko
ko_rep2	/data/ko_rep2.fastq.gz	ko
```

Pass it with `--manifest`:

```bash
rectify run \
    --manifest manifest.tsv \
    --Scer \
    --reference wt \
    -o results/
```

---

## End-to-end workflow

`rectify run --manifest` runs the complete pipeline:

```
For each sample (parallel):
    ├─ Align (if FASTQ) — multi-aligner consensus
    ├─ Correct 3' ends — walk-back + NET-seq refinement
    └─ Write corrected_3ends.tsv + corrected_3ends_index.bed.gz

Combined analysis (after all samples):
    ├─ Pass 1: stream positions → adaptive clustering
    ├─ Pass 2: stream positions → cluster count matrix
    ├─ DESeq2 (gene-level + cluster-level)
    ├─ APA shift analysis
    ├─ GO enrichment
    ├─ Motif discovery (STREME)
    └─ HTML report
```

---

## Running analysis on existing corrected data

If you already have corrected TSVs:

```bash
rectify analyze /dev/null \
    --manifest manifest.tsv \
    --annotation genes.gff.gz \
    --reference wt \
    --run-deseq2 \
    --run-motif \
    -o results/combined/
```

The `manifest.tsv` `path` column should point to the corrected TSVs (or their parent directories).

---

## How the streaming pipeline works

### Memory efficiency

RECTIFY never loads all samples at once. For 21 samples / 150M reads, peak RAM is ~4 GB regardless of depth.

**Three memory tiers:**

1. **Column pruning** — on load, only `chrom`, `strand`, `corrected_position` are retained; all other columns are dropped immediately

2. **Two-pass streaming**:
    - *Pass 1*: each sample's TSV is read sequentially; positions aggregated to unique (chrom, strand, pos) → counts; combined for clustering
    - *Pass 2*: positions streamed in 100k-row chunks; cluster membership looked up via IntervalTree; counts accumulated in a dict

3. **Position index** (`corrected_3ends_index.bed.gz`) — written by `rectify correct` alongside the full TSV; pre-aggregated counts ~300× smaller; when present, both passes use it instead of the full TSV

### Reference condition matching

`--reference` is matched **case-insensitively** against the `condition` column. `--reference wt` matches `WT`, `Wt`, etc.

---

## Output structure

```
results/
├── wt_rep1/
│   ├── corrected_3ends.tsv
│   ├── corrected_3ends_index.bed.gz
│   ├── alignment_features.tsv
│   └── processing_stats.tsv
├── wt_rep2/ ...
├── ko_rep1/ ...
├── ko_rep2/ ...
└── combined/
    ├── cpa_clusters.tsv
    ├── cluster_gene_mapping.tsv
    ├── tables/
    │   ├── deseq2_genes_ko_vs_wt.tsv
    │   ├── deseq2_clusters_ko_vs_wt.tsv
    │   ├── shift_results.tsv
    │   ├── go_enrichment_ko.tsv
    │   └── motif_summary_ko.tsv
    ├── motifs/
    │   └── ko/
    │       ├── streme.html
    │       └── *.meme
    ├── plots/
    │   ├── pca_samples.png
    │   └── heatmap_top_clusters.png
    └── analysis_report.html
```

---

## Generating position indices for existing data

If you corrected samples before the index feature was added, generate indices retroactively:

```python
from pathlib import Path
from rectify.core.bam_processor import write_position_index
from concurrent.futures import ThreadPoolExecutor

samples = ['wt_rep1', 'wt_rep2', 'ko_rep1', 'ko_rep2']
base = Path('results/')

def gen_index(sample):
    tsv = str(base / sample / 'corrected_3ends.tsv')
    write_position_index(tsv, tsv)  # writes alongside the TSV

with ThreadPoolExecutor(max_workers=4) as ex:
    list(ex.map(gen_index, samples))
```

---

## DESeq2: gene-level vs cluster-level

RECTIFY runs DESeq2 at two resolutions:

| Level | What it detects | Example |
|-------|----------------|---------|
| **Gene** | Total expression changes | *HSP82* is 2-fold down in heat shock |
| **Cluster** | CPA site usage changes | *FAS1* shifts from distal to proximal site |

Cluster-level analysis detects isoform-specific changes that gene-level counts would miss (because a proximal-to-distal shift keeps total counts constant while changing isoform usage).

---

## APA shift analysis

`shift_results.tsv` reports genes whose CPA site usage distribution changes significantly between conditions, quantified by Jensen-Shannon divergence. Genes are classified as `proximal_shift` (more reads at upstream CPA sites) or `distal_shift` (more reads at downstream CPA sites).

---

## Limitations of manifest mode

- **Bedgraph generation is skipped** — generate per-sample bigWigs separately:
  ```bash
  rectify export results/wt_rep1/corrected_3ends.tsv --genome genome.fa -o tracks/wt_rep1/
  ```
- **Genomic distribution analysis is skipped** (requires per-read alignment coordinates)
- Both steps can be run from the full per-read TSV if needed
