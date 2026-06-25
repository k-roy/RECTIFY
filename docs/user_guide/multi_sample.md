# Multi-Sample Analysis

RECTIFY's multi-sample mode enables differential expression, APA detection, GO enrichment, and motif discovery across conditions. It uses a streaming pipeline that never loads more than one sample's data into memory, keeping peak RAM low regardless of dataset size.

---

## The manifest file

Provide a tab-separated manifest. `sample_id` and `path` are required;
`condition` is required for differential analysis (DESeq2 / APA / GO):

```tsv
sample_id	path	condition
wt_rep1	/data/wt_rep1.fastq.gz	wt
wt_rep2	/data/wt_rep2.fastq.gz	wt
ko_rep1	/data/ko_rep1.fastq.gz	ko
ko_rep2	/data/ko_rep2.fastq.gz	ko
```

Pass it with `--manifest`:

```bash
rectify run-all \
    --manifest manifest.tsv \
    --Scer \
    --reference wt \
    -o results/
```

---

## End-to-end workflow

`rectify run-all --manifest` runs the complete pipeline:

```text
For each sample (parallel):
    ├─ Align (if FASTQ) — multi-aligner consensus
    ├─ Correct 3' ends — walk-back + NET-seq refinement
    └─ Write corrected_reads.tsv

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

RECTIFY never loads all samples at once — it streams one sample at a time, so peak RAM stays roughly flat as the number of samples grows.

**Three memory tiers:**

1. **Column pruning** — on load, only `chrom`, `strand`, `corrected_3prime` are retained; all other columns are dropped immediately

2. **Two-pass streaming**:
    - *Pass 1*: each sample's TSV is read sequentially; positions aggregated to unique (chrom, strand, pos) → counts; combined for clustering
    - *Pass 2*: positions streamed in chunks; cluster membership looked up via IntervalTree; counts accumulated in a dict

3. **Optional position index** (`<tsv-stem>_index.bed.gz`, e.g. `corrected_reads_index.bed.gz`) — a pre-aggregated count file (~300× smaller). When present beside a sample's TSV, the loader uses it instead of the full TSV; otherwise it falls back to reading the TSV. Generate one with `write_position_index` (see below) if you want the speed-up.

### Reference condition matching

`--reference` is matched **case-insensitively** against the `condition` column. `--reference wt` matches `WT`, `Wt`, etc.

---

## Output structure

```text
results/
├── wt_rep1/
│   ├── corrected_reads.tsv
│   ├── corrected_consensus.bam
│   └── cpa_clusters.tsv
├── wt_rep2/ ...
├── ko_rep1/ ...
├── ko_rep2/ ...
└── combined/
    ├── cpa_clusters.tsv
    ├── tables/
    │   ├── deseq2_genes_ko.tsv
    │   ├── deseq2_clusters_ko.tsv
    │   ├── shift_analysis_ko.tsv
    │   ├── go_enrichment_up_ko.tsv
    │   ├── go_enrichment_down_ko.tsv
    │   └── motif_summary_ko.tsv
    ├── motifs/
    │   └── ko/
    │       ├── streme.html
    │       └── *.meme
    ├── plots/
    │   ├── pca_samples.png
    │   └── sample_heatmap.png
    └── report.html
```

---

## Generating position indices for existing data

If you corrected samples before the index feature was added, generate indices retroactively:

```python
from pathlib import Path
from rectify.core.bam.bam_processor import write_position_index
from concurrent.futures import ThreadPoolExecutor

samples = ['wt_rep1', 'wt_rep2', 'ko_rep1', 'ko_rep2']
base = Path('results/')

def gen_index(sample):
    tsv = str(base / sample / 'corrected_reads.tsv')
    write_position_index(tsv, tsv)  # writes <stem>_index.bed.gz alongside the TSV

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

`tables/shift_analysis_{condition}.tsv` reports genes whose CPA site usage distribution changes significantly between conditions, quantified by Jensen-Shannon divergence. Genes are classified as `proximal_shift` (more reads at upstream CPA sites) or `distal_shift` (more reads at downstream CPA sites).

---

## Limitations of manifest mode

- **Bedgraph generation is skipped** — generate per-sample bigWigs separately:
  ```bash
  rectify export results/wt_rep1/corrected_reads.tsv --genome genome.fa \
      --position-col corrected_3prime -o tracks/wt_rep1/
  ```
- **Genomic distribution analysis is skipped** (requires per-read alignment coordinates)
- Both steps can be run from the full per-read TSV if needed
