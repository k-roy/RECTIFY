# Quick Start

This guide walks through the two most common workflows: single-sample correction and multi-sample analysis with differential expression.

---

## Single sample — bundled yeast data

If you are working with *S. cerevisiae*, RECTIFY bundles everything you need: genome, annotations, GO terms, and WT NET-seq. No external files required.

```bash
pip install rectify-rna

# From FASTQ (RECTIFY aligns + corrects)
rectify run reads.fastq.gz --Scer -o results/

# From pre-aligned BAM (correction only)
rectify correct reads.bam --Scer -o results/corrected.tsv
```

The output directory contains:

```
results/
├── corrected_3ends.tsv          # Per-read corrected 5'/3' positions
├── corrected_3ends_index.bed.gz # Aggregated position counts (300× smaller)
├── alignment_features.tsv       # Per-read alignment metadata
└── processing_stats.tsv         # QC summary
```

---

## Single sample — custom genome

```bash
# Align + correct with custom genome and annotation
rectify run reads.fastq.gz \
    --genome genome.fa.gz \
    --annotation genes.gff.gz \
    -o results/

# Correct only (BAM already aligned)
rectify correct reads.bam \
    --genome genome.fa.gz \
    --annotation genes.gff.gz \
    -o results/corrected.tsv
```

---

## Multi-sample analysis

For differential expression, APA detection, and motif discovery, provide a **manifest TSV** with one row per sample:

```
sample_id	path	condition
wt_rep1	/data/wt_rep1.fastq.gz	wt
wt_rep2	/data/wt_rep2.fastq.gz	wt
ko_rep1	/data/ko_rep1.fastq.gz	ko
ko_rep2	/data/ko_rep2.fastq.gz	ko
```

Run the full pipeline:

```bash
rectify run \
    --manifest manifest.tsv \
    --Scer \
    --reference wt \
    -o results/
```

Output structure:

```
results/
├── wt_rep1/
│   ├── corrected_3ends.tsv
│   └── corrected_3ends_index.bed.gz
├── wt_rep2/ ...
├── ko_rep1/ ...
├── ko_rep2/ ...
└── combined/
    ├── cpa_clusters.tsv
    ├── tables/
    │   ├── deseq2_genes_ko_vs_wt.tsv
    │   ├── deseq2_clusters_ko_vs_wt.tsv
    │   ├── motif_summary_ko.tsv
    │   └── go_enrichment_ko.tsv
    └── plots/
        ├── pca_samples.png
        └── heatmap_top_clusters.png
```

---

## Step-by-step: correction + analysis separately

If you have already corrected each sample, run analysis on the corrected TSVs:

```bash
# Step 1: correct each sample
for sample in wt_rep1 wt_rep2 ko_rep1 ko_rep2; do
    rectify correct raw/${sample}.bam --Scer -o results/${sample}/
done

# Step 2: run combined analysis using manifest
rectify analyze /dev/null \
    --manifest manifest.tsv \
    --Scer \
    --reference wt \
    --run-deseq2 \
    --run-motif \
    -o results/combined/
```

---

## Exporting genome browser tracks

```bash
rectify export results/wt_rep1/corrected_3ends.tsv \
    --genome genome.fa.gz \
    -o tracks/wt_rep1/
```

Produces `wt_rep1.bw` (bigWig) and `wt_rep1.bedgraph` for loading in IGV or UCSC.

---

## Processing NET-seq data

```bash
rectify netseq netseq.bam \
    --genome genome.fa.gz \
    --gff genes.gff.gz \
    -o netseq_output/
```

For *S. cerevisiae*, bundled WT NET-seq is auto-detected by `--Scer`. Provide custom NET-seq with `--netseq-dir`.

---

## HPC / SLURM submission

For large experiments, generate a SLURM array job script:

```bash
rectify run \
    --manifest manifest.tsv \
    --Scer \
    --reference wt \
    -o results/ \
    --profile rectify/slurm_profiles/sherlock_larsms.yaml
```

This generates `results/slurm/rectify_batch_correct.sh` (array job, one task per sample) and `results/slurm/rectify_batch_analyze.sh` (combined analysis, submitted after array). See the [HPC / SLURM guide](user_guide/hpc_slurm.md) for details.

---

## Next steps

- [Input Formats](user_guide/input_formats.md) — supported file types and manifest format
- [Output Formats](user_guide/output_formats.md) — column descriptions for all output files
- [Commands Overview](user_guide/commands/overview.md) — all available subcommands
- [Algorithms](algorithms/overview.md) — how the correction algorithms work
