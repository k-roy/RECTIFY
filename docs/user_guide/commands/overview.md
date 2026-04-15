# Commands Overview

RECTIFY provides 14 subcommands, covering the full pipeline from alignment to visualization.

```
rectify <command> [options]
```

---

## Command summary

| Command | Description |
|:--------|:------------|
| [`rectify correct`](correct.md) | Correct 3' end positions (indel correction, A-tract resolution) |
| [`rectify run`](run.md) | Full pipeline: align (if FASTQ) → correct → analyze |
| [`rectify align`](align.md) | Multi-aligner consensus alignment (DRS-optimized) |
| [`rectify split`](split.md) | Split FASTQ into N equal chunks for parallel SLURM array alignment |
| [`rectify consensus`](consensus.md) | Aligner selection on pre-built per-aligner BAMs (post-merge step) |
| [`rectify analyze`](analyze.md) | Downstream analysis (clustering, DESeq2, GO, motifs) |
| [`rectify export`](export.md) | Export corrected 3' ends to bigWig/bedGraph tracks |
| [`rectify extract`](extract.md) | Extract per-read features from BAM to TSV |
| [`rectify aggregate`](aggregate.md) | Aggregate reads into 3' end, 5' end, and junction datasets |
| [`rectify netseq`](netseq.md) | Process NET-seq BAM files (3' end extraction, deconvolution) |
| [`rectify validate`](validate.md) | Validate corrections against ground truth (NET-seq, annotation) |
| [`rectify train-polya`](train_polya.md) | Train poly(A) tail model from control sites |
| [`rectify batch`](batch.md) | Generate SLURM array job scripts for multi-sample processing |
| [`rectify install-aligners`](install_aligners.md) | Download/install external aligners (deSALT, minimap2, gapmm2, uLTRA) |

---

## Typical workflows

### Single sample (FASTQ → corrected TSV)

```bash
rectify run reads.fastq.gz --Scer -o results/
```

### Single sample (BAM already aligned)

```bash
rectify correct reads.bam --Scer -o results/corrected.tsv
```

### Multi-sample with differential expression

```bash
rectify run --manifest manifest.tsv --Scer --reference wt -o results/
```

### Align only (for custom downstream use)

```bash
rectify align reads.fastq.gz --genome genome.fa --annotation genes.gff -o aligned.bam
```

### Parallel alignment via SLURM array (large datasets)

```bash
# 1. Split into 16 chunks and generate array scripts
rectify split reads.fastq.gz -n 16 -o chunks/ \
    --generate-slurm --aligners minimap2 mapPacBio gapmm2 uLTRA deSALT \
    --genome genome.fa --annotation genes.gff

# 2. Submit array job (16 chunks × 5 aligners = 80 tasks)
sbatch chunks/run_array_align.sh

# 3. After array completes: merge BAMs + run consensus
bash chunks/run_merge_and_consensus.sh
```

### Export bigWig tracks

```bash
rectify export corrected.tsv --genome genome.fa -o tracks/
```

### Validate against NET-seq ground truth

```bash
rectify validate corrected.tsv --netseq-dir netseq/ --Scer -o validation/
```

---

## Global options

All subcommands accept:

```
--help, -h       Show help message and exit
--version        Show RECTIFY version
--verbose, -v    Enable verbose (DEBUG) logging
```
