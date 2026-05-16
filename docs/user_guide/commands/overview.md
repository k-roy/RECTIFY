# Commands Overview

RECTIFY provides subcommands for the full pipeline: alignment, correction
(DRS, dT-primed-cDNA / QuantSeq REV, ONT PCR-cDNA), downstream analysis,
visualization, and utilities.

```
rectify <command> [options]
```

---

## Command summary

### Core pipeline (DRS / dT-primed-cDNA / QuantSeq REV)

| Command | Description |
|:--------|:------------|
| [`rectify correct`](correct.md) | Correct 3' end positions (indel correction, A-tract resolution). Supports DRS, dT-primed-cDNA (QuantSeq), and `--short-read` mode. |
| [`rectify run-all`](run.md) | Full pipeline: align (if FASTQ) -> correct -> analyze |
| [`rectify align`](align.md) | Multi-aligner consensus alignment (long-read default panel; `--short-read` selects bbmap + bwa) |
| [`rectify split`](split.md) | Split FASTQ into N equal chunks for parallel SLURM array alignment |
| [`rectify consensus`](consensus.md) | Aligner selection on pre-built per-aligner BAMs (post-merge step) |
| [`rectify analyze`](analyze.md) | Downstream analysis (clustering, DESeq2, GO, motifs) |

### ONT PCR-cDNA pipeline (PCB114.24)

| Command | Description |
|:--------|:------------|
| [`rectify correct-cdna`](correct_cdna.md) | Stage 1 — UMI extraction, directional clustering, abPOA consensus, pre-trim -> per-cluster FASTQ |
| [`rectify cdna-analyze`](cdna_analyze.md) | Stage 3 — post-align walkback, gene assignment, isoform clustering, T1/T2 pairing |

See the [ONT PCR-cDNA pipeline overview](correct_cdna_overview.md) for the three-stage workflow.

### Poly(A) utilities

| Command | Description |
|:--------|:------------|
| `rectify trim-polya` | DRS Step 0 — trim poly(A) tail + adapter from Dorado-aligned BAMs before alignment |
| `rectify tag-polya` | Annotate BAM reads with poly(A) model scores (`pt`, `ps`); only fills `pt` when Dorado's value is absent |
| `rectify restore-softclip` | DRS Step 4 — re-attach trimmed poly(A) as soft-clip on the corrected BAM (IGV visualization only) |
| [`rectify train-polya`](train_polya.md) | Train poly(A) tail model from control sites |

### Visualization, export, validation

| Command | Description |
|:--------|:------------|
| [`rectify export`](export.md) | Export corrected 3' ends to bigWig/bedGraph tracks |
| [`rectify extract`](extract.md) | Extract per-read features from BAM to TSV |
| [`rectify aggregate`](aggregate.md) | Aggregate reads into 3' end, 5' end, and junction datasets |
| [`rectify netseq`](netseq.md) | Process NET-seq BAM files (3' end extraction, deconvolution) |
| [`rectify validate`](validate.md) | Validate corrections against ground truth (NET-seq, annotation) |

### HPC and installation

| Command | Description |
|:--------|:------------|
| [`rectify batch`](batch.md) | Generate SLURM array job scripts for multi-sample processing |
| [`rectify install-aligners`](install_aligners.md) | Download/install external aligners (deSALT, minimap2, gapmm2, uLTRA) |

---

## Typical workflows

### Single sample (FASTQ → corrected TSV)

```bash
rectify run-all reads.fastq.gz --Scer -o results/
```

### Single sample (BAM already aligned)

```bash
rectify correct reads.bam --Scer -o results/corrected.tsv
```

### Multi-sample with differential expression

```bash
rectify run-all --manifest manifest.tsv --Scer --reference wt -o results/
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
