# rectify run-all

Full end-to-end pipeline: align (if FASTQ) → correct → analyze.

This is the recommended entry point for most workflows. It handles single samples and multi-sample experiments via manifest.

---

## Usage

```bash
# Single sample
rectify run-all <input> [options] -o <output_dir>

# Multi-sample (manifest)
rectify run-all --manifest manifest.tsv [options] -o <output_dir>
```

## Examples

```bash
# Single sample — bundled yeast data
rectify run-all reads.fastq.gz --Scer -o results/

# Single sample — pre-aligned BAM
rectify run-all reads.bam --Scer -o results/

# Custom genome
rectify run-all reads.fastq.gz \
    --genome genome.fa.gz \
    --annotation genes.gff.gz \
    -o results/

# DRS workflow (Dorado direct RNA BAM)
rectify run-all sample_dorado.bam --drs --Scer -o results/

# Multi-sample differential expression
rectify run-all \
    --manifest manifest.tsv \
    --Scer \
    --reference wt \
    -o results/

# Multi-sample with SLURM array jobs
rectify run-all \
    --manifest manifest.tsv \
    --Scer \
    --reference wt \
    -o results/ \
    --profile slurm_profiles/hpc_cpu.yaml
```

---

## Arguments

### Input (one required)

| Argument | Description |
|----------|-------------|
| `input` | Input FASTQ.GZ or BAM file (positional, single sample) |
| `--manifest, -m` | Sample manifest TSV for multi-sample mode |

### Reference data

| Argument | Description |
|----------|-------------|
| `--genome` | Reference genome FASTA (optionally gzipped) |
| `--annotation` | Gene annotation file (GTF or GFF, optionally gzipped) |
| `--Scer` | Use bundled *S. cerevisiae* data |
| `--organism` | Organism name for bundled data (e.g. `yeast`) |

### Required

| Argument | Description |
|----------|-------------|
| `-o, --output-dir` | Output directory |

### Pipeline options

| Argument | Default | Description |
|----------|---------|-------------|
| `--skip-alignment` | off | Skip alignment step (FASTQ input with pre-aligned BAM) |
| `--drs` | off | Enable DRS workflow: adds Step 0 poly(A)+adapter trimming before alignment and Step 4 soft-clip restore after correction (BAM input only) |
| `--no-polya-sequenced` | off | Disable poly(A) trimming and indel correction |
| `--reference` | auto | Reference condition for DESeq2 (case-insensitive) |

### NET-seq

| Argument | Description |
|----------|-------------|
| `--netseq-dir` | Custom NET-seq BigWig directory |

### Aligner options

| Argument | Default | Description |
|----------|---------|-------------|
| `--aligner` | minimap2 | Primary aligner (`minimap2`, `mapPacBio`, `gapmm2`) |
| `--junction-aligners` | — | Optional junction-mode aligners: `uLTRA`, `deSALT` |
| `--parallel-aligners` | off | Run multiple aligners in parallel |
| `--chimeric-consensus` | off | Use chimeric CIGAR assembly (experimental) |
| `--ultra-path` | auto | Path to uLTRA executable |
| `--desalt-path` | auto | Path to deSALT executable |

### Analysis options

| Argument | Description |
|----------|-------------|
| `--go-annotations` | GO annotation file for enrichment analysis |

### Filtering

| Argument | Description |
|----------|-------------|
| `--filter-spikein GENE [GENE ...]` | Remove spike-in reads by gene name |

### Performance

| Argument | Default | Description |
|----------|---------|-------------|
| `--threads` | 4 | Number of threads |
| `--streaming` | off | Streaming mode for large BAMs |
| `--chunk-size` | 10000 | Reads per chunk (streaming mode) |
| `--mapPacBio-chunks` | 1 | Split FASTQ into N chunks for mapPacBio parallelism |

### HPC

| Argument | Description |
|----------|-------------|
| `--profile` | SLURM profile YAML for cluster submission |

---

## Pipeline steps

### Single sample (standard)

```
1. Alignment (skipped if BAM provided)
   ├─ minimap2 (splice-aware, junction-annotated)
   ├─ mapPacBio (PacBio RNA mode)
   ├─ gapmm2 (gap-aware variant)
   └─ Consensus selection by junction scoring

2. Correction (rectify correct — post-consensus, on winning aligner's BAM)
   ├─ Indel artifact correction (walk-back)
   ├─ A-tract ambiguity detection
   ├─ Soft-clip rescue (Module 2G)
   ├─ Junction refinement (Module 2H — requires --aligner-bams)
   ├─ NET-seq refinement (if available)
   └─ Spike-in filtering

3. Analysis (rectify analyze)
   ├─ corrected_3ends.tsv
   ├─ corrected_3ends_index.bed.gz
   └─ processing_stats.tsv
```

### Single sample (DRS — with `--drs`)

```
0. Poly(A) + adapter trimming (rectify trim-polya — 3-pass)
1. Alignment (same as above, on trimmed FASTQ)
2. Correction (rectify correct — post-consensus)
3. Analysis (rectify analyze)
4. Restore poly(A) as soft-clips (rectify restore-softclip — for IGV)
```

### Multi-sample (`--manifest`)

```
For each sample (parallel):
    Steps 0–2 above

Combined analysis:
    ├─ Adaptive CPA clustering (valley-based)
    ├─ Count matrix building (streaming, low RAM)
    ├─ DESeq2 (gene-level + cluster-level)
    ├─ APA shift analysis
    ├─ GO enrichment
    ├─ Motif discovery (STREME)
    └─ HTML report
```

!!! note "DESeq2 requires replicates"
    Differential expression analysis is only run when there are ≥ 2 samples per condition.

---

## Output structure

```
output_dir/
├── <sample>/                   # One directory per sample
│   ├── corrected_3ends.tsv
│   ├── corrected_3ends_index.bed.gz
│   ├── alignment_features.tsv
│   └── processing_stats.tsv
└── combined/                   # Multi-sample only
    ├── cpa_clusters.tsv
    ├── tables/
    ├── plots/
    └── analysis_report.html
```
