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
| `--dT-primed-cDNA` | off | Input is dT-primed cDNA (QuantSeq, etc.) — poly(A) NOT in read. Enables indel artifact correction and poly(A) trimming modules. |
| `--short-read` | off | Input is short-read data (Illumina/Aviti ≤150 bp). Uses bbmap + bwa instead of the long-read aligner panel and disables poly(A)/A-tract/indel modules. |
| `--reference` | auto | Reference condition for DESeq2 (case-insensitive) |

*(deprecated alias retained for backwards compatibility: `--no-polya-sequenced` (synonym for `--dT-primed-cDNA`). Use the current flags above for new pipelines.)*

!!! note "ONT PCR-cDNA (SQK-PCB114.24, SSP + UMI)"
    `run-all` does not yet wire the `--ONT-cDNA` mode end-to-end. For ONT
    PCR-cDNA libraries with UMI adapters, run the cDNA-specific subcommands
    directly: `rectify correct-cdna` (Stage 1 UMI consensus) →
    `rectify align` → `rectify cdna-analyze`. See the
    [ONT PCR-cDNA pipeline guide](correct_cdna_overview.md) for the full
    sequence.

### NET-seq

| Argument | Description |
|----------|-------------|
| `--netseq-dir` | Custom NET-seq BigWig directory |

### Aligner options

| Argument | Default | Description |
|----------|---------|-------------|
| `--base-aligners` | minimap2 mapPacBio gapmm2 (long-read) / bbmap bwa (`--short-read`) | Aligners to include in the consensus pool |
| `--junction-aligners` | uLTRA deSALT (long-read) / [] (short-read) | Junction-aware aligners (requires `--annotation`) |
| `--no-junction-aligners` | — | Explicitly disable uLTRA + deSALT |
| `--parallel-aligners` | off | Run base aligners in parallel (divides `--threads` evenly) |
| `--chimeric-consensus` | on | Per-segment aligner selection (use `--no-chimeric-consensus` to disable) |
| `--ultra-path` | `uLTRA` | Path to uLTRA executable |
| `--desalt-path` | `deSALT` | Path to deSALT executable |

### Chunked alignment (HPC array submission)

| Argument | Default | Description |
|----------|---------|-------------|
| `--chunked-alignment` | off | Generate scheduler array scripts for chunked parallel alignment instead of running inline |
| `--force-no-chunking` | off | DANGER: bypass the chunking requirement for FASTQ inputs (local workstation only) |
| `--target-reads-per-chunk` | 500000 | Target reads per alignment chunk when `--chunked-alignment` is set |
| `--scheduler` | slurm | Target scheduler for generated script headers (`slurm`, `uge`, `pbs`) |
| `--slurm-partition` / `--slurm-account` | — | Header values for SLURM scripts |
| `--uge-queue` / `--uge-pe` | `long.q` / `smp` | Header values for UGE/SGE scripts |
| `--pbs-queue` | `workq` | Header value for PBS/Torque scripts |
| `--python-path` | `sys.executable` | Explicit Python interpreter for generated scripts |
| `--rectify-src` | auto | Path to RECTIFY source checkout for generated scripts |

### Analysis options

| Argument | Description |
|----------|-------------|
| `--go-annotations` | GO annotation file for enrichment analysis |

### Filtering

| Argument | Description |
|----------|-------------|
| `--filter-spikein GENE [GENE ...]` | Remove spike-in reads by gene name |

### Junction scoring pass-through

These flags are forwarded verbatim to `rectify correct`. See the
[rectify correct reference](correct.md#junction-refinement-module-2h) for details.

| Argument | Default | Description |
|----------|---------|-------------|
| `--junction-penalty-table PATH` | (heuristic) | Empirical HP-context penalty table (`penalty_scores.tsv`) |
| `--str-penalty-table PATH` | (none) | STR slippage penalty table (`str_penalty_scores.tsv`) |
| `--junction-overhang-table PATH` | (none) | Empirical junction overhang table produced by `calibrate_junction_overhang.py`. Penalises aligner wins that rely on junctions lacking sufficient flanking overhang for their intron size. Short introns (<500 bp, ≥5 cross-read support) are exempt. |

### BAM output options

| Argument | Default | Description |
|----------|---------|-------------|
| `--write-softclip-bam` | off | Write a soft-clipped corrected BAM alongside the primary hard-clipped BAM. Cat2 soft-clip rescue bases are visible in IGV with "Show soft-clipped bases" enabled. Useful for QC/debugging. |
| `--write-polya-bam` | off | Write a BAM with the original poly(A) tail restored from the DRS trim parquet metadata as a 3' soft clip. Useful for visually validating poly(A) tail handling in IGV. Has no effect without `--drs`. |
| `--bam-dir DIR` | (sample output dir) | Directory to write alignment BAMs (per-aligner and rectified). Useful for inspecting per-aligner BAMs separately from corrected outputs. |
| `--keep-aligner-bams` | off | Retain per-aligner BAMs (minimap2, mapPacBio, gapmm2) after consensus selection. By default these are excluded from the Oak sync to save disk space. |
| `--trust-existing-bams` | off | Reuse pre-existing per-aligner / rectified BAMs without checking provenance sidecars. Use only when BAMs have been manually verified as compatible. |
| `--scratch-dir DIR` | auto | Base directory for intermediate BAM I/O (alignment output, unsorted corrected BAM, sort temp files). Auto-detected from `$SCRATCH`, `$SLURM_TMPDIR`, or `$TMPDIR` inside HPC batch jobs. |

### Performance

| Argument | Default | Description |
|----------|---------|-------------|
| `--threads` | 0 | Number of threads (0 = auto-detect from `SLURM_CPUS_PER_TASK` or CPU count) |
| `--streaming` | off | Streaming mode for large BAMs |
| `--chunk-size` | 10000 | Reads per chunk (streaming mode) |
| `--mapPacBio-chunks` | 1 | Split FASTQ into N chunks for mapPacBio parallelism |
| `--aligner-concurrency` | auto | Run multiple aligner corrections concurrently with a shared worker budget. `auto` disables this on small laptops and reserves 2 CPUs for merge/main-process overhead on cluster nodes. Use `1` to preserve sequential-aligner correction behavior. |
| `--continue-on-error` | off | Continue processing remaining samples if one sample fails (multi-sample mode only). |

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
