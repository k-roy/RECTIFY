# HPC / SLURM

RECTIFY has built-in support for HPC batch schedulers. The `rectify batch` command generates
scheduler scripts and, optionally, submits them automatically.

Supported schedulers: **SLURM**, **UGE/SGE** (Altair Grid Engine), **PBS/Torque**.

---

## Quick start

```bash
# 1. Copy and edit the generic profile for your cluster
cp rectify/slurm_profiles/slurm_generic.yaml my_cluster.yaml
# Edit: set `partition` to your cluster's partition or queue name

# 2. Generate (and optionally submit) SLURM array scripts
rectify batch \
    --manifest manifest.tsv \
    --genome genome.fa \
    --annotation genes.gff \
    --reference wt \
    -o results/ \
    --profile my_cluster.yaml
```

This generates two scheduler scripts:

1. `results/slurm/rectify_batch_correct.sh` — array job (one task per sample)
2. `results/slurm/rectify_batch_analyze.sh` — combined analysis (runs after array)

To submit immediately, add `submit: true` to your profile, or pass `--submit` on the CLI.

---

## Bundled profiles

Three profiles are included in `rectify/slurm_profiles/`:

| Profile | Purpose |
|---------|---------|
| `slurm_generic.yaml` | **Start here.** Generic template with all options documented. |
| `hpc_cpu.yaml` | Generic CPU partition template (use_scratch, streaming on by default). |
| `hpc_gpu.yaml` | Generic GPU partition template. |

For most clusters, start with `slurm_generic.yaml` or `hpc_cpu.yaml`.

---

## Profile reference

All profile fields (with defaults from `slurm_generic.yaml`):

```yaml
partition: null           # REQUIRED: your partition/queue name (null = cluster default)
cpus: 8                   # CPUs per correction task
mem: "32G"                # Memory per correction task
time: "4:00:00"           # Wall-clock limit per correction task
analyze_cpus: 8           # CPUs for combined analysis
analyze_mem: "64G"        # Memory for combined analysis
analyze_time: "8:00:00"   # Wall-clock limit for analysis
job_name: rectify_batch
max_concurrent: 10        # Max simultaneous array tasks
submit: false             # true = submit after generating scripts
use_scratch: false        # true = stage BAMs through $SCRATCH (see below)
streaming: true           # true = stream BAM reads in chunks (recommended)
```

CLI flags override profile values. For example, `--cpus 16` overrides `cpus: 8`.

---

## Scratch staging

When `use_scratch: true` is set and `$SCRATCH` is defined in the job environment,
the generated script:

1. Copies the input BAM to `$SCRATCH` (typically node-local, high-bandwidth)
2. Runs `rectify correct` on the local copy
3. Rsyncs all outputs back to the output directory
4. Cleans up `$SCRATCH`

This avoids I/O contention on shared NFS/Lustre filesystems when many array tasks
run concurrently. On clusters with fast local storage, this typically reduces
correction wall time by 2–3×.

!!! note
    FASTQ inputs are left on shared storage — they are read once sequentially and
    do not benefit from staging.

!!! note
    `$SCRATCH` variable convention varies by cluster. Common names: `$SCRATCH`,
    `$SLURM_TMPDIR`, `$TMPDIR`. RECTIFY checks these in priority order. If none
    is set, I/O falls back to the output directory.

---

## Thread limits

**Batch schedulers will suspend or terminate jobs that use more CPU threads than allocated.**

Python packages (numpy, sklearn, pydeseq2) spawn threads via OpenMP, OpenBLAS, MKL,
and joblib's loky backend. RECTIFY's generated scripts set all limits automatically
using `$RECTIFY_CPUS` (normalised from `SLURM_CPUS_PER_TASK`, `NSLOTS`, or `PBS_NUM_PPN`):

```bash
export OMP_NUM_THREADS=$RECTIFY_CPUS
export OPENBLAS_NUM_THREADS=$RECTIFY_CPUS
export MKL_NUM_THREADS=$RECTIFY_CPUS
export LOKY_MAX_CPU_COUNT=$RECTIFY_CPUS   # critical for pydeseq2/sklearn
```

!!! warning "LOKY_MAX_CPU_COUNT is mandatory"
    Joblib's loky backend ignores `JOBLIB_WORKERS`. You **must** set
    `LOKY_MAX_CPU_COUNT` or pydeseq2 and sklearn will over-subscribe.

If writing custom scripts, set these **before** importing numpy:

```python
import os
os.environ['OMP_NUM_THREADS'] = str(n_cpus)
os.environ['LOKY_MAX_CPU_COUNT'] = str(n_cpus)
# now safe to import numpy, pandas, sklearn
```

---

## Streaming mode

For large BAMs (> 2 GB), use streaming mode to keep peak RAM at ~4–5 GB regardless
of file size:

```bash
rectify correct reads.bam --genome genome.fa --annotation genes.gff --streaming -o results/
```

The `slurm_generic.yaml` profile sets `streaming: true` by default.

Without streaming, all reads are accumulated in RAM before writing. For a 7 GB BAM
with 40M reads this requires ~30–40 GB.

---

## BAM output options

By default, per-aligner BAMs (minimap2, mapPacBio, gapmm2) are discarded after
consensus selection; only the rectified BAM is kept. Use these flags if you want
to inspect the per-aligner outputs:

```bash
rectify run-all sample.fastq.gz --Scer -o results/ \
    --bam-dir /path/to/bams/ \      # Write all alignment BAMs here
    --keep-aligner-bams             # Retain per-aligner BAMs alongside rectified BAM
```

---

## HPC scheduler compatibility

RECTIFY's generated scripts detect the active scheduler automatically using
environment variables and normalise them into `RECTIFY_*` variables:

```bash
# Set at the top of every generated script:
RECTIFY_CPUS         # CPUs allocated to this task
RECTIFY_JOB_ID       # Job identifier
RECTIFY_TASK_ID      # Array task index (0-based)
RECTIFY_SCRATCH_BASE # Fast scratch directory (if available)
```

| Scheduler | CPU source | Job ID | Array task ID |
|-----------|-----------|--------|---------------|
| SLURM | `SLURM_CPUS_PER_TASK` | `SLURM_JOB_ID` | `SLURM_ARRAY_TASK_ID` |
| UGE/SGE | `NSLOTS` | `JOB_ID` | `SGE_TASK_ID` |
| PBS/Torque | `PBS_NUM_PPN` | `PBS_JOBID` | `PBS_ARRAY_INDEX` |

The generated scripts are SLURM `#SBATCH` scripts. For UGE or PBS clusters, use
the script body as a template and replace the directive block with the appropriate
`#$ -` or `#PBS` equivalents.

---

## Python path in generated scripts

Generated scripts activate the conda environment by prepending its `bin/` directory
to `PATH`. Avoid `conda activate` in batch scripts — it can silently fail in
non-interactive shells, leaving the job using system Python.

Recommended pattern:

```bash
# In your SLURM script:
export PATH="/path/to/conda/envs/myenv/bin:$PATH"
python -m rectify correct ...
```

---

## Chunked parallel alignment

For large DRS datasets (> 5 GB FASTQ, multi-hour alignment per aligner), use
`rectify split` to run alignment as a 2D SLURM array: N chunks × M aligners.

### Overview

```
rectify split    →  N chunk FASTQs
SLURM array      →  N × M alignment jobs (one per chunk/aligner pair)
samtools merge   →  one BAM per aligner (across all chunks)
rectify consensus →  best-aligner-per-read rectified BAM
```

### Step 1 — Split and generate scripts

```bash
rectify split reads.fastq.gz \
    -n 16 \
    -o /scratch/chunks/ \
    --generate-slurm \
    --aligners minimap2 mapPacBio gapmm2 uLTRA deSALT \
    --genome /ref/genome.fa.gz \
    --annotation /ref/genes.gff.gz \
    --slurm-partition my-partition \
    --slurm-cpus 16 \
    --slurm-mem 64G \
    --slurm-time 12:00:00
```

This writes to `/scratch/chunks/`:

| File | Purpose |
|------|---------|
| `wt_rep1_chunk_000_of_016.fastq.gz` … `_015_of_016.fastq.gz` | 16 equal chunk FASTQs |
| `chunks_manifest.json` | chunk paths in JSON |
| `run_array_align.sh` | SLURM `--array=0-79` script |
| `run_merge_and_consensus.sh` | post-array merge + consensus script |
| `slurm_logs/` | log directory |

### Step 2 — Submit the array

```bash
sbatch /scratch/chunks/run_array_align.sh
```

Each task decodes its `SLURM_ARRAY_TASK_ID` (0–79) as:

```bash
CHUNK_IDX=$(( SLURM_ARRAY_TASK_ID % 16 ))
ALIGNER_IDX=$(( SLURM_ARRAY_TASK_ID / 16 ))
```

Tasks run `rectify align --no-consensus` — each produces a single-aligner BAM
for that chunk. Thread limits (`OMP_NUM_THREADS`, `LOKY_MAX_CPU_COUNT`, etc.)
are set automatically by the generated script.

### Step 3 — Merge and run consensus

After the array completes:

```bash
bash /scratch/chunks/run_merge_and_consensus.sh
```

This script:

1. `samtools merge` — merges all chunk BAMs per aligner into one sorted BAM
2. `rectify consensus` — runs per-read aligner selection across all merged BAMs

### Scheduler environment variables

The generated scripts use `SLURM_ARRAY_TASK_ID` and `SLURM_CPUS_PER_TASK`.
For other schedulers, set these equivalents before invoking the array script:

| Scheduler | Task ID variable | CPU variable |
|-----------|-----------------|--------------|
| SLURM | `SLURM_ARRAY_TASK_ID` | `SLURM_CPUS_PER_TASK` |
| UGE/SGE | `SGE_TASK_ID` (1-based) | `NSLOTS` |
| PBS/Torque | `PBS_ARRAY_INDEX` | `PBS_NUM_PPN` |

For UGE/SGE, the task IDs are 1-based; adjust the decode logic:

```bash
# UGE wrapper (add before the CHUNK_IDX / ALIGNER_IDX lines):
SLURM_ARRAY_TASK_ID=$(( SGE_TASK_ID - 1 ))
SLURM_CPUS_PER_TASK=$NSLOTS
```

---

## Troubleshooting

**Job fails immediately with "command not found: rectify"**
: The conda environment is not on `PATH`. Set `PATH` explicitly (see above) and
  verify with `which rectify` in an interactive session first.

**Job exceeds memory limit (OOM kill)**
: Enable streaming mode (`streaming: true` in profile or `--streaming` on CLI)
  and increase `mem` to at least 8× the input BAM size.

**Many concurrent tasks all fail with I/O errors**
: NFS contention. Enable scratch staging (`use_scratch: true`) and verify that
  `$SCRATCH` (or equivalent) is defined in your job environment.

**Array task 0 succeeds; others fail**
: Check that array indices map correctly to samples. The generated script uses
  `$RECTIFY_TASK_ID` (0-based) to index into the sample list. Verify that
  `SLURM_ARRAY_TASK_ID` (or scheduler equivalent) is set and exported in your
  job environment.
