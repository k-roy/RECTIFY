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

This generates two scheduler scripts in the output directory:

1. `results/rectify_batch_correct.sh` — array job (one task per sample)
2. `results/rectify_batch_analyze.sh` — combined analysis (runs after array)

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

By default, the per-aligner BAMs (one per aligner in the active panel) are
discarded after consensus selection; only the rectified BAM is kept. Use these
flags if you want to inspect the per-aligner outputs:

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

```text
rectify split     →  N chunk FASTQs + generated stage scripts
alignment arrays  →  per-chunk BAMs (mapPacBio array + "others" array)
merge + prescan   →  one BAM per aligner + variant/junction prescan
correction arrays →  per-aligner, per-chunk correction
chunk-merge + consensus →  best-aligner-per-read consensus BAM
```

### Step 1 — Split and generate scripts

```bash
rectify split reads.fastq.gz \
    -n 16 \
    -o /scratch/chunks/ \
    --generate-slurm \
    --other-aligners minimap2 gapmm2 uLTRA deSALT \
    --genome /ref/genome.fa.gz \
    --annotation /ref/genes.gff.gz \
    --slurm-partition my-partition \
    --slurm-account my-account
```

`--other-aligners` lists the aligners that share the "others" array; mapPacBio
always runs in its own array (omit it with `--skip-map-pacbio`). This writes to
`/scratch/chunks/`:

| File | Purpose |
|------|---------|
| `wt_rep1_chunk_000_of_016.fastq.gz` … `_015_of_016.fastq.gz` | 16 equal chunk FASTQs |
| `chunks_manifest.json` | chunk paths in JSON |
| `run_array_mapPacBio.sh` | mapPacBio alignment array (omitted with `--skip-map-pacbio`) |
| `run_array_others.sh` | minimap2/gapmm2/uLTRA/deSALT alignment array |
| `run_merge_aligners.sh`, `run_prescan.sh` | merge per-aligner BAMs; variant/junction prescan |
| `run_array_correct_<aligner>.sh` | per-aligner correction arrays |
| `run_array_chunk_merge.sh`, `run_final_merge.sh` | per-chunk merge → consensus; final merge |
| `submit_pipeline.sh` | submits the whole DAG with dependencies |
| `logs/` | log directory |

### Step 2 — Submit the pipeline

```bash
bash /scratch/chunks/submit_pipeline.sh
```

`submit_pipeline.sh` submits each generated stage with the correct inter-job
dependencies, following the correct-first ordering (align → merge → prescan →
correct → chunk-merge → consensus). The alignment array tasks run
`rectify align --no-consensus`; thread limits (`OMP_NUM_THREADS`,
`LOKY_MAX_CPU_COUNT`, etc.) are set automatically in each generated script. You
can also submit the individual array scripts by hand (e.g.
`sbatch run_array_others.sh`) for finer control.

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
