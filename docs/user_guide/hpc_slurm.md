# HPC / SLURM

RECTIFY has built-in support for SLURM HPC environments. This page covers the bundled profiles, job submission, and best practices for running on Stanford's Sherlock cluster or any compatible system.

---

## Quick submission

```bash
rectify run \
    --manifest manifest.tsv \
    --Scer \
    --reference wt \
    -o results/ \
    --profile /path/to/rectify/slurm_profiles/sherlock_larsms.yaml
```

This generates and optionally submits two SLURM scripts:

1. `results/slurm/rectify_batch_correct.sh` — SLURM array job (one task per sample)
2. `results/slurm/rectify_batch_analyze.sh` — combined analysis (submitted after array completes)

---

## Bundled profiles

Two profiles are included in `rectify/slurm_profiles/`:

### `sherlock_larsms.yaml` (recommended)

```yaml
partition: larsms,owners
time: "8:00:00"
mem: 128G
cpus: 8
max_concurrent: 20
use_scratch: true   # Stage BAMs through $SCRATCH (75 GB/s)
streaming: true     # Process reads in chunks — peak RSS ~4-5 GB
```

### `sherlock_gpu.yaml`

For GPU-accelerated alignment (e.g. with dorado basecalling):

```yaml
partition: gpu
time: "4:00:00"
mem: 64G
cpus: 8
gpus: 1
```

---

## $SCRATCH staging

**Never run BAM correction I/O directly on Oak in a SLURM array job.**

Oak is NFS-backed shared storage. Running 20 concurrent array tasks each streaming a 7 GB BAM through Oak causes severe I/O contention and inflates wall time 2–3×.

When `use_scratch: true` is set, the generated SLURM script:
1. Copies the input BAM to `$SCRATCH` (node-local, ~75 GB/s)
2. Runs `rectify correct` on the local copy
3. Rsyncs all outputs back to Oak
4. Cleans up `$SCRATCH`

!!! note
    FASTQ inputs are left on Oak — they are read sequentially once and don't benefit from staging.

!!! warning
    `$SCRATCH` on Sherlock is auto-purged after 90 days. Always rsync outputs back to Oak before the purge window.

---

## Thread limits

**SLURM will suspend or ban accounts that spawn more processes than allocated CPUs.**

Python packages (numpy, sklearn, pydeseq2) auto-spawn threads via OpenMP, OpenBLAS, MKL, and joblib's loky backend. RECTIFY's generated SLURM scripts set all thread limits automatically:

```bash
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LOKY_MAX_CPU_COUNT=$SLURM_CPUS_PER_TASK   # critical for pydeseq2/sklearn
```

!!! warning "LOKY_MAX_CPU_COUNT is mandatory"
    Joblib's loky backend ignores `JOBLIB_WORKERS`. You **must** set `LOKY_MAX_CPU_COUNT` or pydeseq2 and sklearn will over-subscribe.

If writing your own scripts, set these **before** importing numpy:

```python
import os
os.environ['OMP_NUM_THREADS'] = '8'
os.environ['LOKY_MAX_CPU_COUNT'] = '8'
# ... now safe to import numpy, pandas, sklearn
import numpy as np
```

---

## Streaming mode

For large BAMs (> 2 GB), use streaming mode to keep peak RAM at ~4–5 GB regardless of file size:

```bash
rectify correct reads.bam --Scer --streaming -o results/
```

The `sherlock_larsms.yaml` profile sets `streaming: true` by default.

Without streaming, the parallel accumulation path loads all results into RAM before writing. For a 7 GB BAM with 40M reads, this requires ~30–40 GB.

---

## Interactive sessions

For quick tests before submitting batch jobs:

```bash
# Request an interactive node (Sherlock)
srun --partition=larsms --account=larsms \
     --cpus-per-task=8 --time=2:00:00 --mem=32G --pty bash

# Set thread limits
export OMP_NUM_THREADS=8
export LOKY_MAX_CPU_COUNT=8

# Test with a small BAM
rectify correct test.bam --Scer --streaming -o test_output/
```

---

## Checking job status

```bash
squeue -u $USER            # Active jobs
sacct -j <jobid>           # Completed job stats
sstat -j <jobid>           # Running job resource usage
sshare -u $USER            # Fairshare score (higher = better priority)
```

---

## Python path in SLURM scripts

Generated scripts use an explicit Python path — never `conda activate`:

```bash
PYTHON="/home/groups/larsms/users/kevinroy/anaconda3/bin/python"
$PYTHON -m rectify correct ...
```

This avoids the common failure where `conda activate` silently falls back to system Python 3.6.

---

## Dry-run validation

Before submitting to a slow partition (e.g., GPU), validate your environment on a fast partition:

```bash
#!/bin/bash
#SBATCH --partition=larsms
#SBATCH --time=0:05:00
#SBATCH --mem=4G

PYTHON="/home/groups/larsms/users/kevinroy/anaconda3/bin/python"

$PYTHON -c "
import rectify
print(f'RECTIFY {rectify.__version__} OK')
import torch; print(f'PyTorch {torch.__version__}, CUDA: {torch.cuda.is_available()}')
"
$PYTHON -m py_compile $(which rectify) && echo 'Syntax OK'
```

Submit this first; if it passes, submit your real job.
