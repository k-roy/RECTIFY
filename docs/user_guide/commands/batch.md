# rectify batch

Generate SLURM array job scripts for multi-sample processing on HPC clusters.

This command creates ready-to-submit SLURM scripts — one array job for per-sample correction, one follow-up job for combined analysis. It is also called internally by `rectify run --profile`.

---

## Usage

```bash
rectify batch [options] -o <output_dir>
```

## Examples

```bash
# Generate SLURM scripts (does not submit)
rectify batch \
    --manifest manifest.tsv \
    --Scer \
    --reference wt \
    -o results/ \
    --profile slurm_profiles/sherlock_larsms.yaml

# Generate and immediately submit
rectify batch \
    --manifest manifest.tsv \
    --Scer \
    -o results/ \
    --profile slurm_profiles/sherlock_larsms.yaml \
    --submit
```

---

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--manifest, -m` | — | Sample manifest TSV |
| `-o, --output-dir` | — | Output directory |
| `--profile` | — | SLURM profile YAML |
| `--Scer` | — | Bundled *S. cerevisiae* data |
| `--genome` | — | Reference genome FASTA |
| `--annotation` | — | Gene annotation GFF/GTF |
| `--reference` | auto | Reference condition for DESeq2 |
| `--submit` | off | Submit generated scripts immediately |

---

## Generated scripts

### `results/slurm/rectify_batch_correct.sh`

A SLURM array job with one task per sample:

```bash
#!/bin/bash
#SBATCH --job-name=rectify_correct
#SBATCH --array=0-3             # 0 to N-1 samples
#SBATCH --partition=larsms,owners
#SBATCH --time=8:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=results/slurm/correct_%A_%a.log

# Thread limits (critical for SLURM compliance)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LOKY_MAX_CPU_COUNT=$SLURM_CPUS_PER_TASK

# Stage BAM to $SCRATCH (if use_scratch: true in profile)
SCRATCH_DIR=$SCRATCH/rectify_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $SCRATCH_DIR
cp /data/sample_${SLURM_ARRAY_TASK_ID}.bam $SCRATCH_DIR/

# Run correction
rectify correct $SCRATCH_DIR/sample.bam --Scer --streaming \
    -o results/sample_${SLURM_ARRAY_TASK_ID}/

# Sync outputs back to Oak
rsync -a results/sample_${SLURM_ARRAY_TASK_ID}/ \
    /oak/results/sample_${SLURM_ARRAY_TASK_ID}/
rm -rf $SCRATCH_DIR
```

### `results/slurm/rectify_batch_analyze.sh`

Submitted with `--dependency=afterok:<array_job_id>`:

```bash
#!/bin/bash
#SBATCH --job-name=rectify_analyze
#SBATCH --partition=larsms,owners
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

rectify analyze /dev/null \
    --manifest manifest.tsv \
    --Scer \
    --reference wt \
    --run-deseq2 \
    --run-motif \
    -o results/combined/
```

---

## SLURM profile format

```yaml
partition: larsms,owners   # SLURM partition(s)
time: "8:00:00"            # Wall time
mem: 128G                  # Memory per task
cpus: 8                    # CPUs per task
max_concurrent: 20         # Max simultaneous array tasks
use_scratch: true          # Stage BAMs through $SCRATCH
streaming: true            # Use streaming mode in rectify correct
```

See [HPC / SLURM](../hpc_slurm.md) for detailed best practices.

---

## Manual submission

```bash
# Submit array job
ARRAY_JOB=$(sbatch --parsable results/slurm/rectify_batch_correct.sh)

# Submit analysis after array completes
sbatch --dependency=afterok:$ARRAY_JOB results/slurm/rectify_batch_analyze.sh
```

---

## Monitoring

```bash
squeue -u $USER                  # Check running jobs
tail -f results/slurm/correct_*.log   # Live log tail
```
