#!/bin/bash
#SBATCH --job-name=cigar_profiler
#SBATCH --partition=larsms,owners
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=/oak/stanford/groups/larsms/Users/kevinroy/common/scripts/nanopore/error_profile_20260422/slurm_%j.log

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LOKY_MAX_CPU_COUNT=$SLURM_CPUS_PER_TASK

PYTHON=/home/groups/larsms/users/kevinroy/anaconda3/bin/python
OUT=/oak/stanford/groups/larsms/Users/kevinroy/common/scripts/nanopore/error_profile_20260422
REF=/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa

echo "Python: $($PYTHON --version)"
echo "Start: $(date)"

$PYTHON \
    /oak/stanford/groups/larsms/Users/kevinroy/common/scripts/nanopore/empirical_cigar_error_profiler.py \
    --run-dir /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/dev_runs/wt_by4742_rep1_chunked_20260412 \
    --reference $REF \
    --output-dir $OUT \
    --isolation-flank 10 \
    --union \
    --str-repeat \
    2>&1 | tee $OUT/profiler.log

echo "End: $(date)"
