#!/bin/bash
#SBATCH --job-name=ds_verify
#SBATCH --partition=larsms,owners
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/scratch/users/kevinroy/rectify_dualsite_verify_20260524/%x_%j.out
#SBATCH --error=/scratch/users/kevinroy/rectify_dualsite_verify_20260524/%x_%j.err

# Single-thread full-pool `correct` on a dense chrI subset, the DIAGD repro.
# Captures wall time, peak RSS (/usr/bin/time -v), output line count + sha256.
# A/B: run once on baseline oak code (bd20f9e), sync M1->oak, run again
# (961c844). Only delta is the dual-site interval-tree pool index.
#
# Params via env (sbatch --export):
#   SUB     input subset BAM (built by the caller)
#   OUT_TSV corrected-reads output path
#   TMO     timeout seconds (default 1200)
#   LABEL   free-text tag for the log

set +u
export PS1=""
source ~/.bashrc
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
set -u
export MKL_DEBUG_CPU_TYPE=5

W=/scratch/users/kevinroy/rectify_wt_by4742_rep1_25846844_0
: "${TMO:=1200}"
RECTIFY_SRC=/oak/stanford/groups/larsms/Users/kevinroy/software/rectify

echo "host=$(hostname) start=$(date) LABEL=${LABEL:-?}"
echo "code marker dual-site: $(grep -c _POOL_FETCH_HALF_WINDOW $RECTIFY_SRC/rectify/core/bam/bam_processor.py) (0=baseline bisect, 1=dual-site)"
echo "subset reads: $(samtools view -c "$SUB")"

/usr/bin/time -v timeout -s ABRT "$TMO" python -X faulthandler -m rectify.cli correct \
    "$SUB" --Scer --skip-junction-refinement \
    --aligner-bams "$W/wt_by4742_rep1.minimap2.bam" \
    --aligner-bams "$W/wt_by4742_rep1.mapPacBio.bam" \
    --aligner-bams "$W/wt_by4742_rep1.gapmm2.bam" \
    --aligner-bams "$W/wt_by4742_rep1.uLTRA.bam" \
    --aligner-bams "$W/wt_by4742_rep1.deSALT.bam" \
    --legacy-single-threaded -j 1 \
    -o "$OUT_TSV" --verbose
RC=$?
echo "exit=$RC end=$(date)"
if [ -f "$OUT_TSV" ]; then
  echo "out_lines=$(wc -l < "$OUT_TSV")"
  echo "out_sha256=$(sha256sum "$OUT_TSV" | cut -d' ' -f1)"
fi
