#!/bin/bash
# =============================================================================
# LongBench RNA004 SIRV spike-in error-structure  (the modern-chemistry arm)
# -----------------------------------------------------------------------------
# Completes the simulator/real three-way with the CURRENT ONT chemistry (RNA004).
# Download the smallest LongBench DRS cell line (H69, ~11.3 GB), align ALL reads
# to the staged SIRV-Set 4 reference (endogenous human won't map → mapped = the
# spike-in subset, absolute truth), and measure the error structure with the SAME
# exonic measure-bam used for yeast / LRGASP-SIRV / SG-NEx (apples-to-apples).
# Spike-ins confirmed present (~0.30% = ~81k reads; SESSION-8). Ref = reuse local
# sirv4.fasta (same SIRV-Set4 the other arms used). RNA004 fastq SEQ is already T.
#
# HEADER: confirm via sherlock-sbatch skill if needed; larsms node, no AVX exclude.
# =============================================================================
#SBATCH --job-name=longbench_h69
#SBATCH --account=larsms
#SBATCH --partition=larsms
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=/scratch/users/kevinroy/longbench_work/lb_h69.%j.out
#SBATCH --error=/scratch/users/kevinroy/longbench_work/lb_h69.%j.err

set -eo pipefail

D=/home/groups/larsms/users/kevinroy/aligner_bench_live
CB=/home/groups/larsms/users/kevinroy/anaconda3
P3=$CB/envs/pbsim3
PY=$CB/envs/rectify/bin/python
W=/scratch/users/kevinroy/longbench_work
SIRV=/scratch/users/kevinroy/sirv_work/sirv4.fasta      # staged SIRV-Set4 ref (reuse)
mkdir -p "$W"; cd "$W"
MM2=$P3/bin/minimap2; SAMTOOLS=$P3/bin/samtools
source "$CB/etc/profile.d/conda.sh"; conda activate "$CB/envs/rectify"
set -u
THREADS=${SLURM_CPUS_PER_TASK:-16}

# ---- 1. download H69 RNA004 DRS reads (idempotent .part -> mv) --------------
URL=https://longbench-data.s3.ap-southeast-2.amazonaws.com/raw/fastq/H69_dRNA_ONT.fastq.gz
if [[ ! -s H69_dRNA_ONT.fastq.gz ]]; then
  echo "[lb] downloading H69 (~11.3 GB) ..."
  wget -q -O H69_dRNA_ONT.fastq.gz.part "$URL"
  mv H69_dRNA_ONT.fastq.gz.part H69_dRNA_ONT.fastq.gz
fi
echo "[lb] reads = $(stat -c%s H69_dRNA_ONT.fastq.gz) bytes"

# ---- 2. align ALL reads to the SIRV ref (mapped = spike-ins) ----------------
# DRS flags + --MD (events_from_bam_read needs MD); U->T fixup kept (harmless on
# already-T RNA004). Endogenous human reads won't map to the 430 KB SIRV ref.
if [[ ! -s h69_sirv.sorted.bam ]]; then
  zcat H69_dRNA_ONT.fastq.gz \
    | awk 'NR%4==2{gsub(/U/,"T");gsub(/u/,"t")}1' \
    | "$MM2" -ax splice -uf -k14 --MD --secondary=no -t "$THREADS" "$SIRV" - \
    | "$SAMTOOLS" sort -@ "$THREADS" -o h69_sirv.sorted.bam -
  "$SAMTOOLS" index h69_sirv.sorted.bam
fi
MAPPED=$("$SAMTOOLS" view -c -F 0x904 h69_sirv.sorted.bam)
echo "[lb] primary-mapped (spike-in) reads = $MAPPED"

# ---- 3. measure (exonic; thin to yeast rate; [600,1000]) --------------------
VAL="$D/scripts/benchmark/error_realism_validate.py"
REPORT="$W/longbench_h69_errstruct.report.txt"
set +e
{
  echo "==== LongBench H69 RNA004 SIRV error structure ===="
  echo "mapped(primary)=$MAPPED  ref=sirv4.fasta (SIRV-Set4)"
  echo "---- measure-bam (exonic, thin 0.0153, [600,1000]) ----"
  PYTHONPATH=$D "$PY" "$VAL" measure-bam h69_sirv.sorted.bam --thin 0.0153 --len-min 600 --len-max 1000
  echo "---- autocorr-bam ----"
  PYTHONPATH=$D "$PY" "$VAL" autocorr-bam h69_sirv.sorted.bam --len-min 600 --len-max 1000 --frac 0.3
} 2>&1 | tee "$REPORT"
RC=${PIPESTATUS[0]}
printf '%s\n' "$RC" > "$W/.lb_rc.tmp" && mv "$W/.lb_rc.tmp" "$W/.lb_rc"
echo "[lb] DONE rc=$RC report=$REPORT"
exit $RC
