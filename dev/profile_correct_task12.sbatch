#!/bin/bash
#SBATCH --job-name=prof_correct_t12
#SBATCH --account=larsms
#SBATCH --partition=larsms
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=00:30:00
#SBATCH --output=prof_correct_t12_%j.out
#SBATCH --error=prof_correct_t12_%j.err
# TRAP 3: AVX-512 SIGILL on AMD Milan nodes — anaconda3 numpy uses AVX-512
#SBATCH --exclude=sh03-07n10,sh03-08n22

# TRAP 2 + 4
export PS1=""
source ~/.bashrc

# Activate conda BEFORE `set -u`: the rectify env's activate.d/java_home.sh hook
# references an unbound JAVA_HOME and dies under `set -u` (Sherlock trap variant).
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify

set -u
set -e
set -o pipefail

cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify

STAGE=/scratch/users/kevinroy/rectify_wt_by4742_rep1_25846844_0
GENOME=rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz
GFF=rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz
OUT=/scratch/users/kevinroy/rectify_perf_task12_$(date +%Y%m%d_%H%M%S)
mkdir -p $OUT
echo "OUTDIR=$OUT"

# All five per-aligner BAMs build the cross-aligner junction pool (same as correct).
ABAMS="$STAGE/wt_by4742_rep1.minimap2.bam $STAGE/wt_by4742_rep1.mapPacBio.bam $STAGE/wt_by4742_rep1.gapmm2.bam $STAGE/wt_by4742_rep1.deSALT.bam"

# Profile three aligners spanning the cost range (mm2 light, mapPacBio mid, gapmm2 heavy tail).
for ALN in minimap2 mapPacBio gapmm2; do
  echo "=================== profiling $ALN ==================="
  python dev/profile_correct_reads.py \
    --bam        $STAGE/wt_by4742_rep1.$ALN.bam \
    --genome     $GENOME \
    --annotation $GFF \
    --aligner-bams $ABAMS \
    --aligner-label $ALN \
    --max-reads  6000 \
    --junction-max-size 10000 \
    --out-prefix $OUT/profile_$ALN
done

echo "ALL DONE. Results in $OUT"
ls -la $OUT
