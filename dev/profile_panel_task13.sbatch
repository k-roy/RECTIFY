#!/bin/bash
#SBATCH --job-name=perf_panel_t13
#SBATCH --account=larsms
#SBATCH --partition=larsms
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --output=perf_panel_t13_%A_%a.out
#SBATCH --error=perf_panel_t13_%A_%a.err
#SBATCH --array=0-4
# TRAP 3: AVX-512 SIGILL on AMD Milan nodes + avoid Kevin's in-flight jobs on
# sh03-07n10 / sh03-08n22 (verification job 26167419 + Sumner array).
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

# Shared output dir keyed on the parent array job id so all 5 aligners land together.
OUT=/scratch/users/kevinroy/rectify_perf_panel_t13_${SLURM_ARRAY_JOB_ID}
mkdir -p $OUT

# All FIVE per-aligner BAMs build the cross-aligner junction pool (exactly as `correct`).
ABAMS="$STAGE/wt_by4742_rep1.minimap2.bam $STAGE/wt_by4742_rep1.mapPacBio.bam $STAGE/wt_by4742_rep1.gapmm2.bam $STAGE/wt_by4742_rep1.deSALT.bam $STAGE/wt_by4742_rep1.uLTRA.bam"

ALIGNERS=(minimap2 mapPacBio gapmm2 deSALT uLTRA)
ALN=${ALIGNERS[$SLURM_ARRAY_TASK_ID]}

# Genome-wide sweep: stride-5 (1 of every 5 primary reads spans all 16 chroms),
# capped at 120k profiled reads/aligner. Stream every read >= 500 ms to the panel
# TSV; accumulate whole-set aggregates instead of buffering all rows.
echo "=================== panel-profiling $ALN ==================="
echo "OUTDIR=$OUT  host=$(hostname)"
python dev/profile_correct_reads.py \
  --bam        $STAGE/wt_by4742_rep1.$ALN.bam \
  --genome     $GENOME \
  --annotation $GFF \
  --aligner-bams $ABAMS \
  --aligner-label $ALN \
  --sample-stride 5 \
  --max-reads  120000 \
  --cost-threshold-ms 500 \
  --junction-max-size 10000 \
  --summary-out $OUT/summary_$ALN.tsv \
  --out-prefix $OUT/profile_$ALN

echo "DONE $ALN. Panel + summary in $OUT"
ls -la $OUT/profile_$ALN.panel.tsv $OUT/summary_$ALN.tsv
