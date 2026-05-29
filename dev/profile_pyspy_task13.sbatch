#!/bin/bash
#SBATCH --job-name=pyspy_t13
#SBATCH --account=larsms
#SBATCH --partition=larsms
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --time=00:40:00
#SBATCH --output=pyspy_t13_%j.out
#SBATCH --error=pyspy_t13_%j.err
# NOTE: larsms has only 2 nodes (sh03-07n10/sh03-08n22); --exclude would leave
# nowhere to run, so it is intentionally omitted. AVX-512 SIGILL is handled by the
# env's MKL_DEBUG_CPU_TYPE=5 (verified: panel jobs run past genome load w/o SIGILL).

export PS1=""
source ~/.bashrc
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
set -u
# NOTE: deliberately NO `set -e` / `set -o pipefail` here. The RID extraction below
# uses `... | sort | head -1`, where head closing the pipe early sends SIGPIPE to sort
# → pipefail would propagate exit 141 → set -e would kill the whole script (this is
# exactly what failed job 26181308 in 4 s). py-spy attach + kill are also expected to
# return nonzero. Each step is guarded explicitly instead.

cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify
STAGE=/scratch/users/kevinroy/rectify_wt_by4742_rep1_25846844_0
GENOME=rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz
GFF=rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz
ABAMS="$STAGE/wt_by4742_rep1.minimap2.bam $STAGE/wt_by4742_rep1.mapPacBio.bam $STAGE/wt_by4742_rep1.gapmm2.bam $STAGE/wt_by4742_rep1.deSALT.bam $STAGE/wt_by4742_rep1.uLTRA.bam"

# PANEL_DIR is the task-13 panel output (passed via --export or edited here).
PANEL_DIR=${PANEL_DIR:-/scratch/users/kevinroy/rectify_perf_panel_t13_PARENTJOBID}
OUT=/scratch/users/kevinroy/rectify_perf_pyspy_t13_$(date +%Y%m%d_%H%M%S)
mkdir -p $OUT
echo "OUTDIR=$OUT  PANEL_DIR=$PANEL_DIR"

# py-spy the single worst read of each aligner (top elapsed_ms row in its panel TSV).
for ALN in minimap2 mapPacBio gapmm2 deSALT uLTRA; do
  PF=$PANEL_DIR/profile_$ALN.panel.tsv
  [ -s "$PF" ] || { echo "no panel for $ALN, skipping"; continue; }
  # column 1 = read_id, column 2 = elapsed_ms; take the max-elapsed row.
  RID=$(tail -n +2 "$PF" | sort -t$'\t' -k2,2 -rn | head -1 | cut -f1)
  MS=$(tail -n +2 "$PF" | sort -t$'\t' -k2,2 -rn | head -1 | cut -f2)
  echo "=============== py-spy on $ALN worst read $RID (${MS} ms) ==============="
  python dev/profile_correct_reads.py \
    --bam $STAGE/wt_by4742_rep1.$ALN.bam \
    --genome $GENOME --annotation $GFF --aligner-bams $ABAMS \
    --aligner-label $ALN --pyspy-read $RID --pyspy-repeat 100000 \
    --out-prefix $OUT/ignore &
  PYPID=$!
  sleep 22
  py-spy record -f raw -o $OUT/stacks_$ALN.folded --pid $PYPID --duration 25 --rate 200 --nonblocking 2>$OUT/pyspy_$ALN.log || true
  kill $PYPID 2>/dev/null || true
  wait $PYPID 2>/dev/null || true
  echo "--- top leaf frames $ALN ---"
  awk '{c=$NF; s=$0; sub(/ [0-9]+$/,"",s); n=split(s,a,";"); sum[a[n]]+=c} END{for(k in sum) print sum[k]"\t"k}' $OUT/stacks_$ALN.folded 2>/dev/null | sort -rn | head -12
done

echo "ALL DONE. Results in $OUT"
ls -la $OUT
