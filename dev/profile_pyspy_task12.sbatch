#!/bin/bash
#SBATCH --job-name=pyspy_t12
#SBATCH --account=larsms
#SBATCH --partition=larsms
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --time=00:30:00
#SBATCH --output=pyspy_t12_%j.out
#SBATCH --error=pyspy_t12_%j.err
#SBATCH --exclude=sh03-07n10,sh03-08n22

export PS1=""
source ~/.bashrc
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
set -u
set -e
set -o pipefail

cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify
STAGE=/scratch/users/kevinroy/rectify_wt_by4742_rep1_25846844_0
GENOME=rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz
GFF=rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz
ABAMS="$STAGE/wt_by4742_rep1.minimap2.bam $STAGE/wt_by4742_rep1.mapPacBio.bam $STAGE/wt_by4742_rep1.gapmm2.bam $STAGE/wt_by4742_rep1.deSALT.bam"
OUT=/scratch/users/kevinroy/rectify_perf_task12_pyspy_$(date +%Y%m%d_%H%M%S)
mkdir -p $OUT
echo "OUTDIR=$OUT"

# --- Part A: py-spy attribution on the single worst read per aligner ---
declare -A WORST=(
  [minimap2]=3706b831-4647-426b-9199-2e9c77e9913d
  [mapPacBio]=b6d5969a-cfac-464c-898c-1c63c22c13f9
  [gapmm2]=a867e7d2-9d58-48a3-b739-e99e3253aa78
)
for ALN in minimap2 mapPacBio gapmm2; do
  RID=${WORST[$ALN]}
  echo "=============== py-spy on $ALN worst read $RID ==============="
  # High repeat ceiling so the loop is still running when py-spy samples; we KILL
  # after recording rather than waiting for all repeats (a single gapmm2 call is
  # ~12 s, so 400 repeats would blow the wall clock — py-spy only needs ~25 s of
  # samples regardless of how many repeats remain).
  python dev/profile_correct_reads.py \
    --bam $STAGE/wt_by4742_rep1.$ALN.bam \
    --genome $GENOME --annotation $GFF --aligner-bams $ABAMS \
    --aligner-label $ALN --pyspy-read $RID --pyspy-repeat 100000 \
    --out-prefix $OUT/ignore &
  PYPID=$!
  # wait for setup (genome load + pool build ~20s) then attach py-spy to the repeat loop
  sleep 22
  py-spy record -f raw -o $OUT/stacks_$ALN.folded --pid $PYPID --duration 25 --rate 200 --nonblocking 2>$OUT/pyspy_$ALN.log || true
  kill $PYPID 2>/dev/null || true
  wait $PYPID 2>/dev/null || true
  echo "--- top leaf frames $ALN ---"
  awk '{c=$NF; s=$0; sub(/ [0-9]+$/,"",s); n=split(s,a,";"); sum[a[n]]+=c} END{for(k in sum) print sum[k]"\t"k}' $OUT/stacks_$ALN.folded 2>/dev/null | sort -rn | head -12
done

# --- Part B: cross-chromosome sampled profile (stride to span all chroms) ---
# Tests whether the "cost no longer candidate-driven" finding holds beyond chrI.
echo "=============== cross-chrom stride profile (gapmm2) ==============="
python dev/profile_correct_reads.py \
  --bam $STAGE/wt_by4742_rep1.gapmm2.bam \
  --genome $GENOME --annotation $GFF --aligner-bams $ABAMS \
  --aligner-label gapmm2_xchrom --max-reads 8000 --sample-stride 40 \
  --junction-max-size 10000 --out-prefix $OUT/profile_gapmm2_xchrom

echo "=============== cross-chrom stride profile (minimap2) ==============="
python dev/profile_correct_reads.py \
  --bam $STAGE/wt_by4742_rep1.minimap2.bam \
  --genome $GENOME --annotation $GFF --aligner-bams $ABAMS \
  --aligner-label minimap2_xchrom --max-reads 8000 --sample-stride 40 \
  --junction-max-size 10000 --out-prefix $OUT/profile_minimap2_xchrom

echo "ALL DONE. Results in $OUT"
ls -la $OUT
