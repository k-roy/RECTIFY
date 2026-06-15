#!/bin/bash
#SBATCH --job-name=cdna_stepCF
#SBATCH --partition=owners
#SBATCH --constraint="CPU_GEN:SKX|CPU_GEN:SPR|CPU_GEN:GEN|CPU_GEN:BGM"
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0/cdna_a549/logs/stepCF.%j.out
# PCR-cDNA Phase-0 retest — Steps F (Round-2 cDNA align, sense-half) + F.0 (lift round-trip) +
# E/F/G (score + verdict). REUSES the directRNA run's 42-cDNA library + blockmaps unchanged
# (same loci, same genome) so only the read protocol differs. Depends on cDNA Step D.
set -eo pipefail
export PS1=""
source ~/.bashrc
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
set -u

P0=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0
RUN=$P0/cdna_a549
SRC=$P0/round2
SCR=$RUN/scripts
RECTIFY=/oak/stanford/groups/larsms/Users/kevinroy/software/rectify
PENDIR=$RECTIFY/rectify/data/genomes/homo_sapiens/penalty_tables
export PYTHONPATH=$RECTIFY:$SCR:${PYTHONPATH:-}
T=$SLURM_CPUS_PER_TASK
REF=$SRC/refs/chr5_6.fa
READS=$RUN/reads/loci_reads.fastq
LIB=$SRC/library          # REUSE directRNA library (42 cDNAs) unchanged

cd $RUN
echo "[$(date)] cDNA Step CF start on $(hostname)"
[ -s "$LIB/library.fa" ] || { echo "FATAL: reused library $LIB/library.fa missing"; exit 1; }
echo "[$(date)] reused library cDNAs: $(grep -c '^>' "$LIB/library.fa")"

# ---- F. Round-2 alignment to the padded cDNA library (sense-half via --for-only) ----
mkdir -p "$RUN/round2_aln"
MM2C=$RUN/round2_aln/cdna_minimap2.bam
if [ ! -s "$MM2C" ]; then
  echo "[$(date)] Round-2 minimap2 (map-ont, no splice, --for-only = RNA-sense half)"
  minimap2 -t $T -ax map-ont --for-only -N 100 -p 0.8 -Y "$LIB/library.fa" "$READS" 2>"$RUN/logs/r2_mm2.log" \
    | samtools sort -@ $T -o "$MM2C" - && samtools index "$MM2C"
fi
echo "[$(date)] Round-2 sense-half alignment count: $(samtools view -c -F 0x900 "$MM2C")"

# ---- F.0 lift-over round-trip on REAL reads (gate before any verdict) ----
echo "[$(date)] F.0 lift-over round-trip"
set +e
python $SCR/liftover_roundtrip.py --blockmaps "$LIB/blockmaps.json" --genome "$REF" \
  --bam $RUN/aligners/minimap2.bam --per-strand 10 2>&1 | tee "$RUN/logs/roundtrip.out"
RT=$?
set -e
if [ $RT -ne 0 ]; then echo "[$(date)] F.0 FAILED (exit $RT) — NOT scoring."; exit 1; fi

# ---- E + F + G score + verdict (patched harness: attribution requires a junction;
#      pre-registered STRONG junction-increasing-win rate is the decision metric) ----
CORARGS="--corrected-bam minimap2=$RUN/corrected/minimap2_corrected.bam \
  --corrected-bam uLTRA=$RUN/corrected/uLTRA_corrected.bam \
  --corrected-bam deSALT=$RUN/corrected/deSALT_corrected.bam"
PENARGS="--penalty-scores $PENDIR/penalty_scores.tsv --str-penalty $PENDIR/str_penalty_scores.tsv"

echo "[$(date)] PRIMARY verdict (minimap2-cDNA, sense-half)"
python $SCR/phase0_score.py --blockmaps "$LIB/blockmaps.json" --genome "$REF" \
  $CORARGS --round2-bam mm2cdna=$MM2C $PENARGS \
  --min-anchor 10 --epsilon 1.0 --out "$RUN/verdict_primary" 2>&1 | tail -45

echo "[$(date)] cDNA Step CF DONE"
cat "$RUN/verdict_primary/verdict.txt"
