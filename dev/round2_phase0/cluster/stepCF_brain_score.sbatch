#!/bin/bash
#SBATCH --job-name=r2p0_brain_CF
#SBATCH --partition=owners
#SBATCH --constraint="CPU_GEN:SKX|CPU_GEN:SPR|CPU_GEN:GEN|CPU_GEN:BGM"
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --output=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0/brain_rna004/logs/stepCF.%j.out
# Round-2 Phase-0 brain RNA004 — Steps C (library) + F (cDNA align) + F.0 (lift round-trip) + E/F/G (score).
# Depends on Step D corrected BAMs in brain_rna004/corrected/.
# Idempotent skip-checks throughout.
set -eo pipefail
export PS1=""
source ~/.bashrc
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
set -u

P0=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0
BR=$P0/brain_rna004
RUN=$P0/round2
SCR=$RUN/scripts
RECTIFY=/oak/stanford/groups/larsms/Users/kevinroy/software/rectify
PENDIR=$RECTIFY/rectify/data/genomes/homo_sapiens/penalty_tables
export PYTHONPATH=$RECTIFY:$SCR:${PYTHONPATH:-}
T=$SLURM_CPUS_PER_TASK
REF=$BR/refs/chr5_6.chr.fa
GTF=$BR/refs/chr5_6.gencode.gtf
READS=$BR/reads/loci_reads.fastq
GENES="DIAPH1 TRIO MYO10 ABLIM3 KIF3A MPC1 COL23A1"

mkdir -p "$BR/library" "$BR/round2_aln" "$BR/verdict" "$BR/logs"
cd $BR
echo "[$(date)] Step CF brain start on $(hostname)"

# ---- C. build library: MANE ∪ Round-1 read-supported gate-passed chains ----
LIB=$BR/library
if [ ! -s "$LIB/library.fa" ]; then
  echo "[$(date)] building cDNA library (MANE ∪ de-novo gate-passed)"
  python $SCR/build_library.py --gtf "$GTF" --genome "$REF" --loci $BR/selected_loci.tsv \
    --genes $GENES --outdir "$LIB" --pad 1500 \
    --round1-bams $BR/corrected/minimap2_corrected.bam $BR/corrected/uLTRA_corrected.bam $BR/corrected/deSALT_corrected.bam \
    --min-anchor 10 --min-chain-support 3 2>&1 | tail -40
fi
echo "[$(date)] library cDNAs: $(grep -c '^>' "$LIB/library.fa")"

# ---- F. Round-2 alignment to padded cDNA library (minimap2 map-ont, no splice) ----
MM2C=$BR/round2_aln/cdna_minimap2.bam
if [ ! -s "$MM2C" ]; then
  echo "[$(date)] Round-2 minimap2 (map-ont, no splice, --for-only)"
  minimap2 -t $T -ax map-ont --for-only -N 100 -p 0.8 -Y "$LIB/library.fa" "$READS" 2>"$BR/logs/r2_mm2.log" \
    | samtools sort -@ $T -o "$MM2C" - && samtools index "$MM2C"
fi
echo "[$(date)] Round-2 alignment counts: $(samtools view -c -F 0x900 "$MM2C")"

# ---- F.0 lift-over round-trip on REAL reads (gate before trusting verdict) ----
echo "[$(date)] F.0 lift-over round-trip"
set +e
python $SCR/liftover_roundtrip.py --blockmaps "$LIB/blockmaps.json" --genome "$REF" \
  --bam $BR/aligners/minimap2.bam --per-strand 10 2>&1 | tee "$BR/logs/roundtrip.out"
RT=$?
set -e
if [ $RT -ne 0 ]; then
  echo "[$(date)] F.0 FAILED (exit $RT) — NOT scoring. Fix liftover.py."
  exit 1
fi

# ---- E + F + G score + verdict ----
CORARGS="--corrected-bam minimap2=$BR/corrected/minimap2_corrected.bam \
  --corrected-bam uLTRA=$BR/corrected/uLTRA_corrected.bam \
  --corrected-bam deSALT=$BR/corrected/deSALT_corrected.bam"
PENARGS="--penalty-scores $PENDIR/penalty_scores.tsv --str-penalty $PENDIR/str_penalty_scores.tsv"

echo "[$(date)] PRIMARY verdict (minimap2-cDNA; brain RNA004 substrate)"
python $SCR/phase0_score.py --blockmaps "$LIB/blockmaps.json" --genome "$REF" \
  $CORARGS --round2-bam mm2cdna=$MM2C $PENARGS \
  --min-anchor 10 --epsilon 1.0 --out "$BR/verdict" 2>&1 | tail -45

echo "[$(date)] Step CF brain DONE"
cat "$BR/verdict/verdict.txt" 2>/dev/null || true
