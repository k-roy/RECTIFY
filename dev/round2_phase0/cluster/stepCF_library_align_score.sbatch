#!/bin/bash
#SBATCH --job-name=r2p0_stepCF
#SBATCH --partition=owners
#SBATCH --constraint="CPU_GEN:SKX|CPU_GEN:SPR|CPU_GEN:GEN|CPU_GEN:BGM"
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --output=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0/round2/logs/stepCF.%j.out
# Round-2 cDNA Phase-0 — Steps C (library) + F (cDNA align) + F.0 (lift round-trip) + E/F/G (score).
# Depends on Step D corrected genome BAMs. Idempotent skip-checks throughout.
set -eo pipefail
export PS1=""
source ~/.bashrc
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
set -u

P0=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0
RUN=$P0/round2
SCR=$RUN/scripts
RECTIFY=/oak/stanford/groups/larsms/Users/kevinroy/software/rectify
PENDIR=$RECTIFY/rectify/data/genomes/homo_sapiens/penalty_tables
export PYTHONPATH=$RECTIFY:$SCR:${PYTHONPATH:-}
T=$SLURM_CPUS_PER_TASK
REF=$RUN/refs/chr5_6.fa
GTF=$RUN/refs/chr5_6.gtf
READS=$RUN/reads/loci_reads.fastq
GENES="TGFBI IK HSD17B4 ANXA6 TRIO XPO5 DIAPH1"
cd $RUN
echo "[$(date)] Step CF start on $(hostname)"

# ---- C. build library: MANE ∪ Round-1 read-supported gate-passed chains ----
LIB=$RUN/library
if [ ! -s "$LIB/library.fa" ]; then
  echo "[$(date)] building cDNA library (MANE ∪ de-novo gate-passed)"
  python $SCR/build_library.py --gtf "$GTF" --genome "$REF" --loci $P0/coverage/selected_loci.tsv \
    --genes $GENES --outdir "$LIB" --pad 1500 \
    --round1-bams $RUN/corrected/minimap2_corrected.bam $RUN/corrected/uLTRA_corrected.bam $RUN/corrected/deSALT_corrected.bam \
    --min-anchor 10 --min-chain-support 3 2>&1 | tail -40
fi
echo "[$(date)] library cDNAs: $(grep -c '^>' "$LIB/library.fa")"

# ---- F. Round-2 alignment of ALL loci reads to the padded cDNA library ----
mkdir -p "$RUN/round2_aln"
MM2C=$RUN/round2_aln/cdna_minimap2.bam
MPBC=$RUN/round2_aln/cdna_mappacbio.bam
if [ ! -s "$MM2C" ]; then
  echo "[$(date)] Round-2 minimap2 (map-ont, no splice, --for-only)"
  minimap2 -t $T -ax map-ont --for-only -N 100 -p 0.8 -Y "$LIB/library.fa" "$READS" 2>"$RUN/logs/r2_mm2.log" \
    | samtools sort -@ $T -o "$MM2C" - && samtools index "$MM2C"
fi
# mapPacBio (BBMapPacBio) DROPPED for human DRS: hard read-length cap ~6019 bp throws
# AssertionError on every full-length transcript read >6 kb (job 29566523), and it games
# HP-ED anyway (feedback_hped_overcall_gaming). minimap2 map-ont is the sole Round-2 aligner.
echo "[$(date)] Round-2 alignment counts:"
for b in "$MM2C" "$MPBC"; do [ -s "$b" ] && echo "  $b : $(samtools view -c -F 0x900 "$b")"; done

# ---- F.0 lift-over round-trip on REAL reads (gate before any verdict) ----
echo "[$(date)] F.0 lift-over round-trip"
set +e
python $SCR/liftover_roundtrip.py --blockmaps "$LIB/blockmaps.json" --genome "$REF" \
  --bam $RUN/aligners/minimap2.bam --per-strand 10 2>&1 | tee "$RUN/logs/roundtrip.out"
RT=$?
set -e
if [ $RT -ne 0 ]; then
  echo "[$(date)] F.0 FAILED (exit $RT) — NOT scoring. Fix liftover.py."
  exit 1
fi

# ---- E + F + G score + verdict ----
CORARGS="--corrected-bam minimap2=$RUN/corrected/minimap2_corrected.bam \
  --corrected-bam uLTRA=$RUN/corrected/uLTRA_corrected.bam \
  --corrected-bam deSALT=$RUN/corrected/deSALT_corrected.bam"
PENARGS="--penalty-scores $PENDIR/penalty_scores.tsv --str-penalty $PENDIR/str_penalty_scores.tsv"

echo "[$(date)] PRIMARY verdict (minimap2-cDNA only; spec PRIMARY aligner)"
python $SCR/phase0_score.py --blockmaps "$LIB/blockmaps.json" --genome "$REF" \
  $CORARGS --round2-bam mm2cdna=$MM2C $PENARGS \
  --min-anchor 10 --epsilon 1.0 --out "$RUN/verdict_primary" 2>&1 | tail -45

if [ -s "$MPBC" ]; then
  echo "[$(date)] SECONDARY verdict (minimap2 + mapPacBio cDNA; mapPacBio games HP-ED — sensitivity check)"
  python $SCR/phase0_score.py --blockmaps "$LIB/blockmaps.json" --genome "$REF" \
    $CORARGS --round2-bam mm2cdna=$MM2C --round2-bam mpbcdna=$MPBC $PENARGS \
    --min-anchor 10 --epsilon 1.0 --out "$RUN/verdict_with_mappacbio" 2>&1 | tail -20
fi

echo "[$(date)] Step CF DONE"
