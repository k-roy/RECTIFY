#!/bin/bash
#SBATCH --job-name=r2p0_brain_D
#SBATCH --partition=owners
#SBATCH --constraint="CPU_GEN:SKX|CPU_GEN:SPR|CPU_GEN:GEN|CPU_GEN:BGM"
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0/brain_rna004/logs/stepD.%j.out
# Round-2 Phase-0 brain RNA004 — Step D: Round-1 genome alignment (3-aligner) + rectify correct.
# Substrate: ENCODE ENCFF043YFO dlpfc_ctrl85 (GENCODE GRCh38, chr-naming, RNA004).
# Refs already built in brain_rna004/refs/ (chr5_6.chr.fa + chr5_6.gencode.gtf).
# Idempotent: every stage skip-checks its output.
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
T=$SLURM_CPUS_PER_TASK

REF=$BR/refs/chr5_6.chr.fa
GTF=$BR/refs/chr5_6.gencode.gtf
POOL=$BR/data/brain_chr5_6.bam
READS=$BR/reads/loci_reads.fastq

mkdir -p "$BR/reads" "$BR/aligners" "$BR/corrected" "$BR/logs"
cd $BR
echo "[$(date)] Step D brain start on $(hostname); cpus=$T"
echo "[$(date)] ref: $(grep -c '^>' $REF) seqs; GTF lines: $(wc -l < $GTF)"

# ---- 1. extract reads at the 7 loci -> fastq (dedup by read id) ----
if [ ! -s "$READS" ]; then
  echo "[$(date)] extracting reads at loci"
  REGS=$(awk 'NR>0{printf "%s:%s-%s ", $2,$3,$4}' "$BR/loci.tsv")
  echo "  regions: $REGS"
  samtools view -b -F 0x900 "$POOL" $REGS \
    | samtools collate -u -O - \
    | samtools fastq - 2>/dev/null \
    | awk 'NR%4==1{id=$1; if(seen[id]++){skip=1}else{skip=0}} skip==0{print}' > "$READS"
  echo "  reads: $(( $(wc -l < "$READS") / 4 ))"
fi

# ---- 2. uLTRA index (chr-prefixed ref + GENCODE GTF) ----
ULTRA_IDX=$BR/refs/ultra_index
if [ ! -s "$ULTRA_IDX/refs_sequences.fa" ]; then
  echo "[$(date)] building uLTRA index"
  rm -rf "$ULTRA_IDX"; mkdir -p "$ULTRA_IDX"
  uLTRA index "$REF" "$GTF" "$ULTRA_IDX" 2>&1 | tail -5
fi

# ---- 3. deSALT index ----
DESALT_IDX=$BR/refs/desalt_index
if [ ! -d "$DESALT_IDX" ] || [ -z "$(ls -A "$DESALT_IDX" 2>/dev/null)" ]; then
  echo "[$(date)] building deSALT index"
  mkdir -p "$DESALT_IDX"
  deSALT index "$REF" "$DESALT_IDX" 2>&1 | tail -3
fi

# ---- 4. alignments ----
MM2=$BR/aligners/minimap2.bam
ULT=$BR/aligners/uLTRA.bam
DES=$BR/aligners/deSALT.bam

if [ ! -s "$MM2" ]; then
  echo "[$(date)] minimap2"
  minimap2 -t $T -ax splice -uf -k14 --max-intron 500000 "$REF" "$READS" 2>"$BR/logs/mm2.log" \
    | samtools sort -@ $T -o "$MM2" - && samtools index "$MM2"
fi

if [ ! -s "$ULT" ]; then
  echo "[$(date)] uLTRA align"
  UOUT=$BR/aligners/ultra_out; mkdir -p "$UOUT"
  uLTRA align --ont --t $T --index "$ULTRA_IDX" --prefix loci "$REF" "$READS" "$UOUT" 2>&1 | tail -5
  samtools sort -@ $T -o "$ULT" "$UOUT/loci.sam" && samtools index "$ULT"
fi

if [ ! -s "$DES" ]; then
  echo "[$(date)] deSALT aln"
  deSALT aln -t $T -x ont1d -T -f "$BR/aligners/desalt.tmp" -o "$BR/aligners/deSALT.sam" "$DESALT_IDX" "$READS" 2>"$BR/logs/desalt.log"
  samtools sort -@ $T -o "$DES" "$BR/aligners/deSALT.sam" && samtools index "$DES"
  rm -f "$BR/aligners/deSALT.sam"
fi

echo "[$(date)] alignment counts (primary):"
for b in "$MM2" "$ULT" "$DES"; do echo "  $b : $(samtools view -c -F 0x900 "$b")"; done

# ---- 5. per-aligner rectify correct (DRS) ----
RECTIFY=/oak/stanford/groups/larsms/Users/kevinroy/software/rectify
export PYTHONPATH=$RECTIFY:$SCR:${PYTHONPATH:-}
for name in minimap2 uLTRA deSALT; do
  IN=$BR/aligners/$name.bam
  OUTTSV=$BR/corrected/${name}_corrected.tsv
  OUTBAM=$BR/corrected/${name}_corrected.bam
  if [ -s "$OUTBAM" ]; then echo "[$(date)] skip correct $name (done)"; continue; fi
  echo "[$(date)] rectify correct $name"
  rectify correct "$IN" --genome "$REF" -o "$OUTTSV" \
    --organism homo_sapiens --write-corrected-bam "$OUTBAM" \
    -j $T 2>"$BR/logs/correct_${name}.log"
  samtools index "$OUTBAM" 2>/dev/null || true
done

echo "[$(date)] Step D brain DONE"
ls -la $BR/corrected/
