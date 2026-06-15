#!/bin/bash
#SBATCH --job-name=r2p0_stepD
#SBATCH --partition=owners
#SBATCH --constraint="CPU_GEN:SKX|CPU_GEN:SPR|CPU_GEN:GEN|CPU_GEN:BGM"
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0/round2/logs/stepD.%j.out
# Round-2 cDNA Phase-0 — Step D: Round-1 genome alignment (3-aligner panel) + per-aligner rectify correct.
# Idempotent: every stage skip-checks its output so a preempted/requeued job resumes cheaply.
set -eo pipefail
export PS1=""
source ~/.bashrc
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
set -u

P0=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0
RUN=$P0/round2
GENO=/oak/stanford/groups/larsms/Projects/split_tap/newvolume/genomes
SRCFA=$GENO/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
SRCGTF=$GENO/Homo_sapiens.GRCh38.109.gtf.gz
POOL=$P0/data/a549_pooled_chr5_6.bam
T=$SLURM_CPUS_PER_TASK

REF=$RUN/refs/chr5_6.fa
GTF=$RUN/refs/chr5_6.gtf
READS=$RUN/reads/loci_reads.fastq

cd $RUN
echo "[$(date)] Step D start on $(hostname); cpus=$T"

# ---- 1. chr5+6 Ensembl reference fasta (stream-extract the two records) ----
if [ ! -s "$REF" ]; then
  echo "[$(date)] building chr5+6 reference fasta (awk stream)"
  zcat "$SRCFA" | awk '/^>/{p=($1==">5"||$1==">6")} p' > "$REF"
  samtools faidx "$REF"
fi
echo "[$(date)] ref seqs:"; grep '^>' "$REF"

# ---- 2. chr5+6 GTF subset ----
if [ ! -s "$GTF" ]; then
  echo "[$(date)] subsetting GTF to chr5+6"
  zcat "$SRCGTF" | awk -F'\t' '$1=="5" || $1=="6" || $1 ~ /^#/' > "$GTF"
fi
echo "[$(date)] GTF lines: $(wc -l < "$GTF")"

# ---- 3. extract reads at the 7 loci -> fastq (dedup by read id) ----
if [ ! -s "$READS" ]; then
  echo "[$(date)] extracting reads at loci"
  REGS=$(awk '{printf "%s:%s-%s ", $2,$3,$4}' "$RUN/loci.tsv")
  samtools view -b -F 0x900 "$POOL" $REGS \
    | samtools collate -u -O - \
    | samtools fastq - 2>/dev/null \
    | awk 'NR%4==1{id=$1; if(seen[id]++){skip=1}else{skip=0}} skip==0{print}' > "$READS"
  echo "  reads: $(( $(wc -l < "$READS") / 4 ))"
fi

# ---- 4. uLTRA index (chr5+6, annotation-guided) ----
ULTRA_IDX=$RUN/refs/ultra_index
if [ ! -s "$ULTRA_IDX/refs_sequences.fa" ]; then
  echo "[$(date)] building uLTRA index"
  rm -rf "$ULTRA_IDX"; mkdir -p "$ULTRA_IDX"
  uLTRA index "$REF" "$GTF" "$ULTRA_IDX" 2>&1 | tail -5
fi

# ---- 5. deSALT index ----
DESALT_IDX=$RUN/refs/desalt_index
if [ ! -d "$DESALT_IDX" ] || [ -z "$(ls -A "$DESALT_IDX" 2>/dev/null)" ]; then
  echo "[$(date)] building deSALT index"
  mkdir -p "$DESALT_IDX"
  deSALT index "$REF" "$DESALT_IDX" 2>&1 | tail -3
fi

# ---- 6. alignments ----
mkdir -p "$RUN/aligners"
MM2=$RUN/aligners/minimap2.bam
ULT=$RUN/aligners/uLTRA.bam
DES=$RUN/aligners/deSALT.bam

if [ ! -s "$MM2" ]; then
  echo "[$(date)] minimap2"
  minimap2 -t $T -ax splice -uf -k14 --max-intron 500000 "$REF" "$READS" 2>"$RUN/logs/mm2.log" \
    | samtools sort -@ $T -o "$MM2" - && samtools index "$MM2"
fi

if [ ! -s "$ULT" ]; then
  echo "[$(date)] uLTRA align"
  UOUT=$RUN/aligners/ultra_out; mkdir -p "$UOUT"
  uLTRA align --ont --t $T --index "$ULTRA_IDX" --prefix loci "$REF" "$READS" "$UOUT" 2>&1 | tail -5
  samtools sort -@ $T -o "$ULT" "$UOUT/loci.sam" && samtools index "$ULT"
fi

if [ ! -s "$DES" ]; then
  echo "[$(date)] deSALT aln"
  deSALT aln -t $T -x ont1d -T -f "$RUN/aligners/desalt.tmp" -o "$RUN/aligners/deSALT.sam" "$DESALT_IDX" "$READS" 2>"$RUN/logs/desalt.log"
  samtools sort -@ $T -o "$DES" "$RUN/aligners/deSALT.sam" && samtools index "$DES"
  rm -f "$RUN/aligners/deSALT.sam"
fi

echo "[$(date)] alignment counts (primary):"
for b in "$MM2" "$ULT" "$DES"; do echo "  $b : $(samtools view -c -F 0x900 "$b")"; done

# ---- 7. per-aligner rectify correct (DRS default; no --annotation on human) ----
mkdir -p "$RUN/corrected"
for name in minimap2 uLTRA deSALT; do
  IN=$RUN/aligners/$name.bam
  OUTTSV=$RUN/corrected/${name}_corrected.tsv
  OUTBAM=$RUN/corrected/${name}_corrected.bam
  if [ -s "$OUTBAM" ]; then echo "[$(date)] skip correct $name (done)"; continue; fi
  echo "[$(date)] rectify correct $name"
  rectify correct "$IN" --genome "$REF" -o "$OUTTSV" \
    --organism homo_sapiens --write-corrected-bam "$OUTBAM" \
    -j $T 2>"$RUN/logs/correct_${name}.log"
  samtools index "$OUTBAM" 2>/dev/null || true
done

echo "[$(date)] Step D DONE"
ls -la $RUN/corrected/
