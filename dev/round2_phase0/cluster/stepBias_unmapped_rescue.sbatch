#!/bin/bash
#SBATCH --job-name=r2p0_bias
#SBATCH --partition=owners
#SBATCH --constraint="CPU_GEN:SKX|CPU_GEN:SPR|CPU_GEN:GEN|CPU_GEN:BGM"
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0/round2/logs/stepBias.%j.out
# Selection-bias check (Kevin): test the cDNA round on GENOME-UNMAPPED reads — the population
# the region-extracted Round-1 test excluded. Extract unmapped reads from the 4 SG-NEx A549
# replicates, align to the 7-locus cDNA library; reads that map there are genome-unmapped but
# cDNA-mappable = the strongest rescue candidates. Confirm with a fresh genome re-align.
set -eo pipefail
export PS1=""
source ~/.bashrc
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
set -u

P0=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0
RUN=$P0/round2
DATA=$P0/data
LIB=$RUN/library
REF=$RUN/refs/chr5_6.fa
T=$SLURM_CPUS_PER_TASK
mkdir -p $RUN/bias/{reads,aln,logs}
cd $RUN
echo "[$(date)] bias check start on $(hostname)"

# ---- 1. extract unmapped primary reads from all replicates -> one fastq ----
UM=$RUN/bias/reads/unmapped.fastq
if [ ! -s "$UM" ]; then
  echo "[$(date)] extracting unmapped reads from 4 replicates"
  : > "$UM"
  for r in 1 4 5 6; do
    b=$DATA/SGNex_A549_directRNA_replicate${r}_run1.bam
    samtools fastq -f 4 -F 0x900 "$b" 2>/dev/null >> "$UM"
  done
  echo "  unmapped reads: $(( $(wc -l < "$UM") / 4 ))"
fi

# ---- 2. align unmapped reads to the cDNA library ----
CDNA=$RUN/bias/aln/unmapped_to_cdna.bam
if [ ! -s "$CDNA" ]; then
  echo "[$(date)] minimap2 unmapped -> cDNA library"
  minimap2 -t $T -ax map-ont --for-only -N 100 -p 0.8 -Y "$LIB/library.fa" "$UM" 2>"$RUN/bias/logs/mm2_cdna.log" \
    | samtools sort -@ $T -o "$CDNA" - && samtools index "$CDNA"
fi
NMAP=$(samtools view -c -F 0x904 "$CDNA" 2>/dev/null)
echo "[$(date)] unmapped reads that map to a cDNA (primary): $NMAP"

# ---- 3. pull those read IDs + their fastq, re-align to GENOME (fresh panel) ----
IDS=$RUN/bias/reads/cdna_hit_ids.txt
HITFQ=$RUN/bias/reads/cdna_hit.fastq
if [ ! -s "$IDS" ]; then
  samtools view -F 0x904 "$CDNA" | cut -f1 | sort -u > "$IDS"
  echo "  unique cDNA-hit read ids: $(wc -l < "$IDS")"
fi
if [ ! -s "$HITFQ" ] && [ -s "$IDS" ]; then
  samtools fastq -F 0x900 "$CDNA" 2>/dev/null | \
    awk 'NR%4==1{h=substr($1,2); getline s; getline p; getline q; if(!seen[h]++){print "@"h"\n"s"\n"p"\n"q}}' > "$HITFQ"
  echo "  cDNA-hit fastq reads: $(( $(wc -l < "$HITFQ") / 4 ))"
fi
GEN=$RUN/bias/aln/cdna_hits_to_genome.bam
if [ ! -s "$GEN" ] && [ -s "$HITFQ" ]; then
  echo "[$(date)] minimap2 cDNA-hit reads -> genome (fresh, my Round-1 params)"
  minimap2 -t $T -ax splice -uf -k14 --max-intron 500000 "$REF" "$HITFQ" 2>"$RUN/bias/logs/mm2_genome.log" \
    | samtools sort -@ $T -o "$GEN" - && samtools index "$GEN"
fi

echo "[$(date)] bias extraction/alignment DONE"
ls -la $RUN/bias/aln/
