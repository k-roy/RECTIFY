#!/bin/bash
# Round-2 cDNA Phase-0 (PCR-cDNA retest) — PREP: remote-extract SG-NEx A549 PCR-cDNA reads
# at the 7 A549 loci into a single deduped FASTQ. RUN ON THE LOGIN NODE (needs outbound HTTPS
# via htslib remote BAM range-reads; compute nodes may lack outbound internet). Idempotent.
#
# Design: same cell line + same loci + same 42-cDNA library as the A549 directRNA Phase-0, so
# the ONLY changed variable is protocol/read-length (PCR-cDNA = full-length vs DRS = 3'-biased).
# This directly tests the opus-review hypothesis that DRS 5'-truncation suppresses the
# rescuable population (reads that terminate near a short-anchored junction).
set -eo pipefail
export PS1=""
source ~/.bashrc
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
set -u

P0=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0
RUN=$P0/cdna_a549
LOCI=$P0/round2/loci.tsv
mkdir -p "$RUN/reads" "$RUN/logs" "$RUN/aligners" "$RUN/corrected" "$RUN/round2_aln"
READS=$RUN/reads/loci_reads.fastq

S3=https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/sequencing_data_ont/bam/genome
REPS=(
  "SGNex_A549_cDNA_replicate1_run2/SGNex_A549_cDNA_replicate1_run2.bam"
  "SGNex_A549_cDNAStranded_replicate3_run3/SGNex_A549_cDNAStranded_replicate3_run3.bam"
  "SGNex_A549_cDNAStranded_replicate5_run2/SGNex_A549_cDNAStranded_replicate5_run2.bam"
)

if [ -s "$READS" ]; then
  echo "[$(date)] reads already extracted: $(( $(wc -l < "$READS") / 4 )) reads — skip"
  exit 0
fi

REGS=$(awk '{printf "%s:%s-%s ", $2,$3,$4}' "$LOCI")
echo "[$(date)] loci regions: $REGS"

TMP=$RUN/reads/_pool.fastq; : > "$TMP"
for rel in "${REPS[@]}"; do
  U=$S3/$rel
  name=$(basename "$rel" .bam)
  echo "[$(date)] remote-extracting $name ..."
  # primary, mapped reads only at the loci -> name-collate -> fastq
  samtools view -b -F 0x904 "$U" $REGS 2>>"$RUN/logs/extract.log" \
    | samtools collate -u -O - 2>>"$RUN/logs/extract.log" \
    | samtools fastq - 2>>"$RUN/logs/extract.log" >> "$TMP"
  echo "  cumulative lines: $(wc -l < "$TMP")"
done

# dedup by read id (a read may appear in >1 rep / multiple loci windows)
awk 'NR%4==1{id=$1; if(seen[id]++){skip=1}else{skip=0}} skip==0{print}' "$TMP" > "$READS"
rm -f "$TMP"
echo "[$(date)] PREP done: $(( $(wc -l < "$READS") / 4 )) unique reads -> $READS"
