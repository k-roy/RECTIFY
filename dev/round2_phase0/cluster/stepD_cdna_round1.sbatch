#!/bin/bash
#SBATCH --job-name=cdna_stepD
#SBATCH --partition=owners
#SBATCH --constraint="CPU_GEN:SKX|CPU_GEN:SPR|CPU_GEN:GEN|CPU_GEN:BGM"
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --output=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0/cdna_a549/logs/stepD.%j.out
# PCR-cDNA Phase-0 retest — Step D: Round-1 genome alignment (3-aligner panel) + per-aligner
# rectify correct, in ONT-cDNA mode. Reuses the EXISTING A549 chr5/6 reference + uLTRA/deSALT
# indexes built for the directRNA run. Consumes the FASTQ produced by prep_extract_cdna_reads.sh
# (run on the login node first). Idempotent skip-checks throughout.
set -eo pipefail
export PS1=""
source ~/.bashrc
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
set -u

P0=/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0
RUN=$P0/cdna_a549
SRC=$P0/round2                 # reuse refs + indexes from the directRNA run
T=$SLURM_CPUS_PER_TASK

REF=$SRC/refs/chr5_6.fa
GTF=$SRC/refs/chr5_6.gtf
ULTRA_IDX=$SRC/refs/ultra_index
DESALT_IDX=$SRC/refs/desalt_index
READS=$RUN/reads/loci_reads.fastq

cd $RUN
echo "[$(date)] cDNA Step D start on $(hostname); cpus=$T"
[ -s "$READS" ] || { echo "FATAL: $READS missing — run prep_extract_cdna_reads.sh on the login node first"; exit 1; }
echo "[$(date)] reads: $(( $(wc -l < "$READS") / 4 ))"
for p in "$REF" "$GTF" "$ULTRA_IDX/database.db" "$DESALT_IDX/ref.seq"; do [ -s "$p" ] || { echo "FATAL: missing reused asset $p"; exit 1; }; done

mkdir -p "$RUN/aligners" "$RUN/logs"
MM2=$RUN/aligners/minimap2.bam
ULT=$RUN/aligners/uLTRA.bam
DES=$RUN/aligners/deSALT.bam

# minimap2 spliced, PCR-cDNA mode: NO -uf (reads come from both strands; minimap2 infers
# transcript strand from splice motifs via the ts tag). This is the key DRS->cDNA change.
if [ ! -s "$MM2" ]; then
  echo "[$(date)] minimap2 -ax splice (cDNA, no -uf)"
  minimap2 -t $T -ax splice -k14 --max-intron 500000 "$REF" "$READS" 2>"$RUN/logs/mm2.log" \
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

# per-aligner rectify correct in ONT-cDNA mode (--ONT-cDNA: SQK-PCB114 PCR-cDNA; poly-A tail
# IS sequenced, is_reverse->minus-strand gene, same convention as DRS).
mkdir -p "$RUN/corrected"
for name in minimap2 uLTRA deSALT; do
  IN=$RUN/aligners/$name.bam
  OUTTSV=$RUN/corrected/${name}_corrected.tsv
  OUTBAM=$RUN/corrected/${name}_corrected.bam
  if [ -s "$OUTBAM" ]; then echo "[$(date)] skip correct $name (done)"; continue; fi
  echo "[$(date)] rectify correct $name (--ONT-cDNA)"
  # cap parallelism at 4: --ONT-cDNA correction duplicates the human chr5+6 genome
  # structures per worker; -j 8 OOM'd at 48G (job 29631826). -j 4 @ 128G is safe and the
  # 4k-read job is tiny so wall-time is unaffected.
  rectify correct "$IN" --genome "$REF" -o "$OUTTSV" \
    --ONT-cDNA --organism homo_sapiens --write-corrected-bam "$OUTBAM" \
    -j 4 2>"$RUN/logs/correct_${name}.log"
  samtools index "$OUTBAM" 2>/dev/null || true
done

echo "[$(date)] cDNA Step D DONE"
ls -la $RUN/corrected/
