#!/bin/bash
# =============================================================================
# WS-1: COMPETITIVE alignment -> fair four-way + biological-excess contrast
# -----------------------------------------------------------------------------
# Settles the RNA004 "hot tail" puzzle by removing the contamination confound:
# the RNA004/LRGASP SIRV numbers came from aligning to a 430 KB SIRV-ONLY ref, so
# human reads that partially match a SIRV contig were KEPT with high error rates.
# Fix: align the RAW reads to the COMBINED hg38+spike-in reference (human maps to
# human; only true spike-ins land on the spike-in contigs), then re-measure.
# Self-contained (no agent loop): alignment + measurement + qscore in one job.
# The heavy/durable part is the alignment; the BAM persists for re-measurement.
# =============================================================================
#SBATCH --job-name=comp_fourway
#SBATCH --account=larsms
#SBATCH --partition=larsms
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/scratch/users/kevinroy/combined_ref/comp.%j.out
#SBATCH --error=/scratch/users/kevinroy/combined_ref/comp.%j.err

set -eo pipefail
D=/home/groups/larsms/users/kevinroy/aligner_bench_live
CB=/home/groups/larsms/users/kevinroy/anaconda3
P3=$CB/envs/pbsim3
PY=$CB/envs/rectify/bin/python
W=/scratch/users/kevinroy/combined_ref
mkdir -p "$W"; cd "$W"
MM2=$P3/bin/minimap2; SAMTOOLS=$P3/bin/samtools
source "$CB/etc/profile.d/conda.sh"; conda activate "$CB/envs/rectify"
set -u
T=${SLURM_CPUS_PER_TASK:-16}
VAL="$D/scripts/benchmark/error_realism_validate.py"

# ---- 1. combined hg38 + spike-in reference (spike-in contigs: chrIS, SIRVomeERCCome)
REF=$W/hg38_sirv.fa
B=https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/annotations/genome_fasta/hg38_sequins_SIRV_ERCCs_longSIRVs.fa
if [[ ! -s $REF ]]; then
  echo "[ws1] downloading combined ref (~2.9 GB) ..."
  wget -q -O $REF.part "$B"; mv $REF.part $REF
fi
[[ -s $REF.fai ]] || "$SAMTOOLS" faidx "$REF"

# ---- 2. competitive align each SIRV-only source to the COMBINED ref ----------
# reads: name -> raw fastq
align() {  # align <tag> <fastq>
  local tag="$1" fq="$2" bam="$W/${1}_combined.sorted.bam"
  [[ -s "$bam" ]] && { echo "[ws1] $tag bam exists"; return 0; }
  echo "[ws1] aligning $tag ..."
  zcat "$fq" | awk 'NR%4==2{gsub(/U/,"T");gsub(/u/,"t")}1' \
    | "$MM2" -ax splice -uf -k14 --MD --secondary=no -t "$T" "$REF" - \
    | "$SAMTOOLS" sort -@ "$T" -o "$bam" -
  "$SAMTOOLS" index "$bam"
}
align h69   /scratch/users/kevinroy/longbench_work/H69_dRNA_ONT.fastq.gz
align encff /scratch/users/kevinroy/sirv_work/ENCFF155CFF.fastq.gz

# ---- 3. measure: competitive SPIKE-IN subset + ENDOGENOUS human subset -------
REPORT=$W/competitive_fourway.report.txt
set +e
{
echo "==== WS-1 COMPETITIVE FOUR-WAY ($(date -u +%FT%TZ)) ===="
for tag in h69 encff; do
  bam=$W/${tag}_combined.sorted.bam
  spk=$("$SAMTOOLS" view -c -F0x904 "$bam" SIRVomeERCCome 2>/dev/null)
  echo "---- $tag : SPIKE-IN (SIRVomeERCCome) competitive, n_primary=$spk ----"
  "$SAMTOOLS" view -b -F0x900 "$bam" SIRVomeERCCome > $W/${tag}_spk.bam 2>/dev/null
  "$SAMTOOLS" index $W/${tag}_spk.bam
  PYTHONPATH=$D "$PY" "$VAL" measure-bam $W/${tag}_spk.bam --thin 0.0153 --len-min 600 --len-max 1000
  PYTHONPATH=$D "$PY" "$VAL" autocorr-bam $W/${tag}_spk.bam --len-min 600 --len-max 1000 --frac 0.3
  echo "---- $tag : ENDOGENOUS human (chr1 subset, biological-excess contrast) ----"
  "$SAMTOOLS" view -b -F0x900 "$bam" chr1 > $W/${tag}_endo_chr1.bam 2>/dev/null
  "$SAMTOOLS" index $W/${tag}_endo_chr1.bam
  endo_n=$("$SAMTOOLS" view -c $W/${tag}_endo_chr1.bam)
  echo "endogenous chr1 n_primary=$endo_n"
  PYTHONPATH=$D "$PY" "$VAL" measure-bam $W/${tag}_endo_chr1.bam --thin 0.0153 --len-min 600 --len-max 1000
done
echo "NOTE: contamination check = compare $tag spike-in n vs the OLD SIRV-only counts"
echo "      (RNA004 old=17461 in longbench_work/h69_sirv.sorted.bam; ENCFF old=4424 in sirv_work/sirv.sorted.bam)."
echo "      A drop = contamination removed. SIRV-vs-endogenous (same sample) = biological excess."
} 2>&1 | tee "$REPORT"

# ---- 4. qscore distribution (reference-FREE) on each raw fastq ----------------
QREPORT=$W/qscore_dist.report.txt
{
for pair in "h69 /scratch/users/kevinroy/longbench_work/H69_dRNA_ONT.fastq.gz" \
            "encff /scratch/users/kevinroy/sirv_work/ENCFF155CFF.fastq.gz"; do
  set -- $pair; tag="$1"; fq="$2"
  echo "---- $tag mean-qscore distribution (first 200k reads) ----"
  zcat "$fq" | "$PY" - "$tag" <<'PY'
import sys
tag=sys.argv[1]; import statistics as st
qs=[]
i=0
for ln in sys.stdin:
    i+=1
    if i%4==0:  # quality line
        q=ln.strip()
        if q: qs.append(sum(ord(c)-33 for c in q)/len(q))
    if len(qs)>=200000: break
if not qs: print("  no reads"); raise SystemExit
qs.sort(); n=len(qs)
import collections
hist=collections.Counter(min(40,int(x)) for x in qs)
mean=sum(qs)/n; med=qs[n//2]
for thr in (7,10,12,15):
    frac=sum(1 for x in qs if x<thr)/n
    print(f"  frac mean-q < {thr}: {frac:.3f}")
print(f"  mean={mean:.2f} median={med:.2f} n={n}")
print("  hist (q-bin:count):", dict(sorted(hist.items())))
PY
done
} 2>&1 | tee "$QREPORT"

RC=$?
printf '%s\n' "$RC" > $W/.comp_rc.tmp && mv $W/.comp_rc.tmp $W/.comp_rc
echo "[ws1] DONE rc=$RC reports: $REPORT ; $QREPORT"
exit $RC
