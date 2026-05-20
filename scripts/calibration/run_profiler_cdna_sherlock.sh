#!/bin/bash
#SBATCH --job-name=cdna_cal
#SBATCH --account=larsms
#SBATCH --partition=larsms
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/users/kevinroy/calibration_2026_05/cdna/slurm_logs/cdna_cal_%j.out
#SBATCH --error=/scratch/users/kevinroy/calibration_2026_05/cdna/slurm_logs/cdna_cal_%j.err
#SBATCH --exclude=sh03-07n10,sh03-08n22

# Sherlock cDNA backup — runs only 2 of 16 chunks (same -n 16 split as H2 so
# the resulting penalties are directly comparable). The script's intent is to
# answer "would 16 chunks tell us anything 2 chunks don't?" by profiling
# chunks 0 and 1 separately and pooled, so per-chunk variation can be eyeballed
# before deciding whether to wait out the H2 16-chunk run.
#
# Traps applied: --account=larsms, --exclude AMD Milan, set -u defer past
# bashrc source, PS1 export so /etc/profile.d doesn't crash on $COLORTERM.

export PS1=""
source ~/.bashrc

. /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify

set -u
set -e
set -o pipefail

# --- Paths ------------------------------------------------------------------
WORK=/scratch/users/kevinroy/calibration_2026_05/cdna
CONS_BAM=$WORK/raw/stage1_consensus.bam
PROFILER=$WORK/scripts/empirical_cigar_error_profiler.py
REF=/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz
CHUNKS_TO_RUN="0 1 2 3 4 5 6 7"   # 8-chunk run for per-UMI-bin convergence

test -f "$REF" || { echo "ERROR: reference not found: $REF" >&2; exit 1; }
test -f "$REF.fai" -a -f "$REF.gzi" || { echo "ERROR: missing .fai/.gzi for $REF" >&2; exit 1; }
test -f "$CONS_BAM" || { echo "ERROR: consensus BAM not found: $CONS_BAM" >&2; exit 1; }

echo "=== cDNA calibration run (Sherlock backup, chunks $CHUNKS_TO_RUN of 16) ==="
echo "host:      $(hostname)"
echo "start:     $(date)"
echo "work:      $WORK"
echo "consensus: $CONS_BAM ($(numfmt --to=iec $(stat -c %s "$CONS_BAM")))"

# --- Step 1: BAM -> FASTQ --------------------------------------------------
FQ="$WORK/wt_rep2_consensus.fastq.gz"
if [ -f "$FQ" ]; then
    echo "  [skip fastq] $FQ exists"
else
    echo "  [fastq] $CONS_BAM -> $FQ"
    samtools fastq -@ 4 "$CONS_BAM" 2>/dev/null | gzip -1 > "$FQ"
fi

# --- Step 2: split into 16 chunks (mirrors H2 layout) ----------------------
CHUNKS="$WORK/chunks"
if compgen -G "$CHUNKS/*chunk*.fastq.gz" >/dev/null; then
    echo "  [skip split] $CHUNKS already populated"
else
    echo "  [split] -n 16"
    mkdir -p "$CHUNKS"
    rectify split "$FQ" -n 16 -o "$CHUNKS"
fi

# --- Step 3: align ONLY the requested chunks ------------------------------
ALIGN_DIR="$CHUNKS/align"
for a in minimap2 mapPacBio gapmm2 uLTRA deSALT; do
    mkdir -p "$ALIGN_DIR/$a"
done
for c in $CHUNKS_TO_RUN; do
    cidx=$(printf '%03d' "$c")
    chunk_fq=$(ls "$CHUNKS"/*chunk_${cidx}_of_*.fastq.gz | head -1)
    chunk_id=$(basename "$chunk_fq" .fastq.gz)
    missing=0
    for a in minimap2 mapPacBio gapmm2 uLTRA deSALT; do
        [ -f "$ALIGN_DIR/$a/${chunk_id}.$a.bam" ] || missing=$((missing+1))
    done
    if [ "$missing" -eq 0 ]; then
        echo "  [skip align] $chunk_id"
        continue
    fi
    echo "  [align] $chunk_id ($missing aligners missing)"
    rectify align "$chunk_fq" \
        --aligners minimap2 mapPacBio gapmm2 \
        --junction-aligners uLTRA deSALT \
        --parallel-aligners \
        --no-consensus \
        --genome "$REF" \
        --annotation /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz \
        -o "$ALIGN_DIR" \
        --threads 8
    for a in minimap2 mapPacBio gapmm2 uLTRA deSALT; do
        for f in "$ALIGN_DIR/${chunk_id}.${a}.bam" "$ALIGN_DIR/${chunk_id}.${a}.bam.bai"; do
            [ -e "$f" ] && mv "$f" "$ALIGN_DIR/$a/"
        done
    done
    rm -f "$ALIGN_DIR/${chunk_id}.rectified.bam" \
          "$ALIGN_DIR/${chunk_id}.rectified.bam.bai" \
          "$ALIGN_DIR/${chunk_id}_junctions.bed"
done

# --- Step 4: profile per chunk AND pooled ---------------------------------
OUT_BASE="$WORK/error_profile_cdna_$(date +%Y%m%d)"
mkdir -p "$OUT_BASE"

# Per-chunk runs so we can compare variation
for c in $CHUNKS_TO_RUN; do
    OUT="$OUT_BASE/chunk_$c"
    mkdir -p "$OUT"
    echo "  [profile] chunk $c -> $OUT"
    python "$PROFILER" \
        --protocol cdna \
        --run-dir "$CHUNKS" \
        --chunks "$c" \
        --reference "$REF" \
        --output-dir "$OUT" \
        --isolation-flank 10 \
        --str-repeat \
        --umi-tag XC \
        --umi-bins 1,2,3+ \
        --umi-source-bam "$CONS_BAM" \
        --workers 1 \
        2>&1 | tee "$OUT/profiler.log"
done

# Pooled run across all CHUNKS_TO_RUN
OUT="$OUT_BASE/pooled"
mkdir -p "$OUT"
echo "  [profile] pooled $CHUNKS_TO_RUN -> $OUT"
python "$PROFILER" \
    --protocol cdna \
    --run-dir "$CHUNKS" \
    --chunks $CHUNKS_TO_RUN \
    --reference "$REF" \
    --output-dir "$OUT" \
    --isolation-flank 10 \
    --str-repeat \
    --umi-tag XC \
    --umi-bins 1,2,3+ \
    --umi-source-bam "$CONS_BAM" \
    --workers $(echo $CHUNKS_TO_RUN | wc -w) \
    2>&1 | tee "$OUT/profiler.log"

echo "=== Done ==="
echo "end:      $(date)"
echo "outputs:  $OUT_BASE/{chunk_$(echo $CHUNKS_TO_RUN | tr ' ' '_'),pooled}/penalty_scores.tsv"
echo
echo "Next: compare chunk_0/penalty_scores.tsv vs chunk_1/penalty_scores.tsv"
echo "vs pooled/penalty_scores.tsv. If per-chunk rates agree to within ~5-10%,"
echo "2 chunks is enough — abort H2 and bundle the Sherlock 'pooled' table."
