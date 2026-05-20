#!/bin/bash
#SBATCH --job-name=qsrev_cal
#SBATCH --account=larsms
#SBATCH --partition=larsms
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/users/kevinroy/calibration_2026_05/qsrev/slurm_logs/qsrev_cal_%j.out
#SBATCH --error=/scratch/users/kevinroy/calibration_2026_05/qsrev/slurm_logs/qsrev_cal_%j.err
#SBATCH --exclude=sh03-07n10,sh03-08n22

# Note: do NOT interleave plain `# comments` inside the #SBATCH block — some
# SLURM versions stop parsing directives at the first non-#SBATCH line, which
# silently dropped the AVX-512 exclude when this script was first submitted
# (job 25262022 landed on sh03-07n10 despite the directive). All explanatory
# comments below the directives, please.
#
# TRAP 1 (--account): without it the job hangs in (AssocGrpCPUMinutesLimit).
# TRAP 3 (AVX-512 SIGILL): sh03-07n10 + sh03-08n22 are AMD Milan; anaconda3
#   numpy uses AVX-512 intrinsics and SIGILLs within ~15 s on those nodes.
# TRAPS 2 + 4: do NOT enable `set -u` before sourcing bashrc; export PS1 first
# (otherwise /etc/profile.d/256term.sh references unbound $COLORTERM → 3 s death).
export PS1=""
source ~/.bashrc

# Sherlock conda activation
. /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify

set -u
set -e
set -o pipefail

# --- Paths ------------------------------------------------------------------
WORK=/scratch/users/kevinroy/calibration_2026_05/qsrev
PROFILER=$WORK/scripts/empirical_cigar_error_profiler.py
RAW=$WORK/raw
# Reference lives in the /oak rectify checkout, not the pip install of the env.
REF=/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz
test -f "$REF" || { echo "ERROR: reference not found: $REF" >&2; exit 1; }
test -f "$REF.fai" -a -f "$REF.gzi" || { echo "ERROR: missing .fai/.gzi for $REF" >&2; exit 1; }

echo "=== QSrev calibration run (Sherlock backup) ==="
echo "host:    $(hostname)"
echo "start:   $(date)"
echo "work:    $WORK"
echo "raw:     $RAW"
echo "profile: $PROFILER"
echo "ref:     $REF"

# --- Step 1: split (skip-if-exists) ----------------------------------------
mkdir -p "$WORK/chunks"
for sample in wt_R1 wt_R2; do
    out="$WORK/chunks/$sample"
    if compgen -G "$out/*chunk*.fastq.gz" >/dev/null; then
        echo "  [skip split] $sample already split"
        continue
    fi
    echo "  [split] $sample"
    mkdir -p "$out"
    rectify split "$RAW/$sample.fastq.gz" -n 8 -o "$out" --short-read --dT-primed-cDNA
done

# --- Step 2: align (bbmap + bwa) -------------------------------------------
for sample in wt_R1 wt_R2; do
    align_dir="$WORK/chunks/$sample/align"
    mkdir -p "$align_dir/bbmap" "$align_dir/bwa"
    for chunk_fq in "$WORK/chunks/$sample"/*chunk*.fastq.gz; do
        chunk_id=$(basename "$chunk_fq" .fastq.gz)
        if [ -f "$align_dir/bbmap/${chunk_id}.bbmap.bam" ] && \
           [ -f "$align_dir/bwa/${chunk_id}.bwa.bam" ]; then
            echo "  [skip align] $sample/$chunk_id"
            continue
        fi
        echo "  [align] $sample/$chunk_id"
        rectify align "$chunk_fq" \
            --short-read \
            --aligners bbmap bwa \
            --no-consensus \
            --genome "$REF" \
            --annotation /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz \
            -o "$align_dir" \
            --threads 4
        for a in bbmap bwa; do
            for f in "$align_dir/${chunk_id}.${a}.bam" "$align_dir/${chunk_id}.${a}.bam.bai"; do
                [ -e "$f" ] && mv "$f" "$align_dir/$a/"
            done
        done
        rm -f "$align_dir/${chunk_id}.rectified.bam" \
              "$align_dir/${chunk_id}.rectified.bam.bai" \
              "$align_dir/${chunk_id}_junctions.bed"
    done
done

# --- Step 3: profile per sample --------------------------------------------
OUT_BASE="$WORK/error_profile_qsrev_$(date +%Y%m%d)"
mkdir -p "$OUT_BASE"
for sample in wt_R1 wt_R2; do
    OUT="$OUT_BASE/$sample"
    mkdir -p "$OUT"
    echo "  [profile] $sample → $OUT"
    python "$PROFILER" \
        --protocol qsrev \
        --run-dir "$WORK/chunks/$sample" \
        --reference "$REF" \
        --output-dir "$OUT" \
        --isolation-flank 10 \
        --str-repeat \
        --workers 8 \
        2>&1 | tee "$OUT/profiler.log"
done

echo "=== Done ==="
echo "end:     $(date)"
echo "outputs: $OUT_BASE/{wt_R1,wt_R2}/error_counts.tsv etc."
