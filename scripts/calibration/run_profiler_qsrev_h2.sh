#!/bin/bash
#$ -N qsrev_cal
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=8G,h_rt=4:00:00
#$ -o /u/scratch/k/kevinroy/calibration_2026_05/qsrev/slurm_logs/qsrev_cal_$JOB_ID.log
#
# Han 2023 QSrev (Illumina QuantSeq REV, 76 bp) penalty-table calibration.
# Driver: scripts/calibration/empirical_cigar_error_profiler.py --protocol qsrev
#
# Wall budget: split <5 min, align ~30 min (8 chunks x 2 aligners x 2 samples,
# serial), profile per sample ~15-30 min each. h_rt=4h for slack.
#
# Outputs: per-sample error_counts.tsv under $OUT_BASE/{wt_R1,wt_R2}/.
# Then bundle_protocol_tables.py merges counts and writes the final
# penalty_scores_qsrev.tsv + str_penalty_scores_qsrev.tsv into rectify/data/.
#
# Submit:   qsub run_profiler_qsrev_h2.sh
# Monitor:  qstat -j <JOB_ID>  ;  tail -f <slurm_logs>/qsrev_cal_<JOB_ID>.log
set -euo pipefail

# --- Env --------------------------------------------------------------------
# H2 batch shells need the modules init explicitly sourced; the default
# /etc/profile.d/modules.sh path does NOT exist on H2 compute nodes.
. /u/local/Modules/default/init/bash
module load conda/23.11.0
source activate rectify

# --- Paths ------------------------------------------------------------------
WORK=/u/scratch/k/kevinroy/calibration_2026_05/qsrev
RAW=/u/project/guillom/shared/raw/SRA/external/han_2023_quantseq_rev/fastq
PROFILER=/u/scratch/k/kevinroy/calibration_2026_05/scripts/empirical_cigar_error_profiler.py
# Bundled reference: H2 ships only .fsa.gz (BGZF) with .fai + .gzi indices.
# pysam.FastaFile reads bgzipped FASTAs transparently when both indices are present.
REF=/u/home/k/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz
test -f "$REF" || { echo "ERROR: reference not found: $REF" >&2; exit 1; }
test -f "$REF.fai" -a -f "$REF.gzi" || { echo "ERROR: missing .fai/.gzi for $REF" >&2; exit 1; }

echo "=== QSrev calibration run ==="
echo "host:     $(hostname)"
echo "start:    $(date)"
echo "work:     $WORK"
echo "raw:      $RAW"
echo "profiler: $PROFILER"
echo "ref:      $REF"

# --- Step 1: split WT FASTQs (8 chunks each) -------------------------------
mkdir -p "$WORK/chunks"
for sample in wt_R1 wt_R2; do
    out="$WORK/chunks/$sample"
    if compgen -G "$out/*chunk*.fastq.gz" >/dev/null; then
        echo "  [skip split] $sample already split"
        continue
    fi
    echo "  [split] $sample"
    mkdir -p "$out"
    rectify split "$RAW/$sample.fastq.gz" \
        -n 8 \
        -o "$out" \
        --short-read \
        --dT-primed-cDNA
done

# --- Step 2: align per chunk (bbmap + bwa) ---------------------------------
# rectify align writes BAMs to a flat dir as <chunk_id>.<aligner>.bam. The
# profiler's find_chunk_bams expects per-aligner subdirs (align/<aligner>/),
# so we move them into place after each chunk completes.
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
            --Scer \
            -o "$align_dir" \
            --threads 4
        # Reorganize: <chunk_id>.<aligner>.bam[.bai] → align/<aligner>/
        for a in bbmap bwa; do
            for f in "$align_dir/${chunk_id}.${a}.bam" "$align_dir/${chunk_id}.${a}.bam.bai"; do
                [ -e "$f" ] && mv "$f" "$align_dir/$a/"
            done
        done
        # Drop the rectified single-aligner sentinel (--no-consensus still
        # produces one; the profiler skips it via the 'rectified' filename
        # filter but cleaning up keeps the dir tidy).
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
echo
echo "Next steps (local M1):"
echo "  1. rsync $OUT_BASE M1:/path/to/rectify/scripts/calibration/runs/"
echo "  2. python scripts/calibration/bundle_protocol_tables.py \\"
echo "         --protocol qsrev \\"
echo "         --inputs runs/error_profile_qsrev_*/wt_R{1,2}"
echo "  3. pytest tests/test_data_bundling.py tests/test_quantseq_rev_walkback.py"
