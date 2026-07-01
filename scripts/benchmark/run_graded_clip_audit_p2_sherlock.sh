#!/bin/bash
#SBATCH --job-name=graded_clip_p2
#SBATCH --account=larsms
#SBATCH --partition=larsms
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/users/kevinroy/graded_clip_audit/slurm_p2_%j.out
#SBATCH --error=/scratch/users/kevinroy/graded_clip_audit/slurm_p2_%j.err
#SBATCH --exclude=sh03-07n10,sh03-08n22

# Option-B 2nd pass (FIX-DEPLOYED) + calibration sweep.
#  - Uses a FRESH rectify tree at fc44ee2 (PYTHONPATH override) so the walkback
#    homopolymer-undercall guard is active — captures the cat2_plus_1-class 9D
#    preservation the stale-tree first pass (job 32183378) undercounted.
#  - Corrects each aligner ONCE → softclipped BAMs, then runs graded_clip_audit.py
#    across a (tail_frac × term_window) GRID (cheap — just re-scores the same BAMs)
#    to calibrate the terminal-window tail estimator.
#
# Traps (Sherlock/larsms): larsms acct/partition, exclude AMD Milan (SIGILL),
# PS1 + bashrc before set -u, conda after.

export PS1=""
source ~/.bashrc
. /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
set -u
set -o pipefail

FRESH=/scratch/users/kevinroy/rectify_fc44ee2          # fc44ee2 code tree (PYTHONPATH)
export PYTHONPATH=$FRESH${PYTHONPATH:+:$PYTHONPATH}
RECT=/oak/stanford/groups/larsms/Users/kevinroy/software/rectify   # data only
GEN=$RECT/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa
GFF=$RECT/rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz
SRC=/scratch/users/kevinroy/rectify_single_22321416_0
BASE=/scratch/users/kevinroy/graded_clip_audit
AUDIT=$BASE/graded_clip_audit.py
WORK=$BASE/p2_${SLURM_JOB_ID:-manual}
RC_FILE=$BASE/.audit_p2_rc
ALIGNERS="minimap2 gapmm2 mapPacBio deSALT uLTRA"

mkdir -p "$WORK/per_aligner" "$WORK/raw"
test -f "$GEN" || { echo "ERROR: genome not found: $GEN" >&2; printf '%s\n' 2 > "$RC_FILE"; exit 2; }
test -f "$AUDIT" || { echo "ERROR: audit script not found: $AUDIT" >&2; printf '%s\n' 2 > "$RC_FILE"; exit 2; }

echo "=== graded clip audit P2 (fix-deployed) — wt_by4742_rep1 ==="
echo "host: $(hostname)  start: $(date)  work: $WORK"
python -c "import rectify,os; from rectify.core.correct import walkback; print('rectify pkg:', os.path.dirname(rectify.__file__)); print('fc44ee2 guard active:', hasattr(walkback,'_MIN_GENOMIC_ANCHOR_3P'))" \
    || { echo "ERROR: fresh rectify import failed" >&2; printf '%s\n' 4 > "$RC_FILE"; exit 4; }

# --- Step 1: correct once per aligner (sort+index raw first) ------------------
for al in $ALIGNERS; do
    raw=$SRC/wt_by4742_rep1.$al.bam
    if [ ! -f "$raw" ]; then echo "  ! missing $raw; skipping $al" >&2; continue; fi
    echo "--> sort+index: $al"
    samtools sort -@ 4 -o "$WORK/raw/$al.bam" "$raw" \
        || { echo "  SORT FAILED $al" >&2; printf '%s\n' 3 > "$RC_FILE"; exit 3; }
    samtools index "$WORK/raw/$al.bam"
    echo "--> rectify correct: $al"
    python -m rectify.cli correct "$WORK/raw/$al.bam" \
        --genome "$GEN" --annotation "$GFF" \
        -o "$WORK/per_aligner/$al.tsv" \
        --write-softclipped-bam "$WORK/per_aligner/$al.softclipped.bam" \
        || { echo "  FAILED $al" >&2; printf '%s\n' 3 > "$RC_FILE"; exit 3; }
    samtools index "$WORK/per_aligner/$al.softclipped.bam"
done

# --- Step 2: calibration grid (re-score the same corrected BAMs) --------------
SUMMARY=$WORK/calibration_summary.tsv
echo -e "tail_frac\tterm_window\tn_flip\tn_tailflip\tn_suspicious" > "$SUMMARY"
RC=0
for tf in 0.5 0.6 0.7 0.8; do
  for tw in 8 12 16; do
    out=$WORK/reshuffle_tf${tf}_tw${tw}.tsv
    log=$(python "$AUDIT" --per-aligner-dir "$WORK/per_aligner" --genome "$GEN" \
            --tail-frac $tf --term-window $tw --out "$out" 2>&1)
    st=$?; [ $st -ne 0 ] && RC=$st
    # parse the "=== N/M reads change ... (X tail-flips, Y SUSPICIOUS) ===" line
    line=$(echo "$log" | grep "reads change winner")
    nflip=$(echo "$line" | sed -E 's#.* ([0-9]+)/[0-9]+ reads change.*#\1#')
    ntail=$(echo "$line" | sed -E 's#.*\(([0-9]+) tail-flips.*#\1#')
    nsus=$(echo "$line"  | sed -E 's#.*, ([0-9]+) SUSPICIOUS.*#\1#')
    echo -e "${tf}\t${tw}\t${nflip}\t${ntail}\t${nsus}" >> "$SUMMARY"
    echo "  tf=$tf tw=$tw -> flip=$nflip tail=$ntail susp=$nsus"
  done
done

echo "=== calibration summary ==="; cat "$SUMMARY"
cp "$SUMMARY" "$BASE/calibration_summary_${SLURM_JOB_ID:-manual}.tsv" 2>/dev/null
# default-params reshuffle for direct 1st-vs-2nd-pass comparison
cp "$WORK/reshuffle_tf0.6_tw12.tsv" "$BASE/graded_clip_reshuffle_p2_${SLURM_JOB_ID:-manual}.tsv" 2>/dev/null
echo "=== done rc=$RC  $(date) ==="
printf '%s\n' "$RC" > "$RC_FILE.tmp" && mv "$RC_FILE.tmp" "$RC_FILE"
exit $RC
