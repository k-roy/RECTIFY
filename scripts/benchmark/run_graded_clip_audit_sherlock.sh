#!/bin/bash
#SBATCH --job-name=graded_clip_audit
#SBATCH --account=larsms
#SBATCH --partition=larsms
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/users/kevinroy/graded_clip_audit/slurm_%j.out
#SBATCH --error=/scratch/users/kevinroy/graded_clip_audit/slurm_%j.err
#SBATCH --exclude=sh03-07n10,sh03-08n22

# Option-B EER-ED winner-reshuffle audit on real DRS data (wt_by4742_rep1).
# Generates corrected per-aligner softclipped BAMs, then runs graded_clip_audit.py
# to compare current flat-clip winners vs graded terminal-window clip-penalty
# winners. See rectify/data/validation/CASE_STUDIES.md cat2_plus_1 + HANDOFF.md.
#
# Traps (Sherlock/larsms): --account/--partition=larsms, --exclude AMD Milan
# (numpy SIGILL), PS1 export + bashrc source before `set -u`, conda after.
#
# FIRST-PASS CAVEAT: this uses the deployed Oak rectify, which may be STALE
# (missing the fc44ee2 walkback fix). The graded-penalty blast radius / false-
# positive rate is largely fix-independent, so the first pass is informative;
# for the fix-synergy pass, deploy fc44ee2's rectify to a fresh tree first and
# point WORK's `python -m rectify.cli` at it via PYTHONPATH.

export PS1=""
source ~/.bashrc
. /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
set -u
set -o pipefail

RECT=/oak/stanford/groups/larsms/Users/kevinroy/software/rectify
GEN=$RECT/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa
GFF=$RECT/rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz
SRC=/scratch/users/kevinroy/rectify_single_22321416_0          # raw per-aligner BAMs
BASE=/scratch/users/kevinroy/graded_clip_audit
AUDIT=$BASE/graded_clip_audit.py                               # deployed alongside this script
WORK=$BASE/run_${SLURM_JOB_ID:-manual}
RC_FILE=$BASE/.audit_rc                                        # sentinel for watchers
ALIGNERS="minimap2 gapmm2 mapPacBio deSALT uLTRA"

mkdir -p "$WORK/per_aligner"
test -f "$GEN" || { echo "ERROR: genome not found: $GEN" >&2; printf '%s\n' 2 > "$RC_FILE"; exit 2; }
test -f "$AUDIT" || { echo "ERROR: audit script not found: $AUDIT" >&2; printf '%s\n' 2 > "$RC_FILE"; exit 2; }

echo "=== graded clip-penalty audit (wt_by4742_rep1) ==="
echo "host: $(hostname)  start: $(date)  work: $WORK"
python -c "from rectify.core.correct import walkback; print('rectify walkback has fc44ee2 guard:', hasattr(walkback,'_MIN_GENOMIC_ANCHOR_3P'))"

# --- Step 1: corrected per-aligner softclipped BAMs (tail bases retained) ----
for al in $ALIGNERS; do
    raw=$SRC/wt_by4742_rep1.$al.bam
    if [ ! -f "$raw" ]; then echo "  ! missing $raw; skipping $al" >&2; continue; fi
    echo "--> rectify correct: $al"
    python -m rectify.cli correct "$raw" \
        --genome "$GEN" --annotation "$GFF" \
        -o "$WORK/per_aligner/$al.tsv" \
        --write-softclipped-bam "$WORK/per_aligner/$al.softclipped.bam" \
        || { echo "  FAILED $al" >&2; printf '%s\n' 3 > "$RC_FILE"; exit 3; }
    samtools index "$WORK/per_aligner/$al.softclipped.bam"
done

# --- Step 2: graded vs flat winner-reshuffle audit ---------------------------
echo "--> graded_clip_audit.py"
python "$AUDIT" \
    --per-aligner-dir "$WORK/per_aligner" \
    --genome "$GEN" \
    --tail-frac 0.6 --term-window 12 \
    --out "$WORK/graded_clip_reshuffle.tsv"
RC=$?

cp "$WORK/graded_clip_reshuffle.tsv" "$BASE/graded_clip_reshuffle_${SLURM_JOB_ID:-manual}.tsv" 2>/dev/null
echo "=== done rc=$RC  $(date) ==="
printf '%s\n' "$RC" > "$RC_FILE.tmp" && mv "$RC_FILE.tmp" "$RC_FILE"   # atomic sentinel
exit $RC
