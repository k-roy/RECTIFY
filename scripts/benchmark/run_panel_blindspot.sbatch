#!/bin/bash
#SBATCH --job-name=nj_panel
#SBATCH --account=larsms
#SBATCH --partition=larsms
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/scratch/users/kevinroy/nj_panel/nj_panel.%j.out
# Full 5-aligner PANEL native novel-junction recovery on the deviation ladder.
# Self-contained + idempotent + sentinel-rc (per the cluster discipline). Runs the
# CURRENT code from the synced dir (NOT the pre-C1 oak checkout).
set -uo pipefail
CODE=/home/groups/larsms/users/kevinroy/nj_panel_code
CORPUS=/scratch/users/kevinroy/nj_panel/corpus
SENT=/scratch/users/kevinroy/nj_panel/.nj_panel_rc
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify
export PYTHONPATH="$CODE"
cd "$CODE"
python scripts/benchmark/run_panel_on_corpus.py --corpus "$CORPUS" \
    --aligners minimap2 gapmm2 deSALT uLTRA --threads 8
RC=$?
printf '%s\n' "$RC" > "${SENT}.tmp" && mv "${SENT}.tmp" "$SENT"
echo "[sbatch] panel run rc=$RC; report at $CORPUS/panel_blindspot.report.txt"
exit $RC
