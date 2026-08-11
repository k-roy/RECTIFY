#!/bin/bash
#SBATCH --job-name=tier2_bench
#SBATCH --account=larsms
#SBATCH --partition=larsms
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --time=00:45:00
#SBATCH --output=tier2_bench_%j.out
#SBATCH --error=tier2_bench_%j.err
# TRAP 3: AVX-512 SIGILL on AMD Milan nodes — anaconda3 numpy uses AVX-512
#SBATCH --exclude=sh03-07n10,sh03-08n22

# TRAP 2 + 4: export PS1 + source bashrc BEFORE set -u (256term.sh references an
# unbound $COLORTERM) so PATH/LD_LIBRARY_PATH propagate from the submission shell.
export PS1=""
source ~/.bashrc
set -u
set -o pipefail
# NOTE: deliberately NOT `set -e` — we want BOTH tasks to attempt and the rc
# sentinel to be written even if one fails (the durable signal for any watcher).

D=/home/groups/larsms/users/kevinroy/aligner_bench_live
CB=/home/groups/larsms/users/kevinroy/anaconda3
P3=$CB/envs/pbsim3
RB=$CB/envs/rectify
cd "$D"
export PYTHONPATH="$D"
mkdir -p "$D/tier2_results"

run_one () {
  # $1 = label (drs|cdna), $2 = errhmm model path
  local label="$1" model="$2"
  local out="$SCRATCH/tier2_${label}_${SLURM_JOB_ID}"
  echo "[tier2.sbatch] === $label === model=$model out=$out"
  "$RB/bin/python" scripts/benchmark/tier2_run.py \
    --gff "$D/refdata/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz" \
    --genome "$D/refdata/S288C_reference_sequence_R64-5-1_20240529.fsa" \
    --errhmm-model "$model" \
    --pbsim-bin "$P3/bin/pbsim" --minimap2 "$P3/bin/minimap2" \
    --samtools "$P3/bin/samtools" \
    --threads 4 --copies 20 --seed 7 --out "$out"
  local task_rc=$?
  if [ -f "$out/tier2_summary.json" ]; then
    cp "$out/tier2_summary.json" "$D/tier2_results/${label}_summary.json"
  fi
  return $task_rc
}

rc=0
run_one drs  "$P3/data/ERRHMM-ONT.model"     || rc=$?
run_one cdna "$P3/data/ERRHMM-ONT-HQ.model"  || rc=$?

# Durable sentinel (rc INSIDE, atomic) — any watcher branches on this.
printf '%s\n' "$rc" > "$D/tier2_results/.tier2_rc.tmp" \
  && mv "$D/tier2_results/.tier2_rc.tmp" "$D/tier2_results/.tier2_rc"
echo "[tier2.sbatch] DONE rc=$rc"
exit $rc
