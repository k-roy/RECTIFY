#!/bin/bash
# Phase D: reanchor optimization timing
# Measures lazy-path wall time after _apply_reanchor_from_clip_len replaces
# the per-base genome scan in apply_corrected_edits_to_read.
#
# Baseline (Phase C, 2026-05-21, job 25662145):
#   Total lazy:  495.1 s
#   mapPacBio HP scoring: 142.3 s (9.3 µs/read, 38% reanchor rate)
#   minimap2  HP scoring:  54.5 s (3.8 µs/read, 10.4% reanchor rate)
#   gapmm2    HP scoring: 102.5 s (7.1 µs/read, 15.5% reanchor rate)
#   BAM write:            136.2 s
#
# Expected improvement: ~25-35% reduction in mapPacBio HP scoring time.

#SBATCH --job-name=rectify_reanchor_timing
#SBATCH --partition=larsms,owners
#SBATCH --account=larsms
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/timing/phase_d_reanchor_timing_%j.out
#SBATCH --error=/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/timing/phase_d_reanchor_timing_%j.err

export MKL_DEBUG_CPU_TYPE=5

export PS1=""
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify

set -euo pipefail

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export LOKY_MAX_CPU_COUNT=4
export NUMEXPR_MAX_THREADS=4

PYTHON="/home/groups/larsms/users/kevinroy/anaconda3/bin/python"
RECTIFY_SRC="/oak/stanford/groups/larsms/Users/kevinroy/software/rectify"
OUTDIR="/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/wt_tfiiib_rep3/chunks"
TIMING_DIR="/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/timing"
CHUNK="008"
GENOME_FSA="${RECTIFY_SRC}/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
PENALTY_TSV="${RECTIFY_SRC}/rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv"
STR_PENALTY_TSV="${RECTIFY_SRC}/rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/str_penalty_scores.tsv"
POOL_PKL="${OUTDIR}/junction_pool.pkl"
SCAN_PKL="${OUTDIR}/rescue_scan.pkl"

mkdir -p "$TIMING_DIR"

echo "Host:   $(hostname)"
echo "Job:    ${SLURM_JOB_ID:-local}"
echo "Chunk:  $CHUNK"
echo "Start:  $(date)"
echo "CPUs:   ${SLURM_CPUS_PER_TASK:-4}"

SCRATCH_WORK="${SCRATCH:-/tmp}/rectify_reanchor_timing_${SLURM_JOB_ID:-$$}"
mkdir -p "$SCRATCH_WORK"
trap "rm -rf '$SCRATCH_WORK'" EXIT

cd "$RECTIFY_SRC"
$PYTHON - << PYEOF
import sys, time, tracemalloc, logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

sys.path.insert(0, '/oak/stanford/groups/larsms/Users/kevinroy/software/rectify')

from rectify.core.consensus.corrected_consensus import merge_corrected_tsvs, write_corrected_consensus_bam
from rectify.utils.genome import load_genome
from rectify.core.splice.hp_penalty import HpPenaltyTable

OUTDIR     = Path('/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/wt_tfiiib_rep3/chunks')
TIMING_DIR = Path('/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/timing')
SCRATCH    = Path('${SCRATCH_WORK}')
CHUNK      = '008'
JOB_ID     = '${SLURM_JOB_ID:-local}'
GENOME_FSA = '${GENOME_FSA}'
PENALTY_TSV    = '${PENALTY_TSV}'
STR_PENALTY_TSV = '${STR_PENALTY_TSV}'
POOL_PKL   = str(OUTDIR / 'junction_pool.pkl')
SCAN_PKL_PATH = OUTDIR / 'rescue_scan.pkl'
SCAN_PKL   = str(SCAN_PKL_PATH) if SCAN_PKL_PATH.exists() else None

ALIGNERS = ['minimap2', 'mapPacBio', 'gapmm2']

per_aligner_tsvs = {}
per_aligner_raw_bams = {}
for a in ALIGNERS:
    chunk_dir = OUTDIR / 'aligner_chunks' / a / f'chunk_{CHUNK}'
    tsv = chunk_dir / 'corrected_reads.tsv'
    raw_bam_matches = list(chunk_dir.glob(f'*.{a}.bam'))
    if tsv.exists():
        per_aligner_tsvs[a] = tsv
    if raw_bam_matches:
        per_aligner_raw_bams[a] = str(raw_bam_matches[0])

log.info('TSVs: %s', list(per_aligner_tsvs.keys()))
log.info('Raw BAMs: %s', list(per_aligner_raw_bams.keys()))

import pandas as pd
for a, tsv in per_aligner_tsvs.items():
    df = pd.read_csv(tsv, sep='\t', low_memory=False)
    n_reanchor = (df['reanchor_clip_len'] > 0).sum() if 'reanchor_clip_len' in df.columns else '?'
    log.info('TSV %s: %d reads, reanchor_clip_len>0: %s (%.1f%%)',
             a, len(df), n_reanchor,
             100 * n_reanchor / len(df) if isinstance(n_reanchor, int) and len(df) else 0)

log.info('Loading genome ...')
t0 = time.perf_counter()
genome = load_genome(GENOME_FSA)
log.info('Genome loaded in %.1fs', time.perf_counter() - t0)

log.info('Loading penalty tables ...')
try:
    penalty_table = HpPenaltyTable.from_tsv(PENALTY_TSV, str_path=STR_PENALTY_TSV)
    log.info('Penalty table loaded (with STR)')
except Exception as e:
    log.warning('STR load failed (%s) — trying without', e)
    try:
        penalty_table = HpPenaltyTable.from_tsv(PENALTY_TSV)
        log.info('Penalty table loaded (without STR)')
    except Exception as e2:
        log.warning('Penalty table load failed (%s) — using None', e2)
        penalty_table = None

# ── Per-aligner HP scoring, individually timed ───────────────────────────
# Time each aligner's HP scoring loop separately by calling score_hp_edits
# directly, then record per-read rate.
log.info('')
log.info('=' * 70)
log.info('PER-ALIGNER HP SCORING (with _apply_reanchor_from_clip_len active)')
log.info('=' * 70)

from rectify.core.consensus.corrected_consensus import _score_hp_distances_for_aligner

aligner_times = {}
for a in ALIGNERS:
    if a not in per_aligner_tsvs or a not in per_aligner_raw_bams:
        log.warning('Skipping %s (missing TSV or raw BAM)', a)
        continue
    log.info('Scoring %s ...', a)
    t_start = time.perf_counter()
    args = (a, per_aligner_raw_bams[a], str(per_aligner_tsvs[a]), genome, penalty_table, True, None)
    _, result = _score_hp_distances_for_aligner(args)
    n_reads = len(result)
    wall = time.perf_counter() - t_start
    us_per_read = 1e6 * wall / n_reads if n_reads else 0
    log.info('  %s: %.1f s, %d reads, %.1f µs/read', a, wall, n_reads, us_per_read)
    aligner_times[a] = (wall, n_reads, us_per_read)

# ── Full lazy path (merge + write consensus BAM) ─────────────────────────
log.info('')
log.info('=' * 70)
log.info('FULL LAZY PATH (merge_corrected_tsvs + write_corrected_consensus_bam)')
log.info('=' * 70)

lazy_out_tsv = SCRATCH / 'lazy_merged_corrected_reads.tsv'
lazy_out_bam = SCRATCH / 'lazy_consensus.bam'

tracemalloc.start()
t_lazy_start = time.perf_counter()

merge_corrected_tsvs(
    per_aligner_tsvs=per_aligner_tsvs,
    output_tsv=lazy_out_tsv,
    per_aligner_raw_bams=per_aligner_raw_bams,
    genome=genome,
    penalty_table=penalty_table,
    lazy_scoring_workers=1,
    strict_lazy_identity=True,
)

t_merge_done = time.perf_counter()
log.info('Merge done in %.1fs', t_merge_done - t_lazy_start)

write_stats = write_corrected_consensus_bam(
    per_aligner_raw_bams=per_aligner_raw_bams,
    per_aligner_tsvs=per_aligner_tsvs,
    merged_tsv=lazy_out_tsv,
    output_bam=lazy_out_bam,
    genome=genome,
    threads=4,
    strict=True,
)

t_lazy_end = time.perf_counter()
_, peak_lazy = tracemalloc.get_traced_memory()
tracemalloc.stop()

wall_lazy = t_lazy_end - t_lazy_start
peak_lazy_mb = peak_lazy / 1e6
wall_bam_write = t_lazy_end - t_merge_done
log.info('BAM write done in %.1fs', wall_bam_write)
log.info('LAZY TOTAL: wall=%.1fs  peak_rss=%.0f MB  stats=%s', wall_lazy, peak_lazy_mb, write_stats)

# ── Results ───────────────────────────────────────────────────────────────
results = f"""
===============================================================
Phase D Reanchor Optimization Timing
Dataset:  wt_tfiiib_rep3 / chunk_008
Commit:   _apply_reanchor_from_clip_len (O(CIGAR-ops) shortcut)
Job:      {JOB_ID}
Date:     2026-05-22

Per-aligner HP scoring:
| Aligner    | Wall (s) | Reads  | µs/read | Baseline µs/read | Speedup |
|------------|----------|--------|---------|-----------------|---------|"""

baseline = {'minimap2': 3.8, 'mapPacBio': 9.3, 'gapmm2': 7.1}
for a in ALIGNERS:
    if a in aligner_times:
        wall, n_reads, us = aligner_times[a]
        base = baseline.get(a, 0)
        spd = base / us if us > 0 else 0
        results += f"\n| {a:<10} | {wall:>8.1f} | {n_reads:>6} | {us:>7.1f} | {base:>15.1f} | {spd:>7.2f}x |"

results += f"""

Full lazy path:
| Stage               | Wall (s) |
|---------------------|----------|
| Merge (HP scoring)  | {t_merge_done - t_lazy_start:>8.1f} |
| BAM write           | {wall_bam_write:>8.1f} |
| Total               | {wall_lazy:>8.1f} |
| Peak RSS            | {peak_lazy_mb:>7.0f} MB |

Phase C baseline (2026-05-21, job 25662145):
  Total lazy: 495.1 s  Peak RSS: 286 MB
===============================================================
"""
print(results)
log.info(results)

result_file = TIMING_DIR / f'phase_d_reanchor_timing_{JOB_ID}.txt'
result_file.write_text(results)
log.info('Results written to %s', result_file)
PYEOF

echo "Done: $(date)"
