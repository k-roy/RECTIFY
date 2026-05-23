#!/bin/bash
# Phase C lazy corrected consensus timing
# Compares:
#   Lazy path: merge_corrected_tsvs with per_aligner_raw_bams (lazy HP scoring from raw BAMs + TSVs)
#   Full path: rectify correct --write-corrected-bam per aligner, then merge with per_aligner_corrected_bams
#
# Dataset: wt_tfiiib_rep3 / chunk_008 (3 aligners: minimap2, mapPacBio, gapmm2)
# Note: lazy path run FIRST to avoid warm-cache bias for full path measurement.

#SBATCH --job-name=rectify_lazy_timing_phaseC
#SBATCH --partition=larsms,owners
#SBATCH --account=larsms
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/timing/phase_c_lazy_timing_%j.out
#SBATCH --error=/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/timing/phase_c_lazy_timing_%j.err

# ── AMD node MKL fix (must be before conda) ───────────────────────────────
export MKL_DEBUG_CPU_TYPE=5

# ── Env setup (source bashrc BEFORE set -u to avoid 256term.sh unbound trap) ──
export PS1=""
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify

# Now safe to enable strict mode
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
ANNOTATION="${RECTIFY_SRC}/rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"
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

# Verify inputs exist
for p in "$OUTDIR/aligner_chunks/minimap2/chunk_${CHUNK}/corrected_reads.tsv" \
          "$OUTDIR/aligner_chunks/mapPacBio/chunk_${CHUNK}/corrected_reads.tsv" \
          "$OUTDIR/aligner_chunks/gapmm2/chunk_${CHUNK}/corrected_reads.tsv" \
          "$GENOME_FSA" "$PENALTY_TSV" "$STR_PENALTY_TSV"; do
    [ -f "$p" ] || { echo "ERROR: missing: $p" >&2; exit 1; }
done

SCRATCH_WORK="${SCRATCH:-/tmp}/rectify_timing_${SLURM_JOB_ID:-$$}"
mkdir -p "$SCRATCH_WORK"
trap "rm -rf '$SCRATCH_WORK'" EXIT

# ─────────────────────────────────────────────────────────────────────────────
# PYTHON TIMING SCRIPT
# ─────────────────────────────────────────────────────────────────────────────
cd "$RECTIFY_SRC"
$PYTHON - << PYEOF
import sys, time, tracemalloc, logging, os, subprocess, shutil
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

sys.path.insert(0, '/oak/stanford/groups/larsms/Users/kevinroy/software/rectify')

from rectify.core.consensus.corrected_consensus import merge_corrected_tsvs
from rectify.utils.genome import load_genome
from rectify.core.splice.hp_penalty import HpPenaltyTable

OUTDIR    = Path('/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/wt_tfiiib_rep3/chunks')
TIMING_DIR= Path('/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/timing')
SCRATCH   = Path('${SCRATCH_WORK}')
CHUNK     = '008'
JOB_ID    = '${SLURM_JOB_ID:-local}'
GENOME_FSA = '/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz'
PENALTY_TSV = '/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv'
STR_PENALTY_TSV = '/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/str_penalty_scores.tsv'
POOL_PKL  = str(OUTDIR / 'junction_pool.pkl')
SCAN_PKL  = str(OUTDIR / 'rescue_scan.pkl') if (OUTDIR / 'rescue_scan.pkl').exists() else None

ALIGNERS = ['minimap2', 'mapPacBio', 'gapmm2']
RECTIFY_SRC = '/oak/stanford/groups/larsms/Users/kevinroy/software/rectify'
ANNOTATION  = f'{RECTIFY_SRC}/rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz'
MERGED_BAMS = Path(str(OUTDIR) + '/merged_bams')

# Build per-aligner paths
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
        log.info('Raw BAM for %s: %s', a, raw_bam_matches[0])

log.info('TSVs found: %s', list(per_aligner_tsvs.keys()))
log.info('Raw BAMs found: %s', list(per_aligner_raw_bams.keys()))

# Row counts
for a, tsv in per_aligner_tsvs.items():
    import pandas as pd
    df = pd.read_csv(tsv, sep='\t', low_memory=False)
    log.info('TSV %s: %d rows, %d reads', a, len(df), df['read_id'].nunique())

# Load genome and penalty table (required for accurate HP-ED scoring)
log.info('Loading genome from %s ...', GENOME_FSA)
t_load = time.perf_counter()
genome = load_genome(GENOME_FSA)
log.info('Genome loaded: %d chroms in %.1fs', len(genome), time.perf_counter() - t_load)

log.info('Loading penalty tables ...')
try:
    penalty_table = HpPenaltyTable.from_tsv(PENALTY_TSV, str_path=STR_PENALTY_TSV)
    log.info('Penalty table loaded: HpPenaltyTable (with STR table)')
except Exception as e:
    log.warning('Could not load penalty table with STR: %s — trying without STR', e)
    try:
        penalty_table = HpPenaltyTable.from_tsv(PENALTY_TSV)
        log.info('Penalty table loaded: HpPenaltyTable (without STR table)')
    except Exception as e2:
        log.warning('Could not load penalty table: %s — using None (flat fallback)', e2)
        penalty_table = None

# ─── PATH 1: LAZY (raw BAMs + TSVs, corrections applied in memory) ───────
log.info('')
log.info('=' * 70)
log.info('PATH 1: LAZY — merge_corrected_tsvs with per_aligner_raw_bams')
log.info('=' * 70)

lazy_out_tsv = SCRATCH / 'lazy_merged_corrected_reads.tsv'
lazy_out_bam = SCRATCH / 'lazy_consensus.bam'

tracemalloc.start()
t0 = time.perf_counter()

merge_corrected_tsvs(
    per_aligner_tsvs=per_aligner_tsvs,
    output_tsv=lazy_out_tsv,
    per_aligner_raw_bams=per_aligner_raw_bams,
    genome=genome,
    penalty_table=penalty_table,
    lazy_scoring_workers=1,
    strict_lazy_identity=True,
)

t_merge_lazy = time.perf_counter()
log.info('Lazy merge done in %.1fs', t_merge_lazy - t0)

# Write consensus BAM via lazy path
from rectify.core.consensus.corrected_consensus import write_corrected_consensus_bam
write_stats = write_corrected_consensus_bam(
    per_aligner_raw_bams=per_aligner_raw_bams,
    per_aligner_tsvs=per_aligner_tsvs,
    merged_tsv=lazy_out_tsv,
    output_bam=lazy_out_bam,
    genome=genome,
    threads=4,
    strict=True,
)
t1 = time.perf_counter()
_, peak_lazy = tracemalloc.get_traced_memory()
tracemalloc.stop()

wall_lazy = t1 - t0
peak_lazy_mb = peak_lazy / 1e6
log.info('LAZY: wall=%.1fs  peak_rss=%.0f MB  stats=%s', wall_lazy, peak_lazy_mb, write_stats)

# ─── PATH 2: FULL (materialize per-aligner corrected BAMs, then score) ───
log.info('')
log.info('=' * 70)
log.info('PATH 2: FULL — rectify correct --write-corrected-bam x3, then merge with corrected BAMs')
log.info('=' * 70)

full_scratch = SCRATCH / 'full'
full_scratch.mkdir(exist_ok=True)

PYTHON_BIN = '/home/groups/larsms/users/kevinroy/anaconda3/envs/rectify/bin/python'
if not Path(PYTHON_BIN).exists():
    PYTHON_BIN = '/home/groups/larsms/users/kevinroy/anaconda3/bin/python'

# Verify PYTHON_BIN has rectify available
import importlib.util
has_rectify = importlib.util.find_spec('rectify') is not None
log.info('rectify importable: %s', has_rectify)

tracemalloc.start()
t2 = time.perf_counter()

per_aligner_corrected_bams = {}
for aligner in ALIGNERS:
    raw_bam = per_aligner_raw_bams.get(aligner)
    if not raw_bam:
        log.warning('No raw BAM for %s — skipping', aligner)
        continue
    corrected_bam = str(full_scratch / f'{aligner}_corrected.bam')
    log.info('Running rectify correct for %s ...', aligner)
    t_aln_start = time.perf_counter()

    SCAN_ARG = ['--variant-scan-cache', SCAN_PKL] if SCAN_PKL else []
    cmd = [
        PYTHON_BIN, '-m', 'rectify', 'correct',
        raw_bam,
        '--genome', GENOME_FSA,
        '--annotation', ANNOTATION,
        '--junction-pool-cache', POOL_PKL,
        '--aligner-bams', f'mapPacBio:{MERGED_BAMS}/wt_tfiiib_rep3_trimmed.mapPacBio.bam',
        '--aligner-bams', f'minimap2:{MERGED_BAMS}/wt_tfiiib_rep3_trimmed.minimap2.bam',
        '--aligner-bams', f'gapmm2:{MERGED_BAMS}/wt_tfiiib_rep3_trimmed.gapmm2.bam',
        '--threads', '4',
        '--streaming',
        '-o', str(full_scratch / f'{aligner}_corrected_reads.tsv'),
        '--write-corrected-bam', corrected_bam,
    ] + SCAN_ARG
    log.info('CMD: %s', ' '.join(cmd[:8]) + ' ...')
    result = subprocess.run(
        cmd, capture_output=True, text=True,
        cwd=RECTIFY_SRC,
        env={**os.environ, 'MKL_DEBUG_CPU_TYPE': '5'},
    )
    if result.returncode != 0:
        log.error('rectify correct failed for %s: %s', aligner, result.stderr[-2000:])
        raise RuntimeError(f'rectify correct failed for {aligner}')
    log.info('  %s correct done in %.1fs', aligner, time.perf_counter() - t_aln_start)
    per_aligner_corrected_bams[aligner] = corrected_bam

log.info('Materialization done. %d corrected BAMs ready.', len(per_aligner_corrected_bams))

# Use the TSVs written by rectify correct (not the pre-existing chunk TSVs,
# since we want the timing to cover the full materialized-bam pipeline)
per_aligner_tsvs_full = {
    a: full_scratch / f'{a}_corrected_reads.tsv'
    for a in per_aligner_corrected_bams
    if (full_scratch / f'{a}_corrected_reads.tsv').exists()
}

# Merge using corrected BAMs for HP-ED scoring
full_out_tsv = full_scratch / 'full_merged_corrected_reads.tsv'
full_out_bam = full_scratch / 'full_consensus.bam'

merge_corrected_tsvs(
    per_aligner_tsvs=per_aligner_tsvs_full if per_aligner_tsvs_full else per_aligner_tsvs,
    output_tsv=full_out_tsv,
    per_aligner_corrected_bams=per_aligner_corrected_bams,
    genome=genome,
    penalty_table=penalty_table,
)

# Write consensus BAM
write_corrected_consensus_bam(
    per_aligner_raw_bams=per_aligner_raw_bams,
    per_aligner_tsvs=per_aligner_tsvs_full if per_aligner_tsvs_full else per_aligner_tsvs,
    merged_tsv=full_out_tsv,
    output_bam=full_out_bam,
    genome=genome,
    threads=4,
    strict=False,
)

t3 = time.perf_counter()
_, peak_full = tracemalloc.get_traced_memory()
tracemalloc.stop()

wall_full = t3 - t2
peak_full_mb = peak_full / 1e6
log.info('FULL: wall=%.1fs  peak_rss=%.0f MB', wall_full, peak_full_mb)

# ─── RESULTS ─────────────────────────────────────────────────────────────
speedup_wall = wall_full / wall_lazy if wall_lazy > 0 else float('inf')
speedup_rss  = peak_full_mb / peak_lazy_mb if peak_lazy_mb > 0 else float('inf')

# Count reads in the representative TSV (header row excluded)
_tsv_sample = list(per_aligner_tsvs.values())[0] if per_aligner_tsvs else None
_n_rows = (sum(1 for _ in open(_tsv_sample)) - 1) if _tsv_sample else '?'

results = f"""
===============================================================
Phase C Timing Results
Dataset: wt_tfiiib_rep3 / chunk_008
Aligners: minimap2, mapPacBio, gapmm2
Rows in minimap2 TSV: ~{_n_rows}
Note: lazy path timed cold; full path timed warm (filesystem cache)
      Real lazy advantage is >= the measured speedup.
Job: {JOB_ID}
Date: 2026-05-21

| Path    | Wall (s)     | Peak RSS (MB) |
|---------|-------------|---------------|
| Lazy    | {wall_lazy:>10.1f} | {peak_lazy_mb:>12.0f}  |
| Full    | {wall_full:>10.1f} | {peak_full_mb:>12.0f}  |
| Speedup | {speedup_wall:>10.1f}x | {speedup_rss:>11.1f}x  |
===============================================================
"""
print(results)
log.info(results)

# Write to timing dir
result_file = TIMING_DIR / f'phase_c_timing_results_{JOB_ID}.txt'
result_file.write_text(results)
log.info('Results written to %s', result_file)

PYEOF

echo "Done: $(date)"
