#!/bin/bash
# Phase E: production-representative lazy timing (lazy_scoring_workers=3)
# + fine-grained BAM write sub-step breakdown.
#
# Phase D exposed two measurement gaps:
#   1. lazy_scoring_workers=1 in Phase C/D is not production config
#      (production uses lazy_scoring_workers=THREADS=8, effectively 3
#       parallel aligner workers).
#   2. The BAM write step (136s, 27% of total) has never been profiled.
#
# This script measures:
#   (a) Full lazy path with lazy_scoring_workers=3 (production config)
#   (b) BAM write sub-steps: TSV load, correction table load, unsorted
#       write loop (3 aligner passes), pysam.sort, pysam.index
#   (c) cProfile top-30 hotspots over the full lazy path
#
# Phase D baseline (workers=1):
#   Merge (HP scoring): 355.6 s
#   BAM write:          136.6 s
#   Total:              492.2 s

#SBATCH --job-name=rectify_phase_e_timing
#SBATCH --partition=larsms,owners
#SBATCH --account=larsms
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/timing/phase_e_parallel_timing_%j.out
#SBATCH --error=/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/timing/phase_e_parallel_timing_%j.err

export MKL_DEBUG_CPU_TYPE=5

export PS1=""
source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
conda activate rectify

set -euo pipefail

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8
export LOKY_MAX_CPU_COUNT=8
export NUMEXPR_MAX_THREADS=8

PYTHON="/home/groups/larsms/users/kevinroy/anaconda3/bin/python"
RECTIFY_SRC="/oak/stanford/groups/larsms/Users/kevinroy/software/rectify"
OUTDIR="/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/wt_tfiiib_rep3/chunks"
TIMING_DIR="/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/timing"
CHUNK="008"
GENOME_FSA="${RECTIFY_SRC}/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
PENALTY_TSV="${RECTIFY_SRC}/rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv"
STR_PENALTY_TSV="${RECTIFY_SRC}/rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/str_penalty_scores.tsv"

mkdir -p "$TIMING_DIR"

echo "Host:   $(hostname)"
echo "Job:    ${SLURM_JOB_ID:-local}"
echo "Chunk:  $CHUNK"
echo "Start:  $(date)"
echo "CPUs:   ${SLURM_CPUS_PER_TASK:-8}"

SCRATCH_WORK="${SCRATCH:-/tmp}/rectify_phase_e_timing_${SLURM_JOB_ID:-$$}"
mkdir -p "$SCRATCH_WORK"
trap "rm -rf '$SCRATCH_WORK'" EXIT

cd "$RECTIFY_SRC"
$PYTHON - << 'PYEOF'
import sys, time, cProfile, pstats, io, tracemalloc, logging
from pathlib import Path
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

import os
SCRATCH    = Path(os.environ.get('SCRATCH_WORK', '/tmp/phase_e'))
SCRATCH.mkdir(parents=True, exist_ok=True)
JOB_ID     = os.environ.get('SLURM_JOB_ID', 'local')
TIMING_DIR = Path('/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/timing')

sys.path.insert(0, '/oak/stanford/groups/larsms/Users/kevinroy/software/rectify')

from rectify.core.consensus.corrected_consensus import (
    merge_corrected_tsvs,
    write_corrected_consensus_bam,
    _normalized_correction_lookup,
    _normalize_bam_read_name,
    _read_corrected_tsv_or_manifest,
)
from rectify.utils.genome import load_genome
from rectify.core.splice.hp_penalty import HpPenaltyTable

OUTDIR     = Path('/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/wt_tfiiib_rep3/chunks')
GENOME_FSA = '/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz'
PENALTY_TSV = '/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv'
STR_PENALTY_TSV = '/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/str_penalty_scores.tsv'
CHUNK = '008'
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

log.info('TSVs: %s  BAMs: %s', list(per_aligner_tsvs.keys()), list(per_aligner_raw_bams.keys()))

log.info('Loading genome ...')
genome = load_genome(GENOME_FSA)
log.info('Loading penalty tables ...')
try:
    penalty_table = HpPenaltyTable.from_tsv(PENALTY_TSV, str_path=STR_PENALTY_TSV)
except Exception:
    penalty_table = HpPenaltyTable.from_tsv(PENALTY_TSV)

# ── (a) Full lazy path with workers=3 ────────────────────────────────────
log.info('')
log.info('=' * 70)
log.info('(a) FULL LAZY PATH — lazy_scoring_workers=3 (production config)')
log.info('=' * 70)

lazy_out_tsv = SCRATCH / 'lazy_merged.tsv'
lazy_out_bam = SCRATCH / 'lazy_consensus.bam'

tracemalloc.start()
t0 = time.perf_counter()

merge_corrected_tsvs(
    per_aligner_tsvs=per_aligner_tsvs,
    output_tsv=lazy_out_tsv,
    per_aligner_raw_bams=per_aligner_raw_bams,
    genome=genome,
    penalty_table=penalty_table,
    lazy_scoring_workers=3,
    strict_lazy_identity=True,
)
t_merge = time.perf_counter()

write_stats = write_corrected_consensus_bam(
    per_aligner_raw_bams=per_aligner_raw_bams,
    per_aligner_tsvs=per_aligner_tsvs,
    merged_tsv=lazy_out_tsv,
    output_bam=lazy_out_bam,
    genome=genome,
    threads=8,
    strict=True,
)
t_done = time.perf_counter()
_, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

wall_merge = t_merge - t0
wall_bam   = t_done - t_merge
wall_total = t_done - t0

log.info('Merge:  %.1f s', wall_merge)
log.info('BAM write: %.1f s', wall_bam)
log.info('Total:  %.1f s  peak_rss=%.0f MB', wall_total, peak / 1e6)

# ── (b) BAM write sub-step breakdown ─────────────────────────────────────
log.info('')
log.info('=' * 70)
log.info('(b) BAM WRITE SUB-STEP BREAKDOWN (second pass, same merged TSV)')
log.info('=' * 70)

import pysam as _pysam
from rectify.core.bam.bam_writer import (
    _load_corrections_from_tsv,
    apply_corrected_edits_to_read,
)

t_b0 = time.perf_counter()
merged_df = _read_corrected_tsv_or_manifest(lazy_out_tsv)
t_b1 = time.perf_counter()
log.info('  TSV load:     %.2f s (%d rows)', t_b1 - t_b0, len(merged_df))

read_to_aligner = {
    _normalize_bam_read_name(str(row.read_id)): str(row.winning_aligner)
    for row in merged_df[['read_id', 'winning_aligner']].drop_duplicates().itertuples(index=False)
}
aligner_to_reads = defaultdict(set)
for rid, al in read_to_aligner.items():
    aligner_to_reads[al].add(rid)
t_b2 = time.perf_counter()
log.info('  Dict build:   %.2f s (%d winners)', t_b2 - t_b1, len(read_to_aligner))

correction_tables = {}
for aligner, tsv_path in per_aligner_tsvs.items():
    corrections = _load_corrections_from_tsv(str(tsv_path))
    correction_tables[aligner] = (corrections, _normalized_correction_lookup(corrections))
t_b3 = time.perf_counter()
log.info('  Corr tables:  %.2f s', t_b3 - t_b2)

header_bam = next(Path(p) for p in per_aligner_raw_bams.values() if Path(p).exists())
unsorted_bam2 = SCRATCH / 'lazy_consensus2.unsorted.bam'
sorted_bam2   = SCRATCH / 'lazy_consensus2.bam'

t_stream_start = time.perf_counter()
aligner_stream_times = {}
written = 0
with _pysam.AlignmentFile(str(header_bam), 'rb') as h:
    header = h.header
with _pysam.AlignmentFile(str(unsorted_bam2), 'wb', header=header) as out:
    for aligner, winning_ids in aligner_to_reads.items():
        raw_bam = per_aligner_raw_bams.get(aligner)
        corrections, norm_corrections = correction_tables[aligner]
        t_al = time.perf_counter()
        n_written = 0
        with _pysam.AlignmentFile(str(raw_bam), 'rb') as bam_in:
            for read in bam_in:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                rid = _normalize_bam_read_name(read.query_name or '')
                if rid not in winning_ids:
                    continue
                correction = corrections.get(read.query_name or '')
                if correction is None:
                    from rectify.core.consensus.corrected_consensus import _lookup_read_correction
                    correction = _lookup_read_correction(read, corrections, norm_corrections)
                if correction is None:
                    continue
                apply_corrected_edits_to_read(read, correction, genome)
                out.write(read)
                n_written += 1
        aligner_stream_times[aligner] = (time.perf_counter() - t_al, n_written)
        written += n_written
t_b4 = time.perf_counter()
log.info('  Unsorted write: %.2f s (%d reads)', t_b4 - t_stream_start, written)
for al, (al_t, al_n) in aligner_stream_times.items():
    log.info('    %s: %.2f s, %d writes, %.1f µs/write', al, al_t, al_n, 1e6 * al_t / al_n if al_n else 0)

_pysam.sort('-@', '8', '-o', str(sorted_bam2), str(unsorted_bam2))
t_b5 = time.perf_counter()
log.info('  pysam.sort:   %.2f s', t_b5 - t_b4)

_pysam.index(str(sorted_bam2))
t_b6 = time.perf_counter()
log.info('  pysam.index:  %.2f s', t_b6 - t_b5)
log.info('  BAM write TOTAL: %.2f s', t_b6 - t_b0)

# ── (c) cProfile top-30 over a fresh merge+write ─────────────────────────
log.info('')
log.info('=' * 70)
log.info('(c) cProfile — top-30 cumulative, full lazy path (workers=3)')
log.info('=' * 70)

lazy_out_tsv2 = SCRATCH / 'lazy_merged_prof.tsv'
lazy_out_bam2 = SCRATCH / 'lazy_consensus_prof.bam'

pr = cProfile.Profile()
pr.enable()

merge_corrected_tsvs(
    per_aligner_tsvs=per_aligner_tsvs,
    output_tsv=lazy_out_tsv2,
    per_aligner_raw_bams=per_aligner_raw_bams,
    genome=genome,
    penalty_table=penalty_table,
    lazy_scoring_workers=3,
    strict_lazy_identity=True,
)
write_corrected_consensus_bam(
    per_aligner_raw_bams=per_aligner_raw_bams,
    per_aligner_tsvs=per_aligner_tsvs,
    merged_tsv=lazy_out_tsv2,
    output_bam=lazy_out_bam2,
    genome=genome,
    threads=8,
    strict=True,
)

pr.disable()

buf = io.StringIO()
ps = pstats.Stats(pr, stream=buf).sort_stats('cumulative')
ps.print_stats(30)
profile_text = buf.getvalue()
log.info('cProfile output:\n%s', profile_text)

# ── Results ───────────────────────────────────────────────────────────────
results = f"""
===============================================================
Phase E: Production-Representative Lazy Timing
Dataset:  wt_tfiiib_rep3 / chunk_008
Config:   lazy_scoring_workers=3, BAM write threads=8
Job:      {JOB_ID}
Date:     2026-05-22

(a) Full lazy path (workers=3):
| Stage      | Wall (s) |
|------------|----------|
| Merge      | {wall_merge:>8.1f} |
| BAM write  | {wall_bam:>8.1f} |
| Total      | {wall_total:>8.1f} |
| Peak RSS   | {peak/1e6:>7.0f} MB |

Phase D baseline (workers=1):  492.2 s total

(b) BAM write sub-steps:
  TSV load:       {t_b1-t_b0:.2f} s
  Dict build:     {t_b2-t_b1:.2f} s
  Corr tables:    {t_b3-t_b2:.2f} s
  Unsorted write: {t_b4-t_stream_start:.2f} s ({written} reads)
    Per-aligner:
{"".join(f"    {al}: {al_t:.2f}s  {al_n} writes  {1e6*al_t/al_n if al_n else 0:.1f} µs/write" + chr(10) for al, (al_t, al_n) in aligner_stream_times.items())}  pysam.sort:     {t_b5-t_b4:.2f} s
  pysam.index:    {t_b6-t_b5:.2f} s
  BAM write TOTAL:{t_b6-t_b0:.2f} s

(c) cProfile top-30 (cumulative):
{profile_text[:3000]}
===============================================================
"""
print(results)

result_file = TIMING_DIR / f'phase_e_parallel_timing_{JOB_ID}.txt'
result_file.write_text(results)
log.info('Results written to %s', result_file)
PYEOF

echo "Done: $(date)"
