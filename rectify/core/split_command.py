"""
Split Command for RECTIFY.

Splits a FASTQ/FASTQ.GZ file into N equal-sized chunks at read boundaries.
Used to parallelize alignment across a scheduler array job (SLURM, UGE/SGE, PBS).

Chunk count is determined automatically from a target reads-per-chunk (default 500 000),
so a 1 M-read dataset gets ~4 chunks while a 50 M-read dataset gets ~100 chunks.
The explicit -n/--n-chunks flag overrides auto-sizing.

Two array scripts are generated (per-aligner resource split):
  run_array_mapPacBio.sh  — N_CHUNKS tasks, 8 cores / 40 GB  (Java -Xmx32g)
  run_array_others.sh     — N_CHUNKS × N_OTHER_ALIGNERS tasks, 2 cores / 16 GB
  run_merge_consensus.sh  — single job, 16 cores / 32 GB
  submit_pipeline.sh      — chains all three as scheduler dependencies

Typical workflow:

  # 1. Auto-split (N computed from read count ÷ 500k)
  rectify split reads.fastq.gz \\
      -o /scratch/chunks \\
      --genome genome.fa.gz --annotation genes.gff.gz \\
      --generate-slurm

  # 2. Submit the chained pipeline
  bash /scratch/chunks/submit_pipeline.sh

Author: Kevin R. Roy
"""

import argparse
import gzip
import json
import logging
import math
import stat
import sys
from pathlib import Path
from typing import List, Optional, Tuple

logger = logging.getLogger(__name__)

# ── Resource specs per aligner (empirically validated on HPC clusters) ──────
# mapPacBio uses -Xmx32g; heap does not shrink with input size.
MPB_CORES   = 8
MPB_MEM_GB  = 40      # 32 GB heap + 8 GB OS/JVM overhead
MPB_TIME    = '6:00:00'

# Other aligners on 1/16-sized chunks are lightweight.
OTHER_CORES  = 2
OTHER_MEM_GB = 16
OTHER_TIME   = '4:00:00'

MERGE_CORES  = 16
MERGE_MEM_GB = 32
MERGE_TIME   = '4:00:00'

# Auto-sizing defaults
DEFAULT_TARGET_READS = 500_000
MIN_CHUNKS = 4
MAX_CHUNKS = 64

# Standard aligners (mapPacBio excluded — handled by separate array)
OTHER_ALIGNERS_DEFAULT = ['minimap2', 'gapmm2', 'uLTRA', 'deSALT']
JUNCTION_ALIGNERS = {'uLTRA', 'deSALT'}


def create_split_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Create split subcommand parser."""
    parser = subparsers.add_parser(
        'split',
        help='Split a FASTQ/FASTQ.GZ into chunks for parallel array-job alignment',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Examples:
  # Auto-size chunks (500k reads/chunk), generate SLURM scripts
  rectify split reads.fastq.gz -o chunks/ --genome g.fa.gz --annotation a.gff.gz --generate-slurm

  # Fixed 20 chunks
  rectify split reads.fastq.gz -n 20 -o chunks/

  # Generate UGE/SGE scripts
  rectify split reads.fastq.gz -o chunks/ --generate-slurm --scheduler uge --uge-queue long.q

  # Dry run — show chunk count without writing
  rectify split reads.fastq.gz --dry-run
        """
    )

    parser.add_argument('reads', type=Path, help='Input FASTQ or FASTQ.GZ file')

    size_group = parser.add_mutually_exclusive_group()
    size_group.add_argument(
        '-n', '--n-chunks',
        type=int, default=None, metavar='N',
        help='Number of chunks (overrides --target-reads-per-chunk)'
    )
    size_group.add_argument(
        '--target-reads-per-chunk',
        type=int, default=DEFAULT_TARGET_READS, metavar='READS',
        help=f'Auto-size chunks to approximately this many reads each '
             f'(min {MIN_CHUNKS}, max {MAX_CHUNKS})'
    )

    parser.add_argument(
        '-o', '--output-dir', type=Path, required=True,
        help='Output directory for chunk FASTQ files and generated scripts'
    )
    parser.add_argument(
        '--prefix', default='',
        help='Prefix for chunk file names (default: derived from input filename)'
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help='Count reads and print chunk sizes without writing files'
    )
    parser.add_argument(
        '--generate-slurm', action='store_true',
        help='Generate scheduler array scripts alongside the chunk files'
    )

    # Script generation options
    script_group = parser.add_argument_group('Script generation options (used with --generate-slurm)')
    script_group.add_argument(
        '--genome', type=Path,
        help='Reference genome path (written into generated scripts)'
    )
    script_group.add_argument(
        '--annotation', type=Path,
        help='Annotation GFF/GTF path (written into generated scripts)'
    )
    script_group.add_argument(
        '--other-aligners', nargs='+',
        default=OTHER_ALIGNERS_DEFAULT,
        metavar='ALIGNER',
        help='Non-mapPacBio aligners to include in the others array'
    )
    script_group.add_argument(
        '--skip-map-pacbio', action='store_true',
        help='Omit the mapPacBio array script (e.g. not installed)'
    )
    script_group.add_argument(
        '--python-path',
        default='python',
        help='Explicit path to Python interpreter (default: python on PATH)'
    )
    script_group.add_argument(
        '--rectify-src',
        default='.',
        help='Path to RECTIFY source checkout (used as working directory; default: current dir)'
    )

    # Scheduler selection
    sched_group = parser.add_argument_group('Scheduler options')
    sched_group.add_argument(
        '--scheduler', choices=['slurm', 'uge', 'pbs'], default='slurm',
        help='Target scheduler for generated script headers'
    )
    # SLURM
    sched_group.add_argument('--slurm-partition', default=None,
                             help='SLURM partition(s)')
    sched_group.add_argument('--slurm-account', default=None,
                             help='SLURM account')
    # UGE/SGE
    sched_group.add_argument('--uge-queue', default='long.q',
                             help='UGE/SGE queue name (-q)')
    sched_group.add_argument('--uge-pe', default='smp',
                             help='UGE/SGE parallel environment (-pe)')
    # PBS/Torque
    sched_group.add_argument('--pbs-queue', default='workq',
                             help='PBS queue name (-q)')

    parser.add_argument('--verbose', action='store_true', help='Verbose logging')
    return parser


def _open_fastq(path: Path):
    """Open a FASTQ or FASTQ.GZ for reading, yielding lines as str."""
    if path.suffix == '.gz' or str(path).endswith('.fastq.gz') or str(path).endswith('.fq.gz'):
        return gzip.open(str(path), 'rt')
    return open(str(path), 'r')


def count_reads(path: Path) -> int:
    """Count reads in a FASTQ/FASTQ.GZ file."""
    n = 0
    with _open_fastq(path) as fh:
        for line in fh:
            n += 1
    if n % 4 != 0:
        raise ValueError(
            f"{path}: line count {n} is not a multiple of 4 — file may be truncated"
        )
    return n // 4


def compute_n_chunks(n_reads: int, target_per_chunk: int = DEFAULT_TARGET_READS) -> int:
    """
    Compute chunk count from target reads-per-chunk.

    Returns max(MIN_CHUNKS, min(MAX_CHUNKS, ceil(n_reads / target_per_chunk))).
    """
    n = math.ceil(n_reads / target_per_chunk)
    return max(MIN_CHUNKS, min(MAX_CHUNKS, n))


def split_fastq(
    input_path: Path,
    output_dir: Path,
    n_chunks: int,
    prefix: str = '',
) -> List[Path]:
    """
    Split a FASTQ/FASTQ.GZ into n_chunks equal files via round-robin interleaving.

    Round-robin gives each chunk an even distribution of read lengths even when
    reads are sorted by genomic position (which correlates with length in RNA-seq).

    Returns list of output chunk paths.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    if not prefix:
        name = input_path.name
        for suffix in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
            if name.endswith(suffix):
                prefix = name[:-len(suffix)]
                break
        else:
            prefix = input_path.stem

    chunk_paths = [
        output_dir / f"{prefix}_chunk_{k:03d}_of_{n_chunks:03d}.fastq.gz"
        for k in range(n_chunks)
    ]
    chunk_files = [gzip.open(str(p), 'wt') for p in chunk_paths]

    try:
        read_idx = 0
        with _open_fastq(input_path) as fh:
            while True:
                header = fh.readline()
                if not header:
                    break
                seq  = fh.readline()
                plus = fh.readline()
                qual = fh.readline()
                if not qual:
                    raise ValueError(
                        f"Truncated FASTQ record at read {read_idx} in {input_path}"
                    )
                dest = chunk_files[read_idx % n_chunks]
                dest.write(header)
                dest.write(seq)
                dest.write(plus)
                dest.write(qual)
                read_idx += 1

        logger.info(
            "Split %d reads into %d chunks (~%d reads/chunk)",
            read_idx, n_chunks, read_idx // n_chunks,
        )
    finally:
        for fh in chunk_files:
            fh.close()

    return chunk_paths


# ── Script header builders ───────────────────────────────────────────────────

def _slurm_headers(
    job_name: str,
    n_tasks: int,
    cores: int,
    mem_gb: int,
    time: str,
    log_pattern: str,
    partition: Optional[str] = None,
    account: Optional[str] = None,
    is_array: bool = True,
) -> str:
    array_line = f'#SBATCH --array=0-{n_tasks - 1}' if is_array else ''
    return f"""\
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --account={account}
#SBATCH --cpus-per-task={cores}
#SBATCH --mem={mem_gb}G
#SBATCH --time={time}
{array_line}
#SBATCH --output={log_pattern}.out
#SBATCH --error={log_pattern}.err"""


def _uge_headers(
    job_name: str,
    n_tasks: int,
    cores: int,
    mem_gb: int,
    time: str,
    log_dir: str,
    queue: str = 'long.q',
    pe: str = 'smp',
    is_array: bool = True,
) -> str:
    mem_per_slot = math.ceil(mem_gb / cores)
    array_line = f'#$ -t 1-{n_tasks}' if is_array else ''   # UGE is 1-based
    return f"""\
#$ -N {job_name}
#$ -q {queue}
#$ -pe {pe} {cores}
#$ -l h_vmem={mem_per_slot}G
#$ -l h_rt={time}
{array_line}
#$ -o {log_dir}/
#$ -e {log_dir}/
#$ -cwd
#$ -j n"""


def _pbs_headers(
    job_name: str,
    n_tasks: int,
    cores: int,
    mem_gb: int,
    time: str,
    log_dir: str,
    queue: str = 'workq',
    is_array: bool = True,
) -> str:
    array_line = f'#PBS -J 0-{n_tasks - 1}' if is_array else ''   # PBS 0-based
    return f"""\
#PBS -N {job_name}
#PBS -q {queue}
#PBS -l select=1:ncpus={cores}:mem={mem_gb}gb
#PBS -l walltime={time}
{array_line}
#PBS -o {log_dir}/
#PBS -e {log_dir}/"""


def _scheduler_headers(
    scheduler: str, job_name: str, n_tasks: int, cores: int, mem_gb: int,
    time: str, log_dir: str, log_pattern: str,
    slurm_partition: Optional[str] = None, slurm_account: Optional[str] = None,
    uge_queue: str = 'long.q', uge_pe: str = 'smp',
    pbs_queue: str = 'workq', is_array: bool = True,
) -> str:
    if scheduler == 'slurm':
        return _slurm_headers(job_name, n_tasks, cores, mem_gb, time,
                              log_pattern, slurm_partition, slurm_account, is_array)
    elif scheduler == 'uge':
        return _uge_headers(job_name, n_tasks, cores, mem_gb, time,
                            log_dir, uge_queue, uge_pe, is_array)
    elif scheduler == 'pbs':
        return _pbs_headers(job_name, n_tasks, cores, mem_gb, time,
                            log_dir, pbs_queue, is_array)
    else:
        raise ValueError(f"Unknown scheduler: {scheduler}")


# ── Task-ID normalisation shim (identical across all generated scripts) ──────
_TASK_ID_SHIM = """\
# ── Scheduler-agnostic task ID and CPU normalisation ──────────────────────
# SLURM: 0-based.  UGE/SGE: 1-based (subtract 1).  PBS: 0-based.
if   [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    RECTIFY_TASK_ID=$SLURM_ARRAY_TASK_ID
    RECTIFY_CPUS=${SLURM_CPUS_PER_TASK:-__CORES__}
elif [ -n "${SGE_TASK_ID:-}" ]; then
    RECTIFY_TASK_ID=$(( SGE_TASK_ID - 1 ))
    RECTIFY_CPUS=${NSLOTS:-__CORES__}
elif [ -n "${PBS_ARRAY_INDEX:-}" ]; then
    RECTIFY_TASK_ID=$(( PBS_ARRAY_INDEX - 1 ))
    RECTIFY_CPUS=${PBS_NUM_PPN:-__CORES__}
else
    echo "ERROR: no scheduler array task variable found" >&2
    exit 1
fi
# ──────────────────────────────────────────────────────────────────────────"""


def _thread_limits_block(cores_var: str = '$RECTIFY_CPUS') -> str:
    return f"""\
export OMP_NUM_THREADS={cores_var}
export OPENBLAS_NUM_THREADS={cores_var}
export MKL_NUM_THREADS={cores_var}
export LOKY_MAX_CPU_COUNT={cores_var}
export NUMEXPR_MAX_THREADS={cores_var}"""


# ── Script body builders ─────────────────────────────────────────────────────

def _mpb_array_body(
    n_chunks: int, sample_prefix: str, output_dir: Path, genome: str,
    annot: str, python_path: str, rectify_src: str,
) -> str:
    """Body of the mapPacBio array script (after headers)."""
    shim = _TASK_ID_SHIM.replace('__CORES__', str(MPB_CORES))
    limits = _thread_limits_block()
    log_dir = output_dir / 'logs'
    return f"""#!/bin/bash
# RECTIFY — mapPacBio chunk-alignment array
# Generated by: rectify split --generate-slurm
# {n_chunks} tasks  ×  1 aligner  =  {n_chunks} total tasks

{_SCHEDULER_HEADER_PLACEHOLDER_MPB}

set -euo pipefail

export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

{shim}

{limits}

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"
GENOME="{genome}"
ANNOT="{annot}"
N_CHUNKS={n_chunks}

CHUNK_PAD=$(printf "%03d" $RECTIFY_TASK_ID)
CHUNK_FASTQ="$OUTDIR/{sample_prefix}_chunk_${{CHUNK_PAD}}_of_{n_chunks:03d}.fastq.gz"
CHUNK_OUTDIR="$OUTDIR/aligner_chunks/mapPacBio/chunk_${{CHUNK_PAD}}"
mkdir -p "$CHUNK_OUTDIR"

echo "Python:  $($PYTHON --version)"
echo "Host:    $(hostname)"
echo "Task:    $RECTIFY_TASK_ID  Chunk: $CHUNK_PAD  CPUs: $RECTIFY_CPUS"
echo "Input:   $CHUNK_FASTQ"
echo "Start:   $(date)"

[ -f "$CHUNK_FASTQ" ] || {{ echo "ERROR: chunk FASTQ not found: $CHUNK_FASTQ" >&2; exit 1; }}

cd "$RECTIFY_SRC"
$PYTHON -m rectify align \\
    "$CHUNK_FASTQ" \\
    --genome "$GENOME" \\
    --annotation "$ANNOT" \\
    --aligners mapPacBio \\
    --no-consensus \\
    --sort --index \\
    -t "$RECTIFY_CPUS" \\
    -o "$CHUNK_OUTDIR" \\
    --verbose

echo "Done: $(date)"
ls -lh "$CHUNK_OUTDIR"/*.bam 2>/dev/null || true
"""


def _others_array_body(
    n_chunks: int, other_aligners: List[str], sample_prefix: str,
    output_dir: Path, genome: str, annot: str, python_path: str, rectify_src: str,
) -> str:
    """Body of the others (minimap2/gapmm2/uLTRA/deSALT) array script."""
    n_other = len(other_aligners)
    n_tasks = n_chunks * n_other
    shim = _TASK_ID_SHIM.replace('__CORES__', str(OTHER_CORES))
    limits = _thread_limits_block()
    aligner_arr = ' '.join(f'"{a}"' for a in other_aligners)
    return f"""#!/bin/bash
# RECTIFY — non-mapPacBio chunk-alignment array
# Generated by: rectify split --generate-slurm
# {n_chunks} chunks × {n_other} aligners = {n_tasks} total tasks
# Aligner order: {', '.join(f'{i}={a}' for i, a in enumerate(other_aligners))}

{_SCHEDULER_HEADER_PLACEHOLDER_OTHERS}

set -euo pipefail

export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

{shim}

{limits}

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"
GENOME="{genome}"
ANNOT="{annot}"
N_CHUNKS={n_chunks}
ALIGNERS=({aligner_arr})

ALIGNER_IDX=$(( RECTIFY_TASK_ID / N_CHUNKS ))
CHUNK_IDX_NUM=$(( RECTIFY_TASK_ID % N_CHUNKS ))
CHUNK_PAD=$(printf "%03d" $CHUNK_IDX_NUM)
ALIGNER="${{ALIGNERS[$ALIGNER_IDX]}}"

CHUNK_FASTQ="$OUTDIR/{sample_prefix}_chunk_${{CHUNK_PAD}}_of_{n_chunks:03d}.fastq.gz"
CHUNK_OUTDIR="$OUTDIR/aligner_chunks/${{ALIGNER}}/chunk_${{CHUNK_PAD}}"
mkdir -p "$CHUNK_OUTDIR"

echo "Python:  $($PYTHON --version)"
echo "Host:    $(hostname)"
echo "Task:    $RECTIFY_TASK_ID  Aligner: $ALIGNER  Chunk: $CHUNK_PAD  CPUs: $RECTIFY_CPUS"
echo "Input:   $CHUNK_FASTQ"
echo "Start:   $(date)"

[ -f "$CHUNK_FASTQ" ] || {{ echo "ERROR: chunk FASTQ not found: $CHUNK_FASTQ" >&2; exit 1; }}

cd "$RECTIFY_SRC"

# Run alignment with retry (deSALT occasionally segfaults on specific read batches;
# a simple retry usually succeeds on the second or third attempt).
MAX_RETRIES=3
for _try in $(seq 1 $MAX_RETRIES); do
    # uLTRA and deSALT are junction-mode aligners
    if [ "$ALIGNER" = "uLTRA" ] || [ "$ALIGNER" = "deSALT" ]; then
        $PYTHON -m rectify align \\
            "$CHUNK_FASTQ" \\
            --genome "$GENOME" \\
            --annotation "$ANNOT" \\
            --aligners none \\
            --junction-aligners "$ALIGNER" \\
            --no-consensus \\
            --sort --index \\
            -t "$RECTIFY_CPUS" \\
            -o "$CHUNK_OUTDIR" \\
            --verbose && break
    else
        $PYTHON -m rectify align \\
            "$CHUNK_FASTQ" \\
            --genome "$GENOME" \\
            --annotation "$ANNOT" \\
            --aligners "$ALIGNER" \\
            --no-consensus \\
            --sort --index \\
            -t "$RECTIFY_CPUS" \\
            -o "$CHUNK_OUTDIR" \\
            --verbose && break
    fi
    echo "Attempt $_try/$MAX_RETRIES failed for $ALIGNER chunk $CHUNK_PAD — retrying in 10s" >&2
    sleep 10
    if [ $_try -eq $MAX_RETRIES ]; then
        echo "ERROR: $ALIGNER chunk $CHUNK_PAD failed after $MAX_RETRIES attempts" >&2
        exit 1
    fi
done

echo "Done: $(date)"
ls -lh "$CHUNK_OUTDIR"/*.bam 2>/dev/null || true
"""


def _merge_consensus_body(
    n_chunks: int, all_aligners: List[str], sample_prefix: str,
    output_dir: Path, genome: str, annot: str, python_path: str, rectify_src: str,
) -> str:
    """Body of the merge + consensus script."""
    limits = _thread_limits_block('$MERGE_CPUS')
    aligner_arr = ' '.join(f'"{a}"' for a in all_aligners)
    return f"""#!/bin/bash
# RECTIFY — merge per-aligner chunk BAMs and run consensus selection
# Run after both array jobs complete.

{_SCHEDULER_HEADER_PLACEHOLDER_MERGE}

set -euo pipefail

export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

MERGE_CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-{MERGE_CORES}}}}}}}

{limits}

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"
GENOME="{genome}"
ANNOT="{annot}"
SAMPLE="{sample_prefix}"
N_CHUNKS={n_chunks}
ALIGNERS=({aligner_arr})

MERGED_DIR="$OUTDIR/merged_bams"
mkdir -p "$MERGED_DIR"

echo "Host:    $(hostname)"
echo "CPUs:    $MERGE_CPUS"
echo "Start:   $(date)"

# ── Merge per-chunk BAMs into one per-aligner BAM ─────────────────────────
BAM_ARGS=""
for ALIGNER in "${{ALIGNERS[@]}}"; do
    echo ""
    echo "=== Merging $ALIGNER ==="
    CHUNK_BAMS=()
    for (( k=0; k<N_CHUNKS; k++ )); do
        CHUNK_PAD=$(printf "%03d" $k)
        CHUNK_OUTDIR="$OUTDIR/aligner_chunks/$ALIGNER/chunk_$CHUNK_PAD"
        BAM=$(ls "$CHUNK_OUTDIR"/*.${{ALIGNER}}.bam 2>/dev/null | head -1 || true)
        if [ -z "$BAM" ]; then
            echo "  WARNING: no BAM in $CHUNK_OUTDIR" >&2
        else
            CHUNK_BAMS+=("$BAM")
        fi
    done

    if [ ${{#CHUNK_BAMS[@]}} -eq 0 ]; then
        echo "  SKIPPING $ALIGNER — no chunk BAMs found" >&2
        continue
    fi

    MERGED_BAM="$MERGED_DIR/${{SAMPLE}}.${{ALIGNER}}.bam"
    echo "  Merging ${{#CHUNK_BAMS[@]}} chunks → $MERGED_BAM"
    samtools merge -f -@ "$MERGE_CPUS" "$MERGED_BAM" "${{CHUNK_BAMS[@]}}"
    samtools sort  -@ "$MERGE_CPUS" -o "${{MERGED_BAM%.bam}}.sorted.bam" "$MERGED_BAM"
    mv "${{MERGED_BAM%.bam}}.sorted.bam" "$MERGED_BAM"
    samtools index -@ "$MERGE_CPUS" "$MERGED_BAM"
    BAM_ARGS="$BAM_ARGS ${{ALIGNER}}:${{MERGED_BAM}}"
    echo "  ✓ $(du -h "$MERGED_BAM" | cut -f1)  $MERGED_BAM"
done

echo ""
echo "=== rectify consensus ==="
CONSENSUS_DIR="$OUTDIR/consensus"
mkdir -p "$CONSENSUS_DIR"

cd "$RECTIFY_SRC"
$PYTHON -m rectify consensus \\
    --genome "$GENOME" \\
    --annotation "$ANNOT" \\
    --prefix "${{SAMPLE}}" \\
    -o "$CONSENSUS_DIR" \\
    --verbose \\
    $BAM_ARGS

echo ""
echo "=== Consensus outputs ==="
ls -lh "$CONSENSUS_DIR/"*.bam 2>/dev/null
echo "Done: $(date)"
"""


def _make_submit_script(
    scheduler: str,
    mpb_script: Optional[Path],
    others_script: Path,
    merge_script: Path,
    output_dir: Path,
) -> str:
    """Generate a submit_pipeline.sh that chains jobs with dependencies."""
    if scheduler == 'slurm':
        lines = ['#!/bin/bash', '# Submit the full chunked-alignment pipeline (SLURM)', '']
        if mpb_script:
            lines += [
                f'MPB_JOB=$(sbatch {mpb_script} | awk \'{{print $4}}\')',
                'echo "mapPacBio array: $MPB_JOB"',
            ]
        lines += [
            f'OTHERS_JOB=$(sbatch {others_script} | awk \'{{print $4}}\')',
            'echo "Others array:    $OTHERS_JOB"',
        ]
        if mpb_script:
            dep = '--dependency=afterok:${MPB_JOB}:${OTHERS_JOB}'
        else:
            dep = '--dependency=afterok:${OTHERS_JOB}'
        lines += [
            f'MERGE_JOB=$(sbatch {dep} {merge_script} | awk \'{{print $4}}\')',
            'echo "Merge+consensus: $MERGE_JOB"',
            'echo ""',
            'echo "Pipeline submitted. Monitor with: squeue -u $USER"',
        ]
    elif scheduler == 'uge':
        lines = ['#!/bin/bash', '# Submit the full chunked-alignment pipeline (UGE/SGE)', '']
        if mpb_script:
            lines += [
                f'MPB_JOB=$(qsub {mpb_script} | awk \'{{print $3}}\')',
                'echo "mapPacBio array: $MPB_JOB"',
            ]
        lines += [
            f'OTHERS_JOB=$(qsub {others_script} | awk \'{{print $3}}\')',
            'echo "Others array:    $OTHERS_JOB"',
        ]
        if mpb_script:
            hold = '-hold_jid ${MPB_JOB},${OTHERS_JOB}'
        else:
            hold = '-hold_jid ${OTHERS_JOB}'
        lines += [
            f'MERGE_JOB=$(qsub {hold} {merge_script} | awk \'{{print $3}}\')',
            'echo "Merge+consensus: $MERGE_JOB"',
            'echo ""',
            'echo "Pipeline submitted. Monitor with: qstat -u $USER"',
        ]
    else:  # pbs
        lines = ['#!/bin/bash', '# Submit the full chunked-alignment pipeline (PBS/Torque)', '']
        if mpb_script:
            lines += [
                f'MPB_JOB=$(qsub {mpb_script})',
                'echo "mapPacBio array: $MPB_JOB"',
            ]
        lines += [
            f'OTHERS_JOB=$(qsub {others_script})',
            'echo "Others array:    $OTHERS_JOB"',
        ]
        if mpb_script:
            dep = '-W depend=afterok:${MPB_JOB}:${OTHERS_JOB}'
        else:
            dep = '-W depend=afterok:${OTHERS_JOB}'
        lines += [
            f'MERGE_JOB=$(qsub {dep} {merge_script})',
            'echo "Merge+consensus: $MERGE_JOB"',
            'echo ""',
            'echo "Pipeline submitted. Monitor with: qstat -u $USER"',
        ]
    return '\n'.join(lines) + '\n'


# Placeholder tokens replaced after header generation
_SCHEDULER_HEADER_PLACEHOLDER_MPB    = '__SCHED_HEADER_MPB__'
_SCHEDULER_HEADER_PLACEHOLDER_OTHERS = '__SCHED_HEADER_OTHERS__'
_SCHEDULER_HEADER_PLACEHOLDER_MERGE  = '__SCHED_HEADER_MERGE__'


def _generate_scripts(
    args: argparse.Namespace,
    n_chunks: int,
    sample_prefix: str,
) -> None:
    """Generate all scheduler scripts for the chunked alignment pipeline."""
    output_dir = args.output_dir.resolve()
    log_dir    = output_dir / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)

    genome_str = str(args.genome.resolve()) if args.genome else '# REQUIRED: set --genome'
    annot_str  = str(args.annotation.resolve()) if args.annotation else '# REQUIRED: set --annotation'
    python_path = getattr(args, 'python_path', 'python')
    rectify_src = getattr(args, 'rectify_src', '.')

    scheduler = getattr(args, 'scheduler', 'slurm')
    slurm_partition = getattr(args, 'slurm_partition', None)
    slurm_account   = getattr(args, 'slurm_account', None)
    uge_queue       = getattr(args, 'uge_queue', 'long.q')
    uge_pe          = getattr(args, 'uge_pe', 'smp')
    pbs_queue       = getattr(args, 'pbs_queue', 'workq')

    other_aligners = list(getattr(args, 'other_aligners', OTHER_ALIGNERS_DEFAULT))
    skip_mpb       = getattr(args, 'skip_map_pacbio', False)

    n_other_tasks = n_chunks * len(other_aligners)
    log_pat = str(log_dir / '%A_%a')   # SLURM pattern; ignored by UGE/PBS headers

    sched_kwargs = dict(
        slurm_partition=slurm_partition, slurm_account=slurm_account,
        uge_queue=uge_queue, uge_pe=uge_pe, pbs_queue=pbs_queue,
    )

    # ── mapPacBio array ────────────────────────────────────────────────────
    mpb_script_path: Optional[Path] = None
    if not skip_mpb:
        mpb_headers = _scheduler_headers(
            scheduler, f'{sample_prefix}_mpb', n_chunks,
            MPB_CORES, MPB_MEM_GB, MPB_TIME,
            str(log_dir), log_pat, is_array=True, **sched_kwargs,
        )
        mpb_body = _mpb_array_body(
            n_chunks, sample_prefix, output_dir, genome_str, annot_str,
            python_path, rectify_src,
        )
        mpb_body = mpb_body.replace(_SCHEDULER_HEADER_PLACEHOLDER_MPB, mpb_headers)
        mpb_script_path = output_dir / 'run_array_mapPacBio.sh'
        mpb_script_path.write_text(mpb_body)
        mpb_script_path.chmod(mpb_script_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)

    # ── Others array ───────────────────────────────────────────────────────
    others_headers = _scheduler_headers(
        scheduler, f'{sample_prefix}_align', n_other_tasks,
        OTHER_CORES, OTHER_MEM_GB, OTHER_TIME,
        str(log_dir), log_pat, is_array=True, **sched_kwargs,
    )
    others_body = _others_array_body(
        n_chunks, other_aligners, sample_prefix, output_dir, genome_str, annot_str,
        python_path, rectify_src,
    )
    others_body = others_body.replace(_SCHEDULER_HEADER_PLACEHOLDER_OTHERS, others_headers)
    others_script_path = output_dir / 'run_array_others.sh'
    others_script_path.write_text(others_body)
    others_script_path.chmod(others_script_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)

    # ── Merge + consensus ─────────────────────────────────────────────────
    all_aligners = (['mapPacBio'] if not skip_mpb else []) + other_aligners
    merge_headers = _scheduler_headers(
        scheduler, f'{sample_prefix}_merge', 1,
        MERGE_CORES, MERGE_MEM_GB, MERGE_TIME,
        str(log_dir), str(log_dir / '%j'), is_array=False, **sched_kwargs,
    )
    merge_body = _merge_consensus_body(
        n_chunks, all_aligners, sample_prefix, output_dir, genome_str, annot_str,
        python_path, rectify_src,
    )
    merge_body = merge_body.replace(_SCHEDULER_HEADER_PLACEHOLDER_MERGE, merge_headers)
    merge_script_path = output_dir / 'run_merge_consensus.sh'
    merge_script_path.write_text(merge_body)
    merge_script_path.chmod(merge_script_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)

    # ── Submission wrapper ─────────────────────────────────────────────────
    submit_body = _make_submit_script(
        scheduler, mpb_script_path, others_script_path, merge_script_path, output_dir,
    )
    submit_path = output_dir / 'submit_pipeline.sh'
    submit_path.write_text(submit_body)
    submit_path.chmod(submit_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)

    # ── Summary ───────────────────────────────────────────────────────────
    logger.info("Generated scheduler scripts (%s):", scheduler.upper())
    if mpb_script_path:
        logger.info("  mapPacBio array (%d tasks, %d cores, %dG):  %s",
                    n_chunks, MPB_CORES, MPB_MEM_GB, mpb_script_path)
    logger.info("  Others array    (%d tasks, %d cores, %dG):  %s",
                n_other_tasks, OTHER_CORES, OTHER_MEM_GB, others_script_path)
    logger.info("  Merge+consensus (1 job,   %d cores, %dG):  %s",
                MERGE_CORES, MERGE_MEM_GB, merge_script_path)
    logger.info("  Submit chain:  bash %s", submit_path)


def generate_alignment_scripts(
    n_chunks: int,
    sample_prefix: str,
    output_dir: Path,
    genome: Optional[Path],
    annotation: Optional[Path],
    python_path: str,
    rectify_src: str,
    scheduler: str = 'slurm',
    other_aligners: Optional[List[str]] = None,
    skip_map_pacbio: bool = False,
    slurm_partition: Optional[str] = None,
    slurm_account: Optional[str] = None,
    uge_queue: str = 'long.q',
    uge_pe: str = 'smp',
    pbs_queue: str = 'workq',
) -> dict:
    """
    Public API: generate alignment array scripts for a given chunk configuration.

    Returns dict with keys: 'mpb' (Path or None), 'others' (Path), 'merge' (Path), 'submit' (Path).
    All scripts are written to output_dir.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    args = argparse.Namespace(
        output_dir=output_dir,
        genome=genome,
        annotation=annotation,
        python_path=python_path,
        rectify_src=rectify_src,
        scheduler=scheduler,
        other_aligners=other_aligners or OTHER_ALIGNERS_DEFAULT,
        skip_map_pacbio=skip_map_pacbio,
        slurm_partition=slurm_partition,
        slurm_account=slurm_account,
        uge_queue=uge_queue,
        uge_pe=uge_pe,
        pbs_queue=pbs_queue,
    )
    _generate_scripts(args, n_chunks, sample_prefix)

    return {
        'mpb':    output_dir / 'run_array_mapPacBio.sh' if not skip_map_pacbio else None,
        'others': output_dir / 'run_array_others.sh',
        'merge':  output_dir / 'run_merge_consensus.sh',
        'submit': output_dir / 'submit_pipeline.sh',
    }


def run_split(args: argparse.Namespace) -> int:
    """Run split command."""
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    if not args.reads.exists():
        logger.error("Reads file not found: %s", args.reads)
        return 1

    # ── Determine chunk count ─────────────────────────────────────────────
    if args.n_chunks is not None:
        n_chunks = args.n_chunks
        if n_chunks < 1:
            logger.error("--n-chunks must be >= 1")
            return 1
        logger.info("Counting reads in %s ...", args.reads)
        n_reads = count_reads(args.reads)
    else:
        target = getattr(args, 'target_reads_per_chunk', DEFAULT_TARGET_READS)
        logger.info("Counting reads in %s ...", args.reads)
        n_reads = count_reads(args.reads)
        n_chunks = compute_n_chunks(n_reads, target)
        logger.info("  Auto-sized: %d reads / %d target = %d chunks "
                    "(min %d, max %d)", n_reads, target, n_chunks, MIN_CHUNKS, MAX_CHUNKS)

    chunk_size = math.ceil(n_reads / n_chunks)
    logger.info("  Total reads:       %s", f"{n_reads:,}")
    logger.info("  Chunks:            %d", n_chunks)
    logger.info("  ~reads per chunk:  %s", f"{chunk_size:,}")

    if args.dry_run:
        return 0

    # ── Derive sample prefix ──────────────────────────────────────────────
    prefix = args.prefix
    if not prefix:
        name = args.reads.name
        for suffix in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
            if name.endswith(suffix):
                prefix = name[:-len(suffix)]
                break
        else:
            prefix = args.reads.stem

    # ── Split FASTQ ───────────────────────────────────────────────────────
    logger.info("Splitting %s into %d chunks → %s", args.reads, n_chunks, args.output_dir)
    chunk_paths = split_fastq(
        input_path=args.reads,
        output_dir=args.output_dir,
        n_chunks=n_chunks,
        prefix=prefix,
    )

    # ── Write manifest ────────────────────────────────────────────────────
    manifest = {
        'input': str(args.reads.resolve()),
        'n_reads': n_reads,
        'n_chunks': n_chunks,
        'reads_per_chunk': chunk_size,
        'chunks': [str(p) for p in chunk_paths],
    }
    manifest_path = args.output_dir / 'chunks_manifest.json'
    with open(str(manifest_path), 'w') as fh:
        json.dump(manifest, fh, indent=2)
    logger.info("Manifest: %s", manifest_path)

    # ── Generate scheduler scripts ────────────────────────────────────────
    if getattr(args, 'generate_slurm', False):
        _generate_scripts(args, n_chunks, prefix)

    return 0


def run(args: argparse.Namespace):
    """Entry point for CLI."""
    sys.exit(run_split(args))
