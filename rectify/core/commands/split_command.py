"""
Split Command for RECTIFY.

Splits a FASTQ/FASTQ.GZ or BAM file into N equal-sized chunks at read boundaries.
Used to parallelize alignment across a scheduler array job (SLURM, UGE/SGE, PBS).
BAM inputs are auto-converted to FASTQ (preserving pt:i: and other aux tags via
`samtools fastq -T '*'`) before chunking, so a Dorado-aligned DRS BAM can go
straight through the chunked alignment pipeline without a manual pre-step.

Chunk count is determined automatically from a target reads-per-chunk (default 50 000),
so a 1 M-read dataset gets ~20 chunks while a 50 M-read dataset gets ~1000 chunks
(capped at MAX_CHUNKS).
The explicit -n/--n-chunks flag overrides auto-sizing.

Correct-first pipeline scripts are generated:
  run_array_mapPacBio.sh          — N_CHUNKS tasks, 8 cores / 40 GB  (Java -Xmx32g)
  run_array_others.sh             — N_CHUNKS × N_OTHER_ALIGNERS tasks, 2 cores / 16 GB
  run_merge_aligners.sh           — single job; merges per-chunk BAMs → per-aligner BAMs (for prescan)
  run_prescan.sh                  — single job; rectify prescan → junction_pool.pkl
  run_array_correct_{aligner}.sh  — one per aligner, N_CHUNKS tasks; uses --junction-pool-cache
  run_array_chunk_merge.sh        — N_CHUNKS tasks; per-chunk TSV merge → consensus BAM → polya BAM
  run_final_merge.sh              — single job; concat TSVs + samtools merge BAMs
  submit_pipeline.sh              — chains all scripts with correct scheduler dependencies

Typical workflow:

  # 1. Auto-split (N computed from read count ÷ 50k)
  rectify split reads.fastq.gz \\
      -o /scratch/chunks \\
      --genome genome.fa.gz --annotation genes.gff.gz \\
      --generate-slurm

  # 2. Submit the chained correct-first pipeline
  bash /scratch/chunks/submit_pipeline.sh

Author: Kevin R. Roy
"""

import argparse
import gzip
import hashlib
import json
import logging
import math
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

from rectify.core.chunking.sidecar import (
    ReadNumSidecarWriter,
    format_fastq_header_with_rn,
    sha256_file,
    split_fastq_header,
)

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
# 50k reads/chunk → ~25 min mapPacBio correction per chunk (validated 2026-04-30)
# Larger chunks (500k) caused 3–6h correction jobs; 50k keeps wall time ≤30 min.
DEFAULT_TARGET_READS = 50_000
MIN_CHUNKS = 4
MAX_CHUNKS = 500

PRESCAN_CORES       = 8
PRESCAN_MEM_GB      = 32
PRESCAN_TIME        = '1:00:00'

CORRECT_CORES       = 8
CORRECT_MEM_GB      = 48
CORRECT_TIME        = '1:30:00'   # 50k reads ~25 min; 1.5h is safe with margin

CHUNK_MERGE_CORES   = 8
CHUNK_MERGE_MEM_GB  = 32
CHUNK_MERGE_TIME    = '1:00:00'

FINAL_MERGE_CORES   = 8
FINAL_MERGE_MEM_GB  = 32
FINAL_MERGE_TIME    = '2:00:00'

# Standard aligners (mapPacBio excluded — handled by separate array)
OTHER_ALIGNERS_DEFAULT = ['minimap2', 'gapmm2', 'uLTRA', 'deSALT']
JUNCTION_ALIGNERS = {'uLTRA', 'deSALT'}

# Short-read aligners (Illumina QuantSeq / cDNA)
SHORT_READ_ALIGNERS_DEFAULT = ['bbmap', 'bwa']

# Short-read resource specs — lighter than long-read (no Java heap, smaller BAMs)
SR_ALIGN_CORRECT_CORES  = 4
SR_ALIGN_CORRECT_MEM_GB = 32
SR_ALIGN_CORRECT_TIME   = '1:30:00'

SR_FINAL_MERGE_CORES  = 8
SR_FINAL_MERGE_MEM_GB = 32
SR_FINAL_MERGE_TIME   = '0:30:00'

# Per-chunk consensus selection (parallelised SLURM array)
# Validated: chunk_000 = 101 min wall for 425k reads at 73 reads/sec
# 22× speedup vs whole-sample sequential consensus
CONSENSUS_PER_CHUNK_CORES   = 8
CONSENSUS_PER_CHUNK_MEM_GB  = 32
CONSENSUS_PER_CHUNK_TIME    = '4:00:00'

# Merge per-chunk consensus BAMs into final sample consensus BAM
MERGE_CONSENSUS_CORES   = 8
MERGE_CONSENSUS_MEM_GB  = 32
MERGE_CONSENSUS_TIME    = '2:00:00'


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

    parser.add_argument(
        'reads', type=Path,
        help=(
            'Input FASTQ, FASTQ.GZ, or BAM file. BAM inputs are converted to '
            'FASTQ before chunking; the derived FASTQ is written into '
            '--output-dir. With --drs, the BAM is run through `trim-polya` '
            'first (yielding a trimmed FASTQ + parquet metadata sidecar so '
            'Step-4 restore-softclip can be invoked later). Without --drs, '
            "the BAM is dumped via `samtools fastq -T '*'` (aux tags kept "
            'as FASTQ comments, no trimming).'
        ),
    )
    parser.add_argument(
        '-2', '--read2', type=Path, default=None,
        help=(
            'Mate-2 (R2) FASTQ for paired-end short-read data. When given, the '
            'positional argument is treated as R1 and chunking is pair-aware: '
            'both mates share one RN and land in the same chunk, emitting paired '
            'chunk FASTQs (_chunk_NNN_of_MMM_R1/_R2). Implies short-read paired '
            'mode (use with --short-read).'
        ),
    )
    parser.add_argument(
        '--drs', action='store_true', default=False,
        help=(
            'Treat a BAM input as a Dorado-aligned DRS BAM and run Step-0 '
            'poly(A)/adapter trim (`trim-polya`) before chunking. The parquet '
            'metadata is preserved at <output-dir>/<stem>_polya_trim_metadata'
            '.parquet for later use by `rectify restore-softclip`. No effect '
            'on FASTQ inputs.'
        ),
    )

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

    from rectify.data import add_organism_args, add_junction_anchor_args
    add_organism_args(script_group)
    add_junction_anchor_args(script_group)
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

    script_group.add_argument(
        '--junction-penalty-table', default='', metavar='TSV',
        help='Empirical HP junction penalty table (forwarded to correction scripts)'
    )
    script_group.add_argument(
        '--str-penalty-table', default='', metavar='TSV',
        help='STR slippage penalty table (forwarded to correction scripts)'
    )
    script_group.add_argument(
        '--junction-overhang-table', default='', metavar='TSV',
        help='Empirical junction overhang table (forwarded to chunk-merge scripts). '
             'Generated by calibrate_junction_overhang.py from a reference run. '
             'When supplied, penalises aligner wins that rely on junctions lacking '
             'sufficient flanking overhang for their intron size.'
    )

    # Short-read mode
    sr_group = parser.add_argument_group(
        'Short-read options (Illumina QuantSeq / dT-primed cDNA)',
        description=(
            'Use --short-read to generate a simplified pipeline for Illumina short reads '
            '(bbmap + bwa aligners, correct-before-consensus, all I/O on $SCRATCH). '
            'Replaces the long-read mapPacBio/minimap2/gapmm2/uLTRA/deSALT pipeline.'
        ),
    )
    sr_group.add_argument(
        '--short-read', action='store_true',
        help=(
            'Generate short-read pipeline (bbmap + bwa). '
            'Disables mapPacBio, prescan, and long-read aligners. '
            'Implies --skip-map-pacbio.'
        ),
    )
    sr_group.add_argument(
        '--dT-primed-cDNA', action='store_true', dest='dt_primed_cdna',
        help=(
            'Pass --dT-primed-cDNA to generated rectify correct calls '
            '(enables AG-mispriming module for oligo-dT cDNA synthesis artifacts). '
            'Auto-enabled when --short-read is set.'
        ),
    )
    sr_group.add_argument(
        '--oak-output-dir', type=Path, default=None, metavar='DIR',
        help=(
            'Destination directory on persistent storage (Oak/NFS) for the final '
            'merged consensus BAM. The final merge script copies the result here '
            'and cleans up $SCRATCH. Defaults to {output-dir}/final/ if not set.'
        ),
    )

    script_group.add_argument(
        '--write-per-aligner-corrected-bams',
        action='store_true',
        help='Also emit per-aligner corrected BAMs for IGV inspection '
             '(default off; lazy consensus path does not need them).',
    )

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
    write_read_num_sidecar: bool = True,
) -> List[Path]:
    """
    Split a FASTQ/FASTQ.GZ into n_chunks equal files via round-robin interleaving.

    Round-robin gives each chunk an even distribution of read lengths even when
    reads are sorted by genomic position (which correlates with length in RNA-seq).
    Each derived FASTQ header carries an ``RN:i:<read_num>`` comment tag. When
    ``write_read_num_sidecar`` is true, a sibling ``<prefix>.read_num_sidecar.parquet``
    maps those read numbers back to the original FASTQ header comment.

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
    sidecar_writer = None
    try:
        if write_read_num_sidecar:
            sidecar_writer = ReadNumSidecarWriter(
                output_dir / f"{prefix}.read_num_sidecar.parquet",
                sample_id=prefix,
                source_fastq=input_path,
                source_fastq_sha256=sha256_file(input_path),
                chunk_count=n_chunks,
            )

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
                chunk_idx = read_idx % n_chunks
                chunk_id = f"chunk_{chunk_idx:03d}_of_{n_chunks:03d}"
                original_qname, fastq_comment = split_fastq_header(header)
                if sidecar_writer is not None:
                    sidecar_writer.add(
                        read_idx, original_qname, fastq_comment, chunk_id, seq, qual
                    )
                dest = chunk_files[chunk_idx]
                dest.write(format_fastq_header_with_rn(
                    original_qname, fastq_comment, read_idx
                ))
                dest.write(seq if seq.endswith('\n') else seq + '\n')
                dest.write(plus if plus.endswith('\n') else plus + '\n')
                dest.write(qual if qual.endswith('\n') else qual + '\n')
                read_idx += 1

        logger.info(
            "Split %d reads into %d chunks (~%d reads/chunk)",
            read_idx, n_chunks, read_idx // n_chunks,
        )
    finally:
        for fh in chunk_files:
            fh.close()
        if sidecar_writer is not None:
            sidecar_writer.close()

    return chunk_paths


def _strip_mate_suffix(qname: str) -> str:
    """Drop a trailing ``/1`` or ``/2`` mate suffix from a FASTQ qname token.

    Paired-end aligners pair mates by IDENTICAL QNAME and distinguish them by
    the R1/R2 input file, so both mates must carry the same bare name. The
    ``@READID 1:N:0:...`` Casava form already shares the ``READID`` token (the
    ``1``/``2`` is in the comment, stripped by ``split_fastq_header``); only the
    older ``@READID/1`` form needs trimming.
    """
    if len(qname) >= 2 and qname[-2] == '/' and qname[-1] in ('1', '2'):
        return qname[:-2]
    return qname


def _derive_paired_prefix(input_r1: Path) -> str:
    """Derive a sample prefix from an R1 filename, stripping the fastq ext and
    a trailing ``_R1`` / ``_1`` mate token (``A549_rep1_R1.fastq.gz`` → ``A549_rep1``)."""
    name = input_r1.name
    for suffix in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    else:
        name = input_r1.stem
    for mate_tok in ('_R1', '_r1', '.R1', '_1'):
        if name.endswith(mate_tok):
            return name[:-len(mate_tok)]
    return name


def split_fastq_paired(
    input_r1: Path,
    input_r2: Path,
    output_dir: Path,
    n_chunks: int,
    prefix: str = '',
    write_read_num_sidecar: bool = True,
) -> Tuple[List[Path], List[Path]]:
    """Split paired-end FASTQs into n_chunks paired files via round-robin over PAIRS.

    RECTIFY's chunking infra is single-end; this is the pair-aware variant. The
    invariants that keep proper-pairs intact across chunking:

      - Round-robin is over the PAIR index, not the record index. Both mates of
        pair ``p`` share ``read_num = p`` and land in the SAME chunk
        (``chunk_idx = p % n_chunks``) at the same local position. The legacy
        single-end formula ``read_num = local_idx * n_chunks + chunk_idx`` is
        therefore unchanged per-pair, so ``reconstruct_posthoc_sidecar_from_chunks``
        works on the R1 chunks without modification.
      - Both mate records in a chunk carry the SAME bare QNAME (``/1``,``/2``
        stripped) plus the SAME ``RN:i:<p>`` tag, so paired aligners pair them by
        identical QNAME and the RN-keyed consensus merge groups them natively.
      - The read-num sidecar stores ONE row per pair (keyed by ``read_num``); its
        seq/qual fingerprint covers R1. R1/R2 pairing integrity is structural
        (shared pair index → same chunk + local position), not fingerprinted.

    A mismatch between mate qnames raises ``ValueError`` — this catches an
    R1/R2 desync (the exact failure mode round-robin would otherwise scatter
    silently) at split time rather than deep in alignment.

    Returns ``(r1_chunk_paths, r2_chunk_paths)``.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    if not prefix:
        prefix = _derive_paired_prefix(input_r1)

    r1_paths = [
        output_dir / f"{prefix}_chunk_{k:03d}_of_{n_chunks:03d}_R1.fastq.gz"
        for k in range(n_chunks)
    ]
    r2_paths = [
        output_dir / f"{prefix}_chunk_{k:03d}_of_{n_chunks:03d}_R2.fastq.gz"
        for k in range(n_chunks)
    ]
    r1_files = [gzip.open(str(p), 'wt') for p in r1_paths]
    r2_files = [gzip.open(str(p), 'wt') for p in r2_paths]
    sidecar_writer = None
    try:
        if write_read_num_sidecar:
            # source_fastq records R1; the combined sha binds both mate files so a
            # swapped/truncated R2 is detectable in provenance.
            combined_sha = hashlib.sha256(
                f"{sha256_file(input_r1)}:{sha256_file(input_r2)}".encode()
            ).hexdigest()
            sidecar_writer = ReadNumSidecarWriter(
                output_dir / f"{prefix}.read_num_sidecar.parquet",
                sample_id=prefix,
                source_fastq=input_r1,
                source_fastq_sha256=combined_sha,
                chunk_count=n_chunks,
            )

        pair_idx = 0
        with _open_fastq(input_r1) as f1, _open_fastq(input_r2) as f2:
            while True:
                h1 = f1.readline()
                h2 = f2.readline()
                if not h1 and not h2:
                    break
                if not h1 or not h2:
                    raise ValueError(
                        f"R1/R2 record count mismatch near pair {pair_idx}: "
                        f"{input_r1} and {input_r2} have unequal lengths"
                    )
                s1, p1, q1 = f1.readline(), f1.readline(), f1.readline()
                s2, p2, q2 = f2.readline(), f2.readline(), f2.readline()
                if not q1 or not q2:
                    raise ValueError(
                        f"Truncated FASTQ record at pair {pair_idx} "
                        f"in {input_r1} / {input_r2}"
                    )

                qn1, comment1 = split_fastq_header(h1)
                qn2, comment2 = split_fastq_header(h2)
                base1 = _strip_mate_suffix(qn1)
                base2 = _strip_mate_suffix(qn2)
                if base1 != base2:
                    raise ValueError(
                        f"Mate qname mismatch at pair {pair_idx}: R1 {qn1!r} vs "
                        f"R2 {qn2!r} (bare {base1!r} != {base2!r}). R1/R2 FASTQs are "
                        f"out of sync; chunking would scatter mates."
                    )

                chunk_idx = pair_idx % n_chunks
                chunk_id = f"chunk_{chunk_idx:03d}_of_{n_chunks:03d}"
                if sidecar_writer is not None:
                    sidecar_writer.add(
                        pair_idx, base1, comment1, chunk_id, s1, q1
                    )

                d1 = r1_files[chunk_idx]
                d1.write(format_fastq_header_with_rn(base1, comment1, pair_idx))
                d1.write(s1 if s1.endswith('\n') else s1 + '\n')
                d1.write(p1 if p1.endswith('\n') else p1 + '\n')
                d1.write(q1 if q1.endswith('\n') else q1 + '\n')

                d2 = r2_files[chunk_idx]
                d2.write(format_fastq_header_with_rn(base1, comment2, pair_idx))
                d2.write(s2 if s2.endswith('\n') else s2 + '\n')
                d2.write(p2 if p2.endswith('\n') else p2 + '\n')
                d2.write(q2 if q2.endswith('\n') else q2 + '\n')

                pair_idx += 1

        logger.info(
            "Split %d read pairs into %d paired chunks (~%d pairs/chunk)",
            pair_idx, n_chunks, pair_idx // n_chunks if n_chunks else 0,
        )
    finally:
        for fh in r1_files:
            fh.close()
        for fh in r2_files:
            fh.close()
        if sidecar_writer is not None:
            sidecar_writer.close()

    return r1_paths, r2_paths


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
#$ -l h_data={mem_per_slot}G
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
export MKL_DEBUG_CPU_TYPE=5
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

PYTHON="{python_path}"
export PATH="$(dirname "$PYTHON"):$HOME/bin:$HOME/.rectify/bin:$PATH"

{shim}

{limits}

RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"
GENOME="{genome}"
ANNOT="{annot}"
N_CHUNKS={n_chunks}

CHUNK_PAD=$(printf "%03d" $RECTIFY_TASK_ID)
CHUNK_FASTQ="$OUTDIR/{sample_prefix}_chunk_${{CHUNK_PAD}}_of_{n_chunks:03d}.fastq.gz"
FINAL_OUTDIR="$OUTDIR/aligner_chunks/mapPacBio/chunk_${{CHUNK_PAD}}"
mkdir -p "$FINAL_OUTDIR"

echo "Python:  $($PYTHON --version)"
echo "Host:    $(hostname)"
echo "Task:    $RECTIFY_TASK_ID  Chunk: $CHUNK_PAD  CPUs: $RECTIFY_CPUS"
echo "Input:   $CHUNK_FASTQ"
echo "Start:   $(date)"

[ -f "$CHUNK_FASTQ" ] || {{ echo "ERROR: chunk FASTQ not found: $CHUNK_FASTQ" >&2; exit 1; }}

# Skip if alignment BAM already exists (safe resume after preemption or requeue)
if compgen -G "$FINAL_OUTDIR/*.bam" > /dev/null 2>&1; then
    echo "Alignment BAM already exists at $FINAL_OUTDIR — skipping: $(date)"
    ls -lh "$FINAL_OUTDIR"/*.bam
    exit 0
fi

# Stage output to fast local scratch to avoid shared-filesystem I/O contention.
# Falls back to /tmp if $SCRATCH is not set (non-HPC or single-node use).
SCRATCH_DIR="${{SCRATCH:-/tmp}}/rectify_mpb_${{SLURM_JOB_ID:-$$}}_${{RECTIFY_TASK_ID}}"
mkdir -p "$SCRATCH_DIR"
trap "rm -rf '$SCRATCH_DIR'" EXIT

cd "$RECTIFY_SRC"
$PYTHON -m rectify align \\
    "$CHUNK_FASTQ" \\
    --genome "$GENOME" \\
    --annotation "$ANNOT" \\
    --aligners mapPacBio \\
    --no-consensus \\
    --sort --index \\
    -t "$RECTIFY_CPUS" \\
    -o "$SCRATCH_DIR" \\
    --verbose

echo "Rsync scratch → output dir: $(date)"
rsync -a "$SCRATCH_DIR/" "$FINAL_OUTDIR/"
echo "Done: $(date)"
ls -lh "$FINAL_OUTDIR"/*.bam 2>/dev/null || true
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

PYTHON="{python_path}"
export PATH="$(dirname "$PYTHON"):$HOME/bin:$HOME/.rectify/bin:$PATH"

{shim}

{limits}

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
FINAL_OUTDIR="$OUTDIR/aligner_chunks/${{ALIGNER}}/chunk_${{CHUNK_PAD}}"
mkdir -p "$FINAL_OUTDIR"

echo "Python:  $($PYTHON --version)"
echo "Host:    $(hostname)"
echo "Task:    $RECTIFY_TASK_ID  Aligner: $ALIGNER  Chunk: $CHUNK_PAD  CPUs: $RECTIFY_CPUS"
echo "Input:   $CHUNK_FASTQ"
echo "Start:   $(date)"

[ -f "$CHUNK_FASTQ" ] || {{ echo "ERROR: chunk FASTQ not found: $CHUNK_FASTQ" >&2; exit 1; }}

# Skip if alignment BAM already exists (safe resume after preemption or requeue)
if compgen -G "$FINAL_OUTDIR/*.bam" > /dev/null 2>&1; then
    echo "Alignment BAM already exists at $FINAL_OUTDIR — skipping: $(date)"
    ls -lh "$FINAL_OUTDIR"/*.bam
    exit 0
fi

# Stage output to fast local scratch to avoid shared-filesystem I/O contention.
# Falls back to /tmp if $SCRATCH is not set (non-HPC or single-node use).
SCRATCH_DIR="${{SCRATCH:-/tmp}}/rectify_${{ALIGNER}}_${{SLURM_JOB_ID:-$$}}_${{RECTIFY_TASK_ID}}"
mkdir -p "$SCRATCH_DIR"
trap "rm -rf '$SCRATCH_DIR'" EXIT

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
            -o "$SCRATCH_DIR" \\
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
            -o "$SCRATCH_DIR" \\
            --verbose && break
    fi
    echo "Attempt $_try/$MAX_RETRIES failed for $ALIGNER chunk $CHUNK_PAD — retrying in 10s" >&2
    sleep 10
    if [ $_try -eq $MAX_RETRIES ]; then
        echo "ERROR: $ALIGNER chunk $CHUNK_PAD failed after $MAX_RETRIES attempts" >&2
        exit 1
    fi
done

echo "Rsync scratch → output dir: $(date)"
rsync -a "$SCRATCH_DIR/" "$FINAL_OUTDIR/"
echo "Done: $(date)"
ls -lh "$FINAL_OUTDIR"/*.bam 2>/dev/null || true
"""


def _merge_aligners_body(
    n_chunks: int, all_aligners: List[str], sample_prefix: str,
    output_dir: Path, python_path: str, rectify_src: str,
) -> str:
    """Merge per-chunk per-aligner BAMs into one sorted BAM per aligner (for prescan)."""
    limits = _thread_limits_block('$MERGE_CPUS')
    aligner_arr = ' '.join(f'"{a}"' for a in all_aligners)
    return f"""#!/bin/bash
# RECTIFY — merge per-chunk per-aligner BAMs into full per-aligner BAMs
# Step: after alignment arrays, before prescan.

{_SCHEDULER_HEADER_PLACEHOLDER_MERGE_ALIGNERS}

set -euo pipefail
PYTHON="{python_path}"
export PATH="$(dirname "$PYTHON"):$HOME/bin:$HOME/.rectify/bin:$PATH"

MERGE_CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-{MERGE_CORES}}}}}}}
{limits}

OUTDIR="{output_dir}"
SAMPLE="{sample_prefix}"
N_CHUNKS={n_chunks}
ALIGNERS=({aligner_arr})

MERGED_DIR="$OUTDIR/merged_bams"
mkdir -p "$MERGED_DIR"

echo "Host:  $(hostname)"
echo "CPUs:  $MERGE_CPUS"
echo "Start: $(date)"

for ALIGNER in "${{ALIGNERS[@]}}"; do
    echo ""
    echo "=== Merging $ALIGNER ==="
    CHUNK_BAMS=()
    for (( k=0; k<N_CHUNKS; k++ )); do
        CHUNK_PAD=$(printf "%03d" $k)
        CHUNK_OUTDIR="$OUTDIR/aligner_chunks/$ALIGNER/chunk_$CHUNK_PAD"
        BAM=$(ls "$CHUNK_OUTDIR"/*.${{ALIGNER}}.bam 2>/dev/null | head -1 || true)
        [ -n "$BAM" ] && CHUNK_BAMS+=("$BAM") || echo "  WARNING: no BAM in $CHUNK_OUTDIR" >&2
    done
    if [ ${{#CHUNK_BAMS[@]}} -eq 0 ]; then
        echo "  SKIPPING $ALIGNER — no chunk BAMs found" >&2; continue
    fi
    MERGED_BAM="$MERGED_DIR/${{SAMPLE}}.${{ALIGNER}}.bam"
    samtools merge -f -@ "$MERGE_CPUS" "$MERGED_BAM" "${{CHUNK_BAMS[@]}}"
    samtools sort  -@ "$MERGE_CPUS" -o "${{MERGED_BAM%.bam}}.sorted.bam" "$MERGED_BAM"
    mv "${{MERGED_BAM%.bam}}.sorted.bam" "$MERGED_BAM"
    samtools index -@ "$MERGE_CPUS" "$MERGED_BAM"
    echo "  ✓ $(du -h "$MERGED_BAM" | cut -f1)  $MERGED_BAM"
done

echo ""
echo "Done: $(date)"
ls -lh "$MERGED_DIR/"
"""


def _prescan_body(
    all_aligners: List[str], sample_prefix: str,
    output_dir: Path, genome: str, annot: str,
    python_path: str, rectify_src: str,
) -> str:
    """Build junction pool pickle and variant scan from merged per-aligner BAMs."""
    limits = _thread_limits_block('$PRESCAN_CPUS')
    aligner_bam_args = '\n'.join(
        f'    --aligner-bams "{a}:$MERGED_DIR/{sample_prefix}.{a}.bam" \\'
        for a in all_aligners
    )
    # Use mapPacBio as the primary BAM for variant scan; fall back to first available aligner.
    primary_aligner = 'mapPacBio' if 'mapPacBio' in all_aligners else all_aligners[0]
    primary_bam = f'$MERGED_DIR/{sample_prefix}.{primary_aligner}.bam'
    return f"""#!/bin/bash
# RECTIFY — prescan: build variant scan + junction pool from merged per-aligner BAMs
# Outputs: rescue_scan.pkl (Module 2D), junction_pool.pkl (Module 2H)
# Both are passed to per-chunk correction jobs via --junction-pool-cache.

{_SCHEDULER_HEADER_PLACEHOLDER_PRESCAN}

set -euo pipefail
PYTHON="{python_path}"
export PATH="$(dirname "$PYTHON"):$HOME/bin:$HOME/.rectify/bin:$PATH"

PRESCAN_CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-{PRESCAN_CORES}}}}}}}
{limits}

RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"
MERGED_DIR="$OUTDIR/merged_bams"

echo "Host:  $(hostname)"
echo "CPUs:  $PRESCAN_CPUS"
echo "Start: $(date)"

cd "$RECTIFY_SRC"
$PYTHON -m rectify prescan \\
    --bam "{primary_bam}" \\
{aligner_bam_args}
    --genome "{genome}" \\
    --annotation "{annot}" \\
    --threads "$PRESCAN_CPUS" \\
    -o "$OUTDIR"

echo "Done: $(date)"
ls -lh "$OUTDIR/rescue_scan.pkl" "$OUTDIR/junction_pool.pkl" 2>/dev/null
"""


def _correct_array_body(
    aligner: str, n_chunks: int, all_aligners: List[str], sample_prefix: str,
    output_dir: Path, genome: str, annot: str,
    python_path: str, rectify_src: str,
    junction_penalty_table: str = '',
    str_penalty_table: str = '',
    write_per_aligner_corrected_bams: bool = False,
) -> str:
    """Per-aligner correction array: one task per chunk. Uses --junction-pool-cache."""
    limits = _thread_limits_block('$CORRECT_CPUS')
    aligner_bam_args = '\n'.join(
        f'    --aligner-bams "{a}:$MERGED_DIR/{sample_prefix}.{a}.bam" \\'
        for a in all_aligners
    )
    jpt_line = f'    --junction-penalty-table "{junction_penalty_table}" \\' if junction_penalty_table else ''
    spt_line = f'    --str-penalty-table "{str_penalty_table}" \\' if str_penalty_table else ''
    write_cb_line = (
        '    --write-corrected-bam "$SCRATCH_WORK/corrected.bam" \\\n'
        if write_per_aligner_corrected_bams else ''
    )
    mv_corrected_bam_block = (
        """\
# Move per-aligner corrected BAM out of scratch to persistent storage for IGV.
samtools sort -@ "$CORRECT_CPUS" \\
    -o "$CHUNK_OUTDIR/corrected.sorted.tmp.bam" \\
    "$SCRATCH_WORK/corrected.bam"
mv "$CHUNK_OUTDIR/corrected.sorted.tmp.bam" "$CHUNK_OUTDIR/corrected.bam"
samtools index "$CHUNK_OUTDIR/corrected.bam"
"""
        if write_per_aligner_corrected_bams else ''
    )
    return f"""#!/bin/bash
# RECTIFY — per-chunk correction for aligner: {aligner}
# Generated by: rectify split --generate-slurm
# {n_chunks} chunks × 1 aligner = {n_chunks} total tasks
# Uses --junction-pool-cache to skip pool rebuild (pre-built by run_prescan.sh)

{_SCHEDULER_HEADER_PLACEHOLDER_CORRECT}

set -euo pipefail
PYTHON="{python_path}"
export PATH="$(dirname "$PYTHON"):$HOME/bin:$HOME/.rectify/bin:$PATH"

if   [ -n "${{SLURM_ARRAY_TASK_ID:-}}" ]; then TASK_ID=$SLURM_ARRAY_TASK_ID; CORRECT_CPUS=${{SLURM_CPUS_PER_TASK:-{CORRECT_CORES}}}
elif [ -n "${{SGE_TASK_ID:-}}" ];          then TASK_ID=$(( SGE_TASK_ID - 1 )); CORRECT_CPUS=${{NSLOTS:-{CORRECT_CORES}}}
elif [ -n "${{PBS_ARRAY_INDEX:-}}" ];      then TASK_ID=$(( PBS_ARRAY_INDEX - 1 )); CORRECT_CPUS=${{PBS_NUM_PPN:-{CORRECT_CORES}}}
else echo "ERROR: no scheduler array task variable" >&2; exit 1; fi

{limits}

RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"
MERGED_DIR="$OUTDIR/merged_bams"
ALIGNER="{aligner}"
CHUNK_PAD=$(printf "%03d" $TASK_ID)
CHUNK_OUTDIR="$OUTDIR/aligner_chunks/$ALIGNER/chunk_$CHUNK_PAD"
IN_BAM=$(ls "$CHUNK_OUTDIR"/*.${{ALIGNER}}.bam 2>/dev/null | head -1)
POOL_PKL="$OUTDIR/junction_pool.pkl"
SCAN_PKL="$OUTDIR/rescue_scan.pkl"

echo "Host:    $(hostname)"
echo "Aligner: $ALIGNER  Chunk: $CHUNK_PAD  CPUs: $CORRECT_CPUS"
echo "Input:   $IN_BAM"
echo "Pool:    $POOL_PKL"
echo "Scan:    $SCAN_PKL"
echo "Start:   $(date)"

[ -f "$IN_BAM" ]   || {{ echo "ERROR: BAM not found: $IN_BAM" >&2; exit 1; }}
[ -f "$POOL_PKL" ] || {{ echo "ERROR: pool pickle not found: $POOL_PKL" >&2; exit 1; }}
[ -f "$SCAN_PKL" ] || {{ echo "WARNING: rescue_scan.pkl not found — variant-aware rescue will re-scan" >&2; }}

# Skip if already complete. Per-aligner corrected BAMs are now debug artifacts;
# corrected TSVs are enough for lazy HP scoring and final consensus BAM writing.
if [ -f "$CHUNK_OUTDIR/corrected_reads.tsv" ]; then
    echo "Chunk $CHUNK_PAD already corrected — skipping: $(date)"
    ls -lh "$CHUNK_OUTDIR/corrected_reads.tsv"
    exit 0
fi

# Checkpoint state must live on the persistent output directory (not $SCRATCH,
# which is ephemeral and node-local — purged on job relaunch).
CHECKPOINT_DIR="$CHUNK_OUTDIR/.checkpoints"
mkdir -p "$CHECKPOINT_DIR"

# Stage INPUT BAM to fast local scratch for random I/O; output BAM also on scratch.
# Falls back to /tmp if $SCRATCH is not set (non-HPC or single-node use).
SCRATCH_WORK="${{SCRATCH:-/tmp}}/rectify_${{SLURM_JOB_ID:-$$}}_${{TASK_ID}}"
mkdir -p "$SCRATCH_WORK"
trap "rm -rf '$SCRATCH_WORK'" EXIT
cp "$IN_BAM" "$SCRATCH_WORK/"
cp "$IN_BAM".bai "$SCRATCH_WORK/" 2>/dev/null || true
LOCAL_BAM="$SCRATCH_WORK/$(basename $IN_BAM)"

SCAN_CACHE_ARG=""
[ -f "$SCAN_PKL" ] && SCAN_CACHE_ARG="--variant-scan-cache $SCAN_PKL"

cd "$RECTIFY_SRC"
$PYTHON -m rectify correct \\
    "$LOCAL_BAM" \\
    --genome "{genome}" \\
    --annotation "{annot}" \\
{jpt_line}
{spt_line}
    $SCAN_CACHE_ARG \\
    --junction-pool-cache "$POOL_PKL" \\
{aligner_bam_args}
{write_cb_line}    --threads "$CORRECT_CPUS" \\
    --streaming \\
    --checkpoint-dir "$CHECKPOINT_DIR" \\
    -o "$CHUNK_OUTDIR/corrected_reads.tsv"

{mv_corrected_bam_block}echo "Done: $(date)"
ls -lh "$CHUNK_OUTDIR/corrected_reads.tsv" 2>/dev/null
"""


def _chunk_merge_body(
    n_chunks: int, all_aligners: List[str], sample_prefix: str,
    output_dir: Path, genome: str, annot: str,
    python_path: str, rectify_src: str,
    junction_overhang_table: str = '',
    min_junction_anchor_bp: int = 0,
) -> str:
    """Per-chunk merge: merge corrected TSVs → consensus BAM → polya BAM."""
    limits = _thread_limits_block('$CHUNK_MERGE_CPUS')
    aligner_arr = ' '.join(f'"{a}"' for a in all_aligners)
    return f"""#!/bin/bash
# RECTIFY — per-chunk merge: corrected TSVs → consensus BAM (→ polya BAM)
# Generated by: rectify split --generate-slurm
# {n_chunks} tasks, one per chunk.

{_SCHEDULER_HEADER_PLACEHOLDER_CHUNK_MERGE}

set -euo pipefail
PYTHON="{python_path}"
export PATH="$(dirname "$PYTHON"):$HOME/bin:$HOME/.rectify/bin:$PATH"

if   [ -n "${{SLURM_ARRAY_TASK_ID:-}}" ]; then TASK_ID=$SLURM_ARRAY_TASK_ID; CHUNK_MERGE_CPUS=${{SLURM_CPUS_PER_TASK:-{CHUNK_MERGE_CORES}}}
elif [ -n "${{SGE_TASK_ID:-}}" ];          then TASK_ID=$(( SGE_TASK_ID - 1 )); CHUNK_MERGE_CPUS=${{NSLOTS:-{CHUNK_MERGE_CORES}}}
elif [ -n "${{PBS_ARRAY_INDEX:-}}" ];      then TASK_ID=$(( PBS_ARRAY_INDEX - 1 )); CHUNK_MERGE_CPUS=${{PBS_NUM_PPN:-{CHUNK_MERGE_CORES}}}
else echo "ERROR: no scheduler array task variable" >&2; exit 1; fi

{limits}

RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"
ALIGNERS=({aligner_arr})
CHUNK_PAD=$(printf "%03d" $TASK_ID)
CHUNK_OUT="$OUTDIR/chunks/chunk_$CHUNK_PAD"
mkdir -p "$CHUNK_OUT"

# PARQUET_PATH: set to per-sample poly-A trim metadata parquet for polya restore.
# Leave empty to skip poly-A restore.
PARQUET_PATH="${{RECTIFY_PARQUET_PATH:-}}"

echo "Host:  $(hostname)"
echo "Chunk: $CHUNK_PAD  CPUs: $CHUNK_MERGE_CPUS"
echo "Start: $(date)"

$PYTHON - << PYEOF
import sys, logging, os
from pathlib import Path
from collections import defaultdict
import pysam

logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
log = logging.getLogger(__name__)

sys.path.insert(0, '{rectify_src}')
from rectify.core.consensus.corrected_consensus import (
    merge_corrected_tsvs,
    write_corrected_consensus_bam,
    _stage_raw_bams,
)
from rectify.core.splice.calibrate_junction_overhang import OverhangTable
from rectify.utils.genome import load_genome
import pandas as pd

OUTDIR   = Path('{output_dir}')
ALIGNERS = {all_aligners!r}
CHUNK    = '$CHUNK_PAD'
CHUNK_OUT = Path('$CHUNK_OUT')
THREADS   = int('$CHUNK_MERGE_CPUS')
PARQUET   = '$PARQUET_PATH'
OVERHANG_TABLE_PATH = '{junction_overhang_table}'
GENOME_PATH = '{genome}'

def _bare(n):
    for sep in ('_pt:i:', ' pt:i:'):
        if sep in n: return n.split(sep)[0]
    return n

# ── Load junction overhang table (optional) ───────────────────────────────
_overhang_table = None
if OVERHANG_TABLE_PATH and Path(OVERHANG_TABLE_PATH).exists():
    try:
        _overhang_table = OverhangTable.from_tsv(OVERHANG_TABLE_PATH)
        log.info('Junction overhang table loaded: %s', OVERHANG_TABLE_PATH)
    except Exception as _e:
        log.warning('Could not load overhang table %s: %s — skipping filter', OVERHANG_TABLE_PATH, _e)

# ── Step 1: merge corrected TSVs ──────────────────────────────────────────
valid_tsvs = {{}}
corrected_bams = {{}}
raw_bams = {{}}
for a in ALIGNERS:
    chunk_dir = OUTDIR / 'aligner_chunks' / a / f'chunk_{{CHUNK}}'
    tsv = chunk_dir / 'corrected_reads.tsv'
    bam = chunk_dir / 'corrected.bam'
    rbam = next(chunk_dir.glob(f'*.{{a}}.bam'), None)
    if tsv.exists(): valid_tsvs[a] = tsv
    if bam.exists(): corrected_bams[a] = bam
    if rbam: raw_bams[a] = rbam

if not valid_tsvs:
    log.error('No corrected TSVs found for chunk %s', CHUNK); sys.exit(1)
log.info('Merging %d TSVs: %s', len(valid_tsvs), list(valid_tsvs.keys()))
genome = load_genome(GENOME_PATH)

# Stage raw BAMs to local NVMe ($L_SCRATCH / $SCRATCH) before HP scoring and BAM
# write to avoid OAK Lustre per-process client cache cold-read overhead.  Workers
# spawned by ProcessPoolExecutor each maintain their own Lustre cache, so without
# staging both HP scoring and write_corrected_consensus_bam hit cold Lustre.
import contextlib as _contextlib
_raw_bams_str = {{a: str(p) for a, p in raw_bams.items()}} if raw_bams else {{}}
_bam_stage_ctx = _contextlib.ExitStack()
_staged = _bam_stage_ctx.enter_context(_stage_raw_bams(_raw_bams_str))

# Prefer lazy HP-edit-distance scoring from raw BAMs + corrected TSVs. Existing
# per-aligner corrected BAMs are still accepted for backward-compatible reruns.
merge_corrected_tsvs(
    valid_tsvs,
    CHUNK_OUT/'corrected_reads.tsv',
    CHUNK_OUT/'aligner_summary.tsv',
    per_aligner_corrected_bams={{a: str(p) for a, p in corrected_bams.items()}} if corrected_bams else None,
    per_aligner_raw_bams=_staged if (raw_bams and not corrected_bams) else None,
    genome=genome,
    overhang_table=_overhang_table,
    lazy_scoring_workers=THREADS,
    min_junction_anchor_bp={min_junction_anchor_bp},
)
merged = pd.read_csv(CHUNK_OUT/'corrected_reads.tsv', sep='\\t', low_memory=False)
log.info('Merged: %d reads', merged['read_id'].nunique())

if 'winning_aligner' in merged.columns:
    for al, cnt in merged.groupby('winning_aligner')['read_id'].nunique().sort_values(ascending=False).items():
        n = merged['read_id'].nunique()
        log.info('  %s: %d (%.1f%%)', al, cnt, cnt/n*100)

# ── Step 2: build corrected consensus BAM ─────────────────────────────────
if raw_bams:
    stats = write_corrected_consensus_bam(
        _staged,
        valid_tsvs,
        CHUNK_OUT/'corrected_reads.tsv',
        CHUNK_OUT/'corrected_consensus.bam',
        genome=genome,
        threads=THREADS,
    )
    log.info('corrected_consensus.bam: %d reads', stats.get('written', 0))
elif corrected_bams:
    read_to_aligner = dict(zip(merged['read_id'], merged.get('winning_aligner', merged.get('aligner', ''))))
    aligner_to_reads = defaultdict(set)
    for rid, al in read_to_aligner.items():
        aligner_to_reads[al].add(rid)

    header_bam = next((b for b in corrected_bams.values() if b.exists()), None)
    if header_bam:
        with pysam.AlignmentFile(str(header_bam), 'rb') as hf:
            header = hf.header.to_dict()
        tmp = str(CHUNK_OUT / 'corrected_consensus.unsorted.bam')
        n_written = 0
        with pysam.AlignmentFile(tmp, 'wb', header=header) as out:
            for al, reads in aligner_to_reads.items():
                b = corrected_bams.get(al)
                if not b or not b.exists(): continue
                with pysam.AlignmentFile(str(b), 'rb') as bam_in:
                    for read in bam_in:
                        if read.is_secondary or read.is_supplementary: continue
                        if _bare(read.query_name) in reads:
                            out.write(read); n_written += 1
        pysam.sort('-o', str(CHUNK_OUT/'corrected_consensus.bam'), '-@', str(THREADS), tmp)
        pysam.index(str(CHUNK_OUT/'corrected_consensus.bam'))
        Path(tmp).unlink(missing_ok=True)
        log.info('corrected_consensus.bam: %d reads', n_written)

# ── Step 3: restore poly-A softclips (if parquet provided) ───────────────
if PARQUET and Path(PARQUET).exists() and raw_bams:
    try:
        from rectify.core.commands.restore_polya_command import restore_polya_softclips
        raw_paths = _staged
        tmp_p = str(CHUNK_OUT / 'corrected_polya.unsorted.bam')
        stats = restore_polya_softclips(
            str(CHUNK_OUT/'corrected_reads.tsv'), raw_paths, PARQUET, tmp_p, threads=THREADS
        )
        pysam.sort('-o', str(CHUNK_OUT/'corrected_polya.bam'), '-@', str(THREADS), tmp_p)
        pysam.index(str(CHUNK_OUT/'corrected_polya.bam'))
        Path(tmp_p).unlink(missing_ok=True)
        log.info('corrected_polya.bam: %s', stats)
    except Exception as e:
        log.warning('poly-A restore failed (non-fatal): %s', e)
elif PARQUET and not Path(PARQUET).exists():
    log.warning('PARQUET_PATH set but file not found: %s — skipping poly-A restore', PARQUET)
else:
    log.info('RECTIFY_PARQUET_PATH not set — skipping poly-A restore')

_bam_stage_ctx.close()
log.info('Chunk %s complete.', CHUNK)
PYEOF

echo "Done: $(date)"
ls -lh "$CHUNK_OUT/"
"""


def _final_merge_body(
    n_chunks: int, all_aligners: List[str], sample_prefix: str,
    output_dir: Path, python_path: str, rectify_src: str,
) -> str:
    """Merge all per-chunk outputs into final per-sample files."""
    limits = _thread_limits_block('$FINAL_MERGE_CPUS')
    return f"""#!/bin/bash
# RECTIFY — final merge: concat per-chunk TSVs and merge per-chunk BAMs
# Run after run_array_chunk_merge.sh completes all tasks.

{_SCHEDULER_HEADER_PLACEHOLDER_FINAL_MERGE}

set -euo pipefail
PYTHON="{python_path}"
export PATH="$(dirname "$PYTHON"):$HOME/bin:$HOME/.rectify/bin:$PATH"

FINAL_MERGE_CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-{FINAL_MERGE_CORES}}}}}}}
{limits}

OUTDIR="{output_dir}"
N_CHUNKS={n_chunks}
CONSENSUS_DIR="$OUTDIR/consensus"
mkdir -p "$CONSENSUS_DIR"

echo "Host:  $(hostname)"
echo "CPUs:  $FINAL_MERGE_CPUS"
echo "Start: $(date)"

# ── Concat corrected_reads.tsv ────────────────────────────────────────────
echo "=== Merging corrected_reads.tsv ==="
HEADER_WRITTEN=0
for (( k=0; k<N_CHUNKS; k++ )); do
    CHUNK_PAD=$(printf "%03d" $k)
    TSV="$OUTDIR/chunks/chunk_$CHUNK_PAD/corrected_reads.tsv"
    if [ ! -f "$TSV" ]; then
        echo "  WARNING: missing $TSV" >&2; continue
    fi
    if [ $HEADER_WRITTEN -eq 0 ]; then
        cat "$TSV" >> "$CONSENSUS_DIR/corrected_reads.tsv"
        HEADER_WRITTEN=1
    else
        tail -n +2 "$TSV" >> "$CONSENSUS_DIR/corrected_reads.tsv"
    fi
done
wc -l "$CONSENSUS_DIR/corrected_reads.tsv"

# ── Merge corrected_consensus.bam ────────────────────────────────────────
echo "=== Merging corrected_consensus.bam ==="
CONSENSUS_BAMS=()
for (( k=0; k<N_CHUNKS; k++ )); do
    CHUNK_PAD=$(printf "%03d" $k)
    BAM="$OUTDIR/chunks/chunk_$CHUNK_PAD/corrected_consensus.bam"
    [ -f "$BAM" ] && CONSENSUS_BAMS+=("$BAM")
done
if [ ${{#CONSENSUS_BAMS[@]}} -gt 0 ]; then
    samtools merge -f -@ "$FINAL_MERGE_CPUS" "$CONSENSUS_DIR/corrected_consensus.unsorted.bam" "${{CONSENSUS_BAMS[@]}}"
    samtools sort -@ "$FINAL_MERGE_CPUS" -o "$CONSENSUS_DIR/corrected_consensus.bam" "$CONSENSUS_DIR/corrected_consensus.unsorted.bam"
    samtools index -@ "$FINAL_MERGE_CPUS" "$CONSENSUS_DIR/corrected_consensus.bam"
    rm -f "$CONSENSUS_DIR/corrected_consensus.unsorted.bam"
    echo "  ✓ $(du -h $CONSENSUS_DIR/corrected_consensus.bam | cut -f1)"
fi

# ── Merge corrected_polya.bam (if present) ───────────────────────────────
POLYA_BAMS=()
for (( k=0; k<N_CHUNKS; k++ )); do
    CHUNK_PAD=$(printf "%03d" $k)
    BAM="$OUTDIR/chunks/chunk_$CHUNK_PAD/corrected_polya.bam"
    [ -f "$BAM" ] && POLYA_BAMS+=("$BAM")
done
if [ ${{#POLYA_BAMS[@]}} -gt 0 ]; then
    echo "=== Merging corrected_polya.bam (${{#POLYA_BAMS[@]}} chunks) ==="
    samtools merge -f -@ "$FINAL_MERGE_CPUS" "$CONSENSUS_DIR/corrected_polya.unsorted.bam" "${{POLYA_BAMS[@]}}"
    samtools sort -@ "$FINAL_MERGE_CPUS" -o "$CONSENSUS_DIR/corrected_polya.bam" "$CONSENSUS_DIR/corrected_polya.unsorted.bam"
    samtools index -@ "$FINAL_MERGE_CPUS" "$CONSENSUS_DIR/corrected_polya.bam"
    rm -f "$CONSENSUS_DIR/corrected_polya.unsorted.bam"
    echo "  ✓ $(du -h $CONSENSUS_DIR/corrected_polya.bam | cut -f1)"
fi

echo ""
echo "=== Final outputs ==="
ls -lh "$CONSENSUS_DIR/"
echo "Done: $(date)"
"""


def _consensus_per_chunk_body(
    n_chunks: int, all_aligners: List[str], sample_prefix: str,
    output_dir: Path, genome: str, annot: str,
    python_path: str, rectify_src: str,
) -> str:
    """Per-chunk consensus selection array (runs after all corrected.bam files exist)."""
    limits = _thread_limits_block('$CPUS')
    aligners_bash = ' '.join(f'"{a}"' for a in all_aligners)
    return f"""#!/bin/bash
# RECTIFY — per-chunk consensus selection array
# Runs rectify consensus independently on each chunk's {len(all_aligners)} corrected BAMs.
# {n_chunks} parallel tasks instead of 1 sequential task over the full read set.
# Validated: ~420k reads/chunk = ~100 min wall; 22× speedup vs sequential.
#
# Run after all run_array_correct_{{aligner}}.sh jobs complete.
# Follow with run_merge_consensus_chunks.sh to merge {n_chunks} BAMs → final consensus.

{_SCHEDULER_HEADER_PLACEHOLDER_CONSENSUS_PER_CHUNK}

set -euo pipefail
PYTHON="{python_path}"
export PATH="$(dirname "$PYTHON"):$HOME/bin:$HOME/.rectify/bin:$PATH"

if   [ -n "${{SLURM_ARRAY_TASK_ID:-}}" ]; then TASK_ID=$SLURM_ARRAY_TASK_ID; CPUS=${{SLURM_CPUS_PER_TASK:-{CONSENSUS_PER_CHUNK_CORES}}}
elif [ -n "${{SGE_TASK_ID:-}}" ];          then TASK_ID=$(( SGE_TASK_ID - 1 )); CPUS=${{NSLOTS:-{CONSENSUS_PER_CHUNK_CORES}}}
elif [ -n "${{PBS_ARRAY_INDEX:-}}" ];      then TASK_ID=$(( PBS_ARRAY_INDEX - 1 )); CPUS=${{PBS_NUM_PPN:-{CONSENSUS_PER_CHUNK_CORES}}}
else echo "ERROR: no scheduler array task variable found" >&2; exit 1; fi

{limits}

RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"
SAMPLE="{sample_prefix}"
GENOME="{genome}"
ANNOT="{annot}"
ALIGNERS=({aligners_bash})

CHUNK_PAD=$(printf "%03d" $TASK_ID)
CHUNK_OUT="$OUTDIR/consensus/chunk_${{CHUNK_PAD}}"
mkdir -p "$CHUNK_OUT"

echo "Host:    $(hostname)"
echo "Chunk:   $CHUNK_PAD  CPUs: $CPUS"
echo "Start:   $(date)"

# Skip if already complete (safe to resubmit after preemption)
CONSENSUS_BAM="$CHUNK_OUT/${{SAMPLE}}_chunk_${{CHUNK_PAD}}.consensus.bam"
if [ -f "$CONSENSUS_BAM" ] && [ -f "${{CONSENSUS_BAM}}.bai" ]; then
    echo "Consensus BAM already exists — skipping: $CONSENSUS_BAM"
    exit 0
fi

# Build aligner:path arguments for this chunk
BAM_ARGS=()
for ALIGNER in "${{ALIGNERS[@]}}"; do
    BAM="$OUTDIR/aligner_chunks/$ALIGNER/chunk_${{CHUNK_PAD}}/corrected.bam"
    if [ ! -f "$BAM" ]; then
        echo "ERROR: corrected BAM not found for $ALIGNER chunk $CHUNK_PAD: $BAM" >&2
        exit 1
    fi
    BAM_ARGS+=("${{ALIGNER}}:${{BAM}}")
done

cd "$RECTIFY_SRC"
$PYTHON -m rectify consensus \\
    "${{BAM_ARGS[@]}}" \\
    --genome "$GENOME" \\
    --annotation "$ANNOT" \\
    --prefix "${{SAMPLE}}_chunk_${{CHUNK_PAD}}" \\
    --no-bedgraph \\
    -o "$CHUNK_OUT"

echo "Done: $(date)"
ls -lh "$CHUNK_OUT/"*.bam 2>/dev/null || true
"""


def _merge_consensus_chunks_body(
    n_chunks: int, all_aligners: List[str], sample_prefix: str,
    output_dir: Path, genome: str,
    python_path: str, rectify_src: str,
) -> str:
    """Merge per-chunk consensus BAMs → final sample consensus BAM."""
    limits = _thread_limits_block('$CPUS')
    return f"""#!/bin/bash
# RECTIFY — merge per-chunk consensus BAMs → final sample consensus BAM
# Aggregates per-chunk aligner stats, coordinate-sorts, indexes, adds MD tags.
# Run after run_array_consensus_per_chunk.sh completes all {n_chunks} tasks.

{_SCHEDULER_HEADER_PLACEHOLDER_MERGE_CONSENSUS}

set -euo pipefail
PYTHON="{python_path}"
export PATH="$(dirname "$PYTHON"):$HOME/bin:$HOME/.rectify/bin:$PATH"

CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-{MERGE_CONSENSUS_CORES}}}}}}}
{limits}

OUTDIR="{output_dir}"
SAMPLE="{sample_prefix}"
N_CHUNKS={n_chunks}
GENOME="{genome}"
CONSENSUS_DIR="$OUTDIR/consensus"
FINAL_BAM="$CONSENSUS_DIR/${{SAMPLE}}.consensus.bam"
FINAL_STATS="$CONSENSUS_DIR/${{SAMPLE}}.consensus_aligner_stats.tsv"

echo "Host:  $(hostname)"
echo "CPUs:  $CPUS"
echo "Start: $(date)"

# Verify all per-chunk consensus BAMs exist before proceeding
echo "=== Verifying per-chunk consensus BAMs ==="
MISSING=0
CHUNK_BAMS=()
for (( k=0; k<N_CHUNKS; k++ )); do
    CHUNK_PAD=$(printf "%03d" $k)
    BAM="$CONSENSUS_DIR/chunk_${{CHUNK_PAD}}/${{SAMPLE}}_chunk_${{CHUNK_PAD}}.consensus.bam"
    if [ ! -f "$BAM" ]; then
        echo "  MISSING: $BAM" >&2
        MISSING=$(( MISSING + 1 ))
    else
        CHUNK_BAMS+=("$BAM")
    fi
done
if [ $MISSING -gt 0 ]; then
    echo "ERROR: $MISSING per-chunk consensus BAMs missing — cannot merge" >&2
    exit 1
fi
echo "  All $N_CHUNKS chunk BAMs present."

# Stage merge to fast local scratch for I/O performance
SCRATCH_DIR="${{SCRATCH:-/tmp}}/rectify_merge_consensus_${{SLURM_JOB_ID:-$$}}"
mkdir -p "$SCRATCH_DIR"
trap "rm -rf '$SCRATCH_DIR'" EXIT

# Merge all per-chunk consensus BAMs
echo ""
echo "=== Merging $N_CHUNKS consensus BAMs ==="
echo "Start merge: $(date)"
samtools merge -@ $CPUS -f "$SCRATCH_DIR/${{SAMPLE}}.consensus.unsorted.bam" "${{CHUNK_BAMS[@]}}"
echo "End merge:   $(date)"

# Coordinate-sort and index
echo "Coordinate-sorting..."
samtools sort -@ $CPUS -o "$SCRATCH_DIR/${{SAMPLE}}.consensus.bam" "$SCRATCH_DIR/${{SAMPLE}}.consensus.unsorted.bam"
samtools index "$SCRATCH_DIR/${{SAMPLE}}.consensus.bam"
rm "$SCRATCH_DIR/${{SAMPLE}}.consensus.unsorted.bam"

# Add MD tags
echo "Adding MD tags..."
samtools calmd -b "$SCRATCH_DIR/${{SAMPLE}}.consensus.bam" "$GENOME" \\
    > "$SCRATCH_DIR/${{SAMPLE}}.consensus.md.bam" 2>/dev/null
if [ -s "$SCRATCH_DIR/${{SAMPLE}}.consensus.md.bam" ]; then
    mv "$SCRATCH_DIR/${{SAMPLE}}.consensus.md.bam" "$SCRATCH_DIR/${{SAMPLE}}.consensus.bam"
    samtools index "$SCRATCH_DIR/${{SAMPLE}}.consensus.bam"
    echo "  MD tags added."
else
    echo "  WARNING: samtools calmd produced empty file — proceeding without MD tags."
    rm -f "$SCRATCH_DIR/${{SAMPLE}}.consensus.md.bam"
fi

# Copy final BAM to persistent storage
cp "$SCRATCH_DIR/${{SAMPLE}}.consensus.bam" "$FINAL_BAM"
cp "$SCRATCH_DIR/${{SAMPLE}}.consensus.bam.bai" "${{FINAL_BAM}}.bai"

echo "Final BAM: $FINAL_BAM ($(ls -lh "$FINAL_BAM" | awk '{{print $5}}'))"

# Aggregate per-chunk aligner stats TSVs
echo ""
echo "=== Aggregating aligner stats ==="
$PYTHON - << 'PYEOF'
import sys
from pathlib import Path
import pandas as pd

OUTDIR    = Path("{output_dir}")
SAMPLE    = "{sample_prefix}"
N_CHUNKS  = {n_chunks}
CONS_DIR  = OUTDIR / "consensus"

dfs = []
for k in range(N_CHUNKS):
    chunk_pad = f"{{k:03d}}"
    tsv = CONS_DIR / f"chunk_{{chunk_pad}}" / f"{{SAMPLE}}_chunk_{{chunk_pad}}.consensus_aligner_stats.tsv"
    if tsv.exists():
        dfs.append(pd.read_csv(tsv, sep='\\t'))
    else:
        print(f"  WARNING: stats TSV missing for chunk {{chunk_pad}}: {{tsv}}", file=sys.stderr)

if not dfs:
    print("ERROR: no stats TSVs found", file=sys.stderr)
    sys.exit(1)

combined = dfs[0].copy()
combined['count'] = sum(df['count'] for df in dfs)
total = combined.loc[combined['metric'] == 'total_reads', 'count'].values[0]
combined['percent'] = (combined['count'] / total * 100).round(2)

out_tsv = CONS_DIR / f"{{SAMPLE}}.consensus_aligner_stats.tsv"
combined.to_csv(out_tsv, sep='\\t', index=False)
print(f"Wrote aggregated stats: {{out_tsv}}")
print(combined[['metric', 'count', 'percent']].to_string(index=False))
PYEOF

echo ""
echo "Done: $(date)"
echo "Output: $FINAL_BAM"
ls -lh "$FINAL_BAM" "${{FINAL_BAM}}.bai" "$FINAL_STATS" 2>/dev/null
"""


def _make_submit_script(
    scheduler: str,
    mpb_script: Optional[Path],
    others_script: Path,
    merge_aligners_script: Path,
    prescan_script: Path,
    correct_scripts: dict,   # {aligner: Path}
    chunk_merge_script: Path,
    final_merge_script: Path,
    consensus_per_chunk_script: Path,
    merge_consensus_script: Path,
    output_dir: Path,
) -> str:
    """Generate a submit_pipeline.sh that chains jobs with correct-first dependencies."""
    if scheduler == 'slurm':
        lines = ['#!/bin/bash', '# Submit the full chunked correct-first pipeline (SLURM)', '']
        if mpb_script:
            lines += [
                f'MPB_JOB=$(sbatch {mpb_script} | awk \'{{print $4}}\')',
                'echo "mapPacBio array: $MPB_JOB"',
            ]
        lines += [
            f'OTHERS_JOB=$(sbatch {others_script} | awk \'{{print $4}}\')',
            'echo "Others array:    $OTHERS_JOB"',
        ]
        align_dep = '--dependency=afterok:${MPB_JOB}:${OTHERS_JOB}' if mpb_script else '--dependency=afterok:${OTHERS_JOB}'
        lines += [
            f'MERGE_ALIGNERS_JOB=$(sbatch {align_dep} {merge_aligners_script} | awk \'{{print $4}}\')',
            'echo "Merge aligners:  $MERGE_ALIGNERS_JOB"',
            f'PRESCAN_JOB=$(sbatch --dependency=afterok:${{MERGE_ALIGNERS_JOB}} {prescan_script} | awk \'{{print $4}}\')',
            'echo "Prescan:         $PRESCAN_JOB"',
        ]
        correct_job_vars = []
        for aligner, script in correct_scripts.items():
            var = f'CORRECT_{aligner.upper()}_JOB'
            lines += [
                f'{var}=$(sbatch --dependency=afterok:${{PRESCAN_JOB}} {script} | awk \'{{print $4}}\')',
                f'echo "Correct {aligner}:   ${{{var}}}"',
            ]
            correct_job_vars.append(f'${{{var}}}')
        correct_dep = ':'.join(correct_job_vars)
        lines += [
            f'CHUNK_MERGE_JOB=$(sbatch --dependency=afterok:{correct_dep} {chunk_merge_script} | awk \'{{print $4}}\')',
            'echo "Chunk merge:     $CHUNK_MERGE_JOB"',
            f'FINAL_MERGE_JOB=$(sbatch --dependency=afterok:${{CHUNK_MERGE_JOB}} {final_merge_script} | awk \'{{print $4}}\')',
            'echo "Final merge:     $FINAL_MERGE_JOB"',
            f'CONSENSUS_PER_CHUNK_JOB=$(sbatch --dependency=afterok:${{FINAL_MERGE_JOB}} {consensus_per_chunk_script} | awk \'{{print $4}}\')',
            'echo "Consensus array: $CONSENSUS_PER_CHUNK_JOB"',
            f'MERGE_CONSENSUS_JOB=$(sbatch --dependency=afterok:${{CONSENSUS_PER_CHUNK_JOB}} {merge_consensus_script} | awk \'{{print $4}}\')',
            'echo "Merge consensus: $MERGE_CONSENSUS_JOB"',
            'echo ""',
            'echo "Full correct-first pipeline submitted. Monitor with: squeue -u $USER"',
        ]
    else:
        lines = [
            '#!/bin/bash',
            f'# Correct-first pipeline submission not yet implemented for {scheduler}.',
            '# Please submit scripts manually in this order:',
            f'# 1. run_array_mapPacBio.sh + run_array_others.sh',
            f'# 2. run_merge_aligners.sh (after alignment)',
            f'# 3. run_prescan.sh (after merge_aligners)',
            f'# 4. run_array_correct_*.sh (after prescan, one per aligner)',
            f'# 5. run_array_chunk_merge.sh (after all correct arrays)',
            f'# 6. run_final_merge.sh (after chunk_merge)',
            f'# 7. run_array_consensus_per_chunk.sh (after final_merge)',
            f'# 8. run_merge_consensus_chunks.sh (after consensus array)',
        ]
    return '\n'.join(lines) + '\n'


# Placeholder tokens replaced after header generation
_SCHEDULER_HEADER_PLACEHOLDER_MPB    = '__SCHED_HEADER_MPB__'
_SCHEDULER_HEADER_PLACEHOLDER_OTHERS = '__SCHED_HEADER_OTHERS__'
_SCHEDULER_HEADER_PLACEHOLDER_MERGE  = '__SCHED_HEADER_MERGE__'

_SCHEDULER_HEADER_PLACEHOLDER_MERGE_ALIGNERS = '__SCHED_HEADER_MERGE_ALIGNERS__'

# Short-read placeholders
_SCHEDULER_HEADER_PLACEHOLDER_SR_ARRAY        = '__SCHED_HEADER_SR_ARRAY__'
_SCHEDULER_HEADER_PLACEHOLDER_SR_FINAL_MERGE  = '__SCHED_HEADER_SR_FINAL_MERGE__'
_SCHEDULER_HEADER_PLACEHOLDER_PRESCAN              = '__SCHED_HEADER_PRESCAN__'
_SCHEDULER_HEADER_PLACEHOLDER_CORRECT              = '__SCHED_HEADER_CORRECT__'
_SCHEDULER_HEADER_PLACEHOLDER_CHUNK_MERGE          = '__SCHED_HEADER_CHUNK_MERGE__'
_SCHEDULER_HEADER_PLACEHOLDER_FINAL_MERGE          = '__SCHED_HEADER_FINAL_MERGE__'
_SCHEDULER_HEADER_PLACEHOLDER_CONSENSUS_PER_CHUNK  = '__SCHED_HEADER_CONSENSUS_PER_CHUNK__'
_SCHEDULER_HEADER_PLACEHOLDER_MERGE_CONSENSUS      = '__SCHED_HEADER_MERGE_CONSENSUS__'


def _generate_scripts(
    args: argparse.Namespace,
    n_chunks: int,
    sample_prefix: str,
    junction_penalty_table: str = '',
    str_penalty_table: str = '',
    junction_overhang_table: str = '',
    write_per_aligner_corrected_bams: bool = False,
) -> None:
    """Generate all scheduler scripts for the chunked correct-first pipeline."""
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
    all_aligners   = (['mapPacBio'] if not skip_mpb else []) + other_aligners

    # Junction-anchor consensus gate threshold baked into the chunk-merge script.
    from rectify.data import resolve_min_junction_anchor_bp
    min_junction_anchor_bp = resolve_min_junction_anchor_bp(args)
    # Co-activation guard: mapPacBio is only safe in the panel while the gate is
    # active. A gate-off run with mapPacBio re-opens the spurious-intron gaming
    # vector the gate exists to close (see memory project-hped-anchor-gate).
    if (not skip_mpb) and min_junction_anchor_bp == 0:
        print(
            "WARNING: mapPacBio is in the aligner panel but the junction-anchor "
            "gate is OFF (min_junction_anchor_bp=0). On human/spliced data this "
            "re-opens the spurious-intron gaming vector mapPacBio was gated for. "
            "Pass --organism human (gate=10), --min-junction-anchor-bp 10, or "
            "--skip-map-pacbio.",
            file=sys.stderr,
        )

    n_other_tasks = n_chunks * len(other_aligners)
    log_pat      = str(log_dir / '%A_%a')   # SLURM pattern; ignored by UGE/PBS headers
    log_pat_job  = str(log_dir / '%j')      # non-array SLURM pattern

    sched_kwargs = dict(
        slurm_partition=slurm_partition, slurm_account=slurm_account,
        uge_queue=uge_queue, uge_pe=uge_pe, pbs_queue=pbs_queue,
    )

    def _write(path: Path, body: str) -> Path:
        path.write_text(body)
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
        return path

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
        mpb_script_path = _write(output_dir / 'run_array_mapPacBio.sh', mpb_body)

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
    others_script_path = _write(output_dir / 'run_array_others.sh', others_body)

    # ── Merge aligners (single job) ────────────────────────────────────────
    merge_aligners_headers = _scheduler_headers(
        scheduler, f'{sample_prefix}_merge_aln', 1,
        MERGE_CORES, MERGE_MEM_GB, MERGE_TIME,
        str(log_dir), log_pat_job, is_array=False, **sched_kwargs,
    )
    merge_aligners_body = _merge_aligners_body(
        n_chunks, all_aligners, sample_prefix, output_dir, python_path, rectify_src,
    )
    merge_aligners_body = merge_aligners_body.replace(
        _SCHEDULER_HEADER_PLACEHOLDER_MERGE_ALIGNERS, merge_aligners_headers
    )
    merge_aligners_script_path = _write(output_dir / 'run_merge_aligners.sh', merge_aligners_body)

    # ── Prescan (single job) ───────────────────────────────────────────────
    prescan_headers = _scheduler_headers(
        scheduler, f'{sample_prefix}_prescan', 1,
        PRESCAN_CORES, PRESCAN_MEM_GB, PRESCAN_TIME,
        str(log_dir), log_pat_job, is_array=False, **sched_kwargs,
    )
    prescan_body_str = _prescan_body(
        all_aligners, sample_prefix, output_dir, genome_str, annot_str,
        python_path, rectify_src,
    )
    prescan_body_str = prescan_body_str.replace(_SCHEDULER_HEADER_PLACEHOLDER_PRESCAN, prescan_headers)
    prescan_script_path = _write(output_dir / 'run_prescan.sh', prescan_body_str)

    # ── Per-aligner correction arrays ─────────────────────────────────────
    correct_script_paths: dict = {}
    for aligner in all_aligners:
        correct_headers = _scheduler_headers(
            scheduler, f'{sample_prefix}_correct_{aligner}', n_chunks,
            CORRECT_CORES, CORRECT_MEM_GB, CORRECT_TIME,
            str(log_dir), log_pat, is_array=True, **sched_kwargs,
        )
        correct_body_str = _correct_array_body(
            aligner, n_chunks, all_aligners, sample_prefix, output_dir,
            genome_str, annot_str, python_path, rectify_src,
            junction_penalty_table=junction_penalty_table,
            str_penalty_table=str_penalty_table,
            write_per_aligner_corrected_bams=write_per_aligner_corrected_bams,
        )
        correct_body_str = correct_body_str.replace(
            _SCHEDULER_HEADER_PLACEHOLDER_CORRECT, correct_headers
        )
        script_path = _write(output_dir / f'run_array_correct_{aligner}.sh', correct_body_str)
        correct_script_paths[aligner] = script_path

    # ── Chunk merge array ─────────────────────────────────────────────────
    chunk_merge_headers = _scheduler_headers(
        scheduler, f'{sample_prefix}_chunk_merge', n_chunks,
        CHUNK_MERGE_CORES, CHUNK_MERGE_MEM_GB, CHUNK_MERGE_TIME,
        str(log_dir), log_pat, is_array=True, **sched_kwargs,
    )
    chunk_merge_body_str = _chunk_merge_body(
        n_chunks, all_aligners, sample_prefix, output_dir, genome_str, annot_str,
        python_path, rectify_src,
        junction_overhang_table=junction_overhang_table,
        min_junction_anchor_bp=min_junction_anchor_bp,
    )
    chunk_merge_body_str = chunk_merge_body_str.replace(
        _SCHEDULER_HEADER_PLACEHOLDER_CHUNK_MERGE, chunk_merge_headers
    )
    chunk_merge_script_path = _write(output_dir / 'run_array_chunk_merge.sh', chunk_merge_body_str)

    # ── Final merge (single job) ──────────────────────────────────────────
    final_merge_headers = _scheduler_headers(
        scheduler, f'{sample_prefix}_final_merge', 1,
        FINAL_MERGE_CORES, FINAL_MERGE_MEM_GB, FINAL_MERGE_TIME,
        str(log_dir), log_pat_job, is_array=False, **sched_kwargs,
    )
    final_merge_body_str = _final_merge_body(
        n_chunks, all_aligners, sample_prefix, output_dir, python_path, rectify_src,
    )
    final_merge_body_str = final_merge_body_str.replace(
        _SCHEDULER_HEADER_PLACEHOLDER_FINAL_MERGE, final_merge_headers
    )
    final_merge_script_path = _write(output_dir / 'run_final_merge.sh', final_merge_body_str)

    # ── Per-chunk consensus array ─────────────────────────────────────────
    consensus_per_chunk_headers = _scheduler_headers(
        scheduler, f'{sample_prefix}_consensus_chunks', n_chunks,
        CONSENSUS_PER_CHUNK_CORES, CONSENSUS_PER_CHUNK_MEM_GB, CONSENSUS_PER_CHUNK_TIME,
        str(log_dir), log_pat, is_array=True, **sched_kwargs,
    )
    consensus_per_chunk_body_str = _consensus_per_chunk_body(
        n_chunks, all_aligners, sample_prefix, output_dir, genome_str, annot_str,
        python_path, rectify_src,
    )
    consensus_per_chunk_body_str = consensus_per_chunk_body_str.replace(
        _SCHEDULER_HEADER_PLACEHOLDER_CONSENSUS_PER_CHUNK, consensus_per_chunk_headers
    )
    consensus_per_chunk_script_path = _write(
        output_dir / 'run_array_consensus_per_chunk.sh', consensus_per_chunk_body_str
    )

    # ── Merge consensus chunks (single job) ──────────────────────────────
    merge_consensus_headers = _scheduler_headers(
        scheduler, f'{sample_prefix}_merge_consensus', 1,
        MERGE_CONSENSUS_CORES, MERGE_CONSENSUS_MEM_GB, MERGE_CONSENSUS_TIME,
        str(log_dir), log_pat_job, is_array=False, **sched_kwargs,
    )
    merge_consensus_body_str = _merge_consensus_chunks_body(
        n_chunks, all_aligners, sample_prefix, output_dir, genome_str,
        python_path, rectify_src,
    )
    merge_consensus_body_str = merge_consensus_body_str.replace(
        _SCHEDULER_HEADER_PLACEHOLDER_MERGE_CONSENSUS, merge_consensus_headers
    )
    merge_consensus_script_path = _write(
        output_dir / 'run_merge_consensus_chunks.sh', merge_consensus_body_str
    )

    # ── Submission wrapper ─────────────────────────────────────────────────
    submit_body = _make_submit_script(
        scheduler,
        mpb_script_path,
        others_script_path,
        merge_aligners_script_path,
        prescan_script_path,
        correct_script_paths,
        chunk_merge_script_path,
        final_merge_script_path,
        consensus_per_chunk_script_path,
        merge_consensus_script_path,
        output_dir,
    )
    submit_path = _write(output_dir / 'submit_pipeline.sh', submit_body)

    # ── Summary ───────────────────────────────────────────────────────────
    logger.info("Generated correct-first scheduler scripts (%s):", scheduler.upper())
    if mpb_script_path:
        logger.info("  mapPacBio array    (%d tasks, %d cores, %dG):  %s",
                    n_chunks, MPB_CORES, MPB_MEM_GB, mpb_script_path)
    logger.info("  Others array       (%d tasks, %d cores, %dG):  %s",
                n_other_tasks, OTHER_CORES, OTHER_MEM_GB, others_script_path)
    logger.info("  Merge aligners     (1 job,   %d cores, %dG):  %s",
                MERGE_CORES, MERGE_MEM_GB, merge_aligners_script_path)
    logger.info("  Prescan            (1 job,   %d cores, %dG):  %s",
                PRESCAN_CORES, PRESCAN_MEM_GB, prescan_script_path)
    for aligner, sp in correct_script_paths.items():
        logger.info("  Correct %-10s (%d tasks, %d cores, %dG):  %s",
                    aligner, n_chunks, CORRECT_CORES, CORRECT_MEM_GB, sp)
    logger.info("  Chunk merge        (%d tasks, %d cores, %dG):  %s",
                n_chunks, CHUNK_MERGE_CORES, CHUNK_MERGE_MEM_GB, chunk_merge_script_path)
    logger.info("  Final merge        (1 job,   %d cores, %dG):  %s",
                FINAL_MERGE_CORES, FINAL_MERGE_MEM_GB, final_merge_script_path)
    logger.info("  Consensus array    (%d tasks, %d cores, %dG):  %s",
                n_chunks, CONSENSUS_PER_CHUNK_CORES, CONSENSUS_PER_CHUNK_MEM_GB,
                consensus_per_chunk_script_path)
    logger.info("  Merge consensus    (1 job,   %d cores, %dG):  %s",
                MERGE_CONSENSUS_CORES, MERGE_CONSENSUS_MEM_GB, merge_consensus_script_path)
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
    junction_penalty_table: str = '',
    str_penalty_table: str = '',
) -> dict:
    """
    Public API: generate correct-first pipeline scripts for a given chunk configuration.

    Returns dict with keys:
      'mpb'                    — Path or None (run_array_mapPacBio.sh)
      'others'                 — Path (run_array_others.sh)
      'merge_aligners'         — Path (run_merge_aligners.sh)
      'prescan'                — Path (run_prescan.sh)
      'correct'                — dict {aligner: Path} (run_array_correct_{aligner}.sh)
      'chunk_merge'            — Path (run_array_chunk_merge.sh)
      'final_merge'            — Path (run_final_merge.sh)
      'consensus_per_chunk'    — Path (run_array_consensus_per_chunk.sh)
      'merge_consensus_chunks' — Path (run_merge_consensus_chunks.sh)
      'submit'                 — Path (submit_pipeline.sh)

    All scripts are written to output_dir.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    _other_aligners = other_aligners or OTHER_ALIGNERS_DEFAULT
    all_aligners = (['mapPacBio'] if not skip_map_pacbio else []) + list(_other_aligners)

    args = argparse.Namespace(
        output_dir=output_dir,
        genome=genome,
        annotation=annotation,
        python_path=python_path,
        rectify_src=rectify_src,
        scheduler=scheduler,
        other_aligners=_other_aligners,
        skip_map_pacbio=skip_map_pacbio,
        slurm_partition=slurm_partition,
        slurm_account=slurm_account,
        uge_queue=uge_queue,
        uge_pe=uge_pe,
        pbs_queue=pbs_queue,
    )
    _generate_scripts(
        args, n_chunks, sample_prefix,
        junction_penalty_table=junction_penalty_table,
        str_penalty_table=str_penalty_table,
    )

    return {
        'mpb':            output_dir / 'run_array_mapPacBio.sh' if not skip_map_pacbio else None,
        'others':         output_dir / 'run_array_others.sh',
        'merge_aligners': output_dir / 'run_merge_aligners.sh',
        'prescan':        output_dir / 'run_prescan.sh',
        'correct':                {a: output_dir / f'run_array_correct_{a}.sh' for a in all_aligners},
        'chunk_merge':            output_dir / 'run_array_chunk_merge.sh',
        'final_merge':            output_dir / 'run_final_merge.sh',
        'consensus_per_chunk':    output_dir / 'run_array_consensus_per_chunk.sh',
        'merge_consensus_chunks': output_dir / 'run_merge_consensus_chunks.sh',
        'submit':                 output_dir / 'submit_pipeline.sh',
    }


# ── Short-read pipeline script generators ────────────────────────────────────

def _short_read_array_body(
    n_chunks: int,
    sample_prefix: str,
    output_dir: Path,
    genome: str,
    annot: str,
    python_path: str,
    rectify_src: str,
    dt_primed_cdna: bool = True,
) -> str:
    """
    Per-chunk short-read pipeline body: align (bbmap+bwa) → correct → consensus.

    I/O policy: all intermediate BAMs live in a per-task $SCRATCH directory.
    Only the per-chunk consensus BAM is copied back to {output_dir}/chunk_outputs/.
    Final Oak copy is handled by run_final_merge_short_read.sh.

    DO NOT change SCRATCH_WORK to point at Oak/NFS — BGZF-compressed BAM writes
    to NFS run at ~1 MB/min; $SCRATCH runs at ~75 GB/s aggregate bandwidth.
    """
    shim = _TASK_ID_SHIM.replace('__CORES__', str(SR_ALIGN_CORRECT_CORES))
    limits = _thread_limits_block()
    dt_flag = '--dT-primed-cDNA' if dt_primed_cdna else ''
    n_pad = f'{n_chunks:03d}'
    return f"""#!/bin/bash
# RECTIFY — short-read per-chunk pipeline (align + correct + consensus)
# Generated by: rectify split --generate-slurm --short-read
# Protocol: Illumina QuantSeq REV / dT-primed cDNA (bbmap + bwa aligners)
# {n_chunks} tasks, one per chunk.
#
# ── I/O POLICY ────────────────────────────────────────────────────────────────
# All intermediate BAMs write to $SCRATCH (per-task scratch dir).
# Only the per-chunk consensus BAM is copied to {output_dir}/chunk_outputs/.
# Final merged BAM is copied to Oak by run_final_merge_short_read.sh.
# DO NOT redirect intermediate BAM writes to Oak/NFS — it is ~75x slower.
# ──────────────────────────────────────────────────────────────────────────────

{_SCHEDULER_HEADER_PLACEHOLDER_SR_ARRAY}

set -euo pipefail
PYTHON="{python_path}"
export PATH="$(dirname "$PYTHON"):$HOME/bin:$HOME/.rectify/bin:$PATH"

{shim}

{limits}

RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"
GENOME="{genome}"
ANNOT="{annot}"
N_CHUNKS={n_chunks}
SAMTOOLS="${{SAMTOOLS:-samtools}}"

CHUNK_PAD=$(printf "%03d" "$RECTIFY_TASK_ID")
PREFIX="{sample_prefix}_chunk_${{CHUNK_PAD}}_of_{n_pad}"
CHUNK_FASTQ="$OUTDIR/${{PREFIX}}.fastq.gz"

# ── Per-task scratch working directory ──────────────────────────────────────
# $SCRATCH is set by the scheduler (Sherlock: /scratch/users/$USER).
# Falls back to /tmp if $SCRATCH is not available.
SCRATCH_WORK="${{SCRATCH:-/tmp}}/rectify_sr_${{SLURM_JOB_ID:-$$}}_${{RECTIFY_TASK_ID}}"
mkdir -p "$SCRATCH_WORK"

BBMAP_BAM="$SCRATCH_WORK/${{PREFIX}}.bbmap.bam"
BWA_BAM="$SCRATCH_WORK/${{PREFIX}}.bwa.bam"
BBMAP_CORR_BAM="$SCRATCH_WORK/${{PREFIX}}.bbmap_corrected.bam"
BWA_CORR_BAM="$SCRATCH_WORK/${{PREFIX}}.bwa_corrected.bam"
BBMAP_CORR_TSV="$SCRATCH_WORK/${{PREFIX}}.bbmap_corrected.tsv"
BWA_CORR_TSV="$SCRATCH_WORK/${{PREFIX}}.bwa_corrected.tsv"
CONSENSUS_BAM="$SCRATCH_WORK/${{PREFIX}}.consensus.bam"

# Per-chunk output directory (on $OUTDIR, which should be on $SCRATCH too)
CHUNK_OUT_DIR="$OUTDIR/chunk_outputs"
mkdir -p "$CHUNK_OUT_DIR"

echo "Python:  $($PYTHON --version 2>&1)"
echo "Host:    $(hostname)"
echo "Task:    $RECTIFY_TASK_ID  Chunk: $CHUNK_PAD  CPUs: $RECTIFY_CPUS"
echo "Scratch: $SCRATCH_WORK"
echo "Input:   $CHUNK_FASTQ"
echo "Start:   $(date)"
echo "======================================"

[ -f "$CHUNK_FASTQ" ] || {{ echo "ERROR: chunk FASTQ not found: $CHUNK_FASTQ" >&2; exit 1; }}

cd "$RECTIFY_SRC"

# ── Step 1: Align (bbmap + bwa) ──────────────────────────────────────────────
echo ""
echo "=== Step 1: rectify align --short-read --no-consensus ==="
time $PYTHON -m rectify align \\
    "$CHUNK_FASTQ" \\
    --short-read \\
    --no-consensus \\
    --genome "$GENOME" \\
    --annotation "$ANNOT" \\
    --sort --index \\
    --prefix "$PREFIX" \\
    -t "$RECTIFY_CPUS" \\
    -o "$SCRATCH_WORK/" \\
    --verbose
echo "  Align done: $(date)"

# Sort+index guard — defensive; rectify align --sort --index should handle this
for BAM in "$BBMAP_BAM" "$BWA_BAM"; do
    if [[ ! -f "${{BAM}}.bai" ]]; then
        echo "  Sorting and indexing $BAM..."
        $SAMTOOLS sort -@ "$RECTIFY_CPUS" -o "${{BAM}}.s.tmp.bam" "$BAM"
        mv "${{BAM}}.s.tmp.bam" "$BAM"
        $SAMTOOLS index "$BAM"
    fi
done

# ── Step 2a: rectify correct --short-read (bbmap) ───────────────────────────
echo ""
echo "=== Step 2a: rectify correct --short-read {dt_flag} (bbmap) ==="
time $PYTHON -m rectify correct "$BBMAP_BAM" \\
    --short-read {dt_flag} \\
    --streaming \\
    --genome "$GENOME" \\
    --annotation "$ANNOT" \\
    -t "$RECTIFY_CPUS" \\
    --write-corrected-bam "$BBMAP_CORR_BAM" \\
    -o "$BBMAP_CORR_TSV"
echo "  bbmap correct done: $(date)"
$SAMTOOLS sort -@ "$RECTIFY_CPUS" -o "${{BBMAP_CORR_BAM}}.s.tmp.bam" "$BBMAP_CORR_BAM"
mv "${{BBMAP_CORR_BAM}}.s.tmp.bam" "$BBMAP_CORR_BAM"
$SAMTOOLS index "$BBMAP_CORR_BAM"

# ── Step 2b: rectify correct --short-read (bwa) ─────────────────────────────
echo ""
echo "=== Step 2b: rectify correct --short-read {dt_flag} (bwa) ==="
time $PYTHON -m rectify correct "$BWA_BAM" \\
    --short-read {dt_flag} \\
    --streaming \\
    --genome "$GENOME" \\
    --annotation "$ANNOT" \\
    -t "$RECTIFY_CPUS" \\
    --write-corrected-bam "$BWA_CORR_BAM" \\
    -o "$BWA_CORR_TSV"
echo "  bwa correct done: $(date)"
$SAMTOOLS sort -@ "$RECTIFY_CPUS" -o "${{BWA_CORR_BAM}}.s.tmp.bam" "$BWA_CORR_BAM"
mv "${{BWA_CORR_BAM}}.s.tmp.bam" "$BWA_CORR_BAM"
$SAMTOOLS index "$BWA_CORR_BAM"

# ── Step 3: rectify consensus ────────────────────────────────────────────────
echo ""
echo "=== Step 3: rectify consensus ==="
time $PYTHON -m rectify consensus \\
    "bbmap:$BBMAP_CORR_BAM" \\
    "bwa:$BWA_CORR_BAM" \\
    --genome "$GENOME" \\
    --annotation "$ANNOT" \\
    --tiebreak compass \\
    --prefix "$PREFIX" \\
    -o "$SCRATCH_WORK/"
echo "  Consensus done: $(date)"

[ -f "$CONSENSUS_BAM" ] || {{ echo "ERROR: consensus BAM not found: $CONSENSUS_BAM" >&2; exit 1; }}
echo "  Consensus BAM: $(du -sh "$CONSENSUS_BAM" | cut -f1)"

# Copy consensus BAM to output directory (the only file that leaves scratch)
cp "$CONSENSUS_BAM" "$CHUNK_OUT_DIR/"

# Remove intermediate BAMs to keep scratch lean; keep correction TSVs for QC
rm -f "$BBMAP_BAM" "${{BBMAP_BAM}}.bai" \\
      "$BWA_BAM" "${{BWA_BAM}}.bai" \\
      "$BBMAP_CORR_BAM" "${{BBMAP_CORR_BAM}}.bai" \\
      "$BWA_CORR_BAM" "${{BWA_CORR_BAM}}.bai" \\
      "$CONSENSUS_BAM"
rmdir "$SCRATCH_WORK" 2>/dev/null || true
echo "  Cleaned up scratch"

echo ""
echo "======================================"
echo "Done: $(date)"
"""


def _short_read_final_merge_body(
    n_chunks: int,
    sample_prefix: str,
    output_dir: Path,
    oak_output_dir: str,
    python_path: str,
    rectify_src: str,
) -> str:
    """
    Final merge body for short-read pipeline: merge per-chunk consensus BAMs,
    sort, index, copy to persistent storage (Oak), clean up scratch.
    """
    limits = _thread_limits_block('$FINAL_MERGE_CPUS')
    n_pad = f'{n_chunks:03d}'
    return f"""#!/bin/bash
# RECTIFY — short-read final merge
# Generated by: rectify split --generate-slurm --short-read
# Merges {n_chunks} per-chunk consensus BAMs → sorted, indexed final BAM → copies to Oak.
# Cleans up per-sample scratch directory after verified copy.

{_SCHEDULER_HEADER_PLACEHOLDER_SR_FINAL_MERGE}

set -euo pipefail
PYTHON="{python_path}"
export PATH="$(dirname "$PYTHON"):$HOME/bin:$HOME/.rectify/bin:$PATH"

FINAL_MERGE_CPUS="${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-{SR_FINAL_MERGE_CORES}}}}}}}
{limits}

SAMTOOLS="${{SAMTOOLS:-samtools}}"
OUTDIR="{output_dir}"
OAK_OUT="{oak_output_dir}"
N_CHUNKS={n_chunks}
MERGED_DIR="$OUTDIR/merged"
mkdir -p "$MERGED_DIR" "$OAK_OUT"

echo "Host:    $(hostname)"
echo "CPUs:    $FINAL_MERGE_CPUS"
echo "Scratch: $OUTDIR"
echo "Oak out: $OAK_OUT"
echo "Start:   $(date)"
echo "======================================"

# ── Collect per-chunk consensus BAMs ─────────────────────────────────────────
echo ""
echo "=== Collecting $N_CHUNKS chunk consensus BAMs ==="
BAM_LIST=()
MISSING=0
for (( k=0; k<N_CHUNKS; k++ )); do
    CHUNK_PAD=$(printf "%03d" $k)
    BAM="$OUTDIR/chunk_outputs/{sample_prefix}_chunk_${{CHUNK_PAD}}_of_{n_pad}.consensus.bam"
    if [[ ! -f "$BAM" ]]; then
        echo "  ERROR: Missing $BAM" >&2
        MISSING=$(( MISSING + 1 ))
    else
        BAM_LIST+=("$BAM")
    fi
done
if [[ $MISSING -gt 0 ]]; then
    echo "ERROR: $MISSING chunk BAM(s) missing — aborting" >&2
    exit 1
fi
echo "  All $N_CHUNKS chunk BAMs found"

# ── Merge ─────────────────────────────────────────────────────────────────────
echo ""
echo "=== Merging ==="
UNSORTED="$MERGED_DIR/{sample_prefix}.consensus.unsorted.bam"
time $SAMTOOLS merge -@ "$FINAL_MERGE_CPUS" -f "$UNSORTED" "${{BAM_LIST[@]}}"
echo "  Merge done: $(date)"

# ── Sort ──────────────────────────────────────────────────────────────────────
echo ""
echo "=== Sorting ==="
FINAL_BAM="$MERGED_DIR/{sample_prefix}.consensus.bam"
time $SAMTOOLS sort -@ "$FINAL_MERGE_CPUS" -o "$FINAL_BAM" "$UNSORTED"
rm "$UNSORTED"
echo "  Sort done: $(date)  ($(du -sh "$FINAL_BAM" | cut -f1))"

# ── Index ─────────────────────────────────────────────────────────────────────
$SAMTOOLS index "$FINAL_BAM"
echo "  Index done: $(date)"

# ── Copy to persistent storage ────────────────────────────────────────────────
echo ""
echo "=== Copying to persistent storage ==="
cp "$FINAL_BAM"     "$OAK_OUT/{sample_prefix}.consensus.bam"
cp "$FINAL_BAM.bai" "$OAK_OUT/{sample_prefix}.consensus.bam.bai"
echo "  Copied: $(date)"
ls -lh "$OAK_OUT/{sample_prefix}.consensus.bam"

# ── Verify copy ───────────────────────────────────────────────────────────────
SRC_SIZE=$(wc -c < "$FINAL_BAM")
DST_SIZE=$(wc -c < "$OAK_OUT/{sample_prefix}.consensus.bam")
if [[ "$SRC_SIZE" -ne "$DST_SIZE" ]]; then
    echo "ERROR: size mismatch — src=$SRC_SIZE dst=$DST_SIZE" >&2
    exit 1
fi
echo "  Size verified: $SRC_SIZE bytes"

# ── Clean up scratch ──────────────────────────────────────────────────────────
echo ""
echo "=== Cleaning up scratch ==="
rm -rf "$OUTDIR"
echo "  Removed: $OUTDIR"

echo ""
echo "======================================"
echo "Done: $(date)"
"""


def _make_short_read_submit_script(
    scheduler: str,
    array_script: Path,
    final_merge_script: Path,
) -> str:
    """Generate submit_pipeline.sh for short-read two-stage pipeline."""
    if scheduler == 'slurm':
        lines = [
            '#!/bin/bash',
            '# Submit the short-read chunked correct-before-consensus pipeline (SLURM)',
            '# Stage 1: per-chunk align+correct+consensus array',
            '# Stage 2: merge all chunks, copy final BAM to Oak',
            '',
            f'ARRAY_JOB=$(sbatch --parsable {array_script})',
            'echo "Array job (align+correct+consensus): $ARRAY_JOB"',
            '',
            f'MERGE_JOB=$(sbatch --parsable --dependency=afterok:${{ARRAY_JOB}} {final_merge_script})',
            'echo "Merge job:                           $MERGE_JOB"',
            '',
            'echo ""',
            'echo "Monitor with: squeue -j $ARRAY_JOB,$MERGE_JOB"',
        ]
    else:
        lines = [
            '#!/bin/bash',
            f'# Short-read pipeline submission not yet implemented for {scheduler}.',
            '# Please submit scripts manually in this order:',
            f'# 1. {array_script}',
            f'# 2. {final_merge_script}  (after all array tasks complete)',
        ]
    return '\n'.join(lines) + '\n'


def _generate_short_read_scripts(
    args: argparse.Namespace,
    n_chunks: int,
    sample_prefix: str,
) -> None:
    """Generate SLURM scripts for the short-read chunked correct-before-consensus pipeline."""
    output_dir = args.output_dir.resolve()
    log_dir    = output_dir / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)

    genome_str = str(args.genome.resolve()) if args.genome else '# REQUIRED: set --genome'
    annot_str  = str(args.annotation.resolve()) if args.annotation else '# REQUIRED: set --annotation'
    python_path = getattr(args, 'python_path', 'python')
    rectify_src = getattr(args, 'rectify_src', '.')
    dt_primed_cdna = getattr(args, 'dt_primed_cdna', True)  # default True for --short-read

    oak_output_dir = getattr(args, 'oak_output_dir', None)
    oak_str = str(oak_output_dir.resolve()) if oak_output_dir else str(output_dir / 'final')

    scheduler = getattr(args, 'scheduler', 'slurm')
    slurm_partition = getattr(args, 'slurm_partition', None)
    slurm_account   = getattr(args, 'slurm_account', None)
    uge_queue       = getattr(args, 'uge_queue', 'long.q')
    uge_pe          = getattr(args, 'uge_pe', 'smp')
    pbs_queue       = getattr(args, 'pbs_queue', 'workq')

    log_pat     = str(log_dir / '%A_%a')
    log_pat_job = str(log_dir / '%j')

    sched_kwargs = dict(
        slurm_partition=slurm_partition, slurm_account=slurm_account,
        uge_queue=uge_queue, uge_pe=uge_pe, pbs_queue=pbs_queue,
    )

    def _write(path: Path, body: str) -> Path:
        path.write_text(body)
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
        return path

    # ── Per-chunk array script ────────────────────────────────────────────────
    array_headers = _scheduler_headers(
        scheduler, f'{sample_prefix}_sr', n_chunks,
        SR_ALIGN_CORRECT_CORES, SR_ALIGN_CORRECT_MEM_GB, SR_ALIGN_CORRECT_TIME,
        str(log_dir), log_pat, is_array=True, **sched_kwargs,
    )
    array_body = _short_read_array_body(
        n_chunks, sample_prefix, output_dir, genome_str, annot_str,
        python_path, rectify_src, dt_primed_cdna=dt_primed_cdna,
    )
    array_body = array_body.replace(_SCHEDULER_HEADER_PLACEHOLDER_SR_ARRAY, array_headers)
    array_script_path = _write(output_dir / 'run_array_short_read.sh', array_body)

    # ── Final merge script ────────────────────────────────────────────────────
    merge_headers = _scheduler_headers(
        scheduler, f'{sample_prefix}_sr_merge', 1,
        SR_FINAL_MERGE_CORES, SR_FINAL_MERGE_MEM_GB, SR_FINAL_MERGE_TIME,
        str(log_dir), log_pat_job, is_array=False, **sched_kwargs,
    )
    merge_body = _short_read_final_merge_body(
        n_chunks, sample_prefix, output_dir, oak_str, python_path, rectify_src,
    )
    merge_body = merge_body.replace(_SCHEDULER_HEADER_PLACEHOLDER_SR_FINAL_MERGE, merge_headers)
    merge_script_path = _write(output_dir / 'run_final_merge_short_read.sh', merge_body)

    # ── Submit wrapper ────────────────────────────────────────────────────────
    submit_body = _make_short_read_submit_script(scheduler, array_script_path, merge_script_path)
    _write(output_dir / 'submit_pipeline.sh', submit_body)

    logger.info("Generated short-read scheduler scripts (%s):", scheduler.upper())
    logger.info("  Array (align+correct+consensus):  %d tasks, %d cores, %dG, %s",
                n_chunks, SR_ALIGN_CORRECT_CORES, SR_ALIGN_CORRECT_MEM_GB, array_script_path)
    logger.info("  Final merge:                      1 job,   %d cores, %dG, %s",
                SR_FINAL_MERGE_CORES, SR_FINAL_MERGE_MEM_GB, merge_script_path)
    logger.info("  Oak output dir:                   %s", oak_str)
    logger.info("  Submit chain:  bash %s", output_dir / 'submit_pipeline.sh')


def _run_split_paired(args: argparse.Namespace) -> int:
    """Paired-end split: chunk R1/R2 by pair, shared RN, paired chunk FASTQs."""
    input_r1 = args.reads
    input_r2 = args.read2
    if not input_r2.exists():
        logger.error("Mate-2 (R2) file not found: %s", input_r2)
        return 1
    if input_r1.suffix == '.bam' or input_r2.suffix == '.bam':
        logger.error("Paired-end mode (--read2) requires FASTQ inputs, not BAM.")
        return 1

    logger.info("Counting read pairs in %s / %s ...", input_r1, input_r2)
    n_r1 = count_reads(input_r1)
    n_r2 = count_reads(input_r2)
    if n_r1 != n_r2:
        logger.error(
            "R1/R2 read-count mismatch: %s has %d, %s has %d — not a valid pair set",
            input_r1, n_r1, input_r2, n_r2,
        )
        return 1
    n_pairs = n_r1

    if args.n_chunks is not None:
        n_chunks = args.n_chunks
        if n_chunks < 1:
            logger.error("--n-chunks must be >= 1")
            return 1
    else:
        target = getattr(args, 'target_reads_per_chunk', DEFAULT_TARGET_READS)
        n_chunks = compute_n_chunks(n_pairs, target)
        logger.info("  Auto-sized: %d pairs / %d target = %d chunks (min %d, max %d)",
                    n_pairs, target, n_chunks, MIN_CHUNKS, MAX_CHUNKS)

    logger.info("  Total pairs:       %s", f"{n_pairs:,}")
    logger.info("  Chunks:            %d", n_chunks)
    logger.info("  ~pairs per chunk:  %s", f"{math.ceil(n_pairs / n_chunks):,}")

    if args.dry_run:
        return 0

    prefix = args.prefix or _derive_paired_prefix(input_r1)

    logger.info("Splitting paired %s + %s into %d chunks → %s",
                input_r1, input_r2, n_chunks, args.output_dir)
    r1_paths, r2_paths = split_fastq_paired(
        input_r1=input_r1,
        input_r2=input_r2,
        output_dir=args.output_dir,
        n_chunks=n_chunks,
        prefix=prefix,
    )

    manifest = {
        'paired': True,
        'input_r1': str(input_r1.resolve()),
        'input_r2': str(input_r2.resolve()),
        'n_pairs': n_pairs,
        'n_chunks': n_chunks,
        'pairs_per_chunk': math.ceil(n_pairs / n_chunks),
        'chunks_r1': [str(p) for p in r1_paths],
        'chunks_r2': [str(p) for p in r2_paths],
        'read_num_sidecar': str(
            (args.output_dir / f'{prefix}.read_num_sidecar.parquet').resolve()
        ),
    }
    manifest_path = args.output_dir / 'chunks_manifest.json'
    with open(str(manifest_path), 'w') as fh:
        json.dump(manifest, fh, indent=2)
    logger.info("Manifest: %s", manifest_path)

    if getattr(args, 'generate_slurm', False):
        # Paired short-read SLURM generation (full COMPASS panel) is wired in P4.
        logger.warning(
            "Paired chunk FASTQs + sidecar written. --generate-slurm for the "
            "paired short-read COMPASS pipeline is not yet wired (P4); chunks are "
            "ready for the array generator."
        )

    return 0


def run_split(args: argparse.Namespace) -> int:
    """Run split command."""
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    if not args.reads.exists():
        logger.error("Reads file not found: %s", args.reads)
        return 1

    # ── Paired-end (short-read) branch ────────────────────────────────────
    # When --read2 is given, the positional `reads` is R1. Pair-aware chunking
    # keeps both mates together (shared RN, same chunk) and emits paired chunk
    # FASTQs. This path is self-contained and returns early.
    if getattr(args, 'read2', None) is not None:
        return _run_split_paired(args)

    # ── BAM input: convert to FASTQ once ──────────────────────────────────
    # For DRS BAMs (--drs), use the same Step-0 trim run-all performs:
    # ``trim_drs_bam_polya`` walks back the poly(A)/adapter tail and emits
    # both a trimmed unaligned BAM and a parquet metadata sidecar (carries
    # original poly-A length + boundary, required by Step-4 restore-softclip).
    # Without --drs we fall back to ``samtools fastq -T '*'`` which preserves
    # all aux tags as FASTQ comments but does NOT trim — appropriate for
    # non-DRS BAMs where the poly-A tail isn't artificial.
    if args.reads.suffix == '.bam':
        args.output_dir.mkdir(parents=True, exist_ok=True)
        stem = args.reads.stem
        derived_fq = args.output_dir / f"{stem}_trimmed.fastq.gz"
        derived_meta = args.output_dir / f"{stem}_polya_trim_metadata.parquet"

        if getattr(args, 'drs', False):
            # DRS path: trim-polya → trimmed BAM → FASTQ (+ parquet metadata)
            if derived_fq.exists() and derived_fq.stat().st_size > 0 and derived_meta.exists():
                logger.info(
                    "DRS BAM: reusing existing trim outputs at %s (+ %s)",
                    derived_fq, derived_meta,
                )
            else:
                if shutil.which('samtools') is None:
                    logger.error("--drs BAM trim requires samtools on PATH.")
                    return 1
                from .drs_trim_command import trim_drs_bam_polya
                trimmed_bam = args.output_dir / f"{stem}_trimmed.bam"
                logger.info(
                    "DRS BAM: trim-polya %s → %s + %s",
                    args.reads, derived_fq, derived_meta,
                )
                trim_stats = trim_drs_bam_polya(
                    input_bam_path=str(args.reads),
                    output_bam_path=str(trimmed_bam),
                    metadata_path=str(derived_meta),
                    threads=max(1, getattr(args, 'threads', 4)),
                )
                logger.info(
                    "  Trimmed %s / %s reads",
                    f"{trim_stats.get('trimmed', 0):,}",
                    f"{trim_stats.get('total', 0):,}",
                )
                cmd = f"samtools fastq -@ {max(1, getattr(args, 'threads', 4) - 1)} -0 {derived_fq} {trimmed_bam}"
                try:
                    subprocess.run(cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    logger.error("samtools fastq failed (exit %d): %s", e.returncode, e)
                    return 1
                logger.info(
                    "  Metadata parquet preserved at %s — pass it to "
                    "`rectify restore-softclip --trim-metadata` after the "
                    "chunked correction pipeline finishes.",
                    derived_meta,
                )
        else:
            # Non-DRS BAM path: dump to FASTQ preserving aux tags as comments
            if shutil.which('samtools') is None:
                logger.error(
                    "BAM input requires samtools on PATH (got %s).",
                    args.reads,
                )
                return 1
            derived_fq = args.output_dir / f"{stem}.fastq.gz"
            if derived_fq.exists() and derived_fq.stat().st_size > 0:
                logger.info("BAM input: reusing existing FASTQ at %s", derived_fq)
            else:
                logger.info(
                    "BAM input: converting %s → %s (samtools fastq -T '*', preserves aux tags as comments)",
                    args.reads, derived_fq,
                )
                cmd = f"samtools fastq -T '*' {args.reads} | gzip > {derived_fq}"
                try:
                    subprocess.run(cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    logger.error("samtools fastq failed (exit %d): %s", e.returncode, e)
                    if derived_fq.exists():
                        derived_fq.unlink()
                    return 1
        args.reads = derived_fq

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
        'read_num_sidecar': str((args.output_dir / f'{prefix}.read_num_sidecar.parquet').resolve()),
    }
    manifest_path = args.output_dir / 'chunks_manifest.json'
    with open(str(manifest_path), 'w') as fh:
        json.dump(manifest, fh, indent=2)
    logger.info("Manifest: %s", manifest_path)

    # ── Generate scheduler scripts ────────────────────────────────────────
    if getattr(args, 'generate_slurm', False):
        if getattr(args, 'short_read', False):
            # Short-read mode: simplified two-stage pipeline (bbmap+bwa, no prescan)
            # --dT-primed-cDNA is implied by --short-read but can also be set explicitly
            if not getattr(args, 'dt_primed_cdna', False):
                args.dt_primed_cdna = True  # default for short-read cDNA protocols
            _generate_short_read_scripts(args, n_chunks, prefix)
        else:
            # Long-read mode: full correct-before-consensus pipeline
            jpt = getattr(args, 'junction_penalty_table', '') or ''
            spt = getattr(args, 'str_penalty_table', '') or ''
            jot = getattr(args, 'junction_overhang_table', '') or ''
            wpacb = getattr(args, 'write_per_aligner_corrected_bams', False)
            _generate_scripts(args, n_chunks, prefix,
                              junction_penalty_table=jpt, str_penalty_table=spt,
                              junction_overhang_table=jot,
                              write_per_aligner_corrected_bams=wpacb)

    return 0


def run(args: argparse.Namespace):
    """Entry point for CLI."""
    sys.exit(run_split(args))
