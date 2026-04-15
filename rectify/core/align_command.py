"""
Align Command for RECTIFY.

Runs multi-aligner alignment pipeline with consensus selection.

Supported aligners:
- minimap2: Fast seed-and-chain baseline
- mapPacBio: BBTools long-read aligner with splice-aware mode
- gapmm2: minimap2 wrapper with terminal exon refinement

Consensus selection:
- Runs all enabled aligners in parallel
- Compares alignments per read
- Selects best based on junction quality (canonical sites, annotation) and 5' rescue
- Outputs single rectified BAM with best alignment per read
- Note: 3' false junctions are handled by walk back correction, not consensus

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def create_align_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Create align subcommand parser."""
    parser = subparsers.add_parser(
        'align',
        help='Align reads using multiple aligners with junction annotation support',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    parser.add_argument(
        'reads',
        type=Path,
        help='Input FASTQ file (or FASTQ.GZ)'
    )

    parser.add_argument(
        '--genome',
        type=Path,
        required=True,
        help='Reference genome FASTA file'
    )

    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory for BAM files'
    )

    # Aligner selection
    aligner_group = parser.add_argument_group('Aligner selection')
    aligner_group.add_argument(
        '--aligners',
        nargs='+',
        choices=['minimap2', 'mapPacBio', 'gapmm2', 'all', 'none'],
        default=['all'],
        help='Default aligners for 3\' end correction. Use "none" to run only --junction-aligners (default: all)'
    )

    aligner_group.add_argument(
        '--junction-aligners',
        nargs='+',
        choices=['uLTRA', 'deSALT'],
        default=[],
        metavar='ALIGNER',
        help=(
            'Opt-in junction-mode aligners to add to the consensus pool '
            '(choices: uLTRA, deSALT). Requires --annotation. '
            'Value is unknown/untested — benchmark before using in production.'
        )
    )

    aligner_group.add_argument(
        '--no-consensus',
        action='store_true',
        help='Skip consensus selection, output separate BAMs per aligner'
    )

    aligner_group.add_argument(
        '--chimeric-consensus',
        action='store_true',
        default=False,
        help=(
            'Use chimeric consensus selection: independently pick the best aligner '
            'for each read segment, then assemble a chimeric CIGAR from the winners. '
            'Experimental — requires further validation before enabling by default.'
        )
    )

    aligner_group.add_argument(
        '--minimap2-path',
        default='minimap2',
        help='Path to minimap2 executable'
    )

    aligner_group.add_argument(
        '--mapPacBio-path',
        default='mapPacBio.sh',
        help='Path to mapPacBio.sh script'
    )

    aligner_group.add_argument(
        '--mapPacBio-chunks',
        type=int,
        default=1,
        metavar='N',
        help=(
            'Split mapPacBio alignment into N independent chunks for parallel '
            'SLURM array execution. Each chunk processes 1/N of the reads '
            '(interleaved, so read-length distribution is even). '
            'Use with --mapPacBio-chunk-idx K to run chunk K, or without '
            '--mapPacBio-chunk-idx to merge existing chunk BAMs.'
        )
    )

    aligner_group.add_argument(
        '--mapPacBio-chunk-idx',
        type=int,
        default=None,
        metavar='K',
        help=(
            '0-based index of the mapPacBio chunk to run (0 to N-1). '
            'Requires --mapPacBio-chunks N. Omit to trigger chunk-merge mode.'
        )
    )

    aligner_group.add_argument(
        '--gapmm2-path',
        default='gapmm2',
        help='Path to gapmm2 executable'
    )

    aligner_group.add_argument(
        '--ultra-path',
        default='uLTRA',
        help='Path to uLTRA executable'
    )

    aligner_group.add_argument(
        '--desalt-path',
        default='deSALT',
        help='Path to deSALT executable'
    )

    # Junction annotation
    junction_group = parser.add_argument_group('Junction annotation')
    junction_group.add_argument(
        '--annotation',
        type=Path,
        help='Gene annotation GFF/GTF for junction hints'
    )

    junction_group.add_argument(
        '--junc-bed',
        type=Path,
        help='Pre-computed junction BED file (overrides --annotation)'
    )

    junction_group.add_argument(
        '--junc-bonus',
        type=int,
        default=9,
        help='Bonus score for annotated junctions (minimap2)'
    )

    # Performance options
    perf_group = parser.add_argument_group('Performance options')
    perf_group.add_argument(
        '-t', '--threads',
        type=int,
        default=8,
        help='Number of threads per aligner'
    )

    perf_group.add_argument(
        '--parallel-aligners',
        action='store_true',
        help='Run aligners in parallel (uses more memory)'
    )

    # Output options
    output_group = parser.add_argument_group('Output options')
    output_group.add_argument(
        '--prefix',
        default='',
        help='Output file prefix'
    )

    output_group.add_argument(
        '--keep-sam',
        action='store_true',
        help='Keep intermediate SAM files'
    )

    output_group.add_argument(
        '--sort',
        action='store_true',
        default=True,
        help='Sort output BAM files by coordinate'
    )

    output_group.add_argument(
        '--index',
        action='store_true',
        default=True,
        help='Index output BAM files'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose logging'
    )

    return parser


def run_align(args: argparse.Namespace) -> int:
    """Run align command."""
    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Validate input
    if not args.reads.exists():
        logger.error(f"Reads file not found: {args.reads}")
        return 1

    if not args.genome.exists():
        logger.error(f"Genome file not found: {args.genome}")
        return 1

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Determine prefix
    prefix = args.prefix if args.prefix else args.reads.stem.replace('.fastq', '').replace('.fq', '')

    # Expand 'all' to list of aligners, then append any junction-mode aligners
    aligners = list(args.aligners)
    if 'none' in aligners:
        aligners = []
    elif 'all' in aligners:
        aligners = ['minimap2', 'mapPacBio', 'gapmm2']

    junction_aligners = getattr(args, 'junction_aligners', []) or []
    for ja in junction_aligners:
        if ja not in aligners:
            aligners.append(ja)

    # Import multi-aligner functions
    from .multi_aligner import (
        run_minimap2,
        run_map_pacbio,
        run_gapmm2,
        run_ultra,
        run_desalt,
        check_aligner_available,
    )

    # Generate junction BED if needed
    junc_bed_path = None
    if args.junc_bed:
        junc_bed_path = str(args.junc_bed)
    elif args.annotation:
        from ..utils.junction_bed import generate_junction_bed
        junc_bed_path = str(args.output_dir / f"{prefix}_junctions.bed")
        generate_junction_bed(str(args.annotation), junc_bed_path)
        logger.info(f"Generated junction BED: {junc_bed_path}")

    # Run aligners
    results = {}

    parallel = getattr(args, 'parallel_aligners', False)

    # Two-phase scheduler: mapPacBio is ~10x slower than minimap2/gapmm2.
    # Phase 1: mapPacBio runs alone with all threads.
    # Phase 2: remaining aligners run in parallel with equal thread shares.
    aligner_thread_counts = {a: args.threads for a in aligners}

    def _run_one_aligner(aligner):
        """Run a single aligner and return (aligner, bam_path_or_None)."""
        if aligner == 'minimap2':
            exec_path = args.minimap2_path
        elif aligner == 'mapPacBio':
            exec_path = args.mapPacBio_path
        elif aligner == 'gapmm2':
            exec_path = args.gapmm2_path
        elif aligner == 'uLTRA':
            exec_path = getattr(args, 'ultra_path', 'uLTRA')
        elif aligner == 'deSALT':
            exec_path = getattr(args, 'desalt_path', 'deSALT')
        else:
            logger.warning(f"Unknown aligner: {aligner}")
            return aligner, None

        if not check_aligner_available(exec_path):
            logger.warning(f"{aligner} not found at {exec_path}, skipping")
            return aligner, None

        import time as _time
        n_threads = aligner_thread_counts[aligner]
        output_bam = args.output_dir / f"{prefix}.{aligner}.bam"

        # Per-aligner BAM checkpoint: skip if final output already exists.
        # This lets the main pipeline skip mapPacBio when chunk BAMs have
        # already been merged by a prior chunk array job.
        if output_bam.exists():
            logger.info(f"{aligner} BAM already exists, skipping: {output_bam}")
            return aligner, str(output_bam)

        logger.info(f"Running {aligner} (threads={n_threads})...")
        _t_aligner = _time.perf_counter()

        try:
            if aligner == 'minimap2':
                run_minimap2(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                    annotation_path=str(args.annotation) if args.annotation else None,
                    junc_bonus=args.junc_bonus,
                    cache_dir=str(args.output_dir),
                )
            elif aligner == 'mapPacBio':
                _n_chunks = getattr(args, 'mapPacBio_chunks', 1) or 1
                _chunk_idx = getattr(args, 'mapPacBio_chunk_idx', None)
                run_map_pacbio(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                    chunk_idx=_chunk_idx,
                    n_chunks=_n_chunks if _n_chunks > 1 else None,
                )
            elif aligner == 'gapmm2':
                run_gapmm2(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                )
            elif aligner == 'uLTRA':
                if not args.annotation:
                    logger.warning("uLTRA requires --annotation; skipping")
                    return aligner, None
                run_ultra(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    annotation_path=str(args.annotation),
                    threads=n_threads,
                )
            elif aligner == 'deSALT':
                run_desalt(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    annotation_path=str(args.annotation) if args.annotation else None,
                    threads=n_threads,
                )

            _elapsed = _time.perf_counter() - _t_aligner
            logger.info(f"{aligner} complete: {output_bam} [{_elapsed:.1f}s]")
            return aligner, str(output_bam)

        except Exception as e:
            _elapsed = _time.perf_counter() - _t_aligner
            logger.error(f"{aligner} failed after {_elapsed:.1f}s: {e}")
            return aligner, None

    if parallel:
        from concurrent.futures import ThreadPoolExecutor, as_completed
        # Phase 1: mapPacBio alone with all threads (~10x slower than others).
        if 'mapPacBio' in aligners:
            logger.info(
                f"Running mapPacBio first with all {args.threads} threads "
                f"(two-phase: mapPacBio → rest in parallel)"
            )
            aligner_name, bam_path = _run_one_aligner('mapPacBio')
            results['mapPacBio'] = bam_path
        # Phase 2: remaining aligners with equal thread shares.
        # deSALT runs sequentially after the parallel pool — it crashes with
        # "double free or corruption" when forked inside a multithreaded process.
        remaining = [a for a in aligners if a != 'mapPacBio']
        if remaining:
            per_thread = max(1, args.threads // len(remaining))
            for a in remaining:
                aligner_thread_counts[a] = per_thread

            alloc_summary = ', '.join(f"{a}={per_thread}" for a in remaining)
            logger.info(
                f"Running {len(remaining)} aligners in parallel "
                f"(threads: {alloc_summary}, total≤{args.threads})"
            )

            parallel_batch = [a for a in remaining if a != 'deSALT']
            sequential_batch = [a for a in remaining if a == 'deSALT']

            with ThreadPoolExecutor(max_workers=max(1, len(parallel_batch))) as pool:
                futures = {pool.submit(_run_one_aligner, a): a for a in parallel_batch}
                for future in as_completed(futures):
                    aligner, bam_path = future.result()
                    results[aligner] = bam_path

            for aligner in sequential_batch:
                aligner_name, bam_path = _run_one_aligner(aligner)
                results[aligner_name] = bam_path
    else:
        for aligner in aligners:
            aligner_name, bam_path = _run_one_aligner(aligner)
            results[aligner_name] = bam_path

    # Summary of alignment step
    logger.info(f"\nAlignment summary:")
    for aligner, bam_path in results.items():
        status = "SUCCESS" if bam_path else "FAILED"
        logger.info(f"  {aligner}: {status}")
        if bam_path:
            logger.info(f"    Output: {bam_path}")

    # Check if we should run consensus selection
    successful_aligners = {k: v for k, v in results.items() if v}

    if not successful_aligners:
        logger.error("No aligners succeeded")
        return 1

    # Skip consensus if only one aligner or --no-consensus
    if len(successful_aligners) == 1 or getattr(args, 'no_consensus', False):
        logger.info("Skipping consensus selection (single aligner or --no-consensus)")
        # Copy single BAM to consensus output path, sort and index
        if len(successful_aligners) == 1:
            single_bam = list(successful_aligners.values())[0]
            rectified_bam = args.output_dir / f"{prefix}.rectified.bam"
            import shutil
            import subprocess as _sp
            threads = getattr(args, 'threads', 1)
            sorted_tmp = str(rectified_bam) + '.sorting_tmp'
            _sp.run(
                ['samtools', 'sort', '-@', str(threads), '-o', str(rectified_bam), str(single_bam)],
                check=True,
            )
            _sp.run(['samtools', 'index', str(rectified_bam)], check=True)
            logger.info(f"Single-aligner output (sorted+indexed): {rectified_bam}")
        return 0

    # Run consensus selection
    import time as _time
    _t_consensus_start = _time.perf_counter()
    logger.info(f"\nRunning consensus selection across {len(successful_aligners)} aligners...")

    from .consensus import run_consensus_selection, load_annotated_junctions

    # Load genome for consensus scoring
    _t_genome_load = _time.perf_counter()
    genome = {}
    import pysam as pysam_lib
    genome_path = str(args.genome)
    try:
        fasta = pysam_lib.FastaFile(genome_path)
    except (OSError, IOError) as e:
        # Handle gzip (not bgzip) compressed genome — auto-convert
        if genome_path.endswith('.gz'):
            import gzip as _gzip
            import subprocess as _sp
            from shutil import which as _which
            logger.warning(f"pysam cannot open {genome_path}: {e}")
            logger.warning("Attempting gzip→bgzip conversion...")
            raw_path = genome_path[:-3]  # strip .gz
            with _gzip.open(genome_path, 'rb') as _fin, open(raw_path, 'wb') as _fout:
                _fout.write(_fin.read())
            import os as _os
            _os.rename(genome_path, genome_path + '.gzip_bak')
            if _which('bgzip'):
                _sp.run(['bgzip', raw_path], check=True)
                _sp.run(['samtools', 'faidx', genome_path], check=True)
                fasta = pysam_lib.FastaFile(genome_path)
                logger.info(f"Successfully converted genome to bgzip: {genome_path}")
            else:
                # No bgzip — use uncompressed
                fasta = pysam_lib.FastaFile(raw_path)
                logger.info(f"Using uncompressed genome: {raw_path}")
        else:
            raise
    for chrom in fasta.references:
        genome[chrom] = fasta.fetch(chrom)
    fasta.close()
    logger.info(f"[TIMING] Genome load: {_time.perf_counter() - _t_genome_load:.1f}s")

    # Load annotated junctions if annotation provided
    annotated_junctions = None
    if args.annotation:
        _t_junc = _time.perf_counter()
        annotated_junctions = load_annotated_junctions(str(args.annotation))
        logger.info(f"[TIMING] Junction load: {_time.perf_counter() - _t_junc:.1f}s")

    # Run aligner selection → rectified BAM
    rectified_bam = args.output_dir / f"{prefix}.rectified.bam"

    try:
        _t_sel = _time.perf_counter()
        use_chimeric = getattr(args, 'chimeric_consensus', False)
        if use_chimeric:
            logger.info("Chimeric consensus enabled (experimental) — segments scored independently")
        stats = run_consensus_selection(
            bam_paths=successful_aligners,
            genome=genome,
            output_bam=str(rectified_bam),
            annotated_junctions=annotated_junctions,
            use_chimeric=use_chimeric,
        )
        logger.info(f"[TIMING] Aligner selection: {_time.perf_counter() - _t_sel:.1f}s")

        logger.info(f"\nRectified BAM: {rectified_bam}")
        logger.info(f"  High confidence: {stats['consensus_high']} reads")
        logger.info(f"  5' rescued: {stats['5prime_rescued']} reads")
        logger.info(f"[TIMING] Aligner selection total (incl. genome/junctions): {_time.perf_counter() - _t_consensus_start:.1f}s")

    except Exception as e:
        logger.error(f"Consensus selection failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Add MD tags via samtools calmd (required for indel correction and
    # alignment identity calculation downstream).
    logger.info("Adding MD tags with samtools calmd...")
    try:
        import subprocess as _sp
        calmd_bam = args.output_dir / f"{prefix}.rectified.md.bam"
        calmd_cmd = [
            'samtools', 'calmd', '-b',
            str(rectified_bam),
            str(args.genome),
        ]
        with open(str(calmd_bam), 'wb') as fh_out:
            result = _sp.run(calmd_cmd, stdout=fh_out, stderr=_sp.PIPE)
        if result.returncode == 0 and calmd_bam.stat().st_size > 0:
            import shutil
            calmd_bam.replace(rectified_bam)
            _sp.run(['samtools', 'index', str(rectified_bam)], check=True)
            logger.info("  MD tags added successfully")
        else:
            logger.warning(
                f"  samtools calmd failed (rc={result.returncode}); "
                "proceeding without MD tags"
            )
            if calmd_bam.exists():
                calmd_bam.unlink()
    except Exception as e:
        logger.warning(f"  samtools calmd error: {e}; proceeding without MD tags")

    return 0


def run(args: argparse.Namespace):
    """Entry point for CLI."""
    sys.exit(run_align(args))
