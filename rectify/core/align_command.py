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
- Outputs single consensus BAM with best alignment per read
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
        choices=['minimap2', 'mapPacBio', 'gapmm2', 'all'],
        default=['all'],
        help='Aligners to run (default: all for consensus)'
    )

    aligner_group.add_argument(
        '--no-consensus',
        action='store_true',
        help='Skip consensus selection, output separate BAMs per aligner'
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
        '--gapmm2-path',
        default='gapmm2',
        help='Path to gapmm2 executable'
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

    # Expand 'all' to list of aligners
    aligners = args.aligners
    if 'all' in aligners:
        aligners = ['minimap2', 'mapPacBio', 'gapmm2']

    # Import multi-aligner functions
    from .multi_aligner import (
        run_minimap2,
        run_map_pacbio,
        run_gapmm2,
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

    # Relative thread weights: minimap2 is ~2x faster than mapPacBio/gapmm2.
    # In parallel mode, allocate proportionally so all three finish at roughly
    # the same time rather than giving minimap2 idle cores it doesn't need.
    # weights sum to 5 → for 16 threads: minimap2=3, mapPacBio=6, gapmm2=6 (=15)
    _THREAD_WEIGHTS = {'minimap2': 1, 'mapPacBio': 2, 'gapmm2': 2}

    if parallel:
        total_weight = sum(_THREAD_WEIGHTS.get(a, 1) for a in aligners)
        aligner_thread_counts = {
            a: max(1, round(args.threads * _THREAD_WEIGHTS.get(a, 1) / total_weight))
            for a in aligners
        }
    else:
        aligner_thread_counts = {a: args.threads for a in aligners}

    def _run_one_aligner(aligner):
        """Run a single aligner and return (aligner, bam_path_or_None)."""
        if aligner == 'minimap2':
            exec_path = args.minimap2_path
        elif aligner == 'mapPacBio':
            exec_path = args.mapPacBio_path
        elif aligner == 'gapmm2':
            exec_path = args.gapmm2_path
        else:
            logger.warning(f"Unknown aligner: {aligner}")
            return aligner, None

        if not check_aligner_available(exec_path):
            logger.warning(f"{aligner} not found at {exec_path}, skipping")
            return aligner, None

        n_threads = aligner_thread_counts[aligner]
        output_bam = args.output_dir / f"{prefix}.{aligner}.bam"
        logger.info(f"Running {aligner} (threads={n_threads})...")

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
                run_map_pacbio(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                )
            elif aligner == 'gapmm2':
                run_gapmm2(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                )

            logger.info(f"{aligner} complete: {output_bam}")
            return aligner, str(output_bam)

        except Exception as e:
            logger.error(f"{aligner} failed: {e}")
            return aligner, None

    if parallel:
        alloc_summary = ', '.join(
            f"{a}={aligner_thread_counts[a]}" for a in aligners
        )
        logger.info(
            f"Running {len(aligners)} aligners in parallel "
            f"(threads: {alloc_summary}, total={args.threads})"
        )
        from concurrent.futures import ThreadPoolExecutor, as_completed
        with ThreadPoolExecutor(max_workers=len(aligners)) as pool:
            futures = {pool.submit(_run_one_aligner, a): a for a in aligners}
            for future in as_completed(futures):
                aligner, bam_path = future.result()
                results[aligner] = bam_path
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
        # Copy single BAM to consensus output path
        if len(successful_aligners) == 1:
            single_bam = list(successful_aligners.values())[0]
            consensus_bam = args.output_dir / f"{prefix}.consensus.bam"
            import shutil
            shutil.copy(single_bam, str(consensus_bam))
            shutil.copy(f"{single_bam}.bai", f"{consensus_bam}.bai")
            logger.info(f"Single-aligner output: {consensus_bam}")
        return 0

    # Run consensus selection
    logger.info(f"\nRunning consensus selection across {len(successful_aligners)} aligners...")

    from .consensus import run_consensus_selection, load_annotated_junctions

    # Load genome for consensus scoring
    genome = {}
    import pysam as pysam_lib
    fasta = pysam_lib.FastaFile(str(args.genome))
    for chrom in fasta.references:
        genome[chrom] = fasta.fetch(chrom)
    fasta.close()

    # Load annotated junctions if annotation provided
    annotated_junctions = None
    if args.annotation:
        annotated_junctions = load_annotated_junctions(str(args.annotation))

    # Run consensus
    consensus_bam = args.output_dir / f"{prefix}.consensus.bam"

    try:
        stats = run_consensus_selection(
            bam_paths=successful_aligners,
            genome=genome,
            output_bam=str(consensus_bam),
            annotated_junctions=annotated_junctions,
        )

        logger.info(f"\nConsensus BAM: {consensus_bam}")
        logger.info(f"  High confidence: {stats['consensus_high']} reads")
        logger.info(f"  5' rescued: {stats['5prime_rescued']} reads")

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
        calmd_bam = args.output_dir / f"{prefix}.consensus.md.bam"
        calmd_cmd = [
            'samtools', 'calmd', '-b',
            str(consensus_bam),
            str(args.genome),
        ]
        with open(str(calmd_bam), 'wb') as fh_out:
            result = _sp.run(calmd_cmd, stdout=fh_out, stderr=_sp.PIPE)
        if result.returncode == 0 and calmd_bam.stat().st_size > 0:
            import shutil
            calmd_bam.replace(consensus_bam)
            _sp.run(['samtools', 'index', str(consensus_bam)], check=True)
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
