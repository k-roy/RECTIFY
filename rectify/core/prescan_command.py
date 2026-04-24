"""
rectify prescan — pre-compute variant scan and junction pool for chunked correction.

Outputs two pickle files used by per-chunk `rectify correct` runs:
  - rescue_scan.pkl    : VariantAwareHomopolymerRescue (Module 2D, Pass 1)
  - junction_pool.pkl  : {'all_junctions': Set, 'annotated_set': Set} (Module 2H)

Usage:
    rectify prescan \\
        --bam merged.bam \\
        --aligner-bams minimap2.bam mapPacBio.bam \\
        --annotation genes.gff.gz \\
        --genome genome.fsa.gz \\
        -o prescan_cache/

Then pass to each chunk:
    rectify correct chunk_001.bam \\
        --variant-scan-cache prescan_cache/rescue_scan.pkl \\
        --junction-pool-cache prescan_cache/junction_pool.pkl \\
        ...
"""

from __future__ import annotations

import argparse
import logging
import pickle
import sys
import time
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def create_prescan_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register `rectify prescan` subparser."""
    parser = subparsers.add_parser(
        'prescan',
        help='Pre-compute variant scan and junction pool for chunked correction',
        description=(
            'Pre-compute the two cross-chunk data structures required by '
            '`rectify correct`: the Module 2D variant frequency table '
            '(rescue_scan.pkl) and the Module 2H junction pool '
            '(junction_pool.pkl). Run once on the merged/full BAM before '
            'fanning out to per-chunk correction jobs.'
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    req = parser.add_argument_group('Required inputs')
    req.add_argument(
        '--bam',
        type=str,
        required=True,
        metavar='BAM',
        help='Merged aligned BAM to scan for variant frequencies (Module 2D). '
             'Should contain all reads across all chunks.',
    )
    req.add_argument(
        '--genome',
        type=str,
        required=True,
        metavar='FASTA',
        help='Reference genome FASTA (.fsa or .fsa.gz).',
    )
    req.add_argument(
        '-o', '--output-dir',
        type=str,
        required=True,
        metavar='DIR',
        help='Directory to write rescue_scan.pkl and junction_pool.pkl.',
    )

    junc = parser.add_argument_group('Junction pool (Module 2H)')
    junc.add_argument(
        '--aligner-bams',
        action='append',
        type=str,
        metavar='BAM',
        default=[],
        help='Per-aligner BAM files for junction pool construction '
             '(repeat for each aligner). Accepts plain paths or '
             '\'aligner:path\' pairs. If omitted, only annotated junctions '
             'are included in the pool.',
    )
    junc.add_argument(
        '--annotation',
        type=str,
        default=None,
        metavar='GFF',
        help='Gene annotation GFF/GFF3 (.gff or .gff.gz). Used to populate '
             'the annotated junction set in the pool.',
    )

    opt = parser.add_argument_group('Options')
    opt.add_argument(
        '--skip-variant-scan',
        action='store_true',
        default=False,
        help='Skip the variant scan step; only produce junction_pool.pkl.',
    )
    opt.add_argument(
        '--skip-junction-pool',
        action='store_true',
        default=False,
        help='Skip the junction pool step; only produce rescue_scan.pkl.',
    )
    opt.add_argument(
        '--threads',
        type=int,
        default=4,
        metavar='N',
        help='Threads for genome loading.',
    )
    opt.add_argument(
        '-v', '--verbose',
        action='store_true',
        default=False,
        help='Verbose logging.',
    )

    return parser


def _strip_aligner_prefix(path_str: str) -> str:
    """Strip 'aligner:' prefix from BAM path (same logic as correct_command)."""
    if ':' in path_str:
        prefix, rest = path_str.split(':', 1)
        if '/' not in prefix:
            return rest
    return path_str


def run(args: argparse.Namespace) -> int:
    """Entry point for `rectify prescan`."""
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s %(levelname)s %(name)s: %(message)s',
        stream=sys.stderr,
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rescue_scan_path = output_dir / 'rescue_scan.pkl'
    junction_pool_path = output_dir / 'junction_pool.pkl'

    t0_total = time.perf_counter()

    # ------------------------------------------------------------------
    # Step 1: Variant scan (Module 2D Pass 1)
    # ------------------------------------------------------------------
    if not args.skip_variant_scan:
        logger.info("=== Step 1: Variant scan (Module 2D) ===")
        logger.info(f"  BAM: {args.bam}")

        from ..utils.genome import load_genome
        t0 = time.perf_counter()
        logger.info("  Loading genome...")
        genome = load_genome(args.genome)
        logger.info(f"  Genome loaded in {time.perf_counter() - t0:.1f}s")

        from .bam_processor import run_variant_aware_scan
        t0 = time.perf_counter()
        logger.info("  Running variant scan...")
        variant_rescue = run_variant_aware_scan(
            bam_path=args.bam,
            genome=genome,
            min_variant_fraction=0.8,
            min_reads_for_variant_call=5,
        )
        elapsed = time.perf_counter() - t0
        logger.info(f"  Variant scan complete in {elapsed:.1f}s")

        with open(rescue_scan_path, 'wb') as fh:
            pickle.dump(variant_rescue, fh, protocol=pickle.HIGHEST_PROTOCOL)
        logger.info(f"  Saved: {rescue_scan_path}")
    else:
        logger.info("Skipping variant scan (--skip-variant-scan)")

    # ------------------------------------------------------------------
    # Step 2: Junction pool (Module 2H)
    # ------------------------------------------------------------------
    if not args.skip_junction_pool:
        logger.info("=== Step 2: Junction pool (Module 2H) ===")

        # Load annotated junctions
        annotated_junctions: set = set()
        if args.annotation:
            from .consensus import load_annotated_junctions
            logger.info(f"  Loading annotated junctions from {args.annotation}")
            annotated_junctions = load_annotated_junctions(args.annotation)
            logger.info(f"  {len(annotated_junctions)} annotated junctions loaded")
        else:
            logger.info("  No --annotation provided; pool will contain aligner BAM junctions only")

        # Strip aligner prefixes from BAM paths
        aligner_bams = [_strip_aligner_prefix(p) for p in (args.aligner_bams or [])]
        if aligner_bams:
            logger.info(f"  Aligner BAMs: {aligner_bams}")
        else:
            logger.info("  No aligner BAMs provided; pool will contain annotated junctions only")

        from .junction_refiner import build_junction_pool
        t0 = time.perf_counter()
        all_junctions, annotated_set = build_junction_pool(aligner_bams, annotated_junctions)
        elapsed = time.perf_counter() - t0
        logger.info(
            f"  Junction pool built in {elapsed:.1f}s: "
            f"{len(all_junctions)} total, {len(annotated_set)} annotated"
        )

        pool_data = {
            'all_junctions': all_junctions,
            'annotated_set': annotated_set,
        }
        with open(junction_pool_path, 'wb') as fh:
            pickle.dump(pool_data, fh, protocol=pickle.HIGHEST_PROTOCOL)
        logger.info(f"  Saved: {junction_pool_path}")
    else:
        logger.info("Skipping junction pool (--skip-junction-pool)")

    logger.info(f"=== prescan complete in {time.perf_counter() - t0_total:.1f}s ===")
    logger.info(f"  Output dir: {output_dir}")
    if not args.skip_variant_scan:
        logger.info(f"  rescue_scan.pkl    -> {rescue_scan_path}")
    if not args.skip_junction_pool:
        logger.info(f"  junction_pool.pkl  -> {junction_pool_path}")

    return 0
