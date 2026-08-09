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
        default=None,
        metavar='FASTA',
        help='Reference genome FASTA (.fsa or .fsa.gz). '
             'Optional when --Scer or --organism is set.',
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
    junc.add_argument(
        '--junction-min-support',
        type=int,
        default=1,
        metavar='N',
        help='Minimum total observed read support across all --aligner-bams '
             'required for non-annotated junctions to enter junction_pool.pkl. '
             'Annotated junctions are always retained. Raising this filters '
             'one-off noisy N-op junctions from large DRS pools.',
    )
    junc.add_argument(
        '--junction-max-size',
        type=int,
        default=None,
        metavar='BP',
        help='Optional maximum observed intron/junction length to include from '
             '--aligner-bams. Annotated junctions are always retained. For '
             'S. cerevisiae, 10000 is a conservative organism-tuned cap.',
    )
    junc.add_argument(
        '--complexity-alpha',
        type=float,
        default=None,
        metavar='ALPHA',
        help='Structural pool-admission gate (planning/648/649): drop a NOVEL '
             'junction when its observed span D exceeds what its worse 15-nt '
             'exonic flank can distinguish from chance '
             '(D * 2^-I_eff > ALPHA). Refuses long-range junctions anchored on '
             'low-complexity/homopolymer flanks — the class that inflates '
             'correct-stage cost panel-wide. Annotated junctions are always '
             'retained. Default: off (no filtering). Suggested: 0.01.',
    )

    from rectify.data import add_organism_args
    add_organism_args(parser)

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

    # Loaded on first use and shared by both steps (the complexity gate in
    # step 2 needs it too, and it is the most expensive thing to load twice).
    _genome_cache = {}

    def _get_genome():
        if 'g' not in _genome_cache:
            from ...utils.genome import load_genome
            _t = time.perf_counter()
            logger.info("  Loading genome...")
            _genome_cache['g'] = load_genome(args.genome)
            logger.info(f"  Genome loaded in {time.perf_counter() - _t:.1f}s")
        return _genome_cache['g']

    # ------------------------------------------------------------------
    # Step 1: Variant scan (Module 2D Pass 1)
    # ------------------------------------------------------------------
    if not args.skip_variant_scan:
        logger.info("=== Step 1: Variant scan (Module 2D) ===")
        logger.info(f"  BAM: {args.bam}")

        genome = _get_genome()

        from ..bam.variant_scan import run_variant_aware_scan
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
            from ..consensus.consensus import load_annotated_junctions
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

        from ..splice.junction_refiner import build_junction_pool
        t0 = time.perf_counter()
        all_junctions, annotated_set = build_junction_pool(
            aligner_bams,
            annotated_junctions,
            min_observed_support=args.junction_min_support,
            max_junction_size=args.junction_max_size,
        )
        elapsed = time.perf_counter() - t0
        logger.info(
            f"  Junction pool built in {elapsed:.1f}s: "
            f"{len(all_junctions)} total, {len(annotated_set)} annotated"
        )

        # Structural complexity gate (planning/649). Off unless requested, so
        # the default pool is byte-identical to before this flag existed.
        if args.complexity_alpha is not None:
            from ..splice.overhang_informativeness import (
                filter_pool_by_flank_complexity,
            )
            t0 = time.perf_counter()
            _before = len(all_junctions)
            all_junctions, _n_refused = filter_pool_by_flank_complexity(
                all_junctions, annotated_set,
                _get_genome(), alpha=args.complexity_alpha,
            )
            logger.info(
                "  Complexity gate (alpha=%g): refused %d of %d novel junctions "
                "in %.1fs; pool %d -> %d",
                args.complexity_alpha, _n_refused,
                _before - len(annotated_set), time.perf_counter() - t0,
                _before, len(all_junctions),
            )

        pool_data = {
            'all_junctions': all_junctions,
            'annotated_set': annotated_set,
            'min_observed_support': args.junction_min_support,
            'max_junction_size': args.junction_max_size,
            'complexity_alpha': args.complexity_alpha,
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
