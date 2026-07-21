"""``rectify umi-dedup`` — collapse a UMI-tagged short-read BAM to one record per molecule.

Runs after alignment + consensus selection (the point at which each read has a
single unambiguous coordinate). Reads the SAM-standard ``RX`` UMI tag (attached
upstream by ``rectify split --umi`` and carried through the COMPASS aligners),
groups reads by (fragment span x UMI) with directional clustering, and writes a
deduplicated BAM in which every downstream count is a MOLECULE count.

Also writes the family-size and within-family UMI edit-distance distributions to
a JSON sidecar — the empirical inputs for calibrating the edit-distance threshold
and deciding whether consensus is worthwhile (design doc S5).
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def create_umi_dedup_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        'umi-dedup',
        help='Collapse a UMI-tagged (RX) short-read BAM to one record per molecule.',
        description=(
            'UMI-aware deduplication for short-read (COMPASS) BAMs. Groups reads by '
            'fragment span x UMI (directional clustering) and writes a BAM with one '
            'record per molecule (tags MI=molecule id, cD=family size). Run AFTER '
            'consensus selection on a coordinate-sorted, RX-tagged BAM.'
        ),
    )
    parser.add_argument('input_bam', type=Path,
                        help='Coordinate-sorted BAM with RX UMI tags.')
    parser.add_argument('output_bam', type=Path,
                        help='Output deduplicated BAM.')
    parser.add_argument('--umi-tag', default='RX',
                        help='BAM tag holding the UMI (default: RX).')
    parser.add_argument('--umi-edit-distance', type=int, default=1,
                        help='Max UMI edit distance within a fragment bucket (default 1).')
    parser.add_argument('--umi-clustering', default='directional',
                        choices=['directional', 'components'],
                        help='UMI clustering method (default: directional).')
    parser.add_argument('--umi-collapse', default='dedup',
                        choices=['dedup', 'consensus'],
                        help='Collapse mode (default: dedup; consensus is a planned extension).')
    parser.add_argument('--antisense', action='store_true',
                        help='Library is ANTISENSE (R1 opposite the mRNA), e.g. TruSeq-dUTP / '
                             'QuantSeq-REV. CORALL is SENSE, so OMIT this for CORALL.')
    parser.add_argument('--stats-json', type=Path, default=None,
                        help='Write the family-size / edit-distance distributions here '
                             '(default: <output_bam>.umi_stats.json).')
    parser.add_argument('--verbose', action='store_true', help='Verbose logging')
    return parser


def run(args: argparse.Namespace) -> int:
    from rectify.core.umi.dedup import dedup_bam

    logging.basicConfig(
        level=logging.DEBUG if getattr(args, 'verbose', False) else logging.INFO,
        format='%(asctime)s %(levelname)s %(message)s',
    )

    if not args.input_bam.exists():
        logger.error("Input BAM not found: %s", args.input_bam)
        return 1
    args.output_bam.parent.mkdir(parents=True, exist_ok=True)

    try:
        stats = dedup_bam(
            str(args.input_bam),
            str(args.output_bam),
            umi_tag=args.umi_tag,
            edit_distance=args.umi_edit_distance,
            clustering=args.umi_clustering,
            r1_sense=not args.antisense,
            collapse=args.umi_collapse,
        )
    except NotImplementedError as e:
        logger.error(str(e))
        return 1

    stats_path = args.stats_json or args.output_bam.with_suffix(
        args.output_bam.suffix + '.umi_stats.json')
    with open(stats_path, 'w') as fh:
        json.dump(stats.as_dict(), fh, indent=2)

    logger.info(
        "Dedup: %d fragments -> %d molecules (%.1f%% duplicates). Stats: %s",
        stats.n_input_fragments, stats.n_molecules,
        100 * stats.duplication_rate, stats_path,
    )
    # index the output for downstream fetch (best-effort)
    try:
        import pysam
        pysam.index(str(args.output_bam))
    except Exception as e:  # noqa: BLE001
        logger.warning("Could not index %s: %s", args.output_bam, e)
    return 0
