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
                        help='dedup (keep the best fragment per molecule) or consensus (one consensus '
                             'fragment per UMI family; see the consensus options below).')
    parser.add_argument('--umi-from-qname', action='store_true',
                        help='When the UMI tag is absent, take the UMI from the QNAME suffix after the last '
                             "'_' (umi_tools extract convention).")
    cg = parser.add_argument_group('consensus options (--umi-collapse consensus)')
    cg.add_argument('--reference', type=Path, default=None,
                    help='Reference FASTA (indexed). Needed to resolve "=" bases (samtools calmd -e) and to '
                         'classify mismatches; without it classification is skipped.')
    cg.add_argument('--span-key', default='unclipped', choices=['unclipped', 'tlen'],
                    help='Molecule key: both soft-clip-corrected 5\' ends (default; robust to primer bases '
                         'soft-clipped differently copy to copy) or R1 5\' + |TLEN| (the dedup key). '
                         'Molecule counts under the other key are reported too. Needs MC tags (samtools fixmate).')
    cg.add_argument('--primer-mate', default='r2', choices=['r1', 'r2'],
                    help='Mate whose first bases are the random RT primer (Lexogen CORALL V2: r2).')
    cg.add_argument('--primer-length', type=int, default=0,
                    help='Random-primer length; 0 disables the primer features (Lexogen CORALL V2: 9).')
    cg.add_argument('--primer-source', default='seq', choices=['seq', 'tag'],
                    help="seq: the primer bases are still in the reads; tag: they were trimmed before "
                         "alignment and carried as pR:Z (sequence, read orientation) and pQ:Z (qualities).")
    cg.add_argument('--primer-action', default='mask', choices=['mask', 'keep'],
                    help='mask (default): set the primer positions of every consensus mate to N/Q2.')
    cg.add_argument('--primer-min-support', type=int, default=2,
                    help='Members per raw UMI needed before a primer difference splits a UMI cluster (default 2).')
    cg.add_argument('--min-posterior', type=float, default=0.9,
                    help='Consensus posterior below which a base becomes N (default 0.9).')
    cg.add_argument('--recurrence-check', action='store_true',
                    help='Pile up every molecule-level mismatch site in the output and report how many are '
                         'shared by most molecules (strain variant or systematic error).')
    cg.add_argument('--threads', type=int, default=1, help='Threads for the samtools sorts.')
    cg.add_argument('--tmp-dir', type=Path, default=None, help='Directory for the intermediate BAMs.')
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

    stats_path = args.stats_json or args.output_bam.with_suffix(
        args.output_bam.suffix + '.umi_stats.json')
    if args.umi_collapse == 'consensus':
        from rectify.core.umi.consensus import consensus_bam, write_stats
        cstats = consensus_bam(
            str(args.input_bam), str(args.output_bam),
            reference=str(args.reference) if args.reference else None,
            umi_tag=args.umi_tag, umi_from_qname=args.umi_from_qname,
            edit_distance=args.umi_edit_distance, clustering=args.umi_clustering,
            r1_sense=not args.antisense, span_key=args.span_key,
            primer_mate=args.primer_mate, primer_length=args.primer_length,
            primer_source=args.primer_source, primer_action=args.primer_action,
            primer_min_support=args.primer_min_support, min_posterior=args.min_posterior,
            recurrence_check=args.recurrence_check, threads=args.threads,
            tmp_dir=str(args.tmp_dir) if args.tmp_dir else None,
        )
        write_stats(cstats, str(stats_path))
        logger.info("Consensus: %d fragments -> %d molecules (%d under the other key). Stats: %s",
                    cstats.n_input_fragments, cstats.n_molecules, cstats.n_molecules_other_key, stats_path)
        return 0

    try:
        stats = dedup_bam(
            str(args.input_bam),
            str(args.output_bam),
            umi_tag=args.umi_tag,
            edit_distance=args.umi_edit_distance,
            clustering=args.umi_clustering,
            r1_sense=not args.antisense,
            collapse=args.umi_collapse,
            umi_from_qname=args.umi_from_qname,
        )
    except NotImplementedError as e:
        logger.error(str(e))
        return 1
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
