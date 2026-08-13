"""
Attribute-reads command for RECTIFY.

Emits a per-READ gene-attribution sidecar TSV alongside an existing corrected-reads
TSV, so a downstream viewer can keep a readthrough molecule with the gene it
INITIATED in instead of silently reassigning it to the gene it terminated in.

Deliberately a post-hoc pass over corrected TSVs, NOT part of ``correct``: an
analyze run and a browser bundle built from *different* corrections have different
3' anchors, which has already blocked one deploy.  Reading the same corrected TSV
the consumer loads makes the read set and the anchors identical by construction.

The per-gene CPA reference is built from CONTROL libraries only.  In a mutant a
gene's own 3'-end distribution already contains the readthrough being measured, so
a mutant-derived reference would be circular.

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path

logger = logging.getLogger(__name__)


def create_attribute_reads_parser(
    subparsers: argparse._SubParsersAction,
) -> argparse.ArgumentParser:
    """Create the ``attribute-reads`` subcommand parser."""
    parser = subparsers.add_parser(
        'attribute-reads',
        help="Emit a per-read gene-attribution sidecar (readthrough-aware) "
             "from existing corrected-reads TSVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        'corrected_tsv',
        type=Path,
        nargs='+',
        help='Corrected-reads TSV(s) to annotate. Must carry alignment_start/'
             'alignment_end — the raw alignment, and the only rescue-proof '
             'source of the 5\' end.',
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Directory for the <sample>.attribution.tsv sidecars',
    )

    ref = parser.add_argument_group('Reference')
    ref.add_argument(
        '--gff',
        type=Path,
        help='Annotation GFF. Omit with --Scer to use the bundled S. cerevisiae '
             'annotation.',
    )
    ref.add_argument(
        '--Scer',
        action='store_true',
        help='Use the bundled S. cerevisiae R64 annotation',
    )

    cal = parser.add_argument_group('CPA reference calibration')
    cal.add_argument(
        '--control-tsv',
        type=Path,
        nargs='+',
        help='Control/WT corrected TSVs used to build the per-gene observed-CPA '
             'reference. Defaults to ALL inputs, which is correct only when none '
             'of them is a readthrough mutant — pass this explicitly for a '
             'CPA-depletion experiment.',
    )
    cal.add_argument(
        '--min-ref-reads',
        type=int,
        default=None,
        help='Reads required before a gene\'s own modal CPA is trusted',
    )

    out = parser.add_argument_group('Output')
    out.add_argument(
        '--unit',
        choices=['reads', 'molecules'],
        default='reads',
        help='UMI-collapsed cDNA TSVs are MOLECULES, not reads — the sidecar '
             'header records this so a downstream count cannot silently change '
             'meaning.',
    )

    return parser


def run(args) -> int:
    from ..analyze import read_attribution as ra

    if args.Scer and not args.gff:
        import rectify
        args.gff = (Path(rectify.__file__).parent / 'data' /
                    'saccharomyces_cerevisiae_R64-5-1_20240529.gff')
    if not args.gff:
        logger.error("no annotation: pass --gff or --Scer")
        return 1
    if not args.gff.exists():
        logger.error("annotation not found: %s", args.gff)
        return 1

    if args.min_ref_reads is not None:
        min_reads = args.min_ref_reads
    else:
        min_reads = ra.MIN_REF_READS

    logger.info("loading annotation %s", args.gff)
    genes = ra.GeneIndex.from_gff(str(args.gff))
    logger.info("genes: %d (%d with CDS)", len(genes.spans), len(genes.cds_end))

    controls = args.control_tsv or args.corrected_tsv
    if not args.control_tsv:
        logger.warning(
            "--control-tsv not given: building the CPA reference from ALL %d "
            "inputs. If any of them is a readthrough mutant this is CIRCULAR — "
            "the gene's own 3'-end distribution already contains the readthrough "
            "being measured.", len(controls))

    ref = ra.build_cpa_reference([str(p) for p in controls], genes,
                                 min_reads=min_reads)
    if not ref:
        logger.error("CPA reference is empty — no gene reached %d reads. "
                     "Check that the control TSVs and the annotation use the "
                     "same chromosome names.", min_reads)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    labels = [Path(p).parent.name or Path(p).stem for p in controls]

    total = 0
    for tsv in args.corrected_tsv:
        sample = tsv.parent.name or tsv.stem
        out = args.output_dir / f"{sample}.attribution.tsv"
        try:
            n = ra.write_attribution_sidecar(
                str(tsv), str(out), genes, ref,
                unit=args.unit, control_labels=labels)
        except ValueError as exc:
            logger.error("%s: %s", tsv, exc)
            return 1
        # A zero-row sidecar joins to nothing downstream while still looking like
        # a delivered file; fail rather than ship it.
        if n == 0:
            logger.error("%s produced 0 rows — refusing to ship an empty sidecar", out)
            return 1
        logger.info("%s -> %s (%d %s)", sample, out.name, n, args.unit)
        total += n

    logger.info("done: %d sidecars, %d %s", len(args.corrected_tsv), total, args.unit)
    return 0
