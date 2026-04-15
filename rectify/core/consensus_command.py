"""
Consensus Command for RECTIFY.

Runs consensus aligner selection on pre-aligned per-aligner BAMs.
Used after merging chunked alignment results.

Accepts BAMs as "aligner:path" pairs, e.g.:
  rectify consensus \\
      minimap2:/path/to/sample.minimap2.sorted.bam \\
      mapPacBio:/path/to/sample.mapPacBio.sorted.bam \\
      gapmm2:/path/to/sample.gapmm2.sorted.bam \\
      uLTRA:/path/to/sample.uLTRA.sorted.bam \\
      deSALT:/path/to/sample.deSALT.sorted.bam \\
      --genome genome.fa \\
      --annotation genes.gff \\
      -o outdir/

Output: <outdir>/<prefix>.consensus.bam (coordinate-sorted, indexed, MD-tagged)

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path

logger = logging.getLogger(__name__)


def create_consensus_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Create consensus subcommand parser."""
    parser = subparsers.add_parser(
        'consensus',
        help='Run consensus aligner selection on pre-aligned per-aligner BAMs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Examples:
  # Select best aligner per read from 5 pre-aligned BAMs
  rectify consensus \\
      minimap2:sample.minimap2.sorted.bam \\
      mapPacBio:sample.mapPacBio.sorted.bam \\
      gapmm2:sample.gapmm2.sorted.bam \\
      uLTRA:sample.uLTRA.sorted.bam \\
      deSALT:sample.deSALT.sorted.bam \\
      --genome genome.fa --annotation genes.gff \\
      -o results/sample/
        """
    )

    parser.add_argument(
        'bams',
        nargs='+',
        metavar='ALIGNER:BAM',
        help=(
            'Per-aligner BAM files in "aligner:path" format. '
            'Accepted aligner names: minimap2, mapPacBio, gapmm2, uLTRA, deSALT.'
        )
    )

    parser.add_argument(
        '--genome',
        type=Path,
        required=True,
        help='Reference genome FASTA file (used for MD-tag recalculation)'
    )

    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory'
    )

    parser.add_argument(
        '--prefix',
        default='',
        help='Output file prefix (default: derived from first BAM filename)'
    )

    parser.add_argument(
        '--annotation',
        type=Path,
        help='Gene annotation GFF/GTF for annotated-junction scoring in consensus selection'
    )

    parser.add_argument(
        '--chimeric',
        action='store_true',
        default=False,
        help=(
            'Use chimeric consensus selection: independently pick the best aligner '
            'for each read segment, then assemble a chimeric CIGAR. Experimental.'
        )
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose logging'
    )

    return parser


def _parse_bam_args(bam_args: list) -> dict:
    """
    Parse ["aligner:path", ...] into {aligner: path, ...}.

    Raises ValueError for malformed or missing entries.
    """
    result = {}
    valid_aligners = {'minimap2', 'mapPacBio', 'gapmm2', 'uLTRA', 'deSALT'}
    for token in bam_args:
        if ':' not in token:
            raise ValueError(
                f"BAM argument must be in 'aligner:path' format, got: {token!r}"
            )
        aligner, path = token.split(':', 1)
        if aligner not in valid_aligners:
            raise ValueError(
                f"Unknown aligner {aligner!r}. Valid names: {sorted(valid_aligners)}"
            )
        bam_path = Path(path)
        if not bam_path.exists():
            raise FileNotFoundError(f"BAM file not found for {aligner}: {bam_path}")
        result[aligner] = str(bam_path)
    return result


def run_consensus(args: argparse.Namespace) -> int:
    """Run consensus selection command."""
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Parse BAM arguments
    try:
        bam_paths = _parse_bam_args(args.bams)
    except (ValueError, FileNotFoundError) as e:
        logger.error(str(e))
        return 1

    if len(bam_paths) < 2:
        logger.error(
            f"Consensus selection requires at least 2 aligner BAMs; "
            f"got {len(bam_paths)}. Use 'rectify align' for single-aligner output."
        )
        return 1

    if not args.genome.exists():
        logger.error(f"Genome file not found: {args.genome}")
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Determine output prefix
    prefix = args.prefix
    if not prefix:
        first_bam = Path(list(bam_paths.values())[0])
        name = first_bam.name
        for aligner in ('minimap2', 'mapPacBio', 'gapmm2', 'uLTRA', 'deSALT'):
            suffix = f'.{aligner}.sorted.bam'
            if name.endswith(suffix):
                prefix = name[:-len(suffix)]
                break
            suffix = f'.{aligner}.bam'
            if name.endswith(suffix):
                prefix = name[:-len(suffix)]
                break
        else:
            prefix = first_bam.stem

    logger.info(
        f"Consensus selection: {len(bam_paths)} aligners "
        f"({', '.join(bam_paths.keys())})"
    )
    logger.info(f"Output prefix: {prefix}")

    import time as _time

    # Load genome
    import pysam as pysam_lib
    logger.info(f"Loading genome: {args.genome}")
    _t = _time.perf_counter()
    genome = {}
    genome_path = str(args.genome)
    try:
        fasta = pysam_lib.FastaFile(genome_path)
    except (OSError, IOError) as e:
        if genome_path.endswith('.gz'):
            import gzip as _gzip
            import subprocess as _sp
            from shutil import which as _which
            logger.warning(f"pysam cannot open {genome_path}: {e}")
            logger.warning("Attempting gzip→bgzip conversion...")
            raw_path = genome_path[:-3]
            with _gzip.open(genome_path, 'rb') as fin, open(raw_path, 'wb') as fout:
                fout.write(fin.read())
            import os as _os
            _os.rename(genome_path, genome_path + '.gzip_bak')
            if _which('bgzip'):
                _sp.run(['bgzip', raw_path], check=True)
                _sp.run(['samtools', 'faidx', genome_path], check=True)
                fasta = pysam_lib.FastaFile(genome_path)
                logger.info(f"Converted genome to bgzip: {genome_path}")
            else:
                fasta = pysam_lib.FastaFile(raw_path)
                logger.info(f"Using uncompressed genome: {raw_path}")
        else:
            raise
    for chrom in fasta.references:
        genome[chrom] = fasta.fetch(chrom)
    fasta.close()
    logger.info(f"[TIMING] Genome load: {_time.perf_counter() - _t:.1f}s")

    # Load annotated junctions
    annotated_junctions = None
    if args.annotation:
        if not args.annotation.exists():
            logger.warning(f"Annotation not found: {args.annotation} — skipping junction scoring")
        else:
            from .consensus import load_annotated_junctions
            _t = _time.perf_counter()
            annotated_junctions = load_annotated_junctions(str(args.annotation))
            logger.info(f"[TIMING] Junction load: {_time.perf_counter() - _t:.1f}s")

    # Run consensus selection
    from .consensus import run_consensus_selection
    output_bam = args.output_dir / f"{prefix}.consensus.bam"

    logger.info(f"Running consensus selection → {output_bam}")
    _t = _time.perf_counter()
    try:
        use_chimeric = getattr(args, 'chimeric', False)
        if use_chimeric:
            logger.info("Chimeric consensus enabled (experimental)")
        stats = run_consensus_selection(
            bam_paths=bam_paths,
            genome=genome,
            output_bam=str(output_bam),
            annotated_junctions=annotated_junctions,
            use_chimeric=use_chimeric,
        )
    except Exception as e:
        logger.error(f"Consensus selection failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

    logger.info(f"[TIMING] Consensus selection: {_time.perf_counter() - _t:.1f}s")
    logger.info(f"  High confidence: {stats.get('consensus_high', 'N/A')}")
    logger.info(f"  5' rescued: {stats.get('5prime_rescued', 'N/A')}")

    # Add MD tags
    logger.info("Adding MD tags with samtools calmd...")
    try:
        import subprocess as _sp
        calmd_bam = args.output_dir / f"{prefix}.consensus.md.bam"
        with open(str(calmd_bam), 'wb') as fh_out:
            result = _sp.run(
                ['samtools', 'calmd', '-b', str(output_bam), str(args.genome)],
                stdout=fh_out, stderr=_sp.PIPE,
            )
        if result.returncode == 0 and calmd_bam.stat().st_size > 0:
            calmd_bam.replace(output_bam)
            _sp.run(['samtools', 'index', str(output_bam)], check=True)
            logger.info("  MD tags added")
        else:
            logger.warning(
                f"  samtools calmd failed (rc={result.returncode}); proceeding without MD tags"
            )
            if calmd_bam.exists():
                calmd_bam.unlink()
    except Exception as e:
        logger.warning(f"  samtools calmd error: {e}; proceeding without MD tags")

    logger.info(f"Output: {output_bam}")
    return 0


def run(args: argparse.Namespace):
    """Entry point for CLI."""
    sys.exit(run_consensus(args))
