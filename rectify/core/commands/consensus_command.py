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

Short-read (Illumina / QuantSeq) correct-before-consensus workflow:
  rectify align --short-read --no-consensus sample.fastq.gz --Scer -o outdir/
  rectify correct --short-read --dT-primed-cDNA outdir/sample.bbmap.bam --write-corrected-bam outdir/sample.bbmap.corrected.bam -o outdir/sample.bbmap.corrected_reads.tsv
  rectify correct --short-read --dT-primed-cDNA outdir/sample.bwa.bam --write-corrected-bam outdir/sample.bwa.corrected.bam -o outdir/sample.bwa.corrected_reads.tsv
  rectify consensus \\
      bbmap:outdir/sample.bbmap.corrected.bam \\
      bwa:outdir/sample.bwa.corrected.bam \\
      --genome genome.fa --annotation genes.gff \\
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
            'Accepted aligner names: minimap2, mapPacBio, gapmm2, uLTRA, deSALT (long-read); bbmap, bwa (short-read).'
        )
    )

    parser.add_argument(
        '--genome',
        type=Path,
        default=None,
        help='Reference genome FASTA file (used for MD-tag recalculation). '
             'Optional when --Scer or --organism is set.'
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

    from rectify.data import add_organism_args
    add_organism_args(parser)

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


    parser.add_argument(
        '--read-num-sidecar',
        type=Path,
        default=None,
        help=(
            'Optional <sample>.read_num_sidecar.parquet from rectify split. '
            'If omitted, RECTIFY looks beside the input BAMs for a matching sidecar.'
        )
    )

    parser.add_argument(
        '--no-bedgraph',
        action='store_true',
        default=False,
        help=(
            'Skip bedGraph/bigWig generation after consensus selection. '
            'By default, strand-specific 3\' and 5\' end coverage files are written '
            'from corrected_reads.tsv if it is present in the output directory.'
        )
    )

    return parser


def _parse_bam_args(bam_args: list) -> dict:
    """
    Parse ["aligner:path", ...] into {aligner: path, ...}.

    Raises ValueError for malformed or missing entries.
    """
    result = {}
    valid_aligners = {'minimap2', 'mapPacBio', 'gapmm2', 'uLTRA', 'deSALT', 'bbmap', 'bwa'}
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


def _autodetect_read_num_sidecar(bam_paths: dict, prefix: str):
    """Find a read-num sidecar beside the input BAMs, if one is unambiguous."""
    candidates = []
    seen = set()
    for path_str in bam_paths.values():
        parent = Path(path_str).parent
        preferred = parent / f"{prefix}.read_num_sidecar.parquet"
        if preferred.exists() and preferred not in seen:
            candidates.append(preferred)
            seen.add(preferred)
        for p in parent.glob("*.read_num_sidecar.parquet"):
            if p not in seen:
                candidates.append(p)
                seen.add(p)
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    logger.warning(
        "Multiple read-num sidecars found; pass --read-num-sidecar explicitly. Candidates: %s",
        ', '.join(str(p) for p in candidates),
    )
    return None


def _generate_bedgraphs(output_dir: Path, prefix: str, genome: dict) -> None:
    """
    Write strand-specific 3' and 5' end bigWig (or bedGraph) files from corrected_reads.tsv.

    Looks for corrected_reads.tsv in output_dir (written by run_final_merge before consensus
    is called). Generates four files per run:
      {prefix}.3prime.plus.bw / {prefix}.3prime.minus.bw
      {prefix}.5prime.plus.bw / {prefix}.5prime.minus.bw
    Falls back to .bedgraph if pyBigWig is not available.
    """
    import pandas as pd
    from .export_command import aggregate_positions, write_bigwig, write_bedgraph, HAS_PYBIGWIG

    tsv_path = output_dir / 'corrected_reads.tsv'
    if not tsv_path.exists():
        logger.warning(
            f"corrected_reads.tsv not found in {output_dir} — skipping bedgraph generation"
        )
        return

    logger.info(f"Generating end-coverage tracks from {tsv_path}")
    df = pd.read_csv(tsv_path, sep='\t', index_col=False)
    logger.info(f"  Loaded {len(df):,} reads")

    chrom_sizes = {chrom: len(seq) for chrom, seq in genome.items()}

    ext = '.bw' if HAS_PYBIGWIG else '.bedgraph'
    write_fn = write_bigwig if HAS_PYBIGWIG else write_bedgraph
    fmt = 'bigWig' if HAS_PYBIGWIG else 'bedGraph'
    logger.info(f"  Output format: {fmt}")

    for end_type, col in [('3prime', 'corrected_3prime'), ('5prime', 'five_prime_position')]:
        if col not in df.columns:
            logger.warning(f"  Column {col!r} not in TSV — skipping {end_type} tracks")
            continue
        # Filter out sentinel / missing values (-1, NaN)
        sub = df[df[col].notna() & (df[col] >= 0)]
        if len(sub) == 0:
            logger.warning(f"  No valid positions for {col} — skipping {end_type} tracks")
            continue
        counts = aggregate_positions(sub, position_col=col)
        for strand, strand_label in [('+', 'plus'), ('-', 'minus')]:
            out_path = output_dir / f"{prefix}.{end_type}.{strand_label}{ext}"
            n = write_fn(counts, out_path, strand, chrom_sizes)
            logger.info(f"  {end_type} {strand_label}: {n:,} positions → {out_path.name}")


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
        for aligner in ('minimap2', 'mapPacBio', 'gapmm2', 'uLTRA', 'deSALT', 'bbmap', 'bwa'):
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

    read_num_sidecar = getattr(args, 'read_num_sidecar', None)
    if read_num_sidecar is None:
        read_num_sidecar = _autodetect_read_num_sidecar(bam_paths, prefix)
    elif not read_num_sidecar.exists():
        logger.error(f"Read-num sidecar not found: {read_num_sidecar}")
        return 1
    if read_num_sidecar is not None:
        logger.info(f"Read-num sidecar: {read_num_sidecar}")

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
            from ..consensus.consensus import load_annotated_junctions
            _t = _time.perf_counter()
            annotated_junctions = load_annotated_junctions(str(args.annotation))
            logger.info(f"[TIMING] Junction load: {_time.perf_counter() - _t:.1f}s")

    # Run consensus selection
    from ..consensus.consensus import run_consensus_selection
    output_bam = args.output_dir / f"{prefix}.consensus.bam"

    # Skip if consensus BAM already exists (safe resume after preemption)
    if output_bam.exists():
        logger.info(f"Consensus BAM already exists — skipping selection: {output_bam}")
        # Still generate stats/report/tracks if they're missing
        stats_tsv = args.output_dir / f"{prefix}.consensus_aligner_stats.tsv"
        if not stats_tsv.exists():
            logger.info("  (stats TSV missing — will regenerate from existing BAM)")
        else:
            if not getattr(args, 'no_bedgraph', False):
                try:
                    _generate_bedgraphs(args.output_dir, prefix, genome)
                except Exception as e:
                    logger.warning(f"Bedgraph generation failed (non-fatal): {e}")
            return 0

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
            read_num_sidecar=str(read_num_sidecar) if read_num_sidecar is not None else None,
        )
    except Exception as e:
        logger.error(f"Consensus selection failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

    logger.info(f"[TIMING] Consensus selection: {_time.perf_counter() - _t:.1f}s")
    logger.info(f"  High confidence: {stats.get('consensus_high', 'N/A')}")
    logger.info(f"  5' rescued: {stats.get('5prime_rescued', 'N/A')}")

    # Write aligner stats TSV and HTML report
    try:
        from ..bam.processing_stats import write_consensus_stats_tsv
        stats_tsv = args.output_dir / f"{prefix}.consensus_aligner_stats.tsv"
        write_consensus_stats_tsv(stats, str(stats_tsv))
    except Exception as _e:
        logger.warning(f"Could not write consensus stats TSV: {_e}")

    try:
        from ..analyze.summary import generate_consensus_html_report
        report_html = args.output_dir / f"{prefix}.consensus_report.html"
        generate_consensus_html_report(stats, str(report_html), sample_name=prefix)
        logger.info(f"Consensus report: {report_html}")
    except Exception as _e:
        logger.warning(f"Could not write consensus HTML report: {_e}")

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

    # Generate strand-specific 3' and 5' end coverage tracks from corrected_reads.tsv
    if not getattr(args, 'no_bedgraph', False):
        try:
            _generate_bedgraphs(args.output_dir, prefix, genome)
        except Exception as e:
            logger.warning(f"Bedgraph generation failed (non-fatal): {e}")

    return 0


def run(args: argparse.Namespace):
    """Entry point for CLI."""
    sys.exit(run_consensus(args))
