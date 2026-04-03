#!/usr/bin/env python3
"""
RECTIFY run-all command — complete pipeline from FASTQ/BAM to final analysis.

Chains every step with provenance tracking and smart step-skipping:

    align → correct → aggregate (3', 5', junctions) → export → analyze

Each step writes a .provenance.json sidecar recording input SHA-256 hashes,
the exact command/flags, and the RECTIFY version.  On re-run, any step whose
inputs AND flags are unchanged is skipped automatically.

Author: Kevin R. Roy
Date: 2026-03-28
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

from .. import __version__
from ..provenance import ProvenanceTracker, generate_header_lines

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_command_string(args: argparse.Namespace) -> str:
    """Reconstruct the command string from parsed args for provenance."""
    parts = ['rectify', 'run-all']
    if hasattr(args, 'bam') and args.bam:
        parts.append(str(args.bam))
    for attr in [
        'genome', 'annotation', 'output_dir', 'organism', 'netseq_dir',
        'aligner', 'manifest', 'go_annotations', 'reference',
    ]:
        val = getattr(args, attr, None)
        if val is not None:
            flag = f"--{attr.replace('_', '-')}"
            parts.append(f"{flag} {val}")
    if getattr(args, 'polya_sequenced', False):
        parts.append('--polya-sequenced')
    if getattr(args, 'force', False):
        parts.append('--force')
    parts.append(f"--threads {getattr(args, 'threads', 4)}")
    return ' '.join(parts)


def _load_genome(genome_path: Path) -> Dict[str, str]:
    """Load genome sequences into memory."""
    import pysam
    genome = {}
    fasta = pysam.FastaFile(str(genome_path))
    for chrom in fasta.references:
        genome[chrom] = fasta.fetch(chrom)
    fasta.close()
    return genome


def _write_tsv_with_header(
    df,
    output_path: Path,
    header_lines: List[str],
    sep: str = '\t',
) -> None:
    """Write a DataFrame to TSV with provenance header comment lines."""
    with open(output_path, 'w') as f:
        for line in header_lines:
            f.write(line + '\n')
        df.to_csv(f, sep=sep, index=False)


def _write_bedgraph_with_header(
    counts_df,
    output_path: Path,
    strand: str,
    header_lines: List[str],
) -> int:
    """Write strand-specific bedGraph with provenance header."""
    strand_df = counts_df[counts_df['strand'] == strand].copy()
    strand_df = strand_df.sort_values(['chrom', 'position'])

    with open(output_path, 'w') as f:
        for line in header_lines:
            f.write(line + '\n')
        for _, row in strand_df.iterrows():
            chrom = row['chrom']
            pos = int(row['position'])
            count = row['count']
            f.write(f"{chrom}\t{pos}\t{pos+1}\t{count}\n")

    return len(strand_df)


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

def step_align(args, tracker: ProvenanceTracker, output_dir: Path) -> Optional[Path]:
    """
    Step 1: Multi-aligner consensus alignment (FASTQ → BAM).

    Only runs if the input is FASTQ.  If input is already BAM, returns it as-is.
    """
    from .preprocess import detect_input_type

    input_path = Path(args.bam)
    input_type = detect_input_type(input_path)

    if input_type == 'bam':
        logger.info("[Step 1/6] Input is BAM — skipping alignment step")
        return input_path

    aligned_bam = output_dir / 'aligned' / 'consensus.bam'
    input_files = [input_path]
    if args.genome:
        input_files.append(Path(args.genome))

    if not args.force and not tracker.step_needed('align', input_files, [aligned_bam]):
        logger.info("[Step 1/6] Alignment is current — skipping")
        tracker.log_skip('align')
        return aligned_bam

    logger.info("[Step 1/6] Running multi-aligner consensus alignment...")
    aligned_dir = output_dir / 'aligned'
    aligned_dir.mkdir(parents=True, exist_ok=True)

    from .align_command import run_alignment
    align_args = argparse.Namespace(
        input=input_path,
        genome=Path(args.genome),
        annotation=Path(args.annotation) if args.annotation else None,
        output=aligned_bam,
        threads=getattr(args, 'threads', 4),
        aligner='consensus',
    )
    run_alignment(align_args)

    tracker.record_step('align', input_files, [aligned_bam])
    logger.info(f"  → {aligned_bam}")
    return aligned_bam


def step_correct(
    args, tracker: ProvenanceTracker, bam_path: Path, output_dir: Path
) -> Path:
    """Step 2: Correct 3' end positions."""
    from . import correct_command

    corrected_tsv = output_dir / 'corrected_3ends.tsv'
    input_files = [bam_path]
    if args.genome:
        input_files.append(Path(args.genome))

    if not args.force and not tracker.step_needed('correct', input_files, [corrected_tsv]):
        logger.info("[Step 2/6] Correction is current — skipping")
        tracker.log_skip('correct')
        return corrected_tsv

    logger.info("[Step 2/6] Correcting 3' end positions...")

    genome_path = Path(args.genome) if args.genome else None
    annotation_path = Path(args.annotation) if args.annotation else None
    organism = getattr(args, 'organism', None)

    # Resolve NET-seq
    netseq_path = None
    custom_netseq = getattr(args, 'netseq_dir', None)
    if custom_netseq:
        netseq_path = Path(custom_netseq)

    correct_args = argparse.Namespace(
        input=bam_path,
        genome=genome_path,
        annotation=annotation_path,
        output=corrected_tsv,
        netseq_dir=netseq_path,
        organism=organism,
        aligner=getattr(args, 'aligner', 'minimap2'),
        polya_sequenced=getattr(args, 'polya_sequenced', False),
        threads=getattr(args, 'threads', 4),
        min_mapq=10,
        skip_secondary=True,
        skip_supplementary=True,
        skip_ag_check=False,
        skip_atract_check=False,
        skip_polya_trim=False,
        skip_indel_correction=False,
        skip_variant_aware=False,
        polya_model=None,
        report=None,
        max_downstream_a=20,
        chunk_size=10000,
        debug=False,
        verbose=False,
        streaming=False,
    )

    correct_command.run(correct_args)

    # Prepend provenance header to the corrected TSV
    _prepend_header(corrected_tsv, tracker.header_lines('correct', input_files))

    tracker.record_step('correct', input_files, [corrected_tsv])
    logger.info(f"  → {corrected_tsv}")
    return corrected_tsv


def step_aggregate(
    args, tracker: ProvenanceTracker, bam_path: Path,
    corrected_tsv: Path, output_dir: Path, genome: Dict[str, str],
) -> Dict[str, Path]:
    """Step 3: Aggregate 3' ends, 5' ends, and junctions."""
    import pandas as pd

    agg_dir = output_dir / 'aggregated'
    agg_dir.mkdir(parents=True, exist_ok=True)

    outputs = {}
    annotation_path = Path(args.annotation) if args.annotation else None
    input_files_common = [bam_path]
    if annotation_path:
        input_files_common.append(annotation_path)

    # -- 3' end aggregation --
    three_prime_tsv = agg_dir / '3prime_clusters.tsv'
    if args.force or tracker.step_needed('aggregate_3prime', input_files_common, [three_prime_tsv]):
        logger.info("[Step 3/6] Aggregating 3' ends with gene attribution...")
        from .aggregate.three_prime import aggregate_3prime_clusters
        from .aggregate_command import load_annotation

        annotation_df = load_annotation(annotation_path) if annotation_path else pd.DataFrame(
            columns=['chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'feature_type']
        )
        clusters_df = aggregate_3prime_clusters(
            bam_path=str(bam_path),
            annotation_df=annotation_df,
            window=10,
            min_reads=1,
        )
        header = tracker.header_lines('aggregate_3prime', input_files_common)
        _write_tsv_with_header(clusters_df, three_prime_tsv, header)
        tracker.record_step('aggregate_3prime', input_files_common, [three_prime_tsv])
        logger.info(f"  → {three_prime_tsv}  ({len(clusters_df)} clusters)")
    else:
        logger.info("[Step 3/6] 3' end aggregation is current — skipping")
        tracker.log_skip('aggregate_3prime')
    outputs['3prime'] = three_prime_tsv

    # -- 5' end aggregation --
    five_prime_tsv = agg_dir / '5prime_clusters.tsv'
    if args.force or tracker.step_needed('aggregate_5prime', input_files_common, [five_prime_tsv]):
        logger.info("[Step 3/6] Aggregating 5' ends with gene attribution...")
        from .aggregate.five_prime import aggregate_5prime_clusters
        from .aggregate_command import load_annotation

        annotation_df = load_annotation(annotation_path) if annotation_path else pd.DataFrame(
            columns=['chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'feature_type']
        )
        clusters_df = aggregate_5prime_clusters(
            bam_path=str(bam_path),
            annotation_df=annotation_df,
            window=10,
            min_reads=1,
        )
        header = tracker.header_lines('aggregate_5prime', input_files_common)
        _write_tsv_with_header(clusters_df, five_prime_tsv, header)
        tracker.record_step('aggregate_5prime', input_files_common, [five_prime_tsv])
        logger.info(f"  → {five_prime_tsv}  ({len(clusters_df)} clusters)")
    else:
        logger.info("[Step 3/6] 5' end aggregation is current — skipping")
        tracker.log_skip('aggregate_5prime')
    outputs['5prime'] = five_prime_tsv

    # -- Junction aggregation --
    junctions_tsv = agg_dir / 'junctions.tsv'
    has_gff = annotation_path and annotation_path.suffix.lower() in ('.gff', '.gff3')

    if has_gff and genome:
        junction_inputs = input_files_common + ([Path(args.genome)] if args.genome else [])
        if args.force or tracker.step_needed('aggregate_junctions', junction_inputs, [junctions_tsv]):
            logger.info("[Step 3/6] Aggregating splice junctions with partial rescue...")
            from .aggregate.junctions import aggregate_junctions, merge_with_partial_evidence, export_junctions
            from .terminal_exon_refiner import load_splice_sites_from_gff, detect_partial_junction_crossings

            junction_df = aggregate_junctions(
                bam_path=str(bam_path), genome=genome, min_reads=1,
            )
            logger.info(f"  Found {len(junction_df)} junctions from CIGAR")

            # Rescue partial evidence
            try:
                splice_index = load_splice_sites_from_gff(str(annotation_path))
                partial_results = detect_partial_junction_crossings(
                    bam_path=str(bam_path), genome=genome,
                    splice_index=splice_index, min_clip_length=1,
                    ambiguous_mode='proportional',
                )
                n_rescued = partial_results['summary']['total_rescued']
                logger.info(f"  Rescued {n_rescued} partial crossings")
                junction_df = merge_with_partial_evidence(
                    junction_df, partial_results, ambiguous_mode='proportional',
                )
            except Exception as e:
                logger.warning(f"  Partial rescue failed: {e}")

            # Write with header
            header = tracker.header_lines('aggregate_junctions', junction_inputs)
            _write_tsv_with_header(junction_df, junctions_tsv, header)
            tracker.record_step('aggregate_junctions', junction_inputs, [junctions_tsv])
            logger.info(f"  → {junctions_tsv}  ({len(junction_df)} junctions)")
        else:
            logger.info("[Step 3/6] Junction aggregation is current — skipping")
            tracker.log_skip('aggregate_junctions')
    else:
        logger.info("[Step 3/6] Skipping junction aggregation (requires GFF annotation + genome)")
        tracker.log_skip('aggregate_junctions')

    outputs['junctions'] = junctions_tsv
    return outputs


def step_export(
    args, tracker: ProvenanceTracker, corrected_tsv: Path, output_dir: Path,
) -> Path:
    """Step 4: Export bedGraph/bigWig tracks (per-replicate and per-condition)."""
    import pandas as pd
    from .export_command import (
        load_corrected_3prime, aggregate_positions,
        process_per_replicate, process_per_condition,
    )

    export_dir = output_dir / 'tracks'
    export_dir.mkdir(parents=True, exist_ok=True)

    # Use a marker file to track whether export is current
    export_marker = export_dir / '.export_complete'
    input_files = [corrected_tsv]

    if not args.force and not tracker.step_needed('export', input_files, [export_marker]):
        logger.info("[Step 4/6] Export is current — skipping")
        tracker.log_skip('export')
        return export_dir

    logger.info("[Step 4/6] Exporting bedGraph/bigWig tracks...")

    # Determine chrom sizes
    chrom_sizes = {}
    if args.genome:
        try:
            from .export_command import get_chrom_sizes_from_genome
            chrom_sizes = get_chrom_sizes_from_genome(Path(args.genome))
        except Exception:
            pass

    if not chrom_sizes:
        # Fallback to bundled S. cerevisiae
        chrom_sizes = {
            'ref|NC_001133|': 230218, 'ref|NC_001134|': 813184,
            'ref|NC_001135|': 316620, 'ref|NC_001136|': 1531933,
            'ref|NC_001137|': 576874, 'ref|NC_001138|': 270161,
            'ref|NC_001139|': 1090940, 'ref|NC_001140|': 562643,
            'ref|NC_001141|': 439888, 'ref|NC_001142|': 745751,
            'ref|NC_001143|': 666816, 'ref|NC_001144|': 1078177,
            'ref|NC_001145|': 924431, 'ref|NC_001146|': 784333,
            'ref|NC_001147|': 1091291, 'ref|NC_001148|': 948066,
            'ref|NC_001224|': 85779,
        }

    # Load corrected data
    df = load_corrected_3prime(corrected_tsv, position_col='position')

    # Determine output format
    try:
        import pyBigWig
        fmt = 'bigwig'
    except ImportError:
        fmt = 'bedgraph'

    # Per-replicate
    rep_dir = export_dir / 'per_replicate'
    rep_dir.mkdir(exist_ok=True)
    process_per_replicate(df, rep_dir, fmt, 'position', chrom_sizes)

    # Per-condition
    cond_dir = export_dir / 'per_condition'
    cond_dir.mkdir(exist_ok=True)
    process_per_condition(df, cond_dir, fmt, 'position', chrom_sizes)

    # Write marker
    export_marker.write_text(f"completed {__version__}\n")

    # Collect all output files for provenance
    output_files = list(export_dir.rglob('*.bw')) + list(export_dir.rglob('*.bedgraph'))
    output_files.append(export_marker)
    tracker.record_step('export', input_files, output_files)

    logger.info(f"  → {export_dir}/")
    return export_dir


def step_analyze(
    args, tracker: ProvenanceTracker, corrected_tsv: Path, output_dir: Path,
) -> Path:
    """Step 5: Downstream analysis (clustering, DESeq2, PCA, motifs)."""
    from .analyze_command import run_analyze

    analysis_dir = output_dir  # analyze writes into output_dir directly
    analysis_marker = output_dir / 'tables' / '.analyze_complete'

    input_files = [corrected_tsv]
    if args.annotation:
        input_files.append(Path(args.annotation))

    if not args.force and not tracker.step_needed('analyze', input_files, [analysis_marker]):
        logger.info("[Step 5/6] Analysis is current — skipping")
        tracker.log_skip('analyze')
        return analysis_dir

    logger.info("[Step 5/6] Running downstream analysis...")

    annotation_path = Path(args.annotation) if args.annotation else None
    genome_path = Path(args.genome) if args.genome else None
    manifest_path = Path(args.manifest) if getattr(args, 'manifest', None) else None
    go_path = Path(args.go_annotations) if getattr(args, 'go_annotations', None) else None

    analyze_args = argparse.Namespace(
        input=corrected_tsv,
        annotation=annotation_path,
        output_dir=output_dir,
        genome=genome_path,
        reference=getattr(args, 'reference', None),
        manifest=manifest_path,
        go_annotations=go_path,
        threads=getattr(args, 'threads', 4),
        sample_column='sample',
        count_column=None,
        cluster_distance=25,
        min_reads=5,
        run_deseq2=True,
        run_motif=True,
        sample_sets=None,
        control=None,
        skip_browser_plots=False,
        max_browser_genes=10,
        top_genes=50,
    )

    try:
        exit_code = run_analyze(analyze_args)
        if exit_code != 0:
            logger.warning(f"Analysis completed with warnings (exit code: {exit_code})")
    except Exception as e:
        logger.error(f"Analysis failed: {e}")

    # Write marker
    tables_dir = output_dir / 'tables'
    tables_dir.mkdir(parents=True, exist_ok=True)
    analysis_marker.write_text(f"completed {__version__}\n")

    output_files = list(tables_dir.glob('*.tsv')) + [analysis_marker]
    tracker.record_step('analyze', input_files, output_files)

    logger.info(f"  → {output_dir}/")
    return analysis_dir


# ---------------------------------------------------------------------------
# Header prepend utility
# ---------------------------------------------------------------------------

def _prepend_header(filepath: Path, header_lines: List[str]) -> None:
    """Prepend provenance header lines to an existing file."""
    filepath = Path(filepath)
    if not filepath.exists():
        return

    original = filepath.read_text()

    # Don't double-prepend
    if original.startswith('# RECTIFY'):
        return

    with open(filepath, 'w') as f:
        for line in header_lines:
            f.write(line + '\n')
        f.write(original)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run(args: argparse.Namespace) -> None:
    """Execute the full run-all pipeline."""
    from ..slurm import set_thread_limits

    # Setup
    set_thread_limits(getattr(args, 'threads', None))
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    command_string = _build_command_string(args)
    tracker = ProvenanceTracker(output_dir, command_string)

    print("=" * 70)
    print(f"RECTIFY v{__version__}: Complete Pipeline (run-all)")
    print("=" * 70)
    print(f"Command: {command_string}")
    print(f"Output:  {output_dir}")
    if getattr(args, 'force', False):
        print("Mode:    FORCE (re-running all steps)")
    else:
        print("Mode:    Incremental (skipping unchanged steps)")
    print()

    # ------------------------------------------------------------------
    # Step 1: Align (if FASTQ input)
    # ------------------------------------------------------------------
    bam_path = step_align(args, tracker, output_dir)

    # ------------------------------------------------------------------
    # Step 2: Correct 3' end positions
    # ------------------------------------------------------------------
    corrected_tsv = step_correct(args, tracker, bam_path, output_dir)

    if not corrected_tsv.exists():
        logger.error(f"Corrected output not found: {corrected_tsv}")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 3: Aggregate (3' ends, 5' ends, junctions)
    # ------------------------------------------------------------------
    genome = {}
    if args.genome:
        try:
            logger.info("Loading genome for aggregation...")
            genome = _load_genome(Path(args.genome))
            logger.info(f"  Loaded {len(genome)} chromosomes")
        except Exception as e:
            logger.warning(f"Could not load genome: {e}")

    agg_outputs = step_aggregate(args, tracker, bam_path, corrected_tsv, output_dir, genome)

    # ------------------------------------------------------------------
    # Step 4: Export bedGraph/bigWig tracks
    # ------------------------------------------------------------------
    step_export(args, tracker, corrected_tsv, output_dir)

    # ------------------------------------------------------------------
    # Step 5: Analyze (clustering, DESeq2, etc.)
    # ------------------------------------------------------------------
    step_analyze(args, tracker, corrected_tsv, output_dir)

    # ------------------------------------------------------------------
    # Step 6: Write pipeline provenance manifest
    # ------------------------------------------------------------------
    logger.info("[Step 6/6] Writing pipeline provenance manifest...")
    manifest_path = tracker.write_manifest()
    logger.info(f"  → {manifest_path}")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    n_completed = sum(1 for s in tracker.steps if s['status'] == 'completed')
    n_skipped = len(tracker.skipped)

    print()
    print("=" * 70)
    print("Pipeline Complete!")
    print("=" * 70)
    print(f"  Steps run:     {n_completed}")
    print(f"  Steps skipped: {n_skipped}  (inputs unchanged)")
    print()
    print(f"Output directory: {output_dir}")
    print(f"Key outputs:")
    print(f"  Corrected 3' ends:  {corrected_tsv}")
    if agg_outputs.get('3prime') and agg_outputs['3prime'].exists():
        print(f"  3' end clusters:    {agg_outputs['3prime']}")
    if agg_outputs.get('5prime') and agg_outputs['5prime'].exists():
        print(f"  5' end clusters:    {agg_outputs['5prime']}")
    if agg_outputs.get('junctions') and agg_outputs['junctions'].exists():
        print(f"  Junctions:          {agg_outputs['junctions']}")
    print(f"  Tracks:             {output_dir / 'tracks'}/")
    print(f"  Analysis:           {output_dir / 'tables'}/")
    print(f"  Provenance:         {manifest_path}")
    print()


# ---------------------------------------------------------------------------
# CLI parser registration
# ---------------------------------------------------------------------------

def create_run_all_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register the run-all subcommand."""
    parser = subparsers.add_parser(
        'run-all',
        help='Complete pipeline: align → correct → aggregate → export → analyze',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Examples:
  # Full pipeline from BAM with bundled yeast data
  rectify run-all reads.bam --genome genome.fa --annotation genes.gff -o results/

  # From FASTQ (triggers multi-aligner consensus first)
  rectify run-all reads.fastq.gz --genome genome.fa --annotation genes.gff -o results/

  # Re-run only changed steps (default — checks input hashes)
  rectify run-all reads.bam --genome genome.fa --annotation genes.gff -o results/

  # Force re-run everything
  rectify run-all reads.bam --genome genome.fa --annotation genes.gff -o results/ --force
        """,
    )

    # Required
    parser.add_argument(
        'bam', type=Path,
        help='Input BAM or FASTQ file (FASTQ triggers alignment step)',
    )
    parser.add_argument(
        '--genome', type=Path, required=True,
        help='Reference genome FASTA file',
    )
    parser.add_argument(
        '--annotation', type=Path, required=True,
        help='Gene annotation file (GTF/GFF)',
    )
    parser.add_argument(
        '-o', '--output-dir', type=Path, required=True,
        help='Output directory for all results',
    )

    # Pipeline control
    ctrl = parser.add_argument_group('Pipeline control')
    ctrl.add_argument(
        '--force', action='store_true',
        help='Force re-run of all steps (ignore provenance cache)',
    )
    ctrl.add_argument(
        '--threads', type=int, default=4,
        help='Number of threads',
    )

    # Correction options
    corr = parser.add_argument_group('Correction options')
    corr.add_argument(
        '--organism', default='yeast',
        help='Organism name (for bundled NET-seq data)',
    )
    corr.add_argument(
        '--netseq-dir', type=Path,
        help='Custom NET-seq directory (overrides bundled data)',
    )
    corr.add_argument(
        '--aligner', choices=['minimap2', 'star', 'bowtie2', 'bwa'],
        default='minimap2', help='Aligner used for BAM file',
    )
    corr.add_argument(
        '--polya-sequenced', action='store_true',
        help='Poly(A) tail was sequenced',
    )

    # Analysis options
    anal = parser.add_argument_group('Analysis options')
    anal.add_argument(
        '--reference', help='Reference condition for DESeq2',
    )
    anal.add_argument(
        '--manifest', type=Path,
        help='Sample metadata TSV (sample, condition columns)',
    )
    anal.add_argument(
        '--go-annotations', type=Path,
        help='GO annotation file for enrichment analysis',
    )

    return parser
