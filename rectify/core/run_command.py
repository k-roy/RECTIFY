#!/usr/bin/env python3
"""
RECTIFY run-all command — complete end-to-end pipeline.

Handles everything from FASTQ (or BAM) to final results:

  Step 0  Alignment    triple-aligner consensus (minimap2 + mapPacBio + gapmm2)
                       Skipped automatically if {sample}.consensus.bam already exists,
                       or if input is already a BAM, or if --skip-alignment is set.

  Step 1  Correction   3' end correction (poly(A) trimming, indel correction,
                       NET-seq refinement, spike-in filtering)

  Step 2  Analysis     Clustering, DESeq2, GO enrichment, motif discovery
                       For multi-sample runs: runs after ALL samples are corrected
                       so DESeq2 has full statistical power across conditions.

Usage:
    # Single sample from FASTQ
    rectify run-all sample.fastq.gz --genome genome.fa --annotation genes.gff -o results/

    # Single sample from BAM (alignment skipped)
    rectify run-all sample.bam --genome genome.fa --annotation genes.gff -o results/

    # Multi-sample (manifest) — parallel correction + combined DESeq2
    rectify run-all --manifest manifest.tsv --genome genome.fa --annotation genes.gff -o results/

    # Multi-sample with SLURM
    rectify run-all --manifest manifest.tsv --genome genome.fa --annotation genes.gff -o results/ \\
        --profile my_cluster.yaml

Author: Kevin R. Roy
Date: 2026-03-28
"""

import argparse
import stat
import sys
import subprocess
import concurrent.futures
from pathlib import Path
from typing import Optional, List, Dict, Tuple


# ---------------------------------------------------------------------------
# Reference path resolution
# ---------------------------------------------------------------------------

def _resolve_reference_paths(args) -> None:
    """
    Resolve genome/annotation/GO annotation paths from explicit args or bundled data.

    Reads args.genome, args.annotation, args.organism.
    Updates args.genome, args.annotation, and args.go_annotations in place.

    Raises SystemExit if a genome cannot be resolved (nothing to align/correct against).
    """
    from rectify.data import ensure_reference_data, get_bundled_go_annotations_path, normalize_organism

    organism = getattr(args, 'organism', None)
    genome_arg = getattr(args, 'genome', None)
    annotation_arg = getattr(args, 'annotation', None)

    if organism is None and genome_arg is None:
        print(
            "ERROR: No reference genome provided.\n"
            "  Supply --genome /path/to/genome.fa,\n"
            "  or use --Scer (S. cerevisiae bundled data),\n"
            "  or use --organism <name> for another supported organism.",
            file=sys.stderr,
        )
        sys.exit(1)

    if organism is not None:
        genome_path, annotation_path, data_source = ensure_reference_data(
            organism=organism,
            custom_genome=genome_arg,
            custom_annotation=annotation_arg,
            verbose=True,
        )
    else:
        # Explicit paths only — validate they exist
        genome_path = Path(genome_arg) if genome_arg else None
        annotation_path = Path(annotation_arg) if annotation_arg else None
        data_source = 'custom' if genome_path else 'none'

    if genome_path is None:
        print(
            f"ERROR: No bundled genome available for organism '{organism}'. "
            "Use --genome to provide a custom reference.",
            file=sys.stderr,
        )
        sys.exit(1)

    args.genome = genome_path
    args.annotation = annotation_path

    # Resolve bundled GO annotations if not already set
    if not getattr(args, 'go_annotations', None) and organism:
        go_path = get_bundled_go_annotations_path(normalize_organism(organism))
        if go_path:
            args.go_annotations = go_path


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def _rectified_bam_path(sample_id: str, sample_output_dir: Path) -> Path:
    """Canonical path for the rectified BAM for a sample."""
    return sample_output_dir / f"{sample_id}.rectified.bam"


_ALIGNER_NAMES = ['minimap2', 'mapPacBio', 'gapmm2', 'uLTRA', 'deSALT']


def _collect_per_aligner_bams(
    sample_id: str,
    sample_output_dir: Path,
) -> Dict[str, Path]:
    """Return per-aligner BAM paths that exist on disk (keyed by aligner name)."""
    bams: Dict[str, Path] = {}
    for aligner in _ALIGNER_NAMES:
        bam = sample_output_dir / f"{sample_id}.{aligner}.bam"
        if bam.exists():
            bams[aligner] = bam
    return bams


def _run_alignment(
    input_path: Path,
    sample_id: str,
    sample_output_dir: Path,
    genome_path: Path,
    annotation_path: Optional[Path],
    threads: int,
    parallel_aligners: bool = False,
    junction_aligners: Optional[List[str]] = None,
    chimeric_consensus: bool = True,
    ultra_path: str = 'uLTRA',
    desalt_path: str = 'deSALT',
    mapPacBio_chunks: int = 1,
) -> Tuple[Dict[str, Path], Path]:
    """
    Run multi-aligner alignment and selection, or return existing rectified.bam.

    Default aligners: minimap2 + mapPacBio + gapmm2 + uLTRA + deSALT (all five),
    with chimeric consensus enabled by default. Pass junction_aligners=[] and
    chimeric_consensus=False to revert to single-best-aligner mode.

    Skips automatically if the rectified.bam already exists — safe to re-run.

    Returns
    -------
    Tuple of (per_aligner_bams, rectified_bam) where per_aligner_bams maps
    aligner name → BAM path for each per-aligner BAM found on disk.
    """
    rectified_bam = _rectified_bam_path(sample_id, sample_output_dir)

    # Backward compatibility: accept old consensus.bam filename from prior runs
    _legacy_bam = sample_output_dir / f"{sample_id}.consensus.bam"
    if not rectified_bam.exists() and _legacy_bam.exists():
        rectified_bam = _legacy_bam

    if rectified_bam.exists():
        print(f"    Skipping alignment — rectified.bam exists: {rectified_bam}")
        per_aligner_bams = _collect_per_aligner_bams(sample_id, sample_output_dir)
        return per_aligner_bams, rectified_bam

    _junction_aligners = junction_aligners or []
    all_aligners = ['minimap2', 'mapPacBio', 'gapmm2'] + _junction_aligners
    aligner_desc = ' + '.join(all_aligners)
    print(f"    Running {len(all_aligners)}-aligner consensus ({aligner_desc})...")
    from .align_command import run_align

    align_args = argparse.Namespace(
        reads=input_path,
        genome=genome_path,
        output_dir=sample_output_dir,
        annotation=annotation_path,
        threads=threads,
        aligners=['all'],
        junction_aligners=_junction_aligners,
        no_consensus=False,
        chimeric_consensus=chimeric_consensus,
        junc_bonus=9,
        junc_bed=None,
        parallel_aligners=parallel_aligners,
        minimap2_path='minimap2',
        mapPacBio_path='mapPacBio.sh',
        gapmm2_path='gapmm2',
        ultra_path=ultra_path,
        desalt_path=desalt_path,
        mapPacBio_chunks=mapPacBio_chunks,
        mapPacBio_chunk_idx=None,  # merge mode: look for existing chunk BAMs
        prefix='',
        keep_sam=False,
        sort=True,
        index=True,
        verbose=False,
    )

    # Temporarily hide scheduler array env vars so consensus.py runs in
    # single-sample mode. In run-all, array indices are for sample-level
    # parallelism; we never want within-sample read-level partitioning here.
    import os as _os
    _array_vars = {k: _os.environ.pop(k) for k in (
        # SLURM
        'SLURM_ARRAY_TASK_ID', 'SLURM_ARRAY_TASK_COUNT',
        'SLURM_ARRAY_TASK_MAX', 'SLURM_ARRAY_TASK_MIN',
        'SLURM_ARRAY_TASK_STEP',
        # UGE/SGE
        'SGE_TASK_ID', 'SGE_TASK_FIRST', 'SGE_TASK_LAST', 'SGE_TASK_STEPSIZE',
        # PBS/Torque
        'PBS_ARRAY_INDEX', 'PBS_ARRAYID',
    ) if k in _os.environ}
    try:
        rc = run_align(align_args)
    finally:
        _os.environ.update(_array_vars)
    if rc != 0:
        raise RuntimeError(f"Triple-aligner failed for {input_path}")

    if not rectified_bam.exists():
        raise RuntimeError(
            f"Alignment completed but rectified.bam not found: {rectified_bam}"
        )

    per_aligner_bams = _collect_per_aligner_bams(sample_id, sample_output_dir)
    return per_aligner_bams, rectified_bam


# ---------------------------------------------------------------------------
# Correction helper
# ---------------------------------------------------------------------------

def _bam_has_md_tags(bam_path: Path) -> bool:
    """Check if a BAM file's reads have MD tags (sample first 10 mapped reads)."""
    import logging
    logger = logging.getLogger(__name__)
    try:
        import pysam
        with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
            checked = 0
            for read in bam:
                if read.is_unmapped or read.is_secondary:
                    continue
                if read.has_tag('MD'):
                    return True
                checked += 1
                if checked >= 10:
                    break
    except FileNotFoundError:
        raise
    except Exception as e:
        logger.warning(f"Could not check MD tags: {e}")
    return False


def _run_correction(
    bam_path: Path,
    output_dir: Path,
    genome_path: Path,
    annotation_path: Optional[Path],
    args,
) -> Path:
    """
    Run rectify correct on a BAM, writing corrected_3ends.tsv into output_dir.
    Returns path to the corrected TSV.
    """
    from . import correct_command

    corrected_tsv = output_dir / 'corrected_3ends.tsv'

    # Indel correction and variant-aware rescue require MD tags in the BAM.
    # Disable them gracefully when MD tags are absent (e.g. rectified BAMs
    # generated by older runs without --MD).
    has_md = _bam_has_md_tags(bam_path)
    # --dT-primed-cDNA enables poly-A trimming / indel-correction modules.
    # Old --no-polya-sequenced is a deprecated alias for NOT passing --dT-primed-cDNA.
    # Default (no flag): direct RNA mode — poly-A IS in the read.
    polya_seq = getattr(args, 'dT_primed_cDNA', False)
    skip_indel = not has_md  # Can't correct indels without MD tags
    skip_variant = not has_md

    if not has_md and polya_seq:
        print(
            "    Note: BAM lacks MD tags — indel correction disabled. "
            "Re-run alignment with current RECTIFY to enable it."
        )

    # Auto-generate report path alongside the corrected TSV
    report_path = corrected_tsv.parent / (corrected_tsv.stem + '_report.html')

    # Derive sample prefix from input BAM stem (strip .rectified / .consensus suffixes)
    _stem = bam_path.stem
    for _sfx in ('.rectified', '.consensus'):
        if _stem.endswith(_sfx):
            _stem = _stem[:-len(_sfx)]
            break
    corrected_bam_path   = output_dir / f"{_stem}.rectified_corrected_3end.bam"
    softclipped_bam_path = output_dir / f"{_stem}.rectified_pA_tail_trimmed.bam"

    correct_args = argparse.Namespace(
        input=bam_path,
        genome=genome_path,
        annotation=annotation_path,
        output=corrected_tsv,
        corrected_bam=corrected_bam_path,
        softclipped_bam=softclipped_bam_path,
        netseq_dir=getattr(args, 'netseq_dir', None),
        organism=getattr(args, 'organism', None),
        aligner=getattr(args, 'aligner', 'minimap2'),
        dT_primed_cDNA=polya_seq,
        polya_sequenced=polya_seq,  # deprecated attr kept for compat
        threads=getattr(args, 'threads', 4),
        filter_spikein=getattr(args, 'filter_spikein', None),
        streaming=getattr(args, 'streaming', False),
        # Defaults
        min_mapq=10,
        skip_secondary=True,
        skip_supplementary=True,
        skip_ag_check=False,
        skip_atract_check=False,
        skip_polya_trim=False,
        skip_indel_correction=skip_indel,
        skip_variant_aware=skip_variant,
        polya_model=None,
        report=report_path,
        max_downstream_a=20,
        chunk_size=getattr(args, 'chunk_size', 10000),
        debug=False,
        verbose=False,
    )

    correct_command.run(correct_args)

    if not corrected_tsv.exists():
        raise RuntimeError(f"Correction did not produce output: {corrected_tsv}")

    return corrected_tsv


def _run_correction_per_aligner(
    per_aligner_bams: Dict[str, Path],
    output_dir: Path,
    genome_path: Path,
    annotation_path: Optional[Path],
    args,
) -> Dict[str, Path]:
    """
    Run ``rectify correct`` on each per-aligner BAM independently.

    Writes per-aligner outputs to ``output_dir/per_aligner_corrected/{aligner}/``.
    Skips any aligner whose TSV already exists (safe to resume).

    Returns
    -------
    Dict mapping aligner name → corrected_3ends.tsv path for each aligner
    whose correction succeeded.
    """
    per_aligner_dir = output_dir / 'per_aligner_corrected'
    per_aligner_dir.mkdir(parents=True, exist_ok=True)

    per_aligner_tsvs: Dict[str, Path] = {}
    for aligner_name, bam_path in per_aligner_bams.items():
        aligner_output_dir = per_aligner_dir / aligner_name
        aligner_output_dir.mkdir(exist_ok=True)
        tsv_path = aligner_output_dir / 'corrected_3ends.tsv'

        if tsv_path.exists():
            print(f"    [{aligner_name}] Skipping — TSV exists: {tsv_path}")
            per_aligner_tsvs[aligner_name] = tsv_path
            continue

        print(f"    [{aligner_name}] Correcting {bam_path.name}...")
        try:
            _run_correction(
                bam_path=bam_path,
                output_dir=aligner_output_dir,
                genome_path=genome_path,
                annotation_path=annotation_path,
                args=args,
            )
        except Exception as exc:
            print(
                f"    WARNING: [{aligner_name}] correction failed: {exc}",
                file=sys.stderr,
            )
            continue

        if tsv_path.exists():
            per_aligner_tsvs[aligner_name] = tsv_path
        else:
            print(
                f"    WARNING: [{aligner_name}] correction produced no output",
                file=sys.stderr,
            )

    return per_aligner_tsvs


# ---------------------------------------------------------------------------
# Combine corrected TSVs
# ---------------------------------------------------------------------------

def _combine_corrected_tsvs(
    samples: List[Dict[str, str]],
    output_dir: Path,
) -> Path:
    """
    Concatenate per-sample corrected TSVs into a single multi-sample TSV.

    Adds a `sample` column equal to sample_id so the analyze command can
    distinguish samples for DESeq2.
    """
    import pandas as pd

    combined_dir = output_dir / 'combined'
    combined_dir.mkdir(parents=True, exist_ok=True)
    combined_tsv = combined_dir / 'corrected_3ends_combined.tsv'

    dfs = []
    missing = []
    for sample in samples:
        tsv_path = output_dir / sample['sample_id'] / 'corrected_3ends.tsv'
        if not tsv_path.exists():
            missing.append(sample['sample_id'])
            print(f"  WARNING: {tsv_path} not found, skipping", file=sys.stderr)
            continue
        df = pd.read_csv(tsv_path, sep='\t')
        df['sample'] = sample['sample_id']
        dfs.append(df)
        print(f"  {sample['sample_id']}: {len(df):,} reads")

    if not dfs:
        raise RuntimeError(
            "No corrected TSVs found to combine. "
            "Ensure the correction step completed successfully."
        )

    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(combined_tsv, sep='\t', index=False)
    print(f"  Combined → {combined_tsv}  ({len(dfs)} samples, {len(combined):,} reads)")

    if missing:
        print(
            f"  WARNING: {len(missing)} samples missing corrected TSV: {', '.join(missing)}",
            file=sys.stderr,
        )

    return combined_tsv


# ---------------------------------------------------------------------------
# Analysis helper
# ---------------------------------------------------------------------------

def _run_analysis(
    corrected_tsv: Path,
    output_dir: Path,
    genome_path: Optional[Path],
    annotation_path: Optional[Path],
    args,
    n_samples: int = 1,
) -> None:
    """Run the analyze command on a (possibly combined) corrected TSV."""
    from .analyze_command import run_analyze

    # Only run DESeq2 when there are multiple samples
    run_deseq2 = n_samples > 1

    analyze_args = argparse.Namespace(
        input=corrected_tsv,
        output=output_dir,
        annotation=annotation_path,
        genome=genome_path,
        reference=getattr(args, 'reference', None),
        go_annotations=getattr(args, 'go_annotations', None),
        threads=getattr(args, 'threads', 4),
        # Clustering
        sample_column='sample',
        count_column=None,
        cluster_distance=25,
        min_reads=5,
        # Analysis flags
        run_deseq2=run_deseq2,
        run_motif=run_deseq2,   # motif discovery only makes sense with DESeq2 results
        sample_sets=None,
        # Filtering
        exclude_mito=True,
        include_mito=False,
        exclude_rdna=True,
        include_rdna=False,
        # Output options
        no_bedgraph=False,
        bedgraph_dir=None,
        no_genomic_distribution=False,
        # Manifest mode (not used in single-file path)
        manifest=None,
        # Motif windows
        motif_upstream=100,
        motif_downstream=50,
    )

    exit_code = run_analyze(analyze_args)
    if exit_code != 0:
        print(f"\nAnalysis completed with warnings (exit code: {exit_code})")


def _run_analysis_manifest(
    manifest_path: Path,
    output_dir: Path,
    genome_path: Optional[Path],
    annotation_path: Optional[Path],
    args,
    n_samples: int = 1,
) -> None:
    """Run the analyze command in manifest mode (memory-efficient multi-sample path)."""
    from .analyze_command import run_analyze

    run_deseq2 = n_samples > 1

    analyze_args = argparse.Namespace(
        # In manifest mode, 'input' is unused but run_analyze checks it after
        # the manifest dispatch, so set it to the manifest path as a fallback.
        input=str(manifest_path),
        manifest=str(manifest_path),
        output=output_dir,
        annotation=annotation_path,
        genome=genome_path,
        reference=getattr(args, 'reference', None),
        go_annotations=getattr(args, 'go_annotations', None),
        threads=getattr(args, 'threads', 4),
        # Clustering
        sample_column='sample',
        count_column=None,
        cluster_distance=25,
        min_reads=5,
        # Analysis flags
        run_deseq2=run_deseq2,
        run_motif=run_deseq2,
        sample_sets=None,
        # Filtering
        exclude_mito=True,
        include_mito=False,
        exclude_rdna=True,
        include_rdna=False,
        # Manifest mode disables bedgraph/genomic distribution (no combined positions_df)
        no_bedgraph=True,
        bedgraph_dir=None,
        no_genomic_distribution=True,
        # Motif windows
        motif_upstream=100,
        motif_downstream=50,
    )

    exit_code = run_analyze(analyze_args)
    if exit_code != 0:
        print(f"\nAnalysis completed with warnings (exit code: {exit_code})")


# ---------------------------------------------------------------------------
# Junction aggregation helper
# ---------------------------------------------------------------------------

def _run_junction_aggregation(
    bam_path: Path,
    genome_path: Path,
    annotation_path: Path,
    output_dir: Path,
    config: Optional[dict] = None,
) -> Optional[Path]:
    """
    Run junction aggregation with partial rescue.
    Only runs if annotation is GFF/GFF3 (has intron features).
    Returns path to junctions TSV, or None if skipped/failed.
    """
    # Handle compressed files: .gff.gz → stem is .gff, suffixes[-1] is .gz
    suffixes = [s.lower() for s in annotation_path.suffixes]
    has_gff = any(s in ('.gff', '.gff3') for s in suffixes)
    if not has_gff:
        print(
            "\n[Junctions] Skipping — requires GFF/GFF3 annotation "
            "(current annotation has no intron features)"
        )
        return None

    print("\n[Junctions] Aggregating splice junctions with partial rescue...")

    from .aggregate.junctions import aggregate_junctions, merge_with_partial_evidence, export_junctions
    from .terminal_exon_refiner import load_splice_sites_from_gff, detect_partial_junction_crossings
    import pysam

    junctions_dir = output_dir / 'junctions'
    junctions_dir.mkdir(parents=True, exist_ok=True)
    junctions_tsv = junctions_dir / 'junctions.tsv'

    try:
        genome = {}
        fasta = pysam.FastaFile(str(genome_path))
        for chrom in fasta.references:
            genome[chrom] = fasta.fetch(chrom)
        fasta.close()

        junction_df = aggregate_junctions(
            bam_path=str(bam_path), genome=genome, min_reads=1,
        )
        print(f"  Found {len(junction_df)} junctions from CIGAR")

        splice_index = load_splice_sites_from_gff(str(annotation_path))
        partial_results = detect_partial_junction_crossings(
            bam_path=str(bam_path),
            genome=genome,
            splice_index=splice_index,
            min_clip_length=1,
            ambiguous_mode='proportional',
        )

        n_rescued = partial_results['stats']['reads_rescued_as_spliced']
        n_ambiguous = partial_results['stats']['reads_ambiguous']
        print(f"  Rescued {n_rescued} partial crossings ({n_ambiguous} ambiguous)")

        junction_df = merge_with_partial_evidence(
            junction_df, partial_results, ambiguous_mode='proportional',
        )
        export_junctions(junction_df, str(junctions_tsv), format='tsv')
        print(f"  Junctions written to {junctions_tsv}")

        # Provenance
        if junctions_tsv.exists():
            from ..utils.provenance import init_provenance
            prov = init_provenance(output_dir, description="RECTIFY junction aggregation", config=config)
            prov.add_output_file(
                junctions_tsv,
                source_files=[bam_path],
                metadata={
                    'n_junctions': len(junction_df),
                    'n_rescued': n_rescued,
                    'n_ambiguous': n_ambiguous,
                },
            )
            prov.save()

        return junctions_tsv

    except Exception as e:
        print(f"\nWarning: Junction aggregation failed: {e}", file=sys.stderr)
        print("Continuing without junction output...")
        return None


# ---------------------------------------------------------------------------
# Per-sample worker (used by ThreadPoolExecutor in multi-sample mode)
# ---------------------------------------------------------------------------

def _process_one_sample(
    sample: Dict[str, str],
    output_dir: Path,
    genome_path: Path,
    annotation_path: Path,
    args,
) -> Tuple[str, int]:
    """
    Process one sample: alignment (if needed) → correction.
    Returns (sample_id, returncode).  Called from threads.
    """
    import re
    import logging
    _logger = logging.getLogger(__name__)
    sample_id = sample['sample_id']
    safe_id = re.sub(r'[^a-zA-Z0-9_\-.]', '_', sample_id)
    if safe_id != sample_id:
        _logger.warning(f"sample_id '{sample_id}' contains invalid characters, using '{safe_id}'")
        sample_id = safe_id
    input_path = Path(sample['path'] if 'path' in sample else sample['bam_path'])
    sample_output = output_dir / sample_id
    sample_output.mkdir(parents=True, exist_ok=True)

    log_file = sample_output / 'rectify_run.log'

    try:
        with open(log_file, 'w') as log:
            # Determine input type
            from .preprocess import detect_input_type
            input_type = detect_input_type(input_path)

            bam_to_correct = input_path

            if input_type in ('fastq', 'fastq.gz') and not getattr(args, 'skip_alignment', False):
                print(f"  [{sample_id}] Aligning…", flush=True)
                # --bam-dir: per-sample subdir within the requested bam_dir
                _bam_dir_arg = getattr(args, 'bam_dir', None)
                if _bam_dir_arg:
                    _align_out = Path(_bam_dir_arg) / sample_id
                    _align_out.mkdir(parents=True, exist_ok=True)
                else:
                    _align_out = sample_output
                try:
                    bam_to_correct = _run_alignment(
                        input_path=input_path,
                        sample_id=sample_id,
                        sample_output_dir=_align_out,
                        genome_path=genome_path,
                        annotation_path=annotation_path,
                        threads=getattr(args, 'threads', 4),
                        parallel_aligners=getattr(args, 'parallel_aligners', False),
                        junction_aligners=getattr(args, 'junction_aligners', []),
                        chimeric_consensus=getattr(args, 'chimeric_consensus', True),
                        ultra_path=getattr(args, 'ultra_path', 'uLTRA'),
                        desalt_path=getattr(args, 'desalt_path', 'deSALT'),
                        mapPacBio_chunks=getattr(args, 'mapPacBio_chunks', 1),
                    )
                    log.write(f"Alignment complete: {bam_to_correct}\n")
                except Exception as e:
                    log.write(f"Alignment failed: {e}\n")
                    return sample_id, 1
            elif input_type in ('fastq', 'fastq.gz'):
                # FASTQ but alignment skipped — look for existing rectified.bam
                rectified_bam = _rectified_bam_path(sample_id, sample_output)
                _legacy_bam = sample_output / f"{sample_id}.consensus.bam"
                if not rectified_bam.exists() and _legacy_bam.exists():
                    rectified_bam = _legacy_bam
                if rectified_bam.exists():
                    bam_to_correct = rectified_bam
                    log.write(f"Using existing rectified.bam: {rectified_bam}\n")
                else:
                    log.write(
                        f"ERROR: --skip-alignment set but no rectified.bam found: {rectified_bam}\n"
                    )
                    return sample_id, 1

            print(f"  [{sample_id}] Correcting 3' ends…", flush=True)
            try:
                _run_correction(
                    bam_path=bam_to_correct,
                    output_dir=sample_output,
                    genome_path=genome_path,
                    annotation_path=annotation_path,
                    args=args,
                )
                log.write("Correction complete\n")
            except Exception as e:
                log.write(f"Correction failed: {e}\n")
                return sample_id, 1

        return sample_id, 0

    except Exception as e:
        with open(log_file, 'a') as log:
            log.write(f"Unexpected error: {e}\n")
        return sample_id, 1


# ---------------------------------------------------------------------------
# Multi-sample pipeline
# ---------------------------------------------------------------------------

def _run_multi_sample(args) -> int:
    """
    Multi-sample pipeline via manifest:
      Stage 1 (parallel): align + correct per sample
      Stage 2:            combine corrected TSVs → add sample column
      Stage 3:            combined analyze (full DESeq2, GO, motifs)
    """
    from .batch_command import parse_manifest, _get_available_cpus

    _resolve_reference_paths(args)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_path = args.genome
    annotation_path = args.annotation

    samples = parse_manifest(Path(args.manifest))
    if not samples:
        print("ERROR: No samples found in manifest.", file=sys.stderr)
        return 1

    # Validate 'path' column exists (manifest uses 'bam_path' or 'path')
    for s in samples:
        if 'path' not in s and 'bam_path' not in s:
            print(
                f"ERROR: manifest row for '{s['sample_id']}' has no path column.",
                file=sys.stderr,
            )
            return 1

    print(f"Found {len(samples)} samples:")
    for s in samples:
        p = s.get('path', s.get('bam_path', '?'))
        cond = f"  [{s['condition']}]" if 'condition' in s else ''
        print(f"  {s['sample_id']}: {p}{cond}")
    print()

    # ── Stage 1: align + correct in parallel ────────────────────────────────
    n_cpus = _get_available_cpus()
    threads_per_sample = getattr(args, 'threads', 4)
    n_workers = max(1, n_cpus // threads_per_sample)

    print(f"[Stage 1/3] Aligning and correcting ({n_workers} parallel workers)...")

    failed = []
    completed = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        future_to_sample = {
            executor.submit(
                _process_one_sample, s, output_dir, genome_path, annotation_path, args
            ): s
            for s in samples
        }

        for future in concurrent.futures.as_completed(future_to_sample):
            s = future_to_sample[future]
            completed += 1
            try:
                sample_id, rc = future.result()
                status = "OK  " if rc == 0 else "FAIL"
                n = len(str(len(samples)))
                print(f"  [{completed:>{n}}/{len(samples)}] [{status}] {sample_id}")
                if rc != 0:
                    failed.append(sample_id)
                    if not getattr(args, 'continue_on_error', False):
                        print(
                            f"\nERROR: {sample_id} failed. "
                            f"Use --continue-on-error to proceed anyway.",
                            file=sys.stderr,
                        )
                        return 1
            except Exception as e:
                completed += 1
                print(f"  [ERROR] {s['sample_id']}: {e}", file=sys.stderr)
                failed.append(s['sample_id'])
                if not getattr(args, 'continue_on_error', False):
                    return 1

    if failed:
        print(f"\n{len(failed)} sample(s) failed: {', '.join(failed)}", file=sys.stderr)
        if len(failed) == len(samples):
            return 1

    # Only combine successfully corrected samples
    successful_samples = [s for s in samples if s['sample_id'] not in failed]

    # ── Stage 2: write manifest pointing to per-sample corrected TSVs ────────
    # No combine step needed — manifest mode in analyze streams each file separately
    print(f"\n[Stage 2/3] Writing sample manifest for combined analysis...")
    import pandas as pd
    combined_dir = output_dir / 'combined'
    combined_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    for s in successful_samples:
        row = {
            'sample_id': s['sample_id'],
            'path': str(output_dir / s['sample_id'] / 'corrected_3ends.tsv'),
        }
        if 'condition' in s:
            row['condition'] = s.get('condition', '')
        manifest_rows.append(row)

    manifest_for_analyze = combined_dir / 'corrected_manifest.tsv'
    pd.DataFrame(manifest_rows).to_csv(manifest_for_analyze, sep='\t', index=False)
    print(f"  Wrote manifest: {manifest_for_analyze} ({len(manifest_rows)} samples)")

    # ── Stage 3: combined analysis (manifest mode) ───────────────────────────
    print(f"\n[Stage 3/3] Running combined analysis (DESeq2, GO, motifs)...")
    try:
        _run_analysis_manifest(
            manifest_path=manifest_for_analyze,
            output_dir=combined_dir,
            genome_path=genome_path,
            annotation_path=annotation_path,
            args=args,
            n_samples=len(successful_samples),
        )
    except Exception as e:
        print(f"ERROR: Combined analysis failed: {e}", file=sys.stderr)
        return 1

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("Pipeline Complete!")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir}")
    print("\nPer-sample outputs:")
    for s in successful_samples:
        print(f"  {output_dir}/{s['sample_id']}/corrected_3ends.tsv")
    print(f"\nCombined analysis: {combined_dir}/")
    print(f"  DESeq2 results:  {combined_dir}/tables/deseq2_genes_*.tsv")
    print(f"  HTML report:     {combined_dir}/report.html")
    print(f"  Provenance:      {combined_dir}/PROVENANCE.json")

    return 0


# ---------------------------------------------------------------------------
# Single-sample pipeline
# ---------------------------------------------------------------------------

def _run_single_sample(args) -> int:
    """
    Single-sample pipeline:
      Step 0 (if FASTQ): multi-aligner alignment → rectified.bam
      Step 1: correction → corrected_3ends.tsv
      Step 2: analysis (no DESeq2 — single sample)
      Step 3 (if GFF): junction aggregation
    """
    _resolve_reference_paths(args)

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_path = args.genome
    annotation_path = args.annotation

    # Pre-flight: verify required input files exist before starting any work
    for path, label in [
        (input_path, 'input'),
        (genome_path, 'genome'),
        (annotation_path, 'annotation GFF'),
    ]:
        if path and not Path(path).exists():
            raise FileNotFoundError(f"{label} file not found: {path}")

    print("=" * 70)
    print("RECTIFY: Single-Sample Pipeline")
    print("=" * 70)
    print(f"\nInput:      {input_path}")
    print(f"Output dir: {output_dir}")

    # Determine input type
    from .preprocess import detect_input_type
    input_type = detect_input_type(input_path)
    print(f"Input type: {input_type}")

    bam_to_correct = input_path
    per_aligner_bams: Dict[str, Path] = {}

    # ── Scratch staging setup ────────────────────────────────────────────────
    # If $SCRATCH is available, run all I/O there and rsync back.
    # This avoids NFS contention across concurrent array tasks.
    # The rectified BAM is always copied back to the output dir even on fresh
    # alignments so it survives $SCRATCH's auto-purge and enables job resumption.
    from ..slurm import make_job_scratch_dir, sync_to_oak
    import shutil as _shutil

    use_scratch = getattr(args, 'use_scratch', True)  # default on when scratch available
    scratch_dir = make_job_scratch_dir('rectify_single') if use_scratch else None
    work_dir = scratch_dir if scratch_dir else output_dir

    # --bam-dir: explicit persistent directory for alignment BAMs.
    # --keep-aligner-bams: sync per-aligner BAMs to Oak (default: discard them).
    _bam_dir_arg = getattr(args, 'bam_dir', None)
    bam_dir: Optional[Path] = Path(_bam_dir_arg) if _bam_dir_arg else None
    if bam_dir:
        bam_dir.mkdir(parents=True, exist_ok=True)
    keep_aligner_bams: bool = getattr(args, 'keep_aligner_bams', False)

    if scratch_dir:
        print(f"\nScratch staging enabled: {scratch_dir}")

    # ── Provenance tracker ───────────────────────────────────────────────────
    # scratch_root/oak_root mapping ensures all sidecar JSON keys use canonical
    # Oak paths — so step_is_current() works on re-runs regardless of $SCRATCH.
    from ..provenance import ProvenanceTracker
    _cmd_string = ' '.join(sys.argv)
    tracker = ProvenanceTracker(
        output_dir=output_dir,
        command_string=_cmd_string,
        scratch_root=scratch_dir,
        oak_root=output_dir if scratch_dir else None,
    )

    # ── Step 0: Alignment ────────────────────────────────────────────────────
    import time as _time
    _pipeline_start = _time.perf_counter()
    if input_type in ('fastq', 'fastq.gz') and not getattr(args, 'skip_alignment', False):
        sample_id = input_path.stem.replace('.fastq', '').replace('.gz', '')
        print(f"\n[Step 1/3] Aligning with triple-aligner...")
        print("-" * 50)
        _t0 = _time.perf_counter()
        # Align to bam_dir if specified; otherwise use work_dir (scratch if available)
        align_output_dir = bam_dir if bam_dir else work_dir
        per_aligner_bams, bam_to_correct = _run_alignment(
            input_path=input_path,
            sample_id=sample_id,
            sample_output_dir=align_output_dir,
            genome_path=genome_path,
            annotation_path=annotation_path,
            threads=getattr(args, 'threads', 4),
            parallel_aligners=getattr(args, 'parallel_aligners', False),
            junction_aligners=getattr(args, 'junction_aligners', []),
            chimeric_consensus=getattr(args, 'chimeric_consensus', True),
            ultra_path=getattr(args, 'ultra_path', 'uLTRA'),
            desalt_path=getattr(args, 'desalt_path', 'deSALT'),
            mapPacBio_chunks=getattr(args, 'mapPacBio_chunks', 1),
        )
        print(f"\nAlignment complete: {bam_to_correct}")
        print(f"[TIMING] Alignment: {_time.perf_counter() - _t0:.1f}s")
        tracker.record_step('align', input_files=[input_path], output_files=[bam_to_correct])
        step_correction = 2
    elif input_type in ('fastq', 'fastq.gz'):
        sample_id = input_path.stem.replace('.fastq', '').replace('.gz', '')
        rectified_bam = _rectified_bam_path(sample_id, output_dir)
        _legacy_bam = output_dir / f"{sample_id}.consensus.bam"
        if not rectified_bam.exists() and _legacy_bam.exists():
            rectified_bam = _legacy_bam
        if rectified_bam.exists():
            print(f"\n[Step 1/3] Alignment — using existing rectified.bam")
            if scratch_dir:
                # Stage existing BAM to scratch
                print(f"  Staging to scratch: {scratch_dir}")
                _shutil.copy2(rectified_bam, scratch_dir / rectified_bam.name)
                bai = Path(str(rectified_bam) + '.bai')
                if bai.exists():
                    _shutil.copy2(bai, scratch_dir / bai.name)
                bam_to_correct = scratch_dir / rectified_bam.name
                tracker.register_staged(bam_to_correct, rectified_bam)
            else:
                bam_to_correct = rectified_bam
        else:
            print(
                f"ERROR: --skip-alignment set but no rectified.bam found: {rectified_bam}",
                file=sys.stderr,
            )
            if scratch_dir:
                _shutil.rmtree(scratch_dir, ignore_errors=True)
            return 1
        step_correction = 2
    else:
        print(f"\nBAM input detected — skipping alignment")
        if scratch_dir:
            # Stage BAM to scratch
            print(f"  Staging to scratch: {scratch_dir}")
            _shutil.copy2(input_path, scratch_dir / input_path.name)
            bai = Path(str(input_path) + '.bai')
            if bai.exists():
                _shutil.copy2(bai, scratch_dir / bai.name)
            bam_to_correct = scratch_dir / input_path.name
            tracker.register_staged(bam_to_correct, input_path)
        step_correction = 1

    # ── Step 1/2: Correction ─────────────────────────────────────────────────
    print(f"\n[Step {step_correction}/3] Correcting 3' end positions...")
    print("-" * 50)
    _t0 = _time.perf_counter()
    try:
        if per_aligner_bams:
            # New workflow: correct each aligner's BAM independently, then merge.
            # The final corrected_3ends.tsv is selected from post-correction features
            # (five_prime_rescued, confidence, 3' agreement) rather than raw alignment
            # features — which are not cross-comparable across aligners.
            from .corrected_consensus import merge_corrected_tsvs, identify_cat5_candidates
            print(f"    Running per-aligner correction ({len(per_aligner_bams)} aligners)...")
            per_aligner_tsvs = _run_correction_per_aligner(
                per_aligner_bams=per_aligner_bams,
                output_dir=work_dir,
                genome_path=genome_path,
                annotation_path=annotation_path,
                args=args,
            )
            if per_aligner_tsvs:
                _per_aligner_dir = work_dir / 'per_aligner_corrected'
                _summary_tsv = _per_aligner_dir / 'comparison_summary.tsv'
                print(f"    Merging {len(per_aligner_tsvs)} per-aligner TSVs...")
                corrected_tsv = merge_corrected_tsvs(
                    per_aligner_tsvs=per_aligner_tsvs,
                    output_tsv=work_dir / 'corrected_3ends.tsv',
                    summary_tsv=_summary_tsv,
                )
                # Identify Cat5 candidates (reads where aligners contribute unique introns)
                _cat5_tsv = _per_aligner_dir / 'cat5_candidates.tsv'
                identify_cat5_candidates(per_aligner_tsvs, output_tsv=_cat5_tsv)
            else:
                # All per-aligner corrections failed — fall back to consensus BAM
                print(
                    "    WARNING: No per-aligner correction succeeded; "
                    "falling back to consensus BAM correction.",
                    file=sys.stderr,
                )
                corrected_tsv = _run_correction(
                    bam_path=bam_to_correct,
                    output_dir=work_dir,
                    genome_path=genome_path,
                    annotation_path=annotation_path,
                    args=args,
                )
        else:
            # No per-aligner BAMs available (BAM input, or aligner BAMs were discarded).
            # Use existing single-BAM correction path.
            corrected_tsv = _run_correction(
                bam_path=bam_to_correct,
                output_dir=work_dir,
                genome_path=genome_path,
                annotation_path=annotation_path,
                args=args,
            )
    except Exception as e:
        print(f"ERROR in correction step: {e}", file=sys.stderr)
        if scratch_dir:
            _shutil.rmtree(scratch_dir, ignore_errors=True)
        return 1
    print(f"\nCorrection complete: {corrected_tsv}")
    print(f"[TIMING] Correction: {_time.perf_counter() - _t0:.1f}s")
    tracker.record_step('correct', input_files=[bam_to_correct], output_files=[corrected_tsv])

    # ── Early partial sync: corrected TSV → Oak ──────────────────────────────
    # After correction, provenance moves back to Oak. Sync the corrected TSV
    # (and its sidecar) immediately so analysis reads from Oak, not ephemeral
    # scratch — matching the canonical paths written into the sidecar.
    if scratch_dir:
        _oak_corrected_tsv = output_dir / corrected_tsv.name
        _shutil.copy2(corrected_tsv, _oak_corrected_tsv)
        _sidecar = corrected_tsv.parent / f".{corrected_tsv.name}.provenance.json"
        if _sidecar.exists():
            _shutil.copy2(_sidecar, output_dir / _sidecar.name)
        corrected_tsv = _oak_corrected_tsv
        work_dir = output_dir  # analysis outputs go directly to Oak

    # Add sample column (needed by analyze even for single sample)
    try:
        import pandas as pd
        df = pd.read_csv(corrected_tsv, sep='\t')
        sample_id = input_path.stem.replace('.fastq', '').replace('.gz', '').replace('.bam', '')
        df['sample'] = sample_id
        df.to_csv(corrected_tsv, sep='\t', index=False)
    except Exception:
        pass  # Not fatal — analyze can work without sample column

    # ── Step 2/3: Analysis ───────────────────────────────────────────────────
    print(f"\n[Step {step_correction + 1}/3] Analyzing results (single-sample)...")
    print("-" * 50)
    print("Note: DESeq2 requires multiple samples — skipped for single-sample run.")
    _t0 = _time.perf_counter()
    try:
        _run_analysis(
            corrected_tsv=corrected_tsv,
            output_dir=work_dir,
            genome_path=genome_path,
            annotation_path=annotation_path,
            args=args,
            n_samples=1,
        )
    except Exception as e:
        print(f"ERROR in analysis step: {e}", file=sys.stderr)
        # Not fatal — correction output is still usable

    # ── Junction aggregation ─────────────────────────────────────────────────
    if genome_path and annotation_path:
        _run_junction_aggregation(
            bam_path=bam_to_correct,
            genome_path=genome_path,
            annotation_path=annotation_path,
            output_dir=work_dir,
            config=vars(args),
        )

    print(f"[TIMING] Analysis:   {_time.perf_counter() - _t0:.1f}s")

    # ── Copy from scratch to Oak ─────────────────────────────────────────────
    if scratch_dir:
        print(f"\nCopying outputs to Oak: {output_dir}")
        _t0 = _time.perf_counter()
        # When bam_dir is set, BAMs went directly there (not scratch), so exclude_aligner_bams
        # only matters when BAMs are on scratch. When keep_aligner_bams is False (default),
        # skip per-aligner BAMs in the scratch→Oak sync to save disk space.
        _exclude_aligner = (not keep_aligner_bams) and (bam_dir is None)
        sync_to_oak(scratch_dir, output_dir, exclude_aligner_bams=_exclude_aligner)
        print(f"[TIMING] Sync to Oak: {_time.perf_counter() - _t0:.1f}s")
        _shutil.rmtree(scratch_dir, ignore_errors=True)

    # ── Provenance manifest ───────────────────────────────────────────────────
    manifest_path = tracker.write_manifest()

    # ── Summary ──────────────────────────────────────────────────────────────
    _pipeline_elapsed = _time.perf_counter() - _pipeline_start
    print("\n" + "=" * 70)
    print("Pipeline Complete!")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir}")
    print(f"  Corrected 3' ends: {corrected_tsv}")
    print(f"  Clusters:          {output_dir / 'cpa_clusters.tsv'}")
    print(f"  HTML report:       {output_dir / 'report.html'}")
    print(f"  Provenance:        {manifest_path}")
    print(f"\n[TIMING] Total pipeline: {_pipeline_elapsed:.1f}s ({_pipeline_elapsed/60:.1f} min)")

    return 0


# ---------------------------------------------------------------------------
# Chunked-alignment script generation
# ---------------------------------------------------------------------------

def _generate_chunked_pipeline(args) -> int:
    """
    Generate a dependency-chained set of scheduler scripts for chunked alignment.

    For FASTQ input: generates split → alignment arrays → merge → correct+analyze chain.
    For BAM input: prints a warning and returns 0 (inline correction will proceed).
    For manifest: generates per-sample chunked alignment chains + correction array + analysis.

    Returns 0 on success, 1 on error.
    """
    from . import split_command

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    scheduler    = getattr(args, 'scheduler', 'slurm')
    python_path  = getattr(args, 'python_path', None) or sys.executable
    rectify_src  = getattr(args, 'rectify_src', None)
    if rectify_src is None:
        # Auto-detect: walk up from this file to the package root
        rectify_src = str(Path(__file__).resolve().parent.parent.parent)

    target_reads = getattr(args, 'target_reads_per_chunk', split_command.DEFAULT_TARGET_READS)
    slurm_partition = getattr(args, 'slurm_partition', None)
    slurm_account   = getattr(args, 'slurm_account', None)
    uge_queue       = getattr(args, 'uge_queue', 'long.q')
    uge_pe          = getattr(args, 'uge_pe', 'smp')
    pbs_queue       = getattr(args, 'pbs_queue', 'workq')

    genome_path     = getattr(args, 'genome', None)
    annotation_path = getattr(args, 'annotation', None)

    sched_kwargs = dict(
        scheduler=scheduler,
        slurm_partition=slurm_partition,
        slurm_account=slurm_account,
        uge_queue=uge_queue,
        uge_pe=uge_pe,
        pbs_queue=pbs_queue,
    )

    def _genome_flag() -> str:
        return f'--genome "{genome_path}"' if genome_path else ''

    def _annot_flag() -> str:
        return f'--annotation "{annotation_path}"' if annotation_path else ''

    def _organism_flag() -> str:
        org = getattr(args, 'organism', None)
        return f'--organism "{org}"' if org else ''

    def _ref_flags() -> str:
        parts = [_genome_flag(), _annot_flag(), _organism_flag()]
        return ' '.join(p for p in parts if p)

    # ── Helper: correction+analysis script header (non-array) ────────────────
    def _make_correct_analyze_header(
        job_name: str, cores: int, mem_gb: int, time: str, log_dir: Path,
    ) -> str:
        log_pat = str(log_dir / '%j')
        return split_command._scheduler_headers(
            scheduler, job_name, 1, cores, mem_gb, time,
            str(log_dir), log_pat,
            is_array=False, **{k: v for k, v in sched_kwargs.items() if k != 'scheduler'},
        )

    # ── Helper: write an executable script ──────────────────────────────────
    def _write_script(path: Path, content: str) -> None:
        path.write_text(content)
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)

    # ── Helper: submission dependency syntax per scheduler ───────────────────
    def _submit_with_dep(script: Path, dep_var: str, dep_type: str = 'afterok') -> str:
        """Return a shell fragment that submits `script` after `dep_var` job IDs complete."""
        if scheduler == 'slurm':
            return (
                f'$(sbatch --dependency={dep_type}:{dep_var} {script} | awk \'{{print $4}}\')'
            )
        elif scheduler == 'uge':
            return f'$(qsub -hold_jid {dep_var} {script} | awk \'{{print $3}}\')'
        else:  # pbs
            return f'$(qsub -W depend=afterok:{dep_var} {script})'

    def _submit_no_dep(script: Path) -> str:
        if scheduler == 'slurm':
            return f'$(sbatch {script} | awk \'{{print $4}}\')'
        elif scheduler == 'uge':
            return f'$(qsub {script} | awk \'{{print $3}}\')'
        else:
            return f'$(qsub {script})'

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # Single-sample path
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    if getattr(args, 'input', None):
        input_path = Path(args.input)
        from .preprocess import detect_input_type
        input_type = detect_input_type(input_path)

        # (BAM inputs are handled in run() before calling this function)
        if input_type not in ('fastq', 'fastq.gz'):
            # Shouldn't normally reach here; run() screens for this.
            return 0

        # Derive sample prefix from FASTQ name
        name = input_path.name
        sample_prefix = name
        for sfx in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
            if name.endswith(sfx):
                sample_prefix = name[:-len(sfx)]
                break

        # Count reads and compute n_chunks
        print(f"Counting reads in {input_path} ...")
        n_reads = split_command.count_reads(input_path)
        n_chunks = split_command.compute_n_chunks(n_reads, target_reads)
        print(f"  {n_reads:,} reads → {n_chunks} chunks (~{n_reads // n_chunks:,} reads/chunk)")

        chunks_dir = output_dir / 'chunks'
        chunks_dir.mkdir(parents=True, exist_ok=True)
        log_dir = output_dir / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)

        # ── 00_split.sh ────────────────────────────────────────────────────
        split_header = split_command._scheduler_headers(
            scheduler, f'{sample_prefix}_split', 1, 4, 16, '2:00:00',
            str(log_dir), str(log_dir / '%j'),
            is_array=False, **{k: v for k, v in sched_kwargs.items() if k != 'scheduler'},
        )
        split_script_content = f"""#!/bin/bash
# RECTIFY — split FASTQ into chunks for parallel alignment
# Generated by: rectify run-all --chunked-alignment

{split_header}

set -euo pipefail

export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

SPLIT_CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-4}}}}}}
export OMP_NUM_THREADS=$SPLIT_CPUS
export OPENBLAS_NUM_THREADS=$SPLIT_CPUS
export MKL_NUM_THREADS=$SPLIT_CPUS
export LOKY_MAX_CPU_COUNT=$SPLIT_CPUS

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"
INPUT="{input_path}"
CHUNKS_DIR="{chunks_dir}"

echo "Python:  $($PYTHON --version)"
echo "Host:    $(hostname)"
echo "Start:   $(date)"
echo "Input:   $INPUT"
echo "Chunks:  $CHUNKS_DIR"

mkdir -p "$CHUNKS_DIR"
cd "$RECTIFY_SRC"

$PYTHON -m rectify split \\
    "$INPUT" \\
    -n {n_chunks} \\
    -o "$CHUNKS_DIR" \\
    --prefix "{sample_prefix}" \\
    --generate-slurm \\
    {_genome_flag()} \\
    {_annot_flag()} \\
    --python-path "{python_path}" \\
    --rectify-src "{rectify_src}" \\
    --scheduler {scheduler} \\
    --slurm-partition {slurm_partition} \\
    --slurm-account {slurm_account} \\
    --uge-queue {uge_queue} \\
    --uge-pe {uge_pe} \\
    --pbs-queue {pbs_queue}

echo "Split complete: $(date)"
ls -lh "$CHUNKS_DIR/"*chunk*.fastq.gz 2>/dev/null | head -5 || true
"""
        split_script = output_dir / '00_split.sh'
        _write_script(split_script, split_script_content)

        # ── Alignment scripts (in chunks_dir) ─────────────────────────────
        # Note: _generate_scripts writes into output_dir (= chunks_dir here)
        # Alignment scripts are generated now with n_chunks already known.
        align_scripts = split_command.generate_alignment_scripts(
            n_chunks=n_chunks,
            sample_prefix=sample_prefix,
            output_dir=chunks_dir,
            genome=genome_path,
            annotation=annotation_path,
            python_path=python_path,
            rectify_src=rectify_src,
            other_aligners=None,  # use defaults
            skip_map_pacbio=False,
            **{k: v for k, v in sched_kwargs.items() if k != 'scheduler'},
            scheduler=scheduler,
        )

        # ── run_correct_analyze.sh ─────────────────────────────────────────
        correct_header = _make_correct_analyze_header(
            f'{sample_prefix}_correct', 8, 32, '4:00:00', log_dir,
        )
        consensus_bam = output_dir / 'consensus' / f'{sample_prefix}.consensus.bam'
        ref_flags = _ref_flags()
        correct_analyze_content = f"""#!/bin/bash
# RECTIFY — correction + analysis after chunked alignment
# Generated by: rectify run-all --chunked-alignment
# Run after run_merge_consensus.sh completes.

{correct_header}

set -euo pipefail

export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

RECTIFY_CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-8}}}}}}
{split_command._thread_limits_block()}

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"
CONSENSUS_BAM="{consensus_bam}"
OUTDIR="{output_dir}"

echo "Python:  $($PYTHON --version)"
echo "Host:    $(hostname)"
echo "Start:   $(date)"
echo "Input:   $CONSENSUS_BAM"

[ -f "$CONSENSUS_BAM" ] || {{ echo "ERROR: consensus BAM not found: $CONSENSUS_BAM" >&2; exit 1; }}

cd "$RECTIFY_SRC"
$PYTHON -m rectify run-all \\
    "$CONSENSUS_BAM" \\
    {ref_flags} \\
    --streaming \\
    --threads $RECTIFY_CPUS \\
    -o "$OUTDIR"

echo "Done: $(date)"
"""
        correct_analyze_script = output_dir / 'run_correct_analyze.sh'
        _write_script(correct_analyze_script, correct_analyze_content)

        # ── submit_pipeline.sh ─────────────────────────────────────────────
        mpb_script  = align_scripts['mpb']
        others_script = align_scripts['others']
        merge_script  = align_scripts['merge']

        if scheduler == 'slurm':
            submit_lines = [
                '#!/bin/bash',
                '# Submit the full chunked-alignment pipeline (SLURM)',
                '# Generated by: rectify run-all --chunked-alignment',
                '',
                f'SPLIT_JOB={_submit_no_dep(split_script)}',
                f'echo "Split job: $SPLIT_JOB"',
            ]
            if mpb_script:
                submit_lines += [
                    f'MPB_JOB=$(sbatch --dependency=afterok:$SPLIT_JOB {mpb_script} | awk \'{{print $4}}\')',
                    'echo "mapPacBio array: $MPB_JOB"',
                ]
            submit_lines += [
                f'OTHERS_JOB=$(sbatch --dependency=afterok:$SPLIT_JOB {others_script} | awk \'{{print $4}}\')',
                'echo "Others array: $OTHERS_JOB"',
            ]
            if mpb_script:
                merge_dep = '--dependency=afterok:$MPB_JOB:$OTHERS_JOB'
            else:
                merge_dep = '--dependency=afterok:$OTHERS_JOB'
            submit_lines += [
                f'MERGE_JOB=$(sbatch {merge_dep} {merge_script} | awk \'{{print $4}}\')',
                'echo "Merge+consensus: $MERGE_JOB"',
                f'CA_JOB=$(sbatch --dependency=afterok:$MERGE_JOB {correct_analyze_script} | awk \'{{print $4}}\')',
                'echo "Correct+analyze: $CA_JOB"',
                'echo ""',
                'echo "Pipeline submitted. Monitor with: squeue -u $USER"',
            ]
        elif scheduler == 'uge':
            submit_lines = [
                '#!/bin/bash',
                '# Submit the full chunked-alignment pipeline (UGE/SGE)',
                '',
                f'SPLIT_JOB={_submit_no_dep(split_script)}',
                f'echo "Split job: $SPLIT_JOB"',
            ]
            if mpb_script:
                submit_lines += [
                    f'MPB_JOB=$(qsub -hold_jid $SPLIT_JOB {mpb_script} | awk \'{{print $3}}\')',
                    'echo "mapPacBio array: $MPB_JOB"',
                ]
            submit_lines += [
                f'OTHERS_JOB=$(qsub -hold_jid $SPLIT_JOB {others_script} | awk \'{{print $3}}\')',
                'echo "Others array: $OTHERS_JOB"',
            ]
            if mpb_script:
                merge_hold = '-hold_jid $MPB_JOB,$OTHERS_JOB'
            else:
                merge_hold = '-hold_jid $OTHERS_JOB'
            submit_lines += [
                f'MERGE_JOB=$(qsub {merge_hold} {merge_script} | awk \'{{print $3}}\')',
                'echo "Merge+consensus: $MERGE_JOB"',
                f'CA_JOB=$(qsub -hold_jid $MERGE_JOB {correct_analyze_script} | awk \'{{print $3}}\')',
                'echo "Correct+analyze: $CA_JOB"',
                'echo ""',
                'echo "Pipeline submitted. Monitor with: qstat -u $USER"',
            ]
        else:  # pbs
            submit_lines = [
                '#!/bin/bash',
                '# Submit the full chunked-alignment pipeline (PBS/Torque)',
                '',
                f'SPLIT_JOB={_submit_no_dep(split_script)}',
                f'echo "Split job: $SPLIT_JOB"',
            ]
            if mpb_script:
                submit_lines += [
                    f'MPB_JOB=$(qsub -W depend=afterok:$SPLIT_JOB {mpb_script})',
                    'echo "mapPacBio array: $MPB_JOB"',
                ]
            submit_lines += [
                f'OTHERS_JOB=$(qsub -W depend=afterok:$SPLIT_JOB {others_script})',
                'echo "Others array: $OTHERS_JOB"',
            ]
            if mpb_script:
                merge_dep = '-W depend=afterok:$MPB_JOB:$OTHERS_JOB'
            else:
                merge_dep = '-W depend=afterok:$OTHERS_JOB'
            submit_lines += [
                f'MERGE_JOB=$(qsub {merge_dep} {merge_script})',
                'echo "Merge+consensus: $MERGE_JOB"',
                f'CA_JOB=$(qsub -W depend=afterok:$MERGE_JOB {correct_analyze_script})',
                'echo "Correct+analyze: $CA_JOB"',
                'echo ""',
                'echo "Pipeline submitted. Monitor with: qstat -u $USER"',
            ]

        submit_pipeline = output_dir / 'submit_pipeline.sh'
        _write_script(submit_pipeline, '\n'.join(submit_lines) + '\n')

        print("\nGenerated chunked-alignment pipeline scripts:")
        print(f"  Split:           {split_script}")
        print(f"  mapPacBio array: {align_scripts['mpb'] or '(skipped)'}")
        print(f"  Others array:    {align_scripts['others']}")
        print(f"  Merge+consensus: {align_scripts['merge']}")
        print(f"  Correct+analyze: {correct_analyze_script}")
        print(f"\nTo launch the pipeline:")
        print(f"  bash {submit_pipeline}")
        return 0

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # Multi-sample path (manifest)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    if getattr(args, 'manifest', None):
        from .batch_command import parse_manifest
        manifest_path = Path(args.manifest)
        samples = parse_manifest(manifest_path)
        if not samples:
            print("ERROR: No samples found in manifest.", file=sys.stderr)
            return 1

        ref_flags = _ref_flags()
        log_dir = output_dir / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)

        fastq_samples = []
        bam_samples   = []

        for s in samples:
            p = Path(s.get('path', s.get('bam_path', '')))
            from .preprocess import detect_input_type
            itype = detect_input_type(p)
            if itype in ('fastq', 'fastq.gz'):
                fastq_samples.append(s)
            else:
                bam_samples.append(s)

        # Per-sample split job IDs (for submit_all.sh)
        # We'll track them as shell variable names: SPLIT_<sample_id>
        sample_split_scripts    = {}
        sample_submit_pipelines = {}

        for s in fastq_samples:
            sample_id  = s['sample_id']
            input_path = Path(s.get('path', s.get('bam_path', '')))
            sample_dir = output_dir / sample_id
            sample_dir.mkdir(parents=True, exist_ok=True)
            chunks_dir = sample_dir / 'chunks'
            chunks_dir.mkdir(parents=True, exist_ok=True)
            slog_dir = sample_dir / 'logs'
            slog_dir.mkdir(parents=True, exist_ok=True)

            # Derive prefix
            name = input_path.name
            prefix = name
            for sfx in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
                if name.endswith(sfx):
                    prefix = name[:-len(sfx)]
                    break

            # Count reads
            print(f"  [{sample_id}] Counting reads ...")
            n_reads  = split_command.count_reads(input_path)
            n_chunks = split_command.compute_n_chunks(n_reads, target_reads)
            print(f"  [{sample_id}] {n_reads:,} reads → {n_chunks} chunks")

            # 00_split.sh for this sample
            split_header = split_command._scheduler_headers(
                scheduler, f'{sample_id}_split', 1, 4, 16, '2:00:00',
                str(slog_dir), str(slog_dir / '%j'),
                is_array=False, **{k: v for k, v in sched_kwargs.items() if k != 'scheduler'},
            )
            split_content = f"""#!/bin/bash
# RECTIFY — split FASTQ for {sample_id}

{split_header}

set -euo pipefail
export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

SPLIT_CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-4}}}}}}
export OMP_NUM_THREADS=$SPLIT_CPUS
export OPENBLAS_NUM_THREADS=$SPLIT_CPUS
export MKL_NUM_THREADS=$SPLIT_CPUS
export LOKY_MAX_CPU_COUNT=$SPLIT_CPUS

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"

mkdir -p "{chunks_dir}"
cd "$RECTIFY_SRC"

echo "Start: $(date)  Host: $(hostname)"

$PYTHON -m rectify split \\
    "{input_path}" \\
    -n {n_chunks} \\
    -o "{chunks_dir}" \\
    --prefix "{prefix}" \\
    --generate-slurm \\
    {_genome_flag()} \\
    {_annot_flag()} \\
    --python-path "{python_path}" \\
    --rectify-src "{rectify_src}" \\
    --scheduler {scheduler} \\
    --slurm-partition {slurm_partition} \\
    --slurm-account {slurm_account} \\
    --uge-queue {uge_queue} \\
    --uge-pe {uge_pe} \\
    --pbs-queue {pbs_queue}

echo "Done: $(date)"
"""
            split_script = sample_dir / '00_split.sh'
            _write_script(split_script, split_content)
            sample_split_scripts[sample_id] = split_script

            # Alignment scripts in chunks_dir
            align_scripts = split_command.generate_alignment_scripts(
                n_chunks=n_chunks,
                sample_prefix=prefix,
                output_dir=chunks_dir,
                genome=genome_path,
                annotation=annotation_path,
                python_path=python_path,
                rectify_src=rectify_src,
                scheduler=scheduler,
                other_aligners=None,
                skip_map_pacbio=False,
                slurm_partition=slurm_partition,
                slurm_account=slurm_account,
                uge_queue=uge_queue,
                uge_pe=uge_pe,
                pbs_queue=pbs_queue,
            )

            # Per-sample submit_pipeline.sh (split → align → merge)
            # Does NOT include correct+analyze — that is a combined array step later
            consensus_bam = chunks_dir / 'consensus' / f'{prefix}.consensus.bam'
            mpb_script    = align_scripts['mpb']
            others_script = align_scripts['others']
            merge_script  = align_scripts['merge']

            if scheduler == 'slurm':
                sub_lines = [
                    '#!/bin/bash',
                    f'# Submit alignment pipeline for {sample_id}',
                    '',
                    f'SPLIT_JOB={_submit_no_dep(split_script)}',
                    f'echo "  [{sample_id}] Split: $SPLIT_JOB"',
                ]
                if mpb_script:
                    sub_lines += [
                        f'MPB_JOB=$(sbatch --dependency=afterok:$SPLIT_JOB {mpb_script} | awk \'{{print $4}}\')',
                    ]
                sub_lines += [
                    f'OTHERS_JOB=$(sbatch --dependency=afterok:$SPLIT_JOB {others_script} | awk \'{{print $4}}\')',
                ]
                if mpb_script:
                    sub_lines += [
                        f'MERGE_JOB=$(sbatch --dependency=afterok:$MPB_JOB:$OTHERS_JOB {merge_script} | awk \'{{print $4}}\')',
                    ]
                else:
                    sub_lines += [
                        f'MERGE_JOB=$(sbatch --dependency=afterok:$OTHERS_JOB {merge_script} | awk \'{{print $4}}\')',
                    ]
                sub_lines += ['echo $MERGE_JOB']  # Print merge job ID so submit_all.sh can capture it
            elif scheduler == 'uge':
                sub_lines = [
                    '#!/bin/bash',
                    f'# Submit alignment pipeline for {sample_id}',
                    '',
                    f'SPLIT_JOB={_submit_no_dep(split_script)}',
                ]
                if mpb_script:
                    sub_lines += [f'MPB_JOB=$(qsub -hold_jid $SPLIT_JOB {mpb_script} | awk \'{{print $3}}\')']
                sub_lines += [f'OTHERS_JOB=$(qsub -hold_jid $SPLIT_JOB {others_script} | awk \'{{print $3}}\')']
                if mpb_script:
                    sub_lines += [f'MERGE_JOB=$(qsub -hold_jid $MPB_JOB,$OTHERS_JOB {merge_script} | awk \'{{print $3}}\')']
                else:
                    sub_lines += [f'MERGE_JOB=$(qsub -hold_jid $OTHERS_JOB {merge_script} | awk \'{{print $3}}\')']
                sub_lines += ['echo $MERGE_JOB']
            else:  # pbs
                sub_lines = [
                    '#!/bin/bash',
                    f'# Submit alignment pipeline for {sample_id}',
                    '',
                    f'SPLIT_JOB={_submit_no_dep(split_script)}',
                ]
                if mpb_script:
                    sub_lines += [f'MPB_JOB=$(qsub -W depend=afterok:$SPLIT_JOB {mpb_script})']
                sub_lines += [f'OTHERS_JOB=$(qsub -W depend=afterok:$SPLIT_JOB {others_script})']
                if mpb_script:
                    sub_lines += [f'MERGE_JOB=$(qsub -W depend=afterok:$MPB_JOB:$OTHERS_JOB {merge_script})']
                else:
                    sub_lines += [f'MERGE_JOB=$(qsub -W depend=afterok:$OTHERS_JOB {merge_script})']
                sub_lines += ['echo $MERGE_JOB']

            sample_submit = sample_dir / 'submit_sample_alignment.sh'
            _write_script(sample_submit, '\n'.join(sub_lines) + '\n')
            sample_submit_pipelines[sample_id] = (sample_submit, consensus_bam)

        # ── Correction array script ────────────────────────────────────────
        # Covers both FASTQ samples (consensus BAM) and BAM samples (direct)
        all_samples_for_correct = []
        for s in fastq_samples:
            _, consensus_bam = sample_submit_pipelines[s['sample_id']]
            all_samples_for_correct.append((s['sample_id'], consensus_bam))
        for s in bam_samples:
            all_samples_for_correct.append(
                (s['sample_id'], Path(s.get('path', s.get('bam_path', ''))))
            )

        n_correct_tasks = len(all_samples_for_correct)
        correct_header = split_command._scheduler_headers(
            scheduler, 'rectify_correct', n_correct_tasks, 8, 32, '4:00:00',
            str(log_dir), str(log_dir / '%A_%a'),
            is_array=True, **{k: v for k, v in sched_kwargs.items() if k != 'scheduler'},
        )

        # Build arrays for bash
        sample_ids_arr   = ' '.join(f'"{sid}"' for sid, _ in all_samples_for_correct)
        consensus_bams_arr = ' '.join(f'"{bam}"' for _, bam in all_samples_for_correct)

        correct_array_content = f"""#!/bin/bash
# RECTIFY — correction array (one task per sample)
# Generated by: rectify run-all --chunked-alignment --manifest

{correct_header}

set -euo pipefail
export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

# ── Scheduler-agnostic task ID ────────────────────────────────────────────
if   [ -n "${{SLURM_ARRAY_TASK_ID:-}}" ]; then
    RECTIFY_TASK_ID=$SLURM_ARRAY_TASK_ID
    RECTIFY_CPUS=${{SLURM_CPUS_PER_TASK:-8}}
elif [ -n "${{SGE_TASK_ID:-}}" ]; then
    RECTIFY_TASK_ID=$(( SGE_TASK_ID - 1 ))
    RECTIFY_CPUS=${{NSLOTS:-8}}
elif [ -n "${{PBS_ARRAY_INDEX:-}}" ]; then
    RECTIFY_TASK_ID=$(( PBS_ARRAY_INDEX - 1 ))
    RECTIFY_CPUS=${{PBS_NUM_PPN:-8}}
else
    echo "ERROR: no scheduler array task variable found" >&2; exit 1
fi
# ─────────────────────────────────────────────────────────────────────────

{split_command._thread_limits_block()}

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"

SAMPLE_IDS=({sample_ids_arr})
CONSENSUS_BAMS=({consensus_bams_arr})

SAMPLE_ID="${{SAMPLE_IDS[$RECTIFY_TASK_ID]}}"
BAM="${{CONSENSUS_BAMS[$RECTIFY_TASK_ID]}}"
SAMPLE_OUTDIR="$OUTDIR/$SAMPLE_ID"
mkdir -p "$SAMPLE_OUTDIR"

echo "Python:  $($PYTHON --version)"
echo "Host:    $(hostname)"
echo "Task:    $RECTIFY_TASK_ID  Sample: $SAMPLE_ID  CPUs: $RECTIFY_CPUS"
echo "Input:   $BAM"
echo "Start:   $(date)"

[ -f "$BAM" ] || {{ echo "ERROR: BAM not found: $BAM" >&2; exit 1; }}

cd "$RECTIFY_SRC"
$PYTHON -m rectify correct \\
    "$BAM" \\
    {_genome_flag()} \\
    {_annot_flag()} \\
    -o "$SAMPLE_OUTDIR/corrected_3ends.tsv" \\
    --streaming \\
    --threads $RECTIFY_CPUS

echo "Done: $(date)"
"""
        correct_array_script = output_dir / 'run_correct_array.sh'
        _write_script(correct_array_script, correct_array_content)

        # ── Analysis script ────────────────────────────────────────────────
        analyze_header = _make_correct_analyze_header(
            'rectify_analyze', 8, 64, '8:00:00', log_dir,
        )
        analyze_manifest_path = output_dir / 'corrected_manifest.tsv'
        analyze_content = f"""#!/bin/bash
# RECTIFY — combined analysis (DESeq2, GO, motifs)
# Generated by: rectify run-all --chunked-alignment --manifest
# Run after run_correct_array.sh completes.

{analyze_header}

set -euo pipefail
export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

ANALYZE_CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-8}}}}}}
{split_command._thread_limits_block('$ANALYZE_CPUS')}

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"
MANIFEST="{analyze_manifest_path}"
COMBINED_DIR="$OUTDIR/combined"

echo "Python:  $($PYTHON --version)"
echo "Host:    $(hostname)"
echo "Start:   $(date)"

mkdir -p "$COMBINED_DIR"

# Write corrected manifest (sample_id, path, condition columns from original manifest)
# This manifest is pre-written by the script-generation step; just verify it exists.
[ -f "$MANIFEST" ] || {{ echo "ERROR: corrected manifest not found: $MANIFEST" >&2; exit 1; }}

cd "$RECTIFY_SRC"
$PYTHON -m rectify analyze /dev/null \\
    --manifest "$MANIFEST" \\
    {_genome_flag()} \\
    {_annot_flag()} \\
    --threads $ANALYZE_CPUS \\
    -o "$COMBINED_DIR"

echo "Done: $(date)"
"""
        analyze_script = output_dir / 'run_analyze.sh'
        _write_script(analyze_script, analyze_content)

        # Write the corrected manifest (paths to per-sample corrected TSVs)
        import csv
        manifest_rows = []
        for s in samples:
            row = {
                'sample_id': s['sample_id'],
                'path': str(output_dir / s['sample_id'] / 'corrected_3ends.tsv'),
            }
            if 'condition' in s:
                row['condition'] = s['condition']
            manifest_rows.append(row)

        fieldnames = ['sample_id', 'path']
        if any('condition' in r for r in manifest_rows):
            fieldnames.append('condition')
        with open(analyze_manifest_path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t',
                                    extrasaction='ignore')
            writer.writeheader()
            writer.writerows(manifest_rows)

        # ── submit_all.sh ──────────────────────────────────────────────────
        if scheduler == 'slurm':
            submit_all_lines = [
                '#!/bin/bash',
                '# Submit all per-sample alignment chains, then correction array, then analysis',
                '# Generated by: rectify run-all --chunked-alignment --manifest',
                '',
                '# ── Per-sample alignment (all run in parallel) ────────────────────',
                'MERGE_JOBS=""',
            ]
            for sid, (sample_submit, _) in sample_submit_pipelines.items():
                var = f'MERGE_{sid.replace("-", "_").replace(".", "_")}'
                submit_all_lines += [
                    f'{var}=$(bash {sample_submit})',
                    f'echo "  [{sid}] merge job: ${{{var}}}"',
                    f'MERGE_JOBS="${{MERGE_JOBS}}:${{{var}}}"',
                ]
            # BAM samples: no alignment needed — correct_array can start immediately
            if bam_samples and not fastq_samples:
                submit_all_lines += [
                    '',
                    '# ── Correction array (no alignment dependency for BAM-only manifest) ─',
                    f'CORRECT_JOB=$(sbatch {correct_array_script} | awk \'{{print $4}}\')',
                ]
            else:
                submit_all_lines += [
                    '',
                    '# ── Correction array (after all merge jobs complete) ─────────────────',
                    'MERGE_JOBS="${MERGE_JOBS#:}"  # strip leading colon',
                    f'CORRECT_JOB=$(sbatch --dependency=afterok:$MERGE_JOBS {correct_array_script} | awk \'{{print $4}}\')',
                ]
            submit_all_lines += [
                'echo "Correction array: $CORRECT_JOB"',
                '',
                '# ── Analysis (after all corrections complete) ────────────────────────',
                f'ANALYZE_JOB=$(sbatch --dependency=afterok:$CORRECT_JOB {analyze_script} | awk \'{{print $4}}\')',
                'echo "Analysis: $ANALYZE_JOB"',
                'echo ""',
                'echo "Full pipeline submitted. Monitor with: squeue -u $USER"',
            ]
        elif scheduler == 'uge':
            submit_all_lines = [
                '#!/bin/bash',
                '# Submit all per-sample alignment chains, then correction array, then analysis',
                '',
                'MERGE_JOBS=""',
            ]
            for sid, (sample_submit, _) in sample_submit_pipelines.items():
                var = f'MERGE_{sid.replace("-", "_").replace(".", "_")}'
                submit_all_lines += [
                    f'{var}=$(bash {sample_submit})',
                    f'MERGE_JOBS="${{MERGE_JOBS}},${{{var}}}"',
                ]
            if bam_samples and not fastq_samples:
                submit_all_lines += [
                    f'CORRECT_JOB=$(qsub {correct_array_script} | awk \'{{print $3}}\')',
                ]
            else:
                submit_all_lines += [
                    'MERGE_JOBS="${MERGE_JOBS#,}"',
                    f'CORRECT_JOB=$(qsub -hold_jid $MERGE_JOBS {correct_array_script} | awk \'{{print $3}}\')',
                ]
            submit_all_lines += [
                'echo "Correction array: $CORRECT_JOB"',
                f'ANALYZE_JOB=$(qsub -hold_jid $CORRECT_JOB {analyze_script} | awk \'{{print $3}}\')',
                'echo "Analysis: $ANALYZE_JOB"',
                'echo "Pipeline submitted. Monitor with: qstat -u $USER"',
            ]
        else:  # pbs
            submit_all_lines = [
                '#!/bin/bash',
                '# Submit all per-sample alignment chains, then correction array, then analysis',
                '',
                'MERGE_JOBS=""',
            ]
            for sid, (sample_submit, _) in sample_submit_pipelines.items():
                var = f'MERGE_{sid.replace("-", "_").replace(".", "_")}'
                submit_all_lines += [
                    f'{var}=$(bash {sample_submit})',
                    f'MERGE_JOBS="${{MERGE_JOBS}}:${{{var}}}"',
                ]
            if bam_samples and not fastq_samples:
                submit_all_lines += [f'CORRECT_JOB=$(qsub {correct_array_script})']
            else:
                submit_all_lines += [
                    'MERGE_JOBS="${MERGE_JOBS#:}"',
                    f'CORRECT_JOB=$(qsub -W depend=afterok:$MERGE_JOBS {correct_array_script})',
                ]
            submit_all_lines += [
                'echo "Correction array: $CORRECT_JOB"',
                f'ANALYZE_JOB=$(qsub -W depend=afterok:$CORRECT_JOB {analyze_script})',
                'echo "Analysis: $ANALYZE_JOB"',
                'echo "Pipeline submitted. Monitor with: qstat -u $USER"',
            ]

        submit_all = output_dir / 'submit_all.sh'
        _write_script(submit_all, '\n'.join(submit_all_lines) + '\n')

        print("\nGenerated multi-sample chunked-alignment pipeline scripts:")
        for sid in [s['sample_id'] for s in fastq_samples]:
            print(f"  [{sid}] {sample_split_scripts[sid]}")
            print(f"  [{sid}] {sample_submit_pipelines[sid][0]}")
        if bam_samples:
            print(f"  {len(bam_samples)} BAM sample(s): direct to correction (no alignment)")
        print(f"  Correction array: {correct_array_script}")
        print(f"  Analysis:         {analyze_script}")
        print(f"\nTo launch the full pipeline:")
        print(f"  bash {submit_all}")
        return 0

    print("ERROR: --chunked-alignment requires either a FASTQ input or --manifest.", file=sys.stderr)
    return 1


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run(args: argparse.Namespace) -> None:
    """Dispatch to single-sample or multi-sample pipeline."""
    from rectify.slurm import set_thread_limits
    set_thread_limits(getattr(args, 'threads', None))

    if getattr(args, 'manifest', None) and getattr(args, 'input', None):
        print(
            "ERROR: Cannot combine a positional input file with --manifest. "
            "Provide one or the other.",
            file=sys.stderr,
        )
        sys.exit(1)

    # --chunked-alignment: generate scripts and exit (don't run inline),
    # UNLESS the input is a BAM (chunked alignment not needed → run inline).
    if getattr(args, 'chunked_alignment', False):
        _resolve_reference_paths(args)
        # Check if it's a BAM — if so _generate_chunked_pipeline prints a warning
        # and returns 0; we then fall through to inline correction.
        input_path_str = getattr(args, 'input', None)
        is_bam_input = False
        if input_path_str:
            from .preprocess import detect_input_type
            itype = detect_input_type(Path(str(input_path_str)))
            is_bam_input = itype not in ('fastq', 'fastq.gz')

        if is_bam_input:
            print(
                f"Warning: --chunked-alignment is only needed for FASTQ inputs. "
                f"Input is a BAM — proceeding with inline BAM correction.",
                file=sys.stderr,
            )
            # Fall through to normal single-sample processing below
        else:
            rc = _generate_chunked_pipeline(args)
            sys.exit(rc)

    if getattr(args, 'manifest', None):
        rc = _run_multi_sample(args)
    elif getattr(args, 'input', None):
        rc = _run_single_sample(args)
    else:
        print(
            "ERROR: Provide either an input file (FASTQ or BAM) or --manifest.",
            file=sys.stderr,
        )
        rc = 1

    sys.exit(rc)
