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
        --profile rectify/slurm_profiles/sherlock_larsms.yaml

Author: Kevin R. Roy
Date: 2026-03-28
"""

import argparse
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

def _consensus_bam_path(sample_id: str, sample_output_dir: Path) -> Path:
    """Canonical path for the consensus BAM for a sample."""
    return sample_output_dir / f"{sample_id}.consensus.bam"


def _run_alignment(
    input_path: Path,
    sample_id: str,
    sample_output_dir: Path,
    genome_path: Path,
    annotation_path: Optional[Path],
    threads: int,
    parallel_aligners: bool = False,
    junction_aligners: Optional[List[str]] = None,
    chimeric_consensus: bool = False,
    ultra_path: str = 'uLTRA',
    desalt_path: str = 'deSALT',
    mapPacBio_chunks: int = 1,
) -> Path:
    """
    Run consensus alignment or return existing consensus.bam.

    Default aligners: minimap2 + mapPacBio + gapmm2.
    Opt-in: uLTRA and/or deSALT via junction_aligners.

    Skips automatically if the consensus.bam already exists — safe to re-run.
    Returns the path to the consensus BAM.
    """
    consensus_bam = _consensus_bam_path(sample_id, sample_output_dir)

    if consensus_bam.exists():
        print(f"    Skipping alignment — consensus.bam exists: {consensus_bam}")
        return consensus_bam

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

    # Temporarily hide SLURM array env vars so consensus.py runs in single-sample
    # mode. In run-all, SLURM array indices are for sample-level parallelism;
    # we never want within-sample read-level partitioning here.
    import os as _os
    _array_vars = {k: _os.environ.pop(k) for k in (
        'SLURM_ARRAY_TASK_ID', 'SLURM_ARRAY_TASK_COUNT',
        'SLURM_ARRAY_TASK_MAX', 'SLURM_ARRAY_TASK_MIN',
        'SLURM_ARRAY_TASK_STEP',
    ) if k in _os.environ}
    try:
        rc = run_align(align_args)
    finally:
        _os.environ.update(_array_vars)
    if rc != 0:
        raise RuntimeError(f"Triple-aligner failed for {input_path}")

    if not consensus_bam.exists():
        raise RuntimeError(
            f"Alignment completed but consensus.bam not found: {consensus_bam}"
        )

    return consensus_bam


# ---------------------------------------------------------------------------
# Correction helper
# ---------------------------------------------------------------------------

def _bam_has_md_tags(bam_path: Path) -> bool:
    """Check if a BAM file's reads have MD tags (sample first 10 mapped reads)."""
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
    except Exception:
        pass
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
    # Disable them gracefully when MD tags are absent (e.g. triple-aligner
    # consensus BAMs generated by older runs without --MD).
    has_md = _bam_has_md_tags(bam_path)
    # --polya-sequenced defaults to True in run-all; --no-polya-sequenced disables it
    polya_seq = not getattr(args, 'no_polya_sequenced', False)
    skip_indel = not has_md  # Can't correct indels without MD tags
    skip_variant = not has_md

    if not has_md and polya_seq:
        print(
            "    Note: BAM lacks MD tags — indel correction disabled. "
            "Re-run alignment with current RECTIFY to enable it."
        )

    # Auto-generate report path alongside the corrected TSV
    report_path = corrected_tsv.parent / (corrected_tsv.stem + '_report.html')

    correct_args = argparse.Namespace(
        input=bam_path,
        genome=genome_path,
        annotation=annotation_path,
        output=corrected_tsv,
        netseq_dir=getattr(args, 'netseq_dir', None),
        organism=getattr(args, 'organism', None),
        aligner=getattr(args, 'aligner', 'minimap2'),
        polya_sequenced=polya_seq,
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
) -> Optional[Path]:
    """
    Run junction aggregation with partial rescue.
    Only runs if annotation is GFF/GFF3 (has intron features).
    Returns path to junctions TSV, or None if skipped/failed.
    """
    has_gff = annotation_path.suffix.lower() in ('.gff', '.gff3')
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

        n_rescued = partial_results['summary']['total_rescued']
        n_ambiguous = partial_results['summary']['total_ambiguous']
        print(f"  Rescued {n_rescued} partial crossings ({n_ambiguous} ambiguous)")

        junction_df = merge_with_partial_evidence(
            junction_df, partial_results, ambiguous_mode='proportional',
        )
        export_junctions(junction_df, str(junctions_tsv), format='tsv')
        print(f"  Junctions written to {junctions_tsv}")

        # Provenance
        if junctions_tsv.exists():
            from ..utils.provenance import init_provenance
            prov = init_provenance(output_dir, description="RECTIFY junction aggregation")
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
    sample_id = sample['sample_id']
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
                try:
                    bam_to_correct = _run_alignment(
                        input_path=input_path,
                        sample_id=sample_id,
                        sample_output_dir=sample_output,
                        genome_path=genome_path,
                        annotation_path=annotation_path,
                        threads=getattr(args, 'threads', 4),
                        parallel_aligners=getattr(args, 'parallel_aligners', False),
                        junction_aligners=getattr(args, 'junction_aligners', []),
                        chimeric_consensus=getattr(args, 'chimeric_consensus', False),
                        ultra_path=getattr(args, 'ultra_path', 'uLTRA'),
                        desalt_path=getattr(args, 'desalt_path', 'deSALT'),
                        mapPacBio_chunks=getattr(args, 'mapPacBio_chunks', 1),
                    )
                    log.write(f"Alignment complete: {bam_to_correct}\n")
                except Exception as e:
                    log.write(f"Alignment failed: {e}\n")
                    return sample_id, 1
            elif input_type in ('fastq', 'fastq.gz'):
                # FASTQ but alignment skipped — look for existing consensus.bam
                consensus_bam = _consensus_bam_path(sample_id, sample_output)
                if consensus_bam.exists():
                    bam_to_correct = consensus_bam
                    log.write(f"Using existing consensus.bam: {consensus_bam}\n")
                else:
                    log.write(
                        f"ERROR: --skip-alignment set but no consensus.bam found: {consensus_bam}\n"
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
      Step 0 (if FASTQ): triple-aligner alignment → consensus.bam
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

    # ── Scratch staging setup ────────────────────────────────────────────────
    # If $SCRATCH is available (Sherlock), run all I/O there and rsync back.
    # This avoids Oak NFS contention across concurrent SLURM array tasks.
    # The consensus BAM is always copied back to Oak even on fresh alignments
    # so it survives $SCRATCH's 90-day auto-purge.
    from ..slurm import make_job_scratch_dir, sync_to_oak
    import shutil as _shutil

    use_scratch = getattr(args, 'use_scratch', True)  # default on when scratch available
    scratch_dir = make_job_scratch_dir('rectify_single') if use_scratch else None
    work_dir = scratch_dir if scratch_dir else output_dir

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
        # Align directly to work_dir (scratch if available, else output_dir)
        bam_to_correct = _run_alignment(
            input_path=input_path,
            sample_id=sample_id,
            sample_output_dir=work_dir,
            genome_path=genome_path,
            annotation_path=annotation_path,
            threads=getattr(args, 'threads', 4),
            parallel_aligners=getattr(args, 'parallel_aligners', False),
            junction_aligners=getattr(args, 'junction_aligners', []),
            chimeric_consensus=getattr(args, 'chimeric_consensus', False),
            ultra_path=getattr(args, 'ultra_path', 'uLTRA'),
            desalt_path=getattr(args, 'desalt_path', 'deSALT'),
        )
        print(f"\nAlignment complete: {bam_to_correct}")
        print(f"[TIMING] Alignment: {_time.perf_counter() - _t0:.1f}s")
        tracker.record_step('align', input_files=[input_path], output_files=[bam_to_correct])
        step_correction = 2
    elif input_type in ('fastq', 'fastq.gz'):
        sample_id = input_path.stem.replace('.fastq', '').replace('.gz', '')
        consensus_bam = _consensus_bam_path(sample_id, output_dir)
        if consensus_bam.exists():
            print(f"\n[Step 1/3] Alignment — using existing consensus.bam")
            if scratch_dir:
                # Stage existing BAM to scratch
                print(f"  Staging to scratch: {scratch_dir}")
                _shutil.copy2(consensus_bam, scratch_dir / consensus_bam.name)
                bai = Path(str(consensus_bam) + '.bai')
                if bai.exists():
                    _shutil.copy2(bai, scratch_dir / bai.name)
                bam_to_correct = scratch_dir / consensus_bam.name
                tracker.register_staged(bam_to_correct, consensus_bam)
            else:
                bam_to_correct = consensus_bam
        else:
            print(
                f"ERROR: --skip-alignment set but no consensus.bam found: {consensus_bam}",
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
        )

    print(f"[TIMING] Analysis:   {_time.perf_counter() - _t0:.1f}s")

    # ── Copy from scratch to Oak ─────────────────────────────────────────────
    if scratch_dir:
        print(f"\nCopying outputs to Oak: {output_dir}")
        _t0 = _time.perf_counter()
        sync_to_oak(scratch_dir, output_dir)
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
# Main entry point
# ---------------------------------------------------------------------------

def run(args: argparse.Namespace) -> None:
    """Dispatch to single-sample or multi-sample pipeline."""
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
