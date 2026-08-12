#!/usr/bin/env python3
"""
RECTIFY Analyze Command

CLI entry point for downstream analysis of corrected 3' end positions.

Thin orchestrator: heavy lifting lives in `rectify/core/analyze/`:
- `loaders.py` — TSV / GTF / position-index loading
- `exclusions.py` — mito / rDNA auto-detection
- `bedgraph.py` — per-condition strand-specific bedgraph emission
- `manifest.py` — memory-efficient manifest-mode runner

Author: Kevin R. Roy
Date: 2026-03-17
"""

import argparse
from pathlib import Path

import pandas as pd

from ..analyze import (
    cluster_cpa_sites,
    build_cluster_count_matrix,
    run_deseq2_gene_level,
    run_deseq2_cluster_level,
    detect_control_samples,
    create_sample_metadata,
    run_pca_analysis,
    plot_pca,
    plot_sample_heatmap,
    run_go_enrichment,
    plot_go_enrichment,
    run_differential_motif_analysis,
    summarize_motif_results,
    analyze_cluster_shifts,
    get_top_shifted_genes,
    plot_gene_browser,
    plot_shift_summary,
    generate_analysis_summary,
    generate_html_report,
    run_3prime_distribution_analysis,
    run_5prime_distribution_analysis,
    run_transcript_body_distribution_analysis,
)
from ..analyze.clustering import annotate_clusters_with_genes, cluster_cpa_sites_adaptive
from ..analyze.deseq2 import extract_condition_from_sample
from ..analyze.loaders import load_corrected_positions, load_annotation
from ..analyze.exclusions import _detect_exclusion_regions
from ..analyze.bedgraph import generate_bedgraphs
from ..analyze.manifest import _run_analyze_manifest, _build_cluster_gene_attribution_table
from ...utils.provenance import init_provenance


def _collect_analyze_inputs(args: argparse.Namespace) -> dict:
    """Build the input-hash dict for skip-check and sidecar, resolving per-sample TSV paths."""
    inputs: dict = {}
    if getattr(args, 'manifest', None):
        inputs['manifest'] = str(args.manifest)
        try:
            _mdf = pd.read_csv(args.manifest, sep='\t')
            for _row in _mdf.itertuples():
                if hasattr(_row, 'sample_id') and hasattr(_row, 'path'):
                    inputs[f'tsv_{_row.sample_id}'] = str(_row.path)
        except Exception:
            pass
    else:
        inputs['input'] = str(args.input)
    if getattr(args, 'annotation', None):
        inputs['annotation'] = str(args.annotation)
    return inputs


def _make_cluster_lookup(clusters_df):
    try:
        from intervaltree import IntervalTree
        trees = {}
        for _, row in clusters_df.iterrows():
            key = (row['chrom'], row['strand'])
            trees.setdefault(key, IntervalTree())[row['start']:row['end'] + 1] = row['cluster_id']

        def lookup(chrom, strand, pos):
            hits = trees.get((chrom, strand), IntervalTree())[pos]
            return hits.pop().data if hits else None
        return lookup
    except ImportError:
        lists = {}
        for _, row in clusters_df.iterrows():
            lists.setdefault((row['chrom'], row['strand']), []).append(
                (row['start'], row['end'], row['cluster_id'])
            )

        def lookup(chrom, strand, pos):
            for start, end, cid in lists.get((chrom, strand), []):
                if start <= pos <= end:
                    return cid
            return None
        return lookup


def run_analyze(args: argparse.Namespace) -> int:
    """
    Run the full analysis pipeline.

    Args:
        args: Parsed command-line arguments

    Returns:
        Exit code (0 for success)
    """
    # Exactly one input channel is required, and manifest mode does not use the
    # positional. It is `nargs='?'` so `--manifest M` alone parses; validate the
    # combination here instead, where the message can name both options.
    # (Passing the manifest path as BOTH input and manifest stays valid — that is
    # what run/multi_sample.py:43 does, and existing callers rely on it.)
    if not getattr(args, 'manifest', None) and not getattr(args, 'input', None):
        # This runs BEFORE the module's logger is set up (it is created inside
        # run_analyze, further down), so use a local one rather than a name
        # that does not exist yet.
        import logging as _log_early
        _log_early.getLogger(__name__).error(
            "rectify analyze needs an input: either a positional TSV of corrected "
            "positions, or --manifest SAMPLES.tsv for multi-sample mode. Got neither."
        )
        return 2

    # Must be called before any numpy/pandas/sklearn/pydeseq2 import side-effects
    # so thread limits take effect before those libraries auto-spawn workers.
    from ...slurm import set_thread_limits
    set_thread_limits(getattr(args, 'threads', None))

    print("=" * 70)
    print("RECTIFY Analysis Pipeline")
    print("=" * 70)

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = output_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)
    tables_dir = output_dir / 'tables'
    tables_dir.mkdir(exist_ok=True)

    import sys as _sys_analyze
    import logging as _log_analyze
    from datetime import datetime as _dt_analyze, timezone as _tz_analyze
    from time import perf_counter as _perf_analyze
    from rectify.core.provenance import can_skip_stage
    from rectify.utils.version import get_rectify_git_sha as _get_sha_analyze

    _sample_id = output_dir.name
    _rectify_sha_analyze = _get_sha_analyze()
    _current_inputs = _collect_analyze_inputs(args)
    _analyze_logger = _log_analyze.getLogger(__name__)

    _skip_decision = can_skip_stage(
        stage='analyze', sample_output=output_dir, sample_id=_sample_id,
        current_inputs=_current_inputs, current_argv=_sys_analyze.argv,
        rectify_git_sha=_rectify_sha_analyze,
        force=getattr(args, 'force_all', False),
        force_stages=set(getattr(args, 'force_stage', '').split(','))
                   if getattr(args, 'force_stage', None) else None,
        accept_prior_provenance=getattr(args, 'accept_prior_provenance', False),
    )
    _analyze_logger.info("[RESUME] sample=%s stage=analyze decision=%s reason=%s",
        _sample_id, 'SKIP' if _skip_decision.skip else 'RUN', _skip_decision.reason)
    if _skip_decision.skip:
        _analyze_logger.info("[RESUME] Skipping analyze stage — prior sidecar valid.")
        return 0
    if getattr(args, 'dry_run_resume', False):
        print(f"[dry-run-resume] stage=analyze decision=RUN reason={_skip_decision.reason}")
        return 0
    _stage_started_at = _dt_analyze.now(_tz_analyze.utc).isoformat()
    _t_analyze_start = _perf_analyze()

    # Dispatch to manifest mode if --manifest is provided
    if getattr(args, 'manifest', None):
        return _run_analyze_manifest(args, output_dir, plots_dir, tables_dir,
                                     _sidecar_inputs=_current_inputs)

    # Initialize provenance tracking
    provenance = init_provenance(
        output_dir,
        description="RECTIFY analysis output (clustering, DESeq2, motifs, etc.)",
        config=vars(args)
    )

    # Load corrected positions
    print(f"\n[1/9] Loading corrected positions from {args.input}...")
    positions_df = load_corrected_positions(args.input, args.sample_column)
    print(f"  Loaded {len(positions_df):,} positions from {positions_df[args.sample_column].nunique()} samples")

    # Load annotation early to detect exclusion regions
    annotation_df = None
    if args.annotation:
        annotation_df = load_annotation(args.annotation)

    # Filter problematic regions (default: exclude mito and rDNA separately)
    exclude_mito = args.exclude_mito and not getattr(args, 'include_mito', False)
    exclude_rdna = getattr(args, 'exclude_rdna', True) and not getattr(args, 'include_rdna', False)

    if exclude_mito or exclude_rdna:
        n_before = len(positions_df)

        # Auto-detect mitochondrial and rDNA from annotation
        mito_chroms, rdna_regions = _detect_exclusion_regions(
            annotation_df, positions_df['chrom'].unique()
        )

        # Report what was auto-detected
        if mito_chroms:
            print(f"  Auto-detected mitochondrial chromosomes: {', '.join(sorted(mito_chroms))}")
        if rdna_regions:
            for chrom, start, end in rdna_regions:
                print(f"  Auto-detected rDNA locus: {chrom}:{start:,}-{end:,} ({(end-start):,} bp)")

        # Mitochondrial filter
        mito_mask = pd.Series(False, index=positions_df.index)
        if exclude_mito and mito_chroms:
            mito_mask = positions_df['chrom'].isin(mito_chroms)

        # rDNA filter (may have multiple regions)
        rdna_mask = pd.Series(False, index=positions_df.index)
        if exclude_rdna and rdna_regions:
            for chrom, start, end in rdna_regions:
                region_mask = (
                    (positions_df['chrom'] == chrom) &
                    (positions_df['corrected_position'] >= start) &
                    (positions_df['corrected_position'] < end)
                )
                rdna_mask = rdna_mask | region_mask

        # Apply combined filter
        exclude_mask = mito_mask | rdna_mask
        positions_df = positions_df[~exclude_mask]

        n_mito = mito_mask.sum()
        n_rdna = rdna_mask.sum()
        n_total_removed = n_before - len(positions_df)

        if n_total_removed > 0:
            print(f"  Excluded {n_total_removed:,} positions:")
            if n_mito > 0:
                print(f"    - Mitochondrial: {n_mito:,} positions")
            if n_rdna > 0:
                print(f"    - rDNA locus: {n_rdna:,} positions")

    # Determine count column (may be auto-set by chunked loading)
    count_col = args.count_column
    if 'count' in positions_df.columns and not count_col:
        count_col = 'count'
        print(f"  Using pre-aggregated counts")
    elif 'fraction' in positions_df.columns and not count_col:
        count_col = 'fraction'
        print(f"  Using fractional assignment weights")

    # Generate bedgraphs (default: enabled, unless --no-bedgraph)
    if not getattr(args, 'no_bedgraph', False):
        bedgraph_dir = Path(args.bedgraph_dir) if args.bedgraph_dir else output_dir / 'bedgraph'
        print(f"\n[Bedgraph] Generating strand-specific bedgraph files...")
        generate_bedgraphs(
            positions_df,
            bedgraph_dir,
            sample_column=args.sample_column,
            position_column='corrected_3prime' if 'corrected_3prime' in positions_df.columns else 'corrected_position',
            normalize_rpm=True,
        )
        print(f"  Saved to {bedgraph_dir}")

    # Genomic distribution analysis (default: enabled if annotation provided)
    if annotation_df is not None and not getattr(args, 'no_genomic_distribution', False):
        position_col = 'corrected_3prime' if 'corrected_3prime' in positions_df.columns else 'corrected_position'

        # 3' end distribution
        print(f"\n[Genomic Distribution] Analyzing 3' end distribution by genomic region...")
        try:
            dist3_files = run_3prime_distribution_analysis(
                positions_df,
                annotation_df,
                output_dir=str(plots_dir),
                sample_column=args.sample_column,
                position_column=position_col,
                count_column=count_col,
            )
            print(f"  Generated {len(dist3_files)} plots/tables")
        except ImportError as e:
            print(f"  Warning: Skipping 3' end distribution (missing dependency: {e})")
        except Exception as e:
            print(f"  Warning: 3' end distribution analysis failed: {e}")

        # Transcript body distribution
        print(f"\n[Genomic Distribution] Analyzing transcript body distribution by RNA biotype...")
        try:
            body_files = run_transcript_body_distribution_analysis(
                positions_df,
                annotation_df,
                output_dir=str(plots_dir),
                sample_column=args.sample_column,
                start_column='alignment_start',
                end_column='alignment_end',
            )
            print(f"  Generated {len(body_files)} plots/tables")
        except ImportError as e:
            print(f"  Warning: Skipping body distribution (missing dependency: {e})")
        except Exception as e:
            print(f"  Warning: Transcript body distribution analysis failed: {e}")

        # 5' end distribution
        print(f"\n[Genomic Distribution] Analyzing 5' end distribution by genomic region...")
        try:
            dist5_files = run_5prime_distribution_analysis(
                positions_df,
                annotation_df,
                output_dir=str(plots_dir),
                sample_column=args.sample_column,
                position_column='five_prime_position',
                count_column=count_col,
            )
            print(f"  Generated {len(dist5_files)} plots/tables")
        except ImportError as e:
            print(f"  Warning: Skipping 5' end distribution (missing dependency: {e})")
        except Exception as e:
            print(f"  Warning: 5' end distribution analysis failed: {e}")

    # Form clusters
    fraction_col = (
        'fraction'
        if 'fraction' in positions_df.columns and count_col != 'fraction'
        else None
    )
    _max_radius = getattr(args, 'max_cluster_radius', 10)
    _min_peak_sep = getattr(args, 'min_peak_sep', 5)
    _min_cluster_samples = getattr(args, 'min_cluster_samples', 2)

    print(f"\n[2/9] Forming CPA clusters (adaptive, radius={_max_radius}bp, "
          f"min_peak_sep={_min_peak_sep}bp, min_samples={_min_cluster_samples})...")
    clusters_df = cluster_cpa_sites_adaptive(
        positions_df,
        max_cluster_radius=_max_radius,
        min_peak_separation=_min_peak_sep,
        min_reads=args.min_reads,
        count_col=count_col,
    )
    print(f"  Formed {len(clusters_df):,} clusters")

    # Annotate with genes (annotation_df already loaded above if provided)
    if args.annotation and annotation_df is not None:
        print(f"\n  Annotating clusters with genes...")
        clusters_df = annotate_clusters_with_genes(clusters_df, annotation_df)
        n_annotated = clusters_df['gene_id'].notna().sum()
        print(f"  Annotated {n_annotated:,} clusters ({100*n_annotated/len(clusters_df):.1f}%)")

    clusters_path = output_dir / 'cpa_clusters.tsv'

    # 5' end (TSS) clustering — mirrors CPA clustering but on five_prime_position
    # with a wider window (75 bp) to account for DRS 5'-end noisiness.
    tss_distance = getattr(args, 'tss_cluster_distance', 75)
    if 'five_prime_position' in positions_df.columns:
        print(f"\n[2b/9] Forming TSS clusters (distance={tss_distance}bp)...")
        try:
            tss_positions_df = positions_df.dropna(subset=['five_prime_position']).copy()
            tss_positions_df['five_prime_position'] = tss_positions_df['five_prime_position'].astype(int)
            tss_clusters_df = cluster_cpa_sites(
                tss_positions_df,
                cluster_distance=tss_distance,
                min_reads=args.min_reads,
                position_col='five_prime_position',
                count_col=count_col,
            )
            if args.annotation and annotation_df is not None:
                tss_clusters_df = annotate_clusters_with_genes(tss_clusters_df, annotation_df)
            tss_clusters_path = output_dir / 'tss_clusters.tsv'
            tss_clusters_df.to_csv(tss_clusters_path, sep='\t', index=False)
            print(f"  Formed {len(tss_clusters_df):,} TSS clusters → {tss_clusters_path}")

            # TSS count matrix
            tss_count_matrix = build_cluster_count_matrix(
                tss_positions_df,
                tss_clusters_df,
                sample_col=args.sample_column,
                count_col=count_col,
                fraction_col=fraction_col,
                position_col='five_prime_position',
            )
            tss_counts_path = output_dir / 'tss_cluster_counts.tsv'
            tss_count_matrix.to_csv(tss_counts_path, sep='\t')
            print(f"  TSS count matrix: {tss_count_matrix.shape[0]:,} clusters × {tss_count_matrix.shape[1]} samples → {tss_counts_path}")
        except Exception as _tss_exc:
            print(f"  Warning: TSS clustering failed: {_tss_exc}")

    # Build count matrix
    print(f"\n[3/9] Building cluster count matrix...")
    count_matrix = build_cluster_count_matrix(
        positions_df,
        clusters_df,
        sample_col=args.sample_column,
        count_col=count_col,
        fraction_col=fraction_col,
    )
    if _min_cluster_samples > 1 and not count_matrix.empty:
        _sample_counts = (count_matrix > 0).sum(axis=1)
        _keep_cluster_ids = set(_sample_counts[_sample_counts >= _min_cluster_samples].index)
        _n_before = len(clusters_df)
        clusters_df = clusters_df[clusters_df['cluster_id'].isin(_keep_cluster_ids)].copy()
        count_matrix = count_matrix.loc[
            [cid for cid in count_matrix.index if cid in _keep_cluster_ids]
        ]
        print(f"  Kept {len(clusters_df):,}/{_n_before:,} clusters present in >= {_min_cluster_samples} samples")
    print(f"  Matrix shape: {count_matrix.shape[0]:,} clusters × {count_matrix.shape[1]} samples")

    sample_names = count_matrix.columns.tolist()
    lookup = _make_cluster_lookup(clusters_df)
    _chrom_format = getattr(args, 'chrom_format', 'passthrough')
    # Build the transcript model once (used by both containment attribution and
    # per-site region classification). None on absent/degraded annotation.
    _transcript_model = None
    if annotation_df is not None:
        from ..analyze.transcript_model import build_transcript_model_for_analyze
        _transcript_model = build_transcript_model_for_analyze(args, annotation_df, _chrom_format)
    cluster_gene_attributions = None
    try:
        _samples_for_attr = [{'sample_id': 'input', 'path': args.input}]
        clusters_df, cluster_gene_attributions = _build_cluster_gene_attribution_table(
            args,
            samples=_samples_for_attr,
            clusters_df=clusters_df,
            lookup=lookup,
            annotation_df=annotation_df,
            chrom_format=_chrom_format,
            model=_transcript_model,
        )
        attr_path = output_dir / 'cluster_gene_attributions.tsv'
        cluster_gene_attributions.to_csv(attr_path, sep='\t', index=False)
        n_attr_clusters = cluster_gene_attributions['cluster_id'].nunique() if not cluster_gene_attributions.empty else 0
        n_attr_genes = cluster_gene_attributions['gene_id'].nunique() if not cluster_gene_attributions.empty else 0
        print(
            f"  Cluster gene attributions: {len(cluster_gene_attributions):,} rows, "
            f"{n_attr_clusters:,} clusters, {n_attr_genes:,} genes → {attr_path}"
        )
    except Exception as _attr_exc:
        mode = getattr(args, 'gene_attribution_mode', 'containment')
        # containment (the default) degrades gracefully to the legacy gene_id
        # already on clusters_df, so the default flip never turns a working
        # --annotation run into a hard failure. body/reference stay hard-fail.
        _soft_modes = {'annotation', 'none', 'legacy-3prime-window',
                       'containment', 'containment-then-annotation'}
        if mode not in _soft_modes:
            print(f"ERROR: cluster gene attribution failed for mode '{mode}': {_attr_exc}", flush=True)
            return 1
        print(f"  WARNING: cluster gene attribution failed for mode '{mode}': {_attr_exc}")
        cluster_gene_attributions = None

    # Per-site region classification (region_class + continuous transcript
    # coordinates + antisense/also_within/intron evidence). Adds columns ONLY;
    # never touches gene_id. Runs regardless of attribution mode; never fatal.
    if _transcript_model is not None:
        try:
            from ..analyze.transcript_model import annotate_clusters_with_transcript_model
            from ...utils.chromosome import normalize_chromosome
            _norm = lambda x: normalize_chromosome(x, _chrom_format)
            clusters_df = annotate_clusters_with_transcript_model(
                clusters_df, _transcript_model,
                samples=[{'sample_id': 'input', 'path': args.input}],
                cluster_lookup=lookup, chrom_normalizer=_norm,
            )
            _n_rc = int(clusters_df['region_class'].notna().sum()) if 'region_class' in clusters_df.columns else 0
            print(f"  Region classification: {_n_rc:,} clusters classified")
        except Exception as _rc_exc:
            print(f"  WARNING: region classification failed: {_rc_exc}")

    # Save clusters after cross-sample filters are applied.
    clusters_df.to_csv(clusters_path, sep='\t', index=False)
    print(f"  Saved clusters to {clusters_path}")

    # Save count matrix
    counts_path = output_dir / 'cluster_counts.tsv'
    count_matrix.to_csv(counts_path, sep='\t')

    # Create sample metadata
    if args.reference:
        reference_condition = args.reference.lower()
    else:
        # Auto-detect reference
        control_samples = detect_control_samples(sample_names)
        if control_samples:
            reference_condition = extract_condition_from_sample(control_samples[0])
            print(f"  Auto-detected reference condition: {reference_condition}")
        else:
            print("  Warning: Could not auto-detect reference. Using first condition.")
            reference_condition = extract_condition_from_sample(sample_names[0])

    sample_metadata = create_sample_metadata(sample_names, control_samples if not args.reference else None)
    sample_metadata.to_csv(output_dir / 'sample_metadata.tsv', sep='\t')

    # Run PCA
    print(f"\n[4/9] Running PCA analysis...")

    # Parse sample sets if provided
    sample_sets = None
    if args.sample_sets:
        import json
        if args.sample_sets.startswith('{'):
            sample_sets = json.loads(args.sample_sets)
        elif Path(args.sample_sets).exists():
            with open(args.sample_sets) as f:
                sample_sets = json.load(f)
        else:
            print(f"  Warning: Could not parse sample_sets: {args.sample_sets}")

    # If sample sets defined, create separate PCA plots
    if sample_sets:
        for set_name, conditions in sample_sets.items():
            # Filter to samples in this set
            set_samples = sample_metadata[
                sample_metadata['condition'].isin(conditions)
            ].index.tolist()
            set_samples = [s for s in set_samples if s in count_matrix.columns]

            if len(set_samples) < 3:
                print(f"  Skipping PCA for {set_name} (only {len(set_samples)} samples)")
                continue

            set_matrix = count_matrix[set_samples]
            set_meta = sample_metadata.loc[sample_metadata.index.isin(set_samples)]

            pca_results = run_pca_analysis(set_matrix)
            if pca_results['pca_coords'] is not None and not pca_results['pca_coords'].empty:
                pca_path = plots_dir / f'pca_{set_name}.png'
                plot_pca(
                    pca_results,
                    sample_metadata=set_meta,
                    color_by='condition',
                    output_path=str(pca_path),
                    title=f'PCA: {set_name} ({", ".join(conditions)})',
                )
                print(f"  Saved PCA plot for {set_name}: {pca_path}")
    else:
        # Default: single PCA with all samples
        pca_results = run_pca_analysis(count_matrix)
        if pca_results['pca_coords'] is not None and not pca_results['pca_coords'].empty:
            pca_path = plots_dir / 'pca_samples.png'
            plot_pca(
                pca_results,
                sample_metadata=sample_metadata,
                color_by='condition',
                output_path=str(pca_path),
                title='Sample PCA (CPA Clusters)',
            )
            print(f"  Saved PCA plot to {pca_path}")

    # Sample heatmap
    print(f"\n[5/9] Creating sample clustering heatmap...")
    if count_matrix.shape[1] < 2:
        print(f"  Skipping heatmap — requires ≥2 samples (have {count_matrix.shape[1]})")
    else:
        heatmap_path = plots_dir / 'sample_heatmap.png'
        plot_sample_heatmap(
            count_matrix,
            sample_metadata=sample_metadata,
            color_by='condition',
            output_path=str(heatmap_path),
        )
        print(f"  Saved heatmap to {heatmap_path}")

    # Run DESeq2
    deseq2_gene_results = {}
    deseq2_cluster_results = {}

    if args.run_deseq2:
        print(f"\n[6/9] Running DESeq2 differential expression...")

        # Guard: DESeq2 requires at least 2 distinct conditions
        _n_conditions = sample_metadata['condition'].nunique()
        if _n_conditions < 2:
            print(f"  Warning: DESeq2 requires ≥2 conditions but only {_n_conditions} found. Skipping.")
        else:
            # Gene-level
            print("  Gene-level analysis...")
            deseq2_gene_results = run_deseq2_gene_level(
                count_matrix,
                clusters_df,
                sample_metadata,
                reference_condition,
                cluster_gene_attributions=cluster_gene_attributions,
                n_cpus=args.threads,
            )
            for condition, result_df in deseq2_gene_results.items():
                result_path = tables_dir / f'deseq2_genes_{condition}.tsv'
                result_df.to_csv(result_path, sep='\t')
                n_sig = (result_df['padj'] < 0.05).sum()
                print(f"    {condition}: {n_sig:,} significant genes")

            # Free memory between gene and cluster analyses
            import gc
            gc.collect()

            # Cluster-level
            print("  Cluster-level analysis...")
            deseq2_cluster_results = run_deseq2_cluster_level(
                count_matrix,
                clusters_df,
                sample_metadata,
                reference_condition,
                n_cpus=args.threads,
            )
            for condition, result_df in deseq2_cluster_results.items():
                result_path = tables_dir / f'deseq2_clusters_{condition}.tsv'
                result_df.to_csv(result_path, sep='\t')
                n_sig = (result_df['padj'] < 0.05).sum()
                print(f"    {condition}: {n_sig:,} significant clusters")

            # Free memory after DESeq2
            gc.collect()
    else:
        print(f"\n[6/9] Skipping DESeq2 (use --run-deseq2 to enable)")

    # GO enrichment
    if args.go_annotations and deseq2_gene_results:
        print(f"\n[7/9] Running GO enrichment analysis...")
        from ..analyze.go_enrichment import load_go_annotations

        go_annotations = load_go_annotations(
            args.go_annotations,
            gene_col='gene_name',
            go_term_col='description',
            category_col='go_category',
        )

        # Build systematic → common name mapping so DESeq2 indices match GO file
        _sys2common = {}
        if 'gene_id' in clusters_df.columns and 'gene_name' in clusters_df.columns:
            for _, _r in clusters_df[['gene_id', 'gene_name']].dropna().drop_duplicates().iterrows():
                _sys2common[_r['gene_id']] = _r['gene_name']
        if cluster_gene_attributions is not None and not cluster_gene_attributions.empty:
            for _, _r in cluster_gene_attributions[['gene_id', 'gene_name']].dropna().drop_duplicates().iterrows():
                _sys2common[_r['gene_id']] = _r['gene_name']

        def _to_common_names(genes):
            return [_sys2common.get(g, g) for g in genes]

        for condition, result_df in deseq2_gene_results.items():
            # Upregulated genes
            up_genes = _to_common_names(result_df[
                (result_df['padj'] < 0.05) & (result_df['log2FoldChange'] > 1)
            ].index.tolist())

            if len(up_genes) >= 10:
                go_up = run_go_enrichment(up_genes, go_annotations)
                if not go_up.empty:
                    go_up.to_csv(tables_dir / f'go_enrichment_up_{condition}.tsv', sep='\t', index=False)
                    plot_go_enrichment(
                        go_up,
                        output_path=str(plots_dir / f'go_enrichment_up_{condition}.png'),
                        title=f'GO Enrichment: Upregulated in {condition}',
                    )

            # Downregulated genes
            down_genes = _to_common_names(result_df[
                (result_df['padj'] < 0.05) & (result_df['log2FoldChange'] < -1)
            ].index.tolist())

            if len(down_genes) >= 10:
                go_down = run_go_enrichment(down_genes, go_annotations)
                if not go_down.empty:
                    go_down.to_csv(tables_dir / f'go_enrichment_down_{condition}.tsv', sep='\t', index=False)
                    plot_go_enrichment(
                        go_down,
                        output_path=str(plots_dir / f'go_enrichment_down_{condition}.png'),
                        title=f'GO Enrichment: Downregulated in {condition}',
                    )
    else:
        print(f"\n[7/9] Skipping GO enrichment (provide --go-annotations)")

    # Motif discovery
    if args.genome and args.run_motif and deseq2_cluster_results:
        print(f"\n[8/9] Running de novo motif discovery...")

        for condition, result_df in deseq2_cluster_results.items():
            enriched = result_df[
                (result_df['padj'] < 0.05) & (result_df['log2FoldChange'] > 1)
            ]
            depleted = result_df[
                (result_df['padj'] < 0.05) & (result_df['log2FoldChange'] < -1)
            ]

            if len(enriched) >= 20 and len(depleted) >= 20:
                motif_results = run_differential_motif_analysis(
                    enriched,
                    depleted,
                    args.genome,
                    str(output_dir / 'motifs' / condition),
                    upstream_window=args.motif_upstream,
                    downstream_window=args.motif_downstream,
                )

                summary = summarize_motif_results(motif_results)
                summary.to_csv(tables_dir / f'motif_summary_{condition}.tsv', sep='\t', index=False)
                print(f"  {condition}: {summary['motifs_found'].sum()} motifs found")
    else:
        print(f"\n[8/9] Skipping motif discovery (provide --genome and --run-motif)")

    # Shift analysis
    shift_results = None
    if len(sample_metadata['condition'].unique()) >= 2:
        print(f"\n[9/9] Running cluster shift analysis...")

        conditions = [c for c in sample_metadata['condition'].unique() if c != reference_condition]

        for condition in conditions:
            shift_df = analyze_cluster_shifts(
                count_matrix,
                clusters_df,
                reference_condition,
                condition,
                sample_metadata,
                cluster_gene_attributions=cluster_gene_attributions,
            )

            if not shift_df.empty:
                shift_df.to_csv(tables_dir / f'shift_analysis_{condition}.tsv', sep='\t', index=False)

                # Plot summary
                plot_shift_summary(
                    shift_df,
                    output_path=str(plots_dir / f'shift_summary_{condition}.png'),
                )

                # Top shifted genes browser plots
                top_shifted = get_top_shifted_genes(shift_df, n_top=20)
                for _, row in top_shifted.head(5).iterrows():
                    plot_gene_browser(
                        row['gene_id'],
                        count_matrix,
                        clusters_df,
                        sample_metadata,
                        [reference_condition, condition],
                        output_path=str(plots_dir / f'browser_{row["gene_name"]}_{condition}.png'),
                    )

                print(f"  {condition}: {len(shift_df)} genes analyzed, "
                      f"{(shift_df['distribution_divergence'] > 0.2).sum()} with large shifts")

                shift_results = shift_df  # Keep last for summary
    else:
        print(f"\n[9/9] Skipping shift analysis (need >=2 conditions)")

    # Generate summary
    print(f"\n[Summary] Generating report...")

    if cluster_gene_attributions is not None and not cluster_gene_attributions.empty:
        n_genes = cluster_gene_attributions['gene_id'].nunique()
    else:
        n_genes = clusters_df['gene_id'].nunique() if 'gene_id' in clusters_df.columns else 0

    summary_df = generate_analysis_summary(
        n_samples=len(sample_names),
        n_clusters=len(clusters_df),
        n_genes=n_genes,
        deseq2_gene_results=deseq2_gene_results,
        deseq2_cluster_results=deseq2_cluster_results,
        reference_condition=reference_condition,
    )

    summary_df.to_csv(output_dir / 'analysis_summary.tsv', sep='\t', index=False)

    # Generate HTML report
    plots_dict = {p.stem: str(p) for p in plots_dir.glob('*.png')}
    tables_dict = {t.stem: str(t) for t in tables_dir.glob('*.tsv')}

    html_path = output_dir / 'report.html'
    generate_html_report(
        summary_df,
        plots_dict,
        tables_dict,
        str(html_path),
        title='RECTIFY Analysis Report',
    )

    # Save provenance
    # Record all output files
    for tsv_file in tables_dir.glob('*.tsv'):
        provenance.add_output_file(tsv_file, source_files=[Path(args.input)])
    for png_file in plots_dir.glob('*.png'):
        provenance.add_output_file(png_file)
    provenance.add_output_file(html_path)
    provenance.add_output_file(output_dir / 'analysis_summary.tsv')
    provenance.save()

    print(f"\n" + "=" * 70)
    print(f"Analysis complete!")
    print(f"  Output directory: {output_dir}")
    print(f"  HTML report: {html_path}")
    print(f"  Provenance: {output_dir / 'PROVENANCE.json'}")
    print("=" * 70)

    _analyze_wall_secs = _perf_analyze() - _t_analyze_start
    try:
        from rectify.core.provenance import ProvenanceRecord, write_stage_sidecar
        _sc_outputs: dict = {}
        if html_path.exists():
            _sc_outputs['html_report'] = str(html_path)
        _summary_tsv = output_dir / 'analysis_summary.tsv'
        if _summary_tsv.exists():
            _sc_outputs['summary'] = str(_summary_tsv)
        _sc_record = ProvenanceRecord.from_components(
            stage='analyze', stage_subtype=None, sample_id=_sample_id,
            sample_output_dir=output_dir, started_at=_stage_started_at,
            completed_at=_dt_analyze.now(_tz_analyze.utc).isoformat(),
            exit_status=0, inputs=_current_inputs, outputs=_sc_outputs,
            stats={'wall_seconds': _analyze_wall_secs},
            skip_check_config={'ignore_argv': ['--threads', '--verbose']},
            argv=_sys_analyze.argv, rectify_git_sha=_rectify_sha_analyze,
        )
        write_stage_sidecar(_sc_record, sample_output=output_dir)
        _analyze_logger.info("[PROVENANCE] Wrote analyze sidecar for %s", _sample_id)
    except Exception as _sc_exc:
        _analyze_logger.warning("[PROVENANCE] Failed to write analyze sidecar (non-fatal): %s", _sc_exc)

    return 0


def create_analyze_parser(subparsers) -> argparse.ArgumentParser:
    """Create argument parser for analyze command."""
    parser = subparsers.add_parser(
        'analyze',
        help='Analyze corrected 3\' end positions',
        description='Perform downstream analysis including clustering, '
                    'differential expression, PCA, GO enrichment, and motif discovery.',
    )

    # Required arguments
    parser.add_argument(
        'input',
        nargs='?',
        default=None,
        help='Input TSV with corrected positions. OMIT when using --manifest: '
             'in manifest mode this value is never read, but argparse made it '
             'required, so a correct manifest-mode invocation died on '
             '"the following arguments are required: input" (reported by the '
             '668-drs-arm session, 4 minutes into a cohort run).',
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory for results',
    )

    # Reference genome and annotation
    parser.add_argument(
        '--genome',
        help='Reference genome FASTA (required for motif discovery)',
    )

    parser.add_argument(
        '--annotation',
        help='Gene annotation file (GTF/GFF or TSV)',
    )

    parser.add_argument(
        '--ncrna-atlas',
        help='Named supplementary ncRNA atlas (CUT/SUT/XUT ...) from the bundled registry '
             '(rectify/data/ncrna_atlases/atlases.yaml), layered on top of --annotation for '
             'the transcript-model classifier. Default: SGD-core only.',
    )

    parser.add_argument(
        '--ncrna-annotations',
        nargs='+',
        metavar='FILE:SOURCE:CLASS',
        help='Ad-hoc supplementary ncRNA track(s) as file.gff:source:class '
             '(e.g. mycuts.gff:MyStudy2024:CUT). Force-tags every feature with the given '
             'source/class regardless of the file column-3 type. Combines with --ncrna-atlas.',
    )

    parser.add_argument(
        '--gene-attribution-window',
        type=int,
        default=100,
        metavar='BP',
        help='For --gene-attribution-mode containment, the downstream window (bp) past a gene '
             '3prime end within which an uncontained cluster is classified downstream_readthrough '
             '(kept equal to the legacy proximity fallback so gene_id and region_class stay '
             'coherent). Default: 100.',
    )

    parser.add_argument(
        '--utr3-proximal-distal-split',
        type=int,
        default=150,
        metavar='BP',
        help='Boundary (nt past the stop codon) between region_class 3primeUTR_proximal and '
             '3primeUTR_distal for the transcript-model classifier. Default: 150.',
    )

    from rectify.data import add_organism_args
    add_organism_args(parser)

    # Sample information
    parser.add_argument(
        '--sample-column',
        default='sample',
        help='Column name for sample identifier (default: sample)',
    )

    parser.add_argument(
        '--count-column',
        help='Column name for read counts (optional)',
    )

    parser.add_argument(
        '--reference',
        help='Reference condition name (auto-detected if not specified)',
    )

    # Clustering parameters
    parser.add_argument(
        '--cluster-distance',
        type=int,
        default=25,
        help='Maximum distance (bp) to merge positions into clusters (default: 25)',
    )

    parser.add_argument(
        '--min-reads',
        type=int,
        default=5,
        help='Minimum reads per cluster (default: 5)',
    )

    parser.add_argument(
        '--max-cluster-radius',
        type=int,
        default=10,
        metavar='BP',
        help='Adaptive clustering: maximum radius (bp) from peak to cluster boundary (default: 10). '
             'Controls how far each CPA peak extends when using valley-based adaptive clustering.',
    )

    parser.add_argument(
        '--min-peak-sep',
        type=int,
        default=5,
        metavar='BP',
        help='Adaptive clustering: minimum separation (bp) between distinct CPA peaks (default: 5). '
             'Peaks closer than this are merged into the dominant peak.',
    )

    parser.add_argument(
        '--min-cluster-samples',
        type=int,
        default=2,
        metavar='N',
        help='Minimum number of samples a cluster must appear in to be retained (default: 2). '
             'Filters singleton-sample clusters that are likely noise.',
    )

    # Analysis options
    parser.add_argument(
        '--run-deseq2',
        action='store_true',
        help='Run DESeq2 differential expression analysis',
    )

    parser.add_argument(
        '--go-annotations',
        help='GO annotation file for enrichment analysis',
    )

    parser.add_argument(
        '--run-motif',
        action='store_true',
        help='Run de novo motif discovery (requires --genome)',
    )

    parser.add_argument(
        '--motif-upstream',
        type=int,
        default=100,
        help='Window upstream of CPA for motif discovery (default: 100bp)',
    )

    parser.add_argument(
        '--motif-downstream',
        type=int,
        default=50,
        help='Window downstream of CPA for motif discovery (default: 50bp)',
    )

    parser.add_argument(
        '--gene-attribution-mode',
        choices=[
            'containment',
            'containment-then-annotation',
            'legacy-3prime-window',
            'annotation',
            'none',
            'body',
            'body-then-annotation',
            'reference',
            'reference-then-annotation',
        ],
        default='containment',
        help='How CPA clusters are assigned to genes for gene-level DESeq2 and shift analysis. '
             'containment (default) attributes each cluster to the gene whose transcribed body '
             'contains its modal position, falling back to the nearest-3prime-end window for '
             'clusters not contained by any gene body (167 Part A). legacy-3prime-window (alias: '
             'annotation) keeps the pure nearest-TES rule. body uses read-body overlap from '
             'corrected TSV alignment spans. reference maps external long-read per-position '
             'attribution TSVs onto current clusters. *-then-annotation fills missing clusters '
             'with the legacy annotation fallback.',
    )

    parser.add_argument(
        '--gene-attributions',
        nargs='+',
        help='External per-position gene attribution TSV file(s) or directories for '
             '--gene-attribution-mode reference/reference-then-annotation. Expected columns: '
             'chrom, position, strand, gene_id, attributed_count.',
    )

    parser.add_argument(
        '--gene-attribution-reference-window',
        type=int,
        default=0,
        metavar='BP',
        help='For --gene-attribution-mode reference/reference-then-annotation, allow a '
             'long-read attributed position to assign to the nearest same-strand current '
             'CPA cluster within BP bases when it does not overlap a cluster exactly. '
             'Exact overlaps are always preferred. Default: 0.',
    )

    parser.add_argument(
        '--body-gene-attribution-mode',
        choices=['origin', 'proportional'],
        default='origin',
        help='Read-body attribution rule for --gene-attribution-mode body/body-then-annotation. '
             'origin assigns each read to the transcription-upstream overlapped gene; '
             'proportional splits by feature overlap length.',
    )

    parser.add_argument(
        '--include-intergenic-attributions',
        action='store_true',
        help='Keep intergenic rows from external attribution inputs. By default they are dropped.',
    )

    # Performance
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads for DESeq2 (default: 4)',
    )

    # Filtering options
    parser.add_argument(
        '--exclude-mito',
        action='store_true',
        default=True,
        help='Exclude mitochondrial chromosome from analysis (default: True)',
    )

    parser.add_argument(
        '--include-mito',
        action='store_true',
        help='Include mitochondrial chromosome (overrides --exclude-mito)',
    )

    parser.add_argument(
        '--exclude-rdna',
        action='store_true',
        default=True,
        help='Exclude rDNA locus from analysis (default: True). Auto-detected from GFF if available.',
    )

    parser.add_argument(
        '--include-rdna',
        action='store_true',
        help='Include rDNA locus (overrides --exclude-rdna)',
    )

    # Bedgraph output (default: enabled)
    parser.add_argument(
        '--no-bedgraph',
        action='store_true',
        help='Disable bedgraph output (bedgraph is generated by default)',
    )

    parser.add_argument(
        '--bedgraph-dir',
        help='Output directory for bedgraph files (default: output_dir/bedgraph)',
    )

    # Genomic distribution analysis (default: enabled if annotation provided)
    parser.add_argument(
        '--no-genomic-distribution',
        action='store_true',
        help="Disable genomic distribution pie chart (shows 3' end distribution across UTR3/CDS/UTR5/intergenic)",
    )

    # Sample set grouping for separate PCA plots
    parser.add_argument(
        '--sample-sets',
        type=str,
        help='JSON string or file defining sample sets for separate PCA plots. '
             'Format: {"set1_name": ["condition1", "condition2"], "set2_name": ["condition3"]}',
    )

    # Manifest mode (memory-efficient multi-sample processing)
    parser.add_argument(
        '--manifest',
        help='Sample manifest TSV (columns: sample_id, path, [condition]). '
             'When provided, processes per-sample TSVs one at a time instead of '
             'a pre-combined TSV.',
    )

    # Chromosome-name normalization (kept off by default for species portability)
    parser.add_argument(
        '--chrom-format',
        choices=['passthrough', 'ucsc', 'ncbi', 'sgd'],
        default='passthrough',
        help='Target chromosome-name convention for output tables. '
             '`passthrough` (default): preserve input names as-is — works for any species. '
             '`ucsc`: chrI / chr1. `ncbi`: ref|NC_001133| / NC_000001.11. `sgd`: BK006935.2. '
             'NOTE: ucsc/ncbi/sgd conversions currently only have lookup tables for S. cerevisiae; '
             'on other species these modes pass through with a warning.',
    )

    from rectify.core.provenance.cli import add_resume_args
    add_resume_args(parser)

    parser.set_defaults(func=run_analyze)

    return parser
