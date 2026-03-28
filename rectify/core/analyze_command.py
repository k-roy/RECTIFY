#!/usr/bin/env python3
"""
RECTIFY Analyze Command

CLI entry point for downstream analysis of corrected 3' end positions.

Author: Kevin R. Roy
Date: 2026-03-17
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd

from .analyze import (
    cluster_cpa_sites,
    build_cluster_count_matrix,
    run_deseq2_gene_level,
    run_deseq2_cluster_level,
    detect_control_samples,
    create_sample_metadata,
    run_pca_analysis,
    plot_pca,
    plot_sample_heatmap,
    plot_cluster_heatmap,
    run_go_enrichment,
    plot_go_enrichment,
    run_differential_motif_analysis,
    summarize_motif_results,
    analyze_cluster_shifts,
    get_top_shifted_genes,
    plot_gene_browser,
    plot_shift_summary,
    generate_summary_tables,
    generate_analysis_summary,
    generate_html_report,
    run_genomic_distribution_analysis,
)
from .analyze.clustering import annotate_clusters_with_genes
from .analyze.deseq2 import extract_condition_from_sample
from ..utils.provenance import init_provenance


def run_analyze(args: argparse.Namespace) -> int:
    """
    Run the full analysis pipeline.

    Args:
        args: Parsed command-line arguments

    Returns:
        Exit code (0 for success)
    """
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
                    (positions_df['corrected_position'] <= end)
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

    # Generate bedgraphs (default: enabled, unless --no-bedgraph)
    if not getattr(args, 'no_bedgraph', False):
        bedgraph_dir = Path(args.bedgraph_dir) if args.bedgraph_dir else output_dir / 'bedgraph'
        print(f"\n[Bedgraph] Generating strand-specific bedgraph files...")
        generate_bedgraphs(
            positions_df,
            bedgraph_dir,
            sample_column=args.sample_column,
            position_column='corrected_3prime' if 'corrected_3prime' in positions_df.columns else 'position',
            normalize_rpm=True,
        )
        print(f"  Saved to {bedgraph_dir}")

    # Genomic distribution analysis (default: enabled if annotation provided)
    if annotation_df is not None and not getattr(args, 'no_genomic_distribution', False):
        print(f"\n[Genomic Distribution] Analyzing 3' end distribution by genomic region...")
        position_col = 'corrected_3prime' if 'corrected_3prime' in positions_df.columns else 'corrected_position'
        try:
            genomic_dist_files = run_genomic_distribution_analysis(
                positions_df,
                annotation_df,
                output_dir=str(plots_dir),
                sample_column=args.sample_column,
                position_column=position_col,
                count_column=count_col,
            )
            print(f"  Generated {len(genomic_dist_files)} plots/tables")
        except ImportError as e:
            print(f"  Warning: Skipping genomic distribution (missing dependency: {e})")
        except Exception as e:
            print(f"  Warning: Genomic distribution analysis failed: {e}")

    # Form clusters
    fraction_col = 'fraction' if 'fraction' in positions_df.columns else None
    print(f"\n[2/9] Forming CPA clusters (distance={args.cluster_distance}bp)...")
    clusters_df = cluster_cpa_sites(
        positions_df,
        cluster_distance=args.cluster_distance,
        min_reads=args.min_reads,
        count_col=count_col,
        fraction_col=fraction_col,
    )
    print(f"  Formed {len(clusters_df):,} clusters")

    # Annotate with genes (annotation_df already loaded above if provided)
    if args.annotation and annotation_df is not None:
        print(f"\n  Annotating clusters with genes...")
        clusters_df = annotate_clusters_with_genes(clusters_df, annotation_df)
        n_annotated = clusters_df['gene_id'].notna().sum()
        print(f"  Annotated {n_annotated:,} clusters ({100*n_annotated/len(clusters_df):.1f}%)")

    # Save clusters
    clusters_path = output_dir / 'cpa_clusters.tsv'
    clusters_df.to_csv(clusters_path, sep='\t', index=False)
    print(f"  Saved clusters to {clusters_path}")

    # Build count matrix
    print(f"\n[3/9] Building cluster count matrix...")
    count_matrix = build_cluster_count_matrix(
        positions_df,
        clusters_df,
        sample_col=args.sample_column,
        count_col=count_col,
        fraction_col=fraction_col,
    )
    print(f"  Matrix shape: {count_matrix.shape[0]:,} clusters × {count_matrix.shape[1]} samples")

    # Save count matrix
    counts_path = output_dir / 'cluster_counts.tsv'
    count_matrix.to_csv(counts_path, sep='\t')

    # Create sample metadata
    sample_names = count_matrix.columns.tolist()

    if args.reference:
        reference_condition = args.reference
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

        # Gene-level
        print("  Gene-level analysis...")
        deseq2_gene_results = run_deseq2_gene_level(
            count_matrix,
            clusters_df,
            sample_metadata,
            reference_condition,
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
        from .analyze.go_enrichment import load_go_annotations

        go_annotations = load_go_annotations(args.go_annotations)

        for condition, result_df in deseq2_gene_results.items():
            # Upregulated genes
            up_genes = result_df[
                (result_df['padj'] < 0.05) & (result_df['log2FoldChange'] > 1)
            ].index.tolist()

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
            down_genes = result_df[
                (result_df['padj'] < 0.05) & (result_df['log2FoldChange'] < -1)
            ].index.tolist()

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

    return 0


def load_corrected_positions(
    filepath: str,
    sample_column: str,
    normalize_chroms: bool = True,
    chrom_format: str = 'ncbi',
    chunk_size: int = 1_000_000,
    max_rows: int = None,
) -> pd.DataFrame:
    """Load corrected positions from TSV file with optional chunked loading.

    Args:
        filepath: Path to TSV file with corrected positions
        sample_column: Column name for sample identifier
        normalize_chroms: Whether to normalize chromosome names
        chrom_format: Target chromosome format ('ncbi', 'ucsc')
        chunk_size: Rows per chunk for large files (default: 1M)
        max_rows: Maximum rows to load (for testing)

    Returns:
        DataFrame with corrected positions
    """
    from ..utils.chromosome import normalize_dataframe_chromosomes
    import os

    # Check file size to determine loading strategy
    file_size = os.path.getsize(filepath)
    file_size_gb = file_size / (1024**3)

    if file_size_gb > 5:
        # Large file - use chunked loading with aggregation
        print(f"  Large file ({file_size_gb:.1f} GB) - using chunked loading...")
        return _load_large_file_chunked(
            filepath, sample_column, normalize_chroms, chrom_format, chunk_size, max_rows
        )

    # Standard loading for smaller files
    df = pd.read_csv(filepath, sep='\t', nrows=max_rows)

    # Handle different position column names
    position_col = None
    for col in ['corrected_position', 'corrected_3prime', 'position']:
        if col in df.columns:
            position_col = col
            break

    if position_col is None:
        raise ValueError("No position column found (tried: corrected_position, corrected_3prime, position)")

    # Standardize to 'corrected_position' if needed
    if position_col != 'corrected_position':
        df['corrected_position'] = df[position_col]
        print(f"  Using '{position_col}' as position column")

    required_cols = ['chrom', 'strand', 'corrected_position', sample_column]
    missing = [c for c in required_cols if c not in df.columns]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Normalize chromosome names for consistent analysis
    if normalize_chroms:
        df = normalize_dataframe_chromosomes(df, 'chrom', chrom_format)
        print(f"  Normalized chromosome names to {chrom_format} format")

    return df


def _load_large_file_chunked(
    filepath: str,
    sample_column: str,
    normalize_chroms: bool,
    chrom_format: str,
    chunk_size: int,
    max_rows: int,
) -> pd.DataFrame:
    """Load large file in chunks, aggregating counts by position.

    Instead of loading all 200M+ rows, aggregate counts by:
    (chrom, strand, position, sample) -> count

    This reduces memory by ~100x for most datasets.
    Uses vectorized pandas operations instead of iterrows for speed.
    """
    from ..utils.chromosome import normalize_chromosome

    print(f"  Aggregating counts by position (chunk_size={chunk_size:,})...")

    # First pass: determine position column name
    header = pd.read_csv(filepath, sep='\t', nrows=0)
    position_col = None
    for col in ['corrected_position', 'corrected_3prime', 'position']:
        if col in header.columns:
            position_col = col
            break

    if position_col is None:
        raise ValueError("No position column found")

    # Check for fraction column (proportional assignment)
    has_fraction = 'fraction' in header.columns

    # Aggregated DataFrames from each chunk
    aggregated_chunks = []

    total_rows = 0
    for chunk_num, chunk in enumerate(pd.read_csv(filepath, sep='\t', chunksize=chunk_size)):
        if max_rows and total_rows >= max_rows:
            break

        # Normalize chromosomes in chunk (vectorized)
        if normalize_chroms:
            chunk['chrom'] = chunk['chrom'].map(
                lambda x: normalize_chromosome(x, chrom_format)
            )

        # VECTORIZED AGGREGATION: Use groupby instead of iterrows
        group_cols = ['chrom', 'strand', position_col, sample_column]

        if has_fraction:
            # Sum fractions for each position/sample combination
            chunk_agg = chunk.groupby(group_cols)['fraction'].sum().reset_index()
            chunk_agg.columns = ['chrom', 'strand', 'corrected_position', sample_column, 'count']
        else:
            # Count occurrences
            chunk_agg = chunk.groupby(group_cols).size().reset_index(name='count')
            chunk_agg.columns = ['chrom', 'strand', 'corrected_position', sample_column, 'count']

        aggregated_chunks.append(chunk_agg)

        total_rows += len(chunk)
        if chunk_num % 10 == 0:
            print(f"    Processed {total_rows:,} rows...", flush=True)

    # Combine all chunk aggregations and re-aggregate
    print(f"  Combining {len(aggregated_chunks)} chunks...")
    combined = pd.concat(aggregated_chunks, ignore_index=True)

    # Final aggregation to merge duplicate keys across chunks
    group_cols = ['chrom', 'strand', 'corrected_position', sample_column]
    final = combined.groupby(group_cols)['count'].sum().reset_index()

    print(f"  Aggregated {total_rows:,} rows into {len(final):,} position/sample combinations")

    return final


def load_annotation(filepath: str, normalize_chroms: bool = True, chrom_format: str = 'ncbi') -> pd.DataFrame:
    """Load gene annotation file (GTF or TSV).

    Args:
        filepath: Path to annotation file
        normalize_chroms: Whether to normalize chromosome names
        chrom_format: Target chromosome format ('ncbi', 'ucsc')

    Returns:
        DataFrame with gene annotations
    """
    from ..utils.chromosome import normalize_dataframe_chromosomes

    if filepath.endswith('.gtf') or filepath.endswith('.gff'):
        # Parse GTF/GFF
        df = _parse_gtf(filepath)
    else:
        # Assume TSV
        df = pd.read_csv(filepath, sep='\t')

    # Normalize chromosome names to match position data
    if normalize_chroms and 'chrom' in df.columns:
        df = normalize_dataframe_chromosomes(df, 'chrom', chrom_format)

    return df


def _detect_exclusion_regions(
    annotation_df: pd.DataFrame,
    data_chroms: list,
) -> tuple:
    """
    Auto-detect mitochondrial chromosomes and rDNA regions from annotation.

    Uses annotation to identify:
    1. Mitochondrial chromosomes (by name or gene content)
    2. rDNA loci (by gene names like RDN*, ribosomal RNA genes)

    Args:
        annotation_df: Gene annotation DataFrame (may be None)
        data_chroms: List of chromosome names in the data

    Returns:
        Tuple of (mito_chroms: set, rdna_regions: list of (chrom, start, end))
    """
    # Default mitochondrial chromosome patterns
    mito_patterns = {'chrM', 'chrMT', 'chrmt', 'MT', 'Mt', 'Mito', 'mitochondrion'}
    mito_ncbi = {'ref|NC_001224|'}  # Yeast mito

    # Find mitochondrial chromosomes in data
    mito_chroms = set()
    for chrom in data_chroms:
        chrom_lower = str(chrom).lower()
        if any(pat.lower() in chrom_lower for pat in mito_patterns):
            mito_chroms.add(chrom)
        if chrom in mito_ncbi:
            mito_chroms.add(chrom)

    # rDNA regions - detect from annotation if available
    rdna_regions = []

    if annotation_df is not None and not annotation_df.empty:
        # Look for mitochondrial genes in annotation
        if 'chrom' in annotation_df.columns:
            for chrom in annotation_df['chrom'].unique():
                chrom_lower = str(chrom).lower()
                if any(pat.lower() in chrom_lower for pat in mito_patterns):
                    mito_chroms.add(chrom)

        # Look for rDNA genes (yeast: RDN5, RDN18, RDN25, RDN37, RDN58, ETS, ITS, NTS)
        # These are the actual ribosomal RNA genes in the rDNA repeat unit
        # Pattern matches gene names starting with RDN, ETS, ITS, NTS (rDNA-specific)
        rdna_gene_patterns = ['^RDN', '^ETS', '^ITS', '^NTS']

        gene_col = 'gene_name' if 'gene_name' in annotation_df.columns else 'gene_id'
        if gene_col in annotation_df.columns:
            import re
            combined_pattern = '|'.join(rdna_gene_patterns)
            rdna_genes = annotation_df[
                annotation_df[gene_col].str.match(combined_pattern, case=False, na=False)
            ]

            if not rdna_genes.empty and 'chrom' in rdna_genes.columns:
                # Group by chromosome and find bounds
                for chrom in rdna_genes['chrom'].unique():
                    chrom_genes = rdna_genes[rdna_genes['chrom'] == chrom]
                    if 'start' in chrom_genes.columns and 'end' in chrom_genes.columns:
                        # Get region bounds with padding
                        start = int(chrom_genes['start'].min()) - 1000
                        end = int(chrom_genes['end'].max()) + 1000
                        start = max(0, start)

                        # Add region if substantial (>5kb suggests rDNA locus)
                        if end - start > 5000:
                            rdna_regions.append((chrom, start, end))

    # Deduplicate and sort rDNA regions
    rdna_regions = list(set(rdna_regions))
    rdna_regions.sort()

    # If no rDNA detected from annotation, use default yeast coordinates
    if not rdna_regions:
        # Yeast rDNA locus on chrXII: ~451,000-469,000 bp
        yeast_rdna_chroms = {'chrXII', 'chr12', 'ref|NC_001144|'}
        for chrom in data_chroms:
            if chrom in yeast_rdna_chroms:
                rdna_regions.append((chrom, 450_000, 490_000))
                break

    return mito_chroms, rdna_regions


def _parse_gtf(filepath: str) -> pd.DataFrame:
    """Parse GTF/GFF3 file to extract gene coordinates.

    Supports both GTF format (key "value") and GFF3 format (key=value).
    For GFF3, extracts:
        - ID -> gene_id
        - Name -> systematic gene name
        - gene -> common gene name (preferred for display)
    """
    genes = []

    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            if fields[2] != 'gene':
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            # Parse attributes - handle both GTF and GFF3 formats
            attrs = {}
            attr_string = fields[8]

            # Detect format: GFF3 uses key=value, GTF uses key "value"
            if '=' in attr_string:
                # GFF3 format: ID=YAL069W;Name=YAL069W;gene=PAU8
                for attr in attr_string.split(';'):
                    attr = attr.strip()
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        # URL decode common escapes
                        value = value.replace('%20', ' ').replace('%3B', ';').replace('%2C', ',')
                        attrs[key] = value
            else:
                # GTF format: gene_id "YAL069W"; gene_name "PAU8"
                for attr in attr_string.split(';'):
                    attr = attr.strip()
                    if ' ' in attr:
                        key, value = attr.split(' ', 1)
                        attrs[key] = value.strip('"')

            # Extract gene_id - try GFF3 keys first, then GTF
            gene_id = attrs.get('ID') or attrs.get('gene_id') or f'{chrom}_{start}'

            # Extract gene_name - prefer common name, fall back to systematic name
            # GFF3: 'gene' has common name, 'Name' has systematic name
            # GTF: 'gene_name' has the name
            gene_name = attrs.get('gene') or attrs.get('gene_name') or attrs.get('Name') or gene_id

            genes.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_id': gene_id,
                'gene_name': gene_name,
            })

    return pd.DataFrame(genes)


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
        help='Input TSV with corrected positions',
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

    parser.set_defaults(func=run_analyze)

    return parser


def generate_bedgraphs(
    positions_df: pd.DataFrame,
    output_dir: Path,
    sample_column: str = 'sample',
    position_column: str = 'corrected_3prime',
    normalize_rpm: bool = True,
) -> None:
    """
    Generate strand-specific bedgraph files per condition from corrected positions.

    Args:
        positions_df: DataFrame with corrected positions
        output_dir: Output directory for bedgraph files
        sample_column: Column name for sample/replicate identifier
        position_column: Column name for corrected position
        normalize_rpm: If True, output RPM-normalized values
    """
    from collections import defaultdict

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Chromosome order for sorted output
    CHROM_ORDER = [
        'chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII',
        'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrmt'
    ]

    # Extract condition from sample name
    conditions = positions_df[sample_column].apply(
        lambda x: extract_condition_from_sample(x)
    ).unique()

    print(f"  Generating bedgraphs for {len(conditions)} conditions...")

    for condition in conditions:
        # Filter to this condition
        condition_mask = positions_df[sample_column].apply(
            lambda x: extract_condition_from_sample(x) == condition
        )
        cond_df = positions_df[condition_mask]

        if len(cond_df) == 0:
            continue

        # Count per chrom/strand/position
        plus_counts = defaultdict(lambda: defaultdict(int))
        minus_counts = defaultdict(lambda: defaultdict(int))

        plus_df = cond_df[cond_df['strand'] == '+']
        minus_df = cond_df[cond_df['strand'] == '-']

        for (chrom, pos), count in plus_df.groupby(['chrom', position_column]).size().items():
            plus_counts[chrom][pos] += count

        for (chrom, pos), count in minus_df.groupby(['chrom', position_column]).size().items():
            minus_counts[chrom][pos] += count

        total_reads = len(cond_df)
        rpm_factor = 1e6 / total_reads if normalize_rpm and total_reads > 0 else 1.0

        # Write bedgraph files
        for strand_name, counts_dict in [('plus', plus_counts), ('minus', minus_counts)]:
            output_path = output_dir / f"{condition}_{strand_name}.bedgraph"

            with open(output_path, 'w') as f:
                f.write(f'track type=bedGraph name="{condition}_{strand_name}" '
                        f'description="RECTIFY 3\' ends ({strand_name} strand)"\n')

                for chrom in CHROM_ORDER:
                    if chrom not in counts_dict:
                        continue

                    for pos in sorted(counts_dict[chrom].keys()):
                        count = counts_dict[chrom][pos]
                        value = count * rpm_factor if normalize_rpm else count
                        start = int(pos) - 1
                        end = int(pos)
                        f.write(f"{chrom}\t{start}\t{end}\t{value:.4f}\n")

        print(f"    {condition}: {total_reads:,} reads")
