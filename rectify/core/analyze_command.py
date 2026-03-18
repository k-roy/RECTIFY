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
)
from .analyze.clustering import annotate_clusters_with_genes
from .analyze.deseq2 import extract_condition_from_sample


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

    # Load corrected positions
    print(f"\n[1/9] Loading corrected positions from {args.input}...")
    positions_df = load_corrected_positions(args.input, args.sample_column)
    print(f"  Loaded {len(positions_df):,} positions from {positions_df[args.sample_column].nunique()} samples")

    # Form clusters
    print(f"\n[2/9] Forming CPA clusters (distance={args.cluster_distance}bp)...")
    clusters_df = cluster_cpa_sites(
        positions_df,
        cluster_distance=args.cluster_distance,
        min_reads=args.min_reads,
        count_col=args.count_column if args.count_column else None,
    )
    print(f"  Formed {len(clusters_df):,} clusters")

    # Annotate with genes
    if args.annotation:
        print(f"\n  Annotating clusters with genes...")
        annotation_df = load_annotation(args.annotation)
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
        count_col=args.count_column if args.count_column else None,
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

    print(f"\n" + "=" * 70)
    print(f"Analysis complete!")
    print(f"  Output directory: {output_dir}")
    print(f"  HTML report: {html_path}")
    print("=" * 70)

    return 0


def load_corrected_positions(
    filepath: str,
    sample_column: str,
    normalize_chroms: bool = True,
    chrom_format: str = 'ncbi',
) -> pd.DataFrame:
    """Load corrected positions from TSV file.

    Args:
        filepath: Path to TSV file with corrected positions
        sample_column: Column name for sample identifier
        normalize_chroms: Whether to normalize chromosome names
        chrom_format: Target chromosome format ('ncbi', 'ucsc')

    Returns:
        DataFrame with corrected positions
    """
    from ..utils.chromosome import normalize_dataframe_chromosomes

    df = pd.read_csv(filepath, sep='\t')

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


def load_annotation(filepath: str) -> pd.DataFrame:
    """Load gene annotation file (GTF or TSV)."""
    if filepath.endswith('.gtf') or filepath.endswith('.gff'):
        # Parse GTF/GFF
        return _parse_gtf(filepath)
    else:
        # Assume TSV
        return pd.read_csv(filepath, sep='\t')


def _parse_gtf(filepath: str) -> pd.DataFrame:
    """Parse GTF file to extract gene coordinates."""
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

            # Parse attributes
            attrs = {}
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if ' ' in attr:
                    key, value = attr.split(' ', 1)
                    attrs[key] = value.strip('"')

            gene_id = attrs.get('gene_id', f'{chrom}_{start}')
            gene_name = attrs.get('gene_name', gene_id)

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

    parser.set_defaults(func=run_analyze)

    return parser
