#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genomic Distribution Analysis Module

Classifies 3' end positions by genomic region (UTR3, CDS, UTR5, intergenic, etc.)
and generates pie charts showing the distribution.

Author: Kevin R. Roy
Date: 2026-03-24
"""

from typing import Dict, List, Optional, Tuple
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# Colors for genomic categories (Wong palette + extensions)
CATEGORY_COLORS = {
    'UTR3': '#009E73',       # Green - expected termination site
    'CDS': '#0072B2',        # Blue
    'UTR5': '#56B4E9',       # Light blue
    'Intergenic': '#E69F00', # Orange
    'ncRNA': '#999999',      # Gray
    'rDNA': '#666666',       # Dark gray
    'Other': '#DDDDDD',      # Light gray
}


def classify_positions_by_region(
    positions_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    position_column: str = 'corrected_position',
) -> pd.DataFrame:
    """
    Classify each 3' end position by genomic region based on gene annotation.

    Categories:
    - UTR3: 3' UTR (downstream of CDS stop, within gene bounds)
    - CDS: Coding sequence
    - UTR5: 5' UTR (upstream of CDS start, within gene bounds)
    - Intergenic: Between genes
    - ncRNA: Non-coding RNA genes (if annotated)

    Args:
        positions_df: DataFrame with 3' end positions (chrom, strand, position)
        annotation_df: Gene annotation DataFrame (chrom, start, end, strand, gene_id)
        position_column: Column name for position

    Returns:
        DataFrame with added 'genomic_region' column
    """
    # Build interval trees for fast lookup
    from intervaltree import IntervalTree

    # Build trees per chromosome/strand
    gene_trees = defaultdict(IntervalTree)

    for _, gene in annotation_df.iterrows():
        chrom = gene['chrom']
        start = gene['start']
        end = gene['end']
        strand = gene.get('strand', '+')

        # Store gene info in interval
        gene_info = {
            'gene_id': gene.get('gene_id', ''),
            'gene_name': gene.get('gene_name', ''),
            'strand': strand,
            'start': start,
            'end': end,
        }

        # Use strand-specific tree
        key = (chrom, strand)
        gene_trees[key][start:end+1] = gene_info

    # Classify each position
    regions = []

    for _, row in positions_df.iterrows():
        chrom = row['chrom']
        pos = row[position_column]
        strand = row.get('strand', '+')

        # Look up in strand-specific tree
        key = (chrom, strand)
        overlaps = gene_trees[key][pos]

        if not overlaps:
            # Check antisense
            anti_key = (chrom, '-' if strand == '+' else '+')
            anti_overlaps = gene_trees[anti_key][pos]

            if anti_overlaps:
                regions.append('Antisense')
            else:
                regions.append('Intergenic')
        else:
            # Inside a gene - determine UTR vs CDS
            # For now, classify as gene body (would need CDS annotations for precise UTR)
            gene_info = list(overlaps)[0].data
            gene_start = gene_info['start']
            gene_end = gene_info['end']
            gene_strand = gene_info['strand']

            # Simple heuristic: first/last 10% is UTR, middle is CDS
            gene_len = gene_end - gene_start
            utr_size = max(50, int(gene_len * 0.1))

            if gene_strand == '+':
                if pos <= gene_start + utr_size:
                    regions.append('UTR5')
                elif pos >= gene_end - utr_size:
                    regions.append('UTR3')
                else:
                    regions.append('CDS')
            else:
                # Minus strand - UTRs are reversed
                if pos >= gene_end - utr_size:
                    regions.append('UTR5')
                elif pos <= gene_start + utr_size:
                    regions.append('UTR3')
                else:
                    regions.append('CDS')

    positions_df = positions_df.copy()
    positions_df['genomic_region'] = regions

    return positions_df


def calculate_genomic_distribution(
    positions_df: pd.DataFrame,
    region_column: str = 'genomic_region',
    count_column: Optional[str] = None,
) -> Dict[str, int]:
    """
    Calculate distribution of 3' ends across genomic regions.

    Args:
        positions_df: DataFrame with genomic_region column
        region_column: Column name for genomic region
        count_column: Optional column for pre-aggregated counts

    Returns:
        Dict mapping region -> count
    """
    if count_column and count_column in positions_df.columns:
        return positions_df.groupby(region_column)[count_column].sum().to_dict()
    else:
        return positions_df[region_column].value_counts().to_dict()


def plot_genomic_distribution_pie(
    distribution: Dict[str, int],
    output_path: str,
    title: str = "Distribution of 3' Ends by Genomic Region",
    min_pct: float = 1.0,
) -> str:
    """
    Create pie chart showing distribution of 3' ends across genomic regions.

    Args:
        distribution: Dict mapping region -> count
        output_path: Output file path
        title: Plot title
        min_pct: Minimum percentage to show as separate slice

    Returns:
        Path to saved plot
    """
    total = sum(distribution.values())

    # Sort by expected biological order
    preferred_order = ['UTR3', 'CDS', 'UTR5', 'Intergenic', 'Antisense', 'ncRNA', 'rDNA', 'Other']

    labels = []
    sizes = []
    colors = []
    other_total = 0

    for cat in preferred_order:
        if cat in distribution:
            count = distribution[cat]
            pct = 100 * count / total

            if pct >= min_pct:
                labels.append(cat)
                sizes.append(count)
                colors.append(CATEGORY_COLORS.get(cat, '#888888'))
            else:
                other_total += count

    # Add any categories not in preferred order
    for cat, count in distribution.items():
        if cat not in preferred_order:
            pct = 100 * count / total
            if pct >= min_pct:
                labels.append(cat)
                sizes.append(count)
                colors.append(CATEGORY_COLORS.get(cat, '#888888'))
            else:
                other_total += count

    # Add "Other" if significant
    if other_total > 0 and (100 * other_total / total) >= 0.5:
        labels.append('Other')
        sizes.append(other_total)
        colors.append('#DDDDDD')

    # Create pie chart
    fig, ax = plt.subplots(figsize=(10, 8))

    def make_autopct(values):
        def my_autopct(pct):
            if pct >= 2:
                return f'{pct:.1f}%'
            elif pct >= 1:
                return f'{pct:.1f}%'
            else:
                return ''
        return my_autopct

    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=labels,
        colors=colors,
        autopct=make_autopct(sizes),
        startangle=90,
        pctdistance=0.75,
        labeldistance=1.1,
        textprops={'fontsize': 11}
    )

    # Make percentage text bold
    for autotext in autotexts:
        autotext.set_fontweight('bold')

    ax.set_title(f"{title}\n(n={total:,} reads)", fontsize=14, fontweight='bold')

    # Add legend with counts
    legend_labels = [f'{lab}: {cnt:,} ({100*cnt/total:.1f}%)'
                     for lab, cnt in zip(labels, sizes)]
    ax.legend(wedges, legend_labels, title="Region", loc="center left",
             bbox_to_anchor=(1, 0, 0.5, 1), fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    return output_path


def plot_genomic_distribution_comparison(
    distributions: Dict[str, Dict[str, int]],
    output_path: str,
    title: str = "Genomic Distribution Comparison",
    condition_labels: Optional[Dict[str, str]] = None,
) -> str:
    """
    Create grouped bar chart comparing genomic distribution across conditions.

    Args:
        distributions: Dict mapping condition -> {region -> count}
        output_path: Output file path
        title: Plot title
        condition_labels: Optional mapping of condition names to display labels

    Returns:
        Path to saved plot
    """
    # Get all regions
    all_regions = set()
    for dist in distributions.values():
        all_regions.update(dist.keys())

    # Define order
    preferred_order = ['UTR3', 'CDS', 'UTR5', 'Intergenic', 'Antisense', 'ncRNA', 'rDNA']
    regions = [r for r in preferred_order if r in all_regions]
    regions.extend([r for r in all_regions if r not in preferred_order])

    # Calculate percentages
    conditions = list(distributions.keys())
    data = np.zeros((len(conditions), len(regions)))

    for i, cond in enumerate(conditions):
        total = sum(distributions[cond].values())
        for j, region in enumerate(regions):
            data[i, j] = 100 * distributions[cond].get(region, 0) / total if total > 0 else 0

    # Create grouped bar chart
    fig, ax = plt.subplots(figsize=(12, 6))

    x = np.arange(len(regions))
    width = 0.8 / len(conditions)

    for i, cond in enumerate(conditions):
        offset = (i - len(conditions)/2 + 0.5) * width
        label = condition_labels.get(cond, cond) if condition_labels else cond
        ax.bar(x + offset, data[i], width, label=label, alpha=0.8)

    ax.set_xlabel('Genomic Region', fontsize=12)
    ax.set_ylabel("% of 3' Ends", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(regions, rotation=45, ha='right')
    ax.legend(title='Condition')
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    return output_path


def run_genomic_distribution_analysis(
    positions_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    output_dir: str,
    sample_column: str = 'sample',
    position_column: str = 'corrected_position',
    count_column: Optional[str] = None,
    condition_labels: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """
    Run full genomic distribution analysis and generate plots.

    Args:
        positions_df: DataFrame with 3' end positions
        annotation_df: Gene annotation DataFrame
        output_dir: Output directory for plots
        sample_column: Column for sample/condition identifier
        position_column: Column for position
        count_column: Optional column for pre-aggregated counts
        condition_labels: Optional display labels for conditions

    Returns:
        Dict mapping plot_name -> output_path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_files = {}

    # Extract condition from sample name
    from ..analyze.deseq2 import extract_condition_from_sample

    positions_df = positions_df.copy()
    positions_df['condition'] = positions_df[sample_column].apply(extract_condition_from_sample)

    conditions = positions_df['condition'].unique()

    print(f"  Running genomic distribution analysis for {len(conditions)} conditions...")

    # Classify positions (this can be slow for large datasets - sample if needed)
    n_positions = len(positions_df)
    if n_positions > 1_000_000:
        print(f"    Sampling {1_000_000:,} positions for classification...")
        sample_df = positions_df.sample(n=1_000_000, random_state=42)
    else:
        sample_df = positions_df

    # Classify by genomic region
    print(f"    Classifying positions by genomic region...")
    classified_df = classify_positions_by_region(
        sample_df, annotation_df, position_column
    )

    # Calculate distribution per condition
    distributions = {}

    for condition in conditions:
        cond_df = classified_df[classified_df['condition'] == condition]
        distributions[condition] = calculate_genomic_distribution(
            cond_df, 'genomic_region', count_column
        )

        # Per-condition pie chart
        total = sum(distributions[condition].values())
        pie_path = output_dir / f'genomic_distribution_pie_{condition}.png'
        label = condition_labels.get(condition, condition) if condition_labels else condition
        plot_genomic_distribution_pie(
            distributions[condition],
            str(pie_path),
            title=f"3' End Distribution: {label}",
        )
        output_files[f'pie_{condition}'] = str(pie_path)
        print(f"    {label}: {total:,} reads classified")

    # Comparison plot if multiple conditions
    if len(conditions) > 1:
        comparison_path = output_dir / 'genomic_distribution_comparison.png'
        plot_genomic_distribution_comparison(
            distributions,
            str(comparison_path),
            title="3' End Genomic Distribution Comparison",
            condition_labels=condition_labels,
        )
        output_files['comparison'] = str(comparison_path)

    # Summary table
    summary_rows = []
    for condition in conditions:
        total = sum(distributions[condition].values())
        for region, count in distributions[condition].items():
            summary_rows.append({
                'condition': condition,
                'region': region,
                'count': count,
                'percent': 100 * count / total if total > 0 else 0,
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_dir / 'genomic_distribution_summary.tsv'
    summary_df.to_csv(summary_path, sep='\t', index=False)
    output_files['summary'] = str(summary_path)

    return output_files
