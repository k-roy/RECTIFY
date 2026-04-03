#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genomic Distribution Analysis Module

Two complementary distribution analyses, each producing a single figure per
run (horizontal bar chart across all conditions + one pie per condition):

─── 3' END DISTRIBUTION ─────────────────────────────────────────────────────
Classifies individual corrected 3' end *positions* by the genomic feature they
fall in.  Mirrors Xu et al. 2009 (Nature 457:1033) panels B/C.

  Priority: UTR3 > snoRNA±300bp > CUT > SUT_XUT > UTR5_CDS > Antisense > Intergenic

  X-axis (barh): % of total reads   Right panel: unique cluster count

─── TRANSCRIPT BODY DISTRIBUTION ────────────────────────────────────────────
Classifies each read by the RNA biotype of the feature it overlaps most (by bp).
Assignment uses the read's full alignment span (alignment_start → alignment_end).

  Categories: protein_coding, CUT, SUT_XUT, snoRNA, tRNA, rRNA, LTR,
              pseudogene, Antisense, Intergenic

  X-axis (barh): % of total reads   Right panel: raw read count per category

─── FIGURE LAYOUT ───────────────────────────────────────────────────────────
Both analyses produce one PNG with:
  Top   : horizontal bar chart (all conditions overlaid, broken x-axis)
  Bottom: pie chart grid (one pie per condition, replicates merged)

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
import matplotlib.ticker as mticker


# Internal category names → display labels
CATEGORY_LABELS = {
    'UTR3':     "3' UTR",
    'snoRNA':   'snoRNA +/- 300 bp',
    'Intergenic': 'intergenic / intronic',
    'CUT':      'CUTs',
    'SUT_XUT':  'SUTs / XUTs',
    'UTR5_CDS': "5' UTR / CDS",
    'Antisense': 'antisense CDS',
}

# Display order (top → bottom in horizontal bar chart, matches figure panels B/C)
CATEGORY_ORDER = ['UTR3', 'snoRNA', 'Intergenic', 'CUT', 'SUT_XUT', 'UTR5_CDS', 'Antisense']

# Colors for genomic categories (Wong palette + extensions)
CATEGORY_COLORS = {
    'UTR3':      '#009E73',   # Green
    'snoRNA':    '#CC79A7',   # Pink/purple
    'Intergenic': '#E69F00',  # Orange
    'CUT':       '#D55E00',   # Vermilion
    'SUT_XUT':   '#F0E442',   # Yellow
    'UTR5_CDS':  '#0072B2',   # Blue
    'Antisense': '#56B4E9',   # Sky blue
    # Legacy names kept for pie chart backward compatibility
    'CDS':       '#0072B2',
    'UTR5':      '#56B4E9',
    'ncRNA':     '#999999',
    'rDNA':      '#666666',
    'Other':     '#DDDDDD',
}

# Flanking window around snoRNA genes for classification (bp)
SNORNA_FLANK_BP = 300

# ── Transcript body category system ──────────────────────────────────────────

BODY_CATEGORY_LABELS = {
    'protein_coding': 'protein-coding',
    'CUT':            'CUTs',
    'SUT_XUT':        'SUTs / XUTs',
    'snoRNA':         'snoRNA',
    'tRNA':           'tRNA',
    'rRNA':           'rRNA',
    'LTR':            'LTR / retrotransposon',
    'pseudogene':     'pseudogene',
    'Antisense':      'antisense',
    'Intergenic':     'intergenic',
}

# Display order (most abundant first for a typical poly(A)+ DRS library)
BODY_CATEGORY_ORDER = [
    'protein_coding', 'CUT', 'SUT_XUT', 'snoRNA',
    'tRNA', 'rRNA', 'LTR', 'pseudogene', 'Antisense', 'Intergenic',
]

BODY_CATEGORY_COLORS = {
    'protein_coding': '#009E73',   # Green
    'CUT':            '#D55E00',   # Vermilion
    'SUT_XUT':        '#F0E442',   # Yellow
    'snoRNA':         '#CC79A7',   # Pink/purple
    'tRNA':           '#56B4E9',   # Sky blue
    'rRNA':           '#E69F00',   # Orange
    'LTR':            '#999999',   # Gray
    'pseudogene':     '#BBBBBB',   # Light gray
    'Antisense':      '#0072B2',   # Blue
    'Intergenic':     '#DDDDDD',   # Very light gray
}

# GFF feature type → body category (lowercase matching)
_BODY_FEATURE_MAP = {
    'mrna': 'protein_coding', 'cds': 'protein_coding',
    'protein_coding': 'protein_coding',
    'cut': 'CUT',
    'sut': 'SUT_XUT', 'xut': 'SUT_XUT',
    'snorna_gene': 'snoRNA', 'snorna': 'snoRNA',
    'trna_gene': 'tRNA', 'trna': 'tRNA',
    'rrna_gene': 'rRNA', 'rrna': 'rRNA', 'rdna_locus': 'rRNA', 'rdna': 'rRNA',
    'ltr_retrotransposon': 'LTR', 'transposable_element_gene': 'LTR',
    'ty_element': 'LTR', 'retrotransposon': 'LTR',
    'pseudogene': 'pseudogene',
    'ncrna_gene': 'Intergenic', 'ncrna': 'Intergenic',  # catch-all ncRNA
}

# Priority order for body category ties (lower index = higher priority)
_BODY_PRIORITY = {cat: i for i, cat in enumerate(BODY_CATEGORY_ORDER)}


def _build_feature_trees(annotation_df: pd.DataFrame) -> dict:
    """
    Build per-chromosome/strand interval trees keyed by feature type.

    Returns a dict: feature_type -> { (chrom, strand) -> IntervalTree }
    Feature types handled: 'gene' (protein-coding), 'mRNA', 'UTR3', 'UTR5',
    'CDS', 'CUT', 'SUT', 'XUT', 'snoRNA' (flanked by SNORNA_FLANK_BP).
    """
    from intervaltree import IntervalTree

    trees: Dict[str, dict] = defaultdict(lambda: defaultdict(IntervalTree))

    ft_col = 'feature_type' if 'feature_type' in annotation_df.columns else 'feature'

    for _, row in annotation_df.iterrows():
        chrom = row['chrom']
        start = int(row['start'])
        end   = int(row['end'])
        strand = row.get('strand', '.')
        ft = str(row.get(ft_col, 'gene')).lower()

        if start >= end:
            continue

        data = {'strand': strand, 'start': start, 'end': end}

        if ft in ('cut',):
            trees['CUT'][(chrom, strand)][start:end] = data
        elif ft in ('sut',):
            trees['SUT'][(chrom, strand)][start:end] = data
        elif ft in ('xut',):
            trees['XUT'][(chrom, strand)][start:end] = data
        elif ft in ('snorna_gene', 'snorna'):
            # Expand by flank on both strands (snoRNAs are processed on both)
            fs = max(0, start - SNORNA_FLANK_BP)
            fe = end + SNORNA_FLANK_BP
            for s in ('+', '-'):
                trees['snoRNA'][(chrom, s)][fs:fe] = data
        elif ft in ('utr3', "3'utr", "three_prime_utr"):
            trees['UTR3'][(chrom, strand)][start:end] = data
        elif ft in ('utr5', "5'utr", "five_prime_utr"):
            trees['UTR5'][(chrom, strand)][start:end] = data
        elif ft in ('cds',):
            trees['CDS'][(chrom, strand)][start:end] = data
            trees['CDS_anti'][(chrom, '+' if strand == '-' else '-')][start:end] = data
        elif ft in ('gene', 'mrna', 'protein_coding', 'transcript'):
            trees['gene'][(chrom, strand)][start:end] = data

    return trees


def classify_positions_by_region(
    positions_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    position_column: str = 'corrected_position',
) -> pd.DataFrame:
    """
    Classify each 3' end position by genomic region.

    Priority (highest → lowest):
      UTR3 → snoRNA → CUT → SUT_XUT → UTR5_CDS → Antisense → Intergenic

    Args:
        positions_df: DataFrame with columns chrom, strand, and position_column
        annotation_df: Gene annotation DataFrame; must have a feature_type column
                       (values: gene, mRNA, CDS, UTR3, UTR5, CUT, SUT, XUT,
                       snoRNA_gene, snoRNA, etc.)
        position_column: Column name for the 3' end position (0-based)

    Returns:
        DataFrame with added 'genomic_region' column using CATEGORY_ORDER labels.
    """
    trees = _build_feature_trees(annotation_df)

    def _hits(tree_key, chrom, strand, pos):
        return bool(trees.get(tree_key, {}).get((chrom, strand), set()).__class__.__mro__) and \
               bool(trees[tree_key][(chrom, strand)][pos]) if tree_key in trees else False

    # Faster lookup helper
    def hits(tree_key, chrom, strand, pos):
        if tree_key not in trees:
            return False
        t = trees[tree_key].get((chrom, strand))
        return bool(t[pos]) if t else False

    regions = []
    for _, row in positions_df.iterrows():
        chrom = row['chrom']
        pos   = int(row[position_column])
        strand = row.get('strand', '+')
        anti  = '-' if strand == '+' else '+'

        if hits('UTR3', chrom, strand, pos):
            cat = 'UTR3'
        elif hits('snoRNA', chrom, strand, pos):
            cat = 'snoRNA'
        elif hits('CUT', chrom, strand, pos):
            cat = 'CUT'
        elif hits('SUT', chrom, strand, pos) or hits('XUT', chrom, strand, pos):
            cat = 'SUT_XUT'
        elif hits('UTR5', chrom, strand, pos) or hits('CDS', chrom, strand, pos):
            cat = 'UTR5_CDS'
        elif hits('CDS_anti', chrom, anti, pos):
            cat = 'Antisense'
        else:
            # Fall back: is the position inside any gene body on same strand?
            if hits('gene', chrom, strand, pos):
                # Gene body hit without explicit UTR/CDS annotation — use 10% heuristic
                gene_data = list(trees['gene'][(chrom, strand)][pos])[0].data
                gs, ge = gene_data['start'], gene_data['end']
                utr_sz = max(50, int((ge - gs) * 0.1))
                if strand == '+':
                    cat = 'UTR3' if pos >= ge - utr_sz else ('UTR5_CDS' if pos <= gs + utr_sz else 'UTR5_CDS')
                else:
                    cat = 'UTR3' if pos <= gs + utr_sz else 'UTR5_CDS'
            else:
                cat = 'Intergenic'

        regions.append(cat)

    out = positions_df.copy()
    out['genomic_region'] = regions
    return out


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


def plot_genomic_distribution_pie_grid(
    distributions: Dict[str, Dict[str, int]],
    output_path: str,
    title: str = "Genomic Distribution by Condition",
    condition_labels: Optional[Dict[str, str]] = None,
    wt_first: bool = True,
    ncols: int = None,
) -> str:
    """
    Create a grid of pie charts showing all conditions side-by-side for easy comparison.

    Args:
        distributions: Dict mapping condition -> {region -> count}
        output_path: Output file path
        title: Overall figure title
        condition_labels: Optional mapping of condition names to display labels
        wt_first: If True, put WT/control condition first (default: True)
        ncols: Number of columns in grid (default: auto based on count)

    Returns:
        Path to saved plot
    """
    conditions = list(distributions.keys())

    # Sort conditions: WT first if requested
    if wt_first:
        wt_names = ['wt', 'WT', 'wild_type', 'wildtype', 'control', 'by4742', 'BY4742']
        wt_conditions = [c for c in conditions if any(w in c.lower() for w in ['wt', 'wild', 'control', 'by4742'])]
        other_conditions = [c for c in conditions if c not in wt_conditions]
        conditions = sorted(wt_conditions) + sorted(other_conditions)

    n_conditions = len(conditions)

    # Determine grid layout
    if ncols is None:
        if n_conditions <= 2:
            ncols = n_conditions
        elif n_conditions <= 4:
            ncols = 2
        elif n_conditions <= 9:
            ncols = 3
        else:
            ncols = 4

    nrows = (n_conditions + ncols - 1) // ncols

    # Create figure
    fig_width = 5 * ncols
    fig_height = 4.5 * nrows
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))

    # Flatten axes for easy iteration
    if n_conditions == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]

    # Preferred order for consistency across all pies
    preferred_order = ['UTR3', 'CDS', 'UTR5', 'Intergenic', 'Antisense', 'ncRNA', 'rDNA', 'Other']

    for idx, condition in enumerate(conditions):
        ax = axes[idx]
        distribution = distributions[condition]
        total = sum(distribution.values())

        # Build labels/sizes in consistent order
        labels = []
        sizes = []
        colors = []
        other_total = 0

        for cat in preferred_order:
            if cat in distribution:
                count = distribution[cat]
                pct = 100 * count / total if total > 0 else 0
                if pct >= 1.0:
                    labels.append(cat)
                    sizes.append(count)
                    colors.append(CATEGORY_COLORS.get(cat, '#888888'))
                else:
                    other_total += count

        # Add categories not in preferred order
        for cat, count in distribution.items():
            if cat not in preferred_order:
                pct = 100 * count / total if total > 0 else 0
                if pct >= 1.0:
                    labels.append(cat)
                    sizes.append(count)
                    colors.append(CATEGORY_COLORS.get(cat, '#888888'))
                else:
                    other_total += count

        if other_total > 0 and (100 * other_total / total) >= 0.5:
            labels.append('Other')
            sizes.append(other_total)
            colors.append('#DDDDDD')

        def make_autopct(pct):
            return f'{pct:.0f}%' if pct >= 5 else ''

        wedges, texts, autotexts = ax.pie(
            sizes,
            labels=None,  # No labels on pie itself
            colors=colors,
            autopct=make_autopct,
            startangle=90,
            pctdistance=0.7,
            textprops={'fontsize': 9, 'fontweight': 'bold'}
        )

        # Get display label
        label = condition_labels.get(condition, condition) if condition_labels else condition
        ax.set_title(f"{label}\n(n={total:,})", fontsize=11, fontweight='bold')

    # Hide empty subplots
    for idx in range(n_conditions, len(axes)):
        axes[idx].set_visible(False)

    # Add shared legend
    # Create legend handles from first plot's data
    first_dist = distributions[conditions[0]]
    total = sum(first_dist.values())
    legend_labels = []
    legend_colors = []

    for cat in preferred_order:
        if any(cat in d for d in distributions.values()):
            legend_labels.append(cat)
            legend_colors.append(CATEGORY_COLORS.get(cat, '#888888'))

    # Add legend at bottom
    legend_handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for c in legend_colors]
    fig.legend(legend_handles, legend_labels, loc='lower center',
               ncol=len(legend_labels), fontsize=10, frameon=False,
               bbox_to_anchor=(0.5, 0.02))

    # Add overall title
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
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


def count_clusters_by_region(
    positions_df: pd.DataFrame,
    region_column: str = 'genomic_region',
    position_column: str = 'corrected_position',
) -> Dict[str, int]:
    """
    Count unique 3' end clusters (unique chrom/strand/position tuples) per region.

    Args:
        positions_df: DataFrame with genomic_region, chrom, strand, and position columns
        region_column: Column holding genomic category
        position_column: Column holding 3' end position

    Returns:
        Dict mapping region -> number of unique cluster positions
    """
    return (
        positions_df
        .drop_duplicates(subset=['chrom', 'strand', position_column])
        .groupby(region_column)
        .size()
        .to_dict()
    )


def plot_genomic_distribution_barh(
    distributions: Dict[str, Dict[str, int]],
    cluster_counts: Dict[str, Dict[str, int]],
    output_path: str,
    title: str = "3' End Distribution by Genomic Category",
    condition_labels: Optional[Dict[str, str]] = None,
    condition_colors: Optional[Dict[str, str]] = None,
    axis_break: float = 10.0,
) -> str:
    """
    Horizontal bar chart matching Xu et al. 2009 figure panels B/C style.

    Shows % of total reads per genomic category per condition, with an optional
    broken x-axis to accommodate both the dominant category (3' UTR ~90%) and
    the minority categories (0–10%).

    Left panel : % of total reads (0 → axis_break), broken to show high end
    Right annotation: cluster count per category per condition

    Args:
        distributions:    condition -> {category -> read_count}
        cluster_counts:   condition -> {category -> unique_cluster_count}
        output_path:      Output file path
        title:            Figure title
        condition_labels: Optional display names for conditions
        condition_colors: Optional bar colors per condition (default: black/gray cycle)
        axis_break:       x value where left panel clips; right panel shows remainder

    Returns:
        Path to saved plot
    """
    conditions = list(distributions.keys())
    n_conds = len(conditions)

    # Default colors: black for first (WT), grays for rest — matches figure style
    default_colors = ['#000000', '#AAAAAA', '#555555', '#CCCCCC', '#333333']
    if condition_colors is None:
        condition_colors = {c: default_colors[i % len(default_colors)]
                            for i, c in enumerate(conditions)}

    # Resolve display labels
    labels_fn = lambda c: (condition_labels or {}).get(c, c)

    # Calculate percentages in CATEGORY_ORDER
    pct: Dict[str, Dict[str, float]] = {}
    for cond, dist in distributions.items():
        total = sum(dist.values()) or 1
        pct[cond] = {cat: 100.0 * dist.get(cat, 0) / total for cat in CATEGORY_ORDER}

    # ---- Figure layout: left (broken axis) + right (cluster counts) --------
    fig = plt.figure(figsize=(10, 4.5))
    # Two axes share the same y: left = small %, right = overflow + cluster counts
    # Use gridspec: [break_left | break_right | spacer | cluster_panel]
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(1, 3, figure=fig, width_ratios=[axis_break, 100 - axis_break, 60],
                  wspace=0.08)
    ax_left  = fig.add_subplot(gs[0])
    ax_right = fig.add_subplot(gs[1], sharey=ax_left)
    ax_clust = fig.add_subplot(gs[2], sharey=ax_left)

    n_cats = len(CATEGORY_ORDER)
    bar_height = 0.72 / max(n_conds, 1)
    y_positions = np.arange(n_cats)

    for i, cond in enumerate(conditions):
        offset = (i - (n_conds - 1) / 2) * bar_height
        color  = condition_colors[cond]
        label  = labels_fn(cond)
        vals   = [pct[cond][cat] for cat in CATEGORY_ORDER]

        ax_left.barh(y_positions + offset, vals, bar_height * 0.9,
                     color=color, label=label, alpha=0.85)
        ax_right.barh(y_positions + offset, vals, bar_height * 0.9,
                      color=color, alpha=0.85)

        # Cluster counts as text on right panel
        clust = cluster_counts.get(cond, {})
        for j, cat in enumerate(CATEGORY_ORDER):
            n_clust = clust.get(cat, 0)
            if n_clust > 0:
                ax_clust.barh(y_positions[j] + offset, n_clust, bar_height * 0.9,
                              color=color, alpha=0.85)

    # ---- Broken-axis formatting ----
    ax_left.set_xlim(0, axis_break)
    ax_right.set_xlim(axis_break, 100)

    # Hide the spines between the two axes
    ax_left.spines['right'].set_visible(False)
    ax_right.spines['left'].set_visible(False)
    ax_left.yaxis.tick_left()
    ax_right.tick_params(left=False)

    # Diagonal break markers
    d = 0.015
    kwargs = dict(transform=ax_left.transAxes, color='k', clip_on=False, lw=1)
    ax_left.plot((1-d, 1+d), (-d, +d), **kwargs)
    ax_left.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
    kwargs.update(transform=ax_right.transAxes)
    ax_right.plot((-d, +d), (-d, +d), **kwargs)
    ax_right.plot((-d, +d), (1-d, 1+d), **kwargs)

    # Y-axis labels (only on left panel)
    ax_left.set_yticks(y_positions)
    ax_left.set_yticklabels([CATEGORY_LABELS[c] for c in CATEGORY_ORDER], fontsize=10)
    ax_left.invert_yaxis()

    # X-axis labels
    ax_left.set_xlabel('')
    ax_right.set_xlabel('percentage of total reads', fontsize=10)
    ax_right.xaxis.set_label_coords(0.0, -0.08)
    ax_left.xaxis.set_major_locator(mticker.MultipleLocator(2))
    ax_right.xaxis.set_major_locator(mticker.MultipleLocator(10))

    # Cluster panel
    ax_clust.set_xlabel('clusters', fontsize=10)
    ax_clust.spines['top'].set_visible(False)
    ax_clust.spines['right'].set_visible(False)
    ax_clust.tick_params(left=False)
    ax_clust.set_yticklabels([])

    # Grids
    for ax in (ax_left, ax_right, ax_clust):
        ax.xaxis.grid(True, alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)

    # Legend
    handles = [plt.Rectangle((0, 0), 1, 1, fc=condition_colors[c], alpha=0.85)
               for c in conditions]
    fig.legend(handles, [labels_fn(c) for c in conditions],
               loc='upper right', fontsize=9, frameon=False,
               bbox_to_anchor=(0.98, 0.98))

    fig.suptitle(title, fontsize=12, fontweight='bold', y=1.02)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def classify_reads_by_transcript_body(
    reads_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    start_column: str = 'alignment_start',
    end_column: str = 'alignment_end',
) -> pd.DataFrame:
    """
    Classify each read by the RNA biotype it overlaps most (majority bp overlap).

    Uses the full alignment span rather than a single position.  When a read
    overlaps multiple feature types the one with the most overlap bp wins;
    ties broken by BODY_CATEGORY_ORDER priority.  Antisense is checked when no
    same-strand feature is found.

    Args:
        reads_df:       DataFrame with chrom, strand, start_column, end_column
        annotation_df:  Gene annotation DataFrame with a feature_type column
        start_column:   0-based alignment start column
        end_column:     0-based alignment end (exclusive)

    Returns:
        DataFrame with added 'body_category' column
    """
    from intervaltree import IntervalTree

    ft_col = 'feature_type' if 'feature_type' in annotation_df.columns else 'feature'

    # Build per-(chrom, strand) interval trees storing body category
    same_trees: Dict[Tuple[str, str], IntervalTree] = defaultdict(IntervalTree)
    anti_trees: Dict[Tuple[str, str], IntervalTree] = defaultdict(IntervalTree)

    for _, row in annotation_df.iterrows():
        chrom  = row['chrom']
        start  = int(row['start'])
        end    = int(row['end'])
        strand = row.get('strand', '.')
        ft     = str(row.get(ft_col, '')).lower()

        body_cat = _BODY_FEATURE_MAP.get(ft)
        if body_cat is None:
            continue
        if start >= end:
            continue

        same_trees[(chrom, strand)][start:end] = body_cat
        anti = '-' if strand == '+' else '+'
        anti_trees[(chrom, anti)][start:end] = body_cat

    categories = []
    for _, row in reads_df.iterrows():
        chrom  = row['chrom']
        strand = row.get('strand', '+')
        r_start = int(row[start_column])
        r_end   = int(row[end_column])
        if r_start >= r_end:
            categories.append('Intergenic')
            continue

        # Accumulate overlap bp per body category on same strand
        overlap_bp: Dict[str, int] = defaultdict(int)
        for iv in same_trees[(chrom, strand)].overlap(r_start, r_end):
            ov = min(iv.end, r_end) - max(iv.begin, r_start)
            if ov > 0:
                overlap_bp[iv.data] += ov

        if overlap_bp:
            # Pick highest-overlap category; break ties by priority
            best = max(overlap_bp,
                       key=lambda c: (overlap_bp[c], -_BODY_PRIORITY.get(c, 999)))
            categories.append(best)
        else:
            # Check antisense
            anti = '-' if strand == '+' else '+'
            anti_hits = any(True for _ in anti_trees[(chrom, anti)].overlap(r_start, r_end))
            categories.append('Antisense' if anti_hits else 'Intergenic')

    out = reads_df.copy()
    out['body_category'] = categories
    return out


def _plot_combined_distribution_figure(
    distributions: Dict[str, Dict[str, int]],
    right_counts: Dict[str, Dict[str, int]],
    output_path: str,
    category_order: List[str],
    category_labels: Dict[str, str],
    category_colors: Dict[str, str],
    figure_title: str,
    barh_title: str,
    right_label: str = 'clusters',
    axis_break: float = 10.0,
    condition_labels: Optional[Dict[str, str]] = None,
    condition_colors: Optional[Dict[str, str]] = None,
) -> str:
    """
    Combined figure: horizontal bar chart (all conditions) on top,
    pie chart grid (one pie per condition) on bottom.  One PNG output.

    Args:
        distributions:    condition → {category → read count}
        right_counts:     condition → {category → count for right panel}
        output_path:      Output PNG path
        category_order:   Ordered list of internal category names
        category_labels:  Internal name → display label
        category_colors:  Internal name → hex color
        figure_title:     Overall figure suptitle
        barh_title:       Title for the bar chart panel
        right_label:      X-axis label for right (counts) panel
        axis_break:       % value where left barh panel clips
        condition_labels: Optional display names per condition
        condition_colors: Optional bar colors per condition

    Returns:
        Path to saved PNG
    """
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    conditions = list(distributions.keys())
    n_conds = len(conditions)
    labels_fn = lambda c: (condition_labels or {}).get(c, c)

    default_colors = ['#000000', '#AAAAAA', '#555555', '#CCCCCC', '#333333',
                      '#888888', '#111111', '#666666']
    if condition_colors is None:
        condition_colors = {c: default_colors[i % len(default_colors)]
                            for i, c in enumerate(conditions)}

    n_cats = len(category_order)
    y_pos  = np.arange(n_cats)

    # ── pct lookup ──────────────────────────────────────────────────────────
    pct: Dict[str, Dict[str, float]] = {}
    for cond, dist in distributions.items():
        total = sum(dist.values()) or 1
        pct[cond] = {cat: 100.0 * dist.get(cat, 0) / total for cat in category_order}

    # ── Pie data (each condition = all reps already merged) ─────────────────
    pie_sizes: Dict[str, List] = {}
    pie_colors: Dict[str, List] = {}
    pie_labels_list: Dict[str, List] = {}
    for cond, dist in distributions.items():
        total = sum(dist.values()) or 1
        sz, cols, lbls = [], [], []
        for cat in category_order:
            v = dist.get(cat, 0)
            if v > 0:
                sz.append(v)
                cols.append(category_colors.get(cat, '#888888'))
                lbls.append(f"{category_labels.get(cat, cat)}\n{100*v/total:.1f}%")
        pie_sizes[cond]  = sz
        pie_colors[cond] = cols
        pie_labels_list[cond] = lbls

    # ── Figure layout ───────────────────────────────────────────────────────
    # Top: barh panel (height ~ n_cats * 0.45 in, min 4)
    # Bottom: pie grid (height ~ 3.5 in per row)
    ncols_pie = min(n_conds, 5)
    nrows_pie = (n_conds + ncols_pie - 1) // ncols_pie
    barh_h = max(4.0, n_cats * 0.48)
    pie_h  = 3.6 * nrows_pie
    fig_h  = barh_h + pie_h + 0.8   # 0.8 spacing
    fig_w  = max(11, ncols_pie * 3.2)

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs_top = GridSpec(2, 1, figure=fig,
                      height_ratios=[barh_h, pie_h],
                      hspace=0.35)

    # ── Barh panel (broken axis) ─────────────────────────────────────────
    gs_barh = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs_top[0],
                                      width_ratios=[axis_break, 100 - axis_break, 55],
                                      wspace=0.06)
    ax_l  = fig.add_subplot(gs_barh[0])
    ax_r  = fig.add_subplot(gs_barh[1], sharey=ax_l)
    ax_rc = fig.add_subplot(gs_barh[2], sharey=ax_l)

    bar_h = 0.72 / max(n_conds, 1)
    for i, cond in enumerate(conditions):
        offset = (i - (n_conds - 1) / 2) * bar_h
        col    = condition_colors[cond]
        lbl    = labels_fn(cond)
        vals   = [pct[cond][cat] for cat in category_order]
        ax_l.barh(y_pos + offset, vals, bar_h * 0.9, color=col, label=lbl, alpha=0.85)
        ax_r.barh(y_pos + offset, vals, bar_h * 0.9, color=col, alpha=0.85)
        rc = right_counts.get(cond, {})
        rvals = [rc.get(cat, 0) for cat in category_order]
        ax_rc.barh(y_pos + offset, rvals, bar_h * 0.9, color=col, alpha=0.85)

    ax_l.set_xlim(0, axis_break)
    ax_r.set_xlim(axis_break, 100)
    ax_l.spines['right'].set_visible(False)
    ax_r.spines['left'].set_visible(False)
    ax_r.tick_params(left=False)
    d = 0.015
    for ax, side in [(ax_l, 'right'), (ax_r, 'left')]:
        kw = dict(transform=ax.transAxes, color='k', clip_on=False, lw=1)
        x0 = 1 if side == 'right' else 0
        ax.plot((x0-d, x0+d), (-d, +d), **kw)
        ax.plot((x0-d, x0+d), (1-d, 1+d), **kw)

    ax_l.set_yticks(y_pos)
    ax_l.set_yticklabels([category_labels.get(c, c) for c in category_order], fontsize=9)
    ax_l.invert_yaxis()
    ax_r.set_xlabel('percentage of total reads', fontsize=9)
    ax_rc.set_xlabel(right_label, fontsize=9)
    ax_l.xaxis.set_major_locator(mticker.MultipleLocator(2))
    ax_r.xaxis.set_major_locator(mticker.MultipleLocator(10))
    ax_rc.tick_params(left=False)
    ax_rc.set_yticklabels([])
    for ax in (ax_l, ax_r, ax_rc):
        ax.xaxis.grid(True, alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    ax_l.set_title(barh_title, fontsize=10, fontweight='bold', loc='left')
    handles = [plt.Rectangle((0,0),1,1, fc=condition_colors[c], alpha=0.85)
               for c in conditions]
    ax_r.legend(handles, [labels_fn(c) for c in conditions],
                fontsize=8, frameon=False, loc='upper right')

    # ── Pie grid ─────────────────────────────────────────────────────────
    gs_pies = GridSpecFromSubplotSpec(nrows_pie, ncols_pie,
                                      subplot_spec=gs_top[1],
                                      hspace=0.15, wspace=0.05)
    for idx, cond in enumerate(conditions):
        row, col_idx = divmod(idx, ncols_pie)
        ax_pie = fig.add_subplot(gs_pies[row, col_idx])
        sz   = pie_sizes[cond]
        cols = pie_colors[cond]
        if sz:
            ax_pie.pie(sz, colors=cols, startangle=90,
                       autopct=lambda p: f'{p:.0f}%' if p >= 5 else '',
                       pctdistance=0.72,
                       textprops={'fontsize': 7, 'fontweight': 'bold'})
        lbl = labels_fn(cond)
        n   = sum(distributions[cond].values())
        ax_pie.set_title(f"{lbl}\n(n={n:,})", fontsize=8, fontweight='bold')

    # Hide unused pie slots
    for idx in range(n_conds, nrows_pie * ncols_pie):
        row, col_idx = divmod(idx, ncols_pie)
        fig.add_subplot(gs_pies[row, col_idx]).set_visible(False)

    # Shared pie legend (bottom of figure)
    legend_patches = [
        plt.Rectangle((0,0),1,1, fc=category_colors.get(c,'#888'), label=category_labels.get(c,c))
        for c in category_order if any(distributions[cond].get(c,0) > 0 for cond in conditions)
    ]
    fig.legend(handles=legend_patches, loc='lower center',
               ncol=min(len(legend_patches), 6), fontsize=7, frameon=False,
               bbox_to_anchor=(0.5, 0.0))

    fig.suptitle(figure_title, fontsize=12, fontweight='bold', y=1.01)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def run_3prime_distribution_analysis(
    positions_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    output_dir: str,
    sample_column: str = 'sample',
    position_column: str = 'corrected_position',
    count_column: Optional[str] = None,
    condition_labels: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """
    3' end distribution analysis.

    Classifies each corrected 3' end position by genomic feature
    (UTR3, snoRNA±300bp, CUT, SUT_XUT, UTR5_CDS, Antisense, Intergenic)
    and writes a combined figure (barh + pie grid) for all conditions.

    Replicates are merged by condition before plotting.

    Args:
        positions_df:     DataFrame with per-read 3' end positions
        annotation_df:    Gene annotation with feature_type column
        output_dir:       Directory for output files
        sample_column:    Column identifying the sample
        position_column:  Column holding corrected 3' end position
        count_column:     Optional pre-aggregated count column
        condition_labels: Optional display names per condition

    Returns:
        Dict mapping output_name → file path
    """
    from ..analyze.deseq2 import extract_condition_from_sample

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_files: Dict[str, str] = {}

    df = positions_df.copy()
    df['condition'] = df[sample_column].apply(extract_condition_from_sample)
    conditions = sorted(df['condition'].unique())

    n = len(df)
    if n > 1_000_000:
        print(f"    Sampling 1,000,000/{n:,} positions for 3' end classification...")
        df = df.sample(n=1_000_000, random_state=42)

    print(f"    Classifying 3' end positions ({len(conditions)} conditions)...")
    classified = classify_positions_by_region(df, annotation_df, position_column)

    distributions: Dict[str, Dict[str, int]] = {}
    cluster_counts: Dict[str, Dict[str, int]] = {}
    for cond in conditions:
        sub = classified[classified['condition'] == cond]
        distributions[cond]  = calculate_genomic_distribution(sub, 'genomic_region', count_column)
        cluster_counts[cond] = count_clusters_by_region(sub, 'genomic_region', position_column)
        total = sum(distributions[cond].values())
        ncl   = sum(cluster_counts[cond].values())
        lbl   = (condition_labels or {}).get(cond, cond)
        print(f"      {lbl}: {total:,} reads, {ncl:,} unique clusters")

    # Combined figure
    fig_path = output_dir / 'genomic_distribution_3prime.png'
    _plot_combined_distribution_figure(
        distributions, cluster_counts, str(fig_path),
        category_order=CATEGORY_ORDER,
        category_labels=CATEGORY_LABELS,
        category_colors=CATEGORY_COLORS,
        figure_title="3\u2019 End Distribution by Genomic Category",
        barh_title="% of reads per category (all conditions)",
        right_label='clusters',
        axis_break=10.0,
        condition_labels=condition_labels,
    )
    output_files['figure'] = str(fig_path)
    print(f"    Combined figure: {fig_path.name}")

    # Summary TSV
    rows = []
    for cond in conditions:
        total = sum(distributions[cond].values()) or 1
        for cat in CATEGORY_ORDER:
            rows.append({
                'condition': cond,
                'category': cat,
                'display_label': CATEGORY_LABELS[cat],
                'reads': distributions[cond].get(cat, 0),
                'reads_pct': 100 * distributions[cond].get(cat, 0) / total,
                'clusters': cluster_counts[cond].get(cat, 0),
            })
    summary_path = output_dir / 'genomic_distribution_3prime_summary.tsv'
    pd.DataFrame(rows).to_csv(summary_path, sep='\t', index=False)
    output_files['summary'] = str(summary_path)
    return output_files


def run_transcript_body_distribution_analysis(
    reads_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    output_dir: str,
    sample_column: str = 'sample',
    start_column: str = 'alignment_start',
    end_column: str = 'alignment_end',
    condition_labels: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """
    Transcript body distribution analysis.

    Assigns each read to an RNA biotype (protein_coding, CUT, SUT_XUT,
    snoRNA, tRNA, rRNA, LTR, pseudogene, Antisense, Intergenic) based on
    majority bp overlap of the alignment span with annotated features.

    Replicates are merged by condition before plotting.

    Args:
        reads_df:         DataFrame with per-read alignment coordinates
        annotation_df:    Gene annotation with feature_type column
        output_dir:       Directory for output files
        sample_column:    Column identifying the sample
        start_column:     0-based alignment start column
        end_column:       0-based alignment end column (exclusive)
        condition_labels: Optional display names per condition

    Returns:
        Dict mapping output_name → file path
    """
    from ..analyze.deseq2 import extract_condition_from_sample

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_files: Dict[str, str] = {}

    df = reads_df.copy()
    df['condition'] = df[sample_column].apply(extract_condition_from_sample)
    conditions = sorted(df['condition'].unique())

    n = len(df)
    if n > 500_000:
        print(f"    Sampling 500,000/{n:,} reads for body classification...")
        df = df.sample(n=500_000, random_state=42)

    print(f"    Classifying transcript bodies ({len(conditions)} conditions)...")
    classified = classify_reads_by_transcript_body(df, annotation_df, start_column, end_column)

    distributions: Dict[str, Dict[str, int]] = {}
    for cond in conditions:
        sub = classified[classified['condition'] == cond]
        distributions[cond] = sub['body_category'].value_counts().to_dict()
        total = sum(distributions[cond].values())
        lbl   = (condition_labels or {}).get(cond, cond)
        print(f"      {lbl}: {total:,} reads classified")

    # Right panel = raw read counts per category (no "cluster" concept here)
    fig_path = output_dir / 'genomic_distribution_body.png'
    _plot_combined_distribution_figure(
        distributions, distributions, str(fig_path),
        category_order=BODY_CATEGORY_ORDER,
        category_labels=BODY_CATEGORY_LABELS,
        category_colors=BODY_CATEGORY_COLORS,
        figure_title="Transcript Body Distribution by RNA Biotype",
        barh_title="% of reads per biotype (all conditions)",
        right_label='reads',
        axis_break=10.0,
        condition_labels=condition_labels,
    )
    output_files['figure'] = str(fig_path)
    print(f"    Combined figure: {fig_path.name}")

    # Summary TSV
    rows = []
    for cond in conditions:
        total = sum(distributions[cond].values()) or 1
        for cat in BODY_CATEGORY_ORDER:
            rows.append({
                'condition': cond,
                'category': cat,
                'display_label': BODY_CATEGORY_LABELS[cat],
                'reads': distributions[cond].get(cat, 0),
                'reads_pct': 100 * distributions[cond].get(cat, 0) / total,
            })
    summary_path = output_dir / 'genomic_distribution_body_summary.tsv'
    pd.DataFrame(rows).to_csv(summary_path, sep='\t', index=False)
    output_files['summary'] = str(summary_path)
    return output_files


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

    # Classify positions (sample if very large)
    n_positions = len(positions_df)
    if n_positions > 1_000_000:
        print(f"    Sampling {1_000_000:,} positions for classification...")
        sample_df = positions_df.sample(n=1_000_000, random_state=42)
    else:
        sample_df = positions_df

    print(f"    Classifying positions by genomic region...")
    classified_df = classify_positions_by_region(
        sample_df, annotation_df, position_column
    )

    # Calculate per-condition read distributions and cluster counts
    distributions: Dict[str, Dict[str, int]] = {}
    all_cluster_counts: Dict[str, Dict[str, int]] = {}

    for condition in conditions:
        cond_df = classified_df[classified_df['condition'] == condition]
        distributions[condition] = calculate_genomic_distribution(
            cond_df, 'genomic_region', count_column
        )
        all_cluster_counts[condition] = count_clusters_by_region(
            cond_df, 'genomic_region', position_column
        )

        total = sum(distributions[condition].values())
        n_clust = sum(all_cluster_counts[condition].values())
        label = condition_labels.get(condition, condition) if condition_labels else condition
        print(f"    {label}: {total:,} reads, {n_clust:,} clusters classified")

    # --- Primary output: horizontal bar chart (panels B/C style) ---
    barh_path = output_dir / 'genomic_distribution.png'
    plot_genomic_distribution_barh(
        distributions,
        all_cluster_counts,
        str(barh_path),
        title="3' End Distribution by Genomic Category",
        condition_labels=condition_labels,
    )
    output_files['barh'] = str(barh_path)
    print(f"    Horizontal bar chart saved: {barh_path.name}")

    # --- Legacy outputs (pie charts + vertical bar) kept for compatibility ---
    for condition in conditions:
        total = sum(distributions[condition].values())
        pie_path = output_dir / f'genomic_distribution_pie_{condition}.png'
        label = condition_labels.get(condition, condition) if condition_labels else condition
        plot_genomic_distribution_pie(
            distributions[condition],
            str(pie_path),
            title=f"3' End Distribution: {label}",
        )
        output_files[f'pie_{condition}'] = str(pie_path)

    if len(conditions) > 1:
        grid_path = output_dir / 'genomic_distribution_pie_grid.png'
        plot_genomic_distribution_pie_grid(
            distributions, str(grid_path),
            title="3' End Distribution by Condition",
            condition_labels=condition_labels, wt_first=True,
        )
        output_files['pie_grid'] = str(grid_path)

        comparison_path = output_dir / 'genomic_distribution_comparison.png'
        plot_genomic_distribution_comparison(
            distributions, str(comparison_path),
            title="3' End Genomic Distribution Comparison",
            condition_labels=condition_labels,
        )
        output_files['comparison'] = str(comparison_path)

    # Summary table (reads + clusters per condition × category)
    summary_rows = []
    for condition in conditions:
        total_reads = sum(distributions[condition].values())
        for cat in CATEGORY_ORDER:
            read_count = distributions[condition].get(cat, 0)
            clust_count = all_cluster_counts[condition].get(cat, 0)
            summary_rows.append({
                'condition': condition,
                'category': cat,
                'display_label': CATEGORY_LABELS[cat],
                'reads': read_count,
                'reads_pct': 100 * read_count / total_reads if total_reads > 0 else 0,
                'clusters': clust_count,
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_dir / 'genomic_distribution_summary.tsv'
    summary_df.to_csv(summary_path, sep='\t', index=False)
    output_files['summary'] = str(summary_path)

    return output_files
