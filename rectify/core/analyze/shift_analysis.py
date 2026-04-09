#!/usr/bin/env python3
"""
Cluster Distribution Shift Analysis Module

Analyzes changes in CPA cluster usage between conditions.

Author: Kevin R. Roy
Date: 2026-03-17
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd

# Try to import matplotlib
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.gridspec import GridSpec
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def analyze_cluster_shifts(
    count_matrix: pd.DataFrame,
    clusters_df: pd.DataFrame,
    condition_a: str,
    condition_b: str,
    sample_metadata: pd.DataFrame,
    min_total_counts: int = 20,
) -> pd.DataFrame:
    """
    Analyze shifts in cluster distribution between two conditions.

    For each gene, compares the distribution of 3' end usage across
    clusters between condition A (reference) and condition B (treatment).

    Args:
        count_matrix: Cluster × Sample count matrix
        clusters_df: Cluster definitions with gene_id
        condition_a: Reference condition name
        condition_b: Treatment condition name
        sample_metadata: Metadata with 'condition' column
        min_total_counts: Minimum total counts across both conditions

    Returns:
        DataFrame with per-gene shift analysis:
            - gene_id, gene_name
            - n_clusters: Number of clusters for this gene
            - major_cluster_a, major_cluster_b: Dominant cluster in each condition
            - shift_bp: Shift in major cluster position (bp)
            - shift_direction: 'upstream', 'downstream', or 'same'
            - distribution_divergence: Jensen-Shannon divergence
    """
    if count_matrix.empty or clusters_df.empty:
        return pd.DataFrame()

    # Get samples for each condition
    samples_a = sample_metadata[
        sample_metadata['condition'] == condition_a
    ].index.tolist()
    samples_b = sample_metadata[
        sample_metadata['condition'] == condition_b
    ].index.tolist()

    samples_a = [s for s in samples_a if s in count_matrix.columns]
    samples_b = [s for s in samples_b if s in count_matrix.columns]

    if not samples_a or not samples_b:
        return pd.DataFrame()

    # Calculate mean counts per condition
    counts_a = count_matrix[samples_a].mean(axis=1)
    counts_b = count_matrix[samples_b].mean(axis=1)

    # Add gene info to clusters
    cluster_info = clusters_df.set_index('cluster_id')

    results = []

    # Group clusters by gene
    for gene_id in clusters_df['gene_id'].dropna().unique():
        gene_clusters = clusters_df[clusters_df['gene_id'] == gene_id]

        if len(gene_clusters) < 1:
            continue

        cluster_ids = gene_clusters['cluster_id'].tolist()

        # Get counts for this gene's clusters
        gene_counts_a = counts_a.loc[counts_a.index.isin(cluster_ids)]
        gene_counts_b = counts_b.loc[counts_b.index.isin(cluster_ids)]

        total_a = gene_counts_a.sum()
        total_b = gene_counts_b.sum()

        # Skip low-count genes
        if total_a + total_b < min_total_counts:
            continue

        # Normalize to fractions
        frac_a = gene_counts_a / total_a if total_a > 0 else gene_counts_a * 0
        frac_b = gene_counts_b / total_b if total_b > 0 else gene_counts_b * 0

        # Find major cluster in each condition
        major_cluster_a = frac_a.idxmax() if total_a > 0 else None
        major_cluster_b = frac_b.idxmax() if total_b > 0 else None

        # Calculate shift
        if major_cluster_a is not None and major_cluster_b is not None:
            pos_a = cluster_info.loc[major_cluster_a, 'modal_position']
            pos_b = cluster_info.loc[major_cluster_b, 'modal_position']
            strand = cluster_info.loc[major_cluster_a, 'strand']

            # Calculate strand-aware shift
            if strand == '+':
                shift_bp = pos_b - pos_a  # Positive = downstream
            else:
                shift_bp = pos_a - pos_b  # For minus strand, lower coord = downstream

            if shift_bp > 2:
                shift_direction = 'downstream'
            elif shift_bp < -2:
                shift_direction = 'upstream'
            else:
                shift_direction = 'same'
        else:
            shift_bp = 0
            shift_direction = 'unknown'

        # Calculate distribution divergence (Jensen-Shannon)
        divergence = _jensen_shannon_divergence(frac_a.values, frac_b.values)

        # Get gene name
        gene_name = gene_clusters['gene_name'].iloc[0] if 'gene_name' in gene_clusters.columns else gene_id

        results.append({
            'gene_id': gene_id,
            'gene_name': gene_name,
            'n_clusters': len(gene_clusters),
            'counts_a': total_a,
            'counts_b': total_b,
            'major_cluster_a': major_cluster_a,
            'major_cluster_b': major_cluster_b,
            'major_frac_a': frac_a.max() if total_a > 0 else 0,
            'major_frac_b': frac_b.max() if total_b > 0 else 0,
            'shift_bp': shift_bp,
            'shift_direction': shift_direction,
            'distribution_divergence': divergence,
            'same_major': major_cluster_a == major_cluster_b,
        })

    return pd.DataFrame(results)


def _jensen_shannon_divergence(p: np.ndarray, q: np.ndarray) -> float:
    """Calculate Jensen-Shannon divergence between two distributions."""
    from scipy.spatial.distance import jensenshannon
    p = np.array(p, dtype=float)
    q = np.array(q, dtype=float)
    # jensenshannon returns sqrt(JSD); square and use base=2 to get JSD in [0, 1]
    return float(jensenshannon(p, q, base=2) ** 2)


def get_top_shifted_genes(
    shift_results: pd.DataFrame,
    n_top: int = 20,
    min_divergence: float = 0.1,
    sort_by: str = 'distribution_divergence',
) -> pd.DataFrame:
    """
    Get top genes with largest distribution shifts.

    Args:
        shift_results: Results from analyze_cluster_shifts()
        n_top: Number of top genes to return
        min_divergence: Minimum JS divergence to include
        sort_by: Column to sort by

    Returns:
        Top shifted genes DataFrame
    """
    if shift_results.empty:
        return shift_results

    # Filter by minimum divergence
    filtered = shift_results[shift_results['distribution_divergence'] >= min_divergence]

    # Sort and take top N
    sorted_df = filtered.sort_values(sort_by, ascending=False)

    return sorted_df.head(n_top)


def plot_gene_browser(
    gene_id: str,
    count_matrix: pd.DataFrame,
    clusters_df: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    conditions: List[str],
    annotation_df: Optional[pd.DataFrame] = None,
    figsize: Tuple[int, int] = (14, 10),
    output_path: Optional[str] = None,
    use_common_names: bool = True,
) -> Optional[plt.Figure]:
    """
    Create genome browser-style plot showing CPA distribution for a gene.

    Follows locus_plotting conventions:
    - Gene track at bottom with CDS boxes as pentagon arrows
    - Strand direction indicated by arrow shape
    - Consistent color scheme

    Args:
        gene_id: Gene identifier
        count_matrix: Cluster × Sample count matrix
        clusters_df: Cluster definitions
        sample_metadata: Sample metadata with 'condition' column
        conditions: List of conditions to show
        annotation_df: Optional gene annotation DataFrame with CDS coordinates
        figsize: Figure size
        output_path: Path to save figure
        use_common_names: Use common gene names when available

    Returns:
        matplotlib Figure
    """
    if not MATPLOTLIB_AVAILABLE:
        return None

    # Get clusters for this gene
    gene_clusters = clusters_df[clusters_df['gene_id'] == gene_id].copy()

    if gene_clusters.empty:
        print(f"No clusters found for gene {gene_id}")
        return None

    gene_clusters = gene_clusters.sort_values('modal_position')
    cluster_ids = gene_clusters['cluster_id'].tolist()
    positions = gene_clusters['modal_position'].values
    strand = gene_clusters['strand'].iloc[0]
    chrom = gene_clusters['chrom'].iloc[0]

    # Determine plot range (extend beyond clusters to show gene context)
    cluster_span = positions.max() - positions.min()
    padding = max(200, cluster_span * 0.2)
    x_min = positions.min() - padding
    x_max = positions.max() + padding

    # Get gene name
    gene_name = gene_clusters['gene_name'].iloc[0] if 'gene_name' in gene_clusters.columns else gene_id
    display_name = gene_name if use_common_names and gene_name else gene_id

    # Set up figure with gene track
    n_condition_panels = len(conditions)
    height_ratios = [3] * n_condition_panels + [1]  # Gene track is shorter

    fig, axes = plt.subplots(
        n_condition_panels + 1, 1,
        figsize=figsize,
        sharex=True,
        gridspec_kw={'height_ratios': height_ratios},
    )

    # Condition colors (consistent with locus_plotting)
    CONDITION_COLORS = {
        0: '#3498db',  # Blue (reference)
        1: '#e74c3c',  # Red (treatment)
        2: '#2ecc71',  # Green
        3: '#9b59b6',  # Purple
        4: '#f39c12',  # Orange
    }

    # Plot each condition's CPA distribution
    for idx, condition in enumerate(conditions):
        ax = axes[idx]

        # Get samples for this condition
        samples = sample_metadata[
            sample_metadata['condition'] == condition
        ].index.tolist()
        samples = [s for s in samples if s in count_matrix.columns]

        if not samples:
            ax.text(0.5, 0.5, f'No samples for {condition}',
                   transform=ax.transAxes, ha='center', va='center')
            ax.set_xlim(x_min, x_max)
            continue

        # Calculate mean counts per cluster
        mean_counts = count_matrix.loc[cluster_ids, samples].mean(axis=1)
        total = mean_counts.sum()
        fractions = mean_counts / total if total > 0 else mean_counts * 0

        # Plot bars at cluster positions
        bar_width = max(10, (x_max - x_min) / 80)
        color = CONDITION_COLORS.get(idx, '#95a5a6')

        ax.bar(
            positions,
            fractions.values,
            width=bar_width,
            color=color,
            alpha=0.85,
            edgecolor='black',
            linewidth=0.5,
            label=condition,
        )

        # Styling
        ax.set_ylabel(f'{condition}\nFraction', fontsize=10)
        max_frac = fractions.max()
        ax.set_ylim(0, max_frac * 1.2 if max_frac > 0 else 1)
        ax.set_xlim(x_min, x_max)

        # Add gridlines
        ax.grid(True, axis='y', alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Add condition label in corner
        ax.text(0.02, 0.95, condition, transform=ax.transAxes,
               fontsize=11, fontweight='bold', va='top',
               color=color)

    # Gene track (bottom panel)
    ax_gene = axes[-1]
    _draw_gene_track(
        ax_gene, gene_id, gene_name, strand, chrom,
        x_min, x_max, positions,
        annotation_df=annotation_df,
        use_common_names=use_common_names,
    )

    # Title
    fig.suptitle(f'{display_name} ({gene_id}) - CPA Site Distribution', fontsize=14, fontweight='bold')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=200, bbox_inches='tight')
        print(f"Saved browser plot to {output_path}")

    return fig


def _draw_gene_track(
    ax,
    gene_id: str,
    gene_name: str,
    strand: str,
    chrom: str,
    region_start: int,
    region_end: int,
    cpa_positions: np.ndarray,
    annotation_df: Optional[pd.DataFrame] = None,
    use_common_names: bool = True,
) -> None:
    """
    Draw gene structure track with CDS/UTR distinction.

    Follows locus_plotting conventions:
    - Pentagon/arrow shapes indicating strand direction
    - CDS shown as filled boxes (full height)
    - UTR shown as thinner lines (half height)
    - Target gene highlighted in red
    - Neighbor genes in gray

    Args:
        ax: Matplotlib axes
        gene_id: Target gene ID
        gene_name: Target gene name
        strand: Strand of target gene
        chrom: Chromosome
        region_start: Start of region to display
        region_end: End of region to display
        cpa_positions: CPA cluster positions (for marking)
        annotation_df: Gene annotation DataFrame
        use_common_names: Use common gene names
    """
    from matplotlib.patches import Polygon, Rectangle

    ax.set_xlim(region_start, region_end)
    ax.set_ylim(-0.5, 1.5)

    # Colors
    TARGET_COLOR = '#e74c3c'  # Red
    TARGET_UTR_COLOR = '#f1948a'  # Light red for UTR
    NEIGHBOR_COLOR = '#95a5a6'  # Gray
    NEIGHBOR_UTR_COLOR = '#d5d8dc'  # Light gray for UTR
    CPA_MARKER_COLOR = '#2c3e50'  # Dark blue-gray

    # Track which genes have been labeled
    genes_labeled = set()
    features_drawn = False

    if annotation_df is not None and not annotation_df.empty:
        # Filter to features in region
        region_features = annotation_df[
            (annotation_df['chrom'] == chrom) &
            (annotation_df['end'] >= region_start) &
            (annotation_df['start'] <= region_end)
        ].copy()

        if 'feature_type' in region_features.columns:
            # Draw UTRs first (thinner, behind CDS)
            utr_types = {'UTR', '5UTR', '3UTR', 'five_prime_UTR', 'three_prime_UTR'}
            utr_features = region_features[region_features['feature_type'].isin(utr_types)]

            for _, feature in utr_features.iterrows():
                start = max(feature['start'], region_start)
                end = min(feature['end'], region_end)
                feat_gene_id = feature.get('gene_id', '')
                feat_gene_name = feature.get('gene_name', feat_gene_id)

                is_target = (feat_gene_id == gene_id or
                            feat_gene_name == gene_name or
                            gene_id in str(feat_gene_id))

                utr_color = TARGET_UTR_COLOR if is_target else NEIGHBOR_UTR_COLOR
                alpha = 0.8 if is_target else 0.4

                # Draw UTR as thin rectangle (half height of CDS)
                rect = Rectangle(
                    (start, 0.35), end - start, 0.3,
                    facecolor=utr_color, edgecolor='black',
                    linewidth=0.5, alpha=alpha
                )
                ax.add_patch(rect)
                features_drawn = True

            # Draw CDS (full height arrows)
            cds_features = region_features[region_features['feature_type'] == 'CDS']

            for _, feature in cds_features.iterrows():
                start = max(feature['start'], region_start)
                end = min(feature['end'], region_end)
                feat_strand = feature.get('strand', strand)
                feat_gene_id = feature.get('gene_id', '')
                feat_gene_name = feature.get('gene_name', feat_gene_id)

                is_target = (feat_gene_id == gene_id or
                            feat_gene_name == gene_name or
                            gene_id in str(feat_gene_id))

                color = TARGET_COLOR if is_target else NEIGHBOR_COLOR
                alpha = 0.9 if is_target else 0.5
                linewidth = 1.5 if is_target else 0.8

                # Draw pentagon arrow for CDS
                _draw_gene_arrow(ax, start, end, feat_strand, 0.5, 0.35,
                               color=color, alpha=alpha, linewidth=linewidth)

                # Label (only once per gene)
                label_key = feat_gene_id or feat_gene_name
                if label_key and label_key not in genes_labeled:
                    label_x = (start + end) / 2
                    display = feat_gene_name if use_common_names and feat_gene_name else feat_gene_id
                    if display:
                        fontsize = 9 if is_target else 7
                        fontweight = 'bold' if is_target else 'normal'
                        ax.text(label_x, 0.5, display[:12], ha='center', va='center',
                               fontsize=fontsize, fontweight=fontweight, color='white')
                        genes_labeled.add(label_key)

                features_drawn = True

        else:
            # No feature_type column - treat all as CDS (legacy behavior)
            for _, feature in region_features.iterrows():
                start = max(feature['start'], region_start)
                end = min(feature['end'], region_end)
                feat_strand = feature.get('strand', strand)
                feat_gene_id = feature.get('gene_id', '')
                feat_gene_name = feature.get('gene_name', feat_gene_id)

                is_target = (feat_gene_id == gene_id or
                            feat_gene_name == gene_name or
                            gene_id in str(feat_gene_id))

                color = TARGET_COLOR if is_target else NEIGHBOR_COLOR
                alpha = 0.9 if is_target else 0.5
                linewidth = 1.5 if is_target else 0.8

                _draw_gene_arrow(ax, start, end, feat_strand, 0.5, 0.35,
                               color=color, alpha=alpha, linewidth=linewidth)

                label_key = feat_gene_id or feat_gene_name
                if label_key and label_key not in genes_labeled:
                    label_x = (start + end) / 2
                    display = feat_gene_name if use_common_names and feat_gene_name else feat_gene_id
                    if display:
                        fontsize = 9 if is_target else 7
                        fontweight = 'bold' if is_target else 'normal'
                        ax.text(label_x, 0.5, display[:12], ha='center', va='center',
                               fontsize=fontsize, fontweight=fontweight, color='white')
                        genes_labeled.add(label_key)

                features_drawn = True

    # If no annotation, draw a simple gene box for the target
    if not features_drawn:
        # Estimate gene position from CPA clusters
        gene_start = cpa_positions.min() - 500
        gene_end = cpa_positions.max() + 100
        _draw_gene_arrow(ax, gene_start, gene_end, strand, 0.5, 0.35,
                        color=TARGET_COLOR, alpha=0.9, linewidth=1.5)
        ax.text((gene_start + gene_end) / 2, 0.5,
               gene_name if gene_name else gene_id,
               ha='center', va='center', fontsize=9, fontweight='bold', color='white')

    # Mark CPA positions with small triangles
    for pos in cpa_positions:
        ax.plot(pos, 1.2, 'v', color=CPA_MARKER_COLOR, markersize=6)

    # Styling
    ax.set_yticks([])
    ax.set_ylabel('Gene', fontsize=10)
    ax.set_xlabel(f'Genomic Position (chr{chrom}, {strand} strand)', fontsize=10)
    ax.axhline(y=0, color='black', linewidth=0.5, alpha=0.3)

    # Add strand direction arrow
    arrow_x = region_start + (region_end - region_start) * 0.95
    if strand == '+':
        ax.annotate('', xy=(arrow_x, -0.3), xytext=(arrow_x - 50, -0.3),
                   arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
        ax.text(arrow_x - 25, -0.45, "5'→3'", ha='center', fontsize=8)
    else:
        ax.annotate('', xy=(arrow_x - 50, -0.3), xytext=(arrow_x, -0.3),
                   arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
        ax.text(arrow_x - 25, -0.45, "3'←5'", ha='center', fontsize=8)


def _draw_gene_arrow(
    ax,
    start: float,
    end: float,
    strand: str,
    y_center: float,
    height: float,
    color: str = '#95a5a6',
    alpha: float = 0.8,
    linewidth: float = 1.0,
) -> None:
    """
    Draw a gene as a pentagon arrow indicating strand direction.

    Args:
        ax: Matplotlib axes
        start: Start coordinate
        end: End coordinate
        strand: '+' or '-'
        y_center: Y position of center
        height: Height of the arrow
        color: Fill color
        alpha: Transparency
        linewidth: Edge line width
    """
    from matplotlib.patches import Polygon

    width = end - start
    arrow_tip = min(width * 0.15, 100)  # Arrow tip size

    y_bot = y_center - height
    y_mid = y_center
    y_top = y_center + height

    if strand == '+':
        vertices = [
            (start, y_bot),
            (start, y_top),
            (end - arrow_tip, y_top),
            (end, y_mid),
            (end - arrow_tip, y_bot),
        ]
    else:
        vertices = [
            (start + arrow_tip, y_top),
            (end, y_top),
            (end, y_bot),
            (start + arrow_tip, y_bot),
            (start, y_mid),
        ]

    arrow_patch = Polygon(
        vertices, closed=True,
        facecolor=color, edgecolor='black',
        linewidth=linewidth, alpha=alpha
    )
    ax.add_patch(arrow_patch)


def plot_shift_summary(
    shift_results: pd.DataFrame,
    figsize: Tuple[int, int] = (14, 10),
    output_path: Optional[str] = None,
) -> Optional[plt.Figure]:
    """
    Create summary plots for shift analysis.

    Creates a 2x2 panel:
    A. Shift direction distribution (bar chart)
    B. Shift magnitude histogram
    C. JS divergence vs shift magnitude scatter
    D. Top shifted genes table

    Args:
        shift_results: Results from analyze_cluster_shifts()
        figsize: Figure size
        output_path: Path to save figure

    Returns:
        matplotlib Figure
    """
    if not MATPLOTLIB_AVAILABLE:
        return None

    if shift_results.empty:
        return None

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 2, figure=fig)

    # A. Shift direction
    ax_a = fig.add_subplot(gs[0, 0])
    direction_counts = shift_results['shift_direction'].value_counts()
    colors = {'upstream': '#2ca02c', 'same': '#7f7f7f', 'downstream': '#d62728'}
    bars = ax_a.bar(
        direction_counts.index,
        direction_counts.values,
        color=[colors.get(d, 'gray') for d in direction_counts.index],
    )
    ax_a.set_ylabel('Number of Genes')
    ax_a.set_title('A. Major Site Shift Direction')

    # B. Shift magnitude histogram
    ax_b = fig.add_subplot(gs[0, 1])
    shifts = shift_results['shift_bp'].dropna()
    ax_b.hist(shifts, bins=30, color='steelblue', edgecolor='black', alpha=0.7)
    ax_b.axvline(0, color='red', linestyle='--', alpha=0.5)
    ax_b.set_xlabel('Shift (bp)')
    ax_b.set_ylabel('Number of Genes')
    ax_b.set_title('B. Shift Magnitude Distribution')

    # C. Divergence vs shift
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.scatter(
        shift_results['shift_bp'],
        shift_results['distribution_divergence'],
        alpha=0.5,
        c=[colors.get(d, 'gray') for d in shift_results['shift_direction']],
    )
    ax_c.set_xlabel('Shift (bp)')
    ax_c.set_ylabel('JS Divergence')
    ax_c.set_title('C. Distribution Change vs Position Shift')

    # D. Top genes table
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.axis('off')

    top_genes = get_top_shifted_genes(shift_results, n_top=10)
    if not top_genes.empty:
        table_data = top_genes[['gene_name', 'shift_bp', 'distribution_divergence']].values
        table = ax_d.table(
            cellText=[[r[0], f'{r[1]:.0f}', f'{r[2]:.3f}'] for r in table_data],
            colLabels=['Gene', 'Shift (bp)', 'JS Divergence'],
            loc='center',
            cellLoc='center',
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.2, 1.5)
    ax_d.set_title('D. Top 10 Shifted Genes')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=200, bbox_inches='tight')

    return fig
