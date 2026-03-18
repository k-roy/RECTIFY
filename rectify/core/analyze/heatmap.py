#!/usr/bin/env python3
"""
Heatmap and Sample Clustering Module

Creates clustered heatmaps for sample QC and visualization.

Author: Kevin R. Roy
Date: 2026-03-17
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd

# Try to import visualization libraries
try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    import seaborn as sns
    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False

try:
    from sklearn.preprocessing import StandardScaler
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

try:
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import pdist
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


def plot_sample_heatmap(
    count_matrix: pd.DataFrame,
    sample_metadata: Optional[pd.DataFrame] = None,
    color_by: Optional[str] = None,
    n_top_features: int = 1000,
    method: str = 'correlation',
    figsize: Tuple[int, int] = (12, 10),
    output_path: Optional[str] = None,
    title: str = 'Sample Correlation Heatmap',
    cmap: str = 'RdBu_r',
) -> Optional[plt.Figure]:
    """
    Create sample-sample correlation heatmap with hierarchical clustering.

    Args:
        count_matrix: Features (rows) × Samples (columns) count matrix
        sample_metadata: Optional metadata for annotation
        color_by: Column in metadata for row/column colors
        n_top_features: Number of top variable features to use
        method: Distance metric ('correlation', 'euclidean')
        figsize: Figure size
        output_path: Path to save figure
        title: Plot title
        cmap: Colormap name

    Returns:
        matplotlib Figure
    """
    if not SEABORN_AVAILABLE:
        print("Warning: seaborn required for clustermap")
        return None

    if count_matrix.empty:
        return None

    # Log transform and select top variable features
    data = np.log2(count_matrix + 1)

    if len(data) > n_top_features:
        variance = data.var(axis=1)
        top_features = variance.nlargest(n_top_features).index
        data = data.loc[top_features]

    # Calculate sample-sample correlation
    # count_matrix is Features × Samples, so data.corr() gives sample correlation
    if method == 'correlation':
        corr_matrix = data.corr()  # Correlation between columns (samples)
    else:
        # For other methods, use distance matrix between samples
        from scipy.spatial.distance import pdist, squareform
        distances = pdist(data.T, metric=method)  # data.T is Samples × Features
        corr_matrix = pd.DataFrame(
            1 - squareform(distances),
            index=data.columns,
            columns=data.columns,
        )

    # Create column colors if metadata provided
    col_colors = None
    if sample_metadata is not None and color_by is not None:
        unique_values = sample_metadata[color_by].unique()
        color_palette = sns.color_palette('tab10', len(unique_values))
        value_to_color = dict(zip(unique_values, color_palette))

        col_colors = sample_metadata.loc[
            corr_matrix.columns, color_by
        ].map(value_to_color)

    # Create clustermap
    g = sns.clustermap(
        corr_matrix,
        cmap=cmap,
        center=None if method != 'correlation' else None,
        vmin=0 if method == 'correlation' else None,
        vmax=1 if method == 'correlation' else None,
        col_colors=col_colors,
        row_colors=col_colors,
        figsize=figsize,
        dendrogram_ratio=0.15,
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
    )

    g.ax_heatmap.set_title(title)

    if output_path:
        plt.savefig(output_path, dpi=200, bbox_inches='tight')
        print(f"Saved sample heatmap to {output_path}")

    return g.fig


def plot_cluster_heatmap(
    count_matrix: pd.DataFrame,
    sample_metadata: Optional[pd.DataFrame] = None,
    color_by: Optional[str] = None,
    n_top_features: int = 500,
    normalize: str = 'zscore',
    cluster_rows: bool = True,
    cluster_cols: bool = False,
    figsize: Tuple[int, int] = (12, 14),
    output_path: Optional[str] = None,
    title: str = 'CPA Cluster Expression',
    cmap: str = 'RdBu_r',
    vmin: float = -2,
    vmax: float = 2,
) -> Optional[plt.Figure]:
    """
    Create expression heatmap of top variable features.

    Args:
        count_matrix: Features (rows) × Samples (columns) count matrix
        sample_metadata: Optional metadata for annotation
        color_by: Column in metadata for column colors
        n_top_features: Number of top variable features to show
        normalize: Normalization method ('zscore', 'log', 'none')
        cluster_rows: Whether to cluster rows (features)
        cluster_cols: Whether to cluster columns (samples)
        figsize: Figure size
        output_path: Path to save figure
        title: Plot title
        cmap: Colormap name
        vmin, vmax: Color scale limits

    Returns:
        matplotlib Figure
    """
    if not SEABORN_AVAILABLE:
        print("Warning: seaborn required for clustermap")
        return None

    if count_matrix.empty:
        return None

    # Prepare data
    data = count_matrix.copy()

    # Log transform first
    data = np.log2(data + 1)

    # Select top variable features
    if len(data) > n_top_features:
        variance = data.var(axis=1)
        top_features = variance.nlargest(n_top_features).index
        data = data.loc[top_features]

    # Normalize
    if normalize == 'zscore':
        if SKLEARN_AVAILABLE:
            scaler = StandardScaler()
            data = pd.DataFrame(
                scaler.fit_transform(data.T).T,
                index=data.index,
                columns=data.columns,
            )
        else:
            # Manual z-score
            data = (data.T - data.mean(axis=1)) / data.std(axis=1)
            data = data.T

    # Create column colors if metadata provided
    col_colors = None
    if sample_metadata is not None and color_by is not None:
        unique_values = sample_metadata[color_by].unique()
        color_palette = sns.color_palette('tab10', len(unique_values))
        value_to_color = dict(zip(unique_values, color_palette))

        col_colors = sample_metadata.loc[
            data.columns, color_by
        ].map(value_to_color)

    # Create clustermap
    g = sns.clustermap(
        data,
        cmap=cmap,
        center=0 if normalize == 'zscore' else None,
        vmin=vmin if normalize == 'zscore' else None,
        vmax=vmax if normalize == 'zscore' else None,
        col_colors=col_colors,
        row_cluster=cluster_rows,
        col_cluster=cluster_cols,
        figsize=figsize,
        dendrogram_ratio=(0.1, 0.15),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        yticklabels=False,  # Hide feature labels
        xticklabels=True,
    )

    g.ax_heatmap.set_title(title)

    if output_path:
        plt.savefig(output_path, dpi=200, bbox_inches='tight')
        print(f"Saved cluster heatmap to {output_path}")

    return g.fig


def plot_condition_dendrogram(
    deseq2_results: Dict[str, pd.DataFrame],
    metric: str = 'correlation',
    figsize: Tuple[int, int] = (10, 6),
    output_path: Optional[str] = None,
) -> Optional[plt.Figure]:
    """
    Plot hierarchical clustering dendrogram of conditions based on fold changes.

    Args:
        deseq2_results: Dict mapping condition -> DESeq2 results
        metric: Distance metric for clustering
        figsize: Figure size
        output_path: Path to save figure

    Returns:
        matplotlib Figure
    """
    if not SCIPY_AVAILABLE or not MATPLOTLIB_AVAILABLE:
        return None

    if len(deseq2_results) < 2:
        print("Warning: Need at least 2 conditions for dendrogram")
        return None

    # Build fold change matrix (genes × conditions)
    # Get common genes
    conditions = list(deseq2_results.keys())
    gene_sets = [set(df.index) for df in deseq2_results.values()]
    common_genes = list(set.intersection(*gene_sets))

    if len(common_genes) < 10:
        print("Warning: Not enough common genes for clustering")
        return None

    # Build matrix
    fc_matrix = np.zeros((len(common_genes), len(conditions)))
    for i, condition in enumerate(conditions):
        df = deseq2_results[condition]
        fc_matrix[:, i] = df.loc[common_genes, 'log2FoldChange'].values

    # Handle NaN
    fc_matrix = np.nan_to_num(fc_matrix, nan=0.0)

    # Calculate distances and linkage
    distances = pdist(fc_matrix.T, metric=metric)
    linkage_matrix = linkage(distances, method='average')

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    dendrogram(linkage_matrix, labels=conditions, ax=ax)

    ax.set_xlabel('Condition')
    ax.set_ylabel(f'{metric.capitalize()} Distance')
    ax.set_title('Condition Clustering by Log2 Fold Changes')

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=200, bbox_inches='tight')

    return fig
