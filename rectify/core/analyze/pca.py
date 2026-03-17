#!/usr/bin/env python3
"""
PCA Analysis Module

Principal Component Analysis for sample clustering and QC.

Author: Kevin R. Roy
Date: 2026-03-17
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd

# Try to import sklearn
try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

# Try to import matplotlib
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def run_pca_analysis(
    count_matrix: pd.DataFrame,
    n_components: int = 10,
    normalize: bool = True,
    log_transform: bool = True,
    min_variance: float = 0.0,
) -> Dict:
    """
    Run PCA on count matrix.

    Args:
        count_matrix: Features (rows) × Samples (columns) count matrix
        n_components: Number of principal components to compute
        normalize: Whether to z-score normalize features
        log_transform: Whether to log2(x+1) transform before PCA
        min_variance: Minimum variance to include a feature

    Returns:
        Dict with:
            - pca_coords: DataFrame with PC coordinates (samples × PCs)
            - variance_explained: Array of variance explained per PC
            - variance_ratio: Array of variance ratio per PC
            - loadings: DataFrame with feature loadings (features × PCs)
            - n_features: Number of features used
    """
    if not SKLEARN_AVAILABLE:
        raise ImportError(
            "scikit-learn is required for PCA analysis. "
            "Install with: pip install scikit-learn"
        )

    if count_matrix.empty:
        return {
            'pca_coords': pd.DataFrame(),
            'variance_explained': np.array([]),
            'variance_ratio': np.array([]),
            'loadings': pd.DataFrame(),
            'n_features': 0,
        }

    # Prepare data matrix (features × samples -> samples × features for PCA)
    data = count_matrix.copy()

    # Log transform
    if log_transform:
        data = np.log2(data + 1)

    # Filter by variance
    if min_variance > 0:
        feature_variance = data.var(axis=1)
        data = data[feature_variance > min_variance]

    if data.empty or len(data) < 2:
        return {
            'pca_coords': pd.DataFrame(),
            'variance_explained': np.array([]),
            'variance_ratio': np.array([]),
            'loadings': pd.DataFrame(),
            'n_features': 0,
        }

    n_features = len(data)

    # Transpose: samples as rows, features as columns
    X = data.T

    # Normalize
    if normalize:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
    else:
        X_scaled = X.values

    # Handle NaN/Inf
    X_scaled = np.nan_to_num(X_scaled, nan=0.0, posinf=0.0, neginf=0.0)

    # Run PCA
    n_components = min(n_components, min(X_scaled.shape))
    pca = PCA(n_components=n_components)
    pca_coords = pca.fit_transform(X_scaled)

    # Create results
    pc_names = [f'PC{i+1}' for i in range(n_components)]

    pca_coords_df = pd.DataFrame(
        pca_coords,
        index=X.index,
        columns=pc_names,
    )

    loadings_df = pd.DataFrame(
        pca.components_.T,
        index=data.index,
        columns=pc_names,
    )

    return {
        'pca_coords': pca_coords_df,
        'variance_explained': pca.explained_variance_,
        'variance_ratio': pca.explained_variance_ratio_,
        'loadings': loadings_df,
        'n_features': n_features,
    }


def plot_pca(
    pca_results: Dict,
    sample_metadata: Optional[pd.DataFrame] = None,
    color_by: Optional[str] = None,
    pc_x: int = 1,
    pc_y: int = 2,
    figsize: Tuple[int, int] = (10, 8),
    output_path: Optional[str] = None,
    title: Optional[str] = None,
    show_labels: bool = True,
    marker_size: int = 200,
    alpha: float = 0.8,
) -> Optional[plt.Figure]:
    """
    Create PCA scatter plot.

    Args:
        pca_results: Results from run_pca_analysis()
        sample_metadata: Optional metadata for coloring/labeling
        color_by: Column in metadata to color by
        pc_x: PC for x-axis (1-indexed)
        pc_y: PC for y-axis (1-indexed)
        figsize: Figure size (width, height)
        output_path: Path to save figure (optional)
        title: Plot title
        show_labels: Whether to show sample labels
        marker_size: Size of scatter markers
        alpha: Transparency of markers

    Returns:
        matplotlib Figure (or None if matplotlib unavailable)
    """
    if not MATPLOTLIB_AVAILABLE:
        print("Warning: matplotlib not available for plotting")
        return None

    pca_coords = pca_results['pca_coords']
    variance_ratio = pca_results['variance_ratio']

    if pca_coords.empty:
        print("Warning: No PCA results to plot")
        return None

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Get PC columns
    pc_x_col = f'PC{pc_x}'
    pc_y_col = f'PC{pc_y}'

    if pc_x_col not in pca_coords.columns or pc_y_col not in pca_coords.columns:
        print(f"Warning: PC{pc_x} or PC{pc_y} not available")
        return None

    x = pca_coords[pc_x_col]
    y = pca_coords[pc_y_col]

    # Determine colors
    if sample_metadata is not None and color_by is not None:
        # Match samples
        samples = pca_coords.index
        colors_series = sample_metadata.loc[
            sample_metadata.index.isin(samples), color_by
        ]

        # Create color mapping
        unique_values = colors_series.unique()
        color_map = plt.cm.get_cmap('tab10')
        value_to_color = {v: color_map(i % 10) for i, v in enumerate(unique_values)}

        colors = [value_to_color.get(colors_series.get(s, None), 'gray')
                  for s in samples]

        # Create scatter with legend
        for value in unique_values:
            mask = [colors_series.get(s, None) == value for s in samples]
            ax.scatter(
                x[mask], y[mask],
                c=[value_to_color[value]],
                label=value,
                s=marker_size,
                alpha=alpha,
                edgecolors='black',
                linewidths=1,
            )
        ax.legend(title=color_by, loc='best')
    else:
        ax.scatter(
            x, y,
            s=marker_size,
            alpha=alpha,
            edgecolors='black',
            linewidths=1,
        )

    # Add sample labels
    if show_labels:
        for sample in pca_coords.index:
            ax.annotate(
                sample,
                (pca_coords.loc[sample, pc_x_col],
                 pca_coords.loc[sample, pc_y_col]),
                fontsize=8,
                alpha=0.7,
            )

    # Labels
    var_x = variance_ratio[pc_x - 1] * 100 if pc_x - 1 < len(variance_ratio) else 0
    var_y = variance_ratio[pc_y - 1] * 100 if pc_y - 1 < len(variance_ratio) else 0

    ax.set_xlabel(f'PC{pc_x} ({var_x:.1f}% variance)')
    ax.set_ylabel(f'PC{pc_y} ({var_y:.1f}% variance)')

    if title:
        ax.set_title(title)
    else:
        ax.set_title(f'PCA: {pca_results["n_features"]:,} features')

    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save if requested
    if output_path:
        plt.savefig(output_path, dpi=200, bbox_inches='tight')
        print(f"Saved PCA plot to {output_path}")

    return fig


def plot_variance_explained(
    pca_results: Dict,
    n_components: int = 10,
    figsize: Tuple[int, int] = (10, 5),
    output_path: Optional[str] = None,
) -> Optional[plt.Figure]:
    """
    Plot variance explained by each PC (scree plot).

    Args:
        pca_results: Results from run_pca_analysis()
        n_components: Number of components to show
        figsize: Figure size
        output_path: Path to save figure

    Returns:
        matplotlib Figure
    """
    if not MATPLOTLIB_AVAILABLE:
        return None

    variance_ratio = pca_results['variance_ratio']

    if len(variance_ratio) == 0:
        return None

    n_components = min(n_components, len(variance_ratio))
    variance_ratio = variance_ratio[:n_components]

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Individual variance
    axes[0].bar(range(1, n_components + 1), variance_ratio * 100)
    axes[0].set_xlabel('Principal Component')
    axes[0].set_ylabel('Variance Explained (%)')
    axes[0].set_title('Individual Variance')

    # Cumulative variance
    cumulative = np.cumsum(variance_ratio) * 100
    axes[1].plot(range(1, n_components + 1), cumulative, 'o-')
    axes[1].axhline(y=80, color='r', linestyle='--', alpha=0.5, label='80%')
    axes[1].axhline(y=90, color='orange', linestyle='--', alpha=0.5, label='90%')
    axes[1].set_xlabel('Number of Components')
    axes[1].set_ylabel('Cumulative Variance (%)')
    axes[1].set_title('Cumulative Variance')
    axes[1].legend()

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=200, bbox_inches='tight')

    return fig
