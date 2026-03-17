#!/usr/bin/env python3
"""
GO Enrichment Analysis Module

Gene Ontology enrichment analysis using Fisher's exact test.

Author: Kevin R. Roy
Date: 2026-03-17
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd

# Try to import scipy
try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

# Try to import matplotlib
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def load_go_annotations(
    filepath: str,
    gene_col: str = 'gene_id',
    go_term_col: str = 'go_term',
    go_id_col: str = 'go_id',
    category_col: str = 'category',
) -> pd.DataFrame:
    """
    Load GO annotations from file.

    Expected format: TSV with columns for gene_id, go_term, go_id, category (C/F/P)

    Args:
        filepath: Path to GO annotation file
        gene_col: Column name for gene identifiers
        go_term_col: Column name for GO term descriptions
        go_id_col: Column name for GO IDs
        category_col: Column name for GO category (C/F/P)

    Returns:
        DataFrame with GO annotations
    """
    df = pd.read_csv(filepath, sep='\t')

    # Standardize column names
    column_mapping = {
        gene_col: 'gene_id',
        go_term_col: 'go_term',
        go_id_col: 'go_id',
        category_col: 'category',
    }

    df = df.rename(columns={k: v for k, v in column_mapping.items() if k in df.columns})

    return df


def run_go_enrichment(
    query_genes: List[str],
    go_annotations: pd.DataFrame,
    background_genes: Optional[List[str]] = None,
    categories: Optional[List[str]] = None,
    min_genes_per_term: int = 3,
    max_genes_per_term: int = 500,
) -> pd.DataFrame:
    """
    Run GO enrichment analysis using Fisher's exact test.

    Args:
        query_genes: List of gene IDs to test for enrichment
        go_annotations: DataFrame with GO annotations
        background_genes: Background gene set (default: all genes in annotations)
        categories: GO categories to include ('C', 'F', 'P')
        min_genes_per_term: Minimum genes per GO term to test
        max_genes_per_term: Maximum genes per GO term to test

    Returns:
        DataFrame with enrichment results:
            - go_id, go_term, category
            - query_count: genes in query with this term
            - background_count: genes in background with this term
            - fold_enrichment: (query_count/query_total) / (bg_count/bg_total)
            - pvalue, padj: Fisher's exact test p-value and FDR-corrected
    """
    if not SCIPY_AVAILABLE:
        raise ImportError(
            "scipy is required for GO enrichment. "
            "Install with: pip install scipy"
        )

    # Determine background
    if background_genes is None:
        background_genes = go_annotations['gene_id'].unique().tolist()

    # Filter to genes in background
    query_genes = [g for g in query_genes if g in background_genes]

    if len(query_genes) == 0:
        return pd.DataFrame()

    # Filter categories
    if categories:
        go_annotations = go_annotations[go_annotations['category'].isin(categories)]

    # Get unique GO terms
    go_terms = go_annotations.groupby(['go_id', 'go_term', 'category'])['gene_id'].apply(set).reset_index()
    go_terms.columns = ['go_id', 'go_term', 'category', 'genes']

    # Filter by term size
    go_terms['n_genes'] = go_terms['genes'].apply(len)
    go_terms = go_terms[
        (go_terms['n_genes'] >= min_genes_per_term) &
        (go_terms['n_genes'] <= max_genes_per_term)
    ]

    if go_terms.empty:
        return pd.DataFrame()

    # Run Fisher's exact test for each term
    query_set = set(query_genes)
    background_set = set(background_genes)
    n_query = len(query_set)
    n_background = len(background_set)

    results = []
    for _, row in go_terms.iterrows():
        term_genes = row['genes'] & background_set  # Only count genes in background

        # Contingency table:
        # | In term | Not in term |
        # Query:      a       b
        # Background: c       d

        a = len(query_set & term_genes)  # Query genes in term
        b = n_query - a                   # Query genes not in term
        c = len(term_genes) - a           # Background genes in term (not in query)
        d = n_background - n_query - c    # Background genes not in term

        # Skip if no overlap
        if a == 0:
            continue

        # Fisher's exact test (one-sided: enrichment)
        _, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative='greater')

        # Fold enrichment
        expected = n_query * len(term_genes) / n_background
        fold_enrichment = a / expected if expected > 0 else 0

        results.append({
            'go_id': row['go_id'],
            'go_term': row['go_term'],
            'category': row['category'],
            'query_count': a,
            'background_count': len(term_genes),
            'query_total': n_query,
            'background_total': n_background,
            'fold_enrichment': fold_enrichment,
            'pvalue': pvalue,
        })

    if not results:
        return pd.DataFrame()

    results_df = pd.DataFrame(results)

    # Multiple testing correction (Benjamini-Hochberg)
    results_df = results_df.sort_values('pvalue')
    n_tests = len(results_df)
    results_df['padj'] = _benjamini_hochberg(results_df['pvalue'].values, n_tests)

    # Sort by padj
    results_df = results_df.sort_values('padj')

    return results_df


def _benjamini_hochberg(pvalues: np.ndarray, n_tests: int) -> np.ndarray:
    """
    Apply Benjamini-Hochberg FDR correction with monotonicity enforcement.
    """
    n = len(pvalues)
    if n == 0:
        return np.array([])

    # Calculate adjusted p-values
    ranks = np.arange(1, n + 1)
    padj = pvalues * n_tests / ranks

    # Enforce monotonicity (running minimum from right to left)
    padj_monotonic = np.minimum.accumulate(padj[::-1])[::-1]

    # Cap at 1.0
    padj_monotonic = np.minimum(padj_monotonic, 1.0)

    return padj_monotonic


def plot_go_enrichment(
    enrichment_df: pd.DataFrame,
    n_terms: int = 20,
    padj_threshold: float = 0.05,
    figsize: Tuple[int, int] = (10, 8),
    output_path: Optional[str] = None,
    title: str = 'GO Enrichment Analysis',
) -> Optional[plt.Figure]:
    """
    Create bar plot of top enriched GO terms.

    Args:
        enrichment_df: Results from run_go_enrichment()
        n_terms: Number of top terms to show
        padj_threshold: Threshold for significance
        figsize: Figure size
        output_path: Path to save figure
        title: Plot title

    Returns:
        matplotlib Figure
    """
    if not MATPLOTLIB_AVAILABLE:
        return None

    if enrichment_df.empty:
        print("Warning: No enrichment results to plot")
        return None

    # Filter significant and take top N
    sig_df = enrichment_df[enrichment_df['padj'] < padj_threshold].copy()

    if sig_df.empty:
        print(f"Warning: No terms significant at padj < {padj_threshold}")
        # Show top terms anyway
        sig_df = enrichment_df.head(n_terms).copy()

    plot_df = sig_df.head(n_terms).copy()

    if plot_df.empty:
        return None

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Color by category
    category_colors = {
        'P': '#2ca02c',  # Green - Biological Process
        'F': '#1f77b4',  # Blue - Molecular Function
        'C': '#ff7f0e',  # Orange - Cellular Component
    }

    colors = [category_colors.get(c, 'gray') for c in plot_df['category']]

    # Plot horizontal bars
    y_pos = np.arange(len(plot_df))
    bars = ax.barh(y_pos, plot_df['fold_enrichment'], color=colors, alpha=0.8)

    # Labels
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['go_term'])
    ax.invert_yaxis()  # Top term at top

    ax.set_xlabel('Fold Enrichment')
    ax.set_title(title)

    # Add gene counts
    for i, (idx, row) in enumerate(plot_df.iterrows()):
        ax.text(
            row['fold_enrichment'] + 0.1,
            i,
            f"({row['query_count']}/{row['background_count']})",
            va='center',
            fontsize=8,
        )

    # Legend for categories
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=category_colors['P'], label='Biological Process'),
        Patch(facecolor=category_colors['F'], label='Molecular Function'),
        Patch(facecolor=category_colors['C'], label='Cellular Component'),
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=200, bbox_inches='tight')
        print(f"Saved GO enrichment plot to {output_path}")

    return fig


def run_enrichment_comparison(
    gene_lists: Dict[str, List[str]],
    go_annotations: pd.DataFrame,
    background_genes: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Run enrichment for multiple gene lists and compare.

    Args:
        gene_lists: Dict mapping list_name -> gene list
        go_annotations: GO annotation DataFrame
        background_genes: Background gene set

    Returns:
        DataFrame with combined enrichment results
    """
    all_results = []

    for list_name, genes in gene_lists.items():
        results = run_go_enrichment(
            genes,
            go_annotations,
            background_genes=background_genes,
        )

        if not results.empty:
            results['gene_list'] = list_name
            all_results.append(results)

    if not all_results:
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)
    return combined
