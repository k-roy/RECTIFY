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
    # n_tests must be the total number of GO terms tested (all size-filtered terms),
    # NOT just the subset with non-zero query overlap. Terms with a==0 were skipped
    # above but they were still hypotheses tested.
    n_tests = len(go_terms)
    results_df = results_df.sort_values('pvalue')
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


def plot_functional_enrichment_with_genes(
    up_genes: List[str],
    down_genes: List[str],
    go_annotations: pd.DataFrame,
    gene_name_mapping: Optional[Dict[str, str]] = None,
    functional_categories: Optional[Dict[str, List[str]]] = None,
    padj_threshold: float = 0.05,
    min_fold_enrichment: float = 1.5,
    max_genes_per_category: int = 15,
    figsize: Tuple[int, int] = (14, 8),
    output_path: Optional[str] = None,
    title: str = 'Functional Category Enrichment',
    condition_name: str = 'condition',
) -> Optional[plt.Figure]:
    """
    Create bidirectional bar plot showing up/down regulated genes by functional category.

    Displays gene names annotated on the bars, similar to the functional_category_enrichment
    style plot.

    Args:
        up_genes: List of upregulated gene IDs
        down_genes: List of downregulated gene IDs
        go_annotations: DataFrame with GO annotations
        gene_name_mapping: Dict mapping systematic names to common names (optional)
        functional_categories: Dict mapping category_name -> [GO terms].
            If None, uses default yeast functional categories.
        padj_threshold: p-value threshold for significance
        min_fold_enrichment: Minimum fold enrichment to include a category
        max_genes_per_category: Max genes to display in annotation
        figsize: Figure size
        output_path: Path to save figure
        title: Plot title
        condition_name: Name of the condition for labels

    Returns:
        matplotlib Figure
    """
    if not MATPLOTLIB_AVAILABLE:
        return None

    # Default functional categories for yeast (curated GO terms)
    if functional_categories is None:
        functional_categories = get_default_yeast_categories()

    # Run enrichment for each category
    all_genes_in_background = go_annotations['gene_id'].unique().tolist()

    category_data = []

    for cat_name, go_terms in functional_categories.items():
        # Get genes in this category
        cat_genes = set(
            go_annotations[go_annotations['go_term'].isin(go_terms)]['gene_id'].unique()
        )

        if not cat_genes:
            continue

        # Count up/down genes in category
        up_in_cat = set(up_genes) & cat_genes
        down_in_cat = set(down_genes) & cat_genes

        # Calculate enrichment (Fisher's exact test)
        if SCIPY_AVAILABLE and (up_in_cat or down_in_cat):
            # Test for upregulated
            up_pval = _fisher_test_category(
                up_in_cat, set(up_genes), cat_genes, set(all_genes_in_background)
            )

            # Test for downregulated
            down_pval = _fisher_test_category(
                down_in_cat, set(down_genes), cat_genes, set(all_genes_in_background)
            )
        else:
            up_pval = 1.0
            down_pval = 1.0

        category_data.append({
            'category': cat_name,
            'up_count': len(up_in_cat),
            'down_count': len(down_in_cat),
            'up_genes': list(up_in_cat),
            'down_genes': list(down_in_cat),
            'up_pval': up_pval,
            'down_pval': down_pval,
            'total_in_cat': len(cat_genes),
        })

    if not category_data:
        print("Warning: No functional categories found with genes")
        return None

    cat_df = pd.DataFrame(category_data)

    # Filter to significant categories (either up or down significant)
    cat_df = cat_df[
        ((cat_df['up_pval'] < padj_threshold) & (cat_df['up_count'] > 0)) |
        ((cat_df['down_pval'] < padj_threshold) & (cat_df['down_count'] > 0))
    ]

    if cat_df.empty:
        # Show top categories by count even if not significant
        cat_df = pd.DataFrame(category_data)
        cat_df['max_count'] = cat_df[['up_count', 'down_count']].max(axis=1)
        cat_df = cat_df.nlargest(10, 'max_count')

    # Sort by total effect size
    cat_df['effect'] = cat_df['up_count'] - cat_df['down_count']
    cat_df = cat_df.sort_values('effect', ascending=True)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    y_positions = np.arange(len(cat_df))
    bar_height = 0.7

    # Colors
    up_color = '#4393c3'    # Blue for upregulated
    down_color = '#d6604d'  # Red/orange for downregulated

    # Plot bars
    for i, (idx, row) in enumerate(cat_df.iterrows()):
        # Downregulated (negative direction)
        if row['down_count'] > 0:
            ax.barh(i, -row['down_count'], height=bar_height, color=down_color, alpha=0.8)

            # Add gene names
            gene_names = _format_gene_names(
                row['down_genes'], gene_name_mapping, max_genes_per_category
            )

            # Significance marker
            sig_marker = '*' if row['down_pval'] < padj_threshold else ''

            # Add count and genes
            ax.text(
                -row['down_count'] / 2, i,
                f"{row['down_count']}{sig_marker}\n{gene_names}",
                ha='center', va='center', fontsize=7, color='white', fontweight='bold'
            )

            # Add p-value
            if row['down_pval'] < padj_threshold:
                ax.text(
                    -row['down_count'] - 0.5, i - 0.25,
                    f"p={row['down_pval']:.1e}",
                    ha='right', va='center', fontsize=6, color='gray'
                )

        # Upregulated (positive direction)
        if row['up_count'] > 0:
            ax.barh(i, row['up_count'], height=bar_height, color=up_color, alpha=0.8)

            # Add gene names
            gene_names = _format_gene_names(
                row['up_genes'], gene_name_mapping, max_genes_per_category
            )

            # Significance marker
            sig_marker = '*' if row['up_pval'] < padj_threshold else ''

            # Add count and genes
            ax.text(
                row['up_count'] / 2, i,
                f"{row['up_count']}{sig_marker}\n{gene_names}",
                ha='center', va='center', fontsize=7, color='white', fontweight='bold'
            )

            # Add p-value
            if row['up_pval'] < padj_threshold:
                ax.text(
                    row['up_count'] + 0.5, i - 0.25,
                    f"p={row['up_pval']:.1e}",
                    ha='left', va='center', fontsize=6, color='gray'
                )

    # Add vertical line at 0
    ax.axvline(x=0, color='black', linewidth=1)

    # Labels
    ax.set_yticks(y_positions)
    ax.set_yticklabels(cat_df['category'], fontsize=10)
    ax.set_xlabel('Number of DE Genes', fontsize=11)
    ax.set_title(f"{title}\n(padj < {padj_threshold}; * = Fisher p < {padj_threshold})", fontsize=12)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=up_color, label=f'UP in {condition_name}'),
        Patch(facecolor=down_color, label=f'DOWN in {condition_name}'),
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    # Adjust x-axis to be symmetric
    max_val = max(cat_df['up_count'].max(), cat_df['down_count'].max())
    ax.set_xlim(-max_val * 1.3, max_val * 1.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=200, bbox_inches='tight')
        print(f"Saved functional enrichment plot to {output_path}")

    return fig


def _fisher_test_category(
    overlap: set,
    query_set: set,
    category_set: set,
    background_set: set,
) -> float:
    """Run Fisher's exact test for category enrichment."""
    if not SCIPY_AVAILABLE:
        return 1.0

    a = len(overlap)
    b = len(query_set) - a
    c = len(category_set & background_set) - a
    d = len(background_set) - len(query_set) - c

    if a == 0:
        return 1.0

    _, pval = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    return pval


def _format_gene_names(
    genes: List[str],
    name_mapping: Optional[Dict[str, str]],
    max_genes: int,
) -> str:
    """Format gene names for display, using common names if available."""
    if not genes:
        return ''

    # Get display names
    if name_mapping:
        display_names = [name_mapping.get(g, g) for g in genes]
    else:
        display_names = genes

    # Truncate if too many
    if len(display_names) > max_genes:
        display_names = display_names[:max_genes] + ['...']

    return ', '.join(display_names)


def get_default_yeast_categories() -> Dict[str, List[str]]:
    """
    Return default functional categories for S. cerevisiae.

    These are curated GO term descriptions that group related processes.
    """
    return {
        'Amino Acid Biosynthesis': [
            'amino acid biosynthetic process',
            'cellular amino acid biosynthetic process',
            'branched-chain amino acid biosynthetic process',
            'aromatic amino acid family biosynthetic process',
            'histidine biosynthetic process',
            'methionine biosynthetic process',
            'lysine biosynthetic process',
            'arginine biosynthetic process',
            'serine family amino acid biosynthetic process',
        ],
        'Purine/Pyrimidine': [
            'purine nucleotide biosynthetic process',
            'pyrimidine nucleotide biosynthetic process',
            'de novo IMP biosynthetic process',
            'purine nucleobase biosynthetic process',
            'nucleobase biosynthetic process',
        ],
        'Heat Shock/Stress': [
            'response to heat',
            'cellular response to heat',
            'protein folding',
            'protein refolding',
            'response to unfolded protein',
            'chaperone-mediated protein folding',
        ],
        'DNA Replication': [
            'DNA replication',
            'DNA replication initiation',
            'DNA synthesis involved in DNA repair',
            'regulation of DNA replication',
        ],
        'Histones/Chromatin': [
            'chromatin organization',
            'nucleosome assembly',
            'chromatin remodeling',
            'histone modification',
            'DNA packaging',
        ],
        'Ribosome Biogenesis': [
            'ribosome biogenesis',
            'rRNA processing',
            'rRNA metabolic process',
            'ribosomal large subunit biogenesis',
            'ribosomal small subunit biogenesis',
            'maturation of LSU-rRNA',
            'maturation of SSU-rRNA',
        ],
        'Translation': [
            'translation',
            'cytoplasmic translation',
            'translational initiation',
            'translational elongation',
        ],
        'mRNA Processing': [
            'mRNA processing',
            'mRNA splicing, via spliceosome',
            'RNA splicing',
            'mRNA polyadenylation',
            'mRNA 3\'-end processing',
        ],
        'Cell Cycle': [
            'cell cycle',
            'mitotic cell cycle',
            'cell division',
            'chromosome segregation',
            'mitotic spindle organization',
        ],
        'Proteasome': [
            'proteasome-mediated ubiquitin-dependent protein catabolic process',
            'proteasomal protein catabolic process',
            'ubiquitin-dependent protein catabolic process',
        ],
        'Glycolysis/TCA': [
            'glycolytic process',
            'tricarboxylic acid cycle',
            'pyruvate metabolic process',
            'gluconeogenesis',
        ],
        'Lipid Metabolism': [
            'lipid biosynthetic process',
            'fatty acid biosynthetic process',
            'sterol biosynthetic process',
            'phospholipid biosynthetic process',
        ],
        'Autophagy': [
            'autophagy',
            'macroautophagy',
            'autophagosome assembly',
        ],
        'Transcription': [
            'transcription by RNA polymerase II',
            'regulation of transcription by RNA polymerase II',
            'transcription initiation from RNA polymerase II promoter',
        ],
    }


def load_gene_name_mapping(
    gff_path: Optional[str] = None,
    mapping_file: Optional[str] = None,
) -> Dict[str, str]:
    """
    Load systematic to common gene name mapping.

    Args:
        gff_path: Path to GFF file (extracts Name and gene attributes)
        mapping_file: Path to TSV with systematic_name, common_name columns

    Returns:
        Dict mapping systematic name -> common name
    """
    mapping = {}

    if mapping_file and Path(mapping_file).exists():
        df = pd.read_csv(mapping_file, sep='\t')
        if 'systematic_name' in df.columns and 'common_name' in df.columns:
            for _, row in df.iterrows():
                if pd.notna(row['common_name']) and row['common_name']:
                    mapping[row['systematic_name']] = row['common_name']

    elif gff_path and Path(gff_path).exists():
        # Parse GFF for gene names
        import re
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.split('\t')
                if len(fields) < 9:
                    continue
                if fields[2] != 'gene':
                    continue

                attrs = fields[8]

                # Extract ID and Name
                id_match = re.search(r'ID=([^;]+)', attrs)
                name_match = re.search(r'Name=([^;]+)', attrs)
                gene_match = re.search(r'gene=([^;]+)', attrs)

                systematic = None
                common = None

                if id_match:
                    systematic = id_match.group(1)
                if name_match:
                    common = name_match.group(1)
                elif gene_match:
                    common = gene_match.group(1)

                if systematic and common and systematic != common:
                    mapping[systematic] = common

    return mapping
