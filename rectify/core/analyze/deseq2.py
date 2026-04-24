#!/usr/bin/env python3
"""
DESeq2 Differential Expression Analysis Module

Performs differential expression analysis at both gene and cluster level
using pyDESeq2.

Author: Kevin R. Roy
Date: 2026-03-17
"""

import logging
from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import re
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Try to import pyDESeq2 (optional dependency)
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    PYDESEQ2_AVAILABLE = True
except ImportError:
    PYDESEQ2_AVAILABLE = False

# Default control sample patterns (case-insensitive)
CONTROL_PATTERNS = [
    r'^wt$', r'^wt_', r'_wt$', r'_wt_',
    r'^wildtype', r'wild_type', r'wild-type',
    r'^control', r'^ctrl$', r'^ctrl_', r'_ctrl$',
    r'^untreated', r'^mock$', r'^dmso$',
    r'^reference', r'^ref$',
]


def detect_control_samples(
    sample_names: List[str],
    patterns: Optional[List[str]] = None,
) -> List[str]:
    """
    Detect control/reference samples from sample names.

    Args:
        sample_names: List of sample names
        patterns: Optional custom regex patterns for control detection

    Returns:
        List of sample names identified as controls
    """
    if patterns is None:
        patterns = CONTROL_PATTERNS

    controls = []
    for sample in sample_names:
        sample_lower = sample.lower()
        for pattern in patterns:
            if re.search(pattern, sample_lower):
                controls.append(sample)
                break

    return controls


def extract_condition_from_sample(sample: str) -> str:
    """
    Extract condition name from sample name.

    Strips replicate suffixes like _rep1, _r1, _1, etc.
    """
    # Try various replicate patterns
    patterns = [
        r'_rep\d+$',  # _rep1, _rep2
        r'_r\d+$',    # _r1, _r2
        r'_\d+$',     # _1, _2
        r'\.rep\d+$', # .rep1, .rep2
        r'\.\d+$',    # .1, .2
    ]

    condition = sample
    for pattern in patterns:
        condition = re.sub(pattern, '', condition, flags=re.IGNORECASE)

    return condition


def create_sample_metadata(
    sample_names: List[str],
    control_samples: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Create sample metadata DataFrame for DESeq2.

    Args:
        sample_names: List of sample names
        control_samples: Optional list of control sample names

    Returns:
        DataFrame with 'sample' and 'condition' columns
    """
    metadata = pd.DataFrame({
        'sample': sample_names,
        'condition': [extract_condition_from_sample(s) for s in sample_names],
    })
    metadata.index = metadata['sample']

    # Detect controls if not provided
    if control_samples is None:
        control_samples = detect_control_samples(sample_names)

    # Add control flag
    metadata['is_control'] = metadata['sample'].apply(
        lambda x: any(extract_condition_from_sample(x) == extract_condition_from_sample(c)
                     for c in control_samples)
    )

    return metadata


def run_deseq2_gene_level(
    cluster_counts: pd.DataFrame,
    clusters_df: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    reference_condition: str,
    min_counts: int = 10,
    min_samples: int = 2,
    n_cpus: int = 4,
) -> Dict[str, pd.DataFrame]:
    """
    Run DESeq2 at gene level (aggregating clusters per gene).

    Args:
        cluster_counts: Cluster × sample count matrix
        clusters_df: Cluster definitions with gene_id column
        sample_metadata: Sample metadata with 'condition' column
        reference_condition: Reference condition for contrasts
        min_counts: Minimum counts per gene for inclusion
        min_samples: Gene must have min_counts in at least this many samples
        n_cpus: Number of CPUs for parallelization

    Returns:
        Dict mapping condition -> DESeq2 results DataFrame
    """
    if not PYDESEQ2_AVAILABLE:
        raise ImportError(
            "pydeseq2 is required for differential expression analysis. "
            "Install with: pip install pydeseq2"
        )

    # Aggregate counts by gene
    gene_counts = _aggregate_counts_by_gene(cluster_counts, clusters_df)

    if gene_counts.empty:
        return {}

    # Filter low-count genes
    gene_counts = _filter_low_counts(gene_counts, min_counts, min_samples)

    if gene_counts.empty:
        return {}

    # Ensure sample order matches metadata
    samples = [s for s in sample_metadata.index if s in gene_counts.columns]
    gene_counts = gene_counts[samples]
    metadata = sample_metadata.loc[samples].copy()

    # Run DESeq2
    results = _run_deseq2(
        gene_counts, metadata, reference_condition, n_cpus
    )

    return results


def run_deseq2_cluster_level(
    cluster_counts: pd.DataFrame,
    clusters_df: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    reference_condition: str,
    min_counts: int = 10,
    n_cpus: int = 4,
) -> Dict[str, pd.DataFrame]:
    """
    Run DESeq2 at cluster level.

    Args:
        cluster_counts: Cluster × sample count matrix
        clusters_df: Cluster definitions
        sample_metadata: Sample metadata with 'condition' column
        reference_condition: Reference condition for contrasts
        min_counts: Minimum total counts per cluster for inclusion
        n_cpus: Number of CPUs for parallelization

    Returns:
        Dict mapping condition -> DESeq2 results DataFrame (includes cluster info)
    """
    if not PYDESEQ2_AVAILABLE:
        raise ImportError(
            "pydeseq2 is required for differential expression analysis. "
            "Install with: pip install pydeseq2"
        )

    if cluster_counts.empty:
        return {}

    # Filter low-count clusters
    total_counts = cluster_counts.sum(axis=1)
    cluster_counts = cluster_counts[total_counts >= min_counts]

    if cluster_counts.empty:
        return {}

    # Ensure sample order matches metadata
    samples = [s for s in sample_metadata.index if s in cluster_counts.columns]
    cluster_counts = cluster_counts[samples]
    metadata = sample_metadata.loc[samples].copy()

    # Run DESeq2
    results = _run_deseq2(
        cluster_counts, metadata, reference_condition, n_cpus
    )

    # Add cluster annotations to results
    for condition, result_df in results.items():
        result_df = result_df.merge(
            clusters_df.set_index('cluster_id'),
            left_index=True,
            right_index=True,
            how='left'
        )
        n_missing = result_df['chrom'].isna().sum() if 'chrom' in result_df.columns else 0
        if n_missing > 0:
            import logging
            logging.warning(
                f"DESeq2 result for '{condition}': {n_missing} clusters have no annotation "
                f"(NaN values in merged output). Check that clusters_df is complete."
            )
        results[condition] = result_df

    return results


def _aggregate_counts_by_gene(
    cluster_counts: pd.DataFrame,
    clusters_df: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate cluster counts to gene level."""
    # Get gene_id for each cluster
    cluster_to_gene = clusters_df.set_index('cluster_id')['gene_id'].to_dict()

    # Filter clusters without gene annotation
    valid_clusters = [c for c in cluster_counts.index if c in cluster_to_gene
                      and cluster_to_gene[c] is not None]
    cluster_counts = cluster_counts.loc[valid_clusters]

    # Add gene_id column
    cluster_counts = cluster_counts.copy()
    cluster_counts['gene_id'] = [cluster_to_gene[c] for c in cluster_counts.index]

    # Aggregate by gene
    gene_counts = cluster_counts.groupby('gene_id').sum()

    return gene_counts


def _filter_low_counts(
    counts: pd.DataFrame,
    min_counts: int,
    min_samples: int,
) -> pd.DataFrame:
    """Filter features with low counts."""
    # Count samples with >= min_counts
    n_samples_above_threshold = (counts >= min_counts).sum(axis=1)

    # Keep features with counts in enough samples
    keep = n_samples_above_threshold >= min_samples

    return counts[keep]


def _run_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    reference_condition: str,
    n_cpus: int,
) -> Dict[str, pd.DataFrame]:
    """
    Run DESeq2 analysis.

    Returns results for each non-reference condition vs reference.
    """
    # Validate reference_condition is present in metadata before proceeding
    available_conditions = metadata['condition'].unique().tolist()
    if reference_condition not in available_conditions:
        raise ValueError(
            f"Reference condition '{reference_condition}' not found in sample metadata. "
            f"Available conditions: {available_conditions}"
        )

    # Convert counts to integers
    counts = counts.astype(int)

    # Create DESeq2 dataset
    # pyDESeq2 expects samples as rows, features as columns
    dds = DeseqDataSet(
        counts=counts.T,  # Transpose: samples × features
        metadata=metadata,
        design_factors="condition",
        n_cpus=n_cpus,
    )

    # Run DESeq2 normalization and dispersion estimation
    dds.deseq2()

    # pydeseq2 converts underscores to hyphens in factor levels (e.g. wt_by4742 → wt-by4742).
    # Build a mapping from original condition names to their transformed versions in dds.obs.
    orig_to_dds = {}
    for orig in metadata['condition'].unique():
        transformed = orig.replace('_', '-')
        if transformed in dds.obs['condition'].values:
            orig_to_dds[orig] = transformed
        else:
            orig_to_dds[orig] = orig  # no transformation occurred

    # Get all conditions
    conditions = metadata['condition'].unique()
    treatment_conditions = [c for c in conditions if c != reference_condition]

    results = {}
    for treatment in treatment_conditions:
        try:
            # Extract results for this contrast (use transformed names for dds.obs lookup)
            stat_res = DeseqStats(
                dds,
                contrast=["condition", orig_to_dds[treatment], orig_to_dds[reference_condition]],
            )
            stat_res.summary()

            # Get results DataFrame
            result_df = stat_res.results_df.copy()
            result_df.index.name = counts.index.name or 'feature_id'

            # Add significance flags
            result_df['significant_padj05'] = result_df['padj'] < 0.05
            result_df['significant_padj01'] = result_df['padj'] < 0.01
            result_df['direction'] = np.where(
                result_df['log2FoldChange'] > 0, 'up', 'down'
            )

            results[treatment] = result_df

        except Exception as e:
            logger.warning(
                "DESeq2 failed for %s vs %s: %s",
                treatment, reference_condition, e,
            )
            continue

    if not results:
        raise RuntimeError(
            f"DESeq2 failed for all contrasts against reference '{reference_condition}'. "
            "Check logs for per-contrast errors."
        )

    return results


def summarize_deseq2_results(
    results: Dict[str, pd.DataFrame],
    padj_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
) -> pd.DataFrame:
    """
    Summarize DESeq2 results across all contrasts.

    Args:
        results: Dict mapping condition -> DESeq2 results
        padj_threshold: Adjusted p-value threshold for significance
        lfc_threshold: Log2 fold-change threshold for biological significance

    Returns:
        Summary DataFrame with counts per condition
    """
    summary_rows = []

    for condition, df in results.items():
        sig = df['padj'] < padj_threshold
        sig_up = sig & (df['log2FoldChange'] > lfc_threshold)
        sig_down = sig & (df['log2FoldChange'] < -lfc_threshold)

        summary_rows.append({
            'condition': condition,
            'total_features': len(df),
            'significant': sig.sum(),
            'significant_up': sig_up.sum(),
            'significant_down': sig_down.sum(),
            'pct_significant': 100 * sig.sum() / len(df) if len(df) > 0 else 0,
        })

    return pd.DataFrame(summary_rows)
