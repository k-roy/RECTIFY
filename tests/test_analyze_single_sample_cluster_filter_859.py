#!/usr/bin/env python3
"""planning/859 regression: a single-sample `run-all --drs` must not write an EMPTY
cpa_clusters.tsv.

Measured on Hoffman2 with `rectify run-all wtaa_rep3_50k.bam --drs --Scer` (unit 859):

    [2/9] Forming CPA clusters (adaptive, radius=10bp, min_peak_sep=5bp, min_samples=2)...
      Formed 1,384 clusters
      Annotated 1,368 clusters (98.8%)
    [3/9] Building cluster count matrix...
      Kept 0/1,384 clusters present in >= 2 samples
      Matrix shape: 0 clusters x 1 samples
      Saved clusters to .../cpa_clusters.tsv          <-- header only

`--min-cluster-samples` is a REPRODUCIBILITY filter; with fewer samples than the
threshold no cluster can satisfy it, so it discarded 100% of them and `analyze`
still exited 0.  The fix clamps the threshold to the number of samples present.
"""
import pandas as pd
import pytest

from rectify.core.analyze.clustering import (
    cluster_cpa_sites_adaptive,
    build_cluster_count_matrix,
)


def _one_sample_positions():
    """One sample, three well-separated CPA sites with ample read support."""
    rows = []
    for pos in (1000, 5000, 9000):
        for i in range(12):
            rows.append({
                'chrom': 'chrI',
                'strand': '+',
                'corrected_3prime': pos + (i % 3),
                'sample': 'wtaa_rep3',
                'fraction': 1.0,
            })
    return pd.DataFrame(rows)


def _apply_presence_filter(clusters_df, count_matrix, min_cluster_samples):
    """The exact guard from analyze_command._run (planning/859 clamp included)."""
    if min_cluster_samples > 1 and not count_matrix.empty:
        n_present = count_matrix.shape[1]
        if n_present < min_cluster_samples:
            min_cluster_samples = n_present
    if min_cluster_samples > 1 and not count_matrix.empty:
        sample_counts = (count_matrix > 0).sum(axis=1)
        keep = set(sample_counts[sample_counts >= min_cluster_samples].index)
        clusters_df = clusters_df[clusters_df['cluster_id'].isin(keep)].copy()
    return clusters_df


def test_single_sample_keeps_its_clusters():
    positions = _one_sample_positions()
    clusters = cluster_cpa_sites_adaptive(
        positions, max_cluster_radius=10, min_peak_separation=5,
        min_reads=5, count_col='fraction',
    )
    assert len(clusters) > 0, "fixture produced no clusters at all"
    matrix = build_cluster_count_matrix(
        positions, clusters, sample_col='sample', count_col='fraction',
    )
    assert matrix.shape[1] == 1, "fixture must be a ONE-sample run"

    kept = _apply_presence_filter(clusters, matrix, min_cluster_samples=2)
    assert len(kept) == len(clusters), (
        "the >=2-samples reproducibility filter wiped out every cluster of a "
        "1-sample run — this is the planning/859 empty-cpa_clusters.tsv defect"
    )


def test_multi_sample_filter_still_bites():
    """The clamp must NOT weaken the filter when the samples are actually there."""
    positions = _one_sample_positions()
    extra = positions[positions['corrected_3prime'] < 2000].copy()
    extra['sample'] = 'wtaa_rep4'
    positions = pd.concat([positions, extra], ignore_index=True)

    clusters = cluster_cpa_sites_adaptive(
        positions, max_cluster_radius=10, min_peak_separation=5,
        min_reads=5, count_col='fraction',
    )
    matrix = build_cluster_count_matrix(
        positions, clusters, sample_col='sample', count_col='fraction',
    )
    assert matrix.shape[1] == 2, "fixture must be a TWO-sample run"

    kept = _apply_presence_filter(clusters, matrix, min_cluster_samples=2)
    assert len(kept) < len(clusters), (
        "with 2 real samples the presence filter must still drop the "
        "single-sample-only clusters"
    )


def test_genomic_distribution_module_has_a_logger():
    """planning/859: `logger` was referenced but never defined in this module, so the
    5'-end guard raised NameError and analyze_command downgraded the whole 5'-end
    distribution to a warning."""
    import logging
    from rectify.core.analyze import genomic_distribution as gd
    assert isinstance(getattr(gd, 'logger', None), logging.Logger)
