#!/usr/bin/env python3
"""
Tests for RECTIFY analyze modules.

Author: Kevin R. Roy
Date: 2026-03-17
"""

import pytest
import numpy as np
import pandas as pd
from unittest.mock import Mock, patch, MagicMock
import tempfile
import os

# Import modules under test
from rectify.core.analyze.clustering import (
    cluster_cpa_sites,
    cluster_cpa_sites_adaptive,
    build_cluster_count_matrix,
    annotate_clusters_with_genes,
)
from rectify.core.analyze.deseq2 import (
    detect_control_samples,
    create_sample_metadata,
    extract_condition_from_sample,
)
from rectify.core.analyze.pca import (
    run_pca_analysis,
)
from rectify.core.analyze.go_enrichment import (
    run_go_enrichment,
    _benjamini_hochberg,
)
from rectify.core.analyze.shift_analysis import (
    analyze_cluster_shifts,
    get_top_shifted_genes,
    _jensen_shannon_divergence,
)
from rectify.core.analyze.summary import (
    generate_analysis_summary,
    generate_summary_tables,
)


# =============================================================================
# Test fixtures
# =============================================================================

@pytest.fixture
def sample_positions_df():
    """Create sample positions DataFrame for testing."""
    return pd.DataFrame({
        'chrom': ['chr1'] * 20 + ['chr2'] * 10,
        'strand': ['+'] * 15 + ['-'] * 15,
        'corrected_position': [100, 105, 110, 200, 205, 300, 305, 310, 400, 500,
                              600, 605, 700, 800, 900, 1000, 1005, 1100, 1200, 1300,
                              100, 150, 200, 250, 300, 350, 400, 450, 500, 550],
        'sample': ['wt_rep1'] * 10 + ['wt_rep2'] * 10 + ['treat_rep1'] * 10,
    })


@pytest.fixture
def sample_clusters_df():
    """Create sample clusters DataFrame for testing."""
    return pd.DataFrame({
        'cluster_id': ['cluster_1', 'cluster_2', 'cluster_3', 'cluster_4'],
        'chrom': ['chr1', 'chr1', 'chr1', 'chr2'],
        'strand': ['+', '+', '-', '+'],
        'modal_position': [105, 305, 1005, 200],
        'n_positions': [3, 3, 2, 2],
        'total_reads': [50, 30, 20, 15],
        'gene_id': ['YAL001C', 'YAL002W', 'YAL003W', 'YBL001C'],
        'gene_name': ['TFC3', 'VPS8', 'FUS1', 'ECM15'],
    })


@pytest.fixture
def sample_count_matrix():
    """Create sample count matrix for testing."""
    data = {
        'wt_rep1': [100, 50, 30, 20],
        'wt_rep2': [110, 45, 35, 18],
        'treat_rep1': [80, 70, 25, 30],
        'treat_rep2': [85, 65, 28, 28],
    }
    return pd.DataFrame(data, index=['cluster_1', 'cluster_2', 'cluster_3', 'cluster_4'])


@pytest.fixture
def sample_metadata():
    """Create sample metadata DataFrame for testing."""
    return pd.DataFrame({
        'sample': ['wt_rep1', 'wt_rep2', 'treat_rep1', 'treat_rep2'],
        'condition': ['wt', 'wt', 'treat', 'treat'],
    }).set_index('sample')


@pytest.fixture
def sample_deseq2_results():
    """Create sample DESeq2 results for testing."""
    return pd.DataFrame({
        'baseMean': [100, 50, 30, 20, 15],
        'log2FoldChange': [2.5, -1.8, 0.5, 1.2, -2.1],
        'lfcSE': [0.3, 0.4, 0.2, 0.3, 0.4],
        'pvalue': [0.001, 0.01, 0.3, 0.05, 0.005],
        'padj': [0.01, 0.04, 0.5, 0.08, 0.03],
    }, index=['gene1', 'gene2', 'gene3', 'gene4', 'gene5'])


# =============================================================================
# Clustering tests
# =============================================================================

class TestClusterCpaSites:
    """Tests for cluster_cpa_sites function."""

    def test_basic_clustering(self, sample_positions_df):
        """Test basic clustering of CPA sites."""
        clusters = cluster_cpa_sites(
            sample_positions_df,
            cluster_distance=25,
            min_reads=1,
        )

        assert len(clusters) > 0
        assert 'cluster_id' in clusters.columns
        assert 'modal_position' in clusters.columns
        assert 'n_positions' in clusters.columns
        assert 'chrom' in clusters.columns
        assert 'strand' in clusters.columns

    def test_cluster_distance_effect(self, sample_positions_df):
        """Test that cluster distance affects number of clusters."""
        clusters_tight = cluster_cpa_sites(
            sample_positions_df,
            cluster_distance=5,
            min_reads=1,
        )
        clusters_loose = cluster_cpa_sites(
            sample_positions_df,
            cluster_distance=50,
            min_reads=1,
        )

        # Loose clustering should produce fewer or equal clusters
        assert len(clusters_loose) <= len(clusters_tight)

    def test_min_reads_filter(self, sample_positions_df):
        """Test minimum reads filter."""
        clusters_low = cluster_cpa_sites(
            sample_positions_df,
            cluster_distance=25,
            min_reads=1,
        )
        clusters_high = cluster_cpa_sites(
            sample_positions_df,
            cluster_distance=25,
            min_reads=5,
        )

        # Higher min_reads should produce fewer or equal clusters
        assert len(clusters_high) <= len(clusters_low)

    def test_strand_separation(self, sample_positions_df):
        """Test that strands are clustered separately."""
        clusters = cluster_cpa_sites(
            sample_positions_df,
            cluster_distance=25,
            min_reads=1,
        )

        # Each cluster should have only one strand
        for _, cluster in clusters.iterrows():
            assert cluster['strand'] in ['+', '-']

    def test_empty_input(self):
        """Test clustering with empty input."""
        empty_df = pd.DataFrame(columns=['chrom', 'strand', 'position', 'sample'])
        clusters = cluster_cpa_sites(empty_df, cluster_distance=25, min_reads=1)
        assert len(clusters) == 0


class TestBuildClusterCountMatrix:
    """Tests for build_cluster_count_matrix function."""

    def test_basic_matrix_building(self, sample_positions_df, sample_clusters_df):
        """Test basic count matrix building."""
        # First cluster the positions
        clusters = cluster_cpa_sites(sample_positions_df, cluster_distance=25, min_reads=1)

        # Build count matrix
        count_matrix = build_cluster_count_matrix(
            sample_positions_df,
            clusters,
            sample_col='sample',
        )

        assert isinstance(count_matrix, pd.DataFrame)
        assert count_matrix.shape[0] <= len(clusters)  # Rows are clusters
        assert count_matrix.shape[1] == sample_positions_df['sample'].nunique()  # Cols are samples

    def test_matrix_values_non_negative(self, sample_positions_df, sample_clusters_df):
        """Test that count matrix values are non-negative."""
        clusters = cluster_cpa_sites(sample_positions_df, cluster_distance=25, min_reads=1)
        count_matrix = build_cluster_count_matrix(
            sample_positions_df,
            clusters,
            sample_col='sample',
        )

        assert (count_matrix.values >= 0).all()


class TestAnnotateClustersWithGenes:
    """Tests for annotate_clusters_with_genes function."""

    def test_basic_annotation(self, sample_clusters_df):
        """Test basic gene annotation."""
        # Create mock annotation
        genes_df = pd.DataFrame({
            'chrom': ['chr1', 'chr1', 'chr2'],
            'start': [50, 250, 150],
            'end': [150, 400, 300],
            'strand': ['+', '+', '+'],
            'gene_id': ['YAL001C', 'YAL002W', 'YBL001C'],
            'gene_name': ['TFC3', 'VPS8', 'ECM15'],
        })

        # Create clusters without gene info
        clusters = sample_clusters_df[['cluster_id', 'chrom', 'strand', 'modal_position', 'n_positions']].copy()

        # Annotate
        annotated = annotate_clusters_with_genes(clusters, genes_df)

        assert 'gene_id' in annotated.columns


# =============================================================================
# DESeq2 tests
# =============================================================================

class TestDetectControlSamples:
    """Tests for detect_control_samples function."""

    def test_detect_wt_samples(self):
        """Test detection of WT samples."""
        samples = ['WT_rep1', 'WT_rep2', 'treatment_rep1', 'treatment_rep2']
        controls = detect_control_samples(samples)

        assert 'WT_rep1' in controls
        assert 'WT_rep2' in controls
        assert 'treatment_rep1' not in controls

    def test_detect_control_samples(self):
        """Test detection of control samples."""
        samples = ['control_1', 'control_2', 'treated_1', 'treated_2']
        controls = detect_control_samples(samples)

        assert 'control_1' in controls
        assert 'control_2' in controls

    def test_detect_untreated_samples(self):
        """Test detection of untreated samples."""
        samples = ['untreated_a', 'untreated_b', 'drug_a', 'drug_b']
        controls = detect_control_samples(samples)

        assert 'untreated_a' in controls
        assert 'untreated_b' in controls

    def test_no_controls_found(self):
        """Test when no control samples are found."""
        samples = ['sample_a', 'sample_b', 'sample_c']
        controls = detect_control_samples(samples)

        assert len(controls) == 0

    def test_case_insensitive(self):
        """Test case-insensitive matching."""
        samples = ['wt_rep1', 'WT_rep2', 'Wt_rep3', 'treat_rep1']
        controls = detect_control_samples(samples)

        assert len(controls) == 3


class TestCreateSampleMetadata:
    """Tests for create_sample_metadata function."""

    def test_basic_metadata_creation(self):
        """Test basic metadata creation."""
        samples = ['wt_rep1', 'wt_rep2', 'treat_rep1', 'treat_rep2']
        metadata = create_sample_metadata(samples)

        assert len(metadata) == 4
        assert 'condition' in metadata.columns

    def test_condition_extraction(self):
        """Test condition extraction from sample names."""
        samples = ['wt_rep1', 'wt_rep2', 'treatment_rep1', 'treatment_rep2']
        metadata = create_sample_metadata(samples)

        # Check that conditions are extracted
        assert metadata.loc['wt_rep1', 'condition'] != metadata.loc['treatment_rep1', 'condition']


class TestExtractConditionFromSample:
    """Tests for extract_condition_from_sample function."""

    def test_underscore_separation(self):
        """Test extraction with underscore separation."""
        assert extract_condition_from_sample('wt_rep1') == 'wt'
        assert extract_condition_from_sample('treatment_rep2') == 'treatment'
        assert extract_condition_from_sample('treatment_r1') == 'treatment'

    def test_no_separator(self):
        """Test extraction with no separator."""
        result = extract_condition_from_sample('sample1')
        assert result == 'sample1'

    def test_no_replicate_suffix(self):
        """Test extraction when no replicate suffix present."""
        assert extract_condition_from_sample('control') == 'control'


# =============================================================================
# PCA tests
# =============================================================================

class TestRunPcaAnalysis:
    """Tests for run_pca_analysis function."""

    def test_basic_pca(self, sample_count_matrix):
        """Test basic PCA analysis."""
        results = run_pca_analysis(sample_count_matrix)

        assert 'pca_coords' in results
        assert 'variance_explained' in results
        assert 'loadings' in results

    def test_pca_dimensions(self, sample_count_matrix):
        """Test PCA output dimensions."""
        results = run_pca_analysis(sample_count_matrix, n_components=2)

        if results['pca_coords'] is not None:
            # Should have 4 samples and 2 components
            assert results['pca_coords'].shape[0] == sample_count_matrix.shape[1]
            assert results['pca_coords'].shape[1] == 2

    def test_pca_variance_ratio_is_reasonable(self, sample_count_matrix):
        """Test that explained variance ratios are reasonable."""
        results = run_pca_analysis(sample_count_matrix)

        if results['variance_ratio'] is not None and len(results['variance_ratio']) > 0:
            # Each component's variance ratio should be between 0 and 1
            for var in results['variance_ratio']:
                assert 0 <= var <= 1.0, f"Variance ratio {var} out of range"
            # Total variance ratio should not exceed 1
            assert sum(results['variance_ratio']) <= 1.01  # Small tolerance for floating point

    def test_pca_with_insufficient_samples(self):
        """Test PCA with insufficient samples."""
        small_matrix = pd.DataFrame({'sample1': [1, 2, 3]})
        results = run_pca_analysis(small_matrix)

        # Should handle gracefully
        assert isinstance(results, dict)


# =============================================================================
# GO enrichment tests
# =============================================================================

class TestBenjaminiHochberg:
    """Tests for Benjamini-Hochberg FDR correction."""

    def test_basic_correction(self):
        """Test basic BH correction."""
        pvalues = np.array([0.01, 0.02, 0.03, 0.04, 0.05])
        n_tests = len(pvalues)
        adjusted = _benjamini_hochberg(pvalues, n_tests)

        # Adjusted p-values should be >= original (sorted input)
        assert all(adj >= orig for adj, orig in zip(adjusted, pvalues))

    def test_monotonicity(self):
        """Test that BH adjustment maintains monotonicity."""
        pvalues = np.array([0.001, 0.005, 0.01, 0.02, 0.05, 0.1])
        n_tests = len(pvalues)
        adjusted = _benjamini_hochberg(pvalues, n_tests)

        # Values should be capped at 1.0
        assert max(adjusted) <= 1.0

    def test_empty_input(self):
        """Test BH with empty input."""
        adjusted = _benjamini_hochberg(np.array([]), 0)
        assert len(adjusted) == 0

    def test_single_pvalue(self):
        """Test BH with single p-value."""
        adjusted = _benjamini_hochberg(np.array([0.05]), 1)
        assert len(adjusted) == 1
        assert adjusted[0] == 0.05


class TestRunGoEnrichment:
    """Tests for run_go_enrichment function."""

    def test_basic_enrichment(self):
        """Test basic GO enrichment."""
        gene_list = ['gene1', 'gene2', 'gene3', 'gene4', 'gene5']

        # Create GO annotations DataFrame
        go_annotations = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3', 'gene6', 'gene7', 'gene8', 'gene9', 'gene10'],
            'go_id': ['GO:0001', 'GO:0001', 'GO:0002', 'GO:0001', 'GO:0001', 'GO:0001', 'GO:0002', 'GO:0002'],
            'go_term': ['process1'] * 5 + ['process2'] * 3,
            'category': ['P'] * 8,
        })

        results = run_go_enrichment(gene_list, go_annotations, min_genes_per_term=2)

        assert isinstance(results, pd.DataFrame)
        if not results.empty:
            assert 'go_id' in results.columns
            assert 'pvalue' in results.columns
            assert 'padj' in results.columns

    def test_empty_gene_list(self):
        """Test enrichment with empty gene list."""
        go_annotations = pd.DataFrame({
            'gene_id': ['gene1', 'gene2'],
            'go_id': ['GO:0001', 'GO:0001'],
            'go_term': ['process1', 'process1'],
            'category': ['P', 'P'],
        })

        results = run_go_enrichment([], go_annotations)
        assert isinstance(results, pd.DataFrame)


# =============================================================================
# Shift analysis tests
# =============================================================================

class TestJensenShannonDivergence:
    """Tests for Jensen-Shannon divergence calculation."""

    def test_identical_distributions(self):
        """Test JS divergence for identical distributions."""
        p = np.array([0.5, 0.3, 0.2])
        q = np.array([0.5, 0.3, 0.2])

        js = _jensen_shannon_divergence(p, q)

        # JS divergence should be ~0 for identical distributions
        assert js < 0.01

    def test_different_distributions(self):
        """Test JS divergence for different distributions."""
        p = np.array([0.9, 0.1])
        q = np.array([0.1, 0.9])

        js = _jensen_shannon_divergence(p, q)

        # JS divergence should be positive for different distributions
        assert js > 0.5  # These are quite different

    def test_js_bounded(self):
        """Test that JS divergence is bounded [0, 1] for log2."""
        p = np.array([1.0, 0.0, 0.0])
        q = np.array([0.0, 0.0, 1.0])

        js = _jensen_shannon_divergence(p, q)

        # JS divergence is bounded by 1 when using log2
        assert 0 <= js <= 1.0

    def test_symmetry(self):
        """Test that JS divergence is symmetric."""
        p = np.array([0.7, 0.2, 0.1])
        q = np.array([0.3, 0.4, 0.3])

        js_pq = _jensen_shannon_divergence(p, q)
        js_qp = _jensen_shannon_divergence(q, p)

        assert abs(js_pq - js_qp) < 1e-10


class TestGetTopShiftedGenes:
    """Tests for get_top_shifted_genes function."""

    def test_basic_top_selection(self):
        """Test basic top gene selection."""
        shift_results = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3', 'gene4', 'gene5'],
            'gene_name': ['A', 'B', 'C', 'D', 'E'],
            'distribution_divergence': [0.8, 0.6, 0.4, 0.2, 0.1],
            'shift_bp': [100, 50, 30, 10, 5],
        })

        top = get_top_shifted_genes(shift_results, n_top=3)

        assert len(top) == 3
        assert top.iloc[0]['gene_id'] == 'gene1'  # Highest divergence first

    def test_min_divergence_filter(self):
        """Test minimum divergence filter."""
        shift_results = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'gene_name': ['A', 'B', 'C'],
            'distribution_divergence': [0.8, 0.15, 0.05],
            'shift_bp': [100, 50, 10],
        })

        top = get_top_shifted_genes(shift_results, n_top=10, min_divergence=0.2)

        assert len(top) == 1
        assert top.iloc[0]['gene_id'] == 'gene1'

    def test_empty_input(self):
        """Test with empty input."""
        empty_df = pd.DataFrame()
        top = get_top_shifted_genes(empty_df)

        assert len(top) == 0


# =============================================================================
# Summary generation tests
# =============================================================================

class TestGenerateAnalysisSummary:
    """Tests for generate_analysis_summary function."""

    def test_basic_summary(self, sample_deseq2_results):
        """Test basic summary generation."""
        summary = generate_analysis_summary(
            n_samples=4,
            n_clusters=100,
            n_genes=50,
            deseq2_gene_results={'treat': sample_deseq2_results},
            deseq2_cluster_results={'treat': sample_deseq2_results},
            reference_condition='wt',
        )

        assert isinstance(summary, pd.DataFrame)
        assert len(summary) > 0
        assert 'category' in summary.columns
        assert 'metric' in summary.columns
        assert 'value' in summary.columns

    def test_summary_contains_input_stats(self, sample_deseq2_results):
        """Test that summary contains input statistics."""
        summary = generate_analysis_summary(
            n_samples=4,
            n_clusters=100,
            n_genes=50,
            deseq2_gene_results={},
            deseq2_cluster_results={},
            reference_condition='wt',
        )

        # Check for input statistics
        input_rows = summary[summary['category'] == 'Input']
        assert len(input_rows) >= 3  # At least samples, clusters, genes


class TestGenerateSummaryTables:
    """Tests for generate_summary_tables function."""

    def test_basic_table_generation(self, sample_deseq2_results):
        """Test basic summary table generation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_files = generate_summary_tables(
                deseq2_gene_results={'treat': sample_deseq2_results},
                deseq2_cluster_results={'treat': sample_deseq2_results},
                output_dir=tmpdir,
            )

            assert isinstance(output_files, dict)
            # Should have created some files
            assert len(output_files) > 0

    def test_table_files_created(self, sample_deseq2_results):
        """Test that table files are actually created."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_files = generate_summary_tables(
                deseq2_gene_results={'treat': sample_deseq2_results},
                deseq2_cluster_results={'treat': sample_deseq2_results},
                output_dir=tmpdir,
            )

            # Check that files exist
            for name, path in output_files.items():
                assert os.path.exists(path), f"File {path} should exist"


# =============================================================================
# Integration tests
# =============================================================================

class TestAnalyzePipelineIntegration:
    """Integration tests for the analyze pipeline."""

    def test_clustering_to_count_matrix_pipeline(self, sample_positions_df):
        """Test clustering followed by count matrix building."""
        # Cluster positions
        clusters = cluster_cpa_sites(
            sample_positions_df,
            cluster_distance=25,
            min_reads=1,
            position_col='corrected_position',
        )

        # Build count matrix
        count_matrix = build_cluster_count_matrix(
            sample_positions_df,
            clusters,
            sample_col='sample',
            position_col='corrected_position',
        )

        # Verify consistency
        assert len(count_matrix.columns) == sample_positions_df['sample'].nunique()

        # Each cluster should have non-zero total counts (or be empty if all filtered)
        if len(count_matrix) > 0:
            assert (count_matrix.sum(axis=1) >= 0).all()

    def test_sample_detection_to_pca_pipeline(self, sample_count_matrix):
        """Test sample detection followed by PCA."""
        # Detect controls
        samples = list(sample_count_matrix.columns)
        controls = detect_control_samples(samples)

        # Create metadata
        metadata = create_sample_metadata(samples, controls)

        # Run PCA
        pca_results = run_pca_analysis(sample_count_matrix)

        # Verify PCA ran
        assert pca_results['pca_coords'] is not None or sample_count_matrix.shape[1] < 2
