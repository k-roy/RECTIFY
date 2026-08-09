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
    _aggregate_counts_by_gene,
)
from rectify.core.analyze.pca import (
    run_pca_analysis,
)
from rectify.core.analyze.bedgraph import generate_bedgraphs
from rectify.core.analyze.loaders import load_annotation, load_corrected_positions
from rectify.core.analyze.genomic_distribution import classify_positions_by_region
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
from rectify.core.analyze.cluster_gene_attribution import (
    build_reference_cluster_lookup,
    reference_position_attributions_to_clusters,
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


class TestClusterCenterOfMass:
    """Tests for the read-weighted cluster_com column (dev/TODO.md)."""

    def _skewed_df(self):
        # Two positions in one cluster, weighted 9:1 toward 100.
        # Read-weighted COM = floor((100*9 + 110*1)/10) = 101.
        # Unweighted median (modal_position) = 105 — deliberately distinct.
        return pd.DataFrame({
            'chrom': ['chr1', 'chr1'],
            'strand': ['+', '+'],
            'corrected_position': [100, 110],
            '_count': [9, 1],
        })

    def test_fixed_distance_com_is_read_weighted(self):
        clusters = cluster_cpa_sites(
            self._skewed_df(), cluster_distance=25, min_reads=1,
            count_col='_count',
        )
        assert len(clusters) == 1
        row = clusters.iloc[0]
        assert 'cluster_com' in clusters.columns
        assert int(row['cluster_com']) == 101
        # COM must differ from the unweighted modal/median position.
        assert int(row['modal_position']) == 105
        assert int(row['n_reads']) == 10

    def test_adaptive_com_is_read_weighted(self):
        # Tight 5-position cluster with the count mass piled on the right edge
        # (104, count 16) but a left tail. Weighted COM = floor(2070/20) = 103,
        # distinct from the peak (modal_position) at 104.
        df = pd.DataFrame({
            'chrom': ['chr1'] * 5,
            'strand': ['+'] * 5,
            'corrected_position': [100, 101, 102, 103, 104],
            '_count': [1, 1, 1, 1, 16],
        })
        clusters = cluster_cpa_sites_adaptive(
            df, max_cluster_radius=25,
            min_peak_separation=10, min_reads=1, count_col='_count',
        )
        assert len(clusters) == 1
        row = clusters.iloc[0]
        assert 'cluster_com' in clusters.columns
        assert int(row['cluster_com']) == 103
        # Adaptive modal_position is the peak (highest-count position) = 104.
        assert int(row['modal_position']) == 104

    def test_empty_input_has_com_column(self):
        empty = pd.DataFrame(columns=['chrom', 'strand', 'corrected_position'])
        assert 'cluster_com' in cluster_cpa_sites(empty, min_reads=1).columns
        assert 'cluster_com' in cluster_cpa_sites_adaptive(empty, min_reads=1).columns


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


class TestWeightedClusterGeneAttribution:
    """Tests for weighted cluster-to-gene attribution in analyze outputs."""

    def test_gene_counts_use_weighted_shared_clusters(self):
        counts = pd.DataFrame(
            {
                'wt_rep1': [100.0, 40.0],
                'wt_rep2': [120.0, 60.0],
                'mut_rep1': [80.0, 30.0],
            },
            index=['cluster_a', 'cluster_b'],
        )
        clusters = pd.DataFrame({
            'cluster_id': ['cluster_a', 'cluster_b'],
            'gene_id': ['fallback_a', 'fallback_b'],
        })
        attr = pd.DataFrame({
            'cluster_id': ['cluster_a', 'cluster_a', 'cluster_b'],
            'gene_id': ['GENE1', 'GENE2', 'GENE2'],
            'gene_name': ['GENE1', 'GENE2', 'GENE2'],
            'raw_attributed_count': [75.0, 25.0, 10.0],
            'attribution_weight': [0.75, 0.25, 1.0],
            'source': ['reference', 'reference', 'reference'],
        })

        gene_counts = _aggregate_counts_by_gene(counts, clusters, attr)

        assert gene_counts.loc['GENE1', 'wt_rep1'] == pytest.approx(75.0)
        assert gene_counts.loc['GENE2', 'wt_rep1'] == pytest.approx(25.0 + 40.0)
        assert gene_counts.loc['GENE2', 'wt_rep2'] == pytest.approx(30.0 + 60.0)

    def test_shift_analysis_uses_weighted_cluster_gene_map(self, sample_metadata):
        counts = pd.DataFrame(
            {
                'wt_rep1': [90.0, 10.0],
                'wt_rep2': [90.0, 10.0],
                'treat_rep1': [10.0, 90.0],
                'treat_rep2': [10.0, 90.0],
            },
            index=['cluster_prox', 'cluster_dist'],
        )
        clusters = pd.DataFrame({
            'cluster_id': ['cluster_prox', 'cluster_dist'],
            'chrom': ['chrI', 'chrI'],
            'strand': ['+', '+'],
            'start': [100, 300],
            'end': [100, 300],
            'modal_position': [100, 300],
            'gene_id': ['nearest_wrong_1', 'nearest_wrong_2'],
            'gene_name': ['nearest_wrong_1', 'nearest_wrong_2'],
        })
        attr = pd.DataFrame({
            'cluster_id': ['cluster_prox', 'cluster_dist'],
            'gene_id': ['GENE_BODY', 'GENE_BODY'],
            'gene_name': ['GENE_BODY', 'GENE_BODY'],
            'raw_attributed_count': [50.0, 50.0],
            'attribution_weight': [1.0, 1.0],
            'source': ['reference', 'reference'],
        })

        result = analyze_cluster_shifts(
            counts,
            clusters,
            'wt',
            'treat',
            sample_metadata,
            cluster_gene_attributions=attr,
            min_total_counts=1,
        )

        assert result['gene_id'].tolist() == ['GENE_BODY']
        assert result.loc[0, 'major_cluster_a'] == 'cluster_prox'
        assert result.loc[0, 'major_cluster_b'] == 'cluster_dist'
        assert result.loc[0, 'shift_bp'] == 200
        assert result.loc[0, 'distribution_divergence'] > 0.2

    def test_reference_position_attributions_map_to_current_clusters(self):
        pos_attr = pd.DataFrame({
            'chrom': ['chrI', 'chrI', 'chrI'],
            'strand': ['+', '+', '+'],
            'position': [101, 101, 305],
            'gene_id': ['GENE1', 'GENE2', 'GENE2'],
            'attributed_count': [3.0, 1.0, 6.0],
        })

        def lookup(chrom, strand, pos):
            if chrom == 'chrI' and strand == '+' and 95 <= pos <= 105:
                return 'cluster_a'
            if chrom == 'chrI' and strand == '+' and 300 <= pos <= 310:
                return 'cluster_b'
            return None

        mapped = reference_position_attributions_to_clusters(pos_attr, lookup, source='drs_origin')

        row_a1 = mapped[(mapped['cluster_id'] == 'cluster_a') & (mapped['gene_id'] == 'GENE1')].iloc[0]
        row_a2 = mapped[(mapped['cluster_id'] == 'cluster_a') & (mapped['gene_id'] == 'GENE2')].iloc[0]
        assert row_a1['attribution_weight'] == pytest.approx(0.75)
        assert row_a2['attribution_weight'] == pytest.approx(0.25)
        assert mapped[mapped['cluster_id'] == 'cluster_b']['gene_id'].tolist() == ['GENE2']

    def test_reference_position_attributions_use_bounded_nearest_cluster_window(self):
        clusters = pd.DataFrame({
            'cluster_id': ['cluster_a', 'cluster_b'],
            'chrom': ['chrI', 'chrI'],
            'strand': ['+', '+'],
            'start': [100, 150],
            'end': [105, 155],
            'modal_position': [102, 152],
        })
        pos_attr = pd.DataFrame({
            'chrom': ['chrI', 'chrI', 'chrI'],
            'strand': ['+', '+', '+'],
            'position': [110, 145, 130],
            'gene_id': ['GENE_A', 'GENE_B', 'GENE_TOO_FAR'],
            'attributed_count': [2.0, 3.0, 10.0],
        })

        exact_lookup = build_reference_cluster_lookup(clusters, max_distance=0)
        exact = reference_position_attributions_to_clusters(pos_attr, exact_lookup)
        assert exact.empty

        window_lookup = build_reference_cluster_lookup(clusters, max_distance=5)
        mapped = reference_position_attributions_to_clusters(pos_attr, window_lookup)

        assert mapped[['cluster_id', 'gene_id']].values.tolist() == [
            ['cluster_a', 'GENE_A'],
            ['cluster_b', 'GENE_B'],
        ]


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


class TestBedgraphCoordinates:
    """Regression: per-condition bedgraph must treat corrected_3prime as 0-based.

    corrected_3prime is stored throughout rectify as 0-based-inclusive (derived
    from pysam reference_end - 1 for is_reverse=True and reference_start for
    is_reverse=False). BED is 0-based half-open, so a single base at 0-based
    pos must emit `start=pos, end=pos+1`. The previous emitter wrote
    `start=pos-1, end=pos`, which shifted every per-condition bedgraph row
    1 bp to the left of the true position — see AGENT_FIXES.md entry
    [2026-05-20] for the CST6 chrIX:287,749 vs 287,750 case.
    """

    def test_single_base_position_emitted_zero_based(self, tmp_path):
        # Construct one read per condition; corrected_3prime is a 0-based
        # coordinate of the last aligned base (matches the convention from
        # protocols/quantseq_rev.py:72 and walkback.py:142).
        positions_df = pd.DataFrame({
            'sample': ['ysh1_rep1'],
            'chrom': ['chrIX'],
            'corrected_3prime': [287749],  # 0-based — IGV 1-based 287,750
            'strand': ['+'],
        })

        generate_bedgraphs(positions_df, tmp_path, normalize_rpm=False)

        out = tmp_path / 'ysh1_plus.bedgraph'
        assert out.exists()
        lines = [l for l in out.read_text().splitlines() if not l.startswith('track')]
        assert len(lines) == 1
        chrom, start, end, value = lines[0].split('\t')
        # A single base at 0-based 287,749 is the half-open interval [287749, 287750).
        assert (chrom, int(start), int(end)) == ('chrIX', 287749, 287750)
        assert float(value) == 1.0

    def test_minus_strand_routed_to_minus_file(self, tmp_path):
        positions_df = pd.DataFrame({
            'sample': ['wt_rep1', 'wt_rep1'],
            'chrom': ['chrI', 'chrI'],
            'corrected_3prime': [1000, 2000],
            'strand': ['+', '-'],
        })
        generate_bedgraphs(positions_df, tmp_path, normalize_rpm=False)

        plus_rows = [l for l in (tmp_path / 'wt_plus.bedgraph').read_text().splitlines()
                     if not l.startswith('track')]
        minus_rows = [l for l in (tmp_path / 'wt_minus.bedgraph').read_text().splitlines()
                      if not l.startswith('track')]
        assert plus_rows == ['chrI\t1000\t1001\t1.0000']
        assert minus_rows == ['chrI\t2000\t2001\t1.0000']

    def test_multiple_reads_at_same_position_aggregate(self, tmp_path):
        positions_df = pd.DataFrame({
            'sample': ['ysh1_rep1'] * 3,
            'chrom': ['chrIX'] * 3,
            'corrected_3prime': [287749, 287749, 287749],
            'strand': ['+', '+', '+'],
        })
        generate_bedgraphs(positions_df, tmp_path, normalize_rpm=False)

        rows = [l for l in (tmp_path / 'ysh1_plus.bedgraph').read_text().splitlines()
                if not l.startswith('track')]
        assert rows == ['chrIX\t287749\t287750\t3.0000']

    def test_fraction_weights_and_non_yeast_contigs_are_emitted(self, tmp_path):
        positions_df = pd.DataFrame({
            'sample': ['sample1', 'sample1'],
            'chrom': ['chrUn_KI270442v1', 'chrUn_KI270442v1'],
            'corrected_3prime': [10, 10],
            'strand': ['+', '+'],
            'fraction': [0.25, 0.75],
        })

        generate_bedgraphs(positions_df, tmp_path, normalize_rpm=False)

        rows = [l for l in (tmp_path / 'sample1_plus.bedgraph').read_text().splitlines()
                if not l.startswith('track')]
        assert rows == ['chrUn_KI270442v1\t10\t11\t1.0000']


class TestAnalyzePositionColumns:
    def test_conflicting_corrected_columns_fail_loudly(self, tmp_path):
        path = tmp_path / 'positions.tsv'
        pd.DataFrame({
            'sample': ['s1'],
            'chrom': ['chrI'],
            'strand': ['+'],
            'corrected_3prime': [100],
            'corrected_position': [101],
        }).to_csv(path, sep='\t', index=False)

        with pytest.raises(ValueError, match='corrected_3prime and corrected_position'):
            load_corrected_positions(str(path), sample_column='sample')


class TestAnnotationFeatureParsing:
    def test_gff_feature_rows_are_preserved_and_half_open(self, tmp_path):
        gff = tmp_path / 'features.gff3'
        gff.write_text(
            '\n'.join([
                '##gff-version 3',
                'chrI\tSGD\tgene\t1\t1000\t.\t+\t.\tID=geneA;Name=YAL001C;gene=TFC3',
                'chrI\tSGD\tfive_prime_UTR\t1\t100\t.\t+\t.\tID=utr5A;Parent=geneA',
                'chrI\tSGD\tCDS\t101\t800\t.\t+\t0\tID=cdsA;Parent=geneA',
                'chrI\tSGD\tthree_prime_UTR\t801\t1000\t.\t+\t.\tID=utr3A;Parent=geneA',
                'chrI\tSGD\tsnoRNA\t1201\t1230\t.\t-\t.\tID=sno1;Name=SNR1',
                'chrI\tSGD\tCUT\t1401\t1500\t.\t+\t.\tID=cut1;Name=CUT001',
                '',
            ])
        )

        ann = load_annotation(str(gff), normalize_chroms=False)

        assert set(['gene', 'UTR5', 'CDS', 'UTR3', 'snoRNA', 'CUT']).issubset(
            set(ann['feature_type'])
        )
        gene = ann[ann['feature_type'] == 'gene'].iloc[0]
        assert int(gene['start']) == 0
        assert int(gene['end']) == 1000
        assert gene['gene_id'] == 'geneA'
        assert gene['gene_name'] == 'TFC3'

        utr3 = ann[ann['feature_type'] == 'UTR3'].iloc[0]
        assert utr3['gene_id'] == 'geneA'
        assert utr3['gene_name'] == 'TFC3'
        assert int(utr3['start']) == 800
        assert int(utr3['end']) == 1000

    def test_feature_rows_drive_genomic_distribution_classification(self, tmp_path):
        gff = tmp_path / 'features.gff3'
        gff.write_text(
            '\n'.join([
                'chrI\tSGD\tgene\t1\t1000\t.\t+\t.\tID=geneA;Name=YAL001C;gene=TFC3',
                'chrI\tSGD\tfive_prime_UTR\t1\t100\t.\t+\t.\tID=utr5A;Parent=geneA',
                'chrI\tSGD\tCDS\t101\t800\t.\t+\t0\tID=cdsA;Parent=geneA',
                'chrI\tSGD\tthree_prime_UTR\t801\t1000\t.\t+\t.\tID=utr3A;Parent=geneA',
                'chrI\tSGD\tsnoRNA\t1201\t1230\t.\t-\t.\tID=sno1;Name=SNR1',
                'chrI\tSGD\tCUT\t1401\t1500\t.\t+\t.\tID=cut1;Name=CUT001',
                '',
            ])
        )
        ann = load_annotation(str(gff), normalize_chroms=False)
        positions = pd.DataFrame({
            'chrom': ['chrI', 'chrI', 'chrI', 'chrI'],
            'strand': ['+', '+', '-', '+'],
            'corrected_position': [850, 200, 1210, 1450],
            'sample': ['s1', 's1', 's1', 's1'],
        })

        classified = classify_positions_by_region(positions, ann)

        assert classified['genomic_region'].tolist() == [
            'UTR3',
            'UTR5_CDS',
            'snoRNA',
            'CUT',
        ]


class TestAdaptiveClusteringBoundaries:
    def test_valley_position_is_assigned_to_a_cluster(self):
        positions_df = pd.DataFrame({
            'chrom': ['chrI'] * 3,
            'strand': ['+'] * 3,
            'corrected_position': [100, 103, 106],
            'count': [10, 1, 10],
        })

        clusters = cluster_cpa_sites_adaptive(
            positions_df,
            max_cluster_radius=10,
            min_peak_separation=5,
            min_reads=1,
            count_col='count',
        )

        assert any(
            row.start <= 103 <= row.end
            for row in clusters.itertuples(index=False)
        )


# =============================================================================
# Transcript-model attribution wiring (planning/167 default-flip regression)
# =============================================================================
import os as _os_tm
from types import SimpleNamespace as _NS

_BUNDLED_GFF_TM = _os_tm.path.join(
    _os_tm.path.dirname(_os_tm.path.dirname(_os_tm.path.abspath(__file__))),
    "rectify", "data", "genomes", "saccharomyces_cerevisiae",
    "saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz",
)


def _tm_cluster(cid, chrom, strand, pos, n=50, r=5):
    return {"cluster_id": cid, "chrom": chrom, "strand": strand,
            "start": pos - r, "end": pos + r, "modal_position": pos, "n_reads": n}


@pytest.mark.skipif(not _os_tm.path.exists(_BUNDLED_GFF_TM), reason="bundled SGD GFF missing")
class TestContainmentAttributionWiring:
    """Exercises _build_cluster_gene_attribution_table's containment mode + the
    default-flip regression that canonical calls must not move (167 §5)."""

    def _dispatch(self, clusters_df, mode):
        from rectify.core.analyze.manifest import _build_cluster_gene_attribution_table
        from rectify.core.commands.analyze_command import _make_cluster_lookup
        ann = load_annotation(_BUNDLED_GFF_TM, normalize_chroms=False)
        # mirror the pipeline: legacy nearest-3' window runs first
        # (analyze_command:333), populating the legacy gene_id that the
        # containment proximity fallback reads.
        clusters_df = annotate_clusters_with_genes(clusters_df.copy(), ann)
        args = _NS(gene_attribution_mode=mode)
        lookup = _make_cluster_lookup(clusters_df)
        return _build_cluster_gene_attribution_table(
            args, samples=[], clusters_df=clusters_df, lookup=lookup,
            annotation_df=ann, chrom_format="passthrough", model=None,
        )

    def test_canonical_gene_id_stable_across_default_flip(self):
        """A canonical 3'UTR cluster (PGK1) attributes to the SAME gene under
        legacy-3prime-window AND containment — the default flip must not move it."""
        clusters = pd.DataFrame([_tm_cluster("d1", "chrIII", "+", 139100)])
        legacy_df, _ = self._dispatch(clusters, "legacy-3prime-window")
        cont_df, _ = self._dispatch(clusters, "containment")
        legacy_gene = legacy_df.set_index("cluster_id").loc["d1", "gene_id"]
        cont_gene = cont_df.set_index("cluster_id").loc["d1", "gene_id"]
        assert legacy_gene == "YCR012W"
        assert cont_gene == "YCR012W"

    def test_containment_rescues_deep_cds_cluster(self):
        """A deep-CDS cluster (RAD3) is intergenic under legacy, rescued under containment."""
        clusters = pd.DataFrame([_tm_cluster("c1", "chrV", "+", 528000)])
        legacy_df, _ = self._dispatch(clusters, "legacy-3prime-window")
        cont_df, _ = self._dispatch(clusters, "containment")
        legacy_gene = legacy_df.set_index("cluster_id").loc["c1", "gene_id"]
        assert legacy_gene is None or pd.isna(legacy_gene)
        assert cont_df.set_index("cluster_id").loc["c1", "gene_id"] == "YER171W"

    def test_readthrough_cluster_keeps_gene_id_under_containment(self):
        """A canonical cluster just past a gene 3' end (not body-contained) keeps
        its gene_id via the proximity-window fallback under the containment default."""
        ann = load_annotation(_BUNDLED_GFF_TM, normalize_chroms=False)
        from rectify.core.analyze.transcript_model import TranscriptModel
        m = TranscriptModel(ann)
        pos = m.genes["YCR012W"].canonical_cpa + 40   # +40 downstream, readthrough
        clusters = pd.DataFrame([_tm_cluster("rt", "chrIII", "+", pos)])
        cont_df, _ = self._dispatch(clusters, "containment")
        assert cont_df.set_index("cluster_id").loc["rt", "gene_id"] == "YCR012W"

    def test_full_column_set_added_by_classifier(self):
        """annotate_clusters_with_transcript_model adds region_class + coords."""
        from rectify.core.analyze.transcript_model import (
            TranscriptModel, annotate_clusters_with_transcript_model)
        ann = load_annotation(_BUNDLED_GFF_TM, normalize_chroms=False)
        m = TranscriptModel(ann)
        clusters = pd.DataFrame([
            _tm_cluster("i1", "chrXVI", "+", 173400),   # intronic
            _tm_cluster("u1", "chrIII", "+", 139100),   # 3'UTR
        ])
        out = annotate_clusters_with_transcript_model(clusters, m)
        assert "region_class" in out.columns and "distance_to_stop_codon" in out.columns
        by = out.set_index("cluster_id")
        assert by.loc["i1", "region_class"] == "intronic"
        assert by.loc["u1", "region_class"].startswith("3primeUTR")
