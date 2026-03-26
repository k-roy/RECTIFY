"""
RECTIFY Analysis Module

Downstream analysis of corrected 3' end positions:
- CPA cluster formation
- Differential expression (DESeq2)
- PCA and sample clustering
- GO enrichment
- De novo motif discovery
- Summary reports
"""

from .clustering import (
    cluster_cpa_sites,
    cluster_cpa_sites_adaptive,
    build_cluster_count_matrix,
)

from .deseq2 import (
    run_deseq2_gene_level,
    run_deseq2_cluster_level,
    detect_control_samples,
    create_sample_metadata,
)

from .pca import (
    run_pca_analysis,
    plot_pca,
)

from .heatmap import (
    plot_sample_heatmap,
    plot_cluster_heatmap,
)

from .go_enrichment import (
    run_go_enrichment,
    plot_go_enrichment,
)

from .motif_discovery import (
    run_motif_discovery,
    extract_sequences_around_clusters,
    run_differential_motif_analysis,
    summarize_motif_results,
)

from .shift_analysis import (
    analyze_cluster_shifts,
    get_top_shifted_genes,
    plot_gene_browser,
    plot_shift_summary,
)

from .summary import (
    generate_summary_tables,
    generate_analysis_summary,
    generate_html_report,
)

from .deconvolution import (
    build_convolution_matrix,
    deconvolve_signal,
    deconvolve_with_confidence,
    identify_peaks,
    verify_boundary_recovery,
    DEFAULT_PSF,
)

from .genomic_distribution import (
    classify_positions_by_region,
    calculate_genomic_distribution,
    plot_genomic_distribution_pie,
    plot_genomic_distribution_pie_grid,
    plot_genomic_distribution_comparison,
    run_genomic_distribution_analysis,
)

from .atract_refiner import (
    ATractRefiner,
    refine_atract_position,
)

from .pan_mutant_refiner import (
    PanMutantRefiner,
)

from .gene_attribution import (
    get_read_5prime_position,
    get_read_3prime_position,
    get_read_body_interval,
    build_cds_interval_tree,
    find_overlapping_genes,
    compute_read_gene_attribution,
    aggregate_attributions_for_3prime_end,
    compute_cluster_attribution,
)

from .junction_validation import (
    JunctionEvidence,
    collect_junction_evidence,
    validate_novel_junctions,
    filter_records_to_validated_junctions,
    check_canonical_splice_motif,
    summarize_junction_validation,
    evidence_to_dataframe,
)

from .junction_analysis import (
    # NOTE: Splice site tolerance matching is disabled pending further testing.
    # The concern is that alternative 3' splice sites can be very close (3-10bp),
    # and tolerance-based snapping could incorrectly merge distinct biological sites.
    # match_junction_with_tolerance,
    # correct_junction_coordinates,
    # correct_all_junction_coordinates,
    # DEFAULT_SPLICE_SITE_TOLERANCE,
    deduplicate_novel_junctions,
    JunctionSignature,
    extract_junction_signatures,
    build_junction_index_from_records,
    build_known_junction_index,
    summarize_junction_analysis,
)

from .apa_detection import (
    APAIsoform,
    GeneAPAProfile,
    assign_reads_to_clusters,
    detect_apa_isoforms,
    build_gene_apa_profiles,
    quantify_apa_usage,
    identify_proximal_distal_tes,
    summarize_apa_detection,
    isoforms_to_dataframe,
)

__all__ = [
    # Clustering
    'cluster_cpa_sites',
    'cluster_cpa_sites_adaptive',
    'build_cluster_count_matrix',
    # DESeq2
    'run_deseq2_gene_level',
    'run_deseq2_cluster_level',
    'detect_control_samples',
    'create_sample_metadata',
    # PCA
    'run_pca_analysis',
    'plot_pca',
    # Heatmap
    'plot_sample_heatmap',
    'plot_cluster_heatmap',
    # GO enrichment
    'run_go_enrichment',
    'plot_go_enrichment',
    # Motif discovery
    'run_motif_discovery',
    'extract_sequences_around_clusters',
    'run_differential_motif_analysis',
    'summarize_motif_results',
    # Shift analysis
    'analyze_cluster_shifts',
    'get_top_shifted_genes',
    'plot_gene_browser',
    'plot_shift_summary',
    # Summary
    'generate_summary_tables',
    'generate_analysis_summary',
    'generate_html_report',
    # Deconvolution
    'build_convolution_matrix',
    'deconvolve_signal',
    'deconvolve_with_confidence',
    'identify_peaks',
    'verify_boundary_recovery',
    'DEFAULT_PSF',
    # Genomic distribution
    'classify_positions_by_region',
    'calculate_genomic_distribution',
    'plot_genomic_distribution_pie',
    'plot_genomic_distribution_pie_grid',
    'plot_genomic_distribution_comparison',
    'run_genomic_distribution_analysis',
    # A-tract refinement
    'ATractRefiner',
    'refine_atract_position',
    'PanMutantRefiner',
    # Gene attribution (5' end analysis)
    'get_read_5prime_position',
    'get_read_3prime_position',
    'get_read_body_interval',
    'build_cds_interval_tree',
    'find_overlapping_genes',
    'compute_read_gene_attribution',
    'aggregate_attributions_for_3prime_end',
    'compute_cluster_attribution',
    # Junction validation
    'JunctionEvidence',
    'collect_junction_evidence',
    'validate_novel_junctions',
    'filter_records_to_validated_junctions',
    'check_canonical_splice_motif',
    'summarize_junction_validation',
    'evidence_to_dataframe',
    # Junction analysis
    # NOTE: Splice site tolerance functions disabled pending testing
    # 'match_junction_with_tolerance',
    # 'correct_junction_coordinates',
    # 'correct_all_junction_coordinates',
    # 'DEFAULT_SPLICE_SITE_TOLERANCE',
    'deduplicate_novel_junctions',
    'JunctionSignature',
    'extract_junction_signatures',
    'build_junction_index_from_records',
    'build_known_junction_index',
    'summarize_junction_analysis',
    # APA detection
    'APAIsoform',
    'GeneAPAProfile',
    'assign_reads_to_clusters',
    'detect_apa_isoforms',
    'build_gene_apa_profiles',
    'quantify_apa_usage',
    'identify_proximal_distal_tes',
    'summarize_apa_detection',
    'isoforms_to_dataframe',
]
