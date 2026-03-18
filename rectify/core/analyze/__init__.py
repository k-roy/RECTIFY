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
]
