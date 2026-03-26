"""
RECTIFY Visualization Module

Provides plotting utilities for RNA 3' end analysis:
- Metagene signal aggregation and plotting
- Consistent loci filtering across multiple data types
- Gene track rendering with box-arrow shapes
- Coverage track visualization
- Multi-track browser-style figures

Install with: pip install rectify-rna[visualize]

Example usage:

    # Quick metagene analysis
    from rectify.visualize import MetagenePipeline, PositionIndex, MetageneConfig

    index = PositionIndex(ends_df, position_col='position')
    pipeline = MetagenePipeline(MetageneConfig(window_upstream=100))
    result = pipeline.compute_profile(loci_df, index)

    # Consistent loci filtering across data types
    # (Use one "priority" dataset to determine which loci to include)
    from rectify.visualize import MetagenePipeline, LociFilter

    pipeline = MetagenePipeline(MetageneConfig(window_upstream=50))

    # Compute filter from priority dataset (e.g., RECTIFY 3' ends)
    loci_filter = pipeline.compute_loci_filter(
        loci_df, rectify_index,
        method='trimmed_mean',  # Exclude top/bottom 10% by signal
        proportion=0.1,
        priority_name='RECTIFY_WT'
    )

    # Apply same filter to all datasets
    results = {
        'RECTIFY': pipeline.compute_profile_filtered(loci_df, rectify_index, loci_filter),
        'NET-seq': pipeline.compute_profile_filtered(loci_df, netseq_index, loci_filter),
        'PAR-CLIP': pipeline.compute_profile_filtered(loci_df, parclip_index, loci_filter),
    }

    # Multi-track browser
    from rectify.visualize import MultiTrackFigure

    fig = (MultiTrackFigure()
        .add_gene_track(gff_features, highlight_gene="ENA1")
        .add_coverage_track("NET-seq", coverage_data)
    )
    fig.save("browser.png", region_start=530000, region_end=535000)

Author: Kevin R. Roy
"""

try:
    import matplotlib
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    raise ImportError(
        "rectify.visualize requires matplotlib. "
        "Install with: pip install rectify-rna[visualize]"
    )

# Config exports
from .config import (
    # Color palettes
    CODON_VARIANT_COLORS,
    GENE_TYPE_COLORS,
    STRAIN_ORIGIN_COLORS,
    CONTROL_TYPE_COLORS,
    EFFECT_DIRECTION_COLORS,
    WONG_COLORS,
    COLUMN_COLORS,
    # Marker shapes
    SUBLIBRARY_MARKERS,
    STRAIN_ORIGIN_MARKERS,
    STRAIN_MARKERS,
    GENE_TYPE_MARKERS,
    SOURCE_MARKERS,
    # Model info
    MODEL_NAMES,
    MODEL_SCORE_COLUMNS,
    PILOT_GENES,
    PILOT_ORF_BY_GENE,
    # Figure configuration
    FigureConfig,
)

# Figure utilities
from .figure_utils import (
    set_publication_style,
    save_multi_format,
    format_genomic_axis,
    add_trt_markers,
    despine,
    add_significance_bracket,
    create_figure_grid,
    # Metagene plotting utilities
    trimmed_mean,
    apply_window_sum_capping,
    plot_ridge_profiles,
    add_metagene_annotations,
    plot_metagene_line,
)

# Metagene
from .metagene import (
    MetageneConfig,
    PositionIndex,
    MetagenePipeline,
    LociFilter,
    build_index_from_bam,
)

# Gene track
from .gene_track import (
    assign_feature_levels,
    draw_gene_track,
    get_genes_in_region,
    draw_gene_arrow,
)

# VEP panels - disabled for now, pending further integration
# from .vep_panels import (
#     extract_position_from_variants,
#     draw_evo2_panel,
#     draw_esm1v_panel,
#     draw_shorkie_pergene_panel,
#     draw_yorzoi_pergene_panel,
#     draw_vep_panel_styled,
# )

# Coverage
from .coverage import (
    extract_coverage_from_bam,
    extract_coverage_from_array,
    draw_coverage_track,
    draw_strand_coverage,
    compare_coverage_tracks,
)

# Multi-track
from .multi_track import (
    MultiTrackFigure,
    create_gene_browser,
)

__all__ = [
    # Flag
    'MATPLOTLIB_AVAILABLE',

    # Config - Color palettes
    'CODON_VARIANT_COLORS',
    'GENE_TYPE_COLORS',
    'STRAIN_ORIGIN_COLORS',
    'CONTROL_TYPE_COLORS',
    'EFFECT_DIRECTION_COLORS',
    'WONG_COLORS',
    'COLUMN_COLORS',

    # Config - Marker shapes
    'SUBLIBRARY_MARKERS',
    'STRAIN_ORIGIN_MARKERS',
    'STRAIN_MARKERS',
    'GENE_TYPE_MARKERS',
    'SOURCE_MARKERS',

    # Config - Model info
    'MODEL_NAMES',
    'MODEL_SCORE_COLUMNS',
    'PILOT_GENES',
    'PILOT_ORF_BY_GENE',

    # Config - Figure settings
    'FigureConfig',

    # Figure utilities
    'set_publication_style',
    'save_multi_format',
    'format_genomic_axis',
    'add_trt_markers',
    'despine',
    'add_significance_bracket',
    'create_figure_grid',
    # Metagene plotting utilities
    'trimmed_mean',
    'apply_window_sum_capping',
    'plot_ridge_profiles',
    'add_metagene_annotations',
    'plot_metagene_line',

    # Metagene
    'MetageneConfig',
    'PositionIndex',
    'MetagenePipeline',
    'LociFilter',
    'build_index_from_bam',

    # Gene track
    'assign_feature_levels',
    'draw_gene_track',
    'get_genes_in_region',
    'draw_gene_arrow',

    # VEP panels - disabled for now
    # 'extract_position_from_variants',
    # 'draw_evo2_panel',
    # 'draw_esm1v_panel',
    # 'draw_shorkie_pergene_panel',
    # 'draw_yorzoi_pergene_panel',
    # 'draw_vep_panel_styled',

    # Coverage
    'extract_coverage_from_bam',
    'extract_coverage_from_array',
    'draw_coverage_track',
    'draw_strand_coverage',
    'compare_coverage_tracks',

    # Multi-track
    'MultiTrackFigure',
    'create_gene_browser',
]
