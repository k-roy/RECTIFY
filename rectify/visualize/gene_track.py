"""
Gene track panel drawing functions for RECTIFY visualization.

Provides functions for drawing gene structure tracks with box-arrow
representations of CDS features.

Features:
- Pentagon/arrow shapes indicating strand direction
- Vertical staggering to avoid overlaps
- Consistent coloring with per-gene expression panels
- Support for common and systematic gene names

Author: Kevin R. Roy
"""

from typing import Dict, List, Optional

try:
    from matplotlib.patches import Polygon
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

from .config import GENE_TYPE_COLORS


def assign_feature_levels(
    features: List[Dict],
    min_gap: int = 500,
) -> List[Dict]:
    """
    Assign y-levels to features to avoid overlaps.

    Uses greedy bin-packing algorithm: sorts features by start position
    and assigns each to the first level with sufficient gap.

    Args:
        features: List of feature dicts with 'start' and 'end' keys
        min_gap: Minimum gap (bp) between features to avoid overlap

    Returns:
        Same list of features with added 'level' key

    Example:
        features = [
            {'start': 100, 'end': 500, 'name': 'gene1'},
            {'start': 450, 'end': 900, 'name': 'gene2'},
            {'start': 1000, 'end': 1500, 'name': 'gene3'},
        ]
        features = assign_feature_levels(features)
        # gene1 and gene3 at level 0, gene2 at level 1
    """
    if not features:
        return features

    # Sort by start position
    sorted_features = sorted(features, key=lambda x: x['start'])

    # Track end position of last feature at each level
    level_ends: List[int] = []

    for feature in sorted_features:
        # Find first level where this feature fits
        placed = False
        for level_idx, level_end in enumerate(level_ends):
            if feature['start'] > level_end + min_gap:
                # Feature fits at this level
                feature['level'] = level_idx
                level_ends[level_idx] = feature['end']
                placed = True
                break

        if not placed:
            # Need a new level
            feature['level'] = len(level_ends)
            level_ends.append(feature['end'])

    return sorted_features


def draw_gene_track(
    ax,
    gene_name: str,
    gff_features: Dict,
    region_start: int,
    region_end: int,
    highlight_gene: bool = True,
    show_all_features: bool = True,
    show_labels: bool = True,
    use_common_names: bool = True,
    neighbor_genes: Optional[Dict[str, str]] = None,
    track_height: float = 0.8,
    arrow_tip_fraction: float = 0.15,
    max_arrow_tip: int = 400,
):
    """
    Draw gene structure track with all CDS annotations in region.

    Features are rendered as pentagonal arrow shapes indicating strand direction.
    Features are staggered vertically to avoid overlaps.

    Args:
        ax: Matplotlib axes to draw on
        gene_name: Target gene name to highlight
        gff_features: Dict of chromosome -> list of GeneFeature objects.
                      Must include '_gene_name_lookup' key for name mapping.
        region_start: Start of genomic region
        region_end: End of genomic region
        highlight_gene: If True, highlight the target gene in red
        show_all_features: If True, show all CDS in region (not just target)
        show_labels: If False, omit gene name labels (for crowded views)
        use_common_names: If True, use common gene names instead of systematic
        neighbor_genes: Optional dict mapping gene_type to gene_name, e.g.,
                       {'upstream_1': 'RSM10', 'downstream_1': 'ENA2'}
                       Used to color genes consistently with expression panels.
        track_height: Height of each gene arrow (fraction of level height)
        arrow_tip_fraction: Arrow tip length as fraction of gene width
        max_arrow_tip: Maximum arrow tip size in bp

    Example:
        from rectify.visualize import draw_gene_track

        # gff_features from rectify GFF parser
        draw_gene_track(
            ax,
            gene_name='ENA1',
            gff_features=gff_features,
            region_start=530000,
            region_end=535000,
            neighbor_genes={'upstream_1': 'RSM10', 'downstream_1': 'ENA2'}
        )
    """
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError("matplotlib is required for draw_gene_track")

    ax.set_xlim(region_start, region_end)

    # Get the gene name lookup table (systematic -> common)
    gene_name_lookup = gff_features.get('_gene_name_lookup', {})

    # Find the chromosome for this gene
    target_chrom = None
    for chrom, features in gff_features.items():
        if chrom.startswith('_'):  # Skip special keys
            continue
        for f in features:
            # Check both common_name and lookup table
            feature_common = getattr(f, 'common_name', None) or gene_name_lookup.get(
                getattr(f, 'gene_id', ''), ''
            )
            if feature_common == gene_name or getattr(f, 'gene_id', '') == gene_name:
                target_chrom = chrom
                break
        if target_chrom:
            break

    if target_chrom is None:
        ax.text(
            0.5, 0.5,
            f'{gene_name}\n(no GFF)',
            ha='center', va='center',
            transform=ax.transAxes,
            fontsize=8,
            color='gray',
        )
        ax.set_yticks([])
        ax.set_ylim(0, 1)
        return

    # Collect all CDS features in the region
    region_features = []
    for f in gff_features.get(target_chrom, []):
        feature_type = getattr(f, 'feature_type', None)
        if feature_type != 'CDS':
            continue

        f_start = getattr(f, 'start', 0)
        f_end = getattr(f, 'end', 0)

        # Check if feature overlaps with region (0-based half-open: [start, end))
        if f_end > region_start and f_start < region_end:
            # Clip to region
            feat_start = max(f_start, region_start)
            feat_end = min(f_end, region_end)

            # Determine names
            systematic_name = getattr(f, 'gene_id', 'unknown') or 'unknown'
            # Look up common name from parent gene feature
            common_name = gene_name_lookup.get(systematic_name, None)
            if not common_name and '_' in systematic_name:
                # Try without _CDS suffix
                parent_name = systematic_name.rsplit('_', 1)[0]
                common_name = gene_name_lookup.get(parent_name, None)

            # Extract just the systematic name
            if systematic_name and '_' in systematic_name:
                systematic_name = systematic_name.split('_')[0]

            display_name = common_name if common_name else systematic_name
            is_target = (
                display_name == gene_name or
                getattr(f, 'gene_id', '') == gene_name or
                common_name == gene_name
            )

            region_features.append({
                'start': feat_start,
                'end': feat_end,
                'strand': getattr(f, 'strand', '+'),
                'name': display_name,
                'common_name': common_name,
                'systematic_name': systematic_name,
                'is_target': is_target,
            })

    if not region_features:
        ax.text(
            0.5, 0.5,
            'No CDS in region',
            ha='center', va='center',
            transform=ax.transAxes,
            fontsize=8,
            color='gray',
        )
        ax.set_yticks([])
        ax.set_ylim(0, 1)
        return

    # Assign levels to avoid overlaps
    region_features = assign_feature_levels(region_features, min_gap=300)
    max_level = max(f['level'] for f in region_features)

    # Set y limits based on number of levels
    ax.set_ylim(-0.5, max_level + 1.5)

    # Build reverse lookup for neighbor genes (gene_name -> gene_type)
    gene_type_lookup = {}
    if neighbor_genes:
        for gene_type, g_name in neighbor_genes.items():
            if g_name:
                gene_type_lookup[g_name] = gene_type
                gene_type_lookup[g_name.upper()] = gene_type

    # Draw each feature as a pentagon/arrow
    for feature in region_features:
        start = feature['start']
        end = feature['end']
        width = end - start
        y_level = feature['level']

        # Determine gene type for coloring
        feature_common = feature.get('common_name', '')
        feature_systematic = feature.get('systematic_name', '')
        gene_type = None

        if feature['is_target']:
            gene_type = 'target'
        elif neighbor_genes:
            gene_type = (
                gene_type_lookup.get(feature_common) or
                gene_type_lookup.get(feature_systematic)
            )

        # Colors based on gene type
        if gene_type and gene_type in GENE_TYPE_COLORS:
            color = GENE_TYPE_COLORS[gene_type]
            alpha = 0.9 if gene_type == 'target' else 0.8
            linewidth = 1.5 if gene_type == 'target' else 1.0
        elif feature['is_target']:
            color = '#e74c3c'  # Red
            alpha = 0.95
            linewidth = 1.5
        else:
            color = '#95a5a6'  # Gray
            alpha = 0.5
            linewidth = 0.8

        # Arrow tip size
        arrow_tip = min(width * arrow_tip_fraction, max_arrow_tip)

        # Y coordinates for pentagon
        half_height = track_height / 2
        y_bot = y_level - half_height
        y_mid = y_level
        y_top = y_level + half_height

        # Create pentagonal arrow shape
        if feature['strand'] == '+':
            vertices = [
                (start, y_bot),
                (start, y_top),
                (end - arrow_tip, y_top),
                (end, y_mid),
                (end - arrow_tip, y_bot),
            ]
        else:
            vertices = [
                (start + arrow_tip, y_top),
                (end, y_top),
                (end, y_bot),
                (start + arrow_tip, y_bot),
                (start, y_mid),
            ]

        arrow_patch = Polygon(
            vertices,
            closed=True,
            facecolor=color,
            edgecolor='black',
            linewidth=linewidth,
            alpha=alpha,
        )
        ax.add_patch(arrow_patch)

        # Label (skip if show_labels is False)
        if show_labels:
            label_x = (start + end) / 2
            label_y = y_level

            # Determine label text
            if use_common_names:
                label_text = feature.get('common_name') or feature.get('systematic_name', feature['name'])
            else:
                label_text = feature.get('systematic_name', feature['name'])

            # Shorter labels for non-target genes
            if feature['is_target']:
                fontsize = 8
                fontweight = 'bold'
                fontcolor = 'white'
            else:
                label_text = label_text[:10] if len(label_text) > 10 else label_text
                fontsize = 6
                fontweight = 'normal'
                fontcolor = 'white'

            ax.text(
                label_x, label_y, label_text,
                ha='center', va='center',
                fontsize=fontsize,
                fontweight=fontweight,
                color=fontcolor,
                alpha=0.95,
            )

    ax.set_yticks([])
    ax.set_ylabel('Genes', fontsize=8)


def get_genes_in_region(
    gff_features: Dict,
    chrom: str,
    region_start: int,
    region_end: int,
) -> Dict[str, str]:
    """
    Get all genes (CDS features) that overlap a genomic region.

    Args:
        gff_features: Dict of chromosome -> list of GeneFeature objects
        chrom: Chromosome name
        region_start: Start of region
        region_end: End of region

    Returns:
        Dict mapping systematic name -> common name for genes in region

    Example:
        genes = get_genes_in_region(gff_features, 'chrIV', 530000, 535000)
        # {'YDR039C': 'ENA2', 'YDR040C': 'ENA1', 'YDR041W': 'RSM10'}
    """
    gene_name_lookup = gff_features.get('_gene_name_lookup', {})
    genes_in_region = {}

    for f in gff_features.get(chrom, []):
        feature_type = getattr(f, 'feature_type', None)
        if feature_type != 'CDS':
            continue

        f_start = getattr(f, 'start', 0)
        f_end = getattr(f, 'end', 0)

        if f_end >= region_start and f_start <= region_end:
            systematic_name = getattr(f, 'gene_id', 'unknown') or 'unknown'
            if '_' in systematic_name:
                systematic_name = systematic_name.split('_')[0]

            common_name = gene_name_lookup.get(systematic_name, None)
            if not common_name:
                gene_id = getattr(f, 'gene_id', '')
                if gene_id and '_' in gene_id:
                    parent_name = gene_id.rsplit('_', 1)[0]
                    common_name = gene_name_lookup.get(parent_name, None)

            genes_in_region[systematic_name] = common_name or systematic_name

    return genes_in_region


def draw_gene_arrow(
    ax,
    start: int,
    end: int,
    strand: str,
    y_level: float = 0,
    height: float = 0.8,
    color: str = '#95a5a6',
    alpha: float = 0.8,
    linewidth: float = 1.0,
    label: Optional[str] = None,
    label_fontsize: int = 7,
):
    """
    Draw a single gene arrow (standalone function).

    Useful for custom gene track compositions.

    Args:
        ax: Matplotlib axes
        start: Start position
        end: End position
        strand: '+' or '-'
        y_level: Y-coordinate for center of arrow
        height: Height of arrow
        color: Fill color
        alpha: Transparency
        linewidth: Edge line width
        label: Optional label text
        label_fontsize: Font size for label
    """
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError("matplotlib is required for draw_gene_arrow")

    width = end - start
    arrow_tip = min(width * 0.15, 400)
    half_height = height / 2

    y_bot = y_level - half_height
    y_mid = y_level
    y_top = y_level + half_height

    if strand == '+':
        vertices = [
            (start, y_bot),
            (start, y_top),
            (end - arrow_tip, y_top),
            (end, y_mid),
            (end - arrow_tip, y_bot),
        ]
    else:
        vertices = [
            (start + arrow_tip, y_top),
            (end, y_top),
            (end, y_bot),
            (start + arrow_tip, y_bot),
            (start, y_mid),
        ]

    patch = Polygon(
        vertices,
        closed=True,
        facecolor=color,
        edgecolor='black',
        linewidth=linewidth,
        alpha=alpha,
    )
    ax.add_patch(patch)

    if label:
        ax.text(
            (start + end) / 2,
            y_level,
            label,
            ha='center',
            va='center',
            fontsize=label_fontsize,
            color='white',
        )
