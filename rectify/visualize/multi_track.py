"""
Multi-track browser-style figure composition for RECTIFY visualization.

Provides a fluent API for building multi-panel genomic browser figures
with gene tracks, coverage, VEP panels, and metagene profiles.

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


class MultiTrackFigure:
    """
    Multi-track browser-style figure composition.

    Provides a fluent API for building complex genomic browser figures
    with multiple panel types.

    Example:
        from rectify.visualize import MultiTrackFigure

        fig = (MultiTrackFigure(figsize=(12, 8))
            .add_gene_track(gff_features, highlight_gene="ENA1")
            .add_coverage_track("NET-seq", netseq_coverage, color="#0072B2")
            .add_coverage_track("3' ends", ends_coverage, color="#D55E00")
            .add_vep_track("evo2", df_evo2)
            .add_vep_track("esm1v", df_esm1v)
        )
        fig.save("browser.png", region_start=530000, region_end=535000)
    """

    def __init__(
        self,
        figsize: Tuple[float, float] = (12, 8),
        panel_height: float = 1.5,
        hspace: float = 0.3,
    ):
        """Initialize MultiTrackFigure.

        Args:
            figsize: Figure size (width, height) in inches
            panel_height: Default height for each panel in inches
            hspace: Vertical spacing between panels
        """
        if not MATPLOTLIB_AVAILABLE:
            raise ImportError("matplotlib is required for MultiTrackFigure")

        self.figsize = figsize
        self.panel_height = panel_height
        self.hspace = hspace
        self.tracks: List[Dict[str, Any]] = []

    def add_gene_track(
        self,
        gff_features: Dict,
        highlight_gene: Optional[str] = None,
        height_ratio: float = 0.7,
        **kwargs,
    ) -> 'MultiTrackFigure':
        """Add gene annotation track.

        Args:
            gff_features: GFF features dict from rectify parser
            highlight_gene: Gene to highlight
            height_ratio: Relative height of this panel
            **kwargs: Additional arguments passed to draw_gene_track

        Returns:
            self for method chaining
        """
        self.tracks.append({
            'type': 'gene',
            'data': gff_features,
            'highlight_gene': highlight_gene,
            'height_ratio': height_ratio,
            'kwargs': kwargs,
        })
        return self

    def add_coverage_track(
        self,
        label: str,
        coverage_data: Union[Tuple[np.ndarray, np.ndarray], Dict],
        color: str = '#0072B2',
        height_ratio: float = 1.0,
        **kwargs,
    ) -> 'MultiTrackFigure':
        """Add coverage plot track.

        Args:
            label: Track label
            coverage_data: Either (positions, depths) tuple or dict with
                           {'positions': array, 'depths': array}
            color: Line/fill color
            height_ratio: Relative height of this panel
            **kwargs: Additional arguments passed to draw_coverage_track

        Returns:
            self for method chaining
        """
        self.tracks.append({
            'type': 'coverage',
            'label': label,
            'data': coverage_data,
            'color': color,
            'height_ratio': height_ratio,
            'kwargs': kwargs,
        })
        return self

    def add_vep_track(
        self,
        panel_type: str,
        df: Any,
        height_ratio: float = 1.0,
        **kwargs,
    ) -> 'MultiTrackFigure':
        """Add VEP panel track.

        Args:
            panel_type: Type of VEP panel ('evo2', 'esm1v', 'shorkie', 'yorzoi')
            df: DataFrame with VEP scores
            height_ratio: Relative height of this panel
            **kwargs: Additional arguments passed to VEP panel function

        Returns:
            self for method chaining
        """
        self.tracks.append({
            'type': 'vep',
            'panel_type': panel_type,
            'data': df,
            'height_ratio': height_ratio,
            'kwargs': kwargs,
        })
        return self

    def add_metagene_track(
        self,
        label: str,
        profile_result: Dict,
        color: str = '#0072B2',
        height_ratio: float = 1.0,
        **kwargs,
    ) -> 'MultiTrackFigure':
        """Add metagene profile track.

        Args:
            label: Track label
            profile_result: Result dict from MetagenePipeline.compute_profile
            color: Line color
            height_ratio: Relative height of this panel
            **kwargs: Additional arguments

        Returns:
            self for method chaining
        """
        self.tracks.append({
            'type': 'metagene',
            'label': label,
            'data': profile_result,
            'color': color,
            'height_ratio': height_ratio,
            'kwargs': kwargs,
        })
        return self

    def add_custom_track(
        self,
        draw_func,
        data: Any = None,
        height_ratio: float = 1.0,
        **kwargs,
    ) -> 'MultiTrackFigure':
        """Add custom track with user-provided drawing function.

        Args:
            draw_func: Function that takes (ax, data, region_start, region_end, **kwargs)
            data: Data to pass to draw_func
            height_ratio: Relative height of this panel
            **kwargs: Additional arguments passed to draw_func

        Returns:
            self for method chaining
        """
        self.tracks.append({
            'type': 'custom',
            'draw_func': draw_func,
            'data': data,
            'height_ratio': height_ratio,
            'kwargs': kwargs,
        })
        return self

    def render(
        self,
        region_start: int,
        region_end: int,
        chrom: Optional[str] = None,
        gene_name: Optional[str] = None,
        title: Optional[str] = None,
        sharex: bool = True,
    ) -> plt.Figure:
        """Render figure with all tracks.

        Args:
            region_start: Start of genomic region
            region_end: End of genomic region
            chrom: Chromosome name (for coverage/VEP panels)
            gene_name: Target gene name (for gene track and per-gene panels)
            title: Optional figure title
            sharex: Whether to share x-axis across panels

        Returns:
            Matplotlib Figure object
        """
        from .gene_track import draw_gene_track
        from .coverage import draw_coverage_track
        from .vep_panels import (
            draw_evo2_panel,
            draw_esm1v_panel,
            draw_shorkie_pergene_panel,
            draw_yorzoi_pergene_panel,
        )
        from .metagene import MetagenePipeline
        from .figure_utils import format_genomic_axis

        n_tracks = len(self.tracks)
        if n_tracks == 0:
            raise ValueError("No tracks added. Use add_* methods first.")

        # Calculate height ratios
        height_ratios = [t['height_ratio'] for t in self.tracks]
        total_height = sum(height_ratios) * self.panel_height

        # Create figure
        fig = plt.figure(figsize=(self.figsize[0], total_height))
        gs = GridSpec(n_tracks, 1, height_ratios=height_ratios, hspace=self.hspace)

        axes = []
        for i, track in enumerate(self.tracks):
            if sharex and i > 0:
                ax = fig.add_subplot(gs[i], sharex=axes[0])
            else:
                ax = fig.add_subplot(gs[i])
            axes.append(ax)

            track_type = track['type']

            if track_type == 'gene':
                target_gene = track.get('highlight_gene') or gene_name
                draw_gene_track(
                    ax,
                    gene_name=target_gene or '',
                    gff_features=track['data'],
                    region_start=region_start,
                    region_end=region_end,
                    **track['kwargs'],
                )

            elif track_type == 'coverage':
                data = track['data']
                if isinstance(data, dict):
                    positions = data['positions']
                    depths = data['depths']
                else:
                    positions, depths = data

                draw_coverage_track(
                    ax,
                    positions,
                    depths,
                    region_start,
                    region_end,
                    color=track['color'],
                    label=track['label'],
                    **track['kwargs'],
                )
                ax.set_ylabel(track['label'], fontsize=8)

            elif track_type == 'vep':
                panel_type = track['panel_type'].lower()
                df = track['data']

                if panel_type == 'evo2':
                    draw_evo2_panel(ax, df, region_start, region_end, **track['kwargs'])
                elif panel_type == 'esm1v':
                    draw_esm1v_panel(ax, df, region_start, region_end, **track['kwargs'])
                elif panel_type == 'shorkie':
                    draw_shorkie_pergene_panel(
                        ax, df, region_start, region_end,
                        gene_name=gene_name or '',
                        **track['kwargs'],
                    )
                elif panel_type == 'yorzoi':
                    draw_yorzoi_pergene_panel(
                        ax, df, region_start, region_end,
                        gene_name=gene_name or '',
                        **track['kwargs'],
                    )
                else:
                    ax.text(0.5, 0.5, f'Unknown VEP type: {panel_type}',
                            ha='center', va='center', transform=ax.transAxes)

            elif track_type == 'metagene':
                result = track['data']
                pipeline = MetagenePipeline()
                pipeline.plot_profile(
                    ax,
                    result,
                    color=track['color'],
                    label=track['label'],
                    **track['kwargs'],
                )

            elif track_type == 'custom':
                track['draw_func'](
                    ax,
                    track['data'],
                    region_start,
                    region_end,
                    **track['kwargs'],
                )

            # Set x-limits for all tracks
            ax.set_xlim(region_start, region_end)

        # Format x-axis on bottom panel only
        for ax in axes[:-1]:
            ax.tick_params(labelbottom=False)

        format_genomic_axis(axes[-1], region_start, region_end)

        if title:
            fig.suptitle(title, fontsize=12, fontweight='bold')

        return fig

    def save(
        self,
        path: Union[str, Path],
        region_start: int,
        region_end: int,
        dpi: int = 150,
        bbox_inches: str = 'tight',
        **render_kwargs,
    ):
        """Render and save figure.

        Args:
            path: Output file path
            region_start: Start of genomic region
            region_end: End of genomic region
            dpi: Output resolution
            bbox_inches: Bounding box setting
            **render_kwargs: Additional arguments passed to render()
        """
        fig = self.render(region_start, region_end, **render_kwargs)

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        fig.savefig(path, dpi=dpi, bbox_inches=bbox_inches)
        plt.close(fig)


def create_gene_browser(
    gff_features: Dict,
    gene_name: str,
    region_start: int,
    region_end: int,
    coverage_tracks: Optional[Dict[str, Tuple]] = None,
    vep_data: Optional[Dict[str, Any]] = None,
    figsize: Tuple[float, float] = (12, 8),
) -> plt.Figure:
    """
    Create a gene browser figure with common defaults.

    Convenience function for quick visualization.

    Args:
        gff_features: GFF features dict
        gene_name: Target gene to highlight
        region_start: Start of region
        region_end: End of region
        coverage_tracks: Optional dict of label -> (positions, depths)
        vep_data: Optional dict of panel_type -> DataFrame
        figsize: Figure size

    Returns:
        Matplotlib Figure

    Example:
        fig = create_gene_browser(
            gff_features,
            gene_name='ENA1',
            region_start=530000, region_end=535000,
            coverage_tracks={'NET-seq': (positions, depths)},
            vep_data={'evo2': df_evo2},
        )
        fig.savefig('browser.png')
    """
    builder = MultiTrackFigure(figsize=figsize)

    # Add gene track
    builder.add_gene_track(gff_features, highlight_gene=gene_name)

    # Add coverage tracks
    if coverage_tracks:
        from .config import WONG_COLORS
        colors = list(WONG_COLORS.values())
        for i, (label, data) in enumerate(coverage_tracks.items()):
            builder.add_coverage_track(label, data, color=colors[i % len(colors)])

    # Add VEP tracks
    if vep_data:
        for panel_type, df in vep_data.items():
            builder.add_vep_track(panel_type, df)

    return builder.render(region_start, region_end, gene_name=gene_name)
