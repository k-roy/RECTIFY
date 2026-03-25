"""
Figure utility functions for RECTIFY visualization.

Provides helper functions for:
- Publication-style figure setup
- Multi-format figure saving
- Genomic axis formatting
- TRT marker annotations
- Axis styling (despine, etc.)

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import List, Dict, Optional, Union

try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def set_publication_style():
    """Set matplotlib rcParams for publication-quality figures.

    Applies consistent styling:
    - Sans-serif fonts
    - Thicker lines and larger markers
    - White background
    - No top/right spines by default
    """
    if not MATPLOTLIB_AVAILABLE:
        return

    plt.rcParams.update({
        # Font settings
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10,

        # Figure settings
        'figure.dpi': 150,
        'figure.facecolor': 'white',
        'figure.edgecolor': 'white',

        # Axes settings
        'axes.linewidth': 1.0,
        'axes.labelsize': 10,
        'axes.titlesize': 12,
        'axes.facecolor': 'white',
        'axes.edgecolor': 'black',
        'axes.spines.top': False,
        'axes.spines.right': False,

        # Tick settings
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'xtick.major.size': 4,
        'ytick.major.size': 4,

        # Line settings
        'lines.linewidth': 1.5,
        'lines.markersize': 6,

        # Legend settings
        'legend.fontsize': 8,
        'legend.framealpha': 0.8,
        'legend.edgecolor': 'gray',

        # Save settings
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.facecolor': 'white',
        'savefig.edgecolor': 'white',
    })


def save_multi_format(
    fig,
    path: Union[str, Path],
    formats: List[str] = None,
    dpi: int = 150,
    bbox_inches: str = 'tight',
    transparent: bool = False,
):
    """Save figure in multiple formats.

    Args:
        fig: Matplotlib figure object
        path: Base path (without extension) or path with extension
        formats: List of formats to save (default: ['png', 'svg', 'pdf'])
        dpi: Resolution for raster formats
        bbox_inches: Bounding box setting ('tight' recommended)
        transparent: Whether to save with transparent background

    Example:
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3])
        save_multi_format(fig, "results/my_plot", formats=['png', 'pdf'])
        # Creates: results/my_plot.png, results/my_plot.pdf
    """
    if formats is None:
        formats = ['png', 'svg', 'pdf']

    path = Path(path)

    # If path has extension, extract base and use that extension
    if path.suffix:
        base_path = path.with_suffix('')
        if path.suffix[1:] not in formats:
            formats = [path.suffix[1:]] + formats
    else:
        base_path = path

    # Ensure parent directory exists
    base_path.parent.mkdir(parents=True, exist_ok=True)

    for fmt in formats:
        output_path = base_path.with_suffix(f'.{fmt}')
        fig.savefig(
            output_path,
            format=fmt,
            dpi=dpi if fmt in ['png', 'jpg', 'jpeg'] else None,
            bbox_inches=bbox_inches,
            transparent=transparent,
        )


def format_genomic_axis(
    ax,
    region_start: int,
    region_end: int,
    show_unit: bool = True,
    n_ticks: int = 5,
):
    """Format x-axis with genomic coordinates.

    Automatically selects appropriate units (bp, kb, Mb) based on region size.

    Args:
        ax: Matplotlib axes object
        region_start: Start of genomic region
        region_end: End of genomic region
        show_unit: Whether to show unit in axis label
        n_ticks: Approximate number of tick marks
    """
    region_size = region_end - region_start

    # Determine appropriate unit
    if region_size >= 1_000_000:
        unit = 'Mb'
        divisor = 1_000_000
    elif region_size >= 1_000:
        unit = 'kb'
        divisor = 1_000
    else:
        unit = 'bp'
        divisor = 1

    # Set axis limits
    ax.set_xlim(region_start, region_end)

    # Create tick formatter
    def format_tick(x, pos):
        return f'{x / divisor:.1f}' if divisor > 1 else f'{int(x)}'

    ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_tick))

    # Set appropriate tick spacing
    tick_spacing = region_size / n_ticks
    # Round to nice numbers
    magnitude = 10 ** int(f'{tick_spacing:.0e}'.split('e')[1])
    tick_spacing = round(tick_spacing / magnitude) * magnitude
    if tick_spacing == 0:
        tick_spacing = magnitude

    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set axis label
    if show_unit:
        ax.set_xlabel(f'Position ({unit})')


def add_trt_markers(
    ax,
    trt_positions: List[Dict],
    y_frac: float = 0.95,
    color: str = 'red',
    linestyle: str = '--',
    linewidth: float = 1.0,
    alpha: float = 0.7,
    show_labels: bool = True,
    label_fontsize: int = 6,
):
    """Add TRT (T-tract) position markers to a plot.

    Draws red dotted vertical lines at TRT positions with T-tract length labels.

    Args:
        ax: Matplotlib axes object
        trt_positions: List of dicts with 'position' and optionally 't_tract_length'
        y_frac: Vertical position for labels (0-1, fraction of y-axis)
        color: Line color
        linestyle: Line style ('--' for dashed, ':' for dotted)
        linewidth: Line width
        alpha: Line transparency
        show_labels: Whether to show T-tract length labels
        label_fontsize: Font size for labels

    Example:
        trts = [
            {'position': 54321, 't_tract_length': 8},
            {'position': 54500, 't_tract_length': 12},
        ]
        add_trt_markers(ax, trts)
    """
    ymin, ymax = ax.get_ylim()
    y_label = ymin + (ymax - ymin) * y_frac

    for trt in trt_positions:
        pos = trt.get('position')
        if pos is None:
            continue

        # Draw vertical line
        ax.axvline(
            x=pos,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            alpha=alpha,
            zorder=1,
        )

        # Add label if T-tract length is provided
        if show_labels and 't_tract_length' in trt:
            t_len = trt['t_tract_length']
            ax.text(
                pos,
                y_label,
                f'T{t_len}',
                fontsize=label_fontsize,
                color=color,
                ha='center',
                va='bottom',
                alpha=alpha,
            )


def despine(ax, which: Optional[List[str]] = None):
    """Remove specified spines from axes.

    Args:
        ax: Matplotlib axes object
        which: List of spines to remove. Default: ['top', 'right']
               Options: 'top', 'right', 'bottom', 'left'

    Example:
        despine(ax)  # Remove top and right spines
        despine(ax, ['top', 'right', 'bottom'])  # Keep only left spine
    """
    if which is None:
        which = ['top', 'right']

    for spine in which:
        ax.spines[spine].set_visible(False)


def add_significance_bracket(
    ax,
    x1: float,
    x2: float,
    y: float,
    text: str,
    bracket_height: float = 0.02,
    fontsize: int = 8,
    color: str = 'black',
):
    """Add significance bracket with text (e.g., p-value annotation).

    Args:
        ax: Matplotlib axes object
        x1: Left x-coordinate
        x2: Right x-coordinate
        y: Y-coordinate for bracket
        text: Text to display (e.g., '***', 'p < 0.001')
        bracket_height: Height of bracket ends as fraction of y-range
        fontsize: Font size for text
        color: Color for bracket and text

    Example:
        add_significance_bracket(ax, 1, 2, 5.5, '***')
    """
    # Get y-axis range for scaling
    ymin, ymax = ax.get_ylim()
    bracket_dy = (ymax - ymin) * bracket_height

    # Draw bracket
    ax.plot(
        [x1, x1, x2, x2],
        [y - bracket_dy, y, y, y - bracket_dy],
        color=color,
        linewidth=1,
    )

    # Add text
    ax.text(
        (x1 + x2) / 2,
        y + bracket_dy / 2,
        text,
        ha='center',
        va='bottom',
        fontsize=fontsize,
        color=color,
    )


def create_figure_grid(
    n_panels: int,
    n_cols: int = 1,
    panel_height: float = 2.0,
    panel_width: float = 8.0,
    height_ratios: Optional[List[float]] = None,
    hspace: float = 0.3,
    sharex: bool = True,
):
    """Create a figure with multiple panel rows.

    Convenient wrapper for creating multi-panel figures with shared x-axis.

    Args:
        n_panels: Number of panels
        n_cols: Number of columns
        panel_height: Height of each panel in inches
        panel_width: Width of figure in inches
        height_ratios: Optional list of relative heights for each row
        hspace: Vertical spacing between panels
        sharex: Whether panels share x-axis

    Returns:
        (fig, axes): Figure and array of axes objects

    Example:
        fig, axes = create_figure_grid(4, panel_height=1.5)
        axes[0].plot(...)  # Gene track
        axes[1].plot(...)  # Coverage
        axes[2].plot(...)  # VEP scores
        axes[3].plot(...)  # Metagene
    """
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError("matplotlib is required for create_figure_grid")

    n_rows = (n_panels + n_cols - 1) // n_cols

    if height_ratios is None:
        height_ratios = [1] * n_rows

    figsize = (panel_width, sum(height_ratios) * panel_height)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=figsize,
        gridspec_kw={'height_ratios': height_ratios, 'hspace': hspace},
        sharex=sharex if n_cols == 1 else False,
        squeeze=False,
    )

    # Flatten for easier access
    axes = axes.flatten()

    # Hide unused axes
    for i in range(n_panels, len(axes)):
        axes[i].set_visible(False)

    return fig, axes[:n_panels]
