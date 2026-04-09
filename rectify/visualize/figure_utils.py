"""
Figure utility functions for RECTIFY visualization.

Provides helper functions for:
- Publication-style figure setup
- Multi-format figure saving
- Genomic axis formatting
- TRT marker annotations
- Axis styling (despine, etc.)
- Metagene plotting utilities (trimmed mean, ridge plots, etc.)

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import List, Dict, Optional, Union, Tuple

import numpy as np

try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


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


# =============================================================================
# Metagene Plotting Utilities
# =============================================================================

def trimmed_mean(arr: np.ndarray, proportion: float = 0.1, axis: int = 0) -> np.ndarray:
    """Calculate trimmed mean along axis, excluding outlier values.

    Uses scipy.stats.trim_mean to exclude the top and bottom proportion of
    values at each position. This is useful for metagene aggregation where
    you want to reduce the influence of outlier loci without excluding them
    entirely.

    Args:
        arr: Input array (typically shape: n_loci x window_size)
        proportion: Proportion to cut from each end (0.1 = 10% from each end,
                    keeping middle 80%). Must be in [0, 0.5).
        axis: Axis along which to compute trimmed mean (default: 0 = across loci)

    Returns:
        Array with trimmed mean values (shape depends on axis)

    Example:
        # Profiles shape: (1000 loci, 101 positions)
        profiles = np.random.randn(1000, 101)

        # 10% trimmed mean at each position
        mean_profile = trimmed_mean(profiles, proportion=0.1, axis=0)
        # Result shape: (101,)
    """
    if not SCIPY_AVAILABLE:
        raise ImportError("scipy is required for trimmed_mean. Install with: pip install scipy")

    return stats.trim_mean(arr, proportion, axis=axis)


def apply_window_sum_capping(
    profiles: np.ndarray,
    percentile: float = 90,
) -> Tuple[np.ndarray, float]:
    """Apply window-sum capping to normalize profiles by scaling outlier loci.

    For each locus, calculates the sum of signal across the window. Loci with
    sums above the specified percentile are scaled down proportionally. This
    preserves the shape of high-signal loci while reducing their contribution
    to the aggregate.

    This is an alternative to trimmed_mean that keeps all loci but reduces
    outlier influence. Useful when you want to retain information from all
    loci but prevent extreme values from dominating the metagene.

    Args:
        profiles: Array of shape (n_loci, window_size)
        percentile: Percentile threshold for capping (default: 90)

    Returns:
        Tuple of (capped_profiles, cap_value) where:
        - capped_profiles: Array with scaled values
        - cap_value: The threshold value used for capping

    Example:
        profiles = np.random.exponential(1, size=(1000, 101))
        capped, cap_val = apply_window_sum_capping(profiles, percentile=90)
        print(f"Capped at window sum = {cap_val:.2f}")
    """
    window_sums = np.sum(profiles, axis=1)
    nonzero = window_sums[window_sums > 0]

    if len(nonzero) == 0:
        return profiles, 0.0

    cap_val = np.percentile(nonzero, percentile)
    scale = np.where(window_sums > cap_val, cap_val / window_sums, 1.0)
    scale = np.where(window_sums == 0, 1.0, scale)

    return profiles * scale[:, np.newaxis], cap_val


def plot_ridge_profiles(
    ax,
    x_positions: np.ndarray,
    profiles: Dict[str, np.ndarray],
    colors: Dict[str, str],
    labels: Dict[str, str],
    order: Optional[List[str]] = None,
    aggregation: str = 'trimmed_mean',
    trim_proportion: float = 0.1,
    offset_scale: float = 0.8,
    fill_alpha: float = 0.6,
    line_alpha: float = 1.0,
    linewidth: float = 1.5,
) -> Dict[str, float]:
    """Plot ridge-style stacked metagene profiles.

    Creates vertically stacked filled profiles where each condition is offset
    vertically. This style is useful for comparing profile shapes across
    conditions while avoiding overlap.

    Args:
        ax: Matplotlib axes object
        x_positions: Array of x-axis positions (e.g., np.arange(-50, 51))
        profiles: Dict mapping condition name -> profile array (n_loci, window_size)
        colors: Dict mapping condition name -> color
        labels: Dict mapping condition name -> display label
        order: List of condition names in bottom-to-top order. If None, uses
               keys from profiles dict.
        aggregation: Aggregation method ('trimmed_mean', 'mean', 'median')
        trim_proportion: Proportion for trimmed_mean (default: 0.1)
        offset_scale: Vertical spacing between conditions as fraction of max signal
        fill_alpha: Alpha for filled region
        line_alpha: Alpha for outline
        linewidth: Line width for outline

    Returns:
        Dict mapping condition name -> y_offset used for that condition

    Example:
        from rectify.visualize import WONG_COLORS, plot_ridge_profiles

        profiles = {
            'wt': wt_profiles,      # shape (n_loci, window_size)
            'mutant': mut_profiles,
        }
        colors = {'wt': WONG_COLORS['blue'], 'mutant': WONG_COLORS['orange']}
        labels = {'wt': 'WT', 'mutant': 'Mutant'}

        fig, ax = plt.subplots()
        offsets = plot_ridge_profiles(
            ax, x_positions, profiles, colors, labels,
            order=['mutant', 'wt'],  # mutant on bottom, wt on top
        )
    """
    if order is None:
        order = list(profiles.keys())

    # Compute aggregated profiles
    agg_profiles = {}
    for key, prof in profiles.items():
        if aggregation == 'trimmed_mean':
            agg_profiles[key] = trimmed_mean(prof, proportion=trim_proportion, axis=0)
        elif aggregation == 'mean':
            agg_profiles[key] = np.mean(prof, axis=0)
        elif aggregation == 'median':
            agg_profiles[key] = np.median(prof, axis=0)
        else:
            raise ValueError(f"Unknown aggregation: {aggregation}. Use 'trimmed_mean', 'mean', or 'median'.")

    # Calculate offset based on max signal
    max_signal = max(agg.max() for agg in agg_profiles.values())
    offset_increment = max_signal * offset_scale

    offsets = {}
    for i, key in enumerate(order):
        y_offset = i * offset_increment
        offsets[key] = y_offset

        mean = agg_profiles[key]
        color = colors.get(key, 'gray')
        label = labels.get(key, key)

        ax.fill_between(
            x_positions,
            y_offset,
            mean + y_offset,
            color=color,
            alpha=fill_alpha,
        )
        ax.plot(
            x_positions,
            mean + y_offset,
            color=color,
            linewidth=linewidth,
            alpha=line_alpha,
            label=label,
        )

    return offsets


def add_metagene_annotations(
    ax,
    x_positions: np.ndarray,
    highlight_region: Optional[Tuple[float, float]] = None,
    highlight_color: str = 'yellow',
    highlight_alpha: float = 0.15,
    grid_range: Optional[Tuple[int, int]] = None,
    grid_color: str = 'gray',
    grid_alpha: float = 0.5,
    zero_line_width: float = 1.0,
    minor_line_width: float = 0.3,
):
    """Add highlight regions and reference grid lines to metagene plots.

    Commonly used to highlight regions of interest (e.g., T-tract region) and
    add vertical reference lines for alignment.

    Args:
        ax: Matplotlib axes object
        x_positions: Array of x-axis positions (used for xlim if not set)
        highlight_region: Tuple (start, end) for highlighted region (e.g., (0, 10))
        highlight_color: Color for highlight region
        highlight_alpha: Alpha for highlight region
        grid_range: Tuple (start, end) for vertical grid lines (e.g., (-10, 10))
        grid_color: Color for grid lines
        grid_alpha: Alpha for grid lines
        zero_line_width: Line width for the line at x=0
        minor_line_width: Line width for other grid lines

    Example:
        # Add T-tract highlight (positions 0-10) and grid lines
        add_metagene_annotations(
            ax, x_positions,
            highlight_region=(0, 10),
            grid_range=(-10, 10),
        )
    """
    # Add highlight region
    if highlight_region is not None:
        ax.axvspan(
            highlight_region[0],
            highlight_region[1],
            color=highlight_color,
            alpha=highlight_alpha,
            zorder=0,
        )

    # Add grid lines
    if grid_range is not None:
        for pos in range(grid_range[0], grid_range[1] + 1):
            lw = zero_line_width if pos == 0 else minor_line_width
            ax.axvline(
                pos,
                color=grid_color,
                linewidth=lw,
                linestyle=':',
                alpha=grid_alpha,
                zorder=1,
            )


def plot_metagene_line(
    ax,
    x_positions: np.ndarray,
    profiles: np.ndarray,
    color: str,
    label: str,
    aggregation: str = 'trimmed_mean',
    trim_proportion: float = 0.1,
    show_sem: bool = True,
    sem_alpha: float = 0.3,
    linewidth: float = 1.5,
) -> Tuple[np.ndarray, np.ndarray]:
    """Plot a single metagene profile with optional SEM shading.

    Useful for overlaid line-style metagene plots where conditions are
    compared on the same y-axis (as opposed to ridge-style stacking).

    Args:
        ax: Matplotlib axes object
        x_positions: Array of x-axis positions
        profiles: Array of shape (n_loci, window_size)
        color: Line and fill color
        label: Legend label
        aggregation: Aggregation method ('trimmed_mean', 'mean')
        trim_proportion: Proportion for trimmed_mean
        show_sem: Whether to show SEM shading
        sem_alpha: Alpha for SEM shading
        linewidth: Line width

    Returns:
        Tuple of (mean_profile, sem_profile)

    Example:
        mean, sem = plot_metagene_line(
            ax, x_positions, wt_profiles,
            color=WONG_COLORS['blue'],
            label='WT (n=3)',
        )
    """
    if aggregation == 'trimmed_mean':
        mean = trimmed_mean(profiles, proportion=trim_proportion, axis=0)
        n = profiles.shape[0]
        k = int(n * trim_proportion)
        n_eff = n - 2 * k
        trimmed_profiles = np.sort(profiles, axis=0)[k:n - k, :]
        sem = np.std(trimmed_profiles, axis=0, ddof=1) / np.sqrt(n_eff)
    elif aggregation == 'mean':
        mean = np.mean(profiles, axis=0)
        n_eff = profiles.shape[0]
        sem = np.std(profiles, axis=0, ddof=1) / np.sqrt(n_eff)
    else:
        raise ValueError(f"Unknown aggregation: {aggregation}")

    if show_sem:
        ax.fill_between(
            x_positions,
            mean - sem,
            mean + sem,
            color=color,
            alpha=sem_alpha,
        )

    ax.plot(
        x_positions,
        mean,
        color=color,
        linewidth=linewidth,
        label=label,
    )

    return mean, sem
