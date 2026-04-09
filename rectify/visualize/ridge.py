"""
Ridge plot utilities for RECTIFY visualization.

Provides functions for creating overlapping ridge-style figures where each
condition/group has its own panel stacked vertically with configurable overlap.
Useful for comparing signal profile shapes across conditions without cluttering
a shared y-axis.

IMPORTANT: Do NOT use constrained_layout=True with overlapping ridge figures —
it fights negative hspace and will un-overlap the panels.

Author: Kevin R. Roy
"""

from typing import List, Optional, Tuple

import numpy as np

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    from scipy.ndimage import gaussian_filter1d
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


def aggregate_replicate_profiles(
    profiles: List[np.ndarray],
    smooth_sigma: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute mean ± SEM from per-replicate 1D profile arrays.

    Optionally applies Gaussian smoothing to each replicate before aggregating.
    Smoothing before aggregation preserves sharp peaks better than smoothing
    the mean profile after aggregation.

    Args:
        profiles: List of 1D arrays (one per replicate), all the same length.
        smooth_sigma: If given, apply Gaussian smoothing with this sigma (in
                      position units) to each replicate before computing mean/SEM.
                      Requires scipy.

    Returns:
        (mean, sem): Two 1D arrays of the same length as each input profile.
                     sem is zeros if only one replicate is provided.

    Raises:
        ValueError: If profiles is empty.
        ImportError: If smooth_sigma is given but scipy is not installed.

    Example:
        # Smooth each replicate at sigma=3, then aggregate
        mean, sem = aggregate_replicate_profiles(raw_arrays, smooth_sigma=3)

        # Or use pre-smoothed arrays
        smoothed = [gaussian_filter1d(y, sigma=3) for y in raw_arrays]
        mean, sem = aggregate_replicate_profiles(smoothed)
    """
    if not profiles:
        raise ValueError("profiles list is empty")

    if smooth_sigma is not None and not SCIPY_AVAILABLE:
        raise ImportError(
            "scipy is required for smooth_sigma. Install with: pip install scipy"
        )

    arrays = []
    for p in profiles:
        arr = np.asarray(p, dtype=float)
        if smooth_sigma is not None:
            arr = gaussian_filter1d(arr, sigma=smooth_sigma)
        arrays.append(arr)

    mat = np.array(arrays)  # shape: (n_replicates, window_size)
    n = len(arrays)
    mean = mat.mean(axis=0)
    sem = mat.std(axis=0, ddof=1) / np.sqrt(n) if n > 1 else np.zeros_like(mean)
    return mean, sem


def plot_ridge_panel(
    ax,
    x: np.ndarray,
    mean: np.ndarray,
    sem: Optional[np.ndarray] = None,
    color: str = '#0072B2',
    fill_alpha: float = 0.55,
    sem_alpha: float = 0.30,
    linewidth: float = 1.2,
    label: Optional[str] = None,
    label_loc: str = 'upper left',
    label_fontsize: int = 8,
    hide_spines: bool = True,
    xlim: Optional[Tuple[float, float]] = None,
) -> None:
    """Draw one overlapping ridge panel: filled area, optional SEM band, mean line.

    Designed to work on axes from create_overlapping_ridge_figure(). The y-axis
    is hidden; panels are compared by visual peak height and shape.

    Args:
        ax: Matplotlib axes (one of the overlapping ridge panels).
        x: x-axis positions array.
        mean: Mean profile (same length as x).
        sem: SEM profile (same length as x). If None, no SEM band is drawn.
        color: Fill and line color.
        fill_alpha: Alpha for main fill (0 to mean).
        sem_alpha: Alpha for SEM band (mean ± sem).
        linewidth: Width of mean line.
        label: Row label string placed inside the panel.
        label_loc: 'upper left' or 'upper right'.
        label_fontsize: Font size for label.
        hide_spines: Remove top, right, left spines (standard for ridge style).
        xlim: (xmin, xmax) limits. If None, uses first and last x values.

    Notes:
        Each panel is drawn with a transparent background so that tall peaks can
        extend visually through the panel above rather than being clipped by its
        white fill. A white horizontal line at y=0 acts as the visual horizon
        separating adjacent panels.

    Example:
        fig, axes, ax_track = create_overlapping_ridge_figure(3, add_gene_track=True)
        x = np.arange(-50, 51)
        mean, sem = aggregate_replicate_profiles(replicate_arrays, smooth_sigma=3)
        plot_ridge_panel(axes[0], x, mean, sem=sem, color='#D55E00', label='WT')
        axes[0].set_xlim(-50, 50)
    """
    x = np.asarray(x)
    mean = np.asarray(mean)

    # Main fill from 0 to mean
    ax.fill_between(x, 0, mean, alpha=fill_alpha, color=color)

    # SEM band: mean ± sem
    if sem is not None:
        sem = np.asarray(sem)
        ax.fill_between(x, mean - sem, mean + sem, alpha=sem_alpha, color=color)

    # Mean line
    ax.plot(x, mean, color=color, linewidth=linewidth)

    # Row label — placed in the LOWER portion of each panel (y=0.20) so it
    # falls in the visible (non-overlap) zone when panels stack with negative
    # hspace.  clip_on=False lets the label render even if the axes bounding
    # box is tight against the data.
    if label is not None:
        ha = 'left' if label_loc == 'upper left' else 'right'
        x_frac = 0.01 if label_loc == 'upper left' else 0.99
        ax.text(x_frac, 0.20, label, transform=ax.transAxes,
                fontsize=label_fontsize, va='bottom', ha=ha, clip_on=False)

    # Transparent background: allows tall peaks to extend visually through the
    # panel above in overlapping ridge layouts instead of being clipped by
    # the upper panel's white background.
    ax.set_facecolor((1, 1, 1, 0))

    # White baseline separator: draws a white horizontal line at y=0 that
    # acts as the visual horizon between overlapping panels, clearly separating
    # each panel's signal from the one below.
    ax.axhline(0, color='white', lw=1.5, zorder=5)

    # Spine and tick styling — hide ALL spines to avoid the bottom-spine
    # cutoff artefact in overlapping ridge layouts
    if hide_spines:
        for spine in ax.spines.values():
            spine.set_visible(False)
    ax.set_yticks([])

    # Clip y-axis at 0: prevents SEM band below zero from stretching ylim
    # downward and creating visual bleed into the panel below
    ax.set_ylim(bottom=0)

    # x-axis limits
    if xlim is not None:
        ax.set_xlim(xlim)
    elif len(x) > 0:
        ax.set_xlim(x[0], x[-1])


def create_overlapping_ridge_figure(
    n_panels: int,
    panel_height: float = 1.6,
    panel_width: float = 9.0,
    hspace: float = -0.40,
    bottom_margin: float = 0.20,
    add_gene_track: bool = False,
    gene_track_height_frac: float = 0.10,
) -> Tuple['plt.Figure', List['plt.Axes'], Optional['plt.Axes']]:
    """Create a figure with n overlapping ridge panels and an optional gene track.

    Panels share the x-axis and overlap vertically via negative hspace,
    producing the classic "ridge plot" waterfall appearance. The y-axis is
    hidden on all panels; only the shape and relative height of the signal
    is meaningful.

    IMPORTANT: Do NOT use constrained_layout=True — it fights the negative
    hspace and will un-overlap the panels.

    Args:
        n_panels: Number of ridge panels (one per condition/group).
        panel_height: Height of each panel in inches.
        panel_width: Width of figure in inches.
        hspace: Vertical overlap between panels. Negative values create
                overlap (default -0.40 ≈ 40% overlap). Use 0 for no overlap.
        bottom_margin: Figure bottom fraction reserved for the gene track
                       and x-axis label.
        add_gene_track: If True, add a gene track axes below the ridge panels.
        gene_track_height_frac: Height of gene track as fraction of total
                                 figure height.

    Returns:
        (fig, axes_list, ax_gene_track)
        - fig: Figure object (do NOT call fig.tight_layout() or set
               constrained_layout=True).
        - axes_list: List of n axes, ordered top-to-bottom. axes_list[0]
                     is the topmost panel.
        - ax_gene_track: Gene track axes if add_gene_track=True, else None.
                         Positioned with fig.add_axes() to avoid constrained-
                         layout conflicts. The caller must set xlim.

    Example:
        from rectify.visualize import (create_overlapping_ridge_figure,
                                        plot_ridge_panel, aggregate_replicate_profiles,
                                        WONG_COLORS, save_multi_format)

        groups = [('WT', '#0072B2'), ('mutant', '#D55E00')]
        fig, axes, ax_track = create_overlapping_ridge_figure(
            n_panels=len(groups), add_gene_track=True)

        x = np.arange(-50, 51)
        for i, (label, color) in enumerate(groups):
            mean, sem = aggregate_replicate_profiles(replicate_data[label], smooth_sigma=3)
            plot_ridge_panel(axes[i], x, mean, sem=sem, color=color, label=label)
            axes[i].set_xlim(-50, 50)

        if ax_track:
            ax_track.set_xlim(-50, 50)
            # ... annotate gene features ...

        save_multi_format(fig, 'plots/my_ridge')
    """
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError(
            "matplotlib is required for create_overlapping_ridge_figure. "
            "Install with: pip install rectify-rna[visualize]"
        )

    total_height = panel_height * n_panels + 2.0
    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(panel_width, total_height),
        sharex=True,
    )

    # Ensure axes is always a list
    if n_panels == 1:
        axes_list = [axes]
    else:
        axes_list = list(axes)

    # Apply negative hspace for overlapping effect
    fig.subplots_adjust(hspace=hspace, bottom=bottom_margin)

    # Reverse z-order so topmost panel (index 0) renders last (in front).
    # Without this, the bottommost panel's white background clips the baseline
    # of all panels above it in the overlap zone.
    for i, ax in enumerate(axes_list):
        ax.set_zorder(n_panels - i)

    ax_gene_track = None
    if add_gene_track:
        # Manually position gene track in the bottom margin using fig.add_axes()
        # so it doesn't interfere with the ridge panels' hspace setting.
        # Coordinates: [left, bottom, width, height] in figure fractions.
        # Match the default left/right margins of subplots (approximately 0.125 and 0.9)
        left = 0.125
        width = 0.775
        track_bottom = 0.03
        track_height = gene_track_height_frac
        ax_gene_track = fig.add_axes([left, track_bottom, width, track_height])

    return fig, axes_list, ax_gene_track
