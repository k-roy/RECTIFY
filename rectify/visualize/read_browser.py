"""
Stacked read browser for RECTIFY visualization.

Efficient genome-browser-style plots of individual nanopore reads,
using matplotlib LineCollection for fast batch rendering.

Key functions:
    assign_rows              — greedy row packing for non-overlapping layout
    parse_junction_strings   — parse 'start-end,start-end' junction columns
    draw_stacked_reads       — low-level: draw pre-classified reads via LineCollection
    plot_stacked_read_panel  — high-level: assign rows + draw in one call

Design principles:
    - Color assignment is the caller's responsibility (domain-specific)
    - LineCollection groups reads by color: ~6 draw calls for 400 reads
    - Junction-spanning reads are supported (exon blocks + thin intron connectors)
    - constrained_layout-friendly: no tight_layout() calls

Performance vs individual ax.plot() calls:
    400 reads × 5 conditions = 2000 Line2D objects → 30 LineCollection objects

Example (minimal):

    from rectify.visualize import assign_rows, plot_stacked_read_panel

    # Pre-assign colors (caller's responsibility)
    reads_df['color'] = reads_df['category'].map(MY_COLORS)

    # Draw
    fig, ax = plt.subplots(constrained_layout=True)
    n_rows = plot_stacked_read_panel(ax, reads_df)
    ax.set_xlim(region_start, region_end)

Author: Kevin R. Roy
"""

from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

try:
    import matplotlib
    from matplotlib.collections import LineCollection
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    raise ImportError(
        "rectify.visualize.read_browser requires matplotlib. "
        "Install with: pip install rectify-rna[visualize]"
    )


# ---------------------------------------------------------------------------
# Row assignment
# ---------------------------------------------------------------------------

def assign_rows(
    starts: Union[Sequence[int], np.ndarray],
    ends: Union[Sequence[int], np.ndarray],
    gap: int = 50,
) -> np.ndarray:
    """
    Greedy interval row packing: assign each read to the lowest row where
    it fits without overlapping any previously placed read.

    Args:
        starts: Alignment start positions (inclusive), any array-like of ints.
        ends:   Alignment end positions (exclusive or inclusive), same length.
        gap:    Minimum gap (bp) required between adjacent reads in the same row.

    Returns:
        1D int array of row indices, same length as starts/ends.
        Row 0 is the bottom row; higher indices are higher rows.

    Example:
        rows = assign_rows(df['alignment_start'], df['alignment_end'], gap=50)
        n_rows = rows.max() + 1
        ax.set_ylim(-0.5, n_rows - 0.5)
    """
    starts = np.asarray(starts, dtype=int)
    ends = np.asarray(ends, dtype=int)
    n = len(starts)
    rows = np.empty(n, dtype=int)
    row_ends: List[int] = []  # current rightmost end in each row

    for i in range(n):
        a, b = starts[i], ends[i]
        placed = False
        for j, re in enumerate(row_ends):
            if a > re + gap:
                row_ends[j] = b
                rows[i] = j
                placed = True
                break
        if not placed:
            rows[i] = len(row_ends)
            row_ends.append(b)

    return rows


# ---------------------------------------------------------------------------
# Junction string parsing
# ---------------------------------------------------------------------------

def parse_junction_strings(
    junction_col,
) -> List[List[Tuple[int, int]]]:
    """
    Parse a column of junction strings into per-read junction lists.

    Junction strings use the format produced by RECTIFY:
        "start-end,start-end,..."
    where each start-end pair is a genomic intron span (0-based, half-open).

    Args:
        junction_col: iterable of junction strings (one per read). NaN / empty
                      strings produce an empty list for that read.

    Returns:
        List of lists. Each inner list contains (intron_start, intron_end) tuples
        for one read. Reads without junctions get an empty list.

    Example:
        reads_df['junction_lists'] = parse_junction_strings(reads_df['junctions'])
    """
    import pandas as pd

    result = []
    for val in junction_col:
        if not val or (hasattr(val, '__class__') and val.__class__.__name__ == 'float'):
            result.append([])
            continue
        s = str(val).strip()
        if not s:
            result.append([])
            continue
        junctions = []
        for part in s.split(','):
            part = part.strip()
            if '-' in part:
                try:
                    js, je = part.split('-', 1)
                    junctions.append((int(js), int(je)))
                except ValueError:
                    pass
        result.append(junctions)
    return result


# ---------------------------------------------------------------------------
# Core drawing
# ---------------------------------------------------------------------------

def draw_stacked_reads(
    ax,
    starts: Union[Sequence[int], np.ndarray],
    ends: Union[Sequence[int], np.ndarray],
    rows: Union[Sequence[int], np.ndarray],
    colors: Union[Sequence[str], str],
    junction_lists: Optional[List[List[Tuple[int, int]]]] = None,
    read_lw: float = 6.0,
    intron_lw: float = 0.8,
    alpha: float = 0.75,
) -> None:
    """
    Draw stacked reads using LineCollection (fast batch rendering).

    Groups reads by color and draws each color group in a single LineCollection,
    reducing matplotlib artist count from O(n_reads) to O(n_colors).

    Exon blocks are drawn at full read_lw. Intron connectors (for spliced reads)
    are drawn at intron_lw with reduced alpha.

    Args:
        ax:             Matplotlib Axes to draw on.
        starts:         Alignment start positions (bp).
        ends:           Alignment end positions (bp).
        rows:           Row assignments from assign_rows().
        colors:         Color string per read, or single color for all reads.
        junction_lists: Per-read list of (intron_start, intron_end) tuples.
                        None or empty list = unspliced read.
                        Use parse_junction_strings() to build from RECTIFY TSV.
        read_lw:        Linewidth for exon/read bars (default: 6).
        intron_lw:      Linewidth for intron connectors (default: 0.8).
        alpha:          Alpha for exon bars (introns drawn at alpha*0.6).

    Example:
        rows = assign_rows(df['alignment_start'], df['alignment_end'])
        draw_stacked_reads(
            ax,
            df['alignment_start'], df['alignment_end'], rows,
            colors=df['color'],
            junction_lists=parse_junction_strings(df['junctions']),
        )
        ax.set_ylim(-0.5, rows.max() + 0.5)
    """
    n = len(starts)
    if n == 0:
        return

    starts = np.asarray(starts, dtype=int)
    ends = np.asarray(ends, dtype=int)
    rows = np.asarray(rows, dtype=int)

    if isinstance(colors, str):
        colors = [colors] * n

    if junction_lists is None:
        junction_lists = [[] for _ in range(n)]

    # Accumulate exon and intron segments grouped by color
    exon_segs: Dict[str, list] = defaultdict(list)
    intron_segs: Dict[str, list] = defaultdict(list)

    for i in range(n):
        s, e, row, color = int(starts[i]), int(ends[i]), int(rows[i]), colors[i]
        junctions = junction_lists[i] if i < len(junction_lists) else []

        if not junctions:
            exon_segs[color].append([(s, row), (e, row)])
        else:
            # Build exon/intron structure from junction spans
            boundaries = [s]
            for js, je in junctions:
                boundaries.append(js)   # exon end
                boundaries.append(je)   # next exon start
            boundaries.append(e)

            # Exon blocks: pairs at even indices
            exon_blocks = list(zip(boundaries[0::2], boundaries[1::2]))
            for es, ee in exon_blocks:
                exon_segs[color].append([(es, row), (ee, row)])

            # Intron connectors: between exon end and next exon start
            for k in range(1, len(boundaries) - 1, 2):
                intron_segs[color].append(
                    [(boundaries[k], row), (boundaries[k + 1], row)]
                )

    # Draw exon blocks (one LineCollection per color)
    for color, segs in exon_segs.items():
        lc = LineCollection(
            segs,
            colors=color,
            linewidths=read_lw,
            capstyle='butt',
            alpha=alpha,
            zorder=2,
        )
        ax.add_collection(lc)

    # Draw intron connectors (one LineCollection per color, thinner)
    for color, segs in intron_segs.items():
        lc = LineCollection(
            segs,
            colors=color,
            linewidths=intron_lw,
            capstyle='butt',
            alpha=alpha * 0.6,
            zorder=1,
        )
        ax.add_collection(lc)


# ---------------------------------------------------------------------------
# High-level panel function
# ---------------------------------------------------------------------------

def plot_stacked_read_panel(
    ax,
    reads_df,
    color_col: str = 'color',
    start_col: str = 'alignment_start',
    end_col: str = 'alignment_end',
    junction_col: Optional[str] = 'junctions',
    gap: int = 50,
    read_lw: float = 6.0,
    intron_lw: float = 0.8,
    alpha: float = 0.75,
) -> int:
    """
    Assign rows and draw stacked reads in one call.

    Caller is responsible for:
      - Pre-sorting reads_df (e.g. by corrected_3prime)
      - Pre-assigning colors in color_col

    Args:
        ax:           Matplotlib Axes.
        reads_df:     DataFrame with at minimum start_col, end_col, color_col.
        color_col:    Column with read color strings. Must exist in reads_df.
        start_col:    Column with alignment start positions.
        end_col:      Column with alignment end positions.
        junction_col: Column with junction strings (RECTIFY format "js-je,js-je,...").
                      Set to None to skip junction parsing (all reads treated as unspliced).
        gap:          Minimum bp gap between reads in the same row.
        read_lw:      Linewidth for read bars.
        intron_lw:    Linewidth for intron connector lines.
        alpha:        Alpha for read bars.

    Returns:
        n_rows: Number of rows used (set ax.set_ylim(-0.5, n_rows - 0.5) after calling).

    Example:
        reads_df['color'] = reads_df['category'].map(CATEGORY_COLORS)
        n_rows = plot_stacked_read_panel(ax, reads_df, gap=50)
        ax.set_xlim(region_start, region_end)
        ax.set_ylim(-0.5, n_rows - 0.5)
    """
    if reads_df.empty:
        return 0

    starts = reads_df[start_col].to_numpy(dtype=int)
    ends = reads_df[end_col].to_numpy(dtype=int)
    rows = assign_rows(starts, ends, gap=gap)
    n_rows = int(rows.max()) + 1

    colors = reads_df[color_col].tolist()

    if junction_col is not None and junction_col in reads_df.columns:
        junctions = parse_junction_strings(reads_df[junction_col])
    else:
        junctions = None

    draw_stacked_reads(
        ax, starts, ends, rows, colors,
        junction_lists=junctions,
        read_lw=read_lw,
        intron_lw=intron_lw,
        alpha=alpha,
    )

    ax.set_ylim(-0.5, n_rows - 0.5)
    return n_rows
