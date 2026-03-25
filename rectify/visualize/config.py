"""
Configuration constants for RECTIFY visualization.

This module consolidates all color schemes, marker shapes, and other
visualization constants used across RNA 3' end plotting.

Includes:
- Codon variant type colors (for VEP panels)
- Gene type colors (for per-gene expression panels)
- Strain and sublibrary markers
- Publication-ready figure configuration

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple

# ============================================================================
# Wong Colorblind-Safe Palette
# ============================================================================
# From: Wong, B. (2011). Points of view: Color blindness. Nature Methods.
# Accessible palette for colorblind viewers

WONG_COLORS: Dict[str, str] = {
    'blue': '#0072B2',
    'orange': '#E69F00',
    'green': '#009E73',
    'pink': '#CC79A7',
    'yellow': '#F0E442',
    'cyan': '#56B4E9',
    'red': '#D55E00',
    'black': '#000000',
}


# ============================================================================
# Codon Variant Type Colors
# ============================================================================
# Used for Evo2, ESM1-v, and variant map panels
# Colors match publication style

CODON_VARIANT_COLORS: Dict[str, str] = {
    # Protein-coding changes
    'missense': '#3498db',           # Blue
    'synonymous': '#2ecc71',         # Green
    'stop_gained': '#e74c3c',        # Red
    'stop_gain': '#e74c3c',          # Red (alias)
    'frameshift': '#9b59b6',         # Purple
    'stop_lost': '#f39c12',          # Orange
    'stop_loss': '#f39c12',          # Orange (alias)
    'inframe_insertion': '#1abc9c',  # Teal
    'inframe_deletion': '#e67e22',   # Dark orange
    'inframe_indel': '#17a2b8',      # Cyan
    'in-frame_indel': '#17a2b8',     # Cyan (alias)
    'start_lost': '#c0392b',         # Dark red
    'start_loss': '#c0392b',         # Dark red (alias)
    'splice_region': '#8e44ad',      # Dark purple
    # Non-coding
    'intergenic': '#95a5a6',         # Gray
    'intronic': '#7f8c8d',           # Dark gray
    'non_coding': '#bdc3c7',         # Light gray
    'non-coding': '#bdc3c7',         # Light gray (alias)
    # Special
    'complete_CDS_deletion': '#2c3e50',  # Dark blue-gray
    'unknown': '#95a5a6',            # Gray
}


# ============================================================================
# Gene Type Colors
# ============================================================================
# Used for Yorzoi/Shorkie per-gene expression panels
# Colors match box-arrow gene track for visual consistency

GENE_TYPE_COLORS: Dict[str, str] = {
    'target': '#e74c3c',             # Red - target gene (most important)
    'upstream_1': '#3498db',         # Blue - immediate upstream
    'upstream_2': '#1abc9c',         # Teal - second upstream
    'downstream_1': '#9b59b6',       # Purple - immediate downstream
    'downstream_2': '#f39c12',       # Orange - second downstream
}


# ============================================================================
# Strain Origin Colors
# ============================================================================
# Used for variant map panel to show strain distribution

STRAIN_ORIGIN_COLORS: Dict[str, str] = {
    'shared_only': '#3498db',        # Blue - only in shared parent
    'shared_plus_left': '#9b59b6',   # Purple - shared + left neighbor
    'shared_plus_right': '#e67e22',  # Orange - shared + right neighbor
    'both_neighbors': '#e74c3c',     # Red - both neighbors, not shared
    'left_only': '#1abc9c',          # Teal - only left neighbor (filtered)
    'right_only': '#f39c12',         # Yellow - only right neighbor (filtered)
    'all_three': '#7f8c8d',          # Gray - in all three (filtered)
    'other_strains': '#bdc3c7',      # Light gray - other combinations
}


# ============================================================================
# Control Type Colors
# ============================================================================
# Used for MAGESTIC control panels (PTC, SYN, DEL)

CONTROL_TYPE_COLORS: Dict[str, str] = {
    'PTC': '#e74c3c',                # Red - premature termination codon
    'SYN': '#2ecc71',                # Green - synonymous
    'DEL': '#2c3e50',                # Dark blue-gray - CDS deletion
}


# ============================================================================
# Effect Direction Colors
# ============================================================================
# Used for trio dissection and other effect-based plots

EFFECT_DIRECTION_COLORS: Dict[str, str] = {
    'positive': '#e74c3c',           # Red - positive effect
    'negative': '#3498db',           # Blue - negative effect
    'neutral': '#95a5a6',            # Gray - no effect
}


# ============================================================================
# Sublibrary Marker Shapes
# ============================================================================
# Used for MAGESTIC screen panels

SUBLIBRARY_MARKERS: Dict[str, str] = {
    'SpG_NGG': 'o',                  # Circle
    'SpG_NGNG': 'v',                 # Triangle down
    'SpG_NGNH': 'D',                 # Diamond
    'SpCas9_NGG': 's',               # Square
    'SpG': 'o',                      # Default SpG
    'SpCas9': 's',                   # Default SpCas9
    'NGG': '^',                      # Triangle up
    'NGNG': 'v',                     # Triangle down
    'NGNH': 'D',                     # Diamond
    'default': 'o',
}


# ============================================================================
# Strain Origin Marker Shapes
# ============================================================================
# Used for variant map panel

STRAIN_ORIGIN_MARKERS: Dict[str, str] = {
    'all_three': 'o',                # Circle - in all 3 strains
    'shared_only': 's',              # Square - only shared parent
    'shared_plus_left': '^',         # Triangle up - shared + left
    'shared_plus_right': 'v',        # Triangle down - shared + right
    'both_neighbors': 'D',           # Diamond - both neighbors, not shared
    'left_only': '<',                # Left arrow - only left neighbor
    'right_only': '>',               # Right arrow - only right neighbor
    'other_strains': 'x',            # X - other combinations
}


# ============================================================================
# Strain Background Marker Shapes
# ============================================================================
# Used for trio dissection panels to distinguish strain backgrounds
# Shape encodes strain while color encodes variant type

STRAIN_MARKERS: Dict[str, str] = {
    # Key strains get the most distinct markers
    "BYa": "o",                      # Circle - MAGESTIC reference
    "RMx": "s",                      # Square
    "YJM981x": "^",                  # Triangle up
    "YJM978x": "v",                  # Triangle down
    "YJM454a": "D",                  # Diamond
    "YJM145x": "P",                  # Plus (filled)
    "273614xa": "*",                 # Star
    "CBS2888a": "X",                 # X (filled)
    "CLIB219x": "h",                 # Hexagon
    "CLIB413a": "p",                 # Pentagon
    "I14a": "<",                     # Triangle left
    "M22": ">",                      # Triangle right
    "PW5a": "8",                     # Octagon
    "Y10x": "H",                     # Hexagon2
    "YPS1009x": "d",                 # Thin diamond
    "YPS163a": "+",                  # Plus
}


# ============================================================================
# Gene Type Marker Shapes
# ============================================================================
# Alternative to strain shapes for per-gene panels

GENE_TYPE_MARKERS: Dict[str, str] = {
    'target': 'o',                   # Circle
    'upstream_1': '^',               # Triangle up
    'upstream_2': 'v',               # Triangle down
    'downstream_1': 's',             # Square
    'downstream_2': 'D',             # Diamond
}


# ============================================================================
# Variant Source Markers (Trio Dissection)
# ============================================================================
# Shape coding for Column 3 (shared parent background)

SOURCE_MARKERS: Dict[str, str] = {
    'insert_left': 'o',              # Circle - from left neighbor
    'insert_right': '^',             # Triangle - from right neighbor
    'insert_shared': 'D',            # Diamond - shared by both neighbors
}


# ============================================================================
# Column Background Colors (Trio Dissection)
# ============================================================================
# Subtle background colors for 4-column layout

COLUMN_COLORS: Dict[str, str] = {
    'magestic': '#f0f8f0',           # Very light green
    'left': '#f0f0f8',               # Very light blue
    'shared': '#f8f0f0',             # Very light red
    'right': '#f8f8f0',              # Very light yellow
}


# ============================================================================
# VEP Model Names and Score Columns
# ============================================================================

MODEL_NAMES: Dict[str, str] = {
    "evo2": "Evo2 40B",
    "yorzoi": "Yorzoi",
    "shorkie": "Shorkie",
    "esm1v": "ESM1-v",
}

MODEL_SCORE_COLUMNS: Dict[str, List[str]] = {
    'evo2': ['evo2_score', 'delta_score', 'raw_score'],
    'esm1v': ['esm1v_delta_logits', 'delta_logits'],
    'yorzoi': ['yorzoi_expression_ratio', 'log2_fold_change'],
    'shorkie': ['shorkie_logSED', 'delta_score'],
}


# ============================================================================
# Pilot Genes
# ============================================================================
# 8 genes for initial trio dissection testing

PILOT_GENES: Dict[str, str] = {
    "YJR066W": "TOR1",
    "YGL167C": "PMR1",
    "YDR040C": "ENA1",
    "YNL085W": "MKT1",
    "YNL083W": "SAL1",
    "YJL005W": "CYR1",
    "YDL035C": "GPR1",
    "YOL081W": "IRA2",
}

# Reverse mapping
PILOT_ORF_BY_GENE: Dict[str, str] = {v: k for k, v in PILOT_GENES.items()}


# ============================================================================
# Figure Configuration
# ============================================================================

@dataclass
class FigureConfig:
    """Configuration for publication-quality figures.

    Usage:
        config = FigureConfig(dpi=300, figsize=(12, 8))
        config.apply()  # Sets matplotlib rcParams
    """
    dpi: int = 150
    figsize: Tuple[float, float] = (10, 6)
    font_family: str = 'sans-serif'
    title_fontsize: int = 12
    label_fontsize: int = 10
    tick_fontsize: int = 8
    legend_fontsize: int = 8
    line_width: float = 1.5
    marker_size: int = 6

    def apply(self):
        """Apply settings to matplotlib rcParams."""
        try:
            import matplotlib.pyplot as plt
            plt.rcParams.update({
                'figure.dpi': self.dpi,
                'figure.figsize': self.figsize,
                'font.family': self.font_family,
                'font.size': self.label_fontsize,
                'axes.titlesize': self.title_fontsize,
                'axes.labelsize': self.label_fontsize,
                'xtick.labelsize': self.tick_fontsize,
                'ytick.labelsize': self.tick_fontsize,
                'legend.fontsize': self.legend_fontsize,
                'lines.linewidth': self.line_width,
                'lines.markersize': self.marker_size,
            })
        except ImportError:
            pass

    @classmethod
    def publication(cls) -> 'FigureConfig':
        """Return publication-ready configuration (300 DPI)."""
        return cls(dpi=300, figsize=(8, 6))

    @classmethod
    def presentation(cls) -> 'FigureConfig':
        """Return presentation-ready configuration (larger fonts)."""
        return cls(
            dpi=150,
            figsize=(12, 8),
            title_fontsize=16,
            label_fontsize=14,
            tick_fontsize=12,
            legend_fontsize=12,
        )
