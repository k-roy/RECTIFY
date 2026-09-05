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
# Used for Evo2, ESM1-v, and variant map panels.
#
# 🔴 REPLACED 2026-09-04 (Chanfreau planning/871). The previous values were the
# flat-UI web palette (#3498db / #2ecc71 / #e74c3c / #9b59b6 ...). They were never
# colour-vision validated, and they assigned RED to `stop_gained` and GREEN to
# `synonymous` — the single most consequential distinction in the panel, drawn as
# exactly the red/green pair that ~8% of male readers cannot separate. It is also
# the known-bad case that palette_check.py plants in its own --self-test.
#
# These values are the eight MAGESTIC semantic tokens, which are already proven
# mutually distinguishable at dE2000 >= 12 under normal vision, deuteranopia and
# protanopia (tokens.json v1.1.0, validated by
# ~/.claude/skills/magestic-figure/tests/palette_check.py). They are REASSIGNED
# here to consequence classes; the hexes are not new and the separation is not a
# fresh claim. Assignment follows severity so hue carries meaning:
#
#   loss of function   crimson family   stop_gained, start_lost   + pink frameshift
#   in-frame altering  amber / blues    missense, inframe_*, stop_lost
#   splice             green            splice_region
#   silent / noncoding neutral ramp     synonymous, intronic, intergenic, unknown
#
# The silent and non-coding tier is deliberately NEUTRAL, not a hue. `synonymous`
# is the null class in these panels and hue is spent only on classes that carry a
# claim (magestic-figure section 1, "structure is never coloured").
#
# ⚠️ KNOWN LIMITATION, recorded rather than papered over. The palette passes the
# CVD checks but NOT the greyscale check ([3b], dL* >= 15 within a legend): amber
# `missense` and sky `inframe_deletion` sit at L* 70.6 and 69.8 and are the same
# grey in a mono printout. An assignment that satisfies greyscale exists but makes
# stop_gained blue and missense purple, i.e. hue stops carrying severity. The
# correct fix is a SECOND CHANNEL, and CODON_VARIANT_MARKERS below supplies it —
# vep_panels.py currently hardcodes marker='o' for the consequence legend, so
# adopting it is the follow-up, not this change.

CODON_VARIANT_COLORS: Dict[str, str] = {
    # --- loss of function -------------------------------------------------
    'stop_gained': '#B02418',            # crimson
    'stop_gain': '#B02418',              # (alias)
    'start_lost': '#762238',             # dark maroon
    'start_loss': '#762238',             # (alias)
    'frameshift': '#CC79A7',             # pink - LoF via indel
    # --- in-frame protein-altering ---------------------------------------
    'missense': '#E69F00',               # amber
    'inframe_insertion': '#0072B2',      # blue
    'inframe_deletion': '#56B4E9',       # sky
    'inframe_indel': '#56B4E9',          # sky (alias of deletion)
    'in-frame_indel': '#56B4E9',         # (alias)
    'stop_lost': '#2400B2',              # deep blue
    'stop_loss': '#2400B2',              # (alias)
    # --- splice -----------------------------------------------------------
    'splice_region': '#009E73',          # green
    # --- silent / non-coding: NEUTRAL, lightness-ordered ------------------
    'synonymous': '#5C6672',             # rule grey - the null class
    'intronic': '#8A97A8',
    'intergenic': '#C3CAD1',
    'non_coding': '#DDE5EC',
    'non-coding': '#DDE5EC',             # (alias)
    # --- special ----------------------------------------------------------
    'complete_CDS_deletion': '#14213D',  # ink - total loss reads darkest
    'unknown': '#C3CAD1',
}


# ============================================================================
# Codon Variant Type Marker Shapes
# ============================================================================
# The non-colour channel for consequence, added 2026-09-04 alongside the palette
# fix above. Shape is what lets the panel survive a mono printout and the tritan
# advisories; gene_type already had GENE_TYPE_MARKERS and consequence did not.
#
# ⚠️ NOT YET CONSUMED. vep_panels.py hardcodes marker='o' at its consequence
# legend and scatter. Wiring it is a rendering change and is deliberately left
# out of this commit, which changes data only.

CODON_VARIANT_MARKERS: Dict[str, str] = {
    'stop_gained': 'X', 'stop_gain': 'X',
    'start_lost': 'P', 'start_loss': 'P',
    'frameshift': 'd',
    'missense': 'o',
    'inframe_insertion': '^', 'inframe_deletion': 'v',
    'inframe_indel': 'D', 'in-frame_indel': 'D',
    'stop_lost': 's', 'stop_loss': 's',
    'splice_region': '*',
    'synonymous': '.', 'intronic': '.', 'intergenic': '.',
    'non_coding': '.', 'non-coding': '.',
    'complete_CDS_deletion': 'H',
    'unknown': '.',
}


# ============================================================================
# Gene Type Colors
# ============================================================================
# Used for Yorzoi/Shorkie per-gene expression panels
# Colors match box-arrow gene track for visual consistency

# 🔴 RE-POINTED 2026-09-05 (rna-figure standard, Chanfreau planning/871): the flat-UI
# hexes were never CVD-validated. Keys unchanged; values are layer-A roles from
# rna_tokens.json -- `target` is FOCAL (the gene the panel is about IS the argument),
# the neighbours take the remaining roles so identity propagates track -> panel.
GENE_TYPE_COLORS: Dict[str, str] = {
    'target': '#D55E00',             # focal
    'upstream_1': '#0072B2',         # reference
    'upstream_2': '#56B4E9',         # stratum_b
    'downstream_1': '#853272',       # control
    'downstream_2': '#E69F00',       # stratum_a
}


# ============================================================================
# Strain Origin Colors
# ============================================================================
# Used for variant map panel to show strain distribution

# Re-pointed 2026-09-05 to validated hexes (rna_tokens.json); keys unchanged.
STRAIN_ORIGIN_COLORS: Dict[str, str] = {
    'shared_only': '#0072B2',        # reference
    'shared_plus_left': '#853272',   # control
    'shared_plus_right': '#E69F00',  # stratum_a
    'both_neighbors': '#D55E00',     # focal
    'left_only': '#009E73',          # match
    'right_only': '#56B4E9',         # stratum_b
    'all_three': '#999999',          # neutral
    'other_strains': '#C7C7C7',      # mute
}


# ============================================================================
# Control Type Colors
# ============================================================================
# Used for MAGESTIC control panels (PTC, SYN, DEL)

# Re-pointed 2026-09-05 to the SAME hexes CODON_VARIANT_COLORS gives these classes
# (stop_gained / synonymous / complete_CDS_deletion), so a control panel and a VEP
# panel agree; the old red/green pair was the CVD failure the validator plants.
CONTROL_TYPE_COLORS: Dict[str, str] = {
    'PTC': '#B02418',                # crimson - premature termination codon
    'SYN': '#5C6672',                # rule grey - synonymous (the null class)
    'DEL': '#14213D',                # ink - CDS deletion
}


# ============================================================================
# Effect Direction Colors
# ============================================================================
# Used for trio dissection and other effect-based plots

# Re-pointed 2026-09-05 (rna_tokens.json); keys unchanged.
EFFECT_DIRECTION_COLORS: Dict[str, str] = {
    'positive': '#D55E00',           # focal
    'negative': '#0072B2',           # reference
    'neutral': '#999999',            # neutral
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
    def house(cls, dpi: int = 300) -> 'FigureConfig':
        """The rna-figure house scale (rna_tokens.json typography.manuscript_pt)."""
        from . import tokens as _TOK
        T = _TOK.typography()
        return cls(dpi=dpi, title_fontsize=T['axis_label'], label_fontsize=T['axis_label'],
                   tick_fontsize=T['tick_label'], legend_fontsize=T['annotation'],
                   line_width=_TOK.stroke()['secondary'])

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
