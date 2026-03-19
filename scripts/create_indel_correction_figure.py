#!/usr/bin/env python3
"""
Create aesthetic indel correction figure for README.

Shows how minimap2 alignments in A-rich regions have indel artifacts
that RECTIFY corrects, with both plus and minus strand examples.

Author: Kevin R. Roy
Date: 2026-03-19
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
import numpy as np

# Set up the figure with high DPI for quality
fig = plt.figure(figsize=(16, 12))
fig.patch.set_facecolor('white')

# Color scheme
COLORS = {
    'A': '#4CAF50',  # Green
    'T': '#2196F3',  # Blue
    'C': '#FFC107',  # Amber
    'G': '#F44336',  # Red
    'match': '#E8F5E9',  # Light green
    'mismatch': '#FFCDD2',  # Light red
    'deletion': '#FFE0B2',  # Light orange
    'softclip': '#E3F2FD',  # Light blue
    'arrow': '#1976D2',  # Dark blue
    'text': '#212121',  # Dark gray
    'plus': '#8BC34A',  # Light green for + strand
    'minus': '#03A9F4',  # Light blue for - strand
}

def draw_base(ax, base, x, y, fontsize=13, bg_color='white', edge_color='#BDBDBD', width=0.85):
    """Draw a single base with colored background."""
    if base == '-':
        # Deletion marker
        rect = FancyBboxPatch((x-width/2, y-0.35), width, 0.7,
                              boxstyle="round,pad=0.03",
                              facecolor=COLORS['deletion'],
                              edgecolor='#E65100', linewidth=1.5)
        ax.add_patch(rect)
        ax.text(x, y, '−', ha='center', va='center', fontsize=fontsize+2,
               fontweight='bold', color='#E65100', family='monospace')
    elif base == ' ':
        pass  # Skip spaces
    else:
        rect = FancyBboxPatch((x-width/2, y-0.35), width, 0.7,
                              boxstyle="round,pad=0.03",
                              facecolor=bg_color,
                              edgecolor=edge_color, linewidth=1)
        ax.add_patch(rect)
        ax.text(x, y, base, ha='center', va='center', fontsize=fontsize,
               fontweight='bold', color=COLORS.get(base, COLORS['text']),
               family='monospace')

def draw_sequence(ax, seq, y, x_start=0, base_colors=None, fontsize=13):
    """Draw a DNA sequence with colored bases."""
    if base_colors is None:
        base_colors = {}
    x = x_start
    for i, base in enumerate(seq):
        if base == ' ':
            x += 1
            continue
        bg = base_colors.get(i, {}).get('bg', 'white')
        edge = base_colors.get(i, {}).get('edge', '#BDBDBD')
        draw_base(ax, base, x, y, fontsize, bg, edge)
        x += 1
    return x

# Create subplots
gs = fig.add_gridspec(3, 1, height_ratios=[1.3, 1.3, 1], hspace=0.4)
ax1 = fig.add_subplot(gs[0])  # Plus strand
ax2 = fig.add_subplot(gs[1])  # Minus strand
ax3 = fig.add_subplot(gs[2])  # Summary

# ============================================================================
# Panel A: Plus Strand Example
# ============================================================================
ax1.set_xlim(-4, 36)
ax1.set_ylim(-1.5, 4)
ax1.axis('off')

# Title with strand indicator
ax1.add_patch(FancyBboxPatch((-3.5, 3.2), 1.5, 0.6,
              boxstyle="round,pad=0.1", facecolor=COLORS['plus'],
              edgecolor='#558B2F', linewidth=2))
ax1.text(-2.75, 3.5, '+', ha='center', va='center', fontsize=16,
         fontweight='bold', color='white')
ax1.text(0, 3.5, 'Plus Strand: Poly(A) tail aligns to downstream genomic A-tract',
         fontsize=13, fontweight='bold', va='center')

# Reference sequence for + strand: CPA site followed by A-tract
ref_plus = "CTGATCAGTC AAAAAAAAAA GTCA"
#           0123456789 0123456789 0123
#                     ^CPA site

# Read with poly(A) tail extending into A-tract
read_plus = "CTGATCAGTC AAAAAAAAAA AAAA"  # Poly(A) continues past genomic

# Color scheme for bases
ref_colors = {}
read_colors = {}

# Mark the A-tract region (positions 10-19)
for i in range(10, 20):
    ref_colors[i] = {'bg': '#C8E6C9', 'edge': '#4CAF50'}

# Mark poly(A) tail region (soft-clipped, positions 20+)
for i in range(20, 24):
    read_colors[i] = {'bg': COLORS['softclip'], 'edge': '#1565C0'}

# Draw reference
ax1.text(-3.5, 2, 'Genome', ha='right', va='center', fontsize=11, fontweight='bold')
draw_sequence(ax1, ref_plus.replace(' ', ''), 2, x_start=0, base_colors=ref_colors)

# Draw aligned read
ax1.text(-3.5, 0.8, 'Read', ha='right', va='center', fontsize=11, fontweight='bold')
draw_sequence(ax1, read_plus.replace(' ', ''), 0.8, x_start=0, base_colors=read_colors)

# Add annotation for CPA site
ax1.annotate('True CPA site', xy=(10, 2.5), xytext=(6, 3.1),
            fontsize=10, ha='center', color='#C62828', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#C62828', lw=1.5))

# Add annotation for A-tract
ax1.annotate('', xy=(10, 2.65), xytext=(20, 2.65),
            arrowprops=dict(arrowstyle='<->', color='#4CAF50', lw=2))
ax1.text(15, 2.9, 'Genomic A-tract', ha='center', fontsize=9, color='#2E7D32')

# Add annotation for soft-clip
ax1.annotate('', xy=(20, 0.35), xytext=(24, 0.35),
            arrowprops=dict(arrowstyle='<->', color='#1565C0', lw=2))
ax1.text(22, 0.05, 'Soft-clip\n(poly-A tail)', ha='center', fontsize=8, color='#1565C0')

# Add correction arrow
ax1.annotate('', xy=(9.5, 0.8), xytext=(20, 0.8),
            arrowprops=dict(arrowstyle='<-', color=COLORS['arrow'], lw=2,
                           connectionstyle="arc3,rad=0.3"))
ax1.text(15, 0.2, 'RECTIFY shifts LEFT\n(upstream to true CPA)', ha='center',
         fontsize=9, color=COLORS['arrow'], fontweight='bold')

# RNA direction indicator
ax1.annotate('', xy=(28, 2), xytext=(25, 2),
            arrowprops=dict(arrowstyle='->', color='#757575', lw=2))
ax1.text(30, 2, "5'→3'\n(RNA)", ha='center', va='center', fontsize=9, color='#757575')

# ============================================================================
# Panel B: Minus Strand Example
# ============================================================================
ax2.set_xlim(-4, 36)
ax2.set_ylim(-1.5, 4.5)
ax2.axis('off')

# Title with strand indicator
ax2.add_patch(FancyBboxPatch((-3.5, 3.2), 1.5, 0.6,
              boxstyle="round,pad=0.1", facecolor=COLORS['minus'],
              edgecolor='#0277BD', linewidth=2))
ax2.text(-2.75, 3.5, '−', ha='center', va='center', fontsize=16,
         fontweight='bold', color='white')
ax2.text(0, 3.5, 'Minus Strand: Poly(A) tail appears as poly(T) in IGV',
         fontsize=13, fontweight='bold', va='center')

# Add explanation box for IGV confusion
ax2.add_patch(FancyBboxPatch((20, 3.0), 13, 1.2,
              boxstyle="round,pad=0.1", facecolor='#FFF3E0',
              edgecolor='#E65100', linewidth=1.5))
ax2.text(26.5, 3.9, '⚠️ IGV shows read sequence', ha='center', fontsize=9, fontweight='bold', color='#E65100')
ax2.text(26.5, 3.4, 'Poly(A) on RNA = T on read', ha='center', fontsize=9, color='#BF360C')

# Reference sequence for - strand: T-tract (complement of A's) followed by CPA site
# On minus strand, poly(A) on RNA = T's on genome template
ref_minus = "TTTT TTTTTTTTTT GACTGATCAG"
#                ^T-tract   ^CPA site

# Read with poly(A) tail (appears as T's in reverse-complement)
read_minus = "TTTT TTTTTTTTTT GACTGATCAG"

# Color scheme for bases
ref_colors2 = {}
read_colors2 = {}

# Mark the T-tract region (positions 0-13) - these are A's on the RNA
for i in range(0, 14):
    ref_colors2[i] = {'bg': '#BBDEFB', 'edge': '#1976D2'}
    if i < 4:
        read_colors2[i] = {'bg': COLORS['softclip'], 'edge': '#1565C0'}

# Draw reference
ax2.text(-3.5, 2, 'Genome', ha='right', va='center', fontsize=11, fontweight='bold')
draw_sequence(ax2, ref_minus.replace(' ', ''), 2, x_start=0, base_colors=ref_colors2)

# Draw aligned read
ax2.text(-3.5, 0.8, 'Read', ha='right', va='center', fontsize=11, fontweight='bold')
draw_sequence(ax2, read_minus.replace(' ', ''), 0.8, x_start=0, base_colors=read_colors2)

# Add annotation for CPA site
ax2.annotate('True CPA site', xy=(14, 2.5), xytext=(18, 3.1),
            fontsize=10, ha='center', color='#C62828', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#C62828', lw=1.5))

# Add annotation for T-tract (A's on RNA)
ax2.annotate('', xy=(0, 2.65), xytext=(14, 2.65),
            arrowprops=dict(arrowstyle='<->', color='#1976D2', lw=2))
ax2.text(7, 2.9, 'T-tract (= A on RNA)', ha='center', fontsize=9, color='#0D47A1')

# Add annotation for soft-clip
ax2.annotate('', xy=(0, 0.35), xytext=(4, 0.35),
            arrowprops=dict(arrowstyle='<->', color='#1565C0', lw=2))
ax2.text(2, 0.0, 'Soft-clip\n(poly-A as T)', ha='center', fontsize=8, color='#1565C0')

# Add correction arrow
ax2.annotate('', xy=(14.5, 0.8), xytext=(4, 0.8),
            arrowprops=dict(arrowstyle='<-', color=COLORS['arrow'], lw=2,
                           connectionstyle="arc3,rad=-0.3"))
ax2.text(9, 0.2, 'RECTIFY shifts RIGHT\n(upstream on − strand)', ha='center',
         fontsize=9, color=COLORS['arrow'], fontweight='bold')

# RNA direction indicator (opposite for minus strand)
ax2.annotate('', xy=(28, 2), xytext=(31, 2),
            arrowprops=dict(arrowstyle='->', color='#757575', lw=2))
ax2.text(25, 2, "3'←5'\n(RNA)", ha='center', va='center', fontsize=9, color='#757575')

# ============================================================================
# Panel C: Summary Box
# ============================================================================
ax3.set_xlim(0, 10)
ax3.set_ylim(0, 3)
ax3.axis('off')

# Create summary box
summary_box = FancyBboxPatch((0.5, 0.3), 9, 2.4,
                              boxstyle="round,pad=0.1",
                              facecolor='#FAFAFA',
                              edgecolor='#424242', linewidth=2)
ax3.add_patch(summary_box)

ax3.text(5, 2.5, 'Summary: Strand-Specific Correction Direction', ha='center',
         fontsize=13, fontweight='bold')

# Plus strand summary
ax3.add_patch(FancyBboxPatch((1, 1.4), 0.4, 0.4,
              boxstyle="round,pad=0.05", facecolor=COLORS['plus'],
              edgecolor='#558B2F'))
ax3.text(1.2, 1.6, '+', ha='center', va='center', fontsize=12,
         fontweight='bold', color='white')
ax3.text(1.8, 1.6, 'Plus strand: Poly(A) extends RIGHT → Correct by shifting LEFT (−)',
         ha='left', va='center', fontsize=11)

# Minus strand summary
ax3.add_patch(FancyBboxPatch((1, 0.7), 0.4, 0.4,
              boxstyle="round,pad=0.05", facecolor=COLORS['minus'],
              edgecolor='#0277BD'))
ax3.text(1.2, 0.9, '−', ha='center', va='center', fontsize=12,
         fontweight='bold', color='white')
ax3.text(1.8, 0.9, 'Minus strand: Poly(A) extends LEFT → Correct by shifting RIGHT (+)',
         ha='left', va='center', fontsize=11)

plt.savefig('/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/docs/figures/indel_correction.png',
            dpi=150, bbox_inches='tight', facecolor='white')
plt.close()

print("Figure saved to docs/figures/indel_correction.png")
