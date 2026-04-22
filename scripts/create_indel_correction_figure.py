#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create aesthetic indel correction figure for README.

Shows how minimap2 alignments in A-rich regions have indel artifacts
that RECTIFY corrects, with both plus and minus strand examples.

Key features shown:
1. Indels where minimap2 forces poly(A) across non-A genomic bases
2. Sequencing errors (T instead of A) in poly(A) tail region
3. Non-A bases in genome where soft-clipping begins (aligner can't extend further)
4. Why we can't just trim poly(A) before alignment

Author: Kevin R. Roy
Date: 2026-03-19
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
import numpy as np

# Set up the figure with high DPI for quality
fig = plt.figure(figsize=(16, 14))
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
    'insertion': '#E1BEE7',  # Light purple
    'softclip': '#E3F2FD',  # Light blue
    'error': '#FFCDD2',  # Light red for sequencing error
    'arrow': '#1976D2',  # Dark blue
    'text': '#212121',  # Dark gray
    'plus': '#8BC34A',  # Light green for + strand
    'minus': '#03A9F4',  # Light blue for - strand
    'nonA': '#FFECB3',  # Light amber for non-A positions where soft-clip starts
}

def draw_base(ax, base, x, y, fontsize=13, bg_color='white', edge_color='#BDBDBD',
              width=0.85, is_deletion=False, is_insertion=False, show_error=False):
    """Draw a single base with colored background."""
    if is_deletion:
        # Deletion marker (gap in read)
        rect = FancyBboxPatch((x-width/2, y-0.35), width, 0.7,
                              boxstyle="round,pad=0.03",
                              facecolor=COLORS['deletion'],
                              edgecolor='#E65100', linewidth=1.5)
        ax.add_patch(rect)
        ax.text(x, y, '-', ha='center', va='center', fontsize=fontsize+2,
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

def draw_sequence_with_features(ax, seq, y, x_start=0, features=None, fontsize=13):
    """
    Draw a DNA sequence with specified features per position.

    features: dict of position -> dict with keys like:
        'bg': background color
        'edge': edge color
        'deletion': True if this is a deletion
        'error': True if sequencing error
    """
    if features is None:
        features = {}
    x = x_start
    for i, base in enumerate(seq):
        if base == ' ':
            x += 1
            continue
        feat = features.get(i, {})
        bg = feat.get('bg', 'white')
        edge = feat.get('edge', '#BDBDBD')
        is_del = feat.get('deletion', False)

        draw_base(ax, base, x, y, fontsize, bg, edge, is_deletion=is_del)
        x += 1
    return x

# Create subplots
gs = fig.add_gridspec(4, 1, height_ratios=[1.8, 1.8, 0.9, 0.9], hspace=0.35)
ax1 = fig.add_subplot(gs[0])  # Plus strand
ax2 = fig.add_subplot(gs[1])  # Minus strand
ax3 = fig.add_subplot(gs[2])  # Why not just trim explanation
ax4 = fig.add_subplot(gs[3])  # Summary

# ============================================================================
# Panel A: Plus Strand Example with indels and sequencing error
# ============================================================================
ax1.set_xlim(-4, 38)
ax1.set_ylim(-2.2, 5.8)
ax1.axis('off')

# Title with strand indicator
ax1.add_patch(FancyBboxPatch((-3.5, 4.8), 1.5, 0.6,
              boxstyle="round,pad=0.1", facecolor=COLORS['plus'],
              edgecolor='#558B2F', linewidth=2))
ax1.text(-2.75, 5.1, '+', ha='center', va='center', fontsize=16,
         fontweight='bold', color='white')
ax1.text(0, 5.1, 'Plus Strand: Poly(A) tail aligns to downstream genomic A-tract',
         fontsize=13, fontweight='bold', va='center')

# Reference sequence for + strand:
# Genome:  C T G A T C A G T C | A A A A A G A A A A A | G T C A G T
#          0 1 2 3 4 5 6 7 8 9  10 11 12 13 14 15 16 17 18 19 20  21 22 23 24 25 26
#                              ^poly(A) start         ^non-A       ^soft-clip boundary

ref_plus = "CTGATCAGTCAAAAAGAAAAAGTCAGT"
#           0         1         2
#           012345678901234567890123456

# Read with poly(A) tail:
# - Position 15 (G in genome) shows as deletion in read (minimap2 skips it)
# - Position 18 shows T (sequencing error in poly(A))
# - Positions 21-26 are soft-clipped (genome has GTCAGT, read has AAAAAA)
read_plus = "CTGATCAGTCAAAAA-AAATAAAAAAAA"  # - is deletion, T is seq error, trailing A's soft-clipped

# Features for reference
ref_features = {}
# A-tract region (positions 10-20) in green
for i in range(10, 21):
    ref_features[i] = {'bg': '#C8E6C9', 'edge': '#4CAF50'}
# Mark the G at position 15 specially
ref_features[15] = {'bg': '#FFCDD2', 'edge': '#C62828'}  # Red - non-A in A-tract
# Non-A positions where soft-clipping occurs (positions 21-26)
for i in range(21, 27):
    ref_features[i] = {'bg': COLORS['nonA'], 'edge': '#FF8F00'}  # Amber - this is where soft-clip starts

# Features for read
read_features = {}
# Aligned A's (before deletion)
for i in range(10, 15):
    read_features[i] = {'bg': '#C8E6C9', 'edge': '#4CAF50'}
# Deletion at position 15
read_features[15] = {'deletion': True}
# More aligned A's
for i in range(16, 18):
    read_features[i] = {'bg': '#C8E6C9', 'edge': '#4CAF50'}
# Sequencing error at position 18 (T instead of A)
read_features[18] = {'bg': COLORS['error'], 'edge': '#C62828'}
# More aligned A's (positions 19-20)
for i in range(19, 21):
    read_features[i] = {'bg': '#C8E6C9', 'edge': '#4CAF50'}
# Soft-clipped region (positions 21-27) - these A's couldn't align to genome's GTCAGT
for i in range(21, 28):
    read_features[i] = {'bg': COLORS['softclip'], 'edge': '#1565C0'}

# Draw reference
ax1.text(-3.5, 2.2, 'Genome', ha='right', va='center', fontsize=11, fontweight='bold')
draw_sequence_with_features(ax1, ref_plus, 2.2, x_start=0, features=ref_features)

# Draw aligned read
ax1.text(-3.5, 1.0, 'Read', ha='right', va='center', fontsize=11, fontweight='bold')
draw_sequence_with_features(ax1, read_plus, 1.0, x_start=0, features=read_features)

# Annotations
# Poly(A) tail start (was "True CPA site")
ax1.annotate('poly(A) tail start', xy=(10, 2.7), xytext=(5, 3.3),
            fontsize=10, ha='center', color='#C62828', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#C62828', lw=1.5))

# A-tract with interruption
ax1.annotate('', xy=(10, 2.85), xytext=(20.5, 2.85),
            arrowprops=dict(arrowstyle='<->', color='#4CAF50', lw=2))
ax1.text(15, 3.1, 'A-tract (with G interruption)', ha='center', fontsize=9, color='#2E7D32')

# Non-A region annotation (where soft-clipping occurs)
ax1.annotate('', xy=(21, 2.85), xytext=(26.5, 2.85),
            arrowprops=dict(arrowstyle='<->', color='#FF8F00', lw=2))
ax1.text(24, 3.1, 'non-A (soft-clip boundary)', ha='center', fontsize=9, color='#E65100')

# Deletion and seq error are visually distinct (orange "-" and red "A") - no labels needed

# Soft-clip annotation
ax1.annotate('', xy=(21, 0.35), xytext=(27.5, 0.35),
            arrowprops=dict(arrowstyle='<->', color='#1565C0', lw=2))
ax1.text(24, -0.1, 'Soft-clip (poly-A tail)', ha='center', fontsize=8, color='#1565C0')

# Correction arrow - just below the bases, arrow shows direction of walking back
# Walking back from last aligned A (pos 20) before soft-clip toward poly(A) start (left)
ax1.annotate('', xy=(9.5, 0.35), xytext=(20, 0.35),
            arrowprops=dict(arrowstyle='->', color=COLORS['arrow'], lw=2.5,
                           connectionstyle="arc3,rad=-0.25"))
ax1.text(15, -0.9, 'RECTIFY walks back, accounting for deletions and seq errors',
         ha='center', va='top', fontsize=9, color=COLORS['arrow'], fontweight='bold', style='italic')

# 5' and 3' labels for plus strand (5' on left, 3' on right) - positioned to not overlap bases
ax1.add_patch(FancyBboxPatch((-3, 1.5), 1.2, 0.5,
              boxstyle="round,pad=0.05", facecolor='#E0E0E0', edgecolor='#757575'))
ax1.text(-2.4, 1.75, "5'", ha='center', va='center', fontsize=11, fontweight='bold', color='#424242')
ax1.add_patch(FancyBboxPatch((28, 1.5), 1.2, 0.5,
              boxstyle="round,pad=0.05", facecolor='#E0E0E0', edgecolor='#757575'))
ax1.text(28.6, 1.75, "3'", ha='center', va='center', fontsize=11, fontweight='bold', color='#424242')

# ============================================================================
# Panel B: Minus Strand Example with indels and sequencing error
# ============================================================================
ax2.set_xlim(-4, 38)
ax2.set_ylim(-2.2, 5.8)
ax2.axis('off')

# Title with strand indicator
ax2.add_patch(FancyBboxPatch((-3.5, 4.8), 1.5, 0.6,
              boxstyle="round,pad=0.1", facecolor=COLORS['minus'],
              edgecolor='#0277BD', linewidth=2))
ax2.text(-2.75, 5.1, '-', ha='center', va='center', fontsize=16,
         fontweight='bold', color='white')
ax2.text(0, 5.1, 'Minus Strand: Poly(A) tail appears as poly(T) in IGV',
         fontsize=13, fontweight='bold', va='center')


# Reference sequence for - strand:
# Genome:  G A C T G A | T T T T T T C G T T T T T T T T | G A C T G A T C
#          0 1 2 3 4 5   6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21   22 23 24 25 26 27 28 29
#          ^non-T region  ^T-tract with CG at 12-13                    ^poly(A) start
# Non-T bases (positions 0-5) are where soft-clipping occurs
# CG at positions 12-13 causes 2 bp deletion (different from + strand's single G)
# Position 17 has T in genome but A in read (seq error)

ref_minus = "GACTGATTTTTTCGTTTTTTTTGACTGATC"
#            0         1         2
#            012345678901234567890123456789

# Read - soft-clipped T's protrude LEFT past genome boundary
# One extra T at position 0 (drawn at x=-1), then positions 1-6 align with genome's GACTGA (soft-clipped)
# 2bp deletion at CG (genome pos 12-13), sequencing error further right
read_minus = "TTTTTTTTTTTTT--TTATTTTGACTGATC"
#             0         1         2          3
#             0123456789012345678901234567890
# 13 T's, then --, then TTATTTTGACTGATC
# Draw starting at x=-1: pos 0->x=-1, pos 13->x=12 (first deletion aligns with genome's C at x=12)

# Features for reference
ref_features2 = {}
# Non-T region at start (positions 0-5) - where soft-clipping occurs
for i in range(0, 6):
    ref_features2[i] = {'bg': COLORS['nonA'], 'edge': '#FF8F00'}
# T-tract region (positions 6-21) in blue
for i in range(6, 22):
    ref_features2[i] = {'bg': '#BBDEFB', 'edge': '#1976D2'}
# Mark the CG at positions 12-13 specially
ref_features2[12] = {'bg': '#FFCDD2', 'edge': '#C62828'}  # Red - C in T-tract
ref_features2[13] = {'bg': '#FFCDD2', 'edge': '#C62828'}  # Red - G in T-tract

# Features for read (drawn starting at x=-1)
read_features2 = {}
# Soft-clipped T's (positions 0-6 in string, x=-1 to x=5) - 1 protrudes left, rest over GACTGA
for i in range(0, 7):
    read_features2[i] = {'bg': COLORS['softclip'], 'edge': '#1565C0'}
# Aligned T's (positions 7-12 in string, x=6 to x=11)
for i in range(7, 13):
    read_features2[i] = {'bg': '#BBDEFB', 'edge': '#1976D2'}
# 2 bp deletion at positions 13-14 (the CG positions - the '--' in the string, x=12-13)
read_features2[13] = {'deletion': True}
read_features2[14] = {'deletion': True}
# Aligned T's (positions 15-16 in string, x=14-15)
for i in range(15, 17):
    read_features2[i] = {'bg': '#BBDEFB', 'edge': '#1976D2'}
# Sequencing error at position 17 (A instead of T, x=16)
read_features2[17] = {'bg': COLORS['error'], 'edge': '#C62828'}
# More aligned T's (positions 18-22, x=17-21)
for i in range(18, 23):
    read_features2[i] = {'bg': '#BBDEFB', 'edge': '#1976D2'}

# Draw reference
ax2.text(-3.5, 2.4, 'Genome', ha='right', va='center', fontsize=11, fontweight='bold')
draw_sequence_with_features(ax2, ref_minus, 2.4, x_start=0, features=ref_features2)

# Draw aligned read - start at x=-1 so one soft-clipped T protrudes left
ax2.text(-3.5, 1.2, 'Read', ha='right', va='center', fontsize=11, fontweight='bold')
draw_sequence_with_features(ax2, read_minus, 1.2, x_start=-1, features=read_features2)

# Annotations
# Non-T region (soft-clip boundary)
ax2.annotate('', xy=(0, 3.0), xytext=(5.5, 3.0),
            arrowprops=dict(arrowstyle='<->', color='#FF8F00', lw=2))
ax2.text(2.75, 3.25, 'non-T (soft-clip boundary)', ha='center', fontsize=9, color='#E65100')

# T-tract with interruption
ax2.annotate('', xy=(6, 3.0), xytext=(21.5, 3.0),
            arrowprops=dict(arrowstyle='<->', color='#1976D2', lw=2))
ax2.text(14, 3.25, 'T-tract (with CG interruption)', ha='center', fontsize=9, color='#0D47A1')

# Poly(A) tail start (was "True CPA site") - positioned to not overlap bases
ax2.annotate('poly(A) tail start', xy=(22, 2.9), xytext=(27, 3.6),
            fontsize=10, ha='center', color='#C62828', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#C62828', lw=1.5))

# Soft-clip annotation (from x=-1 to x=5, one T protrudes left)
ax2.annotate('', xy=(-1, 0.55), xytext=(5.5, 0.55),
            arrowprops=dict(arrowstyle='<->', color='#1565C0', lw=2))
ax2.text(2.25, 0.1, 'Soft-clip (poly-A as T)', ha='center', fontsize=8, color='#1565C0')

# Deletion and seq error are visually distinct (orange "-" and red "A") - no labels needed

# Correction arrow - just below the bases, arrow shows direction of walking back
# Walking back from last aligned T before soft-clip (x=6) toward poly(A) start (x=22)
ax2.annotate('', xy=(22, 0.35), xytext=(6, 0.35),
            arrowprops=dict(arrowstyle='->', color=COLORS['arrow'], lw=2.5,
                           connectionstyle="arc3,rad=0.25"))
ax2.text(14, -1.5, 'RECTIFY walks back, accounting for deletions and seq errors',
         ha='center', va='top', fontsize=9, color=COLORS['arrow'], fontweight='bold', style='italic')

# 5' and 3' labels for minus strand (5' on left, 3' on right - same as plus strand) - positioned to not overlap bases
ax2.add_patch(FancyBboxPatch((-3, 1.7), 1.2, 0.5,
              boxstyle="round,pad=0.05", facecolor='#E0E0E0', edgecolor='#757575'))
ax2.text(-2.4, 1.95, "5'", ha='center', va='center', fontsize=11, fontweight='bold', color='#424242')
ax2.add_patch(FancyBboxPatch((30, 1.7), 1.2, 0.5,
              boxstyle="round,pad=0.05", facecolor='#E0E0E0', edgecolor='#757575'))
ax2.text(30.6, 1.95, "3'", ha='center', va='center', fontsize=11, fontweight='bold', color='#424242')

# ============================================================================
# Panel C: Why we can't just trim poly(A) before alignment
# ============================================================================
ax3.set_xlim(0, 10)
ax3.set_ylim(0, 2.5)
ax3.axis('off')

# Create explanation box
why_box = FancyBboxPatch((0.2, 0.2), 9.6, 2.1,
                          boxstyle="round,pad=0.1",
                          facecolor='#FFFDE7',
                          edgecolor='#F57F17', linewidth=2)
ax3.add_patch(why_box)

ax3.text(5, 2.1, 'Why not just trim poly(A) before alignment?', ha='center',
         fontsize=12, fontweight='bold', color='#E65100')

ax3.text(0.5, 1.5, '1. Genomic A-tracts extend alignment past true poly(A) start',
         ha='left', fontsize=10, color=COLORS['text'])
ax3.text(0.5, 1.1, '2. Indels reveal where aligner forced poly(A) across non-A bases (G, C, T)',
         ha='left', fontsize=10, color=COLORS['text'])
ax3.text(0.5, 0.7, '3. Genome alignment reveals sequencing errors in the poly(A) tail',
         ha='left', fontsize=10, color=COLORS['text'])

ax3.text(5, 0.35, '-> The alignment itself contains information about the poly(A) start!',
         ha='center', fontsize=10, fontweight='bold', color='#2E7D32', style='italic')

# ============================================================================
# Panel D: Summary Box
# ============================================================================
ax4.set_xlim(0, 10)
ax4.set_ylim(0, 2.2)
ax4.axis('off')

# Create summary box
summary_box = FancyBboxPatch((0.5, 0.2), 9, 1.8,
                              boxstyle="round,pad=0.1",
                              facecolor='#FAFAFA',
                              edgecolor='#424242', linewidth=2)
ax4.add_patch(summary_box)

ax4.text(5, 1.85, 'Summary: Strand-Specific Correction Direction', ha='center',
         fontsize=12, fontweight='bold')

# Plus strand summary
ax4.add_patch(FancyBboxPatch((1, 1.1), 0.4, 0.4,
              boxstyle="round,pad=0.05", facecolor=COLORS['plus'],
              edgecolor='#558B2F'))
ax4.text(1.2, 1.3, '+', ha='center', va='center', fontsize=12,
         fontweight='bold', color='white')
ax4.text(1.8, 1.3, 'Plus strand: Poly(A) extends RIGHT -> Correct by shifting LEFT (-)',
         ha='left', va='center', fontsize=10)

# Minus strand summary
ax4.add_patch(FancyBboxPatch((1, 0.5), 0.4, 0.4,
              boxstyle="round,pad=0.05", facecolor=COLORS['minus'],
              edgecolor='#0277BD'))
ax4.text(1.2, 0.7, '-', ha='center', va='center', fontsize=12,
         fontweight='bold', color='white')
ax4.text(1.8, 0.7, 'Minus strand: Poly(A) extends LEFT -> Correct by shifting RIGHT (+)',
         ha='left', va='center', fontsize=10)

plt.savefig(Path(__file__).resolve().parent.parent / 'docs/figures/indel_correction.png',
            dpi=150, bbox_inches='tight', facecolor='white')
plt.close()

print("Figure saved to docs/figures/indel_correction.png")
