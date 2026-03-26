#!/usr/bin/env python3
"""
Generate README figures for RECTIFY.

Creates diagrams for:
1. 5' junction soft-clip rescue
2. 3' false junction (walk back eating through A-tract)

Author: Kevin R. Roy
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Output directory
OUTPUT_DIR = Path(__file__).parent.parent / "docs" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def create_5prime_junction_figure():
    """Create figure showing 5' junction soft-clip rescue."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 6))

    # Colors
    exon_color = '#4CAF50'  # Green
    intron_color = '#E0E0E0'  # Light gray
    softclip_color = '#FF9800'  # Orange
    correct_color = '#2196F3'  # Blue

    # Top panel: The Problem
    ax = axes[0]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 10)
    ax.set_title("THE PROBLEM: 5' soft-clips hide true splice junction", fontsize=12, fontweight='bold')

    # Genome with splice junction
    ax.add_patch(mpatches.Rectangle((5, 7), 30, 1.5, color=exon_color, label='Exon'))
    ax.text(20, 7.75, 'EXON 1', ha='center', va='center', fontsize=10, color='white', fontweight='bold')
    ax.add_patch(mpatches.Rectangle((35, 7), 20, 1.5, color=intron_color))
    ax.text(45, 7.75, 'intron', ha='center', va='center', fontsize=10, style='italic')
    ax.add_patch(mpatches.Rectangle((55, 7), 35, 1.5, color=exon_color))
    ax.text(72, 7.75, 'EXON 2', ha='center', va='center', fontsize=10, color='white', fontweight='bold')
    ax.text(35, 6.2, 'GT', ha='center', va='top', fontsize=9, color='#666')
    ax.text(55, 6.2, 'AG', ha='center', va='top', fontsize=9, color='#666')
    ax.text(2, 7.75, 'Genome:', ha='right', va='center', fontsize=10)

    # Read with soft-clip - soft-clip is ADJACENT to aligned portion (same read!)
    # The soft-clip is at the 5' end (left side of read)
    ax.add_patch(mpatches.Rectangle((40, 4), 15, 1.5, color=softclip_color, label='Soft-clipped'))
    ax.text(47.5, 4.75, 'soft-clip', ha='center', va='center', fontsize=9, color='white')
    ax.add_patch(mpatches.Rectangle((55, 4), 35, 1.5, color=exon_color))
    ax.text(72, 4.75, 'aligned to EXON 2', ha='center', va='center', fontsize=9, color='white')
    ax.text(2, 4.75, 'Read:', ha='right', va='center', fontsize=10)

    # Bracket showing the soft-clip matches EXON 1
    ax.annotate('', xy=(20, 3.2), xytext=(47.5, 3.2),
                arrowprops=dict(arrowstyle='<->', color='red', lw=2))
    ax.text(33, 2.5, 'These bases match EXON 1\nbut aligner soft-clipped them!', ha='center', va='top',
            fontsize=9, color='red')

    ax.axis('off')

    # Bottom panel: The Solution
    ax = axes[1]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 10)
    ax.set_title("RECTIFY'S SOLUTION: Rescue soft-clipped bases through splice junction",
                 fontsize=12, fontweight='bold')

    # Genome (same)
    ax.add_patch(mpatches.Rectangle((5, 7), 30, 1.5, color=exon_color))
    ax.text(20, 7.75, 'EXON 1', ha='center', va='center', fontsize=10, color='white', fontweight='bold')
    ax.add_patch(mpatches.Rectangle((35, 7), 20, 1.5, color=intron_color))
    ax.text(45, 7.75, 'intron', ha='center', va='center', fontsize=10, style='italic')
    ax.add_patch(mpatches.Rectangle((55, 7), 35, 1.5, color=exon_color))
    ax.text(72, 7.75, 'EXON 2', ha='center', va='center', fontsize=10, color='white', fontweight='bold')
    ax.text(35, 6.2, 'GT', ha='center', va='top', fontsize=9, color='#666')
    ax.text(55, 6.2, 'AG', ha='center', va='top', fontsize=9, color='#666')
    ax.text(2, 7.75, 'Genome:', ha='right', va='center', fontsize=10)

    # Corrected read - now aligned through splice junction
    ax.add_patch(mpatches.Rectangle((20, 4), 15, 1.5, color=correct_color, label='Rescued'))
    ax.text(27.5, 4.75, 'EXON 1', ha='center', va='center', fontsize=9, color='white')
    ax.add_patch(mpatches.Rectangle((55, 4), 35, 1.5, color=exon_color))
    ax.text(72, 4.75, 'EXON 2', ha='center', va='center', fontsize=9, color='white')
    ax.text(2, 4.75, 'Corrected:', ha='right', va='center', fontsize=10)

    # Splice junction indicator (N operation in CIGAR)
    ax.plot([35, 35, 55, 55], [4.75, 3, 3, 4.75], 'k-', lw=1.5)
    ax.text(45, 2.5, 'N (splice)', ha='center', va='top', fontsize=9)

    # Success indicator
    ax.text(50, 1, '✓ True read 5\' end recovered • Splice junction rescued', ha='center', va='center',
            fontsize=11, color='green', fontweight='bold')

    ax.axis('off')

    # Legend
    legend_elements = [
        mpatches.Patch(color=exon_color, label='Exon (aligned)'),
        mpatches.Patch(color=intron_color, label='Intron'),
        mpatches.Patch(color=softclip_color, label='Soft-clipped (problem)'),
        mpatches.Patch(color=correct_color, label='Rescued (solution)'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=9)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12)

    output_path = OUTPUT_DIR / "5prime_junction_rescue.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Created: {output_path}")


def create_false_junction_figure():
    """Create figure showing walk back eating through false 3' junction."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 6))

    # Colors
    exon_color = '#4CAF50'  # Green
    atract_color = '#FFC107'  # Amber for A-tract
    polya_color = '#FF5722'  # Deep orange for poly(A)
    false_junc_color = '#F44336'  # Red for false junction
    correct_color = '#2196F3'  # Blue

    # Top panel: The Problem
    ax = axes[0]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 10)
    ax.set_title("THE PROBLEM: Poly(A) tail creates false 3' junction", fontsize=12, fontweight='bold')

    # Genome
    ax.add_patch(mpatches.Rectangle((5, 7), 25, 1.5, color=exon_color))
    ax.text(17.5, 7.75, 'EXON', ha='center', va='center', fontsize=10, color='white', fontweight='bold')
    ax.add_patch(mpatches.Rectangle((30, 7), 15, 1.5, color=atract_color))
    ax.text(37.5, 7.75, 'AAAAA', ha='center', va='center', fontsize=10)
    ax.add_patch(mpatches.Rectangle((45, 7), 15, 1.5, color='#E0E0E0'))
    ax.text(52.5, 7.75, 'gap', ha='center', va='center', fontsize=10, style='italic')
    ax.add_patch(mpatches.Rectangle((60, 7), 15, 1.5, color=atract_color))
    ax.text(67.5, 7.75, 'AAAAA', ha='center', va='center', fontsize=10)
    ax.text(2, 7.75, 'Genome:', ha='right', va='center', fontsize=10)

    # True CPA marker
    ax.annotate('', xy=(30, 9.2), xytext=(30, 8.7),
                arrowprops=dict(arrowstyle='->', color='green', lw=2))
    ax.text(30, 9.5, 'True CPA', ha='center', va='bottom', fontsize=9, color='green', fontweight='bold')

    # Read with false junction - soft-clip is ADJACENT to aligned portion
    ax.add_patch(mpatches.Rectangle((5, 4), 25, 1.5, color=exon_color))
    ax.add_patch(mpatches.Rectangle((30, 4), 35, 1.5, color=polya_color))
    ax.text(47.5, 4.75, 'poly(A) tail aligned', ha='center', va='center', fontsize=9, color='white')
    # Soft-clip directly adjacent to aligned portion
    ax.add_patch(mpatches.Rectangle((65, 4), 15, 1.5, color=polya_color, alpha=0.5))
    ax.text(72.5, 4.75, 'soft-clip', ha='center', va='center', fontsize=9)
    ax.text(2, 4.75, 'Read:', ha='right', va='center', fontsize=10)

    # False junction indicator (in the middle of the poly(A) alignment)
    ax.add_patch(mpatches.Rectangle((44, 3.8), 2, 2, color=false_junc_color))
    ax.annotate('', xy=(45, 2.8), xytext=(45, 3.6),
                arrowprops=dict(arrowstyle='->', color='red', lw=2))
    ax.text(45, 2.3, 'FALSE JUNCTION\n(N operation)', ha='center', va='top',
            fontsize=9, color='red', fontweight='bold')

    ax.axis('off')

    # Bottom panel: The Solution
    ax = axes[1]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 10)
    ax.set_title("RECTIFY'S SOLUTION: Walk back eats through false junction",
                 fontsize=12, fontweight='bold')

    # Genome (same)
    ax.add_patch(mpatches.Rectangle((5, 7), 25, 1.5, color=exon_color))
    ax.text(17.5, 7.75, 'EXON', ha='center', va='center', fontsize=10, color='white', fontweight='bold')
    ax.add_patch(mpatches.Rectangle((30, 7), 15, 1.5, color=atract_color))
    ax.text(37.5, 7.75, 'AAAAA', ha='center', va='center', fontsize=10)
    ax.add_patch(mpatches.Rectangle((45, 7), 15, 1.5, color='#E0E0E0'))
    ax.text(52.5, 7.75, 'gap', ha='center', va='center', fontsize=10, style='italic')
    ax.add_patch(mpatches.Rectangle((60, 7), 15, 1.5, color=atract_color))
    ax.text(67.5, 7.75, 'AAAAA', ha='center', va='center', fontsize=10)
    ax.text(2, 7.75, 'Genome:', ha='right', va='center', fontsize=10)

    # Walk back arrow - points LEFT (walking back from 3' end toward 5' end)
    ax.annotate('', xy=(30, 5.5), xytext=(80, 5.5),
                arrowprops=dict(arrowstyle='->', color=correct_color, lw=3))
    ax.text(55, 6, 'Walk back eats A\'s, discards N', ha='center', va='bottom',
            fontsize=10, color=correct_color, fontweight='bold')

    # Corrected read - just the exon part
    ax.add_patch(mpatches.Rectangle((5, 4), 25, 1.5, color=exon_color))
    ax.text(17.5, 4.75, 'EXON', ha='center', va='center', fontsize=9, color='white', fontweight='bold')
    ax.text(2, 4.75, 'Corrected:', ha='right', va='center', fontsize=10)

    # Correct CPA marker
    ax.annotate('', xy=(30, 3.8), xytext=(30, 3),
                arrowprops=dict(arrowstyle='->', color='green', lw=2))
    ax.text(30, 2.5, 'Correct 3\' end', ha='center', va='top', fontsize=9, color='green', fontweight='bold')

    # Poly(A) length indicator
    ax.add_patch(mpatches.FancyBboxPatch((32, 4.2), 48, 1, boxstyle="round,pad=0.05",
                                          facecolor='none', edgecolor='gray', linestyle='--'))
    ax.text(56, 3.2, 'All these A\'s counted as poly(A) tail', ha='center', va='top',
            fontsize=9, color='gray')

    # Success message
    ax.text(50, 0.8, '✓ False junction ignored • True CPA found • Correct poly(A) length',
            ha='center', va='center', fontsize=11, color='green', fontweight='bold')

    ax.axis('off')

    # Legend
    legend_elements = [
        mpatches.Patch(color=exon_color, label='Exon'),
        mpatches.Patch(color=atract_color, label='Genomic A-tract'),
        mpatches.Patch(color=polya_color, label='Poly(A) tail'),
        mpatches.Patch(color=false_junc_color, label='False junction'),
        mpatches.Patch(color=correct_color, label='Walk back direction'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=5, fontsize=9)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12)

    output_path = OUTPUT_DIR / "false_junction_walkback.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Created: {output_path}")


def create_multi_aligner_consensus_figure():
    """Create figure showing multi-aligner consensus pipeline."""
    fig, axes = plt.subplots(3, 1, figsize=(12, 9))

    # Colors
    minimap2_color = '#2196F3'  # Blue
    mappacbio_color = '#4CAF50'  # Green
    gapmm2_color = '#FF9800'  # Orange
    consensus_color = '#9C27B0'  # Purple
    exon_color = '#E0E0E0'  # Light gray
    softclip_color = '#F44336'  # Red
    rescued_color = '#00BCD4'  # Cyan for rescued

    # Panel 1: The Problem - Different aligners give different results
    ax = axes[0]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 10)
    ax.set_title("STEP 1: Run multiple aligners (each may soft-clip differently)",
                 fontsize=11, fontweight='bold')

    # Genome reference
    ax.add_patch(mpatches.Rectangle((5, 8), 25, 1, color='#4CAF50'))
    ax.text(17.5, 8.5, 'EXON 1', ha='center', va='center', fontsize=8, color='white', fontweight='bold')
    ax.add_patch(mpatches.Rectangle((30, 8), 15, 1, color=exon_color))
    ax.text(37.5, 8.5, 'intron', ha='center', va='center', fontsize=8, style='italic')
    ax.add_patch(mpatches.Rectangle((45, 8), 45, 1, color='#4CAF50'))
    ax.text(67.5, 8.5, 'EXON 2', ha='center', va='center', fontsize=8, color='white', fontweight='bold')
    ax.text(2, 8.5, 'Genome:', ha='right', va='center', fontsize=9)
    ax.text(30, 7.5, 'GT', ha='center', va='top', fontsize=7, color='#666')
    ax.text(45, 7.5, 'AG', ha='center', va='top', fontsize=7, color='#666')

    # minimap2 - soft-clips at 5' end
    ax.add_patch(mpatches.Rectangle((35, 5.5), 10, 0.8, color=softclip_color))
    ax.text(40, 5.9, 'S', ha='center', va='center', fontsize=8, color='white')
    ax.add_patch(mpatches.Rectangle((45, 5.5), 45, 0.8, color=minimap2_color))
    ax.text(2, 5.9, 'minimap2:', ha='right', va='center', fontsize=9, color=minimap2_color)
    ax.text(92, 5.9, '5\' soft-clipped', ha='left', va='center', fontsize=7, color=softclip_color)

    # mapPacBio - finds junction
    ax.add_patch(mpatches.Rectangle((20, 3.8), 10, 0.8, color=mappacbio_color))
    ax.plot([30, 30, 45, 45], [4.2, 3.2, 3.2, 4.2], 'k-', lw=1)
    ax.add_patch(mpatches.Rectangle((45, 3.8), 45, 0.8, color=mappacbio_color))
    ax.text(2, 4.2, 'mapPacBio:', ha='right', va='center', fontsize=9, color=mappacbio_color)
    ax.text(92, 4.2, 'spliced', ha='left', va='center', fontsize=7, color='green')

    # gapmm2 - also soft-clips
    ax.add_patch(mpatches.Rectangle((38, 2.1), 7, 0.8, color=softclip_color))
    ax.text(41.5, 2.5, 'S', ha='center', va='center', fontsize=8, color='white')
    ax.add_patch(mpatches.Rectangle((45, 2.1), 45, 0.8, color=gapmm2_color))
    ax.text(2, 2.5, 'gapmm2:', ha='right', va='center', fontsize=9, color=gapmm2_color)
    ax.text(92, 2.5, '5\' soft-clipped', ha='left', va='center', fontsize=7, color=softclip_color)

    ax.axis('off')

    # Panel 2: Rescue step - try to extend soft-clips through junctions
    ax = axes[1]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 10)
    ax.set_title("STEP 2: Attempt to rescue 5' soft-clips through known junctions",
                 fontsize=11, fontweight='bold')

    # Genome reference (same)
    ax.add_patch(mpatches.Rectangle((5, 8), 25, 1, color='#4CAF50'))
    ax.text(17.5, 8.5, 'EXON 1', ha='center', va='center', fontsize=8, color='white', fontweight='bold')
    ax.add_patch(mpatches.Rectangle((30, 8), 15, 1, color=exon_color))
    ax.text(37.5, 8.5, 'intron', ha='center', va='center', fontsize=8, style='italic')
    ax.add_patch(mpatches.Rectangle((45, 8), 45, 1, color='#4CAF50'))
    ax.text(67.5, 8.5, 'EXON 2', ha='center', va='center', fontsize=8, color='white', fontweight='bold')
    ax.text(2, 8.5, 'Genome:', ha='right', va='center', fontsize=9)

    # minimap2 RESCUED - soft-clip extended through junction
    ax.add_patch(mpatches.Rectangle((20, 5.5), 10, 0.8, color=rescued_color))
    ax.text(25, 5.9, 'rescued', ha='center', va='center', fontsize=7, color='white')
    ax.plot([30, 30, 45, 45], [5.9, 4.9, 4.9, 5.9], 'k-', lw=1)
    ax.add_patch(mpatches.Rectangle((45, 5.5), 45, 0.8, color=minimap2_color))
    ax.text(2, 5.9, 'minimap2:', ha='right', va='center', fontsize=9, color=minimap2_color)
    ax.text(92, 5.9, '✓ rescued', ha='left', va='center', fontsize=7, color='green')

    # mapPacBio - already spliced, no change
    ax.add_patch(mpatches.Rectangle((20, 3.8), 10, 0.8, color=mappacbio_color))
    ax.plot([30, 30, 45, 45], [4.2, 3.2, 3.2, 4.2], 'k-', lw=1)
    ax.add_patch(mpatches.Rectangle((45, 3.8), 45, 0.8, color=mappacbio_color))
    ax.text(2, 4.2, 'mapPacBio:', ha='right', va='center', fontsize=9, color=mappacbio_color)
    ax.text(92, 4.2, '(already spliced)', ha='left', va='center', fontsize=7, color='#666')

    # gapmm2 RESCUED
    ax.add_patch(mpatches.Rectangle((22, 2.1), 8, 0.8, color=rescued_color))
    ax.text(26, 2.5, 'rescued', ha='center', va='center', fontsize=7, color='white')
    ax.plot([30, 30, 45, 45], [2.5, 1.5, 1.5, 2.5], 'k-', lw=1)
    ax.add_patch(mpatches.Rectangle((45, 2.1), 45, 0.8, color=gapmm2_color))
    ax.text(2, 2.5, 'gapmm2:', ha='right', va='center', fontsize=9, color=gapmm2_color)
    ax.text(92, 2.5, '✓ rescued', ha='left', va='center', fontsize=7, color='green')

    ax.axis('off')

    # Panel 3: Score and select best
    ax = axes[2]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 10)
    ax.set_title("STEP 3: Score rescued alignments and select best",
                 fontsize=11, fontweight='bold')

    # Scoring criteria box
    box_x, box_y = 5, 5.5
    ax.add_patch(mpatches.FancyBboxPatch((box_x, box_y), 38, 3.5, boxstyle="round,pad=0.1",
                                          facecolor='#F5F5F5', edgecolor='#333', linewidth=1.5))
    ax.text(box_x + 19, box_y + 2.8, 'Scoring Criteria:', ha='center', va='center',
            fontsize=9, fontweight='bold')
    ax.text(box_x + 2, box_y + 1.8, '• Canonical junctions (GT-AG): +10',
            ha='left', va='center', fontsize=8)
    ax.text(box_x + 2, box_y + 0.9, '• Annotated junctions: +5',
            ha='left', va='center', fontsize=8)
    ax.text(box_x + 2, box_y + 0.0, '• Remaining soft-clip: -2/base',
            ha='left', va='center', fontsize=8)

    # Arrow
    ax.annotate('', xy=(53, 7), xytext=(45, 7),
                arrowprops=dict(arrowstyle='->', color='#333', lw=2))

    # Result box
    ax.add_patch(mpatches.FancyBboxPatch((55, 5.5), 40, 3.5, boxstyle="round,pad=0.1",
                                          facecolor='#E8F5E9', edgecolor=consensus_color, linewidth=2))
    ax.text(75, 8.2, 'BEST ALIGNMENT', ha='center', va='center',
            fontsize=10, fontweight='bold', color=consensus_color)
    ax.text(75, 7, 'Highest score selected', ha='center', va='center', fontsize=8)
    ax.text(75, 6.1, 'Written to consensus BAM', ha='center', va='center', fontsize=8, color='#666')

    # Confidence scoring at bottom
    ax.add_patch(mpatches.Rectangle((10, 1), 80, 3, facecolor='white', edgecolor='#DDD', linewidth=1))
    ax.text(50, 3.4, 'Confidence (junction agreement):', ha='center', va='center', fontsize=9, fontweight='bold')
    ax.text(25, 2.2, 'HIGH', ha='center', va='center', fontsize=8, fontweight='bold', color='green')
    ax.text(25, 1.5, '3/3 agree', ha='center', va='center', fontsize=7)
    ax.text(50, 2.2, 'MEDIUM', ha='center', va='center', fontsize=8, fontweight='bold', color='#FF9800')
    ax.text(50, 1.5, '2/3 agree', ha='center', va='center', fontsize=7)
    ax.text(75, 2.2, 'LOW', ha='center', va='center', fontsize=8, fontweight='bold', color='#F44336')
    ax.text(75, 1.5, '1/3 only', ha='center', va='center', fontsize=7)

    # Dividers
    ax.plot([37.5, 37.5], [1.3, 3.7], 'k-', lw=0.5, alpha=0.3)
    ax.plot([62.5, 62.5], [1.3, 3.7], 'k-', lw=0.5, alpha=0.3)

    ax.axis('off')

    # Legend
    legend_elements = [
        mpatches.Patch(color=minimap2_color, label='minimap2'),
        mpatches.Patch(color=mappacbio_color, label='mapPacBio'),
        mpatches.Patch(color=gapmm2_color, label='gapmm2'),
        mpatches.Patch(color=softclip_color, label='Soft-clipped'),
        mpatches.Patch(color=rescued_color, label='Rescued'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=5, fontsize=9)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.08)

    output_path = OUTPUT_DIR / "multi_aligner_consensus.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Created: {output_path}")


if __name__ == "__main__":
    create_5prime_junction_figure()
    create_false_junction_figure()
    create_multi_aligner_consensus_figure()
    print("\nDone! Figures saved to docs/figures/")
