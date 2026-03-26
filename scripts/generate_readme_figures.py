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

    # Read with soft-clip
    ax.add_patch(mpatches.Rectangle((20, 4), 15, 1.5, color=softclip_color, label='Soft-clipped'))
    ax.text(27.5, 4.75, 'soft-clip', ha='center', va='center', fontsize=9, color='white')
    ax.add_patch(mpatches.Rectangle((55, 4), 35, 1.5, color=exon_color))
    ax.text(72, 4.75, 'aligned to EXON 2', ha='center', va='center', fontsize=9, color='white')
    ax.text(2, 4.75, 'Read:', ha='right', va='center', fontsize=10)

    # Arrow showing the gap
    ax.annotate('', xy=(55, 3.5), xytext=(35, 3.5),
                arrowprops=dict(arrowstyle='<->', color='red', lw=2))
    ax.text(45, 2.8, 'Junction missed!\nAligner soft-clipped instead', ha='center', va='top',
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

    # Corrected read - now spliced through
    ax.add_patch(mpatches.Rectangle((20, 4), 15, 1.5, color=correct_color, label='Rescued'))
    ax.text(27.5, 4.75, 'EXON 1', ha='center', va='center', fontsize=9, color='white')
    ax.add_patch(mpatches.Rectangle((55, 4), 35, 1.5, color=exon_color))
    ax.text(72, 4.75, 'EXON 2', ha='center', va='center', fontsize=9, color='white')
    ax.text(2, 4.75, 'Corrected:', ha='right', va='center', fontsize=10)

    # Splice junction indicator
    ax.plot([35, 35, 55, 55], [4.75, 3, 3, 4.75], 'k-', lw=1.5)
    ax.text(45, 2.5, 'N (splice)', ha='center', va='top', fontsize=9)

    # Success indicator
    ax.text(45, 1, '✓ True 5\' end recovered for TSS analysis', ha='center', va='center',
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

    # Read with false junction
    ax.add_patch(mpatches.Rectangle((5, 4), 25, 1.5, color=exon_color))
    ax.add_patch(mpatches.Rectangle((30, 4), 30, 1.5, color=polya_color))
    ax.text(45, 4.75, 'poly(A) tail aligned', ha='center', va='center', fontsize=9, color='white')
    ax.add_patch(mpatches.Rectangle((75, 4), 15, 1.5, color=polya_color, alpha=0.5))
    ax.text(82.5, 4.75, 'soft-clip', ha='center', va='center', fontsize=9)
    ax.text(2, 4.75, 'Read:', ha='right', va='center', fontsize=10)

    # False junction indicator
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

    # Walk back arrow
    ax.annotate('', xy=(30, 5.5), xytext=(75, 5.5),
                arrowprops=dict(arrowstyle='<-', color=correct_color, lw=3))
    ax.text(52.5, 6, 'Walk back eats A\'s, discards N', ha='center', va='bottom',
            fontsize=10, color=correct_color, fontweight='bold')

    # Corrected read
    ax.add_patch(mpatches.Rectangle((5, 4), 25, 1.5, color=exon_color))
    ax.text(17.5, 4.75, 'EXON', ha='center', va='center', fontsize=9, color='white', fontweight='bold')
    ax.text(2, 4.75, 'Corrected:', ha='right', va='center', fontsize=10)

    # Correct CPA marker
    ax.annotate('', xy=(30, 3.8), xytext=(30, 3),
                arrowprops=dict(arrowstyle='->', color='green', lw=2))
    ax.text(30, 2.5, 'Correct 3\' end', ha='center', va='top', fontsize=9, color='green', fontweight='bold')

    # Poly(A) length indicator
    ax.add_patch(mpatches.FancyBboxPatch((32, 4.2), 55, 1, boxstyle="round,pad=0.05",
                                          facecolor='none', edgecolor='gray', linestyle='--'))
    ax.text(60, 3.2, 'All these A\'s counted as poly(A) tail', ha='center', va='top',
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


if __name__ == "__main__":
    create_5prime_junction_figure()
    create_false_junction_figure()
    print("\nDone! Figures saved to docs/figures/")
