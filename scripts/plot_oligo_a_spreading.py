#!/usr/bin/env python3
"""
Generate visualization of oligo(A) spreading artifact in NET-seq data.

Shows how right-tailed distributions arise from oligo-A tail alignment
to downstream genomic A-tracts, with overlapping tails and soft-clip spillover.

Author: Kevin R. Roy
Date: 2026-03-17
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from scipy.stats import expon
from pathlib import Path

# Output directory
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "docs/figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def generate_right_tailed_counts(positions, peak_pos, amplitude, tail_length=4):
    """
    Generate discrete right-tailed read counts (no left tail).

    Signal is concentrated at peak position with exponential decay to the right only.
    """
    counts = np.zeros(len(positions), dtype=int)

    for i, pos in enumerate(positions):
        if pos == peak_pos:
            # ~54% of signal at true position
            counts[i] = int(amplitude * 0.54)
        elif pos > peak_pos:
            # Exponential decay to the right only
            distance = pos - peak_pos
            if distance <= tail_length * 3:  # Cutoff at 3x tail length
                counts[i] = int(amplitude * 0.46 * np.exp(-distance / tail_length))

    return counts


def create_spreading_visualization():
    """Create visualization with discrete bar counts showing right-tailed spreading."""

    fig, ax = plt.subplots(figsize=(12, 7))

    # Discrete genomic positions
    positions = np.arange(0, 80)

    # A-tract window: 12 bp from position 50-62, boundary at 62
    atract_start = 50
    atract_end = 62

    # True CPA positions - all within 12bp A-tract, blue close to green (4bp upstream)
    # Red at start, Blue 4bp before Green, Green near end
    true_peaks = [
        {"pos": 51, "amp": 28, "color": "#E74C3C", "label": "CPA 1"},
        {"pos": 55, "amp": 40, "color": "#3498DB", "label": "CPA 2"},  # 4bp before green
        {"pos": 59, "amp": 35, "color": "#27AE60", "label": "CPA 3"},  # Near A-tract end
    ]

    # Generate counts for each peak with smooth right-tailed spreading
    total_counts = np.zeros(len(positions), dtype=float)
    individual_counts = []

    for peak in true_peaks:
        counts = np.zeros(len(positions), dtype=float)
        for i, pos in enumerate(positions):
            if pos == peak["pos"]:
                # Main peak (~54% of signal)
                counts[i] = peak["amp"] * 0.54
            elif pos > peak["pos"]:
                # Smooth exponential decay to the right (continuous, no gaps)
                distance = pos - peak["pos"]
                # Use continuous exponential - signal continues beyond A-tract
                counts[i] = peak["amp"] * 0.46 * np.exp(-distance / 3.5)
        individual_counts.append(counts)
        total_counts += counts

    # Convert to integers but keep fractional for smooth appearance
    # Round to nearest integer for display
    display_counts = np.round(total_counts).astype(int)

    # Add minimal background noise outside the main signal area
    np.random.seed(42)
    for i in range(len(positions)):
        if display_counts[i] == 0 and (positions[i] < atract_start - 5 or positions[i] > atract_end + 15):
            display_counts[i] = np.random.choice([0, 0, 0, 1, 1, 2])

    # Ensure continuous signal from green peak through soft-clip region (no zeros)
    green_start = 59
    for i in range(green_start, min(len(positions), green_start + 18)):
        if display_counts[i] == 0:
            # Fill gaps with decaying signal
            distance = i - green_start
            fill_val = max(1, int(35 * 0.46 * np.exp(-distance / 3.5)))
            display_counts[i] = fill_val

    # Draw A-tract region background
    ax.axvspan(atract_start, atract_end, alpha=0.2, color='gold', label='Genomic A-tract (12bp)')
    ax.axvline(atract_end, color='orange', linestyle='-', linewidth=3, alpha=0.9)

    # Plot bars - color by source peak
    bar_colors = []
    for i, pos in enumerate(positions):
        if pos > atract_end:
            bar_colors.append('#CC0000')  # Red for soft-clipped region
        elif display_counts[i] <= 1 and pos < atract_start:
            bar_colors.append('#888888')  # Gray for noise outside A-tract
        else:
            # Find dominant peak contribution
            max_contrib = 0
            dom_color = '#888888'
            for peak, counts in zip(true_peaks, individual_counts):
                if counts[i] > max_contrib:
                    max_contrib = counts[i]
                    dom_color = peak["color"]
            bar_colors.append(dom_color)

    ax.bar(positions, display_counts, width=0.8, color=bar_colors, edgecolor='black', linewidth=0.3)

    # Mark true CPA positions with vertical dashed lines
    for peak in true_peaks:
        ax.axvline(peak["pos"], color=peak["color"], linestyle='--', linewidth=2.5, alpha=0.9,
                   label=f'{peak["label"]} (pos {peak["pos"]})')

    # Label A-tract boundary
    ax.annotate('A-tract\nboundary', xy=(atract_end, max(display_counts) * 0.9), fontsize=14,
                ha='center', color='darkorange', fontweight='bold')

    # Add annotation showing blue shedding into green zone
    ax.annotate('Blue signal\nspills into\ngreen zone',
                xy=(57, 8), xytext=(44, 20),
                fontsize=12, ha='center',
                arrowprops=dict(arrowstyle='->', color='#3498DB', lw=2),
                bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow', edgecolor='#3498DB'))

    # Add annotation for soft-clip region
    ax.annotate('Oligo-A tail from\ngreen peak soft-clipped\nat A-tract end',
                xy=(68, 5), xytext=(72, 18),
                fontsize=12, ha='center',
                arrowprops=dict(arrowstyle='->', color='darkred', lw=2),
                bbox=dict(boxstyle='round,pad=0.4', facecolor='mistyrose', edgecolor='darkred'))

    # Styling - larger text throughout
    ax.set_xlabel('Genomic Position (bp)', fontsize=16, fontweight='bold')
    ax.set_ylabel('Read Counts', fontsize=16, fontweight='bold')
    ax.set_title('Oligo(A) Spreading Artifact in NET-seq', fontsize=18, fontweight='bold')
    ax.tick_params(axis='both', labelsize=14)

    # Legend - larger and more prominent
    ax.legend(loc='upper left', fontsize=14, framealpha=0.95)

    # Set limits - focus on the A-tract region
    ax.set_xlim(35, 80)
    ax.set_ylim(0, max(display_counts) * 1.25)  # Extra space for A-tract annotation

    # Grid
    ax.grid(True, alpha=0.3, linestyle='-', axis='y')
    ax.set_axisbelow(True)

    plt.tight_layout()

    # Save
    output_path = OUTPUT_DIR / "oligo_a_spreading.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")

    alt_path = Path(__file__).resolve().parent.parent / "oligo_a_spreading.png"
    plt.savefig(alt_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Saved: {alt_path}")

    plt.close()

    return output_path


def create_deconvolution_visualization():
    """Show before/after deconvolution comparison - tails "added back" to true peaks."""

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    positions = np.arange(0, 80)

    # A-tract window matching the spreading figure
    atract_start = 50
    atract_end = 62

    # Same peaks as spreading figure - all within 12bp A-tract
    true_peaks = [
        {"pos": 51, "amp": 28, "color": "#E74C3C", "label": "CPA 1"},
        {"pos": 55, "amp": 40, "color": "#3498DB", "label": "CPA 2"},
        {"pos": 59, "amp": 35, "color": "#27AE60", "label": "CPA 3"},
    ]

    np.random.seed(42)

    # Panel A: Observed signal with spreading (same as spreading figure)
    ax1 = axes[0]

    # Generate counts with right-tailed spreading
    total_observed = np.zeros(len(positions), dtype=float)
    individual_counts = []

    for peak in true_peaks:
        counts = np.zeros(len(positions), dtype=float)
        for i, pos in enumerate(positions):
            if pos == peak["pos"]:
                counts[i] = peak["amp"] * 0.54
            elif pos > peak["pos"]:
                distance = pos - peak["pos"]
                counts[i] = peak["amp"] * 0.46 * np.exp(-distance / 3.5)
        individual_counts.append(counts)
        total_observed += counts

    display_observed = np.round(total_observed).astype(int)

    # Ensure continuous signal from green peak through soft-clip region
    green_start = 59
    for i in range(green_start, min(len(positions), green_start + 18)):
        if display_observed[i] == 0:
            distance = i - green_start
            fill_val = max(1, int(35 * 0.46 * np.exp(-distance / 3.5)))
            display_observed[i] = fill_val

    # Color bars by dominant peak, red for soft-clipped
    bar_colors_observed = []
    for i, pos in enumerate(positions):
        if pos > atract_end:
            bar_colors_observed.append('#CC0000')
        else:
            max_contrib = 0
            dom_color = '#888888'
            for peak, counts in zip(true_peaks, individual_counts):
                if counts[i] > max_contrib:
                    max_contrib = counts[i]
                    dom_color = peak["color"]
            bar_colors_observed.append(dom_color)

    # Draw A-tract and bars
    ax1.axvspan(atract_start, atract_end, alpha=0.2, color='gold', label='Genomic A-tract')
    ax1.axvline(atract_end, color='orange', linestyle='-', linewidth=3, alpha=0.9)
    ax1.bar(positions, display_observed, width=0.8, color=bar_colors_observed, edgecolor='black', linewidth=0.3)

    # Mark true CPA positions
    for peak in true_peaks:
        ax1.axvline(peak["pos"], color=peak["color"], linestyle='--', linewidth=2.5, alpha=0.9)

    ax1.set_title('A. Observed Signal (with oligo-A spreading)', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Genomic Position (bp)', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Read Counts', fontsize=14, fontweight='bold')
    ax1.tick_params(axis='both', labelsize=12)

    ax1.annotate('Tails spread\ndownstream',
                xy=(65, 10), xytext=(70, 22),
                fontsize=11, ha='center',
                arrowprops=dict(arrowstyle='->', color='darkred', lw=2),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='mistyrose', edgecolor='darkred'))

    # Panel B: Deconvolved signal - tails "added back" to true peaks
    ax2 = axes[1]

    # Calculate what the true signal should be (tails added back to peaks)
    # Total signal for each peak = sum of all its spread signal
    true_counts = np.zeros(len(positions), dtype=int)
    colors_arr = ['#888888'] * len(positions)

    for peak, counts in zip(true_peaks, individual_counts):
        # Sum all the spread signal and put it back at the true position
        total_peak_signal = int(np.sum(counts))
        true_counts[peak["pos"]] = total_peak_signal
        colors_arr[peak["pos"]] = peak["color"]

    # Draw A-tract and bars
    ax2.axvspan(atract_start, atract_end, alpha=0.2, color='gold', label='Genomic A-tract')
    ax2.axvline(atract_end, color='orange', linestyle='-', linewidth=3, alpha=0.9)
    ax2.bar(positions, true_counts, width=0.8, color=colors_arr, edgecolor='black', linewidth=0.3)

    # Mark true positions
    for peak in true_peaks:
        ax2.axvline(peak["pos"], color=peak["color"], linestyle='--', linewidth=2.5, alpha=0.9,
                   label=f'{peak["label"]} (pos {peak["pos"]})')

    ax2.set_title('B. Deconvolved Signal (oligo-A tails added back)', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Genomic Position (bp)', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Read Counts', fontsize=14, fontweight='bold')
    ax2.tick_params(axis='both', labelsize=12)

    ax2.annotate('All signal\nat true CPA',
                xy=(55, 42), xytext=(45, 35),
                fontsize=11, ha='center',
                arrowprops=dict(arrowstyle='->', color='#3498DB', lw=2),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', edgecolor='#3498DB'))

    # Match y-axis scales and x limits
    max_y = max(display_observed.max(), true_counts.max()) * 1.15
    ax1.set_ylim(0, max_y)
    ax2.set_ylim(0, max_y)
    ax1.set_xlim(35, 80)
    ax2.set_xlim(35, 80)

    for ax in axes:
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_axisbelow(True)
        ax.legend(loc='upper left', fontsize=11, framealpha=0.95)

    plt.tight_layout()

    output_path = OUTPUT_DIR / "oligo_a_deconvolution.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")

    alt_path = Path(__file__).resolve().parent.parent / "oligo_a_deconvolution.png"
    plt.savefig(alt_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Saved: {alt_path}")

    plt.close()

    return output_path


def create_adaptive_clustering_visualization():
    """Create visualization of adaptive valley clustering with discrete bars and cluster stats."""

    fig, ax = plt.subplots(figsize=(14, 8))

    # Discrete positions
    positions = np.arange(0, 100)
    np.random.seed(42)

    # Define clusters with realistic patterns: 1-2 major peaks + minor signals nearby
    # Each cluster has: start, end, peaks as (pos, count) pairs
    clusters = [
        {
            "start": 12, "end": 38, "color": "#FFB3BA", "name": "Cluster 1",
            "peaks": [(22, 28), (24, 35), (25, 18), (27, 8), (29, 4)],  # Major at 24
        },
        {
            "start": 38, "end": 68, "color": "#BAFFC9", "name": "Cluster 2",
            "peaks": [(52, 22), (54, 48), (55, 42), (56, 25), (57, 12), (58, 5), (60, 3)],  # Major at 54-55
        },
        {
            "start": 68, "end": 92, "color": "#BAE1FF", "name": "Cluster 3",
            "peaks": [(78, 15), (80, 32), (81, 18), (82, 8), (84, 3)],  # Major at 80
        },
    ]

    # Generate counts
    counts = np.zeros(len(positions), dtype=int)

    # Add cluster peaks
    for cluster in clusters:
        for pos, count in cluster["peaks"]:
            counts[pos] = count

    # Add background noise (0, 1, 2)
    noise = np.random.choice([0, 0, 0, 0, 0, 1, 1, 2], size=len(positions))
    counts += noise

    # Color bars by cluster
    bar_colors = ['#888888'] * len(positions)
    for cluster in clusters:
        for pos in range(cluster["start"], cluster["end"]):
            if counts[pos] > 2:
                bar_colors[pos] = cluster["color"]

    # Plot bars
    ax.bar(positions, counts, width=0.8, color=bar_colors, edgecolor='black', linewidth=0.3)

    # Draw cluster boundaries and add annotations
    for i, cluster in enumerate(clusters):
        start, end = cluster["start"], cluster["end"]

        # Shade cluster region
        ax.axvspan(start, end, alpha=0.2, color=cluster["color"])

        # Draw boundary lines
        ax.axvline(start, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.axvline(end, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)

        # Calculate cluster statistics
        cluster_counts = counts[start:end]
        peak_signal = max(cluster_counts)
        total_signal = sum(cluster_counts)
        peak_coord = start + np.argmax(cluster_counts)

        # Mark peak with red triangle
        ax.plot(peak_coord, peak_signal + 2, 'rv', markersize=14, markeredgecolor='darkred', markeredgewidth=1.5)

        # Add cluster info box
        info_text = f'{cluster["name"]}\n' \
                   f'Peak: {peak_signal} @ pos {peak_coord}\n' \
                   f'Total: {total_signal} reads\n' \
                   f'Range: [{start}, {end})'

        # Position info box above cluster (raised to avoid valley markers)
        center_x = (start + end) / 2
        ax.text(center_x, max(counts) * 1.12, info_text, ha='center', va='bottom', fontsize=9,
                bbox=dict(facecolor=cluster["color"], alpha=0.85, edgecolor='gray', pad=4),
                family='monospace')

    # Mark valleys with green triangles
    valley_positions = [38, 68]  # Between clusters
    for vpos in valley_positions:
        ax.plot(vpos, counts[vpos] + 1.5, 'g^', markersize=12, markeredgecolor='darkgreen', markeredgewidth=1.5)

    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='v', color='w', markerfacecolor='red',
               markeredgecolor='darkred', markersize=12, label='Peak (local max)'),
        Line2D([0], [0], marker='^', color='w', markerfacecolor='green',
               markeredgecolor='darkgreen', markersize=10, label='Valley (local min)'),
        Line2D([0], [0], color='gray', linestyle='--', linewidth=2, label='Cluster boundary'),
    ]
    ax.legend(handles=legend_elements, loc='center right', fontsize=10)

    # Add algorithm description box
    textstr = '\n'.join([
        'Algorithm:',
        '1. Find peaks (local maxima)',
        '2. Find valleys (local minima)',
        '3. Boundary = midpoint(peak, valley)',
        '   capped at ±10bp from peak'
    ])
    props = dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.9, edgecolor='gray')
    ax.text(0.02, 0.72, textstr, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=props, family='monospace')

    # Styling
    ax.set_xlabel('Genomic Position (bp)', fontsize=12, fontweight='bold')
    ax.set_ylabel("3' End Read Counts", fontsize=12, fontweight='bold')
    ax.set_title('Adaptive Valley Clustering', fontsize=14, fontweight='bold')

    ax.set_xlim(0, 100)
    ax.set_ylim(0, max(counts) * 1.35)  # Extra space for info boxes

    ax.grid(True, alpha=0.3, axis='y')
    ax.set_axisbelow(True)

    plt.tight_layout()

    output_path = OUTPUT_DIR / "adaptive_clustering.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")

    alt_path = Path(__file__).resolve().parent.parent / "adaptive_clustering.png"
    plt.savefig(alt_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Saved: {alt_path}")

    plt.close()

    return output_path


if __name__ == "__main__":
    print("Generating visualizations...")
    create_spreading_visualization()
    create_deconvolution_visualization()
    create_adaptive_clustering_visualization()
    print("Done!")
