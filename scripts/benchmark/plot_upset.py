#!/usr/bin/env python3
"""
plot_upset.py — UpSet plot from per-aligner correction membership data.

Reads upset_data.tsv produced by aligner_contribution_analysis.py
(--per-aligner-correct-dir mode) and generates an UpSet plot showing
which combinations of aligners are capable of correcting each read.

Usage:
    python scripts/benchmark/plot_upset.py \\
        --input  /path/to/upset_data.tsv \\
        --output /path/to/upset_plot.pdf \\
        [--top N] [--min-count N] [--exclude-none] [--title "..."]
        [--width W] [--height H]
"""

import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input", required=True, help="upset_data.tsv from aligner_contribution_analysis.py")
    p.add_argument("--output", default="upset_plot.pdf", help="Output PDF or PNG path")
    p.add_argument("--top", type=int, default=25, help="Keep only top N patterns by count (default: 25)")
    p.add_argument("--min-count", type=int, default=10, help="Minimum read count per pattern (default: 10)")
    p.add_argument("--exclude-none", action="store_true",
                   help="Exclude the '(none_corrected)' row (reads corrected by no aligner)")
    p.add_argument("--title", default="", help="Plot title")
    p.add_argument("--width", type=float, default=0.0, help="Figure width in inches (auto if 0)")
    p.add_argument("--height", type=float, default=6.0, help="Figure height in inches (default: 6)")
    args = p.parse_args()

    try:
        from upsetplot import from_memberships
        from upsetplot import plot as upset_plot
    except ImportError:
        print("upsetplot not installed. Run: pip install upsetplot", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(args.input, sep="\t")

    member_cols = [c for c in df.columns if c.startswith("member_")
                   and c not in ("member_aligners",)]
    aligners = [c.replace("member_", "") for c in member_cols]
    if not aligners:
        print("[ERROR] No member_{aligner} columns found in input.", file=sys.stderr)
        sys.exit(1)

    if args.exclude_none:
        df = df[df["member_aligners"] != "(none_corrected)"]

    df = df[df["count"] >= args.min_count].nlargest(args.top, "count")
    if df.empty:
        print("No patterns remain after filtering.", file=sys.stderr)
        sys.exit(0)

    # Build upsetplot input: one membership tuple per read
    membership_list = []
    for _, row in df.iterrows():
        members = tuple(sorted(a for a in aligners if row.get(f"member_{a}", 0) == 1))
        membership_list.extend([members] * int(row["count"]))

    data = from_memberships(membership_list)

    n_shown = len(df)
    w = args.width if args.width > 0 else max(10.0, n_shown * 0.55)
    fig = plt.figure(figsize=(w, args.height))
    upset_plot(data, fig=fig, show_counts=True, show_setlabels=True,
               totals_plot_elements=3)

    if args.title:
        fig.suptitle(args.title, fontsize=11, y=1.01)

    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"Saved: {args.output}")

    # Console summary
    n_reads = df["count"].sum()
    print(f"\nTop patterns ({n_reads:,} reads across {n_shown} intersections shown):")
    for _, row in df.iterrows():
        pct = 100.0 * row["count"] / n_reads if n_reads else 0.0
        print(f"  {row['member_aligners']:<45} {row['count']:>8,}  ({pct:5.1f}%)")


if __name__ == "__main__":
    main()
