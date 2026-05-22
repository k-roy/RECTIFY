"""
Per-condition strand-specific bedgraph emission from a corrected-positions
DataFrame.

Produces single-base 3'-end pileup tracks (one row per CPA position with
optional RPM normalization). For the manifest-mode equivalent that
streams per-sample TSVs and accumulates per-condition counts without
holding the combined frame in memory, see `manifest._run_analyze_manifest`.
"""

from collections import defaultdict
from pathlib import Path

import pandas as pd

from .deseq2 import extract_condition_from_sample


# Chromosome order for sorted bedgraph output (yeast).
CHROM_ORDER = [
    'chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII',
    'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrmt'
]


def generate_bedgraphs(
    positions_df: pd.DataFrame,
    output_dir: Path,
    sample_column: str = 'sample',
    position_column: str = 'corrected_3prime',
    normalize_rpm: bool = True,
) -> None:
    """
    Generate strand-specific bedgraph files per condition from corrected positions.

    Args:
        positions_df: DataFrame with corrected positions
        output_dir: Output directory for bedgraph files
        sample_column: Column name for sample/replicate identifier
        position_column: Column name for corrected position
        normalize_rpm: If True, output RPM-normalized values
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract condition from sample name
    conditions = positions_df[sample_column].apply(
        lambda x: extract_condition_from_sample(x)
    ).unique()

    print(f"  Generating bedgraphs for {len(conditions)} conditions...")

    for condition in conditions:
        # Filter to this condition
        condition_mask = positions_df[sample_column].apply(
            lambda x: extract_condition_from_sample(x) == condition
        )
        cond_df = positions_df[condition_mask]

        if len(cond_df) == 0:
            continue

        # Count per chrom/strand/position.  Use pre-aggregated or fractional
        # weights when present so bedgraph totals match clustering/count matrices.
        plus_counts = defaultdict(lambda: defaultdict(float))
        minus_counts = defaultdict(lambda: defaultdict(float))
        weight_col = None
        if 'count' in cond_df.columns:
            weight_col = 'count'
        elif 'fraction' in cond_df.columns:
            weight_col = 'fraction'

        plus_df = cond_df[cond_df['strand'] == '+']
        minus_df = cond_df[cond_df['strand'] == '-']

        if weight_col:
            plus_grouped = plus_df.groupby(['chrom', position_column])[weight_col].sum()
            minus_grouped = minus_df.groupby(['chrom', position_column])[weight_col].sum()
            total_reads = float(cond_df[weight_col].sum())
        else:
            plus_grouped = plus_df.groupby(['chrom', position_column]).size()
            minus_grouped = minus_df.groupby(['chrom', position_column]).size()
            total_reads = float(len(cond_df))

        for (chrom, pos), count in plus_grouped.items():
            plus_counts[chrom][pos] += float(count)

        for (chrom, pos), count in minus_grouped.items():
            minus_counts[chrom][pos] += float(count)

        rpm_factor = 1e6 / total_reads if normalize_rpm and total_reads > 0 else 1.0

        # Write bedgraph files atomically: write to a temp file first, then
        # rename into place so that interrupted writes never leave a partial
        # (corrupt) file at the final output path.
        import os as _os
        for strand_name, counts_dict in [('plus', plus_counts), ('minus', minus_counts)]:
            output_path = output_dir / f"{condition}_{strand_name}.bedgraph"
            tmp_path = output_path.with_suffix('.tmp')

            try:
                with open(tmp_path, 'w') as f:
                    f.write(f'track type=bedGraph name="{condition}_{strand_name}" '
                            f'description="RECTIFY 3\' ends ({strand_name} strand)"\n')

                    ordered_chroms = [
                        *[chrom for chrom in CHROM_ORDER if chrom in counts_dict],
                        *sorted(chrom for chrom in counts_dict if chrom not in CHROM_ORDER),
                    ]
                    for chrom in ordered_chroms:
                        if chrom not in counts_dict:
                            continue

                        for pos in sorted(counts_dict[chrom].keys()):
                            count = counts_dict[chrom][pos]
                            value = count * rpm_factor if normalize_rpm else count
                            # corrected_3prime is 0-based-inclusive (derived from
                            # pysam reference_end - 1 or reference_start, both
                            # 0-based). BED is 0-based half-open, so a single
                            # base at 0-based pos is [pos, pos+1).
                            start = int(pos)
                            end = int(pos) + 1
                            f.write(f"{chrom}\t{start}\t{end}\t{value:.4f}\n")

                _os.rename(tmp_path, output_path)
            except Exception:
                if tmp_path.exists():
                    tmp_path.unlink()
                raise
