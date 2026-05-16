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

        # Count per chrom/strand/position
        plus_counts = defaultdict(lambda: defaultdict(int))
        minus_counts = defaultdict(lambda: defaultdict(int))

        plus_df = cond_df[cond_df['strand'] == '+']
        minus_df = cond_df[cond_df['strand'] == '-']

        for (chrom, pos), count in plus_df.groupby(['chrom', position_column]).size().items():
            plus_counts[chrom][pos] += count

        for (chrom, pos), count in minus_df.groupby(['chrom', position_column]).size().items():
            minus_counts[chrom][pos] += count

        total_reads = len(cond_df)
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

                    for chrom in CHROM_ORDER:
                        if chrom not in counts_dict:
                            continue

                        for pos in sorted(counts_dict[chrom].keys()):
                            count = counts_dict[chrom][pos]
                            value = count * rpm_factor if normalize_rpm else count
                            start = int(pos) - 1
                            end = int(pos)
                            f.write(f"{chrom}\t{start}\t{end}\t{value:.4f}\n")

                _os.rename(tmp_path, output_path)
            except Exception:
                if tmp_path.exists():
                    tmp_path.unlink()
                raise
