"""
calibrate_junction_overhang.py — Derive empirical min_overhang(intron_size) curve.

For each intron size bin, computes the minimum exon overhang required to treat
an aligner's junction call as credible.  Calibration uses concordant multi-aligner
reads: junctions reported by ≥ N_CONCORDANT aligners for the same read are treated
as ground-truth, and the minimum overhang seen across those agreeing aligners is
recorded as the evidence floor for that intron size.

The resulting table is used by merge_corrected_tsvs() via the --junction-overhang-table
flag to penalise (not disqualify) aligner wins that rely on junctions with insufficient
flanking evidence relative to their size.  The penalty is soft: a flagged aligner
is sorted after non-flagged aligners, but is still chosen if all aligners are flagged.

Usage (standalone):
    python -m rectify.core.splice.calibrate_junction_overhang \\
        --tsv-dir /path/to/sample/dir \\   # subdirs named by aligner
        --output overhang_table.tsv

    python -m rectify.core.splice.calibrate_junction_overhang \\
        --tsvs mapPacBio:/path/a.tsv deSALT:/path/b.tsv gapmm2:/path/c.tsv \\
        --output overhang_table.tsv

    python -m rectify.core.splice.calibrate_junction_overhang \\
        --merged-tsv corrected_reads.tsv \\   # merged TSV with winning_aligner column
        --output overhang_table.tsv           # (concordance from winners only — less ideal)

Output format (TSV):
    intron_size_max   min_overhang
    100               1
    200               2
    500               5
    ...
    inf               30           # catch-all for sizes beyond calibration data

Author: Kevin R. Roy <kevinrjroy@gmail.com>
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Default aligner subdirectory names to look for when --tsv-dir is used
_DEFAULT_ALIGNERS = ['mapPacBio', 'minimap2', 'gapmm2', 'uLTRA', 'deSALT']

# Number of aligners that must report the same junction for it to be treated
# as concordant ground-truth
N_CONCORDANT_DEFAULT = 4

# Percentile of the minimum-overhang distribution used as the threshold.
# 5th percentile: 95% of real junctions have at least this much overhang.
PERCENTILE_DEFAULT = 5

# Log2 intron-size bin edges (exclusive upper bounds).
# Chosen to give fine resolution over the typical intron range (50–600 bp)
# and coarser resolution for the rare-but-important large-intron region.
_BIN_EDGES = [
    50, 75, 100, 125, 150, 200, 250, 300, 400, 500,
    600, 750, 1000, 1500, 2000, 3000, 5000, float('inf'),
]

# Fallback row added to the end of the table when no observations exist for
# the largest intron sizes.  Forces 30 bp overhang for anything > calibration
# range, preventing chimeric long-range artifacts.
_FALLBACK_LARGE_INTRON_OVERHANG = 30

# Minimum number of observations required per bin to emit a reliable estimate.
_MIN_BIN_OBS = 20


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def _bare_uuid(name: str) -> str:
    """Strip mapPacBio _pt:i:N suffix to recover the canonical bare UUID."""
    if "_pt:i:" in name:
        return name.split("_pt:i:")[0]
    if " pt:i:" in name:
        return name.split(" pt:i:")[0]
    return name


def _parse_junctions(junc_str) -> List[Tuple[int, int]]:
    """Parse semicolon-separated 'start-end' junction string into a list of (start, end) tuples."""
    if pd.isna(junc_str) or not junc_str or str(junc_str) == 'nan':
        return []
    result = []
    for part in str(junc_str).split(';'):
        part = part.strip()
        if not part:
            continue
        fields = part.split('-')
        if len(fields) == 2:
            try:
                result.append((int(fields[0]), int(fields[1])))
            except ValueError:
                pass
    return result


def _load_tsv(path: Path, aligner: str) -> Optional[pd.DataFrame]:
    """Load one per-aligner corrected TSV, retaining only columns needed for calibration.

    Note: some per-aligner TSVs are written with ``index=True`` (minimap2, gapmm2, uLTRA)
    producing one more data column than header columns.  Using ``index_col=False`` prevents
    pandas from using the first column as the row index, keeping all column names aligned.
    """
    required = {'read_id', 'junctions', 'alignment_start', 'alignment_end'}
    try:
        df = pd.read_csv(path, sep='\t', low_memory=False, index_col=False)
    except Exception as exc:
        logger.warning("Cannot read %s: %s", path, exc)
        return None
    missing = required - set(df.columns)
    if missing:
        logger.warning("Skipping %s (missing columns: %s)", path, missing)
        return None
    df = df[['read_id', 'junctions', 'alignment_start', 'alignment_end']].copy()
    df['read_id'] = df['read_id'].astype(str).apply(_bare_uuid)
    df['_aligner'] = aligner
    # Drop rows without junctions (non-spliced reads don't contribute)
    df = df[df['junctions'].notna() & (df['junctions'].astype(str) != 'nan') & (df['junctions'] != '')]
    return df if not df.empty else None


def load_per_aligner_tsvs(
    tsv_dir: Optional[Path] = None,
    tsv_spec: Optional[List[str]] = None,
) -> Dict[str, pd.DataFrame]:
    """Load per-aligner corrected TSVs.

    Parameters
    ----------
    tsv_dir:
        Directory containing subdirectories named by aligner, each with
        corrected_reads.tsv.  If subdirectory does not exist or has no TSV,
        it is silently skipped.
    tsv_spec:
        List of "aligner:path" strings, e.g. ["mapPacBio:/a.tsv", "deSALT:/b.tsv"].

    Returns
    -------
    Dict mapping aligner name → DataFrame.
    """
    dfs: Dict[str, pd.DataFrame] = {}

    if tsv_dir is not None:
        tsv_dir = Path(tsv_dir)
        for aligner in _DEFAULT_ALIGNERS:
            candidate = tsv_dir / aligner / 'corrected_reads.tsv'
            if candidate.exists():
                df = _load_tsv(candidate, aligner)
                if df is not None:
                    dfs[aligner] = df
                    logger.info("Loaded %s: %d spliced reads", aligner, len(df))

    if tsv_spec:
        for spec in tsv_spec:
            if ':' not in spec:
                logger.warning("Skipping malformed --tsvs entry (expected 'aligner:path'): %s", spec)
                continue
            aligner, _, path_str = spec.partition(':')
            path = Path(path_str)
            df = _load_tsv(path, aligner)
            if df is not None:
                dfs[aligner] = df
                logger.info("Loaded %s from %s: %d spliced reads", aligner, path, len(df))

    return dfs


# ---------------------------------------------------------------------------
# Concordance and overhang extraction
# ---------------------------------------------------------------------------

def _explode_junctions(df: pd.DataFrame, aligner: str) -> pd.DataFrame:
    """Vectorized expansion of the junctions column into per-junction rows.

    For each row with a non-empty junction string (format: 'start-end;start-end'),
    creates one output row per (intron_start, intron_end) with left/right overhangs.

    Returns
    -------
    DataFrame with columns: read_id, _aligner, intron_start, intron_end, min_ov.
    Empty DataFrame if no spliced reads with parseable junctions.
    """
    empty = pd.DataFrame(columns=['read_id', '_aligner', 'intron_start', 'intron_end', 'min_ov'])

    # Keep only rows with actual junction coordinate strings (contains '-' separator)
    junc_mask = (
        df['junctions'].notna()
        & df['junctions'].astype(str).str.contains('-', na=False)
    )
    spliced = df[junc_mask].copy()
    if spliced.empty:
        return empty

    spliced = spliced[['read_id', 'alignment_start', 'alignment_end', 'junctions']].copy()
    spliced['alignment_start'] = pd.to_numeric(spliced['alignment_start'], errors='coerce')
    spliced['alignment_end']   = pd.to_numeric(spliced['alignment_end'],   errors='coerce')
    spliced = spliced.dropna(subset=['alignment_start', 'alignment_end'])

    # Parse junction strings into list of (js, je) tuples
    spliced['_junc_list'] = spliced['junctions'].apply(_parse_junctions)
    spliced = spliced[spliced['_junc_list'].apply(len) > 0]
    if spliced.empty:
        return empty

    # Explode: one row per junction
    spliced = spliced.explode('_junc_list')
    spliced = spliced[spliced['_junc_list'].apply(lambda x: isinstance(x, tuple) and len(x) == 2)]
    if spliced.empty:
        return empty
    spliced['intron_start'] = spliced['_junc_list'].apply(lambda x: x[0])
    spliced['intron_end']   = spliced['_junc_list'].apply(lambda x: x[1])

    # Compute per-side overhangs; alignment must enclose the junction
    spliced['left_ov']  = spliced['intron_start'] - spliced['alignment_start'].astype(int)
    spliced['right_ov'] = spliced['alignment_end'].astype(int) - spliced['intron_end']
    valid = spliced[(spliced['left_ov'] >= 0) & (spliced['right_ov'] >= 0)].copy()
    if valid.empty:
        return empty

    valid['min_ov']   = valid[['left_ov', 'right_ov']].min(axis=1)
    valid['_aligner'] = aligner

    return valid[['read_id', '_aligner', 'intron_start', 'intron_end', 'min_ov']].reset_index(drop=True)


def extract_concordant_junction_overhangs(
    dfs: Dict[str, pd.DataFrame],
    n_concordant: int = N_CONCORDANT_DEFAULT,
) -> pd.DataFrame:
    """Extract per-(read, junction) minimum overhang from concordant multi-aligner data.

    For each junction reported by ≥ n_concordant aligners for the same read:
        - Compute left_overhang  = junction_start  - alignment_start
        - Compute right_overhang = alignment_end - junction_end
        - Record intron_size = junction_end - junction_start
        - Record min_overhang = min across ALL concordant aligners reporting this junction

    Uses vectorized pandas operations (explode + groupby) — scales to millions of reads.

    Returns
    -------
    DataFrame with columns: intron_size, min_overhang.
    Each row represents one concordant junction observation (one per read × junction).
    """
    if len(dfs) < n_concordant:
        logger.warning(
            "Only %d aligners loaded; n_concordant=%d — reducing to %d",
            len(dfs), n_concordant, len(dfs),
        )
        n_concordant = max(2, len(dfs))

    # Expand each aligner's spliced reads into per-junction records (vectorized)
    pieces = []
    for aligner, df in dfs.items():
        exp = _explode_junctions(df, aligner)
        if not exp.empty:
            pieces.append(exp)
            logger.info("  %s: %d junction records", aligner, len(exp))

    if not pieces:
        logger.warning("No junction records extracted from input TSVs")
        return pd.DataFrame(columns=['intron_size', 'min_overhang'])

    rec_df = pd.concat(pieces, ignore_index=True)
    logger.info("Total aligner×read×junction triples: %d", len(rec_df))

    # Count how many aligners report each (read, junction) combination
    concordance = (
        rec_df
        .groupby(['read_id', 'intron_start', 'intron_end'])
        .agg(n_aligners=('_aligner', 'nunique'), min_ov=('min_ov', 'min'))
        .reset_index()
    )

    concordant = concordance[concordance['n_aligners'] >= n_concordant].copy()
    logger.info(
        "Concordant junctions (≥%d aligners): %d / %d (%.1f%%)",
        n_concordant, len(concordant), len(concordance),
        100 * len(concordant) / max(1, len(concordance)),
    )

    concordant['intron_size'] = concordant['intron_end'] - concordant['intron_start']
    concordant = concordant[concordant['intron_size'] > 0]

    return concordant[['intron_size', 'min_ov']].rename(columns={'min_ov': 'min_overhang'})


# ---------------------------------------------------------------------------
# Bin, percentile, and isotonic smoothing
# ---------------------------------------------------------------------------

def _isotonic_smooth(values: np.ndarray) -> np.ndarray:
    """Apply pool-adjacent-violators isotonic regression (non-decreasing).

    Ensures that larger introns require at least as much overhang as smaller ones.
    """
    n = len(values)
    if n == 0:
        return values
    result = values.copy().astype(float)
    # Pool-adjacent-violators (PAV) — non-decreasing
    i = 0
    while i < n:
        j = i + 1
        # Find end of current monotone run
        while j < n and result[j] >= result[j - 1]:
            j += 1
        if j == n:
            break
        # Found a violation at j — pool [i, j) to their mean
        start = i
        while start > 0 and result[start - 1] > result[start]:
            start -= 1
        pool_mean = float(np.mean(result[start:j + 1]))
        result[start:j + 1] = pool_mean
        i = start
    return result


def build_overhang_table(
    data: pd.DataFrame,
    percentile: float = PERCENTILE_DEFAULT,
    bin_edges: Optional[List[float]] = None,
    fallback_large_intron: float = _FALLBACK_LARGE_INTRON_OVERHANG,
    min_bin_obs: int = _MIN_BIN_OBS,
) -> pd.DataFrame:
    """Compute the overhang threshold table from concordant junction observations.

    For each intron-size bin:
      1. Collect all min_overhang values.
      2. Compute the specified percentile (default 5th).
      3. Apply isotonic regression to enforce monotone non-decreasing thresholds.

    Parameters
    ----------
    data:
        DataFrame with columns intron_size and min_overhang
        (from extract_concordant_junction_overhangs).
    percentile:
        Percentile of min_overhang to use as threshold (default 5).
    bin_edges:
        Upper-bound edges for intron-size bins (last element must be inf).
    fallback_large_intron:
        min_overhang value to use for bins with too few observations.
    min_bin_obs:
        Minimum observations required to emit a reliable percentile estimate.
        Bins below this count use linear interpolation from neighbours or fallback.

    Returns
    -------
    DataFrame with columns: intron_size_max (float), min_overhang (float).
    The last row has intron_size_max = inf and represents all introns beyond
    the calibration range.
    """
    if bin_edges is None:
        bin_edges = _BIN_EDGES

    assert bin_edges[-1] == float('inf'), "Last bin edge must be float('inf')"

    bin_labels = bin_edges  # upper bounds

    # Assign each observation to a bin
    data = data[data['intron_size'] > 0].copy()
    if data.empty:
        logger.warning("No data for overhang table — returning fallback table")
        return _fallback_table(bin_edges, fallback_large_intron)

    data['_bin'] = pd.cut(
        data['intron_size'],
        bins=[0] + bin_edges[:-1] + [float('inf')],
        labels=bin_labels,
        right=True,
        include_lowest=True,
    ).astype(float)

    # Compute percentile per bin
    rows = []
    for edge in bin_edges:
        bin_data = data[data['_bin'] == edge]['min_overhang']
        n = len(bin_data)
        if n >= min_bin_obs:
            pct = float(np.percentile(bin_data, percentile))
            rows.append({'intron_size_max': float(edge), 'min_overhang': pct, '_n': n})
        else:
            rows.append({'intron_size_max': float(edge), 'min_overhang': None, '_n': n})

    table = pd.DataFrame(rows)
    logger.info(
        "Bin statistics:\n%s",
        table[['intron_size_max', '_n', 'min_overhang']].to_string(index=False),
    )

    # Fill missing bins via forward/backward fill then fallback
    table['min_overhang'] = table['min_overhang'].fillna(method='ffill').fillna(method='bfill')
    table['min_overhang'] = table['min_overhang'].fillna(fallback_large_intron)

    # Enforce minimum of 1 nt (never require 0 overhang)
    table['min_overhang'] = table['min_overhang'].clip(lower=1.0)

    # Isotonic smoothing — ensures larger introns require at least as much as smaller ones
    vals = table['min_overhang'].to_numpy()
    table['min_overhang'] = _isotonic_smooth(vals)

    # Round to 1 decimal place (overhang is in bp, fractional values are from percentile
    # interpolation; round conservatively upward so we err on the side of stringency)
    table['min_overhang'] = np.ceil(table['min_overhang']).astype(float)

    logger.info("Overhang table (smoothed):\n%s",
                table[['intron_size_max', 'min_overhang']].to_string(index=False))

    return table[['intron_size_max', 'min_overhang']]


def _fallback_table(
    bin_edges: List[float],
    fallback_overhang: float,
) -> pd.DataFrame:
    """Return a conservative fallback table when no calibration data is available."""
    rows = []
    for edge in bin_edges:
        if edge <= 100:
            rows.append({'intron_size_max': float(edge), 'min_overhang': 1.0})
        elif edge <= 300:
            rows.append({'intron_size_max': float(edge), 'min_overhang': 3.0})
        elif edge <= 600:
            rows.append({'intron_size_max': float(edge), 'min_overhang': 5.0})
        else:
            rows.append({'intron_size_max': float(edge), 'min_overhang': fallback_overhang})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# OverhangTable — runtime lookup object
# ---------------------------------------------------------------------------

class OverhangTable:
    """Empirical min-overhang-per-intron-size lookup table.

    Loaded from the TSV produced by calibrate_junction_overhang.py.
    Used by merge_corrected_tsvs() to flag aligner×read pairs whose junction
    overhangs are too small for the intron size.

    Parameters
    ----------
    rows:
        Iterable of (intron_size_max, min_overhang) tuples.  The last row
        should have intron_size_max = float('inf') to serve as catch-all.
    """

    def __init__(self, rows: List[Tuple[float, float]]) -> None:
        self._rows: List[Tuple[float, float]] = sorted(rows, key=lambda r: r[0])

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_tsv(cls, path: "Path | str") -> 'OverhangTable':
        """Load from a TSV with columns intron_size_max and min_overhang."""
        df = pd.read_csv(str(path), sep='\t', comment='#')
        df.columns = [c.strip() for c in df.columns]
        if 'intron_size_max' not in df.columns or 'min_overhang' not in df.columns:
            raise ValueError(
                f"Overhang table {path} must have columns "
                f"'intron_size_max' and 'min_overhang'"
            )
        rows = []
        for _, row in df.iterrows():
            smax = row['intron_size_max']
            ov   = float(row['min_overhang'])
            smax = float('inf') if str(smax).strip().lower() in ('inf', 'infinity', 'nan') else float(smax)
            rows.append((smax, ov))
        return cls(rows)

    @classmethod
    def default(cls) -> 'OverhangTable':
        """Conservative default table used when no calibration file is supplied.

        Graduated thresholds for short introns; 22 bp flat from 750 bp to 5 kb;
        50 bp beyond that.  Rationale:

          - Short introns (< 100 bp) are always plausible even with 1 nt overhang
          - Mid-range introns (100–750 bp) need a modest ramp of evidence
          - From 750 bp to 5 kb: 22 bp flat — matches the empirically calibrated
            5th-percentile from concordant reads.  Many exon-1's are only 20–34 nt,
            so a 22 bp threshold avoids penalising reads that start close to the
            intron boundary.  Downstream junction analysis handles real vs. spurious
            junction calling; this filter's goal is the best alignment, not junction
            classification.
          - Beyond 5 kb: 50 bp, and the _add_chimera_flag quality check
            (hp_ed / overhang_bases < 20%) provides the main discrimination;
            very long "introns" in yeast are almost certainly chimeric artifacts
            from gapmm2/deSALT

        Intended as a sane starting point; calibrate from real data for best results.
        """
        return cls([
            (100,          1.0),
            (200,          2.0),
            (300,          3.0),
            (500,          5.0),
            (5000,        22.0),
            (float('inf'), 50.0),
        ])

    # ------------------------------------------------------------------
    # Lookup
    # ------------------------------------------------------------------

    def lookup(self, intron_size: int) -> float:
        """Return the minimum required overhang for an intron of this size."""
        for max_size, min_ov in self._rows:
            if intron_size <= max_size:
                return min_ov
        # Beyond all table rows — use the last row's value
        return self._rows[-1][1] if self._rows else float('inf')

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self._rows, columns=['intron_size_max', 'min_overhang'])

    def to_tsv(self, path: "Path | str") -> None:
        self.to_dataframe().to_csv(str(path), sep='\t', index=False, float_format='%g')

    def __repr__(self) -> str:
        return f"OverhangTable({self._rows!r})"


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def calibrate(
    tsv_dir: Optional[Path] = None,
    tsv_spec: Optional[List[str]] = None,
    output: Path = Path('overhang_table.tsv'),
    n_concordant: int = N_CONCORDANT_DEFAULT,
    percentile: float = PERCENTILE_DEFAULT,
) -> OverhangTable:
    """Run the full calibration pipeline and write the output TSV.

    Returns the OverhangTable object.
    """
    dfs = load_per_aligner_tsvs(tsv_dir=tsv_dir, tsv_spec=tsv_spec)
    if len(dfs) < 2:
        logger.error(
            "At least 2 aligner TSVs required for calibration; got %d. "
            "Returning default table.",
            len(dfs),
        )
        table = OverhangTable.default()
        table.to_tsv(output)
        return table

    logger.info(
        "Calibrating from %d aligners: %s", len(dfs), ', '.join(sorted(dfs))
    )

    # Keep this clamp in sync with the one inside
    # extract_concordant_junction_overhangs (n_concordant = max(2, len(dfs))).
    # Mirror the clamp inside extract_concordant_junction_overhangs so the
    # written header reflects the concordance threshold actually ENFORCED, not
    # the (possibly higher) requested default. Otherwise a 3-aligner run with the
    # default n_concordant=4 prints the impossible "n_aligners: 3  n_concordant: 4".
    effective_n_concordant = max(2, min(n_concordant, len(dfs)))
    obs = extract_concordant_junction_overhangs(dfs, n_concordant=n_concordant)
    if obs.empty:
        logger.warning("No concordant observations found. Returning default table.")
        table = OverhangTable.default()
        table.to_tsv(output)
        return table

    logger.info(
        "Total concordant junction observations: %d  "
        "(intron_size range: %d – %d bp)",
        len(obs), int(obs['intron_size'].min()), int(obs['intron_size'].max()),
    )

    table_df = build_overhang_table(obs, percentile=percentile)
    rows = list(zip(table_df['intron_size_max'].tolist(), table_df['min_overhang'].tolist()))
    table = OverhangTable(rows)

    # Write with a header comment
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w') as fh:
        fh.write(f"# Junction overhang calibration table\n")
        fh.write(f"# Generated by calibrate_junction_overhang.py\n")
        fh.write(f"# n_aligners: {len(dfs)}  n_concordant: {effective_n_concordant}  percentile: {percentile}\n")
        fh.write(f"# n_concordant_observations: {len(obs)}\n")
        fh.write(f"# intron_size_range: {int(obs['intron_size'].min())} – {int(obs['intron_size'].max())} bp\n")
        fh.write("intron_size_max\tmin_overhang\n")
        for max_size, min_ov in table._rows:
            size_str = 'inf' if max_size == float('inf') else str(int(max_size))
            fh.write(f"{size_str}\t{int(min_ov)}\n")

    logger.info("Overhang table written to %s", output)
    return table


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s  %(levelname)-7s  %(message)s',
        datefmt='%H:%M:%S',
    )

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    src = parser.add_mutually_exclusive_group()
    src.add_argument(
        '--tsv-dir', metavar='DIR',
        help='Directory with aligner subdirs (mapPacBio/, deSALT/, ...), '
             'each containing corrected_reads.tsv',
    )
    src.add_argument(
        '--tsvs', nargs='+', metavar='ALIGNER:PATH',
        help='Per-aligner corrected TSVs as "aligner:path" pairs',
    )

    parser.add_argument(
        '-o', '--output', default='overhang_table.tsv', metavar='TSV',
        help='Output TSV path (default: overhang_table.tsv)',
    )
    parser.add_argument(
        '--n-concordant', type=int, default=N_CONCORDANT_DEFAULT,
        metavar='N',
        help=f'Minimum aligners that must agree on a junction for it to be '
             f'treated as ground-truth (default: {N_CONCORDANT_DEFAULT})',
    )
    parser.add_argument(
        '--percentile', type=float, default=PERCENTILE_DEFAULT,
        metavar='PCT',
        help=f'Percentile of min-overhang distribution to use as threshold '
             f'(default: {PERCENTILE_DEFAULT})',
    )

    args = parser.parse_args()

    if not args.tsv_dir and not args.tsvs:
        parser.error('One of --tsv-dir or --tsvs is required')

    calibrate(
        tsv_dir=Path(args.tsv_dir) if args.tsv_dir else None,
        tsv_spec=args.tsvs,
        output=Path(args.output),
        n_concordant=args.n_concordant,
        percentile=args.percentile,
    )


if __name__ == '__main__':
    main()
