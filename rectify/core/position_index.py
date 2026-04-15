#!/usr/bin/env python3
"""
Position index writer for RECTIFY.

Writes a compact, pre-aggregated position count file alongside the full
per-read TSV.  The index is ~300× smaller than the full TSV and is used
by manifest-mode analysis to skip the heavy per-read TSV in both passes.

Index format (gzip-compressed TSV)::

    chrom   corrected_3prime  strand  count
    chrI    12836             +       1.000000
    chrI    41981             +       3.000000

Author: Kevin R. Roy
Date: 2026-03-09
"""

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Union
import gzip
import logging

logger = logging.getLogger(__name__)


def write_position_index(
    results_or_path: Union[List[Dict], str, Path],
    output_tsv_path: str,
) -> None:
    """
    Write a compact position index (chrom, corrected_3prime, strand → count).

    The index path is derived from *output_tsv_path* by replacing ``.tsv`` or
    ``.tsv.gz`` with ``_index.bed.gz``.

    Args:
        results_or_path: List of correction result dicts, **or** a path to an
            existing corrected TSV file (str or Path).
        output_tsv_path: Path to the main corrected TSV; index path is derived
            from this.
    """
    # Derive index path
    _base = str(output_tsv_path)
    if _base.endswith('.tsv.gz'):
        _index_path = _base[:-len('.tsv.gz')] + '_index.bed.gz'
    elif _base.endswith('.tsv'):
        _index_path = _base[:-len('.tsv')] + '_index.bed.gz'
    else:
        _index_path = _base + '_index.bed.gz'

    pos_counts: Dict = defaultdict(float)

    if isinstance(results_or_path, (str, Path)):
        import pandas as pd
        try:
            _df = pd.read_csv(
                str(results_or_path), sep='\t',
                usecols=['chrom', 'corrected_3prime', 'strand', 'fraction'],
            )
            for _, row in _df.iterrows():
                pos_counts[(row['chrom'], row['corrected_3prime'], row['strand'])] += float(row.get('fraction', 1.0))
        except (ValueError, KeyError):
            _df = pd.read_csv(
                str(results_or_path), sep='\t',
                usecols=['chrom', 'corrected_3prime', 'strand'],
            )
            for _, row in _df.iterrows():
                pos_counts[(row['chrom'], row['corrected_3prime'], row['strand'])] += 1.0
    else:
        for r in results_or_path:
            key = (r['chrom'], r['corrected_3prime'], r['strand'])
            pos_counts[key] += float(r.get('fraction', 1.0))

    with gzip.open(_index_path, 'wt') as f:
        f.write('chrom\tcorrected_3prime\tstrand\tcount\n')
        for (chrom, pos, strand), count in sorted(
            pos_counts.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])
        ):
            f.write(f'{chrom}\t{pos}\t{strand}\t{count:.6f}\n')

    logger.info(f"Position index written to {_index_path}")
