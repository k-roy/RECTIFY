#!/usr/bin/env python3
"""
Position index writer for RECTIFY.

Writes a compact, pre-aggregated position count file alongside the full
per-read TSV.  The index is ~300× smaller than the full TSV and is used
by manifest-mode analysis to skip the heavy per-read TSV in both passes.

Index format (gzip-compressed TSV)::

    chrom   corrected_3prime  strand  count   count_ag_rich
    chrI    12836             +       1.000000  0.000000
    chrI    41981             +       3.000000  1.000000

The ``count_ag_rich`` column records the subset of ``count`` whose underlying
reads carry the ``AG_RICH`` qc_flag (oligo-dT internal-priming candidates).
It is emitted only when the source TSV has a ``qc_flags`` column; older
indices without the column are loaded with ``count_ag_rich = 0``.

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

    # Each value is [count_total, count_ag_rich]. count_ag_rich is the subset
    # of count whose underlying reads were flagged AG_RICH (oligo-dT internal-
    # priming candidates). Tracking it here lets manifest-mode analyze
    # propagate AG status to cluster-level annotation without re-streaming
    # the full TSV.
    pos_counts: Dict = defaultdict(lambda: [0.0, 0.0])

    if isinstance(results_or_path, (str, Path)):
        import pandas as pd
        _header = pd.read_csv(str(results_or_path), sep='\t', nrows=0)
        _has_fraction = 'fraction' in _header.columns
        _has_qc = 'qc_flags' in _header.columns
        _usecols = ['chrom', 'corrected_3prime', 'strand']
        if _has_fraction:
            _usecols.append('fraction')
        if _has_qc:
            _usecols.append('qc_flags')
        try:
            _df = pd.read_csv(str(results_or_path), sep='\t', usecols=_usecols)
        except (ValueError, KeyError):
            _df = pd.read_csv(
                str(results_or_path), sep='\t',
                usecols=['chrom', 'corrected_3prime', 'strand'],
            )
            _has_fraction = False
            _has_qc = False
        for _, row in _df.iterrows():
            key = (row['chrom'], row['corrected_3prime'], row['strand'])
            w = float(row.get('fraction', 1.0)) if _has_fraction else 1.0
            pos_counts[key][0] += w
            if _has_qc and 'AG_RICH' in str(row.get('qc_flags', '')):
                pos_counts[key][1] += w
    else:
        for r in results_or_path:
            key = (r['chrom'], r['corrected_3prime'], r['strand'])
            w = float(r.get('fraction', 1.0))
            pos_counts[key][0] += w
            if 'AG_RICH' in str(r.get('qc_flags', '')):
                pos_counts[key][1] += w

    with gzip.open(_index_path, 'wt') as f:
        f.write('chrom\tcorrected_3prime\tstrand\tcount\tcount_ag_rich\n')
        for (chrom, pos, strand), (count, count_ag) in sorted(
            pos_counts.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])
        ):
            f.write(f'{chrom}\t{pos}\t{strand}\t{count:.6f}\t{count_ag:.6f}\n')

    logger.info(f"Position index written to {_index_path}")
