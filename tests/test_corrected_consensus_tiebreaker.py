#!/usr/bin/env python3
"""Tests for paralog-aware tiebreaker in :func:`merge_corrected_tsvs`.

Pins the contract that the fallback sort key (used when ``hp_edit_distance``
is N/A) groups ``_n_agree`` by ``(read_id, chrom, corrected_3prime)`` — not
just ``(read_id, corrected_3prime)``. Without the chrom in the groupby,
paralogous alignments at the same numeric position on different
chromosomes spuriously combine, so a single-aligner cross-chrom outlier
can tie with the multi-aligner same-chrom consensus.

Observed regression (cat1_plus_1 / read 0cb5a111):
  - chrXIV:10610  ← minimap2, gapmm2, uLTRA agree (3 aligners)
  - chrVI:8703    ← deSALT outlier (1 aligner)
The numeric positions are 1907 bp apart, but both loci are subtelomeric
repeats with near-identical sequence. Without `chrom` in `_n_agree`'s
groupby, both rows showed `_n_agree=4` (sum of two unrelated 3-and-1
position counts), and the tail of the fallback sort (span, n_junc) picked
deSALT's chrVI alignment.

After the fix:
  - chrXIV rows show `_n_agree=3` (correct: three aligners at this chrom+pos)
  - chrVI row shows `_n_agree=1` (correct: one aligner at this chrom+pos)
  - The fallback sort therefore prefers chrXIV.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rectify.core.consensus.corrected_consensus import merge_corrected_tsvs


# The minimum schema merge_corrected_tsvs needs in the fallback path
# (no BAMs supplied → no hp_edit_distance, so the legacy sort key fires).
_REQUIRED_COLUMNS = [
    'read_id', 'chrom', 'strand',
    'original_3prime', 'corrected_3prime',
    'five_prime_position', 'five_prime_rescued',
    'alignment_start', 'alignment_end',
    'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
    'polya_length', 'aligned_a_length', 'soft_clip_a_length',
    'junctions', 'n_junctions',
    'five_prime_soft_clip_length', 'three_prime_soft_clip_length',
    'mapq', 'correction_applied', 'confidence', 'qc_flags', 'fraction',
]


def _make_row(
    read_id: str,
    chrom: str,
    corrected_3p: int,
    *,
    strand: str = '+',
    alignment_start: int = 0,
    alignment_end: int = 100,
    n_junctions: int = 0,
    confidence: str = 'high',
) -> dict:
    return {
        'read_id': read_id,
        'chrom': chrom,
        'strand': strand,
        'original_3prime': corrected_3p,
        'corrected_3prime': corrected_3p,
        'five_prime_position': alignment_start,
        'five_prime_rescued': False,
        'alignment_start': alignment_start,
        'alignment_end': alignment_end,
        'ambiguity_min': corrected_3p,
        'ambiguity_max': corrected_3p,
        'ambiguity_range': 0,
        'polya_length': 0,
        'aligned_a_length': 0,
        'soft_clip_a_length': 0,
        'junctions': '',
        'n_junctions': n_junctions,
        'five_prime_soft_clip_length': 0,
        'three_prime_soft_clip_length': 0,
        'mapq': 60,
        'correction_applied': '',
        'confidence': confidence,
        'qc_flags': '',
        'fraction': 1.0,
    }


def _write_tsv(path: Path, rows: list) -> None:
    df = pd.DataFrame(rows, columns=_REQUIRED_COLUMNS)
    df.to_csv(path, sep='\t', index=False)


def test_paralog_tiebreaker_picks_multi_aligner_consensus(tmp_path):
    """A 3-aligner chrXIV consensus must beat a 1-aligner chrVI outlier
    even when both place the read at the same numeric corrected_3prime.
    """
    read_id = 'paralog_read'
    pos = 10610  # the magic position
    # Three aligners agree on chrXIV:10610
    minimap2_tsv = tmp_path / 'minimap2.tsv'
    gapmm2_tsv   = tmp_path / 'gapmm2.tsv'
    ultra_tsv    = tmp_path / 'uLTRA.tsv'
    # One aligner picks chrVI:10610 (paralog)
    desalt_tsv   = tmp_path / 'deSALT.tsv'

    _write_tsv(minimap2_tsv, [_make_row(read_id, 'chrXIV', pos)])
    _write_tsv(gapmm2_tsv,   [_make_row(read_id, 'chrXIV', pos)])
    _write_tsv(ultra_tsv,    [_make_row(read_id, 'chrXIV', pos)])
    _write_tsv(desalt_tsv,   [_make_row(read_id, 'chrVI',  pos)])

    out_tsv = tmp_path / 'merged.tsv'
    merge_corrected_tsvs(
        per_aligner_tsvs={
            'minimap2': minimap2_tsv,
            'gapmm2':   gapmm2_tsv,
            'uLTRA':    ultra_tsv,
            'deSALT':   desalt_tsv,
        },
        output_tsv=out_tsv,
        # No BAMs → fallback sort key path (the one we fixed).
        per_aligner_corrected_bams=None,
    )

    merged = pd.read_csv(out_tsv, sep='\t')
    assert len(merged) == 1, 'one read in → one row out'
    row = merged.iloc[0]
    assert row['chrom'] == 'chrXIV', (
        f'paralog tiebreaker should pick chrXIV (3-aligner consensus), '
        f'got chrom={row["chrom"]!r} winning_aligner={row["winning_aligner"]!r}'
    )
    # winning_aligner can be any of the three chrXIV aligners — they're
    # genuinely tied at this point. What matters is that it's NOT deSALT.
    assert row['winning_aligner'] != 'deSALT', (
        f'paralog tiebreaker must NOT pick the lone cross-chrom outlier'
    )


def test_majority_consensus_picks_chrXIV_even_when_outlier_has_wider_span(tmp_path):
    """The deSALT outlier with a wider span must still lose to the
    same-chrom 3-aligner consensus when _n_agree is computed correctly."""
    read_id = 'paralog_read'
    pos = 10610

    minimap2_tsv = tmp_path / 'minimap2.tsv'
    gapmm2_tsv   = tmp_path / 'gapmm2.tsv'
    ultra_tsv    = tmp_path / 'uLTRA.tsv'
    desalt_tsv   = tmp_path / 'deSALT.tsv'

    # Three aligners: short same-chrom alignment
    _write_tsv(minimap2_tsv, [_make_row(read_id, 'chrXIV', pos,
                                        alignment_start=10500, alignment_end=10650)])
    _write_tsv(gapmm2_tsv,   [_make_row(read_id, 'chrXIV', pos,
                                        alignment_start=10500, alignment_end=10650)])
    _write_tsv(ultra_tsv,    [_make_row(read_id, 'chrXIV', pos,
                                        alignment_start=10500, alignment_end=10650)])
    # deSALT: wider span on the paralog
    _write_tsv(desalt_tsv,   [_make_row(read_id, 'chrVI', pos,
                                        alignment_start=8000, alignment_end=10650)])

    out_tsv = tmp_path / 'merged.tsv'
    merge_corrected_tsvs(
        per_aligner_tsvs={
            'minimap2': minimap2_tsv,
            'gapmm2':   gapmm2_tsv,
            'uLTRA':    ultra_tsv,
            'deSALT':   desalt_tsv,
        },
        output_tsv=out_tsv,
        per_aligner_corrected_bams=None,
    )
    merged = pd.read_csv(out_tsv, sep='\t')
    assert merged.iloc[0]['chrom'] == 'chrXIV', (
        'wider deSALT span must NOT flip the winner when _n_agree groups '
        'by (chrom, pos) correctly'
    )


def test_single_aligner_passthrough_unaffected(tmp_path):
    """Trivial-case path (only one aligner succeeded) doesn't use _n_agree
    at all — must still emit the single row with winning_aligner set."""
    read_id = 'lone_read'
    only_tsv = tmp_path / 'minimap2.tsv'
    _write_tsv(only_tsv, [_make_row(read_id, 'chrI', 1000)])

    out_tsv = tmp_path / 'merged.tsv'
    merge_corrected_tsvs(
        per_aligner_tsvs={'minimap2': only_tsv},
        output_tsv=out_tsv,
        per_aligner_corrected_bams=None,
    )
    merged = pd.read_csv(out_tsv, sep='\t')
    assert len(merged) == 1
    assert merged.iloc[0]['winning_aligner'] == 'minimap2'
    assert merged.iloc[0]['chrom'] == 'chrI'
