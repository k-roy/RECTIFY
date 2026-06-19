#!/usr/bin/env python3
"""Tests for the ``drop_chimeric_winners`` hard filter.

A read whose WINNING alignment is still chimera-flagged (``_effective_chimera_ok==1``)
has no QC-passing alignment across any aligner ("door A") — the soft flag only
re-orders within the read, so the consensus would otherwise emit the least-bad
spurious splice.  The hard filter drops these reads entirely.  Off by default
(K=0 / yeast byte-identical); a no-op when the HP-ED path didn't run (no
``_effective_chimera_ok`` column).

Validated context: A549 chr5 DRS door-A reads (no aligner places a clean junction);
a within-flagged *penalty* cannot fix them (see project-junction-dq-penalty) — only
dropping them does.

Author: Kevin R. Roy
"""
import os
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rectify.core.consensus.corrected_consensus import (
    _chimeric_winner_read_ids,
    merge_corrected_tsvs,
)


def _rep(rows):
    return pd.DataFrame(rows)


def _winners(mapping):
    return pd.DataFrame(
        [{'read_id': r, '_winning_aligner': a} for r, a in mapping.items()]
    )


# ── _chimeric_winner_read_ids ───────────────────────────────────────────────

def test_empty_when_no_effective_chimera_column():
    rep = _rep([{'read_id': 'r1', '_aligner': 'GMAP'}])
    assert _chimeric_winner_read_ids(rep, _winners({'r1': 'GMAP'})) == set()


def test_drops_read_whose_winner_is_flagged():
    """r1 winner flagged (all candidates flagged) → dropped; r2 winner clean → kept."""
    rep = _rep([
        {'read_id': 'r1', '_aligner': 'GMAP',     '_effective_chimera_ok': 1},
        {'read_id': 'r1', '_aligner': 'minimap2', '_effective_chimera_ok': 1},
        {'read_id': 'r2', '_aligner': 'minimap2', '_effective_chimera_ok': 0},
        {'read_id': 'r2', '_aligner': 'GMAP',     '_effective_chimera_ok': 1},
    ])
    winners = _winners({'r1': 'GMAP', 'r2': 'minimap2'})
    assert _chimeric_winner_read_ids(rep, winners) == {'r1'}


def test_only_the_winning_rows_flag_matters():
    """A clean winner is kept even if a non-winning aligner for that read is flagged."""
    rep = _rep([
        {'read_id': 'r3', '_aligner': 'minimap2', '_effective_chimera_ok': 0},  # winner
        {'read_id': 'r3', '_aligner': 'deSALT',   '_effective_chimera_ok': 1},  # loser
    ])
    assert _chimeric_winner_read_ids(rep, _winners({'r3': 'minimap2'})) == set()


def test_all_winners_flagged_drops_all():
    rep = _rep([
        {'read_id': 'r1', '_aligner': 'GMAP', '_effective_chimera_ok': 1},
        {'read_id': 'r2', '_aligner': 'GMAP', '_effective_chimera_ok': 1},
    ])
    assert _chimeric_winner_read_ids(
        rep, _winners({'r1': 'GMAP', 'r2': 'GMAP'})) == {'r1', 'r2'}


# ── Refinement: spare flagged winners that still have a salvageable junction ──

def test_refined_spares_winner_with_one_supported_junction():
    """3-good-1-noisy: a flagged winner with ≥1 supported junction is NOT dropped."""
    rep = _rep([
        {'read_id': 'r1', '_aligner': 'GMAP', '_effective_chimera_ok': 1,
         'junctions': '20-120;200-300'},
    ])
    stats = {(20, 120): (7, 15)}  # J1 supported (cnt 7≥5); J2 absent → unsupported
    assert _chimeric_winner_read_ids(
        rep, _winners({'r1': 'GMAP'}), junction_stats=stats) == set()


def test_refined_drops_winner_when_all_junctions_unsupported():
    rep = _rep([
        {'read_id': 'r1', '_aligner': 'GMAP', '_effective_chimera_ok': 1,
         'junctions': '20-120;200-300'},
    ])
    assert _chimeric_winner_read_ids(
        rep, _winners({'r1': 'GMAP'}), junction_stats={}) == {'r1'}


def test_refined_skips_flagged_winner_with_no_junctions():
    """No junctions (minimal-anchor case) → left to the separate minimal-anchor filter."""
    rep = _rep([
        {'read_id': 'r1', '_aligner': 'GMAP', '_effective_chimera_ok': 1, 'junctions': ''},
    ])
    assert _chimeric_winner_read_ids(
        rep, _winners({'r1': 'GMAP'}), junction_stats={}) == set()


def test_blunt_fallback_without_junction_stats():
    """Without junction_stats → blunt rule (all flagged winners), back-compat."""
    rep = _rep([
        {'read_id': 'r1', '_aligner': 'GMAP', '_effective_chimera_ok': 1,
         'junctions': '20-120'},
    ])
    assert _chimeric_winner_read_ids(rep, _winners({'r1': 'GMAP'})) == {'r1'}


# ── End-to-end: default OFF + no-BAM path is a safe no-op ────────────────────

_COLS = [
    'read_id', 'chrom', 'strand', 'original_3prime', 'corrected_3prime',
    'five_prime_position', 'five_prime_rescued', 'alignment_start', 'alignment_end',
    'ambiguity_min', 'ambiguity_max', 'ambiguity_range', 'polya_length',
    'aligned_a_length', 'soft_clip_a_length', 'junctions', 'n_junctions',
    'five_prime_soft_clip_length', 'three_prime_soft_clip_length', 'mapq',
    'correction_applied', 'confidence', 'qc_flags', 'fraction',
]


def _row(read_id, chrom, c3p):
    d = dict.fromkeys(_COLS, 0)
    d.update(read_id=read_id, chrom=chrom, strand='+', original_3prime=c3p,
             corrected_3prime=c3p, five_prime_rescued=False, alignment_end=100,
             junctions='', mapq=60, correction_applied='', confidence='high',
             qc_flags='', fraction=1.0)
    return d


def _write(path, rows):
    pd.DataFrame(rows, columns=_COLS).to_csv(path, sep='\t', index=False)


@pytest.mark.parametrize("drop_flag", [False, True])
def test_no_bam_path_unaffected(tmp_path, drop_flag):
    """No BAMs → no _effective_chimera_ok column → the filter is inert for BOTH
    drop_chimeric_winners values (the read survives either way)."""
    mm2 = tmp_path / 'minimap2.tsv'
    des = tmp_path / 'deSALT.tsv'
    _write(mm2, [_row('r1', 'chrI', 1000)])
    _write(des, [_row('r1', 'chrI', 1000)])
    out = tmp_path / 'merged.tsv'
    merge_corrected_tsvs(
        per_aligner_tsvs={'minimap2': mm2, 'deSALT': des},
        output_tsv=out,
        per_aligner_corrected_bams=None,
        drop_chimeric_winners=drop_flag,
    )
    merged = pd.read_csv(out, sep='\t')
    assert set(merged['read_id']) == {'r1'}


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
