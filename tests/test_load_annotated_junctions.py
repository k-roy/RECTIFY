#!/usr/bin/env python3
"""Tests for load_annotated_junctions, incl. the exon-derivation fallback.

The fallback fixes a latent bug: an exon-only GTF (GENCODE / Ensembl, and the
SGNEx A549 annotation) carries no explicit `intron` records, so the loader
returned an empty set and silently disabled junction-guided scoring (the
annotated +8 reward in the chimeric score_segment path went dead).

Backward-compat guarantee under test: when ANY `intron` record is present
(e.g. the yeast SGD GTF/GFF), the intron set is returned unchanged and exon
derivation never runs — so intron-bearing annotations stay byte-identical.

Author: Kevin R. Roy
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rectify.core.consensus.consensus import load_annotated_junctions


def _gtf(path, lines):
    with open(path, 'w') as f:
        for ln in lines:
            f.write('\t'.join(str(x) for x in ln) + '\n')


def _exon(chrom, s, e, tid, strand='+'):
    return [chrom, 'src', 'exon', s, e, '.', strand, '.',
            f'transcript_id "{tid}"; gene_id "g_{tid}";']


def _intron(chrom, s, e, strand='+'):
    return [chrom, 'src', 'intron', s, e, '.', strand, '.', 'transcript_id "t";']


# ── backward-compat: intron records present → exon path NOT used ────────────

def test_intron_records_take_priority(tmp_path):
    p = tmp_path / 'a.gtf'
    # intron at 1-based 101..200; plus an exon pair that would derive a
    # DIFFERENT intron — must be ignored because intron records exist.
    _gtf(p, [
        _intron('chrI', 101, 200),
        _exon('chrI', 1, 50, 't2'),
        _exon('chrI', 500, 600, 't2'),
    ])
    js = load_annotated_junctions(str(p))
    # only the explicit intron: (0-based start, exclusive end) = (100, 200)
    assert js == {('chrI', 100, 200, '+')}


# ── exon-only fallback derives introns ──────────────────────────────────────

def test_exon_only_derives_introns(tmp_path):
    p = tmp_path / 'b.gtf'
    # transcript with 3 exons → 2 introns
    _gtf(p, [
        _exon('chrI', 1, 100, 't1'),
        _exon('chrI', 201, 300, 't1'),
        _exon('chrI', 401, 500, 't1'),
    ])
    js = load_annotated_junctions(str(p))
    # intron1: exon1 ends 100 → first intronic base 0-based = 100;
    #          exon2 starts 201 → exclusive end = 200
    # intron2: 300 .. 400
    assert js == {('chrI', 100, 200, '+'), ('chrI', 300, 400, '+')}


def test_derived_matches_equivalent_intron_record(tmp_path):
    """An exon-derived intron must equal what the explicit intron record path
    produces for the same gap — the convention is identical."""
    # exon ends 1-based 100, next starts 1-based 201 → intron spans 1-based
    # 101..200.  As an explicit record that is parts[3]=101, parts[4]=200.
    pe = tmp_path / 'exon.gtf'
    _gtf(pe, [_exon('chrI', 1, 100, 't1'), _exon('chrI', 201, 300, 't1')])
    pi = tmp_path / 'intron.gtf'
    _gtf(pi, [_intron('chrI', 101, 200)])
    assert load_annotated_junctions(str(pe)) == load_annotated_junctions(str(pi))


def test_exon_strand_preserved(tmp_path):
    p = tmp_path / 'c.gtf'
    _gtf(p, [
        _exon('chrI', 1, 100, 't1', strand='-'),
        _exon('chrI', 201, 300, 't1', strand='-'),
    ])
    assert load_annotated_junctions(str(p)) == {('chrI', 100, 200, '-')}


def test_single_exon_transcript_yields_no_intron(tmp_path):
    p = tmp_path / 'd.gtf'
    _gtf(p, [_exon('chrI', 1, 100, 't1')])
    assert load_annotated_junctions(str(p)) == set()


def test_empty_annotation_returns_empty(tmp_path):
    p = tmp_path / 'e.gtf'
    _gtf(p, [['chrI', 'src', 'gene', 1, 100, '.', '+', '.', 'gene_id "g";']])
    assert load_annotated_junctions(str(p)) == set()


if __name__ == '__main__':
    import pytest
    sys.exit(pytest.main([__file__, '-q']))
