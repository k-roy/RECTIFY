"""Tests for FASTQ-comment aux-tag preservation through alignment and consensus.

`rectify trim-cdna-polya` writes `ro:A:` (orientation) and `pl:i:` (poly-A tail length) into the
FASTQ COMMENT. Only minimap2 is invoked with `-y`; the other four aligners do not propagate
comments, and `build_chimeric_read()` rebuilt every record from scratch. Measured on real reads
(planning/561, 563): survival was 40.94% on the 3-aligner default and **0.00%** under chimeric
consensus -- the `run-all` default -- which silently corrupted per-read STRAND (via `ro`), not just
tail length.

Validated by a 60-config matrix (planning/574): all 5 aligners, every subset, both chimeric modes,
plus a reversed-order control -- 60/60 at 100.00%.
"""
import pysam
import pytest


# Aux tags must survive the chimeric consensus
# --------------------------------------------------------------------------

def _template(hdr, tags):
    r = pysam.AlignedSegment(hdr)
    r.query_name = 'r1'
    r.query_sequence = 'A' * 50
    r.query_qualities = pysam.qualitystring_to_array('I' * 50)
    r.flag = 0
    r.reference_id = 0
    r.reference_start = 100
    r.cigar = [(0, 50)]
    r.mapping_quality = 60
    for name, val, typ in tags:
        r.set_tag(name, val, value_type=typ)
    return r


class _CR:
    chimeric_ref_start = 100
    chimeric_cigar = [(0, 50)]
    segment_winners = [(0, 'minimap2', 0, 50)]
    confidence = 'high'
    is_chimeric = False
    n_segments = 1
    n_aligners_used = 1
    five_prime_aligner = 'minimap2'
    three_prime_aligner = 'minimap2'
    interior_aligners = []


def _build(tags):
    from rectify.core.consensus.chimeric_consensus import build_chimeric_read
    hdr = pysam.AlignmentHeader.from_dict(
        {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrI', 'LN': 1000}]})
    t = _template(hdr, tags)
    return build_chimeric_read(template_read=t, ref_start=100,
                               cigar_tuples=[(0, 50)], chimeric_result=_CR(), header=hdr)


def test_cdna_tags_survive_chimeric_consensus():
    """ro/pl/XO/XY/XU carry the ENTIRE ONT-cDNA protocol; losing them voids it.

    Before this fix build_chimeric_read created a fresh AlignedSegment and copied
    no tags at all, so every per-read annotation died at the consensus step even
    for the one aligner that propagates FASTQ comments.
    """
    out = _build([('ro', 'S', 'A'), ('pl', 22, 'i'), ('XO', 'fwd', 'Z'),
                  ('XU', 'ACGTACGT', 'Z'), ('XY', 't1', 'Z')])
    assert out.get_tag('ro') == 'S'
    assert out.get_tag('pl') == 22
    assert out.get_tag('XO') == 'fwd'
    assert out.get_tag('XU') == 'ACGTACGT'
    assert out.get_tag('XY') == 't1'


def test_pl_zero_survives_and_is_not_dropped():
    """`pl:i:0` must round-trip as 0, never vanish.

    Absent means "fall back to the post-alignment measurement"; 0 means "genuinely
    no tail". Conflating them is how a read acquires a confident, plausible, wrong
    tail length.
    """
    out = _build([('pl', 0, 'i'), ('ro', 'U', 'A')])
    assert out.has_tag('pl')
    assert out.get_tag('pl') == 0


def test_positional_tags_are_dropped_as_stale():
    """NM/AS/MD describe the TEMPLATE's alignment; the chimeric CIGAR differs."""
    out = _build([('ro', 'S', 'A'), ('NM', 3, 'i'), ('AS', 99, 'i'), ('MD', '50', 'Z')])
    assert out.get_tag('ro') == 'S'
    for stale in ('NM', 'AS', 'MD'):
        assert not out.has_tag(stale), "%s is stale against the chimeric CIGAR" % stale


def test_chimeric_metadata_tags_still_written():
    """The fix must not displace build_chimeric_read's own provenance tags."""
    out = _build([('ro', 'S', 'A')])
    for t in ('Xa', 'Xc', 'Xz', 'Xg', 'Xm'):
        assert out.has_tag(t)
