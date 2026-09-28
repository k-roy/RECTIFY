"""ISSUE-083 re-split (RECTIFY_2F_RESPLIT): the 5' segment is aligned ACROSS the junction, and the writer draws
the exon-2 head from ``five_prime_exon2_cigar`` instead of a flat ``kM``.

Two synthetic genomes, one per strand, built so that the read's true split puts the last two exon-1 bases where
the aligner had clipped them into an exon-2 prefix: with the fixed split those two bases land on exon 2 and the
exon-1 block needs a deletion beside the N to close (cards 083-2/3/4); the re-split puts them back on exon 1,
the exon-1 block is clean, and the exon-2 head matches the reference. The writer tests pin the surgery's
invariants: the query is conserved, the N runs exactly from the reported donor to the reported acceptor, the
body resumes at the anchor with every downstream base where the aligner had it, and a head that cannot be cut
from the body refuses (read untouched).
"""
import random
from array import array

import pysam
import pytest

from rectify.core.align.local_aligner import (
    _OP_D, _OP_I, _OP_M, _OP_S,
    align_clip_to_exon, cigar_ops_to_str, resplit_across_junction,
)
from rectify.core.bam.read_edits import (
    extend_read_5prime_for_junction_rescue, projected_5prime_rescue_intron_edge,
)

random.seed(83)


def _rand(n):
    return ''.join(random.choice('ACGT') for _ in range(n))


def _rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def _qspan(ops):
    return sum(n for o, n in ops if o in (_OP_M, _OP_I, _OP_S, 7, 8))


def _rspan(ops):
    return sum(n for o, n in ops if o in (_OP_M, _OP_D, 7, 8))


def _genome():
    """exon1 [0,60) · intron [60,260) GT…AG · exon2 [260,320) · tail. The read carries exon1[-20:] + exon2."""
    exon1 = _rand(58) + 'CA'                    # exon 1 ends ...CA
    intron = 'GT' + _rand(196) + 'AG'
    exon2 = 'GG' + _rand(58)                    # exon 2 starts GG
    tail = _rand(40)
    g = exon1 + intron + exon2 + tail
    return g, 60, 260


def _plus_read(g, intron_start, intron_end):
    """Plus strand: read = exon1[-20:] + exon2[:50]. The aligner clipped the 20 exon-1 bases AND the first 2
    exon-2 bases (clip 22), so the body starts 2 bases into exon 2; 2F's fixed split gives k=2."""
    seq = g[intron_start - 20:intron_start] + g[intron_end:intron_end + 50]
    header = pysam.AlignmentHeader.from_references(['chrT'], [len(g)])
    r = pysam.AlignedSegment(header)
    r.query_name = 'plus'
    r.reference_id = 0
    r.flag = 0
    r.reference_start = intron_end + 2
    r.cigartuples = [(4, 22), (0, 48)]
    r.query_sequence = seq
    r.query_qualities = array('B', [30] * len(seq))
    return r


def _minus_read(g, intron_start, intron_end):
    """Minus strand (BAM orientation): the 5' end is the trailing S. Body = exon2 tail-side... mirror of the
    plus case: the record covers exon2' = g[intron_start-50:intron_start] (the read's exon 2 in genomic order,
    lower coordinates) then the clip = the last 2 of it + exon1' = g[intron_end:intron_end+20]."""
    seq = g[intron_start - 50:intron_start] + g[intron_end:intron_end + 20]
    header = pysam.AlignmentHeader.from_references(['chrT'], [len(g)])
    r = pysam.AlignedSegment(header)
    r.query_name = 'minus'
    r.reference_id = 0
    r.flag = 16
    r.reference_start = intron_start - 50
    r.cigartuples = [(0, 48), (4, 22)]
    r.query_sequence = seq
    r.query_qualities = array('B', [30] * len(seq))
    return r


def _minus_genome():
    """For the minus case the donor/acceptor sit the other way round in genomic order: exon2' [0,60) · intron
    [60,260) · exon1' [260,320). The read's exon 1 (its 5' end) is at the HIGHER coordinate."""
    return _genome()


@pytest.mark.parametrize('strand', ['+', '-'])
def test_resplit_puts_the_clipped_exon2_prefix_where_it_belongs(strand):
    if strand == '+':
        g, s, e = _genome()
        read = _plus_read(g, s, e)
    else:
        g, s, e = _minus_genome()
        read = _minus_read(g, s, e)
    clip = 22
    # the fixed split: 2F aligns clip minus the k=2 exon-2 prefix to exon 1 (20 bases)
    align_seq = read.query_sequence[:20] if strand == '+' else read.query_sequence[-20:]
    fixed_ops, _ = align_clip_to_exon(align_seq, g, s, e, strand)
    assert fixed_ops == [(_OP_M, 20)], fixed_ops
    rs = resplit_across_junction(read, clip, len(align_seq), g, s, e, strand)
    assert rs is not None
    assert rs['exon1_ops'] == [(_OP_M, 20)], rs['exon1_ops']
    assert all(o == _OP_M for o, n in rs['exon2_ops']), rs['exon2_ops']
    assert _qspan(rs['exon1_ops']) + _qspan(rs['exon2_ops']) == rs['anchor_query'] if strand == '+' else True
    assert rs['body_replaced_query'] == _qspan(rs['exon1_ops']) + _qspan(rs['exon2_ops']) - clip
    assert rs['body_replaced_query'] >= 12                     # the anchor sits past `extra` body bases
    assert _rspan(rs['exon2_ops']) == rs['body_replaced_ref'] + 2   # the 2 clip bases over exon 2 + the head
    assert len(rs['exon1_seq']) == 20 and rs['exon1_seq'] == (g[s - 20:s] if strand == '+' else g[e:e + 20])


@pytest.mark.parametrize('strand', ['+', '-'])
def test_writer_draws_the_exon2_head_and_conserves_the_read(strand):
    if strand == '+':
        g, s, e = _genome()
        read = _plus_read(g, s, e)
    else:
        g, s, e = _minus_genome()
        read = _minus_read(g, s, e)
    before = read.__copy__()
    rs = resplit_across_junction(read, 22, 20, g, s, e, strand)
    fpp = s - 1 if strand == '+' else e
    edge = projected_5prime_rescue_intron_edge(read, 22, strand, 0, exon2_prefix=2,
                                               exon_cigar_str=cigar_ops_to_str(rs['exon1_ops']),
                                               exon2_cigar_str=cigar_ops_to_str(rs['exon2_ops']))
    assert edge == (e if strand == '+' else s), edge
    assert extend_read_5prime_for_junction_rescue(
        read, fpp, 22, strand, exon_cigar_str=cigar_ops_to_str(rs['exon1_ops']), exon2_prefix=2,
        exon2_cigar_str=cigar_ops_to_str(rs['exon2_ops']))
    assert read.infer_query_length() == before.infer_query_length()
    assert read.query_sequence == before.query_sequence
    ns = []
    pos = read.reference_start
    for o, n in read.cigartuples:
        if o == 3:
            ns.append((pos, pos + n))
        if o in (0, 2, 3, 7, 8):
            pos += n
    assert ns == [(s, e)], (ns, read.cigarstring)
    assert 4 not in {o for o, n in read.cigartuples}            # no clip left: every base is placed
    # every base the aligner had placed keeps its reference position (the body is untouched past the anchor)
    original = dict(before.get_aligned_pairs(matches_only=True))
    for q, r in read.get_aligned_pairs(matches_only=True):
        if q in original:
            assert original[q] == r, (q, r, original[q])
    # and the placed bases are the reference's
    for q, r in read.get_aligned_pairs(matches_only=True):
        assert read.query_sequence[q] == g[r], (q, r)


def test_writer_refuses_a_head_the_body_cannot_give_up():
    g, s, e = _genome()
    read = _plus_read(g, s, e)
    before = read.to_string()
    # an exon-2 head claiming 60 query bases from a 48-base body: the trim fails, the read is untouched
    assert not extend_read_5prime_for_junction_rescue(read, s - 1, 22, '+', exon_cigar_str='20M', exon2_prefix=2,
                                                      exon2_cigar_str='62M')
    assert read.to_string() == before
    assert projected_5prime_rescue_intron_edge(read, 22, '+', 0, exon2_prefix=2, exon_cigar_str='20M',
                                               exon2_cigar_str='62M') is None


def test_no_exon2_cigar_keeps_the_legacy_kM_path():
    g, s, e = _genome()
    read = _plus_read(g, s, e)
    assert extend_read_5prime_for_junction_rescue(read, s - 1, 22, '+', exon_cigar_str='20M', exon2_prefix=2)
    assert read.cigartuples == [(0, 20), (3, e - s), (0, 2), (0, 48)] or read.cigartuples == [(0, 20), (3, e - s), (0, 50)]
