"""ISSUE-089: the 5' rescue's destination guard accepts canonical motifs from either strand (opt-in fix).

`_revert_selfinflicted_noncanonical_n` judged a writer-drawn N with `is_canonical_junction`, which accepts all six
forward-genome pairs, so a CT-AC (the minus strand's GT-AG) drawn on a + strand RNA passed (Reviewer B: RS-28 CT-GC,
RS-29 CT-AC in the re-split census). With RECTIFY_5P_STRANDED_MOTIF=1 and a known RNA strand, only that strand's
motifs count; unset, behaviour is unchanged.
"""
import pysam

from rectify.core.bam.bam_writer import _revert_selfinflicted_noncanonical_n

CHROM = 'chrM89'


def _case():
    g = list('C' * 600)
    g[120:122] = 'CT'          # intron 120..220: forward-genome CT ... AC = a minus-strand GT-AG
    g[218:220] = 'AC'
    genome = {CHROM: ''.join(g)}
    header = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [{'SN': CHROM, 'LN': 600}]})
    r = pysam.AlignedSegment(header)
    r.query_name = 'r'
    r.reference_id = 0
    r.flag = 0
    r.reference_start = 100
    r.cigarstring = '20M100N30M'           # the writer's surgery drew the N
    r.query_sequence = 'C' * 50
    r.query_qualities = pysam.qualitystring_to_array('I' * 50)
    pre_cigar = [(0, 50)]
    return r, genome, pre_cigar


def test_default_keeps_a_wrong_strand_pair(monkeypatch):
    monkeypatch.delenv('RECTIFY_5P_STRANDED_MOTIF', raising=False)
    r, genome, pre = _case()
    assert _revert_selfinflicted_noncanonical_n(r, genome, (), pre, 100, rna_strand='+') is False
    assert r.cigarstring == '20M100N30M'


def test_switch_reverts_a_minus_strand_motif_on_a_plus_rna(monkeypatch):
    monkeypatch.setenv('RECTIFY_5P_STRANDED_MOTIF', '1')
    r, genome, pre = _case()
    assert _revert_selfinflicted_noncanonical_n(r, genome, (), pre, 100, rna_strand='+') is True
    assert r.cigarstring == '50M' and r.reference_start == 100


def test_switch_keeps_it_on_a_minus_rna_and_without_a_strand(monkeypatch):
    monkeypatch.setenv('RECTIFY_5P_STRANDED_MOTIF', '1')
    for strand in ('-', None):
        r, genome, pre = _case()
        assert _revert_selfinflicted_noncanonical_n(r, genome, (), pre, 100, rna_strand=strand) is False, strand
