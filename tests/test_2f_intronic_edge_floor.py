"""ISSUE-017 — the informative floor must cover the truncated comparison.

Mechanism (confirmed by replaying the tester's six chr5 reads with their real
records and candidates): when the alignment starts INSIDE the intron
(``dist == 0``) the rescue loop compares only the intron-mapped bases. A read
1–2 bases inside an annotated intron therefore compares a 1–2-mer, which finds
an ED-0 hit somewhere in the ±15 shift × 11 offset sweep; ``ed_exon == 0``
accepts it, and among ED-0 ties the smallest-|shift| CANONICAL donor wins —
the GTRAGT +5 GT, i.e. a junction 4 nt into the intron, on both strands
(106/160 near-annotated added_nov on the cohort; the annotated exon fit the
clip strictly better in 82/119). The ISSUE-006 floor covered ``five_clip`` and
``len(rescue_seq)`` but never the truncated ``_rlen``.

Fix: a truncated comparison shorter than ``min_informative_clip_bp()`` is not
a sequence search — skip it for that candidate and let the structural Case-4
intronic snap decide (annotated boundary, or "favours intron" = leave the
read alone). The toy intron starts with the GTRAGT consensus so the +4 decoy
exists exactly as on the human introns.
"""
import pysam
import pytest

from rectify.core.splice.overhang_informativeness import COUNTERS, reset_counters
from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation

EXON1_TAIL = 'ACGTTGCATGCAGTCCATG'                       # 19 nt, ends ...CCATG
_GENOME_SEQ = (
    ('T' * (40 - len(EXON1_TAIL) - 1)) + 'A' + EXON1_TAIL  # exon1  [0, 40)
    + 'GTAAGT' + 'N' * 92 + 'AG'                           # intron [40, 140): GTRAGT consensus, decoy GT at 44
    + 'C' * 100                                            # exon2  [140, 240)
)
GENOME = {'chrT': _GENOME_SEQ}
JUNCTION = ('chrT', 40, 140)
DECOY = ('chrT', 44, 140)


@pytest.fixture(autouse=True)
def _fresh(monkeypatch):
    reset_counters()
    monkeypatch.delenv('RECTIFY_2F_NOVEL_GATE', raising=False)
    yield
    reset_counters()


def test_the_locus_is_what_the_test_assumes():
    assert _GENOME_SEQ[40:46] == 'GTAAGT' and _GENOME_SEQ[44:46] == 'GT'
    assert _GENOME_SEQ[138:140] == 'AG'


def _read(start, seq, name):
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'},
                                           'SQ': [{'SN': 'chrT', 'LN': 3_000_000}]})
    r = pysam.AlignedSegment(hdr)
    r.query_name = name
    r.reference_name = 'chrT'
    r.reference_start = start
    r.cigartuples = [(0, len(seq))]
    r.mapping_quality = 60
    r.query_sequence = seq
    return r


def _rescue(read, **kw):
    return rescue_3ss_truncation(read, GENOME, {JUNCTION}, '+',
                                 annotated_junctions={JUNCTION}, **kw)


def test_two_error_bases_never_reach_the_decoy():
    """Alignment starts 2 bases inside the intron with two error bases TA — a
    2-mer absent from the exon-1 tail (no ED-0 at shift 0) but present at
    genome[41:43], which the +4 sweep window reaches: before the fix that ED-0
    hit at shift +4 won and the read was spliced to the decoy donor at 44.
    Now the 2-mer is below the floor, the sequence search is skipped for the
    candidate, and Case 4 snaps the read to the ANNOTATED boundary."""
    assert 'TA' not in EXON1_TAIL[-12:] and _GENOME_SEQ[41:43] == 'TA'
    read = _read(42, 'TA' + _GENOME_SEQ[44:140] + 'C' * 60, 'err_ta')
    res = _rescue(read)
    assert res['rescued_junction'] != DECOY
    assert res['rescued'] and res['rescue_type'] == 'intronic_snap'
    assert res['rescued_junction'] == JUNCTION


def test_two_intronic_bases_that_match_the_intron_are_left_alone():
    """The bases ARE the intron (AA = genome[42:44]): nothing is imperfect, so
    there is no rescue sequence and no junction is created either way."""
    read = _read(42, 'AA' + _GENOME_SEQ[44:140] + 'C' * 60, 'retained2')
    res = _rescue(read)
    assert res['rescued_junction'] != DECOY
    assert not res['rescued']


def test_two_error_bases_inside_an_annotated_intron_snap_to_the_annotated_boundary():
    """Alignment starts 2 bases inside the intron with error bases (CC): the
    2-mer search is skipped; Case 4 cannot say the bases favour the intron, so
    the read is snapped to the ANNOTATED boundary — never to the +4 decoy."""
    read = _read(42, 'CC' + _GENOME_SEQ[44:140] + 'C' * 60, 'err2')
    res = _rescue(read)
    assert res['rescued'] and res['rescue_type'] == 'intronic_snap'
    assert res['rescued_junction'] == JUNCTION
    assert res['landing_annotated'] is True


def test_a_deep_intronic_run_still_uses_the_sequence_search():
    """15 exon-1 bases mapped into the intron carry the floor: the sequence
    rescue runs and lands on the annotated junction."""
    read = _read(140 - 15, _GENOME_SEQ[40 - 15:40] + 'C' * 60, 'deep15')
    res = _rescue(read)
    assert res['rescued'] and res['rescued_junction'] == JUNCTION
    assert COUNTERS.get('intronic_edge_below_floor', 0) == 0


def test_minus_strand_mirror():
    """Same shape on the minus strand: the intron-mapped bases sit at the read's
    right end. A 2-base intronic edge never reaches the shift sweep."""
    # Build a minus-strand locus by reverse-complementing the toy genome.
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    rc = ''.join(comp[b] for b in reversed(_GENOME_SEQ))
    n = len(_GENOME_SEQ)
    genome = {'chrT': rc}
    junction = ('chrT', n - 140, n - 40)          # the same intron on the reverse strand
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'},
                                           'SQ': [{'SN': 'chrT', 'LN': 3_000_000}]})
    r = pysam.AlignedSegment(hdr)
    r.query_name = 'minus_err2'
    r.reference_name = 'chrT'
    r.is_reverse = True
    r.mapping_quality = 60
    # exon-2 body (60 nt) then the 96 intron bases then 2 error bases, ending 2 bases inside the intron
    body = 'G' * 60 + rc[n - 140:n - 42] + 'GG'
    r.reference_start = n - 140 - 60
    r.cigartuples = [(0, len(body))]
    r.query_sequence = body
    res = rescue_3ss_truncation(r, genome, {junction}, '-', annotated_junctions={junction})
    assert res.get('rescued_junction') in (junction, None)
    assert res.get('rescued_junction') != ('chrT', n - 140, n - 44)
