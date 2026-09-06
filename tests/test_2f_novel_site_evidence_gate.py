"""Module 2F novel-site evidence gate (2026-09-05, Sumner cohort; ISSUE-017).

A rescue onto an ANNOTATED junction rests on the annotation's prior and keeps
its historical acceptance (a few matched bases suffice — a read ending a few
bases inside an annotated intron is snapped, as it always was). Onto an
UNANNOTATED candidate (a pool junction, the read's own N-op) the placed exon
segment is judged: ``min_informative_clip_bp()`` matched bases, no gap at the
junction-side end, a bounded indel burden, and a CIGAR at all (fail closed).

What the verdict DOES is a mode (RECTIFY_2F_NOVEL_GATE / `rectify correct
--2f-novel-gate`). ``report`` (the default, arbiter RULING 2): the rescue is
drawn and the verdict is recorded in five_prime_novel_evidence (``pass`` or
the token). ``refuse``: the sequence/snap rescue is refused and the token is
also the five_prime_rescue_refused value. Either way every rescue carries
five_prime_landing_annotated (0/1), so provenance is an offline TSV join.
Both rescue paths that place a segment are judged: the sequence rescue (Cases
1/2 and the terminal peel) and the Case-4 intronic snap.
``annotated_junctions=None`` (a caller without an annotation) is the legacy
behavior: nothing is novel.
"""
import pysam
import pytest

import rectify.core.bam.bam_processor as bp
from rectify.core.bam.output import CORRECTION_TSV_HEADER, correction_result_to_tsv_row
from rectify.core.splice.overhang_informativeness import COUNTERS, reset_counters
from rectify.core.splice.splice_aware_5prime import (
    NOVEL_EXON_REFUSALS,
    NOVEL_GATE_DEFAULT,
    _novel_exon_evidence_refusal,
    min_informative_clip_bp,
    novel_gate_mode,
    rescue_3ss_truncation,
)

EXON1_TAIL = 'ACGTTGCATGCAGTCCATG'          # 19 nt, aperiodic
_GENOME_SEQ = (
    ('T' * (40 - len(EXON1_TAIL) - 1)) + 'A' + EXON1_TAIL   # exon1  [0, 40)
    + 'GT' + 'N' * 96 + 'AG'                                # intron [40, 140)
    + 'C' * 100                                             # exon2  [140, 240)
)
GENOME = {'chrT': _GENOME_SEQ}
JUNCTION = ('chrT', 40, 140)
TOKEN = 'novel_exon_matched_below_floor'


@pytest.fixture(autouse=True)
def _fresh(monkeypatch):
    reset_counters()
    monkeypatch.delenv('RECTIFY_2F_NOVEL_GATE', raising=False)
    yield
    reset_counters()


@pytest.fixture
def refuse(monkeypatch):
    monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', 'refuse')


def _make_read(cigar, seq, start=140, strand='+', name='r'):
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'},
                                           'SQ': [{'SN': 'chrT', 'LN': 3_000_000}]})
    r = pysam.AlignedSegment(hdr)
    r.query_name = name
    r.reference_name = 'chrT'
    r.reference_start = start
    r.cigartuples = cigar
    r.is_reverse = (strand == '-')
    r.mapping_quality = 60
    r.query_sequence = seq
    return r


def _clip_read(clip_len):
    """5' soft clip = the last `clip_len` bases of exon 1 (a perfect rescue)."""
    clip = _GENOME_SEQ[40 - clip_len:40]
    return _make_read([(4, clip_len), (0, 60)], clip + 'C' * 60, name=f'clip{clip_len}')


def _intronic_read(depth):
    """5' end `depth` bases INSIDE the intron, those bases being exon 1's tail:
    the Case-4 intronic-snap shape (the tester's 4M / 3M created junctions)."""
    seq = _GENOME_SEQ[40 - depth:40] + 'C' * 60
    return _make_read([(0, depth + 60)], seq, start=140 - depth, name=f'in{depth}')


def _novel(read, **kw):
    return rescue_3ss_truncation(read, GENOME, {JUNCTION}, '+', annotated_junctions=set(), **kw)


def _annotated(read, **kw):
    return rescue_3ss_truncation(read, GENOME, {JUNCTION}, '+', annotated_junctions={JUNCTION}, **kw)


# ---------------------------------------------------------------------------
# The verdict itself
# ---------------------------------------------------------------------------

def test_refusal_tokens():
    assert min_informative_clip_bp() == 10
    assert _novel_exon_evidence_refusal([(0, 12)], 12, '+', 0.2) == ''
    assert _novel_exon_evidence_refusal([(0, 4)], 4, '+', 0.2) == TOKEN
    # 2M8I2M2I1M1D1M (tester sample): 5 matched, 11 indel bp
    assert _novel_exon_evidence_refusal(
        [(0, 2), (1, 8), (0, 2), (1, 2), (0, 1), (2, 1), (0, 1)], 14, '+', 0.2) == TOKEN
    # enough matched bases but a gap on the junction side: + strand = last op
    assert _novel_exon_evidence_refusal([(0, 12), (1, 2)], 14, '+', 0.2) == 'novel_exon_gap_at_junction'
    assert _novel_exon_evidence_refusal([(1, 2), (0, 12)], 14, '+', 0.2) == ''       # free end
    # - strand = first op abuts the junction
    assert _novel_exon_evidence_refusal([(1, 2), (0, 12)], 14, '-', 0.2) == 'novel_exon_gap_at_junction'
    assert _novel_exon_evidence_refusal([(0, 12), (1, 2)], 14, '-', 0.2) == ''
    # indel burden: 12 matched allows max(5, ceil(0.2*12)=3) = 5 bp of I+D
    assert _novel_exon_evidence_refusal([(0, 6), (1, 5), (0, 6)], 17, '+', 0.2) == ''
    assert _novel_exon_evidence_refusal([(0, 6), (1, 6), (0, 6)], 18, '+', 0.2) == 'novel_exon_indel_burden'
    # the allowance scales with the block
    assert _novel_exon_evidence_refusal([(0, 50), (2, 9), (0, 50)], 100, '+', 0.2) == ''
    # no CIGAR: fail CLOSED (arbiter R1a-1), whatever the placed length
    assert _novel_exon_evidence_refusal(None, 9, '+', 0.2) == 'novel_exon_no_cigar'
    assert _novel_exon_evidence_refusal([], 40, '+', 0.2) == 'novel_exon_no_cigar'
    assert set(NOVEL_EXON_REFUSALS) == {TOKEN, 'novel_exon_gap_at_junction',
                                        'novel_exon_indel_burden', 'novel_exon_no_cigar'}


def test_mode_defaults_to_report_and_rejects_garbage(monkeypatch):
    assert NOVEL_GATE_DEFAULT == 'report' and novel_gate_mode() == 'report'
    monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', 'refuse')
    assert novel_gate_mode() == 'refuse'
    monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', 'whatever')
    assert novel_gate_mode() == 'report'


# ---------------------------------------------------------------------------
# Case-4 intronic snap (the shape of the tiny created junctions)
# ---------------------------------------------------------------------------

def test_annotated_intronic_snap_keeps_the_historical_acceptance():
    res = _annotated(_intronic_read(4))
    assert res['rescued'] and res['rescue_type'] == 'intronic_snap'
    assert res['rescued_junction'] == JUNCTION
    assert res['five_prime_exon_cigar'] == '4M'
    assert res['landing_annotated'] is True and res['novel_evidence'] == ''


def test_legacy_call_without_provenance_is_unchanged():
    res = rescue_3ss_truncation(_intronic_read(4), GENOME, {JUNCTION}, '+')
    assert res['rescued'] and res['rescue_type'] == 'intronic_snap'
    assert res['landing_annotated'] is True and res['novel_evidence'] == ''


def test_report_mode_draws_the_novel_snap_and_records_the_verdict():
    res = _novel(_intronic_read(4))
    assert res['rescued'] and res['rescue_type'] == 'intronic_snap'
    assert res['landing_annotated'] is False
    assert res['novel_evidence'] == TOKEN
    assert not res.get('clip_refused')
    assert COUNTERS[TOKEN] == 1                       # once per read, not per peel depth


def test_refuse_mode_refuses_the_novel_four_base_snap(refuse):
    res = _novel(_intronic_read(4))
    assert not res['rescued']
    assert res['rescue_type'] in ('proximity', 'none')
    assert res.get('clip_refused') == TOKEN
    assert COUNTERS[TOKEN] == 1


def test_novel_landing_site_passes_an_informative_intronic_segment(refuse):
    """15 exon-1 bases inside the intron carry the floor: rescued on a novel
    site even in refuse mode (the reanchor pre-pass turns this shape into a
    soft clip and the sequence rescue takes it; the aligner spells it 14M1I1M —
    15 matched, one insertion inside the allowance, no gap at the junction)."""
    res = _novel(_intronic_read(15))
    assert res['rescued'] and res['rescue_type'] in ('softclip', 'intronic_snap')
    assert res['rescued_junction'] == JUNCTION
    from rectify.core.align.local_aligner import cigar_str_to_ops
    ops = cigar_str_to_ops(res['five_prime_exon_cigar'])
    assert sum(ln for op, ln in ops if op in (0, 7, 8)) >= 15
    assert res['landing_annotated'] is False and res['novel_evidence'] == 'pass'
    assert not res.get('clip_refused')


# ---------------------------------------------------------------------------
# Sequence rescue (Cases 1/2)
# ---------------------------------------------------------------------------

def test_novel_sequence_rescue_passes_an_informative_clip(refuse):
    res = _novel(_clip_read(15))
    assert res['rescued'] and res['rescue_type'] == 'softclip'
    assert res['rescued_junction'] == JUNCTION
    assert res['five_prime_exon_cigar'] == '15M'
    assert res['landing_annotated'] is False and res['novel_evidence'] == 'pass'


def test_novel_sequence_rescue_indel_riddled_placement(monkeypatch):
    """The acceptance test passed (a clean 15-nt clip); the placed exon CIGAR
    is what the gate reads. Force the aligner to spell it I-riddled."""
    def _bad(clip_seq, genome_seq, intron_start, intron_end, strand, **kw):
        return [(0, 3), (1, 4), (0, 1), (1, 3), (0, 1), (1, 3)], intron_start - 5
    import rectify.core.align.local_aligner as la
    monkeypatch.setattr(la, 'align_clip_to_exon', _bad)
    # report (default): the novel verdict alone would have drawn it with the token recorded; ISSUE-028
    # invariant E judges the block as a PLACEMENT decision in both modes (the forged spelling walks the clip
    # against the wrong reference bases: identity under 0.8, -21 bits) and refuses it with its own token
    from rectify.core.splice.splice_aware_5prime import EVIDENCE_REFUSALS
    res = _novel(_clip_read(15))
    assert not res['rescued'] and res.get('clip_refused') in EVIDENCE_REFUSALS, res
    # refuse: the novel verdict's more specific token precedes E's
    monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', 'refuse')
    res = _novel(_clip_read(15))
    assert not res['rescued'] and res.get('clip_refused') == TOKEN
    # ISSUE-026 invariant C-minimal (2026-09-05): the same spelling onto an
    # ANNOTATED site used to be accepted in either mode; the indel-burden SANITY
    # bound now applies to every landing (10 I on 5 M is more gap than half the
    # matches), as a placement decision in both modes, with its own token.
    from rectify.core.splice.splice_aware_5prime import ANNOTATED_INDEL_BURDEN_REFUSAL
    for mode in ('report', 'refuse'):
        monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', mode)
        res2 = _annotated(_clip_read(15))
        assert not res2['rescued'], (mode, res2)
        assert res2['clip_refused'] == ANNOTATED_INDEL_BURDEN_REFUSAL, (mode, res2)


def test_novel_sequence_rescue_gap_at_the_junction(refuse, monkeypatch):
    def _gap(clip_seq, genome_seq, intron_start, intron_end, strand, **kw):
        return [(0, 13), (1, 2)], intron_start - 13          # + strand: last op abuts
    import rectify.core.align.local_aligner as la
    monkeypatch.setattr(la, 'align_clip_to_exon', _gap)
    res = _novel(_clip_read(15))
    assert not res['rescued'] and res['clip_refused'] == 'novel_exon_gap_at_junction'


def test_novel_sequence_rescue_fails_closed_without_a_cigar(monkeypatch):
    def _boom(*a, **kw):
        raise RuntimeError('aligner down')
    import rectify.core.align.local_aligner as la
    monkeypatch.setattr(la, 'align_clip_to_exon', _boom)
    # report: drawn (as any aligner failure always was), verdict says no_cigar
    res = _novel(_clip_read(15))
    assert res['rescued'] and res['five_prime_exon_cigar'] == ''
    assert res['novel_evidence'] == 'novel_exon_no_cigar'
    # refuse: fail closed
    monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', 'refuse')
    res = _novel(_clip_read(15))
    assert not res['rescued'] and res['clip_refused'] == 'novel_exon_no_cigar'
    # Annotated site: still drawn (with an empty exon CIGAR, as before).
    res2 = _annotated(_clip_read(15))
    assert res2['rescued'] and res2['five_prime_exon_cigar'] == '' and res2['novel_evidence'] == ''


# ---------------------------------------------------------------------------
# Wiring through correct_read_3prime and the TSV
# ---------------------------------------------------------------------------

def test_correct_read_3prime_passes_the_annotated_set(monkeypatch):
    seen = {}

    def _spy(read, genome, cands, strand, **kw):
        seen['annotated'] = kw.get('annotated_junctions')
        return {'rescued': False}

    monkeypatch.setattr(bp, '_rescue_3ss', _spy)
    bp.correct_read_3prime(_clip_read(15), GENOME, apply_3ss_rescue=True,
                           annotated_junctions={JUNCTION})
    assert seen['annotated'] == {JUNCTION}
    # No annotation at all -> every candidate is novel, never "legacy".
    bp.correct_read_3prime(_clip_read(15), GENOME, apply_3ss_rescue=True,
                           annotated_junctions=None,
                           pool_chrom_index=bp._build_pool_chrom_index({JUNCTION}))
    assert seen['annotated'] == set()


def _row(read, **kw):
    return bp.correct_read_3prime(read, GENOME, apply_3ss_rescue=True, **kw)[0]


def _cells(row):
    return dict(zip(CORRECTION_TSV_HEADER, correction_result_to_tsv_row(row)))


def test_tsv_columns_report_mode():
    # ISSUE-026 invariant D appended `five_prime_exon2_prefix` after these two;
    # ISSUE-028 invariant E the two block-shape columns after that.
    assert CORRECTION_TSV_HEADER[-5:] == ['five_prime_landing_annotated',
                                          'five_prime_novel_evidence',
                                          'five_prime_exon2_prefix',
                                          'five_prime_exon_identity',
                                          'five_prime_exon_bits']
    pool = bp._build_pool_chrom_index({JUNCTION})
    novel_tiny = _row(_intronic_read(4), annotated_junctions=set(), pool_chrom_index=pool)
    assert novel_tiny['five_prime_rescued'] and novel_tiny['five_prime_rescue_refused'] == ''
    assert novel_tiny['five_prime_landing_annotated'] == 0
    assert novel_tiny['five_prime_novel_evidence'] == TOKEN
    assert _cells(novel_tiny)['five_prime_landing_annotated'] == '0'
    assert _cells(novel_tiny)['five_prime_novel_evidence'] == TOKEN
    novel_ok = _row(_clip_read(15), annotated_junctions=set(), pool_chrom_index=pool)
    assert novel_ok['five_prime_landing_annotated'] == 0
    assert novel_ok['five_prime_novel_evidence'] == 'pass'
    annotated = _row(_intronic_read(4), annotated_junctions={JUNCTION})
    assert annotated['five_prime_rescued'] and annotated['five_prime_landing_annotated'] == 1
    assert annotated['five_prime_novel_evidence'] == ''
    assert _cells(annotated)['five_prime_landing_annotated'] == '1'
    none = _row(_make_read([(0, 60)], 'C' * 60, start=200), annotated_junctions={JUNCTION})
    assert not none['five_prime_rescued']
    assert _cells(none)['five_prime_landing_annotated'] == '' and _cells(none)['five_prime_novel_evidence'] == ''


def test_tsv_columns_refuse_mode(refuse):
    pool = bp._build_pool_chrom_index({JUNCTION})
    refused = _row(_intronic_read(4), annotated_junctions=set(), pool_chrom_index=pool)
    assert not refused['five_prime_rescued']
    assert refused['five_prime_rescue_refused'] == TOKEN
    assert refused['five_prime_landing_annotated'] == 0           # the ATTEMPTED landing was novel
    assert refused['five_prime_novel_evidence'] == ''             # the token lives in the refusal column
    assert _cells(refused)['five_prime_landing_annotated'] == '0'


def test_refuse_then_rerescue_is_traced(refuse, monkeypatch):
    """Refuse mode, defect (c) from the tester's FAST on 34d6852: a read whose
    novel sequence rescue was refused and that a later structural path then
    rescued onto an ANNOTATED intron (hold-out 03c6013e) must say both things:
    landing 1, novel_evidence '<token>>annotated', refusal column empty, the
    token counted once. The body's result is stubbed; this pins the wrapper's
    bookkeeping (the on-read case is the tester's T1 check)."""
    import rectify.core.splice.splice_aware_5prime as sa

    def _body(*a, **kw):
        return {'rescued': True, 'rescue_type': 'intronic_snap', 'five_prime_corrected': 19,
                'rescued_junction': ('chrT', 20, 140), 'edit_distance': -1, 'query_bp': 2,
                'five_prime_exon_cigar': '2M', 'five_prime_upstream_trim': 0,
                'landing_annotated': True, 'novel_evidence': '', 'novel_refused_first': TOKEN}

    monkeypatch.setattr(sa, '_rescue_3ss_truncation_body', _body)
    res = rescue_3ss_truncation(_clip_read(15), GENOME, {JUNCTION}, '+',
                                annotated_junctions={JUNCTION}, terminal_peel=False)
    assert res['rescued'] and res['landing_annotated'] is True
    assert res['novel_evidence'] == TOKEN + '>annotated'
    assert not res.get('clip_refused') and 'novel_refused_first' not in res
    assert COUNTERS[TOKEN] == 1
