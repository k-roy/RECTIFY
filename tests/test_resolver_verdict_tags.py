"""Station A → Module 2F verdict tags (`XW`; Kevin 2026-09-05).

The resolver records what it decided about every assessed terminal clip —
refusals included — and 2F reuses the one sequence-only verdict (`low_info`),
still tries annotation-guided placement after every other refusal, and counts
the disagreements. See ``rectify/core/splice/resolver_verdict.py``.
"""

import random

import pysam
import pytest

import rectify.core.bam.bam_processor as bp
from rectify.core.align.overhang_resolver import (
    ResolverConfig,
    ResolverStats,
    resolve_clip,
    resolve_read,
)
from rectify.core.bam.processing_stats import ProcessingStats, write_stats_tsv
from rectify.core.splice.overhang_informativeness import COUNTERS, reset_counters
from rectify.core.splice.resolver_verdict import (
    SEQUENCE_ONLY_VERDICTS,
    VERDICT_REFUSALS,
    VERDICT_TAG,
    VERDICT_TOKENS,
    format_verdict_entry,
    parse_verdict_tag,
    verdict_for_clip,
)
from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation
from rectify.core.splice.splice_site_index import SpliceSiteIndex


@pytest.fixture(autouse=True)
def _fresh_counters():
    reset_counters()
    yield
    reset_counters()


# ---------------------------------------------------------------------------
# Tag format
# ---------------------------------------------------------------------------

def test_format_and_parse_round_trip():
    value = ';'.join([
        format_verdict_entry('L', 'low_info', 4.62, 0, 30),
        format_verdict_entry('R', 'rejected_edit', 41.0, 5000, 212),
        format_verdict_entry('L', 'repeat', None, 0, 90),      # later L entry wins
    ])
    assert value == 'L:low_info:4.6:0:30;R:rejected_edit:41.0:5000:212;L:repeat:-:0:90'
    parsed = parse_verdict_tag(value)
    assert parsed['R'] == {'token': 'rejected_edit', 'bits': 41.0, 'window': 5000, 'clip_len': 212}
    assert parsed['L'] == {'token': 'repeat', 'bits': None, 'window': 0, 'clip_len': 90}


def test_parse_is_advisory_and_skips_garbage():
    assert parse_verdict_tag(None) == {}
    assert parse_verdict_tag('') == {}
    assert parse_verdict_tag('L:low_info:x:0:30') == {}          # bad float
    assert parse_verdict_tag('M:low_info:1.0:0:30') == {}        # bad side
    assert parse_verdict_tag('L:nonsense:1.0:0:30') == {}        # unknown token
    assert parse_verdict_tag('L:low_info:1.0:0') == {}           # wrong arity
    assert parse_verdict_tag('junk;R:ambiguous:2.0:40:12') == {
        'R': {'token': 'ambiguous', 'bits': 2.0, 'window': 40, 'clip_len': 12}}


def test_vocabulary():
    assert 'resolved' in VERDICT_TOKENS and 'resolved' not in VERDICT_REFUSALS
    assert SEQUENCE_ONLY_VERDICTS == ('low_info',)
    assert VERDICT_TAG == 'XW'


# ---------------------------------------------------------------------------
# Resolver side
# ---------------------------------------------------------------------------

GLEN = 20_000
DON, ACC = 5_000, 5_800          # plus-strand intron [DON, ACC)


def _make_genome():
    rng = random.Random(7)
    seq = [rng.choice('ACGT') for _ in range(GLEN)]
    seq[DON:DON + 2] = list('GT')
    seq[ACC - 2:ACC] = list('AG')
    return ''.join(seq)


GENOME_SEQ = _make_genome()
GENOME = {'chrI': GENOME_SEQ}


@pytest.fixture(scope='module')
def index():
    return SpliceSiteIndex.build(GENOME)


@pytest.fixture()
def cfg():
    return ResolverConfig(alpha=0.01, max_intron=5000)


def _header():
    return pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrI', 'LN': GLEN}],
    })


def _read(name, query, cigar, ref_start, reverse=False):
    r = pysam.AlignedSegment(_header())
    r.query_name = name
    r.query_sequence = query
    r.flag = 16 if reverse else 0
    r.reference_id = 0
    r.reference_start = ref_start
    r.mapping_quality = 60
    r.cigartuples = cigar
    return r


def test_resolve_clip_reports_low_info_verdict(index, cfg):
    vd = {}
    stats = ResolverStats()
    assert resolve_clip(GENOME_SEQ, index, 'chrI', 'L', '+', 'A' * 30,
                        edge=ACC, cfg=cfg, stats=stats, verdict=vd) is None
    assert vd['token'] == 'low_info' and vd['window'] == 0
    assert vd['bits'] is not None and vd['bits'] < 20
    assert stats.refused_low_info == 1


def test_resolve_clip_reports_resolved_verdict(index, cfg):
    vd = {}
    clip = GENOME_SEQ[DON - 30:DON]
    placement = resolve_clip(GENOME_SEQ, index, 'chrI', 'L', '+', clip,
                             edge=ACC, cfg=cfg, stats=ResolverStats(), verdict=vd)
    assert placement is not None and placement.intron_start == DON
    assert vd['token'] == 'resolved' and vd['window'] > 0 and vd['bits'] > 20


def test_resolve_clip_below_min_clip_leaves_verdict_empty(index, cfg):
    vd = {}
    assert resolve_clip(GENOME_SEQ, index, 'chrI', 'L', '+', 'ACG',
                        edge=ACC, cfg=cfg, stats=ResolverStats(), verdict=vd) is None
    assert vd == {}


def test_resolve_read_writes_xw_on_refused_and_resolved_reads(index, cfg):
    stats = ResolverStats()
    # Refused (poly(A) left clip): tagged, alignment untouched, no XJ.
    body = GENOME_SEQ[ACC:ACC + 60]
    polya = _read('polyA', 'A' * 30 + body, [(4, 30), (0, 60)], ACC)
    assert resolve_read(polya, GENOME, index, cfg, stats) is False
    assert not polya.has_tag('XJ')
    entry = verdict_for_clip(polya, 'L', 30)
    assert entry['token'] == 'low_info' and entry['clip_len'] == 30
    assert polya.cigartuples == [(4, 30), (0, 60)]
    # Resolved (exon-1 tail): XJ as before, and XW says so with the ORIGINAL length.
    good = _read('good', GENOME_SEQ[DON - 30:DON] + body, [(4, 30), (0, 60)], ACC)
    assert resolve_read(good, GENOME, index, cfg, stats) is True
    assert good.has_tag('XJ')
    assert parse_verdict_tag(good.get_tag(VERDICT_TAG))['L']['token'] == 'resolved'
    assert stats.extra['verdict_tags_written'] == 2


def test_resolve_read_tags_both_sides_and_right_side_for_reverse_reads(index, cfg):
    body = GENOME_SEQ[ACC:ACC + 60]
    both = _read('both', 'A' * 30 + body + 'T' * 30,
                 [(4, 30), (0, 60), (4, 30)], ACC, reverse=True)
    resolve_read(both, GENOME, index, cfg, ResolverStats())
    parsed = parse_verdict_tag(both.get_tag(VERDICT_TAG))
    assert set(parsed) == {'L', 'R'}
    assert parsed['R']['token'] == 'low_info' and parsed['R']['clip_len'] == 30


def test_resolve_read_no_tag_without_assessed_clips(index, cfg):
    plain = _read('plain', GENOME_SEQ[ACC:ACC + 60], [(0, 60)], ACC)
    resolve_read(plain, GENOME, index, cfg, ResolverStats())
    assert not plain.has_tag(VERDICT_TAG)


# ---------------------------------------------------------------------------
# 2F side — the toy locus from test_2f_min_informative_clip
# ---------------------------------------------------------------------------

EXON1_TAIL = 'ACGTTGCATGCAGTCCATG'
_TOY_SEQ = (
    ('T' * (40 - len(EXON1_TAIL) - 1)) + 'A' + EXON1_TAIL   # exon1  [0, 40)
    + 'GT' + 'N' * 96 + 'AG'                                # intron [40, 140)
    + 'C' * 100                                             # exon2  [140, 240)
)
TOY = {'chrT': _TOY_SEQ}
JUNCTION = {('chrT', 40, 140)}


def _toy_read(clip_len, tag=None):
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'},
                                           'SQ': [{'SN': 'chrT', 'LN': 3_000_000}]})
    r = pysam.AlignedSegment(hdr)
    r.query_name = f'clip{clip_len}'
    r.reference_name = 'chrT'
    r.reference_start = 140
    r.cigartuples = [(4, clip_len), (0, 60)]
    r.is_reverse = False
    r.mapping_quality = 60
    r.query_sequence = _TOY_SEQ[40 - clip_len:40] + 'C' * 60
    if tag is not None:
        r.set_tag(VERDICT_TAG, tag)
    return r


def _rescue(read):
    return rescue_3ss_truncation(read, TOY, JUNCTION, '+')


def test_2f_baseline_rescues_the_informative_clip():
    res = _rescue(_toy_read(35))
    assert res['rescued'] and res['edit_distance'] == 0
    assert res['resolver_verdict'] == '' and res['resolver_relation'] == ''


def test_2f_skips_the_sequence_search_on_low_info():
    res = _rescue(_toy_read(35, tag='L:low_info:3.0:0:35'))
    assert not res['rescued']
    assert res['resolver_verdict'] == 'low_info' and res['resolver_relation'] == 'skip'
    assert COUNTERS['resolver_low_info_skip'] == 1
    assert COUNTERS['assessed'] == 0, "the verdict is reused, never re-derived"


def test_2f_skip_also_applies_below_the_artifact_scope():
    """A 20 nt clip never reaches the >=30 bp short-circuit; the verdict refuses
    the sequence search inside the body instead (the ISSUE-006 floor's shape)."""
    assert _rescue(_toy_read(20))['rescued']
    res = _rescue(_toy_read(20, tag='L:low_info:3.0:0:20'))
    assert not res['rescued'] and res['resolver_relation'] == 'skip'


def test_2f_ignores_a_verdict_for_a_different_clip_length():
    res = _rescue(_toy_read(35, tag='L:low_info:3.0:0:34'))
    assert res['rescued'] and res['resolver_verdict'] == '' and res['resolver_relation'] == ''


def test_2f_ignores_the_other_side():
    res = _rescue(_toy_read(35, tag='R:low_info:3.0:0:35'))
    assert res['rescued'] and res['resolver_relation'] == ''


@pytest.mark.parametrize('token', ['rejected_edit', 'ambiguous', 'no_candidates', 'blowup', 'repeat'])
def test_2f_still_tries_and_counts_the_disagreement(token):
    res = _rescue(_toy_read(35, tag=f'L:{token}:40.0:5000:35'))
    assert res['rescued'] and res['edit_distance'] == 0
    assert res['resolver_verdict'] == token
    assert res['resolver_relation'] == 'rescued_over_refusal'
    assert COUNTERS['resolver_refused_2f_rescued'] == 1


def test_2f_resolved_verdict_is_not_a_disagreement():
    res = _rescue(_toy_read(35, tag='L:resolved:40.0:5000:35'))
    assert res['rescued'] and res['resolver_relation'] == ''


# ---------------------------------------------------------------------------
# Stats plumbing
# ---------------------------------------------------------------------------

def test_correct_read_3prime_carries_the_relation(monkeypatch):
    def _fake(*a, **k):
        return {'rescued': False, 'resolver_verdict': 'low_info', 'resolver_relation': 'skip'}
    monkeypatch.setattr(bp, '_rescue_3ss', _fake)
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrI', 'LN': 2_000_000}]})
    r = pysam.AlignedSegment(hdr)
    r.query_name = 'x'; r.reference_name = 'chrI'; r.reference_start = 1000
    r.cigartuples = [(4, 10), (0, 40)]; r.mapping_quality = 60
    r.query_sequence = 'A' * 10 + 'C' * 40
    row = bp.correct_read_3prime(r, {'chrI': 'C' * 2_000_000}, apply_3ss_rescue=True,
                                 annotated_junctions={('chrI', 980, 1000)})[0]
    assert row['resolver_verdict_5p'] == 'low_info' and row['resolver_2f_relation'] == 'skip'
    s = ProcessingStats()
    s.update_from_result(row)
    assert s.ends_2f_skipped_by_resolver == 1 and s.ends_2f_rescued_over_resolver_refusal == 0


def test_processing_stats_counts_merges_and_writes_rows(tmp_path):
    a = ProcessingStats()
    a.update_from_result({'resolver_2f_relation': 'skip'})
    a.update_from_result({'resolver_2f_relation': 'rescued_over_refusal', 'five_prime_rescued': True})
    a.update_from_result({'resolver_2f_relation': 'skip', 'is_primary_result': False})   # split row: not counted
    a.update_from_result({})
    assert (a.ends_2f_skipped_by_resolver, a.ends_2f_rescued_over_resolver_refusal) == (1, 1)
    b = ProcessingStats(); b.ends_2f_skipped_by_resolver = 4
    a.merge(b)
    assert a.ends_2f_skipped_by_resolver == 5
    assert a.to_dict()['ends_2f_rescued_over_resolver_refusal'] == 1
    a.total_reads_in_bam = 3
    out = tmp_path / 's.tsv'
    write_stats_tsv(a, str(out))
    text = out.read_text()
    assert 'ends_2f_skipped_by_resolver\t5\t' in text
    assert 'ends_2f_rescued_over_resolver_refusal\t1\t' in text
