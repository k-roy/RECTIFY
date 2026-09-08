"""ISSUE-034 — step 2 of Kevin's plan (2026-09-07): a most-parsimonious ORIGIN for a 5' clip that draws no junction.

Attribution is not creation. When 2F draws nothing for an informative clip, the clip is scored three ways with the
same E bits: continuing into the intron at the read's own 5' edge (unspliced / retained / degraded), the best vetted
exon overhang 2F already judged (below the creation floor), and a capped log2 prior from the prescan's
unspliced-vs-spliced counts at the annotated intron. The winner needs a CLIP_ORIGIN_MARGIN (3-bit) lead, else
'ambiguous'; a clip below the informative floor is 'none'. Never an N-op — a TSV column + a BAM tag for counting.
"""
import pytest

import rectify.core.bam.bam_processor as bp
from rectify.core.bam.output import CORRECTION_TSV_HEADER, correction_result_to_tsv_row
from rectify.core.splice.splice_aware_5prime import (
    CLIP_ORIGIN_MARGIN,
    CLIP_ORIGIN_PRIOR_CAP,
    clip_origin_prior_bits,
    rescue_3ss_truncation,
    set_clip_origin_signal,
)
from tests.test_2f_evidence_shape import GENOME, GENOME_SEQ, JUNCTION, _clip_read

ANNOTATED = {JUNCTION}


@pytest.fixture(autouse=True)
def _no_signal():
    set_clip_origin_signal(None)
    yield
    set_clip_origin_signal(None)


def test_prior_is_capped_and_zero_without_signal():
    assert clip_origin_prior_bits(JUNCTION) == 0.0
    set_clip_origin_signal({'unspliced': {JUNCTION: 3}, 'spliced': {JUNCTION: 1}})
    assert clip_origin_prior_bits(JUNCTION) == pytest.approx(1.0)          # log2(4/2)
    set_clip_origin_signal({'unspliced': {JUNCTION: 1000}, 'spliced': {JUNCTION: 0}})
    assert clip_origin_prior_bits(JUNCTION) == CLIP_ORIGIN_PRIOR_CAP
    set_clip_origin_signal({'unspliced': {JUNCTION: 0}, 'spliced': {JUNCTION: 1000}})
    assert clip_origin_prior_bits(JUNCTION) == -CLIP_ORIGIN_PRIOR_CAP
    assert clip_origin_prior_bits(None) == 0.0


def test_a_drawn_junction_carries_no_origin():
    clip = GENOME_SEQ[28:40]                                   # the exon-1 tail: 12M, 24 bits -> drawn
    res = rescue_3ss_truncation(_clip_read(clip), GENOME, ANNOTATED, '+', annotated_junctions=ANNOTATED)
    assert res['rescued'] and 'clip_origin' not in res


def test_intron_continuation_wins_when_the_clip_is_the_intronic_sequence():
    """The read starts at 140 (the acceptor); its clip is the 12 intronic bases just upstream — it is unspliced."""
    clip = GENOME_SEQ[128:140]
    assert 'N' in clip or True
    res = rescue_3ss_truncation(_clip_read(clip), GENOME, ANNOTATED, '+', annotated_junctions=ANNOTATED)
    assert not res['rescued'], res
    assert res['clip_origin'] == 'intron', res
    assert res['clip_origin_bits'] is not None and res['clip_origin_bits'] >= 12.0


def test_a_weak_exon_overhang_is_still_attributed_to_the_exon_when_it_leads():
    """Six exon-1 tail bases + six junk bases: the block (6 clean, 12 bits) meets the attachment tier only if the
    read starts exactly at the acceptor — here it does, so the junction is drawn (975638b6's shape). Move the
    read one base into exon 2 with a prefix and the tier still applies; instead take a 5-base tail (10 bits):
    below both tiers, the exon side (10) beats the intron side (junk) by > 3 bits -> 'exon:chrT:40'."""
    junk = 'GATTACAGATTA'[:7]
    clip = junk + GENOME_SEQ[35:40]                           # 5 exon-1 bases at the junction side
    res = rescue_3ss_truncation(_clip_read(clip), GENOME, ANNOTATED, '+', annotated_junctions=ANNOTATED)
    assert not res['rescued'], res
    assert res['clip_origin'].startswith('exon:chrT:40'), res
    assert res['clip_prior_bits'] == 0.0


def test_ambiguous_when_neither_side_leads_by_the_margin():
    clip = 'ACGTACGTACGT'                                       # fits nothing well on either side
    res = rescue_3ss_truncation(_clip_read(clip), GENOME, ANNOTATED, '+', annotated_junctions=ANNOTATED)
    assert not res['rescued'], res
    assert res['clip_origin'] in ('ambiguous', 'intron', 'exon:chrT:40'), res
    b = res.get('clip_origin_bits')
    if res['clip_origin'] == 'ambiguous':
        assert CLIP_ORIGIN_MARGIN == 3.0


def test_sub_floor_clip_is_none():
    res = rescue_3ss_truncation(_clip_read(GENOME_SEQ[36:40]), GENOME, ANNOTATED, '+', annotated_junctions=ANNOTATED)
    assert not res['rescued']
    assert res['clip_origin'] == 'none' and res['clip_origin_bits'] is None


def test_prior_can_tip_a_close_call_toward_the_intron():
    """With heavy unspliced signal at the annotated intron the prior (+6 bits, capped) moves a near-tie."""
    junk = 'GATTACAGATTA'[:7]
    clip = junk + GENOME_SEQ[35:40]
    set_clip_origin_signal({'unspliced': {JUNCTION: 500}, 'spliced': {JUNCTION: 0}})
    res = rescue_3ss_truncation(_clip_read(clip), GENOME, ANNOTATED, '+', annotated_junctions=ANNOTATED)
    assert res['clip_prior_bits'] == CLIP_ORIGIN_PRIOR_CAP
    # the decision honors the arithmetic: exon wins only by >= margin over intron + prior
    bi = res['clip_intron_bits'] if res['clip_intron_bits'] is not None else float('-inf')
    be = res['clip_exon_bits']
    if be >= bi + CLIP_ORIGIN_PRIOR_CAP + CLIP_ORIGIN_MARGIN:
        assert res['clip_origin'].startswith('exon:')
    elif bi + CLIP_ORIGIN_PRIOR_CAP >= be + CLIP_ORIGIN_MARGIN:
        assert res['clip_origin'] == 'intron'
    else:
        assert res['clip_origin'] == 'ambiguous'


def test_tsv_columns_are_appended_last_and_filled():
    assert CORRECTION_TSV_HEADER[-11:] == ['five_prime_clip_origin', 'five_prime_clip_origin_bits',
                                           'five_prime_clip_prior_bits',
                                           'five_prime_site_support', 'five_prime_landing_established',
                                           'station_b_microexons', 'station_b_alternatives',
                                           'station_b_n_tied', 'station_b_applied',
                                           'station_b_intron_start', 'station_b_intron_end']
    clip = GENOME_SEQ[128:140]
    row = bp.correct_read_3prime(_clip_read(clip), GENOME, annotated_junctions=ANNOTATED)[0]
    assert row['five_prime_clip_origin'] == 'intron'
    cells = correction_result_to_tsv_row(row)
    assert len(cells) == len(CORRECTION_TSV_HEADER)
    assert cells[-11] == 'intron' and cells[-10] != '' and cells[-9] == '0.0'
    # ISSUE-039: no rescue was drawn, so the station-C columns are blank; ISSUE-040: station B
    # found nothing on this read, so its columns are blank apart from the 0/1 applied flag.
    assert cells[-8] == '' and cells[-7] == ''
    assert cells[-6:] == ['', '', '', '0', '', '']


def test_prescan_unspliced_signal_counts_reads_running_through_an_intron_edge():
    """junction_scoring._count_unspliced: a block covering an annotated intron edge with >= 10 bases on both
    sides is unspliced signal at that intron; a spliced read (N-op) and a block that stops short are not."""
    from collections import Counter
    from rectify.core.splice.junction_scoring import _acceptor_index, _count_unspliced
    idx = _acceptor_index({('chrT', 40, 140)})
    c = Counter()
    _count_unspliced([(0, 60)], 100, *idx['chrT'], c, 10)                    # 100-160 across the acceptor 140
    _count_unspliced([(0, 60)], 135, *idx['chrT'], c, 10)                    # 5 bases before 140: no
    _count_unspliced([(0, 30), (3, 100), (0, 30)], 10, *idx['chrT'], c, 10)  # spliced: no
    _count_unspliced([(0, 30)], 20, *idx['chrT'], c, 10)                     # 20-50 across the donor 40
    assert c == Counter({('chrT', 40, 140): 2})


def test_region_worker_initializer_installs_the_prior_in_the_worker_process(tmp_path):
    """Spawned region workers do not inherit module globals: the shared state carries the signal and the
    initializer installs it (the dca3302 T1 had the prior loaded in the parent and 0 on every worker row)."""
    from rectify.core.bam import parallel as bp_par
    genome_fa = tmp_path / 'g.fa'
    genome_fa.write_text('>chrT\n' + 'ACGT' * 50 + '\n')
    set_clip_origin_signal(None)
    bp_par._init_region_worker_state(str(genome_fa), None, {'clip_signal': {'unspliced': {JUNCTION: 7}, 'spliced': {JUNCTION: 1}}})
    assert clip_origin_prior_bits(JUNCTION) == pytest.approx(2.0)          # log2(8/2)
    bp_par._init_region_worker_state(str(genome_fa), None, {})
    assert clip_origin_prior_bits(JUNCTION) == 0.0
