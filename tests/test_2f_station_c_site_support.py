"""ISSUE-039 — station C for the 5' resolver: the population supplies the attachment tier (2026-09-07).

Kevin's rule, from review card R007 (04b17fc6): *when a candidate junction is already carried by an
independent, well-covered set of reads, a read's own clip faces an ATTACHMENT bar, not a CREATION bar.*
The two-tier floor already existed (E_BITS 18 to create a site, E_BITS_ANNOTATED 12 to attach to one);
what station C changes is who supplies the tier — before this, only the annotation did, so a heavily
used NOVEL junction was treated as if each read were inventing it.

Counts are not the evidence and this test says so twice over: a site is established only by reads that
each cross it with a clean 20-base anchor on BOTH flanks (`_junction_anchor_ok`, the arbiter's Ruling-5
truth criterion (i)), the support is the MAX over aligner arms rather than the sum (five arms are five
alignments of the same reads), and cross-library recurrence never enters — the prescan pool is
per-library by construction.

DEFAULT IS REPORT MODE: the columns are emitted, the tier is unchanged, and nothing drawn moves.
`RECTIFY_2F_STATION_C=attach` is the ON arm, so the two arms differ only in the two new columns until
Kevin rules on R007.
"""
import pysam
import pytest

import rectify.core.bam.bam_processor as bp
from rectify.core.bam.output import CORRECTION_TSV_HEADER, correction_result_to_tsv_row
from rectify.core.splice import splice_aware_5prime as S
from rectify.core.splice.junction_scoring import (
    SITE_ESTABLISHED_MIN_READS,
    SITE_SUPPORT_ANCHOR,
    _collect_junction_counts_core,
    build_junction_pool,
)

# chrT: exon 1 [0, 40) ending in an aperiodic tail, intron [40, 140), exon 2 [140, 240).
EXON1_TAIL = 'ACGTTGCATGCAGTCCATG'
GENOME_SEQ = (
    ('T' * (40 - len(EXON1_TAIL) - 1)) + 'A' + EXON1_TAIL
    + 'GT' + 'N' * 96 + 'AG'
    + 'C' * 100
)
GENOME = {'chrT': GENOME_SEQ}
JUNCTION = ('chrT', 40, 140)


@pytest.fixture(autouse=True)
def _fresh(monkeypatch):
    """Every test starts with station C knowing nothing and in report mode."""
    monkeypatch.delenv('RECTIFY_2F_STATION_C', raising=False)
    monkeypatch.delenv('RECTIFY_2F_NOVEL_GATE', raising=False)
    S.set_site_support(None)
    yield
    S.set_site_support(None)


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


def _clip_read(clip_len, name=None):
    """5' soft clip = the last `clip_len` bases of exon 1 — a perfect placement at JUNCTION."""
    clip = GENOME_SEQ[40 - clip_len:40]
    return _make_read([(4, clip_len), (0, 60)], clip + 'C' * 60,
                      name=name or f'clip{clip_len}')


# The pool scan judges anchor complexity on the READ's own bases, never the genome's, so the
# spanning fixture carries an aperiodic sequence of its own: chrT's exon 2 is a C-homopolymer,
# which `_is_low_complexity_anchor` rightly refuses as an anchor.
APERIODIC = 'ACGTTGCATGCAGTCCATGACGGATCTAGCATCG'


def _spanning_bam(tmp_path, name, n_reads, anchor, junction=JUNCTION):
    """A BAM of `n_reads` reads crossing `junction` with `anchor` clean bases each side."""
    path = tmp_path / f'{name}.minimap2.bam'
    hdr = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrT', 'LN': 3_000_000}]}
    _, i_s, i_e = junction
    left = APERIODIC[:anchor]
    right = APERIODIC[-anchor:]
    with pysam.AlignmentFile(str(path), 'wb', header=hdr) as out:
        for k in range(n_reads):
            a = pysam.AlignedSegment(out.header)
            a.query_name = f'{name}_{k}'
            a.reference_name = 'chrT'
            a.reference_start = i_s - anchor
            a.cigartuples = [(0, anchor), (3, i_e - i_s), (0, anchor)]
            a.mapping_quality = 60
            a.query_sequence = left + right
            out.write(a)
    pysam.index(str(path))
    return path


# -------------------------------------------------------------- the signal, in the pool scan
def test_strict_anchor_counts_only_reads_that_clear_the_wider_anchor(tmp_path):
    """20-base anchors count; 12-base ones pass the pool's own floor of 10 and NOT site support."""
    wide = _spanning_bam(tmp_path, 'wide', 4, SITE_SUPPORT_ANCHOR)
    narrow = _spanning_bam(tmp_path, 'narrow', 4, 12)
    from collections import Counter
    for path, expected in ((wide, 4), (narrow, 0)):
        anchor, raw, _uns, strict = _collect_junction_counts_core(
            str(path), unspliced_out=Counter(), strict_anchor=SITE_SUPPORT_ANCHOR)
        assert anchor[JUNCTION] == 4 and raw[JUNCTION] == 4      # both clear the pool's floor of 10
        assert strict[JUNCTION] == expected


def test_site_support_is_the_max_over_arms_never_the_sum(tmp_path):
    """Five aligner arms are five alignments of the SAME reads (the ALIGNER_FAMILY trap)."""
    arms = [str(_spanning_bam(tmp_path, f'arm{i}', 3, SITE_SUPPORT_ANCHOR)) for i in range(4)]
    _all, _annot, signal = build_junction_pool(arms, set(), return_signal=True)
    assert signal['site_support'][JUNCTION] == 3       # not 12


def test_a_pool_without_the_signal_leaves_station_c_inert():
    S.set_site_support(None)
    assert S.site_support_n(JUNCTION) == 0
    assert not S.site_established(JUNCTION)
    # ... and an unestablished novel landing keeps the CREATION tier in both modes.
    assert not S._attachment_tier(JUNCTION, False)


# -------------------------------------------------------------- the tier decision
def test_established_grants_the_attachment_tier_only_in_attach_mode(monkeypatch):
    S.set_site_support({JUNCTION: SITE_ESTABLISHED_MIN_READS})
    assert S.site_established(JUNCTION)
    assert S.station_c_mode() == 'report'
    assert not S._attachment_tier(JUNCTION, False)         # report mode changes no tier
    assert S._attachment_tier(JUNCTION, True)              # annotated is an attachment regardless
    monkeypatch.setenv('RECTIFY_2F_STATION_C', 'attach')
    assert S.station_c_mode() == 'attach'
    assert S._attachment_tier(JUNCTION, False)


def test_support_one_below_the_floor_is_not_established(monkeypatch):
    monkeypatch.setenv('RECTIFY_2F_STATION_C', 'attach')
    S.set_site_support({JUNCTION: SITE_ESTABLISHED_MIN_READS - 1})
    assert not S.site_established(JUNCTION)
    assert not S._attachment_tier(JUNCTION, False)


def test_support_for_a_different_junction_does_not_transfer(monkeypatch):
    monkeypatch.setenv('RECTIFY_2F_STATION_C', 'attach')
    S.set_site_support({('chrT', 41, 140): 50})
    assert not S._attachment_tier(JUNCTION, False)


# -------------------------------------------------------------- end to end
def _row(read, **kw):
    kw.setdefault('annotated_junctions', set())
    kw.setdefault('pool_chrom_index', bp._build_pool_chrom_index({JUNCTION}))
    return bp.correct_read_3prime(read, GENOME, **kw)[0]


def _floors(monkeypatch, creation, attachment):
    """Move the two floors instead of engineering a sequence into the gap between them. The anchored
    aligner is free to slide a block and score a mutated clip back up (it did, twice, while this test
    was being written), so the fixture that survives is a PERFECT clip with the floors placed around
    its score. What is under test is the tier decision; the bit arithmetic is pinned in
    ``test_2f_evidence_shape.py``."""
    monkeypatch.setenv('RECTIFY_2F_EVIDENCE_BITS', str(creation))
    monkeypatch.setenv('RECTIFY_2F_EVIDENCE_BITS_ANNOTATED', str(attachment))


def test_a_novel_clip_between_the_two_floors_is_drawn_only_when_the_site_is_established(monkeypatch):
    """The whole point, end to end: a perfect 12-base clip scores 24 bits. With the creation floor at
    30 and the attachment floor at 12 the read is refused as a site CREATION and drawn as an
    ATTACHMENT — so its fate turns on whether the population has established the junction."""
    _floors(monkeypatch, 30, 12)
    row_off = _row(_clip_read(12, name='between_off'))
    assert not row_off['five_prime_rescued']
    assert row_off['five_prime_exon_bits'] == pytest.approx(24.0)
    assert row_off['five_prime_site_support'] in ('', None)     # nothing drawn -> blank columns

    monkeypatch.setenv('RECTIFY_2F_STATION_C', 'attach')
    S.set_site_support({JUNCTION: SITE_ESTABLISHED_MIN_READS})
    row_on = _row(_clip_read(12, name='between_on'))
    assert row_on['five_prime_rescued']
    assert row_on['five_prime_landing_annotated'] == 0          # still a NOVEL site
    assert row_on['five_prime_landing_established'] == 1
    assert row_on['five_prime_site_support'] == SITE_ESTABLISHED_MIN_READS


def test_attach_mode_does_not_rescue_a_clip_that_fails_the_attachment_floor_too(monkeypatch):
    """Station C lowers the bar to the attachment tier; it does not remove it. The same 24-bit block
    with the attachment floor at 26 stays refused however well established the site is."""
    _floors(monkeypatch, 30, 26)
    monkeypatch.setenv('RECTIFY_2F_STATION_C', 'attach')
    S.set_site_support({JUNCTION: 50})
    row = _row(_clip_read(12, name='under_both'))
    assert not row['five_prime_rescued']


def test_report_mode_still_records_what_the_population_knew(monkeypatch):
    """A rescue that draws anyway carries the support, so the ON arm is predictable from the OFF arm."""
    S.set_site_support({JUNCTION: 9})
    row = _row(_clip_read(20, name='clean20'))                 # 40 bits: over both floors
    assert row['five_prime_rescued']
    assert row['five_prime_site_support'] == 9
    assert row['five_prime_landing_established'] == 1
    assert S.station_c_mode() == 'report'


def test_the_two_columns_are_last_and_blank_without_a_rescue():
    assert CORRECTION_TSV_HEADER[-2:] == ['five_prime_site_support', 'five_prime_landing_established']
    row = _row(_make_read([(0, 60)], 'C' * 60, name='no_clip'))
    cells = correction_result_to_tsv_row(row)
    assert len(cells) == len(CORRECTION_TSV_HEADER)
    assert cells[-2] == '' and cells[-1] == ''
