"""Station B: the correction-regression bypass guard (ISSUE-012).

Every read-evidence signal the classifier uses is LOCAL — junction-proximal
errors in a +/-5 bp window, terminal clip lengths, per-junction annotation
status. A correction can therefore wreck a read globally and still look clean.
Both panel reads that motivated this land every junction ON annotation with
zero clips and junction-proximal error exactly at (not above) the bypass
threshold, so they were labelled ``high_confidence`` and no leg ever offered
them a candidate — while the re-entry arbiter accepts their pre-correction
record on sight (measured, dev/sumner_misplaced_panel_20260904):

    read        cat  hp_ed orig -> corr    gap      gap / ed_orig
    3f941baa    L1       127.25 -> 271.25  144.00       1.13
    45ecf8ed    R5        81.50 -> 292.25  210.75       2.59

The guard's scale is the read's OWN pre-correction hp_ed, so there is no
per-library or per-length constant. ``ratio=1.0`` means "correction more than
doubled the error this read started with". Measured on 1,000 unselected chr5
reads: fires on 1.2% of otherwise-bypassed reads (vs 84.7% for the naive "any
hp_ed regression" rule, which would empty the bypass).
"""

import random

import pysam
import pytest

from rectify.core.consensus.triage import (
    LABEL_HIGH_CONFIDENCE,
    LABEL_TRIAGE,
    REASON_CORRECTION_REGRESSION,
    TriagePolicy,
    classify_bam,
    classify_read,
    triage_realign_bam,
)

GLEN = 3000
DON, ACC = 1000, 1300


def _make_genome():
    rng = random.Random(20260905)
    seq = [rng.choice('ACGT') for _ in range(GLEN)]
    seq[DON:DON + 2] = ['G', 'T']
    seq[ACC - 2:ACC] = ['A', 'G']
    return ''.join(seq)


GSEQ = _make_genome()
GENOME = {'chrA': GSEQ}
HEADER = pysam.AlignmentHeader.from_dict(
    {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrA', 'LN': GLEN}]})

START = 900
EX1, EX2 = 100, 100                       # 100 bp each side of the intron
CLEAN = GSEQ[START:DON] + GSEQ[ACC:ACC + EX2]
J = ('chrA', DON, ACC, '+')
ANNOTATED = {J}
CIGAR = f'{EX1}M{ACC - DON}N{EX2}M'
# The CORRECTED shape carries a 3 bp insertion 50 bp past the acceptor: a real
# `correct` rewrite always edits the CIGAR (the base rewrite comes from indel
# correction, which is a CIGAR edit by construction), and it is far enough from
# the junction to leave every local classifier signal clean. Measured on the
# panel and the 1,000-read hold-out: 0 reads differ from their pre-correction
# record in SEQ alone, so a same-CIGAR fixture would not represent anything.
CIGAR_CORR = f'{EX1}M{ACC - DON}N50M3I47M'


def _read(query, name='r0', cigar=CIGAR, start=START):
    r = pysam.AlignedSegment(header=HEADER)
    r.query_name = name
    r.query_sequence = query
    r.flag = 0
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigarstring = cigar
    return r


def _mutate(query, n):
    """n mismatches placed AWAY from the junction, so every local signal the
    classifier reads stays clean — the shape both panel reads have."""
    q = list(query)
    # exon-interior positions only: >5 bp from either junction edge
    spots = list(range(5, EX1 - 6)) + list(range(EX1 + 6, EX1 + EX2 - 5))
    assert n <= len(spots)
    for i in spots[:n]:
        q[i] = 'A' if q[i] != 'A' else 'C'
    return ''.join(q)


def _hp(read):
    from rectify.core.consensus.corrected_consensus import _cigar_hp_edit_distance
    return _cigar_hp_edit_distance(read, GENOME, None)


# ---------------------------------------------------------------------------
# Fixture sanity: the corrected read must be locally clean, or the test proves
# nothing about the guard (the other signals would have triaged it anyway).
# ---------------------------------------------------------------------------

def test_the_wrecked_read_is_locally_clean_and_would_be_bypassed():
    """Anti-vacuity, and the shape of the real defect: the damage sits 50 bp
    past the acceptor, so zero clips, the junction annotated and exactly
    placed, and no junction-proximal error — every signal the classifier has
    says high-confidence, while hp_ed is 4x the original's."""
    corrected = _read(_mutate(CLEAN, 10), cigar=CIGAR_CORR)
    res = classify_read(corrected, GENOME, {J[:3]}, TriagePolicy())
    assert res.label == LABEL_HIGH_CONFIDENCE
    assert res.reasons == []
    assert res.clip_5p == 0 and res.clip_3p == 0
    assert res.n_unannotated == 0 and res.n_junctions == 1
    assert res.junction_proximal_errors <= 1.0
    assert _hp(_read(_mutate(CLEAN, 10))) == pytest.approx(10.0)
    assert _hp(corrected) == pytest.approx(43.75)


# ---------------------------------------------------------------------------
# The guard
# ---------------------------------------------------------------------------

def test_a_read_correction_more_than_doubled_cannot_be_bypassed():
    """The two panel reads in miniature: hp_ed 10.00 -> 43.75 is a ratio of
    3.375, and every local signal still reads clean."""
    original = _read(_mutate(CLEAN, 10))
    corrected = _read(_mutate(CLEAN, 10), cigar=CIGAR_CORR)
    assert _hp(original) == pytest.approx(10.0)
    assert _hp(corrected) == pytest.approx(43.75)

    res = classify_read(corrected, GENOME, {J[:3]}, TriagePolicy(),
                        original=original)
    assert res.label == LABEL_TRIAGE
    assert REASON_CORRECTION_REGRESSION in res.reasons
    assert res.correction_regression == pytest.approx(3.375)
    assert res.junction_leg, 'must reach the leg where the ISSUE-009 candidate wins'


def test_a_genuinely_high_confidence_read_still_bypasses():
    """The other half of the contract. Correction moved hp_ed 40.00 -> 63.75:
    worse, and the arbiter WOULD take the original back — but by only 0.59 of
    the read's own error burden, which is not evidence of a wrecked alignment.
    This is the case that separates this guard from the naive 'any hp_ed
    regression' rule (84.7% of unselected reads vs 1.2%)."""
    original = _read(_mutate(CLEAN, 40))                     # hp_ed 40.00
    corrected = _read(_mutate(CLEAN, 30), cigar=CIGAR_CORR)  # hp_ed 63.75
    res = classify_read(corrected, GENOME, {J[:3]}, TriagePolicy(),
                        original=original)
    assert res.label == LABEL_HIGH_CONFIDENCE
    assert res.reasons == []
    assert res.correction_regression == pytest.approx(23.75 / 40)

    # ... and this is NOT the trivial case where the original is better anyway:
    # the arbiter WOULD take it back, which is exactly what makes the naive
    # "any hp_ed regression" rule fire on 84.7% of unselected reads.
    from rectify.core.consensus.triage import reentry_accept
    assert reentry_accept(corrected, original, GENOME) is True


def test_a_correction_that_improved_the_read_bypasses():
    original = _read(_mutate(CLEAN, 30), cigar=CIGAR_CORR)   # hp_ed 63.75
    corrected = _read(_mutate(CLEAN, 10))                    # hp_ed 10.00
    res = classify_read(corrected, GENOME, {J[:3]}, TriagePolicy(),
                        original=original)
    assert res.label == LABEL_HIGH_CONFIDENCE
    assert res.correction_regression == pytest.approx(-53.75 / 63.75)


def test_a_perfect_original_made_worse_always_triages():
    """The degenerate input, decided explicitly rather than by accident: with
    hp_ed(original) == 0 the ratio has no scale, so ANY regression triages and
    the recorded ratio is inf."""
    res = classify_read(_read(CLEAN, cigar=CIGAR_CORR), GENOME, {J[:3]},
                        TriagePolicy(), original=_read(CLEAN))
    assert res.label == LABEL_TRIAGE
    assert REASON_CORRECTION_REGRESSION in res.reasons
    assert res.correction_regression == float('inf')

    # ... and a perfect original left perfect is NOT triaged
    same = classify_read(_read(CLEAN), GENOME, {J[:3]}, TriagePolicy(),
                         original=_read(CLEAN))
    assert same.label == LABEL_HIGH_CONFIDENCE
    assert same.correction_regression == 0.0


def test_the_ratio_is_tunable():
    original = _read(_mutate(CLEAN, 40))                     # hp_ed 40.00
    corrected = _read(_mutate(CLEAN, 30), cigar=CIGAR_CORR)  # ratio 0.594
    assert classify_read(corrected, GENOME, {J[:3]},
                         TriagePolicy(max_correction_regression_ratio=0.25),
                         original=original).label == LABEL_TRIAGE
    assert classify_read(corrected, GENOME, {J[:3]},
                         TriagePolicy(max_correction_regression_ratio=1.0),
                         original=original).label == LABEL_HIGH_CONFIDENCE


# ---------------------------------------------------------------------------
# Unreachability without a pre-correction record
# ---------------------------------------------------------------------------

def test_the_guard_cannot_fire_without_an_original(tmp_path):
    """Default-off identity. The signal needs the pre-correction record, so
    with no `original` the classifier is byte-identical to its old behavior —
    including through `classify_bam`, which has no `original` parameter."""
    corrected = _read(_mutate(CLEAN, 10), cigar=CIGAR_CORR)
    res = classify_read(corrected, GENOME, {J[:3]}, TriagePolicy())
    assert res.label == LABEL_HIGH_CONFIDENCE
    assert REASON_CORRECTION_REGRESSION not in res.reasons
    assert res.correction_regression is None and res.hp_ed_original is None

    path = tmp_path / 'corr.bam'
    with pysam.AlignmentFile(str(path), 'wb', header=HEADER) as out:
        out.write(corrected)
    for r in classify_bam(str(path), GENOME, {J}):
        assert REASON_CORRECTION_REGRESSION not in r.reasons


# ---------------------------------------------------------------------------
# End to end: the guard feeds the ISSUE-009 candidate
# ---------------------------------------------------------------------------

def _bam(path, reads):
    with pysam.AlignmentFile(str(path), 'wb', header=HEADER) as out:
        for r in reads:
            out.write(r)
    return str(path)


def test_a_rescued_read_is_reverted_end_to_end(tmp_path):
    """The point of the guard: the read leaves the bypass, reaches the junction
    leg, and the pre-correction candidate wins it back. Without the guard the
    same input is bypassed and nothing happens to it."""
    original = _read(_mutate(CLEAN, 10))                     # hp_ed 10.00
    corrected = _read(_mutate(CLEAN, 10), cigar=CIGAR_CORR)  # hp_ed 43.75
    inp = _bam(tmp_path / 'corrected.bam', [corrected])
    orig = _bam(tmp_path / 'original.bam', [original])

    out = tmp_path / 'triaged.bam'
    rows, stats = triage_realign_bam(
        inp, str(out), genome=GENOME, annotated_junctions=ANNOTATED,
        original_bams=[orig], sort_and_index=False)
    assert stats['correction_regression'] == 1
    assert stats['high_confidence'] == 0 and stats['junction_leg'] == 1
    assert stats['orig_accepted'] == 1
    assert rows[0]['reverted_to_original'] is True
    assert rows[0]['correction_regression'] == 3.375
    with pysam.AlignmentFile(str(out), 'rb', check_sq=False) as bam:
        rec = next(iter(bam.fetch(until_eof=True)))
    assert rec.query_sequence == original.query_sequence

    # Same input, guard disabled by an unreachable ratio -> bypassed, untouched.
    out2 = tmp_path / 'triaged_noguard.bam'
    rows2, stats2 = triage_realign_bam(
        inp, str(out2), genome=GENOME, annotated_junctions=ANNOTATED,
        original_bams=[orig], sort_and_index=False,
        policy=TriagePolicy(max_correction_regression_ratio=1e9))
    assert stats2['correction_regression'] == 0
    assert stats2['high_confidence'] == 1 and stats2['junction_leg'] == 0
    assert stats2['orig_accepted'] == 0
    assert rows2[0]['reverted_to_original'] is False


def test_cli_exposes_the_ratio():
    from rectify.cli import create_parser
    args = create_parser().parse_args(
        ['triage', 'in.bam', '-o', 'out', '--Scer',
         '--max-correction-regression-ratio', '0.5'])
    assert args.max_correction_regression_ratio == 0.5
    assert create_parser().parse_args(
        ['triage', 'in.bam', '-o', 'out', '--Scer']
    ).max_correction_regression_ratio == 1.0
