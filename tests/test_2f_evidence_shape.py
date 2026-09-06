"""ISSUE-028 / invariant E — the placed 5' block must be EVIDENCE for every landing (2026-09-06).

Kevin's read-level review of the 4993253 T0 set found two unchanged controls rescued onto ANNOTATED junctions with
placed blocks of 28 % (26f8fb45: 32 bp, 9=/23X) and 42 % (04b17fc6: 26 bp, 11=/15X) identity — an annotated landing
bypassed the novel-site verdict and the only bound that did apply (invariant C, I+D <= 0.5*matched) ignores
mismatches — and two baseline rescues whose exon CIGAR STARTED with an insertion (`8I6M`, `2M4I2M1I2M`): a leading I
is a soft clip in disguise. The rule (Kevin's ruling, 2026-09-06): identity >= 0.8, no leading I/D (emitted as S),
and an evidence SCORE over the placed block in bits (a match = log2 4 = 2 bits; a mismatch / affine gap at the
anchored aligner's own constants over 2; homopolymer-run errors half) at a cutoff DERIVED from the chance-match
model (docs/algorithms/read_level_review.md § 3; the table is in ISSUE-028). A hard clean-run floor was the wrong
primitive: `6=1X7=` is an excellent 14-nt anchor.

Pinned here, hermetically: the shape measurement (`local_aligner.evidence_shape`: identity, bits, the worked
examples), the leading-indel strip (`strip_leading_indel`, and that the writer emits the stripped bases as S), the
verdict (`_evidence_floor_refusal` with the module constants and the env overrides), the two TSV columns, and the
end-to-end refusal of a low-identity block on an ANNOTATED landing in both gate modes. The replay of the five
reviewed reads lives in ``test_2f_replay_review_4993253.py`` (collaborator data, skip-if-absent).
"""
import pysam
import pytest

import rectify.core.bam.bam_processor as bp
from rectify.core.align.local_aligner import (
    EvidenceShape,
    cigar_str_to_ops,
    evidence_shape,
    strip_leading_indel,
)
from rectify.core.bam.output import CORRECTION_TSV_HEADER, correction_result_to_tsv_row
from rectify.core.bam.read_edits import extend_read_5prime_for_junction_rescue
from rectify.core.splice.splice_aware_5prime import (
    E_BITS,
    E_IDENTITY,
    EVIDENCE_REFUSALS,
    EXON_BITS_REFUSAL,
    EXON_IDENTITY_REFUSAL,
    NOVEL_EXON_REFUSALS,
    PLACEMENT_REFUSALS,
    _evidence_floor_refusal,
    evidence_floor,
    rescue_3ss_truncation,
)

# A synthetic chrT: exon 1 [0, 40) — a 20-T homopolymer then an aperiodic tail — intron [40, 140), exon 2 [140, 240).
EXON1_TAIL = 'ACGTTGCATGCAGTCCATG'
GENOME_SEQ = (
    ('T' * (40 - len(EXON1_TAIL) - 1)) + 'A' + EXON1_TAIL
    + 'GT' + 'N' * 96 + 'AG'
    + 'C' * 100
)
GENOME = {'chrT': GENOME_SEQ}
JUNCTION = ('chrT', 40, 140)
# A minus-strand fixture: exon 1 (RNA 5') at [140, 160) is aperiodic so no homopolymer forgiveness applies.
MINUS_SEQ = 'N' * 140 + 'ACGTTGCATGCAGTCCATGA' + 'T' * 40
# A plus-strand fixture whose exon-1 tail is aperiodic for 30 nt (the worked examples).
PLUS_SEQ = 'A' * 10 + 'TGCATGCAGTCCATGACGGATCTAGCATCG' + 'GT' + 'N' * 96 + 'AG' + 'C' * 100
_FLIP = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}


def _ops(s):
    return cigar_str_to_ops(s)


def _mutate(seq, positions):
    out = list(seq)
    for p in positions:
        out[p] = _FLIP[out[p]]
    return ''.join(out)


# ---------------------------------------------------------------------------------------------- the shape
def test_shape_of_a_clean_plus_block():
    clip = PLUS_SEQ[40 - 15:40]
    sh = evidence_shape(_ops('15M'), clip, PLUS_SEQ, 40, 140, '+')
    assert sh == EvidenceShape(15, 0, 1.0, 30.0, 15, 0)


def test_kevins_worked_examples():
    """The bits the ruling names: `6=1X7=` 24 (pass), 22f609c6's `1I15M` 13=/2X 19.5 (pass), `5M4D2M` 10 (fail),
    `9=/23X` -28 and `11=/15X` -8 (fail) — at ANY cutoff in the model's 14-24 range the verdicts are the same."""
    six_x_seven = evidence_shape(_ops('14M'), _mutate(PLUS_SEQ[26:40], [6]), PLUS_SEQ, 40, 140, '+')
    assert (six_x_seven.matched, six_x_seven.mismatches, six_x_seven.bits) == (13, 1, 24.0)
    assert six_x_seven.junction_clean_run == 7                    # the clean-run floor would have refused it
    ref = MINUS_SEQ[140:155]
    tip = evidence_shape(_ops('1I15M'), 'G' + _mutate(ref, [13, 14]), MINUS_SEQ, 40, 140, '-')
    assert (tip.matched, tip.mismatches, tip.bits) == (13, 2, 19.5)
    assert tip.identity == pytest.approx(13 / 15)
    ref11 = PLUS_SEQ[29:40]
    gap = evidence_shape(_ops('5M4D2M'), ref11[:5] + ref11[9:], PLUS_SEQ, 40, 140, '+')
    assert (gap.matched, gap.mismatches, gap.identity, gap.bits) == (7, 0, 1.0, 10.0)
    # mismatches at block indices 2..24 (genome 10..32, the aperiodic tail — PLUS_SEQ[8:10] is inside the A-run)
    junk32 = evidence_shape(_ops('32M'), _mutate(PLUS_SEQ[8:40], list(range(2, 25))), PLUS_SEQ, 40, 140, '+')
    assert (junk32.matched, junk32.mismatches, junk32.bits) == (9, 23, -28.0)
    junk26 = evidence_shape(_ops('26M'), _mutate(PLUS_SEQ[14:40], list(range(15))), PLUS_SEQ, 40, 140, '+')
    assert (junk26.matched, junk26.mismatches, junk26.bits) == (11, 15, -8.0)
    assert junk26.identity == pytest.approx(11 / 26)
    for cutoff in (14, 16, 18, 19.5):
        assert six_x_seven.bits >= cutoff and tip.bits >= cutoff
        assert gap.bits < cutoff and junk32.bits < cutoff and junk26.bits < cutoff
    assert _evidence_floor_refusal(six_x_seven) == '' and _evidence_floor_refusal(tip) == ''
    assert _evidence_floor_refusal(gap) == EXON_BITS_REFUSAL
    assert _evidence_floor_refusal(junk32) == EXON_IDENTITY_REFUSAL
    assert _evidence_floor_refusal(junk26) == EXON_IDENTITY_REFUSAL


def test_junction_run_is_reported_from_the_junction_side_per_strand():
    # plus: junction = RIGHT end. Mismatch 5 bases in from the junction -> run 4 (bases 16..19 clean)
    clip = _mutate(PLUS_SEQ[40 - 20:40], [15])
    sh = evidence_shape(_ops('20M'), clip, PLUS_SEQ, 40, 140, '+')
    assert (sh.matched, sh.mismatches, sh.junction_clean_run, sh.bits) == (19, 1, 4, 36.0)
    # minus: junction = LEFT end (BAM orientation); the block sits at [intron_end, intron_end + n)
    ref = MINUS_SEQ[140:160]
    sh_m = evidence_shape(_ops('20M'), _mutate(ref, [4]), MINUS_SEQ, 40, 140, '-')
    assert (sh_m.matched, sh_m.mismatches, sh_m.junction_clean_run) == (19, 1, 4)
    sh_m2 = evidence_shape(_ops('20M'), _mutate(ref, [15]), MINUS_SEQ, 40, 140, '-')
    assert sh_m2.junction_clean_run == 15


def test_gaps_are_charged_affinely_in_bits():
    ref = PLUS_SEQ[20:40]                                   # 20 ref bases
    # one 2-base deletion in the middle: 18 matched (36) - (2 + 2*0.5) = 33
    d2 = evidence_shape(_ops('9M2D9M'), ref[:9] + ref[11:], PLUS_SEQ, 40, 140, '+')
    assert (d2.matched, d2.bits) == (18, 33.0)
    # one 3-base insertion: 20 matched (40) - (2 + 1.5) = 36.5
    i3 = evidence_shape(_ops('10M3I10M'), ref[:10] + 'GGG' + ref[10:], PLUS_SEQ, 40, 140, '+')
    assert (i3.matched, i3.bits) == (20, 36.5)


def test_errors_inside_a_homopolymer_run_are_charged_half():
    """GENOME_SEQ's exon 1 opens with a 20-T run (genome homopolymer >= 5)."""
    clip = _mutate(GENOME_SEQ[40 - 30:40], [2])            # index 2 -> genome position 12, inside the T run
    sh = evidence_shape(_ops('30M'), clip, GENOME_SEQ, 40, 140, '+')
    assert (sh.matched, sh.mismatches) == (29, 1)
    assert sh.identity == pytest.approx(29 / 29.5)
    assert sh.bits == 58.0 - 1.0
    # the same single mismatch outside a run costs a whole one
    clip2 = _mutate(GENOME_SEQ[40 - 30:40], [25])          # genome position 35, inside the aperiodic tail
    sh2 = evidence_shape(_ops('30M'), clip2, GENOME_SEQ, 40, 140, '+')
    assert sh2.identity == pytest.approx(29 / 30) and sh2.bits == 58.0 - 2.0
    # a deletion inside the run: (2 + 0.5) / 2 = 1.25
    ref = GENOME_SEQ[10:40]
    d1 = evidence_shape(_ops('5M1D24M'), ref[:5] + ref[6:], GENOME_SEQ, 40, 140, '+')
    assert d1.bits == 58.0 - 1.25
    # forgiveness does not repair the reported run: the anchor breaks on a junction-side mismatch either way
    assert evidence_shape(_ops('30M'), _mutate(GENOME_SEQ[40 - 30:40], [29]), GENOME_SEQ, 40, 140, '+').junction_clean_run == 0


# ---------------------------------------------------------------------------------------------- the strip
def test_leading_insertion_becomes_a_soft_clip():
    """6a6e0b3c's baseline `8I6M` = six placed bases, not fourteen."""
    assert strip_leading_indel(_ops('8I6M'), '+') == (_ops('8S6M'), 8)
    # minus strand: the 5' (free) end is the RIGHT end in BAM orientation
    assert strip_leading_indel(_ops('6M8I'), '-') == (_ops('6M8S'), 8)
    # a leading D consumes no query: dropped; mixed leading runs collapse into one S
    assert strip_leading_indel(_ops('3D6M'), '+') == (_ops('6M'), 0)
    assert strip_leading_indel(_ops('2D3I1D6M'), '+') == (_ops('3S6M'), 3)
    # nothing to strip
    assert strip_leading_indel(_ops('6M2I4M'), '+') == (_ops('6M2I4M'), 0)
    assert strip_leading_indel(_ops('2M4I2M1I2M'), '+') == (_ops('2M4I2M1I2M'), 0)   # 38722d08: no LEADING indel
    # an all-indel CIGAR places nothing
    assert strip_leading_indel(_ops('5I'), '+') == ([], 5)
    assert strip_leading_indel([], '+') == ([], 0)


def test_stripped_block_is_measured_without_the_clip_bases():
    """`8I6M` on a 14-nt clip: the shape is that of the six placed bases (12 bits < E_BITS -> refused)."""
    six = PLUS_SEQ[40 - 6:40]
    clip = 'GATTAC' + 'AG' + six
    ops, unplaced = strip_leading_indel(_ops('8I6M'), '+')
    sh = evidence_shape(ops, clip, PLUS_SEQ, 40, 140, '+')
    assert unplaced == 8 and sh == EvidenceShape(6, 0, 1.0, 12.0, 6, 8)
    assert _evidence_floor_refusal(sh) == EXON_BITS_REFUSAL


def _record(cigar, seq, start, strand):
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrT', 'LN': 3_000_000}]})
    r = pysam.AlignedSegment(hdr)
    r.query_name = 'w'
    r.reference_name = 'chrT'
    r.reference_start = start
    r.cigartuples = cigar
    r.is_reverse = (strand == '-')
    r.query_sequence = seq
    return r


def test_writer_emits_the_stripped_bases_as_s_on_both_strands():
    """The exon ops are prepended (plus) / appended (minus) verbatim and S is query-consuming in the writer's
    span check, so `8S6M` reaches the BAM as `8S 6M N ...` and never as an insertion."""
    six = GENOME_SEQ[34:40]
    r = _record([(4, 14), (0, 60)], 'GATTACAG' + six + 'C' * 60, 140, '+')
    assert extend_read_5prime_for_junction_rescue(r, 39, 14, '+', exon_cigar_str='8S6M')
    assert r.cigartuples[:3] == [(4, 8), (0, 6), (3, 100)] and r.reference_start == 34
    # minus: body over exon 2 ends at intron_start 40, the block starts at intron_end 140
    r2 = _record([(0, 30), (4, 14)], 'C' * 30 + GENOME_SEQ[140:146] + 'GATTACAG', 10, '-')
    assert extend_read_5prime_for_junction_rescue(r2, 140, 14, '-', exon_cigar_str='6M8S')
    assert r2.cigartuples == [(0, 30), (3, 100), (0, 6), (4, 8)]


# ---------------------------------------------------------------------------------------------- the verdict
def test_defaults_and_registered_tokens():
    assert E_IDENTITY == 0.8
    assert 14 <= E_BITS <= 19.5          # the model's range; 22f609c6 (19.5) must pass
    assert evidence_floor() == (E_IDENTITY, E_BITS)
    assert set(EVIDENCE_REFUSALS) == {EXON_IDENTITY_REFUSAL, EXON_BITS_REFUSAL}
    assert set(EVIDENCE_REFUSALS) <= set(PLACEMENT_REFUSALS)
    assert not (set(EVIDENCE_REFUSALS) & set(NOVEL_EXON_REFUSALS))
    assert _evidence_floor_refusal(None) == ''


def test_env_overrides_move_the_operating_point(monkeypatch):
    sh = EvidenceShape(12, 2, 12 / 14, 20.0, 9, 0)             # identity 0.857, 20 bits
    assert _evidence_floor_refusal(sh) == ''
    monkeypatch.setenv('RECTIFY_2F_EVIDENCE_BITS', '22')
    assert evidence_floor() == (E_IDENTITY, 22.0)
    assert _evidence_floor_refusal(sh) == EXON_BITS_REFUSAL
    monkeypatch.setenv('RECTIFY_2F_EVIDENCE_IDENTITY', '0.9')
    assert _evidence_floor_refusal(sh) == EXON_IDENTITY_REFUSAL    # identity checked first
    monkeypatch.setenv('RECTIFY_2F_EVIDENCE_IDENTITY', 'garbage')
    monkeypatch.setenv('RECTIFY_2F_EVIDENCE_BITS', 'garbage')
    assert evidence_floor() == (E_IDENTITY, E_BITS)


# ---------------------------------------------------------------------------------------------- end to end
def _clip_read(clip, name='r'):
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrT', 'LN': 3_000_000}]})
    r = pysam.AlignedSegment(hdr)
    r.query_name = name
    r.reference_name = 'chrT'
    r.mapping_quality = 60
    r.reference_start = 140
    r.cigartuples = [(4, len(clip)), (0, 60)]
    r.query_sequence = clip + 'C' * 60
    return r


def _cells(row):
    return dict(zip(CORRECTION_TSV_HEADER, correction_result_to_tsv_row(row)))


@pytest.mark.parametrize('gate', ['report', 'refuse'])
def test_low_identity_block_on_an_annotated_landing_is_refused_in_both_gate_modes(gate, monkeypatch):
    """04b17fc6's class on chrT: a 20-nt clip with five ISOLATED mismatches (an isolated mismatch, -4, beats a
    1-base gap, -5, so the DP keeps them as X): identity 15/20 = 0.75, over the ANNOTATED junction. A placement
    decision: no rescue, the token in five_prime_rescue_refused, the shape in the TSV."""
    monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', gate)
    good = GENOME_SEQ[20:40]
    noisy = _mutate(good, [0, 3, 6, 9, 12])
    res = rescue_3ss_truncation(_clip_read(noisy), GENOME, {JUNCTION}, '+', annotated_junctions={JUNCTION})
    assert not res['rescued'], res
    assert res.get('clip_refused') == EXON_IDENTITY_REFUSAL, res
    assert res['exon_identity'] == pytest.approx(0.75) and res['exon_bits'] == 30.0 - 10.0
    row = bp.correct_read_3prime(_clip_read(noisy), GENOME, annotated_junctions={JUNCTION})[0]
    assert not row['five_prime_rescued'] and row['five_prime_rescue_refused'] == EXON_IDENTITY_REFUSAL
    assert (40, 140) not in [tuple(j) for j in row['junctions']]
    cells = _cells(row)
    assert cells['five_prime_exon_identity'] == '0.75' and cells['five_prime_exon_bits'] == '20.0'
    # the clean clip on the same junction draws, with a full-marks shape in the columns
    ok = bp.correct_read_3prime(_clip_read(good), GENOME, annotated_junctions={JUNCTION})[0]
    assert ok['five_prime_rescued'] and ok['five_prime_rescue_refused'] == ''
    assert (40, 140) in [tuple(j) for j in ok['junctions']]
    assert _cells(ok)['five_prime_exon_identity'] == '1.00' and _cells(ok)['five_prime_exon_bits'] == '40.0'


def test_short_clean_block_on_an_annotated_landing_draws_no_junction():
    """A block that is clean but carries too few bits (6 matched = 12 bits): no junction reaches the BAM."""
    clip = GENOME_SEQ[34:40]
    res = rescue_3ss_truncation(_clip_read(clip), GENOME, {JUNCTION}, '+', annotated_junctions={JUNCTION})
    assert not res['rescued']
    row = bp.correct_read_3prime(_clip_read(clip), GENOME, annotated_junctions={JUNCTION})[0]
    assert not row['five_prime_rescued'] and (40, 140) not in [tuple(j) for j in row['junctions']]


def test_tsv_columns_are_blank_without_a_placed_block():
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrT', 'LN': 3_000_000}]})
    r = pysam.AlignedSegment(hdr)
    r.query_name = 'plain'
    r.reference_name = 'chrT'
    r.reference_start = 200
    r.cigartuples = [(0, 30)]
    r.query_sequence = 'C' * 30
    r.mapping_quality = 60
    row = bp.correct_read_3prime(r, GENOME, annotated_junctions={JUNCTION})[0]
    cells = _cells(row)
    assert cells['five_prime_exon_identity'] == '' and cells['five_prime_exon_bits'] == ''
    assert CORRECTION_TSV_HEADER[-2:] == ['five_prime_exon_identity', 'five_prime_exon_bits']
