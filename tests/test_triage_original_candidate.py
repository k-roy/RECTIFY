"""Station B: the PRE-CORRECTION alignment as a triage candidate (ISSUE-009).

The junction leg's proposer is ``refine_bam_junctions(motif_blind=True)`` —
Module 2H itself. On a BAM 2H has already refined it re-derives its own fixed
point, so the one placement it structurally cannot propose is the one it moved
away from. Measured on the Sumner RNA004 panel
(dev/sumner_misplaced_panel_20260904 CP3): triage flagged the damage correctly
(R1 16/16, R2 16/16 distressed) and ``reentry_accept`` accepted the ORIGINAL
back on 42/42 reads when handed it, yet only 6/93 damaged reads were reverted —
a missing-candidate failure, not an arbiter failure.

``original_bams`` adds that candidate. It adds NOTHING else: the same strict
hp_ed ``reentry_accept`` decides, so a read whose 2H placement genuinely wins is
never reverted.
"""

import random
from pathlib import Path

import pysam
import pytest

from rectify.core.consensus.triage import triage_realign_bam

GLEN = 3000
DON, ACC = 1000, 1300          # planted GT..AG intron [1000, 1300)
START = 960                    # read starts 40 bp before the donor


def _make_genome():
    rng = random.Random(20260905)
    seq = [rng.choice('ACGT') for _ in range(GLEN)]
    seq[DON:DON + 2] = ['G', 'T']
    seq[ACC - 2:ACC] = ['A', 'G']
    return ''.join(seq)


GSEQ = _make_genome()
GENOME = {'chrA': GSEQ, 'chrB': 'ACGT' * 125}

# Input/corrected BAM header order; the pre-correction BAM deliberately uses
# the OPPOSITE @SQ order below, so a record carried across by reference_id
# would land on the wrong chromosome.
HDR_IN = pysam.AlignmentHeader.from_dict(
    {'HD': {'VN': '1.6'},
     'SQ': [{'SN': 'chrA', 'LN': GLEN}, {'SN': 'chrB', 'LN': 500}]})
HDR_ORIG = pysam.AlignmentHeader.from_dict(
    {'HD': {'VN': '1.6'},
     'SQ': [{'SN': 'chrB', 'LN': 500}, {'SN': 'chrA', 'LN': GLEN}]})

# The read's query is the two exons, so the TRUE placement is a perfect match.
QUERY = GSEQ[START:DON] + GSEQ[ACC:ACC + 40]

# hp_ed 0.00 — donor at the planted GT
CIGAR_TRUE = '40M300N40M'
# hp_ed 1.25 — the exact shape 2H writes for a 1-nt donor slide: the exon
# SEQUENCE has not moved, only the N-op coordinate, paid for with a
# compensating 1I (junction_refiner.py L1300-1309).
CIGAR_SLID = '39M1I301N40M'

J_TRUE = ('chrA', DON, ACC, '+')
J_SLID = ('chrA', DON - 1, ACC, '+')


def _read(header, cigar, name='r0', chrom='chrA'):
    r = pysam.AlignedSegment(header=header)
    r.query_name = name
    r.query_sequence = QUERY
    r.flag = 0
    r.reference_id = header.references.index(chrom)
    r.reference_start = START
    r.mapping_quality = 60
    r.cigarstring = cigar
    return r


def _bam(path, header, reads):
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for r in reads:
            out.write(r)
    return str(path)


def _only(path):
    with pysam.AlignmentFile(str(path), 'rb', check_sq=False) as bam:
        recs = list(bam.fetch(until_eof=True))
    assert len(recs) == 1
    return recs[0]


def _run(tmp_path, in_cigar, orig_cigar, annotated, n_reads=1, **kw):
    """One read: `in_cigar` in the corrected BAM, `orig_cigar` pre-correction."""
    names = [f'r{i}' for i in range(n_reads)]
    inp = _bam(tmp_path / 'corrected.bam', HDR_IN,
               [_read(HDR_IN, in_cigar, n) for n in names])
    orig = _bam(tmp_path / 'original.bam', HDR_ORIG,
                [_read(HDR_ORIG, orig_cigar, n) for n in names])
    out = tmp_path / 'triaged.bam'
    rows, stats = triage_realign_bam(
        inp, str(out), genome=GENOME, annotated_junctions=annotated,
        original_bams=[orig], sort_and_index=False, **kw)
    return rows, stats, out


def test_the_fixture_is_what_the_tests_assume():
    """Anti-vacuity: the slid placement really is the WORSE of the two, and it
    is worse only because of the compensating insertion."""
    from rectify.core.consensus.corrected_consensus import _cigar_hp_edit_distance
    true_ed = _cigar_hp_edit_distance(_read(HDR_IN, CIGAR_TRUE), GENOME, None)
    slid_ed = _cigar_hp_edit_distance(_read(HDR_IN, CIGAR_SLID), GENOME, None)
    assert true_ed == 0.0
    assert slid_ed == pytest.approx(1.25)      # the 1I, nothing else
    assert GSEQ[DON:DON + 2] == 'GT' and GSEQ[ACC - 2:ACC] == 'AG'


def test_distressed_read_is_restored_when_the_original_wins(tmp_path):
    """The 2H-damaged incumbent is reverted to the pre-correction alignment."""
    rows, stats, out = _run(tmp_path, CIGAR_SLID, CIGAR_TRUE,
                            annotated={J_TRUE})
    assert stats['junction_leg'] == 1, 'the read must reach the junction leg'
    assert stats['orig_leg'] == 1 and stats['orig_proposed'] == 1
    assert stats['orig_accepted'] == 1
    assert stats['accepted'] == 1

    rec = _only(out)
    assert rec.cigarstring == CIGAR_TRUE
    assert rec.reference_start == START
    assert rows[0]['reverted_to_original'] is True
    assert rows[0]['accepted'] and rows[0]['realigned']


def test_a_read_whose_proposal_wins_is_not_reverted(tmp_path):
    """The mirror image: the incumbent is the BETTER alignment, so offering the
    original must change nothing. Only the arbiter decides — this leg adds a
    candidate, never an acceptance path."""
    rows, stats, out = _run(tmp_path, CIGAR_TRUE, CIGAR_SLID,
                            annotated={J_SLID})
    assert stats['junction_leg'] == 1
    assert stats['orig_leg'] == 1 and stats['orig_proposed'] == 1
    assert stats['orig_accepted'] == 0

    rec = _only(out)
    assert rec.cigarstring == CIGAR_TRUE, 'a worse original must not be taken'
    assert rows[0]['reverted_to_original'] is False


def test_reverted_record_resolves_the_chromosome_by_name(tmp_path):
    """The pre-correction BAM's @SQ order is the REVERSE of the input's, so a
    record carried across by reference_id would land on chrB. from_dict()
    re-resolves `ref_name` by NAME (the ISSUE-001 failure class)."""
    assert HDR_IN.references != HDR_ORIG.references, 'fixture must differ'
    _rows, stats, out = _run(tmp_path, CIGAR_SLID, CIGAR_TRUE,
                             annotated={J_TRUE})
    assert stats['orig_accepted'] == 1
    rec = _only(out)
    assert rec.reference_name == 'chrA'
    assert rec.reference_start == START


def test_a_chromosome_the_input_header_lacks_is_skipped_and_counted(tmp_path):
    """A pre-correction BAM naming a contig the corrected BAM does not have
    (chr5 vs 5, an extra decoy) must be skipped, never coerced.

    The guard is load-bearing, not belt-and-braces: `from_dict` does NOT raise
    on an unknown `ref_name`. htslib warns ("unrecognized reference name ...;
    treated as unmapped") and hands back a record with reference_id -1 and
    is_unmapped True, which would then be scored and could be written into the
    output BAM. Verified by disabling the guard: this test fails
    (orig_skipped_unknown_chrom 0 != 1) rather than erroring.
    """
    inp = _bam(tmp_path / 'corrected.bam', HDR_IN, [_read(HDR_IN, CIGAR_SLID)])
    other = pysam.AlignmentHeader.from_dict(
        {'HD': {'VN': '1.6'}, 'SQ': [{'SN': '5', 'LN': GLEN}]})
    orig = _bam(tmp_path / 'original.bam', other,
                [_read(other, CIGAR_TRUE, chrom='5')])
    out = tmp_path / 'triaged.bam'
    _rows, stats = triage_realign_bam(
        inp, str(out), genome=GENOME, annotated_junctions={J_TRUE},
        original_bams=[orig], sort_and_index=False)
    assert stats['junction_leg'] == 1, 'the read must reach the leg to be skipped BY IT'
    assert stats['orig_skipped_unknown_chrom'] == 1
    # skipped BEFORE it becomes a candidate — not offered, not proposed
    assert stats['orig_leg'] == 0 and stats['orig_proposed'] == 0
    assert stats['orig_accepted'] == 0
    rec = _only(out)
    assert rec.cigarstring == CIGAR_SLID
    assert rec.reference_name == 'chrA' and not rec.is_unmapped


def test_without_original_bams_nothing_changes(tmp_path):
    """Default-off identity: the counters exist and stay zero, and the damaged
    read is NOT reverted — this is the bug ISSUE-009 reports."""
    inp = _bam(tmp_path / 'corrected.bam', HDR_IN, [_read(HDR_IN, CIGAR_SLID)])
    out = tmp_path / 'triaged.bam'
    rows, stats = triage_realign_bam(
        inp, str(out), genome=GENOME, annotated_junctions={J_TRUE},
        sort_and_index=False)
    assert stats['junction_leg'] == 1
    assert stats['orig_leg'] == stats['orig_proposed'] == stats['orig_accepted'] == 0
    assert _only(out).cigarstring == CIGAR_SLID
    assert rows[0]['reverted_to_original'] is False


def test_only_junction_leg_reads_are_offered_a_candidate(tmp_path):
    """A high-confidence read is bypassed; the leg must not widen the
    distressed minority it acts on."""
    inp = _bam(tmp_path / 'corrected.bam', HDR_IN, [_read(HDR_IN, CIGAR_TRUE)])
    orig = _bam(tmp_path / 'original.bam', HDR_ORIG,
                [_read(HDR_ORIG, CIGAR_SLID)])
    out = tmp_path / 'triaged.bam'
    rows, stats = triage_realign_bam(
        inp, str(out), genome=GENOME, annotated_junctions={J_TRUE},
        original_bams=[orig], sort_and_index=False)
    assert stats['triaged'] == 0 and stats['junction_leg'] == 0
    assert stats['orig_leg'] == 0
    assert _only(out).cigarstring == CIGAR_TRUE
    assert rows[0]['label'] == 'high_confidence'


def test_cli_exposes_original_bams():
    from rectify.cli import create_parser
    args = create_parser().parse_args(
        ['triage', 'in.bam', '-o', 'out', '--Scer',
         '--original-bams', 'a.bam', 'b.bam'])
    assert args.original_bams == ['a.bam', 'b.bam']
