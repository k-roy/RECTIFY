"""ISSUE-020 (arbiter RULING 10 §R37) — Module 2F RANKS 5' rescue candidates with
the ANCHORED placement model.

Mechanism (six-table, `dev/sumner_misplaced_panel_20260904/ISSUES/ISSUE-017_discriminating_bases_six.md`):
the ranking used to sweep `_shift` x `_off` and compare the 5' segment by hp-ED
to a genome window ending `_off` bases BEFORE the junction — an unpenalized
junction-side gap of up to junction_proximity_bp (10). The GTRAGT +5 GT decoy
4 nt into the intron won through that freedom (the six former "+4 by evidence"
reads only at `_off` = 2/10/7/3/3), and the placement (`align_clip_to_exon`,
junction end fixed) then had to spend a `4D` at the junction.

Now every candidate is scored, at its own coordinate and at the best few nearby
shifts, with the SAME Gotoh affine DP the placement runs; hp-ED survives only as
the shift prune; a 5' end inside the intron compares the clip PLUS the
intron-mapped bases (what the placement aligns); dist > 0 trims exon-2 bases from
the READ instead of sliding the genome window.

Consistency invariant (asserted here and, in debug mode
``RECTIFY_2F_CHECK_CONSISTENCY=1``, inside the code):
  (I1) re-scoring the winner's segment at the emitted junction reproduces the
       stored `anchored_deficit`;
  (I2) no anchored-scored candidate had a lower deficit;
  (I3) when the placement aligned the same segment, the emitted exon CIGAR
       scored with the four constants has the same deficit — same model, same
       answer.
"""
import random

import pysam
import pytest

from rectify.core.align import local_aligner as la
from rectify.core.splice import splice_aware_5prime as sa
from rectify.core.splice.overhang_informativeness import COUNTERS, reset_counters
from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation

# --------------------------------------------------------------------------- toy locus
# exon 1 [0, 60) ends ...CAG; intron [60, 200) = GTAAGT + filler + CAG (the +4 GT
# decoy sits at 64, the acceptor AG at 198); exon 2 [200, 300). The intron ENDS
# with CAG — the same three bases that end exon 1 — so an aligner extends the
# read body three bases into the intron: the six-table's "5' inside by 3 nt".
_RNG = random.Random(20260905)


def _filler(n):
    return ''.join(_RNG.choice('ACGT') for _ in range(n))


E1_TAIL = 'GACGTTGCATGCAGTCCATGACCGATTGGACAG'          # 33 nt, ends CAG
EXON1 = ('T' * (60 - len(E1_TAIL))) + E1_TAIL           # [0, 60)
INTRON = 'GTAAGT' + _filler(131) + 'CAG'                # [60, 200), 140 nt
EXON2 = _filler(100)                                    # [200, 300)
GENOME_SEQ = EXON1 + INTRON + EXON2
GENOME = {'chrT': GENOME_SEQ}
ANN = ('chrT', 60, 200)
DECOY = ('chrT', 64, 200)


def _read(name, start, cigar, seq, reverse=False):
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'},
                                           'SQ': [{'SN': 'chrT', 'LN': 3_000_000}]})
    r = pysam.AlignedSegment(hdr)
    r.query_name = name
    r.reference_name = 'chrT'
    r.reference_start = start
    r.cigartuples = cigar
    r.mapping_quality = 60
    r.query_sequence = seq
    r.is_reverse = reverse
    return r


@pytest.fixture(autouse=True)
def _fresh(monkeypatch):
    reset_counters()
    monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', 'report')
    monkeypatch.setenv('RECTIFY_2F_CHECK_CONSISTENCY', '1')   # the debug check raises on violation
    yield
    reset_counters()


def test_the_locus_is_what_the_tests_assume():
    assert GENOME_SEQ[57:60] == 'CAG' and GENOME_SEQ[197:200] == 'CAG'
    assert GENOME_SEQ[60:66] == 'GTAAGT' and GENOME_SEQ[64:66] == 'GT'
    assert GENOME_SEQ[198:200] == 'AG'
    assert len(GENOME_SEQ) == 300


# --------------------------------------------------------------------------- scorer identity
def test_score_right_anchored_equals_the_traceback_alignment_score():
    """(I3) at the unit level: the score-only right-anchored DP equals the score
    of the CIGAR `_align_right_anchored` traces back, on random pairs."""
    rng = random.Random(7)
    for _ in range(200):
        q = ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 30)))
        r = ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 36)))
        ops, ref_skip = la._align_right_anchored(q, r)
        score, state = la.score_right_anchored(q, r)
        assert score == la.affine_cigar_score(ops, q, r[ref_skip:]), (q, r, ops)
        assert state in ('H', 'D', 'I')
        if state == 'D':
            assert ops[-1][0] == la._OP_D


def test_score_left_anchored_equals_the_traceback_alignment_score():
    rng = random.Random(11)
    for _ in range(200):
        q = ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 30)))
        r = ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 36)))
        ops, consumed = la._align_left_anchored(q, r)
        score, j_best = la.score_left_anchored(q, r)
        assert consumed == j_best
        assert score == la.affine_cigar_score(ops, q, r[:consumed]), (q, r, ops)


def test_affine_cigar_score_refuses_ops_that_do_not_consume_the_sequences():
    with pytest.raises(ValueError):
        la.affine_cigar_score([(0, 3)], 'ACGT', 'ACGT')


# --------------------------------------------------------------------------- the decoy shape
def _inside_by_3_plus():
    """Plus-strand read: 15-nt clip = E1[-18:-3], body starts 3 bases INSIDE the
    intron (the CAG that ends the intron = the CAG that ends exon 1)."""
    clip = GENOME_SEQ[42:57]
    body = GENOME_SEQ[197:260]
    return _read('inside3_plus', 197, [(4, 15), (0, len(body))], clip + body)


def test_plus_inside_by_3_lands_annotated_with_no_junction_side_gap():
    read = _inside_by_3_plus()
    res = rescue_3ss_truncation(read, GENOME, {ANN, DECOY}, '+', annotated_junctions={ANN})
    assert res['rescued'] and res['rescue_type'] == 'softclip'
    assert res['rescued_junction'] == ANN
    assert res['landing_annotated'] is True and res['novel_evidence'] == ''
    # The compared segment is clip + the 3 intron-mapped bases = E1[-18:]: a
    # perfect anchored match, so the emitted exon CIGAR is a single M block
    # ending exactly at the donor — no D/I on the junction side.
    assert res['five_prime_exon_cigar'] == '18M'
    assert res['anchored_deficit'] == 0.0
    # (I2): the decoy, scored on the SAME segment with the junction end fixed,
    # is strictly worse — the +4 window can only be reached through a gap.
    seg = (GENOME_SEQ[42:57] + GENOME_SEQ[197:200]).upper()
    assert sa._anchored_deficit(seg, GENOME_SEQ, 64, '+') > res['anchored_deficit']
    assert sa._anchored_deficit(seg, GENOME_SEQ, 60, '+') == res['anchored_deficit']   # (I1)
    assert COUNTERS.get('five_prime_anchored_dps', 0) >= 2   # both candidates were anchored-scored


def test_the_decoy_alone_is_not_placed_through_a_junction_side_gap():
    """Only the +4 candidate on offer: the anchored rank has no window ending at
    the +4 donor that the segment fits (the old sweep found one `_off` = 4 bases
    away, hp-ED 0). Whatever the outcome, no sequence rescue places the segment
    with a gap at the junction, and nothing lands on the decoy by sequence."""
    read = _inside_by_3_plus()
    res = rescue_3ss_truncation(read, GENOME, {DECOY}, '+', annotated_junctions={ANN})
    if res.get('rescue_type') == 'softclip':
        cig = la.cigar_str_to_ops(res['five_prime_exon_cigar'])
        assert cig[-1][0] == la._OP_M
        assert res['novel_evidence'] != 'novel_exon_gap_at_junction'
        assert res['rescued_junction'] != DECOY


def test_minus_strand_mirror_lands_annotated_with_no_junction_side_gap():
    """The same locus reverse-complemented: the read is minus-strand, its clip
    is the trailing S, the intron-mapped bases sit at the segment's junction
    side (its first bases), the anchored DP is left-anchored."""
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rc = ''.join(comp[b] for b in reversed(GENOME_SEQ))
    n = len(GENOME_SEQ)
    genome = {'chrT': rc}
    ann = ('chrT', n - 200, n - 60)
    decoy = ('chrT', n - 200, n - 64)
    # exon-2 body (60 nt, genome-left) + the 3 intron-mapped bases + the clip (E1[-18:-3] in rc)
    body = rc[n - 260:n - 197]
    clip = rc[n - 57:n - 42]
    read = _read('inside3_minus', n - 260, [(0, len(body)), (4, 15)], body + clip, reverse=True)
    res = rescue_3ss_truncation(read, genome, {ann, decoy}, '-', annotated_junctions={ann})
    assert res['rescued'] and res['rescue_type'] == 'softclip'
    assert res['rescued_junction'] == ann
    assert res['five_prime_exon_cigar'] == '18M'
    assert res['anchored_deficit'] == 0.0
    seg = (rc[n - 200:n - 197] + clip).upper()
    assert sa._anchored_deficit(seg, rc, n - 64, '-') > 0.0
    assert sa._anchored_deficit(seg, rc, n - 60, '-') == 0.0


# --------------------------------------------------------------------------- dist > 0 is a trim
def test_dist_gt_0_trims_exon2_bases_from_the_read_not_the_genome_window():
    """The aligner started 2 bases into exon 2 and left those two exon-2 bases in
    the clip: clip = E1[-14:] + exon2[:2]. The ranking trims the 2 junction-side
    bases from the READ (segment = E1[-14:], deficit 0 at the annotated donor);
    the placement still aligns the whole 16-nt clip (the writer's contract), so
    the exon CIGAR spans 16 query bases and the debug check stays silent."""
    clip = GENOME_SEQ[46:60] + GENOME_SEQ[200:202]
    body = GENOME_SEQ[202:262]
    read = _read('dist2_plus', 202, [(4, 16), (0, len(body))], clip + body)
    res = rescue_3ss_truncation(read, GENOME, {ANN, DECOY}, '+', annotated_junctions={ANN})
    assert res['rescued'] and res['rescue_type'] == 'softclip'
    assert res['rescued_junction'] == ANN
    assert res['anchored_deficit'] == 0.0
    ops = la.cigar_str_to_ops(res['five_prime_exon_cigar'])
    assert sum(ln for op, ln in ops if op in (0, 1)) == 16
    assert COUNTERS.get('exon2_trim_consumed_clip', 0) == 0


# --------------------------------------------------------------------------- reranked-between-annotated counter
def test_reranked_between_annotated_counter_fires_when_the_models_disagree():
    """Two ANNOTATED candidates share the acceptor. Window X (donor at 96) is the
    segment with three single-base homopolymer insertions — hp-ED 1.5 (three 0.5
    HP deletions), affine deficit 15 (three gap events); window Y (donor at 60)
    is the segment with two substitutions — hp-ED 2.0, affine deficit 12. The
    pre-020 hp-ED rank picks X, the anchored rank picks Y: the counter counts
    the move, once per read, and the result is Y's junction."""
    rng = random.Random(3)
    seg = 'ACGTTGCAGACTACCGATGAGCTTAGGCAT'              # 30 nt, aperiodic; runs TT(3-4) CC(13-14) TT(22-23) GG(25-26)
    # Y = seg with three substitutions (6 C->T, 17 T->C, 24 A->C; none creates a run):
    # hp-ED 3.0, affine deficit 3 x 6 = 18.
    y = list(seg); y[6] = 'T'; y[17] = 'C'; y[24] = 'C'; y = ''.join(y)
    assert len(y) == len(seg) and sa._hp_edit_distance(seg, y) == 3.0
    # X = seg with two homopolymer DELETIONS (T at 4, T at 23) and two homopolymer
    # INSERTIONS (C after 14 -> CCC, G after 26 -> GGG): same length as seg, so the
    # physical window IS x; hp-ED 4 x 0.5 = 2.0 (cheaper than Y) while the affine
    # DP's optimum against x costs a deficit of 24 (measured; four gap events are
    # 20 on their own and the DP finds nothing cheaper) — dearer than Y's 18.
    # The two models disagree.
    x = seg[:4] + seg[5:15] + 'C' + seg[15:23] + seg[24:27] + 'G' + seg[27:]
    assert len(x) == len(seg) and sa._hp_edit_distance(seg, x) == 2.0
    fill = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    # layout: [0,40) filler + Y (30, ends at 70) | intron2 from 70: GTAAGT + filler to 83 |
    #         X (30) ends at 113 | intron1 from 113: GTAAGT + filler + CAG to 257 | exon 2
    g = fill(40) + y                                     # Y = genome[40:70], donor2 = 70
    g += 'GTAAGT' + fill(7)                             # [70, 83)
    g += x                                              # X = genome[83:113], donor1 = 113
    g += 'GTAAGT' + fill(135) + 'CAG'                   # intron1 [113, 257)
    g += fill(100)                                      # exon 2 [257, 357)
    assert g[40:70] == y and g[83:113] == x and g[255:257] == 'AG'
    genome = {'chrT': g}
    cand_x, cand_y = ('chrT', 113, 257), ('chrT', 70, 257)
    body = g[257:317]
    read = _read('rerank', 257, [(4, 30), (0, len(body))], seg + body)
    res = rescue_3ss_truncation(read, genome, {cand_x, cand_y}, '+',
                                annotated_junctions={cand_x, cand_y})
    assert res['rescued'] and res['rescue_type'] == 'softclip'
    dx = sa._anchored_deficit(seg, g, 113, '+')
    dy = sa._anchored_deficit(seg, g, 70, '+')
    assert dy == 18.0 and dx == 24.0 and dx > dy, (dx, dy)
    assert res['rescued_junction'] == cand_y
    assert res['anchored_deficit'] == 18.0
    assert res['edit_distance'] == 3.0                  # hp-ED of the chosen physical window is still reported
    assert res['reranked_between_annotated'] is True
    assert COUNTERS.get('five_prime_reranked_between_annotated', 0) == 1


def test_reranked_counter_is_silent_when_the_models_agree():
    read = _inside_by_3_plus()
    res = rescue_3ss_truncation(read, GENOME, {ANN, DECOY}, '+', annotated_junctions={ANN})
    assert res['rescued'] and res['reranked_between_annotated'] is False
    assert COUNTERS.get('five_prime_reranked_between_annotated', 0) == 0


# --------------------------------------------------------------------------- the debug check itself
def test_debug_consistency_check_raises_on_a_violated_invariant():
    seg = 'ACGTTGCATGCAGTCCATGACC'
    with pytest.raises(AssertionError):
        sa._check_anchored_consistency(seg.upper(), GENOME_SEQ, ANN, '+', -1.0, seg, None, None, 'x')
