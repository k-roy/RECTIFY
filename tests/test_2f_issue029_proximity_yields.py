"""ISSUE-029 — Case 3 (proximity) must not outrank a SCORED clip (2026-09-06).

A read whose 5' aligned base sits at (or within ``junction_proximity_bp`` of) an annotated exon-2 start used to
fall through to Case 3 whenever the sequence loop emitted nothing: the row was labelled ``proximity``, NAMED the
nearest annotated intron without aligning the clip, and reported ``landing_annotated = False`` for an annotated
coordinate. For an INFORMATIVE clip (>= ``min_informative_clip_bp()``) the sequence loop has already ranked the
clip against every candidate with the anchored placement model; if no landing met the floors (invariants A / C /
E) the read is a no-rescue that carries its refusal token — not a proximity row. Case 3 stays for the zero /
sub-floor clip (no sequence test exists) and for the rare clip that IS evidence at the named donor.

Kevin's pills (2026-09-06, relayed by the coordinator): the +4 placement ``8M3I1M`` for f53d770 5cef5ebb is
DISAPPROVED ("too many mismatches at that overhang"), likewise 975638b6 — so the replay pins the MECHANISM (no
proximity row, nothing drawn, the refusal named), never a +4 coordinate. Every landing still meets invariant E.
"""
import pysam
import pytest

import rectify.core.bam.bam_processor as bp
from rectify.core.splice.splice_aware_5prime import (
    PLACEMENT_REFUSALS,
    _OI_COUNTERS,
    min_informative_clip_bp,
    rescue_3ss_truncation,
)
from tests import _sumner_replay_bundle as SB
from tests.test_2f_evidence_shape import GENOME, GENOME_SEQ, JUNCTION, _clip_read

# Two annotated candidates share the exon-2 start at 140: the real one (40, 140) and a second intron whose donor
# is 30 nt upstream (10, 140) — exon-1 sequence there is the T-homopolymer, so an aperiodic clip cannot fit it.
ANNOTATED = {JUNCTION, ('chrT', 10, 140)}


def _junctions_of(read):
    out, pos = [], read.reference_start
    for op, ln in read.cigartuples or []:
        if op == 3:
            out.append((pos, pos + ln))
        if op in (0, 2, 3, 7, 8):
            pos += ln
    return out


# ------------------------------------------------------------------------------------------------ hermetic
@pytest.mark.parametrize('gate', ['report', 'refuse'])
def test_a_clean_informative_clip_at_the_exon_start_is_a_sequence_rescue_not_proximity(gate, monkeypatch):
    """The read starts exactly at the annotated exon-2 start (dist 0) with a 12-nt clip that is exon-1's tail:
    the anchored ranking places it (``12M``, 24 bits) — ``softclip``, annotated landing, junction drawn."""
    monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', gate)
    clip = GENOME_SEQ[28:40]
    res = rescue_3ss_truncation(_clip_read(clip), GENOME, ANNOTATED, '+', annotated_junctions=ANNOTATED)
    assert res['rescued'] and res['rescue_type'] == 'softclip', res
    assert res['rescued_junction'] == JUNCTION and res['landing_annotated'] is True
    assert res['five_prime_exon_cigar'] == '12M' and res['exon_bits'] == 24.0
    row = bp.correct_read_3prime(_clip_read(clip), GENOME, annotated_junctions=ANNOTATED)[0]
    assert row['five_prime_rescued'] and row['five_prime_rescue_refused'] == ''
    assert (40, 140) in [tuple(j) for j in row['junctions']]


def test_a_sub_floor_clip_at_the_exon_start_still_takes_case_3_and_names_the_annotated_intron():
    """A 6-nt clip is below ``min_informative_clip_bp()``: there is no sequence test, the proximity row stands, and
    ``landing_annotated`` says the named intron IS annotated (it used to be a constant False)."""
    clip = GENOME_SEQ[34:40]
    assert len(clip) < min_informative_clip_bp()
    res = rescue_3ss_truncation(_clip_read(clip), GENOME, ANNOTATED, '+', annotated_junctions=ANNOTATED)
    assert not res['rescued'] and res['rescue_type'] == 'proximity', res
    assert res['rescued_junction'] in ANNOTATED and res['landing_annotated'] is True
    # a pool-only (novel) neighbour is named with landing_annotated False
    res2 = rescue_3ss_truncation(_clip_read(clip), GENOME, {('chrT', 60, 140)}, '+', annotated_junctions=set())
    assert res2['rescue_type'] == 'proximity' and res2['landing_annotated'] is False, res2


@pytest.mark.parametrize('gate', ['report', 'refuse'])
def test_an_informative_clip_that_fits_no_candidate_is_a_refusal_not_a_proximity_row(gate, monkeypatch):
    """The 12-nt clip is the exon-1 tail with four isolated mismatches (identity 0.67): the ranking refuses every
    landing, Case 3 yields (the clip was scored), the row is ``none`` with the refusal token and the shape, and
    the writer draws nothing."""
    monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', gate)
    good = GENOME_SEQ[28:40]
    flip = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}
    noisy = ''.join(flip[c] if i in (1, 4, 7, 10) else c for i, c in enumerate(good))
    before = _OI_COUNTERS.get('five_prime_proximity_yields_to_scored_clip', 0)
    res = rescue_3ss_truncation(_clip_read(noisy), GENOME, ANNOTATED, '+', annotated_junctions=ANNOTATED)
    assert not res['rescued'] and res['rescue_type'] == 'none', res
    assert res['rescued_junction'] is None
    assert res.get('clip_refused') in PLACEMENT_REFUSALS, res
    assert res['exon_bits'] is not None and res['exon_bits'] < 18
    assert _OI_COUNTERS.get('five_prime_proximity_yields_to_scored_clip', 0) > before
    row = bp.correct_read_3prime(_clip_read(noisy), GENOME, annotated_junctions=ANNOTATED)[0]
    assert not row['five_prime_rescued'] and row['five_prime_rescue_refused'] in PLACEMENT_REFUSALS
    assert (40, 140) not in [tuple(j) for j in row['junctions']]


# ------------------------------------------------------------------------------------------------- replay
@pytest.mark.skipif(not SB.bundle_present('f53d770'),
                    reason='Sumner f53d770 replay bundle not present (collaborator data, kept outside the repo)')
@pytest.mark.parametrize('gate', ['report', 'refuse'])
@pytest.mark.parametrize('key', ['5cef5ebb', '5cef5ebb#2', '975638b6'])
def test_5cef5ebb_and_975638b6_are_refusals_not_proximity_rows_and_draw_nothing(key, gate, monkeypatch):
    """5cef5ebb: 12-nt clip GTATGGTGTACA, 5' base 154398497 (annotated exon-2 start + 1). Its best placement in
    ±6 of the annotated donor is the disapproved +4 (`8M3I1M`, 8 matched / 1 mismatch, 10.5 bits); the annotated
    donor gives `2M3I3M4I` (2.5 bits). 975638b6 (−, 10-nt clip TGTTTCGGGG): the annotated donor places `6M4S`
    (12.0 bits). Neither is evidence anywhere: no rescue, the refusal named, no proximity row, record == stock."""
    table = SB.load_bundle('f53d770')
    entry = table[key]
    row, res, rec, stock = SB.replay(entry, monkeypatch, gate)
    assert not res.get('rescued'), res
    assert res.get('rescue_type') != 'proximity', res
    assert res.get('rescued_junction') is None, res
    assert row['five_prime_rescue_refused'] in PLACEMENT_REFUSALS, row['five_prime_rescue_refused']
    assert row['five_prime_exon_bits'] is not None and row['five_prime_exon_bits'] < 18
    assert _junctions_of(rec) == _junctions_of(stock)
    assert not any(s + entry['off'] in (154398398, 154398399) for s, _e in _junctions_of(rec))
