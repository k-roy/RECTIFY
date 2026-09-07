"""T1 replay fixture — the two-tier sha's T1 review bundle (Sumner human RNA004 DRS), 2026-09-07.

Real reads from `i020g_372d6c5` T1 (all chromosomes) replayed through the production path on genome slices
(``tests/_sumner_review_bundle``: GENCODE introns of the slice as the annotated candidates, no pool). Collaborator
data, kept outside the repository (``dev/sumner_misplaced_panel_20260904/holdout/events/372d6c5/``); skipped when absent.

Pinned here — ISSUE-032(c)(i), the annotated-shift margin: on T1 seven baseline-true annotated rescues re-landed on
a NOVEL site 2-5 nt away because the per-candidate shift sweep ranks on the anchored deficit. The unslid annotated
placement now holds unless the shifted winner beats it by ANNOTATED_SHIFT_MARGIN (6) bits:
  7f779873  annotated `2S17M1D1M2D3M` 24.5 bits vs the -3 shift (a GC donor) 26.0  -> annotated
  5d30f4ea  annotated `9M3D13M1I` 30.0 vs +4 33.0                                 -> annotated
  d317bbe2  annotated 24.5 (the sweep's +5 was 22.5)                                -> annotated
  ed3301ff  annotated 30.5 vs +4 `19M2D8M2D3M` 38.0 (margin 7.5 >= 6)               -> the +4 stays (ruling card)
"""
import pytest

from tests import _sumner_review_bundle as RB

pytestmark = pytest.mark.skipif(not RB.bundle_present('372d6c5'),
                                reason='Sumner 372d6c5 T1 review bundle not present (collaborator data, kept outside the repo)')


def _replay(read8, monkeypatch, gate='report'):
    entry = RB.load_bundle('372d6c5')[read8]
    row, res, rec, stock = RB.replay(entry, monkeypatch, gate)
    return row, res, rec, stock, entry['off']


@pytest.mark.parametrize('read8,junction,cigar,bits', [
    ('7f779873', (74203710, 74203917), '2S17M1D1M2D3M', 24.5),
    ('5d30f4ea', (58910081, 58910333), '9M3D13M1I', 30.0),
    ('d317bbe2', (18589683, 18590853), None, 24.5),
])
def test_annotated_placement_holds_within_the_shift_margin(read8, junction, cigar, bits, monkeypatch):
    row, res, rec, stock, off = _replay(read8, monkeypatch)
    assert res.get('rescued') and res.get('landing_annotated') is True, res
    rj = res['rescued_junction']
    assert (rj[1] + off, rj[2] + off) == junction, (rj[1] + off, rj[2] + off)
    assert res.get('exon_bits') == bits, res.get('exon_bits')
    if cigar:
        assert res.get('five_prime_exon_cigar') == cigar, res.get('five_prime_exon_cigar')
    assert row['five_prime_rescue_refused'] == '', row['five_prime_rescue_refused']
    assert (junction[0] - off, junction[1] - off) in RB.nops(rec)


def test_a_shift_that_wins_by_more_than_the_margin_is_kept(monkeypatch):
    """ed3301ff: the +4 donor's block (38.0) beats the annotated one (30.5) by 7.5 bits — above the 6-bit margin —
    so the shifted (novel, GT) placement stands, as 3aea3e5a's clean +4 did. Kevin's queue decides whether a
    +4 block with two 2-base deletions is that kind of evidence."""
    row, res, rec, stock, off = _replay('ed3301ff', monkeypatch)
    assert res.get('rescued') and res.get('landing_annotated') is False, res
    rj = res['rescued_junction']
    assert rj[1] + off == 1223176 and res.get('exon_bits') == 38.0, (rj[1] + off, res.get('exon_bits'))


def test_margin_is_env_tunable(monkeypatch):
    """With a 10-bit margin ed3301ff's +4 (7.5 above the annotated block) is held on the annotated coordinate."""
    from rectify.core.splice.splice_aware_5prime import ANNOTATED_SHIFT_MARGIN, annotated_shift_margin
    assert ANNOTATED_SHIFT_MARGIN == 6.0 and annotated_shift_margin() == 6.0
    monkeypatch.setenv('RECTIFY_2F_ANNOTATED_SHIFT_MARGIN', '10')
    assert annotated_shift_margin() == 10.0
    row, res, rec, stock, off = _replay('ed3301ff', monkeypatch)
    assert res.get('rescued') and res.get('landing_annotated') is True, res
    assert res['rescued_junction'][1] + off == 1223172


def test_a_peel_that_would_shift_the_acceptor_reports_what_is_drawn_or_refuses(monkeypatch):
    """bcd90cad: 11S64M starting exactly at the annotated acceptor 8683201. The terminal peel took four clean
    exon-2 body bases into the exon block (`2S10M2D3M1D`, 12.5 bits); the writer lengthens the N by four, so the
    drawn acceptor would be 8683205 while the TSV named 8683201 (the baseline hid this behind the flat-M
    fallback). The peel now reports the shifted acceptor and judges it at its own tier: a novel, non-canonical
    acceptor at 12.5 bits is refused, the read keeps its record, and no junction is invented (TSV == BAM)."""
    row, res, rec, stock, off = _replay('bcd90cad', monkeypatch)
    assert not res.get('rescued'), res
    assert row['five_prime_rescue_refused'] != ''
    assert RB.nops(rec) == RB.nops(stock)
    assert (8682816 - off, 8683201 - off) not in RB.nops(rec) and (8682816 - off, 8683205 - off) not in RB.nops(rec)
