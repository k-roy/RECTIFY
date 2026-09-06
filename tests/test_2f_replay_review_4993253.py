"""ISSUE-028 / invariant E — replay of the five reads Kevin adjudicated in the 4993253 T0 review (2026-09-06).

The fixer's review bundle (``dev/sumner_misplaced_panel_20260904/holdout/events/4993253/review_2f_4993253/``:
stock.bam + arm BAMs + manifest, hg38 slices + the GENCODE records inside them) replayed through the PRODUCTION
path — ``correct_read_3prime`` with the slice's annotated introns as candidates (no pool file exists for this
bundle; every reviewed rescue landed on an annotated junction and the replay reproduces the T0 arms N-op for
N-op), then the corrected-BAM writer. Collaborator data outside the repository; skips when absent.

Kevin's verdicts (ISSUE-026 § "Kevin's review"; docs/algorithms/read_level_review.md):

  26f8fb45  control newly rescued   WRONG    32-bp block 9=/23X (28 %), MAPQ 1 -> must draw NO rescue
  22f609c6  control newly rescued   RIGHT    1I15M = 13=/2X at the tip, annotated GT-AG -> keeps its junction
  04b17fc6  control newly rescued   WRONG    26-bp block 11=/15X (42 %) -> must draw NO rescue
  6a6e0b3c  vanished rescue_annot   unsupported baseline (`8I6M`: a leading I is a clip) -> stays unrescued
  38722d08  vanished rescue_annot   unsupported baseline (`2M4I2M1I2M`) -> stays unrescued
"""
import pytest

from tests import _sumner_review_bundle as RB

pytestmark = pytest.mark.skipif(not RB.bundle_present('4993253'),
                                reason='Sumner 4993253 review bundle not present (collaborator data, kept outside the repo)')

WRONG_CONTROLS = ['26f8fb45', '04b17fc6']
UNSUPPORTED_BASELINE = ['6a6e0b3c', '38722d08']


def _replay(read8, monkeypatch, gate='report'):
    entry = RB.load_bundle('4993253')[read8]
    row, res, rec, stock = RB.replay(entry, monkeypatch, gate)
    return row, res, rec, stock, entry


@pytest.mark.parametrize('gate', ['report', 'refuse'])
@pytest.mark.parametrize('read8', WRONG_CONTROLS)
def test_the_two_wrong_controls_draw_no_rescue(read8, gate, monkeypatch):
    """The placed block is not evidence (identity 0.28 / 0.42 as drawn; the deepest placement 2F judged is under
    the identity floor too): a placement decision in both gate modes, the token and the shape in the TSV, no new
    N-op, the aligner's own junctions untouched."""
    from rectify.core.splice.splice_aware_5prime import EVIDENCE_REFUSALS, EXON_IDENTITY_REFUSAL
    row, res, rec, stock, entry = _replay(read8, monkeypatch, gate)
    assert not res.get('rescued'), res
    assert RB.nops(rec) == RB.nops(stock), (RB.real(RB.nops(rec), entry['off']), RB.real(RB.nops(stock), entry['off']))
    assert row['five_prime_rescue_refused'] == EXON_IDENTITY_REFUSAL, row['five_prime_rescue_refused']
    assert row['five_prime_rescue_refused'] in EVIDENCE_REFUSALS
    assert row['five_prime_exon_identity'] is not None and row['five_prime_exon_identity'] < 0.8
    assert row['five_prime_exon_bits'] is not None and row['five_prime_exon_bits'] < 12


@pytest.mark.parametrize('gate', ['report', 'refuse'])
def test_22f609c6_keeps_its_annotated_junction(gate, monkeypatch):
    """The right rescue: `1I15M` = 13=/2X at the tip, identity 0.87, 26 - 4 - 2.5 = 19.5 bits — evidence at the
    default floor. The writer draws exactly that junction."""
    row, res, rec, stock, entry = _replay('22f609c6', monkeypatch, gate)
    off = entry['off']
    assert res['rescued'] and res['landing_annotated'] is True
    assert res['rescued_junction'][1:] == (149487098 - off, 149490313 - off), RB.real([res['rescued_junction'][1:]], off)
    assert row['five_prime_rescue_refused'] == ''
    new = [n for n in RB.nops(rec) if n not in RB.nops(stock)]
    assert RB.real(new, off) == [(149487098, 149490313)], RB.real(new, off)
    assert row['five_prime_exon_cigar'] == '1I15M'
    assert row['five_prime_exon_identity'] == pytest.approx(13 / 15)
    assert row['five_prime_exon_bits'] == 19.5


@pytest.mark.parametrize('read8', UNSUPPORTED_BASELINE)
def test_unsupported_baselines_stay_unrescued(read8, monkeypatch):
    row, res, rec, stock, entry = _replay(read8, monkeypatch)
    assert not res.get('rescued'), res
    assert RB.nops(rec) == RB.nops(stock)


def test_6a6e0b3c_is_refused_on_the_bits_after_the_leading_i_is_stripped(monkeypatch):
    """The baseline's `8I6M` was six placed bases. The deepest placement 2F now judges (`8I5M1D` -> `8S5M1D`)
    has 5 matched bases and a deletion at the junction: 10 - 2.5 = 7.5 bits, `exon_bits_below_floor` — not the
    indel-burden token the leading I used to trip."""
    from rectify.core.splice.splice_aware_5prime import EXON_BITS_REFUSAL
    row, res, rec, stock, entry = _replay('6a6e0b3c', monkeypatch)
    assert row['five_prime_rescue_refused'] == EXON_BITS_REFUSAL, row['five_prime_rescue_refused']
    assert row['five_prime_exon_bits'] == 7.5 and row['five_prime_exon_identity'] == 1.0
