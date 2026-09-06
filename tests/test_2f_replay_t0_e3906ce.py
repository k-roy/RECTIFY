"""T0 replay fixture — the tester's e3906ce chrX regression bundle (Sumner human RNA004 DRS), ISSUE-026.

Ten real reads that regressed on `todo/issue-017-off-by-4` @ e3906ce against the
6485226 baseline, replayed through the PRODUCTION path — `correct_read_3prime`
with the tester's pool ∪ GENCODE candidates, then the corrected-BAM writer — on
genome slices (``tests/_sumner_replay_bundle.replay``). The records are
collaborator data and live OUTSIDE the repository
(``dev/sumner_misplaced_panel_20260904/holdout/events/e3906ce/``); the module
skips when they are absent.

Invariants pinned here (ISSUE-026 triage, `T0_E3906CE_TRIAGE_LOG.md`):

  A  a slide off an annotated coordinate is a novel placement and must carry
     novel-site evidence, in BOTH gate modes; an out-of-window shift may not land
     on a non-canonical donor when the unslid coordinate is canonical. The three
     chrX:154398394 reads and dab60caa land on the annotated donor and the writer
     draws it (2586f261 used to become a writer-refused +5 slide).
  B  a 5' rescue may never remove an aligner-called N-op: 771560c4 keeps its own
     15827374-15831447 and gains the annotated 15831568-15845378, never the
     exon-skipping 15827372-15845378.
  C  the indel-burden bound applies to annotated landings too: 1964c591's
     401-nt clip (exon block 255M/148I/23D) is not drawn; the read stays the
     baseline `softclipped_no_junction`.
  -  the zero-clip trio (ec95eb3b / 5f675f8a / 528315b8): a `1M` across an
     intron on a 1-nt acceptor ambiguity is NOT a rescue (orchestrator decision,
     ISSUE-026); no junction is drawn.
"""
import pytest

from tests import _sumner_replay_bundle as SB

pytestmark = pytest.mark.skipif(not SB.bundle_present('e3906ce'),
                                reason='Sumner T0 e3906ce replay bundle not present (collaborator data, kept outside the repo)')

CHRX_LOCUS = ['fd8c2b85', '14be8590', '2586f261']
ZERO_CLIP = ['ec95eb3b', '5f675f8a', '528315b8']


def _replay(read8, monkeypatch, gate='report'):
    """(row, raw 2F result, materialized record, stock record, offset)."""
    entry = SB.load_bundle('e3906ce')[read8]
    row, res, rec, stock = SB.replay(entry, monkeypatch, gate)
    return row, res, rec, stock, entry['off']


def _nops(read):
    return SB.nops(read)


def _real(nops, off):
    return SB.real(nops, off)


# ----------------------------------------------------------------------------------------------- invariant A
@pytest.mark.parametrize('gate', ['report', 'refuse'])
@pytest.mark.parametrize('read8', CHRX_LOCUS)
def test_chrx_locus_lands_on_the_annotated_donor_in_both_gate_modes(read8, gate, monkeypatch):
    """ISSUE-026 (a): the 11-nt clip's best window is 4-5 nt into the intron (the
    GTGAGT +4/+5 GT) on the SAME annotated candidate; that slide is a novel
    placement without evidence (9 matched) or onto a non-canonical donor (+5, TA)
    and is refused as a placement decision — the read lands on 154398394 and the
    writer draws it."""
    row, res, rec, stock, off = _replay(read8, monkeypatch, gate)
    assert res['rescued'], res
    assert res['rescued_junction'][1:] == (154398394 - off, 154398496 - off), (read8, _real([res['rescued_junction'][1:]], off))
    assert res['landing_annotated'] is True and res['novel_evidence'] == ''
    assert row['five_prime_rescue_refused'] == '', row['five_prime_rescue_refused']
    nops = _real(_nops(rec), off)
    assert (154398394, 154398496) in nops, nops
    assert not any(s in (154398398, 154398399) for s, _e in nops), nops   # never the +4 / +5 decoy


def test_dab60caa_keeps_its_annotated_junction(monkeypatch):
    """ISSUE-026 (a), dab60caa: the 52-nt clip's best window slid the annotated
    donor 5 nt LEFT onto an AG (ED 14.5/52); the writer then refused the
    non-canonical destination and the read LOST its annotated junction. The
    slide is refused up front; the annotated 103376613-103377027 stays drawn."""
    row, res, rec, stock, off = _replay('dab60caa', monkeypatch)
    nops = _real(_nops(rec), off)
    assert (103376613, 103377027) in nops, nops
    assert (103377106, 103377509) in nops, nops
    assert row['five_prime_rescue_refused'] == '', row['five_prime_rescue_refused']
    assert res['landing_annotated'] is True


def test_slide_refusal_is_a_placement_decision_not_a_token():
    """The slide refusal never reaches the TSV as a novel-evidence token on the
    fallback placement (the emitted, annotated placement earns ''); the
    constant is registered with the wrapper's trace set so a NON-fallback
    refusal is still attributable."""
    from rectify.core.splice.splice_aware_5prime import (
        ANNOTATED_SLIDE_REFUSAL, NOVEL_EXON_REFUSALS, PLACEMENT_REFUSALS)
    assert ANNOTATED_SLIDE_REFUSAL in PLACEMENT_REFUSALS
    assert ANNOTATED_SLIDE_REFUSAL not in NOVEL_EXON_REFUSALS
    assert set(NOVEL_EXON_REFUSALS) <= set(PLACEMENT_REFUSALS)


# ----------------------------------------------------------------------------------------------- invariant B
@pytest.mark.parametrize('gate', ['report', 'refuse'])
def test_771560c4_keeps_its_own_junction_and_gains_the_annotated_one(gate, monkeypatch):
    """ISSUE-026 (b-i): three annotated candidates share the far site 15845378;
    the 13-nt clip scores the identical window on all three. The nearest near
    site (15831568, 1 nt inside) must win — 15827372's intron contains the
    read's own 15827374-15831447 and its terminal exon, and drawing it deleted
    that junction (exon-skip collapse). Both aligner N-ops survive, the
    annotated 15831568-15845378 is added, the collapse junction never appears."""
    row, res, rec, stock, off = _replay('771560c4', monkeypatch, gate)
    nops = _real(_nops(rec), off)
    assert (15827374, 15831447) in nops, nops          # the aligner's own junction, kept
    assert (15831568, 15845378) in nops, nops          # the annotated one, gained
    assert (15827372, 15845378) not in nops, nops      # never the exon-skip collapse
    assert all(n in nops for n in _real(_nops(stock), off)), (nops, _real(_nops(stock), off))
    assert res['rescued'] and res['rescued_junction'][1:] == (15831568 - off, 15845378 - off)
    assert row['five_prime_rescue_refused'] == ''


def test_a_candidate_containing_the_reads_own_junction_is_never_ranked(monkeypatch):
    """The sequence loop skips 15827372-15845378 (and 15828200-15845378 is not
    the read's junction either): the counter fires and the winner is the
    nearest near site."""
    from rectify.core.splice.splice_aware_5prime import _OI_COUNTERS
    before = _OI_COUNTERS.get('five_prime_candidate_contains_read_nop', 0)
    row, res, rec, stock, off = _replay('771560c4', monkeypatch)
    assert _OI_COUNTERS.get('five_prime_candidate_contains_read_nop', 0) > before
    assert res['rescued_junction'][1:] == (15831568 - off, 15845378 - off)


# ----------------------------------------------------------------------------------------------- invariant C-minimal
@pytest.mark.parametrize('gate', ['report', 'refuse'])
def test_1964c591_annotated_landing_fails_the_indel_burden_bound(gate, monkeypatch):
    """ISSUE-026 (c): the 401-nt clip's anchored placement on the annotated
    314-kb 76693087-77007275 carries an exon block whose I+D burden exceeds the
    allowance (the block overruns the 105-bp exon by 173 nt into the next
    annotated intron). The junction is TRUE — a three-junction 5' extension,
    iteration-4 material — but this block is not evidence for it: the landing is
    refused as a placement in both gate modes, no 76693087-77007275 N-op is
    drawn, the aligner's own junctions survive, and the TSV names the refusal."""
    from rectify.core.splice.splice_aware_5prime import ANNOTATED_INDEL_BURDEN_REFUSAL
    row, res, rec, stock, off = _replay('1964c591', monkeypatch, gate)
    nops = _real(_nops(rec), off)
    assert (76693087, 77007275) not in nops, nops
    assert all(n in nops for n in _real(_nops(stock), off)), (nops, _real(_nops(stock), off))
    tokens = (row['five_prime_rescue_refused'], row.get('five_prime_novel_evidence') or '')
    assert any(ANNOTATED_INDEL_BURDEN_REFUSAL in t for t in tokens), (tokens, res.get('rescue_type'))


def test_annotated_indel_burden_refusal_is_registered():
    from rectify.core.splice.splice_aware_5prime import (
        ANNOTATED_INDEL_BURDEN_REFUSAL, NOVEL_EXON_REFUSALS, PLACEMENT_REFUSALS)
    assert ANNOTATED_INDEL_BURDEN_REFUSAL in PLACEMENT_REFUSALS
    assert ANNOTATED_INDEL_BURDEN_REFUSAL not in NOVEL_EXON_REFUSALS


# ----------------------------------------------------------------------------------------------- decision: zero-clip trio
@pytest.mark.parametrize('read8', ZERO_CLIP)
def test_zero_clip_one_base_acceptor_ambiguity_is_not_a_rescue(read8, monkeypatch):
    """The read's first base (G) is both the intron's last base and exon 1's last
    base; the baseline drew `1M<intron>N` from a 1-mer at ED 0 == ED 0. That is
    not a rescue (ISSUE-026 decision): the record keeps the aligner's N-ops
    exactly and gains none."""
    row, res, rec, stock, off = _replay(read8, monkeypatch)
    assert _nops(rec) == _nops(stock), (_real(_nops(rec), off), _real(_nops(stock), off))
    assert not res.get('rescued'), res
