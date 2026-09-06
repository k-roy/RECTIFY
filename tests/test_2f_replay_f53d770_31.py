"""ISSUE-026 invariant D — "the placement 2F reports is the placement the writer draws".

Replay of the tester's f53d770 T0 bundle (31 chrX reads, Sumner human RNA004 DRS: 11 VANISHED_TP_rescue_annot +
12 VANISHED_FP_added_nov + 4 d4_still_at_decoy_drawn + 4 d4_draws_nothing; three reads listed twice) through the
production path (``tests/_sumner_replay_bundle.replay``). Collaborator data outside the repo; skips when absent.

At f53d770 28/31 violated the invariant: the anchored rank emitted a junction whose donor is non-canonical
(TA-AG, GG-AG, CA-AG, GA-AG, ...; the deficit is the primary key, `_donor_ok` only a tie-breaker) and the writer's
canonical guard reverted it; and on the 15 reads whose alignment starts `dist` bases into exon 2 the writer's
`extend` ran the N to the live edge, so the drawn acceptor sat `dist` past the reported one (drawn WRONG at
6485226 on 9 of them, or pulled back <= 3 nt with a compensating D — the ISSUE-023 shape).

The invariant, asserted for EVERY read with no special case: for an accepted sequence rescue the materialized
record carries exactly one new N-op, equal to `rescued_junction`, with refusal '' and no writer-side acceptor
repair; a read 2F did not rescue gains no N-op. Both gate modes.
"""
import pytest

from tests import _sumner_replay_bundle as SB

pytestmark = pytest.mark.skipif(not SB.bundle_present('f53d770'),
                                reason='Sumner f53d770 31-read replay bundle not present (collaborator data, kept outside the repo)')

_KEYS = sorted(SB.load_bundle('f53d770')) if SB.bundle_present('f53d770') else []


@pytest.mark.parametrize('gate', ['report', 'refuse'])
@pytest.mark.parametrize('key', _KEYS)
def test_drawn_equals_reported_with_no_refusal(key, gate, monkeypatch):
    entry = SB.load_bundle('f53d770')[key]
    row, res, rec, stock = SB.replay(entry, monkeypatch, gate)
    new = SB.new_nops(rec, stock)
    if not res.get('rescued'):
        assert new == [], (key, SB.real(new, entry['off']))
        return
    rj = res['rescued_junction']
    assert row['five_prime_rescue_refused'] == '', (key, row['five_prime_rescue_refused'], SB.real([rj[1:]], entry['off']))
    assert new == [(rj[1], rj[2])], (key, 'reported', SB.real([rj[1:]], entry['off']), 'drawn', SB.real(new, entry['off']))
    assert not SB.body_side_d(rec, new, entry['strand'] if 'strand' in entry else row['strand']), (key, rec.cigarstring)
    # the TSV's junction list agrees with the record (ISSUE-024's invariant, on this read)
    assert (rj[1], rj[2]) in [tuple(j) for j in row['junctions']], (key, row['junctions'])


@pytest.mark.parametrize('key', _KEYS)
def test_no_emitted_junction_is_one_the_writer_reverts(key, monkeypatch):
    """`_writer_would_revert` is the writer's own canonical test asked before the placement is reported."""
    from rectify.core.splice.splice_aware_5prime import _writer_would_revert
    entry = SB.load_bundle('f53d770')[key]
    row, res, rec, stock = SB.replay(entry, monkeypatch)
    if res.get('rescued'):
        rj = res['rescued_junction']
        assert not _writer_would_revert(entry['seq'], rj[1], rj[2]), (key, SB.real([rj[1:]], entry['off']))


def test_the_vanished_reads_draw_again(monkeypatch):
    """The 23 true vanishes of f53d770 (11 TP + 12 FP by the tester's baseline classes) draw a junction again —
    at the coordinates 2F reports. The reads that are NOT drawn are named here so the residue stays visible:
    38722d08 (VANISHED_TP_rescue_annot; 11-nt clip, 1 nt into exon 2) has no canonical shift left once the
    writer-guard is applied to the sweep and ends as a proximity-only row; c5d1c111 (VANISHED_FP_added_nov, listed
    twice) is the +4 decoy 154398398 the baseline DREW on 9 matched bases — invariant A refuses that slide off the
    annotated 154398394 as a placement without novel-site evidence (`novel_exon_matched_below_floor`, both gate
    modes), the unslid coordinate does not pass the exon-vs-intron acceptance, and the read starts exactly at the
    acceptor so Case 4 cannot snap it; ac5225e1 and fb0cdd4e (both VANISHED_FP_added_nov: the baseline drew a wrong
    novel junction 4-5 nt off the annotated one) landed on the annotated site under invariant A with exon blocks of
    `4I6M3D` (7 gap bases on 6 matched) and `3M2I4M5D4M` (7 on 11) — invariant C's sanity bound refuses those
    blocks (`annotated_exon_indel_burden`), so they draw nothing rather than the baseline's false positive.
    Every one of the 11 VANISHED_TP_rescue_annot reads except 38722d08 draws again at its annotated junction."""
    table = SB.load_bundle('f53d770')
    not_drawn, tokens = [], {}
    n_vanished = 0
    for key, entry in sorted(table.items()):
        if not any(c.startswith('VANISHED') for c in entry['classes']):
            continue
        n_vanished += 1
        row, res, rec, stock = SB.replay(entry, monkeypatch)
        if not (res.get('rescued') and row['five_prime_rescue_refused'] == '' and SB.new_nops(rec, stock)):
            not_drawn.append(key)
            tokens[key] = row['five_prime_rescue_refused']
    assert n_vanished >= 23, n_vanished
    assert not_drawn == ['38722d08', 'ac5225e1', 'c5d1c111', 'c5d1c111#2', 'fb0cdd4e'], (not_drawn, tokens)
    assert tokens['c5d1c111'] == 'novel_exon_matched_below_floor', tokens
    assert tokens['ac5225e1'] == 'annotated_exon_indel_burden', tokens
    assert tokens['fb0cdd4e'] == 'annotated_exon_indel_burden', tokens   # a terminal-peel refusal, carried to the TSV
    for k in ('c5d1c111', 'ac5225e1', 'fb0cdd4e'):
        assert 'VANISHED_FP_added_nov' in table[k]['classes'], (k, table[k]['classes'])
    assert 'VANISHED_TP_rescue_annot' in table['38722d08']['classes']
    tp_not_drawn = [k for k in not_drawn if 'VANISHED_TP_rescue_annot' in table[k]['classes']]
    assert tp_not_drawn == ['38722d08'], tp_not_drawn
