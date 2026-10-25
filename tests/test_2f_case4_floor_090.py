"""ISSUE-090: a Case-4 intronic snap of fewer bases than the informative floor skipped the evidence floor (opt-in fix).

Kevin's ruling (read-review section 5): every landing, annotated or novel, meets the evidence floor; a clip below it
gets a most-parsimonious origin, never a junction. The Case-4 snap only applied the floor once at least
`min_informative_clip_bp()` bases were placed, so 083-26..29 (2-4 matched bases) were written at -1.5 to 4.5 bits,
one as a novel AT-AC junction inside an annotated RPS23 exon (master too). RECTIFY_2F_CASE4_FLOOR=1 applies the floor
to every snap; unset keeps today's acceptance.
"""
from tests.test_2f_novel_site_evidence_gate import GENOME, JUNCTION, _intronic_read

from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation


def _snap(depth):
    return rescue_3ss_truncation(_intronic_read(depth), GENOME, {JUNCTION}, '+', annotated_junctions={JUNCTION})


def test_default_keeps_the_short_annotated_snap(monkeypatch):
    monkeypatch.delenv('RECTIFY_2F_CASE4_FLOOR', raising=False)
    res = _snap(4)
    assert res['rescued'] and res['rescue_type'] == 'intronic_snap' and res['five_prime_exon_cigar'] == '4M'


def test_switch_refuses_a_snap_under_the_annotated_floor(monkeypatch):
    monkeypatch.setenv('RECTIFY_2F_CASE4_FLOOR', '1')
    res = _snap(4)                       # 4 matched bases = 8 bits < the 12-bit annotated floor
    assert not (res['rescued'] and res['rescue_type'] == 'intronic_snap')


def test_switch_keeps_a_snap_that_meets_the_floor(monkeypatch):
    monkeypatch.setenv('RECTIFY_2F_CASE4_FLOOR', '1')
    res = _snap(8)                       # 8 matched bases = 16 bits >= 12
    assert res['rescued'] and res['rescue_type'] == 'intronic_snap' and res['five_prime_exon_cigar'] == '8M'
