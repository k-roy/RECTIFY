"""ISSUE-087: a forced 3' tail beats a soft clip only when it fits the genome better than chance (opt-in gate).

Consensus charges an RNA-3' soft clip but not an arm that writes the same bases as mismatches at A positions or
as deletions, so a forced tail always won; on the 2026-09-28 route smoke, tails of 16 nt or more that fit no
better than chance were the one class the landing lost (46 pin-better, 21 landing-better), and one tail was
written as a 5,875-nt deletion. With RECTIFY_TAIL_GATE=filed|written the forced arm is charged the clip arm's 3'
clip; unset, nothing changes.
"""

import random
from dataclasses import dataclass

import pysam
import pytest

from rectify.core.consensus import tail_gate as tg

CHROM = 'chrG'


def _genome(seed=87, n=4000):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


def _read(name, start, cigar, seq, reverse=False):
    header = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [{'SN': CHROM, 'LN': 4000}]})
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.reference_id = 0
    r.reference_start = start
    r.flag = 16 if reverse else 0
    r.mapping_quality = 60
    r.cigarstring = cigar
    r.query_sequence = seq
    r.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
    return r


@dataclass
class _Info:
    effective_three_prime_clip: int
    junction_score: float


def test_off_by_default(monkeypatch):
    monkeypatch.delenv('RECTIFY_TAIL_GATE', raising=False)
    assert tg.gate_mode() is None


def test_a_tail_that_fits_by_chance_is_gated_and_a_real_fit_is_not(monkeypatch):
    g = _genome()
    body = g[1000:1100]
    rng = random.Random(3)
    junk = ''.join(rng.choice('ACGT') for _ in range(24))        # unrelated to the genome after the body
    genome = {CHROM: g}
    clip_arm = _read('r1', 1000, '100M24S', body + junk)
    forced = _read('r1', 1000, '124M', body + junk)              # writes the junk tail onto the genome
    infos = {'minimap2': _Info(10, -20.0), 'uLTRA': _Info(0, 0.0)}
    gated = tg.apply_tail_gate(infos, {'minimap2': clip_arm, 'uLTRA': forced}, genome, mode='filed')
    assert 'uLTRA' in gated and gated['uLTRA'].startswith('chance_fail')
    assert infos['uLTRA'].effective_three_prime_clip == 10 and infos['uLTRA'].junction_score == -10.0
    # The same arms on the real downstream sequence: the tail IS the genome, so the extension stands.
    real = g[1100:1124]
    clip_arm = _read('r2', 1000, '100M24S', body + real)
    exact = _read('r2', 1000, '124M', body + real)
    infos = {'minimap2': _Info(10, -20.0), 'uLTRA': _Info(0, 0.0)}
    for mode in ('filed', 'written'):
        assert tg.apply_tail_gate(infos, {'minimap2': clip_arm, 'uLTRA': exact}, genome, mode=mode) == {}
    assert infos['uLTRA'].junction_score == 0.0


def test_a_large_deletion_inside_the_tail_always_gates():
    g = _genome()
    genome = {CHROM: g}
    body = g[1000:1100]
    tail = g[1100:1115] + g[2600]                                  # 15 nt that match, then one base 1,485 nt away
    clip_arm = _read('r3', 1000, '100M16S', body + tail)
    deleted = _read('r3', 1000, '115M1485D1M', body + tail)
    infos = {'minimap2': _Info(8, -16.0), 'uLTRA': _Info(0, 0.0)}
    gated = tg.apply_tail_gate(infos, {'minimap2': clip_arm, 'uLTRA': deleted}, genome, mode='filed')
    assert gated == {'uLTRA': 'tail_deletion_1485'}


def test_short_or_absent_clips_are_not_tested():
    g = _genome()
    genome = {CHROM: g}
    body = g[1000:1100]
    infos = {'minimap2': _Info(4, -8.0), 'uLTRA': _Info(0, 0.0)}
    clip_arm = _read('r4', 1000, '100M10S', body + 'ACGTACGTAC')
    forced = _read('r4', 1000, '110M', body + 'ACGTACGTAC')
    assert tg.apply_tail_gate(infos, {'minimap2': clip_arm, 'uLTRA': forced}, genome, mode='filed') == {}


def test_minus_strand_tail_is_read_in_rna_orientation():
    g = _genome()
    genome = {CHROM: g}
    body = g[1000:1100]
    rng = random.Random(5)
    junk = ''.join(rng.choice('ACGT') for _ in range(24))
    clip_arm = _read('r5', 1000, '24S100M', junk + body, reverse=True)      # RNA-3' end on the left
    forced = _read('r5', 976, '124M', junk + body, reverse=True)
    infos = {'minimap2': _Info(10, -20.0), 'uLTRA': _Info(0, 0.0)}
    gated = tg.apply_tail_gate(infos, {'minimap2': clip_arm, 'uLTRA': forced}, genome, mode='filed')
    assert 'uLTRA' in gated


def test_decisions_do_not_depend_on_shuffle_order(monkeypatch):
    g = _genome()
    genome = {CHROM: g}
    body = g[1000:1100]
    rng = random.Random(9)
    junk = ''.join(rng.choice('ACGT') for _ in range(30))
    runs = []
    for _ in range(2):
        infos = {'minimap2': _Info(10, -20.0), 'uLTRA': _Info(0, 0.0)}
        runs.append(tg.apply_tail_gate(infos, {'minimap2': _read('r6', 1000, '100M30S', body + junk),
                                               'uLTRA': _read('r6', 1000, '130M', body + junk)}, genome,
                                       mode='filed'))
    assert runs[0] == runs[1]


def test_selection_runs_the_gate_only_when_switched_on(monkeypatch):
    from rectify.core.consensus.extract import extract_alignment_info
    from rectify.core.consensus.select import select_best_alignment
    g = _genome()
    genome = {CHROM: g}
    body = g[1000:1100]
    rng = random.Random(3)
    junk = ''.join(rng.choice('ACGT') for _ in range(24))
    reads = {'minimap2': _read('r7', 1000, '100M24S', body + junk), 'uLTRA': _read('r7', 1000, '124M', body + junk)}
    for mode, expect in (('', {}), ('filed', {'uLTRA'})):
        monkeypatch.setenv('RECTIFY_TAIL_GATE', mode)
        infos = {a: extract_alignment_info(r, a, genome) for a, r in reads.items()}
        res = select_best_alignment(infos, genome, raw_reads=reads)
        assert set(res.tail_gated) == expect, mode
