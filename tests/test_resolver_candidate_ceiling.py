"""Pathological-contig circuit breaker for the overhang resolver.

Candidate enumeration is near_sites x far_sites. On a repetitive or
low-complexity contig that product explodes — 361 candidates/clip measured on
untrimmed cDNA (planning/681 CP6), and a reporter-construct contig ran ~555x
below baseline (planning/673) until an operator set RECTIFY_SKIP_REGIONS by
hand. Skip regions require knowing the bad contig IN ADVANCE, so an
unanticipated custom reference degrades to a silent multi-hour stall.

`max_candidates_per_clip` bounds the per-clip DP work: past the ceiling the
clip is abandoned, the read passes through UNTOUCHED, and
`refused_candidate_blowup` counts it.

The risk this guard must not become: a ceiling set too low silently refuses
legitimate clips — the same false-green class the guard exists to prevent. The
tests pin BOTH directions.

🔴 TWO FIXTURE TRAPS, both hit while writing this file — read before editing it:

1. A splice-site-free genome makes every test here pass VACUOUSLY at
   `candidates_evaluated=0`, never entering the loop under test.
   `test_fixture_actually_enumerates_candidates` exists to fail loudly if that
   returns.
2. A TANDEM-REPEAT genome is not usable either, even though it produces plenty
   of sites. The overhang gate caps W_max at period-1 for periodic sequence, so
   on a period-30 repeat the far window collapses below `min_intron` and
   enumeration yields nothing — silently, with no refusal counter. The fixture
   must not be built from the very sequence class the product defends against.
   (That version passed on a pre-gate branch and failed the moment it met
   master, which is how the trap surfaced.)

Hence the fixture below: an APERIODIC random backbone with splice sites
PLANTED deliberately — one acceptor at the alignment edge, ~300 donors upstream
at irregular offsets inside [min_intron, max_intron]. That yields ~780
candidates for one clip without any periodicity for the gate to catch.
"""

import logging
import random

import pysam

from rectify.core.align import overhang_resolver as ohr
from rectify.core.align.overhang_resolver import (
    ResolverConfig,
    ResolverStats,
    run_overhang_resolver,
)

_SEED = 20260812
_LEN = 12000
_EDGE = 9000  # alignment start; the LEFT clip attaches here


def _planted_genome(tmp_path):
    """Aperiodic backbone + deliberately planted splice sites.

    LEFT clip on the plus strand: the near site is an acceptor (AG) at the
    intron end, the far sites are donors (GT) upstream. Planting one acceptor
    flush with the alignment edge and many donors across the window gives a
    large near x far product — a blow-up — with no repeat structure.
    """
    rng = random.Random(_SEED)
    g = [rng.choice('ACGT') for _ in range(_LEN)]
    g[_EDGE - 2], g[_EDGE - 1] = 'A', 'G'          # acceptor at the edge
    off = 60
    while off < 4800:                               # inside max_intron (5000)
        g[_EDGE - off], g[_EDGE - off + 1] = 'G', 'T'
        off += rng.choice((7, 11, 13, 17, 19, 23))  # irregular => aperiodic
    seq = ''.join(g)
    (tmp_path / 'g.fa').write_text('>chrI\n' + seq + '\n')
    return seq


def _clipped_bam(tmp_path, seq):
    header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'SQ': [{'SN': 'chrI', 'LN': len(seq)}],
    })
    bam = tmp_path / 'in.bam'
    with pysam.AlignmentFile(bam, 'wb', header=header) as out:
        r = pysam.AlignedSegment(header)
        r.query_name = 'r1'
        r.reference_id = 0
        r.reference_start = _EDGE
        r.mapping_quality = 60
        r.query_sequence = seq[_EDGE - 60:_EDGE + 60]
        r.cigarstring = '60S60M'
        r.query_qualities = pysam.qualitystring_to_array('I' * 120)
        r.set_tag('RN', 1, value_type='i')
        out.write(r)
    return bam


def _run(tmp_path, ceiling, tag):
    ohr._BLOWUP_WARNED.clear()
    seq = _planted_genome(tmp_path)
    bam_in = _clipped_bam(tmp_path, seq)
    out = tmp_path / f'out_{tag}.bam'
    run_overhang_resolver(
        str(bam_in), str(tmp_path / 'g.fa'), str(out),
        config=ResolverConfig(max_candidates_per_clip=ceiling),
    )
    return run_overhang_resolver.last_stats, out


def test_fixture_actually_enumerates_candidates(tmp_path):
    """Anti-vacuity guard: without this, every test below could pass trivially."""
    stats, _ = _run(tmp_path, 100_000, 'loose')
    assert stats.clips_assessed == 1
    assert stats.candidates_evaluated > 100, (
        "fixture stopped producing candidates — the ceiling tests would be vacuous"
    )
    assert stats.refused_candidate_blowup == 0


def test_ceiling_bounds_the_work_and_counts_the_refusal(tmp_path):
    """The point of the guard: unbounded work becomes bounded work."""
    loose, _ = _run(tmp_path, 100_000, 'loose')
    tight, _ = _run(tmp_path, 5, 'tight')

    assert tight.refused_candidate_blowup == 1
    # Work stops at the budget (+1, the candidate that trips it) instead of
    # running to completion — that is the stall being converted into a bound.
    assert tight.candidates_evaluated <= 6
    assert tight.candidates_evaluated < loose.candidates_evaluated / 10


def test_abandoned_clip_is_passthrough_not_a_dropped_read(tmp_path):
    """Refusing a clip must never lose a record — the resolver is
    passthrough-or-rewrite over the same read set, and consensus downstream
    keys on RN:i."""
    stats, out = _run(tmp_path, 5, 'tight')
    with pysam.AlignmentFile(out, 'rb', check_sq=False) as f:
        recs = list(f)
    assert [r.query_name for r in recs] == ['r1']
    assert recs[0].get_tag('RN') == 1
    assert stats.resolved == 0  # nothing may be resolved off an abandoned clip


def test_default_ceiling_is_loose_enough_for_real_data():
    """Guard the constant against being tightened into a silent-loss setting.

    Post-trim cDNA runs ~6.7 candidates/clip; the worst legitimate figure ever
    observed is the 361/clip of the untrimmed pathology (planning/681 CP6). A
    default at or below that would refuse real work.
    """
    assert ResolverConfig().max_candidates_per_clip > 361


def test_default_ceiling_does_not_fire_on_this_contig(tmp_path):
    """The other direction: ~780 candidates is dense, but not pathological."""
    stats, _ = _run(tmp_path, ResolverConfig().max_candidates_per_clip, 'default')
    assert stats.refused_candidate_blowup == 0


def test_counter_reaches_the_stats_dict():
    """It must land in stats.json — that file is the beta-ledger record."""
    d = ResolverStats().as_dict()
    assert 'refused_candidate_blowup' in d
    assert d['refused_candidate_blowup'] == 0


def test_warning_is_throttled_to_once_per_contig(caplog):
    """A pathological contig yields one blow-up per clip — potentially millions.
    An unthrottled warning would become the performance problem it reports."""
    ohr._BLOWUP_WARNED.clear()
    caplog.set_level(logging.WARNING)
    for _ in range(50):
        ohr._warn_blowup('chrBad', 2000)
    assert caplog.text.count('candidate blow-up') == 1
    assert 'chrBad' in caplog.text
    ohr._warn_blowup('chrOther', 2000)  # a different contig still gets a line
    assert caplog.text.count('candidate blow-up') == 2
