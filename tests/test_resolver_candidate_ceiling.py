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

    # Mutation-style contrast: the SAME input with the ceiling effectively
    # removed must behave differently, so the assertions below cannot pass for
    # any reason other than the guard firing. (Recommended by the cdna-trim-fix
    # session, which mutation-tested its own pileup no-flip guard the same way.)
    assert loose.refused_candidate_blowup == 0
    assert loose.candidates_evaluated > 100

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

    Measured (planning/681 CP6, real resolver, 2,000 real cDNA molecules):
        pre-trim-fix   17 clips ->  6,137 candidates  =  361 / clip
        post-trim-fix   9 clips ->     60 candidates  =  6.7 / clip

    🔴 361 is the PRE-FIX PATHOLOGY, not a healthy maximum. It is used here as a
    deliberately conservative floor, NOT as a calibration target. Production now
    lives at ~6.7/clip, and a future reader who "tightens toward the post-fix
    mean" would start refusing legitimate clips on repetitive or custom
    references — exactly the case the ceiling exists to survive. Percentiles
    re-derived from pre-fix data would describe a distribution that no longer
    exists.
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
        ohr._warn_blowup('chrBad', 2000, 2001, 5000)
    assert caplog.text.count('candidate blow-up') == 1
    assert 'chrBad' in caplog.text
    ohr._warn_blowup('chrOther', 2000, 2001, 5000)  # a different contig still logs
    assert caplog.text.count('candidate blow-up') == 2


# ---------------------------------------------------------------------------
# ISSUE-010 — the ceiling is a FLOOR that scales with the search window.
# ---------------------------------------------------------------------------

class TestWindowScaledCeiling:
    """`max_candidates_per_clip` was calibrated on yeast at max_intron=5000.
    Candidates are near_sites x far_sites and the far factor is the search
    WINDOW times the reference's site density, so the same healthy clip
    enumerates orders of magnitude more on a mammalian genome at a mammalian
    max_intron — without being pathological.

    Measured on human chr5 (GRCh38, 181.5 Mbp) over the 1,000 real ONT reads of
    `dev/sumner_misplaced_panel_20260904/holdout/chr5_holdout1k.bam`:

        max_intron    median cands/clip     max     abandoned at a FLAT 2000
        5,000 (yeast)              774    2,742     1/50  ( 2%)
        500,000 (human)         31,829  221,086    31/50  (62%)

    and the flat ceiling cost 5 of 7 genuine junctions on that data (removing it
    recovered all 5 at byte-identical coordinates and edit distances). The
    window-scaled ceiling takes that 62% back to 4% while still catching the
    221k-candidate outlier.

    Both directions are pinned: yeast-scale behaviour must be byte-identical,
    human-scale must actually scale."""

    def test_at_or_below_the_reference_window_the_config_value_is_verbatim(self):
        cfg = ResolverConfig()
        for w in (0, 1, 50, 4999, ohr._CEILING_REF_WINDOW):
            assert ohr._candidate_ceiling(cfg, w) == cfg.max_candidates_per_clip, (
                f'yeast-scale behaviour changed at W={w} — every run at or '
                f'below max_intron={ohr._CEILING_REF_WINDOW} must be identical')

    def test_above_the_reference_window_it_scales_linearly(self):
        cfg = ResolverConfig()
        base, ref = cfg.max_candidates_per_clip, ohr._CEILING_REF_WINDOW
        assert ohr._candidate_ceiling(cfg, 10 * ref) == 10 * base
        assert ohr._candidate_ceiling(cfg, 100 * ref) == 100 * base
        # the human case that motivated this: 500 kb window -> 200k candidates,
        # which covers 48 of the 50 measured chr5 clips
        assert ohr._candidate_ceiling(cfg, 500_000) == 200_000

    def test_it_is_a_floor_never_a_reduction(self):
        # a tightened config must not be loosened by the scaling, and a loosened
        # one must not be pulled down toward the default
        tight = ResolverConfig(max_candidates_per_clip=5)
        assert ohr._candidate_ceiling(tight, 1_000_000) >= 5
        loose = ResolverConfig(max_candidates_per_clip=10 ** 9)
        for w in (50, 5000, 500_000):
            assert ohr._candidate_ceiling(loose, w) >= 10 ** 9

    def test_the_guard_still_fires_on_a_genuine_outlier(self):
        # 221,086 candidates in a 500 kb window — the worst real human clip —
        # is still refused. The scaling must not disarm the circuit breaker.
        cfg = ResolverConfig()
        assert 221_086 > ohr._candidate_ceiling(cfg, 500_000)

    def test_yeast_scale_run_is_unchanged_end_to_end(self, tmp_path):
        # the fixture's window is max_intron=5000, so the scaled ceiling must
        # reproduce the flat-ceiling result exactly
        stats, _ = _run(tmp_path, ResolverConfig().max_candidates_per_clip,
                        'scaled_default')
        assert stats.refused_candidate_blowup == 0
        assert stats.candidates_evaluated > 100

    def test_tight_ceiling_still_bounds_work_at_a_large_max_intron(self, tmp_path):
        # scaling is relative to the configured value, so an operator who
        # deliberately tightens the guard keeps a tight guard
        ohr._BLOWUP_WARNED.clear()
        seq = _planted_genome(tmp_path)
        bam_in = _clipped_bam(tmp_path, seq)
        out = tmp_path / 'out_bigwin.bam'
        run_overhang_resolver(
            str(bam_in), str(tmp_path / 'g.fa'), str(out),
            config=ResolverConfig(max_candidates_per_clip=5, max_intron=500_000),
        )
        st = run_overhang_resolver.last_stats
        assert st.refused_candidate_blowup == 1
        assert st.candidates_evaluated <= 5 * (500_000 // ohr._CEILING_REF_WINDOW) + 1

    def test_slow_path_is_announced_once_not_silently_endured(self, caplog, monkeypatch):
        """The scaling trades a silent refusal for real work: on the human
        holdout it took candidates from 65k to 2.08M, which is ~50 min on the
        Python DP versus ~2 min with the numba kernel. A silent multi-hour run
        is the very failure mode this ceiling was added to prevent, so say so."""
        ohr._SLOW_WINDOW_WARNED.clear()
        caplog.set_level(logging.WARNING)
        # 2026-09-05: the kernel is on by default (lazy) — force the pure-Python path
        from rectify.core.splice import overhang_informativeness as _oi
        monkeypatch.setattr(_oi, '_HP_ED_NUMBA_LOADED', True)
        monkeypatch.setattr(_oi, '_hp_ed_bounded_numba', None)
        cfg = ResolverConfig()
        for _ in range(20):
            ohr._candidate_ceiling(cfg, 500_000)
        n = caplog.text.count('RECTIFY_HP_ED_NUMBA=1')
        assert n == 1, f'expected exactly one slow-path warning, got {n}'
        ohr._SLOW_WINDOW_WARNED.clear()

    def test_no_slow_path_warning_at_yeast_scale(self, caplog):
        ohr._SLOW_WINDOW_WARNED.clear()
        caplog.set_level(logging.WARNING)
        cfg = ResolverConfig()
        for w in (50, 5000, 9999):        # < 10x the floor => not the slow case
            ohr._candidate_ceiling(cfg, w)
        assert 'RECTIFY_HP_ED_NUMBA' not in caplog.text

    def test_no_slow_path_warning_when_the_kernel_is_active(self, caplog, monkeypatch):
        from rectify.core.splice import overhang_informativeness as oi
        ohr._SLOW_WINDOW_WARNED.clear()
        caplog.set_level(logging.WARNING)
        monkeypatch.setattr(oi, '_hp_ed_bounded_numba', object(), raising=False)
        ohr._candidate_ceiling(ResolverConfig(), 500_000)
        assert 'RECTIFY_HP_ED_NUMBA' not in caplog.text
        ohr._SLOW_WINDOW_WARNED.clear()

    def test_warning_names_the_contig_count_and_window_not_a_yeast_remedy(self, caplog):
        """ISSUE-010(b): the advice used to declare the contig repetitive and
        tell the operator to add it to RECTIFY_SKIP_REGIONS — on human that
        means discarding a whole chromosome."""
        ohr._BLOWUP_WARNED.clear()
        caplog.set_level(logging.WARNING)
        ohr._warn_blowup('chr5', 200_000, 221_086, 500_000)
        txt = caplog.text
        assert 'chr5' in txt and '221,086' in txt and '500,000' in txt
        assert '200,000' in txt
        assert '--max-intron' in txt and 'RECTIFY_HP_ED_NUMBA' in txt
        # skip-regions may still be MENTIONED, but only as the conditional
        # last resort for a genuinely repetitive contig
        i = txt.find('RECTIFY_SKIP_REGIONS')
        assert i > 0 and 'genuinely repetitive' in txt[:i]
        assert 'repetitive/low-complexity-contig class' not in txt
