"""Composition / interaction tests for the error-realism injector.

These are FAST (no calibration loop — that lives in error_injector.self_check).
They assert (1) NULL params reproduce the marginal/Poisson baseline, (2) each
layer moves its OWN statistic in the right direction, and (3) the BAM adapter's
pure grouping function builds the expected error track. They do NOT assert
realism — realism is SIRV-gated + distributional (see the module docstring).
"""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scripts.benchmark.sim.error_injector import (  # noqa: E402
    InjectorParams, ErrorEvent, inject, simulate_and_measure,
    measure_error_structure, events_from_alignment, _indel_pmf_for_frac_ge2,
    replace,
)

N, L = 1500, 500


def test_null_baseline_is_poisson_single_base():
    """NULL params => dispersion ~1, gap excess ~1, indels ~single-base."""
    p = replace(InjectorParams.null(), base_rate=0.02)
    m = simulate_and_measure(p, random.Random(11), n_reads=N, read_len=L)
    assert 0.8 <= m["dispersion_index"] <= 1.25, m
    assert 0.85 <= m["sub5_gap_excess"] <= 1.15, m
    assert m["indel_run_frac_ge2"] <= 0.02, m


def test_layer1_overdispersion():
    base = replace(InjectorParams.null(), base_rate=0.02)
    m0 = simulate_and_measure(base, random.Random(11), n_reads=N, read_len=L)
    l1 = replace(base, read_rate_dist="gamma", gamma_shape=0.25)
    m1 = simulate_and_measure(l1, random.Random(12), n_reads=N, read_len=L)
    assert m1["dispersion_index"] > 2.0 * m0["dispersion_index"], (m0, m1)


def test_layer2_burst_clustering():
    base = replace(InjectorParams.null(), base_rate=0.02)
    m0 = simulate_and_measure(base, random.Random(11), n_reads=N, read_len=L)
    l2 = replace(base, burst_on=True, p_cold_to_hot=0.05, p_hot_to_cold=0.4, hot_factor=8.0)
    m2 = simulate_and_measure(l2, random.Random(13), n_reads=N, read_len=L)
    assert m2["sub5_gap_excess"] > 1.4 * m0["sub5_gap_excess"], (m0, m2)


def test_layer3_indel_runs():
    base = replace(InjectorParams.null(), base_rate=0.02)
    l3 = replace(base, indel_run_pmf=_indel_pmf_for_frac_ge2(0.4))
    m3 = simulate_and_measure(l3, random.Random(14), n_reads=N, read_len=L)
    assert m3["indel_run_frac_ge2"] > 0.25, m3


def test_burst_does_not_zero_indels_or_subs():
    """Composition: turning on all layers still emits both indels and subs."""
    p = replace(InjectorParams.null(), base_rate=0.02, read_rate_dist="gamma",
                gamma_shape=0.3, burst_on=True, p_cold_to_hot=0.05,
                p_hot_to_cold=0.4, hot_factor=8.0,
                indel_run_pmf=_indel_pmf_for_frac_ge2(0.4))
    m = simulate_and_measure(p, random.Random(15), n_reads=N, read_len=L)
    assert m["marginal_event_rate"] > 0
    assert m["n_indel_events"] > 0


def test_inject_returns_consistent_dirty_seq():
    """Full inject path: dirty length = clean + inserted - deleted bases."""
    clean = "".join(random.Random(3).choice("ACGT") for _ in range(400))
    p = replace(InjectorParams.null(), base_rate=0.05,
                indel_run_pmf=_indel_pmf_for_frac_ge2(0.4))
    dirty, ev = inject(clean, p, random.Random(4))
    ins = sum(e.length for e in ev if e.kind == "ins")
    dele = sum(e.length for e in ev if e.kind == "del")
    assert len(dirty) == len(clean) + ins - dele, (len(dirty), len(clean), ins, dele)
    assert all(e.pos < len(clean) for e in ev)


def test_bam_adapter_pure_grouping():
    """events_from_alignment: I/D ops -> events of their length; mismatches ->
    len-1 sub events; ref span = M/=/X/D/N; clips/insertions don't count."""
    # 10M 2I 5M 3D 4M, mismatches at ABSOLUTE ref positions 102, 107
    # (get_aligned_pairs returns absolute reference coordinates)
    cig = [(4, 6), (0, 10), (1, 2), (0, 5), (2, 3), (0, 4)]
    events, span = events_from_alignment(cig, mismatch_ref_positions=[102, 107], ref_start=100)
    kinds = sorted((e.kind, e.length) for e in events)
    assert ("ins", 2) in kinds
    assert ("del", 3) in kinds
    assert kinds.count(("sub", 1)) == 2
    # ref span = 10 + 5 + 3 + 4 = 22 (soft-clip and insertion excluded)
    assert span == 22, span
    # events sorted by ref position, all >= ref_start
    positions = [e.pos for e in events]
    assert positions == sorted(positions)
    assert min(positions) >= 100


def test_measure_from_explicit_events():
    """measure_error_structure on a hand-built track: a clustered read vs a
    uniform read -> the clustered read raises the sub5 gap excess."""
    uniform = [ErrorEvent(pos=i * 50, kind="sub", length=1, bases="") for i in range(8)]
    clustered = [ErrorEvent(pos=p, kind="sub", length=1, bases="")
                 for p in (10, 11, 12, 13, 200, 201, 202, 203)]
    m = measure_error_structure([uniform, clustered], [400, 400])
    assert m["obs_gap_frac_lt"] > 0  # the clustered read contributes sub-5 gaps
    assert m["n_reads"] == 2


def test_exonic_coords_collapses_introns():
    """The N-in-span fix: on a SPLICED read, exonic_coords excludes the intron
    from the span and maps positions to exon-local coords, so an error gap does
    NOT straddle the intron and the rate is per-exonic-base.
    CIGAR 10M 1D 5M 100N 5M 1D 9M, a substitution at abs ref 120 (2nd exon)."""
    cig = [(0, 10), (2, 1), (0, 5), (3, 100), (0, 5), (2, 1), (0, 9)]
    # --- default (absolute): intron inflates span + the gap straddles it ---
    ev_abs, span_abs = events_from_alignment(cig, [120], ref_start=0)
    assert span_abs == 131, span_abs                 # 10+1+5+100+5+1+9
    pos_abs = sorted(e.pos for e in ev_abs)
    assert pos_abs == [10, 120, 121], pos_abs        # del@10, sub@120, del@121
    assert (120 - 10) == 110                          # gap straddles the 100bp intron
    # --- exonic_coords: intron collapsed out ---
    ev_ex, span_ex = events_from_alignment(cig, [120], ref_start=0, exonic_coords=True)
    assert span_ex == 31, span_ex                     # 10+1+5+5+1+9 (no intron)
    pos_ex = sorted(e.pos for e in ev_ex)
    assert pos_ex == [10, 20, 21], pos_ex             # sub/del shifted left by 100
    # the two adjacent-exon errors are now 1bp apart, not 110 -> no spurious gap
    assert (pos_ex[1] - pos_ex[0]) == 10 and (pos_ex[2] - pos_ex[1]) == 1
    # gap-free read: exonic_coords is a no-op (safe default for callers)
    cig2 = [(0, 10), (2, 1), (0, 9)]
    a0, s0 = events_from_alignment(cig2, [5], ref_start=0)
    a1, s1 = events_from_alignment(cig2, [5], ref_start=0, exonic_coords=True)
    assert s0 == s1 == 20
    assert sorted(e.pos for e in a0) == sorted(e.pos for e in a1)


if __name__ == "__main__":
    import inspect as _inspect
    _g = dict(globals())
    _tests = [(n, f) for n, f in _g.items()
              if n.startswith("test_") and callable(f)
              and not _inspect.signature(f).parameters]
    _tests.sort(key=lambda nf: nf[0])
    for _n, _f in _tests:
        _f()
        print(f"  PASS {_n}")
    print(f"all {len(_tests)} error_injector tests passed")
