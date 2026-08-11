#!/usr/bin/env python3
"""PI #2 MECHANISM probe — does a read's EXONIC error density predict its
JUNCTION-region error? (handoff §PI DESIGN NOTES #2; advisor-gated scope.)

The PI's read-reliability idea: a globally-hot read's junction region is also
error-riddled, so counting it as support for a novel (unannotated) junction
inflates discovery FDR; and a read's hotness is ESTIMABLE FROM ITS EXONIC PORTION
without junction truth → a self-calibrating per-read reliability covariate.

This probe tests the MECHANISM the signal depends on, using the injector's hot
reads. It is DELIBERATELY the qualitative half only:
  * MEASURED here (M1-safe, magnitude-INSENSITIVE): does exonic error density
    correlate with junction-window error across reads, and WHICH layer drives it?
  * NOT measured here (SIRV-gated, magnitude-SENSITIVE — do NOT compute on
    placeholder params): the actual novel-junction FDR LIFT from down-weighting.

Expected, and the design consequence:
  * layer-1 over-dispersion (GLOBAL hotness) -> positive correlation: exonic
    density predicts junction error -> the covariate has signal.
  * layer-2 burst ALONE (LOCAL, regions far apart) -> ~zero correlation: a read
    can be clean in the exon yet bursty at the junction.
  => exonic density predicts junction reliability only PROBABILISTICALLY ->
     SOFT down-weight / posterior input, NOT a hard filter (PI refinement (a)).
"""
from __future__ import annotations

import argparse
import math
import os
import random
import sys
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.dirname(__file__))
from sim.error_injector import (  # noqa: E402
    InjectorParams, _walk_events, replace, _indel_pmf_for_frac_ge2,
)


def _pearson(xs: List[float], ys: List[float]) -> float:
    n = len(xs)
    if n < 2:
        return 0.0
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx <= 0 or syy <= 0:
        return 0.0
    return sxy / math.sqrt(sxx * syy)


def autocorrelation(params: InjectorParams, rng: random.Random,
                    n_reads: int = 4000, read_len: int = 1000,
                    junc_start: int = 850, junc_len: int = 100) -> Dict[str, float]:
    """Exonic predictor = error density OUTSIDE the junction window; outcome =
    junction-window error count. The junction window is placed FAR from the bulk
    exon so a single short burst cannot span both (isolating global vs local)."""
    exon_dens: List[float] = []
    junc_cnt: List[float] = []
    junc_end = junc_start + junc_len
    exon_len = read_len - junc_len
    for _ in range(n_reads):
        _, ev = _walk_events(read_len, params, rng, clean_seq=None)
        jc = sum(1 for e in ev if junc_start <= e.pos < junc_end)
        ec = sum(1 for e in ev if not (junc_start <= e.pos < junc_end))
        exon_dens.append(ec / exon_len)
        junc_cnt.append(float(jc))
    return {
        "pearson_r": _pearson(exon_dens, junc_cnt),
        "mean_junc_err": sum(junc_cnt) / len(junc_cnt),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--reads", type=int, default=4000)
    ap.add_argument("--read-len", type=int, default=1000)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    base_rate = 0.03  # raised so the junction window has measurable variance
    pmf = _indel_pmf_for_frac_ge2(0.39)
    base = replace(InjectorParams.null(), base_rate=base_rate, indel_run_pmf=pmf)
    configs: List[Tuple[str, InjectorParams]] = [
        ("null (iid Poisson)", base),
        ("over-dispersion ONLY (layer 1)",
         replace(base, read_rate_dist="gamma", gamma_shape=0.3)),
        ("burst ONLY (layer 2)",
         replace(base, burst_on=True, p_cold_to_hot=0.05, p_hot_to_cold=0.4, hot_factor=8.0)),
        ("combined (realistic placeholder)",
         replace(base, read_rate_dist="gamma", gamma_shape=0.3, burst_on=True,
                 p_cold_to_hot=0.05, p_hot_to_cold=0.4, hot_factor=8.0)),
    ]

    print("PI #2 MECHANISM probe — exonic error density vs junction-window error")
    print(f"  reads={args.reads} read_len={args.read_len} base_rate={base_rate} "
          f"junction window=[850,950)\n")
    print(f"  {'config':38s} {'pearson_r':>10s} {'mean_junc_err':>14s}")
    print("  " + "-" * 64)
    results = {}
    for name, p in configs:
        m = autocorrelation(p, random.Random(rng.random()),
                            n_reads=args.reads, read_len=args.read_len)
        results[name] = m
        print(f"  {name:38s} {m['pearson_r']:>10.3f} {m['mean_junc_err']:>14.3f}")

    r_null = results["null (iid Poisson)"]["pearson_r"]
    r_disp = results["over-dispersion ONLY (layer 1)"]["pearson_r"]
    r_burst = results["burst ONLY (layer 2)"]["pearson_r"]
    print("\n  INTERPRETATION (mechanism, NOT an FDR-lift claim):")
    disp_drives = r_disp > r_null + 0.05 and r_disp > 0.1
    burst_local = abs(r_burst - r_null) < 0.05
    print(f"   - over-dispersion drives exonic->junction predictability: "
          f"{'YES' if disp_drives else 'NO'} (r {r_null:.3f} -> {r_disp:.3f})")
    print(f"   - burst alone stays LOCAL (no exonic->junction transfer): "
          f"{'YES' if burst_local else 'NO'} (r {r_burst:.3f} vs null {r_null:.3f})")
    print("   => exonic density predicts junction reliability only PROBABILISTICALLY")
    print("      -> SOFT down-weight / posterior input, NOT a hard filter (PI refinement (a)).")
    print("\n  NOTE: the FDR-LIFT quantification on stratum (G) is SIRV-gated "
          "(magnitude-sensitive); not run on placeholder params.")
    ok = disp_drives and burst_local
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
