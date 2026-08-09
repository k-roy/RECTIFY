#!/usr/bin/env python3
"""Sherlock-side error-realism VALIDATION (turnkey; run when Sherlock auth is up).

Completes the FIRST-TASK validation loop (SPEC §ERROR-REALISM): does
pbsim+injector match REAL ONT on the three correlation statistics? Reuses the
SAME measurement code as the injector (``sim.error_injector``) so real / pbsim /
pbsim+injector are scored identically (advisor note #5 — source-agnostic).

Two subcommands:

  measure-bam <bam>     Measure the 3 stats on an aligned BAM via the shared
                        adapter. Needs the MD tag (minimap2 --MD or `calmd`).
                        Optionally Bernoulli-THIN events to a target marginal
                        rate first (the err_corr.py recipe: raw ERRHMM-ONT is
                        ~7x too hot vs modern dorado, so only SHAPE is comparable
                        — thin to the real rate before comparing structure).
                        ** NOT VERIFIED on a real BAM (Sherlock down); the pure
                        grouping `events_from_alignment` IS unit-tested. **

  inject-fastq <in> <out>   Apply the calibrated injector to every read of a
                        FASTQ (= the "+injector" layer on top of pbsim reads),
                        emit a corrupted FASTQ to RE-ALIGN, then `measure-bam`.

VALIDATION PROTOCOL (the RESUME, run on Sherlock):
  1. REAL:   minimap2 -ax splice --MD <real_DRS.fq> | sort > real.bam
             python error_realism_validate.py measure-bam real.bam --thin 0.019
  2. pbsim (before): measure-bam on the pbsim3 round-trip BAM (--thin 0.019)
  3. pbsim+injector (after): inject-fastq pbsim.fq -> pbsim_inj.fq ; re-align ;
             measure-bam pbsim_inj.bam --thin 0.019
  4. CALIBRATE magnitude vs SIRV ABSOLUTE truth: run step (1) on SIRV reads
             aligned to the SIRV reference (read-vs-known-sequence removes the
             alignment-artifact confound) and re-fit calibrate_params to THOSE
             targets, not the read-vs-ref upper bound.
  Compare the 3 stats across {real, pbsim, pbsim+injector}: success = +injector
  closes the pbsim->real gap on dispersion / sub5-gap-excess / indel-run.
"""
from __future__ import annotations

import argparse
import os
import random
import sys
from typing import List

sys.path.insert(0, os.path.dirname(__file__))
from sim.error_injector import (  # noqa: E402
    InjectorParams, ErrorEvent, inject, measure_error_structure,
    events_from_bam_read, placeholder_params,
)


def _thin_events(per_read: List[List[ErrorEvent]], lengths: List[int],
                 target_rate: float, rng: random.Random) -> List[List[ErrorEvent]]:
    """Bernoulli-thin events so the marginal per-position event rate ~= target
    (rate-match before comparing SHAPE; the err_corr.py recipe)."""
    total = sum(len(ev) for ev in per_read)
    tot_len = sum(lengths)
    cur = (total / tot_len) if tot_len else 0.0
    if cur <= target_rate or cur == 0:
        return per_read
    keep = target_rate / cur
    return [[e for e in ev if rng.random() < keep] for ev in per_read]


def cmd_measure_bam(args) -> int:
    import pysam
    bam = pysam.AlignmentFile(args.bam, "rb")
    per_read: List[List[ErrorEvent]] = []
    lengths: List[int] = []
    n = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        ev, span = events_from_bam_read(read, exonic_coords=True)
        if span <= 0:
            continue
        if args.len_min and span < args.len_min:
            continue
        if args.len_max and span > args.len_max:
            continue
        per_read.append(ev)
        lengths.append(span)
        n += 1
        if args.max_reads and n >= args.max_reads:
            break
    if args.thin:
        per_read = _thin_events(per_read, lengths, args.thin, random.Random(args.seed))
    m = measure_error_structure(per_read, lengths)
    win = f"[{args.len_min or 0},{args.len_max or 'inf'}]"
    print(f"[measure-bam] {args.bam}  reads={int(m['n_reads'])} "
          f"len_window={win} thin={args.thin or 'off'}")
    for k in ("mean_ref_len", "marginal_event_rate", "mean_events_per_read",
              "dispersion_index", "overdisp_v", "p90_over_median",
              "sub5_gap_excess", "indel_run_frac_ge2", "n_indel_events"):
        print(f"  {k:24s} {m[k]:.4f}")
    return 0


def _pearson(xs, ys) -> float:
    n = len(xs)
    if n < 2:
        return 0.0
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx <= 0 or syy <= 0:
        return 0.0
    return sxy / (sxx ** 0.5 * syy ** 0.5)


def cmd_autocorr_bam(args) -> int:
    """PI #2(i) on REAL data: does error density in one part of a read predict
    error density in a well-SEPARATED part? Split each read's aligned ref span
    into a head window and a tail window with a gap between them, correlate the
    two densities across reads (length-windowed). A positive r => a GLOBAL hotness
    component (the layer-1 over-dispersion the reliability covariate needs); r~0
    => purely local errors (covariate weak). This real r is the third constraint
    that identifies the global-vs-local split (overdisp_v + gap5x cannot)."""
    import pysam
    bam = pysam.AlignmentFile(args.bam, "rb")
    head_d, tail_d = [], []
    n = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        ev, span = events_from_bam_read(read, exonic_coords=True)
        if span <= 0:
            continue
        if args.len_min and span < args.len_min:
            continue
        if args.len_max and span > args.len_max:
            continue
        # exon-local positions run [reference_start, reference_start+span) with
        # introns collapsed out, so head/tail windows are along the READ axis.
        rs = read.reference_start
        w = int(span * args.frac)               # window length each side
        if w < args.min_win:
            continue
        h0, h1 = rs, rs + w                      # head window
        t0, t1 = rs + span - w, rs + span        # tail window (gap = span - 2w between)
        hc = sum(1 for e in ev if h0 <= e.pos < h1)
        tc = sum(1 for e in ev if t0 <= e.pos < t1)
        head_d.append(hc / w)
        tail_d.append(tc / w)
        n += 1
        if args.max_reads and n >= args.max_reads:
            break
    r = _pearson(head_d, tail_d)
    win = f"[{args.len_min or 0},{args.len_max or 'inf'}]"
    print(f"[autocorr-bam] {args.bam}  reads={n} len_window={win} "
          f"frac={args.frac} (head vs tail {int(args.frac*100)}% windows)")
    print(f"  head-vs-tail error-density pearson_r = {r:.4f}")
    print(f"  INTERPRETATION: r>0 => GLOBAL hotness present (layer-1 over-dispersion;")
    print(f"  the PI-#2 reliability covariate has signal); the magnitude of r is the")
    print(f"  global fraction that pins the injector's global-vs-local split.")
    return 0


def cmd_inject_fastq(args) -> int:
    rng = random.Random(args.seed)
    if args.null:
        params = InjectorParams.null()
        label = "NULL (marginal Poisson)"
    else:
        params = placeholder_params(base_rate=args.base_rate)
        label = (f"PLACEHOLDER (gamma_shape={params.gamma_shape} "
                 f"hot_factor={params.hot_factor} p_cold_to_hot={params.p_cold_to_hot})")
    print(f"[inject-fastq] params: {label}", file=sys.stderr)
    n = 0
    per_read, lengths = [], []     # truth-track accumulation (guard #3)
    with _open(args.infq) as fin, open(args.outfq, "w") as fout:
        while True:
            h = fin.readline()
            if not h:
                break
            seq = fin.readline().strip()
            plus = fin.readline()
            fin.readline()  # original qual (discarded; structure-only validation)
            dirty, ev = inject(seq, params, rng)
            per_read.append(ev)
            lengths.append(len(seq))
            fout.write(h)
            fout.write(dirty + "\n")
            fout.write(plus if plus.startswith("+") else "+\n")
            fout.write("I" * len(dirty) + "\n")
            n += 1
    print(f"[inject-fastq] wrote {n} reads -> {args.outfq}", file=sys.stderr)
    # report the injected-TRUTH structure (no alignment) — compare to measure-bam
    # on the re-aligned output to isolate the ALIGNMENT-ONLY inflation (guard #3).
    m = measure_error_structure(per_read, lengths)
    print("[inject-fastq] injected-TRUTH structure (pre-alignment):")
    for k in ("marginal_event_rate", "overdisp_v", "dispersion_index",
              "sub5_gap_excess", "indel_run_frac_ge2"):
        print(f"  {k:24s} {m[k]:.4f}")
    return 0


def _open(path):
    import gzip
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    mb = sub.add_parser("measure-bam")
    mb.add_argument("bam")
    mb.add_argument("--thin", type=float, default=None,
                    help="Bernoulli-thin events to this marginal rate (e.g. 0.019)")
    mb.add_argument("--len-min", type=int, default=0, dest="len_min",
                    help="length-control: only reads with ref span >= this")
    mb.add_argument("--len-max", type=int, default=0, dest="len_max",
                    help="length-control: only reads with ref span <= this")
    mb.add_argument("--max-reads", type=int, default=0)
    mb.add_argument("--seed", type=int, default=7)
    mb.set_defaults(func=cmd_measure_bam)

    ac = sub.add_parser("autocorr-bam")
    ac.add_argument("bam")
    ac.add_argument("--frac", type=float, default=0.3,
                    help="head/tail window = this fraction of the aligned span each")
    ac.add_argument("--min-win", type=int, default=80, dest="min_win")
    ac.add_argument("--len-min", type=int, default=0, dest="len_min")
    ac.add_argument("--len-max", type=int, default=0, dest="len_max")
    ac.add_argument("--max-reads", type=int, default=0)
    ac.set_defaults(func=cmd_autocorr_bam)

    inj = sub.add_parser("inject-fastq")
    inj.add_argument("infq")
    inj.add_argument("outfq")
    inj.add_argument("--null", action="store_true",
                     help="use NULL marginal params (default: hand-picked PLACEHOLDER)")
    inj.add_argument("--base-rate", type=float, default=0.0153, dest="base_rate")
    inj.add_argument("--seed", type=int, default=7)
    inj.set_defaults(func=cmd_inject_fastq)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
