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
    events_from_bam_read, calibrate_params,
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
        ev, span = events_from_bam_read(read)
        if span <= 0:
            continue
        per_read.append(ev)
        lengths.append(span)
        n += 1
        if args.max_reads and n >= args.max_reads:
            break
    if args.thin:
        per_read = _thin_events(per_read, lengths, args.thin, random.Random(args.seed))
    m = measure_error_structure(per_read, lengths)
    print(f"[measure-bam] {args.bam}  reads={int(m['n_reads'])} "
          f"thin={args.thin or 'off'}")
    for k in ("marginal_event_rate", "dispersion_index", "p90_over_median",
              "sub5_gap_excess", "indel_run_frac_ge2", "n_indel_events"):
        print(f"  {k:24s} {m[k]:.4f}")
    return 0


def cmd_inject_fastq(args) -> int:
    rng = random.Random(args.seed)
    if args.calibrate:
        params, _ = calibrate_params(random.Random(args.seed))
        print(f"[inject-fastq] calibrated PLACEHOLDER params: gamma_shape="
              f"{params.gamma_shape:.3f} p_cold_to_hot={params.p_cold_to_hot:.4f}",
              file=sys.stderr)
    else:
        params = InjectorParams.null()
    n = 0
    with _open(args.infq) as fin, open(args.outfq, "w") as fout:
        while True:
            h = fin.readline()
            if not h:
                break
            seq = fin.readline().strip()
            plus = fin.readline()
            fin.readline()  # original qual (discarded; structure-only validation)
            dirty, _ = inject(seq, params, rng)
            fout.write(h)
            fout.write(dirty + "\n")
            fout.write(plus if plus.startswith("+") else "+\n")
            fout.write("I" * len(dirty) + "\n")
            n += 1
    print(f"[inject-fastq] wrote {n} reads -> {args.outfq}", file=sys.stderr)
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
    mb.add_argument("--max-reads", type=int, default=0)
    mb.add_argument("--seed", type=int, default=7)
    mb.set_defaults(func=cmd_measure_bam)

    inj = sub.add_parser("inject-fastq")
    inj.add_argument("infq")
    inj.add_argument("outfq")
    inj.add_argument("--calibrate", action="store_true",
                     help="use coordinate-descent PLACEHOLDER params (else NULL)")
    inj.add_argument("--seed", type=int, default=7)
    inj.set_defaults(func=cmd_inject_fastq)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
