#!/usr/bin/env python3
"""Cluster driver: run the RECTIFY aligner PANEL on a novel-junction ladder corpus.

Runs each requested aligner via the shipped `run_multi_aligner` wrapper (one at a
time, each in its own try — so one aligner failing on the tiny synthetic contigs does
NOT lose the others), then hands the per-aligner BAMs to panel_blindspot_score for the
per-rung panel-native recovery table. Meant to run inside run_panel_blindspot.sbatch
on Sherlock (rectify env). Writes the report to <corpus>/panel_blindspot.report.txt.

Usage:
  python scripts/benchmark/run_panel_on_corpus.py --corpus <dir> \
      --aligners minimap2 gapmm2 deSALT uLTRA --threads 4
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--aligners", nargs="+",
                    default=["minimap2", "mapPacBio", "deSALT", "uLTRA"])
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    from rectify.core.align.multi_aligner import run_multi_aligner  # noqa: E402

    ref = os.path.join(args.corpus, "ref.fa")
    reads = os.path.join(args.corpus, "reads.fastq")
    outdir = os.path.join(args.corpus, "panel_bams")
    os.makedirs(outdir, exist_ok=True)
    bam_args = []
    for a in args.aligners:
        try:
            print(f"[panel-run] {a} ...", file=sys.stderr, flush=True)
            out = run_multi_aligner(reads, ref, outdir, "nj",
                                    threads=args.threads, aligners=[a])
            bam = out.get(a)
            if bam and os.path.exists(bam):
                bam_args += ["--bam", f"{a}={bam}"]
                print(f"[panel-run] {a} OK -> {bam}", file=sys.stderr, flush=True)
            else:
                print(f"[panel-run] {a} produced no BAM", file=sys.stderr, flush=True)
        except Exception as exc:  # one aligner failing must not sink the rest
            print(f"[panel-run] {a} FAILED: {exc}", file=sys.stderr, flush=True)

    if not bam_args:
        print("[panel-run] NO aligner produced a BAM — abort", file=sys.stderr)
        sys.exit(3)

    report = os.path.join(args.corpus, "panel_blindspot.report.txt")
    scorer = os.path.join(os.path.dirname(__file__), "panel_blindspot_score.py")
    cmd = [sys.executable, scorer, "--corpus", args.corpus, "--out", report] + bam_args
    print("[panel-run] scoring: " + " ".join(cmd), file=sys.stderr, flush=True)
    subprocess.run(cmd, check=True)
    print(f"[panel-run] DONE -> {report}", file=sys.stderr)


if __name__ == "__main__":
    main()
