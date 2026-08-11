#!/usr/bin/env python3
"""SPIKE-IN scorer: build rich per-read truth (truth_schema) from the flat
gen_reads read_truth.tsv + fab_contig_truth.tsv, then score each arm BAM with
scorer.py:score_bam for the ambiguity-aware precision/FDR-by-canonicity metric.

For the FABRICATION panel: only the TRUE canonical junctions are registered as
truth (the drift target is deliberately absent). So a read that arm-B drifts to
the non-canonical drift site scores as a NON-CANONICAL FALSE POSITIVE
(by_canon["noncanonical"]["fp"]) = fabrication — exactly the metric the mission
calls make-or-break. Seed reads (spliced at the drift target) are EXCLUDED from
scoring (they exist only to seed the candidate pool).

Emits per-arm: junction recall, precision, FDR, TP/FP/FN, and the by-canonicity
split (canonical vs non-canonical FP). Stratified per contig (drift distance) if
--per-contig.

Usage:
  spikein_score.py --work-dir fab_sweep --arms A,B,E [--raw aligned.sorted]
"""
from __future__ import annotations
import argparse, os, sys, json
from collections import defaultdict
from typing import Dict, List

_here = os.path.dirname(os.path.abspath(__file__))
_repo = os.path.abspath(os.path.join(_here, "..", "..", ".."))
if _repo not in sys.path:
    sys.path.insert(0, _repo)

from rectify.core.benchmark.scorer import score_bam, load_genome
from rectify.core.benchmark.truth_schema import (
    ReadTruth, JunctionTruth, JunctionClass, SplitTag,
)


def load_fab_truth(work_dir: str) -> Dict[str, dict]:
    """chrom -> {donor, acceptor(true, canonical), drift_acceptor, drift_dist, e5, e3}."""
    p = os.path.join(work_dir, "fab_contig_truth.tsv")
    out = {}
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[idx["chrom"]]] = {
                "donor": int(f[idx["donor"]]),
                "acceptor": int(f[idx["acceptor"]]),
                "drift_acceptor": int(f[idx["drift_acceptor"]]),
                "drift_dist": int(f[idx["drift_dist"]]),
                "e5": int(f[idx["e5"]]), "e3": int(f[idx["e3"]]),
            }
    return out


def build_rich_truth(work_dir: str, genome: Dict[str, str]) -> List[ReadTruth]:
    """One ReadTruth per NON-seed read, carrying ONLY the true canonical junction.
    The drift target is intentionally NOT registered -> a drifted call is a FP."""
    fab = load_fab_truth(work_dir)
    rt_path = os.path.join(work_dir, "read_truth.tsv")
    hdr = open(rt_path).readline().rstrip("\n").split("\t")
    idx = {h: i for i, h in enumerate(hdr)}
    rows: List[ReadTruth] = []
    for line in open(rt_path):
        if line.startswith("read_id\t"):
            continue
        f = line.rstrip("\n").split("\t")
        tid = f[idx["tid"]]
        if tid.endswith("_seed"):
            continue  # seed reads: pool primer only, never scored
        chrom = f[idx["chrom"]]
        info = fab[chrom]
        seq = genome[chrom]
        D, A = info["donor"], info["acceptor"]
        jt = JunctionTruth.from_intron(D, A, "+", JunctionClass.NNC, seq)
        # sanity: the registered truth must be canonical
        rows.append(ReadTruth(
            read_id=f[idx["read_id"]], true_locus=chrom, true_transcript=tid,
            chrom=chrom, strand="+", genome_start=D - info["e5"],
            genome_end=A + info["e3"], junctions=[jt],
            stratum=f"FABD{info['drift_dist']}", split=SplitTag.TEST,
        ))
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--arms", default="A,B,E")
    ap.add_argument("--raw", default="aligned.sorted",
                    help="raw-mm2 BAM stem (baseline, no refinement)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    wd = args.work_dir

    genome = load_genome(os.path.join(wd, "sim_ref.fa"))
    truth = build_rich_truth(wd, genome)
    truth_map = {t.read_id: t for t in truth}
    fab = load_fab_truth(wd)

    arm_stems = [args.raw] + [f"arm_{a}" for a in args.arms.split(",")]
    arm_labels = {args.raw: "raw-mm2"}
    for a in args.arms.split(","):
        arm_labels[f"arm_{a}"] = f"arm-{a}"

    results = {}
    for stem in arm_stems:
        bam = os.path.join(wd, f"{stem}.bam")
        if not os.path.exists(bam):
            continue
        s = score_bam(bam, truth_map, genome, aligner_name=arm_labels[stem])
        summ = s.summary()
        results[arm_labels[stem]] = {
            "recall": summ["junction_recall"],
            "precision": summ["junction_precision"],
            "fdr": summ["junction_fdr"],
            "tp": summ["junction_tp"], "fp": summ["junction_fp"],
            "fn": summ["junction_fn"],
            "fp_noncanonical": s.junction.by_canon.get("noncanonical", {}).get("fp", 0),
            "fp_canonical": s.junction.by_canon.get("canonical", {}).get("fp", 0),
            "reads_scored": summ["reads_scored"],
        }

    # print table
    print(f"# SPIKE-IN scored: {wd}  ({len(truth)} true-junction reads, "
          f"{len(fab)} contigs, drift dists "
          f"{sorted(set(v['drift_dist'] for v in fab.values()))})")
    print(f"{'arm':10s} {'recall':>7s} {'prec':>7s} {'FDR':>7s} {'TP':>5s} "
          f"{'FP':>5s} {'FP_noncanon':>11s} {'FP_canon':>8s} {'FN':>5s}")
    print("-" * 78)
    for lab, r in results.items():
        print(f"{lab:10s} {r['recall']:7.4f} {r['precision']:7.4f} {r['fdr']:7.4f} "
              f"{r['tp']:5d} {r['fp']:5d} {r['fp_noncanonical']:11d} "
              f"{r['fp_canonical']:8d} {r['fn']:5d}")

    out = args.out or os.path.join(wd, "spikein_score.json")
    with open(out, "w") as fh:
        json.dump({"work_dir": wd, "n_true_reads": len(truth),
                   "drift_dists": sorted(set(v["drift_dist"] for v in fab.values())),
                   "arms": results}, fh, indent=2)
    print(f"\n[spikein_score] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
