#!/usr/bin/env python3
"""Panel-native novel-junction recovery per rung — the DECIDING blind-spot number.

Consumes the novel-junction ladder corpus (novel_junction_blindspot.py --emit-corpus)
and the per-aligner BAMs from a cluster 5-aligner PANEL run, and reports, per
deviation rung: EACH aligner's native recovery of the true novel junction AND the
PANEL-UNION recovery (fraction of reads where >=1 aligner calls the true site,
ambiguity-aware). Panel-union recovery is the honest isoform-flattening number the
program needs — minimap2-alone flattening 47-90% (dev/NOVEL_JUNCTION_BLINDSPOT.md) is
only a lower bound; if the WHOLE panel co-fails, the native member is justified.

A read is "recovered" by an aligner iff that aligner's alignment calls the read's
true junction (ambiguity-aware via the shipped normalize_junction) with NO false
junction — same rule as scorer._score_read, applied per aligner per read.

Usage:
  python scripts/benchmark/panel_blindspot_score.py --corpus <dir> \
      --bam minimap2=<...> --bam deSALT=<...> --bam gmap=<...> [--bam uLTRA=<...>] [--bam gapmm2=<...>]
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.scorer import load_genome, extract_junctions  # noqa: E402
from rectify.core.benchmark.truth_schema import read_truth_table  # noqa: E402
from rectify.core.consensus.chimeric_consensus import normalize_junction  # noqa: E402
import pysam  # noqa: E402


def rung_of_read(rid: str) -> str:
    # read_id = nj_<label>_r### (synthetic) or njr_<label>_r### (real-genome);
    # label may contain '_' (e.g. GA-AG_1off). Strip the prefix + the _r### suffix.
    body = rid.split("_", 1)[1] if rid.startswith(("nj_", "njr_")) else rid
    return body.rsplit("_r", 1)[0]


def recovered_ids_from_bam(bam_path, truth_by_id, genome):
    """Return the set of read_ids this aligner RECOVERS (true junction called,
    ambiguity-aware, no FP junction) — mirrors scorer._score_read's junction rule."""
    rec = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            rt = truth_by_id.get(read.query_name)
            if rt is None:
                continue
            gseq = genome.get(read.reference_name, "")
            truth_set = {(j.intron_start, j.intron_end) for j in rt.junctions}
            called = extract_junctions(read.reference_start, read.cigartuples)
            matched, fp = set(), 0
            for (cs, ce) in called:
                ns, ne = normalize_junction(cs, ce, gseq) if gseq else (cs, ce)
                if (ns, ne) in truth_set:
                    matched.add((ns, ne))
                else:
                    fp += 1
            if truth_set and matched == truth_set and fp == 0:
                rec.add(read.query_name)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", required=True, help="dir with ref.fa + truth.tsv")
    ap.add_argument("--bam", action="append", default=[],
                    help="aligner=path (repeatable), e.g. minimap2=/path/mm2.bam")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    genome = load_genome(os.path.join(args.corpus, "ref.fa"))
    truth = read_truth_table(os.path.join(args.corpus, "truth.tsv"))
    truth_by_id = {t.read_id: t for t in truth}
    ids_by_rung = defaultdict(list)
    for t in truth:
        ids_by_rung[rung_of_read(t.read_id)].append(t.read_id)

    aligners = []
    recovered = {}
    for spec in args.bam:
        name, path = spec.split("=", 1)
        if not os.path.exists(path):
            print(f"[panel] SKIP {name}: {path} missing", file=sys.stderr)
            continue
        aligners.append(name)
        recovered[name] = recovered_ids_from_bam(path, truth_by_id, genome)
        print(f"[panel] {name}: recovered {len(recovered[name])} reads", file=sys.stderr)

    # deviation order (must match the ladder)
    ORDER = ["GT-AG_canon", "GC-AG", "AT-AC", "GA-AG_1off", "CT-AC_2off", "CA-TC_deep"]
    rungs = [r for r in ORDER if r in ids_by_rung] + \
            [r for r in ids_by_rung if r not in ORDER]

    lines = []
    hdr = f"{'rung':14s} {'n':>4s} " + " ".join(f"{a[:8]:>8s}" for a in aligners) + f" {'PANEL':>7s} {'BLINDSP':>7s}"
    lines.append("========== PANEL-NATIVE NOVEL-JUNCTION RECOVERY (per rung) ==========")
    lines.append(hdr)
    for r in rungs:
        ids = ids_by_rung[r]
        n = len(ids)
        per = []
        for a in aligners:
            rec = sum(1 for i in ids if i in recovered[a])
            per.append(rec / max(1, n))
        union = sum(1 for i in ids if any(i in recovered[a] for a in aligners))
        u = union / max(1, n)
        lines.append(f"{r:14s} {n:4d} " + " ".join(f"{p:8.3f}" for p in per) +
                     f" {u:7.3f} {1 - u:7.3f}")
    lines.append("")
    lines.append("PANEL = fraction of reads where >=1 aligner recovers the TRUE novel junction (ambiguity-aware)")
    lines.append("BLINDSP = 1 - PANEL = the panel-native isoform-FLATTENING rate (the deciding number)")
    lines.append("READ: canonical rung PANEL ~1.0 = sound. If BLINDSP stays HIGH on non-canonical rungs across")
    lines.append("  the WHOLE panel (not just minimap2) => the panel co-fails / herds => native member JUSTIFIED.")
    lines.append("  If some aligner (deSALT/gmap) recovers what minimap2 flattens => the panel already covers it")
    lines.append("  => the gain is arbitration/consensus, not a new placer.")
    out = "\n".join(lines)
    print(out)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(out + "\n")


if __name__ == "__main__":
    main()
