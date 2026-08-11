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
from rectify.core.benchmark.scorer import (  # noqa: E402
    load_genome, extract_junctions, junction_ambiguity_window,
    _canonical_within_window,
)
from rectify.core.benchmark.truth_schema import read_truth_table  # noqa: E402
from rectify.core.consensus.chimeric_consensus import normalize_junction  # noqa: E402
import pysam  # noqa: E402


def rung_of_read(rid: str) -> str:
    # read_id = nj_<label>_r### (synthetic) or njr_<label>_r### (real-genome);
    # label may contain '_' (e.g. GA-AG_1off). Strip the prefix + the _r### suffix.
    body = rid.split("_", 1)[1] if rid.startswith(("nj_", "njr_")) else rid
    return body.rsplit("_r", 1)[0]


def scan_bam(bam_path, truth_by_id, genome):
    """Per read_id: {rid: (recovered, n_fp, n_fp_noncanon)} for this aligner.

    recovered   = true junction(s) called ambiguity-aware AND no FP junction
                  (mirrors scorer._score_read's junction rule) — the RECALL axis.
    n_fp        = false junctions called (not matching any truth) — the FDR axis:
                  the count that exposes an indiscriminate emitter (mapPacBio's
                  real-data pathology — it "recovers" by calling junctions everywhere).
    n_fp_noncanon = of those FP junctions, how many are NON-canonical (the spurious-
                  novel rate; on an intron-free control read every N-op is fabricated).
    """
    out = {}
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
            matched, fp, fp_nc = set(), 0, 0
            for (cs, ce) in called:
                ns, ne = normalize_junction(cs, ce, gseq) if gseq else (cs, ce)
                if (ns, ne) in truth_set:
                    matched.add((ns, ne))
                else:
                    fp += 1
                    if gseq:
                        la, ra = junction_ambiguity_window(ns, ne, gseq)
                        if not _canonical_within_window(ns, ne, gseq, la, ra):
                            fp_nc += 1
            recovered = bool(truth_set) and matched == truth_set and fp == 0
            out[read.query_name] = (recovered, fp, fp_nc)
    return out


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
    scan = {}
    for spec in args.bam:
        name, path = spec.split("=", 1)
        if not os.path.exists(path):
            print(f"[panel] SKIP {name}: {path} missing", file=sys.stderr)
            continue
        aligners.append(name)
        scan[name] = scan_bam(path, truth_by_id, genome)
        nrec = sum(1 for v in scan[name].values() if v[0])
        print(f"[panel] {name}: recovered {nrec} reads", file=sys.stderr)

    def rec_of(a, i):
        v = scan[a].get(i)
        return bool(v and v[0])

    def fp_of(a, i):
        v = scan[a].get(i)
        return v[1] if v else 0

    def fpnc_of(a, i):
        v = scan[a].get(i)
        return v[2] if v else 0

    # an intron-free control rung carries NO truth junction => recall is N/A there
    def is_intronfree(rid):
        rt = truth_by_id.get(rid)
        return rt is not None and not rt.junctions

    ORDER = ["GT-AG_canon", "GC-AG", "AT-AC", "GA-AG_1off", "CT-AC_2off",
             "CA-TC_deep", "INTRONFREE"]
    rungs = [r for r in ORDER if r in ids_by_rung] + \
            [r for r in ids_by_rung if r not in ORDER]

    lines = []
    lines.append("===== PANEL RECALL (recover true junction) + FDR (false / spurious-noncanonical junctions) =====")
    lines.append("Recall alone is GAMED by indiscriminate emitters (mapPacBio 'recovers' by calling junctions")
    lines.append("everywhere) — read it WITH the FP columns. INTRONFREE control: truth has NO junction, so any")
    lines.append("called junction is FABRICATED (mapPacBio's real-data 97.7%-spurious pathology, reproduced).")
    hdr = f"{'rung':14s} {'n':>4s} "
    for a in aligners:
        hdr += f"{a[:7]+'.rec':>11s} {a[:5]+'.fpnc':>11s} "
    hdr += f"{'PANELrec':>8s}"
    lines.append(hdr)
    for r in rungs:
        ids = ids_by_rung[r]
        n = max(1, len(ids))
        intronfree = all(is_intronfree(i) for i in ids) if ids else False
        row = f"{r:14s} {len(ids):4d} "
        for a in aligners:
            rec = "  n/a " if intronfree else f"{sum(rec_of(a, i) for i in ids)/n:6.3f}"
            fpnc = sum(fpnc_of(a, i) for i in ids) / n   # spurious-noncanonical per read
            row += f"{rec:>11s} {fpnc:11.3f} "
        # panel recall = union of recoverers (only meaningful on junction rungs)
        if intronfree:
            row += f"{'  n/a':>8s}"
        else:
            u = sum(1 for i in ids if any(rec_of(a, i) for a in aligners)) / n
            row += f"{u:8.3f}"
        lines.append(row)
    lines.append("")
    lines.append(".rec  = recall of the TRUE novel junction (ambiguity-aware). .fpnc = spurious NON-canonical")
    lines.append("        junctions PER READ (the FDR axis; on INTRONFREE every one is fabricated from noise).")
    lines.append("PANELrec = union recall (>=1 aligner). READ: a native member is justified when the precise")
    lines.append("  members (mm2/uLTRA/deSALT) FLATTEN non-canonical (low .rec) AND the only recoverer (mapPacBio)")
    lines.append("  does it by FABRICATING (high .fpnc, esp. on INTRONFREE + under error) => no PRECISE novel")
    lines.append("  discoverer exists => the calibrated-empirical member's target. High recall + low .fpnc for a")
    lines.append("  panel member would instead mean the panel already covers it precisely (pivot).")
    out = "\n".join(lines)
    print(out)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(out + "\n")


if __name__ == "__main__":
    main()
