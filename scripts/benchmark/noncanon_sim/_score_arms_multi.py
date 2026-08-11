#!/usr/bin/env python3
"""Compact multi-arm recovery + over-shift table on one panel.

Recovery is ambiguity-aware (chimeric_consensus.normalize_junction), identical to
paired_arm_test.py.  For each requested arm label <L> (reads arm_<L>.bam) we report,
per context cell, the fraction of reads whose refined junction normalizes to the TRUE
junction (recovery), and — vs the arm-B baseline — how many reads MOVED (called junction
!= arm-B's) and how many of those moves LOST recovery that arm-B had (over-shift).

Usage: _score_arms_multi.py <work_dir> B E_m3p0 Ff Ffg [--cells ACC_A_D0,ACC_A_D2]
       (first label is the MOVE baseline; conventionally B)
"""
import os, sys, argparse
from collections import defaultdict
import pysam

_here = os.path.dirname(os.path.abspath(__file__))
_repo = os.path.abspath(os.path.join(_here, "..", "..", ".."))
if _repo not in sys.path:
    sys.path.insert(0, _repo)
from rectify.core.consensus.chimeric_consensus import normalize_junction


def load_genome(fa):
    g, name = {}, None
    for l in open(fa):
        if l.startswith(">"):
            name = l[1:].strip().split()[0]; g[name] = []
        else:
            g[name].append(l.strip())
    return {k: "".join(v) for k, v in g.items()}


def read_junctions(bam):
    out = {}
    b = pysam.AlignmentFile(bam, "rb")
    for r in b:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            out.setdefault(r.query_name, None); continue
        ref = r.reference_start; jn = None
        for op, ln in (r.cigartuples or []):
            if op == 3:
                jn = (ref, ref + ln); break
            if op in (0, 2, 3, 7, 8):
                ref += ln
        out[r.query_name] = jn
    b.close()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work_dir")
    ap.add_argument("labels", nargs="+")
    ap.add_argument("--cells", default="")
    ap.add_argument("--rungs", default="R3,R0flank",
                    help="motif_rungs to score (default cryptic + flank-HP)")
    args = ap.parse_args()
    wd = args.work_dir
    want_cells = set(c for c in args.cells.split(",") if c)
    want_rungs = set(args.rungs.split(","))

    genome = load_genome(os.path.join(wd, "sim_ref.fa"))
    rt = os.path.join(wd, "read_truth.tsv")
    hdr = open(rt).readline().rstrip("\n").split("\t")
    idx = {h: i for i, h in enumerate(hdr)}
    truth = {}
    for line in open(rt):
        if line.startswith("read_id\t"):
            continue
        f = line.rstrip("\n").split("\t")
        truth[f[idx["read_id"]]] = f

    jn = {L: read_junctions(os.path.join(wd, f"arm_{L}.bam")) for L in args.labels}
    base = args.labels[0]

    def rec(rid, L):
        f = truth[rid]; g = genome[f[idx["chrom"]]]
        j = jn[L].get(rid)
        if j is None:
            return None
        tn = normalize_junction(int(f[idx["true_donor"]]), int(f[idx["true_acceptor"]]), g)
        return 1 if normalize_junction(j[0], j[1], g) == tn else 0

    cells = defaultdict(list)
    for rid, f in truth.items():
        if f[idx["has_true_junction"]] != "1":
            continue
        if f[idx["motif_rung"]] not in want_rungs:
            continue
        cells[f[idx["context"]]].append(rid)

    order = sorted(cells)
    print(f"# panel={wd}  baseline(move-ref)=arm_{base}  rungs={sorted(want_rungs)}")
    head = f"{'cell':12s} {'n':>4s}"
    for L in args.labels:
        head += f" | {('rec_'+L):>10s}"
    head += f" || {'moved':>6s} {'overshift':>9s}   (moved/overshift are vs arm_%s, last arm)" % base
    print(head)
    print("-" * len(head))
    test = args.labels[-1]
    for c in order:
        if want_cells and c not in want_cells:
            continue
        rids = cells[c]; n = len(rids)
        line = f"{c:12s} {n:>4d}"
        recs = {}
        for L in args.labels:
            rr = [rec(r, L) for r in rids]
            recs[L] = rr
            got = sum(1 for x in rr if x == 1)
            line += f" | {got/n:10.3f}"
        moved = overshift = 0
        for k, r in enumerate(rids):
            jb = jn[base].get(r); jt = jn[test].get(r)
            if jb != jt:
                moved += 1
                if recs[base][k] == 1 and recs[test][k] != 1:
                    overshift += 1
        line += f" || {moved:>6d} {overshift:>9d}"
        print(line)
    # overall for the scored rungs
    allr = [r for c in order for r in cells[c]]
    print("-" * len(head))
    tot = f"{'OVERALL':12s} {len(allr):>4d}"
    for L in args.labels:
        got = sum(1 for r in allr if rec(r, L) == 1)
        tot += f" | {got/len(allr):10.3f}"
    mv = ov = 0
    for r in allr:
        if jn[base].get(r) != jn[test].get(r):
            mv += 1
            if rec(r, base) == 1 and rec(r, test) != 1:
                ov += 1
    tot += f" || {mv:>6d} {ov:>9d}"
    print(tot)


if __name__ == "__main__":
    main()
