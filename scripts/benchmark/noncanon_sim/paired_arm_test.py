#!/usr/bin/env python3
"""PAIRED arm-B vs arm-C significance test on the cryptic non-canonical reads,
stratified by HP run length (the arm-C dose-response power test).

Because arm-B and arm-C refine the SAME reads, recovery is a PAIRED binary
outcome per read. We test whether the -logP law (arm-C) recovers the true
non-canonical acceptor MORE often than motif-blind (arm-B), with:
  * per-HP-length recovery(B), recovery(C), and delta = C - B;
  * McNemar exact test on discordant pairs (b = B-right&C-wrong, c = B-wrong&C-right);
  * a paired bootstrap 95% CI on the recovery delta;
  * the same for the false-junction (flattened/mis-placed) rate.
Recovery is ambiguity-aware: a read is RECOVERED iff its refined junction
normalizes (chimeric_consensus.normalize_junction, which slides within the
genomic repeat/HP ambiguity) to the SAME coord as the TRUE non-canonical junction
— so 'held anywhere in the A-run' counts as recovered, but 'flattened to the
canonical decoy past the run' does not.

Usage: paired_arm_test.py --work-dir mix_hp_out
"""
from __future__ import annotations
import argparse, os, sys, json, math, random
from collections import defaultdict

import pysam

_here = os.path.dirname(os.path.abspath(__file__))
_repo = os.path.abspath(os.path.join(_here, "..", "..", ".."))
if _repo not in sys.path:
    sys.path.insert(0, _repo)
from rectify.core.consensus.chimeric_consensus import normalize_junction


def load_genome(fa_path):
    g, name = {}, None
    for l in open(fa_path):
        if l.startswith(">"):
            name = l[1:].strip().split()[0]; g[name] = []
        else:
            g[name].append(l.strip())
    return {k: "".join(v) for k, v in g.items()}


def read_junctions(bam_path):
    """read_id -> (donor, acceptor) of the (first) N-op, or None if no N-op/unmapped."""
    out = {}
    b = pysam.AlignmentFile(bam_path, "rb")
    for r in b:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            out.setdefault(r.query_name, None)
            continue
        ref = r.reference_start
        jn = None
        for op, ln in (r.cigartuples or []):
            if op == 3:
                jn = (ref, ref + ln); break
            if op in (0, 2, 3, 7, 8):
                ref += ln
        out[r.query_name] = jn
    b.close()
    return out


def mcnemar_p(b, c):
    """Two-sided exact McNemar p (binomial on discordant pairs)."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    # 2 * sum_{i=0..k} C(n,i) 0.5^n, clipped at 1
    s = sum(math.comb(n, i) for i in range(0, k + 1)) * (0.5 ** n)
    return min(1.0, 2.0 * s)


def boot_ci(pairs, rng, reps=2000):
    """Paired bootstrap 95% CI on mean(C_rec - B_rec) over reads."""
    if not pairs:
        return (0.0, 0.0)
    n = len(pairs)
    deltas = [c - b for (b, c) in pairs]
    out = []
    for _ in range(reps):
        s = 0.0
        for _ in range(n):
            s += deltas[rng.randrange(n)]
        out.append(s / n)
    out.sort()
    return (out[int(0.025 * reps)], out[int(0.975 * reps)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--seed", type=int, default=13)
    args = ap.parse_args()
    wd = args.work_dir
    rng = random.Random(args.seed)

    genome = load_genome(os.path.join(wd, "sim_ref.fa"))
    # read_truth columns
    rt_path = os.path.join(wd, "read_truth.tsv")
    hdr = open(rt_path).readline().rstrip("\n").split("\t")
    idx = {h: i for i, h in enumerate(hdr)}
    truth = {}
    for line in open(rt_path):
        if line.startswith("read_id\t"):
            continue
        f = line.rstrip("\n").split("\t")
        truth[f[idx["read_id"]]] = f

    jn_B = read_junctions(os.path.join(wd, "arm_B.bam"))
    jn_C = read_junctions(os.path.join(wd, "arm_C.bam"))

    # group cryptic reads (R3, has_true_junction=1) by context (HP length label)
    cells = defaultdict(list)   # ctx -> list of (b_rec, c_rec)
    for rid, f in truth.items():
        # scored cells: R3 cryptic non-canonical (hp_power/full) + R0flank canonical
        # flanking-HP (hp_dist). Both carry a true junction whose exact placement is
        # the arm-B-vs-arm-C question.
        if f[idx["motif_rung"]] not in ("R3", "R0flank") or f[idx["has_true_junction"]] != "1":
            continue
        chrom = f[idx["chrom"]]
        td, ta = int(f[idx["true_donor"]]), int(f[idx["true_acceptor"]])
        g = genome[chrom]
        true_norm = normalize_junction(td, ta, g)
        def rec(jn):
            if jn is None:
                return 0
            return 1 if normalize_junction(jn[0], jn[1], g) == true_norm else 0
        b_rec = rec(jn_B.get(rid))
        c_rec = rec(jn_C.get(rid))
        cells[f[idx["context"]]].append((b_rec, c_rec))

    import re
    def hp_key(ctx):
        # numeric sweep position: hp_dist '..._D<dist>' or hp_power 'HP<len>'; else 0
        m = re.search(r'D(\d+)$', ctx) or re.search(r'HP(\d+)', ctx)
        return int(m.group(1)) if m else 0

    def family(ctx):
        # placement family for grouping the dose-response: 'ACC_A','DON_A','BOT_A',
        # 'ACC_T'... (hp_dist, strip the trailing _D<dist>); 'HP' (hp_power); 'plain'.
        m = re.match(r'([A-Z]{3}_[ACGT])_D\d+$', ctx)
        if m:
            return m.group(1)
        return "HP" if ctx.startswith("HP") else ctx

    print(f"{'cell':14s} {'pos':>3s} {'n':>5s} {'recB':>6s} {'recC':>6s} {'dC-B':>7s} "
          f"{'95%CI':>16s} {'McNemar_p':>10s} {'b':>4s} {'c':>4s}")
    print("-" * 92)
    rows = []
    for ctx in sorted(cells, key=lambda c: (family(c), hp_key(c))):
        pairs = cells[ctx]
        n = len(pairs)
        recB = sum(b for b, _ in pairs) / n
        recC = sum(c for _, c in pairs) / n
        delta = recC - recB
        b = sum(1 for x, y in pairs if x == 1 and y == 0)   # B right, C wrong
        c = sum(1 for x, y in pairs if x == 0 and y == 1)   # B wrong, C right
        p = mcnemar_p(b, c)
        lo, hi = boot_ci(pairs, rng)
        L = hp_key(ctx)
        print(f"{ctx:14s} {L:3d} {n:5d} {recB:6.3f} {recC:6.3f} {delta:+7.3f} "
              f"[{lo:+.3f},{hi:+.3f}] {p:10.4f} {b:4d} {c:4d}")
        rows.append({"cell": ctx, "family": family(ctx), "pos": L, "n": n,
                     "recovery_B": round(recB, 4), "recovery_C": round(recC, 4),
                     "delta_C_minus_B": round(delta, 4),
                     "ci95": [round(lo, 4), round(hi, 4)], "mcnemar_p": round(p, 5),
                     "discordant_B_only": b, "discordant_C_only": c})
    out = {"per_cell": rows}
    with open(os.path.join(wd, "paired_arm_test.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    # dose-response summary, PER placement family (the arm-C niche map)
    print("\nDOSE-RESPONSE (arm-C − arm-B recovery vs sweep position, per family):")
    fams = {}
    for r in rows:
        fams.setdefault(r["family"], []).append(r)
    for fam in sorted(fams):
        fr = sorted(fams[fam], key=lambda r: r["pos"])
        if len(fr) < 2:
            continue
        pos = [r["pos"] for r in fr]
        deltas = [r["delta_C_minus_B"] for r in fr]
        peak = max(fr, key=lambda r: r["delta_C_minus_B"])
        shape = ("monotone-up" if all(deltas[i] <= deltas[i + 1] + 1e-9 for i in range(len(deltas) - 1))
                 else "hump/other")
        print(f"  {fam:8s} pos={pos} -> dC-B={deltas}  peak@pos{peak['pos']}"
              f"(Δ{peak['delta_C_minus_B']:+.3f}, p={peak['mcnemar_p']}) [{shape}]")


if __name__ == "__main__":
    main()
