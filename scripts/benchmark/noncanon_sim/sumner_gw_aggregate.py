#!/usr/bin/env python3
"""Aggregate the per-sample genome-wide Sumner discovery outputs -> SMA vs WT.
Reports: per-sample revealed-novel-noncanon yield; RECURRENT revealed junctions (in >=N samples) as the
confident-candidate set; SMA-specific vs WT-specific vs shared. HONEST CAVEAT: on real data with no ground
truth, 'revealed novel non-canonical' is a DISCOVERY YIELD = real flattening-recovery + fabrication; the
FABRICATION RATE is calibrated by the spike-in track, not here. Recurrence + read-support are the only
real-data confidence signals available.
Usage: python sumner_gw_aggregate.py <panel_dir>
"""
import sys, glob, os
from collections import defaultdict

PANEL = sys.argv[1] if len(sys.argv) > 1 else "/scratch/users/kevinroy/sumner_gw/panel"

# per-sample summaries
summ = {}
for f in sorted(glob.glob(os.path.join(PANEL, "*.summary.tsv"))):
    d = {}
    for line in open(f):
        p = line.rstrip("\n").split("\t")
        for i in range(0, len(p) - 1, 2):
            d[p[i]] = p[i + 1]
    summ[d.get("sample", os.path.basename(f))] = d

def grp(s): return "SMA" if s.startswith("SMA") else ("WT" if s.startswith("WT") else "?")

print("=== per-sample revealed novel-non-canonical (genome-wide, 5% downsample) ===")
print(f"{'sample':20s} {'grp':4s} {'sub_reads':>10s} {'revealed':>9s} {'raw_nn':>7s} {'ref_nn':>7s}")
for s in sorted(summ):
    d = summ[s]
    print(f"{s:20s} {grp(s):4s} {int(d.get('sub_reads',0)):10d} {int(d.get('revealed_novel_noncanon',0)):9d} "
          f"{int(d.get('raw_novel_noncanon',0)):7d} {int(d.get('ref_novel_noncanon',0)):7d}")

# recurrence of revealed junctions across samples, split by group
jx_samples = defaultdict(set)   # (chrom,donor,acceptor) -> set(samples)
for f in sorted(glob.glob(os.path.join(PANEL, "*.revealed_noncanon.tsv"))):
    for line in open(f):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            jx_samples[(p[0], int(p[1]), int(p[2]))].add(p[3])

nsma = sum(1 for s in summ if grp(s) == "SMA"); nwt = sum(1 for s in summ if grp(s) == "WT")
print(f"\n=== recurrence ({nsma} SMA, {nwt} WT samples); distinct revealed junctions: {len(jx_samples)} ===")
for minrep in (2, 3, 5):
    rec = {j: ss for j, ss in jx_samples.items() if len(ss) >= minrep}
    sma_only = sum(1 for j, ss in rec.items() if all(grp(x) == "SMA" for x in ss))
    wt_only = sum(1 for j, ss in rec.items() if all(grp(x) == "WT" for x in ss))
    shared = len(rec) - sma_only - wt_only
    print(f"  in >={minrep} samples: {len(rec):6d}  (SMA-only {sma_only}, WT-only {wt_only}, shared {shared})")

# the most-recurrent SMA-specific candidates (biology leads)
print("\n=== top SMA-specific recurrent revealed non-canonical junctions (biology leads) ===")
sma_rec = [(len(ss), j) for j, ss in jx_samples.items()
           if len(ss) >= 3 and all(grp(x) == "SMA" for x in ss)]
sma_rec.sort(reverse=True)
for n, (c, d, a) in sma_rec[:20]:
    print(f"  {c}:{d}-{a}  in {n} SMA samples")
print("\nNOTE: precision (real vs fabricated) is the SPIKE-IN track's job; here recurrence+group-specificity "
      "are the confidence signals. SMN region (chr5:~70.0-70.95Mb) = built-in positive control.")
