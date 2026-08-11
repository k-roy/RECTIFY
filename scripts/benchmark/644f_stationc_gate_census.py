#!/usr/bin/env python
"""644f — Station-C gate census: can pool-level arbitration harvest mapPacBio?

Question (Kevin, 2026-08-10): of the gold junctions the resolver misses, how
many would Station-C triggers flag from a mapPacBio alignment — and how many
of the 5,456 "junk" junctions would come with them? I.e. does judicious
pool-level arbitration (read-support recurrence + motif-class as a SOFT
signal) turn mapPacBio into a high-value probe?

Method: same census + ambiguity-canonicalised matching as 644_t3_score.py
(primary reads, min anchor 8, leftmost-canonical coords; support summed
across exact variants of one canonical junction). Compute mapPacBio's
beyond-mm2 and beyond-(mm2 ∪ resolver_v5) canonical junction sets, classify
each junction {gold / annotated / junk}, then sweep the support gate
(>=1,2,3,5,10) with a canonical-in-class split. This is the
"support-stratified re-census" named in 644's follow-ups.
"""

import json
import re
import sys
from collections import defaultdict

sys.path.insert(0, "/u/scratch/k/kevinroy/644_accept")
import importlib.util
spec = importlib.util.spec_from_file_location(
    "t3score", "/u/scratch/k/kevinroy/644_accept/644_t3_score.py")
t3 = importlib.util.module_from_spec(spec)
# 644_t3_score.py runs main() under __main__ only; safe to import.
spec.loader.exec_module(t3)

BASE = "/u/scratch/k/kevinroy/644_accept/t3/full"
ARMS = {
    "minimap2": f"{BASE}/bams/gold_trimmed.minimap2.bam",
    "resolver_v5": f"{BASE}/resolver.v5.a0.01.bam",
    "mapPacBio": f"{BASE}/bams/gold_trimmed.mapPacBio.bam",
}
GOLD_TSV = "/u/scratch/k/kevinroy/644_accept/617_gold_introns.tsv"
OUT = f"{BASE}/644f_stationc_gate.json"

print("[644f] loading genome + annotation...", flush=True)
G = t3.load_fa(t3.FA)
ann = t3.annotated_introns(G)

# gold catalogue: canonical coords + class column
gold = {}
with open(GOLD_TSV) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    ci = {c: i for i, c in enumerate(header)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        c, s, e = f[ci["chrom"]], int(f[ci["start"]]), int(f[ci["end"]])
        cls = f[ci["source"]] if "source" in ci else "?"
        gold[t3.canon(G, c, s, e)] = cls
print(f"[644f] gold={len(gold)} annotated={len(ann)}", flush=True)

# census each arm; canonical support = sum over exact variants
carm = {}
for name, path in ARMS.items():
    J, support = t3.census_bam(path)
    csup = defaultdict(int)
    for j, n in support.items():
        csup[t3.canon(G, *j)] += n
    carm[name] = csup
    print(f"[644f] {name}: {len(J)} exact / {len(csup)} canonical junctions",
          flush=True)

mm2 = set(carm["minimap2"])
res = set(carm["resolver_v5"])
mpb = carm["mapPacBio"]

def classify(j):
    if j in gold:
        return "gold"
    if j in ann:
        return "annotated"
    return "junk"

def sweep(beyond, label):
    rows = []
    for thr in (1, 2, 3, 5, 10):
        kept = {j: n for j, n in beyond.items() if n >= thr}
        by = defaultdict(list)
        for j, n in kept.items():
            by[classify(j)].append((j, n))
        junk = by.get("junk", [])
        junk_canon = sum(1 for j, _ in junk if t3.canon_in_class(G, *j))
        rows.append({
            "min_support": thr,
            "total": len(kept),
            "gold": len(by.get("gold", [])),
            "annotated": len(by.get("annotated", [])),
            "junk": len(junk),
            "junk_canonical_in_class": junk_canon,
            "gold_list": sorted(
                [[list(j), n, gold[j]] for j, n in by.get("gold", [])],
                key=lambda r: -r[1]),
        })
        print(f"[644f] {label} support>={thr}: total={len(kept)} "
              f"gold={rows[-1]['gold']} ann={rows[-1]['annotated']} "
              f"junk={rows[-1]['junk']} (canon-class {junk_canon})", flush=True)
    return rows

beyond_mm2 = {j: n for j, n in mpb.items() if j not in mm2}
beyond_both = {j: n for j, n in mpb.items() if j not in mm2 and j not in res}
print(f"[644f] beyond-mm2={len(beyond_mm2)}  beyond-(mm2|resolver_v5)={len(beyond_both)}",
      flush=True)

out = {
    "beyond_mm2": sweep(beyond_mm2, "beyond-mm2"),
    "beyond_mm2_and_resolver_v5": sweep(beyond_both, "beyond-both"),
    "arm_canonical_junctions": {k: len(v) for k, v in carm.items()},
}
with open(OUT, "w") as fh:
    json.dump(out, fh, indent=1)
print(f"[644f] wrote {OUT}", flush=True)
