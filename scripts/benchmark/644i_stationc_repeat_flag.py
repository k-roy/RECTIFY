#!/usr/bin/env python
"""644i — Station-C repeat-context flag prototype (the PI's 644h correction).

Kevin (2026-08-11): the copy-misplacement blindness of the overhang-quality
score is not a major problem — Station C should carry a direct repeat flag
(rDNA, tandem arrays with CDS, tandemly arrayed families like CUP1/ENA,
Ty/LTR, subtelomeres). This prototypes the flag two ways and measures the
composite against the 644h census (local 644h.json copy):

1. ANNOTATION flag: GFF repeat-class features +/- MARGIN — rRNA_gene (rDNA),
   LTR_retrotransposon, long_terminal_repeat, transposable_element_gene,
   telomere, telomeric_repeat, X_element, X_element_combinatorial_repeat,
   Y_prime_element, tRNA_gene (dispersed identical copies).
2. VARIANT-MULTIPLICITY flag (label-free, catches CDS tandem arrays like
   CUP1/ENA that annotation cannot): cluster ALL beyond-arm non-canonical
   junctions (no truth labels) at CLUSTER_BP linkage; flag loci carrying
   >= MULTI_K distinct junction variants. Real novel junctions present as
   1-2 variants per locus; systematic repeat misalignment presents as dozens.

Flag semantics per the PI: DEMOTE, never discard — flagged candidates need
orthogonal evidence; unflagged candidates pass on (q, recurrence) alone.

Usage: python 644i_stationc_repeat_flag.py [--json PATH] [--gff PATH]
Writes 644i_stationc_repeat_flag.json next to the input census copy.
"""

import argparse
import gzip
import json
from collections import defaultdict
from pathlib import Path

MARGIN = 500          # bp around an annotated repeat feature
CLUSTER_BP = 1000     # linkage for the label-free variant clusters
MULTI_K = 5           # distinct non-canonical variants at a locus => flagged
RDNA = ("chrXII", 450000, 491000)
REPEAT_TYPES = {
    "rRNA_gene", "LTR_retrotransposon", "long_terminal_repeat",
    "transposable_element_gene", "telomere", "telomeric_repeat",
    "X_element", "X_element_combinatorial_repeat", "Y_prime_element",
    "tRNA_gene",
}
Q_GRID = [0, 40, 60, 80, 90]


def load_repeat_intervals(gff_path):
    iv = defaultdict(list)
    op = gzip.open if str(gff_path).endswith(".gz") else open
    with op(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[2] not in REPEAT_TYPES:
                continue
            iv[f[0]].append((int(f[3]) - 1 - MARGIN, int(f[4]) + MARGIN, f[2]))
    c, s, e = RDNA
    iv[c].append((s, e, "rDNA_region"))
    for chrom in iv:
        iv[chrom].sort()
    return iv


def flag_of(iv, chrom, s, e):
    for a, b, t in iv.get(chrom, ()):
        if s < b and e > a:
            return t
    return None


ap = argparse.ArgumentParser()
ap.add_argument("--json", default=None)
ap.add_argument("--gff", default=None)
ap.add_argument("--selfhom-paf", default=None,
                help="genome self-homology PAF (minimap2 -DP -k19 -w19 -m200 "
                     "G G); intervals +/-200 bp become a second repeat flag")
args = ap.parse_args()
HERE = Path(__file__).resolve()
census_path = Path(args.json) if args.json else HERE.parent / "644h.json"
gff_path = Path(args.gff) if args.gff else (
    HERE.parents[2] / "rectify/data/genomes/saccharomyces_cerevisiae/"
    "saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz")

d = json.load(open(census_path))
iv = load_repeat_intervals(gff_path)
n_feat = sum(len(v) for v in iv.values())
print(f"[644i] repeat intervals: {n_feat} (margin {MARGIN} bp) from {gff_path.name}")

selfhom = defaultdict(list)
if args.selfhom_paf:
    for line in open(args.selfhom_paf):
        f = line.split("\t")
        selfhom[f[0]].append((int(f[2]) - 200, int(f[3]) + 200, "selfhom"))
        selfhom[f[5]].append((int(f[7]) - 200, int(f[8]) + 200, "selfhom"))
    for c in selfhom:
        selfhom[c].sort()
    print(f"[644i] self-homology intervals loaded from {args.selfhom_paf}")

rows = [dict(r, is_gold=True, cls=r["class"]) for r in d["gold_detail"]] + \
       [dict(r, is_gold=False, cls="junk") for r in d["junk_junctions"]]

# label-free variant-multiplicity clusters over ALL non-canonical junctions
nc = sorted((r for r in rows if not r["canon_in_class"]),
            key=lambda r: (r["junction"][0], r["junction"][1]))
clusters = []
for r in nc:
    c, s, e = r["junction"]
    if clusters and clusters[-1]["chrom"] == c and s <= clusters[-1]["end"] + CLUSTER_BP:
        cl = clusters[-1]
        cl["end"] = max(cl["end"], e)
        cl["members"].append(r)
    else:
        clusters.append({"chrom": c, "start": s, "end": e, "members": [r]})
JITTER_BP = 10   # variants with both boundaries within this are ONE junction
# (aligner jitter around a real junction sprays near-duplicate variants; repeat
# misalignment sprays DISTINCT junctions — count after collapsing jitter)


def collapsed_count(members):
    reps = []
    for r in sorted(members, key=lambda r: (r["junction"][1], r["junction"][2])):
        _, s, e = r["junction"]
        for rep in reps:
            if abs(s - rep[0]) <= JITTER_BP and abs(e - rep[1]) <= JITTER_BP:
                break
        else:
            reps.append((s, e))
    return len(reps)


multi_loci = [cl for cl in clusters if collapsed_count(cl["members"]) >= MULTI_K]
print(f"[644i] label-free clusters: {len(clusters)}; "
      f">= {MULTI_K} jitter-collapsed variants: {len(multi_loci)} loci covering "
      f"{sum(len(c['members']) for c in multi_loci)} junctions")

for r in rows:
    c, s, e = r["junction"]
    r["ann_flag"] = flag_of(iv, c, s, e)
    r["multi_flag"] = False
for cl in multi_loci:
    for r in cl["members"]:
        r["multi_flag"] = True

# composite sweep: non-canonical track (the one 644h could not purify)
out_rows = []
print(f"\n[644i] NON-CANONICAL track, phase1-surv>=1 frame "
      f"(644h baseline: q>=80 -> 11/14 gold, 1537 junk):")
sel = [r for r in rows if not r["canon_in_class"] and r["surv_wins"] >= 1]
for qt in Q_GRID:
    base_g = [r for r in sel if r["is_gold"] and r["q_max"] >= qt]
    base_k = [r for r in sel if not r["is_gold"] and r["q_max"] >= qt]
    for flags in ("ann", "ann+selfhom", "ann+multi"):
        def unflagged(r):
            if r["ann_flag"]:
                return False
            if flags == "ann+selfhom" and flag_of(selfhom, *r["junction"]):
                return False
            if flags == "ann+multi" and r["multi_flag"]:
                return False
            return True
        g = [r for r in base_g if unflagged(r)]
        k = [r for r in base_k if unflagged(r)]
        g2 = [r for r in g if r["surv_wins"] >= 2]
        k2 = [r for r in k if r["surv_wins"] >= 2]
        out_rows.append({
            "q_thr": qt, "flags": flags,
            "gold_pass": len(g), "junk_pass": len(k),
            "gold_demoted": len(base_g) - len(g), "junk_demoted": len(base_k) - len(k),
            "gold_pass_recur2": len(g2), "junk_pass_recur2": len(k2),
        })
        print(f"  q>={qt:>2} [{flags:9}]: PASS gold {len(g):>2}/{len(base_g)} "
              f"junk {len(k):>5}/{len(base_k)}   +recur>=2: gold {len(g2)} junk {len(k2)}")

demoted_gold = [
    {"junction": r["junction"], "cls": r["cls"], "ann_flag": r["ann_flag"],
     "multi_flag": r["multi_flag"], "q_max": r["q_max"], "surv": r["surv_wins"]}
    for r in rows if r["is_gold"] and not r["canon_in_class"]
    and (r["ann_flag"] or r["multi_flag"])]
print("\n[644i] demoted gold (needs orthogonal evidence, NOT discarded):")
for g in demoted_gold:
    print(f"  {g['junction']} {g['cls'][:24]} ann={g['ann_flag']} "
          f"multi={g['multi_flag']} q={g['q_max']} surv={g['surv']}")

outp = census_path.parent / "644i_stationc_repeat_flag.json"
json.dump({
    "params": {"MARGIN": MARGIN, "CLUSTER_BP": CLUSTER_BP, "MULTI_K": MULTI_K,
               "REPEAT_TYPES": sorted(REPEAT_TYPES), "RDNA": RDNA},
    "n_repeat_intervals": n_feat,
    "n_multi_loci": len(multi_loci),
    "sweep_noncanonical_surv1": out_rows,
    "demoted_gold": demoted_gold,
}, open(outp, "w"), indent=1)
print(f"\n[644i] wrote {outp}")
