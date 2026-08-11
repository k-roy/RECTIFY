#!/usr/bin/env python
"""644g — phase-1 per-read contest: does consensus refute mapPacBio's junk and
keep its gold? Plus the 5'-clip routing budget.

Kevin's framing (2026-08-10): junk at the ARM level is the wrong denominator —
Rectify phase 1 picks the best alignment per READ, so mapPacBio junk on a read
another aligner aligns better never survives. The questions that matter:

  Q1  Of the reads supporting mapPacBio's beyond-(mm2|resolver) JUNK junctions,
      what fraction does minimap2 beat on a fair per-read score? (junk refuted
      = "of no concern"). Per-junction: how many junk junctions retain >=1 /
      >=2 mpb-winning reads (phase-1 survivors that phase 2 must then judge)?
  Q2  Of the reads supporting the 37 residual GOLD junctions, what fraction
      does mapPacBio WIN? (If gold reads lose the contest, consensus would
      discard mapPacBio's true positives even with it in the panel.)
  Q3  Routing: what fraction of residual-gold-supporting reads are 5'- (or
      any-) soft-clipped >=8/>=15 bp in minimap2's alignment, and what fraction
      of ALL primary reads would such routing send to mapPacBio (compute
      budget)?

Score (flat hp_ed-style proxy, identical for both arms; N-ops free):
  score = NM + soft/hard-clip bases        [NM = subs + indel bases]
Caveat recorded: production phase 1 contests CORRECTED alignments with the
empirical table; this proxy contests raw alignments with flat costs.
"""

import json
import sys
from collections import defaultdict

import importlib.util
spec = importlib.util.spec_from_file_location(
    "t3", "/u/scratch/k/kevinroy/644_accept/644_t3_score.py")
t3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t3)

import pysam

BASE = "/u/scratch/k/kevinroy/644_accept/t3/full"
MM2 = f"{BASE}/bams/gold_trimmed.minimap2.bam"
RES = f"{BASE}/resolver.v5.a0.01.bam"
MPB = f"{BASE}/bams/gold_trimmed.mapPacBio.bam"
GOLD_TSV = "/u/scratch/k/kevinroy/644_accept/617_gold_introns.tsv"
OUT = f"{BASE}/644g_phase1_contest.json"
MIN_ANCHOR = 8

print("[644g] loading genome + catalogues...", flush=True)
G = t3.load_fa(t3.FA)
ann = t3.annotated_introns(G)
gold = {}
with open(GOLD_TSV) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    ci = {c: i for i, c in enumerate(header)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        gold[t3.canon(G, f[ci["chrom"]], int(f[ci["start"]]), int(f[ci["end"]]))] = f[ci["source"]]

def junctions_of(read):
    out = []
    ops = read.cigartuples
    ref = read.reference_start
    for i, (op, ln) in enumerate(ops):
        if op == 3:
            lf = ops[i-1][1] if i >= 1 and ops[i-1][0] in (0, 7, 8) else 0
            rt = ops[i+1][1] if i+1 < len(ops) and ops[i+1][0] in (0, 7, 8) else 0
            if min(lf, rt) >= MIN_ANCHOR:
                out.append((read.reference_name, ref, ref + ln))
        if op in (0, 2, 3, 7, 8):
            ref += ln
    return out

def primary(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary) \
        and read.cigartuples is not None

def flat_score(read):
    nm = read.get_tag("NM") if read.has_tag("NM") else 0
    clips = sum(ln for op, ln in read.cigartuples if op in (4, 5))
    return nm + clips

def clip53(read):
    ct = read.cigartuples
    left = ct[0][1] if ct[0][0] in (4, 5) else 0
    right = ct[-1][1] if ct[-1][0] in (4, 5) else 0
    return (right, left) if read.is_reverse else (left, right)

# pass 1: mm2 + resolver canonical junction sets (for the beyond frame)
def canon_set(path):
    S = set()
    with pysam.AlignmentFile(path, "rb", check_sq=False) as fh:
        for r in fh.fetch(until_eof=True):
            if primary(r):
                for j in junctions_of(r):
                    S.add(t3.canon(G, *j))
    return S

print("[644g] censusing mm2 + resolver junction sets...", flush=True)
mm2_set = canon_set(MM2)
res_set = canon_set(RES)

# pass 2: mapPacBio — supporting reads + mpb score per beyond-both junction
print("[644g] walking mapPacBio for beyond-arm supporting reads...", flush=True)
jclass = {}                       # canonical j -> gold-class/junk
reads_of = defaultdict(set)       # canonical j -> read ids
mpb_read = {}                     # read id -> (score, n_beyond_juncs)
with pysam.AlignmentFile(MPB, "rb", check_sq=False) as fh:
    for r in fh.fetch(until_eof=True):
        if not primary(r):
            continue
        target = []
        for j in junctions_of(r):
            cj = t3.canon(G, *j)
            if cj in mm2_set or cj in res_set:
                continue
            if cj not in jclass:
                jclass[cj] = gold.get(cj) or ("annotated" if cj in ann else "junk")
            target.append(cj)
        if target:
            rid = r.query_name
            for cj in target:
                reads_of[cj].add(rid)
            mpb_read[rid] = (flat_score(r), len(target))
print(f"[644g] beyond-both junctions={len(jclass)} supporting reads={len(mpb_read)}",
      flush=True)

# pass 3: mm2 — contest scores + clip status for the target reads; global clip budget
print("[644g] walking minimap2 for contest + routing...", flush=True)
mm2_read = {}
tot = clip5_8 = clip5_15 = clipany_8 = clipany_15 = 0
with pysam.AlignmentFile(MM2, "rb", check_sq=False) as fh:
    for r in fh.fetch(until_eof=True):
        if not primary(r):
            continue
        c5, c3 = clip53(r)
        tot += 1
        if c5 >= 8: clip5_8 += 1
        if c5 >= 15: clip5_15 += 1
        if max(c5, c3) >= 8: clipany_8 += 1
        if max(c5, c3) >= 15: clipany_15 += 1
        if r.query_name in mpb_read:
            mm2_read[r.query_name] = (flat_score(r), c5, c3)

def contest(rid):
    """Return 'mpb_wins' / 'mm2_wins' / 'tie' / 'mm2_absent'."""
    if rid not in mm2_read:
        return "mpb_wins_absent"
    s_mpb = mpb_read[rid][0]
    s_mm2 = mm2_read[rid][0]
    if s_mpb < s_mm2:
        return "mpb_wins"
    if s_mm2 < s_mpb:
        return "mm2_wins"
    return "tie"

res_rows = {}
for cls in ("junk", "gold"):
    sel = {j: r for j, r in reads_of.items()
           if (jclass[j] != "junk") == (cls == "gold") and jclass[j] != "annotated"}
    read_verdict = defaultdict(int)
    surv1 = surv2 = 0
    routed_reads_8 = routed_reads_15 = read_n = 0
    for j, rids in sel.items():
        wins = 0
        for rid in rids:
            v = contest(rid)
            read_verdict[v] += 1
            read_n += 1
            if v.startswith("mpb_wins") or v == "tie":
                wins += 1
            mm = mm2_read.get(rid)
            if mm is None or mm[1] >= 8:
                routed_reads_8 += 1
            if mm is None or mm[1] >= 15:
                routed_reads_15 += 1
        if wins >= 1: surv1 += 1
        if wins >= 2: surv2 += 1
    res_rows[cls] = {
        "junctions": len(sel),
        "supporting_read_slots": read_n,
        "read_verdicts": dict(read_verdict),
        "junctions_surviving_ge1_win_or_tie": surv1,
        "junctions_surviving_ge2_win_or_tie": surv2,
        "read_slots_routed_5clip_ge8_or_mm2absent": routed_reads_8,
        "read_slots_routed_5clip_ge15_or_mm2absent": routed_reads_15,
    }
    print(f"[644g] {cls}: {res_rows[cls]}", flush=True)

out = {
    "score_def": "NM + soft/hard clip bases (flat; N free); lower wins; raw-alignment proxy",
    "classes": res_rows,
    "routing_budget": {
        "primary_reads": tot,
        "clip5_ge8": clip5_8, "clip5_ge15": clip5_15,
        "clipany_ge8": clipany_8, "clipany_ge15": clipany_15,
    },
    "gold_junction_detail": [
        [list(j), jclass[j], len(reads_of[j]),
         sum(1 for rid in reads_of[j] if contest(rid).startswith("mpb_wins") or contest(rid) == "tie")]
        for j in sorted(reads_of, key=lambda x: -len(reads_of[x])) if jclass[j] not in ("junk", "annotated")
    ],
}
with open(OUT, "w") as fh:
    json.dump(out, fh, indent=1)
print(f"[644g] wrote {OUT}", flush=True)
