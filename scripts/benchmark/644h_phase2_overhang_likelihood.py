#!/usr/bin/env python
"""644h — phase-2 overhang-quality likelihood prototype (Station C's scorer).

The PI's metric (STATIONC_MAPPACBIO_HARVEST_20260810.md §PI-reframe): the
single most important per-junction signal is the longest HIGH-QUALITY overhang
on the SHORT-EXON side — high complexity x length x low error rate. One golden
read with a 40 bp clean, high-complexity short-exon overhang outweighs twenty
reads with 9 bp overhangs in an A-run. This script prototypes that scorer and
asks the 644g question: does it separate the residual gold (37; 14 Gould
non-canonical) from the junk that support-recurrence and motif could not
(5,414 total; 1,779 recurrent non-canonical)?

Per supporting read of each mapPacBio beyond-(mm2 | resolver_v5) junction:
  - short-exon side = the side of the N-op with the SMALLER contiguous aligned
    anchor (query bases to the next N / clip / read end) — the per-read proxy
    for the short exon;
  - overhang = up to CAP=60 junction-adjacent query bases on that side;
  - I_eff   = effective information bits of the overhang (copied verbatim from
    641's overhang_informativeness.effective_information_bits — the resolver's
    own currency: min(HP-discounted composition entropy, order-1 conditional
    entropy x length));
  - errors  = mismatches + inserted bases + deleted reference bases within the
    capped window (alignment vs genome; no MD tag needed);
  - q_read  = max(0, I_eff - ERR_BITS * errors)   [ERR_BITS = 2.0: an errored
    base carries no placement evidence and costs ~2 bits of discrimination].

Junction score q_max = max over supporting reads of q_read.

Calibration positive class: ANNOTATED junctions in the SAME mapPacBio arm
(same error model, same extraction), capped at CAL_READS_PER_J=10 reads per
junction (support-cap noted in the findings doc; gold junctions have <= ~34).

Frames reported (matching 644f/644g): all beyond-both; phase-1 survivors
(>=1 win-or-tie read under the 644g flat contest); recurrent survivors (>=2).
Tracks: canonical-in-class vs non-canonical (the track support cannot purify).

Output: t3/full/644h_phase2_overhang.json + printed sweep tables.
Smoke: --limit-reads N truncates every BAM pass (same code path).
"""

import argparse
import json
import math
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
OUT = f"{BASE}/644h_phase2_overhang.json"

MIN_ANCHOR = 8
CAP = 60                 # junction-adjacent query bases assessed per side
ERR_BITS = 2.0           # bits removed per alignment error in the window
CAL_READS_PER_J = 10     # annotated-junction calibration read cap
HP_DISCOUNT = 0.5
THRESH_GRID = [0, 10, 20, 30, 40, 50, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110]


# --- I_eff: copied verbatim from 641's overhang_informativeness ------------
_ACGT = frozenset("ACGT")


def _entropy_bits(counts, total):
    if total <= 0:
        return 0.0
    h = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            h -= p * math.log2(p)
    return h


def effective_information_bits(seq, hp_discount=HP_DISCOUNT):
    s = [b for b in seq.upper() if b in _ACGT]
    n = len(s)
    if n < 2:
        return 0.0
    eff_counts = {}
    l_eff = 0.0
    i = 0
    while i < n:
        j = i
        while j < n and s[j] == s[i]:
            j += 1
        run_eff = 1.0 + hp_discount * (j - i - 1)
        eff_counts[s[i]] = eff_counts.get(s[i], 0.0) + run_eff
        l_eff += run_eff
        i = j
    i_comp = l_eff * _entropy_bits(eff_counts, l_eff)

    base_counts = {}
    for b in s:
        base_counts[b] = base_counts.get(b, 0.0) + 1.0
    di_counts = {}
    for k in range(n - 1):
        d = s[k] + s[k + 1]
        di_counts[d] = di_counts.get(d, 0.0) + 1.0
    h_cond = max(0.0, _entropy_bits(di_counts, float(n - 1))
                 - _entropy_bits(base_counts, float(n)))
    return min(i_comp, n * h_cond)
# ---------------------------------------------------------------------------


def primary(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary) \
        and read.cigartuples is not None


def junctions_at(read):
    """[(cigar_index, (chrom, start, end)), ...] for anchor->=MIN_ANCHOR N-ops."""
    out = []
    ops = read.cigartuples
    ref = read.reference_start
    for i, (op, ln) in enumerate(ops):
        if op == 3:
            lf = ops[i - 1][1] if i >= 1 and ops[i - 1][0] in (0, 7, 8) else 0
            rt = ops[i + 1][1] if i + 1 < len(ops) and ops[i + 1][0] in (0, 7, 8) else 0
            if min(lf, rt) >= MIN_ANCHOR:
                out.append((i, (read.reference_name, ref, ref + ln)))
        if op in (0, 2, 3, 7, 8):
            ref += ln
    return out


def flat_score(read):
    nm = read.get_tag("NM") if read.has_tag("NM") else 0
    clips = sum(ln for op, ln in read.cigartuples if op in (4, 5))
    return nm + clips


def _offsets(read):
    """Prefix query/ref offsets at each cigar-op boundary."""
    ops = read.cigartuples
    qoff = [0]
    roff = [read.reference_start]
    for op, ln in ops:
        qoff.append(qoff[-1] + (ln if op in (0, 1, 4, 7, 8) else 0))
        roff.append(roff[-1] + (ln if op in (0, 2, 3, 7, 8) else 0))
    return qoff, roff


def _side_features(read, qoff, roff, i_n, chrom_seq, direction):
    """One side of the N-op at cigar index i_n.

    direction=+1: rightward (downstream in reference), -1: leftward.
    Returns dict(anchor_q, seq, errors, first_err) where seq/errors cover the
    capped junction-adjacent window; first_err is the collected-base index of
    the first error (None if the window is clean); anchor_q is the UNCAPPED
    contiguous query anchor to the next N / clip / read end.
    """
    ops = read.cigartuples
    qseq = read.query_sequence
    collected = []
    errors = 0
    first_err = None
    anchor_q = 0

    def note_err(n=1):
        nonlocal errors, first_err
        errors += n
        if first_err is None:
            first_err = len(collected)

    rng = range(i_n + 1, len(ops)) if direction > 0 else range(i_n - 1, -1, -1)
    for k in rng:
        op, ln = ops[k]
        if op in (3, 4, 5):
            break
        if op in (0, 7, 8):
            for t in range(ln):
                if direction > 0:
                    qi, ri = qoff[k] + t, roff[k] + t
                else:
                    qi, ri = qoff[k + 1] - 1 - t, roff[k + 1] - 1 - t
                if len(collected) < CAP:
                    b = qseq[qi].upper()
                    if b == "=":
                        # SAM '=' SEQ compression (mapPacBio): base identical
                        # to the reference by definition.
                        collected.append(chrom_seq[ri] if ri < len(chrom_seq)
                                         else "N")
                    else:
                        collected.append(b)
                        if ri >= len(chrom_seq) or b != chrom_seq[ri]:
                            note_err()
                else:
                    break
            anchor_q += ln
        elif op == 1:
            for t in range(ln):
                if len(collected) < CAP:
                    qi = (qoff[k] + t) if direction > 0 else (qoff[k + 1] - 1 - t)
                    b = qseq[qi].upper()
                    # inserted base: '=' is undefined off-reference — mask it
                    collected.append("N" if b == "=" else b)
                    note_err()
                else:
                    break
            anchor_q += ln
        elif op == 2:
            if len(collected) < CAP:
                note_err(ln)
        if anchor_q >= 10 ** 9:
            break
    seq = "".join(collected)
    if direction < 0:
        seq = seq[::-1]  # junction-adjacent-first is index 0 either way for
        # clean-prefix bookkeeping; reversal is a no-op for I_eff (entropy is
        # orientation-insensitive) — kept for readable per-read dumps.
    return {"anchor_q": anchor_q, "seq": seq, "errors": errors,
            "first_err": first_err}


def read_junction_q(read, qoff, roff, i_n, chrom_seq):
    """Short-exon-side quality for one junction of one read."""
    left = _side_features(read, qoff, roff, i_n, chrom_seq, -1)
    right = _side_features(read, qoff, roff, i_n, chrom_seq, +1)
    short = left if left["anchor_q"] <= right["anchor_q"] else right
    i_eff = effective_information_bits(short["seq"])
    q = max(0.0, i_eff - ERR_BITS * short["errors"])
    clean = short["first_err"] if short["first_err"] is not None \
        else min(short["anchor_q"], CAP)
    return {
        "q": round(q, 2),
        "i_eff": round(i_eff, 2),
        "errors": short["errors"],
        "clean_prefix": clean,
        "anchor_short": short["anchor_q"],
        "anchor_long": max(left["anchor_q"], right["anchor_q"]),
        "side": "L" if short is left else "R",
        "seq": short["seq"],
    }


def canon_set(path, limit):
    S = set()
    with pysam.AlignmentFile(path, "rb", check_sq=False) as fh:
        for n, r in enumerate(fh.fetch(until_eof=True)):
            if limit and n >= limit:
                break
            if primary(r):
                for _, j in junctions_at(r):
                    S.add(t3.canon(G, *j))
    return S


def quantiles(vals, qs=(0.05, 0.25, 0.5, 0.75, 0.95)):
    if not vals:
        return {}
    v = sorted(vals)
    return {f"p{int(q * 100)}": round(v[min(len(v) - 1, int(q * len(v)))], 2)
            for q in qs}


ap = argparse.ArgumentParser()
ap.add_argument("--limit-reads", type=int, default=0,
                help="truncate every BAM pass after N reads (smoke)")
ap.add_argument("--out", default=OUT)
args = ap.parse_args()
LIM = args.limit_reads

print("[644h] loading genome + catalogues...", flush=True)
G = t3.load_fa(t3.FA)
ann = t3.annotated_introns(G)
gold = {}
with open(GOLD_TSV) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    ci = {c: i for i, c in enumerate(header)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        gold[t3.canon(G, f[ci["chrom"]], int(f[ci["start"]]),
                      int(f[ci["end"]]))] = f[ci["source"]]

print("[644h] censusing mm2 + resolver junction sets...", flush=True)
mm2_set = canon_set(MM2, LIM)
res_set = canon_set(RES, LIM)

print("[644h] walking mapPacBio: overhang features for beyond-both "
      "supporting reads + annotated calibration...", flush=True)
jclass = {}                        # cj -> gold-source / junk
jreads = defaultdict(list)        # cj -> [read-level q dicts]
reads_of = defaultdict(set)       # cj -> read ids (for the phase-1 contest)
mpb_read = {}                     # read id -> flat score
ann_qs = defaultdict(list)        # annotated cj -> capped [q_read]
n_noseq = 0
with pysam.AlignmentFile(MPB, "rb", check_sq=False) as fh:
    for n, r in enumerate(fh.fetch(until_eof=True)):
        if LIM and n >= LIM:
            break
        if not primary(r):
            continue
        js = junctions_at(r)
        if not js:
            continue
        offs = None
        for i_n, j in js:
            cj = t3.canon(G, *j)
            chrom_seq = G.get(j[0])
            if chrom_seq is None:
                continue
            if cj in mm2_set or cj in res_set:
                if cj in ann and len(ann_qs[cj]) < CAL_READS_PER_J:
                    if r.query_sequence is None:
                        n_noseq += 1
                        continue
                    if offs is None:
                        offs = _offsets(r)
                    ann_qs[cj].append(
                        read_junction_q(r, offs[0], offs[1], i_n, chrom_seq)["q"])
                continue
            if cj not in jclass:
                jclass[cj] = gold.get(cj) or \
                    ("annotated" if cj in ann else "junk")
            if r.query_sequence is None:
                n_noseq += 1
                continue
            if offs is None:
                offs = _offsets(r)
            feats = read_junction_q(r, offs[0], offs[1], i_n, chrom_seq)
            jreads[cj].append(feats)
            reads_of[cj].add(r.query_name)
            mpb_read[r.query_name] = flat_score(r)
print(f"[644h] beyond-both junctions={len(jclass)} "
      f"supporting reads={len(mpb_read)} ann-cal junctions={len(ann_qs)} "
      f"reads-without-seq={n_noseq}", flush=True)

print("[644h] walking minimap2 for the phase-1 contest...", flush=True)
mm2_score = {}
with pysam.AlignmentFile(MM2, "rb", check_sq=False) as fh:
    for n, r in enumerate(fh.fetch(until_eof=True)):
        if LIM and n >= LIM:
            break
        if primary(r) and r.query_name in mpb_read:
            mm2_score[r.query_name] = flat_score(r)

surv = {}
for cj, rids in reads_of.items():
    wins = sum(1 for rid in rids
               if rid not in mm2_score or mpb_read[rid] <= mm2_score[rid])
    surv[cj] = wins

# --- aggregate + report ----------------------------------------------------
rows = []
for cj, feats in jreads.items():
    best = max(feats, key=lambda f: f["q"])
    qs = sorted((f["q"] for f in feats), reverse=True)
    rows.append({
        "junction": list(cj),
        "class": jclass[cj],
        "is_gold": jclass[cj] not in ("junk", "annotated"),
        "canon_in_class": t3.canon_in_class(G, *cj),
        "n_reads": len(feats),
        "q_max": qs[0],
        "q_2nd": qs[1] if len(qs) > 1 else None,
        "clean_prefix_max": max(f["clean_prefix"] for f in feats),
        "surv_wins": surv.get(cj, 0),
        "best_read": best,
    })

ann_max = [max(v) for v in ann_qs.values() if v]

frames = {
    "all": lambda r: True,
    "phase1_surv1": lambda r: r["surv_wins"] >= 1,
    "phase1_surv2": lambda r: r["surv_wins"] >= 2,
}
tracks = {
    "canonical": lambda r: r["canon_in_class"],
    "noncanonical": lambda r: not r["canon_in_class"],
    "both": lambda r: True,
}

sweep = {}
print(f"\n[644h] annotated calibration (n={len(ann_max)} junctions, "
      f"cap {CAL_READS_PER_J} reads): q_max quantiles {quantiles(ann_max)}",
      flush=True)
for fname, ffn in frames.items():
    for tname, tfn in tracks.items():
        sel = [r for r in rows if ffn(r) and tfn(r)
               and r["class"] != "annotated"]
        g = [r for r in sel if r["is_gold"]]
        k = [r for r in sel if not r["is_gold"]]
        tab = []
        for thr in THRESH_GRID:
            gk = sum(1 for r in g if r["q_max"] >= thr)
            kk = sum(1 for r in k if r["q_max"] >= thr)
            ann_pass = (sum(1 for v in ann_max if v >= thr) / len(ann_max)
                        if ann_max else None)
            tab.append({"q_thr": thr, "gold_kept": gk, "junk_kept": kk,
                        "ann_pass": round(ann_pass, 3) if ann_pass is not None else None})
        # rank AUC (gold vs junk, q_max)
        auc = None
        if g and k:
            gv = sorted(r["q_max"] for r in g)
            kv = sorted(r["q_max"] for r in k)
            import bisect
            wins = ties = 0
            for x in gv:
                lo = bisect.bisect_left(kv, x)
                hi = bisect.bisect_right(kv, x)
                wins += lo
                ties += hi - lo
            auc = round((wins + 0.5 * ties) / (len(gv) * len(kv)), 3)
        sweep[f"{fname}/{tname}"] = {"gold_n": len(g), "junk_n": len(k),
                                     "auc": auc, "table": tab}
        print(f"[644h] frame={fname} track={tname}: gold={len(g)} "
              f"junk={len(k)} AUC={auc}", flush=True)
        for row in tab:
            print(f"        q>={row['q_thr']:>3}: gold {row['gold_kept']:>3} "
                  f"junk {row['junk_kept']:>5} ann_pass {row['ann_pass']}",
                  flush=True)

gold_detail = sorted([r for r in rows if r["is_gold"]],
                     key=lambda r: -r["q_max"])
out = {
    "params": {"MIN_ANCHOR": MIN_ANCHOR, "CAP": CAP, "ERR_BITS": ERR_BITS,
               "CAL_READS_PER_J": CAL_READS_PER_J, "HP_DISCOUNT": HP_DISCOUNT,
               "limit_reads": LIM,
               "score": "q_read = max(0, I_eff(short-side overhang, cap 60bp)"
                        " - 2.0*errors_in_window); junction q_max = max over"
                        " supporting reads"},
    "annotated_calibration": {"n_junctions": len(ann_max),
                              "q_max_quantiles": quantiles(ann_max),
                              "q_max_values": [round(v, 2) for v in ann_max]},
    "sweep": sweep,
    "gold_detail": gold_detail,
    "junk_junctions": [
        {k: r[k] for k in ("junction", "canon_in_class", "n_reads", "q_max",
                           "q_2nd", "clean_prefix_max", "surv_wins")}
        for r in rows if not r["is_gold"] and r["class"] != "annotated"],
}
with open(args.out, "w") as fh:
    json.dump(out, fh, indent=1)
print(f"[644h] wrote {args.out}", flush=True)
