#!/usr/bin/env python3
"""Read-side tri-nucleotide repeat-expansion artifact profiler (organism-agnostic).

Reproduces the yeast wt_by4742_rep1 DRS profile ("Repeat expansion artifacts by
k-mer length") on arbitrary DRS BAMs, to test whether the k=3 / (AAG·CTT)-family
dominance is species-agnostic (Kevin's hunch, 2026-05-26).

Method (recovered from the yeast analysis):
  - For each PRIMARY read, scan three read-side channels for runs >= MIN_LEN bp:
      left soft-clip  (cigar[0]  == S)
      right soft-clip (cigar[-1] == S)
      insertions      (op == I)
  - For each segment, find the minimal-period dominant k-mer (k in 1..5):
    aggregate k-mers by canonical rotation family; the smallest k whose dominant
    canonical family covers >= MIN_FRAC of windows (and the segment spans >= 3 full
    repeats) is the assigned period. This puts poly-A under k=1, (GAA)n under k=3.
  - Report in RNA space: a reverse-strand read's stored sequence is RC of the RNA,
    so the RNA-canonical motif is canonical(rc(dna_motif)) when is_reverse. This
    collapses the +strand AAG and -strand CTT into one RNA motif (same slippage
    event in opposite orientations).
  - EPM = events / total_primary_reads * 1e6.

Motif-agnostic by construction: detects repeat STRUCTURE at any period, never gates
on a specific motif (see no-motifs-unbiased-discovery).

Run on Sherlock (BAMs live there):
    conda activate rectify
    python dev/repeat_expansion_profile.py --bam <sorted.bam> --label <name> [--min-frac 0.75]
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict

import pysam

_S, _I = 4, 1
_COMP = str.maketrans("ACGTN", "TGCAN")


def rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def canonical(kmer: str) -> str:
    """Lexicographically smallest rotation."""
    return min(kmer[i:] + kmer[:i] for i in range(len(kmer)))


def dominant_kmer(seq: str, k: int, min_frac: float):
    """Dominant canonical k-mer family if it covers >= min_frac of windows and the
    segment spans >= 3 full repeats. Returns (canonical_dna_motif, frac) or None."""
    n = len(seq)
    if n < k * 3:
        return None
    n_kmers = n - k + 1
    counts = Counter(seq[i:i + k] for i in range(n_kmers))
    by_canon: Counter = Counter()
    for kmer, c in counts.items():
        if "N" in kmer:
            continue
        by_canon[canonical(kmer)] += c
    if not by_canon:
        return None
    best, cnt = by_canon.most_common(1)[0]
    frac = cnt / n_kmers
    if frac < min_frac:
        return None
    return best, frac


def classify_segment(seq: str, min_frac: float, max_single_base: float = 1.0):
    """Minimal-period dominant repeat. Returns (k, canonical_dna_motif) or None.

    ``max_single_base`` < 1.0 excludes homopolymer-dominated segments from the
    k>=2 repeat classes: a poly-A/poly-T DRS tail with sparse basecall errors can
    dip below the k=1 purity floor and fall into the k=3 ``AAA``/``AAG`` family,
    inflating the apparent (AAG)n count. A genuine (AAG)n is multi-base (no single
    base dominates), so requiring max single-base fraction < ``max_single_base``
    (e.g. 0.75) keeps real expansions and drops poly-A leakage. k=1 (true
    homopolymer / poly-A) is unaffected — it is reported as its own class.
    """
    seq = seq.upper()
    for k in (1, 2, 3, 4, 5):
        r = dominant_kmer(seq, k, min_frac)
        if r:
            if k >= 2 and max_single_base < 1.0 and seq:
                if max(seq.count(b) for b in "ACGT") / len(seq) >= max_single_base:
                    return None  # homopolymer-dominated (poly-A/T leak), not a repeat expansion
            return k, r[0]
    return None


def profile_bam(bam_path: str, min_len: int, min_frac: float, max_reads,
                max_single_base: float = 1.0):
    total_primary = 0
    by_k: Counter = Counter()
    by_k_motif: Counter = Counter()       # (k, rna_motif) -> events
    by_strand_k3: Counter = Counter()     # (strand, rna_motif) -> events for k=3
    by_channel: Counter = Counter()       # channel -> events
    by_k_channel: Counter = Counter()     # (k, channel) -> events

    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            total_primary += 1
            if max_reads and total_primary > max_reads:
                total_primary -= 1
                break
            cig = read.cigartuples
            seq = read.query_sequence
            if not cig or not seq:
                continue
            rev = read.is_reverse

            segments = []  # (channel_bam, seq)
            if cig[0][0] == _S and cig[0][1] >= min_len:
                segments.append(("clip_left", seq[:cig[0][1]]))
            if cig[-1][0] == _S and cig[-1][1] >= min_len:
                segments.append(("clip_right", seq[-cig[-1][1]:]))
            qpos = cig[0][1] if cig[0][0] == _S else 0
            for op, ln in cig:
                if op == _I:
                    if ln >= min_len:
                        segments.append(("insertion", seq[qpos:qpos + ln]))
                    qpos += ln
                elif op in (0, 7, 8):
                    qpos += ln
                # D/N/S(handled)/H consume no further query here (leading S already applied)

            for chan_bam, segseq in segments:
                res = classify_segment(segseq, min_frac, max_single_base)
                if not res:
                    continue
                k, dna_motif = res
                rna_motif = canonical(rc(dna_motif)) if rev else dna_motif
                strand = "-" if rev else "+"
                # RNA 5'/3' channel: forward read left=5', reverse read left=3'
                if chan_bam == "insertion":
                    chan = "insertion"
                elif (chan_bam == "clip_left") != rev:
                    chan = "softclip_5p"
                else:
                    chan = "softclip_3p"
                by_k[k] += 1
                by_k_motif[(k, rna_motif)] += 1
                by_channel[chan] += 1
                by_k_channel[(k, chan)] += 1
                if k == 3:
                    by_strand_k3[(strand, rna_motif)] += 1

    return {
        "total_primary": total_primary,
        "by_k": dict(by_k),
        "by_k_motif": {f"{k}:{m}": c for (k, m), c in by_k_motif.items()},
        "by_strand_k3": {f"{s}:{m}": c for (s, m), c in by_strand_k3.items()},
        "by_channel": dict(by_channel),
        "by_k_channel": {f"{k}:{c}": v for (k, c), v in by_k_channel.items()},
    }


def epm(n, total):
    return n / total * 1e6 if total else 0.0


def report(label, res, min_frac):
    t = res["total_primary"]
    print(f"\n========== {label}  (min_frac={min_frac}, n={t:,} primary reads) ==========")
    print("Total events by repeat period k (EPM):")
    tot_k = sum(res["by_k"].values()) or 1
    for k in (1, 2, 3, 4, 5):
        n = res["by_k"].get(k, 0)
        print(f"  k={k}: {n:>8,}  EPM={epm(n, t):>9.1f}  ({100*n/tot_k:.1f}% of events)")
    k3 = res["by_k"].get(3, 0)
    others = max(1, tot_k - k3)
    print(f"  --> k=3 dominance: {k3} vs {tot_k-k3} others = {k3/others:.1f}x")

    print("Top k=3 RNA motifs (EPM):")
    k3m = sorted(((m, c) for m, c in res["by_k_motif"].items() if m.startswith("3:")),
                 key=lambda x: -x[1])[:8]
    for m, c in k3m:
        print(f"  {m.split(':',1)[1]:>6}: {c:>7,}  EPM={epm(c, t):>8.1f}")

    print("k=3 strand split (top motifs):")
    ss = sorted(res["by_strand_k3"].items(), key=lambda x: -x[1])[:8]
    for sm, c in ss:
        print(f"  {sm:>10}: {c:>7,}")

    print("Channel breakdown (all k, EPM):")
    for chan, n in sorted(res["by_channel"].items(), key=lambda x: -x[1]):
        print(f"  {chan:>12}: {n:>8,}  EPM={epm(n, t):>9.1f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--min-len", type=int, default=30)
    ap.add_argument("--min-frac", type=float, default=0.75)
    ap.add_argument("--max-single-base", type=float, default=1.0,
                    help="exclude homopolymer-dominated (poly-A/T leak) segments from k>=2 "
                         "classes when a single base covers >= this fraction (e.g. 0.75)")
    ap.add_argument("--max-reads", type=int, default=0)
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args()

    print(f"Scanning {args.bam} ...", file=sys.stderr)
    res = profile_bam(args.bam, args.min_len, args.min_frac, args.max_reads or None,
                      args.max_single_base)
    report(args.label, res, args.min_frac)
    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump({"label": args.label, "min_frac": args.min_frac,
                       "min_len": args.min_len, **res}, fh, indent=2)
        print(f"\nWrote {args.json_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
