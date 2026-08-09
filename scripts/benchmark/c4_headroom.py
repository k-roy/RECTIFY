#!/usr/bin/env python3
"""C4 window-selection headroom probe — is the incumbent below ceiling on
*identifiable* paralog reads, and is the gap POOLING-addressable?  (DECODER-FREE;
run this FIRST.)

C4's premise: paralog / multi-copy loci (SMN1/SMN2-style) defeat per-read window
selection — a read maps to the WRONG near-identical copy — and a POA-pooled member
(cluster reads at the locus -> consensus -> align once -> project back) recovers it.

Per the C1/C2/C3 gate discipline, before building any pooling decoder we must
establish the incumbent is BELOW CEILING on *recoverable* reads AND the gap is
POOLING-addressable, NOT the paralog-zero-evidence trap (a fragment whose lone
distinguishing base is destroyed/inverted by noise is informationally identical to
the wrong copy -> NO method recovers it per-read).  This probe measures exactly
that, DECODER-FREE (no pooling member built), so there is zero construction-tuned-
win risk.

The decoder-free gate (mirrors c3_headroom's ceiling - incumbent):

    For each WEAK fragment (covers EXACTLY ONE distinguishing SNP), classify on
    two axes scored vs TRUTH:
      * evidence  : does the read still carry THIS copy's distinguishing base?
            intact      = base == this copy's SNP allele       (IDENTIFIABLE)
            looks_other = base == the OTHER copy's SNP allele   (wrong-evidence)
            third       = base == a non-paralog base            (zero-evidence)
      * placement : did minimap2 map it to its TRUE origin contig?
            correct / WRONG / unplaced  (+ mapq0 split)

    ceiling   = the IDENTIFIABLE universe = intact-base reads (truth recoverable
                per-read because the read still carries decisive evidence)
    incumbent = freq(minimap2 maps an intact-base read to the correct copy)
    HEADROOM  = freq(intact base BUT mis-placed) = the genuine C4-addressable gap

A pooling member can only earn its keep on the HEADROOM cell (identifiable-but-
misplaced).  A below-ceiling number made ONLY of looks_other/third reads is the
pre-committed paralog-zero-evidence NULL — exclude it; no method (incl. C4) recovers
a read whose own distinguishing base is gone.

The probe ALSO audits the smoke (F) pooling-recovery claim two ways:
  * truth-circularity: (F) pools by t.true_transcript (the TRUE copy label), i.e.
    it presupposes the clustering C4 must SOLVE.  We print minimap2's OWN per-copy
    placement distribution to show the incumbent already recovers the pools.
  * minority-collapse: equal reps/copy means "majority base" cannot collapse a
    minority — the membership-preserving fence is untested by (F).

VERDICT (pre-committed):
  HEADROOM ~ 0  => incumbent at ceiling on identifiable reads; the gap is the
    zero-evidence null => C4-as-window-selection REFUTED on this stratum.  Document
    + stop; do NOT build a pooling member to chase a phantom gap.
  HEADROOM > 0  => a genuine identifiable-but-misplaced population exists; THEN
    test whether a truth-FREE clustering + pooling recovers it WITHOUT collapsing
    minorities (the next-cycle ablation).

Why minimap2 alone (not the full panel): C4 is a WHICH-CONTIG (locus-selection)
question, so the single-contig DP arms (align_exon_block_global, the C1/C3 tool)
do NOT apply.  minimap2 is a fair AND conservative proxy: the decisive-SNP logic
is aligner-invariant (any seed-chain aligner wins the right copy by exactly 1
mismatch on an intact read), minimap2/gapmm2/uLTRA share the per-read blind spot
(OVERVIEW), and the shipped scorer/smoke locus metric is itself minimap2-based.

Usage:
  python scripts/benchmark/c4_headroom.py --reps 100 --noise 0.05
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import random
import subprocess
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.scorer import load_genome  # noqa: E402
from scripts.benchmark.sim import controlled  # noqa: E402
import pysam  # noqa: E402


def load_fastq(path):
    seqs = {}
    with pysam.FastxFile(path) as fq:
        for e in fq:
            seqs[e.name] = e.sequence
    return seqs


def run_minimap2(ref_fa, reads_fq, out_bam):
    p = subprocess.run(
        ["minimap2", "-ax", "splice", "-uf", "--eqx", "-k", "14", "-t", "2",
         ref_fa, reads_fq], capture_output=True)
    if p.returncode:
        raise RuntimeError("minimap2 failed: " + p.stderr.decode()[:400])
    s = subprocess.run(["samtools", "sort", "-o", out_bam], input=p.stdout,
                       capture_output=True)
    if s.returncode:
        raise RuntimeError("samtools sort failed: " + s.stderr.decode()[:300])
    subprocess.run(["samtools", "index", out_bam], check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/c4_hr")
    ap.add_argument("--reps", type=int, default=100,
                    help="reads PER copy PER kind (per-locus DEPTH; a C4-pooling "
                         "feature, not diversity — scale --families for diversity)")
    ap.add_argument("--families", type=int, default=3)
    ap.add_argument("--noise", type=float, default=0.05,
                    help="per-base substitution rate (realistic ONT ~0.05-0.10)")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"[c4hr] generating paralog stratum reps={args.reps} families={args.families} "
          f"noise={args.noise} ...", file=sys.stderr)
    rng = random.Random(args.seed)
    refs, reads, truth = controlled.gen_paralog_stratum(
        args.reps, rng, noise=args.noise, n_families=args.families)
    ref_fa = os.path.join(args.out, "ref.fa")
    reads_fq = os.path.join(args.out, "reads.fastq")
    with open(ref_fa, "w") as fh:
        for n, s in refs.items():
            fh.write(f">{n}\n{s}\n")
    with open(reads_fq, "w") as fh:
        for rid, s in reads:
            fh.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")
    truth_map = {t.read_id: t for t in truth}
    genome = load_genome(ref_fa)
    read_seqs = load_fastq(reads_fq)

    bam = os.path.join(args.out, "mm2.bam")
    run_minimap2(ref_fa, reads_fq, bam)
    placed = {}
    with pysam.AlignmentFile(bam, "rb") as b:
        for r in b:
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if r.query_name in truth_map:
                placed[r.query_name] = r

    weak = {rid: t for rid, t in truth_map.items() if "_weak_" in rid}
    span = {rid: t for rid, t in truth_map.items() if "_span_" in rid}

    # ---- per-read 2x2 classification ----
    cells = defaultdict(int)
    pool_truebase = defaultdict(list)          # (locus,copy) -> base@SNP across reads
    pool_mapped = defaultdict(Counter)         # (locus,copy) -> minimap2 mapped-contig counts
    for rid, t in weak.items():
        v = t.variants[0]
        idx = v.pos - t.genome_start
        seq = read_seqs[rid]
        base = seq[idx] if 0 <= idx < len(seq) else "N"
        if base == v.alt_allele:
            evid = "intact"        # this copy's distinguishing base survived -> IDENTIFIABLE
        elif base == v.ref_allele:
            evid = "looks_other"   # flipped to the OTHER copy's base -> wrong-evidence
        else:
            evid = "third"         # flipped to a non-paralog base -> zero-evidence
        r = placed.get(rid)
        if r is None:
            mp, mq0 = "unplaced", False
        else:
            mp = "correct" if r.reference_name == t.chrom else "WRONG"
            mq0 = (r.mapping_quality == 0)
        cells[(evid, mp)] += 1
        if r is not None and mp == "WRONG":
            cells[(evid, "WRONG", "mapq0" if mq0 else "mapq>0")] += 1
        key = (t.true_locus, t.true_transcript)
        pool_truebase[key].append(base)
        if r is not None:
            pool_mapped[key][r.reference_name] += 1

    # span (ceiling control): locus accuracy
    span_correct = sum(1 for rid, t in span.items()
                       if placed.get(rid) is not None
                       and placed[rid].reference_name == t.chrom)
    span_placed = sum(1 for rid in span if rid in placed)
    span_acc = span_correct / span_placed if span_placed else 0.0

    # ---- report ----
    print("\n================ C4 WINDOW-SELECTION HEADROOM (decoder-free) ================")
    print(f"stratum=PARALOG  reps={args.reps}/copy/kind  families={args.families}  "
          f"noise={args.noise}  aligner=minimap2 -ax splice -uf")
    print(f"\nSPAN control (covers ALL spread SNPs, redundant evidence): "
          f"locus_acc={span_acc:.4f}  ({span_correct}/{span_placed} placed)  "
          f"-> at-ceiling check (the metric is not trivially failing)")

    print(f"\nWEAK fragments (cover EXACTLY ONE distinguishing SNP), n={len(weak)}")
    print(f"  {'evidence':12s} {'correct':>8s} {'WRONG':>7s} {'unplaced':>9s} {'n':>6s}")
    universe = {}
    for evid in ("intact", "looks_other", "third"):
        c, w, u = (cells[(evid, "correct")], cells[(evid, "WRONG")],
                   cells[(evid, "unplaced")])
        universe[evid] = c + w + u
        print(f"  {evid:12s} {c:8d} {w:7d} {u:9d} {c + w + u:6d}")

    intact_n = universe["intact"]
    genuine = cells[("intact", "WRONG")]
    genuine_mq0 = cells[("intact", "WRONG", "mapq0")]
    genuine_mqp = cells[("intact", "WRONG", "mapq>0")]
    intact_correct = cells[("intact", "correct")]
    intact_placed = intact_correct + genuine
    incumbent = intact_correct / intact_placed if intact_placed else 0.0
    headroom = genuine / intact_placed if intact_placed else 0.0

    print(f"\nceiling  (IDENTIFIABLE universe = intact-base reads)        = {intact_n}")
    print(f"incumbent(minimap2 correct | intact base, placed)           = {incumbent:.4f} "
          f"({intact_correct}/{intact_placed})")
    print(f"HEADROOM (intact base BUT mis-placed = genuine C4 gap)      = {headroom:.4f} "
          f"({genuine}; mapq0={genuine_mq0}, mapq>0={genuine_mqp})")
    zero_ev = universe["looks_other"] + universe["third"]
    print(f"zero-evidence null (looks_other + third corrupted-base)     = {zero_ev} "
          f"({universe['looks_other']} wrong-evidence + {universe['third']} ambiguous)")
    print(f"intact unplaced (alignment-sensitivity / C5, not C4)        = "
          f"{cells[('intact', 'unplaced')]}")

    # ---- smoke (F) pooling-claim audit ----
    print("\n---- smoke (F) pooling-recovery AUDIT (truth-circular + redundant) ----")
    recovered = redundant = 0
    for key, bases in pool_truebase.items():
        rid0 = next(rid for rid, t in weak.items()
                    if (t.true_locus, t.true_transcript) == key)
        alt = weak[rid0].variants[0].alt_allele
        maj = Counter(bases).most_common(1)[0][0]
        recovered += (maj == alt)
        mc = pool_mapped[key]
        # would minimap2's OWN placement-majority recover this pool, truth-free?
        mm_maj = mc.most_common(1)[0][0] if mc else None
        redundant += (mm_maj == key[1])
        print(f"  {key[1]:10s}: pool-majority-base={maj} true={alt}  "
              f"minimap2 placement={dict(mc)}")
    n_pools = len(pool_truebase)
    print(f"(F) truth-grouped pooling recovery = {recovered}/{n_pools}  "
          f"[but pools are formed FROM truth -> circular]")
    print(f"minimap2 OWN per-read placement-majority already recovers = "
          f"{redundant}/{n_pools} pools  [pooling adds nothing the incumbent didn't deliver]")

    print("\n---- VERDICT (pre-committed) ----")
    if headroom < 0.01:
        print(f"  HEADROOM={headroom:.4f} ~ 0  => minimap2 AT CEILING on identifiable")
        print("  (intact-base) paralog reads. The entire below-ceiling weak gap falls on")
        print("  CORRUPTED-base reads (looks_other = wrong-evidence; third = zero-evidence)")
        print("  -> the pre-committed paralog-zero-evidence NULL. The smoke (F) pooling proof")
        print("  is truth-circular (pools = truth copy label) AND redundant (minimap2's own")
        print("  placement already recovers the pools).")
        print("  => C4-as-window-selection REFUTED on this stratum. Document + stop; do NOT")
        print("     build a pooling member. Named residue: a >=3-copy / mapq0-with-evidence")
        print("     construct (true SMN1/SMN2) is UNEXERCISED here (2 copies + 1 covered SNP")
        print("     => an intact base is decisive by exactly 1 mismatch => mapq>0 always).")
        print("     Defer to a MEASURED real-data trigger; do not synthesize a gap.")
    else:
        print(f"  HEADROOM={headroom:.4f} > 0  => a genuine identifiable-but-misplaced")
        print("  population EXISTS. NEXT: test whether a truth-FREE clustering + pooling")
        print("  recovers these reads WITHOUT collapsing minority copies (abundance-skew),")
        print("  the membership-preserving fence — the next-cycle ablation.")


if __name__ == "__main__":
    main()
