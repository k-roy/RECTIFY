#!/usr/bin/env python3
"""C1 ablation — does a calibrated in-run gap discount improve indel PLACEMENT?

The gate-sanctioned proof for the de-novo aligner's first facet (C1). Three arms,
IDENTICAL except the deletion gap-OPEN cost, ALL passing chrom_ref so the legacy
HP-awareness (homo_mask + homo_mismatch=-2) is active in every arm (the
matched-baseline requirement — without chrom_ref the baseline has NO HP-awareness
and any "win" is a conflated artifact):

  flat    : penalty_table=None              (legacy homo_mismatch only)
  B0      : a CONSTANT in-run gap discount  (no run-length dependence)
  law     : the per-(hp,base) log-odds delta (rate_mean baseline-anchored)

Scored vs TRUTH (never the internal DP score) on the TEST partition:
  * position-exact indel concordance on HP_HARD-noisy   (Claim A: B0/law > flat)
  * false_indel_rate on CLEAN HP/STR runs routed through the arm (must stay ~0)
  * boundary_sub concordance must NOT regress

WHAT THIS CAN AND CANNOT PROVE (see dev/C1_DESIGN.md):
  * Claim A ("an in-run discount improves placement") IS testable here.
  * Claim B ("the table's per-length SHAPE is correct") is NOT — the HP_HARD
    error model (K_DIST) is flat in run length, and length-scaling it from the
    same rate_mean table would make the law win by construction. So law ~ B0 here
    is the EXPECTED, CORRECT result ("length-shape inconclusive on flat-in-L sim
    -> deferred to real SIRV/RNA004"), reported as neither a win nor a failure.

Usage: python scripts/benchmark/c1_ablation.py --out /tmp/c1_abl --reps 400 \
         --penalty-table rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.scorer import (  # noqa: E402
    score_bam, load_genome, cigar_records_to_bam,
)
from rectify.core.benchmark.truth_schema import read_truth_table, SplitTag  # noqa: E402
from rectify.core.align.local_aligner import align_exon_block_global  # noqa: E402
from rectify.core.splice.hp_penalty import HpPenaltyTable  # noqa: E402
from scripts.benchmark.sim import controlled  # noqa: E402
import pysam  # noqa: E402


class ConstantDiscount:
    """B=0 control: a CONSTANT in-run gap-open discount, no run-length shape.
    Same average magnitude as the law arm, so a law-vs-B0 difference isolates the
    run-length SHAPE (which flat-in-L sim cannot validate)."""
    def __init__(self, d_del: float, d_ins: float = 0.0):
        self._d_del, self._d_ins = d_del, d_ins

    def del_open_delta(self, hp, base, lam: float = 1.0):
        return lam * self._d_del

    def ins_open_delta(self, hp, base, lam: float = 1.0):
        return lam * self._d_ins


def load_fastq(path):
    seqs = {}
    with pysam.FastxFile(path) as fq:
        for e in fq:
            seqs[e.name] = e.sequence
    return seqs


def run_arm(read_seqs, truth_subset, genome, ref_fa, out_bam, penalty_table, lam):
    """Matched arm: align each read to its single-contig ref WITH chrom_ref (so
    homo_mask/homo_mismatch are active), optionally with a length-law penalty."""
    records = []
    for rid, t in truth_subset.items():
        seq = read_seqs.get(rid)
        if seq is None:
            continue
        ref = genome[t.chrom]
        cig = align_exon_block_global(seq, ref, chrom_ref=ref, ref_offset=0,
                                      penalty_table=penalty_table, lam=lam)
        records.append((rid, t.chrom, 0, cig, seq))
    cigar_records_to_bam(records, ref_fa, out_bam)
    return score_bam(out_bam, truth_subset, genome, aligner_name="arm")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/c1_abl")
    ap.add_argument("--reps", type=int, default=400)
    ap.add_argument("--penalty-table", required=True)
    ap.add_argument("--str-table", default=None)
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--min-count", type=int, default=100)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"[c1] generating corpus reps={args.reps} ...", file=sys.stderr)
    info = controlled.generate_corpus(args.out, reps=args.reps, seed=7)
    genome = load_genome(info["ref_fa"])
    ref_fa = info["ref_fa"]
    read_seqs = load_fastq(info["reads_fastq"])
    truth = read_truth_table(info["truth_tsv"])

    law_table = HpPenaltyTable.from_tsv(args.penalty_table, args.str_table,
                                        min_count=args.min_count)

    # TEST partition only (held-out region split)
    test = [t for t in truth if t.split == SplitTag.TEST]

    def subset(pred):
        return {t.read_id: t for t in test if pred(t)}

    noisy = subset(lambda t: t.stratum == "HP_HARD" and "_noisy_" in t.read_id)
    bsub = subset(lambda t: t.stratum == "HP_HARD" and "_boundary_sub_" in t.read_id)
    clean = subset(lambda t: t.stratum in ("HP", "STR") and not t.indels)
    print(f"[c1] TEST reads: noisy={len(noisy)} boundary_sub={len(bsub)} clean={len(clean)}",
          file=sys.stderr)

    # B=0 constant = mean of the law deltas over the run (hp, base) actually present
    # in the HP_HARD-noisy set (equal average magnitude, no shape).
    deltas = []
    for t in noisy.values():
        if t.indels and t.indels[0].run_copies:
            ind = t.indels[0]
            deltas.append(law_table.del_open_delta(ind.run_copies,
                                                   ind.run_unit or "A", args.lam))
    const = sum(deltas) / len(deltas) if deltas else 0.0
    print(f"[c1] B=0 constant discount = {const:.4f} (mean law delta, lam={args.lam})",
          file=sys.stderr)

    arms = {
        "flat": (None, 1.0),
        "B0":   (ConstantDiscount(const), 1.0),  # const already includes lam
        "law":  (law_table, args.lam),
    }

    def run_all(sub, tag):
        out = {}
        for name, (pt, lam) in arms.items():
            sc = run_arm(read_seqs, sub, genome, ref_fa,
                         os.path.join(args.out, f"{tag}_{name}.bam"), pt, lam)
            out[name] = sc
        return out

    noisy_s = run_all(noisy, "noisy") if noisy else {}
    bsub_s = run_all(bsub, "bsub") if bsub else {}
    clean_s = run_all(clean, "clean") if clean else {}

    def conc(d, k):
        return d[k].indel.position_exact_concordance if d else float("nan")

    def fir(d, k):
        return d[k].indel.false_indel_rate if d else float("nan")

    print("\n================ C1 ABLATION (TEST split, reps=%d) ================" % args.reps)
    print(f"{'metric':28s} {'flat':>10s} {'B0':>10s} {'law':>10s}")
    print(f"{'HP_HARD-noisy concordance':28s} {conc(noisy_s,'flat'):10.4f} "
          f"{conc(noisy_s,'B0'):10.4f} {conc(noisy_s,'law'):10.4f}")
    print(f"{'boundary_sub concordance':28s} {conc(bsub_s,'flat'):10.4f} "
          f"{conc(bsub_s,'B0'):10.4f} {conc(bsub_s,'law'):10.4f}")
    print(f"{'clean false_indel_rate':28s} {fir(clean_s,'flat'):10.4f} "
          f"{fir(clean_s,'B0'):10.4f} {fir(clean_s,'law'):10.4f}")
    if clean_s:
        print(f"  (clean reads scored per arm: {clean_s['flat'].indel.clean_reads})")

    # ---- go/no-go (Claim A) ----
    print("\n---- VERDICT (Claim A: in-run discount improves placement) ----")
    if noisy_s and clean_s and bsub_s:
        cf, cb, cl = conc(noisy_s,'flat'), conc(noisy_s,'B0'), conc(noisy_s,'law')
        placement = (cl > cf + 1e-9) or (cb > cf + 1e-9)
        no_halluc = max(fir(clean_s,'B0'), fir(clean_s,'law')) <= fir(clean_s,'flat') + 1e-9 \
            and fir(clean_s,'law') < 0.02
        no_regress = conc(bsub_s,'law') >= conc(bsub_s,'flat') - 1e-9
        print(f"  placement improves (B0 or law > flat): {placement} "
              f"(flat={cf:.4f} B0={cb:.4f} law={cl:.4f})")
        print(f"  no hallucinated indels on clean (law fir<=flat, <0.02): {no_halluc} "
              f"(flat={fir(clean_s,'flat'):.4f} law={fir(clean_s,'law'):.4f})")
        print(f"  boundary_sub not regressed (law>=flat): {no_regress} "
              f"(flat={conc(bsub_s,'flat'):.4f} law={conc(bsub_s,'law'):.4f})")
        print(f"  length-SHAPE (law vs B0): law={cl:.4f} B0={cb:.4f} -> "
              f"{'INCONCLUSIVE (expected on flat-in-L sim; defer to real SIRV/RNA004)' if abs(cl-cb)<0.01 else 'DIFFERS (investigate; flat-in-L should NOT separate shape)'}")
        verdict = "PASS (Claim A)" if (placement and no_halluc and no_regress) else "INCONCLUSIVE/FAIL"
        print(f"  ==> {verdict}")
    else:
        print("  insufficient TEST reads in one or more subsets — raise --reps")


if __name__ == "__main__":
    main()
