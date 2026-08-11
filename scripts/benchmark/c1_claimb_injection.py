#!/usr/bin/env python3
"""C1 Claim B — length-stratified placement ablation (the SHAPE axis).

Distinct from Claim A (c1_ablation.py):

  Claim A  = "an in-run gap discount improves indel PLACEMENT"   (PASS, aggregate)
  Claim B  = "the table's per-length SHAPE is correct/transfers" (THIS TRACK)

--------------------------------------------------------------------------------
HONEST SCOPE  (advisor-vetted 2026-06-30)
--------------------------------------------------------------------------------
GENUINE Claim B needs a SIRV-MEASURED per-HP-run-length deletion curve as the
injection rate. That curve was NEVER recorded as numbers (the SPEC SIRV result is
error-STRUCTURE stats — indel_run>=2, overdispersion, autocorr — NOT a
del-rate-by-runlen curve), and reconstructing it means processing the SIRV BAMs
(the deferred multi-night track). So this script does NOT confirm Claim B.

What it DOES do, going one honest step past the flat AGGREGATE of Claim A: reuse
the VETTED HP_HARD `boundary_sub` construction (`scripts/benchmark/sim/controlled.
gen_hp_hard_stratum` — the exact generator that produced Claim A's real
0.00->0.55/0.78 headline) and break the flat/B0/law placement concordance down BY
RUN LENGTH. That per-length breakdown IS the length-SHAPE axis: the law's gap-open
delta EXCEEDS B0 (=the mean delta) at LONG runs and is BELOW it at SHORT runs, so
if the length-shape is doing real work it should show as law>B0 in the LONG-L bin
and law<B0 in the SHORT-L bin, summing to the ~tie/slight-B0-win Claim A reported.

NON-CIRCULARITY: the DP cost is the Scer `penalty_scores.tsv` (rate_mean); NOTHING
is calibrated on this sim (the sim only supplies reads), so scoring all reads is
leak-free. The injection LENGTH is `K_DIST` (flat in L) — so this is the FLAT-in-L
substrate; a real length-correlated (rising) injection is the SIRV-measured piece
that is still deferred. A `--reweight-long` aggregate shows what a rising injection
would do to the headline WITHOUT changing any per-cell number (honest reweighting,
not a new error model).

Report the per-length law-vs-B0 as a MECHANISM finding, never as Claim B. Do NOT
crank lambda to force law>B0 (dev/C1_DESIGN.md lines 75-77).

Usage:
  PYTHONPATH=. python scripts/benchmark/c1_claimb_injection.py \
      --out /tmp/c1_claimb --reps 400 \
      --penalty-table rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv

Author: Kevin R. Roy  (agent C1-ClaimB)
"""
from __future__ import annotations

import argparse
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from rectify.core.benchmark.scorer import (  # noqa: E402
    score_bam, load_genome, cigar_records_to_bam,
)
from rectify.core.align.local_aligner import align_exon_block_global  # noqa: E402
from rectify.core.splice.hp_penalty import HpPenaltyTable  # noqa: E402
from scripts.benchmark.sim import controlled  # noqa: E402


# ---- B=0 control arm (identical to c1_ablation.py) --------------------------
class ConstantDiscount:
    """B=0: a CONSTANT in-run gap-open discount = the MEAN of the law deltas over
    the realized (hp,base) cells. A law-vs-B0 difference then isolates the
    run-length SHAPE — the quantity Claim B is about."""
    def __init__(self, d_del: float, d_ins: float = 0.0):
        self._d_del, self._d_ins = d_del, d_ins

    def del_open_delta(self, hp, base, lam: float = 1.0):
        return lam * self._d_del

    def ins_open_delta(self, hp, base, lam: float = 1.0):
        return lam * self._d_ins


def run_arm(read_seqs, refs, truth_subset, genome, ref_fa, out_bam, penalty_table, lam):
    """Matched arm (mirror c1_ablation.run_arm): align each read to its single-contig
    ref WITH chrom_ref (so homo_mask/homo_mismatch are active in EVERY arm)."""
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
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="/tmp/c1_claimb")
    ap.add_argument("--reps", type=int, default=400)
    ap.add_argument("--penalty-table", required=True)
    ap.add_argument("--str-table", default=None)
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--min-count", type=int, default=100)
    ap.add_argument("--stratum", choices=["boundary_sub", "noisy", "both"],
                    default="boundary_sub",
                    help="which HP_HARD sub-case to score. boundary_sub = pure "
                         "THRESHOLD (flat=0 everywhere); noisy = GRADED (flat near "
                         "ceiling) — the case where a per-length law-vs-B0 signal "
                         "could actually surface (advisor 2026-06-30).")
    ap.add_argument("--reweight-long", action="store_true",
                    help="also print an aggregate that up-weights LONG runs (what a "
                         "rising/SIRV-shaped injection would do to the headline)")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = random.Random(args.seed)

    print(f"[claimB] generating HP_HARD corpus reps={args.reps} (vetted gen_hp_hard_stratum)",
          file=sys.stderr)
    refs, reads_list, truth_list = controlled.gen_hp_hard_stratum(args.reps, rng)

    # write one contig per HP_HARD locus
    ref_fa = os.path.join(args.out, "ref.fa")
    with open(ref_fa, "w") as fh:
        for chrom, seq in refs.items():
            fh.write(f">{chrom}\n{seq}\n")
    import pysam
    pysam.faidx(ref_fa)
    genome = load_genome(ref_fa)

    read_seqs = dict(reads_list)
    truth_by_id = {t.read_id: t for t in truth_list}

    def want(rid):
        bs = "_boundary_sub_" in rid
        ny = "_noisy_" in rid
        if args.stratum == "boundary_sub":
            return bs
        if args.stratum == "noisy":
            return ny
        return bs or ny
    ids = [rid for rid in read_seqs if want(rid)]

    # stratify by run length L (= truth indel.run_copies)
    by_len = defaultdict(list)
    for rid in ids:
        t = truth_by_id[rid]
        L = t.indels[0].run_copies if t.indels else 0
        by_len[L].append(rid)

    law_table = HpPenaltyTable.from_tsv(args.penalty_table, args.str_table,
                                        min_count=args.min_count)

    # B=0 constant = mean law delta over the realized boundary_sub cells
    deltas = []
    for rid in ids:
        t = truth_by_id[rid]
        if t.indels:
            ind = t.indels[0]
            deltas.append(law_table.del_open_delta(ind.run_copies, ind.run_unit or "A", args.lam))
    const = sum(deltas) / len(deltas) if deltas else 0.0
    print(f"[claimB] scored reads={len(ids)}  B0 const discount={const:.4f} "
          f"(mean law delta)", file=sys.stderr)

    arms = {"flat": (None, 1.0),
            "B0":   (ConstantDiscount(const), 1.0),
            "law":  (law_table, args.lam)}

    def score_ids(id_list, tag):
        sub = {rid: truth_by_id[rid] for rid in id_list}
        return {an: run_arm(read_seqs, refs, sub, genome, ref_fa,
                            os.path.join(args.out, f"{tag}_{an}.bam"), pt, lam)
                for an, (pt, lam) in arms.items()}

    def conc(s):
        return s.indel.position_exact_concordance

    print(f"\n============ C1 CLAIM-B: {args.stratum} placement concordance BY RUN LENGTH ============")
    print("(MECHANISM check on FLAT-in-L substrate; NOT Claim B — SIRV length-curve deferred)\n")
    print(f"{'run_len L':10s} {'n':>6s}   {'flat':>8s} {'B0':>8s} {'law':>8s}   {'law-B0':>8s}  verdict")
    per_len = {}
    agg = {an: [0, 0] for an in arms}   # [correct, incorrect] summed across lengths
    for L in sorted(by_len):
        sc = score_ids(by_len[L], f"L{L:02d}")
        for an in arms:
            agg[an][0] += sc[an].indel.correct
            agg[an][1] += sc[an].indel.incorrect
        cf, cb, cl = conc(sc["flat"]), conc(sc["B0"]), conc(sc["law"])
        per_len[L] = (len(by_len[L]), cf, cb, cl)
        v = "law>B0" if cl > cb + 1e-9 else ("law<B0" if cl < cb - 1e-9 else "law~B0")
        print(f"{L:<10d} {len(by_len[L]):6d}   {cf:8.4f} {cb:8.4f} {cl:8.4f}   {cl-cb:+8.4f}  {v}")

    # aggregate (equal-weight reads = flat-in-L injection) — summed from per-length
    def rate(an):
        c, i = agg[an]
        return c / (c + i) if (c + i) else float("nan")
    cf, cb, cl = rate("flat"), rate("B0"), rate("law")
    print(f"\n{'ALL (flat)':10s} {len(ids):6d}   {cf:8.4f} {cb:8.4f} {cl:8.4f}   {cl-cb:+8.4f}  "
          f"{'law>B0' if cl>cb+1e-9 else ('law<B0' if cl<cb-1e-9 else 'law~B0')}")

    # optional: reweight toward long runs (what a rising/SIRV injection would give)
    if args.reweight_long and per_len:
        # weight proportional to L (a simple monotone-rising injection surrogate)
        wf = wb = wl = wn = 0.0
        for L, (n, f, b, l) in per_len.items():
            w = L * n
            wf += w * f; wb += w * b; wl += w * l; wn += w
        print(f"{'reweight-L':10s} {'(w)':>6s}   {wf/wn:8.4f} {wb/wn:8.4f} {wl/wn:8.4f}   "
              f"{(wl-wb)/wn:+8.4f}  {'law>B0' if wl>wb else ('law<B0' if wl<wb else 'law~B0')}")

    # ---- verdict ----
    print("\n---- VERDICT ----")
    print(f"  HARNESS SANITY (aggregate, flat-in-L): B0 {'>=' if cb>=cl-1e-9 else '<'} law "
          f"(B0={cb:.4f} law={cl:.4f}) -- expected B0>=law reproduces Claim A.")
    longLs = [L for L in per_len if L >= 8]
    long_law_wins = [L for L in longLs if per_len[L][3] > per_len[L][2] + 1e-9]
    short_law_loses = [L for L in per_len if L <= 4 and per_len[L][3] < per_len[L][2] - 1e-9]
    if long_law_wins:
        print(f"  MECHANISM SIGNAL: law>B0 at LONG runs L={long_law_wins} "
              f"(and law<=B0 at short L={short_law_loses or 'none'}). This is the "
              f"length-SHAPE doing directional work: bigger in-run discount where runs "
              f"are long. A rising (SIRV-measured) injection would convert this into an "
              f"aggregate win. NOT Claim B yet — the injection here is flat-in-L.")
    else:
        print(f"  NO length-shape signal: law did not exceed B0 at any LONG run "
              f"(L>=8). Honest — report as-is; the flat-in-L placement task does not "
              f"reward the shape. Do NOT crank lambda (C1_DESIGN 75-77).")
    print("\n  GENUINE Claim B remains BLOCKED on the SIRV per-run-length del-rate curve.")


if __name__ == "__main__":
    main()
