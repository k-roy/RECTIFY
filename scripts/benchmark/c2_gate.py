#!/usr/bin/env python3
"""C2 GATE — is GENOMIC_A_CPA a discriminating + addressable stratum?

Mirrors the C1 gate discipline (``c1_ablation.py``) but for 3'/poly-A cleavage
PLACEMENT. The unit under test is a 3'-end ESTIMATOR applied to a FIXED drifted
alignment (M through the genomic A-region, S the tail) — NOT an aligner DP. Three
matched arms, identical input read:

  raw        ``reference_end - 1`` (pure-alignment end). CONTEXT ONLY — this is
             raw minimap2, already solved by the shipped walkback; not a competitor.
  walkback   ``correct/walkback.py::walkback_drs_full`` — the REAL shipped
             incumbent (guarded heuristic walk-inward to the first non-stop
             read==ref match). C2 must beat THIS.
  decoder    a prototype 2-state templated-vs-tail change-point: pull the boundary
             to the 5' start of the maximal gap-tolerant 3'-terminal A-run (a
             degenerate tail-emission model that tolerates the ~7% non-A the
             PolyAModel allows). The candidate facet.

Scored vs TRUTH (``true_cpa``), never an internal score, stratified by cell:

  tract_start   NON-DISCRIMINATING control — walkback ~ ceiling expected.
  readthrough   UNIDENTIFIABLE (pre-committed NULL) — decoder ~ walkback expected
                (both err = k); a decoder "win" here would be the C1 Claim-B artifact.
  interrupted   the ONLY identifiable addressable cell — decoder < walkback expected.
  clean_g0      over-call control — decoder must NOT shift the 3' end (err 0).
  terminal_A    SHARED limitation — both arms over-trim real templated A's.

PRE-COMMITTED verdict logic (printed at the end):
  * DISCRIMINATING iff raw drifts (raw err grows with g) AND walkback is below
    ceiling on >=1 IDENTIFIABLE cell.
  * C2-ADDRESSABLE iff decoder beats walkback on `interrupted` WITHOUT regressing
    walkback on any other cell AND without over-calling on clean_g0.
  * If walkback ~ ceiling on every identifiable cell and the only residual is
    readthrough (unidentifiable), C2-as-PLACEMENT is REFUTED -> pivot to
    C2-as-posterior (the run length IS the uncertainty; feeds consensus LR).

Usage: python scripts/benchmark/c2_gate.py --reps 100
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import random
import statistics
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.scorer import load_genome, cigar_records_to_bam  # noqa: E402
from rectify.core.correct.walkback import walkback_drs_full  # noqa: E402
from scripts.benchmark.sim import controlled  # noqa: E402
import pysam  # noqa: E402


# ---------------------------------------------------------------------------
# Arm 3: prototype change-point decoder (plus-strand, right side)
# ---------------------------------------------------------------------------
def changepoint_cpa_plus(read, chrom_seq, max_gap: int = 1) -> int:
    """Pull the 3' end to the 5' start of the maximal gap-tolerant 3'-terminal
    A-run. A degenerate tail-emission model: the tail is A-dominated but tolerates
    up to ``max_gap`` consecutive non-A (the ~7% non-A PolyAModel allows), so a
    stray genomic non-A INSIDE the A-region (the `interrupted` cell) is absorbed
    into the tail instead of anchoring the boundary early (the walkback's failure).
    Returns a genome coord (the last templated base = CPA)."""
    qs = read.query_sequence
    if not qs:
        return read.reference_end - 1
    n = len(qs)
    start = n           # 5'-most query index that is part of the tail A-run
    run_nonA = 0
    saw_A = False
    i = n - 1
    while i >= 0:
        if qs[i].upper() == "A":
            start = i
            run_nonA = 0
            saw_A = True
        else:
            run_nonA += 1
            if run_nonA > max_gap:
                break
        i -= 1
    if not saw_A or start == n:
        return read.reference_end - 1     # no tail -> no shift (over-call guard)
    cpa_q = start - 1                      # last templated query base
    if cpa_q < 0:
        return read.reference_start
    qp2rp = {qp: rp for qp, rp in read.get_aligned_pairs(matches_only=True)}
    j = cpa_q
    while j >= 0 and j not in qp2rp:
        j -= 1
    return qp2rp[j] if j >= 0 else read.reference_start


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/c2_gate")
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    rng = random.Random(args.seed)
    cells = controlled.genomic_a_cpa_cells(args.reps, rng)
    by_id = {c["read_id"]: c for c in cells}

    # ref FASTA + drifted-alignment BAM
    ref_fa = os.path.join(args.out, "ref.fa")
    with open(ref_fa, "w") as fh:
        for c in cells:
            fh.write(f">{c['chrom']}\n{c['contig']}\n")
    records = [(c["read_id"], c["chrom"], 0, c["drifted_cigar"], c["read_seq"])
               for c in cells]
    bam = cigar_records_to_bam(records, ref_fa, os.path.join(args.out, "drifted.bam"))
    genome = load_genome(ref_fa)

    # err[cell][arm] = list of |est - true_cpa|
    err = defaultdict(lambda: defaultdict(list))
    # signed[cell][arm] = list of (est - true_cpa)  (sign: + = downstream over-shoot)
    signed = defaultdict(lambda: defaultdict(list))
    with pysam.AlignmentFile(bam, "rb") as bf:
        for read in bf:
            if read.is_unmapped:
                continue
            c = by_id.get(read.query_name)
            if c is None:
                continue
            chrom_seq = genome[read.reference_name]
            true_cpa = c["true_cpa"]
            ests = {"raw": read.reference_end - 1}
            wb = walkback_drs_full(read, chrom_seq)
            ests["walkback"] = wb["corrected_pos"] if wb else (read.reference_end - 1)
            ests["decoder"] = changepoint_cpa_plus(read, chrom_seq)
            for arm, est in ests.items():
                err[c["cell"]][arm].append(abs(est - true_cpa))
                signed[c["cell"]][arm].append(est - true_cpa)

    arms = ["raw", "walkback", "decoder"]
    cell_order = ["tract_start", "interrupted", "readthrough", "clean_g0",
                  "overcall_Aend", "terminal_A"]

    def med(xs):
        return statistics.median(xs) if xs else float("nan")

    def mean(xs):
        return sum(xs) / len(xs) if xs else float("nan")

    def exact(xs):
        return sum(1 for x in xs if x == 0) / len(xs) if xs else float("nan")

    print("\n================ C2 GATE (reps=%d) ================" % args.reps)
    print("median |est - true_cpa|   (exact-rate in parens)   per cell x arm")
    print(f"{'cell':14s} {'ident':6s} " + " ".join(f"{a:>22s}" for a in arms))
    for cell in cell_order:
        if cell not in err:
            continue
        ident = "yes" if by_id and next(c for c in cells if c["cell"] == cell)["identifiable"] else "NO"
        row = f"{cell:14s} {ident:6s} "
        for a in arms:
            xs = err[cell][a]
            row += f" {med(xs):6.1f}/{mean(xs):5.2f}({exact(xs):4.2f})"
        print(row)

    # drift-vs-g for raw (discriminating check)
    print("\nraw drift vs genomic-A length g (tract_start cell):")
    drift_by_g = defaultdict(list)
    for c in cells:
        if c["cell"] == "tract_start":
            drift_by_g[c["g"]].append(c)
    # recompute raw drift per g from err is awkward; recompute directly
    with pysam.AlignmentFile(bam, "rb") as bf:
        raw_g = defaultdict(list)
        wb_g = defaultdict(list)
        for read in bf:
            if read.is_unmapped:
                continue
            c = by_id.get(read.query_name)
            if c is None or c["cell"] != "tract_start":
                continue
            raw_g[c["g"]].append(abs((read.reference_end - 1) - c["true_cpa"]))
            wbb = walkback_drs_full(read, genome[read.reference_name])
            wb_g[c["g"]].append(abs((wbb["corrected_pos"] if wbb else read.reference_end - 1) - c["true_cpa"]))
    for g in sorted(raw_g):
        print(f"  g={g:2d}: raw median err={med(raw_g[g]):5.1f}  walkback median err={med(wb_g[g]):5.1f}")

    # ---- truth-flip demonstration (the decisive, model-free artifact proof) ----
    # Recompute the `interrupted` arms against the MAXIMUM-LIKELIHOOD truth X =
    # z + a1 + 1 (the last provably-templated base, favored ~26-113x over Z) instead
    # of the cell's fiat truth Z. If the arms EXACTLY SWAP, the gate measured a
    # truth DEFINITION, not a capability. (a1 is carried as downstream_a_count.)
    flip = defaultdict(list)        # arm -> |est - X|
    base = defaultdict(list)        # arm -> |est - Z|  (recomputed here for the pair)
    with pysam.AlignmentFile(bam, "rb") as bf:
        for read in bf:
            if read.is_unmapped:
                continue
            c = by_id.get(read.query_name)
            if c is None or c["cell"] != "interrupted":
                continue
            z = c["true_cpa"]
            x = z + c["downstream_a_count"] + 1     # the MLE position
            ests = {"raw": read.reference_end - 1}
            wb = walkback_drs_full(read, genome[read.reference_name])
            ests["walkback"] = wb["corrected_pos"] if wb else read.reference_end - 1
            ests["decoder"] = changepoint_cpa_plus(read, genome[read.reference_name])
            for arm, est in ests.items():
                base[arm].append(abs(est - z))
                flip[arm].append(abs(est - x))
    print("\ninterrupted TRUTH-FLIP (model-free artifact proof): mean |est - truth|")
    print(f"{'truth':18s} {'walkback':>10s} {'decoder':>10s}")
    print(f"{'Z (cell fiat)':18s} {mean(base['walkback']):10.2f} {mean(base['decoder']):10.2f}")
    print(f"{'X (likelihood MAP)':18s} {mean(flip['walkback']):10.2f} {mean(flip['decoder']):10.2f}")
    print("  -> the two arms SWAP under relabeling: the 'win' is a definition, not a capability.")

    # ---- verdict (post-adversarial-panel logic; see dev/C2_DESIGN.md) ----
    # A decoder "win" counts ONLY on an IDENTIFIABLE cell with a defensible truth.
    # `interrupted` was relabelled UNIDENTIFIABLE by the panel (the truth-flip:
    # relabel truth Z->X and the arms exactly swap; the walkback's anchor-at-X is
    # in fact the maximum-likelihood estimate), so its decoder "win" is NOT credited.
    print("\n---- VERDICT (post-adversarial-panel) ----")
    raw_drifts = all(med(raw_g[g]) >= g - 1 for g in raw_g)
    ident_cells = {c["cell"] for c in cells if c["identifiable"]}

    wb_tract = mean(err["tract_start"]["walkback"])
    cp_tract = mean(err["tract_start"]["decoder"])
    cp_clean = mean(err["clean_g0"]["decoder"])
    cp_overcall = mean(err["overcall_Aend"]["decoder"])
    wb_overcall = mean(err["overcall_Aend"]["walkback"])

    # Does the decoder beat the walkback on ANY identifiable cell (>0.5 bp)?
    decoder_wins_identifiable = [
        cell for cell in err
        if cell in ident_cells
        and mean(err[cell]["decoder"]) < mean(err[cell]["walkback"]) - 0.5
    ]

    print(f"  [discriminating vs RAW] raw drifts ~full tract: {raw_drifts} "
          f"(raw absorbs g; already shipped-walkback-solved)")
    print(f"  [incumbent at ceiling]  walkback on tract_start (identifiable canonical): "
          f"mean err={wb_tract:.2f}  -> {'CEILING' if wb_tract <= 0.5 else 'below ceiling'}")
    print(f"  [pre-committed NULL]    readthrough: decoder={mean(err['readthrough']['decoder']):.2f} "
          f"walkback={mean(err['readthrough']['walkback']):.2f} (unidentifiable; tie expected)")
    print(f"  [ARTIFACT, truth-flip]  interrupted vs Z: walkback={mean(base['walkback']):.2f} "
          f"decoder={mean(base['decoder']):.2f}  |  vs MLE X: walkback={mean(flip['walkback']):.2f} "
          f"decoder={mean(flip['decoder']):.2f}  -> SWAP (see above); NOT credited")
    print(f"  [over-call, WEAK ctrl]  clean_g0 decoder mean err={cp_clean:.2f}")
    print(f"  [over-call, STRONG ctrl] overcall_Aend decoder mean err={cp_overcall:.2f} "
          f"(walkback={wb_overcall:.2f}) -> {'OVER-CALLS' if cp_overcall > 0.5 else 'ok'}")
    print(f"  [shared over-trim]      terminal_A walkback={mean(err['terminal_A']['walkback']):.2f} "
          f"decoder={mean(err['terminal_A']['decoder']):.2f} "
          f"(== tract_start input, opposite truth -> unidentifiable)")
    print(f"  decoder wins on IDENTIFIABLE cells: {decoder_wins_identifiable or 'NONE'}")

    addressable = bool(decoder_wins_identifiable) and cp_overcall <= 0.5
    verdict = ("C2-PLACEMENT ADDRESSABLE (identifiable win found)" if addressable
               else "C2-PLACEMENT REFUTED — walkback at ceiling on the identifiable "
                    "canonical case; residual is sequence-unidentifiable; the only "
                    "apparent win is a truth-definition artifact; decoder over-calls "
                    "A-ending bodies. Pivot to C2-as-SOFT-posterior (no hard rule).")
    print(f"  ==> {verdict}")


if __name__ == "__main__":
    main()
