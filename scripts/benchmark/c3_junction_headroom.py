#!/usr/bin/env python3
"""C3 junction-arbitration headroom — the 0.09->1.07 locus (LLR-FREE).

The indel probe (c3_headroom.py) refuted C3-as-accuracy on the indel strata: the
shipped hp_edit_distance arbiter is at ceiling.  But the NAMED C3 justification is
the 0.09->1.07 artifact — a JUNCTION-arbitration claim (an aligner's apparent
quality flipped under a pure score re-weighting; the canonical-snap bias).  This
probe tests junction arbitration directly, on M1, with a constructed two-member
family per the SPEC's "fixed placements at fixed loci" case (line 127):

  * mm2   = real minimap2 -ax splice alignment.  On the JUNCTION_DISCOVERY stratum
            it SNAPS a non-canonical truth junction onto a nearby canonical GT-AG
            motif (smoke (G) verified this; the snap costs NM>=1 flanking mismatches).
  * truth = the truth-site placement (the corpus true_cigar) — what an orthogonal
            de-novo placer WOULD produce.  Position-exact by construction.

GOVERNING ARBITER = hp_edit_distance (merge_corrected_tsvs).  In the canonical
correct-first pipeline (run/single_sample.py, split_command.py) the EMITTED record
— junctions included — is the hp_ed winner's (write_corrected_consensus_bam takes
the winner's full corrected record).  hp_ed gives introns FREE (op 3 -> 0 cost),
so it discriminates junction placement ONLY through flanking-exon mismatches — and
the snap INDUCES those.  Prediction (advisor): hp_ed already prefers the truth-site
member -> arbiter at ceiling -> another null.  This MEASURES it.

  junction headroom = freq( truth-site member is position-exact
                            AND the hp_ed arbiter does NOT pick it )

If ~0 => hp_ed arbiter at ceiling on junctions too => C3-as-accuracy fully refuted
(the canonical-snap bias, where it exists, lives in select_best_alignment's
canonical_count TIEBREAKER — a tiebreaker-reweight / Discovery-facet concern, NOT a
calibrated-LLR one).  If the arbiter prefers the SNAP over an available truth
member => a real C3-addressable gap; then decide tiebreaker-reweight vs LLR.

NOTE the C5-vs-C3 boundary: in the REAL panel, if EVERY aligner snaps (no member
places truth), arbitration cannot recover it — that is C5/discovery, not C3.  This
probe INJECTS a truth member to isolate the pure arbitration question (given a
truth placement exists, does the arbiter pick it?).

Usage:
  python scripts/benchmark/c3_junction_headroom.py --out /tmp/c3_jhr --reps 30 \
      --penalty-table rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.scorer import load_genome, cigar_records_to_bam  # noqa: E402
from rectify.core.benchmark.truth_schema import read_truth_table, SplitTag  # noqa: E402
from rectify.core.splice.hp_penalty import HpPenaltyTable  # noqa: E402
from scripts.benchmark.sim import controlled  # noqa: E402
from scripts.benchmark.c3_headroom import (  # noqa: E402
    _cigar_hp_edit_distance, _cigar_aligned_bases, _read_is_position_exact,
    load_fastq, run_minimap2, scan_member,
)
import pysam  # noqa: E402

_CIG_OP = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
_CIG_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_cigar(s: str):
    return [(_CIG_OP[op], int(n)) for n, op in _CIG_RE.findall(s)]


def build_truth_member_bam(read_seqs, truth_subset, ref_fa, out_bam):
    """The truth-site placement member: each read at (chrom, genome_start) with its
    true_cigar — what an orthogonal de-novo placer would emit (position-exact)."""
    records = []
    for rid, t in truth_subset.items():
        seq = read_seqs.get(rid)
        if seq is None or not t.true_cigar:
            continue
        try:
            cig = parse_cigar(t.true_cigar)
        except Exception:
            continue
        # query-consuming ops must sum to len(seq); skip inconsistent rows loudly
        qcount = sum(n for op, n in cig if op in (0, 1, 4, 7, 8))
        if qcount != len(seq):
            continue
        records.append((rid, t.chrom, t.genome_start, cig, seq))
    cigar_records_to_bam(records, ref_fa, out_bam)
    return len(records)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/c3_jhr")
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--penalty-table", required=True)
    ap.add_argument("--str-table", default=None)
    ap.add_argument("--min-count", type=int, default=100)
    ap.add_argument("--strata",
                    default="JUNCTION_DISCOVERY,ANNOTATED,NIC,JUNCTION_AMB",
                    help="junction strata to score")
    ap.add_argument("--all-splits", action="store_true")
    args = ap.parse_args()
    want = {s.strip() for s in args.strata.split(",") if s.strip()}
    os.makedirs(args.out, exist_ok=True)

    print(f"[c3jhr] generating corpus reps={args.reps} ...", file=sys.stderr)
    info = controlled.generate_corpus(args.out, reps=args.reps, seed=7)
    genome = load_genome(info["ref_fa"])
    ref_fa = info["ref_fa"]
    read_seqs = load_fastq(info["reads_fastq"])
    truth = read_truth_table(info["truth_tsv"])
    law_table = HpPenaltyTable.from_tsv(args.penalty_table, args.str_table,
                                        min_count=args.min_count)

    keep = truth if args.all_splits else [t for t in truth if t.split == SplitTag.TEST]
    keep = [t for t in keep if t.stratum in want]
    truth_subset = {t.read_id: t for t in keep}
    print(f"[c3jhr] {len(truth_subset)} junction reads "
          f"({'ALL' if args.all_splits else 'TEST'} splits, strata={sorted(want)})",
          file=sys.stderr)

    # members: real minimap2 (the snap) + truth-site (orthogonal placer)
    mb = os.path.join(args.out, "jmember_mm2.bam")
    tb = os.path.join(args.out, "jmember_truth.bam")
    run_minimap2(ref_fa, info["reads_fastq"], mb)
    n_truth = build_truth_member_bam(read_seqs, truth_subset, ref_fa, tb)
    print(f"[c3jhr] truth-site member: {n_truth} reads placed", file=sys.stderr)
    members = {
        "mm2": scan_member(mb, truth_subset, genome, law_table),
        "truth": scan_member(tb, truth_subset, genome, law_table),
    }

    agg = defaultdict(lambda: dict(n=0, ceiling=0, arb_correct=0, headroom=0,
                                   mm2_snap=0, snap_and_truth=0,
                                   arb_picks_snap_over_truth=0))
    EPS = 1e-9
    for rid, t in truth_subset.items():
        present = [m for m in members if rid in members[m]]
        if "truth" not in present:
            continue  # truth member must be placed to ask the arbitration question
        vals = {m: members[m][rid] for m in present}
        s = agg[t.stratum]
        s["n"] += 1
        truth_exact = vals["truth"][2]
        if any(v[2] for v in vals.values()):
            s["ceiling"] += 1
        # shipped arbiter: argmin hp_ed, span DESC tiebreak
        min_ed = min(v[0] for v in vals.values())
        tied = [m for m in present if vals[m][0] <= min_ed + EPS]
        winner = max(tied, key=lambda m: vals[m][1])
        if vals[winner][2]:
            s["arb_correct"] += 1
        if truth_exact and not vals[winner][2]:
            s["headroom"] += 1
        # snap diagnostics: mm2 present and NOT position-exact == it snapped/missed
        if "mm2" in present and not vals["mm2"][2]:
            s["mm2_snap"] += 1
            if truth_exact:
                s["snap_and_truth"] += 1
                if winner == "mm2":
                    s["arb_picks_snap_over_truth"] += 1

    print("\n========== C3 JUNCTION-ARBITRATION HEADROOM (hp_ed arbiter, LLR-free) ==========")
    print(f"members=[mm2(real snap), truth(orthogonal placer)]  "
          f"split={'ALL' if args.all_splits else 'TEST'}  reps={args.reps}")
    print(f"{'stratum':20s} {'n':>5s} {'ceiling':>8s} {'arbiter':>8s} {'HEADROOM':>9s} "
          f"{'mm2snap':>8s} {'snap&tru':>9s} {'pickSnap':>9s}")
    tot = defaultdict(int)
    for st in sorted(agg):
        s = agg[st]
        n = max(1, s["n"])
        snt = max(1, s["snap_and_truth"])
        for k, v in s.items():
            tot[k] += v
        print(f"{st:20s} {s['n']:5d} {s['ceiling']/n:8.3f} {s['arb_correct']/n:8.3f} "
              f"{s['headroom']/n:9.3f} {s['mm2_snap']/n:8.3f} {s['snap_and_truth']/n:9.3f} "
              f"{s['arb_picks_snap_over_truth']/snt:9.3f}")
    N = max(1, tot["n"]); SNT = max(1, tot["snap_and_truth"])
    print("-" * 92)
    print(f"{'TOTAL':20s} {tot['n']:5d} {tot['ceiling']/N:8.3f} {tot['arb_correct']/N:8.3f} "
          f"{tot['headroom']/N:9.3f} {tot['mm2_snap']/N:8.3f} {tot['snap_and_truth']/N:9.3f} "
          f"{tot['arb_picks_snap_over_truth']/SNT:9.3f}")

    print("\n---- READING ----")
    print("ceiling   = freq a member is position-exact (truth member is exact by construction)")
    print("arbiter   = freq hp_ed arbiter picks a position-exact member")
    print("HEADROOM  = freq truth-site member exact BUT arbiter did NOT pick it")
    print("mm2snap   = freq minimap2 is NOT position-exact (snapped/missed the junction)")
    print("snap&tru  = freq mm2 snapped AND a truth-site member is available (the C3 setup)")
    print("pickSnap  = OF snap&truth reads, freq the hp_ed arbiter picked the SNAP over truth")
    print("            (>0 => a real junction-arbitration gap; =0 => hp_ed prefers truth = at ceiling)")
    print("\n---- VERDICT (pre-committed) ----")
    hr = tot["headroom"] / N
    psnap = tot["arb_picks_snap_over_truth"] / SNT
    if hr < 0.01 and psnap < 0.01:
        print(f"  HEADROOM={hr:.3f}, pickSnap={psnap:.3f} ~ 0 => hp_ed arbiter AT CEILING on")
        print("  junction placement too: given a truth member, it ALREADY beats the snap (introns")
        print("  free + snap induces flanking mismatches). => C3-as-accuracy FULLY REFUTED. The")
        print("  canonical-snap bias, where present, is a TIEBREAKER reweight (select_best_alignment")
        print("  canonical_count) / Discovery-facet concern, NOT a calibrated-LLR one. Do NOT build")
        print("  the LLR merely to host the near-tautological artifact-replay fence.")
    else:
        print(f"  HEADROOM={hr:.3f} pickSnap={psnap:.3f} > 0 => the arbiter prefers the SNAP over an")
        print("  available truth member: a REAL junction-arbitration gap. NEXT: identify what closes")
        print("  it — if it's the canonical_count tiebreaker, a one-line reweight (Discovery) suffices;")
        print("  only a residual coherence defect that a tiebreaker CANNOT fix justifies the LLR.")


if __name__ == "__main__":
    main()
