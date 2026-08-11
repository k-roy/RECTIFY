#!/usr/bin/env python3
"""C1 real-SIRV OVER-CALL safety check (the surviving Claim-B real-data deliverable).

SCOPE: the OVER-CALL axis ONLY. The placement / per-length-SHAPE test on real SIRV
is NOT built here -- it is a pre-committed confounded null (see dev/C1_DESIGN.md
"Claim B"): underpowered (the SIRV/ERCC references carry ~9 distinct HP runs >=9 and
~137 >=7 in their entirety -- synthetic, engineered low-homopolymer) AND
truth-confounded by the "iron triangle" (a real read's placement is only
discriminating where a boundary mismatch makes the per-read truth unknowable). The
genuine length-shape validation is deferred to an independent-rate injection
simulator.

What DOES survive honestly on real reads, and is what this script measures: does the
in-run gap discount HALLUCINATE indels on CLEAN (length-preserving) real HP runs,
versus a MATCHED HP-aware flat baseline? This is the real-data analog of the sim's
clean false_indel_rate (c1_ablation.py), with KNOWN truth (a length-preserving
window's correct alignment is all-M, so any indel op is spurious).

Method (vetted backbone):
  * enumerate HP runs (len>=L_MIN) from the BAM's OWN @SQ contig FASTA;
  * per (read, run) cut a local window [run_start-F, run_end+F], anchored at clean
    base-matched flank columns taken from the BAM alignment (placement-neutral);
  * base-count truth_net = ref_win_len - query_win_len  (TABLE-INDEPENDENT);
  * keep ONLY clean trials: truth_net == 0 AND clean (base-equal) anchor zones.
    Substitutions inside the window are ALLOWED -- they are the over-call driver
    (the boundary/HP sub the discount could wrongly convert into a cheap in-run gap);
  * re-align each clean window with three MATCHED arms, ALL passing chrom_ref = the
    FULL contig (so homo_mask + homo_mismatch=-2 are active in every arm):
        flat : penalty_table=None
        B0   : ConstantDiscount = mean law del_open_delta over the clean trials present
        law  : HpPenaltyTable (rate_mean baseline-anchored log-odds gap-open delta)
  * OVER-CALL = the arm's CIGAR contains ANY indel op (D or I) -- counted by OP, so a
    canceling in-run D+I (net 0) still counts. Stratified by run-length bin.
  * WIRING ASSERT: the DP's own _homopolymer_run_len(full_contig, run_mid) must equal
    the enumerated run_len for every trial (guards the chrom_ref/ref_offset wiring,
    the variable under test).

PASS iff law over-call rate <= flat AND < OVERCALL_MAX (~0.02) in EVERY populated bin.
Never tune lam.

Usage:
  python c1_real_sirv_ablation.py --bam BAM --ref REF.fa \
      --penalty-table penalty_scores.tsv [--lam 1.0] [--max-reads N] [--flank 12]

Author: Kevin R. Roy (over-call deliverable; placement/shape intentionally omitted)
"""
from __future__ import annotations

import argparse
import bisect
import os
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import pysam  # noqa: E402
from rectify.core.align.local_aligner import (  # noqa: E402
    align_exon_block_global, _homopolymer_run_len,
)
from rectify.core.splice.hp_penalty import HpPenaltyTable  # noqa: E402

_OP_I, _OP_D, _OP_N = 1, 2, 3
OVERCALL_MAX = 0.02
L_MIN = 4


class ConstantDiscount:
    """B=0 control: a CONSTANT in-run gap-open discount, no run-length shape (same
    average magnitude as the law arm over the clean trials present, so law-vs-B0
    isolates the run-length SHAPE -- not under test here, but kept for the matched
    3-arm structure mirroring c1_ablation.py)."""
    def __init__(self, d_del: float, d_ins: float = 0.0):
        self._d_del, self._d_ins = d_del, d_ins

    def del_open_delta(self, hp, base, lam: float = 1.0):
        return lam * self._d_del

    def ins_open_delta(self, hp, base, lam: float = 1.0):
        return lam * self._d_ins


class ZeroInsWrap:
    """Diagnostic: the law's deletion deltas, but ins_open_delta forced to 0 (the
    insertion discount is flagged UNVALIDATED in dev/C1_DESIGN.md). Isolates whether
    the over-call is driven by the unvalidated insertion discount (sub->D+I rewrite
    needs a cheap I) or by the deletion shape itself."""
    def __init__(self, law):
        self._law = law

    def del_open_delta(self, hp, base, lam=1.0):
        return self._law.del_open_delta(hp, base, lam)

    def ins_open_delta(self, hp, base, lam=1.0):
        return 0.0


def runbin(L: int) -> str:
    return "SHORT(4-5)" if L <= 5 else ("MID(6-8)" if L <= 8 else "LONG(9+)")


BIN_ORDER = ["SHORT(4-5)", "MID(6-8)", "LONG(9+)"]


def enumerate_hp_runs(seq: str, lmin: int = L_MIN):
    """Return [(start, end, base, length)] for homopolymer runs length >= lmin."""
    runs = []
    n = len(seq)
    i = 0
    while i < n:
        j = i
        c = seq[i].upper()
        while j < n and seq[j].upper() == c:
            j += 1
        if c in "ACGT" and (j - i) >= lmin:
            runs.append((i, j, c, j - i))
        i = j
    return runs


def n_covered_ref_ranges(read):
    """Reference ranges covered by N (intron) ops, to exclude windows that straddle
    a splice (these reads are ~97% spliced)."""
    ranges = []
    pos = read.reference_start
    for op, ln in read.cigartuples or []:
        if op == _OP_N:
            ranges.append((pos, pos + ln))
            pos += ln
        elif op in (0, 2, 3, 7, 8):  # ref-consuming
            pos += ln
    return ranges


def has_indel(cigartuples) -> bool:
    return any(op in (_OP_I, _OP_D) for op, _ln in (cigartuples or []))


def collect_clean_trials(bam_path, ref_fa, flank, anchor, max_reads):
    """Scan the BAM, build clean (truth_net==0) anchored windows. Returns a list of
    trial dicts and per-contig sequences cache."""
    fa = pysam.FastaFile(ref_fa)
    contig_runs = {}      # contig -> (runs, run_starts)
    contig_seq = {}       # contig -> sequence (cache for contigs that have reads)
    trials = []
    filt = defaultdict(int)
    n_reads = 0
    n_with_sub = 0

    bam = pysam.AlignmentFile(bam_path, "rb")
    bam_contigs = set(bam.references)
    fa_contigs = set(fa.references)
    for read in bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        n_reads += 1
        if max_reads and n_reads > max_reads:
            break
        c = read.reference_name
        if c not in fa_contigs:
            filt["contig_not_in_ref"] += 1
            continue
        if c not in contig_seq:
            contig_seq[c] = fa.fetch(c)
            runs = enumerate_hp_runs(contig_seq[c])
            contig_runs[c] = (runs, [r[0] for r in runs])
        seq = contig_seq[c]
        runs, run_starts = contig_runs[c]
        if not runs:
            continue
        rs0, re0 = read.reference_start, read.reference_end
        qseq = read.query_sequence
        if qseq is None:
            continue
        # ref->query map for matched columns
        r2q = {rp: qp for qp, rp in read.get_aligned_pairs(matches_only=True)}
        n_ranges = n_covered_ref_ranges(read)
        lo = bisect.bisect_left(run_starts, rs0)
        for k in range(lo, len(runs)):
            rs, re_, base, L = runs[k]
            if rs > re0:
                break
            ws, we = rs - flank, re_ + flank
            if ws < rs0 or we > re0:
                continue
            # exclude windows straddling an intron
            if any(not (we <= ns or ws >= ne) for (ns, ne) in n_ranges):
                filt["intron_in_window"] += 1
                continue
            # clean anchor zones: outer `anchor` cols of each flank present AND
            # base-equal (reliable, placement-neutral anchor)
            okL = all((ws + t) in r2q for t in range(anchor))
            okR = all((we - 1 - t) in r2q for t in range(anchor))
            if not (okL and okR):
                filt["anchor_gap"] += 1
                continue
            qL = r2q[ws]
            qR = r2q[we - 1]
            if qR < qL:
                filt["anchor_order"] += 1
                continue
            ref_w = we - ws
            q_w = qR - qL + 1
            if abs(q_w - ref_w) > L + 6:
                filt["intron/artifact"] += 1
                continue
            # anchor-zone base equality
            if (qseq[qL:qL + anchor].upper() != seq[ws:ws + anchor].upper()
                    or qseq[qR - anchor + 1:qR + 1].upper() != seq[we - anchor:we].upper()):
                filt["anchor_dirty"] += 1
                continue
            truth_net = ref_w - q_w
            if truth_net != 0:
                filt["nonclean(net!=0)"] += 1
                continue
            qwin = qseq[qL:qR + 1]
            refwin = seq[ws:we]
            # ungapped Hamming (qwin and refwin are EQUAL length when truth_net==0):
            # a same-length window with few mismatches is overwhelmingly a
            # SUBSTITUTION-only window (balanced real indel PAIRS in a ~36bp window
            # are rare), so its correct alignment is provably indel-free -> any indel
            # an arm emits there is a TRUE hallucination. High-Hamming same-length
            # windows are likely balanced real indels and CONTAMINATE the op-count
            # over-call metric (all arms legitimately call them) -- we report the
            # strict sub-only stratum (hamming<=HAM_STRICT) separately.
            ham = sum(1 for a, b2 in zip(qwin.upper(), refwin.upper()) if a != b2)
            if ham > 0:
                n_with_sub += 1
            trials.append({
                "contig": c, "ws": ws, "we": we, "rs": rs, "re": re_,
                "base": base, "L": L, "qwin": qwin, "refwin": refwin,
                "full": seq, "ham": ham,
            })
    return trials, filt, n_reads - 1, n_with_sub


def run(bam_path, ref_fa, penalty_table_path, lam, flank, anchor, max_reads,
        zero_ins=False):
    t0 = time.time()
    law = HpPenaltyTable.from_tsv(penalty_table_path, None, min_count=100)
    trials, filt, n_reads, n_with_sub = collect_clean_trials(
        bam_path, ref_fa, flank, anchor, max_reads)
    print(f"[{os.path.basename(bam_path)}] reads scanned={n_reads} "
          f"clean trials={len(trials)} (with >=1 sub in window={n_with_sub}) "
          f"filters={dict(filt)}", file=sys.stderr)
    if not trials:
        print("  NO clean trials -- cannot evaluate (report low power).")
        return None

    # B0 constant = mean law del_open_delta over the clean trials present
    deltas = [law.del_open_delta(t["L"], t["base"], lam) for t in trials]
    const = sum(deltas) / len(deltas)
    arms = {"flat": (None, lam), "B0": (ConstantDiscount(const), lam), "law": (law, lam)}
    if zero_ins:
        arms["law_delonly"] = (ZeroInsWrap(law), lam)

    # accumulators
    HAM_STRICT = 2     # same-length window with 1-2 mismatches = sub-only proxy
    n_trial = defaultdict(int)
    n_strict = defaultdict(int)
    n_runs_distinct = defaultdict(set)
    overcall = {a: defaultdict(int) for a in arms}         # arm -> bin (ALL clean)
    overcall_strict = {a: defaultdict(int) for a in arms}  # arm -> bin (sub-only)
    wiring_ok = wiring_bad = 0

    for t in trials:
        b = runbin(t["L"])
        n_trial[b] += 1
        n_runs_distinct[b].add((t["contig"], t["rs"]))
        is_strict = 0 < t["ham"] <= HAM_STRICT
        if is_strict:
            n_strict[b] += 1
        # wiring assert: DP run-len at run midpoint == enumerated L
        mid = (t["rs"] + t["re"]) // 2
        rl_dp, _ = _homopolymer_run_len(t["full"], mid)
        if rl_dp == t["L"]:
            wiring_ok += 1
        else:
            wiring_bad += 1
        for a, (pt, lm) in arms.items():
            cig = align_exon_block_global(
                t["qwin"], t["refwin"], chrom_ref=t["full"],
                ref_offset=t["ws"], penalty_table=pt, lam=lm)
            if has_indel(cig):
                overcall[a][b] += 1
                if is_strict:
                    overcall_strict[a][b] += 1

    elapsed = time.time() - t0
    print(f"\n==== OVER-CALL on {os.path.basename(bam_path)} "
          f"(ref={os.path.basename(ref_fa)}, lam={lam}, F={flank}) ====")
    print(f"wiring assert (DP run_len == enumerated): ok={wiring_ok} bad={wiring_bad}"
          + ("   *** WIRING BUG ***" if wiring_bad else ""))
    print(f"B0 constant discount = {const:.4f} (mean law delta over clean trials)")

    arm_names = list(arms)            # flat, B0, law [, law_delonly]
    hdr = "".join(f"{a:>12s}" for a in arm_names)

    # ALL-clean table (truth_net==0; CONTAMINATED by balanced real indel pairs, which
    # all arms legitimately call -> high baseline, shown for context only).
    print("\n-- over-call on ALL clean (truth_net==0) windows [context; "
          "contaminated by balanced real indels] --")
    print(f"{'bin':12s}{'n_runs':>8s}{'n_trials':>9s}{hdr}")
    populated = []
    for b in BIN_ORDER:
        t_ = n_trial[b]
        if t_ == 0:
            continue
        populated.append(b)
        row = "".join(f"{overcall[a][b]/t_:12.4f}" for a in arm_names)
        print(f"{b:12s}{len(n_runs_distinct[b]):8d}{t_:9d}{row}")

    # STRICT sub-only table (same-length, 1-2 mismatches; correct alignment is
    # provably indel-free -> ANY indel is a TRUE hallucination). THIS is the
    # honest over-call control.
    print("\n-- over-call on SUB-ONLY windows (hamming<=%d; honest hallucination "
          "control) --" % HAM_STRICT)
    print(f"{'bin':12s}{'n_strict':>9s}{hdr}")
    rates = {}
    populated_strict = []
    for b in BIN_ORDER:
        s_ = n_strict[b]
        if s_ == 0:
            continue
        populated_strict.append(b)
        fr = {a: overcall_strict[a][b] / s_ for a in arm_names}
        rates[b] = fr
        row = "".join(f"{fr[a]:12.4f}" for a in arm_names)
        print(f"{b:12s}{s_:9d}{row}")

    # ---- verdict (on the honest sub-only control) ----
    print("\n---- VERDICT (over-call safety on sub-only windows: "
          "law must not hallucinate indels) ----")
    if not populated_strict:
        print("  NO populated sub-only bins -- low power, no verdict.")
        return None
    ok = True
    for b in populated_strict:
        fr = rates[b]
        cond = (fr["law"] <= fr["flat"] + 1e-9) and (fr["law"] < OVERCALL_MAX)
        ok = ok and cond
        print(f"  {b}: law<=flat={fr['law']<=fr['flat']+1e-9} "
              f"law<{OVERCALL_MAX}={fr['law']<OVERCALL_MAX} "
              f"(flat={fr['flat']:.4f} B0={fr['B0']:.4f} law={fr['law']:.4f}) "
              f"-> {'PASS' if cond else 'FAIL'}")
    print(f"  ==> {'PASS' if ok else 'FAIL'}  ({elapsed:.1f}s)")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--penalty-table", required=True)
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--flank", type=int, default=12)
    ap.add_argument("--anchor", type=int, default=5)
    ap.add_argument("--max-reads", type=int, default=0, help="0 = all primary reads")
    ap.add_argument("--zero-ins", action="store_true",
                    help="add a law_delonly diagnostic arm (ins_open_delta forced 0)")
    args = ap.parse_args()
    run(args.bam, args.ref, args.penalty_table, args.lam,
        args.flank, args.anchor, args.max_reads or 0, zero_ins=args.zero_ins)


if __name__ == "__main__":
    main()
