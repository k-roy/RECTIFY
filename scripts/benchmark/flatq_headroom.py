#!/usr/bin/env python3
"""Flat-Q quality-axis headroom probe — does per-base quality add ARBITRATION
signal the error-type table (hp_edit_distance) does NOT?  (the C3 flat-Q residue,
SPEC:225).  LLR-FREE in spirit: oracle Q is the UPPER BOUND, not a tuned decoder.

THE PREMISE (confirmed firsthand): RECTIFY is quality-BLIND — `query_qualities` is
propagated-but-unconsumed (chimeric_consensus.py:952, consensus.py:558; the only
reader, terminal_exon_refiner.py:714, slices it as a soft-clip passenger). Every
gate so far ran on FLAT-Q reads (`'I'*len`), so each tested only the error-TYPE
axis. This probe injects INFORMATIVE Q and asks whether it re-ranks the panel
toward truth where the error-type arbiter does NOT.

DESIGN (self-contained — touches NO shared sim file, so the flat-Q corpus is
trivially byte-identical):
  * Take CLEAN controlled-corpus reads (indel strata HP/HP_HARD/STR — where the
    exon-block DP arm is the valid tool, mirroring c3_headroom).
  * Corrupt each with SUBSTITUTION-ONLY background errors from a 2-state burst +
    per-read-multiplier hazard (a faithful sub-only replica of error_injector's
    layer-1/layer-2 process). Sub-only because injected background INDELs trip the
    scorer's has_unexplained gate (SESSION-2 lesson; verified 30/30 clean here).
  * ORACLE Q (advisor): Q_i = -10*log10(hazard_i) — Phred is DEFINED as
    -10log10 P(base wrong) and the hazard IS that probability, so oracle Q is the
    perfectly-calibrated basecaller, parameter-free, the MOST GENEROUS case for the
    Q-hypothesis. If even oracle Q can't beat hp_ed, no noisier real Q can
    (sidesteps the SIRV-gated placeholder-magnitude caveat — an upper bound needs
    no magnitude).

NON-CIRCULARITY FENCE (the arbiter consumes ONLY (alignment, realized Q string) —
never hazard/mult/burst/error-track; that barrier == deployment conditions):
  * read_mult is INERT for same-read arbitration (constant within a read → scales
    every placement equally; Panel A). Only the BURST factor localizes, so only
    burst can re-rank — exactly the channel hp_ed (reference-context-only) lacks.
  * PANEL-A CAVEAT (pre-registered): the injector hazard has NO reference-context
    term, so a positive oracle-Q result is NON-CREDIBLE (oracle Q would be reading
    the generative latent at zero noise, graded against the same latent; real
    basecaller Q is strongly ref-context-correlated, unlike here). Hence this is a
    ONE-SIDED test: a NULL (oracle Q adds 0 headroom) is DECISIVE -> flat-Q caveat
    moot; a WIN demands re-run under degraded Q + a ref-context hazard term before
    it counts. We report both and the per-read mechanism.

The two arbiters, each picks argmin over the panel, winner scored vs TRUTH:
  * hp_ed  : the shipped corrected-consensus cost (X=1, indels via penalty table,
             clips=1, introns free) — the incumbent error-type arbiter.
  * Qcost  : FULL-QUERY emission -log P(read | placement) from oracle Q. Every
             query base carries a cost (matched: -log(1-p); mismatch: -log(p/3);
             clip/insert: background -log(p) emission), so a soft-clip is NOT free
             (advisor full-query rule). This is the quality-aware arbiter.

Usage:
  python scripts/benchmark/flatq_headroom.py --out /tmp/flatq --reps 200 \
      --penalty-table rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import math
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.scorer import (  # noqa: E402
    load_genome, net_indel_in_span, all_indel_positions,
)
from rectify.core.benchmark.truth_schema import (  # noqa: E402
    read_truth_table, SplitTag, IndelKind,
)
from rectify.core.align.local_aligner import align_exon_block_global  # noqa: E402
from rectify.core.splice.hp_penalty import HpPenaltyTable, _hp_run_length  # noqa: E402
from scripts.benchmark.sim import controlled  # noqa: E402
import pysam  # noqa: E402

BASES = "ACGT"


# ---------------------------------------------------------------------------
# Sub-only burst corruption WITH oracle Q (self-contained replica of
# error_injector's layer-1 multiplier + layer-2 2-state burst; INDELS dropped so
# the scorer's has_unexplained gate stays clean). Errors AND Q come from the SAME
# walk => self-consistent. RNG is private to this probe.
# ---------------------------------------------------------------------------
def corrupt_with_oracle_q(clean: str, rng: random.Random,
                          base_rate: float, gamma_shape: float,
                          burst_on: bool, hot_factor: float,
                          p_cold_to_hot: float, p_hot_to_cold: float):
    """Return (dirty_seq, oracle_phred[list], n_err). Sub-only: len preserved, so
    query pos == clean pos and the Q array aligns 1:1 with the read."""
    # layer-1 per-read multiplier (E[m]=1)
    mult = rng.gammavariate(gamma_shape, 1.0 / gamma_shape) if gamma_shape > 0 else 1.0
    # cold factor so stationary-weighted mean == 1 (burst redistributes, not adds)
    if burst_on and (p_cold_to_hot + p_hot_to_cold) > 0:
        ph = p_cold_to_hot / (p_cold_to_hot + p_hot_to_cold)
        pc = 1.0 - ph
        cold_f = max(0.0, (1.0 - ph * hot_factor) / pc) if pc > 0 else 1.0
    else:
        cold_f = 1.0
    out = []
    quals = []
    n_err = 0
    state_hot = False
    for ch in clean:
        if burst_on:
            if state_hot:
                if rng.random() < p_hot_to_cold:
                    state_hot = False
            elif rng.random() < p_cold_to_hot:
                state_hot = True
            factor = hot_factor if state_hot else cold_f
        else:
            factor = 1.0
        hazard = min(1.0, base_rate * mult * factor)
        # ORACLE Phred (clamp hazard away from 0/1 for finite Q in [2, 60])
        h = min(0.999, max(1e-6, hazard))
        q = int(round(min(60.0, max(2.0, -10.0 * math.log10(h)))))
        quals.append(q)
        if rng.random() < hazard:           # substitution event
            out.append(rng.choice([b for b in BASES if b != ch]))
            n_err += 1
        else:
            out.append(ch)
    return "".join(out), quals, n_err


# ---------------------------------------------------------------------------
# Incumbent error-type arbiter cost (VERBATIM hp_edit_distance, as c3_headroom).
# ---------------------------------------------------------------------------
def cigar_hp_edit_distance(read, genome, penalty_table) -> float:
    if read.cigartuples is None:
        return 0.0
    chrom = read.reference_name
    gseq = (genome or {}).get(chrom) if genome else None
    q = read.query_sequence
    rp = read.reference_start
    qp = 0
    total = 0.0
    for op, length in read.cigartuples:
        if op == 7:
            rp += length; qp += length
        elif op == 8:
            total += length; rp += length; qp += length
        elif op == 0:
            if gseq and q:
                rc = gseq[rp:rp + length].upper(); qc = q[qp:qp + length].upper()
                total += sum(a != b for a, b in zip(rc, qc))
            rp += length; qp += length
        elif op == 2:
            for i in range(length):
                r = rp + i
                if penalty_table is not None and gseq and r < len(gseq):
                    total += penalty_table.del_cost(_hp_run_length(gseq, r), gseq[r])
                else:
                    total += 1.0
            rp += length
        elif op == 1:
            if penalty_table is not None and gseq and rp < len(gseq):
                total += length * penalty_table.ins_cost(_hp_run_length(gseq, rp), gseq[rp])
            else:
                total += length * 1.25
            qp += length
        elif op == 3:
            rp += length
        elif op == 4:
            total += length; qp += length
        elif op == 5:
            total += length
    return total


# ---------------------------------------------------------------------------
# Quality-aware arbiter: full-query emission -log P(read | placement) from oracle Q.
# Every query base carries a cost so a soft-clip is NOT free (advisor full-query
# rule). matched base j: -log(1-p_j); mismatch j: -log(p_j/3); inserted/clipped j:
# background -log(p_j) (an unexplained emission). Deletions: a fixed gap penalty
# (no query base) so the two arbiters compare comparable event sets.
# ---------------------------------------------------------------------------
def cigar_q_cost(read, genome, quals, del_gap=4.0) -> float:
    if read.cigartuples is None or quals is None:
        return 0.0
    gseq = (genome or {}).get(read.reference_name)
    q = read.query_sequence
    rp = read.reference_start
    qp = 0
    total = 0.0

    def p_of(j):
        Q = quals[j] if 0 <= j < len(quals) else 30
        return min(0.75, max(1e-6, 10.0 ** (-Q / 10.0)))

    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            for k in range(length):
                j = qp + k
                p = p_of(j)
                if gseq and q and rp + k < len(gseq) and j < len(q):
                    match = gseq[rp + k].upper() == q[j].upper()
                else:
                    match = (op == 7)
                total += -math.log(1.0 - p) if match else -math.log(p / 3.0)
            rp += length; qp += length
        elif op == 2:
            total += del_gap * length
            rp += length
        elif op == 1:
            for k in range(length):
                total += -math.log(p_of(qp + k))   # unexplained emission
            qp += length
        elif op == 3:
            rp += length
        elif op in (4, 5):
            for k in range(length):
                total += -math.log(p_of(qp + k))   # clip is NOT free
            if op == 4:
                qp += length
    return total


def read_is_position_exact(read, rt) -> bool:
    """Indel-strata position-exact rule (no junctions in HP/HP_HARD/STR)."""
    cig = read.cigartuples
    if cig is None:
        return False
    rstart = read.reference_start
    truth_spans = [(ind.eq_start, ind.eq_end) for ind in rt.indels]
    for (ipos, ilen, ikind) in all_indel_positions(rstart, cig):
        covered = (any(s <= ipos <= e for s, e in truth_spans) if ikind == 1
                   else any(s <= ipos < e for s, e in truth_spans))
        if not covered:
            return False
    for ind in rt.indels:
        tn = ind.length if ind.kind == IndelKind.DEL else -ind.length
        in_span, _ = net_indel_in_span(rstart, cig, ind.eq_start, ind.eq_end)
        if in_span != tn:
            return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/flatq")
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--penalty-table", required=True)
    ap.add_argument("--str-table", default=None)
    ap.add_argument("--min-count", type=int, default=100)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--base-rate", type=float, default=0.03)
    ap.add_argument("--gamma-shape", type=float, default=0.5)
    ap.add_argument("--hot-factor", type=float, default=8.0)
    ap.add_argument("--p-cold-to-hot", type=float, default=0.05)
    ap.add_argument("--p-hot-to-cold", type=float, default=0.2)
    ap.add_argument("--strata", default="HP,HP_HARD,STR")
    ap.add_argument("--all-splits", action="store_true", default=True)
    args = ap.parse_args()
    want = {s.strip() for s in args.strata.split(",") if s.strip()}
    os.makedirs(args.out, exist_ok=True)

    print(f"[flatq] corpus reps={args.reps}", file=sys.stderr)
    info = controlled.generate_corpus(args.out, reps=args.reps, seed=args.seed)
    genome = load_genome(info["ref_fa"])
    ref_fa = info["ref_fa"]
    truth = read_truth_table(info["truth_tsv"])
    seqs = {}
    with pysam.FastxFile(info["reads_fastq"]) as fq:
        for e in fq:
            seqs[e.name] = e.sequence

    law = HpPenaltyTable.from_tsv(args.penalty_table, args.str_table, min_count=args.min_count)
    keep = [t for t in truth if t.stratum in want and t.indels]
    print(f"[flatq] {len(keep)} indel-strata reads", file=sys.stderr)

    rng = random.Random(args.seed + 101)  # PRIVATE rng — never the corpus rng

    # corrupt each read; build flat + law member alignments on the DIRTY read;
    # score each member vs truth AND compute hp_ed + Q costs.
    agg = defaultdict(lambda: dict(n=0, ceiling=0, hp_corr=0, q_corr=0,
                                   disagree=0, hp_hr=0, q_hr=0,
                                   q_gain=0, q_break=0, q_gain_dis=0))
    EPS = 1e-9
    n_err_tot = 0
    n_base_tot = 0
    for t in keep:
        clean = seqs.get(t.read_id)
        if clean is None:
            continue
        dirty, quals, n_err = corrupt_with_oracle_q(
            clean, rng, args.base_rate, args.gamma_shape, True,
            args.hot_factor, args.p_cold_to_hot, args.p_hot_to_cold)
        n_err_tot += n_err
        n_base_tot += len(clean)
        ref = genome[t.chrom]
        members = {}
        for name, ptab in (("flat", None), ("law", law)):
            cig = align_exon_block_global(dirty, ref, chrom_ref=ref, ref_offset=0,
                                          penalty_table=ptab, lam=1.0)
            # wrap into a lightweight read-like via a temp AlignedSegment
            seg = pysam.AlignedSegment()
            seg.query_name = t.read_id
            seg.query_sequence = dirty
            seg.reference_id = 0
            seg.reference_start = 0
            seg.cigartuples = cig
            # reference_name needs a header; emulate by stashing chrom on a shim
            members[name] = (cig, seg)

        # score members (need reference_name -> patch via a header-bound seg)
        # Build a minimal header so reference_name resolves to t.chrom.
        hdr = pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.6"}, "SQ": [{"SN": t.chrom, "LN": len(ref)}]})
        vals = {}
        for name, (cig, _seg) in members.items():
            seg = pysam.AlignedSegment(hdr)
            seg.query_name = t.read_id
            seg.query_sequence = dirty
            seg.reference_id = 0
            seg.reference_start = 0
            seg.cigartuples = cig
            exact = read_is_position_exact(seg, t)
            hp = cigar_hp_edit_distance(seg, genome, law)
            qc = cigar_q_cost(seg, genome, quals)
            vals[name] = (hp, qc, exact)

        st = t.stratum
        s = agg[st]
        s["n"] += 1
        any_exact = any(v[2] for v in vals.values())
        if any_exact:
            s["ceiling"] += 1
        n_exact = sum(1 for v in vals.values() if v[2])
        is_dis = 0 < n_exact < len(vals)
        if is_dis:
            s["disagree"] += 1

        # hp_ed arbiter: argmin hp, span tiebreak omitted (2-member, deterministic
        # by min then name) — pick min hp; tie -> first
        hp_win = min(vals, key=lambda m: (vals[m][0], m))
        q_win = min(vals, key=lambda m: (vals[m][1], m))
        if vals[hp_win][2]:
            s["hp_corr"] += 1
        if vals[q_win][2]:
            s["q_corr"] += 1
        if any_exact and not vals[hp_win][2]:
            s["hp_hr"] += 1
        if any_exact and not vals[q_win][2]:
            s["q_hr"] += 1
        # Q GAIN: hp wrong, Q right (the headroom Q fills over the error-type table)
        if any_exact and not vals[hp_win][2] and vals[q_win][2]:
            s["q_gain"] += 1
            if is_dis:
                s["q_gain_dis"] += 1
        # Q BREAK (safety): hp right, Q wrong
        if vals[hp_win][2] and not vals[q_win][2]:
            s["q_break"] += 1

    # ---- report ----
    lines = []
    def P(x): lines.append(x); print(x)
    P("\n========= FLAT-Q QUALITY-AXIS ARBITRATION HEADROOM (oracle Q) =========")
    P(f"panel=[flat,law]  arbiters=[hp_ed, Qcost]  reps={args.reps}  "
      f"meas_sub_rate={n_err_tot/max(1,n_base_tot):.4f}")
    P(f"{'stratum':10s} {'n':>5s} {'ceil':>6s} {'hp_ed':>6s} {'Qcost':>6s} "
      f"{'disagr':>7s} {'Qgain':>6s} {'Qbreak':>7s} {'Qgain|dis':>9s}")
    tot = defaultdict(int)
    for st in sorted(agg):
        s = agg[st]
        n = max(1, s["n"]); nd = max(1, s["disagree"])
        for k in s: tot[k] += s[k]
        P(f"{st:10s} {s['n']:5d} {s['ceiling']/n:6.3f} {s['hp_corr']/n:6.3f} "
          f"{s['q_corr']/n:6.3f} {s['disagree']/n:7.3f} {s['q_gain']/n:6.3f} "
          f"{s['q_break']/n:7.3f} {s['q_gain_dis']/nd:9.3f}")
    N = max(1, tot["n"]); ND = max(1, tot["disagree"])
    P("-" * 78)
    P(f"{'TOTAL':10s} {tot['n']:5d} {tot['ceiling']/N:6.3f} {tot['hp_corr']/N:6.3f} "
      f"{tot['q_corr']/N:6.3f} {tot['disagree']/N:7.3f} {tot['q_gain']/N:6.3f} "
      f"{tot['q_break']/N:7.3f} {tot['q_gain_dis']/ND:9.3f}")
    P("\n---- READING ----")
    P("ceil      = freq >=1 member position-exact (recoverable universe)")
    P("hp_ed     = freq the incumbent error-type arbiter picks a correct member")
    P("Qcost     = freq the oracle-Q arbiter picks a correct member")
    P("Qgain     = freq hp_ed WRONG and Qcost RIGHT (headroom Q fills over error-type)")
    P("Qbreak    = freq hp_ed RIGHT and Qcost WRONG (SAFETY — Q must not break these)")
    P("Qgain|dis = Qgain restricted to member-disagreement reads (the strong test)")
    qg = tot["q_gain"] / N
    P("\n---- VERDICT (pre-committed, ONE-SIDED per Panel A) ----")
    if qg < 0.005:
        P(f"  Qgain={qg:.4f} ~ 0  => even ORACLE Q adds NO arbitration headroom over the")
        P("  error-type table => DECISIVE NULL: the flat-Q caveat is MOOT, the error-type")
        P("  gates stand. (No noisier real Q can beat an upper bound that already fails.)")
    else:
        P(f"  Qgain={qg:.4f} > 0  => oracle Q fills headroom, BUT this is NON-CREDIBLE per")
        P("  Panel A (the injector hazard has no reference-context term; oracle Q reads the")
        P("  generative latent at zero noise). Before counting: re-run under DEGRADED Q +")
        P("  add a ref-context hazard term so hp_ed and Q genuinely compete. Not actionable")
        P("  on this corpus.")

    with open(os.path.join(os.path.dirname(__file__), "flatq_headroom_result.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
