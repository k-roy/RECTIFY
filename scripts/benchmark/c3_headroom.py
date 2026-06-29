#!/usr/bin/env python3
"""C3 headroom probe — is the shipped consensus ARBITER below ceiling on
reads it COULD recover by re-ranking the panel?  (LLR-FREE; run this FIRST.)

C3's premise: replace the integer-max consensus arbiter with a calibrated
likelihood-ratio arbiter.  Per the C1/C2 gate discipline, before building any
LLR decoder we must establish the incumbent ARBITER is BELOW CEILING on
*recoverable* reads AND that the gap is arbitration-addressable.  This probe
measures exactly that, with NO LLR involved (so there is zero risk of a
construction-tuned / hill-climb-into-the-model win):

    arbitration headroom =
        freq( >=1 panel member is position-exact-correct vs TRUTH
              AND the shipped arbiter does NOT select a correct member )

The shipped corrected-consensus arbiter (rectify/core/consensus/corrected_consensus.py,
the use_hp_ed path) picks the winning aligner per read by:
        sort [ _effective_chimera_ok ASC, hp_edit_distance ASC, _span DESC ], .first()
On the controlled single-exon strata the chimera flag is inert, so the arbiter
reduces to argmin(hp_edit_distance) with an alignment-span tiebreak — which is
exactly what this probe replays, scoring the picked member against truth.

If headroom ~ 0  => the arbiter is already at ceiling on recoverable reads
(the C2 outcome: a phantom facet) => REFUTE C3-as-accuracy, document, stop.
If headroom > 0  => characterize it (mis-order vs tie-arbitrary), then and only
then build the coherent LLR and test whether it recovers those reads.

Panel members (the controlled, truth-scored disagreement the gate needs):
  * flat : align_exon_block_global(penalty_table=None)  — the flat-affine family
           (the panel's shared quality-blind bias; C1's matched baseline arm).
  * law  : align_exon_block_global(penalty_table=LAW)   — the C1 native member.
  * mm2  : real minimap2 -ax splice (optional 3rd, --with-mm2) — an INDEPENDENT
           external view (its placements are not produced by any table, the
           hold-out that kills construction-tuned wins).

Fitness = TRUTH (never any internal score), TEST split, ambiguity-aware
(net_indel_in_span / the _score_read position-exact rule), per stratum.

Usage:
  python scripts/benchmark/c3_headroom.py --out /tmp/c3_hr --reps 200 \
      --penalty-table rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.scorer import (  # noqa: E402
    load_genome, cigar_records_to_bam, net_indel_in_span, all_indel_positions,
    extract_junctions,
)
from rectify.core.benchmark.truth_schema import (  # noqa: E402
    read_truth_table, SplitTag, IndelKind,
)
from rectify.core.align.local_aligner import align_exon_block_global  # noqa: E402
from rectify.core.splice.hp_penalty import HpPenaltyTable, _hp_run_length  # noqa: E402
from rectify.core.consensus.chimeric_consensus import normalize_junction  # noqa: E402
from scripts.benchmark.sim import controlled  # noqa: E402
import pysam  # noqa: E402


# --- VERBATIM copies of the shipped corrected-consensus arbiter cost fns
# (rectify/core/consensus/corrected_consensus.py:_cigar_hp_edit_distance /
# _cigar_aligned_bases).  Inlined ONLY because corrected_consensus.py imports
# pandas+Bio at module top, which the thin `pysam` gate env lacks; the cost logic
# below is byte-identical and `_hp_run_length` is the SAME canonical fn (defined
# in hp_penalty.py) the shipped arbiter calls.  Faithful to the real incumbent.
def _cigar_hp_edit_distance(read, genome, penalty_table) -> float:
    if read.cigartuples is None:
        return 0.0
    chrom = read.reference_name
    genome_seq = (genome or {}).get(chrom) if genome else None
    query_seq = read.query_sequence
    ref_pos = read.reference_start
    q_pos = 0
    total = 0.0
    for op, length in read.cigartuples:
        if op == 7:                # =
            ref_pos += length; q_pos += length
        elif op == 8:              # X
            total += length * 1.0; ref_pos += length; q_pos += length
        elif op == 0:              # M (compare bases)
            if genome_seq and query_seq:
                ref_chunk = genome_seq[ref_pos:ref_pos + length].upper()
                q_chunk = query_seq[q_pos:q_pos + length].upper()
                total += sum(r != q for r, q in zip(ref_chunk, q_chunk))
            ref_pos += length; q_pos += length
        elif op == 2:              # D
            for i in range(length):
                rp = ref_pos + i
                if penalty_table is not None and genome_seq and rp < len(genome_seq):
                    hp = _hp_run_length(genome_seq, rp)
                    total += penalty_table.del_cost(hp, genome_seq[rp])
                else:
                    total += 1.0
            ref_pos += length
        elif op == 1:              # I
            if penalty_table is not None and genome_seq and ref_pos < len(genome_seq):
                hp = _hp_run_length(genome_seq, ref_pos)
                total += length * penalty_table.ins_cost(hp, genome_seq[ref_pos])
            else:
                total += length * 1.25
            q_pos += length
        elif op == 3:              # N — free
            ref_pos += length
        elif op == 4:              # S
            total += length * 1.0; q_pos += length
        elif op == 5:              # H
            total += length * 1.0
    return total


def _cigar_aligned_bases(read) -> int:
    if read.cigartuples is None:
        return 0
    return sum(length for op, length in read.cigartuples if op in (0, 7, 8))


def load_fastq(path):
    seqs = {}
    with pysam.FastxFile(path) as fq:
        for e in fq:
            seqs[e.name] = e.sequence
    return seqs


def _read_is_position_exact(read, rt, genome) -> bool:
    """Replicate _score_read's per-read position-exact rule (indels + junctions).

    A member alignment is correct for this read iff:
      * every truth indel's net (D-I) inside its ambiguity span equals the truth
        net AND no indel falls OUTSIDE every truth span (has_unexplained), AND
      * for junction strata: every truth junction is matched (ambiguity-aware via
        normalize_junction) and the member calls NO false-positive junction.
    """
    cig = read.cigartuples
    if cig is None:
        return False
    rstart = read.reference_start
    chrom = read.reference_name
    gseq = genome.get(chrom, "")

    # ---- Junctions (ambiguity-aware) ----
    truth_set = {(j.intron_start, j.intron_end) for j in rt.junctions}
    called = extract_junctions(rstart, cig)
    matched = set()
    for (cs, ce) in called:
        ns, ne = normalize_junction(cs, ce, gseq) if gseq else (cs, ce)
        if (ns, ne) in truth_set:
            matched.add((ns, ne))
        else:
            return False  # any FP junction => not exact
    if truth_set - matched:
        return False      # any missed truth junction => not exact

    # ---- Indels (the framing metric) ----
    truth_spans = [(ind.eq_start, ind.eq_end) for ind in rt.indels]
    has_unexplained = False
    for (ipos, ilen, ikind) in all_indel_positions(rstart, cig):
        if ikind == 1:
            covered = any(s <= ipos <= e for (s, e) in truth_spans)
        else:
            covered = any(s <= ipos < e for (s, e) in truth_spans)
        if not covered:
            has_unexplained = True
            break
    if has_unexplained:
        return False
    for ind in rt.indels:
        truth_net = ind.length if ind.kind == IndelKind.DEL else -ind.length
        in_span, _out = net_indel_in_span(rstart, cig, ind.eq_start, ind.eq_end)
        if in_span != truth_net:
            return False
    return True


def build_member_bam(read_seqs, truth_subset, genome, ref_fa, out_bam,
                     penalty_table, lam):
    """DP-arm member: align each read to its single-contig ref (chrom_ref on so
    homo_mask/homo_mismatch are active in every arm — the matched-baseline rule)."""
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


def run_minimap2(ref_fa, reads_fq, out_bam):
    import subprocess
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


def scan_member(bam_path, truth_subset, genome, penalty_table):
    """Per-read: {read_id: (hp_edit_distance, span, is_position_exact)}."""
    out = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            rid = read.query_name
            rt = truth_subset.get(rid)
            if rt is None:
                continue
            ed = _cigar_hp_edit_distance(read, genome, penalty_table)
            span = _cigar_aligned_bases(read)
            exact = _read_is_position_exact(read, rt, genome)
            out[rid] = (ed, span, exact)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/c3_hr")
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--penalty-table", required=True)
    ap.add_argument("--str-table", default=None)
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--min-count", type=int, default=100)
    ap.add_argument("--with-mm2", action="store_true",
                    help="add a real minimap2 member (external hold-out view)")
    ap.add_argument("--all-splits", action="store_true",
                    help="score ALL reads, not just TEST (more power for the probe)")
    ap.add_argument("--strata", default="HP,HP_HARD,STR",
                    help="comma list of strata to score. Default = the indel strata "
                         "where the exon-block DP arm is the valid tool (sub-vs-indel "
                         "arbitration). Junction/paralog strata need splice aligners, "
                         "not align_exon_block_global, so are excluded by default.")
    args = ap.parse_args()
    want_strata = {s.strip() for s in args.strata.split(",") if s.strip()}
    os.makedirs(args.out, exist_ok=True)

    print(f"[c3hr] generating corpus reps={args.reps} ...", file=sys.stderr)
    info = controlled.generate_corpus(args.out, reps=args.reps, seed=7)
    genome = load_genome(info["ref_fa"])
    ref_fa = info["ref_fa"]
    read_seqs = load_fastq(info["reads_fastq"])
    truth = read_truth_table(info["truth_tsv"])

    law_table = HpPenaltyTable.from_tsv(args.penalty_table, args.str_table,
                                        min_count=args.min_count)

    if args.all_splits:
        keep = truth
    else:
        keep = [t for t in truth if t.split == SplitTag.TEST]
    keep = [t for t in keep if t.stratum in want_strata]
    truth_subset = {t.read_id: t for t in keep}
    print(f"[c3hr] scoring {len(truth_subset)} reads "
          f"({'ALL splits' if args.all_splits else 'TEST split'})", file=sys.stderr)

    # ---- build members ----
    members = {}
    fb = os.path.join(args.out, "member_flat.bam")
    lb = os.path.join(args.out, "member_law.bam")
    build_member_bam(read_seqs, truth_subset, genome, ref_fa, fb, None, 1.0)
    build_member_bam(read_seqs, truth_subset, genome, ref_fa, lb, law_table, args.lam)
    members["flat"] = scan_member(fb, truth_subset, genome, law_table)
    members["law"] = scan_member(lb, truth_subset, genome, law_table)
    if args.with_mm2:
        mb = os.path.join(args.out, "member_mm2.bam")
        try:
            run_minimap2(ref_fa, info["reads_fastq"], mb)
            members["mm2"] = scan_member(mb, truth_subset, genome, law_table)
        except Exception as exc:
            print(f"[c3hr] minimap2 skipped: {exc}", file=sys.stderr)

    member_names = list(members)
    print(f"[c3hr] members: {member_names}", file=sys.stderr)

    # ---- per-read arbitration ----
    # shipped arbiter pick = argmin(hp_edit_distance), tiebreak span DESC.
    # Per stratum accumulate: ceiling (>=1 member exact), arbiter_correct,
    # headroom split into mis-order (strict argmin wrong, a higher-ed member
    # right) vs tie-arbitrary (>=1 member TIED at min-ed right, >=1 tied wrong).
    def _label(t):
        # break out boundary_sub (C1's discriminating cell, law only 0.55) from
        # the HP_HARD aggregate so the indel null is airtight per sub-case.
        if t.stratum == "HP_HARD":
            return "HP_HARD-bsub" if "_boundary_sub_" in t.read_id else "HP_HARD-noisy"
        return t.stratum

    agg = defaultdict(lambda: dict(n=0, ceiling=0, arb_correct=0,
                                   headroom=0, misorder=0, tie_arb=0,
                                   tie_lucky=0, disagree=0, hr_in_disagree=0))

    EPS = 1e-9
    for rid, t in truth_subset.items():
        present = [m for m in member_names if rid in members[m]]
        if not present:
            continue
        vals = {m: members[m][rid] for m in present}  # (ed, span, exact)
        s = agg[_label(t)]
        s["n"] += 1
        any_exact = any(v[2] for v in vals.values())
        if any_exact:
            s["ceiling"] += 1
        # genuine arbitration is only EXERCISED when members differ in correctness
        # (a headroom of 0 over reads where members never disagree is a weak null).
        n_exact = sum(1 for v in vals.values() if v[2])
        is_disagree = (len(present) > 1 and 0 < n_exact < len(present))
        if is_disagree:
            s["disagree"] += 1

        # shipped arbiter: min ed, then max span
        min_ed = min(v[0] for v in vals.values())
        tied = [m for m in present if vals[m][0] <= min_ed + EPS]
        # span tiebreak -> deterministic winner
        winner = max(tied, key=lambda m: vals[m][1])
        if vals[winner][2]:
            s["arb_correct"] += 1

        if any_exact and not vals[winner][2]:
            s["headroom"] += 1
            if is_disagree:
                s["hr_in_disagree"] += 1
            tie_has_correct = any(vals[m][2] for m in tied)
            if len(tied) > 1 and tie_has_correct:
                # the span tiebreak picked a wrong member but a TIED member is
                # correct -> arbitration could recover it for free (tie-arbitrary)
                s["tie_arb"] += 1
            else:
                # the strict argmin is wrong and a STRICTLY higher-ed member is
                # correct -> the scalar mis-ORDERS truth (the coherence defect)
                s["misorder"] += 1
        # diagnostic: arbiter got it right only via the span tiebreak luck
        if vals[winner][2] and len(tied) > 1 and not all(vals[m][2] for m in tied):
            s["tie_lucky"] += 1

    # ---- report ----
    print("\n================ C3 ARBITRATION HEADROOM (LLR-free) ================")
    print(f"panel={member_names}  split={'ALL' if args.all_splits else 'TEST'}  reps={args.reps}")
    print(f"{'stratum':16s} {'n':>5s} {'ceiling':>8s} {'arbiter':>8s} "
          f"{'HEADROOM':>9s} {'misord':>7s} {'tie_arb':>7s} {'disagr':>7s} {'hr|dis':>7s}")
    tot = dict(n=0, ceiling=0, arb_correct=0, headroom=0, misorder=0, tie_arb=0,
               disagree=0, hr_in_disagree=0)
    for st in sorted(agg):
        s = agg[st]
        n = max(1, s["n"])
        nd = max(1, s["disagree"])
        for k in tot:
            tot[k] += s[k]
        print(f"{st:16s} {s['n']:5d} {s['ceiling']/n:8.3f} {s['arb_correct']/n:8.3f} "
              f"{s['headroom']/n:9.3f} {s['misorder']/n:7.3f} {s['tie_arb']/n:7.3f} "
              f"{s['disagree']/n:7.3f} {s['hr_in_disagree']/nd:7.3f}")
    N = max(1, tot["n"])
    ND = max(1, tot["disagree"])
    print("-" * 86)
    print(f"{'TOTAL':16s} {tot['n']:5d} {tot['ceiling']/N:8.3f} {tot['arb_correct']/N:8.3f} "
          f"{tot['headroom']/N:9.3f} {tot['misorder']/N:7.3f} {tot['tie_arb']/N:7.3f} "
          f"{tot['disagree']/N:7.3f} {tot['hr_in_disagree']/ND:7.3f}")

    print("\n---- READING ----")
    print("ceiling   = freq a member IS position-exact (the recoverable universe)")
    print("arbiter   = freq the shipped hp_ed arbiter picks a correct member")
    print("HEADROOM  = ceiling - arbiter = reads recoverable by RE-RANKING the panel")
    print("  misord  = strict argmin(hp_ed) wrong, a higher-ed member correct (coherence defect)")
    print("  tie_arb = hp_ed TIES a correct+wrong member; span tiebreak picked wrong")
    print("disagr    = freq members genuinely DISAGREE in correctness (arbitration exercised)")
    print("hr|dis    = HEADROOM restricted to disagreement reads (the STRONG null: a 0 here")
    print("            means even where members disagree, the arbiter already picks truth)")
    print("\n---- VERDICT (pre-committed) ----")
    hr = tot["headroom"] / N
    if hr < 0.01:
        print(f"  HEADROOM={hr:.3f} ~ 0  => arbiter AT CEILING on recoverable reads.")
        print("  => C3-as-accuracy REFUTED on this corpus (C2 outcome). Document + stop;")
        print("     do NOT build an LLR decoder to chase a phantom gap. The artifact-replay")
        print("     (re-weight -> integer flips, LLR invariant) remains a valid CI fence only.")
    else:
        print(f"  HEADROOM={hr:.3f} > 0  => arbitration gap EXISTS. NEXT: build the coherent")
        print("     -logP LLR and test (matched-arm, TEST split) whether it recovers these")
        print("     reads vs truth, with the zero-evidence guard (truth must be calibrated-")
        print("     distinguishable) and the clean-subset safety control. mis-order reads are")
        print("     the genuine target; tie_arb reads are recoverable by ANY tiebreak (weaker).")


if __name__ == "__main__":
    main()
