#!/usr/bin/env python3
"""Discovery canonical-snap TIEBREAKER probe — the C3 residue / Discovery facet.

C3 refuted calibrated-LLR junction arbitration on BOTH arbiters. For the
``hp_edit_distance`` arbiter (merge_corrected_tsvs) the junction probe
(c3_junction_headroom.py) showed pickSnap=0. C3 then located the residual
canonical-snap bias precisely: it can only survive in the OTHER arbiter,
``select.py::select_best_alignment`` (the ``junction_score`` path), and there
ONLY in the ``_n_annotated``/``canonical_count`` TIEBREAKER (which fires solely
on an exact PRIMARY-score tie). C3's structural spot-check predicted the snap
loses on the PRIMARY ``junction_score`` because ``_count_junction_proximity_errors``
penalizes the snap-induced mismatch within 5bp of the junction → the tiebreaker
never decides → residual FDR = 0.

This probe CLOSES that loose end empirically by driving the ACTUAL shipped
``select_best_alignment(tiebreak='rectify')`` (pulled from the canonical
drs-validation-rebuild branch) on a minimal, conservative 2-member family per
JUNCTION_DISCOVERY read:

  * mm2   = real ``minimap2 -ax splice`` (the snap). On non-canonical discovery
            junctions it snaps the truth site onto a nearby canonical GT-AG.
  * truth = the corpus true_cigar placement (orthogonal de-novo placer;
            position-exact by construction).

The 2-member family is the construction MOST favorable to the snap winning: with
only 2 members the ``_count_3prime_agreement`` tiebreaker is degenerate (both
agree → tie), forcing ``_n_annotated`` then ``canonical_count`` to decide — and
on a ``jd_non_novel`` read the snap has canonical_count=1 vs truth's 0. So a
measured pickSnap=0 here is a STRONG refute; pickSnap>0 is the genuine residual.

Measured (per stratum):
  n                  reads with both members placed
  mm2_snap           freq mm2 is NOT position-exact (snapped/missed)
  snap&truth         freq mm2 snapped AND truth member is position-exact (the setup)
  tie_on_score       OF snap&truth: freq snap and truth TIE on PRIMARY junction_score
  tiebreak_picks_snap OF the ties: freq the rectify TIEBREAKER picked the SNAP over truth
                     (= the residual canonical-snap tiebreaker FDR — C3 predicts ~0)
  winner_snap        OF snap&truth: freq select_best_alignment's winner is the snap
                     (by score OR tiebreak — the total harm, if any)

advisor gate (measurement validity): the annotated set is reconstructed from the
ANNOTATED-stratum truth rows (klass==ANNOTATED). We HARD-ASSERT that on the
ANNOTATED stratum the truth member yields n_annotated>=1 — otherwise the
annotated set is mis-keyed (coord/strand) and _n_annotated silently collapses to
0, under-measuring the tiebreaker.

Usage:
  python scripts/benchmark/discovery_tiebreak_probe.py --out /tmp/disc_tb --reps 40
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.scorer import load_genome  # noqa: E402
from rectify.core.benchmark.truth_schema import (  # noqa: E402
    read_truth_table, SplitTag, JunctionClass,
)
from rectify.core.consensus.extract import extract_alignment_info  # noqa: E402
from rectify.core.consensus.select import select_best_alignment  # noqa: E402
from scripts.benchmark.sim import controlled  # noqa: E402
from scripts.benchmark.c3_headroom import (  # noqa: E402
    load_fastq, run_minimap2, _read_is_position_exact,
)
from scripts.benchmark.c3_junction_headroom import build_truth_member_bam  # noqa: E402
import pysam  # noqa: E402


def _load_bam_by_id(bam_path):
    """{read_id: primary AlignedSegment} (skip unmapped/secondary/supplementary)."""
    out = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            out.setdefault(r.query_name, r)
    return out


def _build_annotated_set(truth_rows):
    """Faithful corpus annotation: every truth junction tagged klass==ANNOTATED,
    keyed (chrom, intron_start, intron_end, strand) — the exact key select.py uses."""
    annot = set()
    for t in truth_rows:
        for j in t.junctions:
            if j.klass == JunctionClass.ANNOTATED:
                annot.add((t.chrom, j.intron_start, j.intron_end, t.strand))
    return annot


def _n_annot(info, annotated):
    return sum(1 for (js, je) in info.junctions
               if (info.chrom, js, je, info.strand) in annotated)


def _ref_edit_distance(read, genome):
    """Aligned edit-distance to reference (I + D + X), encoding-agnostic.

    Mirrors the standard NM aux tag for the ALIGNED portion (soft-clips excluded,
    since junction_score already accounts for clips separately). Works whether the
    aligner emits M (truth member), =/X (minimap2 -ax splice), or mixed. A snap that
    shifts the junction onto a canonical site must absorb the displaced bases as an
    insertion (1I/3I here) -> that misplacement cost surfaces as a strictly higher
    edit distance than the position-exact truth member."""
    cig = read.cigartuples
    if cig is None:
        return 0
    ref = genome.get(read.reference_name, "")
    q = read.query_sequence or ""
    ed = 0
    rpos = read.reference_start
    qpos = 0
    for op, n in cig:
        if op == 1:             # I
            ed += n; qpos += n
        elif op == 2:           # D
            ed += n; rpos += n
        elif op == 8:           # X
            ed += n; qpos += n; rpos += n
        elif op == 7:           # =
            qpos += n; rpos += n
        elif op == 0:           # M: compare bytes vs reference
            for i in range(n):
                if rpos + i < len(ref) and qpos + i < len(q) \
                        and ref[rpos + i].upper() != q[qpos + i].upper():
                    ed += 1
            qpos += n; rpos += n
        elif op == 4:           # S
            qpos += n
        elif op == 3:           # N
            rpos += n
        # H: nothing
    return ed


def _proxerr_fixed_attribution(read, genome, w=5):
    """`_count_junction_proximity_errors` with ONE surgical change: a post-`N`
    insertion is attributed to `intron_end` (first EXON base, in the exon-side
    window) instead of `intron_end-1` (last intron base, the current blind spot).

    Tests advisor/panel-iii option (a): does fixing the attribution make a
    boundary-snapping insertion incur a PRIMARY-score penalty (so the snap loses on
    junction_score and never reaches the tiebreaker)? Everything else mirrors the
    shipped function exactly. Returns weighted error count."""
    from rectify.core.consensus.scoring import _is_homopolymer_position
    cig = read.cigartuples
    if not cig or not read.query_sequence:
        return 0.0
    chrom = read.reference_name
    if not chrom or chrom not in genome:
        return 0.0
    junctions = []
    rp = read.reference_start
    for op, n in cig:
        if op == 3:
            junctions.append((rp, rp + n)); rp += n
        elif op in (0, 2, 7, 8):
            rp += n
    if not junctions:
        return 0.0
    ref = genome[chrom]; q = read.query_sequence
    prox = set()
    for js, je in junctions:
        for r in range(max(0, js - w), js):
            prox.add(r)
        for r in range(je, min(len(ref), je + w)):
            prox.add(r)
    errors = 0.0
    ref_pos = read.reference_start; qpos = 0; prev_rp = None
    for op, n in cig:
        if op == 5:
            pass
        elif op == 4:
            qpos += n
        elif op in (0, 7, 8):
            for i in range(n):
                r = ref_pos + i
                if r in prox and r < len(ref):
                    if ref[r].upper() != 'N' and q[qpos + i].upper() != ref[r].upper():
                        errors += 0.5 if _is_homopolymer_position(ref, r) else 1.0
                prev_rp = r
            qpos += n; ref_pos += n
        elif op == 1:
            if prev_rp is not None and prev_rp in prox:
                errors += 0.5 if _is_homopolymer_position(ref, prev_rp) else 1.0
            qpos += n
        elif op == 2:
            for i in range(n):
                r = ref_pos + i
                if r in prox:
                    errors += 0.5 if _is_homopolymer_position(ref, r) else 1.0
                prev_rp = r
            ref_pos += n
        elif op == 3:
            prev_rp = ref_pos + n          # FIX: intron_end (first exon base), not -1
            ref_pos += n
    return errors


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/disc_tb")
    ap.add_argument("--reps", type=int, default=40)
    ap.add_argument("--strata",
                    default="JUNCTION_DISCOVERY,ANNOTATED,NIC,JUNCTION_AMB",
                    help="junction strata to score")
    ap.add_argument("--all-splits", action="store_true", default=True)
    ap.add_argument("--evidence", default=None,
                    help="optional path to dump per-read evidence TSV")
    args = ap.parse_args()
    want = {s.strip() for s in args.strata.split(",") if s.strip()}
    os.makedirs(args.out, exist_ok=True)

    print(f"[disc_tb] generating corpus reps={args.reps} ...", file=sys.stderr)
    info = controlled.generate_corpus(args.out, reps=args.reps, seed=7)
    genome = load_genome(info["ref_fa"])
    ref_fa = info["ref_fa"]
    read_seqs = load_fastq(info["reads_fastq"])
    truth = read_truth_table(info["truth_tsv"])

    annotated = _build_annotated_set(truth)
    print(f"[disc_tb] annotated set: {len(annotated)} junctions", file=sys.stderr)

    keep = truth if args.all_splits else [t for t in truth if t.split == SplitTag.TEST]
    keep = [t for t in keep if t.stratum in want]
    truth_subset = {t.read_id: t for t in keep}
    print(f"[disc_tb] {len(truth_subset)} junction reads "
          f"(strata={sorted(want)})", file=sys.stderr)

    mb = os.path.join(args.out, "tb_member_mm2.bam")
    tb = os.path.join(args.out, "tb_member_truth.bam")
    run_minimap2(ref_fa, info["reads_fastq"], mb)
    n_truth = build_truth_member_bam(read_seqs, truth_subset, ref_fa, tb)
    print(f"[disc_tb] truth-site member: {n_truth} reads placed", file=sys.stderr)
    mm2_reads = _load_bam_by_id(mb)
    truth_reads = _load_bam_by_id(tb)

    agg = defaultdict(lambda: dict(n=0, mm2_snap=0, snap_and_truth=0,
                                   tie_on_score=0, tiebreak_picks_snap=0,
                                   winner_snap=0))
    # advisor validation gate accumulators
    annot_stratum_truth_nann = []
    ev_rows = []
    _abl_cache = []  # per-read (infos, exact_map, snapped_map, truth_exact_map)

    for rid, t in truth_subset.items():
        if rid not in mm2_reads or rid not in truth_reads:
            continue
        mm2_info = extract_alignment_info(mm2_reads[rid], "mm2", genome)
        truth_info = extract_alignment_info(truth_reads[rid], "truth", genome)
        aligns = {"mm2": mm2_info, "truth": truth_info}
        res = select_best_alignment(aligns, genome, annotated, tiebreak="rectify")

        mm2_score = mm2_info.junction_score
        truth_score = truth_info.junction_score
        tie = len(res.tied_aligners) == 2
        winner = res.best_aligner
        truth_exact = _read_is_position_exact(truth_reads[rid], t, genome)
        mm2_exact = _read_is_position_exact(mm2_reads[rid], t, genome)
        snapped = not mm2_exact
        mm2_nann = _n_annot(mm2_info, annotated)
        truth_nann = _n_annot(truth_info, annotated)
        mm2_info._nann = mm2_nann
        truth_info._nann = truth_nann
        mm2_info._ed = _ref_edit_distance(mm2_reads[rid], genome)
        truth_info._ed = _ref_edit_distance(truth_reads[rid], genome)
        _abl_cache.append((
            {"mm2": mm2_info, "truth": truth_info},
            {"mm2": mm2_exact, "truth": truth_exact},
            {"mm2": snapped, "truth": not truth_exact},
            {"truth": truth_exact, "mm2": mm2_exact},
        ))

        if t.stratum == "ANNOTATED":
            annot_stratum_truth_nann.append(truth_nann)

        s = agg[t.stratum]
        s["n"] += 1
        if snapped:
            s["mm2_snap"] += 1
        if snapped and truth_exact:
            s["snap_and_truth"] += 1
            if tie:
                s["tie_on_score"] += 1
                if winner == "mm2":
                    s["tiebreak_picks_snap"] += 1
            if winner == "mm2":
                s["winner_snap"] += 1

        ev_rows.append((rid, t.stratum, f"{mm2_score:.2f}", f"{truth_score:.2f}",
                        int(tie), mm2_nann, truth_nann,
                        mm2_info.canonical_count, truth_info.canonical_count,
                        winner, int(truth_exact), int(mm2_exact),
                        f"{mm2_info.junction_proximity_errors:.2f}",
                        f"{truth_info.junction_proximity_errors:.2f}",
                        mm2_info._ed, truth_info._ed))

    # ---- REWEIGHT ABLATION (does demoting/dropping canonical_count help?) ----
    # Recompute the winner per read under variant tiebreak tuples, using EXACTLY
    # select.py's components (3'-agreement, n_annotated, canonical_count). Ties
    # within max() resolve to first-in-order; we test BOTH member orders to expose
    # whether any pick is order-arbitrary (no truth-favoring signal) vs. robust.
    # Fitness = winner position-exact (recall) and winner-is-snap (FDR), vs TRUTH.
    def _agree(info, all_c3):
        pos = info.corrected_3prime
        return sum(1 for p in all_c3 if p == pos) if pos is not None else 0

    # ADDITIVE variants test whether a truth-favoring signal (alignment edit
    # distance to reference; lower=cleaner) the shipped tuple OMITS can redirect the
    # snap-vs-truth tie toward truth (steelman iii: the snap pays its junction-shift
    # as an insertion -> higher NM -> a -NM key breaks the tie toward truth).
    variants = {
        "V0_shipped (agree,nann,canon)":     lambda info, ag: (ag, info._nann, info.canonical_count),
        "V1_drop_canon (agree,nann)":        lambda info, ag: (ag, info._nann),
        "V2_only_agree (agree,)":            lambda info, ag: (ag,),
        "V3_add_NM (agree,nann,-ED,canon)":  lambda info, ag: (ag, info._nann, -info._ed, info.canonical_count),
        "V4_NM_first (agree,nann,-ED)":      lambda info, ag: (ag, info._nann, -info._ed),
    }
    abl = {v: dict(win_exact=0, win_snap_harm=0) for v in variants}
    abl_rev = {v: dict(win_exact=0, win_snap_harm=0) for v in variants}  # truth-first order
    for rec in _abl_cache:
        infos, exact_map, snapped_map, truth_exact_map = rec
        all_c3 = [i.corrected_3prime for i in infos.values() if i.corrected_3prime is not None]
        for v, keyfn in variants.items():
            for order, store in ((("mm2", "truth"), abl), (("truth", "mm2"), abl_rev)):
                names = [n for n in order if n in infos]
                # replicate select: max over tied-top-score members by variant key
                top = max(infos[n].junction_score for n in names)
                tied = [n for n in names if infos[n].junction_score == top]
                w = max(tied, key=lambda n: keyfn(infos[n], _agree(infos[n], all_c3)))
                if exact_map[w]:
                    store[v]["win_exact"] += 1
                if (not exact_map[w]) and w == "mm2" and snapped_map["mm2"] and truth_exact_map.get("truth"):
                    store[v]["win_snap_harm"] += 1

    # ---- advisor validation GATE ----
    n_annot_reads = len(annot_stratum_truth_nann)
    n_annot_ok = sum(1 for v in annot_stratum_truth_nann if v >= 1)
    gate_ok = n_annot_reads == 0 or n_annot_ok >= 1
    print(f"\n[GATE] ANNOTATED-stratum reads={n_annot_reads}, "
          f"truth n_annotated>=1 in {n_annot_ok} "
          f"({'OK' if n_annot_ok == n_annot_reads else 'PARTIAL'})", file=sys.stderr)
    if not gate_ok:
        print("[GATE] FAIL: annotated set mis-keyed (coord/strand) — _n_annotated "
              "collapses to 0; measurement INVALID. Aborting.", file=sys.stderr)
        sys.exit(3)

    if args.evidence:
        with open(args.evidence, "w") as fh:
            fh.write("read_id\tstratum\tmm2_score\ttruth_score\ttie\tmm2_nann\t"
                     "truth_nann\tmm2_canon\ttruth_canon\twinner\ttruth_exact\t"
                     "mm2_exact\tmm2_proxerr\ttruth_proxerr\tmm2_ed\ttruth_ed\n")
            for r in ev_rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        print(f"[disc_tb] evidence dumped: {args.evidence} ({len(ev_rows)} rows)",
              file=sys.stderr)

    print("\n===== DISCOVERY CANONICAL-SNAP TIEBREAKER PROBE "
          "(select_best_alignment, tiebreak='rectify') =====")
    print("members=[mm2(real snap), truth(orthogonal placer)]  conservative 2-member family")
    print(f"{'stratum':20s} {'n':>5s} {'mm2snap':>8s} {'snap&tru':>9s} "
          f"{'tie/score':>10s} {'tb->snap':>9s} {'win=snap':>9s}")
    tot = defaultdict(int)
    for st in sorted(agg):
        s = agg[st]
        n = max(1, s["n"]); snt = max(1, s["snap_and_truth"])
        for k, v in s.items():
            tot[k] += v
        print(f"{st:20s} {s['n']:5d} {s['mm2_snap']/n:8.3f} {s['snap_and_truth']/n:9.3f} "
              f"{s['tie_on_score']/snt:10.3f} {s['tiebreak_picks_snap']/snt:9.3f} "
              f"{s['winner_snap']/snt:9.3f}")
    N = max(1, tot["n"]); SNT = max(1, tot["snap_and_truth"])
    print("-" * 78)
    print(f"{'TOTAL':20s} {tot['n']:5d} {tot['mm2_snap']/N:8.3f} {tot['snap_and_truth']/N:9.3f} "
          f"{tot['tie_on_score']/SNT:10.3f} {tot['tiebreak_picks_snap']/SNT:9.3f} "
          f"{tot['winner_snap']/SNT:9.3f}")

    print("\n---- READING ----")
    print("tie/score = OF snap&truth reads, freq snap and truth TIE on PRIMARY junction_score")
    print("tb->snap  = OF snap&truth reads, freq the rectify TIEBREAKER picked the SNAP")
    print("            (= residual canonical-snap tiebreaker FDR; C3 predicts ~0)")
    print("win=snap  = OF snap&truth reads, freq the winner is the snap (score OR tiebreak)")
    print("\n---- REWEIGHT ABLATION (winner recall + harmful snap-FDR, all strata) ----")
    print(f"{'variant':32s} {'win_exact':>10s} {'snap_harm(mm2-1st)':>19s} "
          f"{'snap_harm(truth-1st)':>21s}")
    for v in variants:
        print(f"{v:32s} {abl[v]['win_exact']:10d} {abl[v]['win_snap_harm']:19d} "
              f"{abl_rev[v]['win_snap_harm']:21d}")
    print("  (win_exact: # reads whose winner is position-exact; higher=better recall.")
    print("   snap_harm: # reads where winner is the snap while a truth member was exact.")
    print("   Compare the two order columns: if a variant's snap_harm CHANGES with member")
    print("   order, that pick is ORDER-ARBITRARY — the variant has no truth-favoring signal.)")

    print("\n---- VERDICT (pre-committed) ----")
    tb_snap = tot["tiebreak_picks_snap"] / SNT
    win_snap = tot["winner_snap"] / SNT
    if tot["snap_and_truth"] == 0:
        print("  no snap&truth reads constructed — cannot measure (check minimap2 snapping)")
    elif tb_snap < 1e-9 and win_snap < 1e-9:
        print(f"  tb->snap={tb_snap:.3f}, win=snap={win_snap:.3f} = 0 => the canonical/annotated")
        print("  TIEBREAKER NEVER decides a snap-vs-truth contest: the snap loses on the PRIMARY")
        print("  junction_score (proximity-error penalty), so the tiebreaker is structurally")
        print("  unreachable. TIEBREAKER ALREADY HARMLESS — no reweight needed (C3 confirmed).")
    else:
        print(f"  tb->snap={tb_snap:.3f} win=snap={win_snap:.3f} > 0 => a REAL residual: the")
        print("  tiebreaker (or primary score) prefers the canonical snap over an available truth")
        print("  member. Prototype a reweight (demote canonical_count/_n_annotated) and re-measure")
        print("  snap-FDR + ANNOTATED/NIC recall + the GMAP fence before changing select.py.")


if __name__ == "__main__":
    main()
