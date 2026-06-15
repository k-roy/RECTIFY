#!/usr/bin/env python
"""Round-2 cDNA Phase-0 — Steps E + F + G: weak subset, lift Round-2 cDNA alignments,
score against the CORRECTED genome winner (+ corrected uLTRA), verdict + falsifier gates.

Design (advisor-shaped):
 * Genome side is the CORRECTED per-aligner BAMs (5'/soft-clip/walkback rescue already
   applied) — scoring against RAW genome would be a false-GO generator.
 * cDNA side is the RAW lifted alignment (no walkback) — CONSERVATIVE: any pad/3' soft-clip
   stays as a penalty on the cDNA, making wins harder, so a GO here is robust.
 * Verdict = score_phase0.run_verdict (win-guard, BLOCKER-1 raw-anchor gate, uLTRA-specific,
   trivial-win leak). PLUS two explicit falsifier gates the harness go() does not encode:
     - novelty-suppression: every cDNA win's transferred junctions clear the k>=10 raw
       perfect-anchor gate (already enforced) AND the win introduces no NET new mismatch
       vs the genome read (hp_ed_cdna <= hp_ed_genome by construction; logged per win);
     - CPA-regression: the cDNA's 3' genomic end must not move materially vs the genome call.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
from collections import defaultdict

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from liftover import lift, assert_orientation_invariant, N, S  # noqa: E402
from cdna_models import load_models  # noqa: E402
from score_phase0 import Candidate, run_verdict, decide  # noqa: E402

_REPO_GUESS = None  # rectify importable via PYTHONPATH set by caller
from rectify.core.consensus.corrected_consensus import (  # noqa: E402
    _cigar_hp_edit_distance, _cigar_aligned_bases, _cigar_min_junction_anchor, _NO_JUNCTION_ANCHOR,
)
from rectify.core.splice.hp_penalty import HpPenaltyTable  # noqa: E402


def softclip_bp(read):
    return sum(l for op, l in (read.cigartuples or []) if op == S)


def n_junctions(read):
    return sum(1 for op, l in (read.cigartuples or []) if op == N)


def three_prime_end(read_like, is_reverse):
    """Genomic 3' end of the RNA: reference_start for minus, reference_end-1 for plus."""
    ref = read_like.reference_start
    span = sum(l for op, l in read_like.cigartuples if op in (0, 2, 3, 7, 8))
    return ref if is_reverse else ref + span - 1


def build_genome_candidates(corrected_bams, genome, ptab):
    """Return per-read dict: {read_id: {aligner: Candidate(+meta)}}."""
    per_read = defaultdict(dict)
    meta = defaultdict(dict)  # read_id -> aligner -> (is_reverse, three_prime, ref_name)
    for aligner, path in corrected_bams.items():
        bam = pysam.AlignmentFile(path, "rb")
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            rid = read.query_name
            cand = Candidate.from_read(
                read, genome, ptab, aligner=aligner, read_id=rid,
                five_prime_rescued=0, genome_softclip_bp=softclip_bp(read),
            )
            per_read[rid][aligner] = cand
            meta[rid][aligner] = dict(is_reverse=read.is_reverse,
                                      tp=three_prime_end(read, read.is_reverse),
                                      ref=read.reference_name)
    return per_read, meta


def genome_winner(cands: dict):
    """Best genome candidate: min hp_ed, tie-break max aligned_bases."""
    return min(cands.values(), key=lambda c: (c.hp_ed, -c.aligned_bases))


def build_cdna_candidates(round2_bams, models, genome, ptab):
    """Per read, lift every cDNA alignment and keep the lowest-genome-HP-ED lifted candidate.
    Returns {read_id: (Candidate, meta)} where meta has is_reverse, tp(3'), n_junc, lifted."""
    best = {}
    invariant_fail = 0
    lift_fail = 0
    for aligner, path in round2_bams.items():
        bam = pysam.AlignmentFile(path, "rb")
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_reverse:
                # RNA-sense DRS reads align forward to the RNA-sense cDNA; a reverse hit
                # (e.g. a spurious mapPacBio antisense placement) breaks the lift orientation
                # contract, so drop it rather than mis-lift.
                continue
            cid = read.reference_name
            if cid not in models:
                continue
            m = models[cid]
            if read.query_sequence is None:
                continue
            try:
                lifted = lift(m, read.cigartuples, read.reference_start, read.query_sequence)
            except Exception:
                lift_fail += 1
                continue
            f_ok, l_ok = assert_orientation_invariant(lifted, genome)
            if not (f_ok and l_ok):
                invariant_fail += 1
                continue
            cand = Candidate.from_read(
                lifted, genome, ptab, aligner="cdna_round2", read_id=read.query_name,
                five_prime_rescued=0, genome_softclip_bp=0,
            )
            tp = three_prime_end(lifted, lifted.is_reverse)
            rid = read.query_name
            if rid not in best or cand.hp_ed < best[rid][0].hp_ed:
                best[rid] = (cand, dict(is_reverse=lifted.is_reverse, tp=tp,
                                        n_junc=n_junctions(lifted), cid=cid,
                                        ref_start=lifted.reference_start))
    return best, invariant_fail, lift_fail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blockmaps", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--corrected-bam", action="append", required=True,
                    help="NAME=path for a corrected genome BAM (repeat)")
    ap.add_argument("--round2-bam", action="append", required=True,
                    help="NAME=path for a Round-2 cDNA BAM (repeat)")
    ap.add_argument("--penalty-scores", required=True)
    ap.add_argument("--str-penalty", default=None)
    ap.add_argument("--min-anchor", type=int, default=10)
    ap.add_argument("--epsilon", type=float, default=1.0)
    ap.add_argument("--weak-clip", type=int, default=20)
    ap.add_argument("--cpa-tol", type=int, default=12, help="bp tolerance for CPA regression")
    ap.add_argument("--out", required=True, help="output dir for verdict + per-win TSV")
    args = ap.parse_args()

    cor = dict(kv.split("=", 1) for kv in args.corrected_bam)
    r2 = dict(kv.split("=", 1) for kv in args.round2_bam)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    models = load_models(args.blockmaps)
    fa = pysam.FastaFile(args.genome)
    genome = {c: fa.fetch(c) for c in fa.references}
    ptab = HpPenaltyTable.from_tsv(args.penalty_scores, str_path=args.str_penalty)

    print("[score] building corrected genome candidates ...", flush=True)
    gcands, gmeta = build_genome_candidates(cor, genome, ptab)
    print(f"[score]   {len(gcands)} reads with >=1 corrected genome alignment")

    print("[score] lifting Round-2 cDNA alignments ...", flush=True)
    ccands, inv_fail, lift_fail = build_cdna_candidates(r2, models, genome, ptab)
    print(f"[score]   {len(ccands)} reads with a lifted cDNA candidate "
          f"(orientation-invariant failures: {inv_fail}, lift errors: {lift_fail})")

    # eligible = read has BOTH a genome winner and a lifted cDNA candidate
    triples = []
    rows = []
    eligible = 0
    weak = 0
    ultra_absent = 0
    for rid, cdict in gcands.items():
        if rid not in ccands:
            continue
        gwin = genome_winner(cdict)
        cdna, cmeta = ccands[rid]
        ultra = cdict.get("uLTRA")
        if ultra is None:
            # uLTRA was run on every read; no record => uLTRA failed to place this read,
            # so a cDNA win here genuinely beats uLTRA. Use an inf-HP-ED sentinel.
            ultra = Candidate(read_id=rid, aligner="uLTRA_absent", hp_ed=float("inf"),
                              aligned_bases=0, min_junction_anchor=_NO_JUNCTION_ANCHOR)
            ultra_absent += 1
        eligible += 1
        # Step E weak label (cost-control population, for reporting)
        anchors = gwin.min_junction_anchor
        failed_gate = (anchors < _NO_JUNCTION_ANCHOR and anchors < args.min_anchor)
        hp_disagree = max(c.hp_ed for c in cdict.values()) - min(c.hp_ed for c in cdict.values())
        spliced = gwin.n_junctions >= 1 or cmeta["n_junc"] >= 1 or gwin.genome_softclip_bp >= args.weak_clip
        is_weak = spliced and (gwin.genome_softclip_bp >= args.weak_clip or failed_gate or hp_disagree >= 3.0)
        if is_weak:
            weak += 1
        triples.append((cdna, gwin, ultra))
        rows.append(dict(rid=rid, gwin=gwin, cdna=cdna, ultra=ultra, cmeta=cmeta,
                         gmeta=gmeta[rid].get(gwin.aligner, {}), is_weak=is_weak))

    print(f"[score] eligible reads (genome+cDNA both): {eligible}; weak subset: {weak}")

    v = run_verdict(triples, min_junction_anchor_bp=args.min_anchor, epsilon=args.epsilon)

    # ---- per-win falsifier gates (Step G extra) ----
    win_rows = []
    cpa_regress = 0
    novelty_suspect = 0
    weak_wins = 0
    for row, d in zip(rows, v.decisions):
        if not d.cdna_wins:
            continue
        gwin = row["gwin"]; cdna = row["cdna"]; cmeta = row["cmeta"]; gm = row["gmeta"]
        if row["is_weak"]:
            weak_wins += 1
        # CPA regression: cDNA 3' end vs genome winner 3' end
        cpa_shift = None
        if gm:
            cpa_shift = cmeta["tp"] - gm.get("tp", cmeta["tp"])
            if abs(cpa_shift) > args.cpa_tol:
                cpa_regress += 1
        # novelty-suppression: a win that improved HP-ED while ADDING junctions must carry
        # real perfect anchors (already gated >=min_anchor). Flag if cdna added junctions
        # but the genome winner had a LOWER soft-clip than expected (i.e. nothing to recover).
        added_junc = cmeta["n_junc"] - gwin.n_junctions
        if added_junc > 0 and gwin.genome_softclip_bp < 5 and cdna.min_junction_anchor < args.min_anchor:
            novelty_suspect += 1
        win_rows.append(dict(rid=row["rid"], cause=d.cause, beats_ultra=d.beats_ultra,
                             hp_ed_genome=d.hp_ed_genome, hp_ed_cdna=d.hp_ed_cdna,
                             cdna_anchor=cdna.min_junction_anchor, gwin_clip=gwin.genome_softclip_bp,
                             cdna_njunc=cmeta["n_junc"], gwin_njunc=gwin.n_junctions,
                             cpa_shift=cpa_shift, cid=cmeta["cid"], is_weak=row["is_weak"]))

    # write outputs
    with open(out / "wins.tsv", "w") as fh:
        cols = ["rid", "cause", "beats_ultra", "hp_ed_genome", "hp_ed_cdna", "cdna_anchor",
                "gwin_clip", "cdna_njunc", "gwin_njunc", "cpa_shift", "is_weak", "cid"]
        fh.write("\t".join(cols) + "\n")
        for w in win_rows:
            fh.write("\t".join(str(w.get(c, "")) for c in cols) + "\n")

    summary = v.summary()
    extra = (
        f"  absolute counts            : eligible={v.n_eligible} wins={v.n_cdna_wins} "
        f"weak_wins={weak_wins} beat_ultra={v.n_wins_beating_ultra} no_cause={v.n_wins_no_cause}\n"
        f"  uLTRA-absent eligible reads: {ultra_absent} (cDNA wins among these count as beating uLTRA)\n"
        f"  [gate] CPA-regression wins (>|{args.cpa_tol}|bp 3' shift): {cpa_regress}  [falsifier if >0]\n"
        f"  [gate] novelty-suppression suspects: {novelty_suspect}  [falsifier if >0]\n"
    )
    go_base = v.go()
    go_final = go_base and cpa_regress == 0 and novelty_suspect == 0
    verdict_txt = summary + extra + f"\n[Phase-0 FINAL verdict (incl. extra gates): {'GO' if go_final else 'NO-GO'}]\n"
    print(verdict_txt)
    (out / "verdict.txt").write_text(verdict_txt)
    with open(out / "verdict.json", "w") as fh:
        json.dump(dict(go_base=go_base, go_final=go_final, n_eligible=v.n_eligible,
                       n_cdna_wins=v.n_cdna_wins, weak_wins=weak_wins,
                       n_wins_beating_ultra=v.n_wins_beating_ultra,
                       trivial_win_leak=v.trivial_win_leak,
                       ultra_specific_win_frac=v.ultra_specific_win_frac,
                       net_hp_ed_reduction=v.net_hp_ed_reduction,
                       cpa_regress=cpa_regress, novelty_suspect=novelty_suspect,
                       n_eligible_weak=weak, ultra_absent=ultra_absent), fh, indent=1)


if __name__ == "__main__":
    main()
