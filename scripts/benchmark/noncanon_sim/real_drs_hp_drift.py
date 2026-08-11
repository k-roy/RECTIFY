#!/usr/bin/env python3
"""Real-DRS transfer test for the HP-drift guard (truth = SGD/gencode annotations).

Does the guard's drift-fix — proven on simulated undercalls — hold on REAL ONT DRS
reads, where the homopolymer undercalls are real?  Runs on the cluster where the
BY4742 DRS alignment + genome live.

Design (PI-approved): specificity (drift-fix) on real reads with annotation truth,
plus a do-no-harm check.  Discovery preservation is taken from the sim (zero cost)
and confirmed here only as "the guard is byte-identical away from homopolymers".

Pipeline:
  1. Load the yeast genome + annotated introns (from the GFF).
  2. Flag ANNOTATED junctions that ABUT a homopolymer — a run (>= min_run) reachable
     by a small drift (<= probe bp) of a boundary.  These are the at-risk junctions.
  3. Refine the real minimap2 BAM TWICE: arm-B (motif_blind, no guard) and
     arm-B+guard (motif_blind + hp_drift_margin).  Same pool, same reads.
  4. At the HP-abutting annotated junctions, measure the fraction of spanning reads
     whose refined junction MATCHES the annotation (normalize-aware) vs is drifted
     INTO the homopolymer, per arm.  Guard should raise the match fraction.
  5. DO-NO-HARM: (a) overall annotated-junction match rate must not drop; (b) the two
     arms may differ ONLY at HP-abutting junctions (byte-identical elsewhere).

Usage (on the cluster):
  python real_drs_hp_drift.py \
    --bam  /scratch/.../wt_by4742_rep1.minimap2.namesorted.bam \
    --genome <yeast.fsa> --gff <yeast.gff.gz> \
    --hp-drift-margin 3.0 --out real_drs_hp_drift.json
"""
from __future__ import annotations
import argparse, json, os, sys, tempfile
from collections import defaultdict

import pysam

_here = os.path.dirname(os.path.abspath(__file__))
_repo = os.path.abspath(os.path.join(_here, "..", "..", ".."))
if _repo not in sys.path:
    sys.path.insert(0, _repo)

from rectify.core.splice.junction_refiner import (
    _hp_run_across, refine_bam_junctions, build_junction_pool,
)
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.core.consensus.chimeric_consensus import normalize_junction
from rectify.utils.genome import register_genome_contigs_from_fasta, standardize_chrom_name


def load_genome(fa):
    g, name = {}, None
    opener = __import__("gzip").open if fa.endswith(".gz") else open
    for l in opener(fa, "rt"):
        if l.startswith(">"):
            name = l[1:].split()[0]; g[name] = []
        else:
            g[name].append(l.strip())
    return {standardize_chrom_name(k): "".join(v).upper() for k, v in g.items()}


def hp_abutting(genome, annot, min_run, probe):
    """annotated (chrom, s, e) whose donor or acceptor can drift <=probe bp INTO a run."""
    out = set()
    for t in annot:
        c = standardize_chrom_name(t[0]); s, e = int(t[1]), int(t[2])
        seq = genome.get(c, "")
        if not seq:
            continue
        hit = False
        for pos in (s, e):
            for k in range(1, probe + 1):
                if _hp_run_across(seq, pos + k, min_run) or _hp_run_across(seq, pos - k, min_run):
                    hit = True; break
            if hit:
                break
        if hit:
            out.add((c, s, e))
    return out


def read_junctions(bam):
    """read_id -> list of normalized (donor, acceptor) N-op introns (per primary read)."""
    out = defaultdict(list)
    b = pysam.AlignmentFile(bam, "rb")
    for r in b:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        ref = r.reference_start
        for op, ln in (r.cigartuples or []):
            if op == 3:
                out[(standardize_chrom_name(r.reference_name), r.query_name)].append((ref, ref + ln))
            if op in (0, 2, 3, 7, 8):
                ref += ln
    b.close()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--hp-drift-margin", type=float, default=3.0)
    ap.add_argument("--hp-drift-min-run", type=int, default=4)
    ap.add_argument("--probe", type=int, default=3, help="max drift bp to count a junction HP-abutting")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--out", default="real_drs_hp_drift.json")
    args = ap.parse_args()

    register_genome_contigs_from_fasta(args.genome)
    genome = load_genome(args.genome)
    annot = list(load_annotated_junctions(args.gff))
    annot_norm = {(standardize_chrom_name(t[0]),) + normalize_junction(int(t[1]), int(t[2]),
                  genome.get(standardize_chrom_name(t[0]), "")) for t in annot
                  if genome.get(standardize_chrom_name(t[0]), "")}
    hp_junc = hp_abutting(genome, annot, args.hp_drift_min_run, args.probe)
    hp_junc_norm = {(c,) + normalize_junction(s, e, genome[c]) for (c, s, e) in hp_junc}
    sys.stderr.write(f"[real_drs] {len(annot_norm)} annotated junctions; "
                     f"{len(hp_junc)} HP-abutting (min_run={args.hp_drift_min_run}, probe={args.probe})\n")

    pool, annot_set = build_junction_pool([args.bam], annot)
    arms = {}
    for name, hpd in (("B", 0.0), ("Bguard", args.hp_drift_margin)):
        out_bam = os.path.join(args.outdir, f"real_arm_{name}.bam")
        refine_bam_junctions(
            input_bam=args.bam, output_bam=out_bam, aligner_bams=[args.bam],
            annotated_junctions=annot, genome=genome, penalty_table_path=None,
            prebuilt_junction_pool=pool, prebuilt_annotated_set=annot_set,
            sort_and_index=True, n_workers=args.workers,
            motif_blind=True, hp_drift_margin=hpd, hp_drift_min_run=args.hp_drift_min_run)
        arms[name] = read_junctions(out_bam)

    # --- proximity index: is a junction boundary AT an HP-abutting annotated site? ---
    # (a read that DRIFTED sits a few bp off the annotation, so membership is by proximity,
    #  not exact match — this is what the old `match_stats` got wrong by dividing over ALL reads.)
    from bisect import bisect_left
    hp_don = defaultdict(list); hp_acc = defaultdict(list)
    for (c, s, e) in hp_junc:
        hp_don[c].append(s); hp_acc[c].append(e)
    for d in (hp_don, hp_acc):
        for c in d:
            d[c].sort()

    def _near(sorted_pos, x, W):
        if not sorted_pos:
            return False
        i = bisect_left(sorted_pos, x)
        return any(0 <= j < len(sorted_pos) and abs(sorted_pos[j] - x) <= W
                   for j in (i - 1, i))

    W_MEMBER = 15  # generous: catches drifted placements as "at" the HP-abutting junction

    def at_hp(c, introns):
        return any(_near(hp_don.get(c, []), s, W_MEMBER) or _near(hp_acc.get(c, []), e, W_MEMBER)
                   for (s, e) in introns)

    def matches(c, introns):
        return any((c,) + normalize_junction(s, e, genome.get(c, "")) in annot_norm
                   for (s, e) in introns)

    # per-arm: overall annotated-match rate (per junction) + the match rate RESTRICTED to
    # junctions at HP-abutting sites (the corrected metric).
    arm_stats = {}
    for name in ("B", "Bguard"):
        tot = hit = hp_tot = hp_hit = 0
        for (c, _rid), introns in arms[name].items():
            for (s, e) in introns:
                m = (c,) + normalize_junction(s, e, genome.get(c, "")) in annot_norm
                tot += 1; hit += 1 if m else 0
                if _near(hp_don.get(c, []), s, W_MEMBER) or _near(hp_acc.get(c, []), e, W_MEMBER):
                    hp_tot += 1; hp_hit += 1 if m else 0
        arm_stats[name] = {
            "annotated_match_rate_overall": round(hit / tot, 4) if tot else None,
            "annotated_match_rate_at_hp_abutting": round(hp_hit / hp_tot, 4) if hp_tot else None,
            "n_junctions_scored": tot, "n_junctions_at_hp_abutting": hp_tot,
        }

    # DECISIVE metric: of the reads whose placement the guard CHANGED, did it move the
    # junction TOWARD the annotation (fix) or away (harm)?
    fix = harm = neutral = diffs = diffs_hp = 0
    for k in set(arms["B"]) | set(arms["Bguard"]):
        b = arms["B"].get(k, []); g = arms["Bguard"].get(k, [])
        if b == g:
            continue
        diffs += 1; c = k[0]
        if at_hp(c, b) or at_hp(c, g):
            diffs_hp += 1
        bm, gm = matches(c, b), matches(c, g)
        if gm and not bm:
            fix += 1
        elif bm and not gm:
            harm += 1
        else:
            neutral += 1

    res = {
        "n_annotated": len(annot_norm), "n_hp_abutting": len(hp_junc),
        "hp_drift_margin": args.hp_drift_margin, "hp_drift_min_run": args.hp_drift_min_run,
        "arms": arm_stats,
        "reads_differing_between_arms": diffs,
        "reads_differing_at_hp_abutting": diffs_hp,
        "guard_changes": {"fix_to_annotation": fix, "harm_off_annotation": harm, "neutral": neutral},
    }
    with open(os.path.join(args.outdir, args.out), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))
    print("\nEXPECT: guard_changes fix >> harm (~all fixes, 0 harm); overall match rate not lower "
          "(do-no-harm); match_rate_at_hp_abutting Bguard >= B; reads_differing ~= at_hp_abutting "
          "(guard is HP-specific).")


if __name__ == "__main__":
    main()
