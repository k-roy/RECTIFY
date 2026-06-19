#!/usr/bin/env python3
"""Deliverable B — read-level corroboration of GMAP-only recurrent GT-AG NOVEL junctions.

Turns "recurrence + GT-AG candidates" into a corroboration verdict: for each junction
GMAP places, how many INDEPENDENT (non-GMAP) panel aligners independently anchor the SAME
junction (ambiguity-normalized), at what read depth. Buckets each GMAP novel-canonical
junction into {independently_corroborated / gmap_only_recurrent / gmap_only_singleton}.

Reuses the live rectify library (no re-implementation):
  - collect_junction_counts_from_bam  (junction_scoring.py:473) -> anchor-passing Counter
  - normalize_junction / _canonical_within_window / junction_ambiguity_window
                                       (chimeric_consensus.py:59-155) -> ambiguity-aware
  - load_annotated_junctions          (consensus.py:1222) -> exclude catalogued (novel only)

Run ON Sherlock (data locality + the rectify env has AVX-512). chr5 only.

Author: Kevin R. Roy
"""
import argparse, json, sys
from collections import defaultdict

import pysam  # noqa: F401  (ensures env is the rectify env)
from rectify.core.splice.junction_scoring import collect_junction_counts_from_bam
from rectify.core.consensus.chimeric_consensus import (
    normalize_junction, _canonical_within_window, junction_ambiguity_window,
)
from rectify.core.consensus.consensus import load_annotated_junctions

PANEL = ["minimap2", "uLTRA", "deSALT", "mapPacBio", "GMAP"]
RECURRENCE_MIN = 5          # GMAP anchored reads to call "recurrent"
MIN_ANCHOR = 10             # clean exon overhang on both flanks (DEFAULT_MIN_JUNCTION_ANCHOR)


def load_chrom_seq(fa_path, chrom):
    fa = pysam.FastaFile(fa_path)
    name = chrom if chrom in fa.references else (
        "chr" + chrom if "chr" + chrom in fa.references else chrom.replace("chr", ""))
    seq = fa.fetch(name)
    fa.close()
    return seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam-dir", required=True)
    ap.add_argument("--bam-template", default="a549_chr5_trimmed.{aligner}.bam")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--chrom", default="chr5")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seq = load_chrom_seq(args.genome, args.chrom)
    print(f"[B] {args.chrom} loaded: {len(seq):,} bp", file=sys.stderr)

    annot_raw = load_annotated_junctions(args.gtf)
    # normalize annotated set to leftmost 3-tuple on this chrom (dissolve ambiguity)
    annot = set()
    for t in annot_raw:
        c = t[0]
        if c != args.chrom:
            continue
        s, e = normalize_junction(t[1], t[2], seq)
        annot.add((c, s, e))
    print(f"[B] annotated junctions on {args.chrom}: {len(annot):,}", file=sys.stderr)

    # {normalized_junction: {aligner: anchored_count}}
    stratified = defaultdict(dict)
    for aligner in PANEL:
        bam = f"{args.bam_dir}/{args.bam_template.format(aligner=aligner)}"
        counts = collect_junction_counts_from_bam(
            bam, chrom_filter=args.chrom, min_anchor_overhang=MIN_ANCHOR)
        n = 0
        for (c, s, e), k in counts.items():
            ns, ne = normalize_junction(s, e, seq)
            stratified[(c, ns, ne)][aligner] = stratified[(c, ns, ne)].get(aligner, 0) + k
            n += 1
        print(f"[B] {aligner}: {n:,} anchored junctions -> {len(stratified):,} cum normalized",
              file=sys.stderr)

    # Classify GMAP's junctions
    buckets = {"independently_corroborated": [], "gmap_only_recurrent": [],
               "gmap_only_singleton": [], "gmap_annotated": [], "gmap_noncanonical": []}
    for j, per in stratified.items():
        g = per.get("GMAP", 0)
        if g == 0:
            continue
        c, s, e = j
        l_amb, r_amb = junction_ambiguity_window(s, e, seq)
        canonical = _canonical_within_window(s, e, seq, l_amb, r_amb)
        annotated = j in annot
        others = {a: v for a, v in per.items() if a != "GMAP" and v > 0}
        rec = {"junction": [c, s, e], "gmap_reads": g, "canonical": canonical,
               "annotated": annotated, "independent_aligners": sorted(others),
               "n_independent": len(others), "support_by_aligner": per}
        if annotated:
            buckets["gmap_annotated"].append(rec)
        elif not canonical:
            buckets["gmap_noncanonical"].append(rec)
        elif others:
            buckets["independently_corroborated"].append(rec)
        elif g >= RECURRENCE_MIN:
            buckets["gmap_only_recurrent"].append(rec)
        else:
            buckets["gmap_only_singleton"].append(rec)

    summary = {k: len(v) for k, v in buckets.items()}
    summary["gmap_only_recurrent_GTAG_novel_THE_127"] = len(buckets["gmap_only_recurrent"])
    print(f"[B] SUMMARY: {json.dumps(summary, indent=2)}", file=sys.stderr)
    with open(args.out, "w") as fh:
        json.dump({"summary": summary, "buckets": buckets,
                   "params": {"recurrence_min": RECURRENCE_MIN, "min_anchor": MIN_ANCHOR}}, fh)
    print(f"[B] wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
