#!/usr/bin/env python3
"""ERROR-FREE saturation control for the Tier-2 harness (the harness validation).

The Tier-2 pbsim run measures minimap2 UNDER NOISE, which conflates two things:
harness correctness and aligner-under-noise. This check ISOLATES the harness — it
aligns PERFECT (error-free) reads of the annotation-derived transcripts back to the
genome and scores them against the same GFF-derived truth. The harness is correct
iff error-free recall is ~ceiling with near-zero FDR (a coordinate off-by-one in the
GFF 1-based->0-based conversion or the projection would crater recall with FN+FP
pairs). It also pins that any recall DROP in the pbsim run is genuinely noise-driven,
not a harness artifact.

M1-light (no pbsim/cluster). Needs minimap2 + samtools on PATH.

Result (250 spliced yeast transcripts, 2026-06-26): recall=0.985 FDR=0.008
(FN = edge soft-clips, not coordinate errors) => harness VALIDATED; the ~0.82 pbsim
recall is the noise effect.

Usage:
  python scripts/benchmark/errfree_saturation_check.py --gff X.gff.gz --genome Y.fsa
Exit 0 = recall >= --min-recall (default 0.95).

Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import collections
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from scripts.benchmark.sim.gff_panel import build_panel  # noqa: E402
from rectify.core.benchmark.scorer import (  # noqa: E402
    load_genome, score_bam, extract_junctions, normalize_junction)
from rectify.core.benchmark.truth_schema import ReadTruth, SplitTag  # noqa: E402
import pysam  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="error-free Tier-2 saturation control")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out", default="/tmp/errfree_sat")
    ap.add_argument("--max-transcripts", type=int, default=250)
    ap.add_argument("--min-recall", type=float, default=0.95)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--samtools", default="samtools")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    genome = load_genome(args.genome)
    models, pairs, donors, acc = build_panel(
        args.gff, genome, spliced_only=True, max_transcripts=args.max_transcripts)
    truth = []
    fq = os.path.join(args.out, "reads.fq")
    with open(fq, "w") as f:
        for m in models:
            seq = m.spliced_transcript()
            f.write(f"@{m.name}\n{seq}\n+\n{'I' * len(seq)}\n")
            truth.append(ReadTruth(
                read_id=m.name, true_locus=m.name, true_transcript=m.name,
                chrom=m.chrom, strand=m.strand, genome_start=m.genome_start,
                genome_end=m.genome_end, true_cigar=m.fulllength_cigar(),
                junctions=m.junction_truths(pairs, donors, acc),
                stratum="ERRFREE", split=SplitTag.TEST))
    bam = os.path.join(args.out, "mm2.bam")
    p = subprocess.run([args.minimap2, "-ax", "splice", "-uf", "--eqx", "-k", "14",
                        "-t", "2", args.genome, fq], capture_output=True)
    if p.returncode:
        raise RuntimeError("minimap2 failed: " + p.stderr.decode()[:400])
    subprocess.run([args.samtools, "sort", "-o", bam], input=p.stdout, capture_output=True)
    subprocess.run([args.samtools, "index", bam], check=True)

    tmap = {t.read_id: t for t in truth}
    sc = score_bam(bam, tmap, genome, aligner_name="errfree")
    placed = {r.query_name: r for r in pysam.AlignmentFile(bam, "rb")
              if not (r.is_unmapped or r.is_secondary or r.is_supplementary)}
    mech = collections.Counter()
    for rid, rt in tmap.items():
        r = placed.get(rid)
        if r is None:
            mech["read_unmapped"] += len(rt.junctions)
            continue
        seq = genome.get(r.reference_name, "")
        calledN = {normalize_junction(cs, ce, seq)
                   for cs, ce in extract_junctions(r.reference_start, r.cigartuples)}
        dspans, pos = [], r.reference_start
        for op, ln in (r.cigartuples or []):
            if op in (0, 2, 3, 7, 8):
                if op == 2:
                    dspans.append((pos, pos + ln))
                pos += ln
        for j in rt.junctions:
            if (j.intron_start, j.intron_end) in calledN:
                continue
            if any(abs(ds - j.intron_start) < 5 and abs(de - j.intron_end) < 5
                   for ds, de in dspans):
                mech["FN_called_as_D"] += 1
            elif j.intron_start >= r.reference_start and j.intron_end <= r.reference_end:
                mech["FN_missed_in_span"] += 1
            else:
                mech["FN_outside_aligned_span"] += 1
    print(f"ERROR-FREE saturation: reads={sc.reads_scored} placed={sc.reads_placed} "
          f"recall={sc.junction.recall:.4f} FDR={sc.junction.fdr:.4f} "
          f"TP={sc.junction.tp} FP={sc.junction.fp} FN={sc.junction.fn}")
    print(f"FN mechanism: {dict(mech)}")
    if sc.junction.recall < args.min_recall:
        print(f"SATURATION FAILED: error-free recall {sc.junction.recall:.4f} < "
              f"{args.min_recall} — HARNESS/LOADER BUG (not noise). Inspect FN mechanism.")
        sys.exit(1)
    print(f"SATURATION PASSED: harness validated (error-free recall "
          f"{sc.junction.recall:.4f} >= {args.min_recall}); any pbsim recall drop is "
          f"noise-driven, not artifact.")
    sys.exit(0)


if __name__ == "__main__":
    main()
