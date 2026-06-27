#!/usr/bin/env python3
"""Tier-2 transcriptome benchmark run (BRANCH A: harness scale-up + minimap2 baseline).

Drives the verified pipeline (GFF -> TranscriptModel -> pbsim3 -> MAF->genome
projection -> align -> score) over a REAL annotated transcriptome, at scale. This
is the SPEC's external-validity tier — but read the SCOPE LABELS below before
quoting any number.

WHAT THIS RUN MEASURES (honestly):
  * ANNOTATED junction RECALL + spurious-junction FDR of ONE aligner (minimap2)
    on realistic pbsim3 reads of real annotated transcripts.
  * For yeast this recall is ~ceiling BY DESIGN — yeast is the SATURATION CONTROL
    that validates the harness end-to-end on real coordinates, NOT a discriminator.

WHAT THIS RUN DOES *NOT* MEASURE (do not claim it does):
  * NOT the panel-failure TAIL. ``panel_unplaced`` here = minimap2's own miss rate,
    not "placed by NO aligner". True tail-sizing needs the multi-aligner panel AND
    an injected hard sub-population (elevated error / repeat); a clean run reports
    tail~0, which is meaningless. See BRANCH B in the handoff.
  * NOT novel-junction (NIC/NNC) recall — full-length reads of annotated
    transcripts carry ANNOTATED junctions only. Novel recall needs isoform
    injection (exon-skip -> NIC, novel-site -> NNC). Spurious-junction FDR (an
    N-op not in truth) IS measured.
  * Exact-indel concordance is REPORTED but NOT a Tier-2 metric (pbsim's per-error
    MAF is the wrong truth for placement; it is the Tier-1 controlled metric).

Run on Sherlock. pbsim/minimap2/samtools from the pbsim3 (or rectify) env; the
driver itself needs the rectify env python (numpy/pandas/pysam for ``import
rectify``). NEVER relay the BAM through the M1.

Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from scripts.benchmark.sim.transcript_model import TranscriptModel  # noqa: E402
from scripts.benchmark.sim.gff_panel import build_panel  # noqa: E402
from scripts.benchmark.sim.pbsim3_wrapper import simulate_and_propagate  # noqa: E402
from rectify.core.benchmark.scorer import load_genome, score_bam  # noqa: E402
from rectify.core.benchmark.truth_schema import read_truth_table  # noqa: E402


def replicate(models, copies: int):
    """pbsim templ mode emits ONE read per template RECORD, so replicate each
    transcript into ``copies`` distinct records to get ``copies`` reads/transcript
    (the count mechanism verified in the live round-trip). Replicates share exons +
    genome_seq (same string ref => no memory blow-up)."""
    out = []
    for m in models:
        for c in range(copies):
            out.append(TranscriptModel(name=f"{m.name}__c{c:04d}", chrom=m.chrom,
                                       strand=m.strand, exons=m.exons,
                                       genome_seq=m.genome_seq))
    return out


def minimap2_splice(genome_fa, reads_fq, out_bam, mm2, st, threads):
    p = subprocess.run([mm2, "-ax", "splice", "-uf", "--eqx", "-k", "14",
                        "-t", str(threads), genome_fa, reads_fq], capture_output=True)
    if p.returncode:
        raise RuntimeError("minimap2 failed: " + p.stderr.decode()[:500])
    s = subprocess.run([st, "sort", "-o", out_bam], input=p.stdout, capture_output=True)
    if s.returncode:
        raise RuntimeError("samtools sort failed: " + s.stderr.decode()[:300])
    subprocess.run([st, "index", out_bam], check=True)


def main():
    ap = argparse.ArgumentParser(description="Tier-2 transcriptome benchmark (minimap2 baseline)")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True, help="genome FASTA reads are aligned to")
    ap.add_argument("--errhmm-model", required=True, help="ERRHMM-ONT.model (DRS) / -HQ (cDNA)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--copies", type=int, default=20, help="reads per transcript (templ records)")
    ap.add_argument("--max-transcripts", type=int, default=None, help="subsample for a quick run")
    ap.add_argument("--include-intronless", action="store_true",
                    help="also simulate intronless mRNAs (placement baseline)")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--pbsim-bin", default="pbsim")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print("[tier2] loading genome ...", file=sys.stderr)
    genome = load_genome(args.genome)
    print("[tier2] building panel from GFF ...", file=sys.stderr)
    models, pairs, donors, acceptors = build_panel(
        args.gff, genome, spliced_only=not args.include_intronless,
        max_transcripts=args.max_transcripts, seed=args.seed)
    spliced = [m for m in models if m.introns()]
    print(f"[tier2] panel: {len(models)} transcripts ({len(spliced)} spliced); "
          f"annotated introns={len(pairs)}; x{args.copies} copies -> "
          f"{len(models) * args.copies} template records", file=sys.stderr)

    rep = replicate(models, args.copies)
    print("[tier2] running pbsim3 + MAF->genome projection ...", file=sys.stderr)
    info = simulate_and_propagate(
        rep, args.out, args.errhmm_model, depth=1, pbsim_bin=args.pbsim_bin,
        seed=args.seed, annotated_pairs=pairs, annotated_donors=donors,
        annotated_acceptors=acceptors, stratum="TRANSCRIPTOME")
    n_truth = int(info["n_truth"])
    print(f"[tier2] simulate_and_propagate: {info}", file=sys.stderr)
    if n_truth == 0:
        print("TIER2 FAILED: 0 truth rows (pbsim/MAF parse?)", file=sys.stderr)
        sys.exit(1)

    bam = os.path.join(args.out, "mm2.bam")
    print("[tier2] minimap2 -ax splice ...", file=sys.stderr)
    minimap2_splice(args.genome, info["reads_fastq"], bam, args.minimap2,
                    args.samtools, args.threads)

    truth = read_truth_table(info["truth_tsv"])
    truth_map = {t.read_id: t for t in truth}
    score = score_bam(bam, truth_map, genome, aligner_name="minimap2")
    s = score.summary()
    # honest scope labels travel WITH the numbers in the artifact
    s["_scope"] = {
        "branch": "A (harness scale-up + minimap2 baseline)",
        "junction_metric": "ANNOTATED-recall + spurious-FDR (NOT novel/NIC/NNC recall)",
        "panel_unplaced_meaning": "minimap2 miss rate — NOT the 5-aligner panel tail",
        "indel_concordance": "REPORTED ONLY (Tier-1 metric; pbsim per-error MAF != placement truth)",
        "yeast_recall": "saturation control / harness validation, NOT discrimination",
    }
    out_json = os.path.join(args.out, "tier2_summary.json")
    with open(out_json, "w") as fh:
        json.dump(s, fh, indent=2, default=str)
    print(json.dumps(s, indent=2, default=str), file=sys.stderr)

    failures = []
    if score.reads_missing_contig:
        failures.append(f"reads_missing_contig={score.reads_missing_contig} "
                        "(genome FASTA does not cover all aligned contigs)")
    if score.reads_placed == 0:
        failures.append("no reads placed by minimap2")
    if score.junction.tp == 0:
        failures.append("0 junction TP — projection/strand bug (no real N-op matched truth)")

    print("\n" + "=" * 60, file=sys.stderr)
    if failures:
        print("TIER2 RUN FAILED:", file=sys.stderr)
        for f in failures:
            print("  - " + f, file=sys.stderr)
        sys.exit(1)
    jr, jf = score.junction.recall, score.junction.fdr
    print(f"TIER2 BRANCH-A PASSED — {score.reads_placed}/{score.reads_scored} reads placed, "
          f"0 missing contigs. ANNOTATED junction recall={jr:.3f} spurious-FDR={jf:.3f} "
          f"(TP={score.junction.tp} FP={score.junction.fp} FN={score.junction.fn}). "
          f"Summary+scope labels: {out_json}", file=sys.stderr)
    print(f"  SCOPE: minimap2-only baseline; panel tail + novel-junction recall are BRANCH B "
          f"(needs the multi-aligner panel + isoform/hard-read injection). indel_conc="
          f"{score.indel.position_exact_concordance:.3f} is Tier-1-only (see header).",
          file=sys.stderr)
    sys.exit(0)


if __name__ == "__main__":
    main()
