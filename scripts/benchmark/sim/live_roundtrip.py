#!/usr/bin/env python3
"""LIVE pbsim3 round-trip — the Tier-2 wrapper's only unverified piece.

The truth-correctness-critical step (MAF -> genome projection) is already unit-
verified on a synthetic MAF. This driver closes the remaining MECHANICAL gap:
that real pbsim3 actually runs, its MAF + FASTQ read names line up, and the
projected truth scores against a real minimap2 alignment of the simulated reads.

Self-contained: it builds a few synthetic spliced genes (each on its own contig,
+ AND - strand, with GT-AG introns and HP runs in the exons so both junction and
indel truth are exercised), simulates with the packaged ERRHMM-ONT model, aligns
the reads back to the constructed GENOME, and scores.

Run on Sherlock (pbsim3 env). Tiny (a few transcripts, depth ~30) — no cluster
job needed; the login node is fine. NEVER relay the BAM through the M1.

Usage:
  python scripts/benchmark/sim/live_roundtrip.py \
      --errhmm-model .../envs/pbsim3/data/ERRHMM-ONT.model \
      --out /tmp/pbsim_live --depth 30 --seed 7

Exit 0 = pbsim ran, MAF parsed, reads aligned, truth scored (junction TP>0,
no missing-contig, indel concordance computed).

Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
from scripts.benchmark.sim.transcript_model import TranscriptModel, Exon  # noqa: E402
from scripts.benchmark.sim.pbsim3_wrapper import simulate_and_propagate  # noqa: E402
from rectify.core.benchmark.scorer import load_genome, score_bam, open_fasta  # noqa: E402
from rectify.core.benchmark.truth_schema import read_truth_table  # noqa: E402

BASES = "ACGT"


def _rand(n: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def build_models(rng: random.Random, copies: int = 60):
    """Three 2-exon genes (+, +, -), each its OWN contig, with a GT-AG intron and
    a couple of HP runs embedded in the exons so the projection emits indel truth
    on real pbsim deletions inside those runs.

    pbsim3 templ mode emits ONE read per template RECORD (``--depth`` does not
    multiply count; ``--pass-num`` does, but distinctly-named records give cleaner
    per-read truth). So each gene is replicated into ``copies`` distinct template
    records sharing one contig — pbsim then emits ``copies`` reads per gene, all
    aligning back to that gene's single contig.
    """
    models = []
    specs = [("geneP1", "chrP1", "+"),
             ("geneP2", "chrP2", "+"),
             ("geneM1", "chrM1", "-")]
    for name, chrom, strand in specs:
        # exon1: random + an HP run; intron: GT...AG; exon2: random + HP run
        e1 = _rand(120, rng) + "AAAAAAA" + _rand(40, rng)        # A7 run
        intron = "GT" + _rand(160, rng) + "AG"
        e2 = _rand(50, rng) + "CCCCCCCC" + _rand(110, rng)       # C8 run
        contig = e1 + intron + e2
        e1s, e1e = 0, len(e1)
        i_e = e1e + len(intron)
        e2s, e2e = i_e, i_e + len(e2)
        for c in range(copies):
            models.append(TranscriptModel(
                name=f"{name}__c{c:04d}", chrom=chrom, strand=strand,
                exons=[Exon(e1s, e1e), Exon(e2s, e2e)], genome_seq=contig))
    return models


def write_genome(models, path):
    """One FASTA record per DISTINCT contig (replicate models share a chrom)."""
    seen = {}
    for m in models:
        seen.setdefault(m.chrom, m.genome_seq)
    with open(path, "w") as fh:
        for chrom, seq in seen.items():
            fh.write(f">{chrom}\n{seq}\n")


def minimap2_splice(genome_fa, reads_fq, out_bam, mm2="minimap2", st="samtools"):
    p = subprocess.run([mm2, "-ax", "splice", "-uf", "--eqx", "-k", "14",
                        "-t", "2", genome_fa, reads_fq], capture_output=True)
    if p.returncode:
        raise RuntimeError("minimap2 failed: " + p.stderr.decode()[:500])
    s = subprocess.run([st, "sort", "-o", out_bam], input=p.stdout,
                       capture_output=True)
    if s.returncode:
        raise RuntimeError("samtools sort failed: " + s.stderr.decode()[:300])
    subprocess.run([st, "index", out_bam], check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--errhmm-model", required=True)
    ap.add_argument("--out", default="/tmp/pbsim_live")
    ap.add_argument("--depth", type=int, default=30)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--pbsim-bin", default="pbsim")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--copies", type=int, default=60,
                    help="template records per gene = reads per gene (pbsim templ "
                         "mode emits 1 read/record)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    rng = random.Random(args.seed)
    models = build_models(rng, copies=args.copies)
    genome_fa = os.path.join(args.out, "genome.fa")
    write_genome(models, genome_fa)

    print("[live] running pbsim3 ...", file=sys.stderr)
    info = simulate_and_propagate(
        models, args.out, args.errhmm_model, depth=args.depth,
        pbsim_bin=args.pbsim_bin, seed=args.seed)
    print("[live] simulate_and_propagate:", info, file=sys.stderr)
    n_truth = int(info["n_truth"])
    if n_truth == 0:
        print("LIVE FAILED: pbsim/propagation produced 0 truth rows "
              "(MAF parse or read-name mismatch?)", file=sys.stderr)
        sys.exit(1)

    bam = os.path.join(args.out, "mm2.bam")
    print("[live] minimap2 -ax splice ...", file=sys.stderr)
    minimap2_splice(genome_fa, info["reads_fastq"], bam,
                    mm2=args.minimap2, st=args.samtools)

    genome = load_genome(genome_fa)
    truth = read_truth_table(info["truth_tsv"])
    truth_map = {t.read_id: t for t in truth}
    score = score_bam(bam, truth_map, genome, aligner_name="minimap2")
    s = score.summary()
    print(json.dumps(s, indent=2), file=sys.stderr)

    # This driver verifies the Tier-2 MECHANICAL integration only: pbsim runs, its
    # (gzipped) MAF parses, read names line up with the FASTQ/BAM, the MAF->genome
    # projection yields truth, and the scorer runs on a real minimap2 alignment.
    # The Tier-2 SCIENCE metric is junction recall/FDR (per the SPEC's two-tier
    # split). Exact-indel concordance is a TIER-1 CONTROLLED-TRUTH metric and is
    # deliberately NOT a pass/fail here: pbsim's MAF encodes every stochastic
    # sequencing error as an indel, so the projected per-read "indel truth" is
    # thousands of scattered single-base edits; minimap2 legitimately redistributes
    # them and the scorer's per-read has_unexplained gate (designed for Tier-1's
    # single known indel) then zeroes the read. We REPORT it (with this caveat) but
    # never gate on it.
    failures = []
    if score.reads_missing_contig:
        failures.append(f"reads_missing_contig={score.reads_missing_contig} "
                        "(genome FASTA does not cover all aligned contigs)")
    if score.reads_placed == 0:
        failures.append("no reads placed by minimap2 on the simulated reads")
    if score.junction.tp == 0:
        failures.append("0 junction TP — projected junction truth did not match "
                        "any real minimap2 N-op (projection or strand bug)")
    indel_total = score.indel.correct + score.indel.incorrect
    if indel_total == 0:
        failures.append("projection emitted no IndelTruth at all — pbsim error "
                        "model or MAF parse is broken (expect thousands on noisy reads)")

    print("\n" + "=" * 60, file=sys.stderr)
    if failures:
        print("LIVE ROUND-TRIP FAILED (mechanical integration):", file=sys.stderr)
        for f in failures:
            print("  - " + f, file=sys.stderr)
        sys.exit(1)
    print(f"LIVE ROUND-TRIP PASSED (Tier-2 mechanical integration) — pbsim ran, "
          f"gzipped MAF parsed, {score.reads_placed}/{score.reads_scored} reads "
          f"placed, no missing contigs. Tier-2 science metric: junction recall="
          f"{score.junction.recall:.3f}, FDR={score.junction.fdr:.3f} "
          f"(TP={score.junction.tp} FP={score.junction.fp} FN={score.junction.fn}).",
          file=sys.stderr)
    print(f"  NOTE: indel exact-concordance={score.indel.position_exact_concordance:.3f} "
          f"over {indel_total} rows is EXPECTED-low and NOT a Tier-2 metric — it is a "
          f"TIER-1 controlled-truth measure (pbsim's per-error MAF is the wrong truth "
          f"for placement; see header).", file=sys.stderr)
    sys.exit(0)


if __name__ == "__main__":
    main()
