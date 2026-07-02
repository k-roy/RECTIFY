#!/usr/bin/env python3
"""Novel-junction BLIND-SPOT surface — the discriminating native-aligner measurement.

The Fable Director (dev/DIRECTOR_ASSESSMENT_FABLE.md) identified the ONE measurement
that discriminates the entire "is the native aligner justified?" question and that
the program never made: **how often does the real panel fail to produce the TRUE
novel junction at all** — the isoform-flattening rate — as a function of how far the
true splice motif deviates from the GT-AG consensus. Every prior C-gate tested
ARBITRATION ("given a truth member, does the arbiter pick it?"), a different question.

This gate is deliberately **injector-FREE** (error-free reads): the overview's own
claim is that minimap2 snaps a non-canonical junction onto a nearby canonical motif
*even on error-free reads*. So this measures PURE motif-snapping with ZERO dependence
on the placeholder-magnitude error injector (the program's standing bear case) — and
it transfers to real data immediately.

Construction (generalizes gen_junction_discovery_stratum into a DEVIATION LADDER):
each read = one single-intron transcript, spliced error-free; the intron's donor..
acceptor dinucleotides are set to a rung on the ladder; the flanking exons are random
(so there is no annotation to lean on — pure de-novo discovery). Truth = the true
(novel) junction, verified genuinely non-canonical within its ambiguity window (a snap
is then a real shift, not an ambiguity-equivalent TP).

Metric per rung (scored vs TRUTH, ambiguity-aware, via the shipped scorer):
  native_recovery = junction recall of the TRUE site  (TP / (TP+FN))
  blind_spot      = 1 - native_recovery               (the flattening rate)
  snap            = fp_canonical_snap rate            (canonical FP near the true site)
Guards: the CANONICAL rung is the at-ceiling control (recovery ~1.0 ⇒ the metric is
not trivially broken); addressability is clean by construction (error-free exact-match
exon flanks ⇒ the true site is strictly recoverable by a motif-blind aligner, so any
non-recovery is snapping bias, NOT evidence loss); FDR = the canonical rung must not
manufacture false novel junctions.

Reading: a LARGE blind-spot that grows with deviation ⇒ real isoform-flattening the
panel cannot vote its way out of ⇒ a native evidence-weighing member is justified
(build signal). A ~0 blind-spot ⇒ the panel already recovers novel sites ⇒ pivot the
program to correct-step + real-data transfer. Either way it is decisive.

Usage: python scripts/benchmark/novel_junction_blindspot.py --reps 60
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.scorer import score_bam, load_genome  # noqa: E402
from rectify.core.benchmark.truth_schema import ReadTruth, SplitTag  # noqa: E402
from scripts.benchmark.sim.transcript_model import TranscriptModel, Exon  # noqa: E402
import random  # noqa: E402
import pysam  # noqa: E402

_BASES = "ACGT"

# The deviation ladder: (label, (donor, acceptor), hamming-from-GTAG, is_real_biological_motif)
LADDER = [
    ("GT-AG_canon", ("GT", "AG"), 0, True),   # canonical — at-ceiling control
    ("GC-AG",       ("GC", "AG"), 1, True),    # dev1, real minor motif
    ("AT-AC",       ("AT", "AC"), 2, True),    # dev2, real U12 minor spliceosome
    ("GA-AG_1off",  ("GA", "AG"), 1, False),   # dev1, NON-real (isolates pure motif-snapping)
    ("CT-AC_2off",  ("CT", "AC"), 2, False),   # dev2, non-real
    ("CA-TC_deep",  ("CA", "TC"), 4, False),   # deep non-canonical
]


def _rand_unique(n: int, rng: random.Random) -> str:
    return "".join(rng.choice(_BASES) for _ in range(n))


def build_ladder(reps: int, rng: random.Random, exon: int = 200, core: int = 196):
    """Return (refs, reads[(id,seq)], truth[ReadTruth], rung_of[id])."""
    refs, reads, truth, rung_of = {}, [], [], {}
    li = 0
    for (label, (donor, acceptor), dev, real) in LADDER:
        canonical = (label == "GT-AG_canon")
        for i in range(reps):
            chrom = f"nj_{label}_r{i:03d}"
            model = jt = None
            for _attempt in range(30):
                e1 = _rand_unique(exon, rng)
                e2 = _rand_unique(exon, rng)
                intron = donor + _rand_unique(core, rng) + acceptor
                contig = e1 + intron + e2
                i_start, i_end = exon, exon + len(intron)
                model = TranscriptModel(
                    name=chrom, chrom=chrom, strand="+",
                    exons=[Exon(0, exon), Exon(i_end, len(contig))],
                    genome_seq=contig)
                jt = model.junction_truths(annotated_pairs=None)
                # accept only if the constructed junction's canonicity matches intent
                # (ambiguity-aware — so a "non-canonical" rung is genuinely non-canonical
                # within its slide window and a snap is a real shift, not a tied TP)
                if jt and jt[0].canonical == canonical:
                    break
            else:
                continue
            refs[chrom] = contig
            reads.append((chrom, model.spliced_transcript()))
            truth.append(ReadTruth(
                read_id=chrom, true_locus=f"nj_{label}", true_transcript=chrom,
                chrom=chrom, strand="+", genome_start=model.genome_start,
                genome_end=model.genome_end, true_cigar=model.fulllength_cigar(),
                junctions=jt, true_cpa=model.genome_end - 1,
                stratum="NOVEL_JUNCTION_LADDER", split=SplitTag.TEST, coverage=1))
            rung_of[chrom] = label
            li += 1
    return refs, reads, truth, rung_of


def write_fa(refs, path):
    with open(path, "w") as fh:
        for c, s in refs.items():
            fh.write(f">{c}\n{s}\n")
    pysam.faidx(path)


def write_fq(reads, path):
    with open(path, "w") as fh:
        for rid, seq in reads:
            fh.write(f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n")


def run_minimap2(ref_fa, reads_fq, out_bam):
    p = subprocess.run(["minimap2", "-ax", "splice", "-uf", "--eqx", "-k", "14",
                        "-t", "2", ref_fa, reads_fq], capture_output=True)
    if p.returncode:
        raise RuntimeError("minimap2 failed: " + p.stderr.decode()[:400])
    s = subprocess.run(["samtools", "sort", "-o", out_bam], input=p.stdout,
                       capture_output=True)
    if s.returncode:
        raise RuntimeError("samtools sort failed: " + s.stderr.decode()[:300])
    subprocess.run(["samtools", "index", out_bam], check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/nj_blindspot")
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = random.Random(args.seed)

    print(f"[nj] building deviation ladder reps={args.reps} ...", file=sys.stderr)
    refs, reads, truth, rung_of = build_ladder(args.reps, rng)
    ref_fa = os.path.join(args.out, "ref.fa")
    reads_fq = os.path.join(args.out, "reads.fq")
    bam = os.path.join(args.out, "mm2.bam")
    write_fa(refs, ref_fa)
    write_fq(reads, reads_fq)
    genome = load_genome(ref_fa)
    run_minimap2(ref_fa, reads_fq, bam)

    # score per rung (subset truth by rung, score_bam ambiguity-aware)
    truth_by_rung = {}
    for t in truth:
        truth_by_rung.setdefault(rung_of[t.read_id], []).append(t)

    print("\n============ NOVEL-JUNCTION BLIND-SPOT SURFACE (minimap2, error-FREE) ============")
    print("injector-FREE: pure motif-snapping; transfers to real data. metric = native recovery of the TRUE site.")
    print(f"{'rung':14s} {'dev':>3s} {'real':>4s} {'n':>4s} {'recovery':>9s} {'BLINDSPOT':>10s} "
          f"{'snap':>6s} {'fp/rd':>6s}")
    order = [lab for (lab, _m, _d, _r) in LADDER]
    for lab in order:
        ts = truth_by_rung.get(lab, [])
        if not ts:
            continue
        dev = next(d for (l, _m, d, _r) in LADDER if l == lab)
        real = next(r for (l, _m, _d, r) in LADDER if l == lab)
        sc = score_bam(bam, {t.read_id: t for t in ts}, genome, aligner_name=lab)
        n = len(ts)
        recall = sc.junction.recall
        snaps = sc.junction.fp_canonical_snap
        fp = sc.junction.fp
        print(f"{lab:14s} {dev:3d} {('Y' if real else 'n'):>4s} {n:4d} {recall:9.3f} "
              f"{1 - recall:10.3f} {snaps:6d} {fp / max(1, n):6.2f}")

    print("\n---- READING ----")
    print("recovery  = minimap2 natively calls the TRUE novel junction (ambiguity-aware TP)")
    print("BLINDSPOT = 1 - recovery = the isoform-FLATTENING rate the panel can't vote its way out of")
    print("snap      = canonical FP within the snap window of the true non-canonical site")
    print("dev       = Hamming distance of donor+acceptor from the GT-AG consensus; real = known biological motif")
    print("\n---- VERDICT (pre-committed) ----")
    print("  canonical rung recovery ~1.0 = metric sound (at-ceiling control).")
    print("  If BLINDSPOT grows with deviation (esp. > 0 on real GC-AG/AT-AC) => genuine isoform-flattening")
    print("     the flat-affine panel snaps away; error-free + clean flanks => strictly recoverable =>")
    print("     a native evidence-weighing member is JUSTIFIED (build signal). NEXT: confirm on the full")
    print("     5-aligner panel (cluster) + add exon-size / multi-intron / 5'-terminal rungs + error overlay.")
    print("  If BLINDSPOT ~0 across the ladder => the panel already recovers novel sites => PIVOT the program")
    print("     to correct-step + real-data transfer (the native aligner is not justified on discovery grounds).")


if __name__ == "__main__":
    main()
