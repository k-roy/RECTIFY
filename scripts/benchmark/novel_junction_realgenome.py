#!/usr/bin/env python3
"""Real-genome, multi-intron novel-junction blind-spot ladder — the honest test.

The synthetic ladder (novel_junction_blindspot.py) put each read on its OWN tiny
random contig: no real sequence composition, single intron, no decoy GT-AG sites
elsewhere in a genome. That is the EASY slice, and mapPacBio "broke the herd" there.
This builds the ladder on a REAL GENOME (bundled yeast S288C; --human-genome for a
human port) with:
  * reads aligned to the FULL genome  -> real composition + real DECOY GT-AG sites
    (the context where snapping / mis-mapping actually happens);
  * constructed NOVEL introns at real genomic loci, scanned so the donor..acceptor
    dinucleotides hit a target motif rung (graded deviation from GT-AG);
  * MULTI-INTRON transcripts (--n-introns 2) to test compounding;
  * SHORT exons (--exon 30) for the short-anchor stressor;
  * ERROR overlay (--error-rate, RNA004-bulk; truth-invariant).

Truth = the constructed intron coordinates in GENOME space (novel by construction —
not annotated). A canonical GT-AG rung placed at a NOVEL position is the annotation-
snap control; the non-canonical rungs are the motif-snap discovery targets, now with
real nearby decoy canonical sites to snap TO.

Emits a corpus (ref.fa = the genome; reads.fastq; truth.tsv) for the SAME cluster
panel harness (run_panel_on_corpus.py + panel_blindspot_score.py). Reference is the
whole genome, so minimap2/mapPacBio/deSALT/uLTRA face the real decoy problem.

Usage:
  python scripts/benchmark/novel_junction_realgenome.py --reps 60 --n-introns 2 \
      --emit-corpus <dir>            # yeast S288C by default
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.truth_schema import ReadTruth, SplitTag, write_truth_table  # noqa: E402
from scripts.benchmark.sim.transcript_model import TranscriptModel, Exon  # noqa: E402
import pysam  # noqa: E402

YEAST_GZ = os.path.join(os.path.dirname(__file__), "..", "..", "rectify", "data",
                        "genomes", "saccharomyces_cerevisiae",
                        "S288C_reference_sequence_R64-5-1_20240529.fsa.gz")

# same deviation ladder as the synthetic version
LADDER = [
    ("GT-AG_canon", ("GT", "AG"), 0, True),
    ("AT-AC",       ("AT", "AC"), 2, True),
    ("GA-AG_1off",  ("GA", "AG"), 1, False),
    ("CT-AC_2off",  ("CT", "AC"), 2, False),
    ("CA-TC_deep",  ("CA", "TC"), 4, False),
]


def find_intron(seq, start, donor, acceptor, min_intron=60, max_intron=400,
                max_scan=4000):
    """First [ds, ae) at/after `start` with seq[ds:ds+2]==donor and
    seq[ae-2:ae]==acceptor and intron length in [min,max]. Returns (ds,ae) or None."""
    L = len(seq)
    ds = start
    scanned = 0
    while ds < L - 2 and scanned < max_scan:
        if seq[ds:ds + 2] == donor:
            hi = min(ds + max_intron, L)
            for ae in range(ds + min_intron, hi):
                if seq[ae - 2:ae] == acceptor:
                    return ds, ae
        ds += 1
        scanned += 1
    return None


def build(genome_fa, reps, rng, exon=200, n_introns=1, error_rate=0.0):
    ep = None
    if error_rate > 0:
        from scripts.benchmark.sim.error_injector import inject, placeholder_params
        ep = placeholder_params(base_rate=error_rate)
    fa = pysam.FastaFile(genome_fa)
    chroms = [c for c in fa.references if fa.get_reference_length(c) > 60000]
    seqs = {c: fa.fetch(c).upper() for c in chroms}
    reads, truth, rung_of = [], [], {}
    for (label, (donor, acceptor), dev, real) in LADDER:
        want_canon = (label == "GT-AG_canon")
        made = tries = 0
        while made < reps and tries < reps * 80:
            tries += 1
            chrom = rng.choice(chroms)
            cs = seqs[chrom]
            L = len(cs)
            g0 = rng.randint(1000, L - 6000)
            e_start = g0
            exons, introns, cur, ok = [], [], g0 + exon, True
            for _k in range(n_introns):
                intr = find_intron(cs, cur, donor, acceptor)
                if intr is None:
                    ok = False
                    break
                ds, ae = intr
                exons.append(Exon(e_start, ds))
                introns.append((ds, ae))
                e_start = ae
                cur = ae + exon
            if not ok:
                continue
            g_end = e_start + exon
            if g_end >= L:
                continue
            exons.append(Exon(e_start, g_end))
            rid = f"njr_{label}_r{made:03d}"
            model = TranscriptModel(name=rid, chrom=chrom, strand="+",
                                    exons=exons, genome_seq=cs)
            jt = model.junction_truths(annotated_pairs=None)
            # every constructed intron must classify with the intended canonicity
            # (ambiguity-aware) so a snap is a genuine shift, not a tied TP
            if not jt or len(jt) != n_introns or \
               any(j.canonical != want_canon for j in jt):
                continue
            read = model.spliced_transcript()
            if ep is not None:
                read, _ = inject(read, ep, rng)
            reads.append((rid, read))
            truth.append(ReadTruth(
                read_id=rid, true_locus=f"njr_{label}", true_transcript=rid,
                chrom=chrom, strand="+", genome_start=g0, genome_end=g_end,
                true_cigar=model.fulllength_cigar(), junctions=jt,
                true_cpa=g_end - 1, stratum="NJ_REALGENOME",
                split=SplitTag.TEST, coverage=1))
            rung_of[rid] = label
            made += 1
        print(f"[njr] rung {label}: made {made}/{reps} ({tries} tries)", file=sys.stderr)
    return genome_fa, reads, truth, rung_of


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--exon", type=int, default=200,
                    help="flank exon length; ~30 for the short-anchor stressor")
    ap.add_argument("--n-introns", type=int, default=1)
    ap.add_argument("--error-rate", type=float, default=0.0)
    ap.add_argument("--genome", default=os.path.abspath(YEAST_GZ),
                    help="alignment reference (default bundled yeast S288C)")
    ap.add_argument("--emit-corpus", required=True,
                    help="write ref.fa (=genome) + reads.fastq + truth.tsv + rung_map.tsv")
    args = ap.parse_args()
    rng = random.Random(args.seed)

    print(f"[njr] building REAL-GENOME ladder reps={args.reps} n_introns={args.n_introns} "
          f"exon={args.exon} error={args.error_rate} genome={os.path.basename(args.genome)}",
          file=sys.stderr)
    genome_fa, reads, truth, rung_of = build(
        args.genome, args.reps, rng, exon=args.exon,
        n_introns=args.n_introns, error_rate=args.error_rate)

    d = args.emit_corpus
    os.makedirs(d, exist_ok=True)
    # ref.fa = the full genome (decompress if .gz) so the panel aligns to real
    # sequence with real decoy sites.
    ref_out = os.path.join(d, "ref.fa")
    fa = pysam.FastaFile(genome_fa)
    with open(ref_out, "w") as fh:
        for c in fa.references:
            fh.write(f">{c}\n{fa.fetch(c)}\n")
    pysam.faidx(ref_out)
    with open(os.path.join(d, "reads.fastq"), "w") as fh:
        for rid, seq in reads:
            fh.write(f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n")
    nt = write_truth_table(truth, os.path.join(d, "truth.tsv"))
    with open(os.path.join(d, "rung_map.tsv"), "w") as fh:
        fh.write("read_id\trung\n")
        for rid, lab in rung_of.items():
            fh.write(f"{rid}\t{lab}\n")
    print(f"[njr] corpus -> {d}: ref={os.path.basename(ref_out)} (full genome), "
          f"{len(reads)} reads, {nt} truth rows. Run the panel harness on it.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
