#!/usr/bin/env python3
"""Periodicity earn-its-keep test (human DRS) — Task #7 / HANDOFF_2026-05-25 §7.3.

VERDICT (2026-05-26): NO SEPARATION → drop the periodicity dimension.
  - Current _junction_worst_flank_periodicity (min over 30bp window): 0/6 confirmed
    GAA-bridging artifacts flagged; identical distribution to annotated introns.
  - Even a redesigned boundary-adjacent variant (k=6/10/12) does not separate:
    real annotated splice boundaries are themselves frequently low-complexity
    (39% score >=0.8 at k=6), so no threshold spares real junctions while catching
    the artifacts.
  - Root cause: the artifact's low-complexity lives in the READ (Nanopore GAA
    basecall expansion), not the genomic landing site — the genomic flank carries
    only a short boundary-local repeat, then diverse sequence. Genomic-flank
    periodicity is the wrong probe. Any future "earn its keep" attempt should look
    read-side (soft-clip / read-content periodicity), not genomic.
  Full writeup: dev/handoffs/HANDOFF_2026-05-25_human_drs_validation.md §7.3 addendum.

This script is retained so the question is not re-litigated cold. Run on Sherlock
(genome + data live there):
    conda activate rectify
    python dev/periodicity_sanity_human.py
"""
from __future__ import annotations

import random
import statistics
import sys

from rectify.core.splice.junction_scoring import (
    _junction_worst_flank_periodicity,
    _seq_periodicity,
)

FASTA = "/scratch/users/kevinroy/sumner_lab/references/GRCh38_chr5.fa"
GAMED = "/scratch/users/kevinroy/sumner_lab/mpb_fixed_test/SMA_GSB2394_chr5/gamed_examples.tsv"
ANNOT = "/scratch/users/kevinroy/sumner_lab/alignments/SMA_GSB2394_chr5/annotation.junc.bed"
CHROM = "chr5"
N_ANNOT_SAMPLE = 300
SEED = 7


def load_fasta(path: str) -> dict:
    genome: dict = {}
    name = None
    chunks: list = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    genome[name] = "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        genome[name] = "".join(chunks)
    return genome


def load_gamed(path: str) -> list:
    out = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci, si, ei = header.index("chrom"), header.index("intron_start"), header.index("intron_end")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.append((f[ci], int(f[si]), int(f[ei])))
    return out


def load_annotated(path: str, n: int, seed: int) -> list:
    juncs = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] != CHROM:
                continue
            try:
                js, je = int(f[1]), int(f[2])
            except ValueError:
                continue
            if je > js:
                juncs.append((f[0], js, je))
    juncs = sorted(set(juncs))
    rng = random.Random(seed)
    if len(juncs) > n:
        juncs = rng.sample(juncs, n)
    return juncs


def boundary_max(genome: dict, chrom: str, js: int, je: int, k: int) -> float:
    """Periodicity of the k bases IMMEDIATELY abutting the intron on each flank,
    max over the two flanks. Probes boundary-local repeat (the GAA signature),
    unlike _junction_worst_flank_periodicity which takes the min over a 30bp
    sliding window (rewards any complex anchor within reach)."""
    seq = genome.get(chrom, "")
    left = seq[max(0, js - k):js].upper()
    right = seq[je:je + k].upper()
    return max(_seq_periodicity(left), _seq_periodicity(right))


def stats(label: str, vals: list) -> str:
    if not vals:
        return f"{label}: (none)"
    v = sorted(vals)
    n = len(v)
    return (f"{label}: n={n} med={statistics.median(v):.3f} mean={statistics.mean(v):.3f} "
            f"p90={v[min(n - 1, int(0.9 * n))]:.3f} max={v[-1]:.3f} "
            f">=0.8:{100 * sum(x >= 0.8 for x in v) / n:.0f}%")


def main() -> None:
    print(f"Loading genome {FASTA} ...", file=sys.stderr)
    genome = load_fasta(FASTA)
    print(f"  chr5 len={len(genome.get(CHROM, '')):,}", file=sys.stderr)

    gamed = load_gamed(GAMED)
    annot = load_annotated(ANNOT, N_ANNOT_SAMPLE, SEED)

    # --- Current canonical function: worst-flank, min over 30bp window ---------
    print("\n=== _junction_worst_flank_periodicity (current, reach=30) ===")
    gm = []
    for chrom, js, je in gamed:
        p = _junction_worst_flank_periodicity(genome, chrom, js, je)
        gm.append(p)
        print(f"  {chrom}:{js}-{je}  intron={je - js:>7}bp  periodicity={p:.3f}")
    an = [_junction_worst_flank_periodicity(genome, c, s, e) for c, s, e in annot]
    print("  " + stats("GAMED", gm))
    print("  " + stats("ANNOT", an))

    # --- Boundary-adjacent variant (probe only; not wired anywhere) -----------
    print("\n=== boundary-adjacent periodicity (max over both immediate flanks) ===")
    for k in (6, 10, 12):
        gmk = [boundary_max(genome, c, s, e, k) for c, s, e in gamed]
        ank = [boundary_max(genome, c, s, e, k) for c, s, e in annot]
        print(f"  k={k}: GAMED per-junction {[round(x, 2) for x in gmk]}")
        print("        " + stats("GAMED", gmk))
        print("        " + stats("ANNOT", ank))


if __name__ == "__main__":
    main()
