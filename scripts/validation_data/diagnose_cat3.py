#!/usr/bin/env python
"""Per-read Cat3 diagnostic for Phase 1 Issue 1 (Cat3 5' rescue tiebreaker).

For each Cat3 read in {cat3_plus_1, cat3_plus_2, cat3_minus_1, cat3_minus_2},
report per-aligner:
    - explicit 5' softclip length (read-relative)
    - effective_five_prime_clip (max(softclip, terminal mismatch run))
    - junction_score from score_alignment(...)
    - whether _rescue_5prime_softclip() returns True with the per-read pool
    - junctions list + annotated-junction count
    - canonical splice-site count
And the consensus pick + the three tiebreaker values:
    (a) _count_3prime_agreement
    (b) annotated_junction_matches
    (c) canonical_count
"""

from __future__ import annotations

import sys
from pathlib import Path

import pysam

from rectify.core.consensus.extract import extract_alignment_info
from rectify.core.consensus.scoring import score_alignment, _rescue_5prime_softclip
from rectify.core.consensus.select import select_best_alignment
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.utils.genome import load_genome

REPO = Path(__file__).resolve().parents[2]
VAL = REPO / "rectify" / "data" / "validation"
ALIGNER_DIR = VAL / "aligners"
ALIGNERS = ["minimap2", "mapPacBio", "deSALT", "gapmm2", "uLTRA"]

GENOME_FASTA = REPO / "rectify" / "data" / "genomes" / "saccharomyces_cerevisiae" / \
    "S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
ANNOT_GFF = REPO / "rectify" / "data" / "genomes" / "saccharomyces_cerevisiae" / \
    "saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"

CAT3 = {
    "cat3_plus_1":  "0a28167d",
    "cat3_plus_2":  "79f61403",
    "cat3_minus_1": "ac4db6da",
    "cat3_minus_2": "28ea9379",
}


def collect_reads_by_prefix(bam_path: Path, prefixes: dict[str, str]) -> dict[str, pysam.AlignedSegment]:
    """Return {label: read} for the primary alignment of each prefix in this BAM."""
    out: dict[str, pysam.AlignedSegment] = {}
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for read in bam:
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        for label, prefix in prefixes.items():
            if label in out:
                continue
            if read.query_name.startswith(prefix):
                out[label] = read
                break
        if len(out) == len(prefixes):
            break
    bam.close()
    return out


def main() -> int:
    print(f"genome: {GENOME_FASTA}")
    print(f"gff:    {ANNOT_GFF}")
    print(f"bams:   {ALIGNER_DIR}")
    print()

    genome = load_genome(GENOME_FASTA)
    annotated = load_annotated_junctions(str(ANNOT_GFF))
    print(f"loaded genome ({len(genome)} contigs), annotated junctions ({len(annotated):,})")
    print()

    # Per-aligner: {label: pysam read}
    per_aligner: dict[str, dict[str, pysam.AlignedSegment]] = {}
    for aligner in ALIGNERS:
        bam = ALIGNER_DIR / f"validation_reads.{aligner}.bam"
        per_aligner[aligner] = collect_reads_by_prefix(bam, CAT3)
        missing = set(CAT3) - set(per_aligner[aligner])
        print(f"  {aligner:10s} found {len(per_aligner[aligner])}/{len(CAT3)}"
              f"{'  MISSING ' + ','.join(sorted(missing)) if missing else ''}")
    print()

    for label, prefix in CAT3.items():
        print("=" * 90)
        print(f"{label}  ({prefix})")
        print("=" * 90)

        # Build AlignmentInfo dict for this read across all aligners that have it
        alignments = {}
        for aligner in ALIGNERS:
            read = per_aligner[aligner].get(label)
            if read is None:
                continue
            ai = extract_alignment_info(read, aligner, genome)
            alignments[aligner] = ai

        if not alignments:
            print("  no aligner has this read — skipping")
            continue

        # Build the per-read candidate-junction pool (same logic as select_best_alignment)
        chrom_for_read = next(iter(alignments.values())).chrom
        pool: set = set()
        for j in annotated:
            if j[0] == chrom_for_read:
                pool.add(j)
        for ai in alignments.values():
            for js, je in ai.junctions:
                pool.add((ai.chrom, js, je, ai.strand))

        # Score each + ask whether rescue would fire
        print(f"  chrom={chrom_for_read}")
        print(f"  per-aligner table:")
        print(f"    {'aligner':<10s} {'5pSC':>5s} {'eff5p':>5s} {'score':>7s} "
              f"{'rescued?':>9s} {'#junc':>5s} {'#annot':>6s} {'#canon':>6s} "
              f"{'corr3p':>10s} "
              f"junctions")
        scored: dict[str, float] = {}
        rescued: dict[str, bool] = {}
        for aligner, ai in alignments.items():
            # Score (which internally may flip is_5prime_rescued — fresh ai per call)
            sc = score_alignment(ai, genome, pool)
            scored[aligner] = sc
            # Independent rescue check
            # _rescue_5prime_softclip mutates ai; do it on a fresh extract
            read = per_aligner[aligner][label]
            ai_for_rescue = extract_alignment_info(read, aligner, genome)
            r = _rescue_5prime_softclip(ai_for_rescue, genome, pool)
            rescued[aligner] = bool(r)

            n_annot = sum(
                1 for j in ai.junctions
                if (ai.chrom, j[0], j[1], ai.strand) in annotated
            )
            n_canon = ai.canonical_count
            junc_str = ",".join(f"{a}-{b}" for a, b in ai.junctions[:4])
            if len(ai.junctions) > 4:
                junc_str += f",...(+{len(ai.junctions)-4})"
            print(f"    {aligner:<10s} {ai.five_prime_softclip:>5d} "
                  f"{ai.effective_five_prime_clip:>5d} "
                  f"{sc:>7.2f} "
                  f"{('YES' if r else 'no'):>9s} "
                  f"{len(ai.junctions):>5d} "
                  f"{n_annot:>6d} "
                  f"{n_canon:>6d} "
                  f"{ai.corrected_3prime if ai.corrected_3prime is not None else '-':>10} "
                  f"{junc_str}")

        # Run the real consensus selector
        res = select_best_alignment(alignments, genome, annotated)
        print(f"  → consensus: best={res.best_aligner}  conf={res.confidence}  "
              f"tied={res.tied_aligners}  was_5prime_rescued={res.was_5prime_rescued}")

        # Decompose the three tiebreaker values per aligner
        max_score = max(a.junction_score for a in alignments.values())
        tied = [n for n, a in alignments.items() if a.junction_score == max_score]
        if len(tied) > 1:
            all_c3p = [a.corrected_3prime for a in alignments.values() if a.corrected_3prime is not None]
            print(f"  tiebreaker values (only tied aligners @ max_score={max_score:.2f}):")
            print(f"    {'aligner':<10s} {'(a)agree3p':>11s} {'(b)#annot':>10s} {'(c)#canon':>10s}")
            for name in tied:
                a = alignments[name]
                pos = a.corrected_3prime
                agree = sum(1 for p in all_c3p if p == pos) if pos is not None else 0
                n_annot = sum(
                    1 for j in a.junctions
                    if (a.chrom, j[0], j[1], a.strand) in annotated
                )
                print(f"    {name:<10s} {agree:>11d} {n_annot:>10d} {a.canonical_count:>10d}")
        print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
