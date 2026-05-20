#!/usr/bin/env python
"""Trace rescue_3ss_truncation (Module 2F) for cat3 reads.

Replicates how `process_read` builds the candidate-junction pool and calls
`rescue_3ss_truncation`. Reports the per-junction shift-search outcome and
the function's final pick.
"""

from __future__ import annotations

from pathlib import Path

import pysam

from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation
from rectify.utils.genome import load_genome, standardize_chrom_name

REPO = Path(__file__).resolve().parents[2]
VAL = REPO / "rectify" / "data" / "validation"
GENOME_FA = REPO / "rectify" / "data" / "genomes" / "saccharomyces_cerevisiae" / \
    "S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
ANNOT = REPO / "rectify" / "data" / "genomes" / "saccharomyces_cerevisiae" / \
    "saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"

CAT3 = {
    "cat3_plus_2":  ("79f61403", "chrI",  "+"),
    "cat3_minus_2": ("28ea9379", "chrII", "-"),
}
# Per-aligner BAMs to test (each is what would be fed into bam_processor as input)
ALIGNERS = ["minimap2", "deSALT", "gapmm2", "uLTRA", "mapPacBio"]


def fetch(bam_path: Path, prefix: str) -> pysam.AlignedSegment | None:
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam:
        if r.is_secondary or r.is_supplementary or r.is_unmapped:
            continue
        if r.query_name.startswith(prefix):
            bam.close()
            return r
    bam.close()
    return None


def extract_junctions_from_cigar(read) -> set[tuple[str, int, int]]:
    chrom = standardize_chrom_name(read.reference_name) or read.reference_name
    out: set[tuple[str, int, int]] = set()
    rp = read.reference_start
    for op, length in read.cigartuples or []:
        if op == 4 or op == 5:
            continue
        if op == 3:
            out.add((chrom, rp, rp + length))
            rp += length
        elif op in (0, 7, 8, 2):
            rp += length
    return out


def main() -> int:
    print("Loading genome + annotated junctions …")
    genome = load_genome(GENOME_FA)
    annot4 = load_annotated_junctions(str(ANNOT))
    print(f"  annotated: {len(annot4):,} (4-tuples)")

    for label, (prefix, chrom, strand) in CAT3.items():
        print("\n" + "=" * 92)
        print(f"{label}  ({prefix}, {chrom}{strand})")
        print("=" * 92)
        # 4-tuple annotated junctions on this chrom
        annot_chrom = [j for j in annot4 if j[0] == chrom]
        print(f"  annotated on {chrom}: {len(annot_chrom)}")
        # Show those within 5 bp of expected coords
        if label == "cat3_plus_2":
            target = (142253, 142619)
            print(f"  expected: {chrom}:{target[0]}-{target[1]}  (annotated EFB1 intron)")
        elif label == "cat3_minus_2":
            target = (366502, 366584)
            print(f"  expected: {chrom}:{target[0]}-{target[1]}")
        nearby = [j for j in annot_chrom
                  if abs(j[1] - target[0]) <= 5 or abs(j[2] - target[1]) <= 5]
        print(f"  annotated nearby (±5): {nearby}")

        for aligner in ALIGNERS:
            bam = VAL / "aligners" / f"validation_reads.{aligner}.bam"
            read = fetch(bam, prefix)
            if read is None:
                continue
            ref_start = read.reference_start
            ref_end = read.reference_end
            cigarstr = read.cigarstring
            n_juncs = [j for c, j, k in [(standardize_chrom_name(read.reference_name), 0, 0)]
                       for j in extract_junctions_from_cigar(read)]
            own_juncs = extract_junctions_from_cigar(read)
            print(f"\n  -- {aligner} --")
            print(f"     ref [{ref_start}, {ref_end})  CIGAR head={cigarstr[:60]}{'…' if len(cigarstr)>60 else ''}")
            print(f"     own junctions: {sorted(own_juncs)}")

            # Build _ss_junctions exactly like bam_processor does
            ss_junctions: set = set()
            ss_junctions.update(annot4)               # 4-tuples!
            ss_junctions.update(own_juncs)            # 3-tuples
            # Sanity: are there 4-tuples in the pool?
            arities = {len(j) for j in ss_junctions}
            print(f"     pool arities: {arities}  (size {len(ss_junctions)})")

            # Call rescue_3ss_truncation
            try:
                result = rescue_3ss_truncation(
                    read, genome, ss_junctions, strand=strand,
                )
            except Exception as e:
                print(f"     ERROR rescue_3ss raised: {type(e).__name__}: {e}")
                continue
            print(f"     result: rescued={result.get('rescued')}  "
                  f"type={result.get('rescue_type')}  "
                  f"junction={result.get('rescued_junction')}  "
                  f"five_prime_corrected={result.get('five_prime_corrected')}  "
                  f"edit_distance={result.get('edit_distance')}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
