#!/usr/bin/env python
"""Make-or-break check for the brain RNA004 retest: are the a-priori SRRM4-dependent neural
micro-exons INCLUDED in DLPFC reads? (They are skipped in A549/SRRM4-low — that absence is
why the A549 Phase-0 had no test population.)

For each micro-exon interval, classify every read that SPANS it (ref_start < me_start and
ref_end > me_end) as INCLUDE (the micro-exon midpoint falls in an aligned M/=/X block) or
SKIP (it falls inside an N intron). inclusion_ratio = include / (include + skip).
"""
from __future__ import annotations
import argparse, sys
import pysam

# 1-based inclusive GRCh38 coords from apriori_stress_loci_NOTES.md (chr-named)
MICROEXONS = [
    ("DIAPH1", "chr5", 141571428, 141571436, 9),
    ("MYO10",  "chr5", 16780728, 16780741, 14),
    ("ABLIM3", "chr5", 149252201, 149252208, 8),
    ("KIF3A",  "chr5", 132706451, 132706459, 9),
    ("MPC1",   "chr6", 166370218, 166370221, 4),
    ("TRIO",   "chr5", 14497847, 14497874, 28),
]


def op_at(read, pos0):
    """Return 'M' if genomic pos0 falls in an aligned match block, 'N' if in an intron, else None."""
    ref = read.reference_start
    for op, ln in read.cigartuples or []:
        if op in (0, 7, 8):           # M/=/X
            if ref <= pos0 < ref + ln:
                return "M"
            ref += ln
        elif op == 3:                  # N
            if ref <= pos0 < ref + ln:
                return "N"
            ref += ln
        elif op == 2:                  # D
            ref += ln
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    args = ap.parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")
    print(f"{'gene':8s} {'micro-exon':22s} {'len':4s} {'span':6s} {'incl':6s} {'skip':6s} {'incl_ratio'}")
    for gene, chrom, s1, e1, melen in MICROEXONS:
        s0, e0 = s1 - 1, e1               # 0-based half-open
        mid = (s0 + e0) // 2
        incl = skip = span = 0
        try:
            it = bam.fetch(chrom, s0, e0)
        except (ValueError, KeyError):
            print(f"{gene:8s} {chrom}:{s1}-{e1}  (chrom not in BAM)"); continue
        for r in it:
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.cigartuples is None:
                continue
            if not (r.reference_start < s0 and r.reference_end > e0):
                continue              # must span the micro-exon to be informative
            span += 1
            o = op_at(r, mid)
            if o == "M":
                incl += 1
            elif o == "N":
                skip += 1
        ratio = incl / (incl + skip) if (incl + skip) else 0.0
        print(f"{gene:8s} {chrom}:{s1}-{e1:>9} {melen:4d} {span:6d} {incl:6d} {skip:6d} {ratio:.3f}")


if __name__ == "__main__":
    main()
