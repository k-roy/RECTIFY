#!/usr/bin/env python3
"""Recount proposed destinations on stock BAMs, excluding the moving QNAME.

Input TSV: read_id, chrom, intron_start, intron_end (0-based, half-open).
Additional input columns are retained. Supply ORIGINAL aligner BAMs only:

  python scripts/recount_junction_support.py --targets moves.tsv \
      --stock-bam minimap2.bam --stock-bam ultra.bam --out stock_support.tsv

Support is the union of primary QNAMEs across arms, never their sum. Counts
describe exact genomic junctions without strand inference or normalization.
They are independent of corrections, but are not independent molecules in a
PCR-amplified library. The optional anchor count measures contiguous aligned
length (M/=), not sequence identity, complexity, MAPQ or biological truth.

Streams each BAM once with one reader thread. Memory stores only QNAMEs that
support a requested destination, rather than all BAM records or sequences.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import pysam


REQUIRED_COLUMNS = ("read_id", "chrom", "intron_start", "intron_end")
OUTPUT_COLUMNS = ("stock_support_qnames", "stock_other_qnames", "stock_other_anchored_qnames")


def _junctions(read, min_anchor):
    """Yield exact N coordinates and a length-only two-flank anchor flag."""
    pos = read.reference_start
    cigar = read.cigartuples or []
    for i, (op, length) in enumerate(cigar):
        if op == 3 and length > 0:
            left = cigar[i - 1][1] if i > 0 and cigar[i - 1][0] in (0, 7) else 0
            right = cigar[i + 1][1] if i + 1 < len(cigar) and cigar[i + 1][0] in (0, 7) else 0
            yield (read.reference_name, pos, pos + length), min(left, right) >= min_anchor
        if op in (0, 2, 3, 7, 8):
            pos += length


def recount_support(stock_bams, targets, min_anchor=20):
    """Return one evidence row per target, preserving its original columns."""
    if min_anchor < 1:
        raise ValueError("min_anchor must be positive")
    wanted = {(row["chrom"], int(row["intron_start"]), int(row["intron_end"])) for row in targets}
    if any(start < 0 or end <= start for _, start, end in wanted):
        raise ValueError("target junctions must have 0 <= intron_start < intron_end")
    supporters = defaultdict(set)
    anchored = defaultdict(set)
    for path in stock_bams:
        with pysam.AlignmentFile(str(path), "rb", threads=1) as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary or not read.query_name:
                    continue
                for junction, has_anchor in _junctions(read, min_anchor):
                    if junction in wanted:
                        supporters[junction].add(read.query_name)
                        if has_anchor:
                            anchored[junction].add(read.query_name)
    result = []
    for row in targets:
        junction = (row["chrom"], int(row["intron_start"]), int(row["intron_end"]))
        names, anchored_names = supporters[junction], anchored[junction]
        result.append(dict(row,
                           stock_support_qnames=len(names),
                           stock_other_qnames=len(names) - int(row["read_id"] in names),
                           stock_other_anchored_qnames=len(anchored_names) - int(row["read_id"] in anchored_names)))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--stock-bam", required=True, action="append", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--min-anchor", type=int, default=20)
    args = parser.parse_args(argv)
    if args.out.resolve() in {args.targets.resolve(), *(p.resolve() for p in args.stock_bam)}:
        parser.error("--out must differ from all inputs")
    with args.targets.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        columns = reader.fieldnames or []
        if any(col not in columns for col in REQUIRED_COLUMNS):
            parser.error("--targets must contain " + ", ".join(REQUIRED_COLUMNS))
        targets = list(reader)
    result = recount_support(args.stock_bam, targets, args.min_anchor)
    with args.out.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[c for c in columns if c not in OUTPUT_COLUMNS] + list(OUTPUT_COLUMNS), delimiter="\t")
        writer.writeheader()
        writer.writerows(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
