"""Swap per-aligner BAM records for re-trimmed validation reads.

For each target read prefix, replace the bundle record in each per-aligner BAM
with the corresponding record from a fresh H2 alignment run, preserving any
XV/XG validation tags that were attached to the bundle record.

Run from the rectify repo root:

    python scripts/validation_data/update_validation_aligner_bams.py \
        --new-aligner-dir /tmp/h2_retrim_align/align_out \
        --target-prefixes 0cb5a111,a146838d,...

The script updates ``rectify/data/validation/aligners/validation_reads.<aligner>.bam``
in place, coord-sorts, and rebuilds the .bai index for each aligner.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pysam


ALIGNERS = ["minimap2", "mapPacBio", "gapmm2", "uLTRA", "deSALT"]
VALIDATION_TAGS = ("XV", "XG")


def swap_aligner_bam(bundle_bam: Path, new_bam: Path, target_prefixes: set[str]) -> tuple[int, int]:
    """Swap records for *target_prefixes* in *bundle_bam* with records from *new_bam*.

    Returns (n_replaced, n_skipped_no_new). XV/XG tags from the bundle records
    are carried over to the replacement records. The bundle BAM is rewritten in
    place via a temp file; the .bai is rebuilt at the end.
    """
    # Index the new BAM by read name → record (full QNAME, not just prefix).
    new_records: dict[str, pysam.AlignedSegment] = {}
    with pysam.AlignmentFile(str(new_bam), "rb") as new_fh:
        for read in new_fh:
            qname = read.query_name
            if any(qname.startswith(p) for p in target_prefixes):
                new_records[qname] = read

    # Walk the bundle BAM. For each target-prefix record, look up the new
    # record by full QNAME, copy XV/XG tags from the old record, and write
    # the new one. For non-target records, write unchanged.
    n_replaced = 0
    n_no_match = 0
    tmp_path = bundle_bam.with_suffix(".swap_tmp.bam")
    with pysam.AlignmentFile(str(bundle_bam), "rb") as in_fh, pysam.AlignmentFile(
        str(tmp_path), "wb", header=in_fh.header
    ) as out_fh:
        for read in in_fh:
            qname = read.query_name
            if any(qname.startswith(p) for p in target_prefixes):
                new = new_records.get(qname)
                if new is None:
                    n_no_match += 1
                    print(f"  WARN: {qname} has no matching new record; keeping old", file=sys.stderr)
                    out_fh.write(read)
                    continue
                # Carry XV/XG tags from the bundle record onto the new record.
                for tag in VALIDATION_TAGS:
                    if read.has_tag(tag):
                        new.set_tag(tag, read.get_tag(tag), value_type=read.get_tag(tag, with_value_type=True)[1])
                out_fh.write(new)
                n_replaced += 1
            else:
                out_fh.write(read)

    # Coord-sort the result in place.
    sorted_path = bundle_bam.with_suffix(".sorted_tmp.bam")
    pysam.sort("-o", str(sorted_path), str(tmp_path))
    tmp_path.unlink()
    shutil.move(str(sorted_path), str(bundle_bam))
    pysam.index(str(bundle_bam))

    return n_replaced, n_no_match


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--new-aligner-dir",
        type=Path,
        required=True,
        help="Directory holding fresh per-aligner BAMs (e.g. /tmp/h2_retrim_align/align_out)",
    )
    p.add_argument(
        "--target-prefixes",
        type=str,
        required=True,
        help="Comma-separated read-name prefixes to swap",
    )
    p.add_argument(
        "--bundle-dir",
        type=Path,
        default=Path("rectify/data/validation/aligners"),
        help="Bundle aligners dir (default: rectify/data/validation/aligners)",
    )
    p.add_argument(
        "--new-bam-stem",
        type=str,
        default="retrimmed",
        help="Stem of fresh aligner BAMs (default: retrimmed)",
    )
    args = p.parse_args()

    targets = set(args.target_prefixes.split(","))
    print(f"Target prefixes ({len(targets)}): {sorted(targets)}")

    total_replaced = 0
    for aligner in ALIGNERS:
        bundle_bam = args.bundle_dir / f"validation_reads.{aligner}.bam"
        new_bam = args.new_aligner_dir / f"{args.new_bam_stem}.{aligner}.bam"
        if not bundle_bam.exists():
            print(f"[{aligner}] SKIP: bundle BAM missing at {bundle_bam}")
            continue
        if not new_bam.exists():
            print(f"[{aligner}] SKIP: new BAM missing at {new_bam}")
            continue
        print(f"[{aligner}] swap: {bundle_bam}  <-  {new_bam}")
        n_replaced, n_no_match = swap_aligner_bam(bundle_bam, new_bam, targets)
        print(f"[{aligner}] replaced {n_replaced} records, {n_no_match} kept unchanged (no match)")
        total_replaced += n_replaced

    print(f"\nTotal replaced records across {len(ALIGNERS)} aligners: {total_replaced}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
