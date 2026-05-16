#!/usr/bin/env python3
"""Build a small bundled validation BAM for a RECTIFY pipeline arm.

Reads a TSV manifest of (source_bam, qname_prefix, xv_label, xg_category)
rows, pulls each matching primary alignment from the source BAM, tags it
with XV (label) and XG (category), and writes a sorted+indexed bundled
BAM to the requested output path. The bundled BAM is consumed by the
arm-specific regression suites
(``tests/test_validation_reads_cdna.py`` and
``tests/test_validation_reads_quantseq_rev.py``).

The manifest uses a SHORT QNAME PREFIX so the file stays human-readable;
the script resolves the prefix to a full QNAME at build time via a single
linear scan of the source BAM.

Usage::

    python scripts/validation_data/build_validation_bam.py \\
        --manifest scripts/validation_data/cdna_manifest.tsv \\
        --output rectify/data/validation/validation_reads_cdna.bam

The manifest TSV has 4 columns (tab-separated, with header)::

    source_bam    qname_prefix    xv_label    xg_category
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import pysam


def load_manifest(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"source_bam", "qname_prefix", "xv_label", "xg_category"}
        if not required.issubset(reader.fieldnames or set()):
            raise SystemExit(
                f"Manifest header must contain {required}; got {reader.fieldnames}"
            )
        for r in reader:
            rows.append(r)
    return rows


def find_read(bam_path: Path, qname_prefix: str) -> pysam.AlignedSegment:
    """Return the first primary alignment whose QNAME starts with prefix."""
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if r.query_name.startswith(qname_prefix):
                return r
    raise SystemExit(
        f"No primary alignment with QNAME prefix {qname_prefix!r} found in {bam_path}"
    )


def build_bundle(manifest_path: Path, output_path: Path) -> None:
    rows = load_manifest(manifest_path)

    # Group rows by source BAM to amortize per-file scans.
    by_source: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_source[r["source_bam"]].append(r)

    # Use the first source BAM's header as the bundle header.
    first_bam = Path(next(iter(by_source.keys())))
    if not first_bam.exists():
        raise SystemExit(f"Source BAM not found: {first_bam}")

    with pysam.AlignmentFile(str(first_bam), "rb") as ref_bam:
        header = ref_bam.header.to_dict()

    # Coerce header SQ entries to be consistent across source BAMs (all
    # source BAMs share the same yeast reference).
    out_path_unsorted = output_path.with_suffix(".unsorted.bam")
    out_path_unsorted.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    with pysam.AlignmentFile(str(out_path_unsorted), "wb", header=header) as out_bam:
        for source_str, group in by_source.items():
            source_bam = Path(source_str)
            if not source_bam.exists():
                raise SystemExit(f"Source BAM not found: {source_bam}")
            wanted = {r["qname_prefix"]: r for r in group}
            found: set[str] = set()
            with pysam.AlignmentFile(str(source_bam), "rb") as bam:
                for read in bam:
                    if read.is_secondary or read.is_supplementary or read.is_unmapped:
                        continue
                    qn = read.query_name
                    match_prefix = None
                    for prefix in wanted:
                        if prefix in found:
                            continue
                        if qn.startswith(prefix):
                            match_prefix = prefix
                            break
                    if match_prefix is None:
                        continue
                    row = wanted[match_prefix]
                    # Tag and write.
                    read.set_tag("XV", row["xv_label"], value_type="Z")
                    read.set_tag("XG", row["xg_category"], value_type="Z")
                    out_bam.write(read)
                    found.add(match_prefix)
                    n_written += 1
            missing = set(wanted) - found
            if missing:
                raise SystemExit(
                    f"In {source_bam}: missing QNAME prefixes {sorted(missing)}"
                )

    # Sort + index.
    pysam.sort("-o", str(output_path), str(out_path_unsorted))
    pysam.index(str(output_path))
    out_path_unsorted.unlink()

    print(f"Wrote {n_written} reads to {output_path}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args(argv)
    build_bundle(args.manifest, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
