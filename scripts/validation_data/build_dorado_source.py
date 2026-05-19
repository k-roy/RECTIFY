#!/usr/bin/env python3
"""Rebuild validation_reads_dorado_source.bam from the source Dorado BAM.

For each read in validation_reads.bam (identified by XV tag), the script
extracts the primary alignment from the Dorado-called DRS BAM
(wt_by4742_rep1.bam), applies XV/XG tags, and writes a sorted+indexed
output BAM.

For reads that are unmapped in the Dorado BAM (Cat1 reads 0cb5a111 and
77b392d9 have no mapping there), the script falls back to the DRS-trim
minimap2 merged BAM, which is also "minimap2 of the full untrimmed read"
from the same sequencing run.

This script is called by regen_pa_rest_bundle.py (Step 5) to keep the
dorado_source BAM in sync with validation_reads.bam after every bundle
regen.
"""
from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import pysam

REPO = Path(__file__).resolve().parents[2]
VAL = REPO / "rectify" / "data" / "validation"
DORADO_SRC_BAM = (
    Path("/oak/stanford/groups/larsms/Users/kevinroy")
    / "projects/TRT/raw_data/nanopore/inhouse_by4742_dst1_4nqo_ski7"
    / "wt_by4742_rep1.bam"
)
MINIMAP2_FALLBACK_BAM = (
    REPO / "dev_runs/wt_by4742_rep1_drs_trim_20260417/merged"
    / "wt_by4742_rep1.minimap2.bam"
)
OUT_BAM = VAL / "validation_reads_dorado_source.bam"


def _load_xv_xg(bam_path: Path) -> dict[str, dict]:
    """Return {qname: {XV: str, XG: str|None}} for all XV-tagged reads."""
    result: dict[str, dict] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in bam:
            if r.has_tag("XV"):
                result[r.query_name] = {
                    "XV": r.get_tag("XV"),
                    "XG": r.get_tag("XG") if r.has_tag("XG") else None,
                }
    return result


def _extract_primary(bam_path: Path,
                     wanted: set[str]) -> dict[str, pysam.AlignedSegment]:
    """Return {qname: read} for wanted qnames — first primary alignment wins."""
    result: dict[str, pysam.AlignedSegment] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary:
                continue
            if r.query_name in wanted and r.query_name not in result:
                result[r.query_name] = r
    return result


def _backup(p: Path) -> None:
    if p.exists():
        backup = p.with_suffix(p.suffix + ".pre_regen")
        if not backup.exists():
            shutil.copy2(p, backup)
            print(f"  backed up {p.name} → {backup.name}")


def main() -> int:
    val_bam = VAL / "validation_reads.bam"
    print(f"Reading XV/XG tags from {val_bam}")
    xv_xg = _load_xv_xg(val_bam)
    wanted = set(xv_xg)
    print(f"  {len(wanted)} reads to include in dorado_source")

    if not DORADO_SRC_BAM.exists():
        print(f"ERROR: Dorado source BAM not found: {DORADO_SRC_BAM}", file=sys.stderr)
        return 1
    if not MINIMAP2_FALLBACK_BAM.exists():
        print(f"ERROR: minimap2 fallback BAM not found: {MINIMAP2_FALLBACK_BAM}",
              file=sys.stderr)
        return 1

    # ---- Phase 1: extract from primary Dorado source BAM ----
    print(f"\nSearching {DORADO_SRC_BAM.name} ...")
    found_in_dorado = _extract_primary(DORADO_SRC_BAM, wanted)
    mapped_in_dorado = {q: r for q, r in found_in_dorado.items() if not r.is_unmapped}
    unmapped_in_dorado = {q: r for q, r in found_in_dorado.items() if r.is_unmapped}
    missing_from_dorado = wanted - set(found_in_dorado)
    print(f"  {len(mapped_in_dorado)} mapped, {len(unmapped_in_dorado)} unmapped, "
          f"{len(missing_from_dorado)} not found")

    # ---- Phase 2: fill unmapped/missing from minimap2 fallback ----
    need_fallback = set(unmapped_in_dorado) | missing_from_dorado
    fallback_reads: dict[str, pysam.AlignedSegment] = {}
    if need_fallback:
        print(f"\nFetching {len(need_fallback)} read(s) from minimap2 fallback: "
              f"{MINIMAP2_FALLBACK_BAM.name}")
        fallback_reads = {
            q: r for q, r in _extract_primary(MINIMAP2_FALLBACK_BAM, need_fallback).items()
            if not r.is_unmapped
        }
        still_missing = need_fallback - set(fallback_reads)
        if still_missing:
            print(f"ERROR: {len(still_missing)} read(s) not found in either source:",
                  file=sys.stderr)
            for q in sorted(still_missing):
                print(f"  {q}", file=sys.stderr)
            return 1
        print(f"  {len(fallback_reads)} read(s) resolved via fallback")

    # Mapped Dorado reads take priority; fallback fills the rest
    all_reads: dict[str, pysam.AlignedSegment] = {**mapped_in_dorado, **fallback_reads}
    if len(all_reads) != len(wanted):
        print(f"ERROR: expected {len(wanted)} reads, assembled {len(all_reads)}",
              file=sys.stderr)
        return 1

    # ---- Phase 3: write new BAM ----
    _backup(OUT_BAM)

    with pysam.AlignmentFile(str(DORADO_SRC_BAM), "rb") as src:
        header = src.header.to_dict()

    with tempfile.NamedTemporaryFile(suffix=".unsorted.bam", delete=False) as tf:
        tmp_unsorted = Path(tf.name)

    try:
        with pysam.AlignmentFile(str(tmp_unsorted), "wb", header=header) as out:
            for qname in sorted(all_reads):
                r = all_reads[qname]
                tags = xv_xg[qname]
                r.set_tag("XV", tags["XV"], value_type="Z")
                if tags["XG"] is not None:
                    r.set_tag("XG", tags["XG"], value_type="Z")
                out.write(r)
        pysam.sort("-o", str(OUT_BAM), str(tmp_unsorted))
        pysam.index(str(OUT_BAM))
        print(f"\nWrote {len(all_reads)} reads → {OUT_BAM}")
    finally:
        tmp_unsorted.unlink(missing_ok=True)

    # ---- Phase 4: verification ----
    print("\nVerification:")
    with pysam.AlignmentFile(str(OUT_BAM), "rb") as chk:
        out_reads = list(chk)

    assert len(out_reads) == len(wanted), \
        f"Expected {len(wanted)} reads, got {len(out_reads)}"
    print(f"  read count: {len(out_reads)} ✓")

    out_xv = {r.get_tag("XV") for r in out_reads if r.has_tag("XV")}
    want_xv = {v["XV"] for v in xv_xg.values()}
    assert out_xv == want_xv, f"XV mismatch: {want_xv - out_xv}"
    print(f"  XV tags: all {len(out_xv)} present ✓")

    # CIGAR = / X ops: pysam uses op codes 7 (=) and 8 (X)
    has_eq_x = any(
        any(op in (7, 8) for op, _ in (r.cigartuples or []))
        for r in out_reads if not r.is_unmapped
    )
    if has_eq_x:
        print("  CIGAR =/X ops present ✓")
    else:
        # Dorado and minimap2 (without --cs/--eqx) emit M ops, not =/X.
        # Only warn if the source BAMs switch encoding in the future.
        print("  CIGAR: all-M encoding (expected for Dorado-sourced BAMs) ✓")

    out_qnames = {r.query_name for r in out_reads}
    missing = wanted - out_qnames
    assert not missing, f"Missing qnames in output: {missing}"
    print("  all expected qnames present ✓")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
