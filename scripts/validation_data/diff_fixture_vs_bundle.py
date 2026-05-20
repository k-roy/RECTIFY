#!/usr/bin/env python
"""Compare the in-process fixture pipeline output against the bundled TSV.

For Phase 1 Issue 1 (Cat3): the test fixture's `_run_correct_first_pipeline`
runs per-aligner `rectify correct` × 5 then `merge_corrected_tsvs`. The H2
runall that produced the bundle's `corrected_reads.tsv` is supposed to be
the same path, but the validation suite shows divergence between the two.

This script:
  1. Replicates `_run_correct_first_pipeline` for both --no-2h and --with-2h.
  2. Loads the bundle TSV.
  3. Diffs the four Cat3 read rows side-by-side for a chosen column list.

Output: a text table + writes the in-process TSVs to ./_h2_rebaseline/
for later inspection.

Usage:
    .venv/bin/python scripts/validation_data/diff_fixture_vs_bundle.py
"""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
VAL_DIR = REPO / "rectify" / "data" / "validation"
ALIGNER_DIR = VAL_DIR / "aligners"
GENOME = REPO / "rectify" / "data" / "genomes" / "saccharomyces_cerevisiae" / \
    "S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
ANNOT = REPO / "rectify" / "data" / "genomes" / "saccharomyces_cerevisiae" / \
    "saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"
BUNDLE_TSV = VAL_DIR / "corrected_reads.tsv"

ALIGNERS = ["minimap2", "gapmm2", "mapPacBio", "deSALT", "uLTRA"]
WORK = REPO / "_h2_rebaseline" / "fixture_vs_bundle"

CAT3_FULL_READ_IDS = {
    "cat3_plus_1":  "0a28167d-bb74-4d99-8c12-e4911a01b432",
    "cat3_plus_2":  "79f61403-cf63-4522-b555-569590dc4304",
    "cat3_minus_1": "ac4db6da-...",
    "cat3_minus_2": "28ea9379-1518-42f6-94ac-b00709b37b39",
}
# Will be resolved from the bundle TSV on the fly (cat3_minus_1 prefix-only).

COLS = [
    "winning_aligner",
    "five_prime_position",
    "five_prime_rescued",
    "five_prime_exon_cigar",
    "alignment_start",
    "alignment_end",
    "correction_applied",
    "junctions",
    "n_junctions",
    "corrected_3prime",
    "confidence",
]


def run_correct_first_pipeline(out_dir: Path, *, with_module_2h: bool) -> Path:
    """Mirror tests/test_validation_reads.py::_run_correct_first_pipeline."""
    out_dir.mkdir(parents=True, exist_ok=True)
    per_aligner_tsvs: dict = {}
    per_aligner_corrected_bams: dict = {}

    for aligner in ALIGNERS:
        bam = ALIGNER_DIR / f"validation_reads.{aligner}.bam"
        out_tsv = out_dir / f"{aligner}.tsv"
        out_bam = out_dir / f"{aligner}.corrected.bam"
        cmd = [
            sys.executable, "-m", "rectify.cli", "correct",
            str(bam),
            "--genome", str(GENOME),
            "-o", str(out_tsv),
            "--output-bam", str(out_bam),
            "--annotation", str(ANNOT),
        ]
        if with_module_2h:
            for ab_name in ALIGNERS:
                ab_path = ALIGNER_DIR / f"validation_reads.{ab_name}.bam"
                cmd += ["--aligner-bams", f"{ab_name}:{ab_path}"]
        print(f"  [{aligner}] {' '.join(cmd[2:5])} ...", flush=True)
        res = subprocess.run(cmd, capture_output=True)
        if res.returncode != 0:
            print(f"    FAILED rc={res.returncode}", file=sys.stderr)
            print(res.stderr.decode(errors="replace")[:2000], file=sys.stderr)
            raise SystemExit(2)
        per_aligner_tsvs[aligner] = out_tsv
        if out_bam.exists():
            per_aligner_corrected_bams[aligner] = str(out_bam)

    from rectify.core.consensus.corrected_consensus import merge_corrected_tsvs
    merged = out_dir / "merged_corrected.tsv"
    # Match production (single_sample.py:495 / rebuild_2026_05): omit
    # per_aligner_corrected_bams to use the legacy 5-level sort.
    # _ = per_aligner_corrected_bams  # intentionally unused
    merge_corrected_tsvs(
        per_aligner_tsvs,
        merged,
    )
    return merged


def load_tsv(path: Path) -> dict:
    out: dict = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out.setdefault(row["read_id"], []).append(row)
    return out


def first_row(rows: dict, full_or_prefix: str) -> dict | None:
    if full_or_prefix in rows:
        return rows[full_or_prefix][0]
    for rid, rlist in rows.items():
        if rid.startswith(full_or_prefix):
            return rlist[0]
    return None


def fmt(v: str | None, width: int = 24) -> str:
    if v is None:
        return "-" * min(8, width)
    s = str(v)
    return s if len(s) <= width else s[: width - 1] + "…"


def print_diff_for_read(label: str, prefix: str,
                        bundle: dict, no2h: dict, with2h: dict) -> None:
    print("\n" + "=" * 110)
    print(f"{label}  (read_id ~ {prefix})")
    print("=" * 110)
    b = first_row(bundle, prefix)
    n = first_row(no2h, prefix)
    w = first_row(with2h, prefix)
    print(f"  {'column':<26s} {'bundle':<24s} {'fixture-no2h':<24s} {'fixture-with2h':<24s}")
    print(f"  {'-' * 26:<26s} {'-' * 24:<24s} {'-' * 24:<24s} {'-' * 24:<24s}")
    for col in COLS:
        bv = b.get(col) if b else None
        nv = n.get(col) if n else None
        wv = w.get(col) if w else None
        flag = ""
        if bv != nv:
            flag += " !no2h"
        if bv != wv:
            flag += " !with2h"
        print(f"  {col:<26s} {fmt(bv):<24s} {fmt(nv):<24s} {fmt(wv):<24s}{flag}")


def main() -> int:
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True, exist_ok=True)

    print(f"Running fixture pipeline (NO Module 2H) → {WORK/'no2h'}")
    no2h_merged = run_correct_first_pipeline(WORK / "no2h", with_module_2h=False)
    print(f"  → {no2h_merged}")

    print(f"\nRunning fixture pipeline (WITH Module 2H) → {WORK/'with2h'}")
    with2h_merged = run_correct_first_pipeline(WORK / "with2h", with_module_2h=True)
    print(f"  → {with2h_merged}")

    print(f"\nLoading bundle TSV: {BUNDLE_TSV}")
    bundle = load_tsv(BUNDLE_TSV)
    no2h = load_tsv(no2h_merged)
    with2h = load_tsv(with2h_merged)
    print(f"  bundle:        {len(bundle):5d} read_ids")
    print(f"  fixture-no2h:  {len(no2h):5d} read_ids")
    print(f"  fixture-with2h:{len(with2h):5d} read_ids")

    cat3 = [
        ("cat3_plus_1",  "0a28167d"),
        ("cat3_plus_2",  "79f61403"),
        ("cat3_minus_1", "ac4db6da"),
        ("cat3_minus_2", "28ea9379"),
    ]
    for label, prefix in cat3:
        print_diff_for_read(label, prefix, bundle, no2h, with2h)

    return 0


if __name__ == "__main__":
    sys.exit(main())
