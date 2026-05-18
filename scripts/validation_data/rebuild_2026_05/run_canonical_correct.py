#!/usr/bin/env python3
"""
Run the canonical correct-first pipeline (per-aligner correct + merge) on the
already-aligned + tagged validation BAMs, matching what fixed `run-all --drs`
would produce.

Steps (matching `_run_correction_per_aligner` + `merge_corrected_tsvs`):
  1. For each aligner BAM in aligners_tagged/, run `rectify correct` WITH
     `--aligner-bams` pointing to ALL 5 BAMs (Module 2H enabled).
  2. Merge the 5 per-aligner TSVs into a single corrected_reads.tsv with
     winning_aligner column.

Outputs to ``canonical_corrected/``:
  - per_aligner/<aligner>/corrected_reads.tsv  (×5)
  - per_aligner/<aligner>/*.rectified_corrected_3end.bam  (×5)
  - per_aligner/<aligner>/*.rectified_pA_tail_trimmed.bam  (×5)
  - corrected_reads.tsv  (merged, with winning_aligner column)
"""
import argparse
import subprocess
import sys
from pathlib import Path

BASE = Path('scripts/validation_data/rebuild_2026_05')
TAGGED_DIR = BASE / 'aligners_tagged'
OUT_DIR = BASE / 'canonical_corrected'
ALIGNERS = ['minimap2', 'mapPacBio', 'gapmm2', 'uLTRA', 'deSALT']


def run_per_aligner_correct() -> dict:
    """Run rectify correct on each per-aligner BAM with --aligner-bams pool."""
    all_aligner_bams = [
        f'{a}:{TAGGED_DIR / f"validation_reads.{a}.bam"}'
        for a in ALIGNERS
    ]
    per_aligner_dir = OUT_DIR / 'per_aligner'
    per_aligner_dir.mkdir(parents=True, exist_ok=True)

    tsvs = {}
    corrected_bams = {}
    for aligner in ALIGNERS:
        out_dir = per_aligner_dir / aligner
        out_dir.mkdir(exist_ok=True)
        bam = TAGGED_DIR / f'validation_reads.{aligner}.bam'
        tsv = out_dir / 'corrected_reads.tsv'
        corrected_bam = out_dir / f'validation_reads.{aligner}.rectified_corrected_3end.bam'
        if tsv.exists():
            print(f"  [{aligner}] cached: {tsv}")
            tsvs[aligner] = tsv
            if corrected_bam.exists():
                corrected_bams[aligner] = corrected_bam
            continue

        print(f"  [{aligner}] correcting {bam.name}...")
        cmd = [
            sys.executable, '-m', 'rectify.cli', 'correct',
            str(bam),
            '--Scer',
            '-o', str(tsv),
            '--write-corrected-bam', str(corrected_bam),
            '--write-softclipped-bam', str(out_dir / f'validation_reads.{aligner}.rectified_pA_tail_trimmed.bam'),
        ]
        for a_bam in all_aligner_bams:
            cmd += ['--aligner-bams', a_bam]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ERROR: [{aligner}]\n{result.stderr[-2000:]}", file=sys.stderr)
            continue
        tsvs[aligner] = tsv
        if corrected_bam.exists():
            corrected_bams[aligner] = corrected_bam
    return tsvs, corrected_bams


def merge(tsvs: dict, corrected_bams: dict = None) -> Path:
    from rectify.core.consensus.corrected_consensus import (
        merge_corrected_tsvs, identify_cat5_candidates
    )
    summary_tsv = OUT_DIR / 'per_aligner' / 'comparison_summary.tsv'
    merged = OUT_DIR / 'corrected_reads.tsv'
    print(f"\nMerging {len(tsvs)} per-aligner TSVs → {merged}")
    # Pass per-aligner corrected BAMs when available to activate HP-edit-distance
    # winner selection (correct-first pipeline order per CLAUDE.md). Falls back
    # to legacy 5-key sort transparently when corrected_bams is empty.
    merge_corrected_tsvs(
        per_aligner_tsvs={k: Path(v) for k, v in tsvs.items()},
        output_tsv=merged,
        summary_tsv=summary_tsv,
        per_aligner_corrected_bams=(
            {k: str(v) for k, v in corrected_bams.items()}
            if corrected_bams else None
        ),
        overhang_table=None,
    )
    identify_cat5_candidates(
        {k: Path(v) for k, v in tsvs.items()},
        output_tsv=OUT_DIR / 'per_aligner' / 'cat5_candidates.tsv',
    )
    return merged


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tsvs, corrected_bams = run_per_aligner_correct()
    if not tsvs:
        print("No per-aligner corrections succeeded", file=sys.stderr)
        return 1
    merged = merge(tsvs, corrected_bams)
    print(f"\nDone: {merged}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
