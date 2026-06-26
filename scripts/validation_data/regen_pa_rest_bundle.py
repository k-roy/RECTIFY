#!/usr/bin/env python
"""Regenerate the validation bundle's rectified_pA_tail_{trimmed,soft_clipped}.bam.

Uses `rectify correct` (not the lower-level write_*_bam helpers) so Module 2H
fires on each per-aligner BAM. The 2H pass is required for canonical
GT|AG splice-site refinement on aligners with native N-ops (mapPacBio).

Flow:
  1. For each per-aligner raw BAM: `rectify correct --aligner-bams ...
     --write-corrected-bam X --write-softclipped-bam Y`. Module 2H runs
     because annotation + aligner_bams are supplied.
  2. Re-merge: load each per-aligner TSV and run merge_corrected_tsvs with
     HP-aware edit distance scoring (per_aligner_corrected_bams + genome +
     penalty_table + overhang_table) so the bundle's winner-selection
     reflects the precise empirical-penalty ED metric. Also emits a
     per-aligner comparison summary at
     ``rectified/per_aligner_summary.tsv`` consumed by the renderer's
     end-correction + winner-selection panel.
  3. Filter each per-aligner BAM by winning aligner, concat into the bundle
     paths.
  4. Run restore_polya_from_parquet.py to add poly-A tail to the
     soft-clipped variant.

Backups: existing bundle files copied to *.pre_regen before overwriting (if
no pre_regen already exists).
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam

REPO = Path(__file__).resolve().parents[2]
VAL = REPO / "rectify" / "data" / "validation"
ALIGNER_DIR = VAL / "aligners"
RECTIFIED = VAL / "rectified"
GENOME = REPO / "rectify/data/genomes/saccharomyces_cerevisiae" / \
    "S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
ANNOT = REPO / "rectify/data/genomes/saccharomyces_cerevisiae" / \
    "saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"
CORRECTED_TSV_OUT = VAL / "corrected_reads.tsv"
SUMMARY_TSV_OUT = RECTIFIED / "per_aligner_summary.tsv"
PER_ALIGNER_DIR = RECTIFIED / "per_aligner"
SOFTCLIP_OUT_BAM = RECTIFIED / "rectified_pA_tail_soft_clipped.bam"
TRIMMED_OUT_BAM = RECTIFIED / "rectified_pA_tail_trimmed.bam"
PARQUET_TSV = REPO / "scripts/validation_data/rebuild_2026_05/trimmed/" \
    "validation_reads_polya_trim_metadata.tsv"
PENALTY_DIR = REPO / "rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables"
PENALTY_TSV = PENALTY_DIR / "penalty_scores.tsv"
STR_PENALTY_TSV = PENALTY_DIR / "str_penalty_scores.tsv"
OVERHANG_TSV = PENALTY_DIR / "junction_overhang_table.tsv"

ALIGNERS = ["minimap2", "gapmm2", "mapPacBio", "deSALT", "uLTRA"]


def _backup(p: Path) -> None:
    if p.exists():
        backup = p.with_suffix(p.suffix + ".pre_regen")
        if not backup.exists():
            shutil.copy2(p, backup)
            print(f"  backed up {p.name} → {backup.name}")
        bai = p.with_suffix(".bam.bai")
        if bai.exists() and not bai.with_suffix(".bai.pre_regen").exists():
            shutil.copy2(bai, bai.with_suffix(".bai.pre_regen"))


def main() -> int:
    # Accept both flag spellings: the new --skip-polya-restore and the
    # legacy --skip-polya-splice (kept for backwards compatibility with
    # existing wrappers; "splice" was a misnomer — splicing refers to
    # intron removal, not pA-tail re-attachment).
    skip_polya = ("--skip-polya-restore" in sys.argv
                  or "--skip-polya-splice" in sys.argv)
    print(f"Genome:     {GENOME}")
    print(f"Annotation: {ANNOT}")

    with tempfile.TemporaryDirectory() as td_str:
        td = Path(td_str)
        per_aligner_tsvs: dict[str, Path] = {}
        per_aligner_trimmed: dict[str, Path] = {}
        per_aligner_soft: dict[str, Path] = {}

        # ---- Step 1: per-aligner rectify correct with --aligner-bams ----
        for aligner in ALIGNERS:
            raw = ALIGNER_DIR / f"validation_reads.{aligner}.bam"
            if not raw.exists():
                print(f"  ! missing {raw}; skipping")
                continue
            tsv = td / f"{aligner}.tsv"
            hc_bam = td / f"{aligner}.trimmed.bam"
            sc_bam = td / f"{aligner}.softclipped.bam"
            cmd = [
                sys.executable, "-m", "rectify.cli", "correct",
                str(raw),
                "--genome", str(GENOME),
                "--annotation", str(ANNOT),
                "-o", str(tsv),
                "--write-corrected-bam", str(hc_bam),
                "--write-softclipped-bam", str(sc_bam),
                # Commit B changed `rectify correct`'s default to manifest-only;
                # we need the flat per-aligner TSV at -o for merge_corrected_tsvs.
                "--emit-merged-tsv",
            ]
            for ab in ALIGNERS:
                cmd += ["--aligner-bams",
                        f"{ab}:{ALIGNER_DIR / f'validation_reads.{ab}.bam'}"]
            print(f"\n→ rectify correct: {aligner}")
            res = subprocess.run(cmd, capture_output=True)
            if res.returncode != 0:
                print(f"  FAILED rc={res.returncode}", file=sys.stderr)
                print(res.stderr.decode(errors="replace")[:2000], file=sys.stderr)
                return 2
            per_aligner_tsvs[aligner] = tsv
            if hc_bam.exists():
                per_aligner_trimmed[aligner] = hc_bam
            if sc_bam.exists():
                per_aligner_soft[aligner] = sc_bam

        # ---- Step 2: merge per-aligner TSVs (HP-aware ED mode) ----
        from rectify.core.consensus.corrected_consensus import merge_corrected_tsvs
        merged_tsv = td / "merged.tsv"
        summary_tsv = td / "per_aligner_summary.tsv"

        # Load genome (~12 MB yeast) for HP-aware edit distance.
        print(f"\n→ loading genome for HP-ED scoring: {GENOME.name}")
        genome_dict: dict[str, str] = {}
        with pysam.FastaFile(str(GENOME)) as fa:
            for chrom in fa.references:
                genome_dict[chrom] = fa.fetch(chrom)

        # Optional empirical penalty table (HP-aware del/ins costs).
        penalty_table = None
        if PENALTY_TSV.exists():
            try:
                from rectify.core.splice.junction_refiner import HpPenaltyTable
                str_arg = str(STR_PENALTY_TSV) if STR_PENALTY_TSV.exists() else None
                penalty_table = HpPenaltyTable.from_tsv(
                    str(PENALTY_TSV), str_path=str_arg
                )
                print(f"  loaded penalty table: {PENALTY_TSV.name}"
                      + (f" + {STR_PENALTY_TSV.name}" if str_arg else ""))
            except Exception as exc:
                print(f"  WARNING: penalty table load failed: {exc}; "
                      f"falling back to flat costs", file=sys.stderr)

        # Optional junction overhang filter (chimera pre-filter).
        overhang_table = None
        if OVERHANG_TSV.exists():
            try:
                from rectify.core.splice.calibrate_junction_overhang import OverhangTable
                overhang_table = OverhangTable.from_tsv(OVERHANG_TSV)
                print(f"  loaded overhang table: {OVERHANG_TSV.name}")
            except Exception as exc:
                print(f"  WARNING: overhang table load failed: {exc}; "
                      f"no chimera pre-filter", file=sys.stderr)

        print(f"\n→ merge_corrected_tsvs (HP-ED mode)")
        merge_corrected_tsvs(
            per_aligner_tsvs={k: v for k, v in per_aligner_tsvs.items()},
            output_tsv=merged_tsv,
            summary_tsv=summary_tsv,
            per_aligner_corrected_bams={
                k: str(v) for k, v in per_aligner_trimmed.items()
            },
            genome=genome_dict,
            penalty_table=penalty_table,
            overhang_table=overhang_table,
        )
        # Read winners from the merged TSV
        import csv
        winners: dict[str, str] = {}
        with open(merged_tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                winners[row["read_id"]] = row["winning_aligner"]
        by_aligner: dict[str, set[str]] = {}
        for rid, a in winners.items():
            by_aligner.setdefault(a, set()).add(rid)
        print(f"  winners: " +
              ", ".join(f"{a}={len(by_aligner.get(a, set()))}" for a in ALIGNERS))

        # ---- Step 3: filter each per-aligner BAM by winning aligner, concat ----
        for label, src_map, dst in [
            ("trimmed",     per_aligner_trimmed, TRIMMED_OUT_BAM),
            ("softclipped", per_aligner_soft,    SOFTCLIP_OUT_BAM),
        ]:
            _backup(dst)
            filtered: list[Path] = []
            for aligner in ALIGNERS:
                bp = src_map.get(aligner)
                won = by_aligner.get(aligner, set())
                if bp is None or not won:
                    continue
                out = td / f"{aligner}.{label}.filtered.bam"
                with pysam.AlignmentFile(str(bp), "rb") as bam_in, \
                     pysam.AlignmentFile(str(out), "wb",
                                         header=bam_in.header) as bam_out:
                    for r in bam_in:
                        if r.is_secondary or r.is_supplementary:
                            continue
                        if r.query_name in won:
                            bam_out.write(r)
                filtered.append(out)
            unsorted = td / f"merged.{label}.unsorted.bam"
            pysam.merge("-f", str(unsorted), *[str(p) for p in filtered])
            pysam.sort("-o", str(dst), str(unsorted))
            pysam.index(str(dst))
            print(f"  wrote {dst}")

        # ---- Update the bundle's corrected_reads.tsv too ----
        _backup(CORRECTED_TSV_OUT.with_suffix(".tsv"))  # creates .pre_regen.tsv (best-effort)
        shutil.copy2(merged_tsv, CORRECTED_TSV_OUT)
        print(f"  wrote {CORRECTED_TSV_OUT}")
        # Per-aligner comparison summary — consumed by the renderer's
        # end-correction + winner-selection panel.
        if summary_tsv.exists():
            _backup(SUMMARY_TSV_OUT)
            shutil.copy2(summary_tsv, SUMMARY_TSV_OUT)
            print(f"  wrote {SUMMARY_TSV_OUT}")
        else:
            print(f"  WARNING: summary_tsv not produced by merge_corrected_tsvs",
                  file=sys.stderr)

        # Persist per-aligner corrected BAMs to PER_ALIGNER_DIR/ so the
        # renderer can compute raw edit distance per aligner (HP-aware vs.
        # flat-penalty comparison) and inspect any aligner's CIGAR for any
        # read without re-running the regen pipeline.
        # Also persist per-aligner softclipped variants and splice the pA
        # tail into each so the renderer can show per-aligner pA-rest tracks.
        PER_ALIGNER_DIR.mkdir(parents=True, exist_ok=True)
        for aligner, bam_path in per_aligner_soft.items():
            sc_dst = PER_ALIGNER_DIR / f"{aligner}.softclipped.bam"
            _backup(sc_dst)
            shutil.copy2(bam_path, sc_dst)
            try:
                pysam.index(str(sc_dst))
            except Exception as exc:
                print(f"  WARNING: pysam.index failed on {sc_dst.name}: {exc}",
                      file=sys.stderr)
            # Splice pA tail in-place to create the pA-rest variant.
            pa_dst = PER_ALIGNER_DIR / f"{aligner}.pA_rest.bam"
            if not skip_polya:
                splice_cmd = [
                    sys.executable,
                    str(REPO / "scripts/validation_data/restore_polya_from_parquet.py"),
                    "--in-bam", str(sc_dst),
                    "--parquet", str(PARQUET_TSV),
                    "--out-bam", str(pa_dst),
                ]
                res = subprocess.run(splice_cmd, capture_output=True)
                if res.returncode == 0:
                    try:
                        pysam.index(str(pa_dst))
                    except Exception as exc:
                        print(f"  WARNING: pysam.index failed on {pa_dst.name}: {exc}",
                              file=sys.stderr)
                else:
                    print(f"  WARNING: pA-tail restore failed for {aligner}; "
                          f"falling back to softclipped without pA",
                          file=sys.stderr)
                    if pa_dst.exists():
                        pa_dst.unlink()
                    shutil.copy2(sc_dst, pa_dst)
                    pysam.index(str(pa_dst))
            else:
                shutil.copy2(sc_dst, pa_dst)
                pysam.index(str(pa_dst))
        for aligner, bam_path in per_aligner_trimmed.items():
            dst = PER_ALIGNER_DIR / f"{aligner}.trimmed.bam"
            _backup(dst)
            shutil.copy2(bam_path, dst)
            try:
                pysam.index(str(dst))
            except Exception as exc:
                print(f"  WARNING: pysam.index failed on {dst.name}: {exc}",
                      file=sys.stderr)
        print(f"  wrote {len(per_aligner_trimmed)} per-aligner corrected BAMs "
              f"in {PER_ALIGNER_DIR}")

    # ---- Step 3.5: regenerate 3'-end bedgraphs ----
    # Bundle ships `rectified/corrected_3ends.{plus,minus}.bedgraph` for
    # the plotter / per-category review. After HP-ED winner-selection
    # picks different aligners, the per-read corrected_3prime values
    # change → bedgraph counts at old positions go stale. Regenerate
    # from the freshly-written corrected_reads.tsv.
    try:
        from rectify.core.bam.bedgraph_writers import write_corrected_reads_bedgraph
        _bg_prefix = RECTIFIED / "corrected_3ends"
        print(f"\n→ regenerate 3'-end bedgraphs → {_bg_prefix}.{{plus,minus}}.bedgraph")
        _bg_counts = write_corrected_reads_bedgraph(
            str(CORRECTED_TSV_OUT),
            str(_bg_prefix),
        )
        print(f"  plus={_bg_counts.get('plus', 0)} positions  "
              f"minus={_bg_counts.get('minus', 0)} positions")
    except Exception as _bg_exc:
        print(f"  WARNING: bedgraph regen failed: {_bg_exc}", file=sys.stderr)

    # ---- Step 4b: pre-RECTIFY baseline — splice pA tail into the ORIGINAL
    # minimap2 aligner BAM so the renderer can show what the read would
    # look like without RECTIFY (untrimmed tail + samtools-style 3' end).
    if not skip_polya:
        unrect_dst = PER_ALIGNER_DIR / "minimap2_unrectified.bam"
        _backup(unrect_dst)
        splice_cmd = [
            sys.executable,
            str(REPO / "scripts/validation_data/restore_polya_from_parquet.py"),
            "--in-bam", str(ALIGNER_DIR / "validation_reads.minimap2.bam"),
            "--parquet", str(PARQUET_TSV),
            "--out-bam", str(unrect_dst),
        ]
        res = subprocess.run(splice_cmd, capture_output=True)
        if res.returncode == 0:
            try:
                pysam.index(str(unrect_dst))
            except Exception as exc:
                print(f"  WARNING: pysam.index failed on {unrect_dst.name}: {exc}",
                      file=sys.stderr)
            print(f"  wrote {unrect_dst}")
        else:
            print(f"  WARNING: unrectified-minimap2 pA restore failed", file=sys.stderr)

    # ---- Step 4: parquet poly-A restore on the soft-clipped variant ----
    if not skip_polya:
        print(f"\n→ restore poly-A tail from {PARQUET_TSV}")
        splice_cmd = [
            sys.executable,
            str(REPO / "scripts/validation_data/restore_polya_from_parquet.py"),
            "--in-bam", str(SOFTCLIP_OUT_BAM),
            "--parquet", str(PARQUET_TSV),
            "--out-bam", str(SOFTCLIP_OUT_BAM),
        ]
        res = subprocess.run(splice_cmd, capture_output=True)
        if res.returncode != 0:
            print(res.stderr.decode(errors="replace")[:2000], file=sys.stderr)
            return 2
        print(res.stdout.decode(errors="replace")[-500:])
        # Re-index after rewrite
        pysam.index(str(SOFTCLIP_OUT_BAM))


    # ---- Step 5: rebuild dorado-source BAM to stay in sync with validation_reads.bam ----
    # validation_reads_dorado_source.bam must cover the same 36 reads as
    # validation_reads.bam. When Cat1/Cat2 reads were replaced (v3.2.5) the
    # dorado_source was not updated; this step prevents that from recurring.
    print(f"\n\u2192 rebuild dorado-source BAM")
    dorado_cmd = [
        sys.executable,
        str(REPO / "scripts/validation_data/build_dorado_source.py"),
    ]
    res = subprocess.run(dorado_cmd, capture_output=True)
    print(res.stdout.decode(errors="replace"))
    if res.returncode != 0:
        print(res.stderr.decode(errors="replace")[:2000], file=sys.stderr)
        return 2

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
