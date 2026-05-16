#!/usr/bin/env python3
"""End-to-end chain canary: correct-cdna → align → cdna-analyze.

This is the regression test that protects against ``minimap2 -y`` tag-pass-
through regressions identified in HANDOFF.md. The chain emits per-cluster
SAM-format tags into the FASTQ comment line (correct-cdna), then
``rectify align`` propagates them via ``minimap2 -y`` into the rectified
BAM, and ``cdna-analyze`` re-reads them on the post-align coordinates.

Specifically: if the FASTQ comment is space-separated instead of TAB-
separated, minimap2 collapses everything into one giant Z-string aux tag,
and every downstream consumer breaks. See
``feedback_minimap2_y_tab_sep.md`` in Kevin's memory index.

This test runs only when ALL of:
  - ``pyabpoa`` is importable (needed for correct-cdna's POA consensus),
  - ``minimap2`` and ``samtools`` are on PATH,
  - the chrI fixture files exist on disk.
"""
from __future__ import annotations

import importlib
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Prereq checks
# ---------------------------------------------------------------------------
_CHRI_BAM = Path("/Users/kevinroy/work/ont_cdna/test_data/wt_rep1.chrI.bam")
_GFF = Path("/Users/kevinroy/work/ont_cdna/test_data/saccharomyces_cerevisiae_R64-5-1_20240529.gff")
_FASTA = Path("/Users/kevinroy/work/ont_cdna/test_data/S288C_reference_sequence_R64-5-1_20240529.fsa")


def _have_pyabpoa() -> bool:
    try:
        importlib.import_module("pyabpoa")
        return True
    except ImportError:
        return False


def _have_minimap2() -> bool:
    return shutil.which("minimap2") is not None


def _have_samtools() -> bool:
    return shutil.which("samtools") is not None


def _have_fixture() -> bool:
    return _CHRI_BAM.exists() and _GFF.exists() and _FASTA.exists()


def _missing_prereqs() -> list[str]:
    missing = []
    if not _have_pyabpoa():
        missing.append("pyabpoa (POA consensus dependency)")
    if not _have_minimap2():
        missing.append("minimap2 binary on PATH")
    if not _have_samtools():
        missing.append("samtools binary on PATH")
    if not _have_fixture():
        missing.append(
            f"chrI test fixture (looked at {_CHRI_BAM.parent})"
        )
    return missing


_SKIP_REASON = (
    "chain canary requires: " + ", ".join(_missing_prereqs())
    if _missing_prereqs()
    else None
)


@pytest.mark.slow
@pytest.mark.skipif(_SKIP_REASON is not None, reason=_SKIP_REASON or "prereqs missing")
def test_chain_correct_cdna_to_analyze_preserves_tags(tmp_path):
    """Three-stage chain: correct-cdna → align → cdna-analyze.

    Asserts that the post-align consensus BAM carries every per-cluster
    SAM tag that ``cdna-analyze`` consumes: XU, XO, XC, XR, XT, XY, XF, XB.

    A missing tag in the post-align BAM almost always means minimap2's
    ``-y`` passthrough is collapsing space-separated FASTQ comments into
    one Z-string. The error message points at that specific failure mode.
    """
    # ── correct-cdna ─────────────────────────────────────────────────────
    correct_out = tmp_path / "correct"
    rc = subprocess.run(
        [sys.executable, "-m", "rectify", "correct-cdna",
         str(_CHRI_BAM), "--out", str(correct_out),
         "--reference", str(_FASTA)],
        capture_output=True, text=True, timeout=900,
    )
    assert rc.returncode == 0, (
        f"correct-cdna failed (rc={rc.returncode}):\n"
        f"STDERR: {rc.stderr[-1500:]}"
    )
    fastq = correct_out / "stage1_consensus.fastq.gz"
    assert fastq.exists(), f"correct-cdna did not emit {fastq}"

    # ── align ────────────────────────────────────────────────────────────
    align_out = tmp_path / "align"
    rc = subprocess.run(
        [sys.executable, "-m", "rectify", "align",
         str(fastq), "--genome", str(_FASTA),
         "-o", str(align_out), "--aligners", "minimap2",
         "--prefix", "stage1"],
        capture_output=True, text=True, timeout=1800,
    )
    assert rc.returncode == 0, (
        f"align failed (rc={rc.returncode}):\n"
        f"STDERR: {rc.stderr[-1500:]}"
    )
    aligned_bam = align_out / "stage1.rectified.md.bam"
    if not aligned_bam.exists():
        aligned_bam = align_out / "stage1.rectified.bam"
    assert aligned_bam.exists(), (
        f"no rectified BAM in {align_out}: {sorted(p.name for p in align_out.iterdir())}"
    )

    # ── Tag passthrough check ────────────────────────────────────────────
    # Read the first ~200 mapped records and assert each REQUIRED tag is
    # present on at least one of them. (One record may legitimately be
    # missing one optional tag — e.g. XB for a single-strand cluster —
    # but ANY record carrying the tag set proves the -y plumbing works.)
    import pysam as _pysam
    required = {"XU", "XO", "XC", "XR", "XT", "XY", "XF", "XB"}
    seen: set[str] = set()
    n_checked = 0
    with _pysam.AlignmentFile(str(aligned_bam), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped:
                continue
            n_checked += 1
            have = {t for t, _ in rec.tags}
            seen.update(required & have)
            if seen == required or n_checked >= 200:
                break
    missing = sorted(required - seen)
    assert not missing, (
        "Tag passthrough broken in rectify align (minimap2 -y).\n"
        f"  Records checked: {n_checked}\n"
        f"  Required tags missing from ALL checked records: {missing}\n"
        f"  Tags present: {sorted(seen)}\n"
        f"Most likely cause: the FASTQ comment formatter is using a "
        "SPACE separator. minimap2 -y propagates FASTQ comments verbatim, "
        "but the SAM parser expects TAB-separated tags — a space-separated "
        "comment collapses everything into ONE Z-string aux. Check "
        "write_stage1_fastq() in rectify/core/cdna/io.py: the tag-joiner "
        "MUST use '\\t'.join(...), not ' '.join(...)."
    )

    # ── cdna-analyze ─────────────────────────────────────────────────────
    analyze_out = tmp_path / "analyze"
    rc = subprocess.run(
        [sys.executable, "-m", "rectify", "cdna-analyze",
         str(aligned_bam), "-o", str(analyze_out),
         "--gff", str(_GFF), "--reference", str(_FASTA)],
        capture_output=True, text=True, timeout=600,
    )
    assert rc.returncode == 0, (
        f"cdna-analyze failed (rc={rc.returncode}):\n"
        f"STDERR: {rc.stderr[-1500:]}"
    )
    for f in ("clusters.tsv", "isoforms.tsv", "t1t2_pairs.tsv",
              "consensus_tagged.bam"):
        assert (analyze_out / f).exists(), (
            f"cdna-analyze missing output {f}"
        )

    # cluster count sanity band — per HANDOFF.md item 7
    n_clusters = sum(1 for _ in (analyze_out / "clusters.tsv").open()) - 1
    assert 24000 < n_clusters < 30000, (
        f"clusters.tsv: {n_clusters} rows outside [24000, 30000] band — "
        "may indicate Stage-1 dedup regression"
    )
    n_isoforms = sum(1 for _ in (analyze_out / "isoforms.tsv").open()) - 1
    assert 3000 < n_isoforms < 4500, (
        f"isoforms.tsv: {n_isoforms} rows outside [3000, 4500] band — "
        "may indicate isoform clustering tolerance regression"
    )
