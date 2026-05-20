#!/usr/bin/env python3
"""Integration test for QuantSeq REV walkback via the full 'rectify correct' pipeline.

Validates that the Module 2E walkback wiring introduced in NEW-075 fires correctly
when 'rectify correct --short-read --dT-primed-cDNA' is run on a synthetic
QuantSeq REV BAM.

Genome layout (all on chrI, 0-based):
  [0,   1000): N's             — padding
  [1000, 1010): GCTAGCGTAC    — gene body (non-A)
  [1010, 1020): AAAAAAAAAA    — downstream A-tract (10 A's)
  [1020, 2020): N's             — padding

Reads (all paired with the above):
  Group A (80 reads, is_reverse=True → gene '+'):
    Covers 1000–1019; seq = 'GCTAGCGTAC' + 'AAAAAAAAAA'
    3' end at 1019 (deep in A-tract) → walkback expected → corrected to 1009
  Group B (20 reads, is_reverse=True → gene '+'):
    Covers 1000–1009; seq = 'GCTAGCGTAC'
    3' end at 1009, terminal C matches genome C (non-stop) → Case 1 gate → no walkback
"""
from __future__ import annotations

import csv
import subprocess
import sys
import textwrap
from pathlib import Path

import pysam
import pytest


# ---------------------------------------------------------------------------
# Genome helpers
# ---------------------------------------------------------------------------
GENE_BODY = "GCTAGCGTAC"   # 10 bp, non-A
A_TRACT    = "AAAAAAAAAA"  # 10 bp, pure A
PADDING    = "N" * 1000

GENOME_SEQ = PADDING + GENE_BODY + A_TRACT + PADDING   # 2020 bp total
CHROM = "chrI"
CHROM_LEN = len(GENOME_SEQ)

# Reference positions (0-based)
GENE_BODY_START = 1000
GENE_BODY_END   = 1009   # inclusive, 0-based
A_TRACT_START   = 1010
A_TRACT_END     = 1019   # inclusive, 0-based


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def integration_dir(tmp_path_factory) -> Path:
    """Write genome FASTA, BAM, and return the work directory."""
    d = tmp_path_factory.mktemp("quantseq_rev_integration")

    # --- genome FASTA -------------------------------------------------------
    genome_fa = d / "genome.fa"
    genome_fa.write_text(f">{CHROM}\n{GENOME_SEQ}\n")

    # --- BAM ----------------------------------------------------------------
    unsorted_bam = d / "unsorted.bam"
    sorted_bam   = d / "reads.bam"

    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": CHROM, "LN": CHROM_LEN}],
    })

    with pysam.AlignmentFile(str(unsorted_bam), "wb", header=header) as bam:
        # Group A: poly-A over-extension → walkback expected
        for i in range(80):
            r = pysam.AlignedSegment(header)
            r.query_name           = f"overext_{i:04d}"
            r.reference_name       = CHROM
            r.reference_start      = GENE_BODY_START
            r.cigartuples          = [(0, 20)]
            r.is_reverse           = True
            r.is_unmapped          = False
            r.is_secondary         = False
            r.is_supplementary     = False
            r.mapping_quality      = 60
            r.query_sequence       = GENE_BODY + A_TRACT
            bam.write(r)

        # Group B: clean non-A 3' end → no walkback (Case 1 gate)
        for i in range(20):
            r = pysam.AlignedSegment(header)
            r.query_name           = f"clean_{i:04d}"
            r.reference_name       = CHROM
            r.reference_start      = GENE_BODY_START
            r.cigartuples          = [(0, 10)]
            r.is_reverse           = True
            r.is_unmapped          = False
            r.is_secondary         = False
            r.is_supplementary     = False
            r.mapping_quality      = 60
            r.query_sequence       = GENE_BODY
            bam.write(r)

    pysam.sort("-o", str(sorted_bam), str(unsorted_bam))
    pysam.index(str(sorted_bam))

    return d


@pytest.fixture(scope="module")
def pipeline_output(integration_dir: Path) -> list[dict]:
    """Run 'rectify correct --short-read --dT-primed-cDNA' and return TSV rows."""
    bam_path    = integration_dir / "reads.bam"
    genome_fa   = integration_dir / "genome.fa"
    out_tsv     = integration_dir / "corrected.tsv"

    cmd = [
        sys.executable, "-m", "rectify.cli", "correct",
        str(bam_path),
        "--genome", str(genome_fa),
        "--short-read",
        "--dT-primed-cDNA",
        "-o", str(out_tsv),
        "--emit-merged-tsv",  # Commit B: default is manifest-only; keep merged TSV for test assertions
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(
            f"rectify correct failed (rc={result.returncode}):\n"
            f"STDOUT:\n{result.stdout[-2000:]}\n"
            f"STDERR:\n{result.stderr[-2000:]}"
        )

    rows: list[dict] = []
    with open(out_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
class TestQuantSeqRevIntegration:
    """Full pipeline integration checks for QuantSeq REV (NEW-075)."""

    def test_pipeline_produces_output(self, pipeline_output):
        """Pipeline runs and writes 100 rows (one per synthetic read)."""
        assert len(pipeline_output) == 100

    def test_walkback_fires_for_overextended_reads(self, pipeline_output):
        """Group A (poly-A over-extension): polya_walkback_readgenome fires."""
        overext = [r for r in pipeline_output if r["read_id"].startswith("overext_")]
        assert len(overext) == 80

        walkback_applied = [
            r for r in overext
            if "polya_walkback_readgenome" in r.get("correction_applied", "")
        ]
        # All 80 over-extended reads should fire the walkback
        assert len(walkback_applied) == 80, (
            f"Expected 80 walkback corrections, got {len(walkback_applied)}. "
            f"Sample correction_applied values: "
            f"{[r.get('correction_applied') for r in overext[:5]]}"
        )

    def test_walkback_corrects_position(self, pipeline_output):
        """Group A: corrected 3' end moves from A-tract (1019) back to gene body (1009)."""
        overext = [r for r in pipeline_output if r["read_id"].startswith("overext_")]
        for row in overext:
            if "polya_walkback_readgenome" in row.get("correction_applied", ""):
                corrected = int(row["corrected_3prime"])
                original  = int(row["original_3prime"])
                assert original == A_TRACT_END, f"Expected original=1019, got {original}"
                assert corrected == GENE_BODY_END, (
                    f"Expected corrected=1009, got {corrected}"
                )

    def test_strand_is_gene_strand_not_bam_strand(self, pipeline_output):
        """is_reverse=True reads should report gene '+' strand, not BAM '-' strand."""
        for row in pipeline_output:
            assert row["strand"] == "+", (
                f"Read {row['read_id']}: expected strand='+', got '{row['strand']}'"
            )

    def test_clean_reads_not_walked_back(self, pipeline_output):
        """Group B (clean non-A end): no walkback applied."""
        clean = [r for r in pipeline_output if r["read_id"].startswith("clean_")]
        assert len(clean) == 20

        for row in clean:
            applied = row.get("correction_applied", "")
            assert "polya_walkback_readgenome" not in applied, (
                f"Clean read {row['read_id']} unexpectedly had walkback: {applied}"
            )

    def test_overall_walkback_rate_substantially_above_zero(self, pipeline_output):
        """80% of reads should have walkback applied — far above the legacy ~8% baseline."""
        walkback_count = sum(
            1 for r in pipeline_output
            if "polya_walkback_readgenome" in r.get("correction_applied", "")
        )
        rate = walkback_count / len(pipeline_output)
        assert rate >= 0.70, (
            f"Walkback rate {rate:.1%} is below 70% — expected ~80% for this synthetic dataset. "
            f"This suggests the NEW-075 wiring is not firing."
        )
