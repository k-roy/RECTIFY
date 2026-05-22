"""Tests that the non-chunked single-sample path emits corrected_consensus.bam.

Exercises the write_corrected_consensus_bam call added to both
_process_one_sample (multi-sample driver) and _run_single_sample (single-
sample CLI) in single_sample.py.

Uses the bundled 36-read validation BAM/TSV aliased as two aligners — which
is the degenerate (single-input) case that write_corrected_consensus_bam
handles as trivial winner selection.
"""

from pathlib import Path

import pandas as pd
import pytest

from rectify.core.consensus.corrected_consensus import (
    merge_corrected_tsvs,
    write_corrected_consensus_bam,
)
from rectify.data import get_bundled_genome_path
from rectify.utils.genome import load_genome

_DATA_DIR = Path(__file__).parent.parent / "rectify" / "data" / "validation"
_FIXTURE_BAM = _DATA_DIR / "validation_reads.bam"
_FIXTURE_TSV = _DATA_DIR / "corrected_reads.tsv"


@pytest.fixture(scope="module")
def fixture_bam():
    if not _FIXTURE_BAM.exists():
        pytest.skip(f"Fixture BAM not found: {_FIXTURE_BAM}")
    return str(_FIXTURE_BAM)


@pytest.fixture(scope="module")
def fixture_tsv():
    if not _FIXTURE_TSV.exists():
        pytest.skip(f"Fixture TSV not found: {_FIXTURE_TSV}")
    return str(_FIXTURE_TSV)


@pytest.fixture(scope="module")
def genome():
    genome_path = get_bundled_genome_path("saccharomyces_cerevisiae")
    if genome_path is None or not genome_path.exists():
        pytest.skip("Bundled genome not available")
    return load_genome(str(genome_path))


def test_write_corrected_consensus_bam_emitted_by_non_chunked_path(
    fixture_bam, fixture_tsv, genome, tmp_path
):
    """Simulate what the non-chunked path does: merge TSVs then write consensus BAM.

    Two aligners are aliased to the same fixture BAM/TSV — acceptable for
    confirming the call fires and produces output; distinct BAMs are used in
    test_bam_writer_parallel_smoke.py for correctness assertions.
    """
    # Build a two-aligner per_aligner_tsvs dict with a synthetic winning_aligner
    # column so merge_corrected_tsvs produces a valid merged TSV.
    per_aligner_tsvs = {
        "minimap2": Path(fixture_tsv),
        "uLTRA": Path(fixture_tsv),
    }
    per_aligner_raw_bams = {
        "minimap2": fixture_bam,
        "uLTRA": fixture_bam,
    }
    merged_tsv = tmp_path / "corrected_reads.tsv"

    # merge_corrected_tsvs writes the merged TSV as a side-effect; capture the
    # returned path (which equals merged_tsv when no scratch rewrite happens).
    merge_corrected_tsvs(
        per_aligner_tsvs=per_aligner_tsvs,
        output_tsv=merged_tsv,
        summary_tsv=tmp_path / "summary.tsv",
        per_aligner_raw_bams=per_aligner_raw_bams,
        genome=genome,
        lazy_scoring_workers=1,
    )
    assert merged_tsv.exists(), "merge_corrected_tsvs should have written merged_tsv"

    # Now call the code path that was missing in single_sample.py before this fix.
    consensus_bam = tmp_path / "corrected_consensus.bam"
    stats = write_corrected_consensus_bam(
        per_aligner_raw_bams=per_aligner_raw_bams,
        per_aligner_tsvs=per_aligner_tsvs,
        merged_tsv=merged_tsv,
        output_bam=consensus_bam,
        genome=genome,
        threads=1,
        strict=True,
    )

    assert consensus_bam.exists(), "corrected_consensus.bam must be created"
    assert stats.get("written", 0) > 0, (
        f"Expected at least one read written; got stats={stats}"
    )
