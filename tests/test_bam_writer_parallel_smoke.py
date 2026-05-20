"""Smoke test: write_corrected_bam_parallel produces output equivalent to
the existing single-threaded write_corrected_bam on a small fixture.

Uses the bundled 36-read validation BAM + corrected_reads.tsv.

Run with:
    pytest tests/test_bam_writer_parallel_smoke.py -v
"""

import os
import tempfile
from pathlib import Path

import pytest
import pysam

from rectify.core.bam.bam_writer import write_corrected_bam
from rectify.core.bam.bam_writer_parallel import write_corrected_bam_parallel
from rectify.utils.genome import load_genome
from rectify.data import get_bundled_genome_path

from tests.utils.bam_compare import assert_bams_equivalent

# ---------------------------------------------------------------------------
# Fixture paths (bundled in the repository)
# ---------------------------------------------------------------------------

_DATA_DIR = Path(__file__).parent.parent / "rectify" / "data" / "validation"
_FIXTURE_BAM = _DATA_DIR / "validation_reads.bam"
_FIXTURE_TSV = _DATA_DIR / "corrected_reads.tsv"


# ---------------------------------------------------------------------------
# pytest fixtures
# ---------------------------------------------------------------------------

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
    """Load the bundled S. cerevisiae genome (required for CIGAR surgery)."""
    genome_path = get_bundled_genome_path("saccharomyces_cerevisiae")
    if genome_path is None or not genome_path.exists():
        pytest.skip("Bundled genome not available")
    return load_genome(str(genome_path))


@pytest.fixture(scope="module")
def legacy_bam(fixture_bam, fixture_tsv, genome, tmp_path_factory):
    """Run the single-threaded legacy path once (module scope — reused across parametrize)."""
    tmp = tmp_path_factory.mktemp("legacy")
    out_path = str(tmp / "legacy.bam")
    write_corrected_bam(fixture_bam, fixture_tsv, out_path, genome=genome)
    return out_path


# ---------------------------------------------------------------------------
# Smoke test parametrized over n_threads
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_threads", [1, 2, 4])
def test_parallel_equivalent_to_legacy(
    n_threads, fixture_bam, fixture_tsv, genome, legacy_bam, tmp_path
):
    """write_corrected_bam_parallel with n_threads={1,2,4} must produce output
    QNAME-sort-equivalent to the legacy single-threaded write_corrected_bam.

    Checks:
    1. Output BAM exists.
    2. Read count equals legacy BAM read count.
    3. Every read is equivalent (CIGAR, FLAG, POS, MAPQ, all tags).
    """
    out_path = str(tmp_path / f"parallel_t{n_threads}.bam")

    stats = write_corrected_bam_parallel(
        fixture_bam,
        fixture_tsv,
        out_path,
        n_threads=n_threads,
        genome=genome,
    )

    # Basic sanity on returned stats dict.
    assert "n_regions" in stats
    assert "n_reads_out_total" in stats
    assert stats["n_regions"] >= 1, "Expected at least one region"

    # Output file must exist and be indexed.
    assert os.path.exists(out_path), f"Output BAM not created: {out_path}"
    assert os.path.exists(out_path + ".bai"), "Output BAM index not created"

    # Equivalence check vs legacy.
    assert_bams_equivalent(legacy_bam, out_path)
