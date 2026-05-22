"""Smoke test: write_corrected_bam_parallel produces output equivalent to
the existing single-threaded write_corrected_bam on a small fixture.

Uses the bundled 36-read validation BAM + corrected_reads.tsv.

Run with:
    pytest tests/test_bam_writer_parallel_smoke.py -v
"""

import os
import tempfile
from pathlib import Path
from typing import Dict

import pytest
import pysam
import pandas as pd

from rectify.core.bam.bam_writer import write_corrected_bam
from rectify.core.bam.bam_writer_parallel import write_corrected_bam_parallel
from rectify.core.consensus.corrected_consensus import (
    _read_hp_edit_distances,
    _read_hp_edit_distances_from_raw_bam,
    merge_corrected_tsvs,
    write_corrected_consensus_bam,
)
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


def test_lazy_hp_edit_distance_matches_materialized_corrected_bam(
    fixture_bam, fixture_tsv, genome, legacy_bam
):
    """Lazy raw-BAM scoring must match scoring from a materialized corrected BAM."""
    materialized = _read_hp_edit_distances(legacy_bam, genome=genome)
    lazy = _read_hp_edit_distances_from_raw_bam(
        fixture_bam,
        Path(fixture_tsv),
        genome=genome,
        strict=True,
    )
    assert lazy == materialized


def test_lazy_merge_matches_materialized_corrected_bam_selection(
    fixture_bam, fixture_tsv, genome, legacy_bam, tmp_path
):
    """Raw-BAM lazy HP scoring should choose the same winners as corrected BAM scoring."""
    per_aligner_tsvs = {
        "minimap2": Path(fixture_tsv),
        "uLTRA": Path(fixture_tsv),
    }

    materialized_out = tmp_path / "materialized.tsv"
    materialized_summary = tmp_path / "materialized_summary.tsv"
    merge_corrected_tsvs(
        per_aligner_tsvs,
        materialized_out,
        materialized_summary,
        per_aligner_corrected_bams={
            "minimap2": legacy_bam,
            "uLTRA": legacy_bam,
        },
        genome=genome,
    )

    lazy_out = tmp_path / "lazy.tsv"
    lazy_summary = tmp_path / "lazy_summary.tsv"
    merge_corrected_tsvs(
        per_aligner_tsvs,
        lazy_out,
        lazy_summary,
        per_aligner_raw_bams={
            "minimap2": fixture_bam,
            "uLTRA": fixture_bam,
        },
        genome=genome,
    )

    materialized_df = pd.read_csv(materialized_out, sep="\t")
    lazy_df = pd.read_csv(lazy_out, sep="\t")
    assert list(lazy_df["read_id"]) == list(materialized_df["read_id"])
    assert list(lazy_df["winning_aligner"]) == list(materialized_df["winning_aligner"])


def test_write_corrected_consensus_bam_from_raw_matches_legacy_single_aligner(
    fixture_bam, fixture_tsv, genome, legacy_bam, tmp_path
):
    """Final lazy consensus BAM should match legacy corrected BAM for one winner."""
    merged_tsv = tmp_path / "merged.tsv"
    merged = pd.read_csv(fixture_tsv, sep="\t")
    merged["winning_aligner"] = "minimap2"
    merged.to_csv(merged_tsv, sep="\t", index=False)

    out_bam = tmp_path / "corrected_consensus.bam"
    stats = write_corrected_consensus_bam(
        {"minimap2": fixture_bam},
        {"minimap2": Path(fixture_tsv)},
        merged_tsv,
        out_bam,
        genome=genome,
        threads=1,
    )

    assert stats["written"] == len(merged["read_id"].unique())
    assert out_bam.exists()
    assert Path(str(out_bam) + ".bai").exists()
    assert_bams_equivalent(legacy_bam, str(out_bam))


# ---------------------------------------------------------------------------
# Helpers shared by the two new strengthened tests
# ---------------------------------------------------------------------------

def _build_perturbed_tsv(fixture_tsv_path: str, out_tsv: Path) -> None:
    """Write a perturbed copy of the fixture corrected_reads.tsv.

    For the first 18 reads (half of 36), shift ``corrected_3prime`` toward the
    transcript interior by 50 bp.  The resulting per-aligner corrected BAM will
    have a different hard-clip position on those reads, producing different
    HP-edit-distance values and (crucially) different per-read winner decisions
    when compared to the unperturbed fixture.

    The fixture carries a stale ``winning_aligner`` column from a prior 5-aligner
    run, but merge_corrected_tsvs now drops that column on load (the fix for the
    rename-collision bug), so no stripping is needed here.
    """
    tsv = pd.read_csv(fixture_tsv_path, sep="\t")
    for i in range(18):
        row = tsv.iloc[i]
        if row["strand"] == "+":
            tsv.loc[tsv.index[i], "corrected_3prime"] = int(row["corrected_3prime"]) - 50
        else:
            tsv.loc[tsv.index[i], "corrected_3prime"] = int(row["corrected_3prime"]) + 50
    tsv.to_csv(out_tsv, sep="\t", index=False)


# ---------------------------------------------------------------------------
# P1 fix: non-degenerate two-aligner winner selection
# ---------------------------------------------------------------------------

def test_lazy_merge_distinct_aligners_select_same_winners(
    fixture_bam, fixture_tsv, genome, tmp_path
):
    """Strengthen test_lazy_merge_matches_materialized_corrected_bam_selection.

    The original test passes the SAME BAM and TSV for both aligners, which makes
    every read a perfect tie — winner selection is degenerate and cannot catch
    a materialise-vs-lazy divergence where only one aligner's CIGAR was mutated.

    This test uses TWO DISTINCT (BAM, TSV) pairs:
    - minimap2: the unmodified fixture (original corrected BAM)
    - uLTRA: a perturbed fixture with corrected_3prime shifted -50/+50 bp for
      the first 18 reads (direction inward toward the transcript body)

    The perturbation produces different HP-edit-distances on those 18 reads,
    forcing real per-read winner selection where minimap2 and uLTRA each win
    some reads.

    Assert:
    1. At least one read has each aligner as the winner (non-trivial selection).
    2. Lazy (raw-BAM) and materialized (pre-built corrected BAM) merge agree on
       the per-read winner for every read.
    """
    # Build perturbed TSV and corresponding corrected BAMs.
    perturbed_tsv = tmp_path / "perturbed.tsv"
    original_bam = tmp_path / "original.bam"
    perturbed_bam = tmp_path / "perturbed.bam"

    _build_perturbed_tsv(fixture_tsv, perturbed_tsv)

    # Pass the fixture TSV directly — merge_corrected_tsvs drops the stale
    # winning_aligner column on load, so no pre-stripping is needed.
    write_corrected_bam(fixture_bam, fixture_tsv, str(original_bam), genome=genome)
    write_corrected_bam(fixture_bam, str(perturbed_tsv), str(perturbed_bam), genome=genome)

    per_aligner_tsvs = {
        "minimap2": Path(fixture_tsv),
        "uLTRA": perturbed_tsv,
    }

    # Materialized path: score from pre-built corrected BAMs.
    mat_out = tmp_path / "materialized.tsv"
    merge_corrected_tsvs(
        per_aligner_tsvs,
        mat_out,
        None,
        per_aligner_corrected_bams={
            "minimap2": str(original_bam),
            "uLTRA": str(perturbed_bam),
        },
        genome=genome,
    )

    # Lazy path: apply corrections in memory from raw BAM + TSV.
    lazy_out = tmp_path / "lazy.tsv"
    merge_corrected_tsvs(
        per_aligner_tsvs,
        lazy_out,
        None,
        per_aligner_raw_bams={
            "minimap2": fixture_bam,
            "uLTRA": fixture_bam,
        },
        genome=genome,
    )

    mat_df = pd.read_csv(mat_out, sep="\t")
    lazy_df = pd.read_csv(lazy_out, sep="\t")

    # Guard: both aligners must win at least one read (non-degenerate selection).
    mat_winners = set(mat_df["winning_aligner"].unique())
    assert "minimap2" in mat_winners, (
        "minimap2 never wins — perturbation is too small or fixture changed"
    )
    assert "uLTRA" in mat_winners, (
        "uLTRA never wins — perturbation is too small or fixture changed"
    )

    # Core assertion: lazy and materialized select the same per-read winner.
    assert list(lazy_df["read_id"]) == list(mat_df["read_id"]), (
        "read_id order differs between lazy and materialized"
    )
    assert list(lazy_df["winning_aligner"]) == list(mat_df["winning_aligner"]), (
        "lazy and materialized chose different per-read winners on at least one read"
    )


# ---------------------------------------------------------------------------
# P1 fix: consensus BAM with two distinct aligners
# ---------------------------------------------------------------------------

def test_consensus_bam_per_aligner_writes_correctly(
    fixture_bam, fixture_tsv, genome, tmp_path
):
    """Extend test_write_corrected_consensus_bam_from_raw_matches_legacy_single_aligner.

    Using the same perturbation strategy (minimap2 = original fixture, uLTRA =
    corrected_3prime-shifted fixture), this test verifies that the final consensus
    BAM produced by write_corrected_consensus_bam contains:

    1. Exactly N=36 written reads (all reads accounted for).
    2. For every minimap2-winning read: the consensus CIGAR matches the
       minimap2-specific corrected BAM for that read.
    3. For every uLTRA-winning read: the consensus CIGAR matches the
       uLTRA-specific corrected BAM for that read.
    4. The index file is present alongside the consensus BAM.

    This exercises the two-aligner code path in write_corrected_consensus_bam,
    which iterates aligner_to_reads and must pull the correct raw BAM for each
    winning aligner.
    """
    perturbed_tsv = tmp_path / "perturbed.tsv"
    original_bam = tmp_path / "original.bam"
    perturbed_bam = tmp_path / "perturbed.bam"

    _build_perturbed_tsv(fixture_tsv, perturbed_tsv)

    # Pass the fixture TSV directly — merge_corrected_tsvs drops the stale
    # winning_aligner column on load, so no pre-stripping is needed.
    write_corrected_bam(fixture_bam, fixture_tsv, str(original_bam), genome=genome)
    write_corrected_bam(fixture_bam, str(perturbed_tsv), str(perturbed_bam), genome=genome)

    per_aligner_tsvs: Dict[str, Path] = {
        "minimap2": Path(fixture_tsv),
        "uLTRA": perturbed_tsv,
    }

    # Build merged TSV using materialized BAMs for HP-ED scoring.
    merged_tsv = tmp_path / "merged.tsv"
    merge_corrected_tsvs(
        per_aligner_tsvs,
        merged_tsv,
        None,
        per_aligner_corrected_bams={
            "minimap2": str(original_bam),
            "uLTRA": str(perturbed_bam),
        },
        genome=genome,
    )

    merged_df = pd.read_csv(merged_tsv, sep="\t")
    mm2_wins = set(merged_df.loc[merged_df["winning_aligner"] == "minimap2", "read_id"])
    ultra_wins = set(merged_df.loc[merged_df["winning_aligner"] == "uLTRA", "read_id"])

    # Guard: both aligners must win at least one read.
    assert len(mm2_wins) > 0, "minimap2 must win at least one read"
    assert len(ultra_wins) > 0, "uLTRA must win at least one read"

    # Write the consensus BAM.
    consensus_bam = tmp_path / "consensus.bam"
    stats = write_corrected_consensus_bam(
        {"minimap2": fixture_bam, "uLTRA": fixture_bam},
        per_aligner_tsvs,
        merged_tsv,
        consensus_bam,
        genome=genome,
        threads=1,
    )

    # 1. All 36 reads written.
    n_unique = len(merged_df["read_id"].unique())
    assert stats["written"] == n_unique, (
        f"Expected {n_unique} written reads; got {stats['written']}"
    )
    assert consensus_bam.exists(), "Consensus BAM not created"
    assert Path(str(consensus_bam) + ".bai").exists(), "Consensus BAM index not created"

    # 2 & 3. Per-read CIGAR in the consensus must match the correct per-aligner BAM.
    def _read_cigars(bam_path: str) -> Dict[str, str]:
        cigars = {}
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam:
                if not (read.is_unmapped or read.is_secondary or read.is_supplementary):
                    cigars[read.query_name] = read.cigarstring
        return cigars

    consensus_cigars = _read_cigars(str(consensus_bam))
    mm2_cigars = _read_cigars(str(original_bam))
    ultra_cigars = _read_cigars(str(perturbed_bam))

    # minimap2-winning reads: consensus CIGAR must match original (minimap2) corrected BAM.
    for read_id in mm2_wins:
        if read_id in consensus_cigars and read_id in mm2_cigars:
            assert consensus_cigars[read_id] == mm2_cigars[read_id], (
                f"minimap2-winning read {read_id}: "
                f"consensus={consensus_cigars[read_id]!r} "
                f"mm2={mm2_cigars[read_id]!r}"
            )

    # uLTRA-winning reads: consensus CIGAR must match perturbed (uLTRA) corrected BAM.
    for read_id in ultra_wins:
        if read_id in consensus_cigars and read_id in ultra_cigars:
            assert consensus_cigars[read_id] == ultra_cigars[read_id], (
                f"uLTRA-winning read {read_id}: "
                f"consensus={consensus_cigars[read_id]!r} "
                f"ultra={ultra_cigars[read_id]!r}"
            )


# ---------------------------------------------------------------------------
# Regression: stale winning_aligner column in per-aligner TSV is ignored
# ---------------------------------------------------------------------------

def test_merge_corrected_tsvs_ignores_stale_winning_aligner_column(
    fixture_bam, fixture_tsv, genome, tmp_path
):
    """merge_corrected_tsvs must not let a stale winning_aligner column from a
    per-aligner input TSV shadow the freshly-computed winner.

    The bundled fixture corrected_reads.tsv was generated from a 5-aligner run
    and carries legacy aligner names (e.g. 'deSALT', 'mapPacBio') in its
    winning_aligner column.  Before the fix in _load_tsv, passing that TSV
    directly to merge_corrected_tsvs caused a pandas rename collision:
    all_df already had a 'winning_aligner' column (stale), and renaming
    '_winning_aligner' → 'winning_aligner' produced a duplicate column.
    to_csv wrote both columns under the same header; callers reading the
    output saw the stale values.

    This test calls merge_corrected_tsvs WITHOUT pre-stripping the column and
    asserts the output winners are the freshly-computed HP-ED winners, not the
    stale fixture values.
    """
    # Confirm the fixture actually has the stale winning_aligner column so the
    # test is not vacuously testing a stripped fixture.
    raw_fixture = pd.read_csv(fixture_tsv, sep="\t")
    assert "winning_aligner" in raw_fixture.columns, (
        "Fixture TSV has no winning_aligner column — update this test if the "
        "fixture was regenerated without one."
    )
    stale_values = set(raw_fixture["winning_aligner"].dropna().unique())

    # Build a perturbed second TSV so winner selection is non-degenerate.
    perturbed_tsv = tmp_path / "perturbed.tsv"
    _build_perturbed_tsv(fixture_tsv, perturbed_tsv)

    # Build reference (clean) output using the pre-stripped path that already
    # works, to establish the expected winners.
    clean_tsv = tmp_path / "clean.tsv"
    raw_fixture.drop(columns=["winning_aligner"]).to_csv(clean_tsv, sep="\t", index=False)

    ref_out = tmp_path / "ref.tsv"
    merge_corrected_tsvs(
        {"minimap2": clean_tsv, "uLTRA": perturbed_tsv},
        ref_out,
        per_aligner_raw_bams={
            "minimap2": fixture_bam,
            "uLTRA": fixture_bam,
        },
        genome=genome,
    )
    ref_df = pd.read_csv(ref_out, sep="\t")

    # Now call merge with the ORIGINAL fixture (stale column present) — no stripping.
    test_out = tmp_path / "test.tsv"
    merge_corrected_tsvs(
        {"minimap2": Path(fixture_tsv), "uLTRA": perturbed_tsv},
        test_out,
        per_aligner_raw_bams={
            "minimap2": fixture_bam,
            "uLTRA": fixture_bam,
        },
        genome=genome,
    )
    test_df = pd.read_csv(test_out, sep="\t")

    # The output must have exactly one winning_aligner column (no duplicates).
    assert list(test_df.columns).count("winning_aligner") == 1, (
        "Output TSV has duplicate winning_aligner columns — rename collision not fixed."
    )

    # No stale aligner name should appear in the output.
    output_winners = set(test_df["winning_aligner"].dropna().unique())
    stale_in_output = stale_values - {"minimap2", "uLTRA"}
    assert not (output_winners & stale_in_output), (
        f"Stale winning_aligner values from fixture appeared in output: "
        f"{output_winners & stale_in_output}"
    )

    # The per-read winners must match the reference (stripped) run exactly.
    assert list(test_df["read_id"]) == list(ref_df["read_id"]), (
        "read_id order differs between stale-column run and reference run"
    )
    assert list(test_df["winning_aligner"]) == list(ref_df["winning_aligner"]), (
        "Per-read winners differ between stale-column run and reference run — "
        "stale winning_aligner column is still influencing the output."
    )
