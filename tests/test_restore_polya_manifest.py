"""Tests for restore-softclip corrected TSV input loading."""

from __future__ import annotations  # PEP 563: builtin-generic annotations run on py3.8
from pathlib import Path

import pytest

from rectify.core.commands.restore_polya_command import _load_winning_aligners


def _write_region_tsv(path: Path, rows: list[tuple[str, str]]) -> None:
    with path.open("w") as fh:
        fh.write("read_id\tchrom\tcorrected_3prime\twinning_aligner\n")
        for read_id, aligner in rows:
            fh.write(f"{read_id}\tchrI\t100\t{aligner}\n")


def test_load_winning_aligners_from_single_tsv(tmp_path):
    tsv = tmp_path / "corrected_reads.tsv"
    _write_region_tsv(tsv, [("read1", "minimap2"), ("read2", "deSALT")])

    assert _load_winning_aligners(str(tsv)) == {
        "read1": "minimap2",
        "read2": "deSALT",
    }


def test_load_winning_aligners_from_manifest(tmp_path):
    region_0 = tmp_path / "corrected_reads.region_000.tsv"
    region_1 = tmp_path / "corrected_reads.region_001.tsv"
    _write_region_tsv(region_0, [("read1", "minimap2")])
    _write_region_tsv(region_1, [("read2", "deSALT"), ("read3", "gapmm2")])

    manifest = tmp_path / "corrected_reads.manifest.tsv"
    manifest.write_text(
        "\t".join(["region_id", "chrom", "start", "end", "tsv_path", "n_rows", "sha256"])
        + "\n"
        + f"region_000\tchrI\t0\t100\t{region_0.name}\t1\tunused\n"
        + f"region_001\tchrI\t100\t200\t{region_1.name}\t2\tunused\n"
    )

    assert _load_winning_aligners(str(manifest)) == {
        "read1": "minimap2",
        "read2": "deSALT",
        "read3": "gapmm2",
    }


def test_load_winning_aligners_requires_winning_aligner(tmp_path):
    tsv = tmp_path / "corrected_reads.tsv"
    tsv.write_text("read_id\tchrom\nread1\tchrI\n")

    with pytest.raises(ValueError, match="winning_aligner"):
        _load_winning_aligners(str(tsv))
