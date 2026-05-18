"""Tests for ``scripts/validation_data/splice_polya_from_parquet.py``.

Asserts that the splice script re-attaches the trimmer-determined poly(A)
bases as a 3' soft-clip on the validation bundle's representative
``cat3_plus_2`` read (read_id ``79f61403-cf63-4522-b555-569590dc4304``).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pysam
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SCRIPT_PATH = (
    _REPO_ROOT / "scripts" / "validation_data" / "splice_polya_from_parquet.py"
)
_SOURCE_BAM = (
    _REPO_ROOT
    / "rectify"
    / "data"
    / "validation"
    / "rectified"
    / "rectified_pA_tail_soft_clipped.bam"
)
_METADATA_TSV = (
    _REPO_ROOT
    / "scripts"
    / "validation_data"
    / "rebuild_2026_05"
    / "trimmed"
    / "validation_reads_polya_trim_metadata.tsv"
)
_CAT3_PLUS_2 = "79f61403-cf63-4522-b555-569590dc4304"


def _load_script_module():
    """Import the script as a module from its file path (it lives outside the package)."""
    spec = importlib.util.spec_from_file_location(
        "splice_polya_from_parquet", _SCRIPT_PATH
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


_CIGAR_S = 4
_CIGAR_H = 5


def _strip_trailing_softclip(read: pysam.AlignedSegment, n: int) -> None:
    """Remove a trailing ``n``-base soft-clip in place (mirror of the splice).

    Used in the test fixture to simulate the regenerated-but-poly(A)-dropped
    state, so the test exercises the script regardless of whether the on-disk
    bundle has already been spliced.
    """
    cigar = list(read.cigartuples or [])
    if not cigar:
        return
    trailing_h = 0
    while cigar and cigar[-1][0] == _CIGAR_H:
        trailing_h += cigar.pop()[1]
    if not cigar or cigar[-1][0] != _CIGAR_S or cigar[-1][1] < n:
        return
    op, length = cigar.pop()
    if length > n:
        cigar.append((_CIGAR_S, length - n))
    if trailing_h:
        cigar.append((_CIGAR_H, trailing_h))
    seq = read.query_sequence
    quals = list(read.query_qualities or [])
    has_quals = len(quals) == len(seq)
    read.query_sequence = seq[:-n]
    read.cigartuples = cigar
    if has_quals:
        read.query_qualities = pysam.qualitystring_to_array(
            "".join(chr(min(q + 33, 126)) for q in quals[:-n])
        )


def _extract_single_read_bam(
    src_bam: Path, dst_bam: Path, read_id: str, strip_trailing_n: int = 0
) -> int:
    """Write a 1-read BAM (primary alignment only) for ``read_id``.

    If ``strip_trailing_n`` > 0, strips that many trailing soft-clip bases from
    the plus-strand read before writing (used to simulate the pre-splice state
    even when the bundle on disk has already been spliced).
    """
    count = 0
    with pysam.AlignmentFile(str(src_bam), "rb") as bam_in, pysam.AlignmentFile(
        str(dst_bam), "wb", template=bam_in
    ) as bam_out:
        for read in bam_in:
            if read.query_name != read_id:
                continue
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            if strip_trailing_n > 0 and not read.is_reverse:
                _strip_trailing_softclip(read, strip_trailing_n)
            bam_out.write(read)
            count += 1
    pysam.index(str(dst_bam))
    return count


@pytest.fixture(scope="module")
def parquet_row():
    """Return the trim metadata row for cat3_plus_2 as a dict."""
    assert _METADATA_TSV.exists(), (
        f"Validation-bundle trim metadata TSV not found: {_METADATA_TSV}"
    )
    df = pd.read_csv(_METADATA_TSV, sep="\t")
    rows = df[df["read_id"] == _CAT3_PLUS_2]
    assert len(rows) == 1, (
        f"Expected exactly one trim-metadata row for {_CAT3_PLUS_2}, got {len(rows)}"
    )
    return rows.iloc[0].to_dict()


def test_splice_cat3_plus_2(tmp_path, parquet_row):
    """End-to-end: extract a 1-read BAM, splice, verify the trailing S op and bases."""
    assert _SOURCE_BAM.exists(), f"Source BAM not found: {_SOURCE_BAM}"
    assert _SCRIPT_PATH.exists(), f"Splice script not found: {_SCRIPT_PATH}"

    polya_len = int(parquet_row["polya_len"])
    polya_seq = str(parquet_row["trimmed_3prime_seq"])
    assert len(polya_seq) == polya_len, (
        f"polya_len mismatch in TSV: declared {polya_len}, "
        f"actual seq len {len(polya_seq)}"
    )
    assert polya_len > 0
    # cat3_plus_2 is a plus-strand read per the validation bundle layout
    assert parquet_row["strand"] == "+"

    # Carve out a 1-read BAM. Strip any existing trailing soft-clip of length
    # polya_len so the test simulates the regenerated-but-poly(A)-dropped input
    # state regardless of whether the on-disk bundle has already been spliced.
    in_bam = tmp_path / "cat3_plus_2.in.bam"
    n_extracted = _extract_single_read_bam(
        _SOURCE_BAM, in_bam, _CAT3_PLUS_2, strip_trailing_n=polya_len
    )
    assert n_extracted == 1, (
        f"Expected exactly 1 primary alignment for {_CAT3_PLUS_2}, got {n_extracted}"
    )

    # Capture original (pre-splice) sequence + CIGAR for delta assertions
    with pysam.AlignmentFile(str(in_bam), "rb") as bam:
        orig = next(iter(bam))
    orig_seq = orig.query_sequence
    orig_cigar = list(orig.cigartuples or [])
    # Now the input must NOT carry a trailing polya-sized S op
    assert orig_cigar[-1][0] != 4 or orig_cigar[-1][1] < polya_len, (
        "After stripping, source still has a trailing S op of >= polya_len; "
        f"got CIGAR tail: {orig_cigar[-5:]}"
    )

    # Run the splice script via main()
    out_bam = tmp_path / "cat3_plus_2.out.bam"
    mod = _load_script_module()
    rc = mod.main([
        "--in-bam", str(in_bam),
        "--parquet", str(_METADATA_TSV),
        "--out-bam", str(out_bam),
    ])
    assert rc == 0
    assert out_bam.exists()

    # Verify the spliced read
    with pysam.AlignmentFile(str(out_bam), "rb") as bam:
        reads = [r for r in bam if r.query_name == _CAT3_PLUS_2
                 and not (r.is_secondary or r.is_supplementary)]
    assert len(reads) == 1
    spliced = reads[0]

    # CIGAR: trailing op (ignoring H) must be S of exact size polya_len
    spliced_cigar = list(spliced.cigartuples)
    trailing = next(
        (op_len for op_len in reversed(spliced_cigar) if op_len[0] != 5),
        None,
    )
    assert trailing is not None
    assert trailing[0] == 4, (
        f"Expected trailing soft-clip op (4), got CIGAR tail: {spliced_cigar[-5:]}"
    )
    assert trailing[1] == polya_len, (
        f"Expected trailing S op of length {polya_len}, got {trailing[1]}; "
        f"CIGAR tail: {spliced_cigar[-5:]}"
    )

    # Sequence length grew by exactly polya_len
    assert len(spliced.query_sequence) == len(orig_seq) + polya_len

    # The trailing bases match the parquet poly(A) sequence verbatim
    assert spliced.query_sequence[-polya_len:] == polya_seq, (
        f"Trailing bases mismatch: got {spliced.query_sequence[-polya_len:]!r}, "
        f"expected {polya_seq!r}"
    )

    # Prefix is unchanged
    assert spliced.query_sequence[:-polya_len] == orig_seq

    # reference_start is unchanged
    assert spliced.reference_start == orig.reference_start
