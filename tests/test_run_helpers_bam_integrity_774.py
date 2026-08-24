"""Regression: `_validate_bam_integrity` must not raise NameError.

planning/774 — salvaged from an uncommitted H2 working tree. `helpers.py` called
`_subprocess.run(...)` at the samtools-quickcheck step but never imported
`subprocess`, so EVERY call that got past the `.bai` existence check raised
`NameError: name '_subprocess' is not defined`. The function is called from three
live sites on the `run-all` BAM-reuse path (`run/single_sample.py` :422 and :823,
`run/stages.py` :101), so the bug was reachable in production, not theoretical.
"""
import pathlib

import pytest

from rectify.core.commands.run import helpers
from rectify.core.commands.run.helpers import _validate_bam_integrity


def test_subprocess_is_imported():
    """The module-level name the quickcheck call depends on must exist."""
    assert hasattr(helpers, '_subprocess')
    assert helpers._subprocess.run is not None


def test_missing_bam_is_rejected(tmp_path: pathlib.Path):
    assert _validate_bam_integrity(tmp_path / 'nope.bam') is False


def test_bam_without_index_is_rejected(tmp_path: pathlib.Path):
    bam = tmp_path / 'x.bam'
    bam.write_bytes(b'\x1f\x8b')
    assert _validate_bam_integrity(bam) is False


def test_corrupt_bam_with_index_returns_false_not_nameerror(tmp_path: pathlib.Path):
    """The regression itself: with a .bai present, execution reaches the
    quickcheck call. Before the fix this raised NameError instead of returning."""
    bam = tmp_path / 'x.bam'
    bam.write_bytes(b'\x1f\x8b')          # gzip magic only — not a valid BAM
    (tmp_path / 'x.bam.bai').write_bytes(b'')
    try:
        result = _validate_bam_integrity(bam)
    except NameError as exc:                # pragma: no cover - the bug
        pytest.fail(f'_validate_bam_integrity raised NameError: {exc}')
    except FileNotFoundError:
        pytest.skip('samtools not on PATH')
    assert result is False
