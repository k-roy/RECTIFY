"""Tests for rectify.core.provenance.hashing."""

import hashlib
from pathlib import Path

import pytest

from rectify.core.provenance.hashing import normalized_config_hash, sha256_of_file


# ---------------------------------------------------------------------------
# sha256_of_file
# ---------------------------------------------------------------------------


def test_sha256_tiny_file(tmp_path):
    """SHA-256 of a known small file matches hashlib ground truth."""
    content = b"hello world\n"
    p = tmp_path / "tiny.txt"
    p.write_bytes(content)
    expected = hashlib.sha256(content).hexdigest()
    assert sha256_of_file(p) == expected


def test_sha256_empty_file(tmp_path):
    """Empty file hashes to the empty SHA-256."""
    p = tmp_path / "empty.txt"
    p.write_bytes(b"")
    expected = hashlib.sha256(b"").hexdigest()
    assert sha256_of_file(p) == expected


def test_sha256_large_file_mmap_path(tmp_path):
    """A 20 MB file exercises the mmap code path and matches hashlib."""
    content = b"A" * (20 * 1024 * 1024)  # 20 MiB
    p = tmp_path / "large.bin"
    p.write_bytes(content)
    expected = hashlib.sha256(content).hexdigest()
    assert sha256_of_file(p) == expected


def test_sha256_missing_file(tmp_path):
    """FileNotFoundError on a non-existent path."""
    with pytest.raises(FileNotFoundError):
        sha256_of_file(tmp_path / "does_not_exist.txt")


def test_sha256_accepts_str_path(tmp_path):
    """sha256_of_file accepts a str path as well as Path."""
    p = tmp_path / "str_path.txt"
    p.write_bytes(b"test")
    result_path = sha256_of_file(p)
    result_str = sha256_of_file(str(p))
    assert result_path == result_str


# ---------------------------------------------------------------------------
# normalized_config_hash
# ---------------------------------------------------------------------------


def test_config_hash_stable():
    """Same argv → same hash on repeated calls."""
    argv = ["rectify", "correct", "--sample", "S1", "--bam", "/data/foo.bam"]
    h1 = normalized_config_hash(argv)
    h2 = normalized_config_hash(argv)
    assert h1 == h2


def test_config_hash_starts_with_sha256():
    """Hash result is prefixed with 'sha256:'."""
    h = normalized_config_hash(["foo", "bar"])
    assert h.startswith("sha256:")


def test_config_hash_ignored_flag_unchanged():
    """Adding an ignored flag doesn't change the hash."""
    argv_base = ["rectify", "correct", "--sample", "S1"]
    argv_with = ["rectify", "correct", "--sample", "S1", "--n-threads", "8"]
    h_base = normalized_config_hash(argv_base, ignore={"--n-threads"})
    h_with = normalized_config_hash(argv_with, ignore={"--n-threads"})
    assert h_base == h_with


def test_config_hash_ignored_flag_equals_form():
    """--flag=value form is also filtered by normalized_config_hash."""
    argv_base = ["rectify", "correct", "--sample", "S1"]
    argv_with = ["rectify", "correct", "--sample", "S1", "--n-threads=16"]
    h_base = normalized_config_hash(argv_base, ignore={"--n-threads"})
    h_with = normalized_config_hash(argv_with, ignore={"--n-threads"})
    assert h_base == h_with


def test_config_hash_non_ignored_flag_changes():
    """Adding a non-ignored flag changes the hash."""
    argv_base = ["rectify", "correct", "--sample", "S1"]
    argv_extra = ["rectify", "correct", "--sample", "S1", "--extra-flag", "val"]
    h_base = normalized_config_hash(argv_base)
    h_extra = normalized_config_hash(argv_extra)
    assert h_base != h_extra


def test_config_hash_order_matters():
    """Reordering flags produces a different hash (order is NOT sorted out)."""
    argv_a = ["rectify", "--foo", "--bar"]
    argv_b = ["rectify", "--bar", "--foo"]
    h_a = normalized_config_hash(argv_a)
    h_b = normalized_config_hash(argv_b)
    assert h_a != h_b, (
        "Order should matter: two argvs with same tokens in different order "
        "must produce different hashes."
    )


def test_config_hash_empty_argv():
    """Empty argv is well-defined and produces a stable hash."""
    h = normalized_config_hash([])
    assert h.startswith("sha256:")
    assert h == normalized_config_hash([])
