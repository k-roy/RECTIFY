"""edlib is a CORE dependency, and its absence is loud (2026-09-28).

GitHub CI installs ``.[dev]`` only, and with edlib then absent 13 cDNA tests failed: the fuzzy
SSP/adapter anchoring fell back to exact matching, so an error-bearing SSP made a Type 1 molecule
Type 2 and the walkback misread tail lengths. Every default ``pip install rectify-rna`` shipped
that behaviour without a word.
"""
import logging
import sys
from pathlib import Path

import pytest

from rectify.core.cdna.deps import edlib_available, warn_if_edlib_missing

_PYPROJECT = Path(__file__).resolve().parent.parent / "pyproject.toml"


def _core_dependencies():
    text = _PYPROJECT.read_text()
    block = text[text.index("\ndependencies = ["):]
    return block[:block.index("\n]")]


def test_edlib_is_a_core_dependency():
    if not _PYPROJECT.exists():
        pytest.skip("installed package without its pyproject.toml")
    assert '"edlib' in _core_dependencies()


def test_missing_edlib_warns_and_names_the_consequence(monkeypatch, caplog):
    monkeypatch.setitem(sys.modules, "edlib", None)  # makes `import edlib` raise ImportError
    assert edlib_available() is False
    log = logging.getLogger("test-cdna-deps")
    with caplog.at_level(logging.WARNING, logger="test-cdna-deps"):
        assert warn_if_edlib_missing(log, "correct-cdna") is False
    text = " ".join(r.getMessage() for r in caplog.records)
    assert "edlib is NOT installed" in text and "Type 2" in text and "correct-cdna" in text


def test_available_edlib_is_silent(caplog):
    if not edlib_available():   # not importorskip: pytest 8 re-raises a non-ModuleNotFoundError ImportError
        pytest.skip("edlib is not installed here")
    log = logging.getLogger("test-cdna-deps")
    with caplog.at_level(logging.WARNING, logger="test-cdna-deps"):
        assert warn_if_edlib_missing(log, "cdna-analyze") is True
    assert not caplog.records
