"""Junction-pool provenance and the silent-miss diagnostics (ISSUE-011).

Three ways a junction pool goes wrong without saying anything:

1. a `rectify prescan` pickle written by a build with a coordinate defect is
   reused verbatim (the fix outlives nothing);
2. `--aligner-bams` is the input BAM itself, so the "multi-aligner pool" is the
   input's own error set;
3. the pool is keyed by chromosome names the reads never use, so every candidate
   lookup returns nothing and Module 2H reports "0 refined" as if the alignments
   were already perfect.

Author: Kevin R. Roy (agent S1)
"""

import logging
import sys
from pathlib import Path

import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.splice import junction_scoring as js  # noqa: E402
from rectify.core.commands.correct_command import _warn_self_pool  # noqa: E402


# ---------------------------------------------------------------------------
# 1. cache provenance
# ---------------------------------------------------------------------------

def test_a_freshly_stamped_cache_is_accepted():
    pool_data = {'all_junctions': set(), 'annotated_set': set(),
                 **js.junction_pool_cache_stamp()}
    assert js.junction_pool_cache_problem(pool_data) is None


def test_an_unstamped_cache_is_refused_and_says_why():
    problem = js.junction_pool_cache_problem(
        {'all_junctions': set(), 'annotated_set': set()})
    assert problem and 'cache_format' in problem


def test_a_cache_from_an_older_format_is_refused():
    problem = js.junction_pool_cache_problem(
        {'all_junctions': set(), 'annotated_set': set(), 'cache_format': 1})
    assert problem and '1' in problem


def test_a_cache_from_a_newer_format_is_refused_too():
    problem = js.junction_pool_cache_problem(
        {'cache_format': js.JUNCTION_POOL_CACHE_FORMAT + 1})
    assert problem is not None


def test_non_dict_payload_is_refused_without_raising():
    assert js.junction_pool_cache_problem(['not', 'a', 'dict']) is not None


def test_prescan_writes_the_stamp_this_reader_expects():
    """Writer and reader must not drift: prescan builds pool_data with the stamp."""
    src = (RECTIFY_ROOT / 'rectify/core/commands/prescan_command.py').read_text()
    assert 'junction_pool_cache_stamp()' in src


# ---------------------------------------------------------------------------
# 2. self-pool warning
# ---------------------------------------------------------------------------

def test_self_pool_is_warned(tmp_path, caplog):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    with caplog.at_level(logging.WARNING):
        fired = _warn_self_pool(
            {'bam_path': str(bam), 'aligner_bams': [str(bam)]},
            logging.getLogger("t"),
        )
    assert fired
    assert 'input BAM itself' in caplog.text


def test_distinct_aligner_bams_do_not_warn(tmp_path, caplog):
    bam = tmp_path / "sample.bam"
    other = tmp_path / "minimap2.bam"
    bam.write_bytes(b"")
    other.write_bytes(b"")
    with caplog.at_level(logging.WARNING):
        fired = _warn_self_pool(
            {'bam_path': str(bam), 'aligner_bams': [str(other)]},
            logging.getLogger("t"),
        )
    assert not fired
    assert caplog.text == ""


# ---------------------------------------------------------------------------
# 3. missing-chromosome diagnostic
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _reset_warned():
    js._MISSING_CHROM_WARNED.clear()
    yield
    js._MISSING_CHROM_WARNED.clear()


def test_missing_chrom_key_is_reported_with_the_pool_keys(caplog):
    idx = js._build_junction_index({('chrV', 100, 200)})
    with caplog.at_level(logging.WARNING):
        assert js._candidates_near(idx, 'chr5', 100, 200, radius=50) == []
    assert 'chr5' in caplog.text
    assert 'chrV' in caplog.text


def test_missing_chrom_warns_once_per_chromosome(caplog):
    idx = js._build_junction_index({('chrV', 100, 200)})
    with caplog.at_level(logging.WARNING):
        for _ in range(5):
            js._candidates_near(idx, 'chr5', 100, 200, radius=50)
    assert caplog.text.count('no entry for chromosome') == 1


def test_an_empty_index_does_not_warn(caplog):
    """No pool at all is a different (already logged) situation."""
    with caplog.at_level(logging.WARNING):
        assert js._candidates_near({}, 'chr5', 100, 200, radius=50) == []
    assert caplog.text == ""


def test_a_present_chromosome_with_no_nearby_junction_does_not_warn(caplog):
    idx = js._build_junction_index({('chr5', 10_000, 10_500)})
    with caplog.at_level(logging.WARNING):
        assert js._candidates_near(idx, 'chr5', 100, 200, radius=50) == []
    assert caplog.text == ""
