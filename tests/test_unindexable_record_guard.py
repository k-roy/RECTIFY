"""A single off-contig record must not cost a multi-hour run its output (ISSUE-015).

One 2F rescue placed an 8,835-nt "exon" at chrM:-5843 — 1 record in 145,000 —
and `pysam.index` refused the whole file, so `rectify correct` exited 1 after
1.5 h with the corrected BAM written but unindexed and nothing said about which
read did it.  The admission and writer layers are fixed elsewhere; this is the
floor under them, and it has to hold for a defect nobody has thought of yet.

The contract:
  * the happy path costs nothing — no extra pass over the data;
  * on failure, every offending record is named at ERROR and written to
    `<out>.unindexable_reads.tsv`, the BAM is rewritten without them and
    indexed, and the run continues;
  * `sort_and_index_bam` never raises.

Author: Kevin R. Roy (agent S1)
"""

import sys
from pathlib import Path

import pysam
import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.commands import correct_command as cc  # noqa: E402

CHRM_LEN = 16_569
HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"},
          "SQ": [{"SN": "chrM", "LN": CHRM_LEN}]}


def _mk_read(header, name, pos, cigar, unmapped=False):
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.mapping_quality = 60
    n = sum(l for op, l in cigar if op in (0, 1, 4, 7, 8))
    r.query_sequence = "A" * n
    r.query_qualities = pysam.qualitystring_to_array("I" * n)
    if unmapped:
        r.is_unmapped = True
        r.reference_id = -1
        r.reference_start = -1
        return r
    r.reference_id = 0
    r.cigartuples = cigar
    r.reference_start = pos
    return r


def _write(path, reads_spec):
    header = pysam.AlignmentHeader.from_dict(HEADER)
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for spec in reads_spec:
            out.write(_mk_read(header, *spec))
    return path


# The real shape from ISSUE-015: an 8,134-nt 5' soft clip on the 16,569-nt
# mitochondrial genome, "rescued" into an 8,835-nt exon starting before base 0.
OFFENDER = ("b06afeb3", -5843, [(0, 8835)])
GOOD_A = ("good_a", 100, [(0, 50)])
GOOD_B = ("good_b", 4561, [(4, 8134), (0, 296)])


# ---------------------------------------------------------------------------
# unindexable_reason
# ---------------------------------------------------------------------------

@pytest.fixture
def header():
    return pysam.AlignmentHeader.from_dict(HEADER)


def test_negative_start_is_unindexable(header):
    read = _mk_read(header, *OFFENDER)
    reason = cc.unindexable_reason(read, {"chrM": CHRM_LEN})
    assert reason and 'reference_start=-5843' in reason


def test_start_past_the_contig_end_is_unindexable(header):
    read = _mk_read(header, "past", CHRM_LEN + 10, [(0, 50)])
    reason = cc.unindexable_reason(read, {"chrM": CHRM_LEN})
    assert reason and 'contig length' in reason


def test_end_past_the_contig_end_is_unindexable(header):
    read = _mk_read(header, "over", CHRM_LEN - 10, [(0, 50)])
    reason = cc.unindexable_reason(read, {"chrM": CHRM_LEN})
    assert reason and 'reference_end' in reason


def test_a_normal_record_is_fine(header):
    assert cc.unindexable_reason(_mk_read(header, *GOOD_A),
                                 {"chrM": CHRM_LEN}) is None


def test_an_unmapped_record_is_fine(header):
    """Unplaced reads index fine and must not be quarantined."""
    read = _mk_read(header, "unmapped", 0, [(0, 50)], unmapped=True)
    assert cc.unindexable_reason(read, {"chrM": CHRM_LEN}) is None


def test_a_contig_missing_from_the_header_is_not_judged(header):
    """No length to compare against — do not invent one."""
    assert cc.unindexable_reason(_mk_read(header, *GOOD_A), {}) is None


# ---------------------------------------------------------------------------
# sort_and_index_bam — the failure branch
# ---------------------------------------------------------------------------

@pytest.fixture
def offending_bam(tmp_path):
    return _write(tmp_path / "unsorted.bam", [GOOD_A, OFFENDER, GOOD_B])


def test_the_offending_bam_really_does_break_pysam_index(offending_bam, tmp_path):
    """Guard the fixture: if this ever stops failing, the tests below prove nothing."""
    plain = tmp_path / "plain.bam"
    pysam.sort('-o', str(plain), str(offending_bam))
    with pytest.raises(Exception):
        pysam.index(str(plain))


def test_run_continues_and_the_output_is_indexed(offending_bam, tmp_path):
    out = tmp_path / "corrected.bam"
    report = cc.sort_and_index_bam(offending_bam, out, log=None)

    assert report['indexed'] is True
    assert report['offenders'] == 1
    assert out.exists() and Path(str(out) + '.bai').exists()

    names = {r.query_name for r in pysam.AlignmentFile(str(out))}
    assert names == {"good_a", "good_b"}          # the offender is gone
    assert "b06afeb3" not in names


def test_the_offender_is_written_to_the_tsv_with_diagnostic_fields(
        offending_bam, tmp_path):
    out = tmp_path / "corrected.bam"
    report = cc.sort_and_index_bam(offending_bam, out, log=None)

    tsv = Path(report['offender_tsv'])
    assert tsv.name == "corrected" + cc.UNINDEXABLE_TSV_SUFFIX
    lines = tsv.read_text().strip().split('\n')
    assert lines[0].split('\t') == list(cc._UNINDEXABLE_TSV_COLUMNS)
    assert len(lines) == 2

    row = dict(zip(cc._UNINDEXABLE_TSV_COLUMNS, lines[1].split('\t')))
    assert row['read_name'] == "b06afeb3"
    assert row['contig'] == "chrM"
    assert row['contig_length'] == str(CHRM_LEN)
    assert row['reference_start'] == "-5843"
    assert 'reference_start=-5843' in row['reason']
    assert row['cigar'] == "8835M"


def test_every_offender_is_logged_at_error(offending_bam, tmp_path, caplog):
    import logging
    with caplog.at_level(logging.ERROR):
        cc.sort_and_index_bam(offending_bam, tmp_path / "corrected.bam", log=None)
    assert "b06afeb3" in caplog.text
    assert "chrM" in caplog.text
    assert "UNINDEXABLE RECORD" in caplog.text


# ---------------------------------------------------------------------------
# sort_and_index_bam — the happy path costs nothing
# ---------------------------------------------------------------------------

def test_a_clean_bam_is_sorted_and_indexed_with_no_extra_pass(tmp_path, monkeypatch):
    clean = _write(tmp_path / "clean.unsorted.bam", [GOOD_B, GOOD_A])

    def _must_not_run(*a, **k):
        raise AssertionError("the quarantine scan ran on a clean BAM")

    monkeypatch.setattr(cc, "quarantine_unindexable_records", _must_not_run)

    out = tmp_path / "clean.bam"
    report = cc.sort_and_index_bam(clean, out, log=None)
    assert report == {'indexed': True, 'offenders': 0, 'offender_tsv': None}
    assert Path(str(out) + '.bai').exists()
    assert not list(tmp_path.glob('*' + cc.UNINDEXABLE_TSV_SUFFIX))


# ---------------------------------------------------------------------------
# it never raises
# ---------------------------------------------------------------------------

def test_a_missing_input_does_not_raise(tmp_path):
    report = cc.sort_and_index_bam(tmp_path / "nope.bam", tmp_path / "out.bam",
                                   log=None)
    assert report['indexed'] is False
    assert report['offenders'] == 0


def test_an_index_failure_with_no_off_contig_record_does_not_raise(
        tmp_path, monkeypatch):
    """The other cause: report it, leave the BAM, keep going — and no stray TSV."""
    clean = _write(tmp_path / "clean.unsorted.bam", [GOOD_A])
    real_index = pysam.index

    def _boom(*a, **k):
        raise RuntimeError("simulated htslib failure")

    monkeypatch.setattr(pysam, "index", _boom)
    out = tmp_path / "out.bam"
    report = cc.sort_and_index_bam(clean, out, log=None)
    monkeypatch.setattr(pysam, "index", real_index)

    assert report['indexed'] is False
    assert report['offenders'] == 0
    assert report['offender_tsv'] is None
    assert not list(tmp_path.glob('*' + cc.UNINDEXABLE_TSV_SUFFIX))
