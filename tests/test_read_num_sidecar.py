import importlib.util
import resource

import pytest

from rectify.core.chunking.sidecar import ReadNumSidecar, ReadNumSidecarWriter


def _write_sidecar(path, n):
    with ReadNumSidecarWriter(path, sample_id='sample') as w:
        for i in range(n):
            comment = '' if i == 0 else f'XA:Z:tail{i}\tXC:i:{i}'
            w.add(i, f'read{i}', comment, f'chunk_{i % 4:03d}_of_004', 'ACGT', 'IIII')


def test_round_trip_lookup_and_empty_comment(tmp_path):
    path = tmp_path / 'sample.read_num_sidecar.parquet'
    _write_sidecar(path, 100_000)
    sc = ReadNumSidecar.open(path)
    assert len(sc) == 100_000
    assert sc.lookup(99_999).original_qname == 'read99999'
    assert sc.lookup_by_qname('read17').read_num == 17
    assert sc.lookup(0).fastq_comment == ''


def test_lookup_many_large_batch(tmp_path):
    path = tmp_path / 'sample.read_num_sidecar.parquet'
    _write_sidecar(path, 100_000)
    sc = ReadNumSidecar.open(path)
    rows = sc.lookup_many(range(100_000))
    assert len(rows) == 100_000
    assert rows[42].original_qname == 'read42'


def test_verify_reports_fingerprint_mismatch(tmp_path, caplog):
    path = tmp_path / 'sample.read_num_sidecar.parquet'
    with ReadNumSidecarWriter(path, sample_id='sample') as w:
        w.add(1, 'read1', 'XA:Z:x', 'chunk_000_of_001', 'ACGT', 'IIII')
    sc = ReadNumSidecar.open(path)
    result = sc.verify(1, 'TGCA', 'IIII')
    assert not result.ok
    assert not result.seq_matches
    assert result.qual_matches
    assert 'fingerprint mismatch' in caplog.text


@pytest.mark.skipif(importlib.util.find_spec('pyarrow') is None, reason='pyarrow unavailable')
def test_streaming_write_peak_rss_under_200mb(tmp_path):
    path = tmp_path / 'sample.read_num_sidecar.parquet'
    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    with ReadNumSidecarWriter(path, sample_id='sample') as w:
        for i in range(1_000_000):
            w.add(i, f'read{i}', '', f'chunk_{i % 16:03d}_of_016', 'ACGT', 'IIII')
    after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux reports KiB. Normalize conservatively.
    delta_mb = (after - before) / (1024 * 1024 if after > 10_000_000 else 1024)
    assert delta_mb < 200
