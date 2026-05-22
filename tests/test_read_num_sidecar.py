import importlib.util
import resource

import pytest

from rectify.core.chunking.sidecar import (
    ReadNumSidecar,
    ReadNumSidecarWriter,
    reconstruct_posthoc_sidecar_from_chunks,
)


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


def test_reconstruct_posthoc_sidecar_from_round_robin_chunks(tmp_path):
    chunks = []
    records_by_chunk = {
        0: [
            ('@read0 XA:Z:a\n', 'ACGT\n', '+\n', 'IIII\n'),
            ('@read3\n', 'TTTT\n', '+\n', 'JJJJ\n'),
        ],
        1: [
            ('@read1\n', 'CCCC\n', '+\n', 'KKKK\n'),
            ('@read4\n', 'GGGG\n', '+\n', 'LLLL\n'),
        ],
        2: [
            ('@read2\n', 'AAAA\n', '+\n', 'MMMM\n'),
        ],
    }
    for idx, records in records_by_chunk.items():
        path = tmp_path / f'sample_chunk_{idx:03d}_of_003.fastq'
        with path.open('w') as fh:
            for rec in records:
                fh.writelines(rec)
        chunks.append(path)

    sidecar_path = tmp_path / 'sample.read_num_sidecar.parquet'
    reconstruct_posthoc_sidecar_from_chunks(
        chunks,
        sidecar_path,
        sample_id='sample',
        n_chunks=3,
    )

    sc = ReadNumSidecar.open(sidecar_path)
    assert len(sc) == 5
    assert sc.lookup(0).original_qname == 'read0'
    assert sc.lookup(0).fastq_comment == 'XA:Z:a'
    assert sc.lookup(1).original_qname == 'read1'
    assert sc.lookup(2).original_qname == 'read2'
    assert sc.lookup(3).original_qname == 'read3'
    assert sc.lookup(4).original_qname == 'read4'
    assert (tmp_path / 'sample.POSTHOC_PROVENANCE.json').exists()


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
