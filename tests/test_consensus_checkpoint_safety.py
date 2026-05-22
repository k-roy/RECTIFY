import os

import pysam
import pytest

from rectify.core.consensus import consensus as c


def _write_valid_bam(path):
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'SN': 'chrI', 'LN': 1000}],
    }
    with pysam.AlignmentFile(path, 'wb', header=header) as out:
        read = pysam.AlignedSegment()
        read.query_name = 'read1'
        read.query_sequence = 'ACGT'
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 10
        read.mapping_quality = 60
        read.cigartuples = [(0, 4)]
        read.query_qualities = pysam.qualitystring_to_array('IIII')
        out.write(read)


def test_checkpoint_batch_done_marker_requires_valid_bam(tmp_path):
    bam_path = tmp_path / 'consensus_batch_000000.bam'
    done_path = tmp_path / 'consensus_batch_000000.done'
    bam_path.write_text('not a bam')

    with pytest.raises(RuntimeError):
        c._commit_checkpoint_batch(str(bam_path), str(done_path))

    assert not done_path.exists()


def test_checkpoint_batch_commit_writes_done_after_valid_bam(tmp_path):
    bam_path = tmp_path / 'consensus_batch_000000.bam'
    done_path = tmp_path / 'consensus_batch_000000.done'
    _write_valid_bam(str(bam_path))

    c._commit_checkpoint_batch(str(bam_path), str(done_path))

    assert done_path.read_text() == 'ok\n'


def test_copy2_and_fsync_flushes_destination_and_parent(tmp_path, monkeypatch):
    src = tmp_path / 'src.bam'
    dst = tmp_path / 'dst.bam'
    src.write_bytes(b'bam bytes')
    fsync_calls = []

    def record_fsync(fd):
        fsync_calls.append(fd)

    monkeypatch.setattr(c.os, 'fsync', record_fsync)

    c._copy2_and_fsync(str(src), str(dst))

    assert dst.read_bytes() == b'bam bytes'
    assert len(fsync_calls) >= 2
    assert os.path.exists(dst)
