import subprocess

import pytest

from rectify.core.commands.align_command import _commit_indexed_bam


def test_commit_indexed_bam_indexes_temp_before_replacing_final(tmp_path):
    temp_bam = tmp_path / 'sample.multialigned.md.bam'
    final_bam = tmp_path / 'sample.multialigned.bam'
    temp_bam.write_text('new')
    final_bam.write_text('old')
    calls = []

    def fake_index(cmd, check):
        calls.append(cmd)
        assert cmd[-1] == str(temp_bam)
        (tmp_path / 'sample.multialigned.md.bam.bai').write_text('new-index')

    _commit_indexed_bam(temp_bam, final_bam, fake_index)

    assert calls == [['samtools', 'index', str(temp_bam)]]
    assert final_bam.read_text() == 'new'
    assert (tmp_path / 'sample.multialigned.bam.bai').read_text() == 'new-index'
    assert not temp_bam.exists()


def test_commit_indexed_bam_keeps_final_when_temp_index_fails(tmp_path):
    temp_bam = tmp_path / 'sample.multialigned.md.bam'
    final_bam = tmp_path / 'sample.multialigned.bam'
    temp_bam.write_text('new')
    final_bam.write_text('old')

    def failing_index(cmd, check):
        raise subprocess.CalledProcessError(1, cmd)

    with pytest.raises(subprocess.CalledProcessError):
        _commit_indexed_bam(temp_bam, final_bam, failing_index)

    assert final_bam.read_text() == 'old'
    assert temp_bam.read_text() == 'new'
