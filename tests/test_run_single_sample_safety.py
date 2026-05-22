import subprocess

import pytest

from rectify.core.commands.run import single_sample


def test_run_samtools_fastq_preserves_stderr_on_failure(tmp_path, monkeypatch):
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append((cmd, kwargs))
        return subprocess.CompletedProcess(cmd, 2, stdout='', stderr='bad bam')

    monkeypatch.setattr(subprocess, 'run', fake_run)

    with pytest.raises(RuntimeError, match='bad bam'):
        single_sample._run_samtools_fastq(
            tmp_path / 'trimmed.bam',
            tmp_path / 'trimmed.fastq.gz',
            threads=4,
        )

    assert calls[0][1]['stderr'] is subprocess.PIPE
    assert calls[0][1]['stdout'] is subprocess.PIPE
    assert calls[0][1]['text'] is True
