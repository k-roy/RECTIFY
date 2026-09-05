"""A single-pass mapPacBio run over more than MPB_MONOLITH_MAX_READS reads is refused
with the chunk recipe (2026-08-22: 9 cohort arms lost to ALIGNER_TIMEOUT). Chunk mode
and the explicit override are not refused by the guard."""
import gzip, os
import pytest
from rectify.core.align import multi_aligner as ma


def _fastq(path, n):
    with gzip.open(path, 'wt') as fh:
        for i in range(n):
            fh.write(f'@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n')


@pytest.fixture
def big_fastq(tmp_path, monkeypatch):
    monkeypatch.setattr(ma, 'MPB_MONOLITH_MAX_READS', 50)
    p = tmp_path / 'reads.fastq.gz'; _fastq(p, 60); return p


def test_monolith_refused_with_chunk_recipe(big_fastq, tmp_path, monkeypatch):
    monkeypatch.delenv('RECTIFY_MAPPACBIO_ALLOW_MONOLITH', raising=False)
    with pytest.raises(RuntimeError, match=r'--mapPacBio-chunks \d+ --mapPacBio-chunk-idx K'):
        ma.run_map_pacbio(str(big_fastq), str(tmp_path / 'g.fa'), str(tmp_path / 'out.bam'))


def test_override_env_bypasses_guard(big_fastq, tmp_path, monkeypatch):
    monkeypatch.setenv('RECTIFY_MAPPACBIO_ALLOW_MONOLITH', '1')
    monkeypatch.setattr(ma.shutil, 'which', lambda _n: None)
    with pytest.raises(FileNotFoundError, match='mapPacBio.sh not found'):   # got PAST the guard
        ma.run_map_pacbio(str(big_fastq), str(tmp_path / 'g.fa'), str(tmp_path / 'out.bam'))


def test_chunk_mode_not_refused(big_fastq, tmp_path, monkeypatch):
    monkeypatch.delenv('RECTIFY_MAPPACBIO_ALLOW_MONOLITH', raising=False)
    monkeypatch.setattr(ma.shutil, 'which', lambda _n: None)
    with pytest.raises(FileNotFoundError, match='mapPacBio.sh not found'):
        ma.run_map_pacbio(str(big_fastq), str(tmp_path / 'g.fa'), str(tmp_path / 'out.bam'),
                          chunk_idx=0, n_chunks=4)
