import pytest

from rectify.core.align import multi_aligner as ma


def test_run_multi_aligner_raises_when_requested_minimap2_fails(tmp_path, monkeypatch):
    monkeypatch.setattr(ma, 'run_minimap2', lambda **kwargs: (_ for _ in ()).throw(RuntimeError('boom')))
    monkeypatch.setattr(ma, 'run_map_pacbio', lambda **kwargs: str(tmp_path / 'mpb.bam'))

    with pytest.raises(RuntimeError, match='Required baseline aligner minimap2 failed'):
        ma.run_multi_aligner(
            reads_path='reads.fastq',
            genome_path='genome.fa',
            output_dir=str(tmp_path),
            sample_name='sample',
            aligners=['minimap2', 'mapPacBio'],
        )


def test_run_multi_aligner_allows_non_minimap2_panel(tmp_path, monkeypatch):
    expected = str(tmp_path / 'mpb.bam')
    monkeypatch.setattr(ma, 'run_map_pacbio', lambda **kwargs: expected)

    results = ma.run_multi_aligner(
        reads_path='reads.fastq',
        genome_path='genome.fa',
        output_dir=str(tmp_path),
        sample_name='sample',
        aligners=['mapPacBio'],
    )

    assert results == {'mapPacBio': expected}


def test_run_multi_aligner_raises_when_all_requested_aligners_fail(tmp_path, monkeypatch):
    monkeypatch.setattr(ma, 'run_map_pacbio', lambda **kwargs: (_ for _ in ()).throw(RuntimeError('boom')))

    with pytest.raises(RuntimeError, match='All requested aligners failed'):
        ma.run_multi_aligner(
            reads_path='reads.fastq',
            genome_path='genome.fa',
            output_dir=str(tmp_path),
            sample_name='sample',
            aligners=['mapPacBio'],
        )
