"""ISSUE-022: `correct` must write the same rows at 1 thread as at several, streaming or not.

Root cause (2026-10-13): the Module 2F exclusion detector (rDNA / Pol III regions, cf5ebb9) reached only the
non-streaming multi-thread path. The 1-thread fallback, both streaming paths (production's `--streaming`) and the
single-thread streaming path ran the 3'SS rescue there anyway, while the log said it was skipped. On the upf1d bundle
that turned read ca7a0b66's 5' region, including its aligner-called 128N, into a 62-nt soft clip at 1 thread only
(and moved a 5' end of uLTRA's fa816f03). Every path now receives the detector.
"""

from pathlib import Path

import pytest

from rectify.core.bam import parallel as bp

DATA = Path(__file__).resolve().parent.parent / 'rectify' / 'data'
BAM = DATA / 'validation' / 'aligners' / 'validation_reads_upf1d.mapPacBio.bam'
CA7 = 'ca7a0b66-63bc-4064-a08a-6352edc2ddc3'


@pytest.fixture(scope='module')
def inputs():
    from rectify.data import get_bundled_annotation_path, get_bundled_genome_path
    genome, annot = get_bundled_genome_path('saccharomyces_cerevisiae'), get_bundled_annotation_path('saccharomyces_cerevisiae')
    if not BAM.exists() or genome is None or annot is None:
        pytest.skip('upf1d bundle or bundled S. cerevisiae reference not available')
    from rectify.core.consensus.consensus import load_annotated_junctions
    from rectify.core.exclusion_regions import ExclusionRegionDetector
    det = ExclusionRegionDetector(flanking_bp=100)
    assert det.load_from_gff(Path(str(annot)), exclude_tRNA=True, exclude_snRNA=True, exclude_rDNA=True,
                             exclude_mito=False) > 0
    return dict(genome=str(genome), annotated=load_annotated_junctions(str(annot)), det=det)


def _key(rows):
    out = {}
    for r in rows:
        if r.get('is_primary_result') in (False, 'False', '0'):
            continue
        j = r.get('junctions')
        if isinstance(j, (list, tuple)):
            j = ';'.join(f'{a}-{b}' for a, b in j)
        rescued = str(r.get('five_prime_rescued'))
        rescued = '1' if rescued in ('True', '1') else '0'      # a dict holds bools; a TSV holds 0/1
        out.setdefault(r['read_id'], (str(r['corrected_3prime']), j or '', rescued, str(r.get('five_prime_position'))))
    return out


def _tsv_rows(path):
    import csv
    files = [Path(path)] if Path(path).exists() else sorted(
        f for f in Path(path).parent.glob(Path(path).stem + '.region_*.tsv') if '_stats' not in f.name)
    rows = []
    for f in files:
        with open(f) as fh:
            rows.extend(csv.DictReader(fh, delimiter='\t'))
    return rows


def _in_memory(inputs, n):
    return bp.process_bam_file_parallel(
        bam_path=str(BAM), genome_path=inputs['genome'], n_threads=n, output_path=None, show_progress=False,
        return_stats=False, annotated_junctions=inputs['annotated'], apply_3ss_rescue=True,
        exclusion_detector=inputs['det'])


@pytest.mark.slow
def test_one_thread_matches_two_threads_in_memory(inputs):
    one, two = _key(_in_memory(inputs, 1)), _key(_in_memory(inputs, 2))
    assert one == two
    assert '476768-476896' in one[CA7][1]          # the aligner's 128N survives at 1 thread too


@pytest.mark.slow
def test_streaming_paths_match_the_in_memory_path(inputs, tmp_path):
    ref = _key(_in_memory(inputs, 2))
    single = tmp_path / 'single.tsv'
    bp.process_bam_streaming(bam_path=str(BAM), genome_path=inputs['genome'], output_path=str(single),
                             show_progress=False, annotated_junctions=inputs['annotated'], apply_3ss_rescue=True,
                             exclusion_detector=inputs['det'])
    par = tmp_path / 'par.tsv'
    bp.process_bam_streaming_parallel(bam_path=str(BAM), genome_path=inputs['genome'], output_path=str(par),
                                      n_threads=2, show_progress=False, annotated_junctions=inputs['annotated'],
                                      apply_3ss_rescue=True, exclusion_detector=inputs['det'])
    assert _key(_tsv_rows(single)) == ref
    assert _key(_tsv_rows(par)) == ref


def _tiny_inputs(tmp_path):
    import pysam
    fa = tmp_path / 'g.fa'
    fa.write_text('>chrT\n' + 'ACGTTGCA' * 125 + '\n')
    pysam.faidx(str(fa))
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chrT', 'LN': 1000}]}
    bam = tmp_path / 'in.bam'
    with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
        for i, pos in enumerate((100, 300)):
            r = pysam.AlignedSegment(out.header)
            r.query_name, r.reference_id, r.reference_start, r.flag = f'q{i}', 0, pos, 0
            r.cigarstring = '50M'
            r.query_sequence = ('ACGTTGCA' * 125)[pos:pos + 50]
            r.query_qualities = pysam.qualitystring_to_array('I' * 50)
            r.mapping_quality = 60
            out.write(r)
    pysam.index(str(bam))
    return str(fa), str(bam)


def test_every_in_process_path_hands_the_detector_to_the_per_read_correction(tmp_path, monkeypatch):
    """Fast guard for ISSUE-022: the sequential fallback and the single-thread streaming path pass the detector."""
    fa, bam = _tiny_inputs(tmp_path)
    seen = []
    real = bp.correct_read_3prime

    def spy(read, genome, **kw):
        seen.append(kw.get('exclusion_detector'))
        return real(read, genome, **{k: v for k, v in kw.items() if k != 'exclusion_detector'})

    monkeypatch.setattr(bp, 'correct_read_3prime', spy)
    det = object()
    bp.process_bam_file_parallel(bam_path=bam, genome_path=fa, n_threads=1, output_path=None, show_progress=False,
                                 return_stats=False, apply_3ss_rescue=False, exclusion_detector=det)
    bp.process_bam_streaming(bam_path=bam, genome_path=fa, output_path=str(tmp_path / 's.tsv'), show_progress=False,
                             apply_3ss_rescue=False, exclusion_detector=det)
    assert seen and all(d is det for d in seen), seen
