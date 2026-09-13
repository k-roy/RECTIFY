"""Junction evidence must survive encoding, multi-intron reads, and execution mode."""
from collections import Counter
import copy
import random

import pysam
import pytest

from rectify.core.splice import junction_scoring as J


CHROM = 'chrAudit'
INTRONS = [(100, 200), (300, 400)]


@pytest.fixture
def signal_fixture(tmp_path):
    rng = random.Random(20260913)
    ref = ''.join(rng.choice('ACGT') for _ in range(600))
    fasta = tmp_path / 'reference.fa'
    fasta.write_text(f'>{CHROM}\n{ref}\n')
    pysam.faidx(str(fasta))
    header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': CHROM, 'LN': len(ref)}]})
    read = pysam.AlignedSegment(header)
    read.query_name = 'two_junctions'
    read.reference_id = 0
    read.reference_start = 0
    read.mapping_quality = 60
    read.cigarstring = '100M100N100M100N100M'
    query = list(ref[:100] + ref[200:300] + ref[400:500])
    # Two mismatches at J1, none at J2, one in the body away from BOTH junctions.
    for i in (98, 101, 150):
        query[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[query[i]]
    read.query_sequence = ''.join(query)
    return fasta, ref, read


def _bam(path, read):
    with pysam.AlignmentFile(str(path), 'wb', header=read.header) as out:
        out.write(read)
    return str(path)


def test_mismatch_counts_are_identical_for_literal_and_reference_compressed_seq(signal_fixture):
    fasta, _ref, read = signal_fixture
    encoded = copy.copy(read)
    seq = ['='] * len(read.query_sequence)
    for i in (98, 101, 150):
        seq[i] = read.query_sequence[i]
    encoded.query_sequence = ''.join(seq)
    edges = [p for intron in INTRONS for p in intron]
    with pysam.FastaFile(str(fasta)) as fa:
        literal = J._mismatch_bins(read, fa, edges, J.JUNCTION_MM_WINDOW, J.JUNCTION_MM_VARIANT_REACH)
        compressed = J._mismatch_bins(encoded, fa, edges, J.JUNCTION_MM_WINDOW, J.JUNCTION_MM_VARIANT_REACH)
    assert literal[0] == 2 and literal[2] == 1  # non-vacuous control
    assert compressed == literal


def test_one_bad_junction_does_not_contaminate_its_clean_sibling(signal_fixture, tmp_path):
    fasta, _ref, read = signal_fixture
    bam = _bam(tmp_path / 'arm.minimap2.bam', read)
    *_, mismatch, _positions = J._collect_junction_counts_core(
        bam, unspliced_out=Counter(), fasta_path=str(fasta))
    first, second = [(CHROM, *j) for j in INTRONS]
    assert mismatch[first][0] == 2
    assert mismatch[second][0] == 0
    # Body statistics are shared across the read, but local counts are not.
    assert mismatch[first][2:5] == mismatch[second][2:5] == [1, 250, 1]
    assert mismatch[first][1] == mismatch[second][1] == 25


def test_single_bam_unspliced_signal_is_not_counted_twice(signal_fixture, tmp_path):
    _fasta, ref, read = signal_fixture
    read.reference_start = 50
    read.cigarstring = '100M'
    read.query_sequence = ref[50:150]  # crosses only the first intron's donor
    bam = _bam(tmp_path / 'linear.minimap2.bam', read)
    junction = (CHROM, *INTRONS[0])
    _, _, signal = J.build_junction_pool([bam], {junction}, return_signal=True)
    assert signal['unspliced'][junction] == 1


def test_pool_fallback_preserves_all_signal_arguments(signal_fixture, tmp_path, monkeypatch):
    import concurrent.futures

    fasta, _ref, read = signal_fixture
    paths = [_bam(tmp_path / f'arm{i}.minimap2.bam', read) for i in range(2)]

    def unavailable(*args, **kwargs):
        raise OSError('process pool unavailable for this test')

    monkeypatch.setattr(concurrent.futures, 'ProcessPoolExecutor', unavailable)
    _, _, signal = J.build_junction_pool(paths, set(), return_signal=True, fasta_path=str(fasta))
    first, second = [(CHROM, *j) for j in INTRONS]
    assert signal['site_support'] == {first: 1, second: 1}  # max over arms, never sum
    assert signal['junction_mismatch'][first][4] == 2
    assert signal['junction_mismatch'][second][4] == 2
