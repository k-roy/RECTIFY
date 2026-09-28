"""Full chimeric selection must hand a valid cDNA record to analyze."""
import random
from pathlib import Path

import pysam
import pytest

from rectify.core.consensus.consensus import run_consensus_selection
from rectify.core.commands.cdna_analyze_command import _read_info_from_bam_record
from rectify.utils.genome import register_genome_contigs


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('spliced', [False, True])
@pytest.mark.parametrize('encoded', [False, True])
def test_chimeric_winner_recovers_sibling_comment_block(tmp_path, reverse, spliced, encoded):
    rng = random.Random(510572)
    reference = ''.join(rng.choice('ACGT') for _ in range(1500))
    reference = reference[:245] + 'GT' + reference[247:323] + 'AG' + reference[325:]
    query = reference[200:245] + reference[325:380] if spliced else reference[200:300]
    genome = {'chrMeta': reference}
    register_genome_contigs(genome)
    header = pysam.AlignmentHeader.from_references(['chrMeta'], [len(reference)])
    fasta = tmp_path / 'genome.fa'
    fasta.write_text('>chrMeta\n' + reference + '\n')
    pysam.faidx(str(fasta))
    tags = {'XU': 'CGTA' * 6, 'XO': 'fwd', 'XN': 1, 'XT': 2,
            'XY': 'umi_not_captured', 'XC': 7, 'XF': 0, 'XQ': 11, 'XK': 5}
    paths = {}
    originals = {}
    for arm, start in [('minimap2', 910), ('uLTRA', 200)]:
        read = pysam.AlignedSegment(header)
        read.query_name = 'same_molecule'
        read.flag = 16 if reverse else 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigarstring = '45M80N55M' if arm == 'uLTRA' and spliced else '100M'
        read.query_sequence = query
        read.query_qualities = [35] * 100
        if arm == 'minimap2':
            for tag, value in tags.items():
                read.set_tag(tag, value)
        else:
            read.set_tag('XC', 'NO_SPLICE')
        bam = tmp_path / (arm + '.bam')
        with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
            out.write(read)
        if encoded:
            converted = tmp_path / (arm + '.calmd.bam')
            converted.write_bytes(pysam.calmd('-e', '-b', str(bam), str(fasta)))
            bam = converted
        paths[arm] = str(bam)
        originals[str(bam)] = bam.read_bytes()
    selected = tmp_path / 'selected.bam'
    run_consensus_selection(paths, genome, str(selected), n_workers=1, use_chimeric=True)
    with pysam.AlignmentFile(str(selected)) as handle:
        records = list(handle)
    assert len(records) == 1
    out = records[0]
    assert out.reference_start == 200
    assert out.cigarstring == ('45M80N55M' if spliced else '100M')
    assert out.query_sequence == query
    assert list(out.query_qualities) == [35] * 100
    for tag, value in tags.items():
        assert out.get_tag(tag) == value
    info, count = _read_info_from_bam_record(out, reference)
    assert count == 7
    assert info.orient == ('rev' if reverse else 'fwd')
    assert all(Path(path).read_bytes() == data
               for path, data in originals.items())
