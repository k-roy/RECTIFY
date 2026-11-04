"""Narnavirus house policy (Chanfreau planning/1248, 2026-11-03): every yeast reference carries ScNV-20S / ScNV-23S."""
import gzip
from pathlib import Path

import pytest

from rectify.config import CHROM_SIZES, YEAST_NARNAVIRUS_CONTIGS
from rectify.utils.genome import (looks_like_yeast_reference, missing_yeast_narnavirus_contigs,
                                  reference_contig_names, standardize_chrom_name)

D = Path(__file__).resolve().parents[1] / 'rectify' / 'data' / 'genomes' / 'saccharomyces_cerevisiae'
ROMAN = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']


@pytest.mark.parametrize('names', [
    ['chr' + r for r in ROMAN] + ['chrMito'],
    ROMAN + ['ref|NC_001224|'],
    ['ref|NC_0011%02d|' % (33 + i) for i in range(16)],
])
def test_yeast_reference_without_narnavirus_is_flagged(names):
    assert looks_like_yeast_reference(names)
    assert missing_yeast_narnavirus_contigs(names) == ['ScNV-20S', 'ScNV-23S']
    assert missing_yeast_narnavirus_contigs(names + ['ScNV-20S']) == ['ScNV-23S']
    assert missing_yeast_narnavirus_contigs(names + list(YEAST_NARNAVIRUS_CONTIGS)) == []


def test_non_yeast_reference_is_not_flagged():
    human = ['chr%d' % i for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
    assert not looks_like_yeast_reference(human)
    assert missing_yeast_narnavirus_contigs(human) == []
    # a partial yeast subset (e.g. a single-chromosome test reference) is not a yeast genome reference
    assert missing_yeast_narnavirus_contigs(['chrI', 'chrII']) == []


def test_bundled_yeast_genome_carries_narnavirus_contigs():
    names = reference_contig_names(D / 'S288C_reference_sequence_R64-5-1_20240529.fsa.gz')
    assert names[-2:] == ['ScNV-20S', 'ScNV-23S'] and len(names) == 19
    lengths = dict(l.split('\t')[:2] for l in
                   (D / 'S288C_reference_sequence_R64-5-1_20240529.fsa.gz.fai').read_text().splitlines())
    assert {k: int(lengths[k]) for k in YEAST_NARNAVIRUS_CONTIGS} == YEAST_NARNAVIRUS_CONTIGS
    with gzip.open(D / 'S288C_reference_sequence_R64-5-1_20240529.fsa.ncbi.gz', 'rt') as fh:
        ncbi = [l[1:].split()[0] for l in fh if l.startswith('>')]
    assert ncbi[-2:] == ['ScNV-20S', 'ScNV-23S']


def test_bundled_annotation_models_narnavirus_before_fasta_section():
    with gzip.open(D / 'saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz', 'rt') as fh:
        lines = fh.read().split('\n')
    k = lines.index('##FASTA')
    feats = [l.split('\t') for l in lines[:k] if l.startswith('ScNV-')]
    assert [(f[0], f[2]) for f in feats] == [('ScNV-20S', 'gene'), ('ScNV-20S', 'mRNA'), ('ScNV-20S', 'CDS'),
                                             ('ScNV-23S', 'gene'), ('ScNV-23S', 'mRNA'), ('ScNV-23S', 'CDS')]
    assert '>ScNV-20S' in lines[k:] and '>ScNV-23S' in lines[k:]
    gtf = (D / 'saccharomyces_cerevisiae_R64-5-1_20240529.gtf').read_text()
    assert gtf.count('\nScNV-20S\t') == 3 and gtf.count('\nScNV-23S\t') == 3


def test_narnavirus_contigs_are_valid_coordinates():
    for c, L in YEAST_NARNAVIRUS_CONTIGS.items():
        assert CHROM_SIZES[c] == L
        assert standardize_chrom_name(c) == c


def test_source_fasta_matches_bundled_sequence():
    import pysam
    src = {}
    name = None
    for l in (D / 'narnaviruses_ScNV.fa').read_text().splitlines():
        if l.startswith('>'):
            name = l[1:].split()[0]; src[name] = []
        else:
            src[name].append(l)
    with pysam.FastaFile(str(D / 'S288C_reference_sequence_R64-5-1_20240529.fsa.gz')) as fa:
        for c in YEAST_NARNAVIRUS_CONTIGS:
            assert fa.fetch(c) == ''.join(src[c])
