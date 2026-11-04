#!/usr/bin/env python3
"""Append the S. cerevisiae narnavirus contigs (ScNV-20S, ScNV-23S) to the bundled yeast reference data.

House policy (Kevin Roy, 2026-11-03; Chanfreau planning/1248): every yeast reference RECTIFY aligns to carries the two
narnavirus RNA replicons as extra contigs. They are cytoplasmic, non-genomic and NOT polyadenylated, so oligo-dT libraries
miss them, but any library that tails RNA in vitro (ePAP / E. coli PAP) reads them: in Kevin He's 2026-10-29 ePAP PCR-cDNA,
ScNV-20S was 5-60 % of reads per library (0 % in the untailed libraries of the same RNA). Without a contig those reads are
unmapped and every mapped-fraction gate and per-million denominator is silently wrong.

Source of truth: rectify/data/genomes/saccharomyces_cerevisiae/narnaviruses_ScNV.fa (RefSeq NC_004051.1 -> ScNV-20S,
NC_004050.1 -> ScNV-23S, U->T, 60-col). This script is idempotent: it refuses to append to a file that already has them.

Rewrites, in rectify/data/genomes/saccharomyces_cerevisiae/:
  S288C_reference_sequence_R64-5-1_20240529.fsa.gz (+ .fai, .gzi; bgzip, samtools faidx)
  S288C_reference_sequence_R64-5-1_20240529.fsa.ncbi.gz
  saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz (features before SGD's ##FASTA section, sequences appended to it)
  saccharomyces_cerevisiae_R64-5-1_20240529.gtf

Usage:  python scripts/data/add_yeast_narnavirus_contigs.py   (needs bgzip + samtools on PATH)
"""
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

D = Path(__file__).resolve().parents[2] / 'rectify' / 'data' / 'genomes' / 'saccharomyces_cerevisiae'
SRC = D / 'narnaviruses_ScNV.fa'
# contig, length, CDS start, CDS end (RefSeq, 1-based inclusive), accession, dsRNA form
MODELS = [('ScNV-20S', 2514, 13, 2502, 'NC_004051.1', 'W dsRNA'),
          ('ScNV-23S', 2891, 7, 2832, 'NC_004050.1', 'T dsRNA')]
NAMES = [m[0] for m in MODELS]


def read_fasta(p):
    seqs, name = {}, None
    for line in p.read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = []
        elif line.strip():
            seqs[name].append(line.strip())
    return {k: ''.join(v) for k, v in seqs.items()}


def fasta_block(seqs, header_prefix=''):
    out = []
    for c, L, *_ in MODELS:
        s = seqs[c]
        assert len(s) == L, (c, len(s))
        out.append('>' + header_prefix + c)
        out += [s[i:i + 60] for i in range(0, L, 60)]
    return '\n'.join(out) + '\n'


def gff_features():
    rows = []
    for c, L, cs, ce, acc, ds in MODELS:
        note = ('Saccharomyces%20{}%20RNA%20narnavirus%20({}),%20non-genomic,%20non-polyadenylated;%20RefSeq%20{}'
                .format(c.split('-')[1], ds.replace(' ', '%20'), acc))
        rows += ['\t'.join([c, 'RefSeq', 'gene', '1', str(L), '.', '+', '.',
                            'ID={0};Name={0};gene={0};so_term_name=narnavirus_genome;Note={1}'.format(c, note)]),
                 '\t'.join([c, 'RefSeq', 'mRNA', '1', str(L), '.', '+', '.', 'ID={0}_mRNA;Name={0}_mRNA;Parent={0}'.format(c)]),
                 '\t'.join([c, 'RefSeq', 'CDS', str(cs), str(ce), '.', '+', '0',
                            'Parent={0}_mRNA;Name={0}_CDS;product=RNA-dependent%20RNA%20polymerase'.format(c)])]
    return rows


def gtf_rows():
    # Same feature types and attribute pattern as the SGD rows of the bundled GTF (gene / CDS / mRNA, NO exon lines): that file is the
    # exon-less case COMPASS's pre-flight detects and repairs (tests/test_compass_preflight_866.py); an exon line here would flip it.
    rows = []
    for c, L, cs, ce, *_ in MODELS:
        rows += ['\t'.join([c, 'RefSeq', 'gene', '1', str(L), '.', '+', '.',
                            'gene_id "{0}"; transcript_id "{0}"; Name "{0}";'.format(c)]),
                 '\t'.join([c, 'RefSeq', 'CDS', str(cs), str(ce), '.', '+', '0',
                            'gene_id "{0}_mRNA"; transcript_id "{0}_mRNA"; Name "{0}_CDS";'.format(c)]),
                 '\t'.join([c, 'RefSeq', 'mRNA', '1', str(L), '.', '+', '.',
                            'gene_id "{0}"; transcript_id "{0}"; Name "{0}_mRNA";'.format(c)])]
    return rows


def main():
    seqs = read_fasta(SRC)
    assert sorted(seqs) == sorted(NAMES), sorted(seqs)
    # 1. bgzipped FASTA (+ .fai/.gzi)
    fa_gz = D / 'S288C_reference_sequence_R64-5-1_20240529.fsa.gz'
    body = gzip.decompress(fa_gz.read_bytes()).decode()
    if '>ScNV-20S' in body:
        sys.exit('already augmented: ' + fa_gz.name)
    tmp = D / 'S288C_reference_sequence_R64-5-1_20240529.fsa'
    tmp.write_text(body.rstrip('\n') + '\n' + fasta_block(seqs))
    subprocess.run(['bgzip', '-f', str(tmp)], check=True)          # -> .fsa.gz (replaces the old one)
    for ext in ('.fai', '.gzi'):
        p = Path(str(fa_gz) + ext)
        if p.exists():
            p.unlink()
    subprocess.run(['samtools', 'faidx', str(fa_gz)], check=True)
    # 2. the NCBI-header variant (plain gzip)
    ncbi = D / 'S288C_reference_sequence_R64-5-1_20240529.fsa.ncbi.gz'
    b = gzip.decompress(ncbi.read_bytes()).decode()
    assert '>ScNV-20S' not in b
    ncbi.write_bytes(gzip.compress((b.rstrip('\n') + '\n' + fasta_block(seqs)).encode(), mtime=0))
    # 3. GFF: features before ##FASTA, sequences appended to the embedded FASTA section
    gff = D / 'saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz'
    lines = gzip.decompress(gff.read_bytes()).decode().split('\n')
    assert not any(l.startswith('ScNV-') for l in lines)
    k = lines.index('##FASTA')
    tail = lines[k:]
    while tail and tail[-1] == '':
        tail.pop()
    out = lines[:k] + gff_features() + tail + fasta_block(seqs).rstrip('\n').split('\n')
    gff.write_bytes(gzip.compress(('\n'.join(out) + '\n').encode(), mtime=0))
    # 4. GTF
    gtf = D / 'saccharomyces_cerevisiae_R64-5-1_20240529.gtf'
    g = gtf.read_text()
    assert 'ScNV-20S' not in g
    gtf.write_text(g.rstrip('\n') + '\n' + '\n'.join(gtf_rows()) + '\n')
    # verify
    fai = [l.split('\t')[0] for l in Path(str(fa_gz) + '.fai').read_text().splitlines()]
    assert fai[-2:] == NAMES and len(fai) == 19, fai
    print('bundled yeast reference augmented:', ', '.join(NAMES))


if __name__ == '__main__':
    main()
