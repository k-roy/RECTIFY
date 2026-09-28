"""Actual endpoint loading/index/analysis contracts, without score mocks."""
import argparse
import csv
from contextlib import redirect_stdout, redirect_stderr

import pandas as pd
import pysam
import pytest

from rectify.core.analyze.loaders import load_corrected_positions, _load_large_file_chunked
from rectify.core.analyze.manifest import _tss_position_counts
from rectify.core.bam.bam_processor import correct_read_3prime
from rectify.core.bam.output import write_output_tsv
from rectify.core.commands.analyze_command import create_analyze_parser, run_analyze
from rectify.core.position_index import write_position_index
from rectify.utils.genome import register_genome_contigs


def _write(path, rows):
    with path.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def _rows(chrom, include_tss=True):
    # Alternative CPA rows keep the molecule's5prime position. Another
    # molecule shares its CPA but has a different TSS, exercising joint keys.
    rows = []
    for strand, offset in [('+', 0), ('-', 2000)]:
        for name, five, parts in [('a', 410+offset, [(830+offset, .25), (831+offset, .75)]),
                                  ('b', 514+offset, [(830+offset, 1.)]),
                                  ('c', None, [(830+offset, 1.)])]:
            for end, weight in parts:
                row = dict(read_id=strand+name, sample='WT_1', chrom=chrom, strand=strand,
                           corrected_3prime=end, fraction=weight)
                if include_tss:
                    row['five_prime_position'] = five
                rows.append(row)
    return rows


def _run(tmp_path, input_path, name, manifest=False):
    out = tmp_path/name
    parser = argparse.ArgumentParser()
    create_analyze_parser(parser.add_subparsers(dest='command'))
    argv = ['analyze']
    if manifest:
        path = tmp_path/(name+'.samples.tsv')
        _write(path, [dict(sample_id='WT_1', path=str(input_path), condition='WT')])
        argv += ['--manifest', str(path)]
    else:
        argv += [str(input_path)]
    argv += ['-o', str(out), '--min-reads', '1', '--min-cluster-samples', '1', '--threads', '1',
             '--include-mito', '--no-genomic-distribution', '--gene-attribution-mode', 'none',
             '--sample-sets', '{"test":["WT"]}']
    with (tmp_path/(name+'.log')).open('w') as log, redirect_stdout(log), redirect_stderr(log):
        assert run_analyze(parser.parse_args(argv)) == 0
    return out


def _mass(out, name):
    return pd.read_csv(out/name, sep='\t').drop(columns='cluster_id').sum(numeric_only=True).sum()


@pytest.mark.parametrize('chrom', ['chrV', 'chr5'])
def test_direct_regional_and_chunked_keep_joint_endpoint_mass(tmp_path, chrom):
    rows = _rows(chrom)
    path = tmp_path/'reads.tsv'
    _write(path, rows)
    regional = tmp_path/'regions.tsv'
    _write(regional, [dict(region_id='r0', chrom=chrom, start=0, end=5000,
                          tsv_path='reads.tsv', n_rows=len(rows), sha256='unused_by_loader')])
    small = load_corrected_positions(str(path), 'sample')
    reg = load_corrected_positions(str(regional), 'sample')
    pd.testing.assert_frame_equal(small, reg)
    assert 'five_prime_position' in small
    large = _load_large_file_chunked(str(path), 'sample', True, 'passthrough', 2, None)
    expected = (small.groupby(['chrom', 'strand', 'corrected_position', 'sample', 'five_prime_position'], dropna=False)
                ['fraction'].sum().rename('count').reset_index())
    pd.testing.assert_frame_equal(large, expected)
    assert large['count'].sum() == 6
    assert large.dropna(subset=['five_prime_position'])['count'].sum() == 4
    assert large['five_prime_position'].isna().sum() == 2
    # Explicit pre-aggregated inputs use their count once, including fraction.
    counted = tmp_path/'counted.tsv'
    _write(counted, [dict(r, count=3) for r in rows])
    counted_large = _load_large_file_chunked(str(counted), 'sample', False, 'passthrough', 3, None)
    assert counted_large['count'].sum() == 18
    assert counted_large.dropna(subset=['five_prime_position'])['count'].sum() == 12


@pytest.mark.parametrize('chrom', ['chrV', 'chr5'])
def test_all_entry_paths_keep_weighted_tss_and_cpa_separate(tmp_path, chrom):
    path = tmp_path/'reads.tsv'
    _write(path, _rows(chrom))
    direct = _run(tmp_path, path, 'direct')
    stream = _run(tmp_path, path, 'stream', manifest=True)
    write_position_index(path, str(path))
    indexed = _run(tmp_path, path, 'indexed', manifest=True)
    regional_path=tmp_path/'regions.tsv'
    _write(regional_path,[dict(region_id='r0',chrom=chrom,start=0,end=5000,
                              tsv_path='reads.tsv',n_rows=8,sha256='unused_by_loader')])
    regional=_run(tmp_path,regional_path,'regional')
    chunked_frame=_load_large_file_chunked(str(path),'sample',True,'passthrough',2,None)
    chunked_path=tmp_path/'joint_counts.tsv'
    chunked_frame.to_csv(chunked_path,sep='\t',index=False)
    chunked=_run(tmp_path,chunked_path,'chunked_handoff')
    for out in [direct, stream, indexed, regional, chunked]:
        assert _mass(out, 'cluster_counts.tsv') == 6
        assert _mass(out, 'tss_cluster_counts.tsv') == 4
        assert pd.read_csv(out/'tss_clusters.tsv', sep='\t')['n_reads'].sum() == 4
    for name in ['tss_clusters.tsv', 'tss_cluster_counts.tsv', 'cluster_counts.tsv']:
        pd.testing.assert_frame_equal(pd.read_csv(stream/name, sep='\t'), pd.read_csv(indexed/name, sep='\t'))


@pytest.mark.parametrize('with_tss', [False, True])
def test_legacy_no_fraction_and_missing_tss(tmp_path, with_tss):
    rows = [{k:v for k,v in r.items() if k!='fraction'} for r in _rows('chrV', with_tss)]
    path = tmp_path/'legacy.tsv'
    _write(path, rows)
    loaded = _load_large_file_chunked(str(path), 'sample', False, 'passthrough', 2, None)
    assert loaded['count'].sum() == len(rows)
    assert ('five_prime_position' in loaded) == with_tss
    out = _run(tmp_path, path, 'legacy', manifest=True)
    assert _mass(out, 'cluster_counts.tsv') == len(rows)
    if with_tss:
        assert _mass(out, 'tss_cluster_counts.tsv') == 6
    else:
        assert not (out/'tss_cluster_counts.tsv').exists()


@pytest.mark.parametrize('chrom', ['chrV', 'chr5'])
def test_actual_chimeric_producer_empty_qc_and_ag_rich_control(tmp_path, chrom):
    genome = {chrom: 'CGTC' * 1000}
    register_genome_contigs(genome)
    header = pysam.AlignmentHeader.from_references([chrom], [4000])
    results = []
    for flag, start in [(0, 100), (16, 1900)]:
        r = pysam.AlignedSegment(header)
        r.query_name = f'producer_{flag}'
        r.reference_id = 0
        r.reference_start = start
        r.flag = flag
        r.cigarstring = '3H2S20M37N24M3S5H'
        r.query_sequence = 'CC'+genome[chrom][start:start+20]+genome[chrom][start+57:start+81]+'TTC'
        r.query_qualities = [31] * 49
        r.set_tag('Xz', 1)
        original = r.to_string()
        result = correct_read_3prime(r, genome, apply_atract=False, apply_3ss_rescue=False)
        assert len(result)==1 and result[0]['qc_flags']==''
        assert result[0]['chrom']==chrom
        assert r.to_string()==original
        results += result
    path = tmp_path/'producer.tsv'
    write_output_tsv(results, str(path))
    assert pd.read_csv(path, sep='\t')['qc_flags'].isna().all()
    out = _run(tmp_path, path, 'blank', manifest=True)
    assert _mass(out, 'cluster_counts.tsv') == 2
    assert _mass(out, 'tss_cluster_counts.tsv') == 2
    assert pd.read_csv(out/'cpa_clusters.tsv', sep='\t')['n_reads_ag_rich'].sum()==0
    # Consumer QC policy control is explicitly supplied metadata, not a
    # claim that this synthetic chimeric alignment triggered AG detection.
    rows = _rows(chrom)
    for i,row in enumerate(rows):
        row['qc_flags'] = 'AG_RICH' if i==0 else (None if i%2 else 'PASS')
    mixed = tmp_path/'mixed.tsv';_write(mixed, rows)
    out = _run(tmp_path, mixed, 'mixed', manifest=True)
    assert _mass(out, 'cluster_counts.tsv') == 6
    assert pd.read_csv(out/'cpa_clusters.tsv', sep='\t')['n_reads_ag_rich'].sum()==.25


def test_missing_source_still_fails_with_or_without_index(tmp_path):
    path = tmp_path/'absent.tsv'
    with pytest.raises(FileNotFoundError):
        load_corrected_positions(str(path), 'sample')
    _write(path, _rows('chrV'))
    write_position_index(path, str(path))
    path.unlink()
    with pytest.raises(FileNotFoundError):
        _run(tmp_path, path, 'orphan', manifest=True)


def test_tss_fractional_mass_and_chromosome_exclusion(tmp_path):
    path = tmp_path/'fractional.tsv'
    rows = [dict(chrom='chrV', strand='+', five_prime_position=41, fraction=.3),
            dict(chrom='chrV', strand='+', five_prime_position=41, fraction=.7),
            dict(chrom='chrV', strand='+', five_prime_position=41, fraction=.5),
            dict(chrom='chr5', strand='-', five_prime_position=97, fraction=1.)]
    _write(path, rows)
    assert dict(_tss_position_counts(path, 'passthrough')) == {('chrV','+',41):1.5,('chr5','-',97):1.}
    assert dict(_tss_position_counts(path, 'passthrough', {'chr5'})) == {('chrV','+',41):1.5}


def test_zero_fraction_does_not_create_tss_support(tmp_path):
    rows = _rows('chrV')
    for row in rows:
        if row['read_id']=='+b':
            row['fraction']=0.
    path=tmp_path/'zero.tsv';_write(path,rows)
    stream=_run(tmp_path,path,'zero',manifest=True)
    assert _mass(stream,'cluster_counts.tsv')==5
    assert _mass(stream,'tss_cluster_counts.tsv')==3
    clusters=pd.read_csv(stream/'tss_clusters.tsv',sep='\t')
    assert not ((clusters['strand']=='+') & (clusters['start']==514)).any()


def test_empty_qc_chunk_before_ag_rich_chunk(tmp_path):
    path=tmp_path/'chunked_qc.tsv'
    with path.open('w') as f:
        f.write('read_id\tchrom\tstrand\tcorrected_3prime\tfive_prime_position\tfraction\tqc_flags\n')
        for i in range(100000):
            f.write(f'empty_{i}\tchrV\t+\t831\t410\t0.25\t\n')
        f.write('ag_control\tchrV\t+\t831\t410\t1\tAG_RICH\n')
    out=_run(tmp_path,path,'chunked_qc',manifest=True)
    assert _mass(out,'cluster_counts.tsv')==25001
    assert _mass(out,'tss_cluster_counts.tsv')==25001
    assert pd.read_csv(out/'cpa_clusters.tsv',sep='\t')['n_reads_ag_rich'].sum()==1


def test_fractional_cluster_reports_match_discovery_and_sample_mass(tmp_path):
    path=tmp_path/'partial_assignments.tsv'
    # A bounded subset of proportional assignments has1.5 total mass. A
    # reporting field must not turn that into one whole read.
    rows=[dict(read_id='a',sample='WT_1',chrom='chrV',strand='+',corrected_3prime=831,five_prime_position=410,fraction=1.),
          dict(read_id='b',sample='WT_1',chrom='chrV',strand='+',corrected_3prime=831,five_prime_position=410,fraction=.5)]
    _write(path,rows)
    for name,manifest in [('fractional_direct',False),('fractional_stream',True)]:
        out=_run(tmp_path,path,name,manifest=manifest)
        for endpoint in ['cpa','tss']:
            cluster=pd.read_csv(out/f'{endpoint}_clusters.tsv',sep='\t')
            assert cluster['n_reads'].tolist()==[1.5]
        assert _mass(out,'cluster_counts.tsv')==1.5
        assert _mass(out,'tss_cluster_counts.tsv')==1.5
