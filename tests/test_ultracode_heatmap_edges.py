"""Undefined similarity and singleton plots must not abort valid analysis."""
import argparse
import csv
from contextlib import redirect_stdout, redirect_stderr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from scipy.cluster.hierarchy import leaves_list, linkage

from rectify.core.analyze.heatmap import plot_sample_heatmap, plot_cluster_heatmap
from rectify.core.analyze import heatmap as heatmap_module
from rectify.core.bam.bam_processor import correct_read_3prime
from rectify.core.bam.output import write_output_tsv
from rectify.core.commands.analyze_command import create_analyze_parser, run_analyze
from rectify.utils.genome import register_genome_contigs


def _heatmap_axis(fig, title):
    return next(ax for ax in fig.axes if ax.get_title().startswith(title))


def test_undefined_correlations_remain_missing_in_input_order(tmp_path):
    cases = {
        'constant': pd.DataFrame({'z_last': [2.5]*3, 'a_first': [7.5]*3}),
        'zero': pd.DataFrame({'z_last': [0.]*3, 'a_first': [0.]*3}),
        'mixed': pd.DataFrame({'z_last': [3.]*3, 'a_first': [1.,4.,9.]}),
        'one_feature': pd.DataFrame({'z_last': [3.], 'a_first': [8.]}),
        'one_cell': pd.DataFrame({'only': [3.]}),
    }
    for name, counts in cases.items():
        original = counts.copy(deep=True)
        expected = np.log2(counts+1).corr()
        path = tmp_path/(name+'.png')
        fig = plot_sample_heatmap(counts, output_path=str(path), title=name, figsize=(5,4))
        try:
            ax = _heatmap_axis(fig, name)
            observed = np.ma.asarray(ax.collections[0].get_array()).reshape(expected.shape)
            assert np.array_equal(np.ma.getmaskarray(observed), expected.isna().to_numpy())
            np.testing.assert_allclose(observed.compressed(), expected.to_numpy()[np.isfinite(expected)], rtol=0, atol=0)
            assert [x.get_text() for x in ax.get_xticklabels()] == list(counts.columns)
            assert [x.get_text() for x in ax.get_yticklabels()] == list(counts.columns)
            assert 'Undefined cells' in ax.get_title() and 'input order' in ax.get_title()
            assert path.stat().st_size > 1000
            pd.testing.assert_frame_equal(counts, original)
        finally:
            plt.close(fig)


def test_finite_legacy_clustering_and_values_are_unchanged(tmp_path):
    counts = pd.DataFrame({'z':[1.,8.,3.,19.,7.], 'a':[17.,2.,11.,5.,4.],
                           'q':[2.,7.,4.,16.,9.], 'b':[6.,12.,20.,1.,3.]})
    original = counts.copy(deep=True)
    expected = np.log2(counts+1).corr()
    order = leaves_list(linkage(expected.to_numpy(), method='average', metric='euclidean'))
    fig = plot_sample_heatmap(counts, title='finite', output_path=str(tmp_path/'finite.png'))
    try:
        ax = _heatmap_axis(fig, 'finite')
        assert ax.get_title() == 'finite'
        assert [x.get_text() for x in ax.get_xticklabels()] == list(counts.columns[order])
        assert [x.get_text() for x in ax.get_yticklabels()] == list(counts.columns[order])
        values = np.ma.asarray(ax.collections[0].get_array()).reshape(expected.shape)
        assert not np.ma.getmaskarray(values).any()
        np.testing.assert_allclose(values, expected.to_numpy()[np.ix_(order,order)], rtol=0, atol=0)
        pd.testing.assert_frame_equal(counts, original)
    finally:
        plt.close(fig)


def test_missingness_title_stays_above_group_color_strip(tmp_path):
    counts=pd.DataFrame({'late':[4.]*4,'early':[1.,3.,6.,12.],'middle':[12.,7.,3.,2.]})
    metadata=pd.DataFrame({'group':['one','two','one']},index=counts.columns)
    for size in [(5,4),(6,5),(12,10)]:
        fig=plot_sample_heatmap(counts,metadata,'group',title='metadata',figsize=size,
                                output_path=str(tmp_path/f'metadata_{size[0]}.png'))
        try:
            fig.canvas.draw()
            renderer=fig.canvas.get_renderer()
            ax=_heatmap_axis(fig,'metadata')
            strip=next(a for a in fig.axes if [t.get_text() for t in a.get_yticklabels()]==['group'])
            assert ax.title.get_window_extent(renderer).y0 >= strip.get_window_extent(renderer).y1
            assert [t.get_text() for t in ax.get_xticklabels()]==list(counts.columns)
        finally:
            plt.close(fig)


def test_singleton_sample_and_feature_axes_render(tmp_path):
    cases = {'one_feature':pd.DataFrame({'B':[2.], 'A':[7.]}),
             'one_sample':pd.DataFrame({'only':[1.,4.,9.]}),
             'one_cell':pd.DataFrame({'only':[0.]})}
    for name,counts in cases.items():
        original=counts.copy(deep=True)
        for normalize in ['zscore','log','none']:
            fig=plot_cluster_heatmap(counts, normalize=normalize, cluster_rows=True,
                                     cluster_cols=True, title=name, figsize=(5,4),
                                     output_path=str(tmp_path/(name+'.'+normalize+'.png')))
            try:
                ax=_heatmap_axis(fig,name)
                values=np.ma.asarray(ax.collections[0].get_array())
                assert values.size==counts.size and np.isfinite(values).all()
                pd.testing.assert_frame_equal(counts,original)
            finally:
                plt.close(fig)
    fig=plot_sample_heatmap(cases['one_sample'],title='one_sample')
    try:
        values=np.ma.asarray(_heatmap_axis(fig,'one_sample').collections[0].get_array())
        assert values.size==1 and float(values.reshape(-1)[0])==1.
    finally:
        plt.close(fig)
    assert plot_sample_heatmap(pd.DataFrame()) is None
    assert plot_cluster_heatmap(pd.DataFrame()) is None


def test_single_sample_manual_normalization_matches_available_backend(tmp_path):
    # Exercise the real optional-dependency fallback, not a numeric score mock.
    previous=heatmap_module.SKLEARN_AVAILABLE
    heatmap_module.SKLEARN_AVAILABLE=False
    try:
        counts=pd.DataFrame({'only':[0.,2.,8.]})
        fig=plot_cluster_heatmap(counts,normalize='zscore',cluster_cols=True,
                                 title='manual',output_path=str(tmp_path/'manual.png'))
        try:
            values=np.ma.asarray(_heatmap_axis(fig,'manual').collections[0].get_array())
            np.testing.assert_array_equal(values,np.zeros_like(values))
        finally:
            plt.close(fig)
    finally:
        heatmap_module.SKLEARN_AVAILABLE=previous


def _actual_manifest(tmp_path, chrom, n_features):
    genome={chrom:'CGTC'*600}
    register_genome_contigs(genome)
    header=pysam.AlignmentHeader.from_references([chrom],[2400])
    manifest_rows=[]
    for sample in ['WT_2','WT_1']:
        rows=[]
        with pysam.AlignmentFile(str(tmp_path/(sample+'.bam')),'wb',header=header) as bam:
            for i,(start,flag) in enumerate([(300,0),(900,16),(1500,0)][:n_features]):
                r=pysam.AlignedSegment(header);r.query_name=f'{sample}_{i}';r.reference_id=0
                r.reference_start=start;r.flag=flag;r.cigarstring='40M'
                r.query_sequence=genome[chrom][start:start+40];r.query_qualities=[32]*40;r.set_tag('Xz',1)
                original=r.to_string();result=correct_read_3prime(r,genome,apply_atract=False,apply_3ss_rescue=False)
                assert len(result)==1 and result[0]['chrom']==chrom
                assert r.to_string()==original and result[0]['fraction']==1.
                rows+=result;bam.write(r)
        path=tmp_path/(sample+'.tsv');write_output_tsv(rows,str(path))
        manifest_rows.append(dict(sample_id=sample,path=str(path),condition='WT'))
    manifest=tmp_path/'samples.tsv'
    with manifest.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=['sample_id','path','condition'],delimiter='\t')
        writer.writeheader();writer.writerows(manifest_rows)
    parser=argparse.ArgumentParser();create_analyze_parser(parser.add_subparsers(dest='command'))
    out=tmp_path/'analysis'
    argv=['analyze','--manifest',str(manifest),'-o',str(out),'--min-reads','1',
          '--min-cluster-samples','1','--threads','1','--include-mito','--no-genomic-distribution',
          '--gene-attribution-mode','none','--sample-sets','{"review":["WT"]}']
    with (tmp_path/'analyze.log').open('w') as log,redirect_stdout(log),redirect_stderr(log):
        assert run_analyze(parser.parse_args(argv))==0
    for fn in ['cluster_counts.tsv','tss_cluster_counts.tsv']:
        counts=pd.read_csv(out/fn,sep='\t',index_col=0)
        assert counts.shape==(n_features,2) and (counts==1.).all().all()
        assert counts.to_numpy().sum()==2*n_features
    assert (out/'report.html').is_file()
    assert (out/'plots/sample_heatmap.png').stat().st_size>1000
    plt.close('all')


def test_actual_correction_manifest_completes_for_constant_counts(tmp_path):
    for chrom in ['chrV','chr5']:
        for n_features in [1,3]:
            path=tmp_path/(chrom+'_'+str(n_features));path.mkdir()
            _actual_manifest(path,chrom,n_features)
