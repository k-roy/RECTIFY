"""Actual DESeq2 model/statistics honor the same explicitly limited inference."""
import sys

import joblib
import numpy as np
import pandas as pd

from rectify.core.analyze import deseq2


class StopBeforeExcessWorkers(BaseException):
    """Never let a regression test create an over-budget pool."""


def model_inputs(seed=770924):
    rng=np.random.default_rng(seed)
    samples=['base_'+str(i) for i in range(3)]+['drug_'+str(i) for i in range(3)]+['rescue_'+str(i) for i in range(3)]
    means=np.repeat(rng.uniform(40,260,size=(96,1)),9,axis=1)
    means[:18,3:6]*=2.7
    means[18:32,3:6]*=.45
    means[:18,6:]*=1.6
    means[32:46,6:]*=.5
    values=rng.negative_binomial(12,12/(12+means))+1
    counts=pd.DataFrame(values,index=[f'cpu_feature_{i:03d}' for i in range(96)],columns=samples)
    counts.index.name='synthetic_feature'
    metadata=pd.DataFrame({'condition':['base_line']*3+['drug_exposed']*3+['rescue_arm']*3},index=samples)
    return counts,metadata


def observed_fit(counts,metadata):
    constructors,dispatches=[],[]
    previous=sys.getprofile()
    def profile(frame,event,arg):
        if event=='call' and frame.f_code is joblib.Parallel.__call__.__code__:
            n_jobs=frame.f_locals['self'].n_jobs
            dispatches.append(n_jobs)
            if n_jobs!=1:
                raise StopBeforeExcessWorkers(f'Blocked actual joblib n_jobs={n_jobs}')
        if event=='return' and frame.f_code is deseq2.DeseqStats.__init__.__code__:
            obj=frame.f_locals['self']
            constructors.append(dict(dataset_n_cpus=obj.dds.inference.n_cpus,
                                     stats_n_cpus=obj.inference.n_cpus,
                                     same_inference=obj.inference is obj.dds.inference,
                                     contrast=list(obj.contrast),
                                     fitted_features=len(obj.dds.var_names)))
            if obj.inference.n_cpus!=1:
                raise StopBeforeExcessWorkers('Blocked excess-CPU stats before summary')
    sys.setprofile(profile)
    try:
        result=deseq2._run_deseq2(counts,metadata,'base_line',1)
    finally:
        sys.setprofile(previous)
    return result,constructors,dispatches


def test_actual_fit_and_every_contrast_use_one_inference_cpu():
    counts,metadata=model_inputs()
    before=counts.copy(deep=True)
    results,constructors,dispatches=observed_fit(counts,metadata)
    assert set(results)=={'drug_exposed','rescue_arm'}
    assert len(constructors)==2
    assert all(row['same_inference'] and row['dataset_n_cpus']==row['stats_n_cpus']==1 for row in constructors)
    assert [row['contrast'] for row in constructors]==[['condition','drug-exposed','base-line'],['condition','rescue-arm','base-line']]
    assert dispatches and set(dispatches)=={1}
    pd.testing.assert_frame_equal(counts,before)
    for table in results.values():
        assert table.index.equals(counts.index)
        assert table.index.name=='synthetic_feature'
        assert len(table)==96
        assert np.isfinite(table['log2FoldChange']).all()
        assert np.isfinite(table['stat']).all()
        assert table['pvalue'].notna().all()
        assert {'significant_padj05','significant_padj01','direction'} <= set(table.columns)


def test_invalid_reference_refuses_before_any_model_or_dispatch():
    counts,metadata=model_inputs()
    calls=[];previous=sys.getprofile()
    def profile(frame,event,arg):
        if event=='call' and frame.f_code is joblib.Parallel.__call__.__code__:
            calls.append(frame.f_locals['self'].n_jobs)
            raise StopBeforeExcessWorkers('Unexpected dispatch for missing reference')
    sys.setprofile(profile)
    try:
        try:deseq2._run_deseq2(counts,metadata,'not_present',1)
        except ValueError as exc:assert 'not found' in str(exc)
        else:raise AssertionError('Missing reference must refuse')
    finally:sys.setprofile(previous)
    assert calls==[]
