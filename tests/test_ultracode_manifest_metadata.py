"""The exported manifest design must match the design consumed by analysis."""
import argparse
import csv
from contextlib import redirect_stderr, redirect_stdout

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import pytest

from rectify.core.bam.bam_processor import correct_read_3prime
from rectify.core.bam.output import write_output_tsv
from rectify.core.commands.analyze_command import create_analyze_parser, run_analyze
from rectify.utils.genome import register_genome_contigs


@pytest.mark.parametrize('mode', ['explicit', 'reference', 'fallback'])
def test_native_manifest_exports_final_conditions_and_preserves_counts(tmp_path, mode):
    genome = {'chrV': 'CGTC' * 600}
    register_genome_contigs(genome)
    header = pysam.AlignmentHeader.from_references(['chrV'], [2400])
    samples = ['beta_2', 'alpha_1'] if mode != 'fallback' else ['control_2', 'control_1']
    expected_condition = 'Verified_Group' if mode != 'fallback' else 'control'
    manifest_rows = []
    for sample in samples:
        rows = []
        for i, (start, flag) in enumerate([(300, 0), (900, 16), (1500, 0)]):
            read = pysam.AlignedSegment(header)
            read.query_name = f'{sample}_{i}'
            read.reference_id = 0
            read.reference_start = start
            read.flag = flag
            read.cigarstring = '40M'
            read.query_sequence = genome['chrV'][start:start + 40]
            read.query_qualities = [32] * 40
            read.set_tag('Xz', 1)
            original = read.to_string()
            corrected = correct_read_3prime(read, genome, apply_atract=False, apply_3ss_rescue=False)
            assert read.to_string() == original and len(corrected) == 1
            rows.extend(corrected)
        path = tmp_path / (sample + '.tsv')
        write_output_tsv(rows, str(path))
        row = {'sample_id': sample, 'path': str(path)}
        if mode != 'fallback':
            row['condition'] = expected_condition
        manifest_rows.append(row)
    manifest = tmp_path / 'manifest.tsv'
    pd.DataFrame(manifest_rows).to_csv(manifest, sep='\t', index=False)
    sources = [manifest, *(tmp_path / (s + '.tsv') for s in samples)]
    original_inputs = {p: p.read_bytes() for p in sources}
    out = tmp_path / 'analysis'
    parser = argparse.ArgumentParser()
    create_analyze_parser(parser.add_subparsers(dest='command'))
    argv = ['analyze', '--manifest', str(manifest), '-o', str(out), '--threads', '1',
            '--min-reads', '1', '--min-cluster-samples', '1', '--include-mito',
            '--no-genomic-distribution', '--gene-attribution-mode', 'none']
    if mode == 'reference':
        argv += ['--reference', 'verified_group']
    try:
        with (tmp_path / 'command.log').open('w') as log, redirect_stdout(log), redirect_stderr(log):
            assert run_analyze(parser.parse_args(argv)) == 0
    finally:
        plt.close('all')
    # Read the public TSV explicitly: its historical index and sample column
    # both have the header 'sample', so preserve that wire format here.
    with (out / 'sample_metadata.tsv').open() as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    assert rows[0] == ['sample', 'sample', 'condition', 'is_control']
    assert {r[0]: r[2] for r in rows[1:]} == dict.fromkeys(samples, expected_condition)
    assert all(r[0] == r[1] for r in rows[1:])
    for filename in ['cluster_counts.tsv', 'tss_cluster_counts.tsv']:
        counts = pd.read_csv(out / filename, sep='\t', index_col=0)
        assert counts.shape == (3, 2)
        assert counts.sum().to_dict() == dict.fromkeys(samples, 3.0)
        assert list(counts.columns) == [r[0] for r in rows[1:]]
    summary = pd.read_csv(out / 'analysis_summary.tsv', sep='\t')
    assert any(expected_condition in str(r) for r in summary.to_dict('records'))
    assert '</html>' in (out / 'report.html').read_text().lower()
    assert 'Skipping shift analysis (need >=2 conditions)' in (tmp_path / 'command.log').read_text()
    assert not list(out.rglob('deseq2_*.tsv'))
    assert original_inputs == {p: p.read_bytes() for p in sources}
