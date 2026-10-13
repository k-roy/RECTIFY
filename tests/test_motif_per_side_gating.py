"""`rectify analyze` motif discovery: each side is scanned on its own, against a background set.

Kevin's deliverable spec (TODO_KR_NOTES): de novo motifs 100 nt upstream and 50 nt downstream of the clusters
enriched vs background and depleted vs background, four scans. The callers required BOTH sides to reach 20 clusters,
so one thin side silenced all four (09-11 inventory, gap i). Now a side with >= 20 clusters runs against the
unchanged clusters (subsampled to 5x the side, seed 0), falling back to the other side, then to no background.
"""

from unittest.mock import patch

import pandas as pd

from rectify.core.analyze import motif_discovery as md


def _df(n, start=0):
    return pd.DataFrame({'chrom': ['chrI'] * n, 'strand': ['+'] * n,
                         'modal_position': list(range(start, start + n))}, index=range(start, start + n))


def _fake_extract(df, genome_fasta, upstream_window=100, downstream_window=50, region='both'):
    seqs = [{'cluster_id': i, 'sequence': 'ACGT' * 5} for i in df.index]
    return {'upstream': seqs, 'downstream': seqs}


def _run(tmp_path, enriched, depleted, background, **kw):
    seen = {}

    def _fake_discovery(foreground_sequences, background_sequences=None, output_dir='.', region_name='r', **_):
        seen[region_name] = (len(foreground_sequences), len(background_sequences or []))
        return {'success': True, 'motifs_found': 1, 'n_foreground': len(foreground_sequences),
                'n_background': len(background_sequences or [])}

    with patch.object(md, 'extract_sequences_around_clusters', side_effect=_fake_extract), \
         patch.object(md, 'run_motif_discovery', side_effect=_fake_discovery):
        res = md.run_differential_motif_analysis(enriched, depleted, 'genome.fa', str(tmp_path / 'motifs'),
                                                 n_parallel=1, background_clusters=background, **kw)
    return res, seen


def test_a_thin_side_no_longer_silences_the_other(tmp_path):
    res, seen = _run(tmp_path, _df(25), _df(5, 1000), _df(400, 2000), min_foreground=20)
    assert set(seen) == {'enriched_upstream', 'enriched_downstream'}
    assert seen['enriched_upstream'] == (25, 125)                 # background capped at 5x the side
    assert res['depleted_upstream']['success'] is False
    assert 'skipped' in res['depleted_downstream']['error']


def test_both_sides_run_against_the_background(tmp_path):
    _, seen = _run(tmp_path, _df(30), _df(40, 1000), _df(1000, 3000), min_foreground=20)
    assert seen == {'enriched_upstream': (30, 150), 'enriched_downstream': (30, 150),
                    'depleted_upstream': (40, 200), 'depleted_downstream': (40, 200)}


def test_without_a_background_the_sides_are_contrasted_as_before(tmp_path):
    _, seen = _run(tmp_path, _df(30), _df(40, 1000), None, min_foreground=20)
    assert seen['enriched_upstream'] == (30, 40) and seen['depleted_downstream'] == (40, 30)


def test_a_lone_side_without_any_background_uses_the_shuffled_control(tmp_path):
    _, seen = _run(tmp_path, _df(30), _df(3, 1000), _df(4, 2000), min_foreground=20)
    assert seen == {'enriched_upstream': (30, 0), 'enriched_downstream': (30, 0)}
