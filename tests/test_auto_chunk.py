"""Unit tests for auto-chunk calibration sizing + the UGE header fixes."""
import math

import pytest

from pathlib import Path

from rectify.core.commands.split_command import (
    recommend_chunk_size,
    _uge_headers,
    _run_all_flags_from_split_args,
    _make_submit_script,
)


def _submit(scheduler):
    P = Path
    return _make_submit_script(
        scheduler,
        mpb_script=P('run_array_mapPacBio.sh'),
        others_script=P('run_array_others.sh'),
        merge_aligners_script=P('run_merge_aligners.sh'),
        prescan_script=P('run_prescan.sh'),
        correct_scripts={'minimap2': P('run_array_correct_minimap2.sh'),
                         'deSALT': P('run_array_correct_deSALT.sh')},
        chunk_merge_script=P('run_array_chunk_merge.sh'),
        final_merge_script=P('run_final_merge.sh'),
        consensus_per_chunk_script=P('run_array_consensus_per_chunk.sh'),
        merge_consensus_script=P('run_merge_consensus_chunks.sh'),
        output_dir=P('/tmp/o'),
    )


def test_uge_submit_chain_uses_holdjid_not_sbatch():
    s = _submit('uge')
    assert 'qsub -terse' in s
    assert '-hold_jid' in s
    assert 'sbatch' not in s           # not the SLURM path
    assert 'not yet implemented' not in s
    # merge_aligners holds on BOTH alignment arrays; chunk_merge holds on both correct jobs
    assert '-hold_jid ${MPB_JOB},${OTHERS_JOB}' in s
    assert '${CORRECT_MINIMAP2_JOB},${CORRECT_DESALT_JOB}' in s


def test_slurm_submit_chain_still_uses_sbatch():
    s = _submit('slurm')
    assert 'sbatch' in s and '--dependency=afterok' in s


def _probe(rate_reads_per_min, overhead_s, sizes=(1000, 4000)):
    """Build synthetic probe points for a given marginal rate + fixed overhead."""
    b = 60.0 / rate_reads_per_min          # s per read
    return [(n, overhead_s + b * n) for n in sizes]


def test_drs_like_rate_sizes_to_target():
    # ~375 reads/min (tonight's DRS) -> ~15k reads for a 40-min chunk
    rec = recommend_chunk_size(_probe(375, 30), total_reads=1_000_000, target_minutes=40)
    assert 13_000 <= rec['reads_per_chunk'] <= 17_000
    assert abs(rec['marginal_reads_per_min'] - 375) < 5
    assert 35 <= rec['est_chunk_minutes'] <= 45


def test_cdna_slower_rate_gives_smaller_chunks():
    # ~136 reads/min (tonight's cDNA) -> ~5.4k reads for 40 min — ~3x smaller than DRS
    drs = recommend_chunk_size(_probe(375, 30), 1_000_000, target_minutes=40)
    cdna = recommend_chunk_size(_probe(136, 30), 1_000_000, target_minutes=40)
    assert cdna['reads_per_chunk'] < drs['reads_per_chunk']
    assert 4_000 <= cdna['reads_per_chunk'] <= 7_000


def test_high_overhead_raises_chunk_size_human_index_load():
    # Large fixed overhead (mammalian index load) at a short target -> the
    # overhead-fraction floor pushes chunks LARGER so startup is amortized.
    rec = recommend_chunk_size(_probe(375, overhead_s=300), 1_000_000,
                               target_minutes=10, max_overhead_frac=0.20)
    assert rec['overhead_fraction'] <= 0.21
    assert any('overhead' in n for n in rec['notes'])
    # bigger than the naive marginal-only size at this short target
    naive = (10 * 60) * (375 / 60.0)
    assert rec['reads_per_chunk'] > naive


def test_wall_ceiling_clamp_feasible():
    # Feasible: a huge target is clamped down to the wall cap (chunk count stays
    # well under max_chunks).
    rec = recommend_chunk_size(_probe(375, 30), 1_000_000,
                               target_minutes=10_000, max_chunk_minutes=120)
    assert rec['est_chunk_minutes'] <= 121
    assert rec['n_chunks'] < 500
    assert any('chunk wall' in n for n in rec['notes'])


def test_wall_ceiling_unmeetable_warns():
    # Infeasible: 100M reads can't fit under a 120-min cap within max 500 chunks.
    rec = recommend_chunk_size(_probe(375, 30), 100_000_000,
                               target_minutes=10_000, max_chunk_minutes=120,
                               max_chunks=500)
    assert rec['n_chunks'] == 500
    assert any('cannot keep chunk wall' in n for n in rec['notes'])


def test_chunk_count_clamped():
    rec = recommend_chunk_size(_probe(375, 30), total_reads=5000,
                               target_minutes=40, min_chunks=4)
    assert rec['n_chunks'] >= 4          # min-chunks floor even for tiny inputs


def test_needs_two_points():
    with pytest.raises(ValueError):
        recommend_chunk_size([(1000, 100.0)], total_reads=10000)


def test_uge_header_omits_queue_when_none():
    h = _uge_headers('job', 4, 8, 32, '4:00:00', 'logs', queue=None, pe='shared')
    assert '-q ' not in h
    assert '#$ -pe shared 8' in h
    assert '#$ -t 1-4' in h


def test_uge_header_includes_queue_when_set():
    h = _uge_headers('job', 4, 8, 32, '4:00:00', 'logs', queue='campus.q')
    assert '#$ -q campus.q' in h


def test_flag_reconstruction_drs_5panel():
    import argparse
    ns = argparse.Namespace(
        drs=True, Scer=True, short_read=False,
        other_aligners=['minimap2', 'gapmm2', 'uLTRA', 'deSALT', 'gmap'],
        skip_map_pacbio=False, gmap_db='/db/Scer',
    )
    flags = _run_all_flags_from_split_args(ns)
    assert '--drs' in flags and '--Scer' in flags
    assert '--base-aligners' in flags and 'mapPacBio' in flags and 'minimap2' in flags
    assert '--junction-aligners' in flags and 'uLTRA' in flags and 'gmap' in flags
    assert '--gmap-db' in flags
