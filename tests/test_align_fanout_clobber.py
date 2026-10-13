#!/usr/bin/env python3
"""The per-aligner fan-out must not silently overwrite <prefix>.multialigned.bam (KNOWN_ISSUES, 2026-07-21).

`rectify align` writes its single-aligner copy to <prefix>.multialigned.bam whatever the aligner, so one align job
per aligner with a shared --prefix and -o left that file holding whichever aligner finished last, and nothing
errored (Sumner session C, 2026-07-21; still open on master at the 2026-10-12 audit). A run now refuses when
another aligner's output sits under the same prefix or the align sidecar records a different aligner set.
"""

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import pysam
import pytest

from rectify.core.commands.align_command import _fanout_conflict


def _touch(p):
    Path(p).write_bytes(b'x')


# --------------------------------------------------------------------------- unit: the guard itself

def test_own_outputs_are_no_conflict(tmp_path):
    _touch(tmp_path / 's.minimap2.bam')
    _touch(tmp_path / 's.overhang_resolver.bam')
    assert _fanout_conflict(tmp_path, 's', {'minimap2', 'overhang_resolver'}) is None


def test_another_aligners_bam_is_a_conflict(tmp_path):
    _touch(tmp_path / 's.minimap2.bam')
    _touch(tmp_path / 's.uLTRA.bam')
    msg = _fanout_conflict(tmp_path, 's', {'minimap2', 'overhang_resolver'})
    assert msg and 'uLTRA' in msg


def test_a_kept_provenance_file_counts_as_that_aligners_output(tmp_path):
    # Chains that delete per-aligner BAMs after consensus keep <prefix>.<aligner>.bam.prov.json.
    _touch(tmp_path / 's.deSALT.bam.prov.json')
    msg = _fanout_conflict(tmp_path, 's', {'gapmm2'})
    assert msg and 'deSALT' in msg


def test_other_prefixes_do_not_count(tmp_path):
    _touch(tmp_path / 'other.uLTRA.bam')
    assert _fanout_conflict(tmp_path, 's', {'gapmm2'}) is None


def test_sidecar_recording_another_aligner_set_is_a_conflict(tmp_path):
    (tmp_path / 's.align.provenance.json').write_text(
        json.dumps({'stage': 'align', 'stats': {'aligners': ['deSALT', 'minimap2', 'uLTRA']}}))
    msg = _fanout_conflict(tmp_path, 's', {'minimap2', 'overhang_resolver'})
    assert msg and 'deSALT' in msg and 'uLTRA' in msg


def test_sidecar_of_the_same_aligner_is_a_rerun(tmp_path):
    (tmp_path / 's.align.provenance.json').write_text(json.dumps({'stats': {'aligners': ['gapmm2']}}))
    assert _fanout_conflict(tmp_path, 's', {'gapmm2'}) is None


def test_sidecar_without_an_aligner_list_is_not_judged(tmp_path):
    (tmp_path / 's.align.provenance.json').write_text(json.dumps({'stats': {'wall_seconds': 3.0}}))
    assert _fanout_conflict(tmp_path, 's', {'gapmm2'}) is None


# --------------------------------------------------------------------------- end to end through run_align

def _write_bam(path, name):
    header = {'HD': {'VN': '1.6', 'SO': 'unsorted'}, 'SQ': [{'SN': 'chrI', 'LN': 1000}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        r = pysam.AlignedSegment(out.header)
        r.query_name = name
        r.query_sequence = 'ACGTACGTAC'
        r.flag = 0
        r.reference_id = 0
        r.reference_start = 100
        r.mapping_quality = 60
        r.cigartuples = [(0, 10)]
        r.query_qualities = pysam.qualitystring_to_array('I' * 10)
        out.write(r)


def _args(tmp_path, aligner):
    for f in ('reads.fastq', 'genome.fa', 'anno.gff'):
        (tmp_path / f).touch()
    ns = argparse.Namespace()
    ns.threads = 2
    ns.parallel_aligners = False
    ns.reads = tmp_path / 'reads.fastq'
    ns.genome = tmp_path / 'genome.fa'
    ns.annotation = tmp_path / 'anno.gff'
    ns.output_dir = tmp_path
    ns.prefix = 's'
    ns.minimap2_path = 'minimap2'
    ns.mapPacBio_path = 'mapPacBio.sh'
    ns.gapmm2_path = 'gapmm2'
    ns.ultra_path = 'uLTRA'
    ns.desalt_path = 'deSALT'
    ns.junc_bonus = 9
    ns.no_consensus = False
    ns.junction_aligners = []
    ns.verbose = False
    ns.junc_bed = None
    ns.trust_existing_bams = True
    ns.aligners = [aligner]
    return ns


@pytest.fixture()
def _fake_aligners():
    def _make(name):
        def _runner(*args, output_bam=None, **kwargs):
            _write_bam(output_bam, f'read_from_{name}')
            return output_bam
        return _runner
    patches = [patch(f'rectify.core.align.multi_aligner.{fn}', side_effect=_make(name))
               for fn, name in (('run_gapmm2', 'gapmm2'), ('run_map_pacbio', 'mapPacBio'))]
    patches.append(patch('rectify.core.align.multi_aligner.check_aligner_available', return_value=True))
    for p in patches:
        p.start()
    yield
    for p in patches:
        p.stop()


def _names(bam):
    with pysam.AlignmentFile(str(bam), 'rb', check_sq=False) as fh:
        return [r.query_name for r in fh]


def test_second_aligner_job_on_a_shared_prefix_refuses_and_keeps_the_first_file(tmp_path, _fake_aligners):
    from rectify.core.commands.align_command import run_align
    assert run_align(_args(tmp_path, 'gapmm2')) in (0, None)
    dest = tmp_path / 's.multialigned.bam'
    assert _names(dest) == ['read_from_gapmm2']
    before = dest.read_bytes()
    rc = run_align(_args(tmp_path, 'mapPacBio'))
    assert rc == 1
    assert dest.read_bytes() == before                     # untouched: it still holds gapmm2's records
    assert (tmp_path / 's.mapPacBio.bam').exists()         # the second job's own BAM survives


def test_rerunning_the_same_aligner_still_overwrites(tmp_path, _fake_aligners):
    from rectify.core.commands.align_command import run_align
    assert run_align(_args(tmp_path, 'gapmm2')) in (0, None)
    assert run_align(_args(tmp_path, 'gapmm2')) in (0, None)
    side = json.loads((tmp_path / 's.align.provenance.json').read_text())
    assert side['stats']['aligners'] == ['gapmm2']
