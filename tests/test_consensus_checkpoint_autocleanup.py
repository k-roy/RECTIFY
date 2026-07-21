"""Behavioural tests for the success-path checkpoint auto-cleanup added to
``run_consensus_selection`` (rectify/core/consensus/consensus.py).

Verifies:
  (a) after a SUCCESSFUL consensus with checkpoint_dir set (default
      keep_checkpoints=False), the per-batch checkpoint files are gone;
  (b) with keep_checkpoints=True they remain;
  (c) the final consensus BAM read set is UNCHANGED across
      no-checkpoint / cleanup / keep runs (cleanup must not alter output);
  (d) a consensus whose FINALIZE step fails leaves checkpoints intact
      (resume state preserved — cleanup runs only after a successful sort).
"""

import os
import glob

import pysam
import pytest

from rectify.core.consensus import consensus as c


def _write_multiread_bam(path, n_reads=5):
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'SN': 'chrI', 'LN': 1000}],
    }
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for i in range(n_reads):
            read = pysam.AlignedSegment()
            read.query_name = f'read{i:03d}'
            read.query_sequence = 'ACGT'
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 10 + i * 20
            read.mapping_quality = 60
            read.cigartuples = [(0, 4)]
            read.query_qualities = pysam.qualitystring_to_array('IIII')
            out.write(read)


def _genome():
    return {'chrI': 'ACGT' * 250}  # 1000 bp


def _read_key_set(bam_path):
    keys = set()
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for r in bam.fetch(until_eof=True):
            keys.add((r.query_name, r.reference_start, tuple(r.cigartuples or [])))
    return keys


def _ckpt_files(ckpt_dir):
    batch_bams = glob.glob(os.path.join(ckpt_dir, 'consensus_batch_*.bam'))
    done = glob.glob(os.path.join(ckpt_dir, 'consensus_batch_*.done'))
    js = glob.glob(os.path.join(ckpt_dir, 'consensus_checkpoint.json'))
    return batch_bams, done, js


def _make_inputs(tmp_path):
    bam = tmp_path / 'in.bam'
    _write_multiread_bam(bam, n_reads=5)
    bam_paths = {'minimap2': str(bam), 'uLTRA': str(bam)}
    return bam_paths


def test_success_default_cleans_checkpoints(tmp_path):
    bam_paths = _make_inputs(tmp_path)
    ckpt = tmp_path / 'ckpt_clean'
    out_bam = tmp_path / 'out_clean.bam'

    c.run_consensus_selection(
        bam_paths=bam_paths,
        genome=_genome(),
        output_bam=str(out_bam),
        annotated_junctions=None,
        n_workers=1,
        batch_size=2,           # 5 reads -> 3 batches (2,2,1)
        checkpoint_dir=str(ckpt),
        # keep_checkpoints defaults to False
    )

    assert out_bam.exists(), "final consensus BAM must be written"
    batch_bams, done, js = _ckpt_files(str(ckpt))
    assert batch_bams == [], f"batch BAMs should be removed on success, found {batch_bams}"
    assert done == [], f".done sentinels should be removed on success, found {done}"
    assert js == [], f"consensus_checkpoint.json should be removed on success, found {js}"


def test_keep_checkpoints_retains_resume_state(tmp_path):
    bam_paths = _make_inputs(tmp_path)
    ckpt = tmp_path / 'ckpt_keep'
    out_bam = tmp_path / 'out_keep.bam'

    c.run_consensus_selection(
        bam_paths=bam_paths,
        genome=_genome(),
        output_bam=str(out_bam),
        annotated_junctions=None,
        n_workers=1,
        batch_size=2,
        checkpoint_dir=str(ckpt),
        keep_checkpoints=True,
    )

    assert out_bam.exists()
    batch_bams, done, js = _ckpt_files(str(ckpt))
    assert len(batch_bams) == 3, f"expected 3 retained batch BAMs, found {batch_bams}"
    assert len(done) == 3, f"expected 3 retained .done sentinels, found {done}"
    assert len(js) == 1, f"expected retained consensus_checkpoint.json, found {js}"


def test_output_unchanged_by_cleanup(tmp_path):
    """No-checkpoint, cleanup, and keep runs must all produce the same reads."""
    bam_paths = _make_inputs(tmp_path)

    out_none = tmp_path / 'out_none.bam'
    c.run_consensus_selection(
        bam_paths=bam_paths, genome=_genome(), output_bam=str(out_none),
        annotated_junctions=None, n_workers=1, batch_size=2,
        checkpoint_dir=None,
    )

    out_clean = tmp_path / 'out_clean.bam'
    c.run_consensus_selection(
        bam_paths=bam_paths, genome=_genome(), output_bam=str(out_clean),
        annotated_junctions=None, n_workers=1, batch_size=2,
        checkpoint_dir=str(tmp_path / 'ck_clean'), keep_checkpoints=False,
    )

    out_keep = tmp_path / 'out_keep.bam'
    c.run_consensus_selection(
        bam_paths=bam_paths, genome=_genome(), output_bam=str(out_keep),
        annotated_junctions=None, n_workers=1, batch_size=2,
        checkpoint_dir=str(tmp_path / 'ck_keep'), keep_checkpoints=True,
    )

    k_none = _read_key_set(str(out_none))
    k_clean = _read_key_set(str(out_clean))
    k_keep = _read_key_set(str(out_keep))
    assert k_none == k_clean == k_keep, (
        "cleanup / keep / no-checkpoint outputs must be identical read sets"
    )
    assert len(k_clean) == 5


def test_finalize_failure_preserves_checkpoints(tmp_path, monkeypatch):
    """If the finalize/sort step fails, checkpoints must remain for resume."""
    bam_paths = _make_inputs(tmp_path)
    ckpt = tmp_path / 'ckpt_fail'
    out_bam = tmp_path / 'out_fail.bam'

    def _boom(*args, **kwargs):
        raise RuntimeError("simulated failure during finalize")

    # _copy2_and_fsync is called ONLY in the checkpoint finalize, after all
    # batches are committed and after merge/sort/index succeed, but BEFORE the
    # success-path checkpoint cleanup. Making it raise exercises a finalize
    # failure without disturbing the earlier _ensure_name_sorted (which uses
    # pysam.sort). Checkpoints must survive for resume.
    monkeypatch.setattr(c, '_copy2_and_fsync', _boom)

    with pytest.raises(RuntimeError):
        c.run_consensus_selection(
            bam_paths=bam_paths,
            genome=_genome(),
            output_bam=str(out_bam),
            annotated_junctions=None,
            n_workers=1,
            batch_size=2,
            checkpoint_dir=str(ckpt),
            keep_checkpoints=False,
        )

    batch_bams, done, js = _ckpt_files(str(ckpt))
    assert len(batch_bams) == 3, f"batch BAMs must survive a failed finalize, found {batch_bams}"
    assert len(done) == 3, f".done sentinels must survive a failed finalize, found {done}"
    assert len(js) == 1, f"checkpoint JSON must survive a failed finalize, found {js}"
