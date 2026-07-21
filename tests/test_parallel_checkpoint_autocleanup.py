"""Behavioural tests for the success-path region-checkpoint auto-cleanup added
to ``process_bam_streaming_parallel`` (rectify/core/bam/parallel.py) — the
mechanism that drives the chunked ``rectify correct --checkpoint-dir`` panel.

Verifies:
  (a) after a SUCCESSFUL parallel-streaming run with checkpoint_dir set (default
      keep_checkpoints=False), the per-region checkpoint files are gone;
  (b) with keep_checkpoints=True they remain;
  (c) the output TSV rows are IDENTICAL across cleanup / keep runs.
"""

import glob
import os
from pathlib import Path

import pysam
import pytest

from rectify.core.bam import parallel as p

_CHROM_LEN = 2000


def _make_reference(tmp_path: Path) -> Path:
    seq = list("X" * _CHROM_LEN)
    for i, b in enumerate("ACGTACGTAC"):
        seq[1000 + i] = b
    for i in range(10):
        seq[1010 + i] = "A"
    for i in range(10):
        seq[1020 + i] = "C"
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chrT\n{''.join(seq)}\n")
    pysam.faidx(str(fa))
    return fa


def _make_small_bam(tmp_path: Path) -> Path:
    bam_path = tmp_path / "input.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrT", "LN": _CHROM_LEN}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for i in range(5):
            r = pysam.AlignedSegment(bam.header)
            r.query_name = f"r{i}_test"
            r.query_sequence = "ACGTACGTAC" + "A" * 5
            r.flag = 0
            r.reference_id = 0
            r.reference_start = 1000
            r.mapping_quality = 60
            r.cigartuples = [(0, 15)]
            bam.write(r)
    pysam.index(str(bam_path))
    return bam_path


def _region_files(ckpt_dir):
    tsv = glob.glob(os.path.join(ckpt_dir, 'region_*.tsv'))
    stats = glob.glob(os.path.join(ckpt_dir, 'region_*.stats.json'))
    done = glob.glob(os.path.join(ckpt_dir, 'region_*.done'))
    return tsv, stats, done


def _run(tmp_path, ckpt, out_tsv, keep):
    return p.process_bam_streaming_parallel(
        bam_path=str(_BAM),
        genome_path=str(_REF),
        output_path=str(out_tsv),
        n_threads=2,
        show_progress=False,
        checkpoint_dir=str(ckpt),
        keep_checkpoints=keep,
    )


# Module-scoped inputs so both variants share one BAM/reference.
_REF = None
_BAM = None


@pytest.fixture(autouse=True, scope="module")
def _inputs(tmp_path_factory):
    global _REF, _BAM
    d = tmp_path_factory.mktemp("pcac_inputs")
    _REF = _make_reference(d)
    _BAM = _make_small_bam(d)
    yield


def _read_tsv_rows(path):
    with open(path) as fh:
        return [ln for ln in fh]


def test_parallel_success_default_cleans_region_checkpoints(tmp_path):
    ckpt = tmp_path / "ck_clean"
    out = tmp_path / "out_clean.tsv"
    _run(tmp_path, ckpt, out, keep=False)

    assert out.exists(), "output TSV must be written"
    tsv, stats, done = _region_files(str(ckpt))
    assert tsv == [], f"region_*.tsv should be removed on success, found {tsv}"
    assert stats == [], f"region_*.stats.json should be removed on success, found {stats}"
    assert done == [], f"region_*.done should be removed on success, found {done}"
    assert not (ckpt / "rescue_scan.pkl").exists()


def test_parallel_keep_checkpoints_retains_region_files(tmp_path):
    ckpt = tmp_path / "ck_keep"
    out = tmp_path / "out_keep.tsv"
    _run(tmp_path, ckpt, out, keep=True)

    assert out.exists()
    tsv, stats, done = _region_files(str(ckpt))
    assert len(tsv) >= 1, f"region_*.tsv should be retained with keep, found {tsv}"
    assert len(done) >= 1, f"region_*.done should be retained with keep, found {done}"
    assert len(done) == len(tsv), "one .done per region tsv expected"


def test_parallel_rebuild_failure_preserves_region_checkpoints(tmp_path, monkeypatch):
    """If the final rebuild fails, region checkpoints must be preserved for resume.

    Directly proves the success-only placement: cleanup lives AFTER
    _rebuild_output_from_region_files, so a rebuild that raises must leave the
    committed region files on disk (a requeued/preempted chunk resumes).
    """
    ckpt = tmp_path / "ck_fail"
    out = tmp_path / "out_fail.tsv"

    def _boom(*args, **kwargs):
        raise RuntimeError("simulated rebuild failure")

    monkeypatch.setattr(p, '_rebuild_output_from_region_files', _boom)

    with pytest.raises(RuntimeError):
        p.process_bam_streaming_parallel(
            bam_path=str(_BAM),
            genome_path=str(_REF),
            output_path=str(out),
            n_threads=2,
            show_progress=False,
            checkpoint_dir=str(ckpt),
            keep_checkpoints=False,
        )

    tsv, stats, done = _region_files(str(ckpt))
    assert len(tsv) >= 1, f"region_*.tsv must survive a failed rebuild, found {tsv}"
    assert len(done) >= 1, f"region_*.done must survive a failed rebuild, found {done}"


def test_parallel_output_unchanged_by_cleanup(tmp_path):
    ckpt_c = tmp_path / "ck_c"
    out_c = tmp_path / "oc.tsv"
    _run(tmp_path, ckpt_c, out_c, keep=False)

    ckpt_k = tmp_path / "ck_k"
    out_k = tmp_path / "ok.tsv"
    _run(tmp_path, ckpt_k, out_k, keep=True)

    assert _read_tsv_rows(str(out_c)) == _read_tsv_rows(str(out_k)), (
        "cleanup and keep runs must produce identical output TSV rows"
    )
