"""WP0 remnants of the run-all audits (planning 831 / 832 / 833 / 834 / 835 §4).

Two defects that survived the 836 merge, each of which fails SILENTLY:

* ``_ALIGNER_NAMES`` listed no COMPASS arm, so ``_collect_per_aligner_bams`` returned the bbmap BAM
  only for every ``run-all --short-read`` run while the 5- or 7-way consensus sat unused
  (planning 833 C-2).
* the sample id was re-derived from ``input_path`` AFTER Step 0 had replaced it, renaming every
  ``run-all --ONT-cDNA`` sample to ``stage1_consensus`` (planning 831).

The two already-landed WP0 items are asserted here too, so a future edit cannot quietly undo them:
``--short-read`` must not be a protocol (832 G-0) and ``short_read``/``netseq`` must reach
``correct``'s hand-built Namespace (832 G-1 / 833 C-1 / 834 §6.4).
"""
import argparse

import pytest

from rectify.core.commands.align_command import (
    COMPASS_PE_ALIGNERS,
    COMPASS_SE_ALIGNERS,
)
from rectify.core.commands.run.helpers import _ALIGNER_NAMES, _collect_per_aligner_bams
from rectify.core.commands.run.single_sample import _canonical_sample_id
from rectify.core.commands.run_command import (
    _CORRECT_FORWARDED_FLAGS,
    _MODALITY_FLAGS,
    _PROTOCOL_FLAGS,
    _validate_protocol_flags,
)


# ── 833 C-2: every COMPASS arm must be collectable ──────────────────────────────────────────
@pytest.mark.parametrize("aligner", sorted(set(COMPASS_PE_ALIGNERS) | set(COMPASS_SE_ALIGNERS)))
def test_compass_arms_are_in_aligner_names(aligner):
    assert aligner in _ALIGNER_NAMES, (
        f"{aligner} is in the COMPASS panel but not in _ALIGNER_NAMES, so its BAM is invisible "
        "to correction and the merged TSV (planning 833 C-2)"
    )


def test_collect_per_aligner_bams_finds_a_star_arm(tmp_path):
    """The regression as it actually presented: a STAR BAM on disk, invisible to the collector."""
    for aligner in ('bbmap', 'STAR_default', 'HISAT2_noncanonical'):
        (tmp_path / f"s1.{aligner}.bam").write_bytes(b"x" * 4096)
    found = _collect_per_aligner_bams('s1', tmp_path)
    assert set(found) == {'bbmap', 'STAR_default', 'HISAT2_noncanonical'}


def test_aligner_names_has_no_duplicates():
    assert len(_ALIGNER_NAMES) == len(set(_ALIGNER_NAMES))


# ── 831: the sample id must survive a Step-0 input rewrite ──────────────────────────────────
@pytest.mark.parametrize("name,expected", [
    ("wt_rep3.fastq.gz", "wt_rep3"),
    ("wt_rep3.fq.gz", "wt_rep3"),
    ("wt_rep3.fastq", "wt_rep3"),
    ("sample_dorado.bam", "sample_dorado"),
    ("a.b.fastq.gz", "a.b"),
    ("reads.FASTQ.GZ", "reads"),
])
def test_canonical_sample_id(name, expected):
    assert _canonical_sample_id(f"/tmp/dir/{name}") == expected


def test_canonical_sample_id_is_taken_before_step0_rewrite():
    """The 831 defect in one line: the Step-0 OUTPUT name must never become the sample id."""
    original = "/data/wt_rep3.fastq.gz"
    step0_drs = "/scratch/drs_trim/wt_rep3_trimmed.fastq.gz"
    step0_cdna = "/scratch/cdna/stage1/stage1_consensus.fastq.gz"
    assert _canonical_sample_id(original) == "wt_rep3"
    # Both Step-0 products derive a DIFFERENT id — which is why it must be computed once, up front.
    assert _canonical_sample_id(step0_drs) == "wt_rep3_trimmed"
    assert _canonical_sample_id(step0_cdna) == "stage1_consensus"


def test_single_sample_does_not_rederive_sample_id_after_step0():
    """Guard the call-site contract, not just the helper: one assignment, before Step 0."""
    import inspect

    from rectify.core.commands.run import single_sample

    src = inspect.getsource(single_sample._run_single_sample)
    assert "input_path.stem.replace('.fastq'" not in src, (
        "sample_id re-derived from input_path inside _run_single_sample — Step 0 may already "
        "have replaced it (planning 831)"
    )
    assert src.count("sample_id = _canonical_sample_id(input_path)") == 1


# ── 832 G-0 / G-1, already on master (a9a1943) — locked down here ───────────────────────────
def test_short_read_is_a_modality_not_a_protocol():
    assert 'short_read' not in [dest for _flag, dest in _PROTOCOL_FLAGS]
    assert ('--short-read', 'short_read') in _MODALITY_FLAGS
    args = argparse.Namespace(short_read=True, dT_primed_cDNA=True, drs=False, ONT_cDNA=False)
    assert _validate_protocol_flags(args) is None, (
        "run-all --short-read --dT-primed-cDNA is the documented QuantSeq invocation (832 G-0)"
    )


def test_correct_namespace_forwards_the_modality_and_protocol_flags():
    for dest in ('short_read', 'netseq', 'drs', 'dT_primed_cDNA', 'ONT_cDNA'):
        assert dest in _CORRECT_FORWARDED_FLAGS
