"""GitHub #4 — the region-parallel Stage-1 path must resolve calmd "=" placeholders.

`rectify align` runs `samtools calmd -e`, which rewrites every reference-matching M base in SEQ to
"=".  `restore_eq_seq` puts the real bases back, but ONLY when it is given a reference FASTA.
`_cdna_region_task` never passed one to `write_stage1_fastq`, and that parameter used to default to
`None` — so the omission raised nothing and the region-parallel path emitted consensus records still
full of "=".  Measured on a real PCB114 library (Chanfreau unit 907): 94.9 % of records >50 % "=",
realigning at 0.8 % mapped against 99.1 % for the raw reads, with every stage exiting 0.  The
sequential path passed the reference all along, which is why `--region` runs looked clean.

Layered guards, because no single one of them would have caught it:

  * `test_reference_is_required` — the SIGNATURE.  A caller that forgets `reference` must fail
    loudly.  The silent `= None` default is what let this ship.
  * `test_region_task_passes_its_reference_to_the_writer` — the CALL SITE, so re-adding a default
    elsewhere cannot quietly reintroduce the omission.
  * `test_restore_eq_seq_*` / `test_written_fastq_*` — the BEHAVIOUR, on a calmd-style BAM.
  * `test_align_*` — the second line of defence: `align` refuses a calmd FASTQ from ANY producer,
    so the next path that loses a reference fails loudly instead of mapping 0.8 %.

The first two FAIL on the unfixed tree and pass on the fix (verified 2026-09-08); a guard that
cannot fail is not a guard.

Author: rectify-ce, from the diagnosis in GitHub #4 (reported by the Chanfreau session, whose
workers/absolute-path matrix isolated it to the non-sequential path).
"""
import gzip
import inspect
import random
from pathlib import Path

import pysam
import pytest

from rectify.core.cdna.io import write_stage1_fastq

CHROM = "chrT"
MRNA_START, MRNA_LEN = 1000, 600


def _genome(n=4000):
    rng = random.Random(7)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _calmd_bam(tmp: Path):
    """A mapped read whose SEQ carries calmd '=' for every matching base, plus its reference.

    This is exactly what `rectify align`'s `samtools calmd -e` leaves in the BAM: the aligned block
    is '=' where it matched, and the real base only where it did not.  The soft-clipped flanks keep
    their literal sequence, which is why a naive "does the record have any real bases" check does
    not catch the defect.
    """
    genome = _genome()
    g = genome[MRNA_START:MRNA_START + MRNA_LEN]
    fa = tmp / "ref.fa"
    fa.write_text(f">{CHROM}\n" + "\n".join(genome[i:i + 60] for i in range(0, len(genome), 60)) + "\n")
    pysam.faidx(str(fa))

    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": len(genome)}]})
    a = pysam.AlignedSegment(hdr)
    a.query_name = "calmd_read"
    a.reference_id = 0
    a.reference_start = MRNA_START
    a.mapping_quality = 60
    a.flag = 0
    a.query_sequence = "=" * MRNA_LEN          # every aligned base matched -> all '='
    a.cigartuples = [(0, MRNA_LEN)]
    a.query_qualities = pysam.qualitystring_to_array("?" * MRNA_LEN)
    bam = tmp / "in.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=hdr) as out:
        out.write(a)
    pysam.index(str(bam))
    return bam, fa, g


def test_reference_is_required():
    """`reference` must be a required parameter — no default that can be silently skipped."""
    params = inspect.signature(write_stage1_fastq).parameters
    assert "reference" in params, "write_stage1_fastq must take a reference"
    assert params["reference"].default is inspect.Parameter.empty, (
        "reference must NOT have a default: GitHub #4 shipped because `= None` let "
        "_cdna_region_task omit it without raising, emitting unmappable consensus reads."
    )


def test_region_task_passes_its_reference_to_the_writer():
    """The region task's writer call must forward `reference` (the one-line defect in #4)."""
    from rectify.core.commands import cdna_correct_command as C
    src = inspect.getsource(C._cdna_region_task)
    call = src[src.index("_write("):]
    call = call[:call.index("\n    )") + 6]
    assert "reference=" in call, (
        "_cdna_region_task's write_stage1_fastq call must pass reference=; without it the "
        "region-parallel path emits calmd '=' placeholders (GitHub #4)."
    )


def test_restore_eq_seq_resolves_placeholders_with_a_reference(tmp_path):
    """The restorer itself: '=' in, real reference bases out; and a no-op without a reference."""
    from rectify.core.cdna.consensus import restore_eq_seq

    bam, fa, g = _calmd_bam(tmp_path)
    with pysam.AlignmentFile(str(bam)) as fh:
        seg = next(iter(fh))

    with pysam.FastaFile(str(fa)) as ref:
        restored = restore_eq_seq(seg, ref)
    assert "=" not in restored, "with a reference, no placeholder may survive"
    assert restored == g, "restored bases must equal the reference span"

    # And the failure mode the bug rode in on: no reference -> the '=' string passes through.
    assert "=" in restore_eq_seq(seg, None), (
        "without a reference restore_eq_seq is a no-op — which is why the missing kwarg was "
        "invisible; the guards above are what make that unreachable by accident."
    )


def test_written_fastq_has_no_eq_placeholders(tmp_path):
    """End to end through the real writer: a calmd BAM must not yield '=' in the emitted SEQ."""
    from rectify.core.cdna.io import stream_reads

    bam, fa, g = _calmd_bam(tmp_path)
    reads, _ = stream_reads(bam, None, reference=fa)
    if not reads:
        pytest.skip("synthetic molecule is not UMI-extractable in this build")
    ri = reads[0]
    fq = tmp_path / "stage1.fastq.gz"
    write_stage1_fastq(
        bam, fq, [[ri]],
        umi_canonical={0: ri.umi}, cluster_xf_tier={0: ri.xf_tier},
        cluster_tail_len={0: ri.tail_len}, reference=fa,
    )
    with gzip.open(fq, "rt") as fh:
        lines = fh.read().splitlines()
    seqs = [lines[i] for i in range(1, len(lines), 4)]
    # Count FIRST. "no record contains '=' " is vacuously true of an empty file, so a writer that
    # silently emitted nothing would pass the placeholder check on its own (caught in review by the
    # Chanfreau session). One cluster in -> exactly one record out.
    assert len(seqs) == 1, f"one cluster must yield exactly one record, got {len(seqs)}"
    assert all(seqs), "no emitted record may have an empty SEQ"
    for s in seqs:
        assert "=" not in s, f"emitted SEQ still carries calmd placeholders: {s[:60]}"


# --------------------------------------------------------------------------------------
# The second line of defence: `rectify align` refuses a calmd FASTQ from ANY producer.
# --------------------------------------------------------------------------------------

def _fastq(tmp_path, seqs, name="in.fastq"):
    p = tmp_path / name
    with open(p, "w") as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@r{i}\n{s}\n+\n{'?' * len(s)}\n")
    return p


def test_align_refuses_a_calmd_placeholder_fastq(tmp_path):
    from rectify.core.commands.align_command import _refuse_calmd_placeholder_fastq
    bad = _fastq(tmp_path, ["=" * 100] * 50, "bad.fastq")
    assert _refuse_calmd_placeholder_fastq(bad) == 1, "an all-'=' FASTQ must be refused"


def test_align_accepts_a_normal_fastq(tmp_path):
    from rectify.core.commands.align_command import _refuse_calmd_placeholder_fastq
    good = _fastq(tmp_path, ["ACGT" * 25] * 50, "good.fastq")
    assert _refuse_calmd_placeholder_fastq(good) == 0, "a literal-sequence FASTQ must pass"


def test_align_tolerates_a_stray_equals_sign(tmp_path):
    """A record with a few '=' is not the failure mode — only majority-'=' records are.

    The guard must not fire on incidental characters, or it becomes a thing people work around.
    """
    from rectify.core.commands.align_command import _refuse_calmd_placeholder_fastq
    ok = _fastq(tmp_path, ["ACGT" * 24 + "=ACG"] * 50, "stray.fastq")
    assert _refuse_calmd_placeholder_fastq(ok) == 0


def test_align_guard_survives_an_unreadable_file(tmp_path):
    """An unreadable/binary input is the aligner's error to report, not this pre-check's."""
    from rectify.core.commands.align_command import _refuse_calmd_placeholder_fastq
    assert _refuse_calmd_placeholder_fastq(tmp_path / "does_not_exist.fastq") == 0
