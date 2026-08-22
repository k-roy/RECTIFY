#!/usr/bin/env python3
"""Stage-1 consensus emission must be RNA-sense on every strand/frame (planning/730, A0).

The defect: `rectify correct-cdna` restored the *basecalled* orientation and merely
labelled the molecule with ``XO``, so the pA-first half of a Path-A library reached
`rectify align` antisense. minimap2 ``-uf`` then scored splice motifs on the wrong
strand and manufactured decoy junctions on 15.7 % of those spliced reads (vs 2.8 % on
the sense half; W7 RPS24B / W9 RPS4A in the 730 walkthrough were both this artefact).

Fix: every emitted molecule is reverse-complemented into RNA-sense when its trim frame
is ``rev`` and tagged ``XN:i:1``. These tests enumerate the four (XO × source strand)
cases explicitly — a real-data smoke can contain a skewed mix and pass by accident —
and drive the REAL path (``stream_reads`` -> ``write_stage1_fastq``).
"""
from __future__ import annotations

import gzip
import random
from pathlib import Path

import pysam
import pytest

from rectify.core.cdna._constants import ANCHOR_FWD, SSP_FWD, SSP_RC, UMI_LEN
from rectify.core.cdna.io import stream_reads, write_stage1_fastq

CHROM = "chrR"
MRNA_START, MRNA_LEN = 1000, 600


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _build(tmp: Path, xo: str, is_reverse: bool):
    """One PCB114 molecule whose mRNA body is ``G`` (genome[1000:1600]).

    ``xo='fwd'``: BAM SEQ = [SSP+UMI+GGG] G [polyA+adapter]     (SSP at the LEFT)
    ``xo='rev'``: BAM SEQ = revcomp of that                         (SSP_RC at the RIGHT)
    ``is_reverse`` only sets the flag (which frame stage 1 restores the sequence in).
    Whatever the combination, the molecule's mRNA-sense sequence is ``G``.
    """
    rng = random.Random(7)
    genome = "".join(rng.choice("ACGT") for _ in range(4000))
    umi = "".join(random.Random(13).choice("ACGT") for _ in range(UMI_LEN))
    g = genome[MRNA_START:MRNA_START + MRNA_LEN]
    five = SSP_FWD + umi + "GGG"
    three = "A" * 12 + ANCHOR_FWD

    fa = tmp / "ref.fa"
    fa.write_text(f">{CHROM}\n" + "\n".join(genome[i:i + 60] for i in range(0, len(genome), 60)) + "\n")
    pysam.faidx(str(fa))

    hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": len(genome)}]})
    a = pysam.AlignedSegment(hdr)
    a.query_name = "r_test"
    a.reference_id = 0
    a.reference_start = MRNA_START
    a.mapping_quality = 60
    a.flag = 16 if is_reverse else 0
    if xo == "fwd":
        a.query_sequence = five + g + three
        a.cigartuples = [(4, len(five)), (0, MRNA_LEN), (4, len(three))]
    else:
        a.query_sequence = _revcomp(five + g + three)
        a.cigartuples = [(4, len(three)), (0, MRNA_LEN), (4, len(five))]
    a.query_qualities = pysam.qualitystring_to_array("?" * len(a.query_sequence))
    bam = tmp / "in.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=hdr) as out:
        out.write(a)
    pysam.index(str(bam))
    return bam, fa, g


def _run_stage1(tmp: Path, xo: str, is_reverse: bool):
    bam, fa, g = _build(tmp, xo, is_reverse)
    reads, _ = stream_reads(bam, None, reference=fa)
    assert reads, "the synthetic molecule must be UMI-extractable"
    ri = reads[0]
    assert ri.orient == xo, f"fixture built for XO={xo} but extract_read_info said {ri.orient}"
    fq = tmp / "stage1.fastq.gz"
    stats = write_stage1_fastq(
        bam, fq, [[ri]],
        umi_canonical={0: ri.umi}, cluster_xf_tier={0: ri.xf_tier},
        cluster_tail_len={0: ri.tail_len}, reference=fa,
    )
    with gzip.open(fq, "rt") as fh:
        header, seq = fh.readline().rstrip("\n"), fh.readline().rstrip("\n")
    tags = {t.split(":")[0]: t.split(":", 2)[2] for t in header.split("\t")[1:]}
    return ri, tags, seq, stats, g


@pytest.mark.parametrize("is_reverse", [False, True], ids=["flag0", "flag16"])
@pytest.mark.parametrize("xo", ["fwd", "rev"])
def test_emitted_consensus_is_rna_sense(tmp_path, xo, is_reverse):
    ri, tags, seq, stats, g = _run_stage1(tmp_path, xo, is_reverse)
    assert tags["XO"] == xo, "XO label must be preserved (it is BAM-SEQ-frame and still exact)"
    assert tags.get("XN") == "1", "every emitted molecule must carry the oriented tag XN:i:1"
    # The poly(A) trim may also eat the genomic A's adjacent to the tail (the walkback
    # cannot tell them apart), so the body can come back a few nt short — but it must be
    # the SENSE strand, as a prefix of g, not the antisense strand.
    assert g.startswith(seq) and len(seq) >= MRNA_LEN - 10, (
        f"XO={xo} flag16={is_reverse}: emitted consensus is not the mRNA-sense body "
        f"(starts {seq[:12]!r}, sense starts {g[:12]!r}, antisense starts {_revcomp(g)[:12]!r}; "
        f"len {len(seq)} vs {MRNA_LEN}). Pre-fix this molecule reached minimap2 -uf antisense "
        "and manufactured decoy junctions."
    )
    assert not _revcomp(g).startswith(seq[:40]), "emitted antisense"
    assert SSP_FWD not in seq and SSP_RC not in seq and ANCHOR_FWD not in seq
    # The rev frame is exactly the branch that had to be flipped: (fwd,flag0) and (rev,flag16)
    # are already sense after stage 1's frame restore; (fwd,flag16) and (rev,flag0) are not.
    expected_flip = (xo == "fwd") == is_reverse
    assert stats["reoriented_to_sense"] == int(expected_flip), stats


def test_all_four_cases_emit_the_identical_sequence(tmp_path):
    """The same molecule sequenced four ways must hand the aligner one sequence."""
    out = {}
    for xo in ("fwd", "rev"):
        for is_reverse in (False, True):
            d = tmp_path / f"{xo}_{int(is_reverse)}"
            d.mkdir()
            out[(xo, is_reverse)] = _run_stage1(d, xo, is_reverse)[2]
    assert len(set(out.values())) == 1, {k: v[:15] for k, v in out.items()}
