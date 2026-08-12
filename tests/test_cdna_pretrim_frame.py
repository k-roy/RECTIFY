#!/usr/bin/env python3
"""Stage-1 consensus adapter pretrim must work on BOTH source strands (planning/681).

The defect: ``orient`` is computed by ``extract_read_info`` on ``read.query_sequence``,
which pysam returns in BAM (reference-forward) frame — the ``ReadInfo.orient`` docstring
says so explicitly ("'fwd' = SSP at LEFT of BAM SEQ"). But ``write_stage1_fastq``'s
singleton and ``rep_fallback`` branches reverse-complement the sequence into *basecalled*
frame, and nothing flipped the label with it. So every minus-strand singleton reached
``pretrim_consensus`` labelled ``fwd`` while carrying ``SSP_RC`` at the far end, the exact
23-mer search returned -1, and the trim silently no-op'd at BOTH ends.

Measured cost before the fix (planning/679 CP4-CP10, /681 CP2): 52.0% of molecules had
``XQ:i:0``; the singleton no-ops carried a mean max soft-clip of **145.2 nt** (vs 4.3 nt for
correctly-trimmed ones) of adapter + random 27-nt UMI + poly-A. Because a random UMI is
maximally high-complexity, no complexity gate rejects it, so the junction resolver
enumerated hundreds of GT/AG candidates against sequence that can never align — the
"1-4 reads/s cDNA" pathology.

These tests drive the REAL path (``stream_reads`` -> ``write_stage1_fastq``), not
``pretrim_consensus`` in isolation, because the defect lived in the seam between them.
"""
from __future__ import annotations

import gzip
import random
from pathlib import Path

import pysam
import pytest

from rectify.core.cdna._constants import ANCHOR_FWD, BRIDGE_LEN, SSP_FWD, SSP_RC, UMI_LEN
from rectify.core.cdna.io import stream_reads, write_stage1_fastq

CHROM = "chrR"
MRNA_START, MRNA_LEN = 1000, 600
FIVE_PRIME_LEN = len(SSP_FWD) + UMI_LEN + BRIDGE_LEN  # 23 + 27 + 3 = 53


def _genome() -> str:
    rng = random.Random(7)
    return "".join(rng.choice("ACGT") for _ in range(4000))


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _build(tmp: Path, is_reverse: bool, use_eq: bool):
    """One PCB114 molecule: [SSP+UMI+GGG]softclip [mRNA]aligned [polyA+adapter]softclip.

    ``use_eq`` writes the aligned block as calmd '=' placeholders, which is what
    ``rectify align`` actually stores and what ``restore_eq_seq`` resolves — one of the
    three candidate explanations in the 679 handoff was that its contract failed on that
    path, so both are exercised.
    """
    genome = _genome()
    rng = random.Random(13)
    umi = "".join(rng.choice("ACGT") for _ in range(UMI_LEN))
    mrna = genome[MRNA_START:MRNA_START + MRNA_LEN]

    fa = tmp / "ref.fa"
    fa.write_text(f">{CHROM}\n" + "\n".join(genome[i:i + 60]
                                            for i in range(0, len(genome), 60)) + "\n")
    pysam.faidx(str(fa))

    five = SSP_FWD + umi + "GGG"
    three = "A" * 12 + ANCHOR_FWD

    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": len(genome)}]})
    a = pysam.AlignedSegment(hdr)
    a.query_name = "r_test"
    a.reference_id = 0
    a.reference_start = MRNA_START
    a.mapping_quality = 60
    a.flag = 16 if is_reverse else 0
    a.query_sequence = five + (("=" * MRNA_LEN) if use_eq else mrna) + three
    a.cigartuples = [(4, len(five)), (0, MRNA_LEN), (4, len(three))]
    a.query_qualities = pysam.qualitystring_to_array("?" * len(a.query_sequence))

    bam = tmp / "in.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=hdr) as out:
        out.write(a)
    pysam.index(str(bam))
    return bam, fa, len(five) + MRNA_LEN + len(three)


def _run_stage1(tmp: Path, is_reverse: bool, use_eq: bool):
    bam, fa, untrimmed_len = _build(tmp, is_reverse, use_eq)
    reads, _ = stream_reads(bam, None, reference=fa)
    assert reads, "the synthetic molecule must be UMI-extractable"
    ri = reads[0]
    fq = tmp / "stage1.fastq.gz"
    stats = write_stage1_fastq(
        bam, fq, [[ri]],
        umi_canonical={0: ri.umi}, cluster_xf_tier={0: ri.xf_tier},
        cluster_tail_len={0: ri.tail_len}, reference=fa,
    )
    with gzip.open(fq, "rt") as fh:
        header, seq = fh.readline().rstrip("\n"), fh.readline().rstrip("\n")
    tags = {t.split(":")[0]: t.split(":", 2)[2] for t in header.split("\t")[1:]}
    return ri, tags, seq, stats, untrimmed_len


@pytest.mark.parametrize("use_eq", [False, True], ids=["real_bases", "calmd_eq"])
@pytest.mark.parametrize("is_reverse", [False, True], ids=["plus_strand", "minus_strand"])
def test_adapter_is_stripped_on_both_strands(tmp_path, is_reverse, use_eq):
    ri, tags, seq, _stats, untrimmed_len = _run_stage1(tmp_path, is_reverse, use_eq)

    assert int(tags["XQ"]) == FIVE_PRIME_LEN, (
        f"5' trim no-op on is_reverse={is_reverse}: XQ={tags['XQ']} (expected "
        f"{FIVE_PRIME_LEN}). This is the planning/681 frame bug — orient={ri.orient!r} is a "
        f"BAM-frame label and the emitted consensus was in basecalled frame."
    )
    assert int(tags["XK"]) > 0, f"3' trim no-op on is_reverse={is_reverse}"
    assert SSP_FWD not in seq and SSP_RC not in seq, (
        "adapter survived the trim and will become a long terminal soft-clip"
    )
    assert ANCHOR_FWD not in seq, "PCB114 adapter anchor survived the trim"
    assert len(seq) < untrimmed_len, "nothing was trimmed at all"
    # XQ/XK must account for exactly what was removed, in either frame.
    assert int(tags["XQ"]) + int(tags["XK"]) + len(seq) == untrimmed_len


def test_both_strands_yield_the_same_trim(tmp_path):
    """The molecule is the same molecule; which strand it was sequenced on must not
    change how much adapter is stripped. Pre-fix, minus-strand stripped ZERO."""
    out = {}
    for is_reverse in (False, True):
        d = tmp_path / f"rev{int(is_reverse)}"
        d.mkdir()
        _ri, tags, seq, _s, _n = _run_stage1(d, is_reverse, use_eq=False)
        out[is_reverse] = (int(tags["XQ"]), int(tags["XK"]), len(seq))
    assert out[False] == out[True], (
        f"strand-dependent trim: plus={out[False]} minus={out[True]}"
    )


def _build_multiread(tmp: Path, n_reads: int = 2):
    """A multi-read cluster, every member on the MINUS strand.

    Under the fix these take the ``pileup`` branch, which applies no reverse-complement
    and so stays in BAM (reference) frame — where ``orient`` is already defined. The
    minus strand is what makes this test load-bearing: it is exactly the condition that
    triggers the flip on the singleton branch.
    """
    genome = _genome()
    rng = random.Random(29)
    umi = "".join(rng.choice("ACGT") for _ in range(UMI_LEN))
    mrna = genome[MRNA_START:MRNA_START + MRNA_LEN]

    fa = tmp / "ref.fa"
    fa.write_text(f">{CHROM}\n" + "\n".join(genome[i:i + 60]
                                            for i in range(0, len(genome), 60)) + "\n")
    pysam.faidx(str(fa))

    five, three = SSP_FWD + umi + "GGG", "A" * 12 + ANCHOR_FWD
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": len(genome)}]})
    bam = tmp / "in.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=hdr) as out:
        for i in range(n_reads):
            a = pysam.AlignedSegment(hdr)
            a.query_name = f"r_{i}"
            a.reference_id = 0
            a.reference_start = MRNA_START
            a.mapping_quality = 60
            a.flag = 16                      # every read minus-strand
            a.query_sequence = five + mrna + three
            a.cigartuples = [(4, len(five)), (0, MRNA_LEN), (4, len(three))]
            a.query_qualities = pysam.qualitystring_to_array("?" * len(a.query_sequence))
            out.write(a)
    pysam.index(str(bam))
    return bam, fa, mrna


def test_pileup_branch_is_not_frame_flipped(tmp_path):
    """🔴 The invariant that the 679 brief got BACKWARDS — do not "fix" this branch.

    679's Phase 2 proposed adding the "restore basecalled orientation" RC to the pileup
    branch, on the reasoning that its 94.8% `XQ==0` rate meant it was the most broken.
    Measurement said the opposite (planning/681 CP2): `pileup_consensus` builds from
    ``get_aligned_pairs(matches_only=True)``, so it contains no soft-clip and therefore no
    adapter to strip — its molecules carry 13.7 nt mean max clip vs 145.2 nt for the
    genuinely untrimmed singletons, and the A/B realignment showed it bit-identical before
    and after the fix. It is already in BAM frame, which is the frame `orient` is defined
    in. **Adding the RC here would break a correctly-framed branch.**
    """
    bam, fa, mrna = _build_multiread(tmp_path)
    reads, _ = stream_reads(bam, None, reference=fa)
    assert len(reads) >= 2
    fq = tmp_path / "stage1.fastq.gz"
    stats = write_stage1_fastq(
        bam, fq, [list(reads)],
        umi_canonical={0: reads[0].umi}, cluster_xf_tier={0: reads[0].xf_tier},
        cluster_tail_len={0: reads[0].tail_len}, reference=fa,
    )
    with gzip.open(fq, "rt") as fh:
        header, seq = fh.readline().rstrip("\n"), fh.readline().rstrip("\n")

    assert stats["from_multi_pileup"] == 1, "expected the pileup branch"
    assert "XM:Z:pileup" in header
    assert stats["trim_frame_flipped"] == 0, (
        "the pileup branch applies no reverse-complement, so its frame must NOT be "
        "flipped — see planning/681 CP2 before changing this"
    )
    # The decisive assertion: BAM frame, not basecalled frame, despite is_reverse=True.
    assert seq and seq in mrna, (
        "pileup consensus must stay in BAM (reference) frame; it appears reverse-"
        "complemented, which means a frame flip was wrongly applied to this branch"
    )
    assert seq not in _revcomp(mrna)


def test_frame_counters_are_reported(tmp_path):
    """The fail-loud counter: a silent trim no-op is what let this run for months as an
    unexplained throughput pathology, so it must surface as a number."""
    _ri, _tags, _seq, stats, _n = _run_stage1(tmp_path, is_reverse=True, use_eq=False)
    assert stats["trim_frame_flipped"] == 1, (
        "a minus-strand singleton is RC'd into basecalled frame and must be counted"
    )
    assert stats["trim_noop_5p"] == 0 and stats["trim_noop_3p"] == 0
    # After the flip the label and the sequence agree again, so nothing should look
    # mismatched to pretrim_consensus's own sequence-driven frame detection.
    assert stats["trim_frame_mismatch"] == 0
