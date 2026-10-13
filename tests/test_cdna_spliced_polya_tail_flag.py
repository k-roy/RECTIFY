"""A17 / spliced poly(A) tail: cdna-analyze reports a 3' end placed on a genomic A-run after an intron (report-only).

The audit's pre-registered rule (dev/audit_20261012/AUDIT.md section 6): a terminal block of <= 30 aligned nt next to
an N >= 30 nt, whose genome ends in an A-run (>= 6 of its 8 RNA-3'-most bases A, T in the 'rev' frame). RPS26B's
case is the shape: 32 nt of the mRNA's end then TATTGGTTTTC + A11 placed 609 nt downstream; 744 of 1,675 molecules in
the May cDNA tables, 1 of 4,423 in a current library. Xp:i carries the donor-side end; nothing is re-anchored.
"""

import argparse
import random

import pysam

from rectify.core.cdna.spliced_tail import terminal_arun_junction

RPS26B_BLOCK = "TATTGGTTTTC" + "A" * 11


def _genome_with(pos, seq, length=3000, seed=3):
    rng = random.Random(seed)
    g = list("".join(rng.choice("CG") for _ in range(length)))     # no A/T background: no chance A-runs
    g[pos:pos + len(seq)] = seq
    return "".join(g)


def test_rps26b_shape_is_flagged_with_the_donor_end():
    genome = _genome_with(1709, RPS26B_BLOCK)
    # body 1000..1099 (100M), N 1100..1708 (609N), block 1709..1730 (22M)
    got = terminal_arun_junction([(0, 100), (3, 609), (0, 22)], 1000, genome, rna_reverse=False)
    assert got == (1099, 22, 609)


def test_minus_strand_mirror_is_flagged():
    block_rc = ("T" * 11) + "GAAAACCAATA"            # reverse complement of RPS26B_BLOCK
    genome = _genome_with(1000, block_rc)
    # block 1000..1021 (22M), N 1022..1630, body 1631..1730
    got = terminal_arun_junction([(0, 22), (3, 609), (0, 100)], 1000, genome, rna_reverse=True)
    assert got == (1631, 22, 609)


def test_a_long_last_exon_is_not_flagged():
    genome = _genome_with(1709, "C" * 20 + RPS26B_BLOCK)
    assert terminal_arun_junction([(0, 100), (3, 609), (0, 42)], 1000, genome, rna_reverse=False) is None


def test_a_short_gap_is_not_flagged():
    genome = _genome_with(1120, RPS26B_BLOCK)
    assert terminal_arun_junction([(0, 100), (3, 20), (0, 22)], 1000, genome, rna_reverse=False) is None


def test_a_block_not_on_an_a_run_is_not_flagged():
    genome = _genome_with(1709, "CGCGCGCGCGCGCGCGCGCGCG")
    assert terminal_arun_junction([(0, 100), (3, 609), (0, 22)], 1000, genome, rna_reverse=False) is None


def test_unspliced_records_are_not_flagged():
    genome = _genome_with(1078, RPS26B_BLOCK)
    assert terminal_arun_junction([(0, 100)], 1000, genome, rna_reverse=False) is None


CHROM = "chrP"


def _mol(name, start, cigar, reverse, seq):
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": 3000}]})
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.flag = 16 if reverse else 0
    r.cigartuples = cigar
    r.query_sequence = seq
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    xo = "rev" if reverse else "fwd"
    for tag, val, typ in [("XU", "A" * 27, "Z"), ("XO", xo, "Z"), ("XT", 1, "i"), ("XY", f"umi_captured_{xo}", "Z"),
                          ("XC", 2, "i"), ("XF", 2, "i"), ("XA", 30, "i"), ("XK", 40, "i"), ("XN", 1, "i")]:
        r.set_tag(tag, val, typ)
    return r


def test_cdna_analyze_tags_and_counts_the_spliced_tail(tmp_path, caplog):
    from rectify.core.commands import cdna_analyze_command
    genome = _genome_with(1709, RPS26B_BLOCK)
    fa = tmp_path / "ref.fa"
    fa.write_text(f">{CHROM}\n{genome}\n")
    pysam.faidx(str(fa))
    gff = tmp_path / "ref.gff"
    gff.write_text("##gff-version 3\n"
                   f"{CHROM}\ts\tgene\t901\t1200\t.\t+\t.\tID=g1;Name=g1\n"
                   f"{CHROM}\ts\tgene\t2001\t2400\t.\t+\t.\tID=g2;Name=g2\n")
    spliced = _mol("spliced", 1000, [(0, 100), (3, 609), (0, 22)], False, genome[1000:1100] + RPS26B_BLOCK)
    plain = _mol("plain", 2100, [(0, 150)], False, genome[2100:2250])
    bam = tmp_path / "in.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=spliced.header) as out:
        for r in (spliced, plain):
            out.write(r)
    import logging
    with caplog.at_level(logging.INFO, logger="cdna-analyze"):
        rc = cdna_analyze_command.run(argparse.Namespace(bam=bam, out=tmp_path / "out", gff=gff, reference=fa))
    assert rc == 0
    with pysam.AlignmentFile(str(tmp_path / "out" / "consensus_tagged.bam")) as tagged:
        xp = {r.query_name: (r.get_tag("Xp") if r.has_tag("Xp") else None) for r in tagged.fetch(until_eof=True)}
    assert xp == {"spliced": 1099, "plain": None}
    assert any("1 of 2 molecules end in a block of <= 30 nt on a genomic A-run" in m for m in caplog.messages), \
        caplog.messages
