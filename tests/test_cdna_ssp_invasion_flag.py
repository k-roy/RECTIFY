"""ISSUE-086: cdna-analyze reports molecules whose aligned 5' end sits on a genomic copy of the strand-switching primer.

CDC48 (YDL126C): 233 Ysh1-AA molecules start at chrIV:237,413, where the genome reads GGTGCTGATATTGC in RNA sense,
SSP positions 9-22; direct RNA has no pile there (Chanfreau planning/1036, card C09). The flag is report-only: an
`Xs:i` tag with the match length on the tagged BAM, a log line and a count in the stage sidecar. Nothing is dropped.
"""

import argparse
import random

import pysam

from rectify.core.cdna._constants import SSP_FWD
from rectify.core.cdna.ssp_invasion import SSP_FLAG_MIN_NT, chance_rate, rna_5p_genome, ssp_match_length

CDC48_SITE = "GGTGCTGATATTGC"          # SSP_FWD[8:22]
_RC = str.maketrans("ACGT", "TGCA")


def _rc(s):
    return s.translate(_RC)[::-1]


def test_cdc48_site_matches_fourteen_primer_bases():
    assert SSP_FWD[8:22] == CDC48_SITE
    assert ssp_match_length(CDC48_SITE + "AAACCC") == 14


def test_random_genome_rarely_matches_far():
    rng = random.Random(1)
    lens = [ssp_match_length("".join(rng.choice("ACGT") for _ in range(23))) for _ in range(5000)]
    assert max(lens) < SSP_FLAG_MIN_NT
    assert abs(chance_rate(12) - 12 / 4 ** 12) < 1e-12


def test_window_is_rna_sense_on_both_strands():
    genome = "C" * 50 + CDC48_SITE + "G" * 50
    plus = rna_5p_genome(genome, 50, 120, rna_reverse=False)
    assert plus.startswith(CDC48_SITE)
    genome_m = "C" * 50 + _rc(CDC48_SITE) + "G" * 50          # minus-strand RNA ending at 50 + 14
    minus = rna_5p_genome(genome_m, 10, 64, rna_reverse=True)
    assert minus.startswith(CDC48_SITE)
    assert rna_5p_genome(genome, 990, 1000, rna_reverse=False) is None


CHROM, GLEN = "chrS", 3000


def _genome():
    rng = random.Random(7)
    g = list("".join(rng.choice("ACGT") for _ in range(GLEN)))
    g[1000:1014] = CDC48_SITE                       # + strand molecule starts here
    g[1014] = "A"                                    # ...and the match stops at 14 (SSP_FWD[22] is T)
    g[2086:2100] = _rc(CDC48_SITE)                   # - strand molecule's 5' end is base 2099
    g[2085] = "C"                                    # RNA-sense G after the site: the match stops at 14
    return "".join(g)


def _mol(name, start, length, reverse, genome):
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": GLEN}]})
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.flag = 16 if reverse else 0
    body = genome[start:start + length]
    r.query_sequence = body
    r.cigartuples = [(0, length)]
    r.query_qualities = pysam.qualitystring_to_array("I" * length)
    xo = "rev" if reverse else "fwd"
    for tag, val, typ in [("XU", "A" * 27, "Z"), ("XO", xo, "Z"), ("XT", 1, "i"), ("XY", f"umi_captured_{xo}", "Z"),
                          ("XC", 2, "i"), ("XF", 2, "i"), ("XA", 30, "i"), ("XK", 40, "i"), ("XN", 1, "i")]:
        r.set_tag(tag, val, typ)
    return r


def test_cdna_analyze_tags_and_counts_ssp_starts(tmp_path, caplog):
    from rectify.core.commands import cdna_analyze_command
    genome = _genome()
    fa = tmp_path / "ref.fa"
    fa.write_text(f">{CHROM}\n{genome}\n")
    pysam.faidx(str(fa))
    gff = tmp_path / "ref.gff"
    gff.write_text("##gff-version 3\n"
                   f"{CHROM}\ts\tgene\t901\t1400\t.\t+\t.\tID=g1;Name=g1\n"
                   f"{CHROM}\ts\tgene\t1800\t2200\t.\t-\t.\tID=g2;Name=g2\n"
                   f"{CHROM}\ts\tgene\t301\t700\t.\t+\t.\tID=g3;Name=g3\n")
    recs = [_mol("plain", 400, 200, False, genome),
            _mol("plus_ssp", 1000, 200, False, genome),
            _mol("minus_ssp", 1900, 200, True, genome)]       # aligned 1900..2099, 5' end = 2099; coordinate order
    bam = tmp_path / "in.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=recs[0].header) as out:
        for r in recs:
            out.write(r)
    import logging
    with caplog.at_level(logging.INFO, logger="cdna-analyze"):
        rc = cdna_analyze_command.run(argparse.Namespace(bam=bam, out=tmp_path / "out", gff=gff, reference=fa))
    assert rc == 0
    with pysam.AlignmentFile(str(tmp_path / "out" / "consensus_tagged.bam")) as tagged:
        xs = {r.query_name: (r.get_tag("Xs") if r.has_tag("Xs") else None) for r in tagged.fetch(until_eof=True)}
    assert xs["plus_ssp"] == 14 and xs["minus_ssp"] == 14
    assert xs["plain"] is None
    # The count also goes to the stage sidecar, which refuses paths under $TMPDIR, so the test reads the log line.
    assert any(f"2 of 3 molecules start on a genomic match of >= {SSP_FLAG_MIN_NT} nt" in m for m in caplog.messages), \
        caplog.messages
