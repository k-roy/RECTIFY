"""dorado's per-read poly(A) estimate (`pt:i`) rides Path A into the consensus as XP/XD.

rbrowse request (Kevin, 2026-09-12: "I thought Dorado puts out pt tags for cDNA" — it does,
and the pipeline dropped them). The raw uBAM carries `pt:i` on every read; `correct-cdna`
never read it, and the consensus carried only `XA`, the SEQUENCE-level A-count, which is
~8 nt shorter at the median than dorado's signal estimate and is stripped off the emitted
molecule anyway. Now `ReadInfo.pt` carries the tag, the consensus gets `XP` (mean over
member reads with pt > 0), `XD` (their count) and `XW` (their SD), and both survive the sibling-restore and
CMA whitelists.
"""
import gzip
import random
from pathlib import Path

import pysam

from rectify.core.cdna._constants import ANCHOR_FWD, SSP_FWD, UMI_LEN
from rectify.core.cdna.io import cluster_pt_summary, stream_reads, write_stage1_fastq
from rectify.core.cdna.read_info import ReadInfo, dorado_pt, extract_read_info

CHROM = "chrI"
MRNA_START, MRNA_LEN = 1000, 300


def _genome():
    rng = random.Random(7)
    return "".join(rng.choice("ACGT") for _ in range(4000))


def _bam(tmp: Path, pts):
    """One cluster: n identical Type-1 molecules, one `pt` value each (None = no tag)."""
    genome = _genome()
    rng = random.Random(13)
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
        for i, pt in enumerate(pts):
            a = pysam.AlignedSegment(hdr)
            a.query_name = f"r_{i}"
            a.reference_id = 0
            a.reference_start = MRNA_START
            a.mapping_quality = 60
            a.flag = 0
            a.query_sequence = five + mrna + three
            a.cigartuples = [(4, len(five)), (0, MRNA_LEN), (4, len(three))]
            a.query_qualities = pysam.qualitystring_to_array("?" * len(a.query_sequence))
            if pt is not None:
                a.set_tag("pt", pt, value_type="i")
            out.write(a)
    pysam.index(str(bam))
    return bam, fa


def _stage1_tags(tmp: Path, pts):
    bam, fa = _bam(tmp, pts)
    reads, _ = stream_reads(bam, None, reference=fa)
    assert len(reads) == len(pts)
    fq = tmp / "stage1.fastq.gz"
    write_stage1_fastq(bam, fq, [reads],
                       umi_canonical={0: reads[0].umi},
                       cluster_xf_tier={0: reads[0].xf_tier},
                       cluster_tail_len={0: reads[0].tail_len}, reference=fa)
    with gzip.open(fq, "rt") as fh:
        header = fh.readline().rstrip("\n")
    return reads, {t.split(":")[0]: t.split(":", 2)[2] for t in header.split("\t")[1:]}


class TestReadInfoCarriesPt:
    def test_present_absent_and_negative(self, tmp_path):
        bam, fa = _bam(tmp_path, [40, None, -1, 0])
        reads, _ = stream_reads(bam, None, reference=fa)
        assert [r.pt for r in reads] == [40, None, None, 0]

    def test_dorado_pt_helper(self):
        hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"},
                                               "SQ": [{"SN": CHROM, "LN": 10}]})
        r = pysam.AlignedSegment(hdr)
        assert dorado_pt(r) is None
        r.set_tag("pt", 33, value_type="i")
        assert dorado_pt(r) == 33
        r.set_tag("pt", -1, value_type="i")
        assert dorado_pt(r) is None
        r.set_tag("pt", "x", value_type="Z")
        assert dorado_pt(r) is None

    def test_pt_defaults_to_none_for_hand_built_readinfo(self):
        ri = ReadInfo("r", CHROM, 1, "fwd", "A" * UMI_LEN, False, 2, 12, 0, 1, 1, 0,
                      "umi_captured_fwd")
        assert ri.pt is None


class TestClusterSummary:
    def _r(self, pt):
        return ReadInfo("r", CHROM, 1, "fwd", "", False, 0, 0, 0, 1, 1, 0, "x", pt=pt)

    def test_mean_and_sd_over_positive_values_only(self):
        mean, sd, n = cluster_pt_summary([self._r(40), self._r(50), self._r(-1),
                                          self._r(None), self._r(0)])
        assert (mean, n) == (45.0, 2) and abs(sd - 7.0711) < 1e-3
        mean, sd, n = cluster_pt_summary([self._r(41), self._r(44)])
        assert (mean, n) == (42.5, 2) and abs(sd - 2.1213) < 1e-3
        mean, sd, n = cluster_pt_summary([self._r(118), self._r(20), self._r(35)])
        assert abs(mean - 57.6667) < 1e-3 and n == 3 and abs(sd - 52.7857) < 1e-3

    def test_singleton_has_no_sd(self):
        assert cluster_pt_summary([self._r(118)]) == (118.0, None, 1)

    def test_no_pt_anywhere(self):
        assert cluster_pt_summary([self._r(None), self._r(0)]) == (None, None, 0)

    def test_run_level_correspondence(self):
        from rectify.core.cdna.io import pt_correspondence_summary
        cp = {0: (45.0, 7.07, 2), 1: (50.0, 1.0, 3), 2: (30.0, None, 1), 3: (None, None, 0)}
        c = pt_correspondence_summary(cp)
        assert c['pt_multi_clusters'] == 2
        assert c['pt_sd_median'] in (1.0, 7.07) and c['pt_sd_p90'] == 7.07
        assert 0 < c['pt_cv_median'] < 0.2
        assert pt_correspondence_summary({}) == dict(pt_multi_clusters=0, pt_sd_median=None,
                                                     pt_sd_p90=None, pt_cv_median=None)


class TestConsensusTags:
    def test_xp_and_xd_on_a_cluster_with_pt(self, tmp_path):
        reads, tags = _stage1_tags(tmp_path, [40, 50, None])
        assert tags["XD"] == "2"
        assert tags["XP"] == "45.0"
        assert tags["XW"] == "7.1"
        assert tags["XC"] == "3"
        # XA is the sequence-level count and is untouched by pt
        assert tags["XA"] == str(reads[0].tail_len)

    def test_xd_zero_and_no_xp_when_pt_never_reached_stage1(self, tmp_path):
        _, tags = _stage1_tags(tmp_path, [None, None])
        assert tags["XD"] == "0"
        assert "XP" not in tags and "XW" not in tags

    def test_singleton_with_pt(self, tmp_path):
        _, tags = _stage1_tags(tmp_path, [118])
        assert tags["XD"] == "1" and tags["XP"] == "118.0"
        assert "XW" not in tags                       # no spread from one read


class TestTagsSurviveDownstream:
    def test_in_the_sibling_restore_list_and_the_cma_whitelist(self):
        from rectify.core.consensus.consensus import _CDNA_COMMENT_TAGS
        from rectify.core.multialign.cma_schema import READ_INTRINSIC_TAGS
        for t in ("XP", "XD", "XW"):
            assert t in _CDNA_COMMENT_TAGS
            assert t in READ_INTRINSIC_TAGS

    def test_no_other_writer_uses_xp_or_xd(self):
        import re
        import rectify
        root = Path(rectify.__file__).parent
        pat = re.compile(r"set_tag\(\s*['\"](XP|XD|XW)['\"]")
        hits = [str(py.relative_to(root)) for py in root.rglob("*.py")
                if pat.search(py.read_text())]
        assert hits == [], hits


class TestPretrimKeepsIncomingTags:
    """`run-all --ONT-cDNA` Step 0 rewrites the FASTQ header as bare-UUID + ro + pl; it
    used to drop everything else, so `samtools fastq -T pt` never reached the
    pre-alignment. Well-formed SAM tags are kept, free text is not, ro/pl never doubled."""

    def _trim(self, tmp_path, header):
        from rectify.core.commands.cdna_trim_command import trim_cdna_fastq_polya
        body = "GTTATGTCCTGTCTTTGGTTCAGTTATTGAACCAATGTCACAGGCCTTCCTCCGTGACAAGAAAGTTGTCGGTGTC"
        seq = body + "A" * 25 + "GAAGATAGAGCGACAGGCAAGT"
        fq = tmp_path / "in.fastq.gz"
        with gzip.open(fq, "wt") as fh:
            fh.write(f"@{header}\n{seq}\n+\n{'I' * len(seq)}\n")
        out = tmp_path / "out.fastq.gz"
        trim_cdna_fastq_polya(str(fq), str(out), str(tmp_path / "meta.tsv"),
                              trim_5p_polyt=True)
        with gzip.open(out, "rt") as fh:
            return fh.readline().rstrip("\n")

    def test_pt_and_other_sam_tags_are_carried(self, tmp_path):
        h = self._trim(tmp_path, "uuid-1\tpt:i:40\tqs:f:12.5\tRG:Z:runA")
        fields = h[1:].split("\t")
        assert fields[0] == "uuid-1"
        assert fields[1].startswith("ro:A:") and fields[2].startswith("pl:i:")
        assert fields[3:] == ["pt:i:40", "qs:f:12.5", "RG:Z:runA"]

    def test_free_text_and_our_own_tags_are_not_duplicated(self, tmp_path):
        h = self._trim(tmp_path, "uuid-2 runid=abc ch=7 ro:A:U pl:i:99 pt:i:31")
        fields = h[1:].split("\t")
        assert fields[0] == "uuid-2"
        assert sum(f.startswith("ro:") for f in fields) == 1
        assert sum(f.startswith("pl:") for f in fields) == 1
        assert "pl:i:99" not in fields          # the incoming pl is replaced, not kept
        assert fields[3:] == ["pt:i:31"]

    def test_bare_header_is_unchanged(self, tmp_path):
        h = self._trim(tmp_path, "uuid-3")
        assert len(h[1:].split("\t")) == 3


class TestSilentZeroIsLoud:
    def test_stats_and_warning_when_no_input_read_carries_pt(self, tmp_path, caplog):
        import logging
        bam, fa = _bam(tmp_path, [None, None])
        reads, _ = stream_reads(bam, None, reference=fa)
        with caplog.at_level(logging.WARNING, logger="rectify.core.cdna.io"):
            stats = write_stage1_fastq(bam, tmp_path / "s.fastq.gz", [reads],
                                       umi_canonical={0: reads[0].umi},
                                       cluster_xf_tier={0: reads[0].xf_tier},
                                       cluster_tail_len={0: reads[0].tail_len}, reference=fa)
        assert stats["pt_reads"] == 0 and stats["pt_clusters"] == 0
        assert stats["pt_multi_clusters"] == 0 and stats["pt_sd_median"] is None
        msgs = [r.getMessage() for r in caplog.records]
        assert any("NONE of the 2 clustered reads carries dorado's pt:i tag" in m for m in msgs)
        assert any("samtools fastq -T pt" in m for m in msgs)

    def test_no_warning_when_pt_is_present(self, tmp_path, caplog):
        import logging
        bam, fa = _bam(tmp_path, [40, None])
        reads, _ = stream_reads(bam, None, reference=fa)
        with caplog.at_level(logging.WARNING, logger="rectify.core.cdna.io"):
            stats = write_stage1_fastq(bam, tmp_path / "s.fastq.gz", [reads],
                                       umi_canonical={0: reads[0].umi},
                                       cluster_xf_tier={0: reads[0].xf_tier},
                                       cluster_tail_len={0: reads[0].tail_len}, reference=fa)
        assert stats["pt_reads"] == 1 and stats["pt_clusters"] == 1
        assert not [r for r in caplog.records if "pt:i tag" in r.getMessage()]

    def test_bam_input_with_real_pt_tags_is_the_covered_path(self, tmp_path):
        """The reader is `read.get_tag('pt')` on the input BAM record — a Stage-A BAM
        aligned from `samtools fastq -T pt | minimap2 -y` (Chanfreau 907) qualifies
        directly; no FASTQ comment is involved at correct-cdna time."""
        _, tags = _stage1_tags(tmp_path, [41, 42])
        assert tags["XD"] == "2" and tags["XP"] == "41.5" and tags["XW"] == "0.7"
