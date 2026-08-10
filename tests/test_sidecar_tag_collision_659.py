"""659: sidecar tag restore must OVERWRITE aligner-emitted collisions on the cDNA tag set.

uLTRA emits its own XC:Z:NO_SPLICE / XA:Z: (splice semantics). The old don't-clobber guard in
_restore_sidecar_tags let them shadow the authoritative Stage-1 comment tags, and cdna-analyze's
int(XC) then dropped every uLTRA-won record ("missing required tags", 7.0% of rna15aa_rep1,
biased toward spliced molecules).
"""
from __future__ import annotations

import pysam

from rectify.core.consensus.consensus import _restore_sidecar_tags

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6"},
    "SQ": [{"SN": "chrI", "LN": 100000}],
})


class _Row:
    def __init__(self, comment):
        self.fastq_comment = comment


class _Sidecar:
    def __init__(self, comment):
        self._row = _Row(comment)

    def lookup(self, rn):
        return self._row

    def lookup_by_qname(self, name):
        return self._row


def _ultra_winner():
    r = pysam.AlignedSegment(HEADER)
    r.query_name = "cluster_r000_7"
    r.query_sequence = "ACGT" * 25
    r.flag = 0
    r.reference_id = 0
    r.reference_start = 500
    r.mapping_quality = 60
    r.cigartuples = [(0, 100)]
    r.set_tag("RN", 7, value_type="i")
    r.set_tag("XC", "NO_SPLICE", value_type="Z")   # uLTRA's colliding tag
    r.set_tag("XA", "", value_type="Z")            # uLTRA's colliding tag
    r.set_tag("NM", 3, value_type="i")
    return r


COMMENT = ("XU:Z:TTACACTTCAGGTTGCATTTGCCGAGA\tXO:Z:rev\tXC:i:5\tXR:Z:abc\tXM:Z:rep\t"
           "XF:i:2\tXA:i:12\tXT:i:1\tXY:Z:umi_captured_rev\tXQ:i:0\tXK:i:0\tXB:Z:1/0\t"
           "Zx:Z:keepme_comment")


def test_cdna_tags_overwrite_aligner_collisions():
    read = _ultra_winner()
    _restore_sidecar_tags(read, _Sidecar(COMMENT))
    assert read.get_tag("XC") == 5          # int, not "NO_SPLICE"
    assert int(read.get_tag("XC")) == 5     # cdna-analyze's exact access pattern
    assert read.get_tag("XA") == 12
    assert read.get_tag("XU") == "TTACACTTCAGGTTGCATTTGCCGAGA"
    assert read.get_tag("XT") == 1
    assert read.get_tag("NM") == 3          # aligner-computed tag untouched


def test_non_cdna_tags_keep_dont_clobber_guard():
    read = _ultra_winner()
    read.set_tag("Zx", "keepme_read", value_type="Z")
    _restore_sidecar_tags(read, _Sidecar(COMMENT))
    assert read.get_tag("Zx") == "keepme_read"  # not in _CDNA_COMMENT_TAGS -> guard holds


def test_inject_rn_overwrites_ultra_collisions(tmp_path):
    """The PER-ALIGNER injection choke point (multi_aligner._inject_rn_into_bam) is the site
    that actually bit in production: uLTRA BAM records reached consensus already carrying
    XU (restored) + uLTRA's own XC:Z:NO_SPLICE (guarded, unfixed), so the sibling-restore
    saw XU and early-returned. The comment set must overwrite there too."""
    import gzip

    from rectify.core.align.multi_aligner import _inject_rn_into_bam

    bam_path = tmp_path / "ultra.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=HEADER) as bam:
        bam.write(_ultra_winner())
    fq_path = tmp_path / "stage1.fastq.gz"
    with gzip.open(fq_path, "wt") as fq:
        seq = "ACGT" * 25
        fq.write(f"@cluster_r000_7\t{COMMENT}\n{seq}\n+\n{'?' * len(seq)}\n")

    _inject_rn_into_bam(str(bam_path), qname_to_rn={}, reads_path=str(fq_path))

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        rec = next(iter(bam))
    assert int(rec.get_tag("XC")) == 5      # overwrote XC:Z:NO_SPLICE
    assert rec.get_tag("XA") == 12          # overwrote uLTRA's empty XA:Z
    assert rec.get_tag("XU") == "TTACACTTCAGGTTGCATTTGCCGAGA"
    assert rec.get_tag("NM") == 3           # aligner tag untouched
