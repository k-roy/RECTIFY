"""658: depth-adaptive UMI edit distance + per-library error-profile calibration."""
from __future__ import annotations

import pysam
import pytest

from rectify.core.cdna._constants import SSP_FWD, SSP_RC, UMI_LEN
from rectify.core.cdna.cluster import cluster_reads
from rectify.core.cdna.profile import derive_umi_error_profile
from rectify.core.cdna.read_info import ReadInfo, revcomp

STRUCTURED = "TT" + "ACGA" + "TT" + "CAGC" + "TT" + "GCAC" + "TT" + "AGCA" + "TTT"
assert len(STRUCTURED) == UMI_LEN


def _ri(read_id, anchor, umi, read_type=1):
    return ReadInfo(read_id=read_id, chrom="chrI", anchor=anchor, orient="fwd",
                    umi=umi, is_reverse=False, xf_tier=2, tail_len=10,
                    aln_start=anchor - 500, aln_end=anchor, read_type=read_type,
                    pos5_corrected=anchor - 500, read_subtype="umi_captured_fwd")


def _bucket(anchor, n_hub, with_partner, n_filler, tag):
    hub = "A" * UMI_LEN
    partner = "A" * (UMI_LEN - 2) + "CC"       # Lev 2 from hub
    filler = "G" * UMI_LEN                     # far from everything
    reads = [_ri(f"{tag}_hub{i}", anchor, hub) for i in range(n_hub)]
    if with_partner:
        reads.append(_ri(f"{tag}_partner", anchor, partner))
    reads += [_ri(f"{tag}_fill{i}", anchor, filler) for i in range(n_filler)]
    return reads


class TestAdaptiveEd:
    def test_deep_bucket_drops_to_ed1_shallow_keeps_ed2(self):
        deep = _bucket(1000, n_hub=4, with_partner=True, n_filler=25, tag="d")   # 30 reads
        shallow = _bucket(9000, n_hub=4, with_partner=True, n_filler=0, tag="s")  # 5 reads
        clusters, stats = cluster_reads(deep + shallow, anchor_window=1, max_edit=2,
                                        per_cluster_cap=200, adaptive_threshold=20)
        by_anchor = {}
        for c in clusters:
            by_anchor.setdefault(c[0].anchor, []).append(sorted(r.read_id for r in c))
        # deep bucket (>= threshold) clustered at ed=1: the Lev-2 partner stays separate
        assert len(by_anchor[1000]) == 3
        # shallow bucket clustered at ed=2: partner absorbed by the hub
        assert len(by_anchor[9000]) == 1
        assert stats["adaptive_deep_buckets"] == 1
        assert stats["adaptive_deep_reads"] == 30
        assert stats["adaptive_threshold"] == 20

    def test_threshold_zero_is_prior_behaviour(self):
        deep = _bucket(1000, n_hub=4, with_partner=True, n_filler=25, tag="d")
        clusters, stats = cluster_reads(deep, anchor_window=1, max_edit=2,
                                        per_cluster_cap=200, adaptive_threshold=0)
        assert len(clusters) == 2  # partner merges into the hub at ed=2
        assert stats["adaptive_deep_buckets"] == 0

    def test_type2_buckets_unaffected(self):
        t2 = [_ri(f"t2_{i}", 5000, "", read_type=2) for i in range(50)]
        clusters, stats = cluster_reads(t2, anchor_window=1, max_edit=2,
                                        per_cluster_cap=200, adaptive_threshold=10)
        assert stats["adaptive_deep_buckets"] == 0  # position clustering, not UMI


class TestProfileDerivation:
    def test_synthetic_bam_profiles_both_regimes(self, tmp_path):
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": "chrI", "LN": 100000}],
        })
        fwd_seq = SSP_FWD + STRUCTURED + "GGG" + "C" * 60
        rev_seq = "C" * 40 + revcomp(STRUCTURED) + SSP_RC + "C" * 20
        bam_path = tmp_path / "toy.bam"
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            for name, seq, flag in [("r_fwd", fwd_seq, 0), ("r_rev", rev_seq, 16)]:
                seg = pysam.AlignedSegment(header)
                seg.query_name = name
                seg.query_sequence = seq
                seg.flag = flag
                seg.reference_id = 0
                seg.reference_start = 100
                seg.mapping_quality = 60
                seg.cigartuples = [(0, len(seq))]
                seg.query_qualities = pysam.qualitystring_to_array("?" * len(seq))
                bam.write(seg)
        prof = derive_umi_error_profile(str(bam_path))
        assert prof["n_umi"] == 2
        # fwd/flag0 -> entry; rev(SSP_RC)/flag16 -> UMI at basecalled start -> entry too
        assert set(prof["cells"]) == {"fwd_rev0", "rev_rev1"}
        for cell in prof["cells"].values():
            assert cell["regime"] == "entry"
            assert cell["clean"] == 1.0
        assert prof["regimes"]["entry"]["clean"] == 1.0
        assert prof["regimes"]["entry"]["inframe_sub_per_base"] == 0.0
