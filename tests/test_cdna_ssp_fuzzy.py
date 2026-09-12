"""Type-1 detection tolerates basecall errors in the strand-switching primer.

Regression guard for rbrowse r120_4585 / WT-AA_WT_rep2 (Kevin, 2026-09-11): a full-length
sense molecule whose 23-nt SSP carried ONE miscalled base (the final T read as C) was
classified Type 2 — no UMI, never deduplicated — and kept its whole 143-nt adapter as a 5'
soft clip. In that library 79 of 502 construct reads (16 %) were in this class; 67 of them
carried an intact SSP+UMI with one or two errors. The tier detector and the stage-1 trimmer
were already edlib-tolerant; Type-1 detection was the last exact `find`.

The finder is shared (`read_info.find_ssp_span`) so the classifier and the trimmer can never
disagree about what an SSP is; exact hits keep the historical, unwindowed, byte-identical
path, and the fuzzy fallback is window-gated to the end where the SSP belongs.
"""
import random

import pysam
import pytest

from rectify.core.cdna import consensus as cs
from rectify.core.cdna import read_info as ri
from rectify.core.cdna._constants import BRIDGE_LEN, SSP_FWD, SSP_RC, UMI_LEN
from rectify.core.cdna.read_info import extract_read_info, find_ssp_span, revcomp

# The 143-nt 5' soft clip of r120_4585 as stored (barcode flank, SSP with its last base
# miscalled T->C, UMI, GGG bridge).
R120_4585_CLIP = (
    "GTGTTATGTACATATACTTCGTTCAGTTATGCAGCCAGCCGATATTGGCAAAGCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTG"
    "TTTCTGTTGGTGCTGATATTGCCTTAAGATTCGCATTGCCCTTGGGATTCGGG"
)
R120_4585_UMI = "TTAAGATTCGCATTGCCCTTGGGATTC"  # SSP last base miscalled T->C; the UMI follows it

_ADAPTER_FWD = "GAAGATAGAGCGACAGGCAAGT"


def _body(n=400, seed=1):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


def _read(seq, clip5, clip3, name='r'):
    h = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'},
                                         'SQ': [{'SN': 'chrI', 'LN': 100000}]})
    r = pysam.AlignedSegment(h)
    r.query_name = name
    r.reference_id = 0
    r.reference_start = 5380
    r.mapping_quality = 12
    r.query_sequence = seq
    aligned = len(seq) - clip5 - clip3
    r.cigarstring = f'{clip5}S{aligned}M{clip3}S' if clip3 else f'{clip5}S{aligned}M'
    return r


def _fwd_molecule(ssp, umi, body):
    tail = 'A' * 30 + _ADAPTER_FWD
    return ssp + umi + 'GGG' + body + tail


class TestTheReportedRead:
    def test_r120_4585_is_type1_with_its_umi(self):
        seq = R120_4585_CLIP + _body() + 'A' * 30 + _ADAPTER_FWD
        info = extract_read_info(_read(seq, len(R120_4585_CLIP), 52, 'r120_4585'))
        assert info.read_type == 1
        assert info.read_subtype == 'umi_captured_fwd'
        assert info.orient == 'fwd'
        assert info.umi == R120_4585_UMI

    def test_the_exact_ssp_is_not_in_that_clip(self):
        # documents the cause: one miscalled base, the last T of the SSP
        assert SSP_FWD not in R120_4585_CLIP
        assert SSP_FWD[:-1] in R120_4585_CLIP

    def test_the_trimmer_strips_the_whole_adapter(self):
        seq = R120_4585_CLIP + _body() + 'A' * 30 + _ADAPTER_FWD
        pre = cs.pretrim_consensus(seq, 'fwd', 1)
        # SSP starts at 87 in the clip: 87 + 23 + 27 + 3 = 140 of the 143-nt clip
        assert pre.trim_5p == 143
        assert not pre.seq.startswith(R120_4585_CLIP[:20])


class TestFinder:
    def test_exact_hit_is_the_historical_span(self):
        seq = _fwd_molecule(SSP_FWD, 'C' * UMI_LEN, _body())
        assert find_ssp_span(seq, 'fwd') == (0, len(SSP_FWD))
        rc = revcomp(seq)
        s, e = find_ssp_span(rc, 'rev')
        assert rc[s:e] == SSP_RC

    @pytest.mark.parametrize('mutant', [
        SSP_FWD[:-1] + 'C',                 # substitution at the last base (r120_4585)
        SSP_FWD[:10] + SSP_FWD[11:],        # one deletion
        SSP_FWD[:12] + 'A' + SSP_FWD[12:],  # one insertion
        SSP_FWD[:3] + 'G' + SSP_FWD[4:-2] + 'AC',  # two substitutions
    ])
    def test_error_bearing_ssp_still_yields_type1_and_a_full_umi(self, mutant):
        umi = 'ACGTTGCAACGTTGCAACGTTGCAACG'[:UMI_LEN]
        seq = _fwd_molecule(mutant, umi, _body())
        info = extract_read_info(_read(seq, len(mutant) + UMI_LEN + BRIDGE_LEN, 52))
        assert info.read_type == 1 and info.orient == 'fwd'
        assert len(info.umi) == UMI_LEN
        # an indel inside the SSP shifts the UMI by at most one base
        assert umi[:UMI_LEN - 1] in info.umi or umi[1:] in info.umi or info.umi == umi

    def test_rev_frame_with_an_error_is_type1(self):
        umi = 'ACGTTGCAACGTTGCAACGTTGCAACG'[:UMI_LEN]
        mutant = SSP_FWD[:-1] + 'C'
        seq = revcomp(_fwd_molecule(mutant, umi, _body()))
        info = extract_read_info(_read(seq, 52, len(mutant) + UMI_LEN + BRIDGE_LEN))
        assert info.read_type == 1 and info.orient == 'rev'
        assert len(info.umi) == UMI_LEN

    def test_a_molecule_without_an_ssp_stays_type2(self):
        seq = _body(600, seed=7) + 'A' * 30 + _ADAPTER_FWD
        info = extract_read_info(_read(seq, 0, 52))
        assert info.read_type == 2 and info.umi == ''

    def test_fuzzy_search_is_windowed_to_the_ssp_end(self):
        # an error-bearing SSP buried 500 nt into the read is NOT the 5' adapter
        mutant = SSP_FWD[:-1] + 'C'
        seq = _body(500, seed=3) + mutant + 'C' * UMI_LEN + 'GGG' + _body(300, seed=4)
        assert find_ssp_span(seq, 'fwd') == (-1, -1)
        # ...but an EXACT one anywhere is still honoured (historical behaviour)
        seq_exact = _body(500, seed=3) + SSP_FWD + 'C' * UMI_LEN + 'GGG' + _body(300, seed=4)
        assert find_ssp_span(seq_exact, 'fwd') == (500, 500 + len(SSP_FWD))

    def test_trimmer_and_classifier_share_one_finder(self):
        mutant = SSP_FWD[:-1] + 'C'
        seq = _fwd_molecule(mutant, 'C' * UMI_LEN, _body())
        assert cs._find_ssp(seq, 'fwd') == find_ssp_span(seq, 'fwd')[0] == 0
        assert cs._detect_frame(seq, 'rev', 1) == ('fwd', True)

    def test_without_edlib_exact_only(self, monkeypatch):
        monkeypatch.setattr(ri, 'HAS_EDLIB', False)
        mutant = SSP_FWD[:-1] + 'C'
        assert find_ssp_span(_fwd_molecule(mutant, 'C' * UMI_LEN, _body()), 'fwd') == (-1, -1)
        assert find_ssp_span(_fwd_molecule(SSP_FWD, 'C' * UMI_LEN, _body()), 'fwd') == (0, 23)
