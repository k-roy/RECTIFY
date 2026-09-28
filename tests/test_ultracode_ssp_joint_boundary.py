"""Resolve fuzzy SSP endpoints from existing construct layout, or refuse."""
import random

import pysam
import pytest

from rectify.core.cdna._constants import ANCHOR_FWD, SSP_FWD
from rectify.core.cdna.consensus import pretrim_consensus
from rectify.core.cdna.io import write_stage1_fastq
from rectify.core.cdna.read_info import (
    detect_full_length_tier, extract_read_info, find_ssp_span, revcomp,
)


UMI = 'CAGTCGACTAGTCGACTGACGACTTTT'
TAIL = 'A' * 25 + ANCHOR_FWD + 'CTGCTCGTGC'


def make_case(kind, frame, flag=0, with_tail=True):
    rng = random.Random(519222)
    body = 'CT' + ''.join(rng.choice('ACGT') for _ in range(394)) + 'CGTC'
    primer, umi, bridge, pad = SSP_FWD, UMI, 'GGG', 'GCGTC'
    if kind.startswith('ins22'):
        primer = SSP_FWD[:22] + 'A' + SSP_FWD[22:]
    elif kind.startswith('del22'):
        primer = SSP_FWD[:22] + SSP_FWD[23:]
    elif kind == 'unique_bridge_error':
        primer = SSP_FWD[:7] + 'A' + SSP_FWD[7:]
    elif kind == 'legacy_sub22_umi_G':
        # Exact old synthetic positive from test_cdna_ssp_fuzzy.py. The last
        # SSP substitution ties with a deletion, and UMI-ending G + GGG makes
        # both 27-nt UMI slices bridge-compatible. Neither is uniquely known.
        primer = SSP_FWD[:-1] + 'C'
        umi = 'ACGTTGCAACGTTGCAACGTTGCAACG'
        rng = random.Random(1)
        body = ''.join(rng.choice('ACGT') for _ in range(400))
        pad = ''
    if kind.endswith('bridge_error'):
        bridge = 'GAG'
    if kind == 'ins22_umi_G':
        umi = umi[:-1] + 'G'
    if kind == 'del22_repetitive_umi':
        umi = 'G' * 27
        body = 'G' + body[1:]
    if kind == 'del22_two_loci':
        pad += primer + umi + bridge + 'CGC' * 3
    prefix = pad + primer + umi + bridge
    tail = TAIL if with_tail else ''
    seq = prefix + body + tail
    if frame == 'rev':
        seq = revcomp(seq)
    header = pysam.AlignmentHeader.from_references(['chrT'], [2000])
    read = pysam.AlignedSegment(header)
    read.query_name = f'{kind}_{frame}_{flag}'
    read.reference_id, read.reference_start, read.flag = 0, 500, flag
    read.mapping_quality = 60
    read.query_sequence, read.query_qualities = seq, [33] * len(seq)
    clip3 = f'{len(tail)}S' if tail else ''
    read.cigarstring = (f'{len(prefix)}S400M{clip3}' if frame == 'fwd'
                        else f'{clip3}400M{len(prefix)}S')
    return read, header, body, prefix, umi


@pytest.mark.parametrize('kind', ['ins22', 'del22', 'unique_bridge_error', 'exact_bridge_error'])
@pytest.mark.parametrize('frame', ['fwd', 'rev'])
def test_unique_joint_boundary_emits_exact_rna(tmp_path, kind, frame):
    read, header, body, prefix, umi = make_case(kind, frame)
    info = extract_read_info(read)
    assert info is not None and info.read_type == 1 and info.umi == umi
    pre = pretrim_consensus(read.query_sequence, frame, 1)
    assert pre.trim_5p == len(prefix)
    assert pre.seq == (body if frame == 'fwd' else revcomp(body))
    bam, fastq = tmp_path/'input.bam', tmp_path/'output.fq'
    with pysam.AlignmentFile(str(bam), 'wb', header=header) as handle:
        handle.write(read)
    write_stage1_fastq(bam, fastq, [[info]], {0: info.umi}, {0: info.xf_tier}, {0: 25},
                       reference=None)
    lines = fastq.read_text().splitlines()
    assert lines[1] == body and len(lines[3]) == 400
    assert f'XQ:i:{len(prefix)}' in lines[0] and 'XN:i:1' in lines[0]


@pytest.mark.parametrize('kind', ['ins22_umi_G', 'del22_repetitive_umi',
                                  'del22_two_loci', 'del22_bridge_error',
                                  'legacy_sub22_umi_G'])
@pytest.mark.parametrize('frame', ['fwd', 'rev'])
def test_ambiguous_fuzzy_boundary_refuses_without_deleting_sequence(kind, frame):
    read, _, body, prefix, _ = make_case(kind, frame)
    assert find_ssp_span(read.query_sequence, frame) == (-1, -1)
    info = extract_read_info(read)
    assert info is not None and info.read_type == 2 and info.umi == ''
    assert info.orient == frame
    assert info.read_subtype == 'umi_boundary_ambiguous'
    # Even a caller retaining the old Type1 label must not cut an arbitrary end.
    for read_type in (1, 2):
        pre = pretrim_consensus(read.query_sequence, frame, read_type)
        assert pre.trim_5p == 0
        assert pre.seq == (prefix + body if frame == 'fwd' else revcomp(prefix + body))


def test_type2_equal_tiers_still_choose_forward():
    read, _, body, _, _ = make_case('exact', 'fwd', 16)
    read.query_sequence = 'T' * 12 + body + 'A' * 12
    read.cigarstring = '12S400M12S'
    assert detect_full_length_tier(read.query_sequence, 'fwd') == 1
    assert detect_full_length_tier(read.query_sequence, 'rev') == 1
    info = extract_read_info(read)
    assert info.read_type == 2 and info.orient == 'fwd' and info.umi == ''


def test_exact_repeated_primer_keeps_historical_first_rna_occurrence():
    seq = 'C' * 500 + SSP_FWD + 'G' * 40 + SSP_FWD + 'G' * 40
    assert find_ssp_span(seq, 'fwd') == (500, 500 + len(SSP_FWD))
    n = len(seq)
    assert find_ssp_span(revcomp(seq), 'rev') == (n - 500 - len(SSP_FWD), n - 500)


@pytest.mark.parametrize('kind', ['legacy_sub22_umi_G', 'ins22'])
@pytest.mark.parametrize('frame', ['fwd', 'rev'])
@pytest.mark.parametrize('flag', [0, 16])
def test_no_tail_keeps_unique_ssp_frame_without_inventing_umi(tmp_path, kind, frame, flag):
    read, header, body, prefix, umi = make_case(kind, frame, flag, with_tail=False)
    assert detect_full_length_tier(read.query_sequence, frame) == 0
    info = extract_read_info(read)
    assert info is not None and info.orient == frame and info.xf_tier == 0
    ambiguous = kind == 'legacy_sub22_umi_G'
    assert info.read_type == (2 if ambiguous else 1)
    assert info.umi == ('' if ambiguous else umi)
    if ambiguous:
        assert info.read_subtype == 'umi_boundary_ambiguous'
    pre = pretrim_consensus(read.query_sequence, frame, info.read_type)
    expected = prefix + body if ambiguous else body
    assert pre.seq == (expected if frame == 'fwd' else revcomp(expected))
    assert pre.trim_5p == (0 if ambiguous else len(prefix)) and pre.trim_3p == 0
    bam, fastq = tmp_path/'input.bam', tmp_path/'output.fq'
    with pysam.AlignmentFile(str(bam), 'wb', header=header) as handle:
        handle.write(read)
    write_stage1_fastq(bam, fastq, [[info]], {0: info.umi}, {0: 0}, {0: 0}, reference=None)
    lines = fastq.read_text().splitlines()
    assert lines[1] == expected and len(lines[3]) == len(expected)
    if ambiguous:
        assert 'XY:Z:umi_boundary_ambiguous' in lines[0]
    assert 'XN:i:1' in lines[0] and f'XT:i:{info.read_type}' in lines[0]


def test_no_tail_and_both_ambiguous_frames_still_refuse_unclassified():
    read, _, body, prefix, _ = make_case('legacy_sub22_umi_G', 'fwd', with_tail=False)
    read.query_sequence = prefix + body + revcomp(prefix)
    read.cigarstring = f'{len(prefix)}S400M{len(prefix)}S'
    assert find_ssp_span(read.query_sequence, 'fwd') == (-1, -1)
    assert find_ssp_span(read.query_sequence, 'rev') == (-1, -1)
    assert detect_full_length_tier(read.query_sequence, 'fwd') == 0
    assert detect_full_length_tier(read.query_sequence, 'rev') == 0
    assert extract_read_info(read) is None
