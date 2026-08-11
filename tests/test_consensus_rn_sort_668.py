"""Regression tests for the 668 DRS consensus group-split defect.

RN is assigned in FASTQ record order (``_build_qname_rn_map``), but the
K-way merge inputs were name-sorted. For plain-uuid QNAMEs (ONT DRS)
natural name order is unrelated to FASTQ order, so the RN-keyed merge
desynchronised wherever a read was missing from a subset of aligners and
silently SPLIT (QNAME, RN) groups — ~2.23 primary records per read on the
668 wt_rep1 panel (8,654,356 records from 3,882,433 reads), each subgroup
crowning its own winner (same RN/pos/flag, different Xa).

Fix: sort consensus inputs by the merge key itself (``_ensure_rn_sorted``
when RN-keyed), and fail loud on any out-of-order key in the merge.
"""

import pysam
import pytest

from rectify.core.consensus.consensus import (
    _ensure_rn_sorted,
    _iter_name_grouped_bams,
)


def _header():
    return pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'SQ': [{'SN': 'chrI', 'LN': 1000}],
    })


def _read(header, name, rn, pos=10):
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.reference_id = 0
    r.reference_start = pos
    r.mapping_quality = 60
    r.cigarstring = '4M'
    r.query_sequence = 'ACGT'
    r.query_qualities = pysam.qualitystring_to_array('IIII')
    r.set_tag('RN', rn, value_type='i')
    return r


def _write_bam(path, reads):
    with pysam.AlignmentFile(path, 'wb', header=_header()) as out:
        for r in reads:
            out.write(r)


# The 668 shape: natural QNAME order ('a-…' < 'b-…' < 'c-…') disagrees with
# RN order (b=1, a=2, c=3), and one read ('a-…') is missing from aligner B.
# Under the old name-sorted RN-keyed merge this split read 'b-…' into two
# groups (one per aligner), each of which crowned its own winner.
_A_READS_NAME_ORDER = [('a-2c43-uuid', 2, 10), ('b-4350-uuid', 1, 50), ('c-4469-uuid', 3, 90)]
_B_READS_NAME_ORDER = [('b-4350-uuid', 1, 50), ('c-4469-uuid', 3, 90)]


def test_name_sorted_rn_keyed_input_fails_loud(tmp_path):
    """Name-sorted input with non-monotone RN must raise, not silently split."""
    header = _header()
    a = tmp_path / 'a.bam'
    b = tmp_path / 'b.bam'
    _write_bam(a, [_read(header, n, rn, pos) for n, rn, pos in _A_READS_NAME_ORDER])
    _write_bam(b, [_read(header, n, rn, pos) for n, rn, pos in _B_READS_NAME_ORDER])
    with pytest.raises(RuntimeError, match='not sorted by the K-way merge key'):
        list(_iter_name_grouped_bams({'minimap2': str(a), 'deSALT': str(b)}))


def test_rn_sorted_merge_yields_each_read_once(tmp_path):
    """After _ensure_rn_sorted, every (QNAME, RN) yields exactly one group."""
    header = _header()
    a = tmp_path / 'a.bam'
    b = tmp_path / 'b.bam'
    _write_bam(a, [_read(header, n, rn, pos) for n, rn, pos in _A_READS_NAME_ORDER])
    _write_bam(b, [_read(header, n, rn, pos) for n, rn, pos in _B_READS_NAME_ORDER])

    sorted_paths = {
        'minimap2': _ensure_rn_sorted(str(a)),
        'deSALT': _ensure_rn_sorted(str(b)),
    }
    groups = list(_iter_name_grouped_bams(sorted_paths, use_rn_key=True))

    keys = [k for k, _ in groups]
    assert keys == [1, 2, 3], f"expected one group per RN in order, got {keys}"
    membership = {k: set(reads) for k, reads in groups}
    assert membership[1] == {'minimap2', 'deSALT'}
    assert membership[2] == {'minimap2'}
    assert membership[3] == {'minimap2', 'deSALT'}


def test_ensure_rn_sorted_orders_by_rn(tmp_path):
    header = _header()
    bam = tmp_path / 'x.bam'
    _write_bam(bam, [_read(header, n, rn, pos) for n, rn, pos in _A_READS_NAME_ORDER])
    sorted_path = _ensure_rn_sorted(str(bam))
    assert sorted_path.endswith('.rnsorted.bam')
    with pysam.AlignmentFile(sorted_path, 'rb') as f:
        rns = [r.get_tag('RN') for r in f]
    assert rns == sorted(rns) == [1, 2, 3]
