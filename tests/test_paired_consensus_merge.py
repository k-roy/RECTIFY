"""Paired-end K-way merge: both mates must survive consensus grouping.

Regression test for the ship-blocker the adversarial review found: R1 and R2
share one RN, so the per-aligner `max(candidates, key=n_ops)` collapse used to
keep only the more-spliced mate and silently drop the other — losing ~half the
paired junction evidence. The merge must now emit one group PER MATE.
"""

import pysam

from rectify.core.consensus.consensus import _iter_name_grouped_bams


def _hdr():
    return pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'SQ': [{'SN': 'chr5', 'LN': 1_000_000}],
    })


def _rec(hdr, qname, rn, is_read2, cigar, start=1000):
    r = pysam.AlignedSegment(hdr)
    r.query_name = qname
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigartuples = cigar
    r.query_sequence = 'A' * sum(l for op, l in cigar if op in (0, 1, 4, 7, 8))
    # paired primary, mate flags
    r.is_paired = True
    r.is_read1 = not is_read2
    r.is_read2 = is_read2
    r.is_unmapped = False
    r.is_secondary = False
    r.is_supplementary = False
    r.set_tag('RN', rn, value_type='i')
    return r


def _write_bam(path, records):
    hdr = _hdr()
    with pysam.AlignmentFile(str(path), 'wb', header=hdr) as bam:
        for spec in records:
            bam.write(_rec(hdr, *spec))
    return str(path)


# R1 spans a junction (N op); R2 is ungapped. Same RN=0, same QNAME 'pair0'.
_R1_CIGAR = [(0, 50), (3, 200), (0, 50)]   # 50M 200N 50M  → 1 junction
_R2_CIGAR = [(0, 100)]                      # 100M          → 0 junctions


def test_both_mates_survive_paired_merge(tmp_path):
    # Two aligners, each with both mates of pair0 (R1 then R2, name-sorted).
    a = _write_bam(tmp_path / 'A.bam', [
        ('pair0', 0, False, _R1_CIGAR),
        ('pair0', 0, True, _R2_CIGAR),
    ])
    b = _write_bam(tmp_path / 'B.bam', [
        ('pair0', 0, False, _R1_CIGAR),
        ('pair0', 0, True, _R2_CIGAR),
    ])

    groups = list(_iter_name_grouped_bams({'A': a, 'B': b}))

    # Exactly TWO groups for the pair — one per mate (NOT one collapsed group).
    assert len(groups) == 2, f"expected 2 mate groups, got {len(groups)}"

    keys = [k for k, _ in groups]
    mate_keys = sorted(k[1] for k in keys)
    assert mate_keys == [False, True], f"want R1+R2 mate keys, got {mate_keys}"

    by_mate = {k[1]: grp for k, grp in groups}
    # Each mate group has BOTH aligners (no mate dropped)
    assert set(by_mate[False]) == {'A', 'B'}
    assert set(by_mate[True]) == {'A', 'B'}

    # R1 group keeps the spliced record (junction evidence preserved)
    r1_rec = by_mate[False]['A']
    assert any(op == 3 for op, _ in r1_rec.cigartuples), "R1 junction lost"
    assert r1_rec.is_read1
    # R2 group keeps the ungapped mate (would have been dropped before the fix)
    r2_rec = by_mate[True]['A']
    assert not any(op == 3 for op, _ in r2_rec.cigartuples)
    assert r2_rec.is_read2


def test_unpaired_single_group_unchanged(tmp_path):
    # Single-end (is_paired=False) → one group keyed mate_key=None, both aligners.
    def _se(path):
        hdr = _hdr()
        with pysam.AlignmentFile(str(path), 'wb', header=hdr) as bam:
            r = pysam.AlignedSegment(hdr)
            r.query_name = 'read0'
            r.reference_id = 0
            r.reference_start = 1000
            r.mapping_quality = 60
            r.cigartuples = [(0, 100)]
            r.query_sequence = 'A' * 100
            r.is_paired = False
            r.is_unmapped = False
            r.is_secondary = False
            r.is_supplementary = False
            r.set_tag('RN', 7, value_type='i')
            bam.write(r)
        return str(path)

    groups = list(_iter_name_grouped_bams({'A': _se(tmp_path / 'a.bam'),
                                           'B': _se(tmp_path / 'b.bam')}))
    assert len(groups) == 1
    (key, grp), = groups
    # Unpaired keeps the bare key (no mate tuple) — original contract preserved.
    assert key == 7
    assert set(grp) == {'A', 'B'}
