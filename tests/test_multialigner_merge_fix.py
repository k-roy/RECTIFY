"""Regression tests for the cDNA multi-aligner merge fix (planning/251).

Covers the two independent losses the debug panel (planning/250a-c) pinned in the
cDNA Path-B align2 pipeline, and the defense-in-depth guard:

1. PRIMARY (87% drop): non-globally-unique Stage-1 ``cluster_<cid>`` names from
   per-region ``correct-cdna --region`` output concatenated without a region
   prefix -> the K-way consensus merge keys on the colliding name/RN and
   collapses distinct molecules to one record.
   Fixed by ``_region_cluster_prefix`` (cdna_correct_command). Guarded by
   ``_detect_duplicate_molecule_names`` (consensus) so a regression fails LOUD.

2. SECONDARY (1.53x): a read mapped by minimap2 but WON by a junction aligner
   (uLTRA) that does not propagate the FASTQ comment reaches cdna-analyze without
   the required XU/XO/XT/XY/XC/XF tags (or with colliding uLTRA XC:Z:NO_SPLICE)
   and is dropped. Fixed by ``_restore_comment_tags_from_siblings`` (consensus):
   the winner inherits the authoritative tags from a comment-propagating sibling
   -> the read minimap2 mapped survives. Robust for ANY aligner combination.
"""
from collections import defaultdict

import pysam
import pytest

from rectify.core.consensus import consensus as c
from rectify.core.consensus.consensus import (
    ConsensusResult,
    _detect_duplicate_molecule_names,
    _iter_name_grouped_bams,
    _restore_comment_tags_from_siblings,
)
from rectify.core.commands.cdna_correct_command import _region_cluster_prefix


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _header():
    return pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'SQ': [{'SN': 'chrI', 'LN': 100000}, {'SN': 'chrII', 'LN': 100000}],
    })


def _read(header, name, ref_id=0, pos=10, rn=None, tags=None):
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.reference_id = ref_id
    r.reference_start = pos
    r.mapping_quality = 60
    r.cigarstring = '4M'
    r.query_sequence = 'ACGT'
    r.query_qualities = pysam.qualitystring_to_array('IIII')
    if rn is not None:
        r.set_tag('RN', rn, value_type='i')
    for t, (v, vt) in (tags or {}).items():
        r.set_tag(t, v, value_type=vt)
    return r


def _write_bam(path, reads):
    with pysam.AlignmentFile(path, 'wb', header=_header()) as out:
        for r in reads:
            out.write(r)


# --------------------------------------------------------------------------- #
# PRIMARY fix — globally-unique Stage-1 names
# --------------------------------------------------------------------------- #
def test_region_cluster_prefix_uniquifies_per_region():
    # Bare "cluster" only when unrestricted (whole-BAM single call).
    assert _region_cluster_prefix(None) == "cluster"
    assert _region_cluster_prefix("") == "cluster"
    # Per-region -> region-namespaced, so cluster_0 from chrI != cluster_0 chrII.
    assert _region_cluster_prefix("chrI") == "cluster_chrI"
    assert _region_cluster_prefix("chrII") == "cluster_chrII"
    assert _region_cluster_prefix("chrI") != _region_cluster_prefix("chrII")
    # Sanitized to a normalize-safe token (no ':'/'-' that a later step mangles).
    assert _region_cluster_prefix("chrI:1-1000") == "cluster_chrI_1_1000"
    tok = _region_cluster_prefix("chrI:1-1000")
    assert c._normalize_bam_read_name(tok + "_0") == tok + "_0"


def test_duplicate_molecule_names_guard_raises_on_colliding_names(tmp_path):
    # Simulate per-region output concatenated WITHOUT the prefix fix: two
    # distinct molecules (different chroms/loci) both named 'cluster_0'.
    header = _header()
    a = tmp_path / 'mm2.namesorted.bam'
    b = tmp_path / 'ultra.namesorted.bam'
    reads_a = [_read(header, 'cluster_0', ref_id=0, pos=100, rn=0),
               _read(header, 'cluster_0', ref_id=1, pos=500, rn=0)]
    reads_b = [_read(header, 'cluster_0', ref_id=0, pos=100, rn=0),
               _read(header, 'cluster_0', ref_id=1, pos=500, rn=0)]
    _write_bam(a, reads_a)
    _write_bam(b, reads_b)
    with pytest.raises(RuntimeError, match="Duplicate molecule names"):
        list(_iter_name_grouped_bams({'minimap2': str(a), 'uLTRA': str(b)}))


def test_unique_region_prefixed_names_pass_guard_and_keep_both(tmp_path):
    # WITH the prefix fix: the two molecules carry distinct names and distinct
    # RN -> guard passes and the merge yields TWO groups (nothing collapsed).
    header = _header()
    a = tmp_path / 'mm2.namesorted.bam'
    b = tmp_path / 'ultra.namesorted.bam'
    reads_a = [_read(header, 'cluster_chrI_0', ref_id=0, pos=100, rn=0),
               _read(header, 'cluster_chrII_0', ref_id=1, pos=500, rn=1)]
    reads_b = [_read(header, 'cluster_chrI_0', ref_id=0, pos=100, rn=0),
               _read(header, 'cluster_chrII_0', ref_id=1, pos=500, rn=1)]
    _write_bam(a, reads_a)
    _write_bam(b, reads_b)
    groups = list(_iter_name_grouped_bams({'minimap2': str(a), 'uLTRA': str(b)}))
    assert len(groups) == 2  # both distinct molecules retained
    keys = {g[0] for g in groups}
    assert keys == {0, 1}


def test_guard_exempts_paired_data(tmp_path):
    # Paired reads legitimately share a name (R1/R2); the guard must not fire.
    header = _header()
    a = tmp_path / 'a.bam'
    r1 = _read(header, 'pair0', ref_id=0, pos=100, rn=0)
    r1.is_paired = True
    r1.is_read1 = True
    r2 = _read(header, 'pair0', ref_id=1, pos=900, rn=0)  # mate at a different locus
    r2.is_paired = True
    r2.is_read2 = True
    _write_bam(a, [r1, r2])
    # Should NOT raise despite the shared name at two loci.
    _detect_duplicate_molecule_names({'a': str(a)}, use_rn_key=True)


# --------------------------------------------------------------------------- #
# SECONDARY fix — tag-restore so a minimap2-mapped read is not dropped when a
# junction aligner wins it
# --------------------------------------------------------------------------- #
def _mm2_sibling(header):
    return _read(header, 'cluster_chrI_5', ref_id=0, pos=42, rn=5, tags={
        'XU': ('ACGTUMI', 'Z'), 'XO': ('rev', 'Z'), 'XT': (1, 'i'),
        'XY': ('sub', 'Z'), 'XC': (3, 'i'), 'XF': (2, 'i'),
    })


def test_ultra_won_read_rescued_with_minimap2_tags(tmp_path):
    # uLTRA-won winner: lacks XU/XO/XT/XY/XF and carries a COLLIDING XC:Z:NO_SPLICE.
    header = _header()
    winner = _read(header, 'cluster_chrI_5', ref_id=0, pos=44, rn=5, tags={
        'XC': ('NO_SPLICE', 'Z'), 'XA': ('junk', 'Z'),
    })
    mm = _mm2_sibling(header)
    _restore_comment_tags_from_siblings(winner, {'minimap2': mm, 'uLTRA': winner})
    # Read minimap2 mapped survives: required tags present with authoritative
    # values, and the crash-inducing uLTRA XC:Z:NO_SPLICE is overwritten.
    assert winner.get_tag('XU') == 'ACGTUMI'
    assert winner.get_tag('XC') == 3
    assert winner.get_tag('XF') == 2
    assert winner.get_tag('XO') == 'rev'


def test_ultra_won_read_rescued_any_combo_without_minimap2(tmp_path):
    # Mandate: robust for ANY aligner combination. No minimap2 member; the
    # comment-propagating sibling (a mapPacBio-family record carrying XU) is used.
    header = _header()
    winner = _read(header, 'cluster_chrI_5', ref_id=0, pos=44, rn=5, tags={
        'XC': ('NO_SPLICE', 'Z'),
    })
    pac = _read(header, 'cluster_chrI_5', ref_id=0, pos=42, rn=5, tags={
        'XU': ('GGGGUMI', 'Z'), 'XC': (5, 'i'), 'XF': (1, 'i'),
    })
    _restore_comment_tags_from_siblings(winner, {'uLTRA': winner, 'mapPacBio': pac})
    assert winner.get_tag('XU') == 'GGGGUMI'
    assert winner.get_tag('XC') == 5


def test_authoritative_winner_tags_untouched(tmp_path):
    # A minimap2 winner already carries XU -> nothing is overwritten.
    header = _header()
    winner = _read(header, 'cluster_chrI_5', ref_id=0, pos=42, rn=5, tags={
        'XU': ('TTTTUMI', 'Z'), 'XC': (9, 'i'),
    })
    other = _mm2_sibling(header)
    _restore_comment_tags_from_siblings(winner, {'minimap2': winner, 'uLTRA': other})
    assert winner.get_tag('XU') == 'TTTTUMI'
    assert winner.get_tag('XC') == 9


def test_process_batch_rescues_ultra_won_read_end_to_end(tmp_path, monkeypatch):
    # Integration through _process_and_write_batch: a molecule mapped by BOTH
    # aligners but WON by uLTRA is written WITH the minimap2 tags (not dropped).
    header = _header()
    mm = _mm2_sibling(header)
    ultra = _read(header, 'cluster_chrI_5', ref_id=0, pos=44, rn=5, tags={
        'XC': ('NO_SPLICE', 'Z'),
    })
    out_path = tmp_path / 'out.bam'
    with pysam.AlignmentFile(out_path, 'wb', header=header) as out:
        monkeypatch.setattr(c, 'select_best_alignment', lambda *a, **k: ConsensusResult(
            read_id='cluster_chrI_5', best_aligner='uLTRA', best_alignment=None,
            aligners_compared=['minimap2', 'uLTRA'], confidence='low', n_aligners_agree=1,
        ))
        stats = {
            'consensus_high': 0, 'consensus_medium': 0, 'consensus_low': 0,
            '5prime_rescued': 0, 'tied_score': 0, 'by_aligner': defaultdict(int),
            'by_aligner_combo': defaultdict(int),
        }
        c._process_and_write_batch(
            [('cluster_chrI_5', {'minimap2': object(), 'uLTRA': object()})],
            [('cluster_chrI_5', {'minimap2': mm, 'uLTRA': ultra})],
            {'chrI': 'A' * 100000}, None, out, stats,
        )
    with pysam.AlignmentFile(out_path, 'rb') as inp:
        rec = next(inp)
        # uLTRA won (coords from uLTRA record) but tags rescued from minimap2.
        assert rec.get_tag('XU') == 'ACGTUMI'
        assert rec.get_tag('XC') == 3
