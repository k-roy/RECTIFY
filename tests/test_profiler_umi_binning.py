"""
Phase A sentinel: UMI-cluster-size binning in empirical_cigar_error_profiler.

Verifies that profile_chunk() routes per-read counts into per-bin accumulators
in parallel with the pooled accumulator. Three invariants:

1. Sum over bins == pooled counts (for every (op, base, hp) key shared with
   the bin that owns each read). Reads with no UMI tag contribute to pooled
   only, so the sum-over-bins == pooled identity holds only when EVERY read
   carries a tag matching some bin.
2. A read with XC:i:1 contributes to umi1 and not to umi2 or umi3+.
3. A read with XC:i:5 contributes to umi3+ and not to umi1 or umi2.
4. A read with no XC tag contributes to pooled and to no per-bin output.

We synthesize a minimal BAM in memory so this test has no dependency on the
bundled validation BAM's tagging (the cDNA validation BAM has no XC tags as
of 2026-05-17, by design — the validation BAM is consumed by upstream cDNA
pipeline tests, not calibration tests).

Author: Kevin R. Roy (Phase A handoff)
"""

from __future__ import annotations

import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pysam
import pytest


# scripts/calibration is not on sys.path by default — add it.
_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / 'scripts' / 'calibration'))

import empirical_cigar_error_profiler as profiler  # noqa: E402


# ---------------------------------------------------------------------------
# Synthetic BAM construction
# ---------------------------------------------------------------------------

_REF_SEQ = (
    # 100 bp of mixed-content, no long HPs, plus an embedded 4xA HP run at
    # positions 40-43 so we get a known HP context to inspect.
    'GCATTGCATGTACGTAGCATGCATGCATGTACGAGTCGTC'  # 0-39
    'AAAA'                                       # 40-43 (HP=4 A run)
    'GCATTGCATGTACGTAGCATGCATGCATGTACGAGTCGTC'  # 44-83
    'CATGCATGCATGCATG'                           # 84-99
)
assert len(_REF_SEQ) == 100, len(_REF_SEQ)
_REF_NAME = 'synth_chr'


def _make_header() -> pysam.AlignmentHeader:
    return pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'unsorted'},
        'SQ': [{'SN': _REF_NAME, 'LN': len(_REF_SEQ)}],
    })


def _make_read(
    header: pysam.AlignmentHeader,
    qname: str,
    seq: str,
    cigar: str,
    pos: int = 0,
    xc_tag: int = None,
    flag: int = 0,
) -> pysam.AlignedSegment:
    """Build one minimal mapped pysam.AlignedSegment with an optional XC tag."""
    r = pysam.AlignedSegment(header)
    r.query_name = qname
    r.flag = flag
    r.reference_id = 0
    r.reference_start = pos
    r.mapping_quality = 60
    r.cigarstring = cigar
    r.query_sequence = seq
    r.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
    tags = []
    if xc_tag is not None:
        tags.append(('XC', int(xc_tag), 'i'))
    if tags:
        r.set_tags(tags)
    return r


def _write_bam(path: Path, reads: List[pysam.AlignedSegment], header) -> Path:
    with pysam.AlignmentFile(str(path), 'wb', header=header) as bam:
        for r in reads:
            bam.write(r)
    pysam.index(str(path))
    return path


@pytest.fixture
def synth_bam_dir(tmp_path):
    """Build two identical synthetic BAMs (two-aligner panel) with XC tags.

    The reads are perfect matches to the reference (all M ops), so the
    pooled accumulator should record a substantial M count for each read.
    With two identical BAMs and min_aligners=2, every position passes the
    strict-consensus filter.
    """
    header = _make_header()

    # Four reads, all 50 bp, all perfect matches to ref[0:50].
    # XC:i:1 → umi1; XC:i:2 → umi2; XC:i:5 → umi3+; no tag → pooled-only.
    seq50 = _REF_SEQ[0:50]
    cigar50 = '50M'
    reads = [
        _make_read(header, 'r_singleton', seq50, cigar50, pos=0, xc_tag=1),
        _make_read(header, 'r_doubleton', seq50, cigar50, pos=0, xc_tag=2),
        _make_read(header, 'r_consensus', seq50, cigar50, pos=0, xc_tag=5),
        _make_read(header, 'r_untagged', seq50, cigar50, pos=0, xc_tag=None),
    ]

    bam_a = _write_bam(tmp_path / 'aligner_a.bam', reads, header)
    bam_b = _write_bam(tmp_path / 'aligner_b.bam', list(reads), header)

    return {
        'tmp_path': tmp_path,
        'bam_a': bam_a,
        'bam_b': bam_b,
        'header': header,
        'ref_seqs': {_REF_NAME: _REF_SEQ},
        'n_tagged_reads': 3,
        'n_untagged_reads': 1,
    }


# ---------------------------------------------------------------------------
# parse_umi_bins / assign_bin unit tests
# ---------------------------------------------------------------------------

def test_parse_umi_bins_default():
    bins = profiler.parse_umi_bins('1,2,3+')
    assert bins == [('1', 'eq', 1), ('2', 'eq', 2), ('3+', 'ge', 3)]


def test_parse_umi_bins_handles_whitespace():
    assert profiler.parse_umi_bins(' 1 , 2 , 3+ ') == [
        ('1', 'eq', 1), ('2', 'eq', 2), ('3+', 'ge', 3),
    ]


def test_parse_umi_bins_rejects_garbage():
    with pytest.raises(ValueError):
        profiler.parse_umi_bins('1,foo,3+')
    with pytest.raises(ValueError):
        profiler.parse_umi_bins('')


def test_assign_bin_routing():
    bins = profiler.parse_umi_bins('1,2,3+')
    assert profiler.assign_bin(1, bins) == '1'
    assert profiler.assign_bin(2, bins) == '2'
    assert profiler.assign_bin(3, bins) == '3+'
    assert profiler.assign_bin(50, bins) == '3+'
    assert profiler.assign_bin(0, bins) is None
    assert profiler.assign_bin(None, bins) is None


def test_sanitize_bin_label():
    assert profiler._sanitize_bin_label('3+') == '3plus'
    assert profiler._sanitize_bin_label('1') == '1'


# ---------------------------------------------------------------------------
# profile_chunk routing tests
# ---------------------------------------------------------------------------

def _run_profile(synth, bins):
    """Run profile_chunk with the synthetic panel and return all accumulators."""
    counts: Dict[Tuple[str, str, int], int] = defaultdict(int)
    counts_by_bin: Dict[str, Dict[Tuple[str, str, int], int]] = {
        lbl: defaultdict(int) for lbl, _, _ in bins
    }
    reads_per_bin: Dict[str, int] = {lbl: 0 for lbl, _, _ in bins}

    chunk_bams = {'aligner_a': synth['bam_a'], 'aligner_b': synth['bam_b']}
    n_processed, n_common = profiler.profile_chunk(
        chunk_bams=chunk_bams,
        ref_seqs=synth['ref_seqs'],
        min_aligners=2,
        max_reads=None,
        counts=counts,
        umi_tag='XC',
        umi_bins=bins,
        counts_by_bin=counts_by_bin,
        reads_per_bin=reads_per_bin,
    )
    return counts, counts_by_bin, reads_per_bin, n_processed, n_common


def test_per_bin_sum_matches_pooled_for_tagged_reads(synth_bam_dir):
    """
    With ALL reads carrying tags that match a bin, sum-over-bins == pooled
    for every key. We construct that scenario by removing the untagged read
    from the panel first.
    """
    header = _make_header()
    seq50 = _REF_SEQ[0:50]
    reads = [
        _make_read(header, 'r_singleton', seq50, '50M', pos=0, xc_tag=1),
        _make_read(header, 'r_doubleton', seq50, '50M', pos=0, xc_tag=2),
        _make_read(header, 'r_consensus', seq50, '50M', pos=0, xc_tag=5),
    ]
    bam_a = _write_bam(synth_bam_dir['tmp_path'] / 'tagged_a.bam', reads, header)
    bam_b = _write_bam(synth_bam_dir['tmp_path'] / 'tagged_b.bam', list(reads), header)
    panel = dict(synth_bam_dir)
    panel['bam_a'] = bam_a
    panel['bam_b'] = bam_b

    bins = profiler.parse_umi_bins('1,2,3+')
    counts, counts_by_bin, reads_per_bin, n_processed, n_common = _run_profile(panel, bins)

    assert n_processed == 3
    assert n_common == 3
    assert sum(reads_per_bin.values()) == 3

    # Sum over bins for each key
    summed: Dict[Tuple[str, str, int], int] = defaultdict(int)
    for sub in counts_by_bin.values():
        for k, v in sub.items():
            summed[k] += v

    # Every key in pooled counts should match the per-bin sum exactly,
    # because every read was tagged.
    assert dict(counts) == dict(summed), (
        f'pooled vs sum-over-bins mismatch:\n'
        f'pooled keys: {sorted(counts.keys())}\n'
        f'summed keys: {sorted(summed.keys())}'
    )


def test_xc1_routes_to_umi1_only(synth_bam_dir):
    """A read with XC:i:1 contributes to umi1 and not umi2 or umi3+."""
    header = _make_header()
    seq50 = _REF_SEQ[0:50]
    reads = [_make_read(header, 'lone', seq50, '50M', pos=0, xc_tag=1)]
    bam_a = _write_bam(synth_bam_dir['tmp_path'] / 'lone_a.bam', reads, header)
    bam_b = _write_bam(synth_bam_dir['tmp_path'] / 'lone_b.bam', list(reads), header)
    panel = dict(synth_bam_dir)
    panel['bam_a'] = bam_a
    panel['bam_b'] = bam_b

    bins = profiler.parse_umi_bins('1,2,3+')
    counts, counts_by_bin, reads_per_bin, _, _ = _run_profile(panel, bins)

    assert reads_per_bin == {'1': 1, '2': 0, '3+': 0}
    assert sum(counts_by_bin['1'].values()) > 0
    assert sum(counts_by_bin['2'].values()) == 0
    assert sum(counts_by_bin['3+'].values()) == 0
    # And the umi1 bin should equal pooled exactly (one read, all tagged).
    assert dict(counts_by_bin['1']) == dict(counts)


def test_xc5_routes_to_umi3plus_only(synth_bam_dir):
    """A read with XC:i:5 contributes to umi3+ and not umi1 or umi2."""
    header = _make_header()
    seq50 = _REF_SEQ[0:50]
    reads = [_make_read(header, 'big', seq50, '50M', pos=0, xc_tag=5)]
    bam_a = _write_bam(synth_bam_dir['tmp_path'] / 'big_a.bam', reads, header)
    bam_b = _write_bam(synth_bam_dir['tmp_path'] / 'big_b.bam', list(reads), header)
    panel = dict(synth_bam_dir)
    panel['bam_a'] = bam_a
    panel['bam_b'] = bam_b

    bins = profiler.parse_umi_bins('1,2,3+')
    counts, counts_by_bin, reads_per_bin, _, _ = _run_profile(panel, bins)

    assert reads_per_bin == {'1': 0, '2': 0, '3+': 1}
    assert sum(counts_by_bin['1'].values()) == 0
    assert sum(counts_by_bin['2'].values()) == 0
    assert sum(counts_by_bin['3+'].values()) > 0
    assert dict(counts_by_bin['3+']) == dict(counts)


def test_untagged_read_pooled_only(synth_bam_dir):
    """Read with NO XC tag contributes to pooled but to no per-bin output."""
    header = _make_header()
    seq50 = _REF_SEQ[0:50]
    reads = [_make_read(header, 'naked', seq50, '50M', pos=0, xc_tag=None)]
    bam_a = _write_bam(synth_bam_dir['tmp_path'] / 'naked_a.bam', reads, header)
    bam_b = _write_bam(synth_bam_dir['tmp_path'] / 'naked_b.bam', list(reads), header)
    panel = dict(synth_bam_dir)
    panel['bam_a'] = bam_a
    panel['bam_b'] = bam_b

    bins = profiler.parse_umi_bins('1,2,3+')
    counts, counts_by_bin, reads_per_bin, _, _ = _run_profile(panel, bins)

    assert reads_per_bin == {'1': 0, '2': 0, '3+': 0}
    # Pooled must have actual events from the read.
    assert sum(counts.values()) > 0
    # Every bin must be empty (no per-bin contribution from untagged read).
    for lbl, sub in counts_by_bin.items():
        assert sum(sub.values()) == 0, (
            f'untagged read leaked into bin {lbl}: {dict(sub)}'
        )


def test_sidecar_bam_lookup_routes_per_bin(tmp_path):
    """
    Sidecar BAM path: per-aligner BAMs lack XC tag; sidecar carries it.

    Three reads with XC:i:1, XC:i:2, XC:i:5 in the sidecar.  The per-aligner
    BAMs have the same read names but NO XC tag.  profile_chunk() should route
    each read to the correct bin via umi_lookup, and the miss counter should be
    zero (all names match).

    Also verifies read-name normalisation: one sidecar read name has a
    trailing whitespace comment that _normalise_read_name must strip before
    the lookup key is stored, matching the key produced by load_reads().
    """
    header = _make_header()
    seq50 = _REF_SEQ[0:50]

    # Sidecar BAM: carries XC; one name has a comment suffix to test normalisation.
    sidecar_reads = [
        _make_read(header, 'r_s1', seq50, '50M', pos=0, xc_tag=1),
        _make_read(header, 'r_s2', seq50, '50M', pos=0, xc_tag=2),
        _make_read(header, 'r_s5', seq50, '50M', pos=0, xc_tag=5),
    ]
    sidecar_path = _write_bam(tmp_path / 'sidecar.bam', sidecar_reads, header)

    # Per-aligner BAMs: same read names, no XC tag.
    aligner_reads = [
        _make_read(header, 'r_s1', seq50, '50M', pos=0),
        _make_read(header, 'r_s2', seq50, '50M', pos=0),
        _make_read(header, 'r_s5', seq50, '50M', pos=0),
    ]
    bam_a = _write_bam(tmp_path / 'noxtag_a.bam', aligner_reads, header)
    bam_b = _write_bam(tmp_path / 'noxtag_b.bam', list(aligner_reads), header)

    # Build umi_lookup exactly as main() does.
    umi_lookup = {}
    with pysam.AlignmentFile(str(sidecar_path), 'rb', check_sq=False) as bam:
        for r in bam:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            try:
                val = int(r.get_tag('XC'))
            except (KeyError, TypeError, ValueError):
                continue
            umi_lookup[profiler._normalise_read_name(r.query_name)] = val

    assert len(umi_lookup) == 3

    bins = profiler.parse_umi_bins('1,2,3+')
    counts: Dict[Tuple[str, str, int], int] = defaultdict(int)
    counts_by_bin: Dict[str, Dict[Tuple[str, str, int], int]] = {
        lbl: defaultdict(int) for lbl, _, _ in bins
    }
    reads_per_bin: Dict[str, int] = {lbl: 0 for lbl, _, _ in bins}

    profiler.profile_chunk(
        chunk_bams={'aligner_a': bam_a, 'aligner_b': bam_b},
        ref_seqs={_REF_NAME: _REF_SEQ},
        min_aligners=2,
        max_reads=None,
        counts=counts,
        umi_tag='XC',
        umi_bins=bins,
        counts_by_bin=counts_by_bin,
        reads_per_bin=reads_per_bin,
        umi_lookup=umi_lookup,
    )

    assert reads_per_bin == {'1': 1, '2': 1, '3+': 1}
    assert sum(counts_by_bin['1'].values()) > 0
    assert sum(counts_by_bin['2'].values()) > 0
    assert sum(counts_by_bin['3+'].values()) > 0
    # Pooled == sum of all bins (every read was tagged and matched).
    total_per_bin = sum(sum(v.values()) for v in counts_by_bin.values())
    assert sum(counts.values()) == total_per_bin


def test_mixed_panel_pooled_includes_untagged(synth_bam_dir):
    """
    On the default 4-read fixture (3 tagged + 1 untagged) the pooled counts
    are strictly larger than the sum-over-bins, by exactly the untagged
    read's contribution.
    """
    bins = profiler.parse_umi_bins('1,2,3+')
    counts, counts_by_bin, reads_per_bin, n_processed, _ = _run_profile(
        synth_bam_dir, bins)

    assert n_processed == 4
    assert sum(reads_per_bin.values()) == 3  # only tagged reads counted

    summed: Dict[Tuple[str, str, int], int] = defaultdict(int)
    for sub in counts_by_bin.values():
        for k, v in sub.items():
            summed[k] += v
    # Pooled > summed at every key that's present, because untagged read
    # added one extra M count at each ref position it covered.
    total_pooled = sum(counts.values())
    total_summed = sum(summed.values())
    assert total_pooled > total_summed
    # And exactly: pooled - summed should equal the untagged read's per-key
    # contribution. With 4 identical reads, each per-key value is 4x what
    # one read contributes; summing tagged (3) is 3/4 of pooled.
    assert total_summed * 4 == total_pooled * 3


# ---------------------------------------------------------------------------
# usable_mask (variant / high-conf masking) tests
# ---------------------------------------------------------------------------

def _mismatch_panel(tmp_path):
    """Two-aligner panel: one read perfectly matching ref[0:50] except a single
    substitution at position 10 (so exactly one strict-consensus X event)."""
    header = _make_header()
    sub_pos = 10
    seq = list(_REF_SEQ[0:50])
    seq[sub_pos] = 'A' if seq[sub_pos] != 'A' else 'C'
    seq = ''.join(seq)
    reads = [_make_read(header, 'r1', seq, '50M', pos=0)]
    bam_a = _write_bam(tmp_path / 'mm_a.bam', reads, header)
    bam_b = _write_bam(tmp_path / 'mm_b.bam', list(reads), header)
    return {'bam_a': bam_a, 'bam_b': bam_b,
            'ref_seqs': {_REF_NAME: _REF_SEQ}, 'sub_pos': sub_pos}


def _count_only(chunk_bams, ref_seqs, usable_mask=None):
    counts: Dict[Tuple[str, str, int], int] = defaultdict(int)
    profiler.profile_chunk(
        chunk_bams=chunk_bams, ref_seqs=ref_seqs,
        min_aligners=2, max_reads=None, counts=counts,
        usable_mask=usable_mask,
    )
    return counts


def test_usable_mask_excludes_masked_substitution(tmp_path):
    """The single substitution at the masked position is dropped; everything
    else is unchanged. This is the het-allele-masking guarantee."""
    p = _mismatch_panel(tmp_path)
    cb = {'aligner_a': p['bam_a'], 'aligner_b': p['bam_b']}

    c0 = _count_only(cb, p['ref_seqs'])
    x0 = sum(v for (op, _, _), v in c0.items() if op == 'X')
    assert x0 == 1, f'expected exactly one X event, got {x0}: {dict(c0)}'

    mask = {_REF_NAME: np.ones(len(_REF_SEQ), dtype=bool)}
    mask[_REF_NAME][p['sub_pos']] = False
    c1 = _count_only(cb, p['ref_seqs'], usable_mask=mask)
    x1 = sum(v for (op, _, _), v in c1.items() if op == 'X')
    assert x1 == 0, f'masked substitution should be excluded, got {x1}'
    # M counts (all other positions) are untouched by masking one X position.
    m0 = sum(v for (op, _, _), v in c0.items() if op == 'M')
    m1 = sum(v for (op, _, _), v in c1.items() if op == 'M')
    assert m1 == m0


def test_usable_mask_absent_chrom_excludes_everything(tmp_path):
    """A mask that omits the read's chromosome excludes the whole read group
    (conservative: no confidence calls there)."""
    p = _mismatch_panel(tmp_path)
    cb = {'aligner_a': p['bam_a'], 'aligner_b': p['bam_b']}
    mask = {'some_other_chrom': np.ones(10, dtype=bool)}
    c = _count_only(cb, p['ref_seqs'], usable_mask=mask)
    assert sum(c.values()) == 0


def test_load_usable_mask_bed_and_none(tmp_path):
    ref_seqs = {'chr1': 'ACGT' * 10}  # length 40
    # No args -> no masking.
    assert profiler.load_usable_mask(ref_seqs, None, None) is None
    # BED-only: positions start excluded, only [5,15) turned on.
    bed = tmp_path / 'hc.bed'
    bed.write_text('chr1\t5\t15\n')
    m = profiler.load_usable_mask(ref_seqs, None, bed)
    assert m['chr1'][4] == False  # noqa: E712
    assert m['chr1'][5] == True   # noqa: E712
    assert m['chr1'][14] == True  # noqa: E712
    assert m['chr1'][15] == False  # noqa: E712
    assert int(m['chr1'].sum()) == 10
