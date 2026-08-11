"""Unit tests for the precomputed splice-site index (planning/641 §2.3)."""

import random

import numpy as np
import pytest

from rectify.core.splice.splice_site_index import SpliceSiteIndex


def _rand_genome(n, seed):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


def _brute(seq, dinucs):
    s = seq.upper()
    return sorted(
        p for p in range(len(s) - 1) if s[p:p + 2] in dinucs
    )


@pytest.fixture(scope='module')
def genome():
    return {'chrA': _rand_genome(5000, 42), 'chrB': _rand_genome(1200, 43)}


@pytest.fixture(scope='module')
def index(genome):
    return SpliceSiteIndex.build(genome)


class TestBuild:
    def test_don_gt_plus_matches_brute_force(self, genome, index):
        for chrom, seq in genome.items():
            got = index.sites_in(chrom, 'don_gt_plus', 0, len(seq)).tolist()
            assert got == _brute(seq, {'GT'})

    def test_acc_plus_is_intron_end_exclusive(self, genome, index):
        # acc_plus stores dinuc position + 2 (the intron END coordinate)
        for chrom, seq in genome.items():
            got = index.sites_in(chrom, 'acc_plus', 0, len(seq) + 2).tolist()
            assert got == [p + 2 for p in _brute(seq, {'AG'})]

    def test_don_minus_merges_ac_and_gc(self, genome, index):
        for chrom, seq in genome.items():
            got = index.sites_in(chrom, 'don_minus', 0, len(seq) + 2).tolist()
            assert got == [p + 2 for p in _brute(seq, {'AC', 'GC'})]

    def test_acc_minus_is_intron_start(self, genome, index):
        for chrom, seq in genome.items():
            got = index.sites_in(chrom, 'acc_minus', 0, len(seq)).tolist()
            assert got == _brute(seq, {'CT'})

    def test_don_plus_union_is_sorted(self, genome, index):
        for chrom, seq in genome.items():
            got = index.sites_in(chrom, 'don_plus', 0, len(seq)).tolist()
            assert got == _brute(seq, {'GT', 'GC'})
            assert got == sorted(got)

    def test_lowercase_genome(self):
        idx = SpliceSiteIndex.build({'c': 'ggtaag'.lower()})
        assert idx.sites_in('c', 'don_gt_plus', 0, 10).tolist() == [1]


class TestRangeQuery:
    def test_random_windows_match_brute_force(self, genome, index):
        rng = random.Random(7)
        seq = genome['chrA']
        all_gt = _brute(seq, {'GT'})
        for _ in range(50):
            lo = rng.randint(0, len(seq))
            hi = rng.randint(0, len(seq))
            got = index.sites_in('chrA', 'don_gt_plus', lo, hi).tolist()
            assert got == [p for p in all_gt if lo <= p < hi]

    def test_empty_and_unknown(self, index):
        assert index.sites_in('chrA', 'don_gt_plus', 50, 50).size == 0
        assert index.sites_in('chrNOPE', 'don_gt_plus', 0, 100).size == 0
        assert index.sites_in('chrA', 'acc_plus', -100, 0).size == 0


class TestCache:
    def test_round_trip_and_fingerprint(self, tmp_path, genome):
        fasta = tmp_path / 'g.fa'
        with open(fasta, 'w') as fh:
            for chrom, seq in genome.items():
                fh.write(f'>{chrom}\n{seq}\n')

        idx1 = SpliceSiteIndex.load_or_build(str(fasta), genome)
        cache = tmp_path / 'g.fa.splice_sites.npz'
        assert cache.exists()

        # Reload: fingerprint matches, arrays identical
        idx2 = SpliceSiteIndex.load_or_build(str(fasta), genome)
        for chrom in genome:
            for kind in ('don_gt_plus', 'don_gc_plus', 'acc_plus', 'don_minus', 'acc_minus'):
                np.testing.assert_array_equal(
                    idx1.sites_in(chrom, kind, 0, 10**9),
                    idx2.sites_in(chrom, kind, 0, 10**9),
                )

        # Touch the FASTA (mtime changes) -> stale fingerprint -> rebuild works
        with open(fasta, 'a') as fh:
            fh.write('\n')
        idx3 = SpliceSiteIndex.load_or_build(str(fasta), genome)
        assert idx3.sites_in('chrA', 'don_gt_plus', 0, 10**9).size == \
            idx1.sites_in('chrA', 'don_gt_plus', 0, 10**9).size

    def test_corrupt_cache_rebuilds(self, tmp_path, genome):
        fasta = tmp_path / 'g.fa'
        with open(fasta, 'w') as fh:
            fh.write('>chrA\n' + genome['chrA'] + '\n')
        cache = tmp_path / 'g.fa.splice_sites.npz'
        cache.write_bytes(b'not an npz file')
        idx = SpliceSiteIndex.load_or_build(str(fasta), {'chrA': genome['chrA']})
        assert idx.n_sites('chrA', 'don_gt_plus') > 0
