"""region_skip — the rDNA junction-rescue bypass (planning/644b perf lever)."""

import pytest

from rectify.core.splice.region_skip import (
    YEAST_RDNA_SPEC,
    overlaps_skip_region,
    parse_skip_regions,
    skip_regions_from_env,
)


class TestParse:
    def test_single_region(self):
        r = parse_skip_regions('chrXII:451574-468812')
        assert r == {'chrXII': [(451574, 468812)]}

    def test_yeast_rdna_shorthand(self):
        assert parse_skip_regions('yeast-rdna') == parse_skip_regions(YEAST_RDNA_SPEC)

    def test_mixed_and_multiple(self):
        r = parse_skip_regions('yeast-rdna, chrM:0-100')
        assert set(r) == {'chrXII', 'chrM'}

    def test_colon_in_chrom_name(self):
        # rsplit on ':' keeps exotic contig names with colons intact
        r = parse_skip_regions('HLA:1:5000-6000')
        assert r == {'HLA:1': [(5000, 6000)]}

    @pytest.mark.parametrize('bad', ['chrI', 'chrI:5', 'chrI:9-3', 'chrI:-1-5', ':1-2'])
    def test_malformed_raises(self, bad):
        with pytest.raises(ValueError):
            parse_skip_regions(bad)

    def test_env(self, monkeypatch):
        monkeypatch.setenv('RECTIFY_SKIP_REGIONS', 'yeast-rdna')
        assert 'chrXII' in skip_regions_from_env()
        monkeypatch.delenv('RECTIFY_SKIP_REGIONS')
        assert skip_regions_from_env() == {}


class TestOverlap:
    R = parse_skip_regions('chrXII:451574-468812')

    def test_inside(self):
        assert overlaps_skip_region(self.R, 'chrXII', 455000, 456000)

    def test_straddles_edge(self):
        assert overlaps_skip_region(self.R, 'chrXII', 451000, 451575)

    def test_abuts_half_open(self):
        assert not overlaps_skip_region(self.R, 'chrXII', 468812, 469000)
        assert not overlaps_skip_region(self.R, 'chrXII', 451000, 451574)

    def test_other_chrom(self):
        assert not overlaps_skip_region(self.R, 'chrXI', 455000, 456000)

    def test_empty(self):
        assert not overlaps_skip_region({}, 'chrXII', 455000, 456000)
