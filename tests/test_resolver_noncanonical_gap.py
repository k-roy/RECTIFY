"""The resolver's non-canonical blind spot, pinned as a structural property.

`planning/721` measured, on the same upf1Δ reads through minimap2 and mapPacBio:
the resolver recovers **99.1 %** of the ANNOTATED junctions mapPacBio found and
minimap2 missed, but only **34.5 %** of the NON-CANONICAL ones.

The cause is not a threshold. `SpliceSiteIndex` is built on GT/AG-class
dinucleotides only, so a non-canonical junction has **no entry** and can never be
enumerated as a candidate. This file pins that as a property of the index rather
than re-measuring the 34.5 %, which would need the upf1Δ BAM trio.

WHY THIS TEST EXISTS: `planning/698` described the resolver as doing
"mapPacBio's role, deliberately, with a false-discovery budget". That is true for
annotated junctions and one-third true for non-canonical ones. A future change
that extends the index (a legitimate design choice — it moves both the candidate
space and the `alpha` calibration) should make these tests FAIL, which is the
signal to re-measure the recovery figure and update the docstrings that quote it.
"""

import pytest

from rectify.core.splice.splice_site_index import SpliceSiteIndex


# A synthetic contig carrying one canonical (GT..AG) and one non-canonical
# (AT..AA) intron, with identical flanking exon sequence so the ONLY difference
# between them is the splice dinucleotide.
#
# 🔴 The backbone is constructed to contain NO GT / GC / AG / CT / AC dinucleotide
# anywhere. A first attempt used ordinary random-looking sequence and the test
# failed spuriously: incidental GT/AG pairs inside the exons put indexed sites at
# 29 and 30, within the tolerance of the donor at 31, so the fixture "proved" a
# non-canonical junction was indexed when it had not been. Any edit to these
# strings must preserve the property — `test_fixture_has_no_incidental_sites`
# enforces it.
_EXON_L = 'TCGGATTTTCCATAATTTATTGAAATAAAA'
_INTRON_BODY = ('TAAAATATGAAATCCGATGAAATCGGATATAAATCCCGAAAAAATCCCGATGAATCC'
                'GGATGAATTTC')
_EXON_R = 'TGGGAAATGATGGGATAATCGGATGGGAAT'

CANON = _EXON_L + 'GT' + _INTRON_BODY + 'AG' + _EXON_R
NONCANON = _EXON_L + 'AT' + _INTRON_BODY + 'AA' + _EXON_R

_KINDS = ('don_gt_plus', 'don_gc_plus', 'acc_plus', 'don_minus', 'acc_minus')


def _all_sites(idx, chrom, lo, hi):
    out = set()
    for kind in _KINDS:
        out |= {int(x) for x in idx.sites_in(chrom, kind, lo, hi)}
    return out


def test_fixture_has_no_incidental_splice_dinucleotides():
    """Guards the guard: the backbone must carry no GT/GC/AG/CT/AC pair.

    Without this, an incidental dinucleotide in the exon or intron body puts an
    indexed site next to the junction and the non-canonical assertion below
    passes or fails for the wrong reason (this actually happened on the first
    draft — see the fixture comment).
    """
    forbidden = ('GT', 'GC', 'AG', 'CT', 'AC')
    for name, seq in (('EXON_L', _EXON_L), ('INTRON_BODY', _INTRON_BODY),
                      ('EXON_R', _EXON_R)):
        hits = [seq[i:i + 2] for i in range(len(seq) - 1)
                if seq[i:i + 2] in forbidden]
        assert not hits, f'{name} contains splice dinucleotides {set(hits)}'


def test_canonical_intron_edges_are_indexed():
    """Anti-vacuity: the index must find the CANONICAL junction's edges.

    Without this, the non-canonical test below could pass simply because the
    index was empty or the fixture was malformed.
    """
    idx = SpliceSiteIndex.build({'chrC': CANON})
    donor = len(_EXON_L)                      # position of 'GT'
    acceptor = len(_EXON_L) + 2 + len(_INTRON_BODY) + 2   # intron end (exclusive)
    sites = _all_sites(idx, 'chrC', 0, len(CANON))
    assert sites, 'index produced no sites at all — fixture or build is broken'
    assert any(abs(s - donor) <= 2 for s in sites), (
        f'canonical donor at {donor} not indexed; sites={sorted(sites)}')
    assert any(abs(s - acceptor) <= 2 for s in sites), (
        f'canonical acceptor at {acceptor} not indexed; sites={sorted(sites)}')


def test_noncanonical_intron_edges_are_NOT_indexed():
    """The structural gap: an AT..AA junction has no index entry.

    This is the mechanism behind planning/721's 34.5 % non-canonical recovery.
    If this test starts failing, the index has been extended — re-measure the
    recovery figure and update the docstrings in overhang_resolver.py and
    align_command.py that quote 99 % / 35 %.
    """
    idx = SpliceSiteIndex.build({'chrN': NONCANON})
    donor = len(_EXON_L)
    acceptor = len(_EXON_L) + 2 + len(_INTRON_BODY) + 2
    sites = _all_sites(idx, 'chrN', 0, len(NONCANON))
    assert not any(abs(s - donor) <= 1 for s in sites), (
        'a non-canonical donor is now indexed — planning/721 recovery figures '
        'are stale, re-measure')
    assert not any(abs(s - acceptor) <= 1 for s in sites), (
        'a non-canonical acceptor is now indexed — planning/721 recovery '
        'figures are stale, re-measure')


def test_index_kinds_are_pinned():
    """The site classes themselves encode the (default) limitation.

    2026-08-17 update (planning/722b): this pin fired, as designed, when the
    Prp18-class acceptor extension landed. The index now ALWAYS carries two
    extra acceptor arrays (acc_plus_ext / acc_minus_ext: TG/CG/GG/AT
    transcript-frame acceptors, Roy et al. 2023 NAR gkad968), but querying
    them is OPT-IN via the *_all union kinds + ResolverConfig
    acceptor_classes='prp18'. The DEFAULT query path is still canonical-only,
    so the planning/721 figures (99.1 % annotated / 34.5 % non-canonical)
    continue to describe default behavior. Donors were NOT extended.

    Any further vocabulary growth fails here again — that is the signal to
    re-price (candidate density, alpha) and re-measure before shipping.
    """
    idx = SpliceSiteIndex.build({'chrC': CANON})
    kinds = {k.split('|', 1)[1] for k in idx._arrays}
    # 2026-09-04: fired again, as designed, when the paired AT-AC class landed
    # (don_at_* / acc_ac_*, format v3; ResolverConfig.atac, opt-in). Spelled
    # out literally so that adding a kind to _KINDS cannot satisfy this pin
    # silently — the point is to force a re-price on every vocabulary change.
    assert kinds == {'don_gt_plus', 'don_gc_plus', 'acc_plus', 'don_minus', 'acc_minus',
                     'acc_plus_ext', 'acc_minus_ext',
                     'don_at_plus', 'acc_ac_plus', 'don_at_minus', 'acc_ac_minus'}, (
        f'splice-site kind vocabulary changed ({kinds}) — re-price and '
        f're-measure; see planning/721 + 722b')


def test_default_query_path_is_still_canonical_only():
    """The 721 figures describe DEFAULT behavior; pin that the default
    resolver kind selection never touches the extended arrays."""
    from rectify.core.align.overhang_resolver import (
        _boundary_kinds, _site_kinds,
    )
    canon = set(_KINDS) | {'don_plus'}
    for side in ('L', 'R'):
        for strand in ('+', '-'):
            assert set(_site_kinds(side, strand)) <= canon
    for strand in ('+', '-'):
        assert set(_boundary_kinds(strand)) <= canon


def test_docstring_records_the_gap():
    """The limitation must stay documented where a caller will see it.

    Cheap guard against a future docstring rewrite quietly restoring the
    "replaces mapPacBio" framing that planning/721 refuted.
    """
    from rectify.core.align import overhang_resolver
    doc = overhang_resolver.__doc__ or ''
    assert 'SCOPE LIMIT' in doc, 'the planning/721 scope limit left the docstring'
    assert '34.5' in doc, 'the measured non-canonical recovery figure is missing'
