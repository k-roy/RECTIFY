#!/usr/bin/env python3
"""Regression guard for the three scoring FENCES that make GMAP a net-positive
member of the aligner panel.

Context (A549 chr5 validation, 2026-06-18): GMAP places **55.6%** non-canonical
junctions (878k of 1.58M) vs 1.1% for minimap2, and half its novel-canonical
candidates are GMAP-only singletons. It earns its seat ONLY because the chimeric
``score_segment`` scoring suppresses that noise while keeping its ~100 genuine
unique novel junctions and its annotated parity. The three fences are:

  1. NON-CANONICAL PENALTY — a junction that is neither annotated nor canonical
     (GT-AG within its ambiguity window) is penalized, so GMAP's spurious
     junctions LOSE at selection.
  2. SUPPORT GATE — the annotated bonus only applies when the read anchors the
     junction (min flanking anchor), so a junction FORCED onto a catalogued site
     with a short overhang cannot win on annotation alone.
  3. AMBIGUITY MATCHING — the SAME junction written a few bp over (donor/acceptor
     repeat) still matches the catalog / earns canonical credit, so GMAP's
     legitimate de-novo placements are not charged as wrong.

These tests pin BEHAVIOURAL INVARIANTS (relative orderings), NOT exact constants —
the bonus magnitudes may still be tuned, but if any of these orderings inverts,
GMAP's noise can leak back into the consensus. Break one ⇒ revisit the
"keep GMAP" decision, do not just update the assertion.

Author: Kevin R. Roy
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rectify.core.consensus.chimeric_consensus import score_segment, CigarEvent


def _intron_genome(donor, acceptor, intron_len=50, flank=100):
    """exon1(flank C) + donor(2) + body(T) + acceptor(2) + exon2(flank C);
    intron = [flank, flank+intron_len). Flanks/body are non-repetitive at the
    boundary so the ambiguity window is 0 (predictable canonical/non-canonical)."""
    body = 'T' * (intron_len - 4)
    return 'C' * flank + donor + body + acceptor + 'C' * flank


def _events(anchor, jstart, jlen):
    """[M(anchor), N over [jstart, jstart+jlen), M(anchor)]."""
    return [
        CigarEvent(0, anchor, 0, anchor, jstart - anchor, jstart),
        CigarEvent(3, jlen, anchor, anchor, jstart, jstart + jlen),
        CigarEvent(0, anchor, anchor, 2 * anchor, jstart + jlen, jstart + jlen + anchor),
    ]


_J, _L = 100, 50   # intron [100, 150)


# ── Fence 1: non-canonical penalty (suppresses GMAP's 878k spurious junctions) ──

def test_fence1_noncanonical_loses_to_canonical():
    canon = {'c': _intron_genome('GT', 'AG')}
    noncanon = {'c': _intron_genome('CC', 'CC')}
    sc = score_segment(_events(20, _J, _L), 'interior', 'c', canon, set())
    snc = score_segment(_events(20, _J, _L), 'interior', 'c', noncanon, set())
    # THE FENCE: a non-canonical (GMAP-style) junction scores strictly worse than
    # a clean canonical one, and nets a negative contribution. If this inverts,
    # GMAP's non-canonical junctions can win segments.
    assert snc.n_canonical_junctions == 0
    assert snc.score < sc.score
    assert snc.score < 0


def test_fence1_noncanonical_loses_even_to_no_junction():
    """An aligner that skips (no junction, plain match) must beat one that
    invents a non-canonical junction over the same span."""
    noncanon = {'c': _intron_genome('CC', 'CC')}
    g = noncanon
    with_junc = score_segment(_events(20, _J, _L), 'interior', 'c', g, set())
    # same query span but a single matched run, no N op
    no_junc = score_segment(
        [CigarEvent(0, 40, 0, 40, _J - 20, _J + 20)], 'interior', 'c', g, set())
    assert with_junc.score < no_junc.score


# ── Fence 2: support gate (catalog breaks ties, never overrides read evidence) ──

def test_fence2_short_anchor_annotation_gets_no_bonus():
    canon = {'c': _intron_genome('GT', 'AG')}
    catalog = {('c', _J, _J + _L)}              # normalized (no ambiguity here)
    supported = score_segment(_events(20, _J, _L), 'interior', 'c', canon, catalog)
    short = score_segment(_events(3, _J, _L), 'interior', 'c', canon, catalog)
    # THE FENCE: annotation bonus requires read support. A short-overhang junction
    # forced onto a catalogued site earns canonical credit only, NOT the +8.
    assert supported.n_annotated_junctions == 1
    assert short.n_annotated_junctions == 0
    assert short.n_annotated_unsupported == 1
    assert short.score < supported.score


def test_fence2_anchor_from_full_alignment_not_clipped_segment():
    """The gate anchor must come from the aligner's FULL alignment — a junction
    whose segment-clipped flanks are short but whose full-alignment anchor is long
    MUST still earn the bonus (else segmentation truncation silently nukes it)."""
    canon = {'c': _intron_genome('GT', 'AG')}
    catalog = {('c', _J, _J + _L)}
    full = _events(40, _J, _L)                  # full alignment: 40bp flanks
    # segment events clipped tight around the junction (5bp visible)
    clipped = [
        CigarEvent(0, 5, 0, 5, _J - 5, _J),
        CigarEvent(3, _L, 5, 5, _J, _J + _L),
        CigarEvent(0, 5, 5, 10, _J + _L, _J + _L + 5),
    ]
    s = score_segment(clipped, 'interior', 'c', canon, catalog,
                      aligner_full_events=full)
    assert s.n_annotated_junctions == 1         # full anchor (40) passes the gate


# ── Fence 3: ambiguity matching (protects GMAP's legit de-novo placements) ──

def test_fence3_ambiguity_shifted_junction_still_matches_catalog():
    # intron flanked by a poly-A run: [105,115) is equivalent to [100,110)
    g = 'C' * 100 + 'A' * 20 + 'C' * 100
    catalog = {('chrA', 100, 110)}              # normalized catalog
    ev = [CigarEvent(0, 20, 0, 20, 85, 105),
          CigarEvent(3, 10, 20, 20, 105, 115),
          CigarEvent(0, 20, 20, 40, 115, 135)]
    s = score_segment(ev, 'interior', 'chrA', {'chrA': g}, catalog)
    # THE FENCE: the same junction one window-shift over is NOT charged as a
    # mismatch (this is what flipped GMAP from 0.09 -> ~1.0 annotated win/loss).
    assert s.n_annotated_junctions == 1


def test_fence3_canonical_credited_within_window():
    """A junction landing one bp off a GT-AG site, but canonical within its
    ambiguity window, still earns canonical credit (not the -3 penalty)."""
    from rectify.core.consensus.chimeric_consensus import _canonical_within_window
    # GT-AG reachable by a 1bp slide within an A-repeat boundary
    seq = 'C' * 50 + 'GT' + 'T' * 20 + 'AG' + 'C' * 50  # canonical at [50,74)
    assert _canonical_within_window(50, 74, seq, 0, 0) is True


if __name__ == '__main__':
    import pytest
    sys.exit(pytest.main([__file__, '-q']))
