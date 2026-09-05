"""A D abutting an N is not an alternative junction (the D/N-merge twin).

`junction_scoring._merge_del_into_intron` folds a run of adjacent D/N ops
containing an N into a single N — RECTIFY's convention that `111N 3D` and `114N`
are the same intron.  The junction POOL is built through that normalization, but
`_iter_n_ops` reports the RAW N-op, so a read with `111N 3D` contributes a pool
junction 3 bp longer than the one it is scored at.

The two forms are the same alignment: identical reference span, identical query
bases, only the op labels differ.  So the scorer must not choose between them —
and left to itself it always picks the merged form, because the abutting deletion
stops costing anything once it is relabelled as intron.  On the unselected
hold-out that fake 1.128-unit "improvement" cleared even the >= 1.0
annotated-canonical evidence gate and turned an annotated GT-AG into an
off-annotation GT-TC.

Annotation decides instead — and both directions occur in real data:
    05e2f8d8   annotated 111N -> unannotated 114N   must NOT move
    3178286c   unannotated 109N -> annotated 114N   SHOULD move

Author: Kevin R. Roy (agent S1)
"""

import sys
from pathlib import Path

import pysam
import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.splice import junction_refiner as jr  # noqa: E402
from rectify.core.splice.junction_scoring import _merge_del_into_intron  # noqa: E402

CHROM = "chrT"
GLEN = 600
REF_START = 150
# 50M 100N 3D 47M — exon [150,200), intron [200,300), deletion [300,303)
CIGAR = [(0, 50), (3, 100), (2, 3), (0, 47)]
N_IDX = 1

RAW = (200, 300)
TWIN = (200, 303)
OTHER = (200, 312)      # a genuine alternative, not a twin


def _genome():
    g = list("C" * GLEN)
    g[200:202] = list("GT")
    for e in (300, 303, 312):
        g[e - 2:e] = list("AG")
    return "".join(g)


GENOME = _genome()


def _read(cigar=CIGAR):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": GLEN}]}
    )
    r = pysam.AlignedSegment(header)
    r.query_name = "twin_read"
    r.reference_id = 0
    r.reference_start = REF_START
    r.mapping_quality = 60
    r.cigartuples = cigar
    n = sum(l for op, l in cigar if op in (0, 1, 4, 7, 8))
    r.query_sequence = ("ACGT" * (n // 4 + 1))[:n]
    r.query_qualities = pysam.qualitystring_to_array("I" * n)
    return r


def _winner(monkeypatch, scores, annotated=frozenset()):
    def fake_score(query, q_split, js, je, genome_seq, **kw):
        return scores[(js, je)], 0

    monkeypatch.setattr(jr, "_score_junction", fake_score)
    idx = jr._build_junction_index({(CHROM, s, e) for s, e in scores})
    repl = jr.refine_read_junctions(
        _read(), idx, {(CHROM, s, e) for s, e in annotated}, GENOME, "+",
        boundary_error_window=0,
    )
    return (repl[0][3], repl[0][4]) if repl else None


# ---------------------------------------------------------------------------
# _dn_run_extent
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cigar,n_idx,ns,ne,expected", [
    ([(0, 50), (3, 100), (2, 3), (0, 47)], 1, 200, 300, (200, 303)),   # N then D
    ([(0, 50), (2, 3), (3, 100), (0, 47)], 2, 203, 303, (200, 303)),   # D then N
    ([(0, 50), (2, 3), (3, 100), (2, 2), (0, 47)], 2, 203, 303, (200, 305)),
    ([(0, 50), (3, 100), (0, 47)], 1, 200, 300, (200, 300)),           # no D at all
    ([(0, 50), (3, 100), (2, 4), (3, 60), (0, 47)], 1, 200, 300, (200, 364)),
])
def test_dn_run_extent(cigar, n_idx, ns, ne, expected):
    assert jr._dn_run_extent(cigar, n_idx, ns, ne) == expected


def test_dn_run_extent_agrees_with_the_pool_normalization():
    """The extent must be exactly what `_merge_del_into_intron` would produce."""
    merged = _merge_del_into_intron(CIGAR)
    pos, pool_junctions = REF_START, []
    for op, ln in merged:
        if op == 3:
            pool_junctions.append((pos, pos + ln))
        if op in (0, 2, 3, 7, 8):
            pos += ln
    assert pool_junctions == [jr._dn_run_extent(CIGAR, N_IDX, *RAW)]


# ---------------------------------------------------------------------------
# The policy
# ---------------------------------------------------------------------------

def test_an_annotated_read_is_not_relabelled_onto_an_unannotated_twin(monkeypatch):
    """The 05e2f8d8 case — and the twin's score advantage is large on purpose."""
    assert _winner(monkeypatch, {RAW: 4.7, TWIN: 3.5}, annotated={RAW}) is None


def test_score_alone_never_moves_a_read_to_its_twin(monkeypatch):
    """Neither form annotated: keep the aligner's own, whatever the score says."""
    assert _winner(monkeypatch, {RAW: 4.7, TWIN: 0.0}) is None


def test_an_unannotated_read_is_relabelled_onto_an_annotated_twin(monkeypatch):
    """The 3178286c case — this one is a real correction."""
    assert _winner(monkeypatch, {RAW: 4.7, TWIN: 3.5}, annotated={TWIN}) == TWIN


def test_both_annotated_keeps_the_reads_own_form(monkeypatch):
    assert _winner(monkeypatch, {RAW: 4.7, TWIN: 3.5},
                   annotated={RAW, TWIN}) is None


def test_a_genuine_alternative_is_still_reachable(monkeypatch):
    """Only the twin is special: real candidates are scored as before."""
    assert _winner(monkeypatch, {RAW: 4.7, TWIN: 0.0, OTHER: 0.5},
                   annotated={OTHER}) == OTHER


def test_a_read_without_an_abutting_deletion_is_unaffected(monkeypatch):
    """No D means no twin — (200,303) is then an ordinary candidate."""
    def fake_score(query, q_split, js, je, genome_seq, **kw):
        return {RAW: 4.7, TWIN: 0.0}[(js, je)], 0

    monkeypatch.setattr(jr, "_score_junction", fake_score)
    idx = jr._build_junction_index({(CHROM, s, e) for s, e in (RAW, TWIN)})
    repl = jr.refine_read_junctions(
        _read([(0, 50), (3, 100), (0, 50)]), idx, set(), GENOME, "+",
        boundary_error_window=0,
    )
    assert [(r[3], r[4]) for r in repl] == [TWIN]
