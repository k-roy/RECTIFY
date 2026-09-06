"""ISSUE-021 — Module 2H scored every N-op of a soft-clipped read on a window
shifted ``len(leading S)`` bases into exon 1.

``_iter_n_ops`` reports ``q_split`` counted from the first ALIGNED query base (a
leading soft clip is excluded by design, and that contract is pinned in
``tests/test_junction_multi_intron_coords.py``).  ``refine_read_junctions`` indexes
``read.query_sequence`` — which INCLUDES the soft-clipped bases — so for a read
``5S10M100N10M`` (query ``AAAAA CCCCCCCCCC GGGGGGGGGG``) the scorer received
``q_split == 10`` and read its 30-nt rescue window from ``query_sequence[10:]`` ==
``CCCCCGGGGGGGGGG``: five exon-1 bases first, the true exon-2 start at index 15.
DRS reads are almost always 5' soft-clipped (67 of the 74 to_noncanon rows carry a
leading S; 71–73S clips put the whole window inside the clip), so this affected
most of the N-ops 2H scored on the Sumner cohort.

The fix adds ``read.query_alignment_start`` (= the leading-S length) to ``q_split``
once at the loop head of ``refine_read_junctions``.  These tests pin:

* the window ``_score_junction`` receives starts at the true exon-2 base, for a
  plus-strand read, a minus-strand mirror, and a hard+soft (``3H5S``) prefix;
* the same for ``_positional_signal`` (the drift-gate path reads the same window);
* an UNCLIPPED read hands ``_score_junction`` byte-identical inputs — the same
  ``query_sequence`` object and the raw ``q_split`` — and gets the same output;
* a leading soft clip no longer changes 2H's decision on the real scorer.

Author: Kevin R. Roy (agent, ISSUE-021)
"""

import random
import sys
from pathlib import Path

import pysam
import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.splice import junction_refiner as jr  # noqa: E402

CHROM = "chrT"
GLEN = 600
REF_START = 150

# 5S10M100N10M at ref 150: exon1 [150,160), intron [160,260), exon2 [260,270)
CIGAR_CLIPPED = [(4, 5), (0, 10), (3, 100), (0, 10)]
CIGAR_HS = [(5, 3), (4, 5), (0, 10), (3, 100), (0, 10)]
CIGAR_UNCLIPPED = [(0, 10), (3, 100), (0, 10)]
QUERY = "AAAAA" + "C" * 10 + "G" * 10
EXON2 = "G" * 10
INCUMBENT = (160, 260)
ALTERNATIVE = (160, 263)


def _genome():
    """GT..AG at the incumbent and an AG three bases on for the alternative."""
    g = list("C" * GLEN)
    g[160:162] = list("GT")
    g[258:260] = list("AG")
    g[261:263] = list("AG")
    return "".join(g)


def _read(cigar, query, is_reverse=False):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": GLEN}]}
    )
    r = pysam.AlignedSegment(header)
    r.query_name = "softclip_read"
    r.reference_id = 0
    r.reference_start = REF_START
    r.mapping_quality = 60
    r.is_reverse = is_reverse
    r.cigartuples = cigar
    r.query_sequence = query
    r.query_qualities = pysam.qualitystring_to_array("I" * len(query))
    return r


def _capture_scorer(monkeypatch, read, strand, scores, **kw):
    """Run refine_read_junctions with a fake _score_junction that records what it is
    handed; return the list of (query, q_split, js, je) calls."""
    calls = []

    def fake_score(query, q_split, js, je, genome_seq, **_):
        calls.append((query, q_split, js, je))
        return scores[(js, je)], 0

    monkeypatch.setattr(jr, "_score_junction", fake_score)
    idx = jr._build_junction_index({(CHROM, s, e) for s, e in scores})
    jr.refine_read_junctions(
        read, idx, {(CHROM,) + ALTERNATIVE}, _genome(), strand,
        boundary_error_window=0, **kw,
    )
    assert calls, "the N-op was not scored (pool/filter problem in the test itself)"
    return calls


# ---------------------------------------------------------------------------
# The window handed to _score_junction starts at the true exon-2 base
# ---------------------------------------------------------------------------

def test_iter_n_ops_still_excludes_the_leading_clip():
    """The producer's contract is unchanged: q_split counts aligned bases only."""
    (_, ns, ne, q_split), = list(jr._iter_n_ops(_read(CIGAR_CLIPPED, QUERY)))
    assert (ns, ne) == INCUMBENT
    assert q_split == 10
    # ...which is exactly why indexing query_sequence with it was wrong:
    assert QUERY[q_split:] == "CCCCC" + EXON2


@pytest.mark.parametrize("cigar,clip", [(CIGAR_CLIPPED, 5), (CIGAR_HS, 5)])
def test_plus_strand_window_starts_at_exon2(monkeypatch, cigar, clip):
    read = _read(cigar, QUERY)
    assert read.query_alignment_start == clip
    calls = _capture_scorer(monkeypatch, read, "+", {INCUMBENT: 1.0, ALTERNATIVE: 1.0})
    for query, q_split, _js, _je in calls:
        assert query is read.query_sequence or query == read.query_sequence
        assert q_split == clip + 10
        assert query[q_split:q_split + 30] == EXON2
        assert query[q_split:] == read.query_alignment_sequence[10:]   # the other formulation


def test_minus_strand_mirror_window_starts_at_exon2(monkeypatch):
    """The leading S in CIGAR order is the 3' RNA end of a minus-strand read; the
    shift is the same in CIGAR/genomic orientation, and so is the fix."""
    read = _read(CIGAR_CLIPPED, QUERY, is_reverse=True)
    calls = _capture_scorer(monkeypatch, read, "-", {INCUMBENT: 1.0, ALTERNATIVE: 1.0})
    for query, q_split, _js, _je in calls:
        assert q_split == 15
        assert query[q_split:q_split + 30] == EXON2


def test_positional_signal_reads_the_same_window(monkeypatch):
    """The drift-gate spare path (_positional_signal) indexes q with the same q_split."""
    seen = []
    real = jr._positional_signal

    def spy(genome_seq, q, q_split, ne, new_je, **kw):
        seen.append(q[q_split:q_split + 10])
        return real(genome_seq, q, q_split, ne, new_je, **kw)

    monkeypatch.setattr(jr, "_positional_signal", spy)
    read = _read(CIGAR_CLIPPED, QUERY)
    # alternative wins by 1.0 < hold_margin 2.0 -> veto path -> positional gate consulted
    _capture_scorer(monkeypatch, read, "+", {INCUMBENT: 3.0, ALTERNATIVE: 2.0},
                    hold_margin=2.0, drift_positional_gate=1.0)
    assert seen == [EXON2]


# ---------------------------------------------------------------------------
# Unclipped reads: byte-identical inputs and outputs
# ---------------------------------------------------------------------------

def test_unclipped_read_scorer_inputs_are_unchanged(monkeypatch):
    read = _read(CIGAR_UNCLIPPED, QUERY[5:])
    assert read.query_alignment_start == 0
    (_, _, _, raw_q_split), = list(jr._iter_n_ops(read))
    calls = _capture_scorer(monkeypatch, read, "+", {INCUMBENT: 1.0, ALTERNATIVE: 1.0})
    for query, q_split, _js, _je in calls:
        assert query == read.query_sequence
        assert q_split == raw_q_split == 10


def test_unclipped_read_real_scorer_output_is_unchanged(monkeypatch):
    """Wrap the REAL scorer: every call must receive exactly (query_sequence, raw
    q_split) and return exactly what the real scorer returns for those inputs."""
    read = _read(CIGAR_UNCLIPPED, QUERY[5:])
    (_, _, _, raw_q_split), = list(jr._iter_n_ops(read))
    real = jr._score_junction
    genome = _genome()
    seen = []

    def wrapped(query, q_split, js, je, genome_seq, **kw):
        out = real(query, q_split, js, je, genome_seq, **kw)
        expected = real(read.query_sequence, raw_q_split, js, je, genome, **kw)
        seen.append((query == read.query_sequence, q_split == raw_q_split, out == expected))
        return out

    monkeypatch.setattr(jr, "_score_junction", wrapped)
    idx = jr._build_junction_index({(CHROM,) + INCUMBENT, (CHROM,) + ALTERNATIVE})
    jr.refine_read_junctions(read, idx, {(CHROM,) + ALTERNATIVE}, genome, "+",
                             boundary_error_window=0)
    assert seen and all(all(t) for t in seen)


# ---------------------------------------------------------------------------
# Behavior: a leading soft clip must not change 2H's decision (real scorer)
# ---------------------------------------------------------------------------

def _decision_genome():
    """Seeded random genome; GT donor at 200, a CANONICAL but unannotated incumbent
    acceptor at 300 (CAG) and an annotated canonical alternative at 303 (CAG).
    Both are canonical-class, so no prior separates them: at equal scores the
    incumbent is kept (is_alt), and only the read's exon-2 bases can move it.
    Pre-fix, a 71-nt clip put the whole window inside the clip -> tie -> no move."""
    rng = random.Random(21)
    g = [rng.choice("ACGT") for _ in range(GLEN)]
    g[200:202] = list("GT")
    g[297:300] = list("CAG")
    g[300:303] = list("CAG")
    g[303:306] = list("TTT")      # exon-2 start unlike the CAG the aligner stretched over
    return "".join(g)


def _decision_read(genome, clip):
    """exon1 = genome[150:200]; exon2 = genome[303:343] but ALIGNED from 300 (40M):
    the read's exon-2 bases belong to the alternative acceptor.  Optionally
    prefixed by `clip` junk bases as a leading soft clip."""
    exon1, exon2 = genome[150:200], genome[303:343]
    junk = ("TACG" * 20)[:clip]
    cigar = ([(4, clip)] if clip else []) + [(0, 50), (3, 100), (0, 40)]
    return _read(cigar, junk + exon1 + exon2)


@pytest.mark.parametrize("clip", [5, 20, 71])
def test_leading_clip_does_not_change_the_decision(clip):
    genome = _decision_genome()
    inc, alt = (200, 300), (200, 303)
    idx = jr._build_junction_index({(CHROM,) + inc, (CHROM,) + alt})
    annotated = {(CHROM,) + alt}

    def decide(read):
        repl = jr.refine_read_junctions(read, idx, annotated, genome, "+",
                                        boundary_error_window=0)
        return [(s, e, ns, ne) for (_, s, e, ns, ne) in repl]

    unclipped = decide(_decision_read(genome, 0))
    assert unclipped == [(200, 300, 200, 303)], unclipped   # the read's bases pick the alternative
    assert decide(_decision_read(genome, clip)) == unclipped
