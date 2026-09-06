"""ISSUE-018: Module 2H must not move a canonical junction onto an unannotated
non-canonical site by planting a compensating indel.

Measured on the 16-library Sumner cohort (269ebe7): of 74 ``to_noncanon`` moves,
47 went OFF an annotated/canonical junction ONTO a non-canonical site, and every
one carried the same signature -- the surgery folded a stock ``D`` abutting the
``N`` into the intron (``72M1D1307N2D194M`` -> ``72M1D1309N194M``) or planted a
compensating ``I`` beside the new ``N`` (``1077N4M1I112M`` -> ``1080N3I1M1I112M``).
The read's exon SEQUENCE never moved; only the coordinate did.  The deletion
cost that vanishes when the D is relabeled as intron (1.07-1.41 edit-distance
units with the human tables) is what cleared the annotated-canonical hold, and
the D/N-twin rule never saw the move because the pool is built through
``_merge_del_into_intron``: read A's merged form is read B's "alternative".

Arbiter RULING 3 (b): a junction whose STOCK placement is canonical-class
(``_canonical_tier <= _CANONICAL_TIER_MAX``) may be moved onto an UNANNOTATED
non-canonical-class site only on a preponderance of evidence -- no annotated or
canonical alternative within delta, no I/D abutting the new N on either side, no
rise in indel burden, both anchors at the informative floor.  The gate reads the
post-move anchors from a dry run of the real surgery
(``junction_refiner._noncanon_destination_refusal``).  Moves whose destination
is annotated (the 27 class-B rows) or canonical-class (every drift fix) are not
touched.

Each test builds a hermetic locus: a seeded pseudo-random 600-nt genome (so the
anchors carry real information), a canonical GT..CAG stock intron, and a
destination whose motif is non-canonical.  ``_score_junction`` is faked so the
destination wins by 1.5 (clearing the 1.0 R1 hold) -- the gate is what must
decide.  Pattern: tests/test_canonical_tier_line.py.

Author: Kevin R. Roy (agent rectify-ce, ISSUE-018)
"""

import random
import sys
from collections import Counter
from pathlib import Path

import pysam
import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.splice import junction_refiner as jr  # noqa: E402
from rectify.core.splice.junction_scoring import _CANONICAL_TIER_MAX  # noqa: E402
from rectify.core.splice.splice_aware_5prime import min_informative_clip_bp  # noqa: E402

CHROM = "chrT"
GLEN = 600
REF_START = 150
SEED = 18            # background genome; fixed so every anchor is the same informative sequence

M, I, D, N = 0, 1, 2, 3     # pysam CIGAR op codes

# The three read shapes.  Every read starts at REF_START; the N-op of each is
# located with _iter_n_ops, never hard-coded.
D_MERGE = [(M, 40), (D, 1), (N, 100), (D, 2), (M, 40)]   # 9691a87c's form: 72M1D1307N2D194M
I_SHAPE = [(M, 40), (N, 100), (M, 4), (I, 1), (M, 36)]   # the 1077N4M1I112M form
CLEAN = [(M, 40), (N, 100), (M, 40)]

HIGH, LOW = 2.0, 0.5   # fake scores: the destination wins by 1.5 (> the 1.0 R1 hold)


def _header():
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": GLEN}]}
    )


def _read(cigar, name="issue018"):
    r = pysam.AlignedSegment(_header())
    r.query_name = name
    r.reference_id = 0
    r.reference_start = REF_START
    r.mapping_quality = 60
    r.cigartuples = cigar
    n = sum(l for op, l in cigar if op in (0, 1, 4, 7, 8))
    r.query_sequence = ("ACGT" * (n // 4 + 1))[:n]
    r.query_qualities = pysam.qualitystring_to_array("I" * n)
    return r


def _locus(cigar, exon2):
    """``(read, genome, cigar_idx, ns, ne)`` for a canonical GT..CAG stock intron
    whose downstream exon starts with *exon2* -- the bases that decide the
    destination's acceptor motif and, for the fast path, whether they mirror the
    donor."""
    read = _read(cigar)
    (cigar_idx, ns, ne, _q_split), = list(jr._iter_n_ops(read))
    rng = random.Random(SEED)
    g = [rng.choice("ACGT") for _ in range(GLEN)]

    def plant(pos, s):
        g[pos:pos + len(s)] = list(s)

    plant(ns - 2, "CC")        # exon-1 tail: no accidental donor one base upstream
    plant(ns, "GT")            # canonical donor
    plant(ns + 2, "CCCC")      # intron interior: no accidental donor 1-3 bases downstream
    plant(ne - 3, "CAG")       # canonical YAG acceptor -> tier 0
    plant(ne, exon2)
    return read, "".join(g), cigar_idx, ns, ne


def _tiers(genome, stock, dest):
    return (jr._canonical_tier(*stock, genome, "+"),
            jr._canonical_tier(*dest, genome, "+"))


def _fake_scores(monkeypatch, scores):
    def fake_score(query, q_split, js, je, genome_seq, **kw):
        return scores[(js, je)], 0
    monkeypatch.setattr(jr, "_score_junction", fake_score)


def _refine(read, genome, scores, annotated, *, window=10, counters=None,
            profile=None, **kw):
    idx = jr._build_junction_index({(CHROM, s, e) for s, e in scores})
    return jr.refine_read_junctions(
        read, idx, {(CHROM, s, e) for s, e in annotated}, genome, "+",
        boundary_error_window=window, counters=counters, profile=profile, **kw,
    )


def _apply(read, replacements, genome):
    """``(applied, cigar_string)`` after the real surgery."""
    modified, applied = jr._apply_replacements_to_read(
        read, replacements, genome, "+", 0.25, 15,
    )
    return applied, modified.cigarstring


# ---------------------------------------------------------------------------
# 1. The D-merge shape (the majority signature of the 47)
# ---------------------------------------------------------------------------

class TestDMergeIntoIntron:
    """``40M1D100N2D40M`` -> destination (ns, ne+2).

    The 2D abutting the N folds into the intron: the surgery yields
    ``40M1D102N40M``.  The read's OWN D/N twin is (ns-1, ne+2), because
    ``_dn_run_extent`` folds the D on both sides, so (ns, ne+2) -- another
    read's merged form in a real pool -- is scored as an ordinary alternative
    and the twin rule does not see it.  That is exactly how 9691a87c moved.
    """

    def setup_method(self):
        self.read, self.g, self.idx, self.ns, self.ne = _locus(D_MERGE, "CAA")
        self.stock = (self.ns, self.ne)
        self.dest = (self.ns, self.ne + 2)
        self.scores = {self.stock: HIGH, self.dest: LOW}

    def test_the_shape_is_the_cohort_s(self):
        ts, td = _tiers(self.g, self.stock, self.dest)
        assert ts <= _CANONICAL_TIER_MAX < td                    # GT-CAG -> GT-CA
        twin = jr._dn_run_extent(self.read.cigartuples, self.idx, *self.stock)
        assert twin == (self.ns - 1, self.ne + 2)
        assert self.dest != twin                                 # so it IS scored
        assert jr._has_boundary_error(self.read.cigartuples, *self.stock, REF_START, 10)

    def test_refused_with_adjacent_indel_and_counted(self, monkeypatch):
        _fake_scores(monkeypatch, self.scores)
        counters = Counter()
        profile = jr.JunctionRefineProfile()
        repl = _refine(self.read, self.g, self.scores, [self.stock],
                       counters=counters, profile=profile)
        assert repl == []
        assert counters["noncanon_destination_refused"] == 1
        assert counters["noncanon_destination_refused_adjacent_indel"] == 1
        reasons = [k for k in counters if k.startswith("noncanon_destination_refused_")]
        assert reasons == ["noncanon_destination_refused_adjacent_indel"]
        assert profile.counts["noncanon_destination_refused"] == 1
        assert profile.counts["noncanon_destination_refused_adjacent_indel"] == 1

    def test_the_token_names_the_folded_deletion(self):
        """Direct: the dry run leaves the stock 1D abutting the new N on the left."""
        token = jr._noncanon_destination_refusal(
            self.read, self.idx, *self.stock, *self.dest, self.g, "+", CHROM,
            {(CHROM, *self.stock)}, 0.25, 15,
        )
        assert token == "adjacent_indel"
        assert token in jr.NONCANON_REFUSALS

    def test_the_gate_off_restores_the_cohort_move(self, monkeypatch):
        """A/B baseline: RECTIFY_2H_NONCANON_PREPONDERANCE=0 is the pre-fix 2H."""
        monkeypatch.setattr(jr, "_NONCANON_PREPONDERANCE", False)
        _fake_scores(monkeypatch, self.scores)
        counters = Counter()
        repl = _refine(self.read, self.g, self.scores, [self.stock], counters=counters)
        assert repl == [(self.idx, *self.stock, *self.dest)]
        assert counters["noncanon_destination_refused"] == 0
        # The coordinate moved; the read's bases did not.
        assert _apply(self.read, repl, self.g) == (True, "40M1D102N40M")

    def test_onto_an_annotated_noncanonical_intron_still_moves(self, monkeypatch):
        """Class B (27 rows, e.g. 746a26bd): the destination is a GENCODE intron
        whose own motif is non-canonical.  Annotation decides; the gate must not
        touch it."""
        _fake_scores(monkeypatch, self.scores)
        counters = Counter()
        repl = _refine(self.read, self.g, self.scores, [self.stock, self.dest],
                       counters=counters)
        assert repl == [(self.idx, *self.stock, *self.dest)]
        assert counters["noncanon_destination_refused"] == 0
        assert _apply(self.read, repl, self.g) == (True, "40M1D102N40M")


# ---------------------------------------------------------------------------
# 2. The compensating-insertion shape
# ---------------------------------------------------------------------------

class TestCompensatingInsertion:
    """``40M100N4M1I36M`` -> acceptor +k.  The general path re-labels the
    boundary by planting I(k) right after the new N (the cohort's
    ``1077N4M1I112M`` -> ``1080N3I1M1I112M``).  A +1 shift is caught before
    the dry run: the canonical stock itself is a proper site within delta of
    the destination.  A +3 shift reaches the dry run and is refused on the
    planted insertion.
    """

    @pytest.mark.parametrize("shift,token,moved_cigar", [
        (1, "alternative_within_delta", "40M101N1I3M1I36M"),
        (3, "adjacent_indel", "40M103N3I1M1I36M"),
    ])
    def test_refused_and_restored_with_the_gate_off(self, monkeypatch, shift, token,
                                                    moved_cigar):
        read, g, idx, ns, ne = _locus(I_SHAPE, "CTTA")
        stock, dest = (ns, ne), (ns, ne + shift)
        ts, td = _tiers(g, stock, dest)
        assert ts <= _CANONICAL_TIER_MAX < td
        scores = {stock: HIGH, dest: LOW}
        _fake_scores(monkeypatch, scores)

        counters = Counter()
        assert _refine(read, g, scores, [stock], counters=counters) == []
        assert counters["noncanon_destination_refused"] == 1
        assert counters[f"noncanon_destination_refused_{token}"] == 1

        monkeypatch.setattr(jr, "_NONCANON_PREPONDERANCE", False)
        repl = _refine(read, g, scores, [stock])
        assert repl == [(idx, *stock, *dest)]
        assert _apply(read, repl, g) == (True, moved_cigar)


# ---------------------------------------------------------------------------
# 3. Canonical destinations are outside the gate
# ---------------------------------------------------------------------------

class TestCanonicalDestinationsAreNotGated:
    """The gate's scope is destination-unannotated AND destination
    non-canonical-class.  A canonical -> canonical move (the drift fixes 2H
    exists for) must go through in every annotation configuration."""

    @pytest.mark.parametrize("annotated", [
        pytest.param((), id="novel_to_novel"),
        pytest.param(("dest",), id="drift_fix_onto_annotated"),
        pytest.param(("stock", "dest"), id="isoform_swap"),
        pytest.param(("stock",), id="annotated_to_novel_canonical"),
    ])
    def test_canonical_to_canonical_still_moves(self, monkeypatch, annotated):
        read, g, idx, ns, ne = _locus(CLEAN, "ACGCAG")      # a second CAG 3 nt into exon 2
        stock, dest = (ns, ne), (ns, ne + 6)
        ts, td = _tiers(g, stock, dest)
        assert ts <= _CANONICAL_TIER_MAX and td <= _CANONICAL_TIER_MAX
        scores = {stock: HIGH, dest: LOW}
        _fake_scores(monkeypatch, scores)
        annot = [{"stock": stock, "dest": dest}[k] for k in annotated]
        counters = Counter()
        repl = _refine(read, g, scores, annot, window=0, counters=counters)
        assert repl == [(idx, *stock, *dest)]
        assert counters["noncanon_destination_refused"] == 0
        assert _apply(read, repl, g) == (True, "40M106N6I34M")


def test_motif_blind_disables_the_gate(monkeypatch):
    """Arm B decides on read evidence alone; the gate, like the R1 hold, stands down."""
    read, g, idx, ns, ne = _locus(D_MERGE, "CAA")
    stock, dest = (ns, ne), (ns, ne + 2)
    scores = {stock: HIGH, dest: LOW}
    _fake_scores(monkeypatch, scores)
    counters = Counter()
    repl = _refine(read, g, scores, [stock], counters=counters, motif_blind=True)
    assert repl == [(idx, *stock, *dest)]
    assert counters["noncanon_destination_refused"] == 0


# ---------------------------------------------------------------------------
# 4. Preponderance passing
# ---------------------------------------------------------------------------

class TestPreponderance:
    """``40M100N40M`` -> (ns+2, ne+2), a pure 2-bp slide whose flipping bases
    are identical at both placements (g[ns:ns+2] == g[ne:ne+2] == "GT").  The
    surgery takes the fast path -- ``42M100N38M`` -- so both anchors abut the
    new N as clean M runs, no proper site sits within delta, and the indel
    burden is unchanged.  That is the preponderance the ruling asks for, and
    the move is allowed.
    """

    def setup_method(self):
        self.read, self.g, self.idx, self.ns, self.ne = _locus(CLEAN, "GTA")
        self.stock = (self.ns, self.ne)
        self.dest = (self.ns + 2, self.ne + 2)
        self.scores = {self.stock: HIGH, self.dest: LOW}

    def test_allowed(self, monkeypatch):
        ts, td = _tiers(self.g, self.stock, self.dest)
        assert ts <= _CANONICAL_TIER_MAX < td                    # GT-CAG -> CC-GGT
        assert self.g[self.ns:self.ns + 2] == self.g[self.ne:self.ne + 2]   # fast-path precondition
        _fake_scores(monkeypatch, self.scores)
        counters = Counter()
        repl = _refine(self.read, self.g, self.scores, [self.stock],
                       window=0, counters=counters)
        assert repl == [(self.idx, *self.stock, *self.dest)]
        assert counters["noncanon_destination_refused"] == 0
        assert _apply(self.read, repl, self.g) == (True, "42M100N38M")
        assert min(42, 38) >= min_informative_clip_bp()          # both anchors clear the floor

    def test_a_proper_site_within_delta_defeats_it(self, monkeypatch):
        """Criterion (iv): an annotated junction one base from the destination
        makes the move a near-miss of a real site, not evidence for a new one."""
        _fake_scores(monkeypatch, self.scores)
        counters = Counter()
        near = (self.ns + 3, self.ne + 2)
        repl = _refine(self.read, self.g, self.scores, [self.stock, near],
                       window=0, counters=counters)
        assert repl == []
        assert counters["noncanon_destination_refused"] == 1
        assert counters["noncanon_destination_refused_alternative_within_delta"] == 1

    def test_a_refusal_by_the_surgery_itself_is_not_this_gate_s(self, monkeypatch):
        """When the flipping bases differ, the general path must plant
        D(2)..N..I(2), and the surgery's own both-boundary indel invariant
        refuses that.  The gate does not count it -- the counter is strictly
        this gate's refusals -- and the read stays put."""
        read, g, idx, ns, ne = _locus(CLEAN, "CTA")             # g[ne:ne+2] != "GT"
        stock, dest = (ns, ne), (ns + 2, ne + 2)
        ts, td = _tiers(g, stock, dest)
        assert ts <= _CANONICAL_TIER_MAX < td
        scores = {stock: HIGH, dest: LOW}
        _fake_scores(monkeypatch, scores)
        counters = Counter()
        repl = _refine(read, g, scores, [stock], window=0, counters=counters)
        assert repl == [(idx, *stock, *dest)]
        assert counters["noncanon_destination_refused"] == 0
        assert _apply(read, repl, g) == (False, "40M100N40M")


# ---------------------------------------------------------------------------
# 5. The counter reaches the run's stats line
# ---------------------------------------------------------------------------

class TestCounterPlumbing:
    """The refusal must reach the stats dict whether the driver is sequential
    or the parallel worker -- the READY-note prediction (~34-40 refusals on
    the cohort) is read from the ``noncanon_destination_refused`` log field."""

    @staticmethod
    def _bam(tmp_path, read):
        path = tmp_path / "in.bam"
        with pysam.AlignmentFile(str(path), "wb", header=_header()) as fh:
            fh.write(read)
        return path

    @pytest.mark.parametrize("gate,refused,refined,cigar", [
        pytest.param(True, 1, 0, "40M1D100N2D40M", id="gate_on"),
        pytest.param(False, 0, 1, "40M1D102N40M", id="gate_off"),
    ])
    def test_sequential_driver_stats(self, monkeypatch, tmp_path, gate, refused,
                                     refined, cigar):
        monkeypatch.setattr(jr, "_NONCANON_PREPONDERANCE", gate)
        read, g, idx, ns, ne = _locus(D_MERGE, "CAA")
        stock, dest = (ns, ne), (ns, ne + 2)
        _fake_scores(monkeypatch, {stock: HIGH, dest: LOW})
        out = tmp_path / "out.bam"
        stats = jr.refine_bam_junctions(
            str(self._bam(tmp_path, read)), str(out),
            aligner_bams=[], annotated_junctions=set(),
            genome={CHROM: g},
            prebuilt_junction_pool={(CHROM, *stock), (CHROM, *dest)},
            prebuilt_annotated_set={(CHROM, *stock)},
            boundary_error_window=10, sort_and_index=False, n_workers=1,
        )
        assert stats["n_op_reads"] == 1
        assert stats["noncanon_destination_refused"] == refused     # key always present
        assert stats.get("noncanon_destination_refused_adjacent_indel", 0) == refused
        assert stats["refined"] == refined
        assert stats["unchanged"] == 1 - refined
        with pysam.AlignmentFile(str(out), "rb") as fh:
            assert [r.cigarstring for r in fh] == [cigar]

    def test_worker_returns_a_dict_with_counters(self, monkeypatch):
        read, g, idx, ns, ne = _locus(D_MERGE, "CAA")
        stock, dest = (ns, ne), (ns, ne + 2)
        _fake_scores(monkeypatch, {stock: HIGH, dest: LOW})
        jr._WORKER_POOL_STATE.clear()
        jr._WORKER_POOL_STATE.update({
            "header": _header(),
            "junctions_idx": jr._build_junction_index({(CHROM, *stock), (CHROM, *dest)}),
            "annotated_set": {(CHROM, *stock)},
            "genome": {CHROM: g},
            "kwargs": {"boundary_error_window": 10},
        })
        sam = read.to_string()
        try:
            out = jr._refine_read_batch([sam])
        finally:
            jr._WORKER_POOL_STATE.clear()
        assert isinstance(out, dict)
        assert out["results"] == [(sam, [])]
        assert out["profile"] is None
        assert out["counters"] == {
            "noncanon_destination_refused": 1,
            "noncanon_destination_refused_adjacent_indel": 1,
        }
