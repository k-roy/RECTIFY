"""T0 (chrX+chr1) regressions of the 2H branch at 4f9102d — the four reads the
tester's events tier lost against the 6485226 baseline, replayed and fixed.

Bisected with ``dev/todo_run_20260905/issue018_replay/replay_t0_bundle.py``: every
one flips at ee811b6 (ISSUE-021, the corrected rescue window), none at 64bed7a
(the preponderance gate alone).  None of the four incumbents is in the production
candidate set, so the read must go to the best of the OTHER candidates.

F1 — realizability.  a138b2da (SMA_7.12): with the corrected window an annotated
acceptor 19 nt DOWNSTREAM scored 0.0 (the scorer's free-k prefix matched the
read's own exon-2 bases as "intron tail" — the incumbent alignment relabeled),
the surgery refused that both-boundary relabel as an indel-burden rise, and the
annotated 3-nt donor fix ranked second was discarded with it: the read froze at
the drifted stock.  The ranking is now over candidates the surgery can WRITE
(``_realizable`` dry run); a refused head is skipped and counted
(``unrealizable_winner_skipped``).

F2 — annotated-alternative rule.  527c329f (WT_4.2; canonical incumbent): a GT-GA
+1 beat the annotated GT-AG by 0.125, the gate refused it and the read FROZE
12 nt off.  7628e0dd (GSB_191; GT-AC incumbent): GT-CA needing 3D beat the
annotated GT-AG needing 5D.  A non-canonical unannotated winner now loses to an
annotated canonical candidate that shares one boundary and lies within
``_ANNOT_ALT_DELTA`` on the other (``noncanon_dest_lost_to_annotated_alt``).

F3 — the 2H decision counters reach the stats TSV: ``ProcessingStats.module_2h_counters``
(the ``refine_bam_junctions`` stats dict) is written as ``module_2h_<key>`` rows
right after ``module_2h_failed``, survives the parallel-driver JSON round trip
and merge, and ``correct`` fills it in.

The hermetic tests fake ``_score_junction`` (pattern:
tests/test_2h_noncanon_preponderance.py); the replay tests run the tester's
bundle through the real scorer with the bundled human DRS tables and skip when
the bundle is not on this machine.

Author: Kevin R. Roy (agent rectify-ce, T0 4f9102d triage)
"""

import csv
import random
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam
import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.splice import junction_refiner as jr  # noqa: E402
from rectify.core.splice import junction_scoring as js_mod  # noqa: E402
from rectify.core.splice.junction_scoring import _CANONICAL_TIER_MAX  # noqa: E402
from rectify.utils.genome import standardize_chrom_name  # noqa: E402

CHROM = "chrT"
GLEN = 600
REF_START = 150
SEED = 18
M, I, D, N = 0, 1, 2, 3
CLEAN = [(M, 40), (N, 100), (M, 40)]


# ---------------------------------------------------------------------------
# Hermetic locus (as in test_2h_noncanon_preponderance.py)
# ---------------------------------------------------------------------------

def _header(chrom=CHROM, length=GLEN):
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": length}]}
    )


def _read(cigar, name="t0"):
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


def _locus(cigar, exon2="CTA", plants=()):
    """Canonical GT..CAG stock intron on a seeded genome; *plants* are extra
    ``(pos_relative_to_ne, bases)`` writes (e.g. a second acceptor)."""
    read = _read(cigar)
    (cigar_idx, ns, ne, _q), = list(jr._iter_n_ops(read))
    rng = random.Random(SEED)
    g = [rng.choice("ACGT") for _ in range(GLEN)]

    def plant(pos, s):
        g[pos:pos + len(s)] = list(s)

    plant(ns - 2, "CC")
    plant(ns, "GT")
    plant(ns + 2, "CCCC")
    plant(ne - 3, "CAG")
    plant(ne, exon2)
    for rel, s in plants:
        plant(ne + rel, s)
    return read, "".join(g), cigar_idx, ns, ne


def _fake_scores(monkeypatch, scores):
    def fake_score(query, q_split, js, je, genome_seq, **kw):
        return scores[(js, je)], 0
    monkeypatch.setattr(jr, "_score_junction", fake_score)


def _refine(read, genome, scores, annotated, *, counters=None, **kw):
    idx = jr._build_junction_index({(CHROM, s, e) for s, e in scores})
    return jr.refine_read_junctions(
        read, idx, {(CHROM, s, e) for s, e in annotated}, genome, "+",
        boundary_error_window=0, counters=counters, **kw,
    )


def _apply(read, replacements, genome):
    modified, applied = jr._apply_replacements_to_read(read, replacements, genome, "+", 0.25, 15)
    return applied, modified.cigarstring


# ---------------------------------------------------------------------------
# F1 — the ranking is over candidates the surgery can write
# ---------------------------------------------------------------------------

class TestRealizability:
    """Winner W = (ns+2, ne+2): a both-boundary move whose flipping bases differ,
    so the general path plants D(2)..N..I(2) and the surgery's indel invariant
    refuses it (tests/test_2h_noncanon_preponderance.py pins that refusal).
    Runner-up R = (ns, ne-3): a canonical single-boundary acceptor move the
    surgery writes as 3D."""

    @staticmethod
    def _setup(monkeypatch, *, with_incumbent):
        read, g, idx, ns, ne = _locus(CLEAN, "CTA", plants=[(-6, "CAG")])
        stock, W, R = (ns, ne), (ns + 2, ne + 2), (ns, ne - 3)
        assert jr._canonical_tier(*R, g, "+") <= _CANONICAL_TIER_MAX
        assert not jr._realizable(read, idx, ns, ne, *W, g, "+", 0.25, 15)
        assert jr._realizable(read, idx, ns, ne, *R, g, "+", 0.25, 15)
        scores = {W: 0.0, R: 0.5}
        if with_incumbent:
            scores[stock] = 2.0
        _fake_scores(monkeypatch, scores)
        return read, g, idx, stock, W, R, scores

    @pytest.mark.parametrize("with_incumbent", [False, True],
                             ids=["incumbent_absent_like_production", "incumbent_scored"])
    def test_refused_winner_yields_to_the_runner_up(self, monkeypatch, with_incumbent):
        read, g, idx, stock, W, R, scores = self._setup(monkeypatch, with_incumbent=with_incumbent)
        counters = Counter()
        repl = _refine(read, g, scores, annotated=[R], counters=counters)
        assert repl == [(idx, *stock, *R)]
        assert counters["unrealizable_winner_skipped"] == 1
        assert _apply(read, repl, g) == (True, "40M97N3D40M")

    def test_a_writable_winner_is_untouched(self, monkeypatch):
        """Control: when the head of the ranking is writable nothing changes —
        here R itself wins (W removed from the pool)."""
        read, g, idx, stock, W, R, scores = self._setup(monkeypatch, with_incumbent=True)
        del scores[W]
        counters = Counter()
        repl = _refine(read, g, scores, annotated=[R], counters=counters)
        assert repl == [(idx, *stock, *R)]
        assert counters["unrealizable_winner_skipped"] == 0

    def test_the_incumbent_stops_the_walk(self, monkeypatch):
        """W refused, then the incumbent outranks R: the read stays, one skip counted."""
        read, g, idx, stock, W, R, scores = self._setup(monkeypatch, with_incumbent=True)
        scores[stock] = 0.25            # between W (0.0) and R (0.5)
        counters = Counter()
        assert _refine(read, g, scores, annotated=[R], counters=counters) == []
        assert counters["unrealizable_winner_skipped"] == 1

    def test_nothing_writable_and_no_incumbent_stays(self, monkeypatch):
        read, g, idx, stock, W, R, scores = self._setup(monkeypatch, with_incumbent=False)
        del scores[R]
        counters = Counter()
        assert _refine(read, g, scores, annotated=[], counters=counters) == []
        assert counters["unrealizable_winner_skipped"] == 1

    def test_the_counter_reaches_the_driver_stats(self, monkeypatch, tmp_path):
        read, g, idx, stock, W, R, scores = self._setup(monkeypatch, with_incumbent=False)
        path = tmp_path / "in.bam"
        with pysam.AlignmentFile(str(path), "wb", header=_header()) as fh:
            fh.write(read)
        out = tmp_path / "out.bam"
        stats = jr.refine_bam_junctions(
            str(path), str(out), aligner_bams=[], annotated_junctions=set(),
            genome={CHROM: g},
            prebuilt_junction_pool={(CHROM, *W), (CHROM, *R)},
            prebuilt_annotated_set={(CHROM, *R)},
            boundary_error_window=0, sort_and_index=False, n_workers=1,
        )
        assert stats["refined"] == 1
        assert stats["unrealizable_winner_skipped"] == 1
        with pysam.AlignmentFile(str(out), "rb") as fh:
            assert [r.cigarstring for r in fh] == ["40M97N3D40M"]


# ---------------------------------------------------------------------------
# F2 — a non-canonical winner loses to an annotated canonical alternative within delta
# ---------------------------------------------------------------------------

class TestAnnotatedAlternative:
    """Stock S = (ns, ne) GT..CAG.  A second CAG is planted three bases upstream,
    so A = (ns, ne-3) is canonical (tier 0) and is the ANNOTATED candidate.
    W = (ns, ne-1) is non-canonical (acceptor trinucleotide GCA -> tier 4),
    shares the donor with A and sits 2 nt from it.  Both are single-boundary
    acceptor moves the surgery writes as 1D / 3D.  The gate would refuse W on
    its own (S is a canonical site 1 nt away) and FREEZE the read; the rule
    sends it to A instead."""

    @staticmethod
    def _setup(monkeypatch, *, scores):
        read, g, idx, ns, ne = _locus(CLEAN, "CTA", plants=[(-6, "CAG")])
        S, A, W = (ns, ne), (ns, ne - 3), (ns, ne - 1)
        assert jr._canonical_tier(*S, g, "+") <= _CANONICAL_TIER_MAX
        assert jr._canonical_tier(*A, g, "+") <= _CANONICAL_TIER_MAX
        assert jr._canonical_tier(*W, g, "+") > _CANONICAL_TIER_MAX
        _fake_scores(monkeypatch, {k: v for k, v in scores.items() if k in (S, A, W)}
                     if scores else {})
        return read, g, idx, S, A, W

    def test_the_annotated_alternative_wins(self, monkeypatch):
        read, g, idx, S, A, W = self._setup(monkeypatch, scores=None)
        scores = {W: 0.4, A: 0.5}                       # incumbent absent, as in production
        _fake_scores(monkeypatch, scores)
        counters = Counter()
        repl = _refine(read, g, scores, annotated=[A], counters=counters)
        assert repl == [(idx, *S, *A)]
        assert counters["noncanon_dest_lost_to_annotated_alt"] == 1
        assert counters["noncanon_destination_refused"] == 0   # the gate saw an annotated canonical destination
        assert _apply(read, repl, g) == (True, "40M97N3D40M")

    def test_with_a_scored_incumbent_that_ranks_below(self, monkeypatch):
        read, g, idx, S, A, W = self._setup(monkeypatch, scores=None)
        scores = {W: 0.4, A: 0.5, S: 2.0}
        _fake_scores(monkeypatch, scores)
        counters = Counter()
        assert _refine(read, g, scores, annotated=[A], counters=counters) == [(idx, *S, *A)]
        assert counters["noncanon_dest_lost_to_annotated_alt"] == 1

    def test_a_scored_incumbent_that_outranks_the_alternative_falls_through(self, monkeypatch):
        """No redirect when the incumbent ranks above A: W goes on to the gate,
        which refuses it on its own terms (S is a canonical site 1 nt away)."""
        read, g, idx, S, A, W = self._setup(monkeypatch, scores=None)
        scores = {W: 0.4, S: 0.45, A: 0.5}
        _fake_scores(monkeypatch, scores)
        counters = Counter()
        assert _refine(read, g, scores, annotated=[A], counters=counters) == []
        assert counters["noncanon_dest_lost_to_annotated_alt"] == 0
        assert counters["noncanon_destination_refused_alternative_within_delta"] == 1

    def test_the_incumbent_as_alternative_is_the_gate_s_business(self, monkeypatch):
        """S annotated canonical, W 1 nt away: the only alternative is the
        incumbent -> the rule is silent and the gate refuses W as before, so the
        refusal counter the tester reads and the A/B off-switch keep their meaning
        (tests/test_2h_noncanon_preponderance.py pins both)."""
        read, g, idx, S, A, W = self._setup(monkeypatch, scores=None)
        scores = {W: 0.4, S: 2.0}                        # clears the 1.0 R1 hold: the gate decides
        _fake_scores(monkeypatch, scores)
        counters = Counter()
        assert _refine(read, g, scores, annotated=[S], counters=counters) == []
        assert counters["noncanon_dest_lost_to_annotated_alt"] == 0
        assert counters["noncanon_destination_refused_alternative_within_delta"] == 1

    def test_beyond_delta_the_rule_is_silent(self, monkeypatch):
        """A moved to 4 nt from W (delta is 2): no redirect; with the gate off W
        wins as before, so the test isolates the rule from the gate."""
        monkeypatch.setattr(jr, "_NONCANON_PREPONDERANCE", False)
        read, g, idx, ns, ne = _locus(CLEAN, "CTA", plants=[(-8, "CAG")])
        S, A, W = (ns, ne), (ns, ne - 5), (ns, ne - 1)
        assert jr._canonical_tier(*A, g, "+") <= _CANONICAL_TIER_MAX < jr._canonical_tier(*W, g, "+")
        scores = {W: 0.4, A: 0.5}
        _fake_scores(monkeypatch, scores)
        counters = Counter()
        assert _refine(read, g, scores, annotated=[A], counters=counters) == [(idx, *S, *W)]
        assert counters["noncanon_dest_lost_to_annotated_alt"] == 0

    def test_an_unannotated_canonical_neighbor_does_not_redirect(self, monkeypatch):
        """The alternative must be ANNOTATED: A canonical but unannotated -> no
        redirect (the gate, off here, is what judges such a W)."""
        monkeypatch.setattr(jr, "_NONCANON_PREPONDERANCE", False)
        read, g, idx, S, A, W = self._setup(monkeypatch, scores=None)
        scores = {W: 0.4, A: 0.5}
        _fake_scores(monkeypatch, scores)
        counters = Counter()
        assert _refine(read, g, scores, annotated=[], counters=counters) == [(idx, *S, *W)]
        assert counters["noncanon_dest_lost_to_annotated_alt"] == 0

    def test_motif_blind_disables_the_rule(self, monkeypatch):
        read, g, idx, S, A, W = self._setup(monkeypatch, scores=None)
        scores = {W: 0.4, A: 0.5}
        _fake_scores(monkeypatch, scores)
        counters = Counter()
        assert _refine(read, g, scores, annotated=[A], counters=counters,
                       motif_blind=True) == [(idx, *S, *W)]
        assert counters["noncanon_dest_lost_to_annotated_alt"] == 0

    def test_helper_requires_a_shared_boundary(self):
        """(ns+1, ne-3) is within 2 nt on BOTH boundaries but shares neither
        exactly -> not an alternative under this rule."""
        read, g, idx, ns, ne = _locus(CLEAN, "CTA", plants=[(-6, "CAG")])
        W, A = (ns, ne - 1), (ns, ne - 3)
        tup = lambda s, e: (0.5, 0, 1, 0, 0, s, e, 0)   # noqa: E731
        annotated = {(CHROM, *A), (CHROM, ns + 1, ne - 3)}
        assert jr._annotated_alternative([tup(*A)], CHROM, *W, g, "+", annotated) == tup(*A)
        assert jr._annotated_alternative([tup(ns + 1, ne - 3)], CHROM, *W, g, "+", annotated) is None
        assert jr._annotated_alternative([tup(*W)], CHROM, *W, g, "+", annotated | {(CHROM, *W)}) is None


# ---------------------------------------------------------------------------
# F3 — the 2H counters reach the stats TSV
# ---------------------------------------------------------------------------

class TestModule2hCountersInStats:
    COUNTERS = {
        "n_op_reads": 20, "refined": 7, "noncanon_destination_refused": 3,
        "noncanon_destination_refused_alternative_within_delta": 3,
        "unrealizable_winner_skipped": 1, "noncanon_dest_lost_to_annotated_alt": 2,
    }

    @staticmethod
    def _rows(path):
        rows, order = {}, []
        for line in Path(path).read_text().splitlines():
            if line.startswith("#") or line.startswith("metric\t") or not line.strip():
                continue
            parts = line.split("\t")
            rows[parts[0]] = parts
            order.append(parts[0])
        return rows, order

    def test_rows_follow_module_2h_failed_with_dash_percent(self, tmp_path):
        from rectify.core.bam.processing_stats import ProcessingStats, write_stats_tsv
        stats = ProcessingStats()
        stats.total_reads_in_bam = 10
        stats.module_2h_failed = "RuntimeError: boom"
        stats.module_2h_counters = dict(self.COUNTERS)
        out = tmp_path / "x_stats.tsv"
        write_stats_tsv(stats, str(out))
        rows, order = self._rows(out)
        for k, n in self.COUNTERS.items():
            assert rows[f"module_2h_{k}"][1:3] == [str(n), "-"], k
        keys = [f"module_2h_{k}" for k in sorted(self.COUNTERS)]
        assert order[:1 + len(keys)] == ["module_2h_failed"] + keys
        assert order[1 + len(keys)] == "total_reads_in_bam"

    def test_no_rows_when_2h_did_not_run(self, tmp_path):
        from rectify.core.bam.processing_stats import ProcessingStats, write_stats_tsv
        out = tmp_path / "y_stats.tsv"
        write_stats_tsv(ProcessingStats(), str(out))
        rows, _ = self._rows(out)
        assert not [k for k in rows if k.startswith("module_2h_")]

    def test_merge_sums_and_parallel_round_trip_keeps_types(self):
        import dataclasses
        import json
        from rectify.core.bam import parallel as par
        from rectify.core.bam.processing_stats import ProcessingStats
        a, b = ProcessingStats(), ProcessingStats()
        a.module_2h_counters = {"refined": 1, "unrealizable_winner_skipped": 1}
        b.module_2h_counters = {"refined": 2}
        b.module_2h_failed = "ValueError: x"
        a.merge(b)
        assert a.module_2h_counters == {"refined": 3, "unrealizable_winner_skipped": 1}
        assert a.module_2h_failed == "ValueError: x"
        back = par._stats_from_dict(json.loads(json.dumps(a.to_dict())))
        assert back.module_2h_counters == a.module_2h_counters
        assert back.module_2h_failed == "ValueError: x"        # was int()-coerced before F3
        json.dumps(dataclasses.asdict(a))                        # provenance path still serializes


# ---------------------------------------------------------------------------
# Replay of the tester's bundle (d4 format), skip-if-absent
# ---------------------------------------------------------------------------

BUNDLE =Path.home() / "work/rectify/dev/sumner_misplaced_panel_20260904/holdout/events/4f9102d"
PREFIX = "t0_4f9102d_regress"
TABLES = RECTIFY_ROOT / "rectify/data/genomes/homo_sapiens/penalty_tables"
_OPC = "MIDNSHP=X"
_REF = {0, 2, 3, 7, 8}


def _parse_cigar(s):
    return [(_OPC.index(o), int(n)) for n, o in re.findall(r"(\d+)([MIDNSHP=X])", s)]


def _n_ops(cig, rs):
    out, pos = [], rs
    for op, n in cig:
        if op == 3:
            out.append((pos, pos + n))
        if op in _REF:
            pos += n
    return out


def _load_bundle():
    reads, lib = {}, None
    with open(BUNDLE / f"{PREFIX}.stock.sam") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                head = line[1:].split()
                if head and head[0] != "{src}":
                    lib = head[0]
                continue
            f = line.split("\t")
            reads.setdefault(f[0], {"lib": lib, "arms": []})["arms"].append(f)
    for e in reads.values():
        e["stock"], e["baseline"], e["new"] = e["arms"]
    slices, cur = {}, None
    with open(BUNDLE / f"{PREFIX}.slices.fa") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                m = re.match(r">(\S+):(\d+)-(\d+) offset0based=(\d+)", line)
                cur = (m.group(1), int(m.group(4)), int(m.group(3)))
                slices[cur] = []
            else:
                slices[cur].append(line.upper())
    slices = {k: "".join(v) for k, v in slices.items()}
    pool, annot = defaultdict(set), defaultdict(set)
    with open(BUNDLE / f"{PREFIX}.pool.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["library"], row["read"])
            j = (row["chrom"], int(row["pool_start"]), int(row["pool_end"]))
            pool[key].add(j)
            if row["annotated"] == "1":
                annot[key].add(j)
    with open(BUNDLE / f"{PREFIX}.annotated.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["library"], row["read"])
            j = (row["chrom"], int(row["start"]), int(row["end"]))
            annot[key].add(j)
            pool[key].add(j)
    return reads, slices, pool, annot


def _replay(prefix):
    """``(replay N-ops, baseline N-ops, new N-ops, counters, local annotated set, tier fn)``
    for the read whose name starts with *prefix*, run through THIS tree with
    production parameters and the tester's candidate set."""
    from rectify.core.splice.hp_penalty import HpPenaltyTable

    reads, slices, pool, annot = _load_bundle()
    name, e = next((n, e) for n, e in reads.items() if n.startswith(prefix))
    f = e["stock"]
    chrom, strand = f[2], ("-" if int(f[1]) & 16 else "+")
    key = standardize_chrom_name(chrom)
    cig, rs = _parse_cigar(f[5]), int(f[3]) - 1
    span = sum(n for op, n in cig if op in _REF)
    contig, off = next((dna, o) for (c, o, hi), dna in slices.items()
                       if c == chrom and o <= rs and rs + span <= hi)
    read = pysam.AlignedSegment(_header(key, len(contig)))
    read.query_name, read.flag, read.reference_id = name, int(f[1]), 0
    read.reference_start, read.mapping_quality = rs - off, int(f[4])
    read.cigarstring, read.query_sequence = f[5], f[9]
    read.query_qualities = pysam.qualitystring_to_array(
        f[10] if f[10] != "*" else "I" * len(f[9]))
    pk = (e["lib"], name)
    pool_local = {(key, s - off, t - off) for (_, s, t) in pool[pk]}
    annot_local = {(key, s - off, t - off) for (_, s, t) in annot[pk]}
    table = HpPenaltyTable.from_tsv(str(TABLES / "penalty_scores.tsv"),
                                    str_path=str(TABLES / "str_penalty_scores.tsv"))
    counters = Counter()
    repl = jr.refine_read_junctions(
        read, js_mod._build_junction_index(pool_local), annot_local, contig, strand,
        boundary_error_window=10, max_boundary_shift=50, search_radius=5000,
        hp_pen=0.25, W=15, max_slide=10, penalty_table=table, counters=counters,
    )
    modified, _applied = jr._apply_replacements_to_read(read, repl, contig, strand, 0.25, 15)
    out_n = [(a + off, b + off) for a, b in _n_ops(list(modified.cigartuples), modified.reference_start)]
    base_n = _n_ops(_parse_cigar(e["baseline"][5]), int(e["baseline"][3]) - 1)
    new_n = _n_ops(_parse_cigar(e["new"][5]), int(e["new"][3]) - 1)
    annotated_abs = {(s + off, t + off) for (_, s, t) in annot_local}

    def tier(js, je):
        return js_mod._canonical_tier(js - off, je - off, contig, strand)
    return out_n, base_n, new_n, counters, annotated_abs, tier


needs_bundle = pytest.mark.skipif(
    not (BUNDLE / f"{PREFIX}.stock.sam").exists() or not (TABLES / "penalty_scores.tsv").exists(),
    reason="tester's T0 4f9102d bundle (or the bundled human tables) not on this machine",
)


@needs_bundle
class TestBundleReplay:
    def test_a138b2da_moves_to_the_annotated_donor_again(self):
        """F1: the 19-nt-downstream relabel (225783368-225789258) is skipped as
        unwritable and the annotated 3-nt donor fix is written — the baseline N-op."""
        out_n, base_n, new_n, counters, annotated, tier = _replay("a138b2da")
        assert (225783368, 225789239) in base_n and (225783368, 225789239) in out_n
        assert (225783365, 225789239) in new_n          # 4f9102d had left it at the stock
        assert counters["unrealizable_winner_skipped"] == 1
        assert out_n == base_n

    def test_527c329f_goes_to_the_annotated_acceptor_not_the_freezer(self):
        """F2 (canonical incumbent): GT-GA 103377463 beat the annotated GT-AG
        103377462 by 0.125 under the corrected window; the gate refused it and
        4f9102d left the read 12 nt off at the stock.  Now the annotated
        alternative wins and the gate never fires."""
        out_n, base_n, new_n, counters, annotated, tier = _replay("527c329f")
        assert (103377106, 103377462) in base_n and (103377106, 103377462) in annotated
        assert (103377106, 103377474) in new_n           # 4f9102d: frozen at the stock
        assert (103377106, 103377462) in out_n
        assert counters["noncanon_dest_lost_to_annotated_alt"] == 1
        assert counters["noncanon_destination_refused"] == 0
        assert out_n == base_n

    def test_7628e0dd_goes_to_the_annotated_acceptor_not_gt_ca(self):
        """F2 (non-canonical incumbent GT..AC): GT-CA 154257064 (3D) beat the
        annotated GT-AG 154257062 (5D) under the corrected window."""
        out_n, base_n, new_n, counters, annotated, tier = _replay("7628e0dd")
        assert (154255755, 154257062) in base_n and (154255755, 154257062) in annotated
        assert (154255755, 154257064) in new_n and tier(154255755, 154257064) > _CANONICAL_TIER_MAX
        assert (154255755, 154257062) in out_n
        assert counters["noncanon_dest_lost_to_annotated_alt"] == 1
        assert out_n == base_n

    def test_bbc2d388_is_an_annotated_to_annotated_choice(self):
        """Documented, NOT a fix: both 149197991 (baseline) and 149197995 (4f9102d)
        are GENCODE-annotated canonical donors (minus strand) in the tester's own
        tables, 4 nt apart, and the corrected window prefers 995 (2.265 vs 3.194,
        inside the 1.0 noise floor).  The annotated-alternative rule does not
        apply (the destination is annotated); whether the smaller move should win
        an annotated-vs-annotated near-tie is a ruling question, pinned here so a
        change of behavior is noticed."""
        out_n, base_n, new_n, counters, annotated, tier = _replay("bbc2d388")
        assert (149190531, 149197991) in base_n and (149190531, 149197991) in annotated
        assert (149190531, 149197995) in new_n and (149190531, 149197995) in annotated
        assert tier(149190531, 149197995) <= _CANONICAL_TIER_MAX
        assert (149190531, 149197995) in out_n
        assert counters["noncanon_dest_lost_to_annotated_alt"] == 0
