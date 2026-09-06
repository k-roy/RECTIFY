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
# Replay of the tester's bundle (d4 format), skip-if-absent
# ---------------------------------------------------------------------------

BUNDLE = Path.home() / "work/rectify/dev/sumner_misplaced_panel_20260904/holdout/events/4f9102d"
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
