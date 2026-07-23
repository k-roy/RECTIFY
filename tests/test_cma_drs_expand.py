"""CMA v1 format-core acceptance tests (planning/254 milestone "format-core").

Losslessness is asserted on the load-bearing FIELD VIEW (planning/254 §3.3):
QNAME / is_reverse / RNAME / POS / MAPQ / CIGAR / SEQ / QUAL / NM / MD — never
raw BAM bytes.

The bundled 5-aligner DRS fixture is soft-clip-only, seq-bearing, no RN, 18/36
reverse. The reverse-strand + hard-clip frame (the #1 silent-failure risk,
planning/254 §2.4/§8-1) and the gapmm2 SEQ=None / no-donor cases are NOT in the
fixture and are therefore SYNTHESIZED below.
"""

from array import array
from collections import defaultdict
from pathlib import Path

import pysam
import pytest

import rectify
from rectify.core.consensus.consensus import _normalize_bam_read_name
from rectify.core.consensus.extract import extract_alignment_info
from rectify.core.consensus.select import select_best_alignment
from rectify.core.multialign import (
    build_cma,
    expand,
    load_aligner_records,
    validate_cma,
)
from rectify.core.multialign.cma_schema import TAG_PAYLOAD, decode_eq_seq

FIXTURE = Path(rectify.__file__).parent / "data" / "validation" / "aligners"
ALIGNERS = ["minimap2", "mapPacBio", "gapmm2", "deSALT", "uLTRA"]
PANEL = ALIGNERS


def _bam_paths():
    return {a: str(FIXTURE / f"validation_reads.{a}.bam") for a in ALIGNERS}


@pytest.fixture(scope="module")
def genome():
    from rectify.data import get_bundled_genome_path

    fa = pysam.FastaFile(str(get_bundled_genome_path("saccharomyces_cerevisiae")))
    try:
        return {r: fa.fetch(r) for r in fa.references}
    finally:
        fa.close()


def placement_view(rec, genome=None):
    """The load-bearing, aligner-scratch-free comparison view (QUAL excluded).

    SEQ is compared as DECODED nucleotides (resolving SAM ``=`` match-encoding):
    ``=`` is a placement-relative encoding of the same read bases, so nucleotide
    identity — not the ``=``/explicit byte form — is the losslessness contract.

    QUAL is handled separately (see test_fixture_roundtrip): it is read-intrinsic
    and NOT consumed by any algorithm (planning/254 §2.6), so the CMA stores it
    once from the payload donor and reconstructs it for every placement. The
    fixture exposes per-aligner QUAL anomalies the CMA deliberately canonicalizes
    (deSALT reverse-orients reverse-strand QUAL; uLTRA drops QUAL). Byte-exact
    QUAL reconstruction on consistent inputs is proven by the synthetic tests.
    """
    tags = dict(rec.get_tags())
    return (
        _normalize_bam_read_name(rec.query_name or ""),
        bool(rec.is_reverse),
        rec.reference_name,
        rec.reference_start,
        rec.mapping_quality,
        rec.cigarstring,
        decode_eq_seq(rec, genome),
        tags.get("NM"),
        tags.get("MD"),
    )


# --------------------------------------------------------------------------- #
# Synthetic-record helpers
# --------------------------------------------------------------------------- #
def mk_header(contigs=(("chr1", 100000),)):
    return pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "unsorted"},
            "SQ": [{"SN": n, "LN": l} for n, l in contigs],
        }
    )


def mk_rec(header, qname, flag, rname, pos, mapq, cigar, seq, quals, tags=None):
    r = pysam.AlignedSegment(header)
    r.query_name = qname
    r.flag = flag
    r.reference_id = header.references.index(rname)
    r.reference_start = pos
    r.mapping_quality = mapq
    r.cigarstring = cigar
    if seq is not None:
        r.query_sequence = seq
        if quals is not None:
            r.query_qualities = array("B", quals)
    if tags:
        r.set_tags([(t, v, ty) for t, v, ty in tags])
    return r


# --------------------------------------------------------------------------- #
# 1. Real fixture round-trip (18 reverse-strand frames exercised)
# --------------------------------------------------------------------------- #
def test_fixture_roundtrip(tmp_path, genome):
    bam_paths = _bam_paths()

    truth = defaultdict(dict)
    for a, p in bam_paths.items():
        with pysam.AlignmentFile(p, "rb") as bam:
            for r in bam:
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                truth[_normalize_bam_read_name(r.query_name or "")][a] = placement_view(r, genome)

    with pysam.AlignmentFile(bam_paths["minimap2"], "rb") as bam:
        header = bam.header

    cma = str(tmp_path / "drs.cma.bam")
    stats = build_cma(load_aligner_records(bam_paths), header, cma, PANEL, genome=genome)
    assert stats["reads"] == len(truth) == 36

    assert validate_cma(cma) == []

    seen = set()
    for key, adict in expand(cma):
        seen.add(key)
        assert set(adict) == set(truth[key]), f"aligner set mismatch for {key}"
        fwd_quals = set()
        for a, rec in adict.items():
            # Load-bearing placement fields round-trip EXACTLY (incl. reverse frames).
            assert placement_view(rec, genome) == truth[key][a], f"placement mismatch: {key} / {a}"
            # Every placement gets a QUAL of the right length, frame-sliced.
            q = rec.query_qualities
            assert q is not None and len(q) == len(rec.query_sequence), f"bad QUAL: {key}/{a}"
            fwd_quals.add(tuple(q[::-1]) if rec.is_reverse else tuple(q))
        # QUAL is canonicalized: all aligners of a read agree on ONE forward-read
        # quality (correctly reverse-framed per placement), sourced from the payload.
        assert len(fwd_quals) == 1, f"{key}: QUAL not canonical across placements"
    assert seen == set(truth)


def test_fixture_selection_equivalence(tmp_path, genome):
    """PRIMARY deletion gate (planning/254 §3.3, doc M1 acceptance): the pipeline's
    own consumer — extract_alignment_info + select_best_alignment — picks the
    IDENTICAL winner and corrected 3' end on expand(CMA) vs on the original
    per-aligner records.

    Comparison is against the ``=``-DECODED originals, because that is what
    production BAMs carry (the bundled fixture is ``calmd -e`` ``=``-encoded, an
    artifact of shrinking it). Feeding raw ``=``-SEQ to extract_alignment_info
    mis-scores 5' rescue and changes the winner for 14/36 fixture reads — the CMA
    decodes ``=`` and is reference-correct, so this gate deliberately isolates
    selection fidelity from that (separately-recorded) improvement.
    """
    bam_paths = _bam_paths()
    with pysam.AlignmentFile(bam_paths["minimap2"], "rb") as bam:
        header = bam.header
    orig = {k: ad for k, ad in load_aligner_records(bam_paths)}

    cma = str(tmp_path / "drs.cma.bam")
    build_cma(((k, orig[k]) for k in orig), header, cma, PANEL, genome=genome)
    expanded = {k: ad for k, ad in expand(cma)}
    assert set(expanded) == set(orig)

    def _explicit(reads):
        out = {}
        for a, r in reads.items():
            c = pysam.AlignedSegment.fromstring(r.to_string(), header)
            d = decode_eq_seq(r, genome)
            if d is not None and d != c.query_sequence:
                q = c.query_qualities
                c.query_sequence = d
                c.query_qualities = q
            out[a] = c
        return out

    def _winner(reads):
        al = {a: extract_alignment_info(r, a, genome) for a, r in reads.items()}
        res = select_best_alignment(al, genome, None, tiebreak="rectify")
        c3 = getattr(al.get(res.best_aligner), "corrected_3prime", None)
        return res.best_aligner, c3

    for k in orig:
        assert _winner(expanded[k]) == _winner(_explicit(orig[k])), f"selection diverged for {k}"


def test_fixture_one_seq_copy_per_read(tmp_path):
    """Footprint correctness: exactly one SEQ-bearing (payload) record per read."""
    bam_paths = _bam_paths()
    with pysam.AlignmentFile(bam_paths["minimap2"], "rb") as bam:
        header = bam.header
    cma = str(tmp_path / "drs.cma.bam")
    build_cma(load_aligner_records(bam_paths), header, cma, PANEL)

    seq_bearing = defaultdict(int)
    payload_tag = defaultdict(int)
    with pysam.AlignmentFile(cma, "rb") as f:
        for r in f:
            if r.query_sequence is not None:
                seq_bearing[r.query_name] += 1
            if r.has_tag(TAG_PAYLOAD):
                payload_tag[r.query_name] += 1
    assert seq_bearing and all(v == 1 for v in seq_bearing.values())
    assert all(v == 1 for v in payload_tag.values())


# --------------------------------------------------------------------------- #
# 2. Reverse-strand + hard-clip frame — the load-bearing worked check (§2.4)
# --------------------------------------------------------------------------- #
def test_reverse_strand_hardclip_frame(tmp_path):
    h = mk_header()
    F = "AAAAAAAACCCCGGGT"  # L = 16, payload stored forward
    payload = mk_rec(h, "r1", 0, "chr1", 100, 60, "16M", F, list(range(16)))
    # deSALT reverse, 4H12M: original SEQ is the H-trimmed reverse read.
    des_seq = "GGGGTTTTTTTT"
    des_q = list(range(11, -1, -1))  # [11,10,...,0]
    desalt = mk_rec(h, "r1", 16, "chr1", 200, 50, "4H12M", des_seq, des_q)

    cma = str(tmp_path / "rc.cma.bam")
    build_cma([("r1", {"minimap2": payload, "deSALT": desalt})], h, cma, PANEL)
    assert validate_cma(cma) == []

    out = dict(expand(cma))["r1"]
    assert out["minimap2"].query_sequence == F
    # The exact reconstruction the doc pins: base=revcomp(F)="ACCCGGGGTTTTTTTT",
    # slice [4:16] = "GGGGTTTTTTTT".
    assert out["deSALT"].query_sequence == des_seq
    assert list(out["deSALT"].query_qualities) == des_q
    assert out["deSALT"].is_reverse and out["deSALT"].cigarstring == "4H12M"


# --------------------------------------------------------------------------- #
# 3. gapmm2 SEQ=None variant — SEQ restored from the payload donor
# --------------------------------------------------------------------------- #
def test_gapmm2_seq_none_variant_restored(tmp_path):
    h = mk_header()
    F = "ACGTACGTACGTACGT"
    mm = mk_rec(h, "r2", 0, "chr1", 100, 60, "16M", F, list(range(16)))
    gap = mk_rec(h, "r2", 0, "chr1", 300, 0, "16M", None, None)  # PAF→BAM: SEQ=*

    cma = str(tmp_path / "g.cma.bam")
    build_cma([("r2", {"minimap2": mm, "gapmm2": gap})], h, cma, PANEL)
    assert validate_cma(cma) == []

    out = dict(expand(cma))["r2"]
    assert out["gapmm2"].query_sequence == F  # forward/forward → base=F, slice[0:16]
    assert out["minimap2"].query_sequence == F


# --------------------------------------------------------------------------- #
# 4. Per-aligner MAPQ preserved across a shared placement (lossless MAPQ)
# --------------------------------------------------------------------------- #
def test_per_aligner_mapq_preserved(tmp_path):
    h = mk_header()
    F = "ACGTACGTACGTACGT"
    mm = mk_rec(h, "r3", 0, "chr1", 100, 60, "16M", F, list(range(16)))
    mp = mk_rec(h, "r3", 0, "chr1", 100, 40, "16M", F, list(range(16)))  # same place, diff MAPQ

    cma = str(tmp_path / "mq.cma.bam")
    build_cma([("r3", {"minimap2": mm, "mapPacBio": mp})], h, cma, PANEL)
    assert validate_cma(cma) == []

    out = dict(expand(cma))["r3"]
    assert out["minimap2"].mapping_quality == 60
    assert out["mapPacBio"].mapping_quality == 40
    # one shared placement → one SEQ-bearing record
    with pysam.AlignmentFile(cma, "rb") as f:
        n_seq = sum(1 for r in f if r.query_sequence is not None)
    assert n_seq == 1


# --------------------------------------------------------------------------- #
# 5. RN-keyed grouping — RN (not QNAME) joins aligners; RN rides the payload
# --------------------------------------------------------------------------- #
def test_rn_key_path(tmp_path):
    h = mk_header()
    F1, F2 = "ACGTACGTACGTACGT", "TTTTGGGGCCCCAAAA"

    def write(path, aligner_qsuffix):
        recs = [
            mk_rec(h, f"molA{aligner_qsuffix}", 0, "chr1", 100, 60, "16M", F1,
                   list(range(16)), tags=[("RN", 1, "i")]),
            mk_rec(h, f"molB{aligner_qsuffix}", 16, "chr1", 500, 60, "16M", F2,
                   list(range(16)), tags=[("RN", 2, "i")]),
        ]
        with pysam.AlignmentFile(path, "wb", header=h) as f:
            for r in recs:
                f.write(r)

    p1 = str(tmp_path / "a1.bam")
    p2 = str(tmp_path / "a2.bam")
    write(p1, "_mm")   # distinct QNAMEs per aligner → only RN can join them
    write(p2, "_mp")

    cma = str(tmp_path / "rn.cma.bam")
    stats = build_cma(load_aligner_records({"minimap2": p1, "mapPacBio": p2}), h, cma, PANEL)
    assert stats["reads"] == 2  # joined by RN despite different QNAMEs
    assert validate_cma(cma) == []

    # payload carries RN; both aligners recovered per read
    payload_rns = []
    with pysam.AlignmentFile(cma, "rb") as f:
        for r in f:
            if r.has_tag(TAG_PAYLOAD):
                payload_rns.append(int(r.get_tag("RN")))
    assert sorted(payload_rns) == [1, 2]

    for _key, adict in expand(cma):
        assert set(adict) == {"minimap2", "mapPacBio"}
        seqs = {r.query_sequence for r in adict.values()}
        assert len(seqs) == 1  # both aligners see the same restored SEQ


# --------------------------------------------------------------------------- #
# 6. No-donor fallback — no full-SEQ soft-clip record; lossless, uncompressed
# --------------------------------------------------------------------------- #
def test_no_donor_fallback(tmp_path):
    h = mk_header()
    gap = mk_rec(h, "r4", 0, "chr1", 100, 0, "16M", None, None)          # SEQ=* (PAF)
    des = mk_rec(h, "r4", 16, "chr1", 200, 50, "4H12M", "GGGGTTTTTTTT",
                 list(range(12)))                                          # hard-clipped, partial

    cma = str(tmp_path / "nd.cma.bam")
    build_cma([("r4", {"gapmm2": gap, "deSALT": des})], h, cma, PANEL)
    assert validate_cma(cma) == []

    out = dict(expand(cma))["r4"]
    # gapmm2 had no SEQ and there is no donor → it stays SEQ=* (not fabricated).
    assert out["gapmm2"].query_sequence is None
    # deSALT keeps its own (partial) SEQ verbatim.
    assert out["deSALT"].query_sequence == "GGGGTTTTTTTT"
    assert list(out["deSALT"].query_qualities) == list(range(12))
