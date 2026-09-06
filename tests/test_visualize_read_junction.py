"""rectify.visualize.read_junction -- the per-read junction PNG on a synthetic 3-arm bundle.

A 600-bp slice genome ``chrT`` with a two-exon transcript (exon 1 = [100, 200), intron [200, 400)
GT..AG, exon 2 = [400, 520)) and one read whose 5' end is exon 1's last 20 nt + exon 2:

* ``stock``     soft-clips the 20 nt (``20S120M``) -- no junction;
* ``baseline``  the same record;
* ``candidate`` places the 20 nt on exon 1 across the annotated junction (``20M200N120M``),
  with one mismatch, a 1-nt deletion and a 1-nt insertion planted in exon 2.

The arms disagree at ONE junction, so ``junction_index='auto'`` must give exactly one panel.
"""
from __future__ import annotations

import os
import random

import pytest

pysam = pytest.importorskip("pysam")
pytest.importorskip("matplotlib")

from rectify.visualize import read_junction as RJ  # noqa: E402

try:
    RJ.tokens()
    HAVE_TOKENS = True
except ImportError:
    HAVE_TOKENS = False

pytestmark = pytest.mark.skipif(not HAVE_TOKENS, reason="rna-figure tokens not available on this tree/machine")

EX1, IN_S, IN_E, EX2_E = 100, 200, 400, 520


def _genome(n=600, seed=7):
    rng = random.Random(seed)
    seq = list("".join(rng.choice("ACGT") for _ in range(n)))
    seq[IN_S:IN_S + 2] = "GT"
    seq[IN_E - 2:IN_E] = "AG"
    return "".join(seq)


def _bam(path, header, name, start, cigar, seq, mapq=60):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigarstring = cigar
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    with pysam.AlignmentFile(path, "wb", header=header) as fh:
        fh.write(a)
    pysam.index(path)


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    d = tmp_path_factory.mktemp("rj")
    g = _genome()
    fa = d / "chrT.fa"
    fa.write_text(">chrT\n" + "\n".join(g[i:i + 60] for i in range(0, len(g), 60)) + "\n")
    pysam.faidx(str(fa))
    gtf = d / "t.gtf"
    attrs = 'gene_id "G1"; transcript_id "T1"; gene_name "TEST1";'
    gtf.write_text("\t".join(["chrT", "x", "exon", str(EX1 + 1), str(IN_S), ".", "+", ".", attrs]) + "\n"
                   + "\t".join(["chrT", "x", "exon", str(IN_E + 1), str(EX2_E), ".", "+", ".", attrs]) + "\n")
    exon1_tail = g[IN_S - 20:IN_S]
    exon2 = list(g[IN_E:EX2_E])
    # plant edits in exon 2 as the CANDIDATE and STOCK see them (same read sequence in every arm)
    exon2[10] = "A" if exon2[10] != "A" else "C"          # mismatch at 410
    del exon2[40]                                          # deletion of 440 -> 120 ref bases from 119 read bases
    exon2.insert(80, "G")                                  # insertion before ref 481 (40 read bases past the deletion)
    read = exon1_tail + "".join(exon2)
    n2 = len("".join(exon2))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrT", "LN": len(g)}]}
    name = "aaaa1111-0000-4000-8000-000000000001"
    exon2_cigar = "40M1D40M1I%dM" % (n2 - 40 - 40 - 1)
    _bam(str(d / "stock.bam"), header, name, IN_E, "20S" + exon2_cigar, read)
    _bam(str(d / "baseline.bam"), header, name, IN_E, "20S" + exon2_cigar, read)
    _bam(str(d / "candidate.bam"), header, name, IN_S - 20, "20M%dN" % (IN_E - IN_S) + exon2_cigar, read)
    (d / "manifest.tsv").write_text(
        "# synthetic\nread\tlibrary\tchrom\tstrand\tclass\tstock_junction\tcorrected_junction_baseline\t"
        "corrected_junction_candidate\n" + f"{name}\tLIB\tchrT\t+\tsynthetic\t\t\t{IN_S}-{IN_E}\n")
    return dict(dir=d, fa=fa, gtf=gtf, name=name, genome=g)


def test_data_model_matches_the_planted_edits(bundle):
    G = RJ.Genome(bundle["fa"])
    A = RJ.Annotation(bundle["gtf"])
    recs = RJ.load_arms(bundle["dir"], "aaaa1111")
    assert list(recs) == ["stock", "baseline", "candidate"], "stock first, then the manifest header order"
    v = RJ.arm_view("candidate", recs["candidate"], G, A)
    assert [(i.start, i.end, i.annotated) for i in v.introns] == [(IN_S, IN_E, True)]
    assert v.introns[0].motif == "GT-AG"
    b5 = v.five_block
    assert (b5.rs, b5.re, b5.M, b5.X) == (IN_S - 20, IN_S, 20, 0)
    b3 = v.blocks[1]
    assert (b3.X, b3.I, b3.D) == (1, 1, 1)
    assert v.aligned[IN_E + 40] == "-"                     # the deletion column
    assert IN_E + 81 in v.insertions                       # the insertion boundary (40 M after the 1-nt deletion)
    s = RJ.arm_view("stock", recs["stock"], G, A)
    assert s.five_clip == 20 and not s.introns
    assert s.clip_fit and s.clip_fit["exon_end"] == IN_S and s.clip_fit["matches"] == 20 and s.clip_fit["offset"] == 0
    line = RJ.verdict_line(v, G, RJ.Frame("chrT", "+", (IN_S, IN_E), 16))
    assert "GT-AG" in line and "MQ 60" in line and "run 20|10" in line


def test_junction_selection_auto_is_the_disagreement(bundle):
    G = RJ.Genome(bundle["fa"])
    A = RJ.Annotation(bundle["gtf"])
    views = {n: RJ.arm_view(n, a, G, A) for n, a in RJ.load_arms(bundle["dir"], "aaaa1111").items()}
    assert RJ.select_junctions(views, "auto") == [[(IN_S, IN_E)]]
    assert RJ.select_junctions(views, "all") == [[(IN_S, IN_E)]]
    # no disagreement -> the 5'-most junction
    only = {"stock": views["candidate"], "candidate": views["candidate"]}
    assert RJ.select_junctions(only, "auto") == [[(IN_S, IN_E)]]


def test_render_png_floors_and_panel_count(bundle, tmp_path):
    out = tmp_path / "aaaa1111.png"
    res = RJ.render_read(bundle["dir"], "aaaa1111", bundle["fa"], bundle["gtf"], out)
    assert os.path.exists(out) and os.path.getsize(out) > 10_000
    assert res["panels"] == 1 == len(res["junctions"])       # one disagreement -> one panel
    assert res["junctions"] == [(IN_S, IN_E)]
    assert res["arms"] == ["stock", "baseline", "candidate"]
    assert res["floors_ok"], "P.check_floors must pass"
    assert res["fig"].get_figwidth() == pytest.approx(RJ.tokens()["geometry"]["column_in"]["double"])
    side = tmp_path / "aaaa1111.md"
    assert side.exists() and "GT-AG" in side.read_text()
    # the sidecar guard: a foreign .md beside the PNG is never overwritten
    other = tmp_path / "other.png"
    (tmp_path / "other.md").write_text("# somebody's findings\n")
    with pytest.raises(FileExistsError):
        RJ.render_read(bundle["dir"], "aaaa1111", bundle["fa"], bundle["gtf"], other)


def test_window_refuses_what_the_type_floor_cannot_hold(bundle, tmp_path):
    with pytest.raises(ValueError):
        RJ.render_read(bundle["dir"], "aaaa1111", bundle["fa"], bundle["gtf"], tmp_path / "x.png", window=80)


def test_slice_genome_maps_back_to_genome_coordinates(tmp_path):
    g = _genome()
    fa = tmp_path / "s.fa"
    fa.write_text(">chrT:1000-1600 offset0based=1000\n" + g + "\n")
    pysam.faidx(str(fa))
    G = RJ.Genome(fa)
    assert G.fetch("chrT", 1000 + IN_S, 1000 + IN_S + 2) == "GT"
    assert G.fetch("chrT", 10, 20) is None
