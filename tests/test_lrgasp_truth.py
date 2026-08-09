"""LRGASP NanoSim sim-truth join tests (the NanoSim arm of the three-way).

The load-bearing check: a read's truth junctions are restricted to those it
SPANS WITH ANCHOR (the session-3 fragment bug — attaching the full transcript
junction set to a fragment read deflates recall with false FN). Plus
read_to_isoform parsing (mouse-decoy visibility + ENST filter), verbatim-version
catalogue keying, CPA-at-terminus, and a truth-table round-trip.

Hermetic: synthetic inline GTF + genome (no network).
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scripts.benchmark.sim.lrgasp_truth import (  # noqa: E402
    parse_read_to_isoform, build_truth_catalogue, spanned_anchored_junctions,
    cpa_for_span, read_truth_for_span, MIN_JUNCTION_ANCHOR,
)
from rectify.core.benchmark.truth_schema import (  # noqa: E402
    write_truth_table, read_truth_table, SplitTag,
)

# tx_plus.1: exons [0,20),[40,60),[90,110) -> introns [20,40),[60,90)
GTF = """\
chrT\tt\texon\t1\t20\t.\t+\t.\tgene_id "G1"; transcript_id "tx_plus.1";
chrT\tt\texon\t41\t60\t.\t+\t.\tgene_id "G1"; transcript_id "tx_plus.1";
chrT\tt\texon\t91\t110\t.\t+\t.\tgene_id "G1"; transcript_id "tx_plus.1";
chrT\tt\texon\t201\t260\t.\t+\t.\tgene_id "G2"; transcript_id "mono.2";
"""


def _genome():
    seq = list("C" * 300)
    for i, s in ((20, "GT"), (38, "AG"), (60, "GT"), (88, "AG")):
        for k, ch in enumerate(s):
            seq[i + k] = ch
    return {"chrT": "".join(seq)}


def _write(d, name, text):
    p = os.path.join(d, name)
    with open(p, "w") as fh:
        fh.write(text)
    return p


def test_parse_read_to_isoform_counts_and_filter():
    with tempfile.TemporaryDirectory() as d:
        r2i = _write(d, "r2i.tsv",
                     "r1\tENST001.1\nr2\tENSMUST900.2\nr3\tENST002.3\n# c\nr4\tENST001.1\n")
        all_map, c = parse_read_to_isoform(r2i)
        assert c["total"] == 4 and c["ENST"] == 3 and c["ENSMUST"] == 1, c
        assert all_map["r2"] == "ENSMUST900.2"          # mouse kept when unfiltered
        enst_map, c2 = parse_read_to_isoform(r2i, keep_prefixes=("ENST",))
        assert "r2" not in enst_map and c2["kept"] == 3  # mouse decoy dropped
        assert enst_map["r1"] == "ENST001.1"             # version kept verbatim


def test_catalogue_keys_verbatim_includes_monoexonic():
    with tempfile.TemporaryDirectory() as d:
        gtf = _write(d, "g.gtf", GTF)
        models, pairs, donors, acceptors = build_truth_catalogue(gtf, _genome())
        # spliced_only=False -> monoexonic transcript is present (sim reads come
        # from ALL transcripts), keyed by verbatim transcript_id
        assert set(models) == {"tx_plus.1", "mono.2"}
        assert len(pairs) == 2  # two annotated introns on tx_plus.1


def _cat(d):
    gtf = _write(d, "g.gtf", GTF)
    return build_truth_catalogue(gtf, _genome())


def test_full_span_read_gets_both_junctions_and_cpa():
    with tempfile.TemporaryDirectory() as d:
        models, pairs, donors, acceptors = _cat(d)
        m = models["tx_plus.1"]
        js = spanned_anchored_junctions(m, 0, 109, pairs, donors, acceptors)
        assert [(j.intron_start, j.intron_end) for j in js] == [(20, 40), (60, 90)]
        assert cpa_for_span(m, 0, 109) == 109  # + strand reaches 3' terminus


def test_fragment_read_drops_unspanned_junctions():
    """The session-3 guard: a fragment covering only the FIRST intron must NOT
    carry the second (unspanned) junction as truth, else recall is deflated."""
    with tempfile.TemporaryDirectory() as d:
        models, pairs, donors, acceptors = _cat(d)
        m = models["tx_plus.1"]
        # covers [0,60]: intron1 anchored (cov_hi 60 >= 40+9), intron2 not (needs >=99)
        js = spanned_anchored_junctions(m, 0, 60, pairs, donors, acceptors)
        assert [(j.intron_start, j.intron_end) for j in js] == [(20, 40)]
        assert cpa_for_span(m, 0, 60) is None  # doesn't reach terminus
        # covers [45,109]: intron1 NOT anchored on the left (45 > 20-10), intron2 is
        js2 = spanned_anchored_junctions(m, 45, 109, pairs, donors, acceptors)
        assert [(j.intron_start, j.intron_end) for j in js2] == [(60, 90)]


def test_insufficient_flank_is_not_anchorable():
    with tempfile.TemporaryDirectory() as d:
        models, pairs, donors, acceptors = _cat(d)
        m = models["tx_plus.1"]
        # start at 15 -> only 5bp of exon1 upstream of intron1 (< anchor 10)
        js = spanned_anchored_junctions(m, 15, 109, pairs, donors, acceptors)
        assert (20, 40) not in [(j.intron_start, j.intron_end) for j in js]


def test_truth_table_round_trip():
    with tempfile.TemporaryDirectory() as d:
        models, pairs, donors, acceptors = _cat(d)
        m = models["tx_plus.1"]
        rt = read_truth_for_span("ONT_simulated_read_7", m, 0, 109,
                                 pairs, donors, acceptors, split=SplitTag.TEST)
        p = os.path.join(d, "truth.tsv")
        write_truth_table([rt], p)
        back = read_truth_table(p)
        assert len(back) == 1
        b = back[0]
        assert b.read_id == "ONT_simulated_read_7"
        assert b.true_transcript == "tx_plus.1"
        assert [(j.intron_start, j.intron_end) for j in b.junctions] == [(20, 40), (60, 90)]
        assert b.true_cpa == 109


if __name__ == "__main__":
    test_parse_read_to_isoform_counts_and_filter()
    test_catalogue_keys_verbatim_includes_monoexonic()
    test_full_span_read_gets_both_junctions_and_cpa()
    test_fragment_read_drops_unspanned_junctions()
    test_insufficient_flank_is_not_anchorable()
    test_truth_table_round_trip()
    print(f"all lrgasp_truth tests passed (anchor={MIN_JUNCTION_ANCHOR})")
