"""Station C background-SV flag (planning/730 W6, Kevin 2026-08-20).

Known background-SV regions of the reference (e.g. R64 chrIII SRD1 flank-A,
deleted-and-replaced-by-Ty1 in real strain genomes per the yKR888 T2T
assembly) manufacture canonical-looking junctions from terminal-clip
bridging. A junction overlapping one must demote to
``demote_orthogonal_evidence`` on BOTH tracks — support, motif, and overhang
quality cannot save it, because the reference itself is wrong there.
"""

from pathlib import Path

import pysam
import pytest

from rectify.core.consensus.station_c import (
    find_bundled_background_sv_bed,
    load_background_sv_bed,
    pool_gate,
    write_pool_gate_outputs,
)
from rectify.core.splice.overhang_informativeness import canonical_in_class
from rectify.core.splice.region_skip import (
    YEAST_SRD1_SV_SPEC, parse_skip_regions,
)

# Same high-complexity 60-mers as test_station_c.py.
EX1 = "ACGTGATCCATGCTTACGCTGACTATCGGACTTCAGATCCGTACTGACGATCCATGCATC"
EX2 = "TTCAGCATGGACGTATCACGGATCAGCTACGTTAGCATCGACTGCATCGGATACGCATGC"
FILL = "CATGACTGCATCAGGATCCATCAGCATGACGTATCAGCATGGATCCATCACGGATCATGC"


def _mk_intron(donor: str, acceptor: str, n_fill: int = 76) -> str:
    return donor + (FILL * 3)[:n_fill] + acceptor


@pytest.fixture(scope="module")
def sv(tmp_path_factory):
    """chrA: novel canonical junction (2 clean reads) + a second novel
    canonical junction outside any SV region, + a background-SV BED covering
    only the first."""
    d = tmp_path_factory.mktemp("station_c_bg_sv")
    i_can = _mk_intron("GT", "AG")
    segs = [EX1, i_can, EX2, EX1, i_can, EX2]
    chrA = "".join(segs)
    off = [0]
    for s in segs:
        off.append(off[-1] + len(s))
    j_sv = (off[1], off[2])     # novel canonical, inside the SV BED
    j_ok = (off[4], off[5])     # novel canonical, outside
    assert canonical_in_class(chrA, *j_sv)
    assert canonical_in_class(chrA, *j_ok)

    gff = d / "mini.gff"
    gff.write_text("##gff-version 3\n"
                   f"chrA\tsgd\tgene\t1\t{len(chrA)}\t.\t+\t.\tID=g1\n")
    bed = d / "mini.background_sv.bed"
    bed.write_text("# test background-SV track\n"
                   f"chrA\t{j_sv[0] - 5}\t{j_sv[1] + 5}\ttest_SV\n")

    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrA", "LN": len(chrA)}]}
    bam = d / "mini.bam"

    def read(name, start, intron_len):
        a = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(header))
        a.query_name = name
        a.reference_id = 0
        a.reference_start = start
        a.cigartuples = [(0, len(EX1)), (3, intron_len), (0, len(EX2))]
        a.query_sequence = EX1 + EX2
        a.mapping_quality = 60
        return a

    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for n in range(2):
            out.write(read(f"sv{n}", off[0], len(i_can)))
        for n in range(2):
            out.write(read(f"ok{n}", off[3], len(i_can)))

    return {"dir": d, "gff": gff, "bed": bed, "bam": bam,
            "genome": {"chrA": chrA}, "j_sv": j_sv, "j_ok": j_ok}


def _by_junction(rows):
    return {(r["chrom"], r["start"], r["end"]): r for r in rows}


def test_background_sv_demotes_despite_canonical_support_and_q(sv):
    rows, summary = pool_gate(str(sv["bam"]), sv["genome"], sv["gff"],
                              background_sv_bed=sv["bed"])
    by = _by_junction(rows)

    # Anti-vacuity: the control junction proves the fixture WOULD admit.
    ok = by[("chrA", *sv["j_ok"])]
    assert ok["verdict"] == "admit_candidate"
    assert ok["background_sv_flag"] == ""

    hit = by[("chrA", *sv["j_sv"])]
    assert hit["canonical_in_class"] == 1
    assert hit["support"] == 2
    assert hit["q_max"] >= 40          # everything that would admit it...
    assert hit["background_sv_flag"] == "test_SV"
    assert hit["verdict"] == "demote_orthogonal_evidence"  # ...cannot

    assert summary["background_sv_intervals"] == 1
    assert summary["background_sv_bed"] == str(sv["bed"])


def test_without_bed_the_sv_junction_admits(sv):
    """The demotion comes from the track, not from anything else."""
    rows, summary = pool_gate(str(sv["bam"]), sv["genome"], sv["gff"])
    by = _by_junction(rows)
    assert by[("chrA", *sv["j_sv"])]["verdict"] == "admit_candidate"
    assert summary["background_sv_intervals"] == 0
    assert summary["background_sv_bed"] is None


def test_tsv_carries_background_sv_and_over_max_intron_columns(sv, tmp_path):
    rows, summary = pool_gate(str(sv["bam"]), sv["genome"], sv["gff"],
                              background_sv_bed=sv["bed"])
    tsv, _ = write_pool_gate_outputs(rows, summary, tmp_path / "sample")
    header = tsv.read_text().splitlines()[0].split("\t")
    assert "background_sv_flag" in header
    # over_max_intron was computed since fa7342a but silently dropped by the
    # writer's fixed column list — pinned here so it cannot regress again.
    assert "over_max_intron" in header


def test_load_background_sv_bed_uses_col4_name(sv):
    flags = load_background_sv_bed(sv["bed"])
    assert flags.flag_of("chrA", *sv["j_sv"]) == "test_SV"
    assert flags.flag_of("chrA", *sv["j_ok"]) is None


def test_bundled_bed_covers_the_observed_srd1_junctions(tmp_path):
    """The bundled R64 track exists, is discovered for S288C genomes, and
    covers BOTH junction pairs observed in the wild (planning/730 W6:
    chrIII:148,615-151,517 and 149,817-151,668)."""
    fake_genome = tmp_path / "S288C_reference_sequence_test.fsa"
    fake_genome.write_text(">chrIII\nACGT\n")
    bundled = find_bundled_background_sv_bed(fake_genome)
    assert bundled is not None and bundled.exists()
    flags = load_background_sv_bed(bundled)
    assert flags.flag_of("chrIII", 148615, 151517) is not None
    assert flags.flag_of("chrIII", 149817, 151668) is not None
    # a junction elsewhere on chrIII is untouched
    assert flags.flag_of("chrIII", 200000, 200400) is None


def test_beside_genome_bed_wins_over_bundled(sv, tmp_path):
    g = tmp_path / "S288C_reference_sequence_test2.fsa"
    g.write_text(">chrIII\nACGT\n")
    local = tmp_path / "custom.background_sv.bed"
    local.write_text("chrI\t10\t20\tX\n")
    assert find_bundled_background_sv_bed(g) == local


def test_region_skip_shorthand_matches_the_bundled_interval():
    regions = parse_skip_regions("yeast-srd1-sv")
    assert regions == {"chrIII": [(148500, 151800)]}
    # the shorthand and the bundled BED must not drift apart
    bundled = find_bundled_background_sv_bed(
        Path("S288C_reference_sequence_x.fsa"))
    if bundled:  # bundled data present
        flags = load_background_sv_bed(bundled)
        chrom, span = YEAST_SRD1_SV_SPEC.rsplit(":", 1)
        a, b = (int(x) for x in span.split("-"))
        assert flags.flag_of(chrom, a, a + 1) is not None
        assert flags.flag_of(chrom, b - 1, b) is not None
