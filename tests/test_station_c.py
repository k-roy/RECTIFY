"""Station C v0 (pool-gate) — verdict logic, census, flags, CLI wiring.

Synthetic mini-genome with planted junctions exercises every verdict class;
the bundled validation BAM exercises the real CLI path end-to-end.
"""

import json
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from rectify.core.consensus.station_c import (
    PoolGateConfig,
    find_bundled_selfhom_bed,
    load_annotated_canonical,
    load_repeat_flags,
    load_selfhom_bed,
    pool_gate,
    write_pool_gate_outputs,
)
from rectify.core.splice.overhang_informativeness import canonical_in_class

# Two high-complexity 60-mers (no long homopolymers, no GT/AG at the seams).
EX1 = "ACGTGATCCATGCTTACGCTGACTATCGGACTTCAGATCCGTACTGACGATCCATGCATC"
EX2 = "TTCAGCATGGACGTATCACGGATCAGCTACGTTAGCATCGACTGCATCGGATACGCATGC"
FILL = "CATGACTGCATCAGGATCCATCAGCATGACGTATCAGCATGGATCCATCACGGATCATGC"


def _mk_intron(donor: str, acceptor: str, n_fill: int = 76) -> str:
    return donor + (FILL * 3)[:n_fill] + acceptor


@pytest.fixture(scope="module")
def mini(tmp_path_factory):
    """Genome + GFF + selfhom BED + BAM with one junction per verdict class."""
    d = tmp_path_factory.mktemp("station_c")

    # chrA layout: EX1 | intron GT..AG (annotated) | EX2 | EX1 |
    #              intron GT..AG (novel canonical) | EX2 | EX1 |
    #              intron CC..TT (novel non-canon) | EX2
    i_can = _mk_intron("GT", "AG")
    i_non = _mk_intron("CC", "TT")
    segs = [EX1, i_can, EX2, EX1, i_can, EX2, EX1, i_non, EX2]
    chrA = "".join(segs)
    # chrB: the same non-canonical junction, inside an rRNA_gene flag
    chrB = EX1 + i_non + EX2

    off = [0]
    for s in segs:
        off.append(off[-1] + len(s))
    j_ann = (off[1], off[2])       # annotated canonical
    j_can = (off[4], off[5])       # novel canonical
    j_non = (off[7], off[8])       # novel non-canonical
    j_b = (len(EX1), len(EX1) + len(i_non))

    # the construction must actually give the classes the tests assume
    assert canonical_in_class(chrA, *j_can)
    assert not canonical_in_class(chrA, *j_non)
    assert not canonical_in_class(chrB, *j_b)

    fa = d / "mini.fa"
    fa.write_text(f">chrA\n{chrA}\n>chrB\n{chrB}\n")

    gff = d / "mini.gff"
    gff.write_text(
        "##gff-version 3\n"
        f"chrA\tsgd\tintron\t{j_ann[0] + 1}\t{j_ann[1]}\t.\t+\t.\tID=i1\n"
        f"chrB\tsgd\trRNA_gene\t1\t{len(chrB)}\t.\t+\t.\tID=r1\n"
    )

    bed = d / "mini.selfhomology.bed"
    bed.write_text("# test\n")  # empty track: annotation flag does the work

    header = {"HD": {"VN": "1.6"},
              "SQ": [{"SN": "chrA", "LN": len(chrA)},
                     {"SN": "chrB", "LN": len(chrB)}]}
    bam = d / "mini.bam"

    def read(name, chrom, start, exon1, intron_len, exon2, seq=None):
        a = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(header))
        a.query_name = name
        a.reference_id = 0 if chrom == "chrA" else 1
        a.reference_start = start
        a.cigartuples = [(0, len(exon1)), (3, intron_len), (0, len(exon2))]
        a.query_sequence = seq if seq is not None else exon1 + exon2
        a.mapping_quality = 60
        return a

    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        # annotated junction, 2 reads
        for n in range(2):
            out.write(read(f"ann{n}", "chrA", off[0], EX1, len(i_can), EX2))
        # novel canonical, 2 clean reads -> admit_candidate
        for n in range(2):
            out.write(read(f"can{n}", "chrA", off[3], EX1, len(i_can), EX2))
        # novel non-canonical, 2 clean reads, unflagged -> admit_candidate
        for n in range(2):
            out.write(read(f"non{n}", "chrA", off[6], EX1, len(i_non), EX2))
        # non-canonical inside the rRNA_gene flag -> demote
        out.write(read("flg0", "chrB", 0, EX1, len(i_non), EX2))
        # '='-compressed SEQ variant of the flagged read: must decode, not
        # zero out (the 644h smoke lesson)
        out.write(read("flg1", "chrB", 0, EX1, len(i_non), EX2,
                       seq="=" * (len(EX1) + len(EX2))))

    genome = {"chrA": chrA, "chrB": chrB}
    return {"dir": d, "fa": fa, "gff": gff, "bed": bed, "bam": bam,
            "genome": genome, "j_ann": j_ann, "j_can": j_can,
            "j_non": j_non, "j_b": j_b}


def _rows_by_junction(rows):
    return {(r["chrom"], r["start"], r["end"]): r for r in rows}


def test_pool_gate_verdicts(mini):
    rows, summary = pool_gate(str(mini["bam"]), mini["genome"], mini["gff"],
                              selfhom_bed=mini["bed"])
    by = _rows_by_junction(rows)

    # annotated junction is counted, not reported
    assert summary["n_annotated"] == 1
    assert ("chrA", *mini["j_ann"]) not in by

    can = by[("chrA", *mini["j_can"])]
    assert can["canonical_in_class"] == 1
    assert can["support"] == 2
    assert can["q_max"] >= 40
    assert can["verdict"] == "admit_candidate"

    non = by[("chrA", *mini["j_non"])]
    assert non["canonical_in_class"] == 0
    assert non["support"] == 2
    assert non["q_max"] >= 80  # clean 60-mer overhangs
    assert non["verdict"] == "admit_candidate"

    flg = by[("chrB", *mini["j_b"])]
    assert flg["repeat_flag"] == "rRNA_gene"
    assert flg["verdict"] == "demote_orthogonal_evidence"
    # the '='-compressed read decoded against the reference: q stays high
    assert flg["q_max"] >= 80
    assert flg["support"] == 2


def test_low_support_canonical_is_review(mini):
    cfg = PoolGateConfig(min_support=3)  # push both novel junctions below
    rows, _ = pool_gate(str(mini["bam"]), mini["genome"], mini["gff"],
                        cfg=cfg, selfhom_bed=mini["bed"])
    by = _rows_by_junction(rows)
    assert by[("chrA", *mini["j_can"])]["verdict"] == "review"
    assert by[("chrA", *mini["j_non"])]["verdict"] == "review"
    # flagged stays demoted regardless of support config
    assert by[("chrB", *mini["j_b"])]["verdict"] == "demote_orthogonal_evidence"


def test_selfhom_bed_flags_demote(mini):
    bed = mini["dir"] / "wide.selfhomology.bed"
    j = mini["j_non"]
    bed.write_text(f"chrA\t{j[0] - 10}\t{j[1] + 10}\n")
    rows, summary = pool_gate(str(mini["bam"]), mini["genome"], mini["gff"],
                              selfhom_bed=bed)
    by = _rows_by_junction(rows)
    assert by[("chrA", *mini["j_non"])]["selfhom_flag"] == 1
    assert by[("chrA", *mini["j_non"])]["verdict"] == "demote_orthogonal_evidence"
    # canonical track is never demoted by flags (644i: flags gate the
    # non-canonical track only)
    assert by[("chrA", *mini["j_can"])]["verdict"] == "admit_candidate"
    assert summary["selfhom_intervals"] == 1


def test_outputs_roundtrip(mini, tmp_path):
    rows, summary = pool_gate(str(mini["bam"]), mini["genome"], mini["gff"],
                              selfhom_bed=mini["bed"])
    tsv, js = write_pool_gate_outputs(rows, summary, tmp_path / "s1")
    lines = tsv.read_text().rstrip("\n").split("\n")
    assert len(lines) == 1 + len(rows)
    assert lines[0].split("\t")[0] == "chrom"
    loaded = json.loads(js.read_text())
    assert loaded["verdicts"] == summary["verdicts"]


def test_annotated_loader_counts(mini):
    ann = load_annotated_canonical(mini["gff"], mini["genome"],
                                   PoolGateConfig())
    assert ("chrA", *mini["j_ann"]) in ann


def test_repeat_flag_loader(mini):
    flags = load_repeat_flags(mini["gff"], margin=0)
    assert flags.flag_of("chrB", 5, 25) == "rRNA_gene"
    assert flags.flag_of("chrA", 5, 25) is None


def test_bundled_selfhom_track_exists():
    import rectify
    genome = (Path(rectify.__file__).parent / "data" / "genomes"
              / "saccharomyces_cerevisiae"
              / "S288C_reference_sequence_R64-5-1_20240529.fsa.gz")
    bed = find_bundled_selfhom_bed(genome)
    assert bed is not None and bed.exists()
    flags = load_selfhom_bed(bed)
    assert flags.n_intervals > 100
    # the rDNA array must be covered (the 644i headline case)
    assert flags.flag_of("chrXII", 455000, 465000) == "selfhom"


def test_cli_smoke_on_bundled_validation_bam(tmp_path):
    import rectify
    bam = (Path(rectify.__file__).parent / "data" / "validation" / "aligners"
           / "validation_reads.minimap2.bam")
    out = tmp_path / "pg"
    r = subprocess.run(
        [sys.executable, "-m", "rectify", "pool-gate", str(bam), "--Scer",
         "-o", str(out)],
        capture_output=True, text=True, timeout=600,
    )
    assert r.returncode == 0, r.stderr[-2000:]
    tsv = out.with_suffix(".pool_gate.tsv")
    assert tsv.exists()
    assert "Station C (pool-gate):" in r.stdout


def test_runall_parser_has_station_flags():
    from rectify.cli import create_parser
    parser = create_parser()
    args = parser.parse_args(
        ["run-all", "x.fastq.gz", "--Scer", "-o", "out", "--triage",
         "--triage-clip-legs", "--no-pool-gate"])
    assert args.triage and args.triage_clip_legs and args.no_pool_gate
    args2 = parser.parse_args(["pool-gate", "in.bam", "--Scer", "-o", "p"])
    assert args2.command == "pool-gate"


def test_triage_parser_clip_legs_flag():
    from rectify.cli import create_parser
    parser = create_parser()
    args = parser.parse_args(
        ["triage", "in.bam", "-o", "out", "--Scer", "--clip-legs"])
    assert args.clip_legs
