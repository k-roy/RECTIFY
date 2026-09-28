"""ISSUE-048: placement-relative SEQ must never become a transferable payload."""

from pathlib import Path

import pysam
import pytest

from rectify.core.multialign import build_cma, expand


def records(reverse=False, hardclip=False):
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": "chr1", "LN": 400}],
    })
    donor = pysam.AlignedSegment(header)
    donor.query_name = "physical_ACGA"
    donor.reference_id = 0
    donor.reference_start = 100
    donor.cigarstring = "4M"
    donor.query_sequence = "===="
    donor.query_qualities = [20, 21, 22, 23]
    donor.mapping_quality = 60
    alternate = pysam.AlignedSegment.fromstring(donor.to_string(), header)
    alternate.reference_start = 200
    alternate.flag = 16 if reverse else 0
    alternate.query_sequence = "TCGT" if reverse else "ACGA"
    alternate.query_qualities = [23, 22, 21, 20] if reverse else [20, 21, 22, 23]
    if hardclip:
        seq, qual = alternate.query_sequence, alternate.query_qualities
        alternate.cigarstring = "1H3M"
        alternate.query_sequence = seq[1:]
        alternate.query_qualities = qual[1:]
    genome = {"chr1": "N" * 100 + "ACGA" + "N" * 96 + "ACGG" + "N" * 196}
    return header, {"minimap2": donor, "deSALT": alternate}, genome


@pytest.mark.parametrize("genome", [None, {}, {"chr1": "A" * 102}])
def test_undecodable_payload_refuses_without_output(tmp_path, genome):
    header, arms, _ = records()
    target = tmp_path / "unsafe.cma.bam"
    with pytest.raises(ValueError, match="physical_ACGA"):
        build_cma([("physical_ACGA", arms)], header, str(target), arms, genome=genome)
    assert not target.exists()


@pytest.mark.parametrize("cigar,seq", [("1S3M", "=CGA"), ("1I3M", "=CGA")])
def test_unaligned_equals_cannot_be_decoded_even_with_reference(tmp_path, cigar, seq):
    header, arms, genome = records()
    donor = arms["minimap2"]
    donor.cigarstring, donor.query_sequence = cigar, seq
    with pytest.raises(ValueError, match="physical_ACGA"):
        build_cma([("physical_ACGA", arms)], header, str(tmp_path / "bad.bam"), arms,
                  genome=genome)


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("hardclip", [False, True])
@pytest.mark.parametrize("encoded", [False, True])
def test_expansion_preserves_bases_and_qualities(tmp_path, reverse, hardclip, encoded):
    header, arms, genome = records(reverse, hardclip)
    if not encoded:
        donor = arms["minimap2"]
        donor.query_sequence = "ACGA"
        donor.query_qualities = [20, 21, 22, 23]
    target = tmp_path / "safe.cma.bam"
    build_cma([("physical_ACGA", arms)], header, str(target), arms,
              genome=genome if encoded else None)
    result = dict(expand(str(target)))["physical_ACGA"]
    assert result["minimap2"].query_sequence == "ACGA"
    for name, original in arms.items():
        emitted = result[name]
        assert emitted.query_qualities == original.query_qualities
        assert (emitted.flag, emitted.reference_start, emitted.cigarstring) == (
            original.flag, original.reference_start, original.cigarstring)
    assert result["deSALT"].query_sequence == arms["deSALT"].query_sequence


def test_late_failure_preserves_existing_output_and_cleans_temporary(tmp_path):
    header, unsafe, _ = records()
    _, safe, _ = records()
    safe["minimap2"].query_sequence = "ACGA"
    target = tmp_path / "existing.bam"
    target.write_bytes(b"previous completed output")
    before = set(tmp_path.iterdir())
    with pytest.raises(ValueError):
        build_cma([("a", safe), ("b", unsafe)], header, str(target), unsafe)
    assert target.read_bytes() == b"previous completed output"
    assert set(tmp_path.iterdir()) == before


def test_actual_cli_accepts_reference_and_refuses_unsafe_build(tmp_path, capsys):
    from rectify.cli import main

    header, arms, genome = records()
    inputs = []
    for arm, record in arms.items():
        path = tmp_path / f"{arm}.bam"
        with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
            bam.write(record)
        inputs.append(f"{arm}={path}")
    out = tmp_path / "out.cma.bam"
    argv = ["cma", "build", "--aligner-bams", *inputs, "--out", str(out)]
    with pytest.raises(SystemExit) as failed:
        main(argv)
    assert failed.value.code == 1
    assert not out.exists()
    captured = capsys.readouterr()
    assert "--genome" in captured.err + captured.out
    reference = tmp_path / "genome.fa"
    reference.write_text(">chr1\n" + genome["chr1"] + "\n")
    pysam.faidx(str(reference))
    with pytest.raises(SystemExit) as passed:
        main(argv + ["--genome", str(reference)])
    assert passed.value.code == 0
    assert {r.query_sequence for r in dict(expand(str(out)))["physical_ACGA"].values()} == {"ACGA"}


def test_scer_parser_resolves_bundled_reference():
    from rectify.cli import create_parser
    from rectify.data import resolve_reference_paths

    args = create_parser().parse_args([
        "cma", "build", "--aligner-bams", "minimap2=x.bam", "--out", "x.cma.bam", "--Scer"])
    resolve_reference_paths(args, require_genome=False, verbose=False)
    assert Path(args.genome).is_file()
    assert "saccharomyces_cerevisiae" in str(args.genome)
