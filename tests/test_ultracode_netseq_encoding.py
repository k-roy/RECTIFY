"""Native NET-seq and mapped CPA calls must be invariant to SAM '=' encoding."""
import copy
import csv
import gzip
from dataclasses import asdict

import pysam
import pytest

from rectify.core.bam.netseq_bam_processor import process_netseq_bam, process_netseq_read
from rectify.core.netseq.netseq_rescue import JunctionPool, make_junction, revcomp
from rectify.core.netseq_cpa.pileup import walkback_pileup


def fixture(strand, kind, clipped=True):
    genome = list(("GCTCCGTACGTC" * 60)[:600])
    genome[120:122] = "GT"
    genome[258:270] = "AG" + "GTCTGTCGTC"
    genome[390:400] = "CGTCGAAAAA"
    genome = "".join(genome)
    start, end, clip = (80, 122, "CTGTCGTC") if kind == "rescue" else (360, 400, "AAA" if clipped else "")
    if strand == "-":
        genome = revcomp(genome)
    h = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": len(genome)}]})
    read = pysam.AlignedSegment(h)
    read.query_name = "same_molecule"
    read.reference_id = 0
    read.reference_start = start if strand == "+" else len(genome) - end
    read.mapping_quality = 60
    read.flag = 16 if strand == "+" else 0
    body = genome[read.reference_start:read.reference_start + end - start]
    if strand == "+":
        read.query_sequence = body + clip
        read.cigartuples = [(0, len(body))] + ([(4, len(clip))] if clip else [])
    else:
        read.query_sequence = revcomp(clip) + body
        read.cigartuples = ([(4, len(clip))] if clip else []) + [(0, len(body))]
    read.query_qualities = [31] * len(read.query_sequence)
    read.set_tag("MD", str(len(body)))
    read.set_tag("NM", 0)
    encoded = copy.deepcopy(read)
    bases = list(encoded.query_sequence)
    for q, r in encoded.get_aligned_pairs():
        if q is not None and r is not None:
            assert bases[q] == genome[r]
            bases[q] = "="
    encoded.query_sequence = "".join(bases)
    encoded.query_qualities = read.query_qualities[:]
    j = make_junction("chr1", 120, 260, "+") if strand == "+" else make_junction("chr1", 340, 480, "-")
    return read, encoded, genome, JunctionPool([j])


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("kind,clipped", [("tail", True), ("tail", False), ("rescue", True)])
def test_native_read_and_stream_preserve_equivalent_calls(tmp_path, strand, kind, clipped):
    literal, encoded, genome, pool = fixture(strand, kind, clipped)
    results = []
    for label, read in [("literal", literal), ("encoded", encoded)]:
        original = read.to_string()
        direct = process_netseq_read(read, "chr1", genome={"chr1": genome}, junction_pool=pool)
        assert read.to_string() == original
        bam = tmp_path / f"{label}.bam"
        with pysam.AlignmentFile(str(bam), "wb", header=read.header) as output:
            output.write(read)
        streamed = list(process_netseq_bam(str(bam), genome={"chr1": genome},
                                          junction_pool=pool, show_progress=False))
        assert len(streamed) == 1
        assert asdict(streamed[0]) == asdict(direct)
        results.append(asdict(direct))
    assert results[0] == results[1]
    if kind == "rescue":
        assert results[0]["rescue_k"] == 10
        assert results[0]["rescue_status"] == "spliced_rescued"
    else:
        assert results[0]["tail_walkback"] == (5 if clipped else 0)


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("clipped", [False, True])
def test_mapped_cpa_pileup_preserves_all_counts(tmp_path, strand, clipped):
    literal, encoded, genome, _ = fixture(strand, "tail", clipped)
    ref = tmp_path / "genome.fa"
    ref.write_text(">chr1\n" + genome + "\n")
    pysam.faidx(str(ref))
    results = []
    for label, read in [("literal", literal), ("encoded", encoded)]:
        bam, out = tmp_path / f"{label}.bam", tmp_path / f"{label}.tsv.gz"
        with pysam.AlignmentFile(str(bam), "wb", header=read.header) as output:
            output.write(read)
        before = bam.read_bytes()
        stats = walkback_pileup(bam, ref, out, sample="sample")
        assert bam.read_bytes() == before
        with gzip.open(out, "rt") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert stats["reads_used"] == 1 and len(rows) == 1
        assert int(rows[0]["pos"]) == (394 if strand == "+" else 205)
        results.append((stats, rows))
    assert results[0] == results[1]
