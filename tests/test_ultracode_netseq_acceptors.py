"""Alternative acceptors are scored as a family rather than erased by length."""
from dataclasses import asdict
from itertools import permutations

import pysam
import pytest

from rectify.core.bam.netseq_bam_processor import process_netseq_read
from rectify.core.netseq.netseq_rescue import JunctionPool, make_junction, revcomp
from rectify.utils.genome import register_genome_contigs


@pytest.fixture(autouse=True)
def restore_contig_registry():
    from rectify.utils import genome
    saved = set(genome._KNOWN_CONTIGS)
    try:
        yield
    finally:
        genome._KNOWN_CONTIGS.clear()
        genome._KNOWN_CONTIGS.update(saved)


def case(strand, first="CCCCCCCCCC", second="GACTGTCGTC", clip="GACTGTCGTC"):
    # Production CLI registers the loaded reference before parsing annotation.
    register_genome_contigs(["chr1"])
    g = list(("GCTCCGTACGTC" * 60)[:600])
    g[120:122], g[198:200], g[258:260] = "GT", "AG", "AG"
    g[200:210], g[260:270] = first, second
    g = "".join(g)
    if strand == "-":
        g = revcomp(g)
    h = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": len(g)}]})
    r = pysam.AlignedSegment(h)
    r.query_name, r.reference_id, r.mapping_quality = "alternative", 0, 60
    r.reference_start, r.flag = (80, 16) if strand == "+" else (480, 0)
    body = g[r.reference_start:r.reference_start + 40]
    r.query_sequence = body + clip if strand == "+" else revcomp(clip) + body
    r.cigarstring = f"40M{len(clip)}S" if strand == "+" else f"{len(clip)}S40M"
    r.query_qualities = [33] * len(r.query_sequence)
    coords = [(120, 200), (120, 260)] if strand == "+" else [(400, 480), (340, 480)]
    js = [make_junction("chr1", a, b, strand) for a, b in coords]
    return r, {"chr1": g}, js


@pytest.mark.parametrize("strand", ["+", "-"])
def test_farther_exact_acceptor_survives_annotation_order_and_duplicates(tmp_path, strand):
    read, genome, js = case(strand)
    expected = 269 if strand == "+" else 330
    seen = []
    for order in permutations(js):
        pool = JunctionPool([*order, order[0]])
        assert len(pool) == 2
        row = process_netseq_read(read, "chr1", junction_pool=pool, genome=genome)
        assert row.rescue_status == "spliced_rescued"
        assert row.three_prime_corrected == expected and row.rescue_k == 10
        assert row.rescue_n_candidates == 2
        assert row.to_dict()["rescue_n_candidates"] == 2
        seen.append(asdict(row))
    assert seen[0] == seen[1]
    # Same family through the actual annotation loader.
    gtf = tmp_path / "alternatives.gtf"
    lines = []
    for i, j in enumerate(js):
        for start, end in [(j.intron_start - 40, j.intron_start), (j.intron_end, j.intron_end + 40)]:
            lines.append(f'chr1\ttest\texon\t{start+1}\t{end}\t.\t{strand}\t.\tgene_id "g"; transcript_id "t{i}";')
    gtf.write_text("\n".join(lines) + "\n")
    pool = JunctionPool.from_annotation(gtf)
    assert len(pool) == 2
    assert process_netseq_read(read, "chr1", junction_pool=pool, genome=genome).three_prime_corrected == expected


@pytest.mark.parametrize("strand", ["+", "-"])
def test_two_supported_placements_remain_ambiguous(strand):
    read, genome, js = case(strand, first="GACTGTCGTC")
    row = process_netseq_read(read, "chr1", junction_pool=JunctionPool(js), genome=genome)
    assert row.rescue_status == "ambiguous"
    assert row.three_prime_corrected == row.three_prime_raw
    assert row.rescue_intron_start == row.rescue_intron_end == -1


@pytest.mark.parametrize("strand", ["+", "-"])
def test_larger_family_pays_for_chance_match_search(strand):
    read, genome, js = case(strand, clip="G")
    one = process_netseq_read(read, "chr1", junction_pool=JunctionPool([js[1]]), genome=genome)
    many = process_netseq_read(read, "chr1", junction_pool=JunctionPool(js), genome=genome)
    assert one.rescue_status == "spliced_rescued"  # existing single-candidate policy
    assert many.rescue_status == "exon1_end"       # one base cannot pay for two searches


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("k,expected", [(4, "exon1_end"), (5, "spliced_rescued")])
def test_randomer_floor_and_family_charge_both_apply(strand, k, expected):
    read, genome, js = case(strand, clip="GACTGTCGTC"[:k] + "AAAAAA")
    row = process_netseq_read(read, "chr1", junction_pool=JunctionPool(js), genome=genome, umi_length=6)
    assert row.rescue_status == expected
