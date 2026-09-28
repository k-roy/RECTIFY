"""Fragmented short-read RNA has no general poly(A)-tail trimming premise."""
import csv

import pysam
import pytest

from rectify.cli import main
from rectify.core.bam import bam_writer as bw


@pytest.mark.parametrize("streaming", [False, True])
def test_short_read_cli_preserves_genomic_fragment_ends_in_every_writer(tmp_path, streaming):
    body = "GCTACGTCCGATGCTACGTCGATCGTACCGTCGACTGCTC"
    genome, specs = "N" * 100, []
    for reverse in (False, True):
        for n in (0, 1, 5, 10):
            fragment = "T" * n + body if reverse else body + "A" * n
            start = len(genome)
            genome += fragment + "N" * 100
            for encoded in (False, True):
                specs.append((f"r{int(reverse)}_n{n}_e{int(encoded)}", reverse, fragment, start, encoded))
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chrI\n" + genome + "\n")
    header = pysam.AlignmentHeader.from_dict({"HD": {"SO": "coordinate"},
                                             "SQ": [{"SN": "chrI", "LN": len(genome)}]})
    bam, tsv, out = tmp_path / "input.bam", tmp_path / "corrected.tsv", tmp_path / "corrected.bam"
    truth = {}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as handle:
        for name, reverse, fragment, start, encoded in specs:
            r = pysam.AlignedSegment(header)
            r.query_name, r.reference_id, r.reference_start = name, 0, start
            r.flag, r.mapping_quality = 16 if reverse else 0, 60
            r.cigarstring = f"{len(fragment)}M"
            r.query_sequence = "=" * len(fragment) if encoded else fragment
            r.query_qualities = [35] * len(fragment)
            r.set_tag("MD", str(len(fragment)))
            r.set_tag("NM", 0)
            handle.write(r)
            truth[name] = (start, start + len(fragment), r.cigarstring, fragment, r.query_qualities[:])
    pysam.index(str(bam))
    args = ["correct", str(bam), "--genome", str(fasta), "--short-read", "--threads", "1",
            "-o", str(tsv), "--emit-merged-tsv", "--write-corrected-bam", str(out)]
    if streaming:
        args.append("--streaming")
    main(args)
    with tsv.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 16
    for row in rows:
        assert row["corrected_3prime"] == row["original_3prime"]
        assert row["tail_correction_enabled"] == "0"
        assert "polya_walkback" not in row["correction_applied"]
    loaded = bw._load_corrections_from_tsv(str(tsv))
    assert all(x["tail_correction_enabled"] is False for x in loaded.values())
    soft, dual_hard, dual_soft = [tmp_path / (x + ".bam") for x in ("soft", "dual_hard", "dual_soft")]
    bw.write_softclipped_bam(str(bam), str(tsv), str(soft), {"chrI": genome})
    bw.write_dual_bam(str(bam), str(tsv), str(dual_hard), str(dual_soft), {"chrI": genome})
    for path in (out, soft, dual_hard, dual_soft):
        seen = set()
        with pysam.AlignmentFile(str(path), "rb") as handle:
            for read in handle:
                seen.add(read.query_name)
                assert (read.reference_start, read.reference_end, read.cigarstring,
                        read.query_sequence, read.query_qualities) == truth[read.query_name]
        assert seen == set(truth)


def test_legacy_tsv_retains_existing_tail_behavior(tmp_path):
    path = tmp_path / "old.tsv"
    path.write_text("read_id\tstrand\tcorrected_3prime\nr\t+\t20\n")
    assert bw._load_corrections_from_tsv(str(path))["r"]["tail_correction_enabled"] is True


def test_tail_policy_column_keeps_every_legacy_tsv_position():
    from rectify.core.bam.output import CORRECTION_TSV_HEADER
    # Public schema at accepted e360bb9; positional downstream consumers must remain compatible.
    legacy = ['read_id', 'chrom', 'strand', 'original_3prime', 'corrected_3prime', 'five_prime_position', 'five_prime_rescued', 'five_prime_exon_cigar', 'alignment_start', 'alignment_end', 'ambiguity_min', 'ambiguity_max', 'ambiguity_range', 'polya_length', 'aligned_a_length', 'soft_clip_a_length', 'junctions', 'n_junctions', 'five_prime_soft_clip_length', 'three_prime_soft_clip_length', 'mapq', 'correction_applied', 'confidence', 'qc_flags', 'fraction', 'gene_id', 'pt_tag', 'polya_score', 'polya_source', 'sc_homopolymer_extension', 'sc_rescued_seq', 'sc_original_softclip_len', 'five_prime_intron_clip_pos', 'oc_homopolymer_extension', 'oc_overcall_count', 'oc_terminal_base', 'five_prime_upstream_trim', 'reanchor_clip_len', 'strand_evidence', 'consensus_aligner', 'consensus_confidence', 'consensus_n_agree', 'consensus_tied', 'five_prime_rescue_refused', 'five_prime_landing_annotated', 'five_prime_novel_evidence', 'five_prime_exon2_prefix', 'five_prime_exon_identity', 'five_prime_exon_bits', 'five_prime_clip_origin', 'five_prime_clip_origin_bits', 'five_prime_clip_prior_bits', 'five_prime_site_support', 'five_prime_landing_established', 'station_b_microexons', 'station_b_alternatives', 'station_b_n_tied', 'station_b_applied', 'station_b_intron_start', 'station_b_intron_end']
    # Appended since: five_prime_exon2_cigar (ISSUE-083 re-split, 2026-09-21). Every legacy position holds.
    assert CORRECTION_TSV_HEADER == legacy + ["tail_correction_enabled", "five_prime_exon2_cigar"]
