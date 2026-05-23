"""Tests for poly(A)-artifact false-junction detection.

Regression coverage for NEW-080: the genomic *target* flank used to score the
"junction to genomic A-tract" artifact signal must mirror the read-side flank.
Plus strand reads align poly-A over a genomic A-tract AFTER the N ([end, end+20));
minus strand reads align poly-T over a genomic T-run BEFORE the N ([start-20, start)).
Before the fix, the minus-strand target was fetched from [end, end+20) (the wrong
side of the N) and the artifact was missed whenever the canonical-motif fallback
did not also fire.
"""

import pysam

from rectify.core.splice.false_junction_filter import (
    _GenomeDictReference,
    filter_polya_artifact_junctions,
)


def _make_read(*, chrom="chrI", start=1000, seq="", cigar=((0, 0),), is_reverse=False):
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": 100_000}]}
    )
    r = pysam.AlignedSegment(hdr)
    r.query_name = "test"
    r.reference_name = chrom
    r.reference_start = start
    r.cigartuples = list(cigar)
    r.is_reverse = is_reverse
    r.is_unmapped = False
    r.is_secondary = False
    r.is_supplementary = False
    r.mapping_quality = 60
    r.query_sequence = seq
    return r


def test_minus_strand_artifact_uses_upstream_genomic_flank():
    """NEW-080: minus-strand target flank is [start-20, start), not [end, end+20).

    Genome layout (chrI), all other bases neutral 'C':
      [990, 1010)  = T*20   upstream flank — the genomic T-run the read's poly-T
                            aligns over (the TRUE artifact target)
      [1010, 1012) = CT     minus-strand canonical donor (so has_canonical_motif)
      [1028, 1030) = AC     minus-strand canonical acceptor
      [1030, 1050) = G*20   downstream flank — NOT T-rich (the WRONG target)

    The junction is made canonical on purpose so the `not has_canonical_motif`
    fallback does NOT fire — detection then depends solely on target_a_richness,
    which isolates the flank-selection fix. Pre-fix the target is read from the
    G-rich downstream flank (t_richness 0.0 -> not an artifact); post-fix from the
    T-rich upstream flank (t_richness 1.0 -> artifact).
    """
    chrom = "chrI"
    g = ["C"] * 2000
    for i in range(990, 1010):
        g[i] = "T"
    g[1010], g[1011] = "C", "T"
    g[1028], g[1029] = "A", "C"
    for i in range(1030, 1050):
        g[i] = "G"
    ref = _GenomeDictReference({chrom: "".join(g)})

    # M10 over [1000,1010) (poly-T) | N20 [1010,1030) | M20 [1030,1050) (exon)
    read = _make_read(
        chrom=chrom,
        start=1000,
        seq="T" * 10 + "C" * 20,
        cigar=((0, 10), (3, 20), (0, 20)),
        is_reverse=True,
    )

    real, artifacts = filter_polya_artifact_junctions(read, ref, strand="-")

    assert len(artifacts) == 1, "minus-strand poly-A artifact junction should be detected"
    assert not real
    a = artifacts[0]
    assert a.has_canonical_motif, "fixture is canonical; detection must come from the target flank"
    assert a.target_a_richness >= 0.8, "target must be scored on the T-rich upstream flank, not the G-rich downstream flank"


def test_plus_strand_artifact_still_uses_downstream_flank():
    """Guard: the plus-strand path is unchanged ([end, end+20), A-richness)."""
    chrom = "chrI"
    g = ["C"] * 2000
    g[1020], g[1021] = "G", "T"  # canonical donor
    g[1038], g[1039] = "A", "G"  # canonical acceptor
    for i in range(1040, 1060):
        g[i] = "A"  # genomic A-tract after the N
    ref = _GenomeDictReference({chrom: "".join(g)})

    # M20 over [1000,1020) (exon) | N20 [1020,1040) | M10 [1040,1050) (poly-A)
    read = _make_read(
        chrom=chrom,
        start=1000,
        seq="C" * 20 + "A" * 10,
        cigar=((0, 20), (3, 20), (0, 10)),
        is_reverse=False,
    )

    real, artifacts = filter_polya_artifact_junctions(read, ref, strand="+")

    assert len(artifacts) == 1, "plus-strand poly-A artifact junction should be detected"
    assert artifacts[0].target_a_richness >= 0.8
