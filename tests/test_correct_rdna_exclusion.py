"""Module 2F rRNA/Pol III exclusion: a read mapping into an excluded region
(rDNA, tRNA, SNR6, ...) must SKIP the expensive 3'SS rescue, while an identical
read outside any excluded region still runs it.

These high-coverage non-mRNA loci are dropped downstream by `analyze`; running
the per-read HP-edit-distance rescue over them dominated correct-stage wall time
(2026-05-25). The gate lives in `correct_read_3prime`; the detector is threaded
in from `correct_command` via the parallel worker state.
"""

import pysam

import rectify.core.bam.bam_processor as bp
from rectify.core.exclusion_regions import ExclusionRegionDetector


def _make_read(*, chrom="chrI", start=1000, seq="A" * 10 + "C" * 40,
               cigar=((4, 10), (0, 40))):
    """A 5' soft-clipped read (Gate-A evidence) that reaches the 3'SS rescue."""
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": 2_000_000}]}
    )
    r = pysam.AlignedSegment(hdr)
    r.query_name = "excl_test"
    r.reference_name = chrom
    r.reference_start = start
    r.cigartuples = list(cigar)
    r.is_reverse = False
    r.is_unmapped = False
    r.is_secondary = False
    r.is_supplementary = False
    r.mapping_quality = 60
    r.query_sequence = seq
    return r


def _detector_covering(chrom, start, end, reason="rDNA"):
    det = ExclusionRegionDetector(flanking_bp=0)
    det.add_region(chrom, start, end, reason, "test_locus")
    return det


def _run(read, detector, monkeypatch):
    """Call correct_read_3prime with a spy on the 3'SS rescue; return call count."""
    calls = {"n": 0}

    def _spy(*a, **k):
        calls["n"] += 1
        return {"rescued": False}

    monkeypatch.setattr(bp, "_rescue_3ss", _spy)
    genome = {"chrI": "C" * 2_000_000}
    # annotated_junctions present so the rescue block's data-gate is satisfied;
    # the only thing that should suppress the call is the exclusion gate.
    bp.correct_read_3prime(
        read, genome,
        apply_3ss_rescue=True,
        annotated_junctions={("chrI", 980, 1000)},
        exclusion_detector=detector,
    )
    return calls["n"]


def test_rescue_skipped_inside_excluded_region(monkeypatch):
    read = _make_read(start=1000)
    det = _detector_covering("chrI", 500, 1500)  # read 5' (1000) is inside
    assert _run(read, det, monkeypatch) == 0, \
        "3'SS rescue must be skipped for a read inside an excluded (rDNA/Pol III) region"


def test_rescue_runs_outside_excluded_region(monkeypatch):
    read = _make_read(start=1000)
    det = _detector_covering("chrI", 50_000, 60_000)  # far from the read
    assert _run(read, det, monkeypatch) >= 1, \
        "3'SS rescue must still run for a read outside any excluded region"


def test_no_detector_runs_rescue(monkeypatch):
    read = _make_read(start=1000)
    assert _run(read, None, monkeypatch) >= 1, \
        "with no exclusion detector, rescue runs as before (default behaviour)"
