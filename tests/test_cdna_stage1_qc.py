"""Stage-1 cDNA QC must be emitted on EVERY run, and must not depend on --workers.

Regression guard for a real production failure: ``_cdna_region_task`` computed the
read-type / XF / tail metrics per region and then discarded them, and the parent summed
only a fixed key list that excluded them. So the serial path printed a full QC block
while the parallel path — the path any real run takes, since it is the fast one — printed
almost none of it. Production libraries shipped with no read-type breakdown, and the
missing numbers had to be recomputed by hand from the output FASTQ.

The invariant these tests pin: for the same input, serial and parallel QC agree.
"""
from __future__ import annotations

import json

import pytest

from rectify.core.cdna.qc import (SHALLOW_READ_THRESHOLD, TYPE1_READ_FRAC_WARN,
                                  collect_qc, render_qc)


def _stats(**kw):
    base = dict(type1_reads=800, type2_reads=200, type1_clusters=400,
                type2_clusters=190, buckets_dropped_polyA_pileup=0,
                reads_in_dropped_buckets=0,
                subtype_reads={"umi_captured_fwd": 400, "umi_captured_rev": 400,
                               "umi_not_captured": 200})
    base.update(kw)
    return base


def _fq(**kw):
    base = dict(written=590, from_singletons=500, from_multi_pileup=90,
                from_multi_fallback=0, trim_frame_flipped=295,
                trim_frame_mismatch=0, trim_noop_5p=0, trim_noop_3p=60)
    base.update(kw)
    return base


def test_read_and_molecule_fractions_are_reported_separately():
    """The read-level and molecule-level Type-1 fractions are different quantities.

    Conflating them produced a bad QC threshold in real use: a molecule-level fraction was
    compared against the read-level ~82 % expectation and read as a failure.
    """
    qc = collect_qc(fastq_stats=_fq(), stats=_stats(), n_clusters=590,
                    tier_counts={0: 90, 1: 100, 2: 400}, tail_hist=[0] * 10,
                    n_input_reads=1000, workers=1)
    rt = qc["read_type"]
    assert rt["type1_read_frac"] == pytest.approx(0.8)
    assert rt["type1_cluster_frac"] == pytest.approx(400 / 590)
    assert rt["type1_read_frac"] != rt["type1_cluster_frac"]
    out = render_qc(qc)
    assert "READS" in out and "MOLECULES" in out


def test_shallow_library_suppresses_the_low_duplication_flag():
    """A near-zero duplication rate at low depth is expected, not a dedup failure.

    Without this, every small/failed-prep library trips a false alarm.
    """
    shallow = collect_qc(fastq_stats=_fq(written=990), stats=_stats(),
                         n_clusters=990, tier_counts={2: 990}, tail_hist=[0] * 10,
                         n_input_reads=1000, workers=1)
    assert shallow["shallow_library"] is True
    assert shallow["qc_depth_caveat"]
    assert not any(f.startswith("LOW_UMI_DUPLICATION") for f in shallow["flags"])

    deep_n = SHALLOW_READ_THRESHOLD * 2
    deep = collect_qc(fastq_stats=_fq(written=int(deep_n * 0.99)), stats=_stats(),
                      n_clusters=int(deep_n * 0.99), tier_counts={2: 10},
                      tail_hist=[0] * 10, n_input_reads=deep_n, workers=1)
    assert deep["shallow_library"] is False
    assert any(f.startswith("LOW_UMI_DUPLICATION") for f in deep["flags"])


def test_low_type1_read_fraction_is_flagged():
    qc = collect_qc(fastq_stats=_fq(), stats=_stats(type1_reads=600, type2_reads=400),
                    n_clusters=590, tier_counts={2: 590}, tail_hist=[0] * 10,
                    n_input_reads=1000, workers=1)
    assert qc["read_type"]["type1_read_frac"] < TYPE1_READ_FRAC_WARN
    assert any(f.startswith("LOW_TYPE1_READ_FRAC") for f in qc["flags"])


def test_qc_is_json_serialisable():
    """The sidecar is what a cross-sample PCA/heatmap pass consumes; console text is not."""
    qc = collect_qc(fastq_stats=_fq(), stats=_stats(), n_clusters=590,
                    tier_counts={2: 590}, tail_hist=[0] * 10,
                    n_input_reads=1000, workers=8)
    json.loads(json.dumps(qc))


def test_serial_and_parallel_shapes_agree():
    """Same numbers in, same QC out — whichever path assembled them.

    Serial passes per-cluster maps; parallel passes pre-aggregated counts. This is the
    exact divergence that let the parallel path ship without QC.
    """
    xf = {i: 2 for i in range(400)}
    xf.update({400 + i: 1 for i in range(100)})
    xf.update({500 + i: 0 for i in range(90)})
    tails = {i: 25 for i in range(400)}

    serial = collect_qc(fastq_stats=_fq(), stats=_stats(), n_clusters=590,
                        cluster_xf_tier=xf, cluster_tail_len=tails,
                        n_input_reads=1000, workers=1)
    parallel = collect_qc(fastq_stats=_fq(), stats=_stats(), n_clusters=590,
                          tier_counts=serial_tiers(serial), tail_hist=list(
                              serial["polya_tail_hist_xf2"].values()),
                          n_input_reads=1000, workers=8)
    for k in ("read_type", "subtype_reads", "umi_duplication_rate",
              "xf_tier_clusters", "polya_tail_hist_xf2", "input_reads",
              "output_molecules"):
        assert serial[k] == parallel[k], k


def serial_tiers(qc):
    x = qc["xf_tier_clusters"]
    return {0: x["xf0_not_detected"], 1: x["xf1_unanchored_medium"],
            2: x["xf2_anchored_high"]}
