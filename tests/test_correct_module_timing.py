"""Per-module wall-time counter in `correct` (Kevin, 2026-09-05).

Station A (the resolver) and Module 2F both ask "does this overhang belong
across a junction?"; before deciding whether a resolver-verdict skip in 2F is
worth a cross-arm cache, the 2F / 2H / walkback / indel share of the correct
stage has to be MEASURED. These tests pin the plumbing:

* `correct_read_3prime` returns `module_seconds` on its primary row — one entry
  per per-read module plus `read_total` — and the map is NOT a TSV column.
* `ProcessingStats` sums the map over primary rows only, merges across regions,
  and round-trips through the region-checkpoint JSON beside the int counters.
* The stats TSV and the text report surface the rows.
"""

import json
import time

import pysam

import rectify.core.bam.bam_processor as bp
from rectify.core.bam import parallel as bam_parallel
from rectify.core.bam.output import CORRECTION_TSV_HEADER, correction_result_to_tsv_row
from rectify.core.bam.processing_stats import (
    MODULE_DESCRIPTIONS,
    PER_READ_MODULES,
    ProcessingStats,
    format_module_seconds,
    generate_stats_report,
    module_seconds_rows,
    write_stats_tsv,
)


def _make_read(*, chrom="chrI", start=1000, seq="A" * 10 + "C" * 40,
               cigar=((4, 10), (0, 40))):
    """A 5' soft-clipped read that reaches the 3'SS rescue (Module 2F)."""
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": 2_000_000}]}
    )
    r = pysam.AlignedSegment(hdr)
    r.query_name = "timing_test"
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


_GENOME = {"chrI": "C" * 2_000_000}


def _run(read, **kw):
    rows = bp.correct_read_3prime(
        read, _GENOME,
        apply_3ss_rescue=True,
        annotated_junctions={("chrI", 980, 1000)},
        **kw,
    )
    assert len(rows) == 1
    return rows[0]


def test_primary_row_carries_every_module_timer():
    row = _run(_make_read())
    ms = row['module_seconds']
    for name in PER_READ_MODULES + ('read_total',):
        assert name in ms, name
        assert ms[name] >= 0.0, name
    # Every module timer is a slice of the whole-function clock.
    assert ms['read_total'] >= max(v for k, v in ms.items() if k != 'read_total')
    assert set(ms) <= set(MODULE_DESCRIPTIONS)


def test_2f_timer_measures_the_rescue_call(monkeypatch):
    """The 2F entry is the wall time INSIDE `_rescue_3ss`, not an estimate."""
    def _slow(*a, **k):
        time.sleep(0.02)
        return {"rescued": False}

    monkeypatch.setattr(bp, "_rescue_3ss", _slow)
    ms = _run(_make_read())['module_seconds']
    assert ms['five_prime_rescue_2f'] >= 0.015
    assert ms['read_total'] >= ms['five_prime_rescue_2f']


def test_2f_timer_is_zero_when_the_module_is_skipped():
    row = bp.correct_read_3prime(_make_read(), _GENOME, apply_3ss_rescue=False)[0]
    assert row['module_seconds']['five_prime_rescue_2f'] == 0.0


def test_module_seconds_is_not_a_tsv_column():
    """The corrected-TSV schema (44 columns, checkpoint-guarded) is untouched."""
    assert 'module_seconds' not in CORRECTION_TSV_HEADER
    row = _run(_make_read())
    assert len(correction_result_to_tsv_row(row)) == len(CORRECTION_TSV_HEADER)


def test_processing_stats_sums_primary_rows_only():
    s = ProcessingStats()
    s.update_from_result({'module_seconds': {'five_prime_rescue_2f': 1.0, 'read_total': 2.0}})
    s.update_from_result({'module_seconds': {'five_prime_rescue_2f': 5.0, 'read_total': 9.0},
                          'is_primary_result': False})   # NET-seq split row
    s.update_from_result({'module_seconds': {'five_prime_rescue_2f': 0.5, 'read_total': 1.0}})
    s.update_from_result({})                              # chimeric stub: no map
    assert s.module_seconds == {'five_prime_rescue_2f': 1.5, 'read_total': 3.0}
    assert s.reads_processed == 3


def test_processing_stats_merge_and_add():
    a = ProcessingStats(); a.add_module_seconds('junction_refinement_2h', 4.0)
    b = ProcessingStats(); b.add_module_seconds('junction_refinement_2h', 1.0)
    b.add_module_seconds('bam_processing', 7.0)
    a.merge(b)
    assert a.module_seconds == {'junction_refinement_2h': 5.0, 'bam_processing': 7.0}
    # to_dict() stays an int-only map (the checkpoint reader int-casts it).
    assert 'module_seconds' not in a.to_dict()


def test_region_checkpoint_json_round_trips_module_seconds():
    s = ProcessingStats()
    s.reads_processed = 3
    s.add_module_seconds('five_prime_rescue_2f', 0.25)
    s.add_module_seconds('read_total', 1.5)
    text = json.dumps(bam_parallel._stats_to_checkpoint_dict(s), sort_keys=True)
    back = bam_parallel._stats_from_dict(json.loads(text))
    assert back.reads_processed == 3
    assert back.module_seconds == {'five_prime_rescue_2f': 0.25, 'read_total': 1.5}
    # A pre-timing checkpoint (no map) still loads.
    old = bam_parallel._stats_from_dict({'reads_processed': 2})
    assert old.reads_processed == 2 and old.module_seconds == {}


def test_rows_report_per_read_share_of_read_total():
    rows = module_seconds_rows({
        'read_total': 10.0, 'five_prime_rescue_2f': 2.5, 'junction_refinement_2h': 4.0,
    })
    by_name = {name: (secs, pct) for name, secs, pct, _ in rows}
    assert by_name['five_prime_rescue_2f'] == (2.5, '25.00')
    assert by_name['junction_refinement_2h'] == (4.0, '-')      # stage wall, no share
    assert by_name['read_total'] == (10.0, '-')
    assert 'five_prime_rescue_2f 2.5s (25.00%)' in format_module_seconds(
        {'read_total': 10.0, 'five_prime_rescue_2f': 2.5})


def test_stats_tsv_and_report_surface_module_seconds(tmp_path):
    s = ProcessingStats()
    s.total_reads_in_bam = 1
    s.reads_processed = 1
    s.confidence_high = 1
    s.add_module_seconds('read_total', 4.0)
    s.add_module_seconds('five_prime_rescue_2f', 1.0)
    s.add_module_seconds('junction_refinement_2h', 2.0)
    out = tmp_path / 'stats.tsv'
    write_stats_tsv(s, str(out))
    lines = [l for l in out.read_text().splitlines() if l.startswith('module_seconds_')]
    assert 'module_seconds_five_prime_rescue_2f\t1.000\t25.00\t' + \
        MODULE_DESCRIPTIONS['five_prime_rescue_2f'] in lines
    assert any(l.startswith('module_seconds_junction_refinement_2h\t2.000\t-\t') for l in lines)
    report = generate_stats_report(s)
    assert 'Module wall time' in report and 'five_prime_rescue_2f' in report
    # No rows and no section when nothing was recorded.
    empty = tmp_path / 'empty.tsv'
    write_stats_tsv(ProcessingStats(), str(empty))
    assert 'module_seconds_' not in empty.read_text()
    assert 'Module wall time' not in generate_stats_report(ProcessingStats())
