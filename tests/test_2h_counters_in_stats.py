#!/usr/bin/env python3
"""Module 2H's decision counters reach the stats TSV (T0 4f9102d F3, ported 2026-10-13 from the stranded 018 branch).

The refiner's stats dict (n_op_reads, refined, unchanged, unrealizable_winner_skipped, ...) was a log line only. It is now
copied into ProcessingStats.module_2h_counters and written as one `module_2h_<key>` row per key right after the
module_2h_failed marker (ISSUE-025), and the parallel driver's region-JSON round trip keeps both fields' types (it
int()-coerced every field, so a non-empty module_2h_failed would have raised).
"""
from pathlib import Path


class TestModule2hCountersInStats:
    COUNTERS = {
        "n_op_reads": 20, "refined": 7, "noncanon_destination_refused": 3,
        "noncanon_destination_refused_alternative_within_delta": 3,
        "unrealizable_winner_skipped": 1, "noncanon_dest_lost_to_annotated_alt": 2,
    }

    @staticmethod
    def _rows(path):
        rows, order = {}, []
        for line in Path(path).read_text().splitlines():
            if line.startswith("#") or line.startswith("metric\t") or not line.strip():
                continue
            parts = line.split("\t")
            rows[parts[0]] = parts
            order.append(parts[0])
        return rows, order

    def test_rows_follow_module_2h_failed_with_dash_percent(self, tmp_path):
        from rectify.core.bam.processing_stats import ProcessingStats, write_stats_tsv
        stats = ProcessingStats()
        stats.total_reads_in_bam = 10
        stats.module_2h_failed = "RuntimeError: boom"
        stats.module_2h_counters = dict(self.COUNTERS)
        out = tmp_path / "x_stats.tsv"
        write_stats_tsv(stats, str(out))
        rows, order = self._rows(out)
        for k, n in self.COUNTERS.items():
            assert rows[f"module_2h_{k}"][1:3] == [str(n), "-"], k
        keys = [f"module_2h_{k}" for k in sorted(self.COUNTERS)]
        assert order[:1 + len(keys)] == ["module_2h_failed"] + keys
        assert order[1 + len(keys)] == "total_reads_in_bam"

    def test_no_rows_when_2h_did_not_run(self, tmp_path):
        from rectify.core.bam.processing_stats import ProcessingStats, write_stats_tsv
        out = tmp_path / "y_stats.tsv"
        write_stats_tsv(ProcessingStats(), str(out))
        rows, _ = self._rows(out)
        assert not [k for k in rows if k.startswith("module_2h_")]

    def test_merge_sums_and_parallel_round_trip_keeps_types(self):
        import dataclasses
        import json
        from rectify.core.bam import parallel as par
        from rectify.core.bam.processing_stats import ProcessingStats
        a, b = ProcessingStats(), ProcessingStats()
        a.module_2h_counters = {"refined": 1, "unrealizable_winner_skipped": 1}
        b.module_2h_counters = {"refined": 2}
        b.module_2h_failed = "ValueError: x"
        a.merge(b)
        assert a.module_2h_counters == {"refined": 3, "unrealizable_winner_skipped": 1}
        assert a.module_2h_failed == "ValueError: x"
        back = par._stats_from_dict(json.loads(json.dumps(a.to_dict())))
        assert back.module_2h_counters == a.module_2h_counters
        assert back.module_2h_failed == "ValueError: x"        # was int()-coerced before F3
        json.dumps(dataclasses.asdict(a))                        # provenance path still serializes
