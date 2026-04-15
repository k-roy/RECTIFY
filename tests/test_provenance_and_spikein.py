#!/usr/bin/env python3
"""
Tests for provenance tracking and spike-in filter integration.

Covers:
- ProvenanceTracker saves PROVENANCE.json and README.md correctly
- ProcessingStats is JSON-serializable via dataclasses.asdict()
- SpikeInFilter ENO2 known signature is loadable
- correct_command config includes filter_spikein and filter_spikein propagates
- CLI exposes --filter-spikein on 'correct' and 'run-all' is present
"""

import dataclasses
import json
import sys
import tempfile
import argparse
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))


# ─────────────────────────────────────────────────────────────────────────────
# Provenance
# ─────────────────────────────────────────────────────────────────────────────

class TestProvenanceTracker:

    def test_save_creates_json_and_readme(self, tmp_path):
        from rectify.utils.provenance import ProvenanceTracker
        tracker = ProvenanceTracker(tmp_path, description="test run")
        tracker.set_command(["rectify", "correct", "foo.bam", "-o", "out.tsv"])
        tracker.set_config({"bam_path": "foo.bam", "threads": 4})
        tracker.save()

        assert (tmp_path / "PROVENANCE.json").exists()
        assert (tmp_path / "README.md").exists()

    def test_json_is_valid(self, tmp_path):
        from rectify.utils.provenance import ProvenanceTracker
        tracker = ProvenanceTracker(tmp_path)
        tracker.set_command(["rectify"])
        tracker.save()

        data = json.loads((tmp_path / "PROVENANCE.json").read_text())
        assert "runs" in data
        assert len(data["runs"]) == 1
        assert data["runs"][0]["command"] == "rectify"

    def test_processing_stats_serializable_via_asdict(self):
        """The bug: ProcessingStats passed directly to json.dumps fails.
        The fix: dataclasses.asdict() must succeed."""
        from rectify.core.processing_stats import ProcessingStats
        stats = ProcessingStats()
        stats.spikein_reads_filtered = 12345

        # Direct serialization must fail (confirms the bug existed)
        with pytest.raises(TypeError):
            json.dumps(stats)

        # asdict serialization must succeed (confirms the fix works)
        d = dataclasses.asdict(stats)
        serialized = json.dumps(d)
        assert '"spikein_reads_filtered": 12345' in serialized

    def test_add_output_file_with_stats_metadata(self, tmp_path):
        """Provenance.save() must not raise when stats metadata is a dict."""
        from rectify.utils.provenance import ProvenanceTracker
        from rectify.core.processing_stats import ProcessingStats

        stats = ProcessingStats()
        stats.reads_processed = 9999

        tracker = ProvenanceTracker(tmp_path)
        # Write a dummy file so add_output_file doesn't skip it
        dummy = tmp_path / "corrected_3ends.tsv"
        dummy.write_text("read_id\tcorrected_3prime\n")

        tracker.add_output_file(
            dummy,
            source_files=[Path("input.bam")],
            metadata={"stats": dataclasses.asdict(stats)},  # the fix
        )
        tracker.save()  # must not raise

        data = json.loads((tmp_path / "PROVENANCE.json").read_text())
        assert data["runs"][0]["outputs"][0]["metadata"]["stats"]["reads_processed"] == 9999

    def test_incremental_runs_appended(self, tmp_path):
        from rectify.utils.provenance import ProvenanceTracker
        for i in range(3):
            t = ProvenanceTracker(tmp_path, description="multi-run")
            t.set_command([f"run-{i}"])
            t.save()

        data = json.loads((tmp_path / "PROVENANCE.json").read_text())
        assert len(data["runs"]) == 3


# ─────────────────────────────────────────────────────────────────────────────
# Spike-in filter
# ─────────────────────────────────────────────────────────────────────────────

class TestSpikeInFilter:

    def test_eno2_known_signature_exists(self):
        from rectify.core.spikein_filter import KNOWN_SPIKEIN_SIGNATURES
        assert "ENO2" in KNOWN_SPIKEIN_SIGNATURES
        sig = KNOWN_SPIKEIN_SIGNATURES["ENO2"]
        assert sig.get("systematic_name") == "YHR174W"

    def test_add_known_signature(self):
        from rectify.core.spikein_filter import SpikeInFilter
        f = SpikeInFilter()
        f.add_known_signature("ENO2")
        assert len(f.signatures) == 1
        assert f.signatures[0]["gene"] == "ENO2"

    def test_unknown_gene_not_added(self):
        from rectify.core.spikein_filter import SpikeInFilter
        f = SpikeInFilter()
        f.add_known_signature("NOTREAL")
        assert len(f.signatures) == 0


# ─────────────────────────────────────────────────────────────────────────────
# correct_command config
# ─────────────────────────────────────────────────────────────────────────────

class TestCorrectCommandConfig:

    def _make_args(self, **kwargs):
        defaults = dict(
            input=Path("/dev/null"),
            genome=None,
            output=None,
            annotation=None,
            dT_primed_cDNA=False,
            skip_atract_check=False,
            skip_ag_check=False,
            skip_polya_trim=False,
            skip_indel_correction=False,
            skip_variant_aware=False,
            polya_model=None,
            organism=None,
            netseq_dir=None,
            netseq_samples=None,
            threads=4,
            verbose=False,
            filter_spikein=None,
        )
        defaults.update(kwargs)
        return argparse.Namespace(**defaults)

    def test_filter_spikein_defaults_none(self):
        from rectify.core.correct_command import validate_inputs
        args = self._make_args(input=Path(__file__))
        config = validate_inputs(args)
        assert config["filter_spikein"] is None

    def test_filter_spikein_propagated(self):
        from rectify.core.correct_command import validate_inputs
        args = self._make_args(input=Path(__file__), filter_spikein=["ENO2"])
        config = validate_inputs(args)
        assert config["filter_spikein"] == ["ENO2"]


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

class TestCLI:

    def test_run_all_command_registered(self):
        from rectify.cli import create_parser
        parser = create_parser()
        choices = parser._subparsers._group_actions[0].choices
        assert "run-all" in choices

    def test_run_command_removed(self):
        from rectify.cli import create_parser
        parser = create_parser()
        choices = parser._subparsers._group_actions[0].choices
        assert "run" not in choices

    def test_correct_has_filter_spikein(self):
        from rectify.cli import create_parser
        parser = create_parser()
        correct_parser = parser._subparsers._group_actions[0].choices["correct"]
        actions = {a.dest for a in correct_parser._actions}
        assert "filter_spikein" in actions

    def test_filter_spikein_parsed(self):
        from rectify.cli import create_parser
        parser = create_parser()
        with tempfile.NamedTemporaryFile(suffix=".bam") as bam, \
             tempfile.NamedTemporaryFile(suffix=".fa") as genome:
            args = parser.parse_args([
                "correct", bam.name,
                "--genome", genome.name,
                "-o", "/tmp/out.tsv",
                "--filter-spikein", "ENO2",
            ])
        assert args.filter_spikein == ["ENO2"]

    def test_filter_spikein_multiple_genes(self):
        from rectify.cli import create_parser
        parser = create_parser()
        with tempfile.NamedTemporaryFile(suffix=".bam") as bam, \
             tempfile.NamedTemporaryFile(suffix=".fa") as genome:
            args = parser.parse_args([
                "correct", bam.name,
                "--genome", genome.name,
                "-o", "/tmp/out.tsv",
                "--filter-spikein", "ENO2", "ACT1",
            ])
        assert args.filter_spikein == ["ENO2", "ACT1"]
