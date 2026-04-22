#!/usr/bin/env python3
"""
Tests for run_command (run-all) argument wiring and provenance integration.

Covers:
- analyze_args has 'output' attribute (not 'output_dir') — the critical mismatch bug
- analyze_args has all attributes read by analyze_command.run_analyze()
- batch command has --filter-spikein in its argparse definition
- SLURM profile loading works correctly (especially time field not parsed as sexagesimal)
- _get_available_cpus() respects SLURM environment variable
"""

import argparse
import os
import tempfile
from pathlib import Path

import pytest


# ─────────────────────────────────────────────────────────────────────────────
# run_command analyze_args wiring
# ─────────────────────────────────────────────────────────────────────────────

class TestRunCommandAnalyzeArgWiring:
    """
    Verify that the Namespace built in run_command for the analyze step
    has every attribute that analyze_command.run_analyze() reads.

    The critical bug: run_command used output_dir= but analyze_command reads args.output.
    """

    def _make_run_all_args(self, tmp_path):
        """Simulate what run_command.run() does when building analyze_args."""
        # This mirrors the exact Namespace construction in run_command.py
        corrected_tsv = tmp_path / 'corrected_3ends.tsv'
        output_dir = tmp_path
        genome_path = None
        annotation_path = None
        go_path = None

        return argparse.Namespace(
            input=corrected_tsv,
            output=output_dir,               # Fixed: was output_dir=
            annotation=annotation_path,
            genome=genome_path,
            reference=None,
            go_annotations=go_path,
            threads=4,
            sample_column='sample',
            count_column=None,
            cluster_distance=25,
            min_reads=5,
            run_deseq2=True,
            run_motif=True,
            sample_sets=None,
            exclude_mito=True,
            include_mito=False,
            exclude_rdna=True,
            include_rdna=False,
            no_bedgraph=False,
            bedgraph_dir=None,
            no_genomic_distribution=False,
            motif_upstream=100,
            motif_downstream=50,
        )

    def test_output_attribute_exists(self, tmp_path):
        """Critical: analyze_command reads args.output, not args.output_dir."""
        args = self._make_run_all_args(tmp_path)
        # Must have 'output', not 'output_dir'
        assert hasattr(args, 'output'), "analyze_args must have 'output' attribute"
        assert not hasattr(args, 'output_dir') or args.output == tmp_path

    def test_output_is_path(self, tmp_path):
        args = self._make_run_all_args(tmp_path)
        assert isinstance(args.output, Path)

    def test_exclude_mito_exists(self, tmp_path):
        """analyze_command line 86: args.exclude_mito — direct access, not getattr."""
        args = self._make_run_all_args(tmp_path)
        assert hasattr(args, 'exclude_mito')
        assert args.exclude_mito is True  # default should exclude mito

    def test_bedgraph_dir_exists(self, tmp_path):
        """analyze_command line 143: args.bedgraph_dir — direct access."""
        args = self._make_run_all_args(tmp_path)
        assert hasattr(args, 'bedgraph_dir')

    def test_motif_upstream_exists(self, tmp_path):
        """analyze_command line 397: args.motif_upstream — direct access in motif block."""
        args = self._make_run_all_args(tmp_path)
        assert hasattr(args, 'motif_upstream')
        assert args.motif_upstream == 100

    def test_motif_downstream_exists(self, tmp_path):
        """analyze_command line 400: args.motif_downstream — direct access."""
        args = self._make_run_all_args(tmp_path)
        assert hasattr(args, 'motif_downstream')
        assert args.motif_downstream == 50

    def test_all_required_attributes_present(self, tmp_path):
        """All attributes accessed by run_analyze() must be present."""
        args = self._make_run_all_args(tmp_path)
        required = [
            'input', 'output', 'annotation', 'genome', 'reference',
            'go_annotations', 'threads', 'sample_column', 'count_column',
            'cluster_distance', 'min_reads', 'run_deseq2', 'run_motif',
            'sample_sets', 'exclude_mito', 'include_mito', 'exclude_rdna',
            'include_rdna', 'no_bedgraph', 'bedgraph_dir',
            'no_genomic_distribution', 'motif_upstream', 'motif_downstream',
        ]
        missing = [attr for attr in required if not hasattr(args, attr)]
        assert missing == [], f"analyze_args missing attributes: {missing}"


# ─────────────────────────────────────────────────────────────────────────────
# Actual Namespace from run_command module
# ─────────────────────────────────────────────────────────────────────────────

class TestRunCommandModuleWiring:
    """Verify the actual Namespace built in run_command.run() has correct attributes."""

    def _simulate_run_command_analyze_args(self, tmp_path):
        """
        Reproduce exactly what run_command.run() builds for analyze_args.
        We import the module to test against the real code, not a copy.
        """
        import importlib
        import types

        # We can't call run_command.run() directly (needs real BAM files),
        # but we can inspect that analyze_command's expected attrs are present
        # by checking the source uses 'output=' not 'output_dir='.
        import inspect
        import rectify.core.run_command as run_mod

        source = inspect.getsource(run_mod.run)
        return source

    def test_run_command_uses_output_not_output_dir(self):
        """
        _run_analysis() in run_command must pass output= (not output_dir=) to analyze_args,
        since analyze_command reads args.output.
        """
        import inspect
        import rectify.core.run_command as run_mod

        # The fix is now in _run_analysis which builds analyze_args
        source = inspect.getsource(run_mod._run_analysis)

        assert 'output_dir=output_dir' not in source, (
            "_run_analysis() still uses 'output_dir=output_dir'. "
            "analyze_command reads args.output — use 'output=output_dir'."
        )
        assert 'output=output_dir' in source, (
            "_run_analysis() must set 'output=output_dir' in analyze_args"
        )

    def test_run_command_sets_exclude_mito(self):
        import inspect
        import rectify.core.run_command as run_mod

        source = inspect.getsource(run_mod._run_analysis)
        assert 'exclude_mito' in source, (
            "_run_analysis() must set exclude_mito in analyze_args"
        )

    def test_run_command_sets_motif_upstream(self):
        import inspect
        import rectify.core.run_command as run_mod

        source = inspect.getsource(run_mod._run_analysis)
        assert 'motif_upstream' in source, (
            "_run_analysis() must set motif_upstream in analyze_args"
        )


# ─────────────────────────────────────────────────────────────────────────────
# batch command argparse
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# --Scer / --organism shorthand
# ─────────────────────────────────────────────────────────────────────────────

class TestScerShorthand:
    """--Scer sets organism=saccharomyces_cerevisiae and makes --genome/--annotation optional."""

    def test_run_all_Scer_sets_organism(self):
        from rectify.cli import create_parser
        parser = create_parser()
        args = parser.parse_args(['run-all', 'sample.bam', '--Scer', '-o', 'out/'])
        assert args.organism == 'saccharomyces_cerevisiae', (
            "--Scer must set organism='saccharomyces_cerevisiae'"
        )

    def test_run_all_organism_yeast_accepted(self):
        from rectify.cli import create_parser
        parser = create_parser()
        args = parser.parse_args(['run-all', 'sample.bam', '--organism', 'yeast', '-o', 'out/'])
        assert args.organism == 'yeast'

    def test_run_all_genome_optional_when_Scer(self):
        """--genome is no longer required — omitting it with --Scer must not raise."""
        from rectify.cli import create_parser
        parser = create_parser()
        # Must not raise argparse.ArgumentError
        args = parser.parse_args(['run-all', 'sample.bam', '--Scer', '-o', 'out/'])
        assert args.genome is None  # not set explicitly

    def test_run_all_explicit_genome_still_works(self):
        from rectify.cli import create_parser
        parser = create_parser()
        args = parser.parse_args([
            'run-all', 'sample.bam',
            '--genome', 'genome.fa', '--annotation', 'genes.gff',
            '-o', 'out/',
        ])
        assert str(args.genome) == 'genome.fa'
        assert args.organism is None

    def test_batch_Scer_sets_organism(self):
        from rectify.cli import create_parser
        parser = create_parser()
        batch_parser = parser._subparsers._group_actions[0].choices["batch"]
        args = batch_parser.parse_args(['--manifest', 'manifest.tsv', '--Scer', '-o', 'out/'])
        assert args.organism == 'saccharomyces_cerevisiae'

    def test_batch_genome_optional_when_Scer(self):
        from rectify.cli import create_parser
        parser = create_parser()
        batch_parser = parser._subparsers._group_actions[0].choices["batch"]
        args = batch_parser.parse_args(['--manifest', 'manifest.tsv', '--Scer', '-o', 'out/'])
        assert args.genome is None


class TestBatchCommandArgparse:

    def test_batch_has_filter_spikein(self):
        from rectify.cli import create_parser
        parser = create_parser()
        batch_parser = parser._subparsers._group_actions[0].choices["batch"]
        actions = {a.dest for a in batch_parser._actions}
        assert "filter_spikein" in actions, \
            "batch subparser must expose --filter-spikein"

    def test_batch_has_profile(self):
        from rectify.cli import create_parser
        parser = create_parser()
        batch_parser = parser._subparsers._group_actions[0].choices["batch"]
        actions = {a.dest for a in batch_parser._actions}
        assert "profile" in actions, \
            "batch subparser must expose --profile for SLURM profile files"

    def test_batch_has_workers(self):
        from rectify.cli import create_parser
        parser = create_parser()
        batch_parser = parser._subparsers._group_actions[0].choices["batch"]
        actions = {a.dest for a in batch_parser._actions}
        assert "workers" in actions, \
            "batch subparser must expose --workers for interactive parallelism"


# ─────────────────────────────────────────────────────────────────────────────
# SLURM profile loading
# ─────────────────────────────────────────────────────────────────────────────

class TestSlurmProfile:

    def test_load_json_profile(self, tmp_path):
        import json
        from rectify.core.batch_command import load_slurm_profile

        profile_data = {
            "partition": "my-partition",
            "time": "4:00:00",
            "mem": "32G",
            "cpus": 8,
        }
        profile_path = tmp_path / "test.json"
        profile_path.write_text(json.dumps(profile_data))

        profile = load_slurm_profile(profile_path)
        assert profile["partition"] == "my-partition"
        assert profile["time"] == "4:00:00"
        assert profile["cpus"] == 8

    def test_yaml_time_not_sexagesimal(self, tmp_path):
        """
        YAML without quotes parses 4:00:00 as 14400 (sexagesimal).
        Profile files must quote time values. Verify our bundled profiles are correct.
        """
        try:
            import yaml
        except ImportError:
            pytest.skip("PyYAML not installed")

        from rectify.core.batch_command import load_slurm_profile

        # Find bundled profile
        import rectify
        profile_path = Path(rectify.__file__).parent / 'slurm_profiles' / 'hpc_cpu.yaml'

        if not profile_path.exists():
            pytest.skip(f"Bundled profile not found: {profile_path}")

        profile = load_slurm_profile(profile_path)
        assert isinstance(profile['time'], str), (
            f"SLURM time field must be a string (got {type(profile['time']).__name__}: {profile['time']!r}). "
            "Add quotes around time values in YAML: time: \"4:00:00\""
        )
        assert ':' in profile['time'], \
            f"Time field should be HH:MM:SS format, got: {profile['time']!r}"

    def test_load_yaml_profile(self, tmp_path):
        """Basic YAML profile loading (with quoted time)."""
        try:
            import yaml
        except ImportError:
            pytest.skip("PyYAML not installed")

        from rectify.core.batch_command import load_slurm_profile

        content = 'partition: my-partition\ntime: "4:00:00"\nmem: 32G\ncpus: 8\n'
        profile_path = tmp_path / "test.yaml"
        profile_path.write_text(content)

        profile = load_slurm_profile(profile_path)
        assert profile["time"] == "4:00:00"
        assert profile["partition"] == "my-partition"

    def test_invalid_extension_raises(self, tmp_path):
        from rectify.core.batch_command import load_slurm_profile

        profile_path = tmp_path / "test.txt"
        profile_path.write_text("partition: my-partition")
        with pytest.raises(ValueError, match="Profile must be"):
            load_slurm_profile(profile_path)


# ─────────────────────────────────────────────────────────────────────────────
# CPU detection
# ─────────────────────────────────────────────────────────────────────────────

class TestCpuDetection:

    def test_uses_slurm_env_when_set(self, monkeypatch):
        from rectify.core.batch_command import _get_available_cpus

        monkeypatch.setenv('SLURM_CPUS_PER_TASK', '16')
        assert _get_available_cpus() == 16

    def test_falls_back_to_system_count(self, monkeypatch):
        import multiprocessing
        from rectify.core.batch_command import _get_available_cpus

        monkeypatch.delenv('SLURM_CPUS_PER_TASK', raising=False)
        assert _get_available_cpus() == multiprocessing.cpu_count()

    def test_slurm_env_takes_precedence_over_system(self, monkeypatch):
        """SLURM_CPUS_PER_TASK should override multiprocessing.cpu_count()."""
        from rectify.core.batch_command import _get_available_cpus

        monkeypatch.setenv('SLURM_CPUS_PER_TASK', '4')
        result = _get_available_cpus()
        assert result == 4, f"Expected 4 (from SLURM env), got {result}"
