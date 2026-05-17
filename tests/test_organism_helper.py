"""
Tests for the shared --organism / --Scer helper.

Covers:
  - add_organism_args(parser_or_group) attaches both flags identically on a
    parser or argument group
  - --Scer is a store_const shorthand for --organism saccharomyces_cerevisiae
  - --Scer is wired into every subcommand that operates on a reference
    (align, analyze, aggregate, extract, validate, netseq, consensus,
     train-polya, prescan, split, export, correct, run-all, batch)
  - resolve_reference_paths(args) fills args.genome / args.annotation /
    args.go_annotations from the bundled organism when missing
  - resolve_reference_paths is a no-op when args lacks the organism slot
  - resolve_reference_paths preserves explicit --genome / --annotation
  - CLI hook in rectify.cli.main fills args.genome before dispatch for
    every non-bespoke subcommand
"""

import argparse
from pathlib import Path

import pytest


# Subcommands that should accept --Scer / --organism after Wave 3
SCER_COMMANDS = [
    ('align',       ['reads.fq', '-o', 'out/']),
    ('analyze',     ['in.tsv', '-o', 'out/']),
    ('aggregate',   ['in.bam', '-o', 'out/']),
    ('extract',     ['in.bam', '-o', 'out.tsv']),
    ('validate',    ['in.tsv', '-o', 'out.tsv']),
    ('netseq',      ['in.bam', '-o', 'out/']),
    ('consensus',   ['minimap2:a.bam', '-o', 'out/']),
    ('train-polya', ['in.bam', '--control-sites', 'c.tsv', '-o', 'm.json']),
    ('prescan',     ['--bam', 'in.bam', '-o', 'out/']),
    ('split',       ['reads.fq.gz', '-o', 'chunks/']),
    ('export',      ['in.tsv', '-o', 'out/']),
    ('correct',     ['in.bam', '-o', 'out.tsv']),
    ('run-all',     ['sample.bam', '-o', 'out/']),
    ('batch',       ['--manifest', 'm.tsv', '-o', 'out/']),
]


# ─────────────────────────────────────────────────────────────────────────────
# add_organism_args helper
# ─────────────────────────────────────────────────────────────────────────────

class TestAddOrganismArgs:

    def test_attaches_organism_and_scer_to_parser(self):
        from rectify.data import add_organism_args
        parser = argparse.ArgumentParser()
        add_organism_args(parser)
        dests = {a.dest for a in parser._actions}
        assert 'organism' in dests

    def test_scer_is_store_const_shorthand(self):
        from rectify.data import add_organism_args
        parser = argparse.ArgumentParser()
        add_organism_args(parser)
        args = parser.parse_args(['--Scer'])
        assert args.organism == 'saccharomyces_cerevisiae'

    def test_organism_default_none(self):
        from rectify.data import add_organism_args
        parser = argparse.ArgumentParser()
        add_organism_args(parser)
        args = parser.parse_args([])
        assert args.organism is None

    def test_organism_explicit_value_preserved(self):
        from rectify.data import add_organism_args
        parser = argparse.ArgumentParser()
        add_organism_args(parser)
        args = parser.parse_args(['--organism', 'yeast'])
        assert args.organism == 'yeast'

    def test_attaches_to_argument_group(self):
        """--Scer must work the same way when added to an argument group."""
        from rectify.data import add_organism_args
        parser = argparse.ArgumentParser()
        group = parser.add_argument_group('Reference')
        add_organism_args(group)
        args = parser.parse_args(['--Scer'])
        assert args.organism == 'saccharomyces_cerevisiae'


# ─────────────────────────────────────────────────────────────────────────────
# Per-command --Scer wiring (parametrized across the full inventory)
# ─────────────────────────────────────────────────────────────────────────────

class TestScerOnEveryCommand:

    @pytest.mark.parametrize('command,extra_argv', SCER_COMMANDS)
    def test_scer_sets_organism(self, command, extra_argv):
        from rectify.cli import create_parser
        parser = create_parser()
        args = parser.parse_args([command, '--Scer'] + extra_argv)
        assert args.organism == 'saccharomyces_cerevisiae', (
            f"--Scer must set organism on '{command}'"
        )

    @pytest.mark.parametrize('command,extra_argv', SCER_COMMANDS)
    def test_organism_yeast_accepted(self, command, extra_argv):
        from rectify.cli import create_parser
        parser = create_parser()
        args = parser.parse_args([command, '--organism', 'yeast'] + extra_argv)
        assert args.organism == 'yeast'


# ─────────────────────────────────────────────────────────────────────────────
# resolve_reference_paths
# ─────────────────────────────────────────────────────────────────────────────

class TestResolveReferencePaths:

    def test_no_organism_attr_is_noop(self):
        """A namespace without `organism` must not raise or mutate."""
        from rectify.data import resolve_reference_paths
        ns = argparse.Namespace(foo='bar')
        resolve_reference_paths(ns, require_genome=False)
        # No new attributes appeared
        assert not hasattr(ns, 'genome')
        assert not hasattr(ns, 'annotation')

    def test_fills_genome_and_annotation_from_organism(self):
        """--Scer alone should resolve to the bundled FASTA + GFF."""
        from rectify.data import resolve_reference_paths
        ns = argparse.Namespace(
            organism='saccharomyces_cerevisiae',
            genome=None,
            annotation=None,
        )
        resolve_reference_paths(ns, require_genome=False, verbose=False)
        assert ns.genome is not None
        assert Path(ns.genome).exists(), f"bundled genome missing: {ns.genome}"
        assert ns.annotation is not None
        assert Path(ns.annotation).exists(), f"bundled annotation missing: {ns.annotation}"

    def test_explicit_paths_take_priority(self, tmp_path):
        """If --genome is given, organism resolution must not overwrite it."""
        from rectify.data import resolve_reference_paths
        custom = tmp_path / 'custom.fa'
        custom.write_text('>chrI\nACGT\n')
        ns = argparse.Namespace(
            organism='saccharomyces_cerevisiae',
            genome=custom,
            annotation=None,
        )
        resolve_reference_paths(ns, require_genome=False, verbose=False)
        assert Path(ns.genome) == custom

    def test_fills_go_annotations_when_slot_exists(self):
        """When the command has a --go-annotations slot, organism fills it."""
        from rectify.data import resolve_reference_paths
        ns = argparse.Namespace(
            organism='saccharomyces_cerevisiae',
            genome=None,
            annotation=None,
            go_annotations=None,
        )
        resolve_reference_paths(ns, require_genome=False, verbose=False)
        assert ns.go_annotations is not None
        assert Path(ns.go_annotations).exists()

    def test_no_go_annotations_slot_ignored(self):
        """If args has no go_annotations attr, function must not set it."""
        from rectify.data import resolve_reference_paths
        ns = argparse.Namespace(
            organism='saccharomyces_cerevisiae',
            genome=None,
            annotation=None,
        )
        resolve_reference_paths(ns, require_genome=False, verbose=False)
        assert not hasattr(ns, 'go_annotations')

    def test_no_organism_no_genome_no_slot_returns(self):
        """With require_genome=False and no genome slot, the function returns silently."""
        from rectify.data import resolve_reference_paths
        ns = argparse.Namespace(organism=None)
        # Must not raise
        resolve_reference_paths(ns, require_genome=False, verbose=False)


# ─────────────────────────────────────────────────────────────────────────────
# CLI main() hook
# ─────────────────────────────────────────────────────────────────────────────

class TestCliMainHook:
    """
    The hook in rectify.cli.main() must fill args.genome from --Scer for every
    non-bespoke subcommand. We exercise it by calling parse_args and running
    the resolver the same way main() does, then asserting on the result.
    """

    def test_align_genome_filled_from_scer(self):
        from rectify.cli import create_parser
        from rectify.data import resolve_reference_paths
        parser = create_parser()
        args = parser.parse_args(['align', 'reads.fq', '--Scer', '-o', 'out/'])
        resolve_reference_paths(args, require_genome=False, verbose=False)
        assert args.genome is not None
        assert Path(args.genome).exists()

    @pytest.mark.parametrize('command,extra_argv', [
        ('align',     ['reads.fq', '-o', 'out/']),
        ('analyze',   ['in.tsv', '-o', 'out/']),
        ('aggregate', ['in.bam', '-o', 'out/']),
        ('consensus', ['minimap2:a.bam', '-o', 'out/']),
        ('extract',   ['in.bam', '-o', 'out.tsv']),
        ('netseq',    ['in.bam', '-o', 'out/']),
    ])
    def test_genome_resolved_for_non_bespoke_commands(self, command, extra_argv):
        from rectify.cli import create_parser
        from rectify.data import resolve_reference_paths
        parser = create_parser()
        args = parser.parse_args([command, '--Scer'] + extra_argv)
        resolve_reference_paths(args, require_genome=False, verbose=False)
        assert args.genome is not None, f"{command} should auto-fill --genome from --Scer"
        assert Path(args.genome).exists()
