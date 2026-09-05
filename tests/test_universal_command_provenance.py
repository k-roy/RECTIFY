#!/usr/bin/env python3
"""EVERY rectify subcommand must leave a record of which code ran.

Audit 2026-08-02: only 5 of 26 command modules wrote a provenance sidecar, so
most rectify operations produced outputs with no version record at all. The fix
hooks the single CLI dispatch point rather than hand-wiring each module, so a
command added tomorrow is covered without anyone remembering to do it. These
tests guard that property.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pytest

from rectify.core.provenance import auto_sidecar as aus


class TestOutputDirResolution:
    @pytest.mark.parametrize("attr", ["output_dir", "outdir", "output", "out"])
    def test_common_output_arg_names(self, attr, tmp_path):
        ns = argparse.Namespace(**{attr: str(tmp_path)})
        assert aus.resolve_output_dir(ns) == tmp_path

    def test_a_file_path_resolves_to_its_parent(self, tmp_path):
        ns = argparse.Namespace(output=str(tmp_path / "calls.tsv"))
        assert aus.resolve_output_dir(ns) == tmp_path

    def test_no_output_arg_returns_none(self):
        assert aus.resolve_output_dir(argparse.Namespace(verbose=True)) is None


class TestCommandSidecar:
    def test_writes_and_pins_the_version(self, tmp_path):
        ns = argparse.Namespace(output_dir=str(tmp_path), threads=8, verbose=False)
        p = aus.write_command_sidecar("split", ns, "2026-08-02T00:00:00Z", 1.5, 0)
        assert p is not None and p.exists()
        d = json.loads(p.read_text())
        assert d["stage"] == "split"
        assert d["exit_status"] == 0
        # THE point of the whole exercise:
        assert d["rectify_git_sha"] != "unknown"
        assert d["rectify_git_sha"].startswith(("git:", "pkghash:"))
        # the knobs that shaped the run are captured
        assert d["args"]["threads"] == 8

    def test_records_a_failure_too(self, tmp_path):
        ns = argparse.Namespace(output_dir=str(tmp_path))
        p = aus.write_command_sidecar("aggregate", ns, "2026-08-02T00:00:00Z", 0.1, 2)
        d = json.loads(p.read_text())
        assert d["exit_status"] == 2, "a failed run must still be attributable"

    def test_non_data_commands_are_skipped(self, tmp_path):
        ns = argparse.Namespace(output_dir=str(tmp_path))
        for cmd in sorted(aus.SKIP_COMMANDS):
            assert aus.write_command_sidecar(cmd, ns, "t", 0.0, 0) is None
        assert not list(tmp_path.glob("*.json"))

    def test_skip_list_stays_small(self):
        """Guard against the skip list quietly becoming the escape hatch."""
        assert aus.SKIP_COMMANDS <= {"test", "install-aligners", "verify"}

    def test_never_raises_on_a_bad_output_dir(self, tmp_path, monkeypatch):
        # The unwritable output_dir forces the Path.cwd() fallback (see
        # auto_sidecar.write_command_sidecar) — chdir into tmp_path first so
        # that fallback write lands in the test's own temp dir instead of
        # wherever pytest happened to be invoked from (e.g. the repo root).
        monkeypatch.chdir(tmp_path)
        ns = argparse.Namespace(output_dir="/proc/nonexistent/cannot/create")
        # Must degrade, never propagate — provenance may not break a command.
        p = aus.write_command_sidecar("export", ns, "t", 0.0, 0)
        assert p == tmp_path / "export.command.provenance.json"
        assert p.exists()

    def test_write_is_atomic_no_tmp_left_behind(self, tmp_path):
        ns = argparse.Namespace(output_dir=str(tmp_path))
        aus.write_command_sidecar("prescan", ns, "t", 0.0, 0)
        assert not list(tmp_path.glob("*.tmp")), "partial file left on disk"


class TestCliHook:
    def test_main_emits_a_sidecar_for_a_module_with_no_provenance_of_its_own(
        self, tmp_path, monkeypatch
    ):
        """The hook is at the dispatch point, so coverage does not depend on the
        command module knowing anything about provenance."""
        import rectify.cli as cli

        called = {}

        def fake_dispatch(args, parser):
            called["yes"] = True          # a command that writes no sidecar

        monkeypatch.setattr(cli, "_dispatch", fake_dispatch)
        monkeypatch.setattr(
            cli, "create_parser",
            lambda: _parser_returning(argparse.Namespace(
                command="split", output_dir=str(tmp_path))),
        )
        cli.main([])
        assert called.get("yes")
        sidecars = list(tmp_path.glob("split.command.provenance.json"))
        assert sidecars, "dispatch hook did not write a command sidecar"
        assert json.loads(sidecars[0].read_text())["rectify_git_sha"] != "unknown"

    def test_nonzero_exit_is_still_recorded(self, tmp_path, monkeypatch):
        import rectify.cli as cli

        def fake_dispatch(args, parser):
            raise SystemExit(3)

        monkeypatch.setattr(cli, "_dispatch", fake_dispatch)
        monkeypatch.setattr(
            cli, "create_parser",
            lambda: _parser_returning(argparse.Namespace(
                command="export", output_dir=str(tmp_path))),
        )
        with pytest.raises(SystemExit):
            cli.main([])
        d = json.loads((tmp_path / "export.command.provenance.json").read_text())
        assert d["exit_status"] == 3


def _parser_returning(ns):
    class _P:
        def parse_args(self, argv=None):
            return ns

        def print_help(self):
            pass
    return _P()
