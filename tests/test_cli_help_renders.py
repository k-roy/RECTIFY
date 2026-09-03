"""Every subcommand's ``--help`` must render.

argparse %-formats every help string, so a literal percent sign has to be written ``%%``. A bare
``%`` in `rectify align`'s aligner help raised ``TypeError: %o format: an integer is required, not
dict`` -- i.e. `rectify align --help` was UNREACHABLE, and three independent audits
(planning 830 G2, 831, 832) each rediscovered it. One parametrized test closes the whole class.
"""
from __future__ import annotations

import argparse

import pytest

from rectify.cli import create_parser


def _subparsers():
    parser = create_parser()
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action.choices
    raise AssertionError("rectify CLI has no subparsers")


SUBCOMMANDS = sorted(_subparsers())


def test_there_are_subcommands():
    assert len(SUBCOMMANDS) > 10


@pytest.mark.parametrize("name", SUBCOMMANDS)
def test_subcommand_help_renders(name):
    _subparsers()[name].format_help()


def test_top_level_help_renders():
    create_parser().format_help()
