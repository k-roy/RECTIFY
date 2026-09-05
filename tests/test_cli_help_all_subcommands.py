"""`--help` must not crash, for any subcommand (ISSUE-014).

argparse %-expands every ``help=`` string against the action's own attributes,
so a bare ``%`` in prose is a format conversion. ``'99% of'`` is read as the
``'% o'`` octal conversion and ``rectify align --help`` died with

    TypeError: %o format: an integer is required, not dict

That is invisible until someone actually asks for help — no import, no test and
no run touches the formatter — and it hid the resolver post-pass flags from the
one person trying to discover them. Percent signs in help text must be written
``%%``.

The loop below is the point: a single broken string anywhere in the CLI fails
this file, so the class of bug is pinned rather than the one instance.
"""

import argparse

import pytest

from rectify.cli import create_parser


def _subparsers(parser, path=()):
    """Yield (path, parser) for the root parser and every subparser, at any
    nesting depth."""
    yield path, parser
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            seen = set()
            for name, sub in action.choices.items():
                if id(sub) in seen:      # aliases share one parser object
                    continue
                seen.add(id(sub))
                yield from _subparsers(sub, path + (name,))


ALL = list(_subparsers(create_parser()))
IDS = [' '.join(p) or '(root)' for p, _ in ALL]


def test_the_parser_actually_has_subcommands():
    """Anti-vacuity guard: if create_parser() ever stops registering
    subcommands, the parametrised test below would pass with nothing to check."""
    assert len(ALL) > 5, f'only {len(ALL)} parsers discovered — the sweep is vacuous'
    assert any('align' in p for p, _ in ALL)


@pytest.mark.parametrize('path,parser', ALL, ids=IDS)
def test_format_help_does_not_raise(path, parser):
    text = parser.format_help()
    assert text, f'{" ".join(path)} produced empty help'


@pytest.mark.parametrize('path,parser', ALL, ids=IDS)
def test_format_usage_does_not_raise(path, parser):
    assert parser.format_usage()


def test_align_help_renders_the_percentages_literally():
    """The specific regression: the escaped signs must survive as one '%'
    each, not disappear or double."""
    align = next(p for path, p in ALL if path and path[-1] == 'align')
    text = align.format_help()
    assert '99% of' in text.replace('\n', ' ').replace('  ', ' ') or '99%' in text
    assert '%%' not in text, 'an escape leaked into the rendered help'


def test_help_exits_zero_not_with_a_traceback():
    """End to end, the way a user hits it."""
    parser = create_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(['align', '--help'])
    assert exc.value.code == 0
