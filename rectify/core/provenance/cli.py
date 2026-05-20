"""Shared argparse helpers for resume flags. Used by every stage parser
that wants to honor ``--force-*`` / ``--accept-prior-provenance`` /
``--dry-run-resume``.
"""

import argparse


def add_resume_args(parser: argparse.ArgumentParser) -> None:
    """Add the resume/sidecar flags to ``parser``.

    Flags added:

    - ``--force-all`` — ignore all sidecars; rerun everything.
    - ``--force-stage NAME[,NAME...]`` — ignore sidecars for these stages only.
    - ``--accept-prior-provenance`` — treat git_sha mismatch as non-blocking.
    - ``--dry-run-resume`` — print SKIP/RUN decisions and exit.
    """
    g = parser.add_argument_group("Resume / sidecar")
    g.add_argument(
        "--force-all",
        action="store_true",
        default=False,
        help="Ignore all sidecars; rerun every stage.",
    )
    g.add_argument(
        "--force-stage",
        type=str,
        default=None,
        metavar="NAME[,NAME...]",
        help=(
            "Comma-separated list of stage names to force-rerun "
            "(e.g. analyze,correct). Stages downstream of a forced stage "
            "are also forced because their input sha256 will change."
        ),
    )
    g.add_argument(
        "--accept-prior-provenance",
        action="store_true",
        default=False,
        help=(
            "Treat rectify git_sha mismatch between prior and current run "
            "as non-blocking. Default behavior is to rerun any stage whose "
            "prior run was at a different git_sha (safer; conservative). "
            "Use this for cosmetic commits that don't affect output bytes "
            "(docs, comments)."
        ),
    )
    g.add_argument(
        "--dry-run-resume",
        action="store_true",
        default=False,
        help=(
            "Print the SKIP/RUN decision for each stage with full diff "
            "(which input or argv changed) and exit. No work done."
        ),
    )
