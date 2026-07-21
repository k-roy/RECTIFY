"""Shared ``--umi*`` argument group for RECTIFY commands.

Mirrors the ``add_organism_args(parser)`` precedent (``align_command.py``): a
single reusable arg-adder so the standard UMI vocabulary is identical across
``split`` / ``align`` / any future command, rather than each rolling its own.

Two of these flags (``--umi-edit-distance``, ``--umi-clustering``) already exist
verbatim on ``correct-cdna``; keeping the same names here makes the short-read
option a LIFT of an existing vocabulary, not a competing one.

Defaults encode the design decisions (see planning/299_corall_umi_design.md); all
three deviations from a naive reading of the request are deliberate and each is
exposed as a parameter so it can be settled empirically:

  * ``--umi-edit-distance`` default **1** (not 2): five independent tools default
    to 1; raising a 12-nt UMI to 2 grows the false-merge neighbourhood 17x
    (37 -> 631 of 4^12) to recover ~0.007% of UMIs at Q30.
  * ``--umi-position-tolerance`` default **0** (exact): no mainstream tool uses a
    tolerance window; soft-clip correction handles real end-wobble deterministically.
  * ``--umi-collapse`` default **dedup** (not consensus): pending the empirical
    family-size distribution; consensus is available for deep families.
"""
from __future__ import annotations

import argparse

from .extract import CORALL_UMI_LENGTH, DEFAULT_UMI_SEPARATOR, UMI_LOCATIONS

UMI_CLUSTERING_CHOICES = ("directional", "components")
UMI_COLLAPSE_CHOICES = ("dedup", "consensus")


def add_umi_args(parser: argparse.ArgumentParser) -> argparse._ArgumentGroup:
    """Add the standard UMI option group to ``parser``. Returns the group."""
    g = parser.add_argument_group(
        "UMI options (standard fixed-length N-mer, e.g. Lexogen CORALL V2)"
    )
    g.add_argument(
        "--umi", action="store_true", dest="umi",
        help="Enable UMI handling for short-read data. Off by default.",
    )
    g.add_argument(
        "--umi-location", dest="umi_location", default="read-id",
        choices=list(UMI_LOCATIONS),
        help="Where the UMI is: 'read-id' (already umi_tools-extracted, in the "
             "read name; the default and the path for pre-extracted data), "
             "'r1-start' (raw: first --umi-length bases of Read 1, the CORALL "
             "layout), or 'r2-start'. Default: read-id.",
    )
    g.add_argument(
        "--umi-length", dest="umi_length", type=int, default=None,
        help=f"UMI length in nt. Required for r1-start/r2-start; for CORALL V2 use "
             f"{CORALL_UMI_LENGTH}. For read-id it is an optional validation check.",
    )
    g.add_argument(
        "--umi-separator", dest="umi_separator", default=DEFAULT_UMI_SEPARATOR,
        help="Separator between read name and UMI for --umi-location=read-id "
             "(umi_tools uses '_'). Default: '_'.",
    )
    g.add_argument(
        "--umi-edit-distance", dest="umi_edit_distance", type=int, default=1,
        help="Max edit distance for merging UMIs within a fragment bucket. "
             "Default 1 (the field standard; see design doc for why 2 is worse "
             "for a 12-nt UMI).",
    )
    g.add_argument(
        "--umi-clustering", dest="umi_clustering", default="directional",
        choices=list(UMI_CLUSTERING_CHOICES),
        help="UMI graph clustering: 'directional' (umi_tools 2x count-gradient "
             "rule; default) or 'components' (single-linkage; comparison only).",
    )
    g.add_argument(
        "--umi-position-tolerance", dest="umi_position_tolerance", type=int, default=0,
        help="Bp tolerance when matching fragment 5' positions. Default 0 (exact, "
             "soft-clip corrected). Raising this AND --umi-edit-distance together "
             "compounds false merges; change only with empirical justification.",
    )
    g.add_argument(
        "--umi-collapse", dest="umi_collapse", default="dedup",
        choices=list(UMI_COLLAPSE_CHOICES),
        help="How to collapse a UMI family: 'dedup' (keep best-quality "
             "representative; default) or 'consensus' (quality-aware per-base "
             "consensus; worthwhile only for deep families).",
    )
    return g


def validate_umi_args(args: argparse.Namespace) -> None:
    """Fail fast on incoherent UMI arg combinations, with actionable messages."""
    if not getattr(args, "umi", False):
        return
    location = getattr(args, "umi_location", "read-id")
    length = getattr(args, "umi_length", None)
    if location in ("r1-start", "r2-start") and not length:
        raise ValueError(
            f"--umi-length is required with --umi-location={location} "
            f"(the number of bases to slice; CORALL V2 = {CORALL_UMI_LENGTH})."
        )
    if getattr(args, "umi_edit_distance", 1) < 0:
        raise ValueError("--umi-edit-distance must be >= 0.")
    if getattr(args, "umi_position_tolerance", 0) < 0:
        raise ValueError("--umi-position-tolerance must be >= 0.")
