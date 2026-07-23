"""`rectify cma` migration command — build / validate / verify / expand.

Exercises the streaming build path (coord-sorted inputs → name-sort →
_iter_name_grouped_bams → build_cma) and the verify gate on the bundled
5-aligner DRS fixture.
"""

import argparse
from pathlib import Path

import pysam

import rectify
from rectify.data import get_bundled_genome_path
from rectify.core.commands import cma_command

FIXTURE = Path(rectify.__file__).parent / "data" / "validation" / "aligners"
ALIGNERS = ["minimap2", "mapPacBio", "gapmm2", "deSALT", "uLTRA"]


def _bams():
    return [str(FIXTURE / f"validation_reads.{a}.bam") for a in ALIGNERS]


def _genome():
    return str(get_bundled_genome_path("saccharomyces_cerevisiae"))


def _args(**kw):
    return argparse.Namespace(**kw)


def test_cli_parser_wires_cma():
    from rectify.cli import create_parser

    args = create_parser().parse_args(["cma", "validate", "--cma", "x.bam"])
    assert args.command == "cma" and args.cma_subcommand == "validate"


def test_build_validate_verify_expand(tmp_path):
    bams = _bams()
    genome = _genome()
    cma = str(tmp_path / "drs.cma.bam")

    # build (bare paths → aligner inferred from *.<aligner>.bam)
    assert cma_command._run_build(
        _args(aligner_bams=bams, out=cma, panel=",".join(ALIGNERS), genome=genome)
    ) == 0
    assert Path(cma).exists()

    # validate
    assert cma_command._run_validate(_args(cma=cma)) == 0

    # verify — BOTH deletion gates must pass on the fixture
    assert cma_command._run_verify(
        _args(cma=cma, aligner_bams=bams, max_reads=0, genome=genome)
    ) == 0

    # expand all aligners → 36 reads × 5 aligners
    out_all = str(tmp_path / "expand_all.bam")
    assert cma_command._run_expand(_args(cma=cma, out=out_all, aligner=None, genome=genome)) == 0
    with pysam.AlignmentFile(out_all, "rb") as f:
        assert sum(1 for _ in f) == 180

    # expand a single aligner → 36 reads
    out_mm = str(tmp_path / "expand_mm.bam")
    assert cma_command._run_expand(_args(cma=cma, out=out_mm, aligner="minimap2", genome=genome)) == 0
    with pysam.AlignmentFile(out_mm, "rb") as f:
        assert sum(1 for _ in f) == 36


def test_aligner_name_inference():
    parsed = cma_command._parse_aligner_bams(
        ["/x/s.minimap2.bam", "uLTRA=/y/z.bam"]
    )
    assert parsed == {"minimap2": "/x/s.minimap2.bam", "uLTRA": "/y/z.bam"}


def test_build_cma_from_bams_helper(tmp_path):
    """The reusable builder used by both `cma build` and `align --emit-cma`."""
    from rectify.core.multialign import build_cma_from_bams, validate_cma

    aligner_bams = {a: str(FIXTURE / f"validation_reads.{a}.bam") for a in ALIGNERS}
    cma = str(tmp_path / "h.cma.bam")
    stats = build_cma_from_bams(aligner_bams, cma, panel=ALIGNERS, genome=None)
    assert stats["reads"] == 36
    assert validate_cma(cma) == []


def test_align_parser_has_emit_cma():
    import argparse

    from rectify.cli import create_parser

    p = create_parser()
    sub = next(a for a in p._actions if isinstance(a, argparse._SubParsersAction))
    align = sub.choices["align"]
    opts = {s for act in align._actions for s in act.option_strings}
    assert "--emit-cma" in opts
