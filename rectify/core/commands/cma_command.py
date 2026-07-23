"""``rectify cma`` — compressed-multialign (CMA) migration & inspection.

Subcommands:
  build     per-aligner BAMs → one CMA BAM (SEQ/QUAL stored once).
  expand    CMA → per-aligner BAM(s) (materialize the reconstructed placements).
  validate  structural invariants of a CMA.
  verify    run BOTH deletion gates (planning/254 §3.3) comparing expand(CMA) to
            the original per-aligner BAMs — the check that licenses deleting them.

The migration ingest for existing runs (doc M2): point ``build`` at a run's
merged per-aligner BAMs, ``verify`` on a subsample, then (separately, with an
explicit file list) delete the originals. This command NEVER deletes anything.
"""

from __future__ import annotations

import argparse
import collections.abc
import logging
import os
import tempfile

import pysam

from ..multialign import build_cma, expand, load_aligner_records, validate_cma
from ..multialign.cma_schema import decode_eq_seq

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
class _LazyGenome(collections.abc.Mapping):
    """Fetch-on-demand ``{reference_name: seq}`` mapping over a FASTA.

    A read-only Mapping so it is a drop-in for the ``Dict[str, str]`` that
    ``extract_alignment_info`` / ``select_best_alignment`` / ``decode_eq_seq``
    expect (``genome[chrom]``, ``chrom in genome``, ``genome.get(chrom)``).
    Contigs are fetched+cached on first access, so explicit-SEQ (production)
    BAMs — where ``decode_eq_seq`` never runs — load nothing.
    """

    def __init__(self, fa_path):
        self._fa = pysam.FastaFile(fa_path)
        self._refs = set(self._fa.references)
        self._cache = {}

    def __getitem__(self, name):
        if name not in self._refs:
            raise KeyError(name)
        if name not in self._cache:
            self._cache[name] = self._fa.fetch(name)
        return self._cache[name]

    def __iter__(self):
        return iter(self._fa.references)

    def __len__(self):
        return len(self._refs)

    def close(self):
        self._fa.close()


def _parse_aligner_bams(items):
    """``aligner=path`` or bare ``path`` (aligner inferred from ``*.<aligner>.bam``)."""
    out = {}
    for it in items:
        if "=" in it:
            name, path = it.split("=", 1)
        else:
            path = it
            base = os.path.basename(path)
            stem = base[:-4] if base.endswith(".bam") else base
            name = stem.split(".")[-1]
        if name in out:
            raise ValueError(f"duplicate aligner name '{name}' in --aligner-bams")
        out[name] = path
    return out


def _ensure_name_sorted(path, tmpdir):
    """Return a name(queryname)-sorted copy of ``path`` (the K-way merge needs it)."""
    with pysam.AlignmentFile(path, "rb") as f:
        so = f.header.to_dict().get("HD", {}).get("SO")
    if so == "queryname":
        return path
    out = os.path.join(tmpdir, os.path.basename(path) + ".nsort.bam")
    pysam.sort("-n", "-o", out, path)
    return out


def _genome_from_args(args):
    g = getattr(args, "genome", None)
    return _LazyGenome(g) if g else None


# --------------------------------------------------------------------------- #
# build
# --------------------------------------------------------------------------- #
def _run_build(args):
    aligner_bams = _parse_aligner_bams(args.aligner_bams)
    panel = args.panel.split(",") if args.panel else list(aligner_bams)
    genome = _genome_from_args(args)
    from ..consensus.consensus import _iter_name_grouped_bams

    with tempfile.TemporaryDirectory(prefix="cma_build_") as td:
        sorted_paths = {a: _ensure_name_sorted(p, td) for a, p in aligner_bams.items()}
        first = next(iter(sorted_paths.values()))
        with pysam.AlignmentFile(first, "rb") as f:
            header = f.header
        logger.info("Building CMA from %d aligner BAMs → %s", len(sorted_paths), args.out)
        stats = build_cma(
            _iter_name_grouped_bams(sorted_paths), header, args.out, panel, genome=genome
        )
    if genome:
        genome.close()
    problems = validate_cma(args.out)
    print(f"[cma build] {stats['reads']} reads, {stats['records']} records → {args.out}")
    if problems:
        print(f"[cma build] VALIDATION FAILED ({len(problems)} problems):")
        for p in problems[:20]:
            print(f"   - {p}")
        return 1
    print("[cma build] validation OK")
    return 0


# --------------------------------------------------------------------------- #
# expand
# --------------------------------------------------------------------------- #
def _run_expand(args):
    genome = _genome_from_args(args)
    with pysam.AlignmentFile(args.cma, "rb") as cma:
        header = cma.header
    want = args.aligner
    n = 0
    with pysam.AlignmentFile(args.out, "wb", header=header) as out:
        for _key, adict in expand(args.cma, genome):
            for aligner, rec in adict.items():
                if want and aligner != want:
                    continue
                out.write(rec)
                n += 1
    if genome:
        genome.close()
    print(f"[cma expand] wrote {n} records → {args.out}"
          + (f" (aligner={want})" if want else ""))
    return 0


# --------------------------------------------------------------------------- #
# validate
# --------------------------------------------------------------------------- #
def _run_validate(args):
    problems = validate_cma(args.cma)
    if problems:
        print(f"[cma validate] {len(problems)} problems:")
        for p in problems[:50]:
            print(f"   - {p}")
        return 1
    print("[cma validate] OK — invariants hold")
    return 0


# --------------------------------------------------------------------------- #
# verify — both deletion gates
# --------------------------------------------------------------------------- #
def _placement_view(rec, genome):
    tags = dict(rec.get_tags())
    return (
        bool(rec.is_reverse),
        rec.reference_name,
        rec.reference_start,
        rec.mapping_quality,
        rec.cigarstring,
        decode_eq_seq(rec, genome),
        tags.get("NM"),
        tags.get("MD"),
    )


def _explicit_clone(rec, header, genome):
    c = pysam.AlignedSegment.fromstring(rec.to_string(), header)
    d = decode_eq_seq(rec, genome)
    if d is not None and d != c.query_sequence:
        q = c.query_qualities
        c.query_sequence = d
        c.query_qualities = q
    return c


def _winner(reads, genome):
    from ..consensus.extract import extract_alignment_info
    from ..consensus.select import select_best_alignment

    al = {a: extract_alignment_info(r, a, genome) for a, r in reads.items()}
    res = select_best_alignment(al, genome, None, tiebreak="rectify")
    c3 = getattr(al.get(res.best_aligner), "corrected_3prime", None)
    return res.best_aligner, c3


def _run_verify(args):
    aligner_bams = _parse_aligner_bams(args.aligner_bams)
    genome = _genome_from_args(args)
    with pysam.AlignmentFile(args.cma, "rb") as cma:
        header = cma.header

    # Originals (bounded to a validation sample if --max-reads).
    orig = {}
    for k, adict in load_aligner_records(aligner_bams):
        orig[k] = adict
        if args.max_reads and len(orig) >= args.max_reads:
            break
    expanded = {k: adict for k, adict in expand(args.cma, genome)}

    sec_fail = prim_fail = checked = 0
    for k in orig:
        if k not in expanded:
            sec_fail += 1
            continue
        checked += 1
        o, e = orig[k], expanded[k]
        if set(o) != set(e):
            sec_fail += 1
            continue
        # Secondary gate: load-bearing field view per aligner.
        if any(_placement_view(e[a], genome) != _placement_view(o[a], genome) for a in o):
            sec_fail += 1
        # Primary gate: selection winner + corrected 3' end (vs =-decoded originals).
        o_expl = {a: _explicit_clone(o[a], header, genome) for a in o}
        if _winner(e, genome) != _winner(o_expl, genome):
            prim_fail += 1
    if genome:
        genome.close()

    print(f"[cma verify] checked {checked} reads (of {len(orig)} sampled)")
    print(f"[cma verify] SECONDARY (field view) failures: {sec_fail}")
    print(f"[cma verify] PRIMARY   (selection winner+3') failures: {prim_fail}")
    ok = sec_fail == 0 and prim_fail == 0
    print(f"[cma verify] {'PASS — safe to delete the per-aligner BAMs' if ok else 'FAIL — do NOT delete'}")
    return 0 if ok else 1


# --------------------------------------------------------------------------- #
# parser + dispatch
# --------------------------------------------------------------------------- #
def create_cma_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    from ...data import add_organism_args

    p = subparsers.add_parser("cma", help="Compressed-multialign build/expand/validate/verify")
    sub = p.add_subparsers(dest="cma_subcommand")

    b = sub.add_parser("build", help="Build a CMA from per-aligner BAMs")
    b.add_argument("--aligner-bams", nargs="+", required=True,
                   help="aligner=path ... (or bare path; aligner inferred from *.<aligner>.bam)")
    b.add_argument("--out", required=True, help="output CMA BAM (name/RN-sorted)")
    b.add_argument("--panel", default=None, help="comma-list panel/aligner order")
    add_organism_args(b)

    e = sub.add_parser("expand", help="Materialize per-aligner records from a CMA")
    e.add_argument("--cma", required=True)
    e.add_argument("--out", required=True)
    e.add_argument("--aligner", default=None, help="only this aligner (default: all)")
    add_organism_args(e)

    v = sub.add_parser("validate", help="Structural validation of a CMA")
    v.add_argument("--cma", required=True)

    vf = sub.add_parser("verify", help="Run both deletion gates vs the original per-aligner BAMs")
    vf.add_argument("--cma", required=True)
    vf.add_argument("--aligner-bams", nargs="+", required=True)
    vf.add_argument("--max-reads", type=int, default=0,
                    help="verify only the first N reads (validation sample; 0=all)")
    add_organism_args(vf)
    return p


def run(args):
    sc = getattr(args, "cma_subcommand", None)
    if sc == "build":
        return _run_build(args)
    if sc == "expand":
        return _run_expand(args)
    if sc == "validate":
        return _run_validate(args)
    if sc == "verify":
        return _run_verify(args)
    print("usage: rectify cma {build,expand,validate,verify} ...")
    return 1
