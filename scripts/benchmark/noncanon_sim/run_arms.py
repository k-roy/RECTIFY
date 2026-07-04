#!/usr/bin/env python3
"""Component C — ARM-A/B/C ALIGN+REFINE runner for the yeast non-canonical
long-read splice-simulator vetting harness.

See ``scripts/benchmark/noncanon_sim/SPEC.md`` for the coordination contract.

Pipeline (this component):

    reads.fastq + sim_ref.fa
        --> minimap2 -ax splice  (spliced ONT alignment)
        --> samtools sort + index         (aligned.sorted.bam)
        --> build_junction_pool( [aligned.sorted.bam], annotated_from_GFF )
        --> refine_bam_junctions THREE ways (n_workers=1), each sorted+indexed:
              arm_A.bam : motif_blind=False                      (incumbent)
              arm_B.bam : motif_blind=True                       (motif-blind)
              arm_C.bam : motif_blind=True + penalty_table_path  (-logP law)

Only ``arm_{A,B,C}.bam`` are the contract outputs (SPEC item 6); the sorted
alignment BAM and an ``arms_stats.json`` sidecar are written for the
integrator / component-D scorer.

INPUTS consumed
    reads.fastq   (SPEC item 4) — read id ``<tid>_r<NNN>``
    sim_ref.fa    (SPEC item 2) — the genomic reference the reads align to;
                  truth junction coords are in this frame.

OUTPUTS produced (in --outdir)
    aligned.sorted.bam{,.bai}   — minimap2 alignment (intermediate)
    arm_A.bam{,.bai}            — refined, motif_blind=False
    arm_B.bam{,.bai}            — refined, motif_blind=True
    arm_C.bam{,.bai}            — refined, motif_blind=True + penalty table
    arms_stats.json             — per-arm refine stats + exact commands (aux)

The refined BAMs are sorted+indexed and contain ALL reads (N-op reads with
refined junction boundaries; non-N-op reads passed through unchanged), in the
sim_ref.fa coordinate frame, ready for the ambiguity-aware junction scorer.

Author: Kevin R. Roy (harness component C)
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

# ---------------------------------------------------------------------------
# Resolve the rectify package (repo root is three levels up from this file:
# scripts/benchmark/noncanon_sim/run_arms.py -> repo root).
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import pysam  # noqa: E402

from rectify.core.splice.junction_refiner import (  # noqa: E402
    build_junction_pool,
    refine_bam_junctions,
)
from rectify.core.consensus.consensus import load_annotated_junctions  # noqa: E402
from rectify.utils.genome import (  # noqa: E402
    register_genome_contigs,
    standardize_chrom_name,
)

# ---------------------------------------------------------------------------
# Bundled defaults (yeast). SPEC: genome/annotation/penalty table are bundled.
# ---------------------------------------------------------------------------
_DATA = _REPO_ROOT / "rectify" / "data" / "genomes" / "saccharomyces_cerevisiae"
DEFAULT_GFF = _DATA / "saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"
# Arm-C penalty table: SPEC defines C as "B + penalty_table_path=<yeast
# penalty_scores.tsv>". Use ONLY the DRS penalty_scores.tsv — deliberately NOT
# the STR table (production passes str_penalty_table too, but SPEC does not, and
# adding it would make C not-quite-C).
DEFAULT_PENALTY_TABLE = _DATA / "penalty_tables" / "penalty_scores.tsv"


# ---------------------------------------------------------------------------
# minimap2 alignment
# ---------------------------------------------------------------------------
def _require_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(f"'{name}' not found on PATH (required by run_arms.py)")
    return path


def align_reads(
    reads_fastq: Path,
    sim_ref: Path,
    out_bam: Path,
    *,
    u_mode: str = "f",
    max_intron: int = 5000,
    kmer: int = 14,
    threads: int = 1,
    extra: Optional[List[str]] = None,
) -> List[str]:
    """Align ONT-like long reads with a spliced aligner (minimap2 -ax splice),
    then coordinate-sort + index.  Returns the minimap2 argv (for provenance).

    Alignment is deliberately ANNOTATION-BLIND (no --junc-bed): the experiment
    measures the *refiner's* motif handling, so biasing the initial placement
    toward annotated/canonical sites would confound the snap-vs-hold signal.

    -uf (u_mode='f') = discover canonical GT-AG on the transcript-forward
    strand.  Correct for sense-oriented direct-RNA / cDNA reads and handles
    BOTH genomic strands (a minus-strand gene's sense read maps reverse).  Flip
    to u_mode='b' only if the read generator emits mixed-orientation reads.
    """
    mm2 = _require_tool("minimap2")
    samtools = _require_tool("samtools")

    mm2_cmd = [
        mm2,
        "-ax", "splice",
        f"-u{u_mode}",
        f"-k{kmer}",
        "-G", str(max_intron),
        "--splice-flank=no",
        "--secondary=no",
        "--MD",
        "-t", str(threads),
    ]
    if extra:
        mm2_cmd.extend(extra)
    mm2_cmd.extend([str(sim_ref), str(reads_fastq)])

    out_bam.parent.mkdir(parents=True, exist_ok=True)
    # minimap2 (SAM on stdout) | samtools sort -> out_bam
    proc_mm2 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE)
    sort_cmd = [samtools, "sort", "-@", str(threads), "-o", str(out_bam), "-"]
    proc_sort = subprocess.run(sort_cmd, stdin=proc_mm2.stdout)
    proc_mm2.stdout.close()
    mm2_rc = proc_mm2.wait()
    if mm2_rc != 0:
        raise RuntimeError(f"minimap2 failed (rc={mm2_rc}): {' '.join(mm2_cmd)}")
    if proc_sort.returncode != 0:
        raise RuntimeError(f"samtools sort failed (rc={proc_sort.returncode})")
    subprocess.run([samtools, "index", str(out_bam)], check=True)
    return mm2_cmd


# ---------------------------------------------------------------------------
# Genome dict (sim_ref frame) + chrom-resolution guard
# ---------------------------------------------------------------------------
def load_sim_genome(sim_ref: Path) -> Dict[str, str]:
    """Load the sim_ref.fa into a {chrom: seq} dict in EXACTLY the frame the
    reads were aligned to, and make ``standardize_chrom_name`` resolve those
    contigs to a present key.

    Silent-failure guard (advisor-flagged): the refiner does
    ``genome.get(standardize_chrom_name(read.reference_name), '')`` — a miss
    yields ``''`` and the read is passed through UNCHANGED, so all three arms
    would come out byte-identical and the benchmark would be garbage while
    "succeeding".  We (a) register the contigs so standardize returns them
    verbatim, and (b) additionally alias each key under its standardized name
    (covers the NCBI-header case where standardize maps e.g. NC_001133 -> chrI
    before consulting the registry).  The per-read assertion in ``main`` is the
    final backstop.
    """
    fa = pysam.FastaFile(str(sim_ref))
    genome: Dict[str, str] = {name: fa.fetch(name) for name in fa.references}
    register_genome_contigs(genome.keys())
    # Alias under standardized names so genome.get(standardize(ref)) always hits.
    for name in list(genome.keys()):
        std = standardize_chrom_name(name)
        if std != name and std not in genome:
            genome[std] = genome[name]
    return genome


def _assert_chroms_resolve(bam_path: Path, genome: Dict[str, str], n_check: int = 50) -> None:
    """Confirm aligned reads' chroms resolve to a non-empty genome sequence.

    Raises RuntimeError (not a silent pass-through) if the sim_ref frame and the
    refiner's chrom standardization disagree.
    """
    seen = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.reference_name is None:
                continue
            std = standardize_chrom_name(read.reference_name)
            seq = genome.get(std, "")
            if not seq:
                raise RuntimeError(
                    "CHROM RESOLUTION FAILURE: aligned read on contig "
                    f"'{read.reference_name}' -> standardize -> '{std}' is not a "
                    "non-empty key in the sim_ref genome dict. The refiner would "
                    "silently pass every read through unchanged (all arms "
                    "identical). Genome keys: "
                    f"{sorted(genome.keys())[:10]}{'...' if len(genome) > 10 else ''}"
                )
            seen += 1
            if seen >= n_check:
                break
    if seen == 0:
        raise RuntimeError(
            f"No mapped reads found in {bam_path}; minimap2 produced no "
            "alignments. Check reads.fastq / sim_ref.fa orientation (-u mode)."
        )


# ---------------------------------------------------------------------------
# Arm runner
# ---------------------------------------------------------------------------
def run_arms(
    aligned_bam: Path,
    genome: Dict[str, str],
    annotated_junctions,
    outdir: Path,
    *,
    penalty_table_path: Path,
    threads: int = 1,
) -> Dict[str, dict]:
    """Build the junction pool once, then run the three refiner arms.

    The SAME prebuilt pool + annotated set feed all three arms, so the only
    difference across arms is ``motif_blind`` (A vs B) and the penalty table
    (C) — never the candidate pool.
    """
    aligner_bams = [str(aligned_bam)]

    # Build the pool once (advisor: identical pool across arms).
    all_junctions, annot_set = build_junction_pool(aligner_bams, annotated_junctions)

    arm_specs = [
        # (name, motif_blind, penalty_table_path)
        ("A", False, None),
        ("B", True, None),
        ("C", True, str(penalty_table_path)),
    ]
    results: Dict[str, dict] = {}
    for name, motif_blind, pen_path in arm_specs:
        out_bam = outdir / f"arm_{name}.bam"
        stats = refine_bam_junctions(
            input_bam=str(aligned_bam),
            output_bam=str(out_bam),
            aligner_bams=aligner_bams,
            annotated_junctions=annotated_junctions,
            genome=genome,
            penalty_table_path=pen_path,
            prebuilt_junction_pool=all_junctions,
            prebuilt_annotated_set=annot_set,
            sort_and_index=True,
            sort_threads=threads,
            n_workers=1,
            motif_blind=motif_blind,
        )
        results[name] = {
            "output_bam": str(out_bam),
            "motif_blind": motif_blind,
            "penalty_table_path": pen_path,
            "stats": stats,
        }
    return {
        "pool_junctions": len(all_junctions),
        "pool_annotated": len(annot_set),
        "arms": results,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Component C: align (minimap2 -ax splice) + refine 3 arms "
                    "(A=motif_blind off, B=motif_blind on, C=B+penalty table).",
    )
    p.add_argument("--reads", required=True, type=Path,
                   help="reads.fastq (SPEC item 4)")
    p.add_argument("--sim-ref", required=True, type=Path,
                   help="sim_ref.fa the reads align to (SPEC item 2)")
    p.add_argument("--outdir", required=True, type=Path,
                   help="output directory for arm_{A,B,C}.bam + aligned.sorted.bam")
    p.add_argument("--gff", type=Path, default=DEFAULT_GFF,
                   help=f"annotation GFF for annotated junctions (default bundled yeast: {DEFAULT_GFF})")
    p.add_argument("--penalty-table", type=Path, default=DEFAULT_PENALTY_TABLE,
                   help=f"arm-C penalty_scores.tsv (default bundled yeast DRS: {DEFAULT_PENALTY_TABLE})")
    p.add_argument("--u-mode", default="f", choices=["f", "b", "n"],
                   help="minimap2 -u splice-strand mode (default 'f' for "
                        "sense-oriented direct-RNA/cDNA reads; use 'b' if the "
                        "read generator emits mixed orientation)")
    p.add_argument("--max-intron", type=int, default=5000,
                   help="minimap2 -G max intron (default 5000 = yeast)")
    p.add_argument("--kmer", type=int, default=14,
                   help="minimap2 -k (default 14 for noisy ONT)")
    p.add_argument("--threads", type=int, default=1,
                   help="threads for minimap2/samtools (refiner is n_workers=1)")
    p.add_argument("--mm2-extra", default=None,
                   help="extra minimap2 flags (space-separated string), appended verbatim")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    for pth, label in [(args.reads, "reads"), (args.sim_ref, "sim_ref"), (args.gff, "gff")]:
        if not Path(pth).exists():
            raise FileNotFoundError(f"--{label} not found: {pth}")
    if not Path(args.penalty_table).exists():
        raise FileNotFoundError(f"--penalty-table not found: {args.penalty_table}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    aligned_bam = outdir / "aligned.sorted.bam"

    mm2_extra = args.mm2_extra.split() if args.mm2_extra else None

    # 1. Align.
    mm2_cmd = align_reads(
        args.reads, args.sim_ref, aligned_bam,
        u_mode=args.u_mode, max_intron=args.max_intron, kmer=args.kmer,
        threads=args.threads, extra=mm2_extra,
    )

    # 2. Genome dict in sim_ref frame + chrom-resolution guard.
    genome = load_sim_genome(args.sim_ref)
    _assert_chroms_resolve(aligned_bam, genome)

    # 3. Annotated junctions from the (yeast) GFF.
    annotated_junctions = load_annotated_junctions(str(args.gff))

    # 4. Refine three arms.
    summary = run_arms(
        aligned_bam, genome, annotated_junctions, outdir,
        penalty_table_path=args.penalty_table, threads=args.threads,
    )

    manifest = {
        "reads": str(args.reads),
        "sim_ref": str(args.sim_ref),
        "gff": str(args.gff),
        "aligned_bam": str(aligned_bam),
        "minimap2_cmd": mm2_cmd,
        "u_mode": args.u_mode,
        "n_annotated_junctions_gff": len(annotated_junctions),
        **summary,
    }
    stats_path = outdir / "arms_stats.json"
    with open(stats_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    # Machine-readable summary to stdout (advisor: surface per-arm stats so the
    # integrator can tell "chrom didn't resolve" from "no N-ops placed").
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
