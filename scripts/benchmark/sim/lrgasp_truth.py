#!/usr/bin/env python3
"""LRGASP NanoSim sim-truth join — the NanoSim arm of the simulator three-way.

The benchmark's external-validity plan triangulates three read sources with the
SAME scorer: our **pbsim3** vs LRGASP's **NanoSim** vs **real SIRV** reads. The
SIRV arm (absolute truth) and the pbsim3 arm (per-read MAF truth) already exist.
This module builds the NanoSim arm.

LRGASP ground truth is ANNOTATION-LEVEL, not per-read-alignment: a 2-column
``read_to_isoform.tsv`` (``read_id <TAB> transcript_id``, no header), where the
transcript_id is a versioned Ensembl id (human ``ENST`` + ~20% interleaved mouse
``ENSMUST`` artificial-novel decoys). So unlike pbsim3 there is NO per-read edit
script — we know only which transcript a read came from. The truth we CAN build
is the transcript-of-origin's annotated junctions + 3' terminus.

CRITICAL (the session-3 fragment bug, re-encountered): the scorer counts EVERY
unmatched truth junction as FN with no span intersection (``scorer.py`` ~L428).
NanoSim reads are fragments, so attaching the FULL transcript junction set to a
fragment read deflates recall with false FN. Since ``read_to_isoform`` carries no
per-read coordinates, the truth junctions must be restricted to the read's
**spanned + anchored** junctions using the read's ALIGNED genomic span at score
time — mirroring ``pbsim3_wrapper.project_maf_record`` (same ``MIN_JUNCTION_ANCHOR``
gate) so the two simulator arms are scored identically. The one caveat vs pbsim3:
pbsim3 anchors on the TRUE (MAF) span; here we only have the read's aligned span,
so a badly-clipping aligner could shrink its own truth (mild circularity). Recall
within the confidently-aligned core and spurious-junction FP are robust to it;
this is the documented annotation-level limitation, consistent with the SPEC's
"the three-way is DISTRIBUTIONAL, not locus-matched".

Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from typing import Dict, Iterator, List, Optional, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
from rectify.core.benchmark.truth_schema import (  # noqa: E402
    ReadTruth, JunctionTruth, SplitTag, write_truth_table,
)
from scripts.benchmark.sim.transcript_model import TranscriptModel  # noqa: E402
from scripts.benchmark.sim.gff_panel import build_panel_from_gtf  # noqa: E402

# Mirror pbsim3_wrapper.MIN_JUNCTION_ANCHOR (the 10bp anchored-junction gate in
# junction_scoring) so the NanoSim and pbsim3 arms apply IDENTICAL recall rules.
MIN_JUNCTION_ANCHOR = 10


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


# ---------------------------------------------------------------------------
# read_to_isoform parsing
# ---------------------------------------------------------------------------
def parse_read_to_isoform(path: str,
                          keep_prefixes: Optional[Tuple[str, ...]] = None
                          ) -> Tuple[Dict[str, str], Dict[str, int]]:
    """Parse the 2-column ``read_id <TAB> transcript_id`` truth join.

    No header. ``keep_prefixes`` filters by transcript_id prefix (e.g.
    ``("ENST",)`` to drop the mouse ``ENSMUST`` decoys, or ``None`` to keep all).
    Returns ``(read_to_tid, counts)`` where ``counts`` reports total / kept /
    by-prefix so the mouse-decoy fraction (~20%) is visible — a silent drop of
    those reads is the documented "human-only GTF fails ~20% of reads" trap.
    """
    read_to_tid: Dict[str, str] = {}
    counts = {"total": 0, "kept": 0, "ENST": 0, "ENSMUST": 0, "other": 0}
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rid, tid = parts[0], parts[1]
            counts["total"] += 1
            if tid.startswith("ENST"):
                counts["ENST"] += 1
            elif tid.startswith("ENSMUST"):
                counts["ENSMUST"] += 1
            else:
                counts["other"] += 1
            if keep_prefixes and not tid.startswith(keep_prefixes):
                continue
            read_to_tid[rid] = tid
            counts["kept"] += 1
    return read_to_tid, counts


# ---------------------------------------------------------------------------
# transcript truth catalogue
# ---------------------------------------------------------------------------
def build_truth_catalogue(gtf_path: str, genome: Dict[str, str], seed: int = 7):
    """Build the per-transcript truth lookup from an exon-feature GTF.

    Returns ``(models_by_tid, pairs, donors, acceptors)``. ``spliced_only=False``
    because sim reads originate from ALL transcripts (mono- AND multi-exon), not
    just the junction-bearing ones. Keyed by VERBATIM ``transcript_id`` (version
    included) so the ``read_to_isoform`` join matches exactly.
    """
    models, pairs, donors, acceptors = build_panel_from_gtf(
        gtf_path, genome, spliced_only=False, seed=seed)
    models_by_tid = {m.name: m for m in models}
    return models_by_tid, pairs, donors, acceptors


# ---------------------------------------------------------------------------
# spanned + anchored truth (the recall-correctness core)
# ---------------------------------------------------------------------------
def spanned_anchored_junctions(model: TranscriptModel, cov_lo: int, cov_hi: int,
                               pairs, donors, acceptors,
                               anchor: int = MIN_JUNCTION_ANCHOR
                               ) -> List[JunctionTruth]:
    """Truth junctions a read covering genome ``[cov_lo, cov_hi]`` (inclusive hi)
    can actually anchor: ``cov_lo <= intron_start - anchor`` AND
    ``cov_hi >= intron_end + anchor - 1`` (>=anchor exon bp on BOTH flanks).
    Identical predicate to ``pbsim3_wrapper.project_maf_record``."""
    out = []
    for j in model.junction_truths(pairs, donors, acceptors):
        if cov_lo <= j.intron_start - anchor and cov_hi >= j.intron_end + anchor - 1:
            out.append(j)
    return out


def cpa_for_span(model: TranscriptModel, cov_lo: int, cov_hi: int) -> Optional[int]:
    """The true 3' cleavage coord IF the read reaches the transcript terminus,
    else None (mirrors pbsim3_wrapper). + strand: ``genome_end-1`` when
    ``cov_hi >= genome_end-1``; - strand: ``genome_start`` when
    ``cov_lo <= genome_start``."""
    if model.strand == "+":
        return (model.genome_end - 1) if cov_hi >= model.genome_end - 1 else None
    return model.genome_start if cov_lo <= model.genome_start else None


def read_truth_for_span(read_id: str, model: TranscriptModel,
                        cov_lo: int, cov_hi: int,
                        pairs, donors, acceptors,
                        split: SplitTag = SplitTag.TEST,
                        stratum: str = "LRGASP_NANOSIM") -> ReadTruth:
    """Build a ``ReadTruth`` for one NanoSim read given its aligned genomic span
    ``[cov_lo, cov_hi]`` (inclusive hi, e.g. ``read.reference_start`` /
    ``read.reference_end - 1``). Junctions are spanned+anchored; ``indels`` is
    empty (NanoSim gives no per-read edit script — annotation-level truth only);
    ``true_cigar`` left blank for the same reason."""
    js = spanned_anchored_junctions(model, cov_lo, cov_hi, pairs, donors, acceptors)
    return ReadTruth(
        read_id=read_id, true_locus=model.name, true_transcript=model.name,
        chrom=model.chrom, strand=model.strand,
        genome_start=cov_lo, genome_end=cov_hi + 1,
        true_cigar="", junctions=js, indels=[],
        true_cpa=cpa_for_span(model, cov_lo, cov_hi),
        stratum=stratum, split=split)


# ---------------------------------------------------------------------------
# BAM-driven driver (Sherlock side — aligned NanoSim reads -> truth table)
# ---------------------------------------------------------------------------
def build_truths_from_bam(bam_path: str, read_to_tid: Dict[str, str],
                          models_by_tid: Dict[str, TranscriptModel],
                          pairs, donors, acceptors,
                          stratum: str = "LRGASP_NANOSIM"
                          ) -> Tuple[List[ReadTruth], Dict[str, int]]:
    """For each primary-mapped aligned read, look up its transcript-of-origin and
    emit a span-restricted ``ReadTruth``. Returns ``(rows, stats)``. ``stats``
    surfaces reads with no transcript in the catalogue (the human-only-GTF trap)
    and reads whose contig differs from the transcript's (chimeric / mismapped)."""
    import pysam
    rows: List[ReadTruth] = []
    stats = {"aligned": 0, "no_truth_tid": 0, "tid_not_in_panel": 0,
             "chrom_mismatch": 0, "emitted": 0}
    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        stats["aligned"] += 1
        tid = read_to_tid.get(read.query_name)
        if tid is None:
            stats["no_truth_tid"] += 1
            continue
        model = models_by_tid.get(tid)
        if model is None:
            stats["tid_not_in_panel"] += 1
            continue
        if read.reference_name != model.chrom:
            stats["chrom_mismatch"] += 1
            continue
        cov_lo = read.reference_start
        cov_hi = read.reference_end - 1
        rows.append(read_truth_for_span(read.query_name, model, cov_lo, cov_hi,
                                        pairs, donors, acceptors, stratum=stratum))
        stats["emitted"] += 1
    return rows, stats


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="LRGASP NanoSim sim-truth: read_to_isoform + GTF -> truth table")
    ap.add_argument("--read-to-isoform", required=True)
    ap.add_argument("--gtf", required=True, help="exon-feature GTF (GENCODE / SIRV)")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--bam", help="aligned NanoSim BAM; if omitted, AUDIT only")
    ap.add_argument("--out-truth", help="write truth table TSV (requires --bam)")
    ap.add_argument("--enst-only", action="store_true",
                    help="drop mouse ENSMUST decoys (else keep all)")
    args = ap.parse_args()
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
    from rectify.core.benchmark.scorer import load_genome

    keep = ("ENST",) if args.enst_only else None
    read_to_tid, counts = parse_read_to_isoform(args.read_to_isoform, keep)
    print(f"[lrgasp] read_to_isoform: total={counts['total']} kept={counts['kept']} "
          f"ENST={counts['ENST']} ENSMUST(mouse decoy)={counts['ENSMUST']} "
          f"other={counts['other']}")

    genome = load_genome(args.genome)
    models_by_tid, pairs, donors, acceptors = build_truth_catalogue(args.gtf, genome)
    print(f"[lrgasp] catalogue: transcripts={len(models_by_tid)} "
          f"annotated_pairs={len(pairs)}")
    in_panel = sum(1 for t in set(read_to_tid.values()) if t in models_by_tid)
    print(f"[lrgasp] distinct truth transcripts present in catalogue: {in_panel} "
          f"/ {len(set(read_to_tid.values()))}")

    if not args.bam:
        print("[lrgasp] AUDIT only (no --bam); not building per-read truth.")
        return
    rows, stats = build_truths_from_bam(args.bam, read_to_tid, models_by_tid,
                                        pairs, donors, acceptors)
    print(f"[lrgasp] from BAM: {stats}")
    if args.out_truth:
        n = write_truth_table(rows, args.out_truth)
        print(f"[lrgasp] wrote {n} truth rows -> {args.out_truth}")


if __name__ == "__main__":
    main()
