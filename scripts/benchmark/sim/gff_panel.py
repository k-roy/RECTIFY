#!/usr/bin/env python3
"""Real-GFF -> TranscriptModel loader for the Tier-2 transcriptome benchmark.

The Tier-1 controlled generators (``controlled.py``) build SYNTHETIC mini-loci.
Tier-2 instead simulates reads from a REAL annotated transcriptome so junction
recall/FDR is measured on real intron biology (the SPEC's external-validity tier).
This module turns a GFF + genome FASTA into ``TranscriptModel``s whose junction
truth is DERIVED FROM THE ANNOTATION (so every annotated intron a full-length read
crosses is an ANNOTATED-class truth junction).

Yeast SGD GFF specifics (R64-5-1): there is NO plain ``exon`` feature — the file
carries ``mRNA`` spans + explicit ``intron`` features (whose ``Parent=`` is a
comma-separated list of mRNA IDs). Exon blocks are therefore reconstructed by
cutting each mRNA's introns out of its span. tRNA introns (Parent ending in a
non-mRNA id) are skipped. GFF coords are 1-based inclusive -> converted to
0-based half-open here (the convention the scorer + TranscriptModel use).

NIC/NNC NOTE (see SIMULATION_BENCHMARK_SPEC.md): full-length reads of annotated
transcripts yield ANNOTATED junctions only. Novel-junction RECALL needs isoform
injection (exon-skip -> NIC, novel-site -> NNC); spurious-junction FDR IS
measurable from annotated-only truth. So a Tier-2 run on this loader measures
ANNOTATED-recall + spurious-FDR, NOT novel-junction recall. Inject isoforms (a
later increment) to measure novel recall.

Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import gzip
import os
import random
import sys
from typing import Dict, List, Optional, Set, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
from scripts.benchmark.sim.transcript_model import TranscriptModel, Exon  # noqa: E402


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _attrs(field9: str) -> Dict[str, str]:
    out = {}
    for kv in field9.rstrip(";").split(";"):
        kv = kv.strip()
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def parse_gff(gff_path: str) -> Tuple[Dict[str, Tuple[str, str, int, int]],
                                      Dict[str, List[Tuple[int, int]]]]:
    """Parse mRNA + intron features.

    Returns ``(mrnas, introns_by_mrna)`` where ``mrnas[id] = (chrom, strand,
    start0, end0)`` (0-based half-open genome span) and ``introns_by_mrna[id] =
    [(istart0, iend0), ...]`` (0-based half-open intron spans linked via the
    intron feature's comma-separated ``Parent=``). Only introns whose parent is a
    known mRNA id are kept (tRNA introns dropped)."""
    mrnas: Dict[str, Tuple[str, str, int, int]] = {}
    raw_introns: List[Tuple[str, int, int, List[str]]] = []
    with _open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _src, ftype, start, end, _sc, strand, _frame, attr = f[:9]
            if ftype == "mRNA":
                a = _attrs(attr)
                mid = a.get("ID")
                if mid:
                    mrnas[mid] = (chrom, strand, int(start) - 1, int(end))
            elif ftype == "intron":
                a = _attrs(attr)
                parents = a.get("Parent", "").split(",") if a.get("Parent") else []
                raw_introns.append((chrom, int(start) - 1, int(end), parents))
    introns_by_mrna: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, s0, e0, parents in raw_introns:
        for p in parents:
            if p in mrnas and mrnas[p][0] == chrom:
                introns_by_mrna.setdefault(p, []).append((s0, e0))
    return mrnas, introns_by_mrna


def _exons_from_span(start0: int, end0: int,
                     introns0: List[Tuple[int, int]]) -> List[Exon]:
    """Exon blocks = the mRNA span minus its (sorted) introns."""
    exons: List[Exon] = []
    cur = start0
    for (is0, ie0) in sorted(introns0):
        if is0 < cur or ie0 > end0:
            continue  # intron outside the mRNA span (defensive)
        if is0 > cur:
            exons.append(Exon(cur, is0))
        cur = ie0
    if cur < end0:
        exons.append(Exon(cur, end0))
    return exons


def annotated_sets(mrnas, introns_by_mrna
                   ) -> Tuple[Set[Tuple[str, int, int]],
                              Set[Tuple[str, int]], Set[Tuple[str, int]]]:
    """The catalogue used to classify a read's junctions ANNOTATED/NIC/NNC: every
    annotated intron, donor, and acceptor across ALL mRNAs."""
    pairs: Set[Tuple[str, int, int]] = set()
    donors: Set[Tuple[str, int]] = set()
    acceptors: Set[Tuple[str, int]] = set()
    for mid, (chrom, _strand, _s, _e) in mrnas.items():
        for (is0, ie0) in introns_by_mrna.get(mid, []):
            pairs.add((chrom, is0, ie0))
            donors.add((chrom, is0))
            acceptors.add((chrom, ie0))
    return pairs, donors, acceptors


def build_panel(gff_path: str, genome: Dict[str, str],
                spliced_only: bool = True, max_transcripts: Optional[int] = None,
                seed: int = 7) -> Tuple[List[TranscriptModel],
                                        Set[Tuple[str, int, int]],
                                        Set[Tuple[str, int]],
                                        Set[Tuple[str, int]]]:
    """Build ``TranscriptModel``s from the GFF + an in-memory genome dict.

    ``spliced_only`` keeps only multi-exon transcripts (the junction-recall set);
    set False to include intronless mRNAs (placement baseline). ``max_transcripts``
    randomly subsamples (seeded) for a quick run. Returns the models + the
    annotation sets for junction classification.
    """
    mrnas, introns_by_mrna = parse_gff(gff_path)
    pairs, donors, acceptors = annotated_sets(mrnas, introns_by_mrna)
    models: List[TranscriptModel] = []
    skipped_no_contig = 0
    for mid, (chrom, strand, s0, e0) in mrnas.items():
        introns0 = introns_by_mrna.get(mid, [])
        if spliced_only and not introns0:
            continue
        if chrom not in genome:
            skipped_no_contig += 1
            continue
        exons = _exons_from_span(s0, e0, introns0)
        if not exons:
            continue
        models.append(TranscriptModel(name=mid, chrom=chrom, strand=strand,
                                      exons=exons, genome_seq=genome[chrom]))
    if skipped_no_contig:
        print(f"[gff_panel] WARNING: skipped {skipped_no_contig} mRNAs whose contig "
              f"is absent from the genome FASTA", file=sys.stderr)
    models.sort(key=lambda m: m.name)
    if max_transcripts is not None and len(models) > max_transcripts:
        random.Random(seed).shuffle(models)
        models = sorted(models[:max_transcripts], key=lambda m: m.name)
    return models, pairs, donors, acceptors


def main():
    ap = argparse.ArgumentParser(description="GFF->TranscriptModel panel loader (audit/dry-run)")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--max-transcripts", type=int, default=None)
    ap.add_argument("--include-intronless", action="store_true")
    args = ap.parse_args()
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
    from rectify.core.benchmark.scorer import load_genome
    genome = load_genome(args.genome)
    models, pairs, donors, acceptors = build_panel(
        args.gff, genome, spliced_only=not args.include_intronless,
        max_transcripts=args.max_transcripts)
    n_introns = sum(len(m.introns()) for m in models)
    print(f"models={len(models)} total_introns_in_models={n_introns} "
          f"annotated_pairs={len(pairs)} donors={len(donors)} acceptors={len(acceptors)}")
    # sanity: every model intron classifies ANNOTATED against the catalogue
    n_anno = n_other = 0
    for m in models[:200]:
        for j in m.junction_truths(pairs, donors, acceptors):
            if j.klass.value == "ANNOTATED":
                n_anno += 1
            else:
                n_other += 1
    print(f"[sample of <=200 models] junctions ANNOTATED={n_anno} non-ANNOTATED={n_other} "
          f"(expect non-ANNOTATED=0 for annotation-derived models)")


if __name__ == "__main__":
    main()
