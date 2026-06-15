#!/usr/bin/env python
"""Shared helpers for Round-2 cDNA Phase-0 lift/score: load CdnaModels from blockmaps.json,
genome->transcript projection (inverse of lift, for the F.0 round-trip), and read->cdna
candidate construction.
"""
from __future__ import annotations
import json
from pathlib import Path
from typing import Dict, List, Tuple
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from liftover import Block, CdnaModel, revcomp, M, I, D, N, S, H, EQ, X  # noqa: E402

_REF_CONSUMING = frozenset({M, D, N, EQ, X})
_QUERY_CONSUMING = frozenset({M, I, S, EQ, X})


def load_models(blockmaps_json) -> Dict[str, CdnaModel]:
    with open(blockmaps_json) as fh:
        data = json.load(fh)
    models = {}
    for cid, d in data.items():
        blocks = [Block(t_start=b["t_start"], g_start=b["g_start"], g_end=b["g_end"],
                        is_pad=b["is_pad"], kind=b["kind"]) for b in d["blocks"]]
        models[cid] = CdnaModel(cdna_id=cid, chrom=d["chrom"], strand=d["strand"],
                                blocks=blocks, length=d["length"])
    return models


def g2t(model: CdnaModel, g: int):
    """Inverse of CdnaModel.t2g: plus-strand genome coord -> transcript coord, or None
    if g falls outside the model's blocks (i.e. inside an intron / off-cDNA)."""
    for b in model.blocks:
        if b.g_start <= g < b.g_end:
            off = g - b.g_start
            return b.t_start + off if model.strand == "+" else b.t_start + (b.length - 1 - off)
    return None


def project_genome_to_cdna(model: CdnaModel, cigartuples, ref_start, query_seq):
    """Project a genome alignment of a read that splices cleanly along `model` into a
    cDNA-local (cigar, cdna_ref_start, cdna_query). Returns None unless the read's exon
    structure EXACTLY matches `model`'s block tiling — i.e. every ref-consuming base maps
    into an exon block AND consecutive ref bases map to consecutive transcript coords
    (so a read using a non-MANE junction, where its intron length differs from the block
    gap, is rejected rather than silently remapped onto MANE coords).

    The cDNA-local CIGAR drops genome N-ops (they become block gaps) and merges the
    flanking match runs; lift() should re-insert exactly those N-ops.
    """
    step = 1 if model.strand == "+" else -1
    g = ref_start
    genome_ops: List[Tuple[int, int]] = []
    first_ref_g = None
    last_ref_g = None
    prev_t = None  # transcript coord of the last ref base emitted (contiguity tracker)
    for op, ln in cigartuples:
        if op == N:
            g += ln
            continue
        if op in (M, EQ, X, D):
            t0 = g2t(model, g)
            t1 = g2t(model, g + ln - 1)
            if t0 is None or t1 is None:
                return None
            # run must stay within one exon block: end is start +/- (ln-1)
            if t1 != t0 + step * (ln - 1):
                return None
            # contiguity with the previous ref base (across N's / I's): no transcript gap
            if prev_t is not None and t0 != prev_t + step:
                return None
            prev_t = t1
            if first_ref_g is None:
                first_ref_g = g
            last_ref_g = g + ln - 1
            genome_ops.append((op, ln))
            g += ln
        elif op in (I, S, H):
            genome_ops.append((op, ln))
        else:
            return None
    if first_ref_g is None:
        return None

    # 2. merge adjacent identical ops (N removal fused exon matches)
    merged = _merge_keep_clip(genome_ops)

    if model.strand == "+":
        cdna_cigar = merged
        cdna_ref_start = g2t(model, first_ref_g)
        cdna_query = query_seq
    else:
        # transcript order is reverse of genome order; revcomp query, reverse ops,
        # swap leading/trailing clips (handled by reversal)
        cdna_cigar = list(reversed(merged))
        cdna_ref_start = g2t(model, last_ref_g)   # transcript-5' = highest genome coord
        cdna_query = revcomp(query_seq)
    return cdna_cigar, cdna_ref_start, cdna_query


def _merge_keep_clip(cig):
    out = []
    for op, ln in cig:
        if ln == 0:
            continue
        if out and out[-1][0] == op and op != N:
            out[-1] = (op, out[-1][1] + ln)
        else:
            out.append((op, ln))
    return out


def norm_cigar(cig):
    """Canonical form for comparison: merge adjacent identical (non-N) ops, drop zero-len."""
    return _merge_keep_clip(list(cig))
