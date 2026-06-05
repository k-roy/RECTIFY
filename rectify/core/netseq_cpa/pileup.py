"""Tier 1 — genome-aware, junction-anchored poly(A) pileup on MAPPED NET-seq reads.

The non-circular discovery + concordance tier. For each mapped read, walk the
RNA 3' end back through the poly-A tract and measure the non-templated poly(A)
length two ways, both genome-aware:

  * **walkback non-templated** (``n_nt``): terminal read A's aligned over
    NON-A genome (A-tract-absorbed tail the aligner placed as matches/mismatches).
  * **genome-aware soft-clip** (``softclip_run_nt``): the junction-adjacent A-run
    in the RNA-3' soft-clip MINUS the genomic A-run the aligner clipped at that
    junction — only NON-templated clipped A's count.

  oaNT = n_nt + softclip_run_nt  (the decontaminated, internal-priming-immune tail).

NET-seq is antisense; geometry is identical to QuantSeq REV:

  + gene -> read maps REVERSE; RNA 3' end (CPA) = reference_end (RIGHT);
            poly-A = A-run on the RIGHT (stop_base 'A').
  - gene -> read maps FORWARD; RNA 3' end (CPA) = reference_start (LEFT);
            poly-A = T-run on the LEFT (stop_base 'T').

``walkback_and_count`` is a faithful port of the RECTIFY 3' walkback that
additionally exposes the templated/non-templated terminal counts the canonical
:func:`rectify.core.correct.protocols.netseq.walkback_netseq` wrapper does not
return. ``corrected_position_matches_walkback_netseq`` is a parity check used by
the test suite to lock this port to the canonical wrapper's corrected CPA.

Output (one row per (chrom, corrected 0-based CPA, gene strand)):
  chrom, pos, strand, n_reads,
  n_oa_ge1/2/3, n_tot_ge3, sum_oa            (LEGACY genome-blind soft-clip),
  n_oaNT_ge1/2/3, sum_oaNT, sample           (GENOME-AWARE — the validated metric).

Validated 2026-06-03: register peak at offset 0, ~2x over matched gene-body
background; ``n_oaNT_ge2`` is the headline metric. See
``analyses/netseq_oligoA_cpa_20260603/``.
"""
from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Mapping, Optional, Tuple

import pysam

RIGHT, LEFT = "right", "left"
MAXLEN_DEFAULT = 25

PILEUP_COLUMNS = [
    "chrom", "pos", "strand", "n_reads",
    "n_oa_ge1", "n_oa_ge2", "n_oa_ge3", "n_tot_ge3", "sum_oa",
    "n_oaNT_ge1", "n_oaNT_ge2", "n_oaNT_ge3", "sum_oaNT", "sample",
]


def geometry_for_read(read: pysam.AlignedSegment) -> Tuple[str, str, str]:
    """Return ``(gene_strand, side, stop_base)`` for an antisense NET-seq read."""
    if read.is_reverse:                 # + gene (antisense)
        return "+", RIGHT, "A"
    return "-", LEFT, "T"               # - gene


def walkback_and_count(
    read: pysam.AlignedSegment, chrom_seq: str, side: str, stop_base: str
) -> Tuple[int, int, int, int]:
    """Walk the 3' end back through the ``stop_base`` tract.

    Returns ``(original, corrected, n_templated, n_nontemplated)``. Faithful to
    the RECTIFY 3' walkback: scan aligned (read, ref) pairs inward from the 3'
    end; count terminal read ``stop_base`` over genomic ``stop_base`` as
    templated and over non-``stop_base`` genome as non-templated; stop at the
    first position where read==genome and genome != stop_base (the true CPA).

    The single gate that skips walkback: terminal base is NOT ``stop_base`` AND
    matches the genome -> already correctly anchored (returns original).
    """
    original = read.reference_end - 1 if side == RIGHT else read.reference_start
    pairs = [
        (qp, rp)
        for qp, rp in read.get_aligned_pairs(matches_only=False)
        if qp is not None and rp is not None
    ]
    if not pairs:
        return original, original, 0, 0
    scan = list(reversed(pairs)) if side == RIGHT else pairs
    qs = read.query_sequence
    if qs is None:
        return original, original, 0, 0
    fqp, frp = scan[0]
    if qs[fqp].upper() != stop_base and qs[fqp].upper() == chrom_seq[frp].upper():
        return original, original, 0, 0
    corrected = original
    n_tmpl = n_nt = 0
    for qp, rp in scan:
        rb = qs[qp].upper()
        gb = chrom_seq[rp].upper()
        if rb == gb and rb != stop_base:
            corrected = rp
            break
        if rb == stop_base:
            if gb == stop_base:
                n_tmpl += 1
            else:
                n_nt += 1
    return original, corrected, n_tmpl, n_nt


def softclip_run(read: pysam.AlignedSegment, side: str, stop_base: str) -> int:
    """Longest run of ``stop_base`` in the RNA-3' soft-clip. GENOME-BLIND
    (legacy): counts clipped genomic A-tracts as poly-A. Kept for comparison."""
    ct = read.cigartuples
    seq = read.query_sequence
    if not ct or not seq:
        return 0
    seq = seq.upper()
    if side == RIGHT:
        op, ln = ct[-1]
        clip = seq[len(seq) - ln:] if op == 4 else ""
    else:
        op, ln = ct[0]
        clip = seq[:ln] if op == 4 else ""
    if not clip:
        return 0
    return max((len(m.group()) for m in re.finditer(stop_base + "+", clip)), default=0)


def softclip_run_nt(
    read: pysam.AlignedSegment, chrom_seq: str, side: str, stop_base: str
) -> int:
    """GENOME-AWARE soft-clip poly-A: junction-adjacent ``stop_base`` run in the
    RNA-3' clip MINUS the genomic ``stop_base`` run the aligner clipped at that
    junction. Only NON-templated clipped bases count."""
    ct = read.cigartuples
    seq = read.query_sequence
    if not ct or not seq:
        return 0
    seq = seq.upper()
    if side == RIGHT:
        op, ln = ct[-1]
        if op != 4:
            return 0
        clip = seq[len(seq) - ln:]
        m = re.match(stop_base + "+", clip)
        clip_run = len(m.group()) if m else 0
        g = read.reference_end
        gr = 0
        while g < len(chrom_seq) and chrom_seq[g].upper() == stop_base:
            gr += 1
            g += 1
    else:
        op, ln = ct[0]
        if op != 4:
            return 0
        clip = seq[:ln]
        m = re.search(stop_base + "+$", clip)
        clip_run = len(m.group()) if m else 0
        g = read.reference_start - 1
        gr = 0
        while g >= 0 and chrom_seq[g].upper() == stop_base:
            gr += 1
            g -= 1
    return max(0, clip_run - gr)


def corrected_position_matches_walkback_netseq(
    read: pysam.AlignedSegment, chrom_seq: str
) -> bool:
    """Parity check: this module's ``walkback_and_count`` corrected CPA equals
    the canonical :func:`...protocols.netseq.walkback_netseq` corrected CPA.

    Used by the test suite to guard against the port drifting from the wrapper.
    """
    from rectify.core.correct.protocols.netseq import walkback_netseq

    gstrand, side, stop = geometry_for_read(read)
    _, corr_here, _, _ = walkback_and_count(read, chrom_seq, side, stop)
    _orig, corr_wrap, _applied, _gs = walkback_netseq(read, chrom_seq)
    return corr_here == corr_wrap and gstrand == _gs


def walkback_pileup(
    bam_path: str | Path,
    genome_path: str | Path,
    out_path: str | Path,
    sample: str,
    *,
    chrom_map: Optional[Mapping[str, str]] = None,
    limit: int = 0,
    maxlen: int = MAXLEN_DEFAULT,
    reads_parquet=None,
    cpa_set: Optional[set] = None,
) -> Dict[str, int]:
    """Tier-1 genome-aware poly(A) pileup over the mapped reads of ``bam_path``.

    Writes a gzipped TSV (``PILEUP_COLUMNS``) to ``out_path``; returns summary
    stats ``{reads_seen, reads_used, positions, oaNT_ge2_reads, sum_oaNT}``.

    ``chrom_map`` (e.g. :data:`...netseq_cpa.NCBI_TO_CHROM`) translates BAM
    reference names to standard chrI-style names so the pileup chrom column
    matches the DRS cluster map. ``None`` -> identity (names already standard).

    ``reads_parquet`` (a :class:`...reads_parquet.ParquetReadWriter`), when given,
    receives one ``tier="mapped"`` per-read record for every used read — emitted
    just before the ``reads_used`` increment so the mapped-row count equals
    ``reads_used`` exactly. ``cpa_set`` (from :func:`...concordance.load_cpa_set`)
    sets each record's ``at_cpa``; pass ``None`` to emit ``at_cpa=None`` (no DRS
    map -> unknowable, distinct from not-at-CPA). The per-read emit is purely
    additive: the TSV aggregation below is unchanged, preserving pileup parity.
    """
    chrom_map = dict(chrom_map or {})
    fa = pysam.FastaFile(str(genome_path))
    seq_cache: Dict[str, str] = {}
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    agg: Dict[Tuple[str, int, str], list] = defaultdict(lambda: [0] * 10)
    n_seen = n_used = oaNT2 = sum_oaNT = 0
    try:
        for read in bam.fetch(until_eof=True):
            n_seen += 1
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.reference_end is None):
                continue
            short = read.reference_name.split()[0]
            chrom = chrom_map.get(short, short)
            gstrand, side, stop = geometry_for_read(read)
            if short not in seq_cache:
                seq_cache[short] = fa.fetch(short)
            cs = seq_cache[short]
            try:
                orig, corr, n_tmpl, n_nt = walkback_and_count(read, cs, side, stop)
            except (IndexError, ValueError):
                continue
            scrun = softclip_run(read, side, stop)             # genome-blind (legacy)
            scrun_nt = softclip_run_nt(read, cs, side, stop)   # genome-aware
            oa = min(n_nt + scrun, maxlen)                     # legacy tail
            oa_nt = min(n_nt + scrun_nt, maxlen)               # decontaminated tail
            tot = min(n_tmpl + n_nt + scrun, maxlen)           # incl templated A-run
            c = agg[(chrom, corr, gstrand)]
            c[0] += 1
            if oa >= 1:
                c[1] += 1
            if oa >= 2:
                c[2] += 1
            if oa >= 3:
                c[3] += 1
            if tot >= 3:
                c[4] += 1
            c[5] += oa
            if oa_nt >= 1:
                c[6] += 1
            if oa_nt >= 2:
                c[7] += 1
                oaNT2 += 1
            if oa_nt >= 3:
                c[8] += 1
            c[9] += oa_nt
            sum_oaNT += oa_nt
            if reads_parquet is not None:
                at_cpa = ((chrom, gstrand, corr) in cpa_set
                          if cpa_set is not None else None)
                reads_parquet.add(
                    read_id=read.query_name.split()[0],
                    chrom=chrom, cpa_pos=int(corr), gene_strand=gstrand,
                    oaNT_tail_len=int(oa_nt), at_cpa=at_cpa, tier="mapped",
                    mapq=int(read.mapping_quality or 0),
                )
            n_used += 1
            if limit and n_seen >= limit:
                break
    finally:
        bam.close()
        fa.close()

    n_emit = 0
    with gzip.open(out_path, "wt") as fh:
        fh.write("\t".join(PILEUP_COLUMNS) + "\n")
        for (chrom, pos, strand), c in agg.items():
            fh.write(
                f"{chrom}\t{pos}\t{strand}\t"
                + "\t".join(str(x) for x in c)
                + f"\t{sample}\n"
            )
            n_emit += 1
    return {
        "reads_seen": n_seen,
        "reads_used": n_used,
        "positions": n_emit,
        "oaNT_ge2_reads": oaNT2,
        "sum_oaNT": sum_oaNT,
    }
