"""Tier 1 — single-pass PARE census over MAPPED reads: 5'P cut-site pileup,
read-length distributions, and genome-aware junction-anchored poly(A) capture.

PARE is SENSE (read base 1 = the RNA 5'-monophosphate; the read maps in the
RNA's own orientation) — the mirror of the NET-seq/QuantSeq-REV antisense
geometry:

  + gene -> read maps FORWARD; 5'P = reference_start (LEFT);
            RNA 3' end (potential pA junction) = reference_end (RIGHT);
            poly-A = A-run on the RIGHT (stop_base 'A').
  - gene -> read maps REVERSE; 5'P = reference_end - 1 (RIGHT);
            RNA 3' end = reference_start (LEFT);
            poly-A = T-run on the LEFT (stop_base 'T').

The poly(A) metrics (``walkback_and_count``, ``softclip_run``,
``softclip_run_nt``) are IMPORTED from :mod:`rectify.core.netseq_cpa.pileup`
unchanged — they are (side, stop_base)-parameterized; only the geometry mapping
is PARE-specific. The CPA pileup therefore shares ``PILEUP_COLUMNS`` with the
NET-seq arm and feeds :mod:`rectify.core.netseq_cpa.concordance` as-is.

The 5'P pileup enforces the workspace-mandatory strand-aware 5'-clip gate: a
read whose RNA-5' terminus is soft/hard-clipped has an alignment-artifact 5'
end, so it must NOT contribute a 5'P call (default ``max_fivep_clip=0``); the
excluded fraction is always reported. The CPA side takes all mapped reads —
the pA junction is validated by its own genome-aware non-templated A evidence,
and the per-read parquet carries ``five_p_clip`` so downstream analyses can
gate further.

Outputs (all in ``out_dir``):
  cpa_pileup.tsv.gz        one row per (chrom, corrected pA junction, strand);
                           columns = netseq_cpa PILEUP_COLUMNS
  fivep_pileup.tsv.gz      chrom, pos, strand, n_reads, sample  (clip-gated)
  <sample>.5p_plus.bedgraph / <sample>.5p_minus.bedgraph
                           chrom, 0-based start, 1-based end, count — one row
                           per 5'P position (END pileup, never coverage), RNA
                           strand = mapped strand (sense protocol), sorted by
                           BAM header reference order then position
  read_lengths.tsv.gz      metric, value, count  (query_len, aligned_len,
                           clip5, clip3, oaNT over mapped primary reads)
"""
from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, Mapping, Optional, Tuple

import pysam

from ..netseq_cpa.pileup import (
    LEFT,
    MAXLEN_DEFAULT,
    PILEUP_COLUMNS,
    RIGHT,
    softclip_run,
    softclip_run_nt,
    walkback_and_count,
)

FIVEP_COLUMNS = ["chrom", "pos", "strand", "n_reads", "sample"]
# ``unmapped_len`` and ``mapq`` exist so short-read loss and ambiguous placement
# are VISIBLE rather than silent: closely spaced cleavage sites are resolved only
# by short, uniquely placed reads, so a pipeline that quietly drops 20-nt tags —
# or keeps them at MAPQ 0 — would understate exactly the sites we are hunting.
LENGTH_METRICS = ("query_len", "aligned_len", "clip5", "clip3", "oaNT",
                  "mapq", "unmapped_len")
MAX_FIVEP_CLIP_DEFAULT = 0
AMBIGUOUS_MAPQ = 3          # below this a placement is effectively multi-mapping


def geometry_for_read(read: pysam.AlignedSegment) -> Tuple[str, str, str]:
    """Return ``(gene_strand, pA_side, stop_base)`` for a SENSE PARE read."""
    if read.is_reverse:                 # - gene (sense)
        return "-", LEFT, "T"
    return "+", RIGHT, "A"              # + gene


def five_prime_of(read: pysam.AlignedSegment) -> Tuple[int, int]:
    """Strand-aware ``(0-based 5'P position, RNA-5'-side clip length)``.

    The clip length counts soft AND hard clip (ops 4/5) at the RNA 5' terminus;
    a non-zero value means the read's base 1 was not placed by the aligner, so
    its apparent 5'P is an alignment artifact, not a molecule end.
    """
    ct = read.cigartuples or []
    if read.is_reverse:
        pos = read.reference_end - 1
        op, ln = ct[-1] if ct else (0, 0)
    else:
        pos = read.reference_start
        op, ln = ct[0] if ct else (0, 0)
    return pos, (ln if op in (4, 5) else 0)


def _rna3_clip_len(read: pysam.AlignedSegment, side: str) -> int:
    """Soft/hard clip length at the RNA 3' side (the pA-junction side)."""
    ct = read.cigartuples or []
    if not ct:
        return 0
    op, ln = ct[-1] if side == RIGHT else ct[0]
    return ln if op in (4, 5) else 0


def _write_bedgraphs(
    agg_5p: Dict[Tuple[str, int, str], int],
    out_dir: Path,
    sample: str,
    ref_order: Dict[str, int],
) -> Dict[str, int]:
    """5'P bedgraphs, one row per position: chrom, pos, pos+1, count."""
    rows = {"+": [], "-": []}
    for (chrom, pos, strand), n in agg_5p.items():
        rows[strand].append((ref_order.get(chrom, len(ref_order)), chrom, pos, n))
    stats = {}
    for strand, suffix in (("+", "plus"), ("-", "minus")):
        path = out_dir / f"{sample}.5p_{suffix}.bedgraph"
        rows[strand].sort()
        with open(path, "w") as fh:
            for _o, chrom, pos, n in rows[strand]:
                fh.write(f"{chrom}\t{pos}\t{pos + 1}\t{n}\n")
        stats[f"bedgraph_{suffix}_rows"] = len(rows[strand])
    return stats


def pare_pileup(
    bam_path: str | Path,
    genome_path: str | Path,
    out_dir: str | Path,
    sample: str,
    *,
    chrom_map: Optional[Mapping[str, str]] = None,
    limit: int = 0,
    maxlen: int = MAXLEN_DEFAULT,
    max_fivep_clip: int = MAX_FIVEP_CLIP_DEFAULT,
    reads_parquet=None,
    cpa_set: Optional[set] = None,
) -> Dict[str, object]:
    """Single Tier-1 pass over the mapped primary reads of ``bam_path``.

    Returns summary stats: ``reads_seen``, ``reads_used``, ``fivep_used``,
    ``fivep_clip_excluded``, ``fivep_clip_fraction``, ``read3p_positions``,
    ``fivep_positions``, ``oaNT_ge2_reads``, ``sum_oaNT`` + bedgraph row counts.

    ``reads_parquet`` (a :class:`rectify.core.pare.reads_parquet.PareReadWriter`)
    receives one ``tier="mapped"`` record per used read — including reads with
    ``oaNT_tail_len == 0`` — carrying ``five_p_pos``/``five_p_clip`` so the
    per-read (5'P, pA-junction) pairing is queryable without re-walking the BAM.
    ``cpa_set`` (from :func:`rectify.core.netseq_cpa.concordance.load_cpa_set`)
    sets ``at_cpa``; ``None`` -> ``at_cpa=None`` (no DRS map = unknowable).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    chrom_map = dict(chrom_map or {})
    fa = pysam.FastaFile(str(genome_path))
    seq_cache: Dict[str, str] = {}
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    ref_order = {chrom_map.get(r.split()[0], r.split()[0]): i
                 for i, r in enumerate(bam.references)}

    agg_cpa: Dict[Tuple[str, int, str], list] = defaultdict(lambda: [0] * 10)
    agg_5p: Dict[Tuple[str, int, str], int] = defaultdict(int)
    hist: Dict[str, Dict[int, int]] = {m: defaultdict(int) for m in LENGTH_METRICS}
    n_seen = n_used = n_5p = n_5p_excl = oaNT2 = sum_oaNT = n_unmapped = 0
    n_ambig = 0
    try:
        for read in bam.fetch(until_eof=True):
            n_seen += 1
            if read.is_secondary or read.is_supplementary:
                continue
            if read.is_unmapped or read.reference_end is None:
                # census the losses: a short-read floor that silently discards
                # 20-nt tags would erase closely spaced cleavage sites
                hist["unmapped_len"][len(read.query_sequence or "")] += 1
                n_unmapped += 1
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
            oa = min(n_nt + scrun, maxlen)
            oa_nt = min(n_nt + scrun_nt, maxlen)
            tot = min(n_tmpl + n_nt + scrun, maxlen)

            # ------------------------------------------------ CPA (pA junction)
            c = agg_cpa[(chrom, corr, gstrand)]
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

            # ------------------------------------------------ 5'P (clip-gated)
            fp_pos, fp_clip = five_prime_of(read)
            if fp_clip <= max_fivep_clip:
                agg_5p[(chrom, fp_pos, gstrand)] += 1
                n_5p += 1
            else:
                n_5p_excl += 1

            # ------------------------------------------------ read-length census
            qlen = read.query_length or (len(read.query_sequence or ""))
            hist["query_len"][qlen] += 1
            hist["aligned_len"][read.query_alignment_length] += 1
            hist["clip5"][fp_clip] += 1
            hist["clip3"][_rna3_clip_len(read, side)] += 1
            hist["oaNT"][oa_nt] += 1
            mq = int(read.mapping_quality or 0)
            hist["mapq"][mq] += 1
            if mq < AMBIGUOUS_MAPQ:
                n_ambig += 1

            if reads_parquet is not None:
                at_cpa = ((chrom, gstrand, corr) in cpa_set
                          if cpa_set is not None else None)
                reads_parquet.add(
                    read_id=read.query_name.split()[0],
                    chrom=chrom, cpa_pos=int(corr), gene_strand=gstrand,
                    oaNT_tail_len=int(oa_nt), at_cpa=at_cpa, tier="mapped",
                    mapq=int(read.mapping_quality or 0),
                    five_p_pos=int(fp_pos), five_p_clip=int(fp_clip),
                )
            n_used += 1
            if limit and n_seen >= limit:
                break
    finally:
        bam.close()
        fa.close()

    with gzip.open(out_dir / "cpa_pileup.tsv.gz", "wt") as fh:
        fh.write("\t".join(PILEUP_COLUMNS) + "\n")
        for (chrom, pos, strand), c in agg_cpa.items():
            fh.write(f"{chrom}\t{pos}\t{strand}\t"
                     + "\t".join(str(x) for x in c) + f"\t{sample}\n")

    with gzip.open(out_dir / "fivep_pileup.tsv.gz", "wt") as fh:
        fh.write("\t".join(FIVEP_COLUMNS) + "\n")
        for (chrom, pos, strand), n in sorted(
                agg_5p.items(),
                key=lambda kv: (ref_order.get(kv[0][0], len(ref_order)), kv[0][1])):
            fh.write(f"{chrom}\t{pos}\t{strand}\t{n}\t{sample}\n")

    with gzip.open(out_dir / "read_lengths.tsv.gz", "wt") as fh:
        fh.write("metric\tvalue\tcount\n")
        for metric in LENGTH_METRICS:
            for value in sorted(hist[metric]):
                fh.write(f"{metric}\t{value}\t{hist[metric][value]}\n")

    stats: Dict[str, object] = {
        "reads_seen": n_seen,
        "reads_used": n_used,
        "reads_unmapped": n_unmapped,
        "ambiguous_reads": n_ambig,
        "ambiguous_fraction": (n_ambig / n_used) if n_used else 0.0,
        "fivep_used": n_5p,
        "fivep_clip_excluded": n_5p_excl,
        "fivep_clip_fraction": (n_5p_excl / n_used) if n_used else 0.0,
        # NOT a count of CPA calls: this is every distinct walkback-corrected
        # READ 3' END. In a library whose 3' end is set by random priming or a
        # fixed-length enzymatic cut it is not a 3'-end measurement at all. The
        # CPA-evidence number is oaNT_ge2_reads (non-templated poly(A) support).
        "read3p_positions": len(agg_cpa),
        "fivep_positions": len(agg_5p),
        "oaNT_ge2_reads": oaNT2,
        "sum_oaNT": sum_oaNT,
    }
    stats.update(_write_bedgraphs(agg_5p, out_dir, sample, ref_order))
    return stats
