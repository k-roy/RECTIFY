"""Rescue tier — recover poly-A-dominated UNMAPPED PARE reads (sense).

A PARE fragment whose body-to-tail ratio is low fails the aligner's min-ratio
filter and lands in the unmapped pile, exactly like the NET-seq case — but in
SENSE orientation: the poly-A run is at the read's 3' END (the NET-seq rescue
trims a 5' poly-T instead). This tier:

  1. :func:`trim_unmapped_polya` — for each unmapped read, strip trailing
     last-cycle ``N``s, then the 3' poly-A run, and emit the genomic anchor as
     FASTQ with the poly-A length encoded in the read name (``_paN``). Keeps
     only reads with an anchor >= ``min_anchor`` and poly-A >= ``min_pa``.
  2. (the caller re-maps the anchors with bbmap)
  3. :func:`quantify_rescue` — for each re-mapped anchor, the poly-A-side end
     is the pA junction (CPA candidate) and the OTHER end is the 5'P cut site,
     so each rescued read yields a single-molecule (5'P, pA-junction) pair.
     Reports the at-DRS-CPA fraction when a cluster map is supplied.

Anchor-side geometry after re-mapping (anchor is the read minus its 3' tail):
  FORWARD anchor -> + gene: pA junction = reference_end - 1, 5'P = reference_start.
  REVERSE anchor -> - gene: pA junction = reference_start, 5'P = reference_end - 1.

``min_anchor`` defaults LOWER than the NET-seq arm's 20: PARE inserts are short
(51-nt reads with heavy 3'-adapter read-through), so a tailed read's body is
often < 20 nt; ambiguous placements are handled by the ``min_mapq`` gate at
quantify time instead. Report-and-gate, never silently drop.
"""
from __future__ import annotations

import gzip
import re
import statistics as st
from collections import defaultdict
from pathlib import Path
from typing import Dict, Mapping, Optional

import pysam

from ..netseq_cpa.concordance import REG_DEFAULT

MIN_PA_DEFAULT = 8
MIN_ANCHOR_DEFAULT = 15
PALEN_CAP = 60


def trim_unmapped_polya(
    bam_path: str | Path,
    out_fq: str | Path,
    *,
    min_pa: int = MIN_PA_DEFAULT,
    min_anchor: int = MIN_ANCHOR_DEFAULT,
) -> Dict[str, int]:
    """Trim the 3' poly-A run from unmapped reads; emit genomic anchors as FASTQ.

    Returns ``{unmapped, rescued}``. The poly-A length is encoded as ``_pa<N>``
    appended to the read name so :func:`quantify_rescue` can read it back after
    re-mapping.
    """
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    n_um = n_kept = 0
    with gzip.open(out_fq, "wt") as out:
        for read in bam.fetch(until_eof=True):
            if not read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            n_um += 1
            s = read.query_sequence
            if not s:
                continue
            s = s.rstrip("N")                    # strip trailing last-cycle N's
            m = re.search("A+$", s)
            pa = len(m.group()) if m else 0
            if pa < min_pa:
                continue
            anchor = s[:len(s) - pa]
            if len(anchor) < min_anchor:
                continue
            q = read.query_qualities
            qa = ("".join(chr(x + 33) for x in q[:len(anchor)])
                  if q is not None else "I" * len(anchor))
            name = read.query_name.split()[0]    # QNAME may carry SRA description
            out.write(f"@{name}_pa{pa}\n{anchor}\n+\n{qa}\n")
            n_kept += 1
    bam.close()
    return {"unmapped": n_um, "rescued": n_kept}


def quantify_rescue(
    bam_path: str | Path,
    drs_clusters: Optional[str | Path] = None,
    *,
    min_mapq: int = 3,
    reg: int = REG_DEFAULT,
    chrom_map: Optional[Mapping[str, str]] = None,
    reads_parquet=None,
) -> Dict[str, object]:
    """Quantify re-mapped rescue anchors; each yields a (5'P, pA-junction) pair.

    Returns ``mapped``, ``at_cpa``, ``frac_at_cpa`` (0.0 with no DRS map) and
    poly-A length summaries. ``reads_parquet`` receives one ``tier="rescued"``
    record per anchor clearing ``min_mapq`` (``at_cpa`` null without a map;
    the parquet ``mapq`` column FLOORS at ``min_mapq`` for rescued rows).
    """
    chrom_map = dict(chrom_map or {})
    cpa: Optional[Dict[tuple, set]] = None
    if drs_clusters is not None:
        import csv
        cpa = defaultdict(set)
        with open(drs_clusters) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                c, s, m = r["chrom"], r["strand"], int(float(r["modal_position"]))
                for g in range(m - reg, m + reg + 1):
                    cpa[(c, s)].add(g)

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    n_map = n_atcpa = 0
    palen, palen_atcpa = [], []
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        mm = re.search(r"_pa(\d+)$", read.query_name)
        if not mm:
            continue
        L = int(mm.group(1))
        short = read.reference_name.split()[0]
        chrom = chrom_map.get(short, short)
        if read.is_reverse:                       # - gene (sense anchor)
            strand = "-"
            pos, five_p = read.reference_start, read.reference_end - 1
        else:                                     # + gene
            strand = "+"
            pos, five_p = read.reference_end - 1, read.reference_start
        n_map += 1
        palen.append(min(L, PALEN_CAP))
        at_cpa: Optional[bool] = None
        if cpa is not None:
            at_cpa = pos in cpa.get((chrom, strand), ())
            if at_cpa:
                n_atcpa += 1
                palen_atcpa.append(min(L, PALEN_CAP))
        if reads_parquet is not None:
            reads_parquet.add(
                read_id=re.sub(r"_pa\d+$", "", read.query_name.split()[0]),
                chrom=chrom, cpa_pos=int(pos), gene_strand=strand,
                oaNT_tail_len=int(min(L, PALEN_CAP)), at_cpa=at_cpa,
                tier="rescued", mapq=int(read.mapping_quality or 0),
                five_p_pos=int(five_p), five_p_clip=0,
            )
    bam.close()

    return {
        "mapped": n_map,
        "at_cpa": n_atcpa,
        "frac_at_cpa": (n_atcpa / n_map) if n_map else 0.0,
        "palen_median": st.median(palen) if palen else 0,
        "palen_mean": st.mean(palen) if palen else 0.0,
        "palen_max": max(palen) if palen else 0,
        "palen_atcpa_median": st.median(palen_atcpa) if palen_atcpa else 0,
        "palen_atcpa_mean": st.mean(palen_atcpa) if palen_atcpa else 0.0,
    }
