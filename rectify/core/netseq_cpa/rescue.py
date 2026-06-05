"""Rescue tier — recover poly-A-dominated UNMAPPED NET-seq reads.

Standard NET-seq pipelines lose the CPA intermediates whose poly-A tail is a
large fraction of the read: the aligner's min-ratio filter rejects them (aligned
fraction too low) and dumps them in the unmapped pile (27-98x poly-A-enriched).
This tier recovers them:

  1. :func:`trim_unmapped_polyt` — for each unmapped read, strip a leading
     first-cycle ``N`` (Xi etc.), then the 5' poly-T run (= RNA poly-A in the
     antisense first-strand cDNA), and emit the genomic anchor as FASTQ with the
     poly-A length encoded in the read name (``_paN``). Keeps only reads with a
     clean anchor >= ``min_anchor`` and poly-T >= ``min_pt`` (poly-A-only junk
     with no anchor is dropped).
  2. (the caller re-maps the anchors with bbmap/bowtie2)
  3. :func:`quantify_rescue` — for each re-mapped anchor, the poly-T-side end is
     the RNA 3' end = CPA (antisense). Report how many land on known DRS CPAs
     (+- ``reg``) and the poly-A length distribution.

Validated 2026-06-03 across 4 labs / 3 protocols: 78-89% of recovered anchors
land on DRS CPAs (vs ~5% genome null), poly-A length median 12-13 nt.
"""
from __future__ import annotations

import csv
import gzip
import re
import statistics as st
from collections import defaultdict
from pathlib import Path
from typing import Dict, Mapping, Optional

import pysam

MIN_PT_DEFAULT = 8
MIN_ANCHOR_DEFAULT = 20
REG_DEFAULT = 5


def trim_unmapped_polyt(
    bam_path: str | Path,
    out_fq: str | Path,
    *,
    min_pt: int = MIN_PT_DEFAULT,
    min_anchor: int = MIN_ANCHOR_DEFAULT,
) -> Dict[str, int]:
    """Trim 5' poly-T from unmapped reads; emit genomic anchors as FASTQ.

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
            ln = len(s) - len(s.lstrip("N"))     # strip a leading first-cycle N
            s = s[ln:]
            m = re.match("T+", s)
            pt = len(m.group()) if m else 0
            if pt < min_pt:
                continue
            anchor = s[pt:]
            if len(anchor) < min_anchor:
                continue
            q = read.query_qualities
            off = ln + pt
            qa = "".join(chr(x + 33) for x in q[off:]) if q is not None else "I" * len(anchor)
            name = read.query_name.split()[0]    # QNAME may carry SRA description (spaces)
            out.write(f"@{name}_pa{pt}\n{anchor}\n+\n{qa}\n")
            n_kept += 1
    bam.close()
    return {"unmapped": n_um, "rescued": n_kept}


def quantify_rescue(
    bam_path: str | Path,
    drs_clusters: str | Path,
    *,
    min_mapq: int = 3,
    reg: int = REG_DEFAULT,
    chrom_map: Optional[Mapping[str, str]] = None,
    reads_parquet=None,
) -> Dict[str, object]:
    """Quantify re-mapped rescue anchors against the DRS CPA map.

    Returns a dict with ``mapped``, ``at_cpa``, ``frac_at_cpa``, and poly-A
    length summaries (``palen_median``/``mean``/``max`` over all recovered, and
    ``palen_atcpa_median``/``mean`` at DRS CPAs).

    ``reads_parquet`` (a :class:`...reads_parquet.ParquetReadWriter`), when given,
    receives one ``tier="rescued"`` record per re-mapped anchor that clears
    ``min_mapq`` — so the rescued-row count equals the returned ``mapped``, and
    the parquet ``mapq`` column FLOORS at ``min_mapq`` for rescued rows.
    """
    chrom_map = dict(chrom_map or {})
    cpa: Dict[tuple, set] = defaultdict(set)
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
        if read.is_reverse:
            strand, pos = "+", read.reference_end - 1
        else:
            strand, pos = "-", read.reference_start
        n_map += 1
        palen.append(min(L, 60))
        at_cpa = pos in cpa.get((chrom, strand), ())
        if at_cpa:
            n_atcpa += 1
            palen_atcpa.append(min(L, 60))
        if reads_parquet is not None:
            reads_parquet.add(
                read_id=re.sub(r"_pa\d+$", "", read.query_name.split()[0]),
                chrom=chrom, cpa_pos=int(pos), gene_strand=strand,
                oaNT_tail_len=int(min(L, 60)), at_cpa=at_cpa, tier="rescued",
                mapq=int(read.mapping_quality or 0),
            )
    bam.close()

    out: Dict[str, object] = {
        "mapped": n_map,
        "at_cpa": n_atcpa,
        "frac_at_cpa": (n_atcpa / n_map) if n_map else 0.0,
        "palen_median": st.median(palen) if palen else 0,
        "palen_mean": st.mean(palen) if palen else 0.0,
        "palen_max": max(palen) if palen else 0,
        "palen_atcpa_median": st.median(palen_atcpa) if palen_atcpa else 0,
        "palen_atcpa_mean": st.mean(palen_atcpa) if palen_atcpa else 0.0,
    }
    return out
