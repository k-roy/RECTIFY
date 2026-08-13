#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-read gene attribution sidecar (rbrowse spec `702`, design `703`).

WHY THIS EXISTS
---------------
Downstream viewers decide which 3' end "belongs to" which gene from the 3' end
alone.  For a readthrough molecule -- one that initiates in gene A and terminates
inside gene B -- that counts the molecule FOR B and subtracts it FROM A.  In a
CPA-depletion experiment readthrough rises by construction, so the error is
directional and correlated with the genotype being compared: it inflates the
apparent decrease of exactly the genes the experiment is about.

This module emits a per-READ sidecar that carries both ends, so the consumer can
keep a read with its initiating gene and still see where it terminated.

TWO CONSTRAINTS FROM THE SPEC, BOTH LOAD-BEARING
------------------------------------------------
1. **Never key output on an analyze ``cluster_id``.**  rbrowse re-clusters live
   in the browser from user-dragged sliders, so an analyze cluster boundary has
   no counterpart there and any join on it is unresolvable.  Output is per-read,
   joined on ``read_id`` (DRS) or ``umikey`` (UMI-collapsed cDNA) -- which
   survives re-clustering at any window by construction.

2. **Do not re-run ``correct``.**  An analyze run and a bundle built from
   *different* corrections have different 3' anchors, which has already blocked
   one deploy.  Everything here is computed from the corrected-reads TSV that the
   consumer itself loads, so the read set and the 3' anchors are identical by
   construction rather than by convention.  No BAM is opened.

THE 5' COORDINATE TRAP
----------------------
``five_prime_position`` is NOT always the observed 5' end.  The Cat3 rescue
(``extend_read_5prime_for_junction_rescue``) EXTENDS a read's 5' end upstream
across an annotated junction, and ``rescue_3ss_truncation`` writes the moved
coordinate into that column -- so for ``five_prime_rescued=1`` reads it is a
model-derived coordinate.  The spec's central claim ("the pore truncates from the
5' end, so a 5' end upstream of X PROVES the molecule began upstream of X") only
holds for an observed coordinate.  We therefore derive the 5' end from
``alignment_start``/``alignment_end`` (the raw input alignment, untouched by
correction) and report the rescue flag in ``origin5_evidence`` rather than
silently laundering it.

THE TERMINATION BOUNDARY (measured, not assumed)
------------------------------------------------
Deciding "did this read terminate inside gene G?" by containment in G's annotated
interval is wrong, and wrong in the dangerous direction: it INVENTS readthrough in
the wild type.  Measured on the anchor-away DRS cohort, containment-only takes
PLB3's WT escape rate from 0/845 to ~204/845, because a gene's dominant CPA
frequently sits a little past its annotated end.

A fixed downstream window does not fix it either.  Genome-wide over 5,607 genes
with >=20 WT reads (`planning/706`), the per-gene modal 3' end relative to the
annotated transcript end is distributed p50 -29, p75 0, p90 +114, p95 +281 --
but ~4% of genes sit far outside any fixed window (``YGL007C-A``: n=2,136 with
46% of reads at one base **+3,249 bp** past the annotated end).  Any constant
truncates exactly those genes.

So the boundary is per-gene and data-derived (`planning/707` rule B)::

    boundary(G) = max(annotated_end(G), observed_modal_CPA(G)) + GENOME_TAIL_BP

``observed_modal_CPA`` comes from the CONTROL libraries only: in a mutant, a
gene's own 3'-end distribution already contains the readthrough being measured, so
a mutant-derived reference would be circular.  ``GENOME_TAIL_BP`` is the measured
genome-wide p95, not a guess.

Rule B was selected by measurement against an independently-derived browser-side
number (WT 0/845, rna15 94/201, ysh1 143/355); rules using the gene's own q95/q99
inflate WT because a sharp modal peak underestimates the real spread.

Author: Kevin R. Roy
"""

from __future__ import annotations

import bisect
import csv
import gzip
import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)

#: Genome-wide p95 of the per-gene modal 3'-end offset past the annotated
#: transcript end, measured on the WT anchor-away DRS libraries
#: (`planning/706`, 5,607 genes with >=20 reads).  Used as the tolerance added on
#: top of each gene's OWN observed CPA -- it is a measured quantile, not a
#: hand-chosen window, and it is only ever additive to per-gene evidence.
GENOME_TAIL_BP = 281

#: Reads required before a gene's own modal CPA is trusted; below this the
#: annotated transcript end carries the boundary alone.
MIN_REF_READS = 20

#: Fixed-window comparison rule, retained ONLY so the consumer can A/B it against
#: the per-gene rule (Kevin, 2026-08-13: ship both, don't hard-code one).  This is
#: the CONTRAST column, never the default.
FIXED_WINDOW_BP = 300

SIDECAR_HEADER = [
    "join_key",
    "attr_genes",
    "attr_primary",
    "attr_rule",
    "origin5",
    "origin5_evidence",
    "escapes_gene",
    "escapes_gene_cpa",
    "region_class",
    "dist_to_stop",
    "ambiguity_range",
    "strand_evidence",
]

#: Columns read from the corrected-reads TSV.  ``alignment_start``/``alignment_end``
#: are required: they are the raw alignment and the only rescue-proof 5' source.
REQUIRED_COLUMNS = ("chrom", "strand", "alignment_start", "alignment_end")


# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------

class GeneIndex:
    """Position -> gene lookup over 0-based half-open intervals.

    Built from GFF ``gene`` features.  Note these are UTR-INCLUSIVE in the SGD
    R64 GFF -- ``gene`` is the union of the per-condition ``mRNA`` transcript
    models, so it already contains the 3'UTR and does NOT stop at the CDS.  The
    separate ``CDS`` extent is kept so ``dist_to_stop`` measures from the real
    stop codon rather than from the transcript end.
    """

    def __init__(self) -> None:
        self.spans: Dict[str, Tuple[str, int, int, str]] = {}
        self.cds_end: Dict[str, int] = {}
        self._idx: Dict[Tuple[str, str], Tuple[List[int], List[Tuple[int, int, str]]]] = {}

    @classmethod
    def from_gff(cls, path: str) -> "GeneIndex":
        self = cls()
        cds: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        opener = gzip.open if str(path).endswith(".gz") else open
        # errors='replace': SGD GFFs carry non-ASCII in free-text attributes and
        # a strict decode aborts the whole load over a description field.
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                if f[2] == "gene":
                    gid = _attr(f[8], "ID=")
                    if gid:
                        self.spans[gid] = (f[0], int(f[3]) - 1, int(f[4]), f[6])
                elif f[2] == "CDS":
                    for parent in _attr(f[8], "Parent=").split(","):
                        base = parent.split("_")[0]
                        if base:
                            cds[base].append((int(f[3]) - 1, int(f[4])))
        for gid, ivs in cds.items():
            if gid not in self.spans:
                continue
            strand = self.spans[gid][3]
            # 3'-most CDS base in transcription direction = the stop codon side
            self.cds_end[gid] = max(e for _, e in ivs) - 1 if strand == "+" \
                else min(s for s, _ in ivs)
        self._build()
        return self

    def _build(self) -> None:
        by: Dict[Tuple[str, str], List[Tuple[int, int, str]]] = defaultdict(list)
        for gid, (chrom, s, e, strand) in self.spans.items():
            if e > s:
                by[(chrom, strand)].append((s, e, gid))
        for key, items in by.items():
            items.sort()
            self._idx[key] = ([it[0] for it in items], items)

    def at(self, chrom: str, strand: str, pos: int) -> Optional[str]:
        """Gene containing ``pos`` on ``strand``, else None."""
        entry = self._idx.get((chrom, strand))
        if entry is None:
            return None
        starts, items = entry
        i = bisect.bisect_right(starts, pos) - 1
        # Walk back past genes that start earlier but still span pos; yeast genes
        # rarely nest, but overlapping annotations do exist.
        while i >= 0 and starts[i] > pos - _MAX_GENE_SPAN:
            s, e, gid = items[i]
            if s <= pos < e:
                return gid
            i -= 1
        return None

    def overlapping(self, chrom: str, strand: str, start: int, end: int) -> List[str]:
        """Gene ids overlapping [start, end) on ``strand``, largest overlap first."""
        entry = self._idx.get((chrom, strand))
        if entry is None:
            return []
        starts, items = entry
        i = bisect.bisect_right(starts, end) - 1
        hits = []
        while i >= 0 and starts[i] > start - _MAX_GENE_SPAN:
            s, e, gid = items[i]
            ov = min(end, e) - max(start, s)
            if ov > 0:
                hits.append((ov, gid))
            i -= 1
        hits.sort(key=lambda x: (-x[0], x[1]))
        return [gid for _, gid in hits]

    def annotated_end_offset_base(self, gid: str) -> int:
        """Genomic coordinate of the gene's annotated 3'-most base."""
        chrom, s, e, strand = self.spans[gid]
        return (e - 1) if strand == "+" else s


_MAX_GENE_SPAN = 20_000


def _attr(field: str, key: str) -> str:
    for kv in field.split(";"):
        if kv.startswith(key):
            return kv[len(key):]
    return ""


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def observed_5prime(strand: str, alignment_start: int, alignment_end: int) -> int:
    """Observed 5' end from the RAW alignment.

    Deliberately not ``five_prime_position``: that column is moved upstream by the
    Cat3 junction rescue, which would manufacture the very "began upstream"
    evidence the caller is trying to prove.  See the module docstring.
    """
    return alignment_start if strand == "+" else alignment_end - 1


def signed_offset_past(strand: str, pos: int, boundary_base: int) -> int:
    """Distance of ``pos`` past ``boundary_base``, in the transcription direction.

    Positive means downstream of the boundary for both strands.
    """
    return (pos - boundary_base) if strand == "+" else (boundary_base - pos)


# ---------------------------------------------------------------------------
# Reference CPA profile (control libraries only -- see module docstring)
# ---------------------------------------------------------------------------

def build_cpa_reference(
    control_tsvs: Sequence[str],
    genes: GeneIndex,
    min_reads: int = MIN_REF_READS,
) -> Dict[str, Dict[str, int]]:
    """Per-gene observed 3'-end profile, from CONTROL libraries only.

    For each read whose OBSERVED 5' end lies inside gene G, record how far past
    G's annotated 3' end the corrected 3' end sits.  The modal value is G's
    dominant CPA.  Returns ``{gene_id: {n, mode, q50, q95}}`` for genes clearing
    ``min_reads``.

    Control-only is not a preference: in a mutant, G's own 3'-end distribution
    already contains the readthrough being measured, so a mutant-derived reference
    would be circular.
    """
    offsets: Dict[str, List[int]] = defaultdict(list)
    for tsv in control_tsvs:
        for rec in _iter_corrected(tsv):
            gid = genes.at(rec["chrom"], rec["strand"], rec["five_p"])
            if gid is None:
                continue
            offsets[gid].append(
                signed_offset_past(rec["strand"], rec["three_p"],
                                   genes.annotated_end_offset_base(gid))
            )

    ref: Dict[str, Dict[str, int]] = {}
    for gid, vals in offsets.items():
        if len(vals) < min_reads:
            continue
        vals.sort()
        # Mode via a coarse bin then the exact peak inside it: raw single-base
        # modes are noisy at low depth, and DRS 3' ends spread over a few bases.
        binned = Counter((v // 10) * 10 for v in vals)
        best_bin = binned.most_common(1)[0][0]
        mode = Counter(v for v in vals if best_bin <= v < best_bin + 10).most_common(1)[0][0]
        ref[gid] = {
            "n": len(vals),
            "mode": mode,
            "q50": vals[len(vals) // 2],
            "q95": vals[min(len(vals) - 1, int(0.95 * len(vals)))],
        }
    logger.info("CPA reference: %d genes profiled from %d control libraries",
                len(ref), len(control_tsvs))
    return ref


def termination_boundary(
    gid: str,
    genes: GeneIndex,
    ref: Dict[str, Dict[str, int]],
) -> int:
    """Offset past the annotated end beyond which a 3' end has left the gene.

    Rule B (`planning/707`): ``max(annotated_end, observed_modal_CPA) + tail``.
    Falls back to the annotated end plus the same measured tail when the gene has
    too few control reads to profile -- which is itself a measured fallback, since
    the annotated end tracks the observed CPA to a median of -29 bp.
    """
    prof = ref.get(gid)
    base = max(0, prof["mode"]) if prof else 0
    return base + GENOME_TAIL_BP


# ---------------------------------------------------------------------------
# TSV streaming
# ---------------------------------------------------------------------------

def _iter_corrected(path: str):
    """Stream a corrected-reads TSV, yielding only what attribution needs.

    csv.reader rather than pandas: these files reach ~200 MB x many libraries and
    a DataFrame of all 38 columns does not fit the usual per-task memory.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        idx = {name: i for i, name in enumerate(header)}
        missing = [c for c in REQUIRED_COLUMNS if c not in idx]
        if missing:
            raise ValueError(
                f"{path}: corrected TSV lacks {missing}. alignment_start/"
                "alignment_end are required -- they are the raw alignment and the "
                "only rescue-proof source of the 5' end."
            )
        # DRS keys on read_id; the UMI-collapsed cDNA arm keys on umikey and has
        # no read_id at all.  Whichever is present becomes join_key.
        key_col = next((c for c in ("read_id", "umikey") if c in idx), None)
        if key_col is None:
            raise ValueError(f"{path}: no read_id or umikey column to join on")
        pos_col = "corrected_3prime" if "corrected_3prime" in idx else "corrected_position"
        if pos_col not in idx:
            raise ValueError(f"{path}: no corrected 3' end column")

        # Optional columns: absent from the 11-column UMI-collapsed cDNA schema.
        i_rescued = idx.get("five_prime_rescued")
        i_strand_ev = idx.get("strand_evidence")
        # Passed through, not used for attribution: `ambiguity_range` is the
        # per-read window in which the true CPA cannot be localised, so its
        # median over reads in view IS the assay's 3'-end resolution -- which is
        # the defensible way to choose a clustering window, as opposed to picking
        # the centre of the widest plateau in a parameter sweep (that plateau
        # sits where the parameter has already destroyed the structure).
        i_ambig = idx.get("ambiguity_range")
        for row in rd:
            try:
                strand = row[idx["strand"]]
                a_s = int(float(row[idx["alignment_start"]]))
                a_e = int(float(row[idx["alignment_end"]]))
                three_p = int(float(row[idx[pos_col]]))
            except (ValueError, IndexError):
                continue
            if strand not in ("+", "-"):
                continue
            yield {
                "key": row[idx[key_col]],
                "chrom": row[idx["chrom"]],
                "strand": strand,
                "aln_start": a_s,
                "aln_end": a_e,
                "three_p": three_p,
                "five_p": observed_5prime(strand, a_s, a_e),
                "rescued": i_rescued is not None and row[i_rescued] == "1",
                "strand_evidence": row[i_strand_ev] if i_strand_ev is not None else "",
                "ambiguity_range": row[i_ambig] if i_ambig is not None else "",
            }


# ---------------------------------------------------------------------------
# Per-read classification
# ---------------------------------------------------------------------------

def classify_read(rec: dict, genes: GeneIndex, ref: Dict[str, Dict[str, int]]) -> dict:
    """Classify one read into the sidecar columns.

    ``escapes_gene`` (fixed window) and ``escapes_gene_cpa`` (per-gene observed
    CPA) are emitted side by side on purpose: the fixed rule is the contrast the
    consumer A/Bs against, not a second opinion to be averaged in.
    """
    chrom, strand = rec["chrom"], rec["strand"]
    attr_genes = genes.overlapping(chrom, strand, rec["aln_start"], rec["aln_end"])
    attr_primary = attr_genes[0] if attr_genes else ""
    attr_rule = "body_max_overlap" if attr_primary else "none"

    init_gene = genes.at(chrom, strand, rec["five_p"])

    # origin5.  `initiating` requires positive 5'-end evidence (cDNA full-length
    # tier 2), which no schema we currently ship carries -- so it is structurally
    # unreachable here rather than guessed, which is what keeps a DRS read from
    # ever being called `initiating`.
    if not attr_primary:
        origin5, evidence = "unknown", "none"
    elif init_gene == attr_primary:
        origin5, evidence = "internal", "none"
    elif init_gene is None and _starts_upstream(rec, genes, attr_primary):
        origin5, evidence = "upstream_origin", "drs_5p_upstream"
    elif init_gene is not None and init_gene != attr_primary:
        origin5, evidence = "upstream_origin", "drs_5p_upstream"
    else:
        origin5, evidence = "unknown", "none"

    # A rescued 5' coordinate is model-derived; say so instead of laundering it.
    # The classification above already avoids five_prime_position, but the
    # consumer still needs to know the read carries a moved 5' end.
    if rec["rescued"] and evidence == "drs_5p_upstream":
        evidence = "five_prime_rescued"

    escapes = escapes_cpa = ""
    dist_to_stop = ""
    if init_gene:
        end_base = genes.annotated_end_offset_base(init_gene)
        past = signed_offset_past(strand, rec["three_p"], end_base)
        if past > FIXED_WINDOW_BP:
            escapes = init_gene
        if past > termination_boundary(init_gene, genes, ref):
            escapes_cpa = init_gene

    if attr_primary and attr_primary in genes.cds_end:
        dist_to_stop = signed_offset_past(
            strand, rec["three_p"], genes.cds_end[attr_primary])

    return {
        "join_key": rec["key"],
        "attr_genes": "|".join(attr_genes),
        "attr_primary": attr_primary,
        "attr_rule": attr_rule,
        "origin5": origin5,
        "origin5_evidence": evidence,
        "escapes_gene": escapes,
        "escapes_gene_cpa": escapes_cpa,
        # region_class is deliberately empty in v1: the vocabulary lives in
        # cluster-keyed cpa_clusters.tsv, and the spec forbids joining on
        # cluster_id.  Reimplementing it here would create a second source of
        # truth that could drift from rectify's own.
        "region_class": "",
        "dist_to_stop": dist_to_stop,
        "ambiguity_range": rec["ambiguity_range"],
        "strand_evidence": rec["strand_evidence"],
    }


def _starts_upstream(rec: dict, genes: GeneIndex, gid: str) -> bool:
    """True when the observed 5' end lies upstream of the gene it overlaps.

    Provable in DRS: the pore truncates from the 5' end, so truncation can only
    move a 5' end DOWNSTREAM.  A 5' end upstream of the gene therefore cannot be
    a truncation artifact.
    """
    chrom, s, e, strand = genes.spans[gid]
    start_base = s if strand == "+" else e - 1
    return signed_offset_past(strand, rec["five_p"], start_base) < 0


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def write_attribution_sidecar(
    corrected_tsv: str,
    output_path: str,
    genes: GeneIndex,
    ref: Dict[str, Dict[str, int]],
    unit: str = "reads",
    control_labels: Optional[Sequence[str]] = None,
) -> int:
    """Write the per-read sidecar for one library.  Returns rows written.

    The header comments are not decoration: the consumer joins this file blind, so
    the unit (reads vs UMI-collapsed molecules) and the provenance of the CPA
    reference have to travel WITH the data or a count silently changes meaning.
    """
    n = 0
    with open(output_path, "w", newline="") as fh:
        fh.write(f"# unit: {unit}\n")
        fh.write(f"# join_key: read_id (DRS) or umikey (UMI-collapsed cDNA)\n")
        fh.write(f"# termination_rule: max(annotated_end, observed_modal_CPA) "
                 f"+ {GENOME_TAIL_BP} bp (measured genome-wide p95)\n")
        fh.write(f"# escapes_gene uses a fixed {FIXED_WINDOW_BP} bp window and is the "
                 f"CONTRAST column; escapes_gene_cpa is the per-gene rule\n")
        if control_labels:
            fh.write(f"# cpa_reference_from: {','.join(control_labels)}\n")
        fh.write(f"# cpa_reference_genes: {len(ref)}\n")
        w = csv.DictWriter(fh, fieldnames=SIDECAR_HEADER, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for rec in _iter_corrected(corrected_tsv):
            w.writerow(classify_read(rec, genes, ref))
            n += 1
    logger.info("wrote %s (%d rows)", output_path, n)
    return n
