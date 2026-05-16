"""Per-read splice classification and per-gene aggregation.

Classification (**strict 0-bp matching** on intron coordinates, after optional
boundary snapping):

- ``unspliced``: read spans an annotated intron but has no N-cigar in that
  range — intron retention.
- ``annotated_spliced``: every N-cigar in the spanned region exactly matches
  an annotated intron coordinate pair.
- ``alternative_spliced``: at least one N-cigar's donor (5'SS) OR acceptor
  (3'SS) matches an annotated splice site, but the other end does not.
  Captures alt-5'SS / alt-3'SS / exon-skipping events.
- ``novel_spliced``: at least one N-cigar has *neither* end matching any
  annotated splice site of the host gene.
- ``no_intron_span``: read does not fully span any annotated intron (e.g.,
  truncated read or read aligned to a different transcript boundary). Useful
  to filter out from per-gene fractions.

The "strict 0-bp" rule is applied *after* boundary snapping: an aligner that
places an N-cigar 1–2 bp off the true donor/acceptor due to basecaller error
in a homopolymer near the splice site would otherwise be miscategorised as
"novel". The :mod:`rectify.core.splice.junction_refiner` module handles this in the
correction pass (commit ``de27543`` and earlier). For the audit script below
we use a lightweight inline snapping based on canonical signal proximity —
this is the same conceptual model as ``junction_refiner`` but without the
full HP-aware DP scoring.

See ``.github/ISSUE_walkback_integration.md`` and the run-all audit notes.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

import pysam


# CIGAR op code for N (intron skip) is 3.
_CIGAR_N = 3
_CIGAR_M = 0
_CIGAR_D = 2
_CIGAR_EQ = 7
_CIGAR_X = 8
# These three advance reference position.
_CIGAR_REFCONSUMERS = {_CIGAR_M, _CIGAR_D, _CIGAR_N, _CIGAR_EQ, _CIGAR_X}


# Outcome strings — public constants so callers can branch on them.
CLASS_UNSPLICED = "unspliced"
CLASS_ANNOTATED = "annotated_spliced"
CLASS_ALTERNATIVE = "alternative_spliced"
CLASS_NOVEL = "novel_spliced"
CLASS_NO_SPAN = "no_intron_span"


@dataclass(frozen=True)
class GeneIntronSet:
    """Annotated introns for a single gene, plus precomputed splice-site sets.

    All coordinates are 0-based half-open in reference coordinates: an intron
    [start, end) means the splice donor lies at the exon position immediately
    before ``start`` and the splice acceptor lies at the exon position
    starting at ``end``. For BAM N-cigars at reference position ``R`` with
    length ``L``, the N-op represents an intron at ``[R, R+L)``, so the
    intron tuple is ``(R, R+L)`` — directly comparable to the annotated set.
    """

    gene_id: str
    chrom: str
    strand: str
    introns: FrozenSet[Tuple[int, int]]
    donors: FrozenSet[int] = field(init=False)
    acceptors: FrozenSet[int] = field(init=False)

    def __post_init__(self):
        # Donor and acceptor sets — strand-agnostic at the coord level.
        # For gene on + strand: donor at intron_start, acceptor at intron_end.
        # For gene on - strand: roles swap but the *coordinates* still come
        # from the (start, end) tuple; classification only cares about which
        # set of integers an N-op end appears in, so we don't care which
        # genomic end is "donor" vs "acceptor" — we just need both sides
        # tracked.
        object.__setattr__(
            self, "donors", frozenset(start for start, _ in self.introns)
        )
        object.__setattr__(
            self, "acceptors", frozenset(end for _, end in self.introns)
        )


def extract_n_ops(read: pysam.AlignedSegment) -> List[Tuple[int, int]]:
    """Return list of ``(intron_start, intron_end)`` for every N-cigar.

    Empty list if the read has no N-ops (i.e., no introns crossed).
    """
    if read.is_unmapped or read.cigartuples is None:
        return []
    out: List[Tuple[int, int]] = []
    ref_pos = read.reference_start
    for op, length in read.cigartuples:
        if op == _CIGAR_N:
            out.append((ref_pos, ref_pos + length))
        if op in _CIGAR_REFCONSUMERS:
            ref_pos += length
    return out


def read_spans_intron(
    read: pysam.AlignedSegment, intron: Tuple[int, int], margin: int = 5
) -> bool:
    """True iff the read's alignment fully covers ``intron`` plus ``margin``
    bp of flanking exon on each side. The margin prevents counting reads
    that barely overlap an intron edge.
    """
    if read.is_unmapped:
        return False
    start, end = intron
    return read.reference_start <= start - margin and read.reference_end >= end + margin


def snap_n_op_to_annotated(
    n_op: Tuple[int, int],
    introns: FrozenSet[Tuple[int, int]],
    tolerance: int,
) -> Tuple[int, int]:
    """Snap an N-op to the nearest annotated intron if both endpoints lie
    within ``tolerance`` bp.

    This is a deliberately conservative inline snapper — strictly tighter
    than :mod:`junction_refiner`. We only snap when *both* the donor and
    acceptor sides are within tolerance of an *annotated* intron; that
    prevents a basecaller HP slip near one splice site from being
    "rescued" into the nearest annotated junction when the other end is
    far away.

    Returns the snapped tuple if a match is found, else the original.
    """
    if tolerance <= 0:
        return n_op
    start, end = n_op
    best: Optional[Tuple[int, int]] = None
    best_d = tolerance + 1
    for a_start, a_end in introns:
        ds = abs(a_start - start)
        de = abs(a_end - end)
        if ds <= tolerance and de <= tolerance:
            d = ds + de
            if d < best_d:
                best_d = d
                best = (a_start, a_end)
    return best if best is not None else n_op


def classify_read(
    read: pysam.AlignedSegment,
    gene: GeneIntronSet,
    snap_tolerance: int = 0,
    span_margin: int = 5,
) -> str:
    """Classify one read against one gene's annotated introns.

    Parameters
    ----------
    snap_tolerance
        If ``> 0``, each N-op is snapped to the nearest annotated intron
        when both endpoints are within ``snap_tolerance`` bp (see
        :func:`snap_n_op_to_annotated`). After snapping, comparison is
        **strict 0-bp**. Default 0 means no snapping — pure strict match.
    span_margin
        Read must overlap each side of the intron by at least this many bp
        to count as "spanning" (used for the unspliced detection).
    """
    n_ops = extract_n_ops(read)
    introns = gene.introns

    # Restrict N-ops to those falling inside the gene's annotated region.
    # An N-op outside the gene's intron range is from a different splice
    # event (e.g., a different transcript) and shouldn't influence the
    # classification of *this* gene's splice pattern.
    if not introns:
        return CLASS_NO_SPAN
    gene_min = min(s for s, _ in introns)
    gene_max = max(e for _, e in introns)
    n_ops_in_gene = [
        op for op in n_ops
        if op[0] >= gene_min - 50 and op[1] <= gene_max + 50
    ]

    # If snap_tolerance > 0, snap each N-op.
    if snap_tolerance > 0:
        n_ops_in_gene = [
            snap_n_op_to_annotated(op, introns, snap_tolerance)
            for op in n_ops_in_gene
        ]

    spans_any = any(read_spans_intron(read, i, span_margin) for i in introns)

    if not spans_any:
        return CLASS_NO_SPAN

    if not n_ops_in_gene:
        # Read spans an annotated intron but has no N-cigar → intron retention.
        return CLASS_UNSPLICED

    # At least one N-op present. Classify.
    annotated_set = introns
    donors = gene.donors
    acceptors = gene.acceptors

    any_novel = False
    any_alt = False
    all_match = True
    for op in n_ops_in_gene:
        if op in annotated_set:
            continue
        all_match = False
        donor_hit = op[0] in donors
        acceptor_hit = op[1] in acceptors
        if donor_hit or acceptor_hit:
            any_alt = True
        else:
            any_novel = True

    if all_match:
        return CLASS_ANNOTATED
    if any_novel:
        return CLASS_NOVEL
    if any_alt:
        return CLASS_ALTERNATIVE
    return CLASS_ANNOTATED  # unreachable


@dataclass
class GeneSpliceStats:
    """Per-gene per-sample aggregate."""

    gene_id: str
    sample: str
    counts: Counter = field(default_factory=Counter)

    @property
    def total_spanning(self) -> int:
        """Reads that span at least one annotated intron — denominator for
        meaningful percentages."""
        return sum(
            v for k, v in self.counts.items() if k != CLASS_NO_SPAN
        )

    def pct(self, cls: str) -> float:
        denom = self.total_spanning
        if denom == 0:
            return float("nan")
        return 100.0 * self.counts.get(cls, 0) / denom

    def as_row(self) -> Dict[str, object]:
        return {
            "gene_id": self.gene_id,
            "sample": self.sample,
            "n_spanning": self.total_spanning,
            "n_unspliced": self.counts.get(CLASS_UNSPLICED, 0),
            "n_annotated": self.counts.get(CLASS_ANNOTATED, 0),
            "n_alternative": self.counts.get(CLASS_ALTERNATIVE, 0),
            "n_novel": self.counts.get(CLASS_NOVEL, 0),
            "n_no_span": self.counts.get(CLASS_NO_SPAN, 0),
            "pct_unspliced": round(self.pct(CLASS_UNSPLICED), 2),
            "pct_annotated": round(self.pct(CLASS_ANNOTATED), 2),
            "pct_alternative": round(self.pct(CLASS_ALTERNATIVE), 2),
            "pct_novel": round(self.pct(CLASS_NOVEL), 2),
        }


def classify_bam_for_genes(
    bam_path: str,
    genes: Sequence[GeneIntronSet],
    sample_name: str,
    snap_tolerance: int = 0,
    span_margin: int = 5,
    min_mapq: int = 0,
) -> List[GeneSpliceStats]:
    """Walk a BAM, fetching reads per gene region and classifying each.

    Returns one ``GeneSpliceStats`` per input gene.
    """
    stats = [GeneSpliceStats(g.gene_id, sample_name) for g in genes]
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for g, s in zip(genes, stats):
            try:
                fetch_iter = bam.fetch(g.chrom, min(d for d, _ in g.introns) - 1000, max(e for _, e in g.introns) + 1000)
            except ValueError:
                # Region not in BAM header — skip.
                continue
            for read in fetch_iter:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < min_mapq:
                    continue
                cls = classify_read(read, g, snap_tolerance=snap_tolerance, span_margin=span_margin)
                s.counts[cls] += 1
    return stats
