"""UMI-aware deduplication of paired-end short-read BAMs.

Collapses PCR/optical duplicates -- reads originating from the SAME molecule --
into one record per molecule, using (fragment span x UMI) as the molecule
identity and umi_tools-style directional clustering to absorb sequencing errors
in the UMI. The output is a clean BAM that every downstream analysis
(junction aggregation, coverage, ...) then reads as MOLECULE counts for free.

WHY THIS MATTERS FOR SPLICING (the whole point)
-----------------------------------------------
Novel/unannotated junctions are low-count calls, and low-count calls are exactly
what PCR jackpotting inflates: a real novel junction is supported by MANY DISTINCT
molecules, while an amplification artifact is ONE molecule copied many times.
Raw read count cannot tell them apart; molecule count can. So the headline metric
this feeds is deduped-molecules-per-junction, and the diagnostic is raw >> deduped.

APPROACH (umi_tools-compatible, memory-light, single pass over a sorted BAM)
---------------------------------------------------------------------------
The molecule key is taken from the R1 read alone, so we never need the mate in
hand (the classic umi_tools ``--paired`` trick):

    (contig, strand, R1 5'-unclipped, |TLEN|, spliced) x UMI

  * ``R1 5'-unclipped`` -- soft-clip-corrected 5' position (see grouping.py);
    two copies of one molecule share it even if one was soft-clipped.
  * ``|TLEN|`` -- template length = the full fragment span (outer coordinates),
    so R2's end is captured without reading R2. PCR duplicates share it by
    construction; different fragment lengths at the same 5' separate.
  * ``spliced`` -- keeps a spliced and an unspliced read at one position distinct
    (umi_tools ``--spliced-is-unique``); different CIGARs also give different
    |TLEN|, so this is belt-and-braces.

Within each such bucket the distinct UMIs are clustered directionally; each
cluster is one molecule. The representative (dedup) is the highest-quality
fragment; both its mates are written. The keep/drop decision is propagated to R2
by QNAME in a second pass, so mates never separate.

Consensus collapse (``collapse='consensus'``) is a documented extension point;
it is gated on the empirically-measured family-size distribution (most Illumina
UMI families are 1-2 reads, where consensus adds nothing) and is not implemented
in this first cut -- dedup is the dominant win. The stats this module emits ARE
that family-size measurement.
"""
from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, NamedTuple, Optional, Tuple

from .clustering import canonical_umi, umi_components, umi_components_directional
from .grouping import five_prime_unclipped, is_spliced

logger = logging.getLogger(__name__)

_DEFAULT_UMI_TAG = "RX"
_MOLECULE_TAG = "MI"       # SAM-standard molecular identifier (fgbio convention)
_FAMILY_SIZE_TAG = "cD"    # reads contributing to this molecule (fgbio convention)


class FragmentRecord(NamedTuple):
    """One fragment's dedup-relevant fields, extracted from its R1 read.

    BAM-free so the clustering core is unit-testable with literal values.
    """
    qname: str
    contig: str
    strand: str
    five_prime: int
    span: int          # |TLEN| (0 when unavailable, e.g. mate unmapped)
    spliced: bool
    umi: str
    score: int         # higher = better representative (e.g. summed mapping+base Q)


class BucketKey(NamedTuple):
    contig: str
    strand: str
    five_prime: int
    span: int
    spliced: bool


@dataclass
class UmiDedupStats:
    """Everything the empirical calibration (design doc S5) needs, produced as a
    side effect of a real dedup run rather than a separate analysis."""
    n_input_fragments: int = 0
    n_molecules: int = 0
    n_no_umi: int = 0
    # family_size -> number of molecules with that many reads
    family_size_hist: Dict[int, int] = field(default_factory=lambda: defaultdict(int))
    # within-family UMI Hamming-distance -> count of (member, canonical) pairs
    within_family_edit_hist: Dict[int, int] = field(default_factory=lambda: defaultdict(int))

    @property
    def duplication_rate(self) -> float:
        """Fraction of input fragments that were duplicates (1 - molecules/reads)."""
        if self.n_input_fragments == 0:
            return 0.0
        return 1.0 - self.n_molecules / self.n_input_fragments

    def as_dict(self) -> dict:
        return {
            "n_input_fragments": self.n_input_fragments,
            "n_molecules": self.n_molecules,
            "n_no_umi": self.n_no_umi,
            "duplication_rate": round(self.duplication_rate, 6),
            "family_size_hist": dict(sorted(self.family_size_hist.items())),
            "within_family_edit_hist": dict(sorted(self.within_family_edit_hist.items())),
        }


def _hamming(a: str, b: str) -> int:
    """Hamming distance for equal-length strings; length diff counts as large."""
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for x, y in zip(a, b) if x != y)


def select_molecules(
    fragments: List[FragmentRecord],
    edit_distance: int = 1,
    clustering: str = "directional",
    stats: Optional[UmiDedupStats] = None,
) -> Tuple[Dict[str, str], Dict[str, int]]:
    """Cluster fragments into molecules.

    Returns ``(keepers, family_size)`` where ``keepers`` maps each kept QNAME (one
    per molecule -- the highest-score fragment) to a stable ``molecule_id``
    (``<contig>:<strand>:<5'>:<span>:<n>``, suitable for the MI tag), and
    ``family_size`` maps ``molecule_id`` to the number of member fragments (the cD
    tag). Pure function over :class:`FragmentRecord` -- no BAM.

    ``stats`` (when provided) is populated with the family-size and within-family
    edit-distance distributions -- the empirical calibration inputs (design S5).
    """
    cluster_fn = (umi_components_directional if clustering == "directional"
                  else umi_components)

    # Bucket fragments by (position span), then cluster the UMIs inside each bucket.
    buckets: Dict[BucketKey, List[int]] = defaultdict(list)
    for i, fr in enumerate(fragments):
        buckets[BucketKey(fr.contig, fr.strand, fr.five_prime, fr.span, fr.spliced)].append(i)

    keepers: Dict[str, str] = {}
    family_size: Dict[str, int] = {}
    for bkey, idxs in buckets.items():
        umis = [fragments[i].umi for i in idxs]
        clusters = cluster_fn(umis, edit_distance)
        for c_no, member_local in enumerate(clusters):
            member_frags = [fragments[idxs[m]] for m in member_local]
            # Canonical UMI = the cluster's most-frequent UMI (its inferred true
            # sequence); error-variant UMIs are rarer. Used both to prefer the
            # representative from the dominant UMI and to measure within-family error.
            canon = canonical_umi([f.umi for f in member_frags])
            # Representative = highest score, then prefer a read carrying the
            # canonical UMI (not a rare error variant), then smallest qname.
            rep = max(member_frags,
                      key=lambda f: (f.score, f.umi == canon, _neg_qname(f.qname)))
            mol_id = f"{bkey.contig}:{bkey.strand}:{bkey.five_prime}:{bkey.span}:{c_no}"
            keepers[rep.qname] = mol_id
            family_size[mol_id] = len(member_frags)
            if stats is not None:
                stats.n_molecules += 1
                stats.family_size_hist[len(member_frags)] += 1
                for f in member_frags:
                    stats.within_family_edit_hist[_hamming(f.umi, canon)] += 1
    return keepers, family_size


def _neg_qname(qname: str) -> Tuple[int, ...]:
    """Sort key that makes 'lexicographically smallest' win under max()."""
    return tuple(-ord(c) for c in qname)


def _read_score(read) -> int:
    """Representative-selection score: mapping quality + mean base quality.

    Higher-quality reads make better representatives (dedup mode) -- this is the
    cheap analogue of consensus for the common shallow-family case.
    """
    mq = read.mapping_quality or 0
    quals = read.query_qualities
    mean_bq = int(sum(quals) / len(quals)) if quals else 0
    return int(mq) + mean_bq


def extract_fragments(
    bam,
    umi_tag: str = _DEFAULT_UMI_TAG,
    r1_sense: bool = True,
    stats: Optional[UmiDedupStats] = None,
) -> List[FragmentRecord]:
    """Pull one :class:`FragmentRecord` per usable R1 (or unpaired) primary read.

    ``bam`` is an open ``pysam.AlignmentFile``. Reads without the UMI tag are
    counted in ``stats.n_no_umi`` and skipped (they cannot be deduped safely).
    """
    frags: List[FragmentRecord] = []
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        # One record per fragment: use R1 for pairs, the read itself if unpaired.
        if read.is_paired and not read.is_read1:
            continue
        try:
            umi = read.get_tag(umi_tag)
        except KeyError:
            if stats is not None:
                stats.n_no_umi += 1
            continue
        if not umi:
            if stats is not None:
                stats.n_no_umi += 1
            continue
        fp = five_prime_unclipped(read)
        if fp is None:
            continue
        # Gene/RNA strand: R1 orientation, inverted for antisense chemistries.
        read_is_forward = not read.is_reverse
        sense = read_is_forward if r1_sense else (not read_is_forward)
        strand = "+" if sense else "-"
        span = abs(read.template_length or 0)
        frags.append(FragmentRecord(
            qname=read.query_name,
            contig=read.reference_name,
            strand=strand,
            five_prime=fp,
            span=span,
            spliced=is_spliced(read),
            umi=umi,
            score=_read_score(read),
        ))
    if stats is not None:
        stats.n_input_fragments = len(frags)
    return frags


def dedup_bam(
    input_bam: str,
    output_bam: str,
    umi_tag: str = _DEFAULT_UMI_TAG,
    edit_distance: int = 1,
    clustering: str = "directional",
    r1_sense: bool = True,
    collapse: str = "dedup",
) -> UmiDedupStats:
    """Write a UMI-deduplicated copy of ``input_bam`` and return run statistics.

    ``input_bam`` must be coordinate-sorted and carry the UMI on ``umi_tag``
    (``RX`` by default) on at least the R1 reads. One molecule -> its
    highest-quality fragment survives; both mates of the survivor are written and
    tagged with ``MI`` (molecule id) and ``cD`` (family size).
    """
    import pysam

    if collapse not in ("dedup", "consensus"):
        raise ValueError(f"unknown collapse mode {collapse!r}")
    if collapse == "consensus":
        raise NotImplementedError(
            "collapse='consensus' is a planned extension (fgbio-style read-space "
            "posterior); it is gated on the measured family-size distribution "
            "(design doc S5.1). Use collapse='dedup' for now, and inspect the "
            "returned UmiDedupStats.family_size_hist to decide if consensus is "
            "worth it (most Illumina UMI families are 1-2 reads)."
        )

    stats = UmiDedupStats()

    # Pass 1: extract fragments (R1 only), cluster into molecules, pick keepers.
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        fragments = extract_fragments(bam_in, umi_tag=umi_tag, r1_sense=r1_sense, stats=stats)
    keepers, family_size = select_molecules(
        fragments, edit_distance=edit_distance, clustering=clustering, stats=stats
    )

    # Pass 2: stream input, write only the kept molecules' reads (both mates),
    # tagging MI (molecule id) + cD (family size).
    n_written = 0
    with pysam.AlignmentFile(input_bam, "rb") as bam_in, \
         pysam.AlignmentFile(output_bam, "wb", header=bam_in.header) as bam_out:
        for read in bam_in:
            mid = keepers.get(read.query_name)
            if mid is None:
                continue
            read.set_tag(_MOLECULE_TAG, mid, value_type="Z")
            cd = family_size.get(mid)
            if cd is not None:
                read.set_tag(_FAMILY_SIZE_TAG, int(cd), value_type="i")
            bam_out.write(read)
            n_written += 1

    logger.info(
        "UMI dedup %s: %d fragments -> %d molecules (%.1f%% duplicates); "
        "%d reads written; %d reads lacked a %s tag",
        input_bam, stats.n_input_fragments, stats.n_molecules,
        100 * stats.duplication_rate, n_written, stats.n_no_umi, umi_tag,
    )
    return stats
