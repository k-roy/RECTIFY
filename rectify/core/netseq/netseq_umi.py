"""UMI-aware molecule counting for NET-seq (the RECTIFY NET-seq arm's standardized UMI path).

This is the NET-seq member of RECTIFY's per-chemistry UMI contract:

    | chemistry             | UMI length | UMI location    | dedup key                            |
    | ONT cDNA (PCB114.24)  | structured | adapter         | corrected 3' end x UMI (+ cluster)   |
    | Corall V2 TruSeq      | 12 nt      | read 1 5'       | fragment span x UMI (paired-end)     |
    | **NET-seq (here)**    | **6 nt**   | **read 5'**     | **corrected 3' end x strand x UMI**  |

WHY THE KEY IS NOT ``grouping.five_prime_unclipped()``  (the single most important point)
----------------------------------------------------------------------------------------
``five_prime_unclipped`` adds the 5' soft clip back to recover where a molecule *started*. That is right for
a chemistry whose 5' clip is genomic sequence the aligner declined to place. **NET-seq is not such a
chemistry:** the read 5' end IS the RNA 3' end, so the 5' clip holds the randomer UMI and the non-templated
tail -- material that is non-genomic BY DEFINITION. Keying on it would displace each read by its own tail
length, splitting one molecule across several positions in proportion to how long its tail was, i.e. biasing
the free-vs-tailed ratio precisely where the ratio is the measurement. The key here is therefore the
RECTIFY-**corrected** 3' end (post-walkback), and dedup runs AFTER correction, never before: two PCR copies
of one molecule can be soft-clipped differently (A-tract ambiguity, an error in the tail), so their raw ends
differ while their corrected ends agree.

ORIENTATION -- MEASURED, and it is the inverse of what a sense-read chemistry would give
---------------------------------------------------------------------------------------
In Churchman-style NET-seq the sequenced read is revcomp(RNA), so:

    gene strand  = '+' if read.is_reverse else '-'          (BAM strand is the INVERSE of the gene strand)
    RNA 3' end   = reference_end - 1 for a '+' gene, reference_start for a '-' gene   (the read's 5' terminus)

Verified three independent ways on GSE25107 (2011) and GSE159603 (2022) (planning ``478a`` CP1):
  * BAM strand == annotated gene strand in only 18-19% of reads over genes -- i.e. ~81% inverted;
  * the randomer-free 2011 library's read-5' soft clip is 85.6% T (first base 84.5% T) -- an RNA poly(A)
    tail seen through the reverse complement -- while its read-3' clip is compositionally flat;
  * the read-3' A-richness in 2022 is the linker's constant A (cutadapt: "bases preceding removed adapters:
    A 99.8%"), not a tail.
This also reproduces ``408c_clipdump.py``'s documented convention and the 697/697 intron self-test.

THE RANDOMER GEOMETRY -- three cases, and only blind length-slicing handles all three
------------------------------------------------------------------------------------
With a randomer of length ``N`` at the read 5' terminus and an observed 5'-terminal soft clip of length ``L``:

    L >  N   randomer + (L-N) nt of non-templated tail; the alignment terminus is the true RNA 3' end.
    L == N   randomer only, tail 0;                     the alignment terminus is the true RNA 3' end.
    L <  N   **(N-L) randomer bases were aligned to the genome by chance.** The tail is 0, and the alignment
             terminus OVERSHOOTS the true RNA 3' end by (N-L) bases -- corrected by
             :func:`randomer_overshoot`. The UMI is still ``seq[:N]``: clipping is an alignment fact, not a
             library fact. (Measured at ~26% of reads in GSE159603, so ignoring this case is not an edge case.)

WHY CLUSTERING DEFAULTS TO EXACT (edit distance 0)
--------------------------------------------------
A 6-nt UMI spans U = 4^6 = 4096 barcodes. At a deep NET-seq 3'-end position the barcode space fills:
k/U = 0.55 was measured at a single chrI position carrying 20k reads. In that regime the distance-1 graph is
dense (every 6-mer has 18 neighbours at distance 1), so ``umi_tools``-style directional clustering collapses
giant components and DESTROYS molecules instead of correcting errors. Directional clustering is calibrated
for 10-12-nt UMIs and does not transfer. ``edit_distance > 0`` is exposed so the effect can be MEASURED, not
because it is advisable here.

SATURATION IS A FIRST-CLASS OUTPUT, NOT A FOOTNOTE
--------------------------------------------------
Distinct UMIs observed at a position with n molecules: E[k] = U(1 - e^(-n/U)); the ML inverse is
``m_hat = ln(1 - k/U) / ln(1 - 1/U)`` (:func:`occupancy_corrected_molecules`). Naive molecule counts are
therefore compressed exactly at the highest-signal positions -- which are the biologically interesting ones.
Callers must (1) report reads and molecules side by side, (2) apply the occupancy correction, and
(3) fall back to read counts wherever ``k/U > SATURATION_KU_THRESHOLD`` (:func:`is_saturated`).
"""
from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Dict, List, NamedTuple, Optional, Tuple

from ..umi.clustering import canonical_umi, umi_components, umi_components_directional
from ..umi.dedup import UmiDedupStats
from ..umi.extract import extract_umi_from_read_id, is_valid_umi
from ..umi.grouping import is_spliced

#: Where the NET-seq UMI can be found.
NETSEQ_UMI_SOURCES = ("read5p", "name")

#: The Churchman/Couvillion NET-seq randomer length (measured, GSE159603).
NETSEQ_UMI_LENGTH = 6

#: Above this occupancy the barcode space is exhausted and molecule counts are not usable.
SATURATION_KU_THRESHOLD = 0.5

_CIGAR_SOFT_CLIP = 4
_CIGAR_REF_SKIP = 3


class NetseqBucketKey(NamedTuple):
    """The molecule-identifying key, excluding the UMI.

    Deliberately NOT ``core.umi.grouping.FragmentKey``: ``span`` is meaningless for a single-end library
    (there is no TLEN) and ``five_prime`` is the wrong anchor here (see the module docstring). ``spliced``
    is kept -- a NET-seq read crossing a junction is a different molecular species at the same 3' end.
    """
    contig: str
    strand: str          # gene/RNA strand, NOT the BAM strand
    corrected_3p: int    # post-walkback RNA 3' end, 0-based
    spliced: bool


class NetseqFragment(NamedTuple):
    """One read's dedup-relevant fields. BAM-free so the selection core is unit-testable."""
    qname: str
    contig: str
    strand: str
    corrected_3p: int
    spliced: bool
    umi: str
    score: int


def netseq_gene_strand(read: Any) -> str:
    """The gene/RNA strand of a NET-seq read -- the INVERSE of the BAM strand.

    See the module docstring: the read is revcomp(RNA), so a read from a '+' gene aligns reverse.
    """
    return "+" if read.is_reverse else "-"


def five_prime_terminal_clip(read: Any) -> int:
    """Soft clip at the READ's 5' terminus (== the RNA 3' terminus in this chemistry).

    ``cigartuples`` is in reference orientation, so for a reverse-mapped read the read's 5' terminus is the
    LAST cigar operation.
    """
    ct = read.cigartuples
    if not ct:
        return 0
    op, length = ct[-1] if read.is_reverse else ct[0]
    return length if op == _CIGAR_SOFT_CLIP else 0


def randomer_overshoot(clip_length: int, umi_length: int) -> int:
    """How many bases the alignment terminus overshoots the true RNA 3' end.

    Non-zero only in the ``L < N`` case, where part of the randomer aligned to the genome by chance.
    Returns a magnitude; the caller applies it in the RNA 3'-ward direction (see
    :func:`netseq_rna_three_prime`).
    """
    return max(0, umi_length - clip_length)


def netseq_rna_three_prime(read: Any, umi_length: int = 0) -> Tuple[int, str]:
    """``(rna_3p_position, gene_strand)`` for a NET-seq read, before walkback.

    Corrects the ``L < umi_length`` randomer-overshoot case. ``umi_length=0`` (a randomer-free library)
    disables that correction, which is then a no-op by construction.
    """
    strand = netseq_gene_strand(read)
    pos = (read.reference_end - 1) if strand == "+" else read.reference_start
    if umi_length:
        over = randomer_overshoot(five_prime_terminal_clip(read), umi_length)
        if over:
            # The RNA 3'-ward direction is towards higher coordinates on '+', lower on '-'.
            pos = pos - over if strand == "+" else pos + over
    return pos, strand


def netseq_nontemplated_tail(read: Any, umi_length: int = 0) -> str:
    """The non-templated 3'-terminal bases, in READ orientation (so an RNA poly(A) tail reads as ``T*``).

    Empty when the randomer partially aligned (``L < umi_length``): a molecule with a tail could not have
    had its randomer aligned, so that case carries no tail by construction.
    """
    clip = five_prime_terminal_clip(read)
    if clip <= umi_length:
        return ""
    seq = read.get_forward_sequence()
    if not seq:
        return ""
    return seq[umi_length:clip]


def extract_netseq_umi(
    read: Any,
    umi_length: int,
    umi_source: str = "read5p",
    separator: str = "_",
) -> Optional[str]:
    """Pull the UMI off a NET-seq read.

    ``read5p`` slices the first ``umi_length`` bases of the read AS SEQUENCED
    (``get_forward_sequence``) -- correct regardless of how much of the randomer the aligner happened to
    place. ``name`` parses a ``umi_tools extract``-style suffix, for libraries whose randomer was moved into
    the read name before alignment (the only option once the randomer has been trimmed out of SEQ).

    Returns None when no valid UMI is available; callers count that, they do not guess.
    """
    if umi_source not in NETSEQ_UMI_SOURCES:
        raise ValueError(f"unknown netseq umi source {umi_source!r}; expected one of {NETSEQ_UMI_SOURCES}")
    if umi_length <= 0:
        raise ValueError(f"umi_length must be positive, got {umi_length}")
    if umi_source == "name":
        return extract_umi_from_read_id(read.query_name, expected_length=umi_length, separator=separator)
    seq = read.get_forward_sequence()
    if not seq or len(seq) < umi_length:
        return None
    umi = seq[:umi_length]
    if not is_valid_umi(umi, expected_length=umi_length, allow_n=False):
        return None
    return umi


def occupancy_corrected_molecules(k: int, umi_length: int) -> float:
    """ML estimate of the true molecule count from ``k`` distinct UMIs seen at one position.

    Inverts E[k] = U(1 - e^(-m/U)) exactly: ``m_hat = ln(1 - k/U) / ln(1 - 1/U)``. Returns ``inf`` when the
    barcode space is completely exhausted (``k == U``), where the data carry no information about ``m``.
    """
    u = 4 ** umi_length
    if k <= 0:
        return 0.0
    if k >= u:
        return float("inf")
    return math.log(1.0 - k / u) / math.log(1.0 - 1.0 / u)


def is_saturated(k: int, umi_length: int, threshold: float = SATURATION_KU_THRESHOLD) -> bool:
    """True when occupancy is high enough that molecule counting must fall back to read counting."""
    return (k / (4 ** umi_length)) > threshold


def _neg_qname(qname: str) -> Tuple[int, ...]:
    return tuple(-ord(c) for c in qname)


def read_score(read: Any) -> int:
    """Representative-selection score: mapping quality + mean base quality (same rule as ``core/umi/dedup``)."""
    mq = read.mapping_quality or 0
    quals = read.query_qualities
    mean_bq = int(sum(quals) / len(quals)) if quals else 0
    return int(mq) + mean_bq


def select_netseq_molecules(
    fragments: List[NetseqFragment],
    edit_distance: int = 0,
    clustering: str = "exact",
    stats: Optional[UmiDedupStats] = None,
) -> Tuple[Dict[str, str], Dict[str, int]]:
    """Cluster NET-seq reads into molecules.

    Mirrors ``core.umi.dedup.select_molecules`` -- same representative rule, same ``UmiDedupStats``
    population -- differing ONLY in the bucket key (see :class:`NetseqBucketKey`). Returns
    ``(keepers, family_size)``: ``keepers`` maps the kept QNAME (one per molecule, the highest-score read) to
    a stable molecule id ``<contig>:<strand>:<corrected_3p>:<n>``, and ``family_size`` maps that id to its
    member count.

    ``clustering='exact'`` (the default, and the right one for a 6-nt UMI) groups by identical UMI. The
    graph clusterers are reachable for measurement via ``clustering='directional'|'components'``.
    """
    buckets: Dict[NetseqBucketKey, List[int]] = defaultdict(list)
    for i, fr in enumerate(fragments):
        buckets[NetseqBucketKey(fr.contig, fr.strand, fr.corrected_3p, fr.spliced)].append(i)

    keepers: Dict[str, str] = {}
    family_size: Dict[str, int] = {}
    for bkey, idxs in buckets.items():
        umis = [fragments[i].umi for i in idxs]
        if clustering == "exact" or edit_distance <= 0:
            groups: Dict[str, List[int]] = defaultdict(list)
            for local, u in enumerate(umis):
                groups[u].append(local)
            clusters = list(groups.values())
        elif clustering == "directional":
            clusters = umi_components_directional(umis, edit_distance)
        elif clustering == "components":
            clusters = umi_components(umis, edit_distance)
        else:
            raise ValueError(f"unknown clustering {clustering!r}")

        for c_no, member_local in enumerate(clusters):
            member_frags = [fragments[idxs[m]] for m in member_local]
            canon = canonical_umi([f.umi for f in member_frags])
            rep = max(member_frags, key=lambda f: (f.score, f.umi == canon, _neg_qname(f.qname)))
            mol_id = f"{bkey.contig}:{bkey.strand}:{bkey.corrected_3p}:{c_no}"
            keepers[rep.qname] = mol_id
            family_size[mol_id] = len(member_frags)
            if stats is not None:
                stats.n_molecules += 1
                stats.family_size_hist[len(member_frags)] += 1
                for f in member_frags:
                    d = (sum(1 for x, y in zip(f.umi, canon) if x != y)
                         if len(f.umi) == len(canon) else max(len(f.umi), len(canon)))
                    stats.within_family_edit_hist[d] += 1
    if stats is not None:
        stats.n_input_fragments = len(fragments)
    return keepers, family_size


class PositionMoleculeCount(NamedTuple):
    """Per-position molecule accounting, with the saturation escape hatch attached to the datum."""
    contig: str
    strand: str
    position: int
    reads: int
    distinct_umis: int          # k
    molecules_corrected: float  # m_hat
    saturated: bool


def summarize_positions(
    fragments: List[NetseqFragment],
    umi_length: int,
    threshold: float = SATURATION_KU_THRESHOLD,
) -> List[PositionMoleculeCount]:
    """Reads, distinct UMIs, occupancy-corrected molecules and the saturation flag, per 3'-end position.

    This is the object callers should emit alongside a read-level track -- never a molecule count alone.
    """
    per_pos: Dict[Tuple[str, str, int], Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for fr in fragments:
        per_pos[(fr.contig, fr.strand, fr.corrected_3p)][fr.umi] += 1
    out: List[PositionMoleculeCount] = []
    for (contig, strand, pos), umis in per_pos.items():
        k = len(umis)
        out.append(PositionMoleculeCount(
            contig=contig, strand=strand, position=pos,
            reads=sum(umis.values()), distinct_umis=k,
            molecules_corrected=occupancy_corrected_molecules(k, umi_length),
            saturated=is_saturated(k, umi_length, threshold),
        ))
    out.sort(key=lambda r: (r.contig, r.strand, r.position))
    return out


def stream_netseq_positions(
    fragments,
    umi_length: int,
    edit_distance: int = 0,
    clustering: str = "exact",
    stats: Optional[UmiDedupStats] = None,
    flush_pad: int = 2000,
    threshold: float = SATURATION_KU_THRESHOLD,
):
    """Streaming counterpart of :func:`summarize_positions` + :func:`select_netseq_molecules`.

    ⚠️ **This, not the list-based pair, is what a real library must go through.** A NET-seq library is
    ~50M reads over ~4M distinct 3'-end positions; materialising one ``NetseqFragment`` per read and one
    bucket per position costs tens of GB *on top of* the ~23 GB the read-level NET-seq pass already uses.
    The list-based functions stay for unit tests and small inputs, where being pure is worth more.

    Yields :class:`PositionMoleculeCount` in streaming order and populates ``stats`` as it goes, holding
    only a coordinate window in memory: a bucket is finalised once the input cursor has moved more than
    ``flush_pad`` past it. That is safe because the caller feeds fragments in BAM coordinate order and a
    read's corrected 3' end is within one read length of its ``reference_start`` -- ``flush_pad`` is an
    order of magnitude larger than that jitter. A contig change flushes everything unconditionally.
    """
    buckets: Dict[NetseqBucketKey, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    cur_contig: Optional[str] = None
    since_check = 0

    def _finalise(key: NetseqBucketKey) -> PositionMoleculeCount:
        umis = buckets.pop(key)
        reads = sum(umis.values())
        if clustering == "exact" or edit_distance <= 0:
            groups = [[u] for u in umis]
        elif clustering == "directional":
            keys = list(umis)
            groups = [[keys[i] for i in comp] for comp in umi_components_directional(keys, edit_distance)]
        elif clustering == "components":
            keys = list(umis)
            groups = [[keys[i] for i in comp] for comp in umi_components(keys, edit_distance)]
        else:
            raise ValueError(f"unknown clustering {clustering!r}")
        if stats is not None:
            stats.n_input_fragments += reads
            for members in groups:
                stats.n_molecules += 1
                stats.family_size_hist[sum(umis[u] for u in members)] += 1
        k = len(groups)
        return PositionMoleculeCount(
            contig=key.contig, strand=key.strand, position=key.corrected_3p, reads=reads,
            distinct_umis=k, molecules_corrected=occupancy_corrected_molecules(k, umi_length),
            saturated=is_saturated(k, umi_length, threshold),
        )

    for fr in fragments:
        if fr.contig != cur_contig:
            for key in list(buckets):
                yield _finalise(key)
            cur_contig = fr.contig
        buckets[NetseqBucketKey(fr.contig, fr.strand, fr.corrected_3p, fr.spliced)][fr.umi] += 1
        since_check += 1
        if since_check >= 100000:
            since_check = 0
            cutoff = fr.corrected_3p - flush_pad
            for key in [x for x in buckets if x.corrected_3p < cutoff]:
                yield _finalise(key)
    for key in list(buckets):
        yield _finalise(key)


def iter_netseq_fragments(
    bam_path: str,
    umi_length: int,
    umi_source: str = "read5p",
    exclusion_detector: Any = None,
    min_mapq: int = 0,
    rna3p_at: str = "read5p",
    max_reads: Optional[int] = None,
    stats: Optional[UmiDedupStats] = None,
):
    """Yield one :class:`NetseqFragment` per usable NET-seq read.

    Filtering (unmapped / secondary / supplementary / MAPQ / exclusion regions), chromosome
    standardisation and the corrected 3'-end position all come from
    ``core.bam.netseq_bam_processor``, so the molecule track and the read track are guaranteed to sit on the
    SAME positions -- a molecule count anchored on a different position than its read count would be
    unusable side by side, which is exactly how it must be reported.

    The randomer-overshoot correction (``L < umi_length``) is applied on top of the processor's position.
    """
    import pysam

    from ..bam.netseq_bam_processor import (
        get_netseq_3prime_position,
        standardize_netseq_chrom,
    )

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        chrom_map = {sq["SN"]: standardize_netseq_chrom(sq["SN"]) for sq in bam.header["SQ"]}
        n = 0
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            pos, strand, _ = get_netseq_3prime_position(read, rna3p_at=rna3p_at)
            over = randomer_overshoot(five_prime_terminal_clip(read), umi_length)
            if over:
                pos = pos - over if strand == "+" else pos + over
            chrom = chrom_map.get(read.reference_name, read.reference_name)
            if exclusion_detector is not None and exclusion_detector.is_excluded(chrom, pos):
                continue
            umi = extract_netseq_umi(read, umi_length, umi_source)
            if umi is None:
                if stats is not None:
                    stats.n_no_umi += 1
                continue
            yield NetseqFragment(
                qname=read.query_name, contig=chrom, strand=strand, corrected_3p=pos,
                spliced=is_spliced(read), umi=umi, score=read_score(read),
            )
            n += 1
            if max_reads and n >= max_reads:
                break
    finally:
        bam.close()


def netseq_fragment_from_read(
    read: Any,
    contig: str,
    umi_length: int,
    umi_source: str = "read5p",
    corrected_3p: Optional[int] = None,
) -> Optional[NetseqFragment]:
    """Build a :class:`NetseqFragment`, or None when the read carries no usable UMI.

    ``corrected_3p`` lets the caller pass the position produced by the full RECTIFY walkback; when omitted
    the pre-walkback RNA 3' end is used (correct for the randomer-overshoot case, but without the A-tract
    walkback).
    """
    umi = extract_netseq_umi(read, umi_length, umi_source)
    if umi is None:
        return None
    pos, strand = netseq_rna_three_prime(read, umi_length)
    return NetseqFragment(
        qname=read.query_name,
        contig=contig,
        strand=strand,
        corrected_3p=pos if corrected_3p is None else corrected_3p,
        spliced=is_spliced(read),
        umi=umi,
        score=read_score(read),
    )
