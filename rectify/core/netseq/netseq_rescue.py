"""
Donor-side junction rescue and non-templated 3'-tail calling for NET-seq reads.

Two corrections that a NET-seq 3'-end track needs and that nothing in RECTIFY performed before
(audited in planning ``834`` §7 / ``835`` §2.1 -- ``rectify correct``'s only junction-aware
clip re-placer, Module 2F ``rescue_3ss_truncation``, moves the transcript **5'** coordinate):

**1. Donor-side ("5'SS") junction rescue.**  Under Churchman NET-seq geometry the RNA 3' end is the
read's 5' terminus, so a nascent RNA whose 3' end sits only 1-10 nt into exon 2 gives an alignment
whose RNA-3'-most block stops at the exon-1 side of the intron: the exon-2 overhang is too short to
anchor across the junction, so a local aligner soft-clips it (STAR ``--alignEndsType Local``, bbmap)
or mis-extends the alignment a few nt into the intron.  Either way the read is called at the 5'
splice site and becomes a **false splicing intermediate** -- at exactly the coordinate where the real
splicing-intermediate signal lives.  :func:`rescue_read` matches the clip (plus any mis-extended
intronic bases) against the start of exon 2 of a pooled junction and moves the 3' end there.

**2. Non-templated tail calling.**  The same clip also carries the poly(A)/oligo(A) tail, and in
randomer libraries it carries a 6-nt randomer *distal* to the tail (planning ``478`` §1: the randomer
and the tail occupy the SAME terminus).  :func:`call_tail` strips the declared randomer from the
distal end, takes the A-run adjacent to the alignment, and adds the aligned A's recovered by a
genome-aware walkback (Chanfreau invariant 7) -- never the 0.8 A-fraction-over-the-whole-clip test
that :func:`~rectify.core.bam.netseq_bam_processor.detect_oligo_a_in_softclip` uses, which fails on
every randomer-bearing read.

Coordinates
-----------
Everything is expressed on the **gene strand**, which for ``--rna3p-at read5p`` is the INVERSE of the
BAM strand (:func:`~rectify.core.bam.netseq_bam_processor.netseq_gene_strand`).  For an intron at
0-based half-open ``[intron_start, intron_end)``:

===============  =========================  ==========================
term             ``+`` gene                 ``-`` gene
===============  =========================  ==========================
``donor``        ``intron_start``           ``intron_end - 1``
``exon1_last``   ``intron_start - 1``       ``intron_end``
``exon2_first``  ``intron_end``             ``intron_start - 1``
5'-ward in RNA   decreasing coordinate      increasing coordinate
===============  =========================  ==========================

``donor`` is the FIRST INTRONIC base and ``exon2_first`` the FIRST EXON-2 base, both on the gene
strand.  The external ``--junction-pool`` TSV uses the same convention (see :func:`load_junction_tsv`).

Author: Kevin R. Roy
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pysam

from ..unified_record import UnifiedReadRecord


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    """Reverse complement of an ACGTN string (case preserved by the translation table)."""
    return seq.translate(_COMPLEMENT)[::-1]


def longest_common_prefix(a: str, b: str) -> int:
    """Length of the longest common prefix of two strings (case-insensitive)."""
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


# ---------------------------------------------------------------------------------------------
# Junction pool
# ---------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class NetseqJunction:
    """One pooled intron, pre-resolved to gene-strand donor / exon-boundary coordinates."""

    chrom: str
    strand: str
    intron_start: int   # 0-based genomic, inclusive
    intron_end: int     # 0-based genomic, exclusive
    donor: int          # first intronic base, gene strand
    exon1_last: int     # last exon-1 base, gene strand
    exon2_first: int    # first exon-2 base, gene strand

    @property
    def intron_length(self) -> int:
        return self.intron_end - self.intron_start


def make_junction(chrom: str, intron_start: int, intron_end: int, strand: str) -> NetseqJunction:
    """Build a :class:`NetseqJunction` from a genomic intron interval and its GENE strand."""
    if strand == '+':
        donor, exon1_last, exon2_first = intron_start, intron_start - 1, intron_end
    else:
        donor, exon1_last, exon2_first = intron_end - 1, intron_end, intron_start - 1
    return NetseqJunction(
        chrom=chrom, strand=strand,
        intron_start=intron_start, intron_end=intron_end,
        donor=donor, exon1_last=exon1_last, exon2_first=exon2_first,
    )


class JunctionPool:
    """Donor-indexed pool of introns, queried by the RNA 3' end of a read.

    The index key is ``exon1_last`` (the last exon-1 base on the gene strand), because a read that
    needs rescuing has its RNA-3'-most aligned base AT that coordinate (clean soft clip) or up to
    ``max_intronic`` nt gene-strand-ward of it (aligner mis-extension into the intron).
    """

    def __init__(self, junctions: Iterable[NetseqJunction] = ()):
        self._by_exon1: Dict[Tuple[str, str], Dict[int, NetseqJunction]] = {}
        self.n_junctions = 0
        #: kept / dropped intron counts by PARENT FEATURE TYPE (mRNA, tRNA, snoRNA, ...)
        self.by_parent_type: Dict[str, int] = {}
        self.dropped_by_parent_type: Dict[str, int] = {}
        for j in junctions:
            self.add(j)

    def add(self, junction: NetseqJunction) -> None:
        """Add one junction. Ties on ``exon1_last`` keep the SHORTEST intron (the nearest acceptor)."""
        key = (junction.chrom, junction.strand)
        slot = self._by_exon1.setdefault(key, {})
        prev = slot.get(junction.exon1_last)
        if prev is None:
            slot[junction.exon1_last] = junction
            self.n_junctions += 1
        elif junction.intron_length < prev.intron_length:
            slot[junction.exon1_last] = junction

    def __len__(self) -> int:
        return self.n_junctions

    def lookup(
        self,
        chrom: str,
        strand: str,
        p: int,
        max_intronic: int,
    ) -> Optional[Tuple[int, NetseqJunction]]:
        """Nearest pooled donor for an RNA 3' end at ``p``.

        Returns ``(n_intronic, junction)`` where ``n_intronic`` is the number of aligned bases at or
        past the donor (0 == the aligned end sits exactly on ``exon1_last``), or ``None`` when no
        pooled junction has its exon-1 boundary within ``max_intronic`` nt gene-strand-upstream of
        ``p``.  The SMALLEST ``n_intronic`` wins, i.e. the donor nearest the aligned end.
        """
        slot = self._by_exon1.get((chrom, strand))
        if not slot:
            return None
        step = -1 if strand == '+' else 1  # gene-strand-ward from p back to exon1_last
        for n_intronic in range(0, max_intronic + 1):
            exon1_last = p + step * n_intronic
            junction = slot.get(exon1_last)
            if junction is not None:
                return n_intronic, junction
        return None

    # -- constructors -------------------------------------------------------------------------

    @classmethod
    def from_annotation(cls, annotation_path, include_trna: bool = False,
                        include_organellar: bool = False) -> "JunctionPool":
        """Build the pool from a GFF/GTF.

        Introns are read with :func:`parse_gff_introns`, which also resolves each intron's PARENT
        FEATURE TYPE, so tRNA introns can be excluded -- see ``include_trna``. When the annotation
        carries no ``intron`` feature at all (an exon-only GTF), it falls back to
        :func:`rectify.core.consensus.consensus.load_annotated_junctions`, which derives introns
        from ``exon`` records; parent types are unknown on that path and nothing is filtered.

        Args:
            annotation_path: GFF3/GTF, optionally gzipped.
            include_organellar: keep mitochondrial/organellar introns in the pool. **OFF by
                default** — they are group I/II SELF-SPLICING introns on a genome Pol II does not
                transcribe, and their parent type is ``mRNA`` so the tRNA filter misses them
                (measured on wt_rep3: 94 of 580 rescues, 16 %, were fabricated on chrMito).
            include_trna: keep tRNA introns in the pool. **OFF by default.** tRNA splicing is a
                different machinery (a tRNA endonuclease/ligase pathway, not the spliceosome) at
                Pol III loci that ``rectify netseq`` excludes from its tracks anyway, so a tRNA
                intron in the pool can only manufacture rescues -- measured on wt_rep3, the single
                tRNA locus ``YNCO0031W_tRNA`` collected 46 of them.

        The pool records ``by_parent_type`` (kept) and ``dropped_by_parent_type``, which
        ``rectify netseq`` prints so the composition is visible rather than assumed.
        """
        pool = cls()
        rows = parse_gff_introns(annotation_path)
        if rows:
            for chrom, start, end, strand, parent_type, parent_id in rows:
                if end <= start or strand not in ('+', '-'):
                    continue
                if not include_trna and is_trna_parent(parent_type, parent_id):
                    pool.dropped_by_parent_type[parent_type] = \
                        pool.dropped_by_parent_type.get(parent_type, 0) + 1
                    continue
                if not include_organellar and is_organellar_contig(chrom):
                    key = f"{parent_type} (organellar)"
                    pool.dropped_by_parent_type[key] = \
                        pool.dropped_by_parent_type.get(key, 0) + 1
                    continue
                pool.by_parent_type[parent_type] = pool.by_parent_type.get(parent_type, 0) + 1
                pool.add(make_junction(chrom, start, end, strand))
            return pool

        from ..consensus.consensus import load_annotated_junctions

        for chrom, start, end, strand in load_annotated_junctions(str(annotation_path)):
            if end <= start or strand not in ('+', '-'):
                continue
            if not include_organellar and is_organellar_contig(chrom):
                pool.dropped_by_parent_type['unknown (organellar)'] = \
                    pool.dropped_by_parent_type.get('unknown (organellar)', 0) + 1
                continue
            pool.by_parent_type['unknown'] = pool.by_parent_type.get('unknown', 0) + 1
            pool.add(make_junction(chrom, start, end, strand))
        return pool

    @classmethod
    def from_tsv(cls, tsv_path) -> "JunctionPool":
        """Build the pool from an external junction table (see :func:`load_junction_tsv`)."""
        pool = cls()
        for junction in load_junction_tsv(tsv_path):
            pool.add(junction)
        return pool


#: Parent feature types that mean "this intron is spliced by the tRNA endonuclease, not the
#: spliceosome". SGD writes the parent of a tRNA intron as a ``tRNA`` feature whose own parent is a
#: ``tRNA_gene``; the id also carries the ``_tRNA`` suffix, and both are checked because an
#: annotation that omits the ``tRNA`` feature line still names its parent that way.
_TRNA_PARENT_TYPES = frozenset({'trna', 'trna_gene'})


def is_organellar_contig(chrom: str) -> bool:
    """Whether ``chrom`` is a mitochondrial (or otherwise organellar) contig.

    Mitochondrial introns are **group I and group II self-splicing** introns — in S. cerevisiae,
    ~32 of them across COX1, COB and 21S rRNA. They are not spliceosomal, they carry no
    GT..AG-style donor/acceptor grammar this rescue models, and Pol II does not transcribe the
    mitochondrial genome at all, so a nascent-Pol-II 3' end cannot legitimately sit at one. Their
    parent feature type is ``mRNA``, exactly like a nuclear intron, so :func:`is_trna_parent`
    does not catch them — the CONTIG is the signal.

    Measured on wt_rep3 (planning 856, 2026-09-03): with mito introns in the pool, **94 of 580
    rescues — 16 %% — were on chrMito**, i.e. fabricated 3'-end re-placements on a genome the assay
    does not measure. Same argument, and the same measurement, that dropped tRNA introns.

    The patterns are shared with :class:`rectify.core.exclusion_regions.ExclusionRegionDetector`.
    """
    from ..exclusion_regions import ExclusionRegionDetector

    lowered = (chrom or '').lower()
    if any(pat.lower() in lowered for pat in ExclusionRegionDetector.MITO_PATTERNS):
        return True
    return chrom in ExclusionRegionDetector.MITO_NCBI


def is_trna_parent(parent_type: str, parent_id: str) -> bool:
    """Whether an intron's parent identifies it as a tRNA intron."""
    if (parent_type or '').lower() in _TRNA_PARENT_TYPES:
        return True
    return (parent_id or '').endswith('_tRNA')


def parse_gff_introns(annotation_path):
    """Read ``intron``-class features from a GFF3, resolving each one's PARENT FEATURE TYPE.

    Returns a list of ``(chrom, intron_start, intron_end, strand, parent_type, parent_id)`` with
    ``intron_start`` 0-based and ``intron_end`` exclusive, chrom standardized the same way
    :func:`rectify.core.consensus.consensus.load_annotated_junctions` standardizes it (so pool
    lookups match the read chrom names).

    The parent type is what makes a tRNA intron distinguishable from an mRNA one: SGD's GFF gives
    an intron ``Parent=YNCO0031W_tRNA``, and a separate line declares ``ID=YNCO0031W_tRNA`` with
    feature type ``tRNA``. Any subtype containing "intron" is collected
    (``intron``, ``five_prime_UTR_intron``, ...). Returns ``[]`` for an annotation with no intron
    features, which is the caller's signal to fall back to exon-derived junctions.
    """
    import gzip as _gzip

    from ...utils.genome import standardize_chrom_name

    path = str(annotation_path)
    _open = _gzip.open if path.endswith('.gz') else open

    id_type: Dict[str, str] = {}
    raw: List[Tuple[str, int, int, str, str]] = []
    with _open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            attrs = {}
            for kv in parts[8].split(';'):
                if '=' in kv:
                    key, value = kv.split('=', 1)
                    attrs[key.strip()] = value.strip()
            feature = parts[2]
            if 'ID' in attrs:
                id_type[attrs['ID']] = feature
            if 'intron' in feature.lower():
                parent = (attrs.get('Parent') or '').split(',')[0]
                strand = parts[6] if parts[6] in ('+', '-') else '+'
                raw.append((standardize_chrom_name(parts[0]),
                            int(parts[3]) - 1, int(parts[4]), strand, parent))

    return [(chrom, start, end, strand, id_type.get(parent, 'unknown'), parent)
            for chrom, start, end, strand, parent in raw]


def load_junction_tsv(tsv_path) -> List[NetseqJunction]:
    """Read an external junction table.

    Header-driven TSV with columns ``chrom``, ``donor``, ``acceptor``, ``strand``.  **Both sites are
    INTRONIC and on the gene strand**: ``donor`` = the FIRST intronic base (the base immediately
    3'-ward of the last exon-1 base), ``acceptor`` = the LAST intronic base (the base immediately
    5'-ward of the first exon-2 base), 0-based.  For a ``-`` gene ``donor > acceptor``.  This is the
    convention a long-read-derived pool should be written in; anything else shifts every rescued
    3' end by one nucleotide.
    """
    junctions: List[NetseqJunction] = []
    path = Path(tsv_path)
    with open(path) as fh:
        header: Optional[List[str]] = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if header is None:
                header = [f.strip().lower() for f in fields]
                missing = {'chrom', 'donor', 'acceptor', 'strand'} - set(header)
                if missing:
                    raise ValueError(
                        f"{path}: junction TSV needs columns chrom/donor/acceptor/strand; "
                        f"missing {sorted(missing)}"
                    )
                continue
            row = dict(zip(header, fields))
            strand = row['strand'].strip()
            if strand not in ('+', '-'):
                continue
            chrom = row['chrom'].strip()
            donor = int(row['donor'])
            acceptor = int(row['acceptor'])
            if strand == '+':
                intron_start, intron_end = donor, acceptor + 1
            else:
                intron_start, intron_end = acceptor, donor + 1
            if intron_end <= intron_start:
                continue
            junctions.append(make_junction(chrom, intron_start, intron_end, strand))
    return junctions


# ---------------------------------------------------------------------------------------------
# Read geometry helpers
# ---------------------------------------------------------------------------------------------

def rna_clip(read: pysam.AlignedSegment, strand: str) -> str:
    """The read-5' soft clip in RNA orientation, first base adjacent to the RNA 3' end.

    Under ``--rna3p-at read5p`` the gene strand is the inverse of the BAM strand, so a ``+`` gene
    read maps reverse and its read-5' terminus is the RIGHT end of the reference-forward
    ``query_sequence`` (returned as-is, because RNA == genome forward for a ``+`` gene); a ``-`` gene
    read maps forward and its read-5' terminus is the LEFT end (returned reverse-complemented, so
    index 0 is again the base adjacent to the aligned RNA 3' end).
    """
    if not read.cigartuples or not read.query_sequence:
        return ""
    seq = read.query_sequence.upper()
    if strand == '+':
        op, length = read.cigartuples[-1]
        return seq[len(seq) - length:] if op == 4 else ""
    op, length = read.cigartuples[0]
    return revcomp(seq[:length]) if op == 4 else ""


def _ref_to_query(read: pysam.AlignedSegment) -> Dict[int, int]:
    """Map reference position -> query index for the read's aligned (M/=/X) bases."""
    return {r: q for q, r in read.get_aligned_pairs(matches_only=True)}


def aligned_intronic_bases(
    read: pysam.AlignedSegment,
    strand: str,
    p: int,
    n_intronic: int,
    ref_to_query: Optional[Dict[int, int]] = None,
) -> str:
    """The ``n_intronic`` aligned READ bases lying at/past the donor, in RNA orientation.

    Ordered starting at the donor and running toward the RNA 3' end, i.e. the same direction as
    :func:`rna_clip`, so ``aligned_intronic_bases(...) + rna_clip(...)`` is contiguous RNA.  READ
    bases, not genome bases -- the aligner forced them onto intron sequence and it is the read that
    must be matched to exon 2.
    """
    if n_intronic <= 0:
        return ""
    seq = (read.query_sequence or "").upper()
    if ref_to_query is None:
        ref_to_query = _ref_to_query(read)
    # gene-strand positions donor .. p
    positions = [p - (n_intronic - 1) + i for i in range(n_intronic)] if strand == '+' \
        else [p + (n_intronic - 1) - i for i in range(n_intronic)]
    out = []
    for ref_pos in positions:
        q = ref_to_query.get(ref_pos)
        if q is None:
            return ""  # an indel/N inside the intronic stub: refuse rather than guess
        base = seq[q]
        out.append(base if strand == '+' else base.translate(_COMPLEMENT))
    return "".join(out)


def genome_gene_strand(genome_seq: str, start: int, length: int, strand: str) -> str:
    """``length`` genome bases starting at ``start`` and running gene-strand-ward, RNA orientation."""
    if length <= 0:
        return ""
    if strand == '+':
        return genome_seq[start:start + length].upper()
    lo = max(0, start - length + 1)
    return revcomp(genome_seq[lo:start + 1].upper())


# ---------------------------------------------------------------------------------------------
# Non-templated tail
# ---------------------------------------------------------------------------------------------

@dataclass
class TailCall:
    """Result of :func:`call_tail`."""

    tail_len: int = 0            # clip_a_run + walkback
    clip_a_run: int = 0          # A-run of the clip adjacent to the alignment (randomer stripped)
    walkback: int = 0            # aligned A's over genomic A's recovered 5'-ward from p
    tail_class: str = "none"     # none | A_only | A_rich | other
    clip_rna: str = ""           # the whole read-5' clip in RNA orientation
    tail_region: str = ""        # clip_rna with the distal randomer removed
    randomer: str = ""           # the distal umi_length nt of clip_rna (declared, not detected)


def call_tail(
    read: pysam.AlignedSegment,
    strand: str,
    p: int,
    genome_seq: Optional[str],
    umi_length: int = 0,
    clip_rna: Optional[str] = None,
    ref_to_query: Optional[Dict[int, int]] = None,
    max_walkback: int = 60,
    require_clip_evidence: bool = False,
) -> TailCall:
    """Call the non-templated 3' tail of one NET-seq read.

    The tail sits between the last genomically-encoded base and the randomer: in RNA orientation the
    clip is ``[tail A's][randomer]``, so the randomer is the DISTAL part and is stripped first --
    with ``umi_length=6`` a 6-nt clip is randomer-only and contributes no tail, which is exactly the
    ~30 % class of planning ``829`` §4 that the legacy 0.8-A-fraction test mis-handles.

    The walkback (Chanfreau invariant 7) then recovers tail A's the aligner placed on genomic A's:
    walking 5'-ward in RNA from ``p`` while the read base AND the genome base are both ``A``.  It is
    gated on the clip's tail region being entirely A -- if a non-A base sits 3' of the A-run, the RNA
    3' terminus is bounded there and the aligned A's cannot belong to the same tail.

    Args:
        read: pysam AlignedSegment.
        strand: GENE strand ('+' / '-').
        p: RNA 3' end in reference coordinates (already X-mismatch trimmed if that is wanted).
        genome_seq: chromosome sequence; ``None`` disables the walkback (clip A-run only).
        umi_length: declared randomer length; 0 = no randomer.
        clip_rna: precomputed :func:`rna_clip` result (avoids recomputation).
        ref_to_query: precomputed reference->query map.
        max_walkback: safety bound on the walkback.
        require_clip_evidence: only walk back when the clip carries at least one non-templated A
            adjacent to the alignment. OFF by default, which is invariant 7 / rectify's own
            ``walkback_3prime`` (a terminal read A over a genomic A is deliberately NOT skipped).
            Turn it ON for a nascent-RNA library where most 3' ends carry no tail at all: MEASURED
            on wt_rep3 chrI+chrII, 41,711 of 42,644 walkbacks had no clip evidence, and at RPL32
            (whose exon 1 ends ...AAAA) 24 of the 33 reads sitting on the exon-1 3' end -- 22 of
            them classified as splicing intermediates -- were walked 4 nt off it, erasing the peak.
    """
    clip = rna_clip(read, strand) if clip_rna is None else clip_rna
    randomer = ""
    if umi_length > 0 and clip:
        randomer = clip[max(0, len(clip) - umi_length):]
        tail_region = clip[:max(0, len(clip) - umi_length)]
    else:
        tail_region = clip

    a_run = 0
    while a_run < len(tail_region) and tail_region[a_run] == 'A':
        a_run += 1

    walkback = 0
    _evidence_ok = (a_run > 0) or not require_clip_evidence
    if genome_seq is not None and a_run == len(tail_region) and _evidence_ok:
        # The whole non-randomer clip region is A (possibly empty) -> the tail may continue into
        # bases the aligner placed on the genome.
        seq = (read.query_sequence or "").upper()
        if ref_to_query is None:
            ref_to_query = _ref_to_query(read)
        step = -1 if strand == '+' else 1
        read_a, genome_a = ('A', 'A') if strand == '+' else ('T', 'T')
        pos = p
        while walkback < max_walkback and 0 <= pos < len(genome_seq):
            q = ref_to_query.get(pos)
            if q is None:
                break
            if seq[q] != read_a or genome_seq[pos].upper() != genome_a:
                break
            walkback += 1
            pos += step

    tail_len = a_run + walkback
    if tail_region:
        n_a = tail_region.count('A')
        if n_a == len(tail_region):
            tail_class = "A_only"
        elif n_a / len(tail_region) >= 0.5:
            tail_class = "A_rich"
        else:
            tail_class = "other"
    else:
        tail_class = "A_only" if walkback > 0 else "none"

    return TailCall(
        tail_len=tail_len, clip_a_run=a_run, walkback=walkback, tail_class=tail_class,
        clip_rna=clip, tail_region=tail_region, randomer=randomer,
    )


# ---------------------------------------------------------------------------------------------
# Donor-side junction rescue
# ---------------------------------------------------------------------------------------------

#: Read classes emitted by :func:`rescue_read`.
RESCUE_CLASSES = (
    "none",             # no pooled donor within reach of this read's RNA 3' end
    "exon1_end",        # aligned end exactly at the exon-1 3' end, nothing matches exon 2
    "spliced_rescued",  # k exon-2 bases recovered; 3' end moved into exon 2
    "ambiguous",        # k >= min_k but the non-templated remainder r is not an allowed length
    "intronic_end",     # alignment runs into the intron and nothing matches exon 2
)


@dataclass
class RescueCall:
    """Result of :func:`rescue_read`."""

    status: str = "none"
    k: int = 0                       # exon-2 bases recovered
    r: int = 0                       # non-templated remainder, len(S) - k
    n_intronic: int = 0              # aligned bases at/past the donor
    position: Optional[int] = None   # rescued 3' end (reference coordinate); None unless rescued
    decoy_k: int = 0                 # LCP against a decoy acceptor -- the chance-match null
    decoy_would_rescue: bool = False  # the SAME acceptance rule applied to the decoy acceptor
    junction: Optional[NetseqJunction] = None
    s_seq: str = ""                  # the matched string S, RNA orientation


def allowed_remainders(umi_length: int) -> Tuple[int, ...]:
    """Non-templated remainder lengths a rescue may leave over.

    ``0`` always (a clean overhang).  With a declared randomer of ``N`` nt, also ``N-1``, ``N``,
    ``N+1``: the randomer is non-templated, and planning ``478`` §1 measured that its terminal 1-2 nt
    align by chance in a 1 : 0.28 : 0.09 ratio, so a 6-nt randomer presents as 5, 6 or 7 leftover nt.
    """
    if umi_length <= 0:
        return (0,)
    return tuple(sorted({0, umi_length - 1, umi_length, umi_length + 1}))


def rescue_read(
    read: pysam.AlignedSegment,
    strand: str,
    chrom: str,
    p: int,
    pool: Optional[JunctionPool],
    genome_seq: Optional[str],
    umi_length: int = 0,
    max_intronic: int = 10,
    min_k: int = 1,
    min_k_with_remainder: int = 4,
    clip_rna: Optional[str] = None,
    ref_to_query: Optional[Dict[int, int]] = None,
    decoy_offset: int = 50,
) -> RescueCall:
    """Match a NET-seq read's RNA-3' clip to the start of exon 2 across a pooled junction.

    ``S`` = the aligned bases at/past the donor (RNA orientation, donor first) + the read-5' soft
    clip in RNA orientation.  ``k`` = the longest prefix of ``S`` equal to the genome from
    ``exon2_first`` gene-strand-ward; ``r = len(S) - k`` is the non-templated remainder.  A rescue is
    accepted when ``r`` is in :func:`allowed_remainders` and ``k`` clears a **remainder-aware**
    floor, and puts the 3' end ``k - 1`` nt into exon 2.

    **The floor is remainder-aware because the chance channel IS the remainder.** With ``r == 0``
    every base of the clip is exon-2 sequence, so even a 1-nt overhang with nothing beyond it is
    real evidence — ``min_k`` applies.  With ``r > 0`` a randomer is being invoked to explain the
    rest of the clip, and a randomer's first base matches exon 2 a quarter of the time; the decoy
    null bears this out (decoy ``k=1`` 70 against 61 observed, ``k=2`` 24 against 4, ``k=3`` 10
    against 3, on a 504-read candidate set), so ``min_k_with_remainder`` applies instead.

    ``decoy_k`` is the same LCP computed against a decoy acceptor ``decoy_offset`` nt further into
    exon 2, and ``decoy_would_rescue`` runs the identical acceptance rule against it.  That is the
    chance-match floor for the rescue count, and is why the accepted-``k`` distribution can be read
    as evidence rather than assertion.
    """
    if pool is None or genome_seq is None:
        return RescueCall()

    hit = pool.lookup(chrom, strand, p, max_intronic)
    if hit is None:
        return RescueCall()
    n_intronic, junction = hit

    intronic = aligned_intronic_bases(read, strand, p, n_intronic, ref_to_query=ref_to_query)
    if n_intronic > 0 and not intronic:
        # indel inside the intronic stub -- do not guess
        return RescueCall(status="none", n_intronic=n_intronic, junction=junction)
    clip = rna_clip(read, strand) if clip_rna is None else clip_rna
    s_seq = intronic + clip

    exon2 = genome_gene_strand(genome_seq, junction.exon2_first, len(s_seq), strand)
    k = longest_common_prefix(s_seq, exon2)
    r = len(s_seq) - k

    decoy_start = junction.exon2_first + (decoy_offset if strand == '+' else -decoy_offset)
    decoy = genome_gene_strand(genome_seq, decoy_start, len(s_seq), strand)
    decoy_k = longest_common_prefix(s_seq, decoy)
    allowed = allowed_remainders(umi_length)

    def _floor(rr: int) -> int:
        """The k floor that applies to a remainder of ``rr``."""
        return min_k if rr == 0 else min_k_with_remainder

    def _accept(kk: int, rr: int) -> bool:
        """The acceptance rule, applied identically to the real and the decoy acceptor."""
        return kk >= _floor(rr) and rr in allowed

    decoy_r = len(s_seq) - decoy_k
    decoy_would_rescue = _accept(decoy_k, decoy_r)
    common = dict(k=k, r=r, n_intronic=n_intronic, decoy_k=decoy_k,
                  decoy_would_rescue=decoy_would_rescue, junction=junction, s_seq=s_seq)

    # ORDER MATTERS: the k FLOOR is tested first, the remainder second.
    #
    # A read at a donor with a 6-nt randomer clip whose first 2 bases match exon 2 by chance has
    # k=2, r=4 -- and 4 is not an allowed remainder, so a remainder-first test called it
    # `ambiguous` (measured on the full wt_rep3 library: 66 reads at k=2, 12 at k=3). Those are
    # genuine exon-1 3' ends carrying a randomer, and hiding them in `ambiguous` removes them from
    # the splicing-intermediate tally. Below the floor the match is not evidence of ANYTHING, so
    # the read simply keeps the class its geometry implies.
    if k < _floor(r):
        if n_intronic == 0:
            return RescueCall(status="exon1_end", **common)
        return RescueCall(status="intronic_end", **common)
    if r not in allowed:
        # Cleared the floor, but the leftover is not a length any randomer can explain --
        # genuinely undecidable, and kept out of BOTH the rescued and the exon1_end tallies.
        return RescueCall(status="ambiguous", **common)
    step = 1 if strand == '+' else -1
    return RescueCall(status="spliced_rescued",
                      position=junction.exon2_first + step * (k - 1), **common)


# ---------------------------------------------------------------------------------------------
# Per-read record
# ---------------------------------------------------------------------------------------------

@dataclass
class NetseqReadRecord(UnifiedReadRecord):
    """:class:`UnifiedReadRecord` plus the NET-seq rescue / tail fields.

    A subclass rather than new fields on the shared record so that the parquet schema of the DRS,
    cDNA and short-read pipelines is untouched.
    """

    rescue_status: str = "none"
    rescue_k: int = 0
    rescue_r: int = 0
    rescue_n_intronic: int = 0
    rescue_decoy_k: int = 0
    rescue_decoy_would_rescue: bool = False
    rescue_intron_start: int = -1
    rescue_intron_end: int = -1
    tail_len: int = 0
    tail_clip_a_run: int = 0
    tail_walkback: int = 0
    tail_class: str = "none"
    clip_rna: str = ""
    three_prime_shift: int = 0   # corrected - raw, signed in reference coordinates


# ---------------------------------------------------------------------------------------------
# Per-sample summary
# ---------------------------------------------------------------------------------------------

@dataclass
class NetseqCorrectionSummary:
    """Compact per-sample counters written to ``<sample>.netseq_summary.json``."""

    reads: int = 0
    junction_pool_size: int = 0
    reads_near_donor: int = 0
    rescued: int = 0
    rescued_mis_extended: int = 0        # rescued reads that had aligned bases inside the intron
    rescued_also_tailed: int = 0         # rescued reads whose clip also looked like a tail
    exon1_end: int = 0
    ambiguous: int = 0
    intronic_end: int = 0
    rescued_by_k: Dict[int, int] = field(default_factory=dict)
    rescued_by_k_clean: Dict[int, int] = field(default_factory=dict)   # r == 0  (min_k floor)
    rescued_by_k_randomer: Dict[int, int] = field(default_factory=dict)  # r > 0 (min_k_with_remainder floor)
    near_donor_k_clean: Dict[int, int] = field(default_factory=dict)   # ALL candidates, r == 0
    near_donor_k_randomer: Dict[int, int] = field(default_factory=dict)  # ALL candidates, r > 0
    decoy_k_at_donor: Dict[int, int] = field(default_factory=dict)     # chance-match null
    decoy_rescued: int = 0               # the SAME rule applied to a decoy acceptor 50 nt away
    rescued_by_r: Dict[int, int] = field(default_factory=dict)
    tailed_ge1: int = 0
    tailed_ge3: int = 0
    tail_walkback_gt0: int = 0
    tail_clip_evidence: int = 0          # >=1 non-templated A in the clip adjacent to the alignment
    tail_walkback_only: int = 0          # tail called from aligned A's ALONE -- no clip evidence
    tail_len_hist: Dict[int, int] = field(default_factory=dict)
    tail_class_hist: Dict[str, int] = field(default_factory=dict)
    xtrim_reads: int = 0                 # reads with a terminal oligo(A) X/MD mismatch trim
    end_moved: int = 0                   # corrected end != raw end
    end_shift_hist: Dict[int, int] = field(default_factory=dict)  # |shift| (gene-strand nt)

    def observe(self, record) -> None:
        """Fold one finished record into the counters.

        Reads the NET-seq fields with ``getattr`` defaults so a plain
        :class:`~rectify.core.unified_record.UnifiedReadRecord` (rescue and tail both disabled)
        still aggregates instead of raising.
        """
        self.reads += 1
        record = _RecordView(record)
        st = record.rescue_status
        if st != "none":
            self.reads_near_donor += 1
        if st == "spliced_rescued":
            self.rescued += 1
            _bump(self.rescued_by_k, record.rescue_k)
            _bump(self.rescued_by_r, record.rescue_r)
            if record.rescue_r == 0:
                _bump(self.rescued_by_k_clean, record.rescue_k)
            else:
                _bump(self.rescued_by_k_randomer, record.rescue_k)
            if record.rescue_n_intronic > 0:
                self.rescued_mis_extended += 1
            if record.tail_len > 0:
                self.rescued_also_tailed += 1
        elif st == "exon1_end":
            self.exon1_end += 1
        elif st == "ambiguous":
            self.ambiguous += 1
        elif st == "intronic_end":
            self.intronic_end += 1
        if st != "none":
            _bump(self.decoy_k_at_donor, record.rescue_decoy_k)
            if record.rescue_k >= 1:
                if record.rescue_r == 0:
                    _bump(self.near_donor_k_clean, record.rescue_k)
                else:
                    _bump(self.near_donor_k_randomer, record.rescue_k)
            if record.rescue_decoy_would_rescue:
                self.decoy_rescued += 1

        if record.tail_len >= 1:
            self.tailed_ge1 += 1
        if record.tail_len >= 3:
            self.tailed_ge3 += 1
        if record.tail_walkback > 0:
            self.tail_walkback_gt0 += 1
        if getattr(record, 'tail_clip_a_run', 0) > 0:
            self.tail_clip_evidence += 1
        elif record.tail_walkback > 0:
            self.tail_walkback_only += 1
        _bump(self.tail_len_hist, record.tail_len)
        self.tail_class_hist[record.tail_class] = self.tail_class_hist.get(record.tail_class, 0) + 1
        if record.aligned_a_length > 0:
            self.xtrim_reads += 1
        if record.three_prime_shift != 0:
            self.end_moved += 1
            _bump(self.end_shift_hist, abs(record.three_prime_shift))

    def as_dict(self) -> Dict:
        """JSON-ready dict (histogram keys become strings, as JSON requires)."""
        out: Dict = {}
        for key, value in self.__dict__.items():
            if isinstance(value, dict):
                out[key] = {str(k): v for k, v in sorted(value.items(), key=lambda kv: str(kv[0]))}
            else:
                out[key] = value
        frac = (lambda n: round(n / self.reads, 6)) if self.reads else (lambda n: None)
        out['frac_tailed_ge1'] = frac(self.tailed_ge1)
        out['frac_tailed_ge3'] = frac(self.tailed_ge3)
        out['frac_end_moved'] = frac(self.end_moved)
        return out


def _bump(hist: Dict[int, int], key: int) -> None:
    hist[key] = hist.get(key, 0) + 1


class _RecordView:
    """Read-only view giving a record the NET-seq fields' defaults when they are absent."""

    _DEFAULTS = {
        'rescue_status': 'none', 'rescue_k': 0, 'rescue_r': 0, 'rescue_n_intronic': 0,
        'rescue_decoy_k': 0, 'rescue_decoy_would_rescue': False,
        'tail_len': 0, 'tail_clip_a_run': 0, 'tail_walkback': 0, 'tail_class': 'none',
        'aligned_a_length': 0, 'three_prime_shift': 0,
    }

    __slots__ = ('_r',)

    def __init__(self, record):
        self._r = record

    def __getattr__(self, name):
        return getattr(self._r, name, self._DEFAULTS.get(name, 0))
