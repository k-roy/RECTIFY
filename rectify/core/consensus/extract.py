"""
Alignment-info dataclasses and CIGAR → AlignmentInfo extraction.

The type system shared by ``scoring.py``, ``select.py``, and
``consensus.py``. Pure CIGAR/genome inspection — no aligner scoring logic.

Includes the small CIGAR-extract helpers ``extract_junctions_from_cigar``,
``get_softclip_lengths``, and ``check_canonical_splice_sites``. The
latter sits here (rather than in ``scoring.py`` as the original handoff
brief suggested) so that ``extract_alignment_info`` can call it without
creating a circular import: scoring depends on AlignmentInfo, so the
dependency graph must run ``extract → (nothing) → scoring → select``.

The module-level ``_atract_cache`` is keyed by ``(chrom, raw_3prime, strand)``;
many reads share their 3' endpoint, so caching avoids repeated A-tract
genome sequence lookups during ``extract_alignment_info``.
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pysam

logger = logging.getLogger(__name__)


# Module-level A-tract result cache: (chrom, pos, strand) → atract dict
# Many reads share the same 3' endpoint; this avoids repeated sequence lookups.
_atract_cache: Dict[Tuple[str, int, str], dict] = {}

# Canonical splice site dinucleotides
CANONICAL_5SS = {'GT', 'GC'}  # 5' splice site (donor)
CANONICAL_3SS = {'AG'}  # 3' splice site (acceptor)

_RC_TABLE = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def _revcomp(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1].upper()


@dataclass
class AlignmentInfo:
    """Stores key info about a read's alignment from one aligner."""
    read_id: str
    aligner: str
    chrom: str
    strand: str
    reference_start: int
    reference_end: int
    cigar_string: str

    # Junction information
    junctions: List[Tuple[int, int]] = field(default_factory=list)  # (start, end) pairs

    # Soft-clip info (5' and 3')
    five_prime_softclip: int = 0
    three_prime_softclip: int = 0

    # 3' end A-tract correction (computed pre-consensus using genome sequence only;
    # full indel correction requiring MD tags is applied post-consensus)
    corrected_3prime: Optional[int] = None   # estimated true CPA position
    three_prime_atract_depth: int = 0        # A's downstream of raw 3' end (0 = clean landing)

    # Soft-clipped sequence at the 5' end (for sequence-based rescue scoring)
    five_prime_softclip_seq: str = ""

    # Query bases in the terminal mismatch/indel region at the 5' end.
    # mapPacBio forces mismatches/indels at splice junction boundaries instead of
    # soft-clipping; those forced-error bases are structurally equivalent to a
    # soft-clip and should be used for sequence-based rescue if five_prime_softclip_seq
    # is empty.  Length = effective_five_prime_clip - five_prime_softclip.
    effective_five_prime_clip_seq: str = ""

    # Quality scores
    junction_score: float = 0.0
    canonical_count: int = 0
    non_canonical_count: int = 0

    # Effective 5' clip for scoring: max(explicit_soft_clip, terminal_error_region).
    # Some aligners (mapPacBio) substitute forced mismatches/indels for soft-clips at
    # splice junction boundaries — identical structural situation, different CIGAR encoding.
    # This field ensures fair scoring regardless of the aligner's clipping policy.
    effective_five_prime_clip: int = 0

    # Effective 3' clip for scoring: non-poly(A) terminal errors at the 3' end.
    # Aligners should only soft-clip the poly(A) tail; clipping or force-mismatching
    # real exon sequence means the aligner stopped before the true 3' end.
    # Complementary to the A-tract depth penalty (which catches going too far INTO poly(A)).
    effective_three_prime_clip: int = 0

    # Mismatches + indels within junction_window_bp of any splice junction.
    # Cleaner aligners (mapPacBio) score 0; aligners with forced errors near
    # junctions (minimap2, gapmm2) accumulate a penalty here.
    # Float because homopolymer-context errors are weighted 0.5 (nanopore DRS
    # commonly undercalls homopolymer lengths, making such errors less diagnostic
    # of misalignment than non-homopolymer mismatches).
    junction_proximity_errors: float = 0.0

    # Flags
    has_false_3prime_junction: bool = False
    is_best: bool = False


@dataclass
class ConsensusResult:
    """Result of consensus selection for a read."""
    read_id: str
    best_aligner: str
    best_alignment: Optional[AlignmentInfo]
    aligners_compared: List[str]

    # Consensus metrics
    n_aligners_agree: int = 0  # Number of aligners with same junctions
    n_tied_score: int = 1      # Number of aligners with equal top junction score
    tied_aligners: List[str] = field(default_factory=list)  # All aligners tied for top score
    confidence: str = ""  # 'high', 'medium', 'low'

    # Rescue info
    was_5prime_rescued: bool = False  # Did we pick an alignment that spliced vs clipped
    false_junction_removed: bool = False

    # All alignments for this read (for debugging)
    all_alignments: Dict[str, AlignmentInfo] = field(default_factory=dict)


def extract_junctions_from_cigar(read: pysam.AlignedSegment) -> List[Tuple[int, int]]:
    """
    Extract splice junctions from CIGAR N operations.

    Returns list of (intron_start, intron_end) in 0-based coords.
    intron_end is exclusive.
    """
    if not read.cigartuples:
        return []

    junctions = []
    ref_pos = read.reference_start

    for op, length in read.cigartuples:
        if op == 3:  # N = skipped region (intron)
            junctions.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op in (0, 2, 7, 8):  # M, D, =, X consume reference
            ref_pos += length

    return junctions


def get_softclip_lengths(read: pysam.AlignedSegment) -> Tuple[int, int]:
    """
    Get 5' and 3' soft-clip lengths for a read.

    For plus strand: 5' is left, 3' is right
    For minus strand: 5' is right, 3' is left (in aligned orientation)

    Returns: (five_prime_clip, three_prime_clip)
    """
    if not read.cigartuples:
        return (0, 0)

    left_clip = 0
    right_clip = 0

    # Left soft-clip
    if read.cigartuples[0][0] == 4:  # S = soft-clip
        left_clip = read.cigartuples[0][1]

    # Right soft-clip
    if read.cigartuples[-1][0] == 4:
        right_clip = read.cigartuples[-1][1]

    # Adjust for strand
    if read.is_reverse:
        return (right_clip, left_clip)  # 5' is on right for minus strand
    else:
        return (left_clip, right_clip)  # 5' is on left for plus strand


def check_canonical_splice_sites(
    junctions: List[Tuple[int, int]],
    chrom: str,
    genome: Dict[str, str],
    strand: str = '+',
) -> Tuple[int, int]:
    """
    Count canonical vs non-canonical splice sites for junctions.

    Returns: (canonical_count, non_canonical_count)
    """
    if chrom not in genome:
        return (0, 0)

    seq = genome[chrom]
    canonical = 0
    non_canonical = 0

    for start, end in junctions:
        if start < 0 or end > len(seq):
            continue

        if strand == '-':
            # Minus-strand transcript orientation: genomic right edge is the
            # donor and genomic left edge is the acceptor.
            five_ss = _revcomp(seq[end - 2:end].upper())
            three_ss = _revcomp(seq[start:start + 2].upper())
        else:
            # Plus-strand transcript orientation.
            five_ss = seq[start:start + 2].upper()
            three_ss = seq[end - 2:end].upper()

        if five_ss in CANONICAL_5SS and three_ss in CANONICAL_3SS:
            canonical += 1
        else:
            non_canonical += 1

    return (canonical, non_canonical)


def extract_alignment_info(
    read: pysam.AlignedSegment,
    aligner: str,
    genome: Dict[str, str],
) -> AlignmentInfo:
    """Extract alignment info from a pysam read.

    Computes corrected_3prime pre-consensus using A-tract detection (genome-only,
    no MD tags required). Full indel correction (MD-dependent) is applied
    post-consensus as a refinement step.
    """
    from ..polya.atract_detector import calculate_atract_ambiguity
    from .scoring import _get_effective_5prime_clip, _get_effective_3prime_clip, _count_junction_proximity_errors

    junctions = extract_junctions_from_cigar(read)
    five_clip, three_clip = get_softclip_lengths(read)

    # Effective 5' clip for scoring: includes terminal mismatch/indel region.
    # Must be computed before the AlignmentInfo is constructed.
    effective_five_clip = _get_effective_5prime_clip(read, genome)

    # Effective 3' clip for scoring: non-poly(A) terminal errors at 3' end.
    effective_three_clip = _get_effective_3prime_clip(read, genome)

    # Extract soft-clipped bases at the 5' end for sequence-based rescue
    five_prime_seq = ""
    if five_clip > 0 and read.query_sequence:
        if not read.is_reverse:
            five_prime_seq = read.query_sequence[:five_clip]
        else:
            five_prime_seq = read.query_sequence[-five_clip:]

    # Extract the terminal mismatch/indel region for mapPacBio rescue.
    # MPB forces mismatches/indels instead of soft-clipping at splice boundaries;
    # effective_five_clip > five_clip means there are terminal alignment errors.
    # We extract those query bases so _rescue_5prime_softclip can compare them
    # against the upstream exon end — same test as for an explicit soft-clip.
    effective_five_prime_seq = ""
    terminal_error_len = effective_five_clip - five_clip
    if terminal_error_len > 0 and read.query_sequence:
        qlen = len(read.query_sequence)
        if not read.is_reverse:
            # Plus strand: terminal errors are the first `terminal_error_len` aligned bases
            # (immediately after the explicit soft-clip region, if any)
            region_start = five_clip
            region_end = five_clip + terminal_error_len
            if region_end <= qlen:
                effective_five_prime_seq = read.query_sequence[region_start:region_end]
        else:
            # Minus strand: terminal errors are at the high end of query_sequence
            # (immediately before the explicit soft-clip region, if any)
            region_end = qlen - five_clip
            region_start = region_end - terminal_error_len
            if region_start >= 0:
                effective_five_prime_seq = read.query_sequence[region_start:region_end]

    chrom = read.reference_name
    strand = '-' if read.is_reverse else '+'
    canonical, non_canonical = check_canonical_splice_sites(junctions, chrom, genome, strand)

    # Estimate corrected 3' end using A-tract ambiguity detection.
    # Raw 3' end: reference_end - 1 for + strand, reference_start for - strand.
    raw_3prime = (read.reference_end - 1) if strand == '+' else read.reference_start
    corrected_3prime = raw_3prime
    atract_depth = 0

    chrom_std = chrom
    if genome.get(chrom_std) is None:
        # Try standardized name
        from ...utils.genome import standardize_chrom_name
        chrom_std = standardize_chrom_name(chrom) or chrom

    if genome.get(chrom_std) is not None:
        try:
            # Cache A-tract results by (chrom, pos, strand) — many reads share 3' ends
            _cache_key = (chrom_std, raw_3prime, strand)
            if _cache_key not in _atract_cache:
                _atract_cache[_cache_key] = calculate_atract_ambiguity(
                    genome, chrom_std, raw_3prime, strand, downstream_bp=10
                )
            atract = _atract_cache[_cache_key]
            atract_depth = atract.get('downstream_a_count') or 0
            # Best-guess corrected position: ambiguity_min for +, ambiguity_max for -
            if strand == '+':
                corrected_3prime = atract.get('ambiguity_min', raw_3prime)
            else:
                corrected_3prime = atract.get('ambiguity_max', raw_3prime)
        except Exception:
            pass  # Non-fatal; raw position used

    # Count mismatches/indels within 5 bp of each junction boundary.
    # Used by score_alignment() to prefer aligners with clean junction handling.
    junction_prox_errors = _count_junction_proximity_errors(read, genome)

    return AlignmentInfo(
        read_id=read.query_name,
        aligner=aligner,
        chrom=chrom,
        strand=strand,
        reference_start=read.reference_start,
        reference_end=read.reference_end,
        cigar_string=read.cigarstring or "",
        junctions=junctions,
        five_prime_softclip=five_clip,
        three_prime_softclip=three_clip,
        five_prime_softclip_seq=five_prime_seq,
        effective_five_prime_clip=effective_five_clip,
        effective_five_prime_clip_seq=effective_five_prime_seq,
        effective_three_prime_clip=effective_three_clip,
        corrected_3prime=corrected_3prime,
        three_prime_atract_depth=atract_depth,
        canonical_count=canonical,
        non_canonical_count=non_canonical,
        junction_proximity_errors=junction_prox_errors,
        has_false_3prime_junction=False,  # Not used; walk back handles this
    )
