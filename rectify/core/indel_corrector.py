#!/usr/bin/env python3
"""
Indel artifact correction for RECTIFY.

This module implements detection and correction of insertion/deletion artifacts
at A-tracts. These artifacts occur when aligners try to force-align poly(A) tails
to genomic A-tracts, creating spurious indels.

Problem: Aligners (minimap2, BWA, STAR) create deletion artifacts when aligning
poly(A) tails to genomic sequences. Example:
    Read:    ATTAAAAAA
    Genomic: ATTTAAAAA
    Aligned: ATT--A--AAAA (aligner removes genomic TT)

This shifts the apparent 3' end position and must be corrected.

Three correction approaches are provided:
1. `detect_indel_artifacts()` + `correct_position_for_indels()`:
   Detects small indels near 3' end and corrects for them

2. `find_polya_boundary()`:
   Genome-aware algorithm that walks backwards from mapped 3' end,
   comparing read to genome, to find the TRUE CPA (cleavage and
   polyadenylation) site where they agree on a non-A base.

3. `rescue_softclip_at_homopolymer()`:
   Rescues soft-clipped bases adjacent to homopolymers. When basecaller
   under-calls homopolymer length, aligner may soft-clip the adjacent
   non-homopolymer base(s) instead of placing an indel. This function
   extends the 3' end to include matching soft-clipped bases.
   Use for LOCAL aligners (minimap2 default, STAR, etc.).

4. `rescue_mismatch_inside_homopolymer()`:
   Rescues misaligned bases inside homopolymers. When using GLOBAL alignment
   (no soft-clipping), the aligner forces mismatches inside homopolymer tracts.
   This function finds non-homopolymer bases in the read that were force-aligned
   and treats them as the true 3' end.
   Use for GLOBAL aligners (BBMap, minimap2 with local=false).

The `find_polya_boundary()` approach is recommended for nanopore data as it
handles complex poly(A) alignment artifacts including multiple deletions,
mismatches (base-calling errors), and genomic A-tracts.

For comprehensive coverage, use `correct_3prime_position()` which applies
all strategies and picks the best correction.

Author: Kevin R. Roy
Date: 2026-03-09
Updated: 2026-03-16 - Added genome-aware poly-A boundary detection
Updated: 2026-03-26 - Added soft-clip rescue at homopolymer boundaries
Updated: 2026-03-26 - Added variant-aware two-pass homopolymer rescue
"""

from typing import Dict, List, Optional, Tuple
from collections import defaultdict, Counter
from dataclasses import dataclass, field
import pysam

from ..config import (
    INDEL_MAX_SIZE,
    INDEL_SEARCH_WINDOW,
    INDEL_FLANK_A_THRESHOLD,
    INDEL_MIN_FLANK_LENGTH,
)
from ..utils.alignment import extract_deletions, extract_insertions
from ..utils.genome import fetch_genomic_sequence
from ..config import CHROM_TO_GENOME


# =============================================================================
# Genome-aware Poly-A Boundary Detection
# =============================================================================

def find_polya_boundary(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Dict[str, str],
    min_polya_len: int = 5
) -> Optional[Dict]:
    """
    Find the true 3' end by comparing mRNA and genome sequences.

    This is the recommended approach for nanopore data. It handles complex
    poly(A) alignment artifacts by walking backwards from the mapped 3' end
    and finding the first position where the read and genome AGREE on a non-A base.

    Algorithm:
    1. Build aligned positions from CIGAR (read_base, genome_base pairs)
    2. Start at the mapped 3' end
    3. Walk BACKWARDS toward the mRNA body
    4. Find the FIRST non-A base where genome and mRNA AGREE
    5. That's the true CPA

    Example (+ strand):
        mRNA:    CTGACGATGAAAAAGAAAATAAACAA-AAAAAAA
        genome:  CTGACGATGAAAAAGAAAATAAA-AAGAAAAAAA
                                   ^   ^
                                   T   C (mismatch: genome=A, read=C -> seq error)

        Walking backwards from 3' end:
        - Positions in poly-A region: A=A (match), A=A (match), etc.
        - Position 23: C vs A (MISMATCH) -> sequencing error, keep scanning
        - Position 19: T vs T (MATCH on non-A) -> TRUE CPA!

    Args:
        read: pysam AlignedSegment
        strand: '+' or '-'
        genome: Dict of chromosome sequences (NCBI format keys)
        min_polya_len: Minimum poly-A length to consider for correction (default 5)

    Returns:
        Dict with:
            - corrected_pos: True CPA position (0-based)
            - original_pos: Original mapped 3' position
            - polya_aligned_bp: Number of bp in the aligned poly-A region
            - correction_bp: How much the position was corrected
        Or None if no poly-A alignment artifacts found
    """
    cigar = read.cigartuples
    seq = read.query_sequence
    if not cigar or not seq:
        return None

    # Get genomic sequence
    chrom = read.reference_name
    genome_chrom = CHROM_TO_GENOME.get(chrom, chrom)
    genome_seq = genome.get(genome_chrom, '')
    if not genome_seq:
        return None

    # CIGAR ops: M=0, I=1, D=2, N=3, S=4, H=5

    if strand == '+':
        # For + strand, 3' end is at right side
        # Walk BACKWARDS from 3' end to find first non-A agreement

        ref_pos = read.reference_start
        read_pos = 0

        # Skip 5' soft-clip
        if cigar[0][0] == 4:
            read_pos = cigar[0][1]

        # Build list of aligned positions with both read and ref bases
        # (read_pos, ref_pos, read_base, genome_base)
        aligned_positions = []

        for op, length in cigar:
            if op == 4:  # Soft-clip
                continue
            elif op == 0 or op == 7 or op == 8:  # M, =, X - match/mismatch
                for i in range(length):
                    if ref_pos < len(genome_seq):
                        read_base = seq[read_pos].upper() if read_pos < len(seq) else 'N'
                        genome_base = genome_seq[ref_pos].upper()
                        aligned_positions.append((read_pos, ref_pos, read_base, genome_base))
                    read_pos += 1
                    ref_pos += 1
            elif op == 1:  # Insertion - consumes query only
                read_pos += length
            elif op == 2:  # Deletion - consumes reference only
                # Mark deletions (read has no base, genome has base)
                for i in range(length):
                    if ref_pos < len(genome_seq):
                        genome_base = genome_seq[ref_pos].upper()
                        aligned_positions.append((None, ref_pos, None, genome_base))
                    ref_pos += 1

        if not aligned_positions:
            return None

        # Walk BACKWARDS from the 3' end (last aligned position)
        # Find the FIRST non-A base where genome and read AGREE
        raw_3prime = read.reference_end - 1
        true_cpa_ref_pos = None

        for i in range(len(aligned_positions) - 1, -1, -1):
            rp, refp, rb, gb = aligned_positions[i]
            if rp is None:  # Skip deletions
                continue
            # Check: do genome and read agree on a non-A base?
            if rb == gb and gb != 'A':
                true_cpa_ref_pos = refp
                break

        if true_cpa_ref_pos is not None:
            polya_len = raw_3prime - true_cpa_ref_pos
            # Only return correction if we found significant poly-A
            if polya_len >= min_polya_len:
                return {
                    'corrected_pos': true_cpa_ref_pos,
                    'original_pos': raw_3prime,
                    'polya_aligned_bp': polya_len,
                    'correction_bp': raw_3prime - true_cpa_ref_pos,
                }

    else:  # strand == '-'
        # For - strand, 3' end is at left side (reference_start)
        # Poly-A appears as T's in genomic coordinates
        # Walk FORWARD from 3' end to find first non-T agreement

        ref_pos = read.reference_start
        read_pos = 0

        # Skip 3' soft-clip (at left for - strand)
        if cigar[0][0] == 4:
            read_pos = cigar[0][1]

        # Build aligned positions
        aligned_positions = []

        for op, length in cigar:
            if op == 4:  # Soft-clip
                continue
            elif op == 0 or op == 7 or op == 8:  # M, =, X
                for i in range(length):
                    if ref_pos < len(genome_seq):
                        read_base = seq[read_pos].upper() if read_pos < len(seq) else 'N'
                        genome_base = genome_seq[ref_pos].upper()
                        aligned_positions.append((read_pos, ref_pos, read_base, genome_base))
                    read_pos += 1
                    ref_pos += 1
            elif op == 1:  # Insertion
                read_pos += length
            elif op == 2:  # Deletion
                for i in range(length):
                    if ref_pos < len(genome_seq):
                        genome_base = genome_seq[ref_pos].upper()
                        aligned_positions.append((None, ref_pos, None, genome_base))
                    ref_pos += 1

        if not aligned_positions:
            return None

        # Walk FORWARD from the 3' end (first aligned position)
        # Find the FIRST non-T base where genome and read AGREE
        raw_3prime = read.reference_start
        true_cpa_ref_pos = None

        for i in range(len(aligned_positions)):
            rp, refp, rb, gb = aligned_positions[i]
            if rp is None:  # Skip deletions
                continue
            # Check: do genome and read agree on a non-T base?
            if rb == gb and gb != 'T':
                true_cpa_ref_pos = refp
                break

        if true_cpa_ref_pos is not None:
            polya_len = true_cpa_ref_pos - raw_3prime
            # Only return correction if we found significant poly-A
            if polya_len >= min_polya_len:
                return {
                    'corrected_pos': true_cpa_ref_pos,
                    'original_pos': raw_3prime,
                    'polya_aligned_bp': polya_len,
                    'correction_bp': true_cpa_ref_pos - raw_3prime,
                }

    return None


# =============================================================================
# Poly-A Prefix Rescue: extend past soft-clipped genomic A/T bases
# =============================================================================

def rescue_polya_prefix_in_softclip(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Dict[str, str],
    current_pos: Optional[int] = None,
    min_softclip_len: int = 5,
    min_polya_fraction: float = 0.7,
    max_rescue_bases: int = 5,
) -> Optional[Dict]:
    """
    Extend 3' end by rescuing genomic A/T bases incorrectly grouped with poly-A soft-clip.

    Problem: minimap2 cannot distinguish a genomic 'A' immediately adjacent to a
    T-tract from the poly-A tail. When the base at position P+1 (just past the
    alignment boundary) is 'A' and the poly-A tail also begins with 'A', minimap2
    soft-clips P+1 together with the poly-A tail. This shifts corrected_3prime 1 bp
    upstream (from the true CPA site at P+1 to the apparent site at P).

    Example (+ strand, center = first T of T-tract):
        Genome: ...[X][A][T][T][T][T][T]   <-- 'A' at center-1, T-tract at center+
        Read:   ...[X][A][A][A][A][A][A]   <-- poly-A tail

        minimap2 aligns X at center-2, then soft-clips the rest:
            [X] | [A][T][T][T][T][T][A][A][A][A][A]  (soft-clipped)
             raw_pos=center-2

        True CPA is center-1, but appears as center-2.

    Fix: if genome[P+1] == softclip[0] == 'A', extend corrected_3prime by 1
    (continuing as long as both are 'A', up to max_rescue_bases).

    For minus strand reads (poly-A → T's in BAM soft-clip):
        Extends leftward when genome[P-1] == softclip[-1] == 'T'.
        This is a no-op for TRT loci (center+1 is always non-A) but handles
        symmetric cases on other minus strand loci.

    Args:
        read: pysam AlignedSegment
        strand: '+' or '-'
        genome: Dict of chromosome sequences (NCBI format keys)
        current_pos: Current corrected position to extend from (default: BAM raw position)
        min_softclip_len: Minimum soft-clip length to consider (default 5)
        min_polya_fraction: Minimum fraction of A's (+ strand) or T's (- strand) in
            soft-clip to confirm poly-A tail (default 0.7)
        max_rescue_bases: Maximum bases to rescue (default 5)

    Returns:
        Dict with:
            - corrected_pos: Extended position (0-based)
            - original_pos: Starting position (= current_pos or BAM raw)
            - rescued_bases: Number of bases rescued
            - rescued_seq: Rescued genomic sequence
        Or None if no rescue applicable
    """
    cigar = read.cigartuples
    seq = read.query_sequence
    if not cigar or not seq:
        return None

    chrom = read.reference_name
    genome_chrom = CHROM_TO_GENOME.get(chrom, chrom)
    genome_seq = genome.get(genome_chrom, '')
    if not genome_seq:
        return None

    if strand == '+':
        # 3' end is at right side; check right soft-clip
        if cigar[-1][0] != 4:
            return None

        softclip_len = cigar[-1][1]
        if softclip_len < min_softclip_len:
            return None

        softclip_seq = seq[-softclip_len:].upper()
        if softclip_seq.count('A') / softclip_len < min_polya_fraction:
            return None

        raw_pos = current_pos if current_pos is not None else (read.reference_end - 1)

        # Extend rightward while soft-clip base == genomic base == 'A'
        rescued = 0
        rescued_chars: List[str] = []
        for i in range(min(softclip_len, max_rescue_bases)):
            genome_pos = raw_pos + 1 + i
            if genome_pos >= len(genome_seq):
                break
            genome_base = genome_seq[genome_pos].upper()
            softclip_base = softclip_seq[i]
            if genome_base == softclip_base == 'A':
                rescued += 1
                rescued_chars.append('A')
            else:
                break

        if rescued == 0:
            return None

        return {
            'corrected_pos': raw_pos + rescued,
            'original_pos': raw_pos,
            'rescued_bases': rescued,
            'rescued_seq': ''.join(rescued_chars),
        }

    else:  # strand == '-'
        # 3' end is at left side; check left soft-clip
        # For minus strand reads in BAM, query_sequence is reverse-complemented,
        # so the poly-A tail (A's in RNA) appears as T's at the START of the sequence.
        if cigar[0][0] != 4:
            return None

        softclip_len = cigar[0][1]
        if softclip_len < min_softclip_len:
            return None

        softclip_seq = seq[:softclip_len].upper()
        if softclip_seq.count('T') / softclip_len < min_polya_fraction:
            return None

        raw_pos = current_pos if current_pos is not None else read.reference_start

        # Extend leftward while (soft-clip base from right end) == genomic base == 'T'
        rescued = 0
        rescued_chars: List[str] = []
        for i in range(min(softclip_len, max_rescue_bases)):
            genome_pos = raw_pos - 1 - i
            if genome_pos < 0:
                break
            genome_base = genome_seq[genome_pos].upper()
            # Walk from rightmost soft-clip base (closest to alignment boundary)
            softclip_base = softclip_seq[softclip_len - 1 - i]
            if genome_base == softclip_base == 'T':
                rescued += 1
                rescued_chars.append('T')
            else:
                break

        if rescued == 0:
            return None

        return {
            'corrected_pos': raw_pos - rescued,
            'original_pos': raw_pos,
            'rescued_bases': rescued,
            'rescued_seq': ''.join(reversed(rescued_chars)),
        }


# =============================================================================
# Soft-clip Rescue at Homopolymer Boundaries
# =============================================================================

def rescue_softclip_at_homopolymer(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Dict[str, str],
    min_homopolymer_len: int = 3,
    end: str = '3prime'
) -> Optional[Dict]:
    """
    Rescue soft-clipped bases at homopolymer boundaries.

    Problem: When the basecaller under-calls homopolymer length (e.g., calls 8 U's
    instead of 9), the aligner may soft-clip the adjacent non-homopolymer base(s)
    instead of placing an indel within the homopolymer. This causes the apparent
    RNA end to appear INSIDE the homopolymer tract rather than at its boundary.

    Example (+ strand, poly-A at 3' end):
        True RNA:   ...GCTTTTTTTTC-AAAAAAAAA  (C followed by poly-A)
        Basecalled: ...GCTTTTTTTCAAAAAAAA     (shorter poly-U -> shorter poly-A)
        Aligned:    ...GCTTTTTTTT|c           (C soft-clipped, 3' end inside T-tract)
                                ^-- apparent 3' end (wrong!)
                                 ^-- true 3' end (the C)

    This also applies to 5' ends when there are homopolymers near the TSS.

    This function:
    1. Detects soft-clip adjacent to homopolymer (any base: A, T, G, C)
    2. Compares soft-clipped bases to reference sequence beyond alignment
    3. Extends RNA end to include matching bases until mismatch or homopolymer entry

    Args:
        read: pysam AlignedSegment
        strand: '+' or '-'
        genome: Dict of chromosome sequences (NCBI format keys)
        min_homopolymer_len: Minimum homopolymer run to trigger rescue (default 3)
        end: '3prime' or '5prime' - which end of the RNA to rescue

    Returns:
        Dict with:
            - corrected_pos: Extended position (0-based)
            - original_pos: Original mapped position
            - rescued_bases: Number of soft-clipped bases rescued
            - rescued_seq: Sequence of rescued bases
            - adjacent_homopolymer: The homopolymer base at boundary (A/T/G/C)
            - end: Which end was rescued ('3prime' or '5prime')
        Or None if no rescue applicable
    """
    cigar = read.cigartuples
    seq = read.query_sequence
    if not cigar or not seq:
        return None

    # Get genomic sequence
    chrom = read.reference_name
    genome_chrom = CHROM_TO_GENOME.get(chrom, chrom)
    genome_seq = genome.get(genome_chrom, '')
    if not genome_seq:
        return None

    # CIGAR ops: M=0, I=1, D=2, N=3, S=4, H=5

    # Determine which side of the alignment to examine based on strand and end:
    # - 3' end: + strand = right side, - strand = left side
    # - 5' end: + strand = left side, - strand = right side
    use_right_side = (strand == '+' and end == '3prime') or (strand == '-' and end == '5prime')

    if use_right_side:
        # Right side: soft-clip at END of CIGAR, position at reference_end - 1
        if cigar[-1][0] != 4:  # No soft-clip at this end
            return None

        softclip_len = cigar[-1][1]
        if softclip_len == 0:
            return None

        # Get the soft-clipped sequence (rightmost bases of read)
        softclip_seq = seq[-softclip_len:].upper()

        # Get the aligned position (just before soft-clip)
        raw_pos = read.reference_end - 1  # 0-based, last aligned position

        # Check if there's a homopolymer at the alignment boundary
        # Look at the last few aligned positions in the genome
        boundary_start = max(0, raw_pos - min_homopolymer_len + 1)
        boundary_seq = genome_seq[boundary_start:raw_pos + 1].upper()

        if len(boundary_seq) < min_homopolymer_len:
            return None

        # Check if all bases in boundary region are the same (homopolymer)
        boundary_base = boundary_seq[-1]  # Base at boundary
        if not all(b == boundary_base for b in boundary_seq):
            return None  # Not a homopolymer at boundary

        # Now check soft-clipped bases against reference beyond alignment
        # KEY INSIGHT: The reference homopolymer may extend beyond the alignment.
        # We need to skip over remaining reference homopolymer bases to find where
        # the soft-clipped bases match.
        #
        # Example:
        #   Read aligned to: ...TTTTTTTT (ending at pos 369,394)
        #   Reference:       ...TTTTTTTTTTG (T-tract extends to 369,396, G at 369,397)
        #   Soft-clip:       GTTC...
        #   The G matches ref position 369,397, so rescue 3 bases (skip 2 T's + G)

        # First, find where the reference homopolymer ends
        homopolymer_extension = 0
        for offset in range(1, min(softclip_len + 10, len(genome_seq) - raw_pos)):
            ref_pos = raw_pos + offset
            if ref_pos >= len(genome_seq):
                break
            ref_base = genome_seq[ref_pos].upper()
            if ref_base != boundary_base:
                # Found end of homopolymer
                break
            homopolymer_extension += 1

        # Now try to match soft-clip to reference starting after homopolymer
        rescued_count = 0
        rescued_bases = []
        match_start = raw_pos + homopolymer_extension + 1

        for i in range(softclip_len):
            ref_pos = match_start + i
            if ref_pos >= len(genome_seq):
                break

            ref_base = genome_seq[ref_pos].upper()
            softclip_base = softclip_seq[i]

            # Check if soft-clipped base matches reference
            if softclip_base == ref_base:
                rescued_count += 1
                rescued_bases.append(softclip_base)
            else:
                # Mismatch - stop rescue
                # But only accept rescue if we got at least 1 base
                break

        if rescued_count == 0:
            return None

        # The corrected position is after all rescued bases
        # Include the skipped homopolymer bases in the total shift
        total_shift = homopolymer_extension + rescued_count

        return {
            'corrected_pos': raw_pos + total_shift,
            'original_pos': raw_pos,
            'rescued_bases': total_shift,  # Total bases rescued (homopolymer + matched)
            'rescued_seq': ''.join(rescued_bases),
            'adjacent_homopolymer': boundary_base,
            'homopolymer_extension': homopolymer_extension,  # How many T's were under-called
            'end': end,
        }

    else:  # use_left_side
        # Left side: soft-clip at START of CIGAR, position at reference_start
        if cigar[0][0] != 4:  # No soft-clip at this end
            return None

        softclip_len = cigar[0][1]
        if softclip_len == 0:
            return None

        # Get the soft-clipped sequence (leftmost bases of read)
        softclip_seq = seq[:softclip_len].upper()

        # Get the aligned position (just after soft-clip)
        raw_pos = read.reference_start  # 0-based, first aligned position

        # Check if there's a homopolymer at the alignment boundary
        # Look at the first few aligned positions in the genome
        boundary_end = min(len(genome_seq), raw_pos + min_homopolymer_len)
        boundary_seq = genome_seq[raw_pos:boundary_end].upper()

        if len(boundary_seq) < min_homopolymer_len:
            return None

        # Check if all bases in boundary region are the same (homopolymer)
        boundary_base = boundary_seq[0]  # Base at boundary
        if not all(b == boundary_base for b in boundary_seq):
            return None  # Not a homopolymer at boundary

        # Now check soft-clipped bases against reference beyond alignment
        # Similar to right side: we need to skip over remaining homopolymer bases
        # in the reference before matching soft-clipped bases.
        #
        # Walk through reference LEFTWARD from alignment start to find
        # where homopolymer ends

        # First, find where the reference homopolymer ends (leftward)
        homopolymer_extension = 0
        for offset in range(1, min(softclip_len + 10, raw_pos + 1)):
            ref_pos = raw_pos - offset
            if ref_pos < 0:
                break
            ref_base = genome_seq[ref_pos].upper()
            if ref_base != boundary_base:
                # Found end of homopolymer
                break
            homopolymer_extension += 1

        # Now try to match soft-clip to reference starting before homopolymer
        rescued_count = 0
        rescued_bases = []
        match_start = raw_pos - homopolymer_extension - 1

        # Walk through soft-clip BACKWARDS (from alignment boundary outward)
        for i in range(softclip_len):
            ref_pos = match_start - i
            if ref_pos < 0:
                break

            ref_base = genome_seq[ref_pos].upper()
            # Soft-clip bases: walk from right to left (alignment boundary outward)
            softclip_base = softclip_seq[softclip_len - 1 - i]

            # Check if soft-clipped base matches reference
            if softclip_base == ref_base:
                rescued_count += 1
                rescued_bases.append(softclip_base)
            else:
                # Mismatch - stop rescue
                break

        if rescued_count == 0:
            return None

        # The corrected position is before all rescued bases
        total_shift = homopolymer_extension + rescued_count

        return {
            'corrected_pos': raw_pos - total_shift,
            'original_pos': raw_pos,
            'rescued_bases': total_shift,  # Total bases rescued
            'rescued_seq': ''.join(reversed(rescued_bases)),  # Restore original order
            'adjacent_homopolymer': boundary_base,
            'homopolymer_extension': homopolymer_extension,
            'end': end,
        }


def rescue_mismatch_inside_homopolymer(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Dict[str, str],
    min_homopolymer_len: int = 3,
    end: str = '3prime'
) -> Optional[Dict]:
    """
    Rescue misaligned bases inside homopolymers from global alignment.

    Problem: When using global alignment (no soft-clipping, e.g., BBMap or minimap2
    with local=false), if the basecaller under-calls homopolymer length, the aligner
    places mismatches inside the homopolymer tract. The apparent RNA end is inside
    the homopolymer when it should be at its boundary.

    Example (- strand, T-tract = poly-A in RNA at 3' end):
        True RNA:       ...C-AAAAAAAAAA 3'  (ends with C before poly-A)
        Basecalled:     ...CAAAAAAAA       (under-called poly-U -> shorter poly-A)
        Reference:      ...TTTTTTTTTTG...  (T-tract; G is complement of C)
        Global aligned: ...TTT[G]TTTTTTT   (G forced to align inside T-tract)
                              ^ mismatch: read=G, ref=T

        The G misaligned inside the T-tract is actually the true 3' end (the C).

    This also applies to 5' ends when there are homopolymers near the TSS.

    This function:
    1. Checks if RNA end is inside a homopolymer tract
    2. Walks INTO alignment looking for non-homopolymer bases in the READ
    3. Extends RNA end to include those bases as the true boundary

    Note: This is complementary to `rescue_softclip_at_homopolymer()` which handles
    aligners that soft-clip. Use both for comprehensive coverage.

    Args:
        read: pysam AlignedSegment
        strand: '+' or '-'
        genome: Dict of chromosome sequences (NCBI format keys)
        min_homopolymer_len: Minimum homopolymer run to trigger rescue (default 3)
        end: '3prime' or '5prime' - which end of the RNA to rescue

    Returns:
        Dict with:
            - corrected_pos: Extended position (0-based)
            - original_pos: Original mapped position
            - rescued_bases: Number of misaligned bases rescued
            - rescued_seq: Sequence of rescued bases
            - homopolymer_base: The homopolymer base (A/T/G/C)
            - end: Which end was rescued ('3prime' or '5prime')
        Or None if no rescue applicable
    """
    cigar = read.cigartuples
    seq = read.query_sequence
    if not cigar or not seq:
        return None

    # Get genomic sequence
    chrom = read.reference_name
    genome_chrom = CHROM_TO_GENOME.get(chrom, chrom)
    genome_seq = genome.get(genome_chrom, '')
    if not genome_seq:
        return None

    # CIGAR ops: M=0, I=1, D=2, N=3, S=4, H=5

    # Determine which side of the alignment to examine based on strand and end:
    # - 3' end: + strand = right side, - strand = left side
    # - 5' end: + strand = left side, - strand = right side
    use_right_side = (strand == '+' and end == '3prime') or (strand == '-' and end == '5prime')

    # Build aligned positions (same for both sides)
    ref_pos = read.reference_start
    read_pos = 0

    # Skip soft-clip at start if present
    if cigar[0][0] == 4:
        read_pos = cigar[0][1]

    aligned_positions = []
    for op, length in cigar:
        if op == 4:  # Soft-clip
            continue
        elif op == 0 or op == 7 or op == 8:  # M, =, X
            for i in range(length):
                if ref_pos < len(genome_seq) and read_pos < len(seq):
                    read_base = seq[read_pos].upper()
                    genome_base = genome_seq[ref_pos].upper()
                    aligned_positions.append((read_pos, ref_pos, read_base, genome_base))
                read_pos += 1
                ref_pos += 1
        elif op == 1:  # Insertion
            read_pos += length
        elif op == 2:  # Deletion
            for i in range(length):
                if ref_pos < len(genome_seq):
                    genome_base = genome_seq[ref_pos].upper()
                    aligned_positions.append((None, ref_pos, None, genome_base))
                ref_pos += 1

    if not aligned_positions:
        return None

    if use_right_side:
        # Right side: position at reference_end - 1, walk backwards
        raw_pos = read.reference_end - 1  # 0-based, last aligned position

        # Check if we're inside a homopolymer at the boundary
        boundary_end = min(len(genome_seq), raw_pos + min_homopolymer_len + 1)
        boundary_seq = genome_seq[raw_pos:boundary_end].upper()

        if len(boundary_seq) < min_homopolymer_len:
            return None

        # Check if the position is inside a homopolymer
        homopolymer_base = boundary_seq[0]
        if not all(b == homopolymer_base for b in boundary_seq[:min_homopolymer_len]):
            return None  # Not inside a homopolymer

        # Walk BACKWARDS from end looking for non-homopolymer bases in the READ
        rescued_count = 0
        rescued_bases = []
        new_pos = raw_pos

        for i in range(len(aligned_positions) - 1, -1, -1):
            rp, refp, rb, gb = aligned_positions[i]
            if rp is None:  # Skip deletions
                continue

            # If read base is NOT the homopolymer base, it's likely the true end
            if rb != homopolymer_base:
                new_pos = refp
                rescued_count = raw_pos - refp
                rescued_bases.append(rb)
                break

            # If we've moved past the homopolymer region in the genome, stop
            if gb != homopolymer_base:
                break

        if rescued_count == 0:
            return None

        return {
            'corrected_pos': new_pos,
            'original_pos': raw_pos,
            'rescued_bases': rescued_count,
            'rescued_seq': ''.join(rescued_bases),
            'homopolymer_base': homopolymer_base,
            'end': end,
        }

    else:  # use_left_side
        # Left side: position at reference_start, walk forwards
        raw_pos = read.reference_start  # 0-based, first aligned position

        # Check if we're inside a homopolymer at the boundary
        boundary_start = max(0, raw_pos - min_homopolymer_len)
        boundary_seq = genome_seq[boundary_start:raw_pos + 1].upper()

        if len(boundary_seq) < min_homopolymer_len:
            return None

        # Check if the position is inside a homopolymer
        homopolymer_base = boundary_seq[-1]
        if not all(b == homopolymer_base for b in boundary_seq[-min_homopolymer_len:]):
            return None  # Not inside a homopolymer

        # Walk FORWARD from end looking for non-homopolymer bases in the READ
        rescued_count = 0
        rescued_bases = []
        new_pos = raw_pos

        for i in range(len(aligned_positions)):
            rp, refp, rb, gb = aligned_positions[i]
            if rp is None:  # Skip deletions
                continue

            # If read base is NOT the homopolymer base, it's likely the true end
            if rb != homopolymer_base:
                new_pos = refp
                rescued_count = refp - raw_pos
                rescued_bases.append(rb)
                break

            # If we've moved past the homopolymer region in the genome, stop
            if gb != homopolymer_base:
                break

        if rescued_count == 0:
            return None

        return {
            'corrected_pos': new_pos,
            'original_pos': raw_pos,
            'rescued_bases': rescued_count,
            'rescued_seq': ''.join(rescued_bases),
            'homopolymer_base': homopolymer_base,
            'end': end,
        }


# =============================================================================
# Variant-Aware Homopolymer Rescue (Two-Pass Approach)
# =============================================================================

@dataclass
class MismatchFrequency:
    """Track mismatch frequencies at a genomic position."""
    chrom: str
    position: int
    ref_base: str
    homopolymer_base: str
    total_reads: int = 0
    mismatch_counts: Dict[str, int] = field(default_factory=lambda: defaultdict(int))

    @property
    def mismatch_fraction(self) -> float:
        """Fraction of reads with any non-homopolymer base at this position."""
        if self.total_reads == 0:
            return 0.0
        total_mismatches = sum(self.mismatch_counts.values())
        return total_mismatches / self.total_reads

    @property
    def dominant_mismatch(self) -> Tuple[str, float]:
        """Return the most common mismatch base and its frequency."""
        if not self.mismatch_counts or self.total_reads == 0:
            return ('', 0.0)
        most_common = max(self.mismatch_counts.items(), key=lambda x: x[1])
        return (most_common[0], most_common[1] / self.total_reads)


class VariantAwareHomopolymerRescue:
    """
    Two-pass variant-aware homopolymer rescue.

    Problem: When rescuing mismatches inside homopolymers, we might incorrectly
    "rescue" true SNPs as basecalling artifacts. A real T>C SNP 1bp upstream of
    a T-tract would look identical to a basecalling error.

    Solution: Use read frequency to distinguish:
    - Basecalling errors: Stochastic, affect ~5-15% of reads at a position
    - True SNPs: Consistent, affect >80% of reads at the position
    - Heterozygous variants: ~50% of reads

    Two-pass approach:
    1. First pass: Scan all reads to build per-position mismatch frequencies
    2. Second pass: Only rescue at positions where mismatch frequency is LOW
       (indicating basecalling error, not true variant)

    Usage:
        rescue = VariantAwareHomopolymerRescue(min_variant_fraction=0.8)

        # Pass 1: Build frequency map
        for read in bam:
            rescue.scan_read(read, strand, genome)
        rescue.finalize_scan()

        # Pass 2: Apply rescue with variant filtering
        for read in bam:
            result = rescue.rescue_with_variant_filter(read, strand, genome)
    """

    def __init__(
        self,
        min_variant_fraction: float = 0.8,
        min_reads_for_variant_call: int = 5,
        min_homopolymer_len: int = 4,  # Increased from 3 for safety
        max_rescue_bases: int = 3,      # Limit rescue to 3bp max
    ):
        """
        Initialize variant-aware rescue.

        Args:
            min_variant_fraction: If mismatch at position exceeds this fraction,
                treat as true variant, not basecalling error (default 0.8 = 80%)
            min_reads_for_variant_call: Minimum reads at position to make variant
                call (default 5). Below this, flag as AMBIGUOUS.
            min_homopolymer_len: Minimum homopolymer length to trigger rescue
                (default 4, increased from 3 for safety)
            max_rescue_bases: Maximum bases to rescue in one operation (default 3,
                nanopore under-calling rarely exceeds this)
        """
        self.min_variant_fraction = min_variant_fraction
        self.min_reads_for_variant_call = min_reads_for_variant_call
        self.min_homopolymer_len = min_homopolymer_len
        self.max_rescue_bases = max_rescue_bases

        # Position frequency tracking: (chrom, position) -> MismatchFrequency
        self._mismatch_freq: Dict[Tuple[str, int], MismatchFrequency] = {}

        # Track potential variants for reporting
        self._potential_variants: List[MismatchFrequency] = []

        # Scan state
        self._scan_complete = False

    def scan_read(
        self,
        read: pysam.AlignedSegment,
        strand: str,
        genome: Dict[str, str],
        end: str = '3prime'
    ) -> None:
        """
        First pass: Scan a read to update mismatch frequencies.

        Call this for every read in the BAM before applying rescue.

        Args:
            read: pysam AlignedSegment
            strand: '+' or '-'
            genome: Dict of chromosome sequences
            end: '3prime' or '5prime' - which end to analyze
        """
        cigar = read.cigartuples
        seq = read.query_sequence
        if not cigar or not seq:
            return

        # Get genomic sequence
        chrom = read.reference_name
        genome_chrom = CHROM_TO_GENOME.get(chrom, chrom)
        genome_seq = genome.get(genome_chrom, '')
        if not genome_seq:
            return

        # Determine which side to examine
        use_right_side = (strand == '+' and end == '3prime') or (strand == '-' and end == '5prime')

        # Build aligned positions
        ref_pos = read.reference_start
        read_pos = 0
        if cigar[0][0] == 4:  # Skip soft-clip at start
            read_pos = cigar[0][1]

        aligned_positions = []
        for op, length in cigar:
            if op == 4:  # Soft-clip
                continue
            elif op == 0 or op == 7 or op == 8:  # M, =, X
                for i in range(length):
                    if ref_pos < len(genome_seq) and read_pos < len(seq):
                        read_base = seq[read_pos].upper()
                        genome_base = genome_seq[ref_pos].upper()
                        aligned_positions.append((read_pos, ref_pos, read_base, genome_base))
                    read_pos += 1
                    ref_pos += 1
            elif op == 1:  # Insertion
                read_pos += length
            elif op == 2:  # Deletion
                ref_pos += length

        if not aligned_positions:
            return

        # Find positions inside homopolymers where read has non-homopolymer base
        if use_right_side:
            raw_pos = read.reference_end - 1

            # Check if we're at/near a homopolymer
            boundary_end = min(len(genome_seq), raw_pos + self.min_homopolymer_len + 1)
            boundary_seq = genome_seq[raw_pos:boundary_end].upper()

            if len(boundary_seq) < self.min_homopolymer_len:
                return

            homopolymer_base = boundary_seq[0]
            if not all(b == homopolymer_base for b in boundary_seq[:self.min_homopolymer_len]):
                return  # Not inside a homopolymer

            # Walk backwards looking for mismatches
            for i in range(len(aligned_positions) - 1, -1, -1):
                rp, refp, rb, gb = aligned_positions[i]
                if rp is None:
                    continue

                # If read base differs from homopolymer, record it
                if rb != homopolymer_base:
                    key = (chrom, refp)
                    if key not in self._mismatch_freq:
                        self._mismatch_freq[key] = MismatchFrequency(
                            chrom=chrom,
                            position=refp,
                            ref_base=gb,
                            homopolymer_base=homopolymer_base
                        )

                    self._mismatch_freq[key].total_reads += 1
                    self._mismatch_freq[key].mismatch_counts[rb] += 1
                    break  # Only record first mismatch

                # Stop if we've moved past homopolymer region
                if gb != homopolymer_base:
                    break

        else:  # left side
            raw_pos = read.reference_start

            boundary_start = max(0, raw_pos - self.min_homopolymer_len)
            boundary_seq = genome_seq[boundary_start:raw_pos + 1].upper()

            if len(boundary_seq) < self.min_homopolymer_len:
                return

            homopolymer_base = boundary_seq[-1]
            if not all(b == homopolymer_base for b in boundary_seq[-self.min_homopolymer_len:]):
                return

            # Walk forward looking for mismatches
            for i in range(len(aligned_positions)):
                rp, refp, rb, gb = aligned_positions[i]
                if rp is None:
                    continue

                if rb != homopolymer_base:
                    key = (chrom, refp)
                    if key not in self._mismatch_freq:
                        self._mismatch_freq[key] = MismatchFrequency(
                            chrom=chrom,
                            position=refp,
                            ref_base=gb,
                            homopolymer_base=homopolymer_base
                        )

                    self._mismatch_freq[key].total_reads += 1
                    self._mismatch_freq[key].mismatch_counts[rb] += 1
                    break

                if gb != homopolymer_base:
                    break

    def finalize_scan(self) -> None:
        """
        Complete the scanning phase and identify potential variants.

        Call this after scanning all reads, before applying rescue.
        """
        self._potential_variants = []

        for key, freq in self._mismatch_freq.items():
            if freq.total_reads >= self.min_reads_for_variant_call:
                if freq.mismatch_fraction >= self.min_variant_fraction:
                    self._potential_variants.append(freq)

        self._scan_complete = True

    def is_likely_variant(self, chrom: str, position: int) -> Tuple[bool, str]:
        """
        Check if a position is likely a true variant (not basecalling error).

        Args:
            chrom: Chromosome name
            position: 0-based genomic position

        Returns:
            Tuple of (is_variant, confidence):
            - is_variant: True if likely a true variant
            - confidence: 'HIGH' (>80% of reads), 'MEDIUM' (15-80%),
                         'LOW' (<15% - likely basecalling error),
                         'UNKNOWN' (insufficient reads)
        """
        if not self._scan_complete:
            raise RuntimeError("Must call finalize_scan() before checking variants")

        key = (chrom, position)
        if key not in self._mismatch_freq:
            return (False, 'UNKNOWN')

        freq = self._mismatch_freq[key]

        if freq.total_reads < self.min_reads_for_variant_call:
            return (False, 'UNKNOWN')

        frac = freq.mismatch_fraction

        if frac >= self.min_variant_fraction:
            return (True, 'HIGH')  # Likely true variant
        elif frac >= 0.15:
            return (False, 'MEDIUM')  # Ambiguous
        else:
            return (False, 'LOW')  # Likely basecalling error - OK to rescue

    def rescue_with_variant_filter(
        self,
        read: pysam.AlignedSegment,
        strand: str,
        genome: Dict[str, str],
        end: str = '3prime'
    ) -> Optional[Dict]:
        """
        Second pass: Apply rescue with variant filtering.

        Like rescue_mismatch_inside_homopolymer(), but skips positions
        that are likely true variants based on population frequency.

        Args:
            read: pysam AlignedSegment
            strand: '+' or '-'
            genome: Dict of chromosome sequences
            end: '3prime' or '5prime'

        Returns:
            Dict with rescue result (same as rescue_mismatch_inside_homopolymer)
            plus additional fields:
            - variant_check: 'SKIPPED_LIKELY_VARIANT', 'RESCUED', or 'NO_RESCUE_NEEDED'
            - variant_confidence: confidence level from is_likely_variant()
            Or None if no rescue applicable
        """
        if not self._scan_complete:
            raise RuntimeError("Must call finalize_scan() before rescue")

        # Get the basic rescue result
        result = rescue_mismatch_inside_homopolymer(
            read, strand, genome, self.min_homopolymer_len, end
        )

        if result is None:
            return None

        # Check if rescued position is likely a variant
        chrom = read.reference_name
        rescued_pos = result['corrected_pos']

        is_variant, confidence = self.is_likely_variant(chrom, rescued_pos)

        if is_variant:
            # This looks like a true variant - don't rescue!
            return {
                **result,
                'variant_check': 'SKIPPED_LIKELY_VARIANT',
                'variant_confidence': confidence,
                'corrected_pos': result['original_pos'],  # Revert to original
                'rescued_bases': 0,
            }

        # Check rescue size limit
        if result['rescued_bases'] > self.max_rescue_bases:
            return {
                **result,
                'variant_check': 'RESCUE_TOO_LARGE',
                'variant_confidence': confidence,
                'corrected_pos': result['original_pos'],  # Revert
                'rescued_bases': 0,
            }

        # Safe to rescue
        return {
            **result,
            'variant_check': 'RESCUED',
            'variant_confidence': confidence,
        }

    def get_potential_variants(self) -> List[Dict]:
        """
        Get list of positions that look like true variants.

        These should be reported to the user as potential reference issues
        or strain-specific SNPs that should NOT be rescued.

        Returns:
            List of dicts with variant information
        """
        if not self._scan_complete:
            raise RuntimeError("Must call finalize_scan() first")

        variants = []
        for freq in self._potential_variants:
            dom_base, dom_frac = freq.dominant_mismatch
            variants.append({
                'chrom': freq.chrom,
                'position': freq.position,
                'ref_base': freq.ref_base,
                'homopolymer_base': freq.homopolymer_base,
                'total_reads': freq.total_reads,
                'mismatch_fraction': freq.mismatch_fraction,
                'dominant_mismatch_base': dom_base,
                'dominant_mismatch_fraction': dom_frac,
            })

        return sorted(variants, key=lambda x: (x['chrom'], x['position']))

    def get_statistics(self) -> Dict:
        """Get summary statistics from the scan."""
        if not self._scan_complete:
            raise RuntimeError("Must call finalize_scan() first")

        total_positions = len(self._mismatch_freq)
        positions_with_coverage = sum(
            1 for f in self._mismatch_freq.values()
            if f.total_reads >= self.min_reads_for_variant_call
        )

        return {
            'total_positions_scanned': total_positions,
            'positions_with_sufficient_coverage': positions_with_coverage,
            'potential_variants_detected': len(self._potential_variants),
            'min_variant_fraction_threshold': self.min_variant_fraction,
            'min_reads_for_variant_call': self.min_reads_for_variant_call,
        }


def process_bam_with_variant_aware_rescue(
    bam_path: str,
    genome: Dict[str, str],
    output_variants_path: Optional[str] = None,
    min_variant_fraction: float = 0.8,
    min_reads_for_variant_call: int = 5,
    region: Optional[str] = None,
) -> Tuple[VariantAwareHomopolymerRescue, Dict]:
    """
    Process a BAM file with variant-aware homopolymer rescue.

    This is a convenience function that handles the two-pass approach:
    1. Scans all reads to build mismatch frequency map
    2. Identifies potential variants (positions with high mismatch frequency)
    3. Returns the rescue object ready for applying corrections

    Args:
        bam_path: Path to BAM file
        genome: Dict of chromosome sequences
        output_variants_path: Optional path to write potential variants TSV
        min_variant_fraction: Threshold for calling a position as variant (default 0.8)
        min_reads_for_variant_call: Minimum reads at position (default 5)
        region: Optional region string (e.g., "chrI:1000-2000") to limit scan

    Returns:
        Tuple of (VariantAwareHomopolymerRescue object, statistics dict)

    Example:
        >>> rescue, stats = process_bam_with_variant_aware_rescue(
        ...     "sample.bam", genome, output_variants_path="potential_variants.tsv"
        ... )
        >>> print(f"Found {stats['potential_variants_detected']} potential variants")
        >>>
        >>> # Now use rescue object in second pass
        >>> with pysam.AlignmentFile("sample.bam") as bam:
        ...     for read in bam:
        ...         strand = '-' if read.is_reverse else '+'
        ...         result = rescue.rescue_with_variant_filter(read, strand, genome)
    """
    rescue = VariantAwareHomopolymerRescue(
        min_variant_fraction=min_variant_fraction,
        min_reads_for_variant_call=min_reads_for_variant_call,
    )

    # Pass 1: Scan all reads
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        fetch_kwargs = {}
        if region:
            fetch_kwargs['region'] = region

        for read in bam.fetch(**fetch_kwargs):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            strand = '-' if read.is_reverse else '+'
            rescue.scan_read(read, strand, genome, end='3prime')

    # Finalize scan
    rescue.finalize_scan()

    # Get statistics
    stats = rescue.get_statistics()

    # Write potential variants if requested
    if output_variants_path:
        variants = rescue.get_potential_variants()
        if variants:
            with open(output_variants_path, 'w') as f:
                header = ['chrom', 'position', 'ref_base', 'homopolymer_base',
                          'total_reads', 'mismatch_fraction',
                          'dominant_mismatch_base', 'dominant_mismatch_fraction']
                f.write('\t'.join(header) + '\n')
                for v in variants:
                    row = [
                        v['chrom'],
                        str(v['position']),
                        v['ref_base'],
                        v['homopolymer_base'],
                        str(v['total_reads']),
                        f"{v['mismatch_fraction']:.3f}",
                        v['dominant_mismatch_base'],
                        f"{v['dominant_mismatch_fraction']:.3f}",
                    ]
                    f.write('\t'.join(row) + '\n')

    return rescue, stats


def correct_rna_end_position(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Dict[str, str],
    end: str = '3prime',
    min_polya_len: int = 5,
    min_homopolymer_len: int = 3,
    variant_aware_rescue: Optional['VariantAwareHomopolymerRescue'] = None,
) -> Dict:
    """
    Comprehensive RNA end correction combining multiple rescue strategies.

    This function applies correction strategies to handle different aligner behaviors:
    - For 3' ends: poly-A boundary detection + soft-clip rescue + mismatch rescue
    - For 5' ends: soft-clip rescue + mismatch rescue (no poly-A detection)

    The strategies handle different alignment scenarios:
    - Soft-clip rescue: For aligners that soft-clip at homopolymer boundaries
    - Mismatch rescue: For global aligners that force-align, placing mismatches

    Args:
        read: pysam AlignedSegment
        strand: '+' or '-'
        genome: Dict of chromosome sequences (NCBI format keys)
        end: '3prime' or '5prime' - which end to correct
        min_polya_len: Minimum poly-A length to trigger boundary correction (default 5)
        min_homopolymer_len: Minimum homopolymer run to trigger rescue (default 3)
        variant_aware_rescue: Optional VariantAwareHomopolymerRescue object for
            filtering out likely true variants during mismatch rescue. If provided,
            positions with high mismatch frequency (>80% of reads) will NOT be
            rescued, preventing false correction of real SNPs.

    Returns:
        Dict with:
            - original_pos: Original mapped position (0-based)
            - corrected_pos: Final corrected position
            - end: Which end was corrected ('3prime' or '5prime')
            - polya_correction: Dict from find_polya_boundary() or None (3' only)
            - softclip_rescue: Dict from rescue_softclip_at_homopolymer() or None
            - mismatch_rescue: Dict from rescue_mismatch_inside_homopolymer() or None
            - total_correction_bp: Total position shift
            - variant_skipped: True if mismatch rescue was skipped due to likely variant
    """
    # Determine which side based on strand and end
    use_right_side = (strand == '+' and end == '3prime') or (strand == '-' and end == '5prime')

    # Get original position
    if use_right_side:
        original_pos = read.reference_end - 1  # 0-based inclusive
    else:
        original_pos = read.reference_start  # 0-based

    current_pos = original_pos
    polya_result = None

    # Apply poly-A boundary detection (only for 3' ends)
    if end == '3prime':
        polya_result = find_polya_boundary(read, strand, genome, min_polya_len)
        if polya_result:
            current_pos = polya_result['corrected_pos']

    # Apply soft-clip rescue (for local aligners that soft-clip)
    softclip_result = rescue_softclip_at_homopolymer(read, strand, genome, min_homopolymer_len, end)

    # Apply mismatch rescue (for global aligners that force-align)
    # Use variant-aware rescue if provided
    variant_skipped = False
    if variant_aware_rescue is not None:
        mismatch_result = variant_aware_rescue.rescue_with_variant_filter(
            read, strand, genome, end
        )
        if mismatch_result and mismatch_result.get('variant_check') == 'SKIPPED_LIKELY_VARIANT':
            variant_skipped = True
            mismatch_result = None  # Don't use this rescue
    else:
        mismatch_result = rescue_mismatch_inside_homopolymer(read, strand, genome, min_homopolymer_len, end)

    # Combine all corrections - use the one that moves end furthest from homopolymer
    # For 3' + strand or 5' - strand (right side): further = higher position
    # For 3' - strand or 5' + strand (left side): further = lower position
    candidate_positions = [current_pos]

    if softclip_result:
        candidate_positions.append(softclip_result['corrected_pos'])
    if mismatch_result:
        candidate_positions.append(mismatch_result['corrected_pos'])

    if use_right_side:
        final_pos = max(candidate_positions)
    else:
        final_pos = min(candidate_positions)

    return {
        'original_pos': original_pos,
        'corrected_pos': final_pos,
        'end': end,
        'polya_correction': polya_result,
        'softclip_rescue': softclip_result,
        'mismatch_rescue': mismatch_result,
        'total_correction_bp': abs(final_pos - original_pos),
        'variant_skipped': variant_skipped,
    }


def correct_3prime_position(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Dict[str, str],
    min_polya_len: int = 5,
    min_homopolymer_len: int = 3
) -> Dict:
    """
    Comprehensive 3' end correction (wrapper for correct_rna_end_position).

    See correct_rna_end_position() for full documentation.

    Returns:
        Dict with original_3prime, corrected_3prime, and correction details.
    """
    result = correct_rna_end_position(
        read, strand, genome, '3prime', min_polya_len, min_homopolymer_len
    )
    # Rename keys for backwards compatibility
    return {
        'original_3prime': result['original_pos'],
        'corrected_3prime': result['corrected_pos'],
        'polya_correction': result['polya_correction'],
        'softclip_rescue': result['softclip_rescue'],
        'mismatch_rescue': result['mismatch_rescue'],
        'total_correction_bp': result['total_correction_bp'],
    }


def correct_5prime_position(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Dict[str, str],
    min_homopolymer_len: int = 3
) -> Dict:
    """
    Comprehensive 5' end correction (wrapper for correct_rna_end_position).

    See correct_rna_end_position() for full documentation.

    Returns:
        Dict with original_5prime, corrected_5prime, and correction details.
    """
    result = correct_rna_end_position(
        read, strand, genome, '5prime', min_polya_len=0, min_homopolymer_len=min_homopolymer_len
    )
    return {
        'original_5prime': result['original_pos'],
        'corrected_5prime': result['corrected_pos'],
        'softclip_rescue': result['softclip_rescue'],
        'mismatch_rescue': result['mismatch_rescue'],
        'total_correction_bp': result['total_correction_bp'],
    }


# =============================================================================
# Small Indel Detection (Original Approach)
# =============================================================================

def is_atract_deletion(
    deletion: Dict,
    genome_seq: str,
    read_seq: str,
    strand: str
) -> bool:
    """
    Determine if a deletion is an A-tract alignment artifact.

    Artifact criteria (all aligners):
    1. Deletion size ≤ INDEL_MAX_SIZE (default 3bp)
    2. Deleted bases are NOT A's (+ strand) or T's (- strand)
    3. Flanking regions in read are A-rich (>INDEL_FLANK_A_THRESHOLD)

    Args:
        deletion: Deletion dict with 'pos', 'length', 'ref_seq' keys
        genome_seq: Genomic sequence at deletion site
        read_seq: Read sequence flanking deletion
        strand: Gene strand ('+' or '-')

    Returns:
        True if deletion is likely an A-tract artifact
    """
    # Check deletion size
    if deletion['length'] > INDEL_MAX_SIZE:
        return False

    # Get deleted bases
    deleted_bases = deletion.get('ref_seq', '').upper()
    if not deleted_bases:
        return False

    # Check if deleted bases are A's/T's (these are expected in A-tracts)
    target_base = 'A' if strand == '+' else 'T'
    if all(b == target_base for b in deleted_bases):
        # Deletion is all A's/T's - this is expected, not an artifact
        return False

    # Check flanking regions are A/T-rich (depending on strand)
    # This indicates the read has poly(A) sequence at this position
    # For minus strand reads, poly(A) appears as poly(T) in genomic orientation
    if len(read_seq) < INDEL_MIN_FLANK_LENGTH * 2:
        return False

    # Get flanks around deletion position in read
    del_pos = deletion['read_pos']
    left_start = max(0, del_pos - INDEL_MIN_FLANK_LENGTH)
    left_flank = read_seq[left_start:del_pos]
    right_flank = read_seq[del_pos:del_pos + INDEL_MIN_FLANK_LENGTH]

    # Check richness of flanks (A for + strand, T for - strand)
    flank_base = 'A' if strand == '+' else 'T'

    if len(left_flank) >= INDEL_MIN_FLANK_LENGTH:
        left_richness = left_flank.count(flank_base) / len(left_flank)
        if left_richness < INDEL_FLANK_A_THRESHOLD:
            return False

    if len(right_flank) >= INDEL_MIN_FLANK_LENGTH:
        right_richness = right_flank.count(flank_base) / len(right_flank)
        if right_richness < INDEL_FLANK_A_THRESHOLD:
            return False

    return True


def detect_indel_artifacts(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Optional[Dict[str, str]] = None
) -> List[Dict]:
    """
    Detect indel artifacts near the 3' end of read.

    Scans CIGAR for deletions/insertions within INDEL_SEARCH_WINDOW of 3' end.
    Classifies each as artifact or legitimate variation.

    Args:
        read: pysam AlignedSegment
        strand: Gene strand ('+' or '-')
        genome: Optional genome dict for sequence validation

    Returns:
        List of artifact dicts with:
            - type: 'deletion' or 'insertion'
            - pos: Position in read (0-based)
            - ref_pos: Position in reference
            - length: Indel size
            - is_artifact: True if classified as artifact
    """
    artifacts = []

    # Determine 3' end position in read
    read_len = read.query_length
    if strand == '+':
        # 3' end is at right (end of read)
        three_prime_pos = read_len
        search_start = max(0, three_prime_pos - INDEL_SEARCH_WINDOW)
        search_end = three_prime_pos
    else:
        # 3' end is at left (start of read)
        three_prime_pos = 0
        search_start = 0
        search_end = min(read_len, INDEL_SEARCH_WINDOW)

    # Extract deletions
    deletions = extract_deletions(read)
    for deletion in deletions:
        # Check if deletion is near 3' end
        del_pos = deletion['read_pos']
        if search_start <= del_pos <= search_end:
            # Check if it's an artifact
            read_seq = read.query_sequence
            genome_seq = ''
            if genome:
                # Get genomic sequence for validation
                chrom = read.reference_name
                ref_pos = deletion['ref_pos']
                genome_seq = fetch_genomic_sequence(
                    genome, chrom, ref_pos - 10, ref_pos + 10
                )

            is_artifact = is_atract_deletion(
                deletion, genome_seq, read_seq, strand
            )

            artifacts.append({
                'type': 'deletion',
                'pos': del_pos,
                'ref_pos': deletion['ref_pos'],
                'length': deletion['length'],
                'ref_seq': deletion.get('ref_seq', ''),
                'is_artifact': is_artifact,
            })

    # Extract insertions (less common but can occur)
    insertions = extract_insertions(read)
    for insertion in insertions:
        ins_pos = insertion['read_pos']
        if search_start <= ins_pos <= search_end:
            # Insertions near 3' end are suspicious if small
            is_artifact = insertion['length'] <= INDEL_MAX_SIZE

            artifacts.append({
                'type': 'insertion',
                'pos': ins_pos,
                'ref_pos': insertion['ref_pos'],
                'length': insertion['length'],
                'seq': insertion.get('seq', ''),
                'is_artifact': is_artifact,
            })

    return artifacts


def correct_position_for_indels(
    original_pos: int,
    indel_artifacts: List[Dict],
    strand: str
) -> Dict:
    """
    Correct 3' end position by accounting for DELETION artifacts only.

    IMPORTANT: Only deletions affect reference coordinates. Insertions add bases
    to the read but do NOT shift the reference position. Therefore, only deletion
    artifacts require position correction.

    For deletions in homopolymer regions (nanopore error):
    - The aligner skips genomic bases that should be in the alignment
    - This shifts reference_end to a higher coordinate than the true 3' end
    - Correction: move 3' position UPSTREAM by deletion length

    Args:
        original_pos: Original 3' end position (0-based)
        indel_artifacts: List of artifact dicts from detect_indel_artifacts()
        strand: Gene strand ('+' or '-')

    Returns:
        Dict with:
            - corrected_position: Position after removing artifacts
            - correction_bp: Number of bp corrected
            - n_deletions: Number of deletion artifacts corrected
            - n_insertions: Number of insertion artifacts (flagged only, no correction)
    """
    correction_bp = 0
    n_deletions = 0
    n_insertions = 0

    for artifact in indel_artifacts:
        if not artifact['is_artifact']:
            continue

        if artifact['type'] == 'deletion':
            # Deletion artifact: aligner skipped genomic bases
            # The reference coordinates are shifted downstream by deletion length
            # Correct by moving position UPSTREAM:
            #   + strand: upstream = leftward = subtract
            #   - strand: upstream = rightward = add
            if strand == '+':
                correction_bp -= artifact['length']
            else:
                correction_bp += artifact['length']
            n_deletions += 1

        elif artifact['type'] == 'insertion':
            # Insertion artifact: extra bases in read that aren't in reference
            # IMPORTANT: Insertions do NOT affect reference coordinates!
            # The reference_start and reference_end are unchanged by insertions.
            # We count them for QC but do NOT apply position correction.
            n_insertions += 1
            # No correction_bp change for insertions

    return {
        'corrected_position': original_pos + correction_bp,
        'correction_bp': correction_bp,
        'n_deletions': n_deletions,
        'n_insertions': n_insertions,
    }


def correct_indels_from_read(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Optional[Dict[str, str]] = None
) -> Dict:
    """
    Full indel artifact correction for a single read.

    This function:
    1. Detects indel artifacts near 3' end
    2. Calculates corrected 3' end position
    3. Returns summary statistics

    Args:
        read: pysam AlignedSegment
        strand: Gene strand ('+' or '-')
        genome: Optional genome dict for validation

    Returns:
        Dict with:
            - original_3prime: Original 3' end position (0-based)
            - corrected_3prime: Corrected position
            - correction_bp: Position shift (bp)
            - has_artifacts: True if artifacts detected
            - artifacts: List of artifact dicts
    """
    # Get original 3' end position
    if strand == '+':
        original_3prime = read.reference_end - 1  # 0-based inclusive
    else:
        original_3prime = read.reference_start  # 0-based

    # Detect artifacts
    artifacts = detect_indel_artifacts(read, strand, genome)

    # Count artifacts
    artifact_list = [a for a in artifacts if a['is_artifact']]
    has_artifacts = len(artifact_list) > 0

    # Correct position
    if has_artifacts:
        correction = correct_position_for_indels(
            original_3prime, artifact_list, strand
        )
        corrected_3prime = correction['corrected_position']
        correction_bp = correction['correction_bp']
    else:
        corrected_3prime = original_3prime
        correction_bp = 0

    return {
        'original_3prime': original_3prime,
        'corrected_3prime': corrected_3prime,
        'correction_bp': correction_bp,
        'has_artifacts': has_artifacts,
        'artifacts': artifact_list,
    }


def calculate_indel_statistics(indel_results: list) -> Dict:
    """
    Calculate summary statistics for indel correction.

    Args:
        indel_results: List of result dicts from correct_indels_from_read()

    Returns:
        Dict with summary statistics
    """
    if not indel_results:
        return {
            'total': 0,
            'with_artifacts': 0,
            'artifact_rate': 0.0,
            'mean_correction': 0.0,
        }

    import numpy as np

    with_artifacts = [r for r in indel_results if r['has_artifacts']]
    corrections = [abs(r['correction_bp']) for r in with_artifacts]

    return {
        'total': len(indel_results),
        'with_artifacts': len(with_artifacts),
        'artifact_rate': len(with_artifacts) / len(indel_results),
        'mean_correction': np.mean(corrections) if corrections else 0.0,
        'median_correction': np.median(corrections) if corrections else 0.0,
        'max_correction': max(corrections) if corrections else 0,
    }


def format_indel_report(stats: Dict) -> str:
    """
    Format indel correction statistics as human-readable report.

    Args:
        stats: Statistics dict from calculate_indel_statistics()

    Returns:
        Formatted report string
    """
    report = []
    report.append("=" * 60)
    report.append("Indel Artifact Correction Summary")
    report.append("=" * 60)
    report.append("")

    report.append("Overall:")
    report.append(f"  Total reads:            {stats['total']:,}")
    report.append(f"  With artifacts:         {stats['with_artifacts']:,} ({stats['artifact_rate']:.1%})")
    report.append("")

    if stats['with_artifacts'] > 0:
        report.append("Position Corrections:")
        report.append(f"  Mean:                   {stats['mean_correction']:.1f} bp")
        report.append(f"  Median:                 {stats['median_correction']:.1f} bp")
        report.append(f"  Maximum:                {stats['max_correction']} bp")

    report.append("")
    report.append("=" * 60)

    return "\n".join(report)
