#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
False Junction Filter for RECTIFY

This module detects and filters false splice junctions that are artifacts of
poly(A) tail mapping. Aligners sometimes create CIGAR N operations (junctions)
to map poly(A) tails to distant genomic A-tracts.

Problem:
    When a poly(A) tail is partially sequenced, aligners may:
    1. Map the poly(A) to a nearby genomic A-tract
    2. Create a false "intron" (N operation) to bridge the gap
    3. This results in spurious splice junctions near the 3' end

Detection Criteria:
    A junction is flagged as a likely poly(A) artifact if:
    1. It is near the read's 3' end (within POLYA_JUNCTION_WINDOW bp)
    2. The downstream region (after junction) is highly A-rich
    3. The junction target connects to a genomic A-tract
    4. The "intron" has atypical splice site motifs (not GT-AG)

Integration:
    - Run this filter BEFORE the standard junction extraction
    - Artifact junctions are treated like large deletions
    - The poly(A) length calculation absorbs the artifact
    - Position is corrected to the upstream exon boundary

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-24
"""

import logging
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
import pysam

logger = logging.getLogger(__name__)

from ..config import POLYA_RICHNESS_THRESHOLD
POLYPOLYA_RICHNESS_THRESHOLD = POLYA_RICHNESS_THRESHOLD  # legacy alias


# =============================================================================
# Genome reference adapter
# =============================================================================

class _GenomeDictReference:
    """Thin adapter so a genome dict can be passed where pysam.FastaFile is expected."""

    def __init__(self, genome: Dict[str, str]):
        self._genome = genome

    def fetch(self, chrom: str, start: int, end: int) -> str:
        seq = self._genome.get(chrom, '')
        return seq[max(0, start):end]

    def get_reference_length(self, chrom: str) -> int:
        return len(self._genome.get(chrom, ''))


# =============================================================================
# Configuration
# =============================================================================

POLYA_JUNCTION_WINDOW = 50      # bp from 3' end to consider for artifacts
# POLYA_RICHNESS_THRESHOLD sourced from config.POLYPOLYA_RICHNESS_THRESHOLD (0.8)
GENOMIC_ATRACT_THRESHOLD = 0.8  # A-richness in junction target


# =============================================================================
# Helper Functions
# =============================================================================

def calculate_a_richness(sequence: str) -> float:
    """
    Calculate A-richness in sequence.

    Args:
        sequence: DNA sequence (case-insensitive)

    Returns:
        A-richness (0.0-1.0)
    """
    if not sequence:
        return 0.0
    sequence = sequence.upper()
    return sequence.count('A') / len(sequence)


def calculate_t_richness(sequence: str) -> float:
    """
    Calculate T-richness in sequence (for minus strand).

    Args:
        sequence: DNA sequence (case-insensitive)

    Returns:
        T-richness (0.0-1.0)
    """
    if not sequence:
        return 0.0
    sequence = sequence.upper()
    return sequence.count('T') / len(sequence)


def get_read_3prime_position(read: pysam.AlignedSegment, strand: str) -> int:
    """
    Get the 3' end genomic position from a read.

    Args:
        read: pysam AlignedSegment
        strand: '+' or '-'

    Returns:
        3' end position (0-based)
    """
    if strand == '+':
        return read.reference_end - 1
    else:
        return read.reference_start


def extract_junctions_from_cigar(read: pysam.AlignedSegment) -> List[Tuple[int, int, int]]:
    """
    Extract splice junctions (N operations) from CIGAR.

    Args:
        read: pysam AlignedSegment

    Returns:
        List of (intron_start, intron_end, cigar_index) tuples
    """
    if read.cigartuples is None:
        return []

    junctions = []
    ref_pos = read.reference_start

    for i, (op, length) in enumerate(read.cigartuples):
        if op == 3:  # N = skipped region (intron)
            junctions.append((ref_pos, ref_pos + length, i))
            ref_pos += length
        elif op in (0, 2, 7, 8):  # M, D, =, X consume reference
            ref_pos += length
        # I, S, H don't consume reference

    return junctions


def get_sequence_after_junction(
    read: pysam.AlignedSegment,
    junction_cigar_index: int,
    length: int = 20,
) -> str:
    """
    Get the read sequence after a junction (rightward in cigar / genomic order).

    Args:
        read: pysam AlignedSegment
        junction_cigar_index: Index of junction in CIGAR
        length: Number of bases to extract

    Returns:
        Sequence after the junction
    """
    if read.query_sequence is None or read.cigartuples is None:
        return ""

    # Find query position after junction
    query_pos = 0
    for i, (op, op_len) in enumerate(read.cigartuples):
        if i > junction_cigar_index:
            break
        if op in (0, 1, 4, 7, 8):  # M, I, S, =, X consume query
            query_pos += op_len

    end_pos = min(query_pos + length, len(read.query_sequence))
    return read.query_sequence[query_pos:end_pos]


def get_sequence_before_junction(
    read: pysam.AlignedSegment,
    junction_cigar_index: int,
    length: int = 20,
) -> str:
    """
    Get the read sequence immediately before a junction (leftward in cigar / genomic order).

    For minus-strand reads stored in BAM, the 3' poly-A tail is reverse-complemented
    to poly-T and sits at the start of query_sequence (left side of the alignment).
    A false N near the 3' end thus has poly-T *before* it in cigar order, not after.

    Args:
        read: pysam AlignedSegment
        junction_cigar_index: Index of junction in CIGAR
        length: Number of bases to extract

    Returns:
        Sequence before the junction (up to `length` bases, ending at the junction)
    """
    if read.query_sequence is None or read.cigartuples is None:
        return ""

    # Accumulate query length up to (but not including) the N op
    query_pos = 0
    for i, (op, op_len) in enumerate(read.cigartuples):
        if i >= junction_cigar_index:
            break
        if op in (0, 1, 4, 7, 8):  # M, I, S, =, X consume query
            query_pos += op_len

    start_pos = max(0, query_pos - length)
    return read.query_sequence[start_pos:query_pos]


# =============================================================================
# Core Detection Functions
# =============================================================================

@dataclass
class JunctionAnalysis:
    """Results of junction artifact analysis."""
    junction_start: int
    junction_end: int
    junction_size: int
    distance_to_3prime: int
    downstream_a_richness: float
    target_a_richness: float
    has_canonical_motif: bool
    is_artifact: bool
    artifact_reason: str


def analyze_junction_for_artifact(
    junction: Tuple[int, int, int],
    read: pysam.AlignedSegment,
    reference: Optional[pysam.FastaFile],
    strand: str,
) -> JunctionAnalysis:
    """
    Analyze a single junction to determine if it's a poly(A) artifact.

    Args:
        junction: (start, end, cigar_index) tuple
        read: pysam AlignedSegment
        reference: Reference genome FastaFile (optional, for motif checking)
        strand: '+' or '-'

    Returns:
        JunctionAnalysis with artifact determination
    """
    start, end, cigar_idx = junction

    # Get 3' end position
    three_prime = get_read_3prime_position(read, strand)

    # Calculate distance to 3' end
    if strand == '+':
        dist_to_3prime = three_prime - end
    else:
        dist_to_3prime = start - three_prime

    # Initialize analysis
    analysis = JunctionAnalysis(
        junction_start=start,
        junction_end=end,
        junction_size=end - start,
        distance_to_3prime=dist_to_3prime,
        downstream_a_richness=0.0,
        target_a_richness=0.0,
        has_canonical_motif=False,
        is_artifact=False,
        artifact_reason="",
    )

    # Check if junction is too far from 3' end
    if dist_to_3prime > POLYA_JUNCTION_WINDOW:
        return analysis

    # Get the poly-A-side sequence: for plus strand it is after the N (rightward toward 3'),
    # for minus strand it is before the N (leftward toward 3').  In the BAM the minus-strand
    # poly-A is stored as poly-T at the start of query_sequence (reverse-complemented), so
    # it sits to the LEFT of the false N in cigar order.
    if strand == '+':
        polya_side_seq = get_sequence_after_junction(read, cigar_idx, 20)
        analysis.downstream_a_richness = calculate_a_richness(polya_side_seq)
    else:
        polya_side_seq = get_sequence_before_junction(read, cigar_idx, 20)
        analysis.downstream_a_richness = calculate_t_richness(polya_side_seq)

    # Check if downstream region is A-rich
    if analysis.downstream_a_richness < POLYA_RICHNESS_THRESHOLD:
        return analysis

    # Check junction target in genome (if reference available)
    if reference is not None:
        try:
            target_seq = reference.fetch(
                read.reference_name,
                end,
                min(end + 20, reference.get_reference_length(read.reference_name))
            )
            if strand == '+':
                analysis.target_a_richness = calculate_a_richness(target_seq)
            else:
                analysis.target_a_richness = calculate_t_richness(target_seq)

            # Check splice site motifs
            if start >= 2:
                donor_seq = reference.fetch(read.reference_name, start, start + 2)
                acceptor_seq = reference.fetch(read.reference_name, end - 2, end)

                # Canonical: GT-AG (plus strand) or CT-AC in genomic coords (minus strand)
                if strand == '+':
                    analysis.has_canonical_motif = (
                        donor_seq.upper() == 'GT' and acceptor_seq.upper() == 'AG'
                    )
                else:
                    # Minus strand: RC of GT-AG = CT-AC in genomic coordinates
                    # donor is at junction_start (left boundary), acceptor at junction_end-2
                    analysis.has_canonical_motif = (
                        donor_seq.upper() == 'CT' and acceptor_seq.upper() == 'AC'
                    )

        except Exception as e:
            logger.debug(f"Junction motif analysis failed for read '{read.query_name}': {e}")

    # Final artifact determination
    if analysis.downstream_a_richness >= POLYA_RICHNESS_THRESHOLD:
        if analysis.target_a_richness >= GENOMIC_ATRACT_THRESHOLD:
            analysis.is_artifact = True
            analysis.artifact_reason = "Junction to genomic A-tract with A-rich downstream"
        elif not analysis.has_canonical_motif:
            # Non-canonical junction near 3' end with A-rich downstream
            analysis.is_artifact = True
            analysis.artifact_reason = "Non-canonical junction near 3' end with A-rich downstream"

    return analysis


def filter_polya_artifact_junctions(
    read: pysam.AlignedSegment,
    reference: Optional[pysam.FastaFile] = None,
    strand: Optional[str] = None,
) -> Tuple[List[Tuple[int, int]], List[JunctionAnalysis]]:
    """
    Filter junctions, separating real from poly(A) artifacts.

    Args:
        read: pysam AlignedSegment
        reference: Reference genome FastaFile (optional)
        strand: Strand override, or None to infer from read

    Returns:
        (real_junctions, artifact_analyses)
        - real_junctions: List of (start, end) tuples for valid junctions
        - artifact_analyses: List of JunctionAnalysis for artifacts
    """
    if strand is None:
        strand = '-' if read.is_reverse else '+'

    # Extract all junctions
    all_junctions = extract_junctions_from_cigar(read)

    real_junctions = []
    artifact_analyses = []

    for junction in all_junctions:
        start, end, _ = junction

        analysis = analyze_junction_for_artifact(
            junction, read, reference, strand
        )

        if analysis.is_artifact:
            artifact_analyses.append(analysis)
        else:
            real_junctions.append((start, end))

    return real_junctions, artifact_analyses


def correct_3prime_for_artifact_junctions(
    read: pysam.AlignedSegment,
    artifacts: List[JunctionAnalysis],
    strand: str,
) -> Dict:
    """
    Correct the 3' end position accounting for artifact junctions.

    When a junction is identified as an artifact, the "true" 3' end
    should be at the upstream boundary of the artifact, not where
    the poly(A) was mapped.

    Args:
        read: pysam AlignedSegment
        artifacts: List of JunctionAnalysis for artifacts
        strand: '+' or '-'

    Returns:
        Dict with:
            'corrected_3prime': Corrected 3' end position
            'artifact_length': Total length absorbed by artifacts
            'n_artifacts': Number of artifact junctions
    """
    if not artifacts:
        three_prime = get_read_3prime_position(read, strand)
        return {
            'corrected_3prime': three_prime,
            'artifact_length': 0,
            'n_artifacts': 0,
        }

    # Get the artifact closest to 3' end
    artifacts_sorted = sorted(
        artifacts,
        key=lambda a: a.distance_to_3prime
    )
    closest_artifact = artifacts_sorted[0]

    # The corrected 3' end is the upstream boundary of the artifact
    if strand == '+':
        corrected_3prime = closest_artifact.junction_start - 1
        artifact_length = get_read_3prime_position(read, strand) - corrected_3prime
    else:
        corrected_3prime = closest_artifact.junction_end
        artifact_length = corrected_3prime - get_read_3prime_position(read, strand)

    return {
        'corrected_3prime': corrected_3prime,
        'artifact_length': artifact_length,
        'n_artifacts': len(artifacts),
    }


# =============================================================================
# Integration Functions
# =============================================================================

def process_read_junctions(
    read: pysam.AlignedSegment,
    reference: Optional[pysam.FastaFile] = None,
    strand: Optional[str] = None,
) -> Dict:
    """
    Process read junctions with artifact filtering.

    This is the main entry point for junction processing.

    Args:
        read: pysam AlignedSegment
        reference: Reference genome FastaFile (optional but recommended)
        strand: Strand override, or None to infer from read

    Returns:
        Dict with:
            'real_junctions': List of (start, end) for valid junctions
            'n_real_junctions': Count of valid junctions
            'artifact_junctions': List of JunctionAnalysis for artifacts
            'n_artifact_junctions': Count of artifact junctions
            'corrected_3prime': 3' end position after artifact correction
            'artifact_length': Total length absorbed by artifacts
    """
    if strand is None:
        strand = '-' if read.is_reverse else '+'

    # Filter junctions
    real_junctions, artifact_analyses = filter_polya_artifact_junctions(
        read, reference, strand
    )

    # Correct 3' end for artifacts
    correction = correct_3prime_for_artifact_junctions(
        read, artifact_analyses, strand
    )

    return {
        'real_junctions': real_junctions,
        'n_real_junctions': len(real_junctions),
        'artifact_junctions': artifact_analyses,
        'n_artifact_junctions': len(artifact_analyses),
        'corrected_3prime': correction['corrected_3prime'],
        'artifact_length': correction['artifact_length'],
    }


# =============================================================================
# Validation
# =============================================================================

def validate_junction_filter():
    """
    Print validation examples for junction filtering.
    """
    print("=" * 60)
    print("FALSE JUNCTION FILTER VALIDATION")
    print("=" * 60)

    # Create mock read for testing
    class MockRead:
        def __init__(self, start, end, is_reverse, cigartuples, sequence):
            self.reference_start = start
            self.reference_end = end
            self.is_reverse = is_reverse
            self.cigartuples = cigartuples
            self.query_sequence = sequence
            self.reference_name = 'chrI'

    # Test case 1: Real junction (internal, non-A-rich downstream)
    print("\nTest 1: Real junction (internal, non-A-rich downstream)")
    read1 = MockRead(
        start=1000,
        end=2100,
        is_reverse=False,
        cigartuples=[(0, 500), (3, 100), (0, 500)],  # 500M100N500M = 1000 query bases
        sequence="ATGCATGCATGC" * 84  # 84 * 12 = 1008 chars (>1000)
    )
    result1 = filter_polya_artifact_junctions(read1, None, '+')
    print(f"  Real junctions: {result1[0]}")
    print(f"  Artifact analyses: {len(result1[1])}")
    assert len(result1[0]) == 1, "Should have 1 real junction"
    assert len(result1[1]) == 0, "Should have 0 artifact junctions"
    print("  PASSED: Internal junction correctly classified as real")

    # Test case 2: Artifact junction (near 3' end, A-rich downstream)
    print("\nTest 2: Artifact junction (near 3' end, A-rich downstream)")
    # CIGAR: 980M100N20M = 1000 query bases
    # Junction at position 1000+980=1980 to 2080, with only 20bp after
    # Downstream sequence is A-rich
    read2 = MockRead(
        start=1000,
        end=2100,  # 1000 + 980 + 100 + 20 = 2100
        is_reverse=False,
        cigartuples=[(0, 980), (3, 100), (0, 20)],  # 980M100N20M
        sequence="ATGCATGCATGC" * 82 + "A" * 16  # 82*12=984, +16 = 1000 chars
    )
    result2 = filter_polya_artifact_junctions(read2, None, '+')
    print(f"  Real junctions: {result2[0]}")
    print(f"  Artifact analyses: {len(result2[1])}")
    if result2[1]:
        print(f"  Artifact reason: {result2[1][0].artifact_reason}")
        print(f"  Downstream A-richness: {result2[1][0].downstream_a_richness:.2f}")
        print("  PASSED: Near-3' junction with A-rich downstream classified as artifact")
    else:
        # Debug: print what the downstream sequence looks like
        junctions = extract_junctions_from_cigar(read2)
        if junctions:
            downstream = get_sequence_after_junction(read2, junctions[0][2], 20)
            print(f"  DEBUG: downstream sequence = '{downstream}'")
            print(f"  DEBUG: downstream A-richness = {calculate_a_richness(downstream):.2f}")
            # Check distance to 3' end
            three_prime = get_read_3prime_position(read2, '+')
            dist = three_prime - junctions[0][1]
            print(f"  DEBUG: distance to 3' end = {dist}")

    # Test case 3: Real junction (near 3' end but NOT A-rich downstream)
    print("\nTest 3: Real junction (near 3' end but NOT A-rich downstream)")
    read3 = MockRead(
        start=1000,
        end=2100,
        is_reverse=False,
        cigartuples=[(0, 980), (3, 100), (0, 20)],  # Same CIGAR as test 2
        sequence="ATGCATGCATGC" * 82 + "GCTAGCTAGCTAGCTA"  # 984 + 16 = 1000, NOT A-rich
    )
    result3 = filter_polya_artifact_junctions(read3, None, '+')
    print(f"  Real junctions: {result3[0]}")
    print(f"  Artifact analyses: {len(result3[1])}")
    if len(result3[0]) == 1 and len(result3[1]) == 0:
        print("  PASSED: Near-3' junction with non-A-rich downstream classified as real")
    else:
        junctions = extract_junctions_from_cigar(read3)
        if junctions:
            downstream = get_sequence_after_junction(read3, junctions[0][2], 20)
            print(f"  DEBUG: downstream sequence = '{downstream}'")
            print(f"  DEBUG: downstream A-richness = {calculate_a_richness(downstream):.2f}")

    print("\n" + "=" * 60)
    print("VALIDATION COMPLETE")
    print("=" * 60)


if __name__ == '__main__':
    validate_junction_filter()
