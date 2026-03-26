#!/usr/bin/env python3
"""
Full-Length Read Classifier for RECTIFY

Implements Bambu-style classification to distinguish full-length reads from
5'-truncated reads in nanopore Direct RNA Sequencing (DRS) data.

DRS Context:
- Sequencing starts at the 3' end (polyA tail) and extends 5'-ward
- Many reads are truncated before reaching the true transcription start site
- Full-length reads provide more confident 5' end positions

Key Features:
- Heuristic classification based on soft-clips and alignment fraction
- Optional annotation-based distance-to-TSS checking
- Confidence weighting for downstream 5' end analysis

Reference:
    Bambu: Context-aware transcript quantification from long-read RNA-seq data
    (Nature Methods, 2023)

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-25
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import numpy as np

# Import UnifiedReadRecord for type hints
from rectify.core.unified_record import UnifiedReadRecord


# Default thresholds
DEFAULT_MAX_5PRIME_SOFTCLIP = 50      # bp
DEFAULT_MIN_ALIGNMENT_FRACTION = 0.8   # fraction of read aligned
DEFAULT_MAX_TSS_DISTANCE = 500         # bp from annotated TSS


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class FullLengthClassification:
    """
    Classification of read completeness (full-length vs truncated).

    Attributes:
        read_id: Unique read identifier
        is_full_length: True if read is classified as full-length
        confidence: Confidence score 0-1 (1 = high confidence full-length)
        truncation_reasons: List of reasons for truncation classification
        five_prime_softclip: Length of 5' soft-clip
        alignment_fraction: Fraction of read that aligned
        distance_to_tss: Distance to nearest annotated TSS (if available)
    """
    read_id: str
    is_full_length: bool
    confidence: float
    truncation_reasons: List[str] = field(default_factory=list)
    five_prime_softclip: int = 0
    alignment_fraction: float = 1.0
    distance_to_tss: Optional[int] = None

    def to_dict(self) -> Dict:
        """Convert to dictionary for DataFrame construction."""
        return {
            'read_id': self.read_id,
            'is_full_length': self.is_full_length,
            'confidence': self.confidence,
            'truncation_reason': '; '.join(self.truncation_reasons) if self.truncation_reasons else None,
            'five_prime_softclip': self.five_prime_softclip,
            'alignment_fraction': self.alignment_fraction,
            'distance_to_tss': self.distance_to_tss,
        }


# =============================================================================
# Classification Functions
# =============================================================================

def classify_full_length_heuristic(
    record: UnifiedReadRecord,
    max_5prime_softclip: int = DEFAULT_MAX_5PRIME_SOFTCLIP,
    min_alignment_fraction: float = DEFAULT_MIN_ALIGNMENT_FRACTION,
    max_tss_distance: int = DEFAULT_MAX_TSS_DISTANCE,
    tss_positions: Optional[Dict[str, List[int]]] = None,
) -> FullLengthClassification:
    """
    Heuristically classify a read as full-length or truncated.

    Indicators of truncation:
    1. Large 5' soft-clip (sequence extends beyond alignment)
    2. Low alignment fraction (much of read didn't align)
    3. 5' end far from annotated TSS (if annotation provided)

    DRS context: Reads start at 3' (polyA) and extend 5'-ward.
    Truncated reads fail to reach the true 5' end (TSS).

    Args:
        record: UnifiedReadRecord to classify
        max_5prime_softclip: Maximum 5' soft-clip for full-length
        min_alignment_fraction: Minimum alignment fraction for full-length
        max_tss_distance: Maximum distance from annotated TSS
        tss_positions: Dict mapping chrom to list of TSS positions

    Returns:
        FullLengthClassification with assessment
    """
    reasons = []
    confidence_penalties = []

    # Calculate values
    five_prime_softclip = record.five_prime_soft_clip_length

    # Calculate alignment fraction
    query_length = record.alignment_end - record.alignment_start
    if record.five_prime_soft_clip_length or record.three_prime_soft_clip_length:
        total_length = (
            query_length +
            record.five_prime_soft_clip_length +
            record.three_prime_soft_clip_length
        )
    else:
        total_length = query_length

    alignment_fraction = query_length / total_length if total_length > 0 else 1.0

    # Check 5' soft-clip length
    if five_prime_softclip > max_5prime_softclip:
        reasons.append(f"5prime_softclip_{five_prime_softclip}bp")
        # Penalty proportional to excess
        penalty = min(0.4, (five_prime_softclip - max_5prime_softclip) / 100)
        confidence_penalties.append(penalty)

    # Check alignment fraction
    if alignment_fraction < min_alignment_fraction:
        reasons.append(f"low_alignment_{alignment_fraction:.2f}")
        # Penalty proportional to deficit
        penalty = min(0.4, (min_alignment_fraction - alignment_fraction) / 0.3)
        confidence_penalties.append(penalty)

    # Check distance to annotated TSS (if available)
    distance_to_tss = None
    if tss_positions and record.chrom in tss_positions:
        # Get 5' end position
        if record.strand == '+':
            five_prime_pos = record.five_prime_corrected
        else:
            five_prime_pos = record.five_prime_corrected

        # Find nearest TSS
        tss_list = tss_positions[record.chrom]
        if tss_list:
            distances = [abs(five_prime_pos - tss) for tss in tss_list]
            distance_to_tss = min(distances)

            if distance_to_tss > max_tss_distance:
                reasons.append(f"far_from_tss_{distance_to_tss}bp")
                # Penalty proportional to excess distance
                penalty = min(0.3, (distance_to_tss - max_tss_distance) / 1000)
                confidence_penalties.append(penalty)

    # Compute final classification
    is_full_length = len(reasons) == 0
    confidence = max(0.0, 1.0 - sum(confidence_penalties))

    return FullLengthClassification(
        read_id=record.read_id,
        is_full_length=is_full_length,
        confidence=confidence,
        truncation_reasons=reasons,
        five_prime_softclip=five_prime_softclip,
        alignment_fraction=alignment_fraction,
        distance_to_tss=distance_to_tss,
    )


def classify_all_reads(
    records: List[UnifiedReadRecord],
    max_5prime_softclip: int = DEFAULT_MAX_5PRIME_SOFTCLIP,
    min_alignment_fraction: float = DEFAULT_MIN_ALIGNMENT_FRACTION,
    max_tss_distance: int = DEFAULT_MAX_TSS_DISTANCE,
    tss_positions: Optional[Dict[str, List[int]]] = None,
) -> Dict[str, FullLengthClassification]:
    """
    Classify all reads as full-length or truncated.

    Args:
        records: List of UnifiedReadRecord objects
        max_5prime_softclip: Maximum 5' soft-clip for full-length
        min_alignment_fraction: Minimum alignment fraction for full-length
        max_tss_distance: Maximum distance from annotated TSS
        tss_positions: Dict mapping chrom to list of TSS positions

    Returns:
        Dictionary mapping read_id to FullLengthClassification
    """
    classifications = {}

    for record in records:
        classification = classify_full_length_heuristic(
            record,
            max_5prime_softclip=max_5prime_softclip,
            min_alignment_fraction=min_alignment_fraction,
            max_tss_distance=max_tss_distance,
            tss_positions=tss_positions,
        )
        classifications[record.read_id] = classification

    return classifications


# =============================================================================
# 5' End Weighting
# =============================================================================

def weight_5prime_ends(
    records: List[UnifiedReadRecord],
    classifications: Dict[str, FullLengthClassification],
) -> Dict[Tuple[str, str, int], float]:
    """
    Weight 5' end positions by full-length confidence.

    Full-length reads contribute weight=1.0 to their 5' position.
    Truncated reads contribute weight=confidence (0-1).

    Args:
        records: List of UnifiedReadRecord objects
        classifications: Dictionary of read classifications

    Returns:
        Dictionary mapping (chrom, strand, position) to total weight
    """
    position_weights = defaultdict(float)

    for record in records:
        # Get 5' end position
        five_prime_pos = record.five_prime_corrected

        # Get classification
        classification = classifications.get(record.read_id)

        if classification:
            if classification.is_full_length:
                weight = 1.0
            else:
                weight = classification.confidence
        else:
            weight = 1.0  # Default if not classified

        key = (record.chrom, record.strand, five_prime_pos)
        position_weights[key] += weight * record.count

    return dict(position_weights)


def get_weighted_5prime_distribution(
    records: List[UnifiedReadRecord],
    classifications: Dict[str, FullLengthClassification],
    gene_start: int,
    gene_end: int,
    chrom: str,
    strand: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get weighted 5' end distribution for a genomic region.

    Args:
        records: List of UnifiedReadRecord objects
        classifications: Dictionary of read classifications
        gene_start: Start of region
        gene_end: End of region
        chrom: Chromosome
        strand: Strand

    Returns:
        Tuple of (positions array, weights array)
    """
    positions = []
    weights = []

    for record in records:
        if record.chrom != chrom or record.strand != strand:
            continue

        five_prime_pos = record.five_prime_corrected

        if gene_start <= five_prime_pos <= gene_end:
            classification = classifications.get(record.read_id)

            if classification:
                weight = 1.0 if classification.is_full_length else classification.confidence
            else:
                weight = 1.0

            positions.append(five_prime_pos)
            weights.append(weight * record.count)

    return np.array(positions), np.array(weights)


# =============================================================================
# TSS Annotation Loading
# =============================================================================

def load_tss_from_gff(
    gff_path: str,
    feature_type: str = 'gene',
) -> Dict[str, List[int]]:
    """
    Load transcription start site (TSS) positions from GFF annotation.

    TSS is the 5' end of the gene:
    - Plus strand: gene start (leftmost coordinate)
    - Minus strand: gene end (rightmost coordinate)

    Args:
        gff_path: Path to GFF/GTF file
        feature_type: Feature type to extract TSS from

    Returns:
        Dictionary mapping chrom to list of TSS positions
    """
    tss_positions = defaultdict(list)

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            if parts[2] == feature_type:
                chrom = parts[0]
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])
                strand = parts[6]

                # TSS is 5' end
                if strand == '+':
                    tss = start
                else:
                    tss = end - 1

                tss_positions[chrom].append(tss)

    # Sort positions for efficient lookup
    for chrom in tss_positions:
        tss_positions[chrom].sort()

    return dict(tss_positions)


# =============================================================================
# Summary Statistics
# =============================================================================

def summarize_full_length_classification(
    classifications: Dict[str, FullLengthClassification],
) -> Dict:
    """
    Generate summary statistics for full-length classification.

    Args:
        classifications: Dictionary of read classifications

    Returns:
        Dictionary with summary statistics
    """
    n_total = len(classifications)
    n_full_length = sum(1 for c in classifications.values() if c.is_full_length)
    n_truncated = n_total - n_full_length

    confidences = [c.confidence for c in classifications.values()]
    softclips = [c.five_prime_softclip for c in classifications.values()]
    alignment_fractions = [c.alignment_fraction for c in classifications.values()]

    # Count truncation reasons
    reason_counts = defaultdict(int)
    for c in classifications.values():
        for reason in c.truncation_reasons:
            # Extract reason type (before the underscore with value)
            reason_type = reason.split('_')[0] + '_' + reason.split('_')[1]
            reason_counts[reason_type] += 1

    return {
        'total_reads': n_total,
        'full_length_reads': n_full_length,
        'truncated_reads': n_truncated,
        'full_length_fraction': n_full_length / n_total if n_total > 0 else 0,
        'mean_confidence': np.mean(confidences) if confidences else 0,
        'median_confidence': np.median(confidences) if confidences else 0,
        'mean_5prime_softclip': np.mean(softclips) if softclips else 0,
        'median_5prime_softclip': np.median(softclips) if softclips else 0,
        'mean_alignment_fraction': np.mean(alignment_fractions) if alignment_fractions else 0,
        'truncation_reasons': dict(reason_counts),
    }


# =============================================================================
# DataFrame Export
# =============================================================================

def classifications_to_dataframe(
    classifications: Dict[str, FullLengthClassification],
) -> 'pd.DataFrame':
    """
    Convert classifications to pandas DataFrame.

    Args:
        classifications: Dictionary of read classifications

    Returns:
        DataFrame with classification information
    """
    import pandas as pd

    rows = [c.to_dict() for c in classifications.values()]
    if not rows:
        return pd.DataFrame(columns=[
            'read_id', 'is_full_length', 'confidence',
            'truncation_reason', 'five_prime_softclip',
            'alignment_fraction', 'distance_to_tss'
        ])

    return pd.DataFrame(rows)


if __name__ == '__main__':
    # Quick test
    print("Testing full_length_classifier module...")

    # Create test records
    test_records = [
        UnifiedReadRecord(
            read_id="full_length_001",
            chrom="chrI",
            strand="+",
            five_prime_raw=1000,
            five_prime_corrected=1000,
            first_exon_start=1000,
            starts_in_intron=False,
            three_prime_raw=3000,
            three_prime_corrected=2995,
            alignment_start=1000,
            alignment_end=3000,
            five_prime_soft_clip_length=10,  # Small soft-clip
            three_prime_soft_clip_length=5,
            mapq=60,
        ),
        UnifiedReadRecord(
            read_id="truncated_001",
            chrom="chrI",
            strand="+",
            five_prime_raw=1500,
            five_prime_corrected=1500,
            first_exon_start=1500,
            starts_in_intron=False,
            three_prime_raw=3000,
            three_prime_corrected=2995,
            alignment_start=1500,
            alignment_end=3000,
            five_prime_soft_clip_length=100,  # Large soft-clip
            three_prime_soft_clip_length=5,
            mapq=50,
        ),
        UnifiedReadRecord(
            read_id="truncated_002",
            chrom="chrI",
            strand="+",
            five_prime_raw=2000,
            five_prime_corrected=2000,
            first_exon_start=2000,
            starts_in_intron=False,
            three_prime_raw=3000,
            three_prime_corrected=2995,
            alignment_start=2000,
            alignment_end=3000,
            five_prime_soft_clip_length=20,
            three_prime_soft_clip_length=200,  # Low alignment fraction
            mapq=45,
        ),
    ]

    # Classify
    classifications = classify_all_reads(test_records)

    print("\nClassifications:")
    for read_id, classification in classifications.items():
        status = "Full-length" if classification.is_full_length else "Truncated"
        print(f"  {read_id}: {status} (confidence: {classification.confidence:.2f})")
        if classification.truncation_reasons:
            print(f"    Reasons: {', '.join(classification.truncation_reasons)}")

    # Summary
    summary = summarize_full_length_classification(classifications)
    print(f"\nSummary:")
    print(f"  Full-length: {summary['full_length_reads']}/{summary['total_reads']}")
    print(f"  Mean confidence: {summary['mean_confidence']:.2f}")

    # Weight 5' ends
    weights = weight_5prime_ends(test_records, classifications)
    print(f"\n5' end weights:")
    for (chrom, strand, pos), weight in sorted(weights.items()):
        print(f"  {chrom}:{strand}:{pos} -> weight={weight:.2f}")

    print("\nfull_length_classifier module ready!")
