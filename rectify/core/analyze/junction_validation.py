#!/usr/bin/env python3
"""
Junction Validation Module for RECTIFY

Implements ESPRESSO-style multi-read junction validation to reduce false positive
novel junction calls. Novel junctions (not in annotation) require multiple
supporting reads before being considered valid.

Key Features:
- Multi-read evidence requirement (default: ≥2 reads)
- Canonical splice motif checking (GT-AG, GC-AG)
- Integration with existing UnifiedReadRecord structure

Reference:
    ESPRESSO: Robust discovery and quantification of transcript isoforms
    from error-prone long-read RNA-seq data (Science Advances, 2023)

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-25
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import numpy as np

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False

# Import UnifiedReadRecord for type hints
from rectify.core.unified_record import UnifiedReadRecord


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class JunctionEvidence:
    """
    Evidence supporting a splice junction.

    Attributes:
        junction: (donor, acceptor) genomic coordinates
        chrom: Chromosome name
        strand: Strand (+/-)
        supporting_reads: List of read IDs supporting this junction
        has_canonical_motif: True if GT-AG or GC-AG splice motif
        mean_mapping_quality: Average MAPQ of supporting reads
        is_annotated: True if junction is in annotation
    """
    junction: Tuple[int, int]
    chrom: str
    strand: str
    supporting_reads: List[str] = field(default_factory=list)
    has_canonical_motif: bool = True
    mean_mapping_quality: float = 0.0
    is_annotated: bool = False

    @property
    def donor(self) -> int:
        """Donor (5' splice site) position."""
        return self.junction[0]

    @property
    def acceptor(self) -> int:
        """Acceptor (3' splice site) position."""
        return self.junction[1]

    @property
    def n_supporting(self) -> int:
        """Number of reads supporting this junction."""
        return len(self.supporting_reads)

    @property
    def intron_length(self) -> int:
        """Length of the intron."""
        return self.acceptor - self.donor

    def is_validated(self, min_reads: int = 2) -> bool:
        """
        Check if junction passes validation criteria.

        ESPRESSO criterion: Novel junctions require ≥min_reads supporting reads.
        Annotated junctions are always valid.
        """
        if self.is_annotated:
            return True
        return self.n_supporting >= min_reads

    def to_dict(self) -> Dict:
        """Convert to dictionary for DataFrame construction."""
        return {
            'chrom': self.chrom,
            'donor': self.donor,
            'acceptor': self.acceptor,
            'strand': self.strand,
            'n_reads': self.n_supporting,
            'has_canonical_motif': self.has_canonical_motif,
            'mean_mapq': self.mean_mapping_quality,
            'is_annotated': self.is_annotated,
            'intron_length': self.intron_length,
        }


# =============================================================================
# Canonical Splice Motif Checking
# =============================================================================

CANONICAL_DONOR_MOTIFS = {'GT', 'GC'}  # 5' splice site
CANONICAL_ACCEPTOR_MOTIFS = {'AG'}      # 3' splice site


def check_canonical_splice_motif(
    junction: Tuple[int, int],
    genome: 'pysam.FastaFile',
    chrom: str,
    strand: str = '+',
) -> bool:
    """
    Check if junction has canonical splice motif.

    Canonical motifs:
    - GT-AG (most common, ~98% of introns)
    - GC-AG (less common, ~1% of introns)

    Args:
        junction: (donor, acceptor) positions
        genome: pysam.FastaFile for sequence lookup
        chrom: Chromosome name
        strand: Strand (+/-) for motif orientation

    Returns:
        True if junction has canonical splice motif
    """
    if not HAS_PYSAM:
        return True  # Assume canonical if can't check

    donor, acceptor = junction

    try:
        # Get 2bp at donor site (intron start)
        donor_seq = genome.fetch(chrom, donor, donor + 2).upper()
        # Get 2bp before acceptor site (intron end)
        acceptor_seq = genome.fetch(chrom, acceptor - 2, acceptor).upper()

        if strand == '-':
            # Reverse complement for minus strand.
            # On the minus strand the canonical donor dinucleotide is at the
            # donor genomic coordinate and the acceptor dinucleotide is at the
            # acceptor genomic coordinate; both must be reverse-complemented.
            donor_seq = _reverse_complement(
                genome.fetch(chrom, donor, donor + 2).upper()
            )
            acceptor_seq = _reverse_complement(
                genome.fetch(chrom, acceptor - 2, acceptor).upper()
            )

        return donor_seq in CANONICAL_DONOR_MOTIFS and acceptor_seq in CANONICAL_ACCEPTOR_MOTIFS

    except Exception:
        return True  # Assume canonical on error


def _reverse_complement(seq: str) -> str:
    """Get reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


# =============================================================================
# Main Validation Functions
# =============================================================================

def collect_junction_evidence(
    records: List[UnifiedReadRecord],
    known_junctions: Optional[Set[Tuple[str, int, int, str]]] = None,
) -> Dict[Tuple[str, int, int, str], JunctionEvidence]:
    """
    Collect evidence for all junctions from a set of read records.

    Groups reads by junction and collects supporting evidence including
    read IDs, mapping quality, and annotation status.

    Args:
        records: List of UnifiedReadRecord objects
        known_junctions: Set of (chrom, donor, acceptor, strand) tuples
            representing annotated junctions. If None, all junctions
            are treated as novel.

    Returns:
        Dictionary mapping (chrom, donor, acceptor, strand) to JunctionEvidence
    """
    if known_junctions is None:
        known_junctions = set()

    # Collect junction evidence
    junction_data = defaultdict(lambda: {
        'reads': [],
        'mapqs': [],
    })

    for record in records:
        if not record.junctions:
            continue

        for donor, acceptor in record.junctions:
            key = (record.chrom, donor, acceptor, record.strand)
            junction_data[key]['reads'].append(record.read_id)
            junction_data[key]['mapqs'].append(record.mapq)

    # Create JunctionEvidence objects
    evidence = {}
    for (chrom, donor, acceptor, strand), data in junction_data.items():
        is_annotated = (chrom, donor, acceptor, strand) in known_junctions

        evidence[(chrom, donor, acceptor, strand)] = JunctionEvidence(
            junction=(donor, acceptor),
            chrom=chrom,
            strand=strand,
            supporting_reads=data['reads'],
            mean_mapping_quality=np.mean(data['mapqs']) if data['mapqs'] else 0.0,
            is_annotated=is_annotated,
        )

    return evidence


def validate_novel_junctions(
    records: List[UnifiedReadRecord],
    known_junctions: Optional[Set[Tuple[str, int, int, str]]] = None,
    min_supporting_reads: int = 2,
    require_canonical_motif: bool = True,
    genome: Optional['pysam.FastaFile'] = None,
) -> Dict[Tuple[str, int, int, str], JunctionEvidence]:
    """
    Validate novel junctions using ESPRESSO-style multi-read evidence.

    Novel junctions (not in annotation) require:
    1. ≥min_supporting_reads independent reads
    2. Optionally: canonical splice motif (GT-AG or GC-AG)

    Args:
        records: List of UnifiedReadRecord objects
        known_junctions: Set of annotated junctions as (chrom, donor, acceptor, strand)
        min_supporting_reads: Minimum reads required for novel junctions
        require_canonical_motif: If True, novel junctions must have canonical motif
        genome: pysam.FastaFile for motif checking (optional)

    Returns:
        Dictionary of validated junctions with their evidence
    """
    if known_junctions is None:
        known_junctions = set()

    # Collect all junction evidence
    all_evidence = collect_junction_evidence(records, known_junctions)

    # Check canonical motifs if requested and genome available
    if require_canonical_motif and genome is not None:
        for key, evidence in all_evidence.items():
            chrom, donor, acceptor, strand = key
            evidence.has_canonical_motif = check_canonical_splice_motif(
                (donor, acceptor), genome, chrom, strand
            )

    # Filter to validated junctions
    validated = {}
    for key, evidence in all_evidence.items():
        # Annotated junctions are always valid
        if evidence.is_annotated:
            validated[key] = evidence
            continue

        # Novel junctions need sufficient support
        if evidence.n_supporting < min_supporting_reads:
            continue

        # Check canonical motif if required
        if require_canonical_motif and not evidence.has_canonical_motif:
            continue

        validated[key] = evidence

    return validated


def filter_records_to_validated_junctions(
    records: List[UnifiedReadRecord],
    validated_junctions: Dict[Tuple[str, int, int, str], JunctionEvidence],
) -> List[UnifiedReadRecord]:
    """
    Filter read records to only include validated junctions.

    Removes junctions that failed validation from each record.
    Records with no remaining junctions are still included.

    Args:
        records: List of UnifiedReadRecord objects
        validated_junctions: Dictionary of validated junctions

    Returns:
        List of UnifiedReadRecord with filtered junctions
    """
    filtered_records = []

    for record in records:
        if not record.junctions:
            filtered_records.append(record)
            continue

        # Filter junctions to only validated ones
        valid_junctions = []
        for donor, acceptor in record.junctions:
            key = (record.chrom, donor, acceptor, record.strand)
            if key in validated_junctions:
                valid_junctions.append((donor, acceptor))

        # Create new record with filtered junctions
        # Using dataclass replace pattern
        filtered_record = UnifiedReadRecord(
            read_id=record.read_id,
            chrom=record.chrom,
            strand=record.strand,
            five_prime_raw=record.five_prime_raw,
            five_prime_corrected=record.five_prime_corrected,
            first_exon_start=record.first_exon_start,
            starts_in_intron=record.starts_in_intron,
            three_prime_raw=record.three_prime_raw,
            three_prime_corrected=record.three_prime_corrected,
            alignment_start=record.alignment_start,
            alignment_end=record.alignment_end,
            junctions=valid_junctions,
            n_junctions=len(valid_junctions),
            junctions_filtered=record.n_junctions - len(valid_junctions),
            five_prime_soft_clip_length=record.five_prime_soft_clip_length,
            three_prime_soft_clip_length=record.three_prime_soft_clip_length,
            five_prime_soft_clip_seq=record.five_prime_soft_clip_seq,
            three_prime_soft_clip_seq=record.three_prime_soft_clip_seq,
            polya_length=record.polya_length,
            aligned_a_length=record.aligned_a_length,
            soft_clip_a_length=record.soft_clip_a_length,
            mapq=record.mapq,
            cigar_summary=record.cigar_summary,
            alignment_identity=record.alignment_identity,
            best_aligner=record.best_aligner,
            consensus_confidence=record.consensus_confidence,
            n_aligners_support=record.n_aligners_support,
            count=record.count,
        )
        filtered_records.append(filtered_record)

    return filtered_records


# =============================================================================
# Annotation Loading
# =============================================================================

def load_known_junctions_from_gff(
    gff_path: str,
    feature_type: str = 'intron',
) -> Set[Tuple[str, int, int, str]]:
    """
    Load known junctions from a GFF/GTF annotation file.

    Args:
        gff_path: Path to GFF/GTF file
        feature_type: Feature type to extract junctions from

    Returns:
        Set of (chrom, donor, acceptor, strand) tuples
    """
    known_junctions = set()

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

                known_junctions.add((chrom, start, end, strand))

    return known_junctions


def load_known_junctions_from_bed(
    bed_path: str,
) -> Set[Tuple[str, int, int, str]]:
    """
    Load known junctions from a BED file (e.g., from minimap2 --junc-bed).

    Expects BED6 format: chrom, start, end, name, score, strand

    Args:
        bed_path: Path to BED file

    Returns:
        Set of (chrom, donor, acceptor, strand) tuples
    """
    known_junctions = set()

    with open(bed_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue

            chrom = parts[0]
            start = int(parts[1])  # Already 0-based
            end = int(parts[2])
            strand = parts[5] if len(parts) > 5 else '+'

            known_junctions.add((chrom, start, end, strand))

    return known_junctions


# =============================================================================
# Summary Statistics
# =============================================================================

def summarize_junction_validation(
    all_evidence: Dict[Tuple[str, int, int, str], JunctionEvidence],
    validated: Dict[Tuple[str, int, int, str], JunctionEvidence],
) -> Dict:
    """
    Generate summary statistics for junction validation.

    Args:
        all_evidence: All junction evidence before filtering
        validated: Validated junctions after filtering

    Returns:
        Dictionary with summary statistics
    """
    n_total = len(all_evidence)
    n_validated = len(validated)
    n_filtered = n_total - n_validated

    # Count by annotation status
    n_annotated_total = sum(1 for e in all_evidence.values() if e.is_annotated)
    n_novel_total = n_total - n_annotated_total

    n_annotated_validated = sum(1 for e in validated.values() if e.is_annotated)
    n_novel_validated = n_validated - n_annotated_validated

    # Count canonical motifs in validated
    n_canonical = sum(1 for e in validated.values() if e.has_canonical_motif)

    # Calculate statistics
    return {
        'total_junctions': n_total,
        'validated_junctions': n_validated,
        'filtered_junctions': n_filtered,
        'annotated_total': n_annotated_total,
        'novel_total': n_novel_total,
        'annotated_validated': n_annotated_validated,
        'novel_validated': n_novel_validated,
        'novel_filtered': n_novel_total - n_novel_validated,
        'canonical_motif_count': n_canonical,
        'validation_rate': n_validated / n_total if n_total > 0 else 0,
        'novel_validation_rate': n_novel_validated / n_novel_total if n_novel_total > 0 else 0,
    }


# =============================================================================
# DataFrame Export
# =============================================================================

def evidence_to_dataframe(
    evidence: Dict[Tuple[str, int, int, str], JunctionEvidence],
) -> 'pd.DataFrame':
    """
    Convert junction evidence to pandas DataFrame.

    Args:
        evidence: Dictionary of junction evidence

    Returns:
        DataFrame with junction information
    """
    import pandas as pd

    rows = [e.to_dict() for e in evidence.values()]
    if not rows:
        return pd.DataFrame(columns=[
            'chrom', 'donor', 'acceptor', 'strand', 'n_reads',
            'has_canonical_motif', 'mean_mapq', 'is_annotated', 'intron_length'
        ])

    df = pd.DataFrame(rows)
    return df.sort_values(['chrom', 'donor']).reset_index(drop=True)


if __name__ == '__main__':
    # Quick test
    print("Testing junction_validation module...")

    # Create test records
    test_records = [
        UnifiedReadRecord(
            read_id="read_001",
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
            junctions=[(1200, 1400), (2000, 2200)],  # Two junctions
            n_junctions=2,
            mapq=60,
        ),
        UnifiedReadRecord(
            read_id="read_002",
            chrom="chrI",
            strand="+",
            five_prime_raw=1050,
            five_prime_corrected=1050,
            first_exon_start=1050,
            starts_in_intron=False,
            three_prime_raw=2950,
            three_prime_corrected=2945,
            alignment_start=1050,
            alignment_end=2950,
            junctions=[(1200, 1400)],  # Only first junction
            n_junctions=1,
            mapq=55,
        ),
        UnifiedReadRecord(
            read_id="read_003",
            chrom="chrI",
            strand="+",
            five_prime_raw=1100,
            five_prime_corrected=1100,
            first_exon_start=1100,
            starts_in_intron=False,
            three_prime_raw=2900,
            three_prime_corrected=2895,
            alignment_start=1100,
            alignment_end=2900,
            junctions=[(1200, 1400), (2500, 2700)],  # First + novel junction
            n_junctions=2,
            mapq=50,
        ),
    ]

    # Collect evidence
    evidence = collect_junction_evidence(test_records)
    print(f"\nTotal junctions found: {len(evidence)}")

    for key, e in evidence.items():
        print(f"  {key}: {e.n_supporting} reads")

    # Validate with min_reads=2
    validated = validate_novel_junctions(
        test_records,
        min_supporting_reads=2,
        require_canonical_motif=False,  # No genome for test
    )
    print(f"\nValidated junctions (≥2 reads): {len(validated)}")

    # Summarize
    summary = summarize_junction_validation(evidence, validated)
    print(f"\nSummary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")

    print("\njunction_validation module ready!")
