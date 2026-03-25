"""
Terminal exon refinement for soft-clipped or misaligned read ends.

Many direct RNA-seq reads have soft-clipped 5' ends that represent
unmapped portions of the first exon. This module attempts to:

1. Detect soft-clipped ends near potential splice sites
2. Search for the soft-clipped sequence in nearby upstream/downstream exons
3. Refine the alignment by inserting appropriate splice junctions

This is applied BEFORE COMPASS scoring to improve junction detection.

Design principle:
- Uses known splice sites from annotation to guide refinement
- Does NOT bias scoring - just improves alignment accuracy
- Falls back to soft-clip if no good alignment found

Author: Kevin R. Roy
"""

import logging
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
from pathlib import Path

import pysam

logger = logging.getLogger(__name__)


# Splice site motifs
CANONICAL_5SS = 'GT'  # 5' splice site (donor) - intron start
CANONICAL_3SS = 'AG'  # 3' splice site (acceptor) - intron end

# Minimum soft-clip length to attempt refinement
MIN_SOFTCLIP_LENGTH = 10

# Maximum distance to search for splice site
MAX_SPLICE_SITE_DISTANCE = 50

# Minimum alignment score for refined exon (as fraction of perfect match)
MIN_ALIGNMENT_SCORE_FRACTION = 0.7


@dataclass
class SpliceSite:
    """Represents a known splice site from annotation."""
    chrom: str
    position: int  # 0-based position of GT (5'SS) or AG (3'SS)
    strand: str
    site_type: str  # '5ss' or '3ss'
    gene_name: str = ''
    intron_id: str = ''


@dataclass
class RefinementResult:
    """Result of attempting to refine a soft-clipped end."""
    read_id: str
    original_cigar: str
    refined_cigar: Optional[str] = None
    refined_start: Optional[int] = None
    junction_added: Optional[Tuple[int, int]] = None  # (intron_start, intron_end)
    soft_clip_length: int = 0
    soft_clip_sequence: str = ''
    matched_splice_site: Optional[SpliceSite] = None
    alignment_score: float = 0.0
    success: bool = False
    reason: str = ''


@dataclass
class SpliceSiteIndex:
    """Index of splice sites for fast lookup."""
    # Indexed by chrom -> position -> SpliceSite
    five_ss: Dict[str, Dict[int, SpliceSite]] = field(default_factory=dict)
    three_ss: Dict[str, Dict[int, SpliceSite]] = field(default_factory=dict)

    def add_site(self, site: SpliceSite):
        """Add a splice site to the index."""
        if site.site_type == '5ss':
            if site.chrom not in self.five_ss:
                self.five_ss[site.chrom] = {}
            self.five_ss[site.chrom][site.position] = site
        else:
            if site.chrom not in self.three_ss:
                self.three_ss[site.chrom] = {}
            self.three_ss[site.chrom][site.position] = site

    def get_nearby_5ss(self, chrom: str, position: int, max_distance: int = 50) -> List[SpliceSite]:
        """Get 5' splice sites within distance of position."""
        sites = []
        if chrom in self.five_ss:
            for pos, site in self.five_ss[chrom].items():
                if abs(pos - position) <= max_distance:
                    sites.append(site)
        return sorted(sites, key=lambda s: abs(s.position - position))

    def get_nearby_3ss(self, chrom: str, position: int, max_distance: int = 50) -> List[SpliceSite]:
        """Get 3' splice sites within distance of position."""
        sites = []
        if chrom in self.three_ss:
            for pos, site in self.three_ss[chrom].items():
                if abs(pos - position) <= max_distance:
                    sites.append(site)
        return sorted(sites, key=lambda s: abs(s.position - position))


def extract_junctions_from_bam(bam_path: str, min_reads: int = 1) -> SpliceSiteIndex:
    """Extract splice junctions from BAM file N operations.

    This captures NOVEL junctions detected by aligners that may not be
    in the annotation. Combined with annotated junctions, this provides
    the most complete picture for terminal exon refinement.

    Args:
        bam_path: Path to BAM file
        min_reads: Minimum reads supporting a junction

    Returns:
        SpliceSiteIndex with detected 5'SS and 3'SS positions
    """
    # Count junction occurrences
    junction_counts: Dict[Tuple[str, int, int, str], int] = {}

    bam = pysam.AlignmentFile(bam_path, 'rb')

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        if not read.cigartuples:
            continue

        strand = '-' if read.is_reverse else '+'
        chrom = read.reference_name

        # Extract N operations (splice junctions)
        ref_pos = read.reference_start
        for op, length in read.cigartuples:
            if op == 3:  # N = skipped region (intron)
                intron_start = ref_pos
                intron_end = ref_pos + length
                key = (chrom, intron_start, intron_end, strand)
                junction_counts[key] = junction_counts.get(key, 0) + 1

            # Update reference position
            if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
                ref_pos += length

    bam.close()

    # Build splice site index from junctions with sufficient support
    index = SpliceSiteIndex()

    for (chrom, intron_start, intron_end, strand), count in junction_counts.items():
        if count < min_reads:
            continue

        # Create splice site entries
        # For both strands, intron_start is always the 5'SS genomic position
        # and intron_end-2 is the 3'SS genomic position (AG dinucleotide start)
        intron_id = f"novel_{chrom}_{intron_start}_{intron_end}"

        if strand == '+':
            # + strand: 5'SS at intron start (GT), 3'SS at intron end (AG)
            five_ss = SpliceSite(
                chrom=chrom, position=intron_start, strand=strand,
                site_type='5ss', intron_id=intron_id
            )
            three_ss = SpliceSite(
                chrom=chrom, position=intron_end - 2, strand=strand,
                site_type='3ss', intron_id=intron_id
            )
        else:
            # - strand: genomic coordinates are reversed relative to RNA
            # But in CIGAR, N operations are always in genomic left-to-right order
            # So intron_start is still the leftmost position
            five_ss = SpliceSite(
                chrom=chrom, position=intron_end - 2, strand=strand,
                site_type='5ss', intron_id=intron_id
            )
            three_ss = SpliceSite(
                chrom=chrom, position=intron_start, strand=strand,
                site_type='3ss', intron_id=intron_id
            )

        index.add_site(five_ss)
        index.add_site(three_ss)

    n_junctions = len(junction_counts)
    n_above_threshold = sum(1 for c in junction_counts.values() if c >= min_reads)
    logger.info(f"Extracted {n_junctions} junctions from BAM, {n_above_threshold} with >= {min_reads} reads")

    return index


def merge_splice_indices(*indices: SpliceSiteIndex) -> SpliceSiteIndex:
    """Merge multiple splice site indices into one.

    When the same position appears in multiple indices, the first
    occurrence is kept (annotation takes precedence over novel).

    Args:
        *indices: Variable number of SpliceSiteIndex objects

    Returns:
        Merged SpliceSiteIndex
    """
    merged = SpliceSiteIndex()

    for index in indices:
        for chrom, sites in index.five_ss.items():
            for pos, site in sites.items():
                if chrom not in merged.five_ss:
                    merged.five_ss[chrom] = {}
                if pos not in merged.five_ss[chrom]:
                    merged.five_ss[chrom][pos] = site

        for chrom, sites in index.three_ss.items():
            for pos, site in sites.items():
                if chrom not in merged.three_ss:
                    merged.three_ss[chrom] = {}
                if pos not in merged.three_ss[chrom]:
                    merged.three_ss[chrom][pos] = site

    n_5ss = sum(len(v) for v in merged.five_ss.values())
    n_3ss = sum(len(v) for v in merged.three_ss.values())
    logger.info(f"Merged index has {n_5ss} 5'SS and {n_3ss} 3'SS")

    return merged


def load_splice_sites_from_gff(gff_path: str) -> SpliceSiteIndex:
    """Load splice sites from GFF intron features.

    Args:
        gff_path: Path to GFF annotation file

    Returns:
        SpliceSiteIndex with 5'SS and 3'SS positions
    """
    index = SpliceSiteIndex()

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attrs = fields

            if 'intron' not in feature.lower():
                continue

            # Parse attributes for name
            name = ''
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('Name='):
                    name = attr.split('=')[1]
                    break

            # Convert to 0-based
            start_0 = int(start) - 1
            end_0 = int(end)  # Half-open

            # For + strand: 5'SS is at intron start, 3'SS is at intron end
            # For - strand: 5'SS is at intron end, 3'SS is at intron start
            if strand == '+':
                five_ss = SpliceSite(
                    chrom=chrom, position=start_0, strand=strand,
                    site_type='5ss', intron_id=name
                )
                three_ss = SpliceSite(
                    chrom=chrom, position=end_0 - 2, strand=strand,  # AG is 2bp before end
                    site_type='3ss', intron_id=name
                )
            else:
                # Minus strand: intron coordinates are still left-to-right in genome
                # but 5'SS (GT) is at the right (end) and 3'SS (AG) is at left (start)
                five_ss = SpliceSite(
                    chrom=chrom, position=end_0 - 2, strand=strand,
                    site_type='5ss', intron_id=name
                )
                three_ss = SpliceSite(
                    chrom=chrom, position=start_0, strand=strand,
                    site_type='3ss', intron_id=name
                )

            index.add_site(five_ss)
            index.add_site(three_ss)

    n_5ss = sum(len(v) for v in index.five_ss.values())
    n_3ss = sum(len(v) for v in index.three_ss.values())
    logger.info(f"Loaded {n_5ss} 5'SS and {n_3ss} 3'SS from {gff_path}")

    return index


@dataclass
class MismatchCluster:
    """A cluster of mismatches/indels that may indicate a missed splice junction."""
    chrom: str
    start: int  # Genomic start of problematic region
    end: int    # Genomic end of problematic region
    strand: str
    n_mismatches: int
    n_insertions: int
    n_deletions: int
    total_events: int
    density: float  # Events per bp
    read_id: str
    position_in_read: int  # Query position


def detect_mismatch_clusters(
    read: pysam.AlignedSegment,
    window_size: int = 15,
    min_events: int = 3,
    min_density: float = 0.2
) -> List[MismatchCluster]:
    """Detect clusters of mismatches/indels that may indicate missed junctions.

    mapPacBio doesn't soft-clip, so failed junction crossings manifest as
    regions of high mismatch/indel density where the aligner is forcing
    alignment through what should be an intron.

    Args:
        read: pysam AlignedSegment
        window_size: Size of sliding window for density calculation
        min_events: Minimum events (mismatches + indels) in window
        min_density: Minimum events per bp to flag as cluster

    Returns:
        List of MismatchCluster objects
    """
    if not read.cigartuples:
        return []

    strand = '-' if read.is_reverse else '+'
    chrom = read.reference_name

    # Build list of mismatch/indel positions
    events = []  # List of (ref_pos, query_pos, event_type, length)

    ref_pos = read.reference_start
    query_pos = 0

    # Get MD tag for mismatch positions if available
    try:
        md_tag = read.get_tag('MD')
        has_md = True
    except KeyError:
        has_md = False

    for op, length in read.cigartuples:
        if op == 0:  # M (match/mismatch)
            # Without MD tag, we can't distinguish match from mismatch
            # We'll rely on insertions/deletions for now
            ref_pos += length
            query_pos += length
        elif op == 1:  # I (insertion)
            events.append((ref_pos, query_pos, 'I', length))
            query_pos += length
        elif op == 2:  # D (deletion)
            events.append((ref_pos, query_pos, 'D', length))
            ref_pos += length
        elif op == 3:  # N (skip/intron) - ignore, this is a proper junction
            ref_pos += length
        elif op == 4:  # S (soft-clip)
            query_pos += length
        elif op == 7:  # = (sequence match)
            ref_pos += length
            query_pos += length
        elif op == 8:  # X (sequence mismatch)
            events.append((ref_pos, query_pos, 'X', length))
            ref_pos += length
            query_pos += length

    if not events:
        return []

    # Sliding window to find clusters
    clusters = []
    i = 0

    while i < len(events):
        # Find all events within window_size of this event
        window_start = events[i][0]
        window_events = []

        j = i
        while j < len(events) and events[j][0] - window_start <= window_size:
            window_events.append(events[j])
            j += 1

        # Calculate density
        if len(window_events) >= min_events:
            window_end = window_events[-1][0] + window_events[-1][3]
            window_span = max(1, window_end - window_start)
            density = len(window_events) / window_span

            if density >= min_density:
                n_mismatches = sum(1 for e in window_events if e[2] == 'X')
                n_insertions = sum(1 for e in window_events if e[2] == 'I')
                n_deletions = sum(1 for e in window_events if e[2] == 'D')

                cluster = MismatchCluster(
                    chrom=chrom,
                    start=window_start,
                    end=window_end,
                    strand=strand,
                    n_mismatches=n_mismatches,
                    n_insertions=n_insertions,
                    n_deletions=n_deletions,
                    total_events=len(window_events),
                    density=density,
                    read_id=read.query_name,
                    position_in_read=window_events[0][1]
                )
                clusters.append(cluster)

                # Skip past this cluster
                i = j
                continue

        i += 1

    return clusters


def detect_large_deletions_as_potential_junctions(
    read: pysam.AlignedSegment,
    min_deletion_size: int = 30,
    splice_index: Optional[SpliceSiteIndex] = None,
    max_distance_from_ss: int = 10
) -> List[Dict]:
    """Detect large deletions that may represent missed splice junctions.

    mapPacBio optimizes alignments by shifting indels around. A "missed"
    splice junction often manifests as:
    - A large deletion (size similar to intron)
    - Potentially with flanking mismatches where junction refinement failed

    This is different from minimap2's soft-clipping approach.

    Args:
        read: pysam AlignedSegment
        min_deletion_size: Minimum deletion size to consider (bp)
        splice_index: Optional splice site index for validation
        max_distance_from_ss: Max distance from known SS to validate

    Returns:
        List of dicts with potential missed junction info
    """
    if not read.cigartuples:
        return []

    strand = '-' if read.is_reverse else '+'
    chrom = read.reference_name

    potential_junctions = []
    ref_pos = read.reference_start

    for i, (op, length) in enumerate(read.cigartuples):
        if op == 2 and length >= min_deletion_size:  # D = deletion
            # This large deletion could be a missed intron
            deletion_start = ref_pos
            deletion_end = ref_pos + length

            result = {
                'chrom': chrom,
                'strand': strand,
                'deletion_start': deletion_start,
                'deletion_end': deletion_end,
                'deletion_size': length,
                'read_id': read.query_name,
                'cigar_index': i,
                'near_annotated_ss': False,
                'matching_intron': None
            }

            # Check if this deletion aligns with known splice sites
            if splice_index:
                # Look for 5'SS near deletion start and 3'SS near deletion end
                nearby_5ss = splice_index.get_nearby_5ss(
                    chrom, deletion_start, max_distance_from_ss
                )
                nearby_3ss = splice_index.get_nearby_3ss(
                    chrom, deletion_end - 2, max_distance_from_ss  # AG is 2bp before end
                )

                for ss5 in nearby_5ss:
                    for ss3 in nearby_3ss:
                        if ss5.intron_id == ss3.intron_id and ss5.strand == strand:
                            result['near_annotated_ss'] = True
                            result['matching_intron'] = ss5.intron_id
                            break

            potential_junctions.append(result)

        # Update reference position
        if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
            ref_pos += length

    return potential_junctions


def detect_junction_boundary_artifacts(
    read: pysam.AlignedSegment,
    splice_index: SpliceSiteIndex,
    window_around_junction: int = 10
) -> List[Dict]:
    """Detect alignment artifacts at junction boundaries.

    mapPacBio shifts indels and mismatches to optimize alignment scores.
    This can create characteristic patterns around missed junctions:
    - Mismatches clustered just before/after where a junction should be
    - Small indels that "absorb" the junction boundary

    Args:
        read: pysam AlignedSegment
        splice_index: Index of known splice sites
        window_around_junction: Window size to check around known SS

    Returns:
        List of potential refinement opportunities
    """
    if not read.cigartuples:
        return []

    strand = '-' if read.is_reverse else '+'
    chrom = read.reference_name

    # Get all known splice site positions in this read's region
    read_start = read.reference_start
    read_end = read.reference_end

    # Collect 3'SS and 5'SS positions that fall within the read
    relevant_sites = []

    for pos, site in splice_index.three_ss.get(chrom, {}).items():
        if read_start <= pos <= read_end and site.strand == strand:
            relevant_sites.append(('3ss', pos, site))

    for pos, site in splice_index.five_ss.get(chrom, {}).items():
        if read_start <= pos <= read_end and site.strand == strand:
            relevant_sites.append(('5ss', pos, site))

    if not relevant_sites:
        return []

    # Build a map of ref_pos -> CIGAR operations for quick lookup
    # and track where indels/mismatches occur
    event_positions = []  # (ref_pos, event_type, length)

    ref_pos = read.reference_start
    for op, length in read.cigartuples:
        if op == 1:  # Insertion
            event_positions.append((ref_pos, 'I', length))
        elif op == 2:  # Deletion
            event_positions.append((ref_pos, 'D', length))
        elif op == 8:  # Mismatch (X)
            event_positions.append((ref_pos, 'X', length))

        if op in (0, 2, 3, 7, 8):
            ref_pos += length

    # Check each splice site for nearby artifacts
    artifacts = []

    for site_type, pos, site in relevant_sites:
        # Count events in window around this splice site
        nearby_events = [
            e for e in event_positions
            if abs(e[0] - pos) <= window_around_junction
        ]

        if nearby_events:
            total_event_size = sum(e[2] for e in nearby_events)
            n_deletions = sum(1 for e in nearby_events if e[1] == 'D')
            total_deletion_size = sum(e[2] for e in nearby_events if e[1] == 'D')

            # Flag if there's significant artifact near this SS
            if total_event_size >= 3 or total_deletion_size >= 10:
                artifacts.append({
                    'chrom': chrom,
                    'strand': strand,
                    'splice_site_type': site_type,
                    'splice_site_pos': pos,
                    'intron_id': site.intron_id,
                    'n_events_nearby': len(nearby_events),
                    'total_event_size': total_event_size,
                    'n_deletions': n_deletions,
                    'total_deletion_size': total_deletion_size,
                    'read_id': read.query_name,
                    'events': nearby_events
                })

    return artifacts


def analyze_mapPacBio_mismatch_patterns(
    bam_path: str,
    splice_index: SpliceSiteIndex,
    max_distance: int = 20
) -> Dict[str, int]:
    """Analyze mismatch clusters in mapPacBio alignments near splice sites.

    This diagnostic function helps understand how often mapPacBio
    introduces mismatch clusters where splice junctions should be.

    Args:
        bam_path: Path to mapPacBio BAM file
        splice_index: Index of known splice sites
        max_distance: Maximum distance from splice site

    Returns:
        Dict with statistics
    """
    stats = {
        'total_reads': 0,
        'reads_with_mismatch_clusters': 0,
        'clusters_near_splice_sites': 0,
        'clusters_not_near_splice_sites': 0,
        'total_clusters': 0
    }

    bam = pysam.AlignmentFile(bam_path, 'rb')

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        stats['total_reads'] += 1

        clusters = detect_mismatch_clusters(read)

        if clusters:
            stats['reads_with_mismatch_clusters'] += 1
            stats['total_clusters'] += len(clusters)

            for cluster in clusters:
                # Check if cluster is near a splice site
                nearby_5ss = splice_index.get_nearby_5ss(
                    cluster.chrom, cluster.start, max_distance
                )
                nearby_3ss = splice_index.get_nearby_3ss(
                    cluster.chrom, cluster.start, max_distance
                )

                if nearby_5ss or nearby_3ss:
                    stats['clusters_near_splice_sites'] += 1
                else:
                    stats['clusters_not_near_splice_sites'] += 1

    bam.close()

    return stats


def get_soft_clip_info(read: pysam.AlignedSegment) -> Dict:
    """Extract soft-clip information from a read.

    Args:
        read: pysam AlignedSegment

    Returns:
        Dict with 'left' and 'right' soft-clip info
    """
    info = {
        'left': {'length': 0, 'sequence': '', 'quality': None},
        'right': {'length': 0, 'sequence': '', 'quality': None}
    }

    if not read.cigartuples:
        return info

    seq = read.query_sequence
    qual = read.query_qualities

    # Left soft-clip (5' for + strand)
    if read.cigartuples[0][0] == 4:  # S = soft clip
        length = read.cigartuples[0][1]
        info['left'] = {
            'length': length,
            'sequence': seq[:length],
            'quality': qual[:length] if qual else None
        }

    # Right soft-clip (3' for + strand)
    if read.cigartuples[-1][0] == 4:
        length = read.cigartuples[-1][1]
        info['right'] = {
            'length': length,
            'sequence': seq[-length:],
            'quality': qual[-length:] if qual else None
        }

    return info


def simple_align(query: str, target: str, max_mismatches: int = None) -> Tuple[int, int, int]:
    """Simple alignment of query to target allowing mismatches.

    Args:
        query: Query sequence
        target: Target sequence to search in
        max_mismatches: Maximum allowed mismatches (default: 30% of query length)

    Returns:
        Tuple of (best_position, best_score, n_mismatches)
        Position is 0-based index in target where best alignment starts
        Returns (-1, 0, len(query)) if no good alignment found
    """
    if not query or not target:
        return -1, 0, len(query) if query else 0

    if max_mismatches is None:
        max_mismatches = int(len(query) * 0.3)

    best_pos = -1
    best_score = 0
    best_mismatches = len(query)

    # Slide query along target
    for i in range(len(target) - len(query) + 1):
        mismatches = 0
        for j, (q, t) in enumerate(zip(query, target[i:i+len(query)])):
            if q.upper() != t.upper():
                mismatches += 1
                if mismatches > max_mismatches:
                    break

        score = len(query) - mismatches
        if score > best_score:
            best_score = score
            best_pos = i
            best_mismatches = mismatches

    return best_pos, best_score, best_mismatches


def attempt_5prime_refinement(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    splice_index: SpliceSiteIndex,
    min_softclip: int = MIN_SOFTCLIP_LENGTH,
    max_search_distance: int = MAX_SPLICE_SITE_DISTANCE
) -> RefinementResult:
    """Attempt to refine the 5' end of a read by finding upstream exon.

    For reads with soft-clipped 5' ends, this function:
    1. Checks if there's a known 3'SS (AG) near the aligned 5' end
    2. Searches upstream for where the soft-clipped sequence might align
    3. If found, proposes a refined alignment with splice junction

    Args:
        read: pysam AlignedSegment
        genome: Dict mapping chrom -> sequence
        splice_index: Index of known splice sites
        min_softclip: Minimum soft-clip length to attempt refinement
        max_search_distance: Maximum distance to search for splice site

    Returns:
        RefinementResult with success status and refined CIGAR if successful
    """
    result = RefinementResult(
        read_id=read.query_name,
        original_cigar=read.cigarstring
    )

    strand = '-' if read.is_reverse else '+'
    chrom = read.reference_name

    # Get soft-clip info
    clip_info = get_soft_clip_info(read)

    # For + strand, 5' soft-clip is on the left
    # For - strand, 5' soft-clip is on the right
    if strand == '+':
        clip = clip_info['left']
        aligned_5prime = read.reference_start  # Position where alignment begins
    else:
        clip = clip_info['right']
        aligned_5prime = read.reference_end - 1  # Rightmost aligned position

    result.soft_clip_length = clip['length']
    result.soft_clip_sequence = clip['sequence']

    # Check if soft-clip is long enough
    if clip['length'] < min_softclip:
        result.reason = f"Soft-clip too short ({clip['length']} < {min_softclip})"
        return result

    # Look for nearby 3'SS (the acceptor site at the 5' end of an exon)
    # For + strand: 3'SS should be just before the aligned start
    # For - strand: 3'SS should be just after the aligned end
    nearby_3ss = splice_index.get_nearby_3ss(chrom, aligned_5prime, max_search_distance)

    if not nearby_3ss:
        result.reason = f"No known 3'SS within {max_search_distance}bp"
        return result

    # Try each nearby 3'SS
    if chrom not in genome:
        result.reason = f"Chromosome {chrom} not in genome"
        return result

    genome_seq = genome[chrom]

    for ss in nearby_3ss:
        # Only consider sites on the same strand
        if ss.strand != strand:
            continue

        # The soft-clipped sequence should align to the upstream exon
        # (before the intron start = 5'SS)
        # We need to find the corresponding 5'SS for this 3'SS
        intron_id = ss.intron_id

        # Search for matching 5'SS (same intron)
        matching_5ss = None
        for pos, site in splice_index.five_ss.get(chrom, {}).items():
            if site.intron_id == intron_id:
                matching_5ss = site
                break

        if not matching_5ss:
            continue

        # Calculate the upstream exon region to search
        if strand == '+':
            # Upstream exon is to the left of the 5'SS
            search_start = max(0, matching_5ss.position - 500)
            search_end = matching_5ss.position
            search_region = genome_seq[search_start:search_end]

            # Try to align soft-clipped sequence to end of upstream exon
            # The soft-clip should align to the 3' end of the upstream exon
            query_rc = clip['sequence']  # For + strand, use as-is
        else:
            # For - strand, upstream exon is to the right of 5'SS
            search_start = matching_5ss.position + 2  # After GT
            search_end = min(len(genome_seq), search_start + 500)
            search_region = genome_seq[search_start:search_end]

            # Need reverse complement of soft-clip
            query_rc = reverse_complement(clip['sequence'])

        # Attempt alignment
        pos, score, mismatches = simple_align(query_rc, search_region)

        if pos < 0:
            continue

        # Calculate alignment score as fraction of perfect match
        score_fraction = score / len(clip['sequence'])

        if score_fraction >= MIN_ALIGNMENT_SCORE_FRACTION:
            result.success = True
            result.matched_splice_site = ss
            result.alignment_score = score_fraction

            if strand == '+':
                # New 5' end is in the upstream exon
                new_start = search_start + pos
                intron_start = search_start + pos + len(clip['sequence'])
                intron_end = aligned_5prime
            else:
                new_start = search_start + pos
                intron_start = aligned_5prime + 1
                intron_end = new_start

            result.junction_added = (intron_start, intron_end)
            result.refined_start = new_start
            result.reason = f"Found alignment to upstream exon (score={score_fraction:.2f})"

            # Build refined CIGAR
            # Original CIGAR with soft-clip replaced by match + intron + match
            # This is simplified - full implementation would rebuild entire CIGAR
            intron_length = abs(intron_end - intron_start)
            result.refined_cigar = f"{len(clip['sequence'])}M{intron_length}N{read.cigarstring.replace(f'{clip['length']}S', '')}"

            return result

    result.reason = "No good alignment found to any upstream exon"
    return result


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def refine_terminal_exons(
    bam_path: str,
    genome: Dict[str, str],
    splice_index: SpliceSiteIndex,
    output_path: Optional[str] = None,
    min_softclip: int = MIN_SOFTCLIP_LENGTH
) -> Dict[str, RefinementResult]:
    """Refine terminal exons for all reads in a BAM file.

    Args:
        bam_path: Path to input BAM file
        genome: Dict mapping chrom -> sequence
        splice_index: Index of known splice sites
        output_path: Optional path for refined BAM output
        min_softclip: Minimum soft-clip length to attempt refinement

    Returns:
        Dict mapping read_id to RefinementResult
    """
    results = {}

    bam = pysam.AlignmentFile(bam_path, 'rb')

    n_total = 0
    n_with_softclip = 0
    n_refined = 0

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        n_total += 1

        # Check if read has significant soft-clip at 5' end
        clip_info = get_soft_clip_info(read)
        strand = '-' if read.is_reverse else '+'

        five_prime_clip = clip_info['left'] if strand == '+' else clip_info['right']

        if five_prime_clip['length'] >= min_softclip:
            n_with_softclip += 1

            result = attempt_5prime_refinement(
                read, genome, splice_index, min_softclip
            )
            results[read.query_name] = result

            if result.success:
                n_refined += 1

    bam.close()

    logger.info(f"Processed {n_total} reads")
    logger.info(f"  With soft-clip >= {min_softclip}bp: {n_with_softclip}")
    logger.info(f"  Successfully refined: {n_refined}")

    return results


def build_comprehensive_splice_index(
    gff_path: Optional[str] = None,
    bam_paths: Optional[List[str]] = None,
    min_reads_per_junction: int = 2
) -> SpliceSiteIndex:
    """Build a comprehensive splice site index from annotation AND alignments.

    This combines:
    1. Known junctions from GFF annotation (high confidence)
    2. Novel junctions detected by aligners (may be alternative/cryptic)

    The combination ensures we can refine terminal exons using both
    well-characterized splice sites and newly discovered ones.

    Args:
        gff_path: Path to GFF annotation file (optional)
        bam_paths: List of BAM file paths to extract novel junctions (optional)
        min_reads_per_junction: Minimum reads for novel junction inclusion

    Returns:
        Merged SpliceSiteIndex with all splice sites
    """
    indices = []

    # Load annotated junctions first (takes precedence)
    if gff_path:
        logger.info(f"Loading annotated junctions from {gff_path}")
        annotated_index = load_splice_sites_from_gff(gff_path)
        indices.append(annotated_index)

    # Extract novel junctions from BAM files
    if bam_paths:
        for bam_path in bam_paths:
            logger.info(f"Extracting junctions from {bam_path}")
            bam_index = extract_junctions_from_bam(bam_path, min_reads_per_junction)
            indices.append(bam_index)

    if not indices:
        logger.warning("No splice sites loaded - provide GFF and/or BAM files")
        return SpliceSiteIndex()

    # Merge all indices
    return merge_splice_indices(*indices)


def analyze_soft_clips_near_junctions(
    bam_path: str,
    splice_index: SpliceSiteIndex,
    max_distance: int = 20
) -> Dict[str, int]:
    """Analyze how many soft-clips occur near known splice sites.

    This diagnostic function helps understand the potential for refinement.

    Args:
        bam_path: Path to BAM file
        splice_index: Index of known splice sites
        max_distance: Maximum distance from splice site

    Returns:
        Dict with statistics
    """
    stats = {
        'total_reads': 0,
        'reads_with_5prime_softclip': 0,
        'softclips_near_3ss': 0,
        'softclips_not_near_3ss': 0,
        'softclip_lengths': []
    }

    bam = pysam.AlignmentFile(bam_path, 'rb')

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        stats['total_reads'] += 1

        clip_info = get_soft_clip_info(read)
        strand = '-' if read.is_reverse else '+'
        chrom = read.reference_name

        five_prime_clip = clip_info['left'] if strand == '+' else clip_info['right']

        if five_prime_clip['length'] >= MIN_SOFTCLIP_LENGTH:
            stats['reads_with_5prime_softclip'] += 1
            stats['softclip_lengths'].append(five_prime_clip['length'])

            # Check if near a 3'SS
            aligned_5prime = read.reference_start if strand == '+' else read.reference_end - 1
            nearby_3ss = splice_index.get_nearby_3ss(chrom, aligned_5prime, max_distance)

            if nearby_3ss:
                stats['softclips_near_3ss'] += 1
            else:
                stats['softclips_not_near_3ss'] += 1

    bam.close()

    return stats
