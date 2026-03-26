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
from typing import Dict, List, Optional, Set, Tuple
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


def detect_junction_truncated_reads(
    bam_path: str,
    splice_index: SpliceSiteIndex,
    tolerance: int = 10
) -> Dict[str, List[Dict]]:
    """Detect reads that are truncated at junction boundaries.

    These are reads that end precisely at a splice site without crossing
    the junction. Common in mapPacBio alignments where soft-clipping is rare.

    For minus strand genes:
    - Reads ending at 3'SS (acceptor) without upstream exon coverage
    - These represent 5' truncated transcripts

    Args:
        bam_path: Path to BAM file
        splice_index: Index of known splice sites
        tolerance: Position tolerance for boundary detection (bp)

    Returns:
        Dict with 'truncated_at_3ss' and 'truncated_at_5ss' lists
    """
    results = {
        'truncated_at_3ss': [],
        'truncated_at_5ss': [],
        'stats': {
            'total_reads': 0,
            'truncated_at_3ss': 0,
            'truncated_at_5ss': 0,
            'properly_spliced': 0
        }
    }

    bam = pysam.AlignmentFile(bam_path, 'rb')

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        results['stats']['total_reads'] += 1

        chrom = read.reference_name
        strand = '-' if read.is_reverse else '+'
        read_start = read.reference_start
        read_end = read.reference_end

        # Check if read has any splice junctions
        has_junction = any(op == 3 for op, _ in (read.cigartuples or []))

        # Check if read ends near a known 3'SS
        for pos, site in splice_index.three_ss.get(chrom, {}).items():
            if site.strand != strand:
                continue

            # For + strand: 3'SS is at higher position, read ends there
            # For - strand: 3'SS is at lower position, read starts there
            if strand == '+':
                if abs(read_end - pos - 2) <= tolerance and not has_junction:
                    results['truncated_at_3ss'].append({
                        'read_id': read.query_name,
                        'chrom': chrom,
                        'strand': strand,
                        'read_start': read_start,
                        'read_end': read_end,
                        'splice_site_pos': pos,
                        'intron_id': site.intron_id,
                        'distance_from_ss': read_end - pos - 2
                    })
                    results['stats']['truncated_at_3ss'] += 1
                    break
            else:  # minus strand
                if abs(read_end - pos - 2) <= tolerance and not has_junction:
                    results['truncated_at_3ss'].append({
                        'read_id': read.query_name,
                        'chrom': chrom,
                        'strand': strand,
                        'read_start': read_start,
                        'read_end': read_end,
                        'splice_site_pos': pos,
                        'intron_id': site.intron_id,
                        'distance_from_ss': read_end - pos - 2
                    })
                    results['stats']['truncated_at_3ss'] += 1
                    break

        # Check if read starts near a known 5'SS (truncated at 5' end)
        for pos, site in splice_index.five_ss.get(chrom, {}).items():
            if site.strand != strand:
                continue

            if strand == '+':
                if abs(read_start - pos) <= tolerance and not has_junction:
                    results['truncated_at_5ss'].append({
                        'read_id': read.query_name,
                        'chrom': chrom,
                        'strand': strand,
                        'read_start': read_start,
                        'read_end': read_end,
                        'splice_site_pos': pos,
                        'intron_id': site.intron_id,
                        'distance_from_ss': read_start - pos
                    })
                    results['stats']['truncated_at_5ss'] += 1
                    break
            else:
                if abs(read_start - pos) <= tolerance and not has_junction:
                    results['truncated_at_5ss'].append({
                        'read_id': read.query_name,
                        'chrom': chrom,
                        'strand': strand,
                        'read_start': read_start,
                        'read_end': read_end,
                        'splice_site_pos': pos,
                        'intron_id': site.intron_id,
                        'distance_from_ss': read_start - pos
                    })
                    results['stats']['truncated_at_5ss'] += 1
                    break

        if has_junction:
            results['stats']['properly_spliced'] += 1

    bam.close()

    logger.info(f"Total reads: {results['stats']['total_reads']:,}")
    logger.info(f"Properly spliced: {results['stats']['properly_spliced']:,}")
    logger.info(f"Truncated at 3'SS: {results['stats']['truncated_at_3ss']:,}")
    logger.info(f"Truncated at 5'SS: {results['stats']['truncated_at_5ss']:,}")

    return results


@dataclass
class PartialJunctionEvidence:
    """Evidence of a partial junction crossing from soft-clipped sequence.

    Even a few nucleotides of soft-clipped sequence at a splice site boundary
    can provide evidence of a splicing event if they match the expected
    upstream/downstream exon sequence.

    Ambiguity handling:
    - ambiguous=True: Multiple 5'SS have identical upstream sequences
    - alternative_introns: List of other introns that match equally well
    - For ambiguous reads, use proportional attribution or hold out
    """
    read_id: str
    chrom: str
    strand: str
    read_start: int
    read_end: int
    splice_site_pos: int
    splice_site_type: str  # '3ss' or '5ss'
    intron_id: str
    soft_clip_length: int
    soft_clip_sequence: str
    expected_exon_sequence: str
    matches: int
    mismatches: int
    match_fraction: float
    inferred_junction: Tuple[int, int]  # (intron_start, intron_end)
    confidence: str  # 'high', 'medium', 'low'
    ambiguous: bool = False  # True if multiple equally-good matches exist
    alternative_introns: List[str] = field(default_factory=list)  # Other introns with same upstream seq


def detect_partial_junction_crossings(
    bam_path: str,
    genome: Dict[str, str],
    splice_index: SpliceSiteIndex,
    boundary_tolerance: int = 5,
    min_clip_length: int = 1,
    min_match_fraction: float = 0.6,
    ambiguous_mode: str = 'proportional'
) -> Dict[str, any]:
    """Detect reads with soft-clips at splice sites that provide junction evidence.

    This function identifies reads that:
    1. End precisely at a splice site boundary (within tolerance)
    2. Have soft-clipped sequence at the 5' end
    3. The soft-clipped sequence matches the expected upstream exon

    Even a single nucleotide match provides some evidence of splicing.
    This rescues the ~71% of reads with 11-20bp soft-clips that would
    otherwise be counted as unspliced.

    Ambiguity handling:
    When multiple 5'SS have identical upstream sequences that match the soft-clip,
    the read is marked as ambiguous. Options:
    - 'proportional': Attribute reads in same proportion as confident reads
    - 'holdout': Exclude ambiguous reads from counts

    Args:
        bam_path: Path to BAM file
        genome: Dict mapping chrom -> sequence
        splice_index: Index of known splice sites
        boundary_tolerance: Max distance from splice site boundary (bp)
        min_clip_length: Minimum soft-clip length to consider (default: 1)
        min_match_fraction: Minimum fraction of matches required (0.0-1.0)
        ambiguous_mode: How to handle ambiguous reads ('proportional' or 'holdout')

    Returns:
        Dict with:
            - 'evidence': List of PartialJunctionEvidence objects
            - 'stats': Statistics dict
            - 'by_intron': Dict mapping intron_id to list of evidence
    """
    results = {
        'evidence': [],
        'stats': {
            'total_reads': 0,
            'reads_at_splice_boundary': 0,
            'reads_with_softclip_at_boundary': 0,
            'reads_with_matching_sequence': 0,
            'reads_rescued_as_spliced': 0,
            'reads_ambiguous': 0,
            'by_clip_length': {},
            'by_confidence': {'high': 0, 'medium': 0, 'low': 0}
        },
        'by_intron': {},
        'ambiguous_mode': ambiguous_mode
    }

    # Pre-compute upstream exon sequences for all 5'SS to detect ambiguity
    # Key: (chrom, strand, sequence) -> list of intron_ids with this upstream seq
    upstream_seq_to_introns: Dict[Tuple[str, str, str], List[str]] = {}

    for chrom, sites in splice_index.five_ss.items():
        if chrom not in genome:
            continue
        genome_seq = genome[chrom]

        for pos, site in sites.items():
            # Get upstream exon sequence (first 30bp for comparison)
            if site.strand == '-':
                # Upstream exon is at higher coords (past 5'SS)
                exon_start = pos + 2
                exon_end = min(exon_start + 30, len(genome_seq))
                upstream_seq = genome_seq[exon_start:exon_end]
            else:
                # Upstream exon is at lower coords (before 5'SS)
                exon_end = pos
                exon_start = max(0, exon_end - 30)
                upstream_seq = genome_seq[exon_start:exon_end]

            key = (chrom, site.strand, upstream_seq[:20])  # Use first 20bp for comparison
            if key not in upstream_seq_to_introns:
                upstream_seq_to_introns[key] = []
            upstream_seq_to_introns[key].append(site.intron_id)

    # Log potential ambiguity
    n_ambiguous_seqs = sum(1 for introns in upstream_seq_to_introns.values() if len(introns) > 1)
    if n_ambiguous_seqs > 0:
        logger.info(f"Found {n_ambiguous_seqs} upstream exon sequences shared by multiple introns")

    bam = pysam.AlignmentFile(bam_path, 'rb')

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        results['stats']['total_reads'] += 1

        chrom = read.reference_name
        strand = '-' if read.is_reverse else '+'
        read_start = read.reference_start
        read_end = read.reference_end

        # Skip if already has a junction (properly spliced)
        has_junction = any(op == 3 for op, _ in (read.cigartuples or []))
        if has_junction:
            continue

        # Get soft-clip info
        clip_info = get_soft_clip_info(read)

        # For - strand (direct RNA): right clip is at 5' end of RNA (high genomic coord)
        # The read ends at ref_end, and we're looking for 3'SS near ref_end
        if strand == '-':
            five_prime_clip = clip_info['right']
            boundary_pos = read_end  # Position where alignment ends
        else:
            five_prime_clip = clip_info['left']
            boundary_pos = read_start  # Position where alignment starts

        # Check if read ends near a 3'SS (acceptor site at intron end)
        # For - strand: 3'SS is at the lower coordinate of the intron
        for pos, site in splice_index.three_ss.get(chrom, {}).items():
            if site.strand != strand:
                continue

            # Check if boundary is near this 3'SS
            # For - strand, the 3'SS position marks where the intron ends
            # Reads ending at the 3'SS have ref_end near pos
            distance = abs(boundary_pos - pos)

            if distance <= boundary_tolerance:
                results['stats']['reads_at_splice_boundary'] += 1

                # Check if there's a soft-clip
                if five_prime_clip['length'] < min_clip_length:
                    continue

                results['stats']['reads_with_softclip_at_boundary'] += 1

                # Track clip length distribution
                clip_len_bin = f"{(five_prime_clip['length'] // 5) * 5}-{(five_prime_clip['length'] // 5) * 5 + 4}"
                results['stats']['by_clip_length'][clip_len_bin] = \
                    results['stats']['by_clip_length'].get(clip_len_bin, 0) + 1

                # Find the corresponding 5'SS for this intron
                intron_id = site.intron_id
                matching_5ss = None
                for fpos, fsite in splice_index.five_ss.get(chrom, {}).items():
                    if fsite.intron_id == intron_id and fsite.strand == strand:
                        matching_5ss = fsite
                        break

                if not matching_5ss:
                    continue

                # Get the expected upstream exon sequence
                # For - strand: upstream exon is at HIGHER coordinates (past the 5'SS)
                # The soft-clipped sequence should match the end of the upstream exon
                if chrom not in genome:
                    continue

                genome_seq = genome[chrom]
                clip_seq = five_prime_clip['sequence']

                if strand == '-':
                    # For - strand, upstream exon is at higher coords (past 5'SS)
                    # The 5'SS position is the rightmost position of the intron
                    # Upstream exon starts just after 5'SS
                    exon_start = matching_5ss.position + 2  # Just past the 5'SS
                    exon_end = min(exon_start + len(clip_seq) + 20, len(genome_seq))
                    expected_seq = genome_seq[exon_start:exon_end]

                    # IMPORTANT: pysam stores soft-clipped sequence in GENOMIC orientation
                    # (already reverse-complemented for - strand alignments)
                    # So we compare the clip directly to the genomic exon sequence
                    query_seq = clip_seq
                else:
                    # For + strand, upstream exon is at lower coords (before 5'SS)
                    exon_end = matching_5ss.position
                    exon_start = max(0, exon_end - len(clip_seq) - 20)
                    expected_seq = genome_seq[exon_start:exon_end]
                    query_seq = clip_seq

                # Align the soft-clipped sequence to the expected exon
                best_pos, best_score, n_mismatches = simple_align(
                    query_seq, expected_seq, max_mismatches=int(len(query_seq) * (1 - min_match_fraction))
                )

                if best_pos >= 0:
                    n_matches = len(query_seq) - n_mismatches
                    match_fraction = n_matches / len(query_seq) if len(query_seq) > 0 else 0

                    if match_fraction >= min_match_fraction:
                        results['stats']['reads_with_matching_sequence'] += 1

                        # Check for ambiguity: are there other introns with same upstream sequence?
                        # Use the first 20bp of the matching region for comparison
                        compare_len = min(20, len(query_seq))
                        matched_portion = expected_seq[best_pos:best_pos + compare_len]
                        ambiguity_key = (chrom, strand, matched_portion)

                        is_ambiguous = False
                        alternative_introns = []

                        if ambiguity_key in upstream_seq_to_introns:
                            matching_introns = upstream_seq_to_introns[ambiguity_key]
                            if len(matching_introns) > 1:
                                # Multiple introns have identical upstream sequences
                                is_ambiguous = True
                                alternative_introns = [i for i in matching_introns if i != intron_id]
                                results['stats']['reads_ambiguous'] += 1

                        # Determine confidence based on match quality and length
                        if match_fraction >= 0.9 and len(query_seq) >= 10:
                            confidence = 'high'
                        elif match_fraction >= 0.75 and len(query_seq) >= 5:
                            confidence = 'medium'
                        else:
                            confidence = 'low'

                        results['stats']['by_confidence'][confidence] += 1
                        results['stats']['reads_rescued_as_spliced'] += 1

                        # Determine inferred junction coordinates
                        # The 3'SS position marks where the AG dinucleotide starts
                        # The 5'SS position marks where the GT dinucleotide starts
                        # CIGAR intron coords use half-open: [intron_start, intron_end)
                        # where intron_end is 2bp past the AG (exclusive end)
                        if strand == '-':
                            # Minus strand: 3'SS at lower coord, 5'SS at higher coord
                            # intron_start = 3'SS position
                            # intron_end = 5'SS position + 2 (past the GT)
                            inferred_junction = (pos, matching_5ss.position + 2)
                        else:
                            # Plus strand: 5'SS at lower coord, 3'SS at higher coord
                            # intron_start = 5'SS position
                            # intron_end = 3'SS position + 2 (past the AG)
                            inferred_junction = (matching_5ss.position, pos + 2)

                        evidence = PartialJunctionEvidence(
                            read_id=read.query_name,
                            chrom=chrom,
                            strand=strand,
                            read_start=read_start,
                            read_end=read_end,
                            splice_site_pos=pos,
                            splice_site_type='3ss',
                            intron_id=intron_id,
                            soft_clip_length=five_prime_clip['length'],
                            soft_clip_sequence=clip_seq,
                            expected_exon_sequence=expected_seq[:len(query_seq)+5],
                            matches=n_matches,
                            mismatches=n_mismatches,
                            match_fraction=match_fraction,
                            inferred_junction=inferred_junction,
                            confidence=confidence,
                            ambiguous=is_ambiguous,
                            alternative_introns=alternative_introns
                        )

                        results['evidence'].append(evidence)

                        if intron_id not in results['by_intron']:
                            results['by_intron'][intron_id] = []
                        results['by_intron'][intron_id].append(evidence)

                break  # Only process one splice site per read

    bam.close()

    logger.info(f"Partial junction detection complete:")
    logger.info(f"  Total reads: {results['stats']['total_reads']:,}")
    logger.info(f"  At splice boundary: {results['stats']['reads_at_splice_boundary']:,}")
    logger.info(f"  With soft-clip at boundary: {results['stats']['reads_with_softclip_at_boundary']:,}")
    logger.info(f"  With matching sequence: {results['stats']['reads_with_matching_sequence']:,}")
    logger.info(f"  Rescued as spliced: {results['stats']['reads_rescued_as_spliced']:,}")
    logger.info(f"  Ambiguous (multiple matching introns): {results['stats']['reads_ambiguous']:,}")
    logger.info(f"  Ambiguous mode: {ambiguous_mode}")

    return results


def compute_junction_counts_with_ambiguity(
    detection_results: Dict,
    mode: str = 'proportional'
) -> Dict[str, Dict]:
    """Compute junction counts with ambiguity handling.

    For unambiguous reads: count directly toward their intron
    For ambiguous reads:
    - 'proportional': Attribute in same proportion as confident reads
    - 'holdout': Exclude from counts

    Args:
        detection_results: Output from detect_partial_junction_crossings()
        mode: 'proportional' or 'holdout'

    Returns:
        Dict mapping intron_id to {
            'confident_count': int,  # Unambiguous reads
            'ambiguous_count': int,  # Ambiguous reads (raw)
            'attributed_count': float,  # After proportional attribution
            'total_count': float,  # Final count (confident + attributed)
            'alternative_introns': set,  # Introns that share upstream seq
        }
    """
    counts = {}

    # First pass: count confident (unambiguous) reads per intron
    for evidence in detection_results['evidence']:
        intron_id = evidence.intron_id

        if intron_id not in counts:
            counts[intron_id] = {
                'confident_count': 0,
                'ambiguous_count': 0,
                'ambiguous_reads': [],  # Store for proportional attribution
                'attributed_count': 0.0,
                'total_count': 0.0,
                'alternative_introns': set()
            }

        if evidence.ambiguous:
            counts[intron_id]['ambiguous_count'] += 1
            counts[intron_id]['ambiguous_reads'].append(evidence)
            counts[intron_id]['alternative_introns'].update(evidence.alternative_introns)
        else:
            counts[intron_id]['confident_count'] += 1

    if mode == 'holdout':
        # Simply use confident counts only
        for intron_id in counts:
            counts[intron_id]['attributed_count'] = 0.0
            counts[intron_id]['total_count'] = float(counts[intron_id]['confident_count'])

    elif mode == 'proportional':
        # For each group of introns sharing upstream sequences,
        # distribute ambiguous reads proportionally based on confident counts

        # Build groups of introns that share upstream sequences
        processed_groups: Set[frozenset] = set()

        for intron_id, data in counts.items():
            if data['ambiguous_count'] == 0:
                # No ambiguous reads, just use confident count
                counts[intron_id]['total_count'] = float(data['confident_count'])
                continue

            # Build the group of introns that share this sequence
            group = {intron_id} | data['alternative_introns']
            group_key = frozenset(group)

            if group_key in processed_groups:
                continue
            processed_groups.add(group_key)

            # Get confident counts for all introns in the group
            group_confident_counts = {}
            for g_intron in group:
                if g_intron in counts:
                    group_confident_counts[g_intron] = counts[g_intron]['confident_count']
                else:
                    group_confident_counts[g_intron] = 0

            total_confident = sum(group_confident_counts.values())

            if total_confident == 0:
                # No confident reads in any intron - distribute equally
                for g_intron in group:
                    if g_intron in counts:
                        n_ambiguous = counts[g_intron]['ambiguous_count']
                        fraction = 1.0 / len(group)
                        counts[g_intron]['attributed_count'] = n_ambiguous * fraction
            else:
                # Distribute proportionally based on confident counts
                for g_intron in group:
                    if g_intron in counts:
                        n_ambiguous = counts[g_intron]['ambiguous_count']
                        fraction = group_confident_counts[g_intron] / total_confident
                        counts[g_intron]['attributed_count'] = n_ambiguous * fraction

        # Compute total counts
        for intron_id in counts:
            counts[intron_id]['total_count'] = (
                counts[intron_id]['confident_count'] +
                counts[intron_id]['attributed_count']
            )

    # Remove internal tracking data before returning
    for intron_id in counts:
        del counts[intron_id]['ambiguous_reads']
        counts[intron_id]['alternative_introns'] = list(counts[intron_id]['alternative_introns'])

    return counts


def summarize_ambiguity_report(
    detection_results: Dict,
    counts: Dict[str, Dict]
) -> str:
    """Generate a summary report of ambiguity handling.

    Args:
        detection_results: Output from detect_partial_junction_crossings()
        counts: Output from compute_junction_counts_with_ambiguity()

    Returns:
        Formatted string report
    """
    lines = []
    lines.append("=" * 60)
    lines.append("AMBIGUITY HANDLING REPORT")
    lines.append("=" * 60)

    mode = detection_results.get('ambiguous_mode', 'unknown')
    lines.append(f"\nMode: {mode}")

    stats = detection_results['stats']
    total_rescued = stats['reads_rescued_as_spliced']
    total_ambiguous = stats['reads_ambiguous']

    lines.append(f"\nTotal rescued reads: {total_rescued:,}")
    lines.append(f"Ambiguous reads: {total_ambiguous:,} ({100*total_ambiguous/max(1,total_rescued):.1f}%)")
    lines.append(f"Unambiguous reads: {total_rescued - total_ambiguous:,}")

    # Find introns with ambiguity
    introns_with_ambiguity = [
        (intron_id, data) for intron_id, data in counts.items()
        if data['ambiguous_count'] > 0
    ]

    if introns_with_ambiguity:
        lines.append(f"\nIntrons with ambiguous reads: {len(introns_with_ambiguity)}")
        lines.append("-" * 60)

        # Sort by ambiguous count
        introns_with_ambiguity.sort(key=lambda x: x[1]['ambiguous_count'], reverse=True)

        for intron_id, data in introns_with_ambiguity[:20]:  # Top 20
            lines.append(f"\n{intron_id}:")
            lines.append(f"  Confident: {data['confident_count']}")
            lines.append(f"  Ambiguous: {data['ambiguous_count']}")
            lines.append(f"  Attributed: {data['attributed_count']:.1f}")
            lines.append(f"  Total: {data['total_count']:.1f}")
            if data['alternative_introns']:
                lines.append(f"  Shares sequence with: {', '.join(data['alternative_introns'][:3])}")

    lines.append("\n" + "=" * 60)

    return "\n".join(lines)


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
