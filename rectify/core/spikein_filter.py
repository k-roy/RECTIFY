#!/usr/bin/env python3
"""
Spike-in RNA detection and filtering for RECTIFY.

Provides two main capabilities:
1. INFER spike-in signatures from BAM data (for any gene with spike-in controls)
2. FILTER reads using identified spike-in signatures

Spike-in constructs typically contain:
- Real gene CDS (identical to genomic)
- Synthetic 5' UTR (different from genomic)
- Synthetic 3' UTR (different from genomic)

The inference algorithm:
1. Extract reads from a suspected spike-in region
2. Look at sequences following the CDS end (3' UTR region)
3. Compare to genomic 3' UTR sequence
4. Identify the synthetic 3' UTR as the dominant non-genomic pattern

Author: Kevin R. Roy
Date: 2026-03-13
"""

import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import Counter
import pysam


def edit_distance(s1: str, s2: str) -> int:
    """
    Compute Levenshtein edit distance between two strings.

    Uses dynamic programming with O(min(len(s1), len(s2))) space.
    """
    if len(s1) < len(s2):
        return edit_distance(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]


def find_best_match_position(seq: str, pattern: str, window_size: int = None) -> Tuple[int, int, int]:
    """
    Find the position in seq where pattern has the minimum edit distance.

    Returns:
        Tuple of (best_position, min_edit_distance, subsequence_length)
    """
    if window_size is None:
        window_size = len(pattern)

    best_pos = -1
    best_dist = len(pattern) + 1  # Max possible distance
    best_len = window_size

    # Slide window across sequence
    for i in range(len(seq) - window_size + 1):
        subseq = seq[i:i + window_size]
        dist = edit_distance(subseq, pattern)
        if dist < best_dist:
            best_dist = dist
            best_pos = i
            best_len = window_size

    return best_pos, best_dist, best_len


# Pre-defined spike-in signatures (can be extended)
# Uses multiple 6bp cores requiring 2+ matches (handles nanopore ~5-10% errors)
KNOWN_SPIKEIN_SIGNATURES = {
    'ENO2': {
        'gene': 'ENO2',
        'systematic_name': 'YHR174W',
        'chrom': 'ref|NC_001140|',
        'chrom_alt': 'chrVIII',
        'strand': '+',  # PLUS strand gene
        'start': 451297,
        'end': 452951,
        'cds_end': 452640,
        # Synthetic 3' UTR pattern (found in spike-in construct)
        # This appears in the soft-clipped region before poly-A
        'synthetic_3utr_full': 'CGTCGTGAGTAGTGAACCGTAAGC',
        # Multiple 6bp cores - require 2+ matches for classification
        'synthetic_3utr_cores': [
            'CGTCGT', 'TCGTGA', 'GTGAGT', 'GAGTAG',
            'AGTAGT', 'TAGTGA', 'GTGAAC', 'GAACCG',
        ],
        # Genomic 3' UTR pattern (from reference genome R64-5-1)
        # Starts immediately after stop codon (TAA at 452638-452640)
        # Reference: ref|NC_001140|:452641-452670 = AGTGCTTTTAACTAAGAATTATTAGTCTTT
        'genomic_3utr_full': 'AGTGCTTTTAACTAAGAATTATTAGTCTTT',
        # Multiple 6bp cores for genomic pattern (corrected)
        'genomic_3utr_cores': [
            'AGTGCT', 'TGCTTT', 'GCTTTT', 'CTTTTA',
            'TTTAAC', 'TAACTA', 'AACTAA', 'ACTAAG',
        ],
        # Minimum cores required for match (default 2)
        'min_cores_required': 2,
    }
}


class SpikeInInference:
    """
    Infer spike-in signatures from BAM data.

    For genes used as spike-in controls, this class identifies the
    synthetic UTR sequences that distinguish spike-in from endogenous reads.
    """

    def __init__(self, bam_path: str, genome_path: str):
        """
        Initialize spike-in inference.

        Args:
            bam_path: Path to BAM file
            genome_path: Path to reference genome FASTA
        """
        self.bam_path = bam_path
        self.genome_path = genome_path

    def get_genomic_sequence(self, chrom: str, start: int, end: int) -> str:
        """Extract genomic sequence from reference."""
        try:
            cmd = f"samtools faidx {self.genome_path} '{chrom}:{start}-{end}' 2>/dev/null"
            result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, timeout=30)
            lines = result.stdout.decode().strip().split('\n')
            return ''.join(lines[1:]).upper()
        except (subprocess.TimeoutExpired, subprocess.SubprocessError, UnicodeDecodeError) as e:
            # Silent failure is intentional - caller handles empty sequence
            return ''

    def infer_spikein_signature(
        self,
        gene_name: str,
        chrom: str,
        gene_start: int,
        gene_end: int,
        cds_end: int,
        strand: str = '+',
        max_reads: int = 50000,
        min_pattern_count: int = 100
    ) -> Optional[Dict]:
        """
        Infer spike-in signature for a gene.

        Strategy:
        1. Get genomic 3' UTR sequence (after CDS end)
        2. Extract reads from gene region
        3. Find sequences following CDS end in reads
        4. Compare to genomic - the dominant non-genomic pattern is synthetic

        Args:
            gene_name: Gene name (e.g., 'ENO2')
            chrom: Chromosome name
            gene_start: Gene start position
            gene_end: Gene end position
            cds_end: CDS end position
            strand: Strand ('+' or '-')
            max_reads: Maximum reads to analyze
            min_pattern_count: Minimum count to consider a pattern significant

        Returns:
            Dict with spike-in signature, or None if not found
        """
        print(f"Inferring spike-in signature for {gene_name}...")

        # Get genomic 3' UTR sequence (20bp after CDS end)
        if strand == '+':
            genomic_3utr = self.get_genomic_sequence(chrom, cds_end, cds_end + 30)
        else:
            genomic_3utr = self.get_genomic_sequence(chrom, cds_end - 30, cds_end)
            # Reverse complement for minus strand
            genomic_3utr = self._revcomp(genomic_3utr)

        if not genomic_3utr:
            print(f"  ERROR: Could not get genomic sequence")
            return None

        print(f"  Genomic 3' UTR: {genomic_3utr[:20]}...")

        # Get CDS end motif (for finding position in reads)
        if strand == '+':
            cds_end_motif = self.get_genomic_sequence(chrom, cds_end - 10, cds_end)
        else:
            cds_end_motif = self.get_genomic_sequence(chrom, cds_end, cds_end + 10)
            cds_end_motif = self._revcomp(cds_end_motif)

        print(f"  CDS end motif: {cds_end_motif}")

        # Extract patterns following CDS end in reads
        patterns_with_genomic = []
        patterns_without_genomic = []
        total_reads = 0

        bam = pysam.AlignmentFile(self.bam_path, 'rb')

        for read in bam.fetch(chrom, gene_start, gene_end):
            if total_reads >= max_reads:
                break
            total_reads += 1

            seq = read.query_sequence
            if not seq:
                continue

            # Find CDS end motif and extract following sequence
            # Use a shorter CDS end motif for more matches
            short_motif = cds_end_motif[-8:] if len(cds_end_motif) >= 8 else cds_end_motif
            match = re.search(rf'{short_motif}([ACGT]{{30,50}})', seq)

            if match:
                following = match.group(1)[:30]

                # Check if it matches genomic
                if genomic_3utr[:15] in following:
                    patterns_with_genomic.append(following)
                else:
                    patterns_without_genomic.append(following)

        bam.close()

        print(f"  Analyzed {total_reads:,} reads")
        print(f"  Patterns with genomic 3' UTR: {len(patterns_with_genomic):,}")
        print(f"  Patterns without genomic 3' UTR: {len(patterns_without_genomic):,}")

        if len(patterns_without_genomic) < min_pattern_count:
            print(f"  Not enough non-genomic patterns (need {min_pattern_count})")
            return None

        # Find the most common non-genomic pattern (synthetic 3' UTR)
        pattern_counts = Counter(patterns_without_genomic)
        most_common = pattern_counts.most_common(1)

        if not most_common:
            return None

        synthetic_pattern, count = most_common[0]

        # Extract a core signature (middle 20bp for robustness)
        if len(synthetic_pattern) >= 25:
            synthetic_core = synthetic_pattern[5:25]
        else:
            synthetic_core = synthetic_pattern[:20]

        print(f"  Synthetic 3' UTR pattern: {synthetic_pattern}")
        print(f"  Synthetic core (for filtering): {synthetic_core}")
        print(f"  Count: {count:,} ({100*count/len(patterns_without_genomic):.1f}% of non-genomic)")

        signature = {
            'gene': gene_name,
            'chrom': chrom,
            'start': gene_start,
            'end': gene_end,
            'cds_end': cds_end,
            'strand': strand,
            'genomic_3utr': genomic_3utr[:20],
            'synthetic_3utr_full': synthetic_pattern,
            'synthetic_3utr_core': synthetic_core,
            'total_reads': total_reads,
            'spikein_reads': len(patterns_without_genomic),
            'endogenous_reads': len(patterns_with_genomic),
            'confidence': count / len(patterns_without_genomic) if patterns_without_genomic else 0
        }

        return signature

    def _revcomp(self, seq: str) -> str:
        """Reverse complement a DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(b, 'N') for b in reversed(seq.upper()))


class SpikeInFilter:
    """
    Filter spike-in reads from BAM files.

    Uses sequence-based filtering to identify spike-in reads by their
    synthetic UTR signatures. This correctly preserves endogenous reads
    even when they map to the same genomic region as the spike-in.

    Uses fast multi-core exact matching by default, with optional
    edit distance fallback for ambiguous reads.
    """

    def __init__(
        self,
        signatures: Optional[List[Dict]] = None,
        use_edit_distance: bool = False
    ):
        """
        Initialize spike-in filter.

        Args:
            signatures: List of spike-in signature dicts
            use_edit_distance: Use slow edit distance matching (default: False)
        """
        self.signatures = signatures or []
        self.use_edit_distance = use_edit_distance

    def add_known_signature(self, gene_name: str):
        """Add a pre-defined spike-in signature by gene name."""
        if gene_name in KNOWN_SPIKEIN_SIGNATURES:
            self.signatures.append(KNOWN_SPIKEIN_SIGNATURES[gene_name])

    def add_signature(self, signature: Dict):
        """Add a custom spike-in signature."""
        self.signatures.append(signature)

    def extract_3utr_from_read(self, read: pysam.AlignedSegment) -> Optional[str]:
        """
        Extract the 3' UTR region from a read's soft-clipped sequence.

        For direct RNA-seq of plus strand genes, the 3' UTR appears in
        the RIGHT soft-clip (before poly-A tail).

        Args:
            read: pysam AlignedSegment

        Returns:
            3' UTR sequence (portion before poly-A), or None if not found
        """
        seq = read.query_sequence
        if not seq:
            return None

        cigar = read.cigartuples
        if not cigar:
            return None

        # Check for right soft-clip (operation 4 = S)
        if cigar[-1][0] == 4:  # S = soft clip
            right_clip_len = cigar[-1][1]
            if right_clip_len >= 20:
                right_clip = seq[-right_clip_len:]

                # Check if it contains poly-A
                a_count = right_clip.count('A')
                if a_count / len(right_clip) > 0.5:
                    # Find where poly-A starts (consecutive A's)
                    polya_start = len(right_clip)
                    for i in range(len(right_clip)):
                        # Count consecutive A's from this position
                        consecutive = 0
                        for j in range(i, len(right_clip)):
                            if right_clip[j] == 'A':
                                consecutive += 1
                            else:
                                break
                        if consecutive >= 8:
                            polya_start = i
                            break

                    # Extract 3' UTR (before poly-A)
                    utr = right_clip[:polya_start]
                    if len(utr) >= 10:
                        return utr

        return None

    def classify_read_fast(self, seq: str) -> Tuple[str, str]:
        """
        Fast classification using multiple short exact matches.

        Classification requires at least min_cores (default 2) matches to avoid
        false positives from random 6bp matches. A single 6bp core can match
        by chance, especially near poly(A) tails.

        Classification logic:
        1. Strong match: ≥min_cores of one type, <min_cores of other → classify
        2. Both strong: use majority vote
        3. Neither strong: ambiguous (avoid false positives)

        Args:
            seq: Read sequence string

        Returns:
            Tuple of (classification, gene_name)
            - classification: 'spikein', 'endogenous', or 'ambiguous'
            - gene_name: Name of matched gene (or '')
        """
        if not seq:
            return 'ambiguous', ''

        for sig in self.signatures:
            gene_name = sig.get('gene', 'unknown')
            min_cores = sig.get('min_cores_required', 2)

            # Get cores (use multiple short patterns for robustness)
            syn_cores = sig.get('synthetic_3utr_cores', [])
            gen_cores = sig.get('genomic_3utr_cores', [])

            # Fallback: generate 6bp cores from full sequence
            if not syn_cores and sig.get('synthetic_3utr_full'):
                full = sig['synthetic_3utr_full']
                syn_cores = [full[i:i+6] for i in range(0, len(full)-5, 3)]
            if not gen_cores and sig.get('genomic_3utr_full'):
                full = sig['genomic_3utr_full']
                gen_cores = [full[i:i+6] for i in range(0, len(full)-5, 3)]

            # Count matching cores
            syn_matches = sum(1 for core in syn_cores if core in seq)
            gen_matches = sum(1 for core in gen_cores if core in seq)

            # Require at least min_cores matches to classify
            # Single 6bp matches can occur by chance - don't use them
            has_strong_synthetic = syn_matches >= min_cores
            has_strong_genomic = gen_matches >= min_cores

            if has_strong_synthetic and not has_strong_genomic:
                return 'spikein', gene_name
            elif has_strong_genomic and not has_strong_synthetic:
                return 'endogenous', gene_name
            elif has_strong_synthetic and has_strong_genomic:
                # Both strong - use majority vote
                if syn_matches > gen_matches:
                    return 'spikein', gene_name
                elif gen_matches > syn_matches:
                    return 'endogenous', gene_name
                else:
                    return 'ambiguous', gene_name

            # Neither has enough matches - ambiguous
            # Do NOT classify based on single core matches (too many false positives)

        return 'ambiguous', ''

    def classify_read_by_edit_distance(
        self,
        read: pysam.AlignedSegment
    ) -> Tuple[str, str, float]:
        """
        Slow but precise classification using edit distance.

        Use this for edge cases or validation. For bulk filtering,
        use classify_read_fast() instead.

        Args:
            read: pysam AlignedSegment

        Returns:
            Tuple of (classification, gene_name, confidence)
        """
        seq = read.query_sequence
        if not seq:
            return 'ambiguous', '', 0.0

        for sig in self.signatures:
            synthetic = sig.get('synthetic_3utr_full', '')
            genomic = sig.get('genomic_3utr_full', '')
            gene_name = sig.get('gene', 'unknown')

            if not synthetic or not genomic:
                continue

            # Find best match positions
            _, dist_synthetic, _ = find_best_match_position(seq, synthetic)
            _, dist_genomic, _ = find_best_match_position(seq, genomic)

            # Normalize
            pattern_len = len(synthetic)
            norm_synthetic = dist_synthetic / pattern_len
            norm_genomic = dist_genomic / pattern_len

            # Classify
            if min(norm_synthetic, norm_genomic) > 0.5:
                continue

            confidence = abs(norm_genomic - norm_synthetic)

            if norm_synthetic < norm_genomic:
                return 'spikein', gene_name, confidence
            elif norm_genomic < norm_synthetic:
                return 'endogenous', gene_name, confidence
            else:
                return 'ambiguous', gene_name, 0.0

        return 'ambiguous', '', 0.0

    def is_spikein_read(self, read: pysam.AlignedSegment) -> Tuple[bool, str]:
        """
        Check if a read is from a spike-in control RNA.

        Uses the 3' UTR region (soft-clipped sequence before poly-A) for
        accurate classification. Falls back to full sequence if no soft-clip.

        Args:
            read: pysam AlignedSegment

        Returns:
            Tuple of (is_spikein, spikein_name)
        """
        classification, gene_name = self.classify_read_by_3utr(read)
        return classification == 'spikein', gene_name

    def classify_read(self, read: pysam.AlignedSegment) -> str:
        """
        Classify a read as spike-in, endogenous, or ambiguous.

        Uses the 3' UTR region (soft-clipped sequence before poly-A) for
        accurate classification. Falls back to full sequence if no soft-clip.

        Args:
            read: pysam AlignedSegment

        Returns:
            Classification: 'spikein', 'endogenous', or 'ambiguous'
        """
        classification, _ = self.classify_read_by_3utr(read)
        return classification

    def classify_read_by_3utr(self, read: pysam.AlignedSegment) -> Tuple[str, str]:
        """
        Classify a read using soft-clip analysis.

        Classification logic:
        1. Long soft-clip with synthetic 3' UTR + poly-A → spike-in
        2. Long soft-clip with genomic 3' UTR + poly-A → endogenous
        3. Short soft-clip with only poly-A → endogenous (3' UTR aligned to reference)
        4. No poly-A soft-clip → ambiguous (internal read)

        Args:
            read: pysam AlignedSegment

        Returns:
            Tuple of (classification, gene_name)
        """
        seq = read.query_sequence
        if not seq:
            return 'ambiguous', ''

        cigar = read.cigartuples
        if not cigar:
            return 'ambiguous', ''

        # Check for right soft-clip
        if cigar[-1][0] != 4:  # No soft-clip
            return 'ambiguous', ''

        right_clip_len = cigar[-1][1]
        if right_clip_len < 5:
            return 'ambiguous', ''

        right_clip = seq[-right_clip_len:]

        # Check A content (poly-A indicator)
        a_frac = right_clip.count('A') / len(right_clip)

        if a_frac < 0.5:
            # Not a poly-A clip, likely internal read
            return 'ambiguous', ''

        # Try to extract 3' UTR (sequence before poly-A)
        utr_seq = self.extract_3utr_from_read(read)

        if utr_seq and len(utr_seq) >= 15:
            # Long soft-clip with 3' UTR → classify based on 3' UTR sequence
            return self.classify_read_fast(utr_seq)
        elif right_clip_len <= 20 and a_frac > 0.7:
            # Short soft-clip with mostly poly-A → endogenous
            # (3' UTR aligned to reference, only poly-A is soft-clipped)
            # Return first gene signature as context
            gene_name = self.signatures[0].get('gene', 'unknown') if self.signatures else ''
            return 'endogenous', gene_name
        else:
            # Intermediate case - check for synthetic cores in the clip
            # Require >=2 matches to avoid false positives from random 6bp matches
            min_cores = 2
            for sig in self.signatures:
                gene_name = sig.get('gene', 'unknown')
                syn_cores = sig.get('synthetic_3utr_cores', [])
                gen_cores = sig.get('genomic_3utr_cores', [])

                syn_matches = sum(1 for core in syn_cores if core in right_clip)
                gen_matches = sum(1 for core in gen_cores if core in right_clip)

                # Only classify as spike-in if strong synthetic signal and weak genomic
                if syn_matches >= min_cores and gen_matches < min_cores:
                    return 'spikein', gene_name
            # No strong synthetic pattern found - classify as endogenous
            gene_name = self.signatures[0].get('gene', 'unknown') if self.signatures else ''
            return 'endogenous', gene_name

    def filter_bam(
        self,
        input_bam: str,
        output_bam: str,
        report_path: Optional[str] = None,
        keep_ambiguous: bool = True
    ) -> Dict:
        """
        Filter spike-in reads from a BAM file.

        Args:
            input_bam: Input BAM file path
            output_bam: Output BAM file path (filtered)
            report_path: Optional path for filtering report
            keep_ambiguous: If True (default), keep ambiguous reads; only confirmed spike-ins are removed

        Returns:
            Dict with filtering statistics
        """
        stats = {
            'total_reads': 0,
            'spikein_reads': 0,
            'endogenous_reads': 0,
            'ambiguous_reads': 0,
            'kept_reads': 0,
            'by_gene': {}
        }

        for sig in self.signatures:
            stats['by_gene'][sig.get('gene', 'unknown')] = 0

        bam_in = pysam.AlignmentFile(input_bam, 'rb')
        bam_out = pysam.AlignmentFile(output_bam, 'wb', template=bam_in)

        for read in bam_in:
            stats['total_reads'] += 1

            seq = read.query_sequence
            classification, gene_name = self.classify_read_fast(seq)

            if classification == 'spikein':
                stats['spikein_reads'] += 1
                if gene_name in stats['by_gene']:
                    stats['by_gene'][gene_name] += 1
            elif classification == 'endogenous':
                stats['endogenous_reads'] += 1
                bam_out.write(read)
                stats['kept_reads'] += 1
            else:  # ambiguous
                stats['ambiguous_reads'] += 1
                if keep_ambiguous:
                    bam_out.write(read)
                    stats['kept_reads'] += 1

            if stats['total_reads'] % 100000 == 0:
                print(f"  Processed {stats['total_reads']:,} reads, "
                      f"filtered {stats['spikein_reads']:,} spike-in, "
                      f"{stats['ambiguous_reads']:,} ambiguous")

        bam_in.close()
        bam_out.close()

        # Index output BAM
        pysam.index(output_bam)

        # Write report
        if report_path:
            self._write_report(stats, input_bam, output_bam, report_path)

        return stats

    def _write_report(
        self,
        stats: Dict,
        input_bam: str,
        output_bam: str,
        report_path: str
    ):
        """Write filtering report to file."""
        with open(report_path, 'w') as f:
            f.write("=== Spike-in Filtering Report ===\n\n")
            f.write(f"Input BAM: {input_bam}\n")
            f.write(f"Output BAM: {output_bam}\n\n")

            f.write("Signatures used:\n")
            for sig in self.signatures:
                cores = sig.get('synthetic_3utr_cores', [])
                core_str = sig.get('synthetic_3utr_full', ', '.join(cores[:3]) if cores else 'N/A')
                f.write(f"  {sig.get('gene', 'unknown')}: {core_str}\n")

            f.write(f"\nTotal reads: {stats['total_reads']:,}\n")
            if stats['total_reads'] > 0:
                total = stats['total_reads']
                f.write(f"  Spike-in (removed):   {stats['spikein_reads']:,} "
                        f"({100*stats['spikein_reads']/total:.1f}%)\n")
                f.write(f"  Endogenous (kept):    {stats['endogenous_reads']:,} "
                        f"({100*stats['endogenous_reads']/total:.1f}%)\n")
                f.write(f"  Ambiguous:            {stats['ambiguous_reads']:,} "
                        f"({100*stats['ambiguous_reads']/total:.1f}%)\n")
                f.write(f"\nReads in output BAM: {stats['kept_reads']:,} "
                        f"({100*stats['kept_reads']/total:.1f}%)\n")

            f.write("\nSpike-in reads by gene:\n")
            for gene, count in stats['by_gene'].items():
                f.write(f"  {gene}: {count:,}\n")


def infer_spikein_from_bam(
    bam_path: str,
    genome_path: str,
    gene_name: str,
    chrom: str,
    gene_start: int,
    gene_end: int,
    cds_end: int,
    strand: str = '+',
    output_json: Optional[str] = None
) -> Optional[Dict]:
    """
    Convenience function to infer spike-in signature from BAM.

    Args:
        bam_path: Path to BAM file
        genome_path: Path to reference genome FASTA
        gene_name: Gene name (e.g., 'ENO2')
        chrom: Chromosome name
        gene_start: Gene start position
        gene_end: Gene end position
        cds_end: CDS end position
        strand: Strand ('+' or '-')
        output_json: Optional path to save signature as JSON

    Returns:
        Dict with spike-in signature, or None if not found
    """
    inferrer = SpikeInInference(bam_path, genome_path)
    signature = inferrer.infer_spikein_signature(
        gene_name=gene_name,
        chrom=chrom,
        gene_start=gene_start,
        gene_end=gene_end,
        cds_end=cds_end,
        strand=strand
    )

    if signature and output_json:
        import json
        with open(output_json, 'w') as f:
            json.dump(signature, f, indent=2)
        print(f"Saved signature to: {output_json}")

    return signature


def filter_spikein_reads(
    input_bam: str,
    output_bam: str,
    signatures: Optional[List[Dict]] = None,
    known_genes: Optional[List[str]] = None,
    report_path: Optional[str] = None,
    keep_ambiguous: bool = True
) -> Dict:
    """
    Convenience function to filter spike-in reads from BAM.

    Args:
        input_bam: Input BAM file path
        output_bam: Output BAM file path
        signatures: List of spike-in signature dicts
        known_genes: List of known spike-in gene names (e.g., ['ENO2'])
        report_path: Optional path for filtering report
        keep_ambiguous: If True (default), keep ambiguous reads; only confirmed spike-ins are removed

    Returns:
        Dict with filtering statistics
    """
    filt = SpikeInFilter()

    # Add known signatures
    if known_genes:
        for gene in known_genes:
            filt.add_known_signature(gene)

    # Add custom signatures
    if signatures:
        for sig in signatures:
            filt.add_signature(sig)

    return filt.filter_bam(input_bam, output_bam, report_path, keep_ambiguous=keep_ambiguous)
