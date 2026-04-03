"""
Splice site motif scoring and non-canonical splice site detection.

Ported from nanocompass src/sequence_analysis.py and src/noncanonical_ss_detection.py.

Author: Kevin R. Roy
"""

from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import logging

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# SpliceMotifScorer
# ---------------------------------------------------------------------------

class SpliceMotifScorer:
    """Score splice site sequences against consensus motifs.

    Supports organism-specific consensus sequences and position-specific
    penalty matrices.
    """

    # Yeast (S. cerevisiae) defaults — highly constrained splice sites
    YEAST_DEFAULTS = dict(
        five_ss_consensus="GTATGT",
        five_ss_penalties=[4.0, 3.0, 1.0, 1.0, 2.0, 1.0],
        three_ss_consensus="YYYAG",
        three_ss_penalties=[1.0, 1.0, 1.0, 3.0, 4.0],
    )

    # Metazoan (human/mouse) defaults — less constrained
    METAZOAN_DEFAULTS = dict(
        five_ss_consensus="GTAAGT",
        five_ss_penalties=[4.0, 3.0, 1.0, 1.0, 1.0, 1.0],
        three_ss_consensus="YYYYAG",
        three_ss_penalties=[1.0, 1.0, 1.0, 1.0, 3.0, 4.0],
    )

    # IUPAC degenerate base codes
    IUPAC_CODES: Dict[str, List[str]] = {
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T'],
    }

    def __init__(
        self,
        five_ss_consensus: str,
        five_ss_penalties: List[float],
        three_ss_consensus: str,
        three_ss_penalties: List[float],
    ):
        self.five_ss_consensus = five_ss_consensus.upper()
        self.five_ss_penalties = five_ss_penalties
        self.three_ss_consensus = three_ss_consensus.upper()
        self.three_ss_penalties = three_ss_penalties

    @classmethod
    def for_yeast(cls) -> 'SpliceMotifScorer':
        return cls(**cls.YEAST_DEFAULTS)

    @classmethod
    def for_metazoan(cls) -> 'SpliceMotifScorer':
        return cls(**cls.METAZOAN_DEFAULTS)

    def _matches(self, observed: str, consensus: str) -> bool:
        if consensus in self.IUPAC_CODES:
            return observed.upper() in self.IUPAC_CODES[consensus]
        return observed.upper() == consensus.upper()

    def score_five_ss(self, sequence: str) -> float:
        """Score 5' splice site sequence (lower = better match, 0 = perfect)."""
        if len(sequence) != len(self.five_ss_consensus):
            return 100.0
        return sum(
            pen
            for (obs, cons), pen in zip(
                zip(sequence.upper(), self.five_ss_consensus),
                self.five_ss_penalties,
            )
            if not self._matches(obs, cons)
        )

    def score_three_ss(self, sequence: str) -> float:
        """Score 3' splice site sequence (lower = better match, 0 = perfect)."""
        if len(sequence) != len(self.three_ss_consensus):
            return 100.0
        return sum(
            pen
            for (obs, cons), pen in zip(
                zip(sequence.upper(), self.three_ss_consensus),
                self.three_ss_penalties,
            )
            if not self._matches(obs, cons)
        )

    def classify_splice_site_type(self, five_dinuc: str, three_dinuc: str) -> str:
        """Classify a splice site pair by dinucleotide type."""
        five = five_dinuc.upper()
        three = three_dinuc.upper()
        if five == 'GT' and three == 'AG':
            return 'GT-AG (canonical)'
        if five == 'GC' and three == 'AG':
            return 'GC-AG (non-canonical)'
        if five == 'AT' and three == 'AC':
            return 'AT-AC (minor spliceosome)'
        if three != 'AG':
            if three in ('AA', 'GA'):
                return f'{five}-{three} (non-canonical, RAG-type)'
            if three in ('CG', 'GG', 'TG'):
                return f'{five}-{three} (non-canonical, BG-type)'
            if three in ('AT', 'CT', 'TT'):
                return f'{five}-{three} (non-canonical, HAU-type)'
        return f'{five}-{three} (non-canonical)'


# ---------------------------------------------------------------------------
# Splice site sequence extraction (genome-dict based, no repeated file opens)
# ---------------------------------------------------------------------------

_RC_TABLE = str.maketrans('ACGTacgt', 'TGCAtgca')


def _rc(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def get_splice_site_dinucleotides(
    genome: Dict[str, str],
    chrom: str,
    intron_start: int,
    intron_end: int,
    strand: str,
) -> Tuple[str, str]:
    """Return (five_dinuc, three_dinuc) for a junction in transcript orientation.

    Uses 0-based half-open coordinates: intron spans [intron_start, intron_end).

    Plus strand:  five_dinuc = genome[intron_start : intron_start+2]
                  three_dinuc = genome[intron_end-2 : intron_end]
    Minus strand: five_dinuc = RC(genome[intron_end-2 : intron_end])
                  three_dinuc = RC(genome[intron_start : intron_start+2])
    """
    seq = genome.get(chrom, '')
    if not seq or intron_start < 0 or intron_end > len(seq) or intron_start >= intron_end:
        return ('', '')

    if strand == '+':
        five_dinuc = seq[intron_start : intron_start + 2].upper()
        three_dinuc = seq[intron_end - 2 : intron_end].upper()
    else:
        five_dinuc = _rc(seq[intron_end - 2 : intron_end]).upper()
        three_dinuc = _rc(seq[intron_start : intron_start + 2]).upper()

    return five_dinuc, three_dinuc


def get_splice_site_sequences(
    genome: Dict[str, str],
    chrom: str,
    intron_start: int,
    intron_end: int,
    strand: str,
    context: int = 6,
) -> Tuple[str, str]:
    """Return (five_ss_seq, three_ss_seq) in transcript orientation.

    five_ss_seq:  context bases starting at the intron donor (5'SS)
    three_ss_seq: context bases ending at the intron acceptor (3'SS)
    """
    seq = genome.get(chrom, '')
    if not seq:
        return ('', '')

    if strand == '+':
        five_ss_seq = seq[intron_start : intron_start + context].upper()
        three_ss_seq = seq[max(0, intron_end - context) : intron_end].upper()
    else:
        raw_five = seq[max(0, intron_end - context) : intron_end].upper()
        raw_three = seq[intron_start : intron_start + context].upper()
        five_ss_seq = _rc(raw_five)
        three_ss_seq = _rc(raw_three)

    return five_ss_seq, three_ss_seq


# ---------------------------------------------------------------------------
# NonCanonicalSSDetector
# ---------------------------------------------------------------------------

class NonCanonicalSSDetector:
    """Detect and classify non-canonical splice site dinucleotide pairs."""

    CANONICAL_PAIRS = {
        ('GT', 'AG'),  # U2-type major (~99 %)
        ('GC', 'AG'),  # U2-type minor (~0.5–1 %)
        ('AT', 'AC'),  # U12-type (~0.01 %)
    }

    RARE_NONCANONICAL = {
        ('GT', 'CG'),
        ('GC', 'CG'),
        ('GT', 'TG'),
    }

    def __init__(self, min_support: int = 3, max_error_rate: float = 0.05):
        self.min_support = min_support
        self.max_error_rate = max_error_rate
        self._dinuc_counts: Dict[Tuple[str, str], int] = defaultdict(int)
        self._total = 0

    def is_canonical(self, five_dinuc: str, three_dinuc: str) -> bool:
        return (five_dinuc.upper(), three_dinuc.upper()) in self.CANONICAL_PAIRS

    def classify(self, five_dinuc: str, three_dinuc: str) -> Dict:
        """Classify a dinucleotide pair. Returns a dict with type, subtype, likely_artifact."""
        five = five_dinuc.upper()
        three = three_dinuc.upper()
        pair = (five, three)

        if pair == ('GT', 'AG'):
            return {'type': 'canonical', 'subtype': 'GT-AG (U2-major)',
                    'expected_frequency': 0.99, 'likely_artifact': False}
        if pair == ('GC', 'AG'):
            return {'type': 'canonical', 'subtype': 'GC-AG (U2-minor)',
                    'expected_frequency': 0.008, 'likely_artifact': False}
        if pair == ('AT', 'AC'):
            return {'type': 'canonical', 'subtype': 'AT-AC (U12-type)',
                    'expected_frequency': 0.0001, 'likely_artifact': False}
        if pair in self.RARE_NONCANONICAL:
            return {'type': 'rare_noncanonical', 'subtype': f'{five}-{three}',
                    'expected_frequency': 0.00001, 'likely_artifact': False}
        if 'N' in five or 'N' in three:
            return {'type': 'unknown', 'subtype': f'{five}-{three}',
                    'expected_frequency': 0, 'likely_artifact': True}
        return {'type': 'noncanonical', 'subtype': f'{five}-{three}',
                'expected_frequency': 0, 'likely_artifact': True}

    def record(self, five_dinuc: str, three_dinuc: str) -> None:
        """Record a junction's dinucleotide pair for summary statistics."""
        self._dinuc_counts[(five_dinuc.upper(), three_dinuc.upper())] += 1
        self._total += 1

    def summary(self) -> List[Dict]:
        """Return sorted list of dinucleotide pair frequencies."""
        rows = []
        for (five, three), cnt in sorted(
            self._dinuc_counts.items(), key=lambda x: -x[1]
        ):
            freq = cnt / self._total if self._total else 0
            info = self.classify(five, three)
            rows.append({
                'five_dinuc': five,
                'three_dinuc': three,
                'pair': f'{five}-{three}',
                'count': cnt,
                'frequency': freq,
                **info,
            })
        return rows

    @property
    def total_junctions(self) -> int:
        return self._total
