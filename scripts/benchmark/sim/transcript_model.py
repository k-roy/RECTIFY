#!/usr/bin/env python3
"""Transcript<->genome model: the truth-propagation backbone (component 1 core).

Junction / CPA / NIC-NNC truth is DERIVED FROM THE CONSTRUCTION, not from the
simulator (the simulator only reports read->transcript-of-origin). A
``TranscriptModel`` carries the exon structure relative to the genome, so:

* the spliced transcript sequence is built from the genome exons;
* every intron yields a ``JunctionTruth`` (donor/acceptor, canonical?, motif),
  classified ANNOTATED / NIC / NNC against the supplied annotation;
* a per-base map ``transcript_to_genome`` lets the pbsim3 MAF (read<->transcript)
  be COMPOSED into an exact read<->genome alignment — the per-read ground-truth
  alignment the framing metric (exact-indel concordance) requires.

This module is the shared backbone of both benchmark tiers:
* Tier-1 controlled generators (``controlled.py``) build small ``TranscriptModel``s
  with engineered ambiguity/HP/STR strata.
* Tier-2 pbsim3 wrapper (``pbsim3_wrapper.py``) simulates reads from the model's
  spliced FASTA and projects the MAF back through ``transcript_to_genome``.

Author: Kevin R. Roy
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
from rectify.core.benchmark.truth_schema import (  # noqa: E402
    JunctionTruth,
    JunctionClass,
)

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


@dataclass
class Exon:
    start: int   # 0-based, half-open, genome coords (plus-strand orientation)
    end: int


@dataclass
class TranscriptModel:
    """A transcript laid out on the genome by its exon blocks.

    Exons are stored in GENOME (plus-strand) coordinate order regardless of
    strand. The spliced transcript runs 5'->3' (so for a '-' strand transcript
    the spliced sequence is the reverse complement of the concatenated genome
    exons). Junction/indel truth is always reported in GENOME coordinates (the
    coordinate system the aligner BAM uses), matching the scorer.
    """

    name: str
    chrom: str
    strand: str
    exons: List[Exon]          # genome-sorted, ascending
    genome_seq: str            # full contig sequence (for motif/ambiguity calc)

    def __post_init__(self):
        self.exons = sorted(self.exons, key=lambda e: e.start)

    # ---- sequence -------------------------------------------------------
    def genomic_spliced(self) -> str:
        """Spliced sequence in GENOME (plus-strand) orientation (exons joined)."""
        return "".join(self.genome_seq[e.start:e.end] for e in self.exons)

    def spliced_transcript(self) -> str:
        """Mature transcript 5'->3' (reverse-complemented for '-' strand)."""
        g = self.genomic_spliced()
        return g if self.strand == "+" else revcomp(g)

    @property
    def genome_start(self) -> int:
        return self.exons[0].start

    @property
    def genome_end(self) -> int:
        return self.exons[-1].end

    # ---- coordinate maps ------------------------------------------------
    def genomic_offset_to_genome(self, off: int) -> int:
        """Map an offset into the *genomic-spliced* string to a genome position."""
        rem = off
        for e in self.exons:
            ln = e.end - e.start
            if rem < ln:
                return e.start + rem
            rem -= ln
        return self.exons[-1].end  # clamp

    def transcript_pos_to_genomic_offset(self, tpos: int) -> int:
        """Map a transcript (5'->3') position to a genomic-spliced offset."""
        total = sum(e.end - e.start for e in self.exons)
        if self.strand == "+":
            return tpos
        return total - 1 - tpos  # '-' strand: 5'->3' reverses the genomic order

    def transcript_pos_to_genome(self, tpos: int) -> int:
        return self.genomic_offset_to_genome(self.transcript_pos_to_genomic_offset(tpos))

    # ---- introns / junctions -------------------------------------------
    def introns(self) -> List[Tuple[int, int]]:
        """Intron spans ``[start, end)`` in genome coords (between adjacent exons)."""
        out = []
        for a, b in zip(self.exons, self.exons[1:]):
            if b.start > a.end:
                out.append((a.end, b.start))
        return out

    def junction_truths(self,
                        annotated_pairs: Optional[Set[Tuple[str, int, int]]] = None,
                        annotated_donors: Optional[Set[Tuple[str, int]]] = None,
                        annotated_acceptors: Optional[Set[Tuple[str, int]]] = None,
                        ) -> List[JunctionTruth]:
        """Build ``JunctionTruth`` per intron, classified ANNOTATED/NIC/NNC.

        * ANNOTATED: the (chrom, start, end) intron is in ``annotated_pairs``.
        * NIC: BOTH the donor (chrom, start) AND acceptor (chrom, end) are
          individually catalogued (in ``annotated_donors``/``annotated_acceptors``)
          but the pairing is not — a novel isoform from known sites.
        * NNC: at least one site is uncatalogued.

        Donor = intron_start (genome), acceptor = intron_end (genome). On '-'
        strand the biological donor is at intron_end, but the catalogue/coords
        are genome-based so membership tests stay genome-based; the NIC/NNC logic
        is symmetric in (start, end).
        """
        annotated_pairs = annotated_pairs or set()
        annotated_donors = annotated_donors or set()
        annotated_acceptors = annotated_acceptors or set()
        out = []
        for (s, e) in self.introns():
            if (self.chrom, s, e) in annotated_pairs:
                klass = JunctionClass.ANNOTATED
            else:
                d_known = (self.chrom, s) in annotated_donors
                a_known = (self.chrom, e) in annotated_acceptors
                klass = JunctionClass.NIC if (d_known and a_known) else JunctionClass.NNC
            out.append(JunctionTruth.from_intron(s, e, self.strand, klass, self.genome_seq))
        return out

    # ---- full-length read truth (no internal indel; junction round-trip) -
    def fulllength_cigar(self) -> str:
        """SAM CIGAR for a perfect full-length read of this transcript aligned to
        the genome: ``<exon>M (<intron>N <exon>M)*``."""
        ops = []
        for i, e in enumerate(self.exons):
            ops.append(f"{e.end - e.start}M")
            if i + 1 < len(self.exons):
                gap = self.exons[i + 1].start - e.end
                if gap > 0:
                    ops.append(f"{gap}N")
        return "".join(ops)
