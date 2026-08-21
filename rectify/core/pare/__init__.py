"""RECTIFY PARE (5'-monophosphate degradome) arm.

PARE / degradome sequencing captures RNA 5'-monophosphate ends: **read base 1 =
the 5'P cut site** (endonucleolytic cleavage, decapping, or an exonuclease
stall), and the read extends 3'-ward INTO the captured fragment in the RNA's
own orientation — PARE is a SENSE protocol, the 5'-end mirror of NET-seq.

Because PARE inserts are short (the whole fragment usually fits in one read), a
fragment whose 3' end is POLYADENYLATED leaves a non-templated A-run between
the genomic anchor and the 3' adapter. End-to-end aligners (bwa aln) reject
those reads outright and standard 5'P pipelines never see them; a local aligner
soft-clips the A-run instead. This arm recovers the evidence two ways, exactly
as the validated NET-seq CPA arm does (:mod:`rectify.core.netseq_cpa`):

- **Tier 1** (:mod:`rectify.core.pare.pileup`): one pass over the mapped reads
  producing (a) the strand-aware 5'P cut-site pileup with the MANDATORY 5'-clip
  gate + clipped-fraction report, (b) a genome-aware junction-anchored poly(A)
  pileup at the RNA 3' side (walkback + soft-clip metrics imported unchanged
  from ``netseq_cpa.pileup`` — only the sense geometry differs), and (c) the
  mapped read-length census.
- **Rescue** (:mod:`rectify.core.pare.rescue`): poly-A-dominated UNMAPPED reads
  — trim the trailing 3' poly-A (sense), re-map the genomic anchor, read the
  CPA junction off the poly-A-side end. The anchor's OTHER end is the 5'P cut
  site, so every rescued read yields a single-molecule (5'P, pA-junction) pair
  — the direct evidence that an endonucleolytic cut and a polyadenylated 3'
  end coexisted on one short RNA.

Concordance against the orthogonal DRS CPA map reuses
:mod:`rectify.core.netseq_cpa.concordance` unchanged (shared pileup columns).
"""
from __future__ import annotations

# Chromosome-name translation + genome-null constant are shared with the
# NET-seq CPA arm; import from there so the two arms can never drift.
from ..netseq_cpa import NCBI_TO_CHROM, SCER_NUCLEAR_BP

__all__ = ["NCBI_TO_CHROM", "SCER_NUCLEAR_BP"]
