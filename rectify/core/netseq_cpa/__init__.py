"""RECTIFY NET-seq CPA-intermediate capture arm.

Recovers cleavage-and-polyadenylation (CPA) intermediates from NET-seq /
mNET-seq libraries — Pol II-associated nascent RNA 3' ends that map to poly-A
sites and carry a non-templated poly(A) tail. Standard NET-seq pipelines
silently discard the evidence: poly-A-dominated reads fail the aligner's
min-ratio filter and land in the unmapped pile (27-98x poly-A-enriched). This
arm recovers them via a two-tier strategy:

- **Tier 1** (:mod:`rectify.core.netseq_cpa.pileup`): genome-aware,
  junction-anchored poly-A metric on the *mapped* reads — the non-circular
  discovery + concordance tier. NET-seq antisense geometry (poly-A on the read
  5' side) is identical to QuantSeq REV; the per-read walkback delegates to the
  validated :func:`rectify.core.correct.protocols.netseq.walkback_netseq`.
- **Rescue** (:mod:`rectify.core.netseq_cpa.rescue`): the poly-A-dominated
  *unmapped* reads — trim the 5' poly-T (= RNA poly-A in the antisense cDNA),
  re-map the genomic anchor, read the CPA off the poly-T-side end. 78-89% of
  recovered anchors land on known DRS CPAs (vs ~5% genome null) across 4 labs.
- **Tier 2** (:mod:`rectify.core.netseq_cpa.cpa_window_rescue`): residual
  sub-uniqueness anchors aligned against a CPA-window mini-reference only.
  CIRCULAR for site discovery — for poly-A LENGTH characterization at
  already-established CPAs ONLY; tagged distinctly, never pooled into
  discovery/concordance numbers.

Paired-end mNET-seq is reduced to read 2 (the nascent 3' end read) upstream by
:mod:`rectify.core.netseq_cpa.pe_reduce`.

Concordance against the orthogonal DRS CPA map is the validation; always
reported against the genome null. See ``analyses/netseq_oligoA_cpa_20260603/``.
"""
from __future__ import annotations

# Yeast S. cerevisiae R64 RefSeq accession -> chrI-style name (lab convention).
# Pass to pileup/rescue/concordance when the BAM was aligned to the NCBI-named
# genome; omit (None -> identity) when reference names are already chrI-style.
NCBI_TO_CHROM = {
    f"ref|NC_{i:06d}|": f"chr{r}"
    for i, r in zip(
        range(1133, 1149),
        ["I", "II", "III", "IV", "V", "VI", "VII", "VIII",
         "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"],
    )
}
NCBI_TO_CHROM["ref|NC_001224|"] = "chrMito"

# S. cerevisiae R64 nuclear genome length (excl. mito); genome-null denominator
# = NUCLEAR_BP * 2 strands. Override for other organisms.
SCER_NUCLEAR_BP = 12_071_326

__all__ = ["NCBI_TO_CHROM", "SCER_NUCLEAR_BP"]
