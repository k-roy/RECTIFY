"""ONT PCR-cDNA protocol wrapper for the read-vs-reference walkback.

ONT PCR-cDNA (SQK-PCB114 and relatives) is built by strand-switching PCR, so the
sequenced molecules are double-stranded and **reads arrive in BOTH orientations**:
roughly a third of the reads are the mRNA sense strand and roughly a third its
reverse complement (the rest carry no tail evidence at all).  This is the one
structural difference from DRS, where the read *is* the RNA and BAM strand is
therefore RNA strand.

Why a fixed rule cannot work
----------------------------
For a gene on the ``+`` strand:

    sense read      -> aligns forward  -> is_reverse=False -> RNA 3' end at reference_end
    antisense read  -> aligns reverse  -> is_reverse=True  -> RNA 3' end at reference_end

and for a gene on the ``-`` strand both orientations put the RNA 3' end at
``reference_start``.  So ``is_reverse`` alone does **not** determine the RNA
strand: gene strand is ``+`` iff ``read_is_sense XOR is_reverse``.  Until
2026-08-01 ``rectify correct --ONT-cDNA`` applied the DRS mapping
(``is_reverse`` -> gene strand) unconditionally, which is right for sense reads
and wrong for antisense ones — measured at ~2/3 of a real library.  See
``planning/541_ont_cdna_strand_fix.md`` and contract C8 in
``planning/518_ms2_metric_contract.md``.

Relationship to the ``correct-cdna`` (UMI-consensus) path
---------------------------------------------------------
RECTIFY has TWO cDNA routes and they must not be confused:

* **Path A (UMI-aware):** ``correct-cdna`` -> ``align -y`` -> ``cdna-analyze``.
  Stage 1 emits one consensus per molecule and tags it ``XO:Z:fwd|rev``
  (plus ``XY:Z:umi_captured_fwd|rev``); ``align -y`` carries those into the BAM
  and ``cdna-analyze`` maps ``{fwd: '+', rev: '-'}`` to get the gene strand.
* **Path B (pre-UMI-collapse reads):** ``trim-cdna-polya`` -> aligner ->
  ``correct --ONT-cDNA``.  This is the route used for internal-poly(A) work,
  where per-read tail/site *distributions* are wanted and UMI collapse is
  deliberately skipped (contract C9).

⚠️ **Path A does NOT normalise orientation** — despite what
``docs/figures/generate_cdna_umi_v3.py`` claimed until 2026-08-01, stage 1
restores the *basecalled* orientation and merely LABELS the molecule with
``XO``; a ``rev`` molecule stays antisense in the FASTQ.  So reads in BOTH
orientations reach the aligner on either path, and ``is_reverse`` is never the
gene strand.  ``read_info.py`` says so explicitly: *"Direction of 'polyA side'
is determined by **orient**, not is_reverse."*

``XO`` is defined on BAM SEQ (reference orientation), whereas ``ro`` below is
defined on the basecalled read; the two are equivalent under
``gene_strand = '+' iff (read_is_sense XOR is_reverse)``.  Both are consumed
here so this function is correct on a BAM from either path.

Resolution order (documented precedence, C8 + 535i section 3.2)
---------------------------------------------------------------
1. ``XO:Z:fwd|rev`` (or ``XY:Z:umi_captured_fwd|rev``) — the canonical
   orientation label written by ``rectify correct-cdna``.  Defined on BAM SEQ:
   ``fwd`` = poly-A at the RIGHT of BAM SEQ = gene '+'; ``rev`` = poly-T at the
   LEFT = gene '-'.  Preferred because it is the established convention and is
   backed by per-molecule consensus rather than a single read.
2. ``ro:A:`` BAM tag — read-intrinsic tail evidence written by
   ``rectify trim-cdna-polya`` and carried through alignment by ``minimap2 -y``.
   ``S`` = 3' poly-A (sense), ``A`` = 5' poly-T (antisense), ``B`` = both tails
   (resolved as **sense**: measured 98.3% agreement with the annotated gene
   strand, indistinguishable from the pure sense class), ``U`` = unresolved.
   Needs no annotation and is independent of alignment.
3. **Maximally-overlapping annotated gene.**  Queried on BOTH strands, so it is
   *not* the circular strand-matched lookup that ``compute_read_gene_attribution``
   performs (that one filters to the alignment strand and can only ever return a
   strand-matched gene).
4. Otherwise ``None`` -> ``unassigned``.  Never guessed.

Once the RNA strand is known, side and stop-base follow the same mapping DRS
uses, because in alignment orientation ``pysam`` has already put the read into
reference orientation: ``+`` -> 3' end at ``reference_end``, poly-A on the RIGHT;
``-`` -> 3' end at ``reference_start``, poly-T on the LEFT.  Only the *strand
determination* differs from DRS, so the walkback core is shared unchanged.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import pysam

from rectify.core.correct.walkback import (
    THREE_PRIME_SIDE_LEFT,
    THREE_PRIME_SIDE_RIGHT,
    _classify_artifact_nops,
    walkback_3prime_guarded,
)

#: BAM aux tag carrying the per-read orientation label written by
#: ``rectify trim-cdna-polya`` (see ``cdna_trim_command``).
ORIENTATION_TAG = "ro"

#: Canonical per-molecule orientation tag written by ``rectify correct-cdna``
#: stage 1 and carried into the BAM by ``rectify align -y`` / ``minimap2 -y``.
#: Defined on BAM SEQ, so it maps straight to gene strand.  ``XY`` carries the
#: same information as ``umi_captured_fwd`` / ``umi_captured_rev``.
CDNA_ORIENT_TAG = "XO"
CDNA_SUBTYPE_TAG = "XY"

#: ``rectify cdna_analyze_command`` uses exactly this mapping; kept identical so
#: the two consumers of ``XO`` cannot drift apart.
ORIENT_TO_STRAND = {"fwd": "+", "rev": "-"}

#: ``strand_evidence`` values emitted into the corrected-3'-ends TSV.
EVIDENCE_XO_ORIENT = "XO_orient"
EVIDENCE_POLYA_3P = "polyA_3p"
EVIDENCE_POLYT_5P = "polyT_5p"
EVIDENCE_GENE_OVERLAP = "gene_overlap"
EVIDENCE_UNASSIGNED = "unassigned"


def _drs_rule_strand(read: pysam.AlignedSegment) -> str:
    """Gene strand under the DRS convention (correct for *sense* cDNA reads)."""
    return "-" if read.is_reverse else "+"


def _flip(strand: str) -> str:
    return "-" if strand == "+" else "+"


def max_overlap_gene_strand(
    read: pysam.AlignedSegment,
    gene_interval_trees: Dict,
    chrom: Optional[str] = None,
) -> Optional[str]:
    """Strand of the annotated gene with the largest overlap of the read body.

    Both strands are queried and the larger overlap wins; ties resolve to the
    DRS-rule strand so the behaviour is deterministic.  Returns ``None`` when the
    read overlaps no annotated feature on either strand.
    """
    if not gene_interval_trees:
        return None
    if chrom is None:
        chrom = read.reference_name
    start, end = read.reference_start, read.reference_end
    if start is None or end is None or end <= start:
        return None

    from rectify.core.analyze.gene_attribution import find_overlapping_genes

    best_strand, best_bp = None, 0
    # Query the DRS-rule strand second so an exact tie keeps it.
    for strand in (_flip(_drs_rule_strand(read)), _drs_rule_strand(read)):
        hits = find_overlapping_genes(chrom, start, end, strand, gene_interval_trees)
        bp = max((h["overlap_bp"] for h in hits), default=0)
        if bp > best_bp:
            best_strand, best_bp = strand, bp
    return best_strand


def resolve_rna_strand(
    read: pysam.AlignedSegment,
    gene_interval_trees: Optional[Dict] = None,
    chrom: Optional[str] = None,
) -> Tuple[Optional[str], str]:
    """Resolve the RNA (gene) strand of one ONT PCR-cDNA read.

    Returns ``(strand, evidence)`` where *strand* is ``'+'``/``'-'`` or ``None``
    when unresolvable, and *evidence* is one of the ``EVIDENCE_*`` constants.
    See the module docstring for the precedence.
    """
    # --- 1. canonical correct-cdna orientation label (XO, else XY) ---
    # Defined on BAM SEQ, so it gives the gene strand directly with no
    # is_reverse arithmetic.  Same mapping cdna_analyze_command.py uses.
    orient = None
    try:
        orient = read.get_tag(CDNA_ORIENT_TAG)
    except KeyError:
        try:
            # XY is 'umi_captured_fwd' / 'umi_captured_rev' / 'umi_not_captured'
            sub = read.get_tag(CDNA_SUBTYPE_TAG)
            if isinstance(sub, str) and sub.startswith("umi_captured_"):
                orient = sub[len("umi_captured_"):]
        except KeyError:
            pass
    if orient in ORIENT_TO_STRAND:
        return ORIENT_TO_STRAND[orient], EVIDENCE_XO_ORIENT

    # --- 2. read-intrinsic tail evidence from the pre-align trim ---
    try:
        label = read.get_tag(ORIENTATION_TAG)
    except KeyError:
        label = None
    if label in ("S", "B"):
        return _drs_rule_strand(read), EVIDENCE_POLYA_3P
    if label == "A":
        return _flip(_drs_rule_strand(read)), EVIDENCE_POLYT_5P

    # --- 3. maximally-overlapping annotated gene (both strands) ---
    if gene_interval_trees:
        gs = max_overlap_gene_strand(read, gene_interval_trees, chrom=chrom)
        if gs is not None:
            return gs, EVIDENCE_GENE_OVERLAP

    # --- 4. never guess ---
    return None, EVIDENCE_UNASSIGNED


def three_prime_position(read: pysam.AlignedSegment, rna_strand: str) -> Optional[int]:
    """RNA 3'-end position implied by *rna_strand*, in alignment coordinates."""
    if rna_strand == "+":
        return (read.reference_end - 1) if read.reference_end is not None else None
    return read.reference_start


def walkback_ont_cdna(
    read: pysam.AlignedSegment,
    chrom_seq: str,
    rna_strand: str,
    *,
    artifact_analyses=None,
    large_del_min_bp: int = 5,
    poly_noise_window: int = 50,
    tail_context_k: int = 4,
    tail_context_far: int = 16,
    max_scan_depth: int = 1000,
    early_exit_min_homopolymer_len: int = 4,
) -> Optional[dict]:
    """Poly-A walkback for ONT PCR-cDNA, driven by an already-resolved strand.

    Identical to :func:`rectify.core.correct.walkback.walkback_drs_full` except
    that side and stop-base come from *rna_strand* (resolved per read by
    :func:`resolve_rna_strand`) instead of from ``is_reverse``.  Returns the same
    dict shape — ``{'corrected_pos', 'original_pos', 'polya_aligned_bp',
    'correction_bp'}`` — or ``None`` when no correction is needed.
    """
    if rna_strand == "-":
        side = THREE_PRIME_SIDE_LEFT
        stop_base = "T"
    else:
        side = THREE_PRIME_SIDE_RIGHT
        stop_base = "A"
    artifact_starts = _classify_artifact_nops(
        read, chrom_seq, rna_strand, artifact_analyses
    )
    return walkback_3prime_guarded(
        read,
        chrom_seq,
        side,
        stop_base=stop_base,
        artifact_n_ref_starts=artifact_starts,
        large_del_min_bp=large_del_min_bp,
        poly_noise_window=poly_noise_window,
        tail_context_k=tail_context_k,
        tail_context_far=tail_context_far,
        max_scan_depth=max_scan_depth,
        early_exit_min_homopolymer_len=early_exit_min_homopolymer_len,
    )
