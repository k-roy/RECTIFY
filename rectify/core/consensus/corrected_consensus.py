"""
corrected_consensus.py — Select winning aligner per read across per-aligner corrected outputs.

For each read, selects the winning aligner's full corrected result from N per-aligner
corrected_reads.tsv files.  The winning aligner's corrected BAM record is then used
downstream for the consensus corrected BAM — meaning the entire read is taken from the
winner (5' junction rescue, intron corrections, and 3' CPA position).

Selection uses HP-aware edit distance computed from the final corrected CIGAR:

  Primary:   hp_edit_distance (ascending) — sum of HP-aware penalties over the
             corrected CIGAR.  Soft-clips count as 1.0/base (flat); intron N-ops
             are free (0); deletions and insertions use the empirical penalty table
             when provided.  Aligners whose corrected reads fit the genome better
             (lower edit distance) win.  5' Cat3 junction rescues naturally emerge
             as winning signals because soft-clipped exon1 bases are converted to
             aligned matches (penalty ≈ 0) by the rescue.

  Tiebreak:  alignment span (descending) — wider alignment preferred when edit
             distances are identical.

Also identifies Cat5 candidates: reads where ≥2 aligners each uniquely contribute
at least one correctly-rescued intron not present in the other aligner's result.
"""

from __future__ import annotations

import contextlib
import logging
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple, Any

import pandas as pd

try:
    import pysam as _pysam
except ImportError:  # pragma: no cover
    _pysam = None  # type: ignore

from rectify.core.splice.calibrate_junction_overhang import OverhangTable, _parse_junctions as _parse_junctions_list
from rectify.core.consensus.consensus import (
    _normalize_bam_read_name,
    _UNDERSCORE_COMMENT_RE,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# HP-aware edit distance from corrected CIGAR
# ---------------------------------------------------------------------------


def _cigar_hp_edit_distance(read, genome: Optional[Dict[str, str]], penalty_table) -> float:
    """Compute HP-aware edit distance for one corrected pysam read.

    CIGAR op costs:
        = (sequence match):   0
        X (mismatch):         1.0 flat
        M (alignment match):  0 if query base == ref base, else 1.0 flat
                              (requires genome for base comparison; falls back
                              to 0 when genome is absent)
        D (deletion):         HP-aware del_cost from penalty_table; fallback 1.0
        I (insertion):        HP-aware ins_cost from penalty_table; fallback 1.25
        N (intron):           0  — free pass
        S (soft-clip):        1.0 per base (flat)
        H (hard-clip):        1.0 per base (flat)

    Soft- and hard-clip are penalised symmetrically. Hard-clips in corrected
    BAMs come from rectify's CIGAR surgery (``clip_read_to_corrected_3prime``)
    when walkback or rescue shifts the 3' end past originally-aligned bases.
    Without this penalty, an aligner that resolved an HP-undercall by
    walkback+hard-clip (e.g. minimap2 ``2S 17M`` → ``5H 14M``) would beat an
    aligner that represented the same artefact in-line as a deletion (e.g.
    uLTRA ``2= 1D 17=``) — even though the latter is more informative and
    places the 3' end at the biologically correct position. Counting H
    against the aligner forces winner-selection to value preserved evidence
    over discarded evidence.
    """
    if read.cigartuples is None:
        return 0.0

    chrom = read.reference_name
    genome_seq: Optional[str] = (genome or {}).get(chrom) if genome else None
    query_seq: Optional[str] = read.query_sequence  # None for unmapped

    ref_pos = read.reference_start
    q_pos = 0
    total = 0.0

    for op, length in read.cigartuples:
        if op == 7:                # = (sequence match)
            ref_pos += length
            q_pos += length
        elif op == 8:              # X (mismatch)
            total += length * 1.0
            ref_pos += length
            q_pos += length
        elif op == 0:              # M (match or mismatch — compare bases)
            if genome_seq and query_seq:
                ref_chunk = genome_seq[ref_pos:ref_pos + length].upper()
                q_chunk   = query_seq[q_pos:q_pos + length].upper()
                # numpy byte comparison — fast for the long M blocks typical in
                # nanopore reads.  Falls back to zip() for short chunks.
                try:
                    import numpy as _np
                    ra = _np.frombuffer(ref_chunk.encode('ascii'), dtype=_np.uint8)
                    qa = _np.frombuffer(q_chunk.encode('ascii'),   dtype=_np.uint8)
                    total += float(_np.count_nonzero(ra != qa))
                except Exception:
                    total += sum(r != q for r, q in zip(ref_chunk, q_chunk))
            ref_pos += length
            q_pos += length
        elif op == 2:              # D (deletion)
            for i in range(length):
                rp = ref_pos + i
                if penalty_table is not None and genome_seq and rp < len(genome_seq):
                    from rectify.core.splice.junction_refiner import _hp_run_length
                    hp = _hp_run_length(genome_seq, rp)
                    total += penalty_table.del_cost(hp, genome_seq[rp])
                else:
                    total += 1.0
            ref_pos += length
        elif op == 1:              # I (insertion — ref_pos does not advance)
            if penalty_table is not None and genome_seq and ref_pos < len(genome_seq):
                from rectify.core.splice.junction_refiner import _hp_run_length
                hp = _hp_run_length(genome_seq, ref_pos)
                total += length * penalty_table.ins_cost(hp, genome_seq[ref_pos])
            else:
                total += length * 1.25
            q_pos += length
        elif op == 3:              # N (intron skip)
            ref_pos += length      # free pass
        elif op == 4:              # S (soft-clip)
            total += length * 1.0
            q_pos += length
        elif op == 5:              # H (hard-clip — query_seq already excludes these)
            total += length * 1.0
        # P (6): no penalty, no advance
    return total


def _cigar_aligned_bases(read) -> int:
    """Count M/=/X CIGAR operations — bases that are truly aligned to the reference.

    Excludes soft-clips (S), hard-clips (H), deletions (D), insertions (I), and
    intron skips (N).  The result is the total span of query bases that are mapped
    to genomic sequence, which is used as a proxy for "decent-complexity mapped
    sequence."  Reads with very few aligned bases (e.g. 40S20M4S → 20 bp) are
    minimal-anchor alignment artifacts common in poly-A/T regions.
    """
    if read.cigartuples is None:
        return 0
    # op codes: 0=M, 7==, 8=X
    return sum(length for op, length in read.cigartuples if op in (0, 7, 8))


_NO_JUNCTION_ANCHOR = 10 ** 9  # sentinel for reads with no N-ops (never gated)


def _cigar_min_junction_anchor(read, genome: Optional[Dict[str, str]]) -> int:
    """Min over all N-ops of the junction-local perfect-match anchor (both flanks).

    For each intron (N-op), walk outward from the junction on each flank and count
    consecutive read bases that exactly match the genome, stopping at the first
    mismatch or non-aligned op.  Returns the minimum of (left, right) anchors across
    every N-op in the read — the gate is only as strong as the weakest flank of the
    weakest junction.  Reads with no N-ops (or no genome) return a large sentinel so
    they are never flagged by the anchor gate.

    A read base of ``=`` is the SAM match-encoding (BBMap/mapPacBio write ``=``/``X``
    CIGAR ops with an ``=``-encoded SEQ); it counts as a match.  For aligners that
    emit literal bases under ``M`` ops, the base is compared to the genome directly.
    Genome is required to resolve literal-base matches; without it the read is not
    gated.
    """
    ct = read.cigartuples
    if not ct:
        return _NO_JUNCTION_ANCHOR
    chrom = read.reference_name
    gseq: Optional[str] = (genome or {}).get(chrom) if genome else None
    query: Optional[str] = read.query_sequence
    if gseq is None or query is None:
        return _NO_JUNCTION_ANCHOR

    def _match(qb: str, gb: str) -> bool:
        return qb == '=' or qb.upper() == gb.upper()

    ref_pos = read.reference_start
    q_pos = 0
    best = _NO_JUNCTION_ANCHOR
    G = len(gseq)
    L = len(query)
    for op, length in ct:
        if op == 3:  # N — measure anchors at (q_pos, ref_pos) | (q_pos, ref_pos+length)
            # left flank: walk back from the base just 5' of the intron
            la = 0
            q = q_pos - 1
            r = ref_pos - 1
            while q >= 0 and r >= 0 and _match(query[q], gseq[r]):
                la += 1; q -= 1; r -= 1
            # right flank: walk forward from the base just 3' of the intron
            ra = 0
            q = q_pos
            r = ref_pos + length
            while q < L and r < G and _match(query[q], gseq[r]):
                ra += 1; q += 1; r += 1
            best = min(best, la, ra)
        # advance positions (ref-consuming: M,D,N,=,X ; query-consuming: M,I,S,=,X)
        if op in (0, 2, 3, 7, 8):
            ref_pos += length
        if op in (0, 1, 4, 7, 8):
            q_pos += length
    return best


def _read_hp_edit_distances(
    bam_path: str,
    genome: Optional[Dict[str, str]] = None,
    penalty_table=None,
) -> Dict[str, Tuple[float, int, int]]:
    """Read a corrected BAM and return {bare_read_id: (hp_edit_distance, aligned_bases)}.

    ``hp_edit_distance`` is the HP-aware edit distance from ``_cigar_hp_edit_distance``.
    ``aligned_bases`` is the count of M/=/X CIGAR bases from ``_cigar_aligned_bases``.
    Both values are returned together to avoid scanning the BAM twice.
    """
    if _pysam is None:
        raise RuntimeError("pysam is required to compute HP edit distances from BAM")
    results: Dict[str, Tuple[float, int, int]] = {}
    try:
        with _pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam:
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                rid = _normalize_bam_read_name(read.query_name or "")
                results[rid] = (
                    _cigar_hp_edit_distance(read, genome, penalty_table),
                    _cigar_aligned_bases(read),
                    _cigar_min_junction_anchor(read, genome),
                )
    except Exception as exc:
        logger.warning("Failed to compute HP edit distances from %s: %s", bam_path, exc)
    return results


def _normalized_correction_lookup(
    corrections: Dict[str, dict],
) -> Dict[str, dict]:
    """Build a normalized-QNAME correction lookup, excluding collisions."""
    norm_counts = Counter(_normalize_bam_read_name(str(read_id)) for read_id in corrections)
    return {
        norm_id: corrections[raw_id]
        for raw_id in corrections
        for norm_id in [_normalize_bam_read_name(str(raw_id))]
        if norm_counts[norm_id] == 1
    }


def _lookup_read_correction(
    read,
    corrections: Dict[str, dict],
    normalized_corrections: Dict[str, dict],
) -> Optional[dict]:
    raw_name = read.query_name or ""
    correction = corrections.get(raw_name)
    if correction is not None:
        return correction
    return normalized_corrections.get(_normalize_bam_read_name(raw_name))


def _as_int(value, default: int = 0) -> int:
    try:
        if pd.isna(value):
            return default
        return int(value)
    except (TypeError, ValueError):
        return default


def _as_str(value) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return str(value)



def _read_hp_edit_distances_from_raw_bam(
    bam_path: str,
    corrected_tsv_path: Path,
    genome: Optional[Dict[str, str]] = None,
    penalty_table=None,
    *,
    strict: bool = True,
    read_id_filter: Optional[Set[str]] = None,
) -> Dict[str, Tuple[float, int, int]]:
    """Compute HP edit distances after applying TSV corrections in memory.

    This is the lazy equivalent of reading a materialized per-aligner corrected
    BAM.  It streams the raw aligner BAM, applies the same corrected-BAM edits
    used by ``write_corrected_bam()``, scores the mutated record, then discards
    it.
    """
    if _pysam is None:
        raise RuntimeError("pysam is required to compute HP edit distances from BAM")

    from rectify.core.bam.bam_writer import (
        _load_corrections_from_tsv,
        apply_corrected_edits_to_read,
    )

    corrections = _load_corrections_from_tsv(str(corrected_tsv_path))
    normalized_corrections = _normalized_correction_lookup(corrections)
    normalized_filter = (
        {_normalize_bam_read_name(str(read_id)) for read_id in read_id_filter}
        if read_id_filter is not None
        else None
    )
    seen_corrected_ids: Set[str] = set()
    results: Dict[str, Tuple[float, int, int]] = {}

    try:
        with _pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam:
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                read_id = _normalize_bam_read_name(read.query_name or "")
                if normalized_filter is not None and read_id not in normalized_filter:
                    continue
                correction = _lookup_read_correction(read, corrections, normalized_corrections)
                if correction is None:
                    continue
                seen_corrected_ids.add(_normalize_bam_read_name(str(correction.get('read_id', read_id))))
                apply_corrected_edits_to_read(read, correction, genome)
                results[read_id] = (
                    _cigar_hp_edit_distance(read, genome, penalty_table),
                    _cigar_aligned_bases(read),
                    _cigar_min_junction_anchor(read, genome),
                )
    except Exception:
        if strict:
            raise
        logger.exception("Failed lazy HP edit-distance scoring from %s", bam_path)

    if strict:
        expected_ids = {
            _normalize_bam_read_name(str(correction.get('read_id', read_id)))
            for read_id, correction in corrections.items()
        }
        if normalized_filter is not None:
            expected_ids &= normalized_filter
        missing = expected_ids - seen_corrected_ids
        if missing:
            sample = ', '.join(sorted(missing)[:5])
            raise RuntimeError(
                f"Lazy HP scoring missed {len(missing)} corrected reads in {bam_path}; "
                f"examples: {sample}"
            )

    logger.info(
        "Lazy HP scoring %s: scored=%d",
        bam_path,
        len(results),
    )
    return results


def _score_hp_distances_for_aligner(args: Tuple[str, str, Optional[str], Optional[Dict[str, str]], Any, bool, Optional[Set[str]]]) -> Tuple[str, Dict[str, Tuple[float, int, int]]]:
    """Process-pool worker for one aligner's HP edit-distance scan."""
    aligner_name, bam_path, tsv_path, genome, penalty_table, strict, read_id_filter = args
    if tsv_path is None:
        if strict:
            raise FileNotFoundError(f"Missing corrected TSV for lazy HP scoring: {aligner_name}")
        logger.warning("Missing corrected TSV for lazy HP scoring, skipping %s", aligner_name)
        return aligner_name, {}
    return aligner_name, _read_hp_edit_distances_from_raw_bam(
        str(bam_path),
        Path(tsv_path),
        genome,
        penalty_table,
        strict=strict,
        read_id_filter=read_id_filter,
    )


def _score_corrected_bam_for_aligner(args: Tuple[str, str, Optional[Dict[str, str]], Any]) -> Tuple[str, Dict[str, Tuple[float, int, int]]]:
    """Process-pool worker for one corrected BAM HP edit-distance scan."""
    aligner_name, bam_path, genome, penalty_table = args
    return aligner_name, _read_hp_edit_distances(str(bam_path), genome, penalty_table)


def write_corrected_consensus_bam(
    per_aligner_raw_bams: Dict[str, str],
    per_aligner_tsvs: Dict[str, Path],
    merged_tsv: Path,
    output_bam: Path,
    genome: Optional[Dict[str, str]] = None,
    *,
    threads: int = 1,
    strict: bool = True,
) -> Dict[str, int]:
    """Write the final corrected consensus BAM directly from raw BAMs + TSVs.

    The merged TSV supplies ``read_id -> winning_aligner``.  For each winning
    raw BAM record, this function applies the aligner's corrected TSV row in
    memory and writes only that corrected winning record to the final BAM.
    """
    if _pysam is None:
        raise RuntimeError("pysam is required to write corrected consensus BAM")

    from rectify.core.bam.bam_writer import (
        _load_corrections_from_tsv,
        apply_corrected_edits_to_read,
    )

    merged_df = _read_corrected_tsv_or_manifest(Path(merged_tsv))
    if 'winning_aligner' not in merged_df.columns:
        raise ValueError(f"{merged_tsv} has no winning_aligner column")
    if 'read_id' not in merged_df.columns:
        raise ValueError(f"{merged_tsv} has no read_id column")

    read_to_aligner = {
        _normalize_bam_read_name(str(row.read_id)): str(row.winning_aligner)
        for row in merged_df[['read_id', 'winning_aligner']].drop_duplicates().itertuples(index=False)
    }
    aligner_to_reads: Dict[str, Set[str]] = defaultdict(set)
    for read_id, aligner in read_to_aligner.items():
        aligner_to_reads[aligner].add(read_id)

    correction_tables: Dict[str, Tuple[Dict[str, dict], Dict[str, dict]]] = {}
    for aligner, tsv_path in per_aligner_tsvs.items():
        corrections = _load_corrections_from_tsv(str(tsv_path))
        correction_tables[aligner] = (
            corrections,
            _normalized_correction_lookup(corrections),
        )

    header_bam = next((Path(p) for p in per_aligner_raw_bams.values() if Path(p).exists()), None)
    if header_bam is None:
        raise FileNotFoundError("No raw BAMs available for corrected consensus output")

    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    unsorted_bam = output_bam.with_name(output_bam.name + '.unsorted.bam')

    written_ids: Set[str] = set()
    stats = {
        'requested': len(read_to_aligner),
        'written': 0,
        'skipped_no_raw_bam': 0,
        'skipped_no_correction': 0,
    }

    with _pysam.AlignmentFile(str(header_bam), 'rb') as header_in:
        header = header_in.header
    with _pysam.AlignmentFile(str(unsorted_bam), 'wb', header=header) as out_bam:
        for aligner, winning_ids in aligner_to_reads.items():
            raw_bam = per_aligner_raw_bams.get(aligner)
            if not raw_bam or not Path(raw_bam).exists():
                stats['skipped_no_raw_bam'] += len(winning_ids)
                if strict:
                    raise FileNotFoundError(f"Missing raw BAM for winning aligner {aligner}: {raw_bam}")
                continue
            if aligner not in correction_tables:
                stats['skipped_no_correction'] += len(winning_ids)
                if strict:
                    raise FileNotFoundError(f"Missing corrected TSV for winning aligner {aligner}")
                continue
            corrections, normalized_corrections = correction_tables[aligner]
            with _pysam.AlignmentFile(str(raw_bam), 'rb') as bam_in:
                for read in bam_in:
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    read_id = _normalize_bam_read_name(read.query_name or "")
                    if read_id not in winning_ids:
                        continue
                    correction = _lookup_read_correction(read, corrections, normalized_corrections)
                    if correction is None:
                        stats['skipped_no_correction'] += 1
                        if strict:
                            raise RuntimeError(
                                f"Missing correction row for winning read {read_id} "
                                f"from aligner {aligner}"
                            )
                        continue
                    apply_corrected_edits_to_read(read, correction, genome)
                    out_bam.write(read)
                    written_ids.add(read_id)
                    stats['written'] += 1

    _pysam.sort('-@', str(max(1, threads)), '-o', str(output_bam), str(unsorted_bam))
    _pysam.index(str(output_bam))
    try:
        unsorted_bam.unlink()
    except OSError:
        pass

    if strict:
        missing = set(read_to_aligner) - written_ids
        if missing:
            sample = ', '.join(sorted(missing)[:5])
            raise RuntimeError(
                f"Corrected consensus BAM missed {len(missing)} winning reads; "
                f"examples: {sample}"
            )

    logger.info(
        "write_corrected_consensus_bam: requested=%d written=%d -> %s",
        stats['requested'],
        stats['written'],
        output_bam,
    )
    return stats


# ---------------------------------------------------------------------------
# UpSet plot helpers
# ---------------------------------------------------------------------------

def _plot_aligner_agreement_upset(
    summary_df: "pd.DataFrame",
    output_path: Path,
) -> None:
    """Write an UpSet plot showing which aligners tied for the lowest HP edit distance per read.

    For each read, aligner X is in the "tied" set if its hp_edit_distance equals
    the winner's — i.e., it achieved the same minimum edit distance before the
    alignment-span tiebreaker resolved the winner.

    Parameters
    ----------
    summary_df:
        DataFrame with columns read_id, _aligner, _is_winner, hp_edit_distance
        (from merge_corrected_tsvs with per_aligner_corrected_bams provided).
        Falls back to a simpler indicator plot if hp_edit_distance is absent.
    output_path:
        Destination .png path.
    """
    try:
        import upsetplot
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("upsetplot or matplotlib not available — skipping UpSet plot")
        return

    if 'hp_edit_distance' in summary_df.columns:
        # HP edit distance mode: tied = same minimum edit distance as winner
        winner_ed = (
            summary_df[summary_df['_is_winner']][['read_id', 'hp_edit_distance']]
            .drop_duplicates('read_id')
            .rename(columns={'hp_edit_distance': '_w_ed'})
        )
        merged = summary_df.merge(winner_ed, on='read_id', how='inner')
        merged['_tied'] = merged['hp_edit_distance'] == merged['_w_ed']
        subtitle = "tied on HP-aware edit distance; winner resolved by alignment span"
    else:
        # Legacy fallback: tied on _five_rescued + _conf_rank
        winner_scores = (
            summary_df[summary_df['_is_winner']][['read_id', '_five_rescued', '_conf_rank']]
            .drop_duplicates('read_id')
            .rename(columns={'_five_rescued': '_w_rescued', '_conf_rank': '_w_conf'})
        )
        merged = summary_df.merge(winner_scores, on='read_id', how='inner')
        merged['_tied'] = (
            (merged['_five_rescued'] == merged['_w_rescued']) &
            (merged['_conf_rank']    == merged['_w_conf'])
        )
        subtitle = "tied on five_prime_rescued + confidence; winner resolved by span/junctions"

    aligners = sorted(merged['_aligner'].unique())

    tied_wide = (
        merged.pivot_table(index='read_id', columns='_aligner',
                           values='_tied', aggfunc='max')
        .reindex(columns=aligners)
        .fillna(False)
        .astype(bool)
    )

    counts = upsetplot.from_indicators(tied_wide)

    fig = plt.figure(figsize=(max(8, len(aligners) * 1.5), 5))
    upsetplot.plot(counts, fig=fig, show_counts=True, subset_size='count',
                   totals_plot_elements=2, sort_by='cardinality',
                   sort_categories_by='-cardinality')
    fig.suptitle(
        f"Aligners tied for best score per read\n({subtitle})",
        y=1.02,
    )
    fig.savefig(str(output_path), dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info("UpSet plot written to %s", output_path)

_CONFIDENCE_RANK: Dict[str, int] = {'high': 2, 'medium': 1, 'low': 0}


# ---------------------------------------------------------------------------
# Chimeric-junction overhang filter
# ---------------------------------------------------------------------------

# Parameters for the short-intron relaxation (user-tunable but sane defaults):
#   introns < _SHORT_INTRON_BP: relax to 1 nt overhang if the junction is
#   "well-supported" (≥ _MIN_JUNCTION_SUPPORT aligner×read pairs AND at least
#   one of those has ≥ _GOOD_OVERHANG_BP overhang).
_SHORT_INTRON_BP       = 500
_MIN_JUNCTION_SUPPORT  = 5     # minimum aligner×read pairs reporting the junction
_GOOD_OVERHANG_BP      = 10    # at least one observation must have this much overhang

# Parameters for the long-intron quality check.
#   introns > _LONG_INTRON_BP: apply a higher min-overhang threshold AND require
#   that the overhang alignment quality is reasonable.  The quality metric is
#   hp_edit_distance / (left_overhang + right_overhang): the denominator is the
#   total exon-aligned bases for the junction (excluding the intron gap), which
#   is the natural per-base error rate for the alignment evidence supporting the
#   junction.  False introns invented by gapmm2/deSALT that bridge unrelated
#   sequences tend to have either short overhangs OR elevated per-base error.
_LONG_INTRON_BP        = 5000   # introns longer than this trigger the quality check
_LONG_INTRON_OVERHANG  = 50     # min overhang required for long introns
_MAX_OVERHANG_ERROR    = 0.20   # max allowed hp_ed / (left_ov + right_ov)

# Minimum aligned (M/=/X) bases required — reads below this threshold are
# minimal-anchor alignment artifacts common in poly-A/T regions (e.g. CIGAR
# 40S20M4S has only 20 aligned bases despite MAPQ 60, NM=0).  They are
# flagged even when they have no junctions (the chimera loop otherwise exits
# early with flag=0 for junction-free reads).
_MIN_ALIGNED_BASES     = 50

# Junction-local perfect-match anchor gate (off by default — 0 disables it).
# An N-op alignment is flagged suspicious when the exact-match run of read bases
# immediately abutting the intron (on its weaker flank) is shorter than this,
# UNLESS every junction in the read is rescued by the cross-read support
# relaxation (recurrent real junctions are spared).  Empirically (A549 chr5 DRS,
# 2026-06-13) real introns have a ~13 bp median min-flank perfect-match anchor vs
# ~2 bp for spurious mapPacBio splits; K≈10 rejects ~94% of spurious while the
# support relaxation rescues real-but-noisy-anchor junctions.  Default 0 keeps the
# pipeline (and the valued yeast mapPacBio panel) byte-identical; human runs set 10.
# See memory project-hped-anchor-gate.
_MIN_JUNCTION_ANCHOR   = 0


def _compute_junction_stats(
    all_df: "pd.DataFrame",
    min_junction_anchor_bp: int = 0,
) -> Dict[Tuple[int, int], Tuple[int, int]]:
    """Build a per-junction support dictionary from the concatenated all_df.

    For each (intron_start, intron_end) pair, returns:
        (n_aligner_read_pairs, max_min_overhang_seen)

    where:
        n_aligner_read_pairs = how many (read_id, aligner) rows reported this junction
        max_min_overhang_seen = maximum of min(left_ov, right_ov) across all those rows

    This is used by the short-intron relaxation: if a short intron is reported
    by many reads and at least one has good overhang, it is almost certainly real.

    When min_junction_anchor_bp > 0 and the 'min_junction_anchor' column is present,
    only rows whose per-read perfect-match anchor meets the threshold contribute to
    the support count.  This ensures that the rescue escape hatch for low-anchor
    alignments (e.g. GMAP with median 4 bp anchor) requires at least one K-quality
    sequence observation to vouch for the junction — not just a high observation count
    from consistently low-quality reads.
    """
    junc_stats: Dict[Tuple[int, int], Tuple[int, int]] = {}  # {(js,je): (count, max_min_ov)}

    aln_start_col  = all_df.get('alignment_start', pd.Series(0, index=all_df.index)).fillna(0)
    aln_end_col    = all_df.get('alignment_end',   pd.Series(0, index=all_df.index)).fillna(0)
    junc_col       = all_df.get('junctions',       pd.Series('',  index=all_df.index))
    anchor_gate_on = min_junction_anchor_bp > 0 and 'min_junction_anchor' in all_df.columns
    anchor_col     = (
        all_df['min_junction_anchor'].fillna(_NO_JUNCTION_ANCHOR).astype('int64')
        if anchor_gate_on else None
    )

    for idx in all_df.index:
        if anchor_gate_on and int(anchor_col.iat[idx]) < min_junction_anchor_bp:
            continue  # this row's junctions don't meet K-quality; skip its votes
        junc_str  = junc_col.iat[idx] if idx < len(junc_col) else ''
        aln_start = int(aln_start_col.iat[idx])
        aln_end   = int(aln_end_col.iat[idx])

        for js, je in _parse_junctions_list(junc_str):
            left_ov  = js - aln_start
            right_ov = aln_end - je
            if left_ov < 0 or right_ov < 0:
                continue
            min_ov = min(left_ov, right_ov)
            key = (js, je)
            cnt, max_ov = junc_stats.get(key, (0, 0))
            junc_stats[key] = (cnt + 1, max(max_ov, min_ov))

    return junc_stats


def _add_chimera_flag(
    rep_df: "pd.DataFrame",
    overhang_table: OverhangTable,
    junction_stats: Dict[Tuple[int, int], Tuple[int, int]],
    short_intron_bp: int = _SHORT_INTRON_BP,
    min_junction_support: int = _MIN_JUNCTION_SUPPORT,
    good_overhang_bp: int = _GOOD_OVERHANG_BP,
    long_intron_bp: int = _LONG_INTRON_BP,
    long_intron_overhang: int = _LONG_INTRON_OVERHANG,
    max_overhang_error: float = _MAX_OVERHANG_ERROR,
    min_aligned_bases: int = _MIN_ALIGNED_BASES,
    min_junction_anchor_bp: int = _MIN_JUNCTION_ANCHOR,
) -> "pd.Series":
    """Return a 0/1 Series indicating whether each row in rep_df has a suspicious intron.

    A row is flagged (1) when ANY of the following conditions holds:

    Minimal-anchor filter (applied to ALL reads, with or without junctions):
        ``aligned_bases < min_aligned_bases`` — the read has fewer than
        ``min_aligned_bases`` (default 50) M/=/X CIGAR bases.  These are
        alignment artifacts common in poly-A/T regions (e.g. CIGAR 40S20M4S
        → only 20 aligned bases despite MAPQ 60, NM=0).  The ``aligned_bases``
        column must be present in rep_df (joined from per-aligner corrected
        BAMs via ``_read_hp_edit_distances``).  When the column is absent the
        check is skipped.

    Junction overhang filter (applied to reads that have N-ops):
        ANY junction fails the overhang threshold for its intron size AND does
        not qualify for the short-intron relaxation.

    Short-intron relaxation (intron_size < short_intron_bp):
        A junction qualifies for relaxed treatment (min_overhang ≥ 1 is sufficient)
        when the junction appears in ≥ min_junction_support aligner×read pairs and
        at least one of those pairs has ≥ good_overhang_bp overhang.  This handles:
          - Annotated splice sites (high read support guarantees they are real)
          - Novel alternative junctions with convincing cross-read evidence
        Rationale: a junction seen across many reads with good overhang in at least
        one cannot be an alignment artifact; requiring large overhang in every read
        would wrongly penalise aligners that clip close to a real short intron.

    Long-intron quality check (intron_size > long_intron_bp):
        In addition to the stricter min-overhang threshold (long_intron_overhang),
        the alignment quality of the overhang region must be acceptable.  Quality is
        measured as hp_edit_distance / (left_overhang + right_overhang): the per-base
        HP-aware error rate over the exon-aligned bases flanking the intron.  False
        introns invented by gapmm2/deSALT that bridge unrelated sequences tend to have
        either short overhangs OR elevated per-base error in those overhangs.  A real
        exon-1 (≥ long_intron_overhang bp) aligning cleanly validates the junction;
        a forced alignment of unrelated sequence fails the quality gate.

    Parameters
    ----------
    rep_df:
        Representative row per (read_id, aligner) from merge_corrected_tsvs.
        Must contain hp_edit_distance (float, may be inf) for the quality check.
        May contain aligned_bases (int) for the minimal-anchor check.
    overhang_table:
        Calibrated or default OverhangTable.
    junction_stats:
        From _compute_junction_stats(all_df).
    short_intron_bp, min_junction_support, good_overhang_bp:
        Parameters for the short-intron relaxation condition.
    long_intron_bp, long_intron_overhang, max_overhang_error:
        Parameters for the long-intron quality check.
    min_aligned_bases:
        Minimum number of M/=/X CIGAR bases required.  Reads below this threshold
        are flagged regardless of junction content.

    Returns
    -------
    Integer Series (0 = no chimeric junction found, 1 = suspicious junction found).
    Indexed to match rep_df.
    """
    aln_start_col   = rep_df.get('alignment_start', pd.Series(0, index=rep_df.index)).fillna(0)
    aln_end_col     = rep_df.get('alignment_end',   pd.Series(0, index=rep_df.index)).fillna(0)
    junc_col        = rep_df.get('junctions',       pd.Series('',  index=rep_df.index))
    hp_ed_col       = rep_df.get('hp_edit_distance', pd.Series(float('inf'), index=rep_df.index)).fillna(float('inf'))
    # aligned_bases is optional — only present when per-aligner corrected BAMs were supplied
    has_aligned_col = 'aligned_bases' in rep_df.columns
    aln_bases_col   = rep_df['aligned_bases'].fillna(0).astype(int) if has_aligned_col else None
    # min_junction_anchor is optional — only present when HP-ED scoring ran; the gate
    # is off unless min_junction_anchor_bp > 0 (default keeps behavior byte-identical)
    anchor_gate_on  = min_junction_anchor_bp > 0 and 'min_junction_anchor' in rep_df.columns
    anchor_col      = (
        rep_df['min_junction_anchor'].fillna(_NO_JUNCTION_ANCHOR).astype('int64')
        if anchor_gate_on else None
    )

    flags: List[int] = []
    for i, idx in enumerate(rep_df.index):
        junc_str  = junc_col.iat[i]
        aln_start = int(aln_start_col.iat[i])
        aln_end   = int(aln_end_col.iat[i])
        hp_ed     = float(hp_ed_col.iat[i])
        juncs     = _parse_junctions_list(junc_str)

        # ── Minimal-anchor filter (no-junction reads included) ────────────
        if has_aligned_col:
            ab = int(aln_bases_col.iat[i])
            if ab < min_aligned_bases:
                flags.append(1)
                continue

        if not juncs:
            flags.append(0)
            continue

        row_is_chimeric = False
        for js, je in juncs:
            intron_size = je - js
            if intron_size <= 0:
                continue

            left_ov  = js - aln_start
            right_ov = aln_end - je
            if left_ov < 0 or right_ov < 0:
                # Alignment coordinates don't enclose the junction — suspicious
                row_is_chimeric = True
                break
            min_ov   = min(left_ov, right_ov)

            # Long-intron path: stricter overhang + quality check
            if intron_size > long_intron_bp:
                required = max(overhang_table.lookup(intron_size), long_intron_overhang)
                if min_ov < required:
                    row_is_chimeric = True
                    break
                # Quality check: hp_ed per exon-aligned base must be below threshold
                exon_bases = left_ov + right_ov
                if exon_bases > 0 and hp_ed < float('inf'):
                    per_base_err = hp_ed / exon_bases
                    if per_base_err >= max_overhang_error:
                        row_is_chimeric = True
                        break
                elif hp_ed == float('inf'):
                    # No HP edit distance available — treat as suspicious for long introns
                    row_is_chimeric = True
                    break
                continue  # passes long-intron checks

            required = overhang_table.lookup(intron_size)

            if min_ov >= required:
                continue  # passes strict threshold

            # Below strict threshold — check short-intron relaxation
            if intron_size <= short_intron_bp and min_ov >= 1:
                key = (js, je)
                cnt, max_ov = junction_stats.get(key, (0, 0))
                if cnt >= min_junction_support and max_ov >= good_overhang_bp:
                    # Well-supported short intron: allow even 1 nt overhang
                    continue

            # Failed both strict and relaxed checks
            row_is_chimeric = True
            break

        # ── Junction-local perfect-match anchor gate ─────────────────────
        # An N-op may only win the consensus by splitting a noisy read if the
        # exact-match run abutting the intron is long enough.  Spurious
        # mapPacBio splits leave a ~2 bp anchor; real introns ~13 bp.  Real
        # but noisy-anchor junctions are recurrent → spared when EVERY junction
        # in the read is rescued by the cross-read support relaxation.  This is
        # the per-aligner discriminator read-support alone cannot provide
        # (the harmful wins are sole-wins, support=1).
        if not row_is_chimeric and anchor_gate_on:
            anchor = int(anchor_col.iat[i])
            if anchor < min_junction_anchor_bp:
                all_supported = True
                for js, je in juncs:
                    cnt, max_ov = junction_stats.get((js, je), (0, 0))
                    if not (cnt >= min_junction_support and max_ov >= good_overhang_bp):
                        all_supported = False
                        break
                if not all_supported:
                    row_is_chimeric = True

        flags.append(1 if row_is_chimeric else 0)

    return pd.Series(flags, index=rep_df.index, dtype=int, name='_chimera_ok')


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_junctions(junc_str) -> frozenset:
    """Parse semicolon-separated 'start-end' junction string into a frozenset."""
    if pd.isna(junc_str) or not junc_str or str(junc_str) == 'nan':
        return frozenset()
    result = set()
    for part in str(junc_str).split(';'):
        part = part.strip()
        if not part:
            continue
        fields = part.split('-')
        if len(fields) == 2:
            try:
                result.add((int(fields[0]), int(fields[1])))
            except ValueError:
                pass
    return frozenset(result)


def _normalize_read_id(read_id_series: "pd.Series") -> "pd.Series":
    """Vectorized application of ``_normalize_bam_read_name`` to a TSV read_id column.

    Covers every suffix shape ``_normalize_bam_read_name`` handles: mapPacBio
    ``_pt:i:N``, BBmap ``_<N>_length=<N>``, generic SAM aux leaks
    (``_XA:Z:``, ``_XC:i:``, etc.), and Dorado metadata keys
    (``_runid=``, ``_ch=``, ``_start_time=``, ...).

    Equivalent to ``read_id_series.apply(_normalize_bam_read_name)`` but
    vectorized for the TSV-merge hot path on large per-aligner outputs.
    """
    # Strip post-whitespace (mirrors the first half of _normalize_bam_read_name).
    out = read_id_series.str.split(' ', n=1).str[0].str.split('\t', n=1).str[0]
    # Strip known underscore-encoded suffixes using the canonical regex
    # imported from consensus.py so the two normalizers cannot drift.
    out = out.str.replace(_UNDERSCORE_COMMENT_RE.pattern, '', regex=True)
    return out


_CORRECTED_TSV_MANIFEST_HEADER = [
    'region_id', 'chrom', 'start', 'end', 'tsv_path', 'n_rows', 'sha256',
]


def _read_single_corrected_tsv(tsv_path: Path) -> pd.DataFrame:
    import warnings as _warnings
    with _warnings.catch_warnings():
        _warnings.simplefilter('ignore')
        return pd.read_csv(tsv_path, sep='\t', index_col=False)


def _read_corrected_tsv_or_manifest(tsv_path: Path) -> pd.DataFrame:
    with tsv_path.open() as fh:
        header = fh.readline().rstrip('\n').split('\t')

    if header != _CORRECTED_TSV_MANIFEST_HEADER:
        return _read_single_corrected_tsv(tsv_path)

    from ..bam.tsv_partition import load_manifest

    dfs = []
    for entry in load_manifest(tsv_path):
        region_tsv = entry['tsv_path']
        region_df = _read_single_corrected_tsv(region_tsv)
        if not region_df.empty:
            dfs.append(region_df)
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def _load_tsv(aligner_name: str, tsv_path: Path) -> Optional[pd.DataFrame]:
    """Load one per-aligner TSV, returning None on failure.

    Some per-aligner TSVs (minimap2, gapmm2, uLTRA) are written with
    ``index=True`` (32 header columns, 33 data columns).  Without
    ``index_col=False``, pandas auto-uses the first data column (UUID) as
    the row index and shifts all column names one position right, mapping
    ``read_id`` → chromosome, ``chrom`` → strand, etc.  Using
    ``index_col=False`` prevents this by treating every column as data and
    dropping the trailing extra field instead.
    """
    if not tsv_path.exists():
        logger.warning("Per-aligner TSV not found, skipping: %s", tsv_path)
        return None
    try:
        df = _read_corrected_tsv_or_manifest(tsv_path)
        if df.empty:
            logger.warning("Empty per-aligner TSV, skipping: %s", tsv_path)
            return None
        # Drop any pre-existing winning_aligner column.  Per-aligner TSVs may
        # carry a stale winning_aligner from a prior merge run (e.g. the
        # validation fixture was generated from a 5-aligner merge and has
        # legacy aligner names in that column).  Keeping it causes a pandas
        # rename collision at line 1354: after the merge step all_df already
        # contains "winning_aligner" (stale) and the fresh "_winning_aligner"
        # from winner_df; renaming _winning_aligner → winning_aligner produces
        # two columns with the same name and to_csv writes both — the stale
        # one appears first, so callers that read the output see wrong winners.
        # merge_corrected_tsvs is the canonical source of winning_aligner;
        # whatever the per-aligner input says is always superseded.
        if 'winning_aligner' in df.columns:
            df = df.drop(columns=['winning_aligner'])
        df['_aligner'] = aligner_name
        if 'read_id' in df.columns:
            df['read_id'] = _normalize_read_id(df['read_id'])
        return df
    except Exception as exc:
        logger.warning("Failed to load %s: %s", tsv_path, exc)
        return None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def _stage_raw_bams(
    per_aligner_raw_bams: Dict[str, str],
    scratch_root: Optional[str] = None,
) -> Iterator[Dict[str, str]]:
    """Copy raw BAMs to local scratch before merge + write to avoid Lustre cold-read overhead.

    ProcessPoolExecutor workers maintain per-process Lustre client caches. BAMs read
    inside workers don't warm the main process's cache, so both HP scoring and the
    subsequent BAM write hit cold Lustre. Staging to local NVMe ($L_SCRATCH / $SCRATCH)
    once in the main process eliminates this for all readers.

    Yields a {aligner: staged_path_str} dict. Cleans up the staging directory on exit.
    """
    if not per_aligner_raw_bams:
        yield per_aligner_raw_bams
        return
    if scratch_root is None:
        scratch_root = os.environ.get('L_SCRATCH', os.environ.get('SCRATCH')) or None
    stage_dir = Path(tempfile.mkdtemp(prefix='rectify_bam_stage_', dir=scratch_root))
    try:
        staged: Dict[str, str] = {}
        for aligner, bam_path in per_aligner_raw_bams.items():
            dest = stage_dir / Path(str(bam_path)).name
            shutil.copy2(str(bam_path), str(dest))
            staged[aligner] = str(dest)
            logger.info(
                'Staged %s BAM to local scratch (%.0f MB): %s',
                aligner, Path(str(bam_path)).stat().st_size / 1e6, dest,
            )
        yield staged
    finally:
        shutil.rmtree(stage_dir, ignore_errors=True)


def merge_corrected_tsvs(
    per_aligner_tsvs: Dict[str, Path],
    output_tsv: Path,
    summary_tsv: Optional[Path] = None,
    per_aligner_corrected_bams: Optional[Dict[str, str]] = None,
    per_aligner_raw_bams: Optional[Dict[str, str]] = None,
    genome: Optional[Dict[str, str]] = None,
    penalty_table=None,
    overhang_table: Optional[OverhangTable] = None,
    strict_lazy_identity: bool = True,
    lazy_scoring_workers: int = 1,
    min_junction_anchor_bp: int = _MIN_JUNCTION_ANCHOR,
) -> Path:
    """
    Merge N per-aligner corrected TSVs into a single corrected_reads.tsv.

    Winner selection (when *per_aligner_corrected_bams* is provided):

    0. ``_chimera_ok`` (0 = plausible, 1 = suspicious) — computed from junction
       overhang vs. intron size using *overhang_table*.  Aligners whose introns
       lack sufficient flanking evidence are sorted after plausible ones.  When
       all aligners are suspicious, the HP edit distance resolves the tie.
       This pre-filter prevents chimeric deSALT/gapmm2 long-range false junctions
       from winning just because their enormous spans yield large HP edit distances
       relative to the true alignment.  Pass ``--junction-overhang-table`` to enable.
    1. ``hp_edit_distance`` — HP-aware edit distance from the corrected CIGAR.
       Soft-clips count as 1.0/base; intron N-ops are free; deletions and
       insertions use *penalty_table* when provided.  Lower = better.
    2. alignment span (tiebreak) — wider alignment preferred.

    Winner selection (fallback when no BAMs supplied):

    1. ``five_prime_rescued``, 2. ``confidence``, 3. 3-prime position agreement,
    4. alignment span, 5. junction count  (legacy five-level sort).

    Multi-peak NET-seq reads (multiple rows per read_id with different
    ``fraction`` values) are handled correctly: the winning aligner is
    determined from a single representative row per read, then *all* rows
    from that aligner for that read_id are written to the output.

    Parameters
    ----------
    per_aligner_tsvs:
        Mapping from aligner name to corrected_reads.tsv path.
    output_tsv:
        Destination path for the merged corrected_reads.tsv.
    summary_tsv:
        Optional path for a per-read aligner comparison table.
    per_aligner_corrected_bams:
        Optional mapping from aligner name to corrected.bam path.  When
        provided, HP-aware edit distances are computed from the final corrected
        CIGARs and used as the primary winner criterion.
    per_aligner_raw_bams:
        Optional mapping from aligner name to raw/pre-correction BAM path.  When
        provided without ``per_aligner_corrected_bams``, corrections are applied
        in memory from each aligner's TSV and HP-aware edit distances are scored
        lazily without materializing per-aligner corrected BAMs.
    genome:
        Optional genome dict {chrom: seq} for HP run-length computation during
        edit distance scoring.  If absent, flat penalties are used.
    penalty_table:
        Optional ``HpPenaltyTable`` for HP-context-aware del/ins costs.  If
        absent, fixed fallback costs are used (del=1.0, ins=1.25).
    overhang_table:
        Optional ``OverhangTable`` from ``calibrate_junction_overhang.py``.
        When provided, aligners with junctions whose overhang is insufficient
        for their intron size are penalised (sorted last) before HP edit
        distance is applied.  Short introns (< 500 bp) are treated with
        relaxed thresholds when backed by strong cross-read evidence.
        Load with ``OverhangTable.from_tsv(path)`` or use
        ``OverhangTable.default()`` for the built-in conservative defaults.
    strict_lazy_identity:
        If true, lazy HP scoring fails fast when a TSV read cannot be matched
        back to the raw BAM by raw or normalized QNAME.
    lazy_scoring_workers:
        Number of aligner-level workers for HP edit-distance scoring. Each
        worker scans one BAM and applies one aligner's corrections in memory.

    Returns
    -------
    Path to the written output_tsv.
    """
    dfs: Dict[str, pd.DataFrame] = {}
    for aligner_name, tsv_path in per_aligner_tsvs.items():
        df = _load_tsv(aligner_name, tsv_path)
        if df is not None:
            dfs[aligner_name] = df

    if not dfs:
        raise ValueError("No valid per-aligner corrected TSVs found")

    # Trivial case — only one aligner succeeded
    if len(dfs) == 1:
        aligner_name, df = next(iter(dfs.items()))
        out_df = df.drop(columns=['_aligner'])
        # winning_aligner is required by restore_polya_softclips and
        # build_corrected_consensus_bam — add it even in the single-aligner case.
        out_df['winning_aligner'] = aligner_name
        out_df.to_csv(output_tsv, sep='\t', index=False)
        logger.info("Single aligner (%s) — wrote %d rows to %s",
                    aligner_name, len(out_df), output_tsv)
        return output_tsv

    all_df = pd.concat(dfs.values(), ignore_index=True)

    # ── Span (always needed as tiebreaker) ────────────────────────────────
    aln_start = all_df.get('alignment_start', pd.Series(0, index=all_df.index)).fillna(0)
    aln_end   = all_df.get('alignment_end',   pd.Series(0, index=all_df.index)).fillna(0)
    all_df['_span'] = (aln_end - aln_start).astype(int)

    # five_prime_rescued is needed by both the HP sort and the lazy scoring
    # prefilter below.
    rescued_col = all_df.get('five_prime_rescued', pd.Series(0, index=all_df.index))
    all_df['_five_rescued'] = (
        pd.to_numeric(rescued_col, errors='coerce').fillna(0).astype(int)
    )

    # Before paying for HP-aware CIGAR scoring, identify aligner/read pairs
    # that can still win after the chimera-overhang tier.  Candidates with a
    # worse effective chimera flag cannot win because HP edit distance is only
    # sorted after that flag.  We still score the sole surviving candidate for
    # each read so the minimal-anchor hard filter keeps using corrected CIGAR
    # aligned-bases for the eventual winner.
    lazy_read_filters: Optional[Dict[str, Set[str]]] = None
    if per_aligner_raw_bams and not per_aligner_corrected_bams:
        rep_preview = (
            all_df.sort_values(
                ['read_id', '_aligner', 'fraction'],
                ascending=[True, True, False],
                na_position='last',
            )
            .groupby(['read_id', '_aligner'], sort=False)
            .first()
            .reset_index()
        )
        _preview_ov_table = overhang_table if overhang_table is not None else OverhangTable.default()
        preview_junc_stats = _compute_junction_stats(all_df, min_junction_anchor_bp=min_junction_anchor_bp)
        preview_chimera_flags = _add_chimera_flag(rep_preview, _preview_ov_table, preview_junc_stats)
        rep_preview['_chimera_ok'] = preview_chimera_flags
        rep_preview['_effective_chimera_ok'] = (
            rep_preview['_chimera_ok'] & (rep_preview['_five_rescued'] == 0)
        ).astype(int)
        keep_pairs = []
        for _rid, _grp in rep_preview.groupby('read_id', sort=False):
            _best_flag = _grp['_effective_chimera_ok'].min()
            keep_pairs.append(_grp.loc[_grp['_effective_chimera_ok'] == _best_flag, ['read_id', '_aligner']])
        if keep_pairs:
            lazy_read_filters = {}
            keep_df = pd.concat(keep_pairs, ignore_index=True)
            for _aligner, _grp in keep_df.groupby('_aligner', sort=False):
                lazy_read_filters[str(_aligner)] = {
                    _normalize_bam_read_name(str(_rid)) for _rid in _grp['read_id']
                }
            kept_pairs = sum(len(v) for v in lazy_read_filters.values())
            logger.info(
                "Lazy HP scoring prefilter retained %d / %d aligner-read pairs",
                kept_pairs,
                len(rep_preview),
            )

    # ── HP-aware edit distance (primary criterion when BAMs provided) ──────
    use_hp_ed = bool(per_aligner_corrected_bams or per_aligner_raw_bams)
    if use_hp_ed:
        if per_aligner_raw_bams and not per_aligner_corrected_bams and strict_lazy_identity:
            missing_raw = sorted(set(dfs) - set(per_aligner_raw_bams))
            if missing_raw:
                raise FileNotFoundError(
                    "Missing raw BAMs for lazy HP scoring: " + ', '.join(missing_raw)
                )
        if per_aligner_corrected_bams:
            logger.info("Computing HP-aware edit distances from corrected BAMs...")
        else:
            logger.info("Computing HP-aware edit distances lazily from raw BAMs + TSVs...")
        ed_rows = []
        score_aligners = per_aligner_corrected_bams or per_aligner_raw_bams or {}
        n_workers = max(1, min(int(lazy_scoring_workers or 1), len(score_aligners)))
        if n_workers > 1:
            if per_aligner_corrected_bams:
                worker_args = [
                    (aligner_name, str(bam_path), genome, penalty_table)
                    for aligner_name, bam_path in score_aligners.items()
                ]
                worker_fn = _score_corrected_bam_for_aligner
            else:
                worker_args = [
                    (
                        aligner_name,
                        str(bam_path),
                        str(per_aligner_tsvs[aligner_name])
                        if aligner_name in per_aligner_tsvs
                        else None,
                        genome,
                        penalty_table,
                        strict_lazy_identity,
                        lazy_read_filters.get(aligner_name, set()) if lazy_read_filters else None,
                    )
                    for aligner_name, bam_path in score_aligners.items()
                ]
                worker_fn = _score_hp_distances_for_aligner
            logger.info("HP edit-distance scoring using %d aligner workers", n_workers)
            pool = ProcessPoolExecutor(max_workers=n_workers)
            try:
                futures = [pool.submit(worker_fn, args) for args in worker_args]
                for future in as_completed(futures):
                    aligner_name, ed_dict = future.result()
                    logger.info(
                        "HP edit distances computed for %s: %d reads",
                        aligner_name,
                        len(ed_dict),
                    )
                    for rid, (ed, ab, anchor) in ed_dict.items():
                        ed_rows.append({
                            'read_id': rid,
                            '_aligner': aligner_name,
                            'hp_edit_distance': ed,
                            'aligned_bases': ab,
                            'min_junction_anchor': anchor,
                        })
            finally:
                # All results consumed; release workers in background so the
                # 22-s queue-drain overhead overlaps with post-scoring merge work.
                pool.shutdown(wait=False)
        else:
            for aligner_name, bam_path in score_aligners.items():
                logger.info("HP edit-distance scoring for %s", aligner_name)
                if per_aligner_corrected_bams:
                    ed_dict = _read_hp_edit_distances(str(bam_path), genome, penalty_table)
                else:
                    tsv_path = per_aligner_tsvs.get(aligner_name)
                    if tsv_path is None:
                        if strict_lazy_identity:
                            raise FileNotFoundError(
                                f"Missing corrected TSV for lazy HP scoring: {aligner_name}"
                            )
                        logger.warning(
                            "Missing corrected TSV for lazy HP scoring, skipping %s",
                            aligner_name,
                        )
                        continue
                    ed_dict = _read_hp_edit_distances_from_raw_bam(
                        str(bam_path),
                        Path(tsv_path),
                        genome,
                        penalty_table,
                        strict=strict_lazy_identity,
                        read_id_filter=(
                            lazy_read_filters.get(aligner_name, set())
                            if lazy_read_filters
                            else None
                        ),
                    )
                logger.info(
                    "HP edit distances computed for %s: %d reads",
                    aligner_name,
                    len(ed_dict),
                )
                for rid, (ed, ab, anchor) in ed_dict.items():
                    ed_rows.append({
                        'read_id': rid,
                        '_aligner': aligner_name,
                        'hp_edit_distance': ed,
                        'aligned_bases': ab,
                        'min_junction_anchor': anchor,
                    })
        if ed_rows:
            ed_df = pd.DataFrame(ed_rows)
            all_df = all_df.merge(ed_df, on=['read_id', '_aligner'], how='left')
            all_df['hp_edit_distance'] = all_df['hp_edit_distance'].fillna(float('inf'))
            all_df['aligned_bases']    = all_df['aligned_bases'].fillna(0).astype(int)
            # reads not scored (no HP-ED) get the no-junction sentinel so the
            # anchor gate never fires on them
            all_df['min_junction_anchor'] = (
                all_df['min_junction_anchor'].fillna(_NO_JUNCTION_ANCHOR).astype('int64')
            )
            logger.info("HP edit distances joined for %d aligner×read pairs", len(ed_rows))
        else:
            logger.warning("No HP edit distances computed — falling back to legacy sort")
            use_hp_ed = False

    # ── five_prime_rescued flag (always computed — needed for sort in both paths) ──
    rescued_col = all_df.get('five_prime_rescued', pd.Series(0, index=all_df.index))
    all_df['_five_rescued'] = (
        pd.to_numeric(rescued_col, errors='coerce').fillna(0).astype(int)
    )

    # ── Legacy scoring columns (fallback) ─────────────────────────────────
    if not use_hp_ed:
        all_df['_conf_rank'] = (
            all_df['confidence'].map(_CONFIDENCE_RANK).fillna(0).astype(int)
        )
        all_df['_n_junc'] = all_df.get('n_junctions', 0).fillna(0).astype(int)
        # Group on (read_id, chrom, corrected_3prime) — NOT just (read_id,
        # corrected_3prime). Paralogous loci can place a read at the same
        # numeric position on different chromosomes (e.g. cat1_plus_1: deSALT
        # → chrVI:8703, minimap2/gapmm2/uLTRA → chrXIV:10610). Grouping on
        # position alone made a 1-aligner cross-chrom outlier collapse into
        # the 3-aligner same-chrom consensus and tie in `_n_agree`, letting
        # the tail of the sort key (span / n_junc) pick the paralog.
        pos_counts = (
            all_df.groupby(['read_id', 'chrom', 'corrected_3prime'])
            .size()
            .reset_index(name='_n_agree')
        )
        all_df = all_df.merge(
            pos_counts, on=['read_id', 'chrom', 'corrected_3prime'], how='left'
        )

    # ── Representative row per (read_id, aligner) for ranking ────────────
    rep_df = (
        all_df.sort_values(
            ['read_id', '_aligner', 'fraction'],
            ascending=[True, True, False],
            na_position='last',
        )
        .groupby(['read_id', '_aligner'], sort=False)
        .first()
        .reset_index()
    )

    # ── Chimeric-junction overhang filter (always active) ────────────────
    # Use the calibrated table when supplied; fall back to OverhangTable.default()
    # so the filter is never silently skipped.  The default 22 bp threshold is
    # sufficient to catch e.g. a 3681N intron with only 6M overhang on one side.
    _ov_table = overhang_table if overhang_table is not None else OverhangTable.default()
    if overhang_table is None:
        logger.info("No junction-overhang table supplied — using OverhangTable.default()")
    logger.info("Applying chimeric-junction overhang filter ...")
    junc_stats = _compute_junction_stats(all_df, min_junction_anchor_bp=min_junction_anchor_bp)
    chimera_flags = _add_chimera_flag(
        rep_df, _ov_table, junc_stats, min_junction_anchor_bp=min_junction_anchor_bp,
    )
    rep_df['_chimera_ok'] = chimera_flags
    n_flagged = int(chimera_flags.sum())
    n_total   = len(chimera_flags)
    logger.info(
        "Chimeric flag: %d / %d aligner×read pairs flagged (%.1f%%)",
        n_flagged, n_total, 100 * n_flagged / max(1, n_total),
    )

    # Sort so the best candidate per read comes first.
    # _effective_chimera_ok applies the chimera penalty only for non-rescued reads:
    #   - _five_rescued=0 (normal): _effective_chimera_ok = _chimera_ok (0=good, 1=suspect)
    #   - _five_rescued=1 (Cat3):   _effective_chimera_ok = 0 always (exempt from penalty)
    # Cat3 reads have intentionally short 5' overhangs (5–20 bp) that trigger _chimera_ok=1,
    # but they are valid corrections.  Sorting _five_rescued unconditionally BEFORE
    # _chimera_ok made rescued reads unconditional winners even when their hp_edit_distance
    # was worse than a non-rescued aligner's.  The exemption approach lets quality
    # (hp_edit_distance) decide among rescued reads while still protecting them from the
    # chimera penalty.
    if use_hp_ed:
        rep_df['_effective_chimera_ok'] = (
            rep_df['_chimera_ok'] & (rep_df['_five_rescued'] == 0)
        ).astype(int)
        rep_df = rep_df.sort_values(
            ['read_id', '_effective_chimera_ok', 'hp_edit_distance', '_span'],
            ascending=[True, True, True, False],
        )
    else:
        rep_df = rep_df.sort_values(
            ['read_id', '_five_rescued', '_chimera_ok', '_conf_rank', '_n_agree', '_span', '_n_junc'],
            ascending=[True, False, True, False, False, False, False],
        )

    # Pick winning aligner per read
    winner_cols = ['read_id', '_aligner']
    winner_df = (
        rep_df.groupby('read_id', sort=False)
        .first()
        .reset_index()[winner_cols]
        .rename(columns={'_aligner': '_winning_aligner'})
    )

    # ── Hard-filter: drop reads whose winning alignment is a minimal-anchor artifact ──
    # Soft _chimera_ok=1 only re-orders within-read aligners; if ALL aligners for a read
    # produce a short-anchor alignment the "winner" still passes.  Hard-filter here so
    # these reads are absent from the output TSV and consensus BAM entirely.
    if 'aligned_bases' in rep_df.columns:
        winner_ab = (
            rep_df[['read_id', '_aligner', 'aligned_bases']]
            .merge(winner_df, on='read_id', how='inner')
        )
        winner_ab = winner_ab[winner_ab['_aligner'] == winner_ab['_winning_aligner']]
        minimal_ids = set(
            winner_ab.loc[
                winner_ab['aligned_bases'] < _MIN_ALIGNED_BASES, 'read_id'
            ]
        )
        if minimal_ids:
            logger.info(
                "Hard-filtering %d reads with winning aligned_bases < %d "
                "(minimal-anchor artifacts)",
                len(minimal_ids), _MIN_ALIGNED_BASES,
            )
            winner_df = winner_df[~winner_df['read_id'].isin(minimal_ids)]

    # ── Comparison summary (always built; written to file only if requested) ─
    summary = rep_df.merge(winner_df, on='read_id', how='left')
    summary['_is_winner'] = summary['_aligner'] == summary['_winning_aligner']

    # ── Per-aligner effective-utility cluster annotation ──
    #
    # Within each read, group aligners by their (chrom, strand,
    # corrected_3prime, corrected_junctions) tuple — aligners landing at the same biological
    # answer are functionally interchangeable for RECTIFY's outputs, even
    # if their CIGARs differ in mismatch/indel placement.  This makes the
    # HP-ED tiebreak within a cluster informative without losing the signal
    # that multiple aligners are "right enough."
    #
    # The new ``effective_group`` column assigns each row a letter (A, B,
    # C, …) keyed by its cluster within the read.  Cluster letters are
    # ordered by descending size (most-popular cluster = A).  Ties are
    # broken alphabetically by aligner name.
    #
    # ``effectively_matched_winner`` is True when the row's cluster is the
    # same as the winning aligner's cluster.
    def _eff_key(_row):
        _juncs = _row.get('junctions', '')
        # pandas .get() returns NaN (a truthy float) when the column is
        # missing, so `or ''` does NOT catch it — coerce to str first.
        if not isinstance(_juncs, str):
            _juncs = ''
        # Normalize: sort the donor-acceptor tuples so semicolon-order
        # differences across aligners don't shadow real equivalence.
        _parts = tuple(sorted(p for p in _juncs.split(';') if p))
        try:
            _c3p = int(_row['corrected_3prime'])
        except (TypeError, ValueError, KeyError):
            _c3p = -1
        return (_row.get('chrom', ''), _row.get('strand', ''), _c3p, _parts)

    summary['_eff_key'] = summary.apply(_eff_key, axis=1)
    summary['effective_group'] = ''
    summary['effectively_matched_winner'] = False
    for _rid, _grp in summary.groupby('read_id', sort=False):
        # Cluster sizes within this read (descending), then alphabetical
        # aligner name as tiebreak for deterministic labeling.
        _key_to_aligners: Dict[tuple, list] = {}
        for _idx, _row in _grp.iterrows():
            _key_to_aligners.setdefault(_row['_eff_key'], []).append(_row['_aligner'])
        _sorted_keys = sorted(
            _key_to_aligners.items(),
            key=lambda kv: (-len(kv[1]), sorted(kv[1])[0]),
        )
        _key_to_letter = {_k: chr(ord('A') + _i) for _i, (_k, _) in enumerate(_sorted_keys)}
        # Winner's effective key — match column at winner row
        _winner_row = _grp[_grp['_is_winner']]
        _winner_key = (
            _winner_row.iloc[0]['_eff_key']
            if not _winner_row.empty
            else None
        )
        for _idx, _row in _grp.iterrows():
            summary.at[_idx, 'effective_group'] = _key_to_letter[_row['_eff_key']]
            summary.at[_idx, 'effectively_matched_winner'] = (
                _winner_key is not None and _row['_eff_key'] == _winner_key
            )

    if summary_tsv:
        keep = [
            'read_id', '_aligner', '_is_winner', '_chimera_ok', 'hp_edit_distance',
            'aligned_bases', 'corrected_3prime', 'five_prime_position',
            'junctions', '_conf_rank', '_n_agree', '_n_junc', '_five_rescued',
            'effective_group', 'effectively_matched_winner',
            'correction_applied', 'confidence',
        ]
        keep = [c for c in keep if c in summary.columns]
        summary[keep].sort_values(['read_id', '_aligner']).to_csv(
            summary_tsv, sep='\t', index=False
        )
        logger.info("Comparison summary → %s", summary_tsv)

        # Sample-wide effective-utility rollup: per aligner, how often did
        # it carry the same biological answer as the consensus winner?
        # An aligner can lose HP-ED within the winning cluster (cluster A)
        # and still be in the right cluster — this captures that signal.
        try:
            _n_reads = summary['read_id'].nunique()
            _rollup = (
                summary.groupby('_aligner')['effectively_matched_winner']
                .sum()
                .sort_values(ascending=False)
            )
            logger.info("Sample-wide effective utility (out of %d reads):", _n_reads)
            for _aln, _matched in _rollup.items():
                _pct = 100.0 * _matched / max(1, _n_reads)
                logger.info("  %-10s matched winner cluster in %5d / %5d reads (%5.1f%%)",
                            _aln, int(_matched), _n_reads, _pct)
        except Exception as _roll_exc:
            logger.warning("effective-utility rollup failed: %s", _roll_exc)

    # ── Select all rows from winning aligner (handles multi-peak reads) ──
    result_df = all_df.merge(winner_df, on='read_id', how='left')
    result_df = result_df[result_df['_aligner'] == result_df['_winning_aligner']].copy()

    # Expose winning_aligner as a public column (needed by restore-softclip
    # to pull reads from the correct raw aligner BAM).
    if '_winning_aligner' in result_df.columns:
        result_df.rename(columns={'_winning_aligner': 'winning_aligner'}, inplace=True)

    drop_cols = [
        '_aligner', '_conf_rank', '_five_rescued',
        '_span', '_n_junc', '_n_agree', 'hp_edit_distance', '_chimera_ok',
        '_effective_chimera_ok', 'aligned_bases',
    ]
    result_df.drop(
        columns=[c for c in drop_cols if c in result_df.columns],
        inplace=True,
    )
    result_df.reset_index(drop=True, inplace=True)

    result_df.to_csv(output_tsv, sep='\t', index=False)
    logger.info(
        "Merged %d aligners → %d rows → %s",
        len(dfs), len(result_df), output_tsv,
    )

    # UpSet plot: aligners tied for best score per read
    upset_path = Path(output_tsv).with_name(
        Path(output_tsv).stem + '_aligner_upset.png'
    )
    _plot_aligner_agreement_upset(summary, upset_path)

    return output_tsv


def identify_cat5_candidates(
    per_aligner_tsvs: Dict[str, Path],
    output_tsv: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Identify reads where ≥2 aligners each contribute a unique post-correction junction.

    A read is a Cat5 candidate when aligner A has junction J_A and aligner B has
    junction J_B, where J_A ∉ B's junctions AND J_B ∉ A's junctions — meaning
    neither aligner alone produces the complete splice pattern.

    Junctions here are taken from the ``junctions`` column of each
    corrected_reads.tsv (post-correction, so Cat3-rescued junctions are included
    and false N-ops absorbed by Cat4 are excluded).

    Parameters
    ----------
    per_aligner_tsvs:
        Mapping from aligner name to corrected_reads.tsv path.
    output_tsv:
        Optional path to write the candidate table.

    Returns
    -------
    DataFrame with columns: read_id, aligner_a, aligner_b,
    unique_junctions_a, unique_junctions_b, chrom, strand.
    An empty DataFrame is returned when no candidates exist.
    """
    per_aligner_data: Dict[str, pd.DataFrame] = {}
    for aligner_name, tsv_path in per_aligner_tsvs.items():
        df = _load_tsv(aligner_name, tsv_path)
        if df is None:
            continue
        # One row per read_id for junction comparison
        rep = df.groupby('read_id').first().reset_index()
        rep['_junc_set'] = rep['junctions'].apply(_parse_junctions)
        per_aligner_data[aligner_name] = rep[
            ['read_id', 'chrom', 'strand', '_junc_set']
        ].copy()

    if len(per_aligner_data) < 2:
        empty = pd.DataFrame(columns=[
            'read_id', 'aligner_a', 'aligner_b',
            'unique_junctions_a', 'unique_junctions_b', 'chrom', 'strand',
        ])
        if output_tsv:
            empty.to_csv(output_tsv, sep='\t', index=False)
        return empty

    aligner_names = list(per_aligner_data.keys())

    # Build per-aligner lookup dicts for fast access
    junc_lookup: Dict[str, Dict[str, frozenset]] = {
        a: dict(zip(df['read_id'], df['_junc_set']))
        for a, df in per_aligner_data.items()
    }
    chrom_lookup: Dict[str, Dict[str, str]] = {
        a: dict(zip(df['read_id'], df['chrom']))
        for a, df in per_aligner_data.items()
    }
    strand_lookup: Dict[str, Dict[str, str]] = {
        a: dict(zip(df['read_id'], df['strand']))
        for a, df in per_aligner_data.items()
    }

    # Reads present in ≥2 aligners
    presence: Dict[str, int] = {}
    for reads in junc_lookup.values():
        for rid in reads:
            presence[rid] = presence.get(rid, 0) + 1
    multi_reads = {rid for rid, cnt in presence.items() if cnt >= 2}

    candidates = []
    for i, aligner_a in enumerate(aligner_names):
        for aligner_b in aligner_names[i + 1:]:
            shared = (
                set(junc_lookup[aligner_a])
                & set(junc_lookup[aligner_b])
                & multi_reads
            )
            for read_id in shared:
                juncs_a = junc_lookup[aligner_a].get(read_id, frozenset())
                juncs_b = junc_lookup[aligner_b].get(read_id, frozenset())
                if not juncs_a and not juncs_b:
                    continue
                unique_a = juncs_a - juncs_b
                unique_b = juncs_b - juncs_a
                if unique_a and unique_b:
                    candidates.append({
                        'read_id': read_id,
                        'aligner_a': aligner_a,
                        'aligner_b': aligner_b,
                        'unique_junctions_a': ';'.join(
                            f'{s}-{e}' for s, e in sorted(unique_a)
                        ),
                        'unique_junctions_b': ';'.join(
                            f'{s}-{e}' for s, e in sorted(unique_b)
                        ),
                        'chrom': chrom_lookup[aligner_a].get(read_id, ''),
                        'strand': strand_lookup[aligner_a].get(read_id, ''),
                    })

    result = pd.DataFrame(candidates)
    if not result.empty:
        result = result.drop_duplicates(subset=['read_id', 'aligner_a', 'aligner_b'])
        result = result.sort_values(['chrom', 'read_id', 'aligner_a'])

    if output_tsv:
        result.to_csv(output_tsv, sep='\t', index=False)
        logger.info("Cat5 candidates: %d read-pair combinations → %s",
                    len(result), output_tsv)

    return result
