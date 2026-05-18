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

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

try:
    import pysam as _pysam
except ImportError:  # pragma: no cover
    _pysam = None  # type: ignore

from rectify.core.splice.calibrate_junction_overhang import OverhangTable, _parse_junctions as _parse_junctions_list

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# HP-aware edit distance from corrected CIGAR
# ---------------------------------------------------------------------------

def _bare_uuid(name: str) -> str:
    """Strip mapPacBio _pt:i:N or ' pt:i:N' suffix to recover the canonical read ID."""
    if "_pt:i:" in name:
        return name.split("_pt:i:")[0]
    if " pt:i:" in name:
        return name.split(" pt:i:")[0]
    return name


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


def _read_hp_edit_distances(
    bam_path: str,
    genome: Optional[Dict[str, str]] = None,
    penalty_table=None,
) -> Dict[str, Tuple[float, int]]:
    """Read a corrected BAM and return {bare_read_id: (hp_edit_distance, aligned_bases)}.

    ``hp_edit_distance`` is the HP-aware edit distance from ``_cigar_hp_edit_distance``.
    ``aligned_bases`` is the count of M/=/X CIGAR bases from ``_cigar_aligned_bases``.
    Both values are returned together to avoid scanning the BAM twice.
    """
    if _pysam is None:
        raise RuntimeError("pysam is required to compute HP edit distances from BAM")
    results: Dict[str, Tuple[float, int]] = {}
    try:
        with _pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam:
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                rid = _bare_uuid(read.query_name or "")
                results[rid] = (
                    _cigar_hp_edit_distance(read, genome, penalty_table),
                    _cigar_aligned_bases(read),
                )
    except Exception as exc:
        logger.warning("Failed to compute HP edit distances from %s: %s", bam_path, exc)
    return results


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


def _compute_junction_stats(all_df: "pd.DataFrame") -> Dict[Tuple[int, int], Tuple[int, int]]:
    """Build a per-junction support dictionary from the concatenated all_df.

    For each (intron_start, intron_end) pair, returns:
        (n_aligner_read_pairs, max_min_overhang_seen)

    where:
        n_aligner_read_pairs = how many (read_id, aligner) rows reported this junction
        max_min_overhang_seen = maximum of min(left_ov, right_ov) across all those rows

    This is used by the short-intron relaxation: if a short intron is reported
    by many reads and at least one has good overhang, it is almost certainly real.
    """
    junc_stats: Dict[Tuple[int, int], Tuple[int, int]] = {}  # {(js,je): (count, max_min_ov)}

    aln_start_col = all_df.get('alignment_start', pd.Series(0, index=all_df.index)).fillna(0)
    aln_end_col   = all_df.get('alignment_end',   pd.Series(0, index=all_df.index)).fillna(0)
    junc_col      = all_df.get('junctions',       pd.Series('',  index=all_df.index))

    for idx in all_df.index:
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
    """Strip mapPacBio's pt:i:N suffix from read IDs so all aligners share a common key.

    mapPacBio embeds the FASTQ header's auxiliary tag in the BAM read name.
    The separator is space before ``samtools sort`` and underscore after:
    - Pre-sort (live correction output): 'UUID pt:i:25'
    - Post-sort / merged BAMs: 'UUID_pt:i:25'

    All other aligners use just 'UUID'.  Without normalization, merge_corrected_tsvs
    treats these as different reads, producing ~50% logical duplicates in the merged
    output.
    """
    read_id_series = read_id_series.copy()
    # Space form (pre-samtools-sort): 'UUID pt:i:N'
    space_mask = read_id_series.str.contains(' pt:i:', na=False, regex=False)
    if space_mask.any():
        read_id_series[space_mask] = read_id_series[space_mask].str.split(' pt:i:').str[0]
    # Underscore form (post-samtools-sort): 'UUID_pt:i:N'
    under_mask = read_id_series.str.contains('_pt:i:', na=False, regex=False)
    if under_mask.any():
        read_id_series[under_mask] = read_id_series[under_mask].str.split('_pt:i:').str[0]
    return read_id_series


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
        import warnings as _warnings
        with _warnings.catch_warnings():
            _warnings.simplefilter('ignore')
            df = pd.read_csv(tsv_path, sep='\t', index_col=False)
        if df.empty:
            logger.warning("Empty per-aligner TSV, skipping: %s", tsv_path)
            return None
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

def merge_corrected_tsvs(
    per_aligner_tsvs: Dict[str, Path],
    output_tsv: Path,
    summary_tsv: Optional[Path] = None,
    per_aligner_corrected_bams: Optional[Dict[str, str]] = None,
    genome: Optional[Dict[str, str]] = None,
    penalty_table=None,
    overhang_table: Optional[OverhangTable] = None,
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

    # ── HP-aware edit distance (primary criterion when BAMs provided) ──────
    use_hp_ed = bool(per_aligner_corrected_bams)
    if use_hp_ed:
        logger.info("Computing HP-aware edit distances from corrected BAMs...")
        ed_rows = []
        for aligner_name, bam_path in (per_aligner_corrected_bams or {}).items():
            ed_dict = _read_hp_edit_distances(str(bam_path), genome, penalty_table)
            for rid, (ed, ab) in ed_dict.items():
                ed_rows.append({'read_id': rid, '_aligner': aligner_name,
                                 'hp_edit_distance': ed, 'aligned_bases': ab})
        if ed_rows:
            ed_df = pd.DataFrame(ed_rows)
            all_df = all_df.merge(ed_df, on=['read_id', '_aligner'], how='left')
            all_df['hp_edit_distance'] = all_df['hp_edit_distance'].fillna(float('inf'))
            all_df['aligned_bases']    = all_df['aligned_bases'].fillna(0).astype(int)
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
    junc_stats = _compute_junction_stats(all_df)
    chimera_flags = _add_chimera_flag(rep_df, _ov_table, junc_stats)
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
    # Within each read, group aligners by their (corrected_3prime,
    # corrected_junctions) tuple — aligners landing at the same biological
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
        _juncs = _row.get('junctions', '') or ''
        # Normalize: sort the donor-acceptor tuples so semicolon-order
        # differences across aligners don't shadow real equivalence.
        _parts = tuple(sorted(p for p in _juncs.split(';') if p))
        try:
            _c3p = int(_row['corrected_3prime'])
        except (TypeError, ValueError, KeyError):
            _c3p = -1
        return (_c3p, _parts)

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
