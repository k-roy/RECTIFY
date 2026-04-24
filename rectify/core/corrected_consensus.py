"""
corrected_consensus.py — Select best corrected 3' end positions across per-aligner outputs.

For each read, selects the winning aligner's corrected output from N per-aligner
corrected_3ends.tsv files.  The selection uses post-correction features
(five_prime_rescued, confidence, corrected_3prime agreement, alignment span,
junction count) rather than raw alignment features such as AS tag or soft-clip
length — which are not cross-comparable across aligners.

Also identifies Cat5 candidates: reads where ≥2 aligners each uniquely contribute
at least one correctly-rescued intron not present in the other aligner's result.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

logger = logging.getLogger(__name__)

_CONFIDENCE_RANK: Dict[str, int] = {'high': 2, 'medium': 1, 'low': 0}


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

    mapPacBio embeds the FASTQ header's auxiliary tag verbatim into the BAM read
    name as a space-separated suffix: 'UUID pt:i:25'.  All other aligners strip
    the suffix and use just 'UUID'.  Without normalization, merge_corrected_tsvs
    treats these as different reads, producing ~50% logical duplicates in the merged
    output.

    Note: the separator is a space character, not an underscore.
    """
    mask = read_id_series.str.contains(' pt:i:', na=False, regex=False)
    if mask.any():
        read_id_series = read_id_series.copy()
        read_id_series[mask] = read_id_series[mask].str.split(' pt:i:').str[0]
    return read_id_series


def _load_tsv(aligner_name: str, tsv_path: Path) -> Optional[pd.DataFrame]:
    """Load one per-aligner TSV, returning None on failure."""
    if not tsv_path.exists():
        logger.warning("Per-aligner TSV not found, skipping: %s", tsv_path)
        return None
    try:
        df = pd.read_csv(tsv_path, sep='\t')
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
) -> Path:
    """
    Merge N per-aligner corrected TSVs into a single corrected_3ends.tsv.

    For each read_id the winning aligner is chosen by (in priority order):

    1. ``five_prime_rescued`` — prefer aligner where Cat3 rescue fired (1 > 0)
    2. ``confidence``         — prefer 'high' > 'medium' > 'low'
    3. corrected_3prime agreement — prefer position agreed on by most aligners
    4. alignment span         — prefer wider alignment (alignment_end - alignment_start)
    5. ``n_junctions``        — prefer more junctions (more completely spliced)

    Multi-peak NET-seq reads (multiple rows per read_id with different
    ``fraction`` values) are handled correctly: the winning aligner is
    determined from a single representative row per read, then *all* rows
    from that aligner for that read_id are written to the output.

    Parameters
    ----------
    per_aligner_tsvs:
        Mapping from aligner name to corrected_3ends.tsv path.
    output_tsv:
        Destination path for the merged corrected_3ends.tsv.
    summary_tsv:
        Optional path for a per-read aligner comparison table (useful for
        Cat5 candidate inspection and pipeline debugging).

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
        out_df.to_csv(output_tsv, sep='\t', index=False)
        logger.info("Single aligner (%s) — wrote %d rows to %s",
                    aligner_name, len(out_df), output_tsv)
        return output_tsv

    all_df = pd.concat(dfs.values(), ignore_index=True)

    # ── Scoring columns ───────────────────────────────────────────────────
    all_df['_conf_rank'] = (
        all_df['confidence'].map(_CONFIDENCE_RANK).fillna(0).astype(int)
    )
    rescued_col = all_df.get('five_prime_rescued', pd.Series(0, index=all_df.index))
    all_df['_five_rescued'] = rescued_col.fillna(0).astype(int)

    aln_start = all_df.get('alignment_start', pd.Series(0, index=all_df.index)).fillna(0)
    aln_end   = all_df.get('alignment_end',   pd.Series(0, index=all_df.index)).fillna(0)
    all_df['_span'] = (aln_end - aln_start).astype(int)

    all_df['_n_junc'] = all_df.get('n_junctions', 0).fillna(0).astype(int)

    # Count how many aligners agree on each corrected_3prime per read
    pos_counts = (
        all_df.groupby(['read_id', 'corrected_3prime'])
        .size()
        .reset_index(name='_n_agree')
    )
    all_df = all_df.merge(pos_counts, on=['read_id', 'corrected_3prime'], how='left')

    # ── Representative row per (read_id, aligner) for ranking ────────────
    # For NET-seq multi-peak reads, use the row with the highest fraction as
    # the representative (it carries the primary corrected position).
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

    # Sort so the best candidate per read comes first
    rep_df = rep_df.sort_values(
        ['read_id', '_five_rescued', '_conf_rank', '_n_agree', '_span', '_n_junc'],
        ascending=[True, False, False, False, False, False],
    )

    # Pick winning aligner per read
    winner_cols = ['read_id', '_aligner']
    winner_df = (
        rep_df.groupby('read_id', sort=False)
        .first()
        .reset_index()[winner_cols]
        .rename(columns={'_aligner': '_winning_aligner'})
    )

    # ── Optional comparison summary ───────────────────────────────────────
    if summary_tsv:
        summary = rep_df.merge(
            winner_df, on='read_id', how='left'
        )
        summary['_is_winner'] = summary['_aligner'] == summary['_winning_aligner']
        keep = [
            'read_id', '_aligner', '_is_winner',
            'corrected_3prime', 'five_prime_position',
            'junctions', '_conf_rank', '_n_agree', '_n_junc', '_five_rescued',
            'correction_applied', 'confidence',
        ]
        keep = [c for c in keep if c in summary.columns]
        summary[keep].sort_values(['read_id', '_aligner']).to_csv(
            summary_tsv, sep='\t', index=False
        )
        logger.info("Comparison summary → %s", summary_tsv)

    # ── Select all rows from winning aligner (handles multi-peak reads) ──
    result_df = all_df.merge(winner_df, on='read_id', how='left')
    result_df = result_df[result_df['_aligner'] == result_df['_winning_aligner']].copy()

    drop_cols = [
        '_aligner', '_winning_aligner', '_conf_rank', '_five_rescued',
        '_span', '_n_junc', '_n_agree',
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
    corrected_3ends.tsv (post-correction, so Cat3-rescued junctions are included
    and false N-ops absorbed by Cat4 are excluded).

    Parameters
    ----------
    per_aligner_tsvs:
        Mapping from aligner name to corrected_3ends.tsv path.
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
