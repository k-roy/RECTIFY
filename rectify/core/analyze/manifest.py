"""
Manifest-mode runner for the RECTIFY analyze pipeline.

Memory-efficient multi-sample analysis: streams per-sample TSVs (via a
manifest mapping sample_id → TSV path), aggregates positions across all
samples, clusters once, then re-streams each sample to build the
clusters × samples count matrix without ever holding the combined
positions frame in memory. Downstream steps (PCA, DESeq2, GO, motifs,
shift, HTML report) run identically to the standard pipeline.
"""

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd

from .clustering import (
    cluster_cpa_sites,
    cluster_cpa_sites_adaptive,
    annotate_clusters_with_genes,
)
from .deseq2 import (
    extract_condition_from_sample,
    detect_control_samples,
    create_sample_metadata,
    run_deseq2_gene_level,
    run_deseq2_cluster_level,
)
from .pca import run_pca_analysis, plot_pca
from .heatmap import plot_sample_heatmap
from .go_enrichment import run_go_enrichment, plot_go_enrichment
from .motif_discovery import (
    run_differential_motif_analysis,
    summarize_motif_results,
)
from .shift_analysis import (
    analyze_cluster_shifts,
    get_top_shifted_genes,
    plot_gene_browser,
    plot_shift_summary,
)
from .summary import generate_analysis_summary, generate_html_report
from .genomic_distribution import (
    run_3prime_distribution_analysis,
    run_5prime_distribution_analysis,
    run_transcript_body_distribution_analysis,
)
from .loaders import load_annotation, load_position_index
from .exclusions import _detect_exclusion_regions


def _run_analyze_manifest(
    args: argparse.Namespace,
    output_dir: Path,
    plots_dir: Path,
    tables_dir: Path,
) -> int:
    """
    Memory-efficient manifest-mode analysis pipeline.

    Instead of loading a large combined TSV, processes each sample's corrected
    TSV individually using a two-pass approach:

    Pass 1 (clustering): Stream through each sample to aggregate position counts,
        then cluster the combined positions.

    Pass 2 (count matrix): Stream through each sample again, look up cluster IDs,
        and accumulate counts per cluster/sample.

    Downstream steps (PCA, DESeq2, GO, motifs, HTML report) run identically to
    the standard pipeline.

    Bedgraphs and genomic distribution are generated in a lightweight Pass 3 that
    streams only the needed columns (chrom, strand, position) one sample at a time,
    accumulating per-condition counts without ever holding the full combined dataset
    in memory.

    5' TSS clustering mirrors 3' CPA clustering: a parallel Pass 1b aggregates
    five_prime_position counts, then clusters them with a wider window (75 bp default).
    """
    from ...utils.chromosome import normalize_chromosome

    print("=" * 70)
    print("RECTIFY Analysis Pipeline (Manifest Mode)")
    print("=" * 70)

    # Load manifest
    manifest_df = pd.read_csv(args.manifest, sep='\t')
    required_manifest_cols = ['sample_id', 'path']
    missing_manifest_cols = [c for c in required_manifest_cols if c not in manifest_df.columns]
    if missing_manifest_cols:
        print(f"ERROR: Manifest missing columns: {missing_manifest_cols}", flush=True)
        return 1

    samples = manifest_df.to_dict('records')
    print(f"\nManifest: {len(samples)} samples")
    for s in samples:
        cond_str = f" [{s['condition']}]" if 'condition' in s else ''
        print(f"  {s['sample_id']}: {s['path']}{cond_str}")

    # Chromosome format settings: passthrough-by-default for species portability.
    # Set by --chrom-format CLI flag.
    chrom_format = getattr(args, 'chrom_format', 'passthrough')
    sample_column = getattr(args, 'sample_column', 'sample')

    # Load annotation early; thread chrom_format so its chrom names match
    # the cluster table (which is also normalized to chrom_format).
    annotation_df = None
    if args.annotation:
        annotation_df = load_annotation(args.annotation, chrom_format=chrom_format)

    # ─────────────────────────────────────────────────────────────────────────
    # Pass 1: Aggregate positions across all samples for clustering
    # ─────────────────────────────────────────────────────────────────────────
    print(f"\n[1/9] Aggregating positions across {len(samples)} samples (Pass 1)...")

    all_pos_dfs = []
    for s in samples:
        sample_id = s['sample_id']
        tsv_path = s['path']

        # Try loading from compact index first
        idx_df = load_position_index(tsv_path, sample_id, normalize_chroms=True, chrom_format=chrom_format)
        if idx_df is not None:
            n_ag = int(idx_df['count_ag_rich'].sum()) if 'count_ag_rich' in idx_df.columns else 0
            print(f"  {sample_id}: loaded index ({len(idx_df):,} positions, {n_ag:,} AG-rich reads tracked)")
            # Index has corrected_position, chrom, strand, count, count_ag_rich, sample
            all_pos_dfs.append(
                idx_df[['chrom', 'strand', 'corrected_position', 'count', 'count_ag_rich', 'sample']]
            )
            continue

        # Fall back: stream the full TSV with minimal columns
        print(f"  {sample_id}: streaming TSV (no index found)...")
        _header = pd.read_csv(tsv_path, sep='\t', nrows=0)
        _pos_col = None
        for _c in ['corrected_position', 'corrected_3prime', 'position']:
            if _c in _header.columns:
                _pos_col = _c
                break
        if _pos_col is None:
            print(f"  WARNING: No position column in {tsv_path}, skipping.")
            continue

        _has_qc = 'qc_flags' in _header.columns
        _usecols = ['chrom', 'strand', _pos_col]
        if _has_qc:
            _usecols.append('qc_flags')
        # Each value is [count_total, count_ag_rich]. AG_RICH reads are KEPT
        # in the cluster-discovery pool (v0.9.2+): the WT-vs-Ysh1AA contrast
        # showed condition-dependent AG enrichment reflecting real CPA-readthrough
        # biology, not protocol artefact. AG counts are tracked separately so
        # cluster-level annotation can carry an `ag_rich_fraction` flag for
        # downstream filtering (e.g. DRS-based orthogonal voting).
        _agg = defaultdict(lambda: [0.0, 0.0])
        _n_ag_rich = 0
        for _chunk in pd.read_csv(tsv_path, sep='\t', chunksize=100_000, usecols=_usecols):
            if _has_qc:
                _ag_mask = _chunk['qc_flags'].str.contains('AG_RICH', na=False)
                _n_ag_rich += int(_ag_mask.sum())
            else:
                _ag_mask = None
            _chunk['chrom'] = _chunk['chrom'].map(lambda x: normalize_chromosome(x, chrom_format))
            # Vectorized groupby (~20-50x faster than itertuples+dict for large chunks)
            for (chrom, strand, pos), cnt in (
                _chunk.groupby(['chrom', 'strand', _pos_col], sort=False).size().items()
            ):
                _agg[(chrom, strand, pos)][0] += float(cnt)
            if _ag_mask is not None and _ag_mask.any():
                for (chrom, strand, pos), cnt in (
                    _chunk.loc[_ag_mask]
                    .groupby(['chrom', 'strand', _pos_col], sort=False)
                    .size()
                    .items()
                ):
                    _agg[(chrom, strand, pos)][1] += float(cnt)

        if _has_qc and _n_ag_rich > 0:
            print(f"  {sample_id}: tracked {_n_ag_rich:,} AG-rich reads (kept in cluster pool)")

        if not _agg:
            print(f"  WARNING: No positions found in {tsv_path}, skipping.")
            continue

        _rows = [
            {
                'chrom': c, 'strand': st, 'corrected_position': pos,
                'count': cnt, 'count_ag_rich': cnt_ag, 'sample': sample_id,
            }
            for (c, st, pos), (cnt, cnt_ag) in _agg.items()
        ]
        all_pos_dfs.append(pd.DataFrame(_rows))
        print(f"  {sample_id}: aggregated {len(_rows):,} positions")

    if not all_pos_dfs:
        print("ERROR: No positions loaded from any sample.", flush=True)
        return 1

    positions_agg = pd.concat(all_pos_dfs, ignore_index=True)
    print(f"  Total aggregated positions: {len(positions_agg):,} across {positions_agg['sample'].nunique()} samples")

    # Filter exclusion regions
    exclude_mito = args.exclude_mito and not getattr(args, 'include_mito', False)
    exclude_rdna = getattr(args, 'exclude_rdna', True) and not getattr(args, 'include_rdna', False)
    mito_chroms: set = set()
    if exclude_mito or exclude_rdna:
        n_before = len(positions_agg)
        mito_chroms, rdna_regions = _detect_exclusion_regions(annotation_df, positions_agg['chrom'].unique())
        if mito_chroms:
            print(f"  Auto-detected mitochondrial chromosomes: {', '.join(sorted(mito_chroms))}")
        mito_mask = pd.Series(False, index=positions_agg.index)
        if exclude_mito and mito_chroms:
            mito_mask = positions_agg['chrom'].isin(mito_chroms)
        rdna_mask = pd.Series(False, index=positions_agg.index)
        if exclude_rdna and rdna_regions:
            for chrom, start, end in rdna_regions:
                rdna_mask = rdna_mask | (
                    (positions_agg['chrom'] == chrom) &
                    (positions_agg['corrected_position'] >= start) &
                    (positions_agg['corrected_position'] <= end)
                )
        positions_agg = positions_agg[~(mito_mask | rdna_mask)]
        n_removed = n_before - len(positions_agg)
        if n_removed > 0:
            print(f"  Excluded {n_removed:,} positions (mito/rDNA)")

    # Cluster CPA sites using aggregated data
    count_col = 'count'
    _max_radius = getattr(args, 'max_cluster_radius', 10)
    _min_peak_sep = getattr(args, 'min_peak_sep', 5)
    _min_cluster_samples = getattr(args, 'min_cluster_samples', 2)

    print(f"\n[2/9] Forming CPA clusters (adaptive, radius={_max_radius}bp, "
          f"min_peak_sep={_min_peak_sep}bp, min_samples={_min_cluster_samples})...")
    clusters_df = cluster_cpa_sites_adaptive(
        positions_agg,
        max_cluster_radius=_max_radius,
        min_peak_separation=_min_peak_sep,
        min_reads=args.min_reads,
        count_col=count_col,
    )
    print(f"  Formed {len(clusters_df):,} clusters")

    # Annotate clusters with AG-richness aggregates. Post-hoc join keeps the
    # clusterer signature stable; runs in O(positions + clusters) via interval
    # lookups grouped by (chrom, strand).
    if 'count_ag_rich' in positions_agg.columns and len(clusters_df) > 0:
        _pos_totals = (
            positions_agg.groupby(
                ['chrom', 'strand', 'corrected_position'], sort=False
            )[['count', 'count_ag_rich']]
            .sum()
            .reset_index()
        )
        # Per-cluster aggregate: sum count_ag_rich and count over positions in
        # [start, end] on the same (chrom, strand).
        _n_ag = []
        _n_tot = []
        _modal_ag = []
        # Group positions by (chrom, strand) once for cheap repeated lookups
        # via the position index inside each per-(chrom,strand) frame.
        _by_cs = {
            k: v.set_index('corrected_position')[['count', 'count_ag_rich']]
            for k, v in _pos_totals.groupby(['chrom', 'strand'])
        }
        for _, c in clusters_df.iterrows():
            sub = _by_cs.get((c['chrom'], c['strand']))
            if sub is None:
                _n_ag.append(0.0)
                _n_tot.append(0.0)
                _modal_ag.append(False)
                continue
            in_range = sub.loc[(sub.index >= c['start']) & (sub.index <= c['end'])]
            _n_ag.append(float(in_range['count_ag_rich'].sum()))
            _n_tot.append(float(in_range['count'].sum()))
            # modal_pos_ag_rich: at the modal_position, is the AG fraction > 0.5?
            modal_row = sub.loc[sub.index == c['modal_position']]
            if len(modal_row) > 0 and float(modal_row['count'].iloc[0]) > 0:
                _modal_ag.append(
                    float(modal_row['count_ag_rich'].iloc[0])
                    / float(modal_row['count'].iloc[0])
                    > 0.5
                )
            else:
                _modal_ag.append(False)
        clusters_df['n_reads_ag_rich'] = _n_ag
        clusters_df['ag_rich_fraction'] = [
            (a / t) if t > 0 else 0.0 for a, t in zip(_n_ag, _n_tot)
        ]
        clusters_df['modal_pos_ag_rich'] = _modal_ag
        _ag_threshold = getattr(args, 'cluster_ag_fraction_threshold', 0.5)
        clusters_df['cluster_ag_flag'] = (
            (clusters_df['ag_rich_fraction'] > _ag_threshold)
            | clusters_df['modal_pos_ag_rich']
        )
        n_ag_flagged = int(clusters_df['cluster_ag_flag'].sum())
        print(
            f"  AG-rich annotation: {n_ag_flagged:,} clusters flagged "
            f"(modal_pos AG_RICH OR cluster ag_rich_fraction > {_ag_threshold})"
        )

    if args.annotation and annotation_df is not None:
        print(f"  Annotating clusters with genes...")
        clusters_df = annotate_clusters_with_genes(clusters_df, annotation_df)
        n_annotated = clusters_df['gene_id'].notna().sum()
        print(f"  Annotated {n_annotated:,} clusters ({100*n_annotated/len(clusters_df):.1f}%)")

    clusters_path = output_dir / 'cpa_clusters.tsv'
    clusters_df.to_csv(clusters_path, sep='\t', index=False)
    print(f"  Saved clusters to {clusters_path}")

    # ─────────────────────────────────────────────────────────────────────────
    # Pass 1b: Aggregate 5' end (TSS) positions and cluster them
    # ─────────────────────────────────────────────────────────────────────────
    tss_distance = getattr(args, 'tss_cluster_distance', 75)
    tss_clusters_df = None
    tss_lookup = None
    print(f"\n[2b/9] Aggregating 5' TSS positions (Pass 1b)...")
    tss_agg_dict: defaultdict = defaultdict(float)
    for s in samples:
        sample_id = s['sample_id']
        tsv_path = s['path']
        _header = pd.read_csv(tsv_path, sep='\t', nrows=0)
        if 'five_prime_position' not in _header.columns:
            continue
        _usecols_tss = ['chrom', 'strand', 'five_prime_position']
        try:
            for _chunk in pd.read_csv(tsv_path, sep='\t', chunksize=100_000, usecols=_usecols_tss):
                _chunk = _chunk.dropna(subset=['five_prime_position'])
                if _chunk.empty:
                    continue
                _chunk['chrom'] = _chunk['chrom'].map(lambda x: normalize_chromosome(x, chrom_format))
                _chunk['five_prime_position'] = _chunk['five_prime_position'].astype(int)
                # Apply mito/rDNA exclusion if available
                if exclude_mito and mito_chroms:
                    _chunk = _chunk[~_chunk['chrom'].isin(mito_chroms)]
                for (chrom, strand, pos), cnt in (
                    _chunk.groupby(['chrom', 'strand', 'five_prime_position'], sort=False)
                    .size().items()
                ):
                    tss_agg_dict[(chrom, strand, pos)] += float(cnt)
        except Exception as _tss_agg_exc:
            print(f"  WARNING: TSS aggregation failed for {sample_id}: {_tss_agg_exc}")

    if tss_agg_dict:
        _tss_rows = [
            {'chrom': c, 'strand': st, 'five_prime_position': pos, 'count': cnt}
            for (c, st, pos), cnt in tss_agg_dict.items()
        ]
        tss_agg_df = pd.DataFrame(_tss_rows)
        print(f"  TSS: {len(tss_agg_df):,} unique positions aggregated across {len(samples)} samples")
        try:
            tss_clusters_df = cluster_cpa_sites(
                tss_agg_df,
                cluster_distance=tss_distance,
                min_reads=args.min_reads,
                position_col='five_prime_position',
                count_col='count',
            )
            if args.annotation and annotation_df is not None:
                tss_clusters_df = annotate_clusters_with_genes(tss_clusters_df, annotation_df)
            tss_clusters_path = output_dir / 'tss_clusters.tsv'
            tss_clusters_df.to_csv(tss_clusters_path, sep='\t', index=False)
            print(f"  Formed {len(tss_clusters_df):,} TSS clusters → {tss_clusters_path}")

            # Build TSS cluster lookup (same pattern as CPA lookup below)
            try:
                from intervaltree import IntervalTree as _IT
                _tss_trees: dict = {}
                for _, _row in tss_clusters_df.iterrows():
                    _k = (_row['chrom'], _row['strand'])
                    if _k not in _tss_trees:
                        _tss_trees[_k] = _IT()
                    _tss_trees[_k][_row['start']:_row['end'] + 1] = _row['cluster_id']

                def tss_lookup(chrom, strand, pos):
                    _hits = _tss_trees.get((chrom, strand), _IT())[pos]
                    return _hits.pop().data if _hits else None
            except ImportError:
                _tss_lists: dict = {}
                for _, _row in tss_clusters_df.iterrows():
                    _tss_lists.setdefault((_row['chrom'], _row['strand']), []).append(
                        (_row['start'], _row['end'], _row['cluster_id'])
                    )

                def tss_lookup(chrom, strand, pos):
                    for _s, _e, _cid in _tss_lists.get((chrom, strand), []):
                        if _s <= pos <= _e:
                            return _cid
                    return None
        except Exception as _tss_cl_exc:
            print(f"  Warning: TSS clustering failed: {_tss_cl_exc}")
    else:
        print(f"  Note: five_prime_position column not found in any sample TSV — TSS clustering skipped.")

    # Build cluster lookup structure
    try:
        from intervaltree import IntervalTree
        trees = {}
        for _, row in clusters_df.iterrows():
            key = (row['chrom'], row['strand'])
            if key not in trees:
                trees[key] = IntervalTree()
            trees[key][row['start']:row['end'] + 1] = row['cluster_id']

        def lookup(chrom, strand, pos):
            hits = trees.get((chrom, strand), IntervalTree())[pos]
            return hits.pop().data if hits else None
    except ImportError:
        from bisect import bisect_left  # noqa: F401
        _lists = {}
        for _, row in clusters_df.iterrows():
            key = (row['chrom'], row['strand'])
            _lists.setdefault(key, []).append((row['start'], row['end'], row['cluster_id']))
        for k in _lists:
            _lists[k].sort()

        def lookup(chrom, strand, pos):
            for s, e, cid in _lists.get((chrom, strand), []):
                if s <= pos <= e:
                    return cid
            return None

    # ─────────────────────────────────────────────────────────────────────────
    # Pass 2: Build CPA count matrix + TSS count matrix by streaming each sample
    # Also accumulate per-condition bedgraph position counts.
    # ─────────────────────────────────────────────────────────────────────────
    print(f"\n[3/9] Building cluster count matrix (Pass 2)...")

    count_accumulator: defaultdict = defaultdict(lambda: defaultdict(float))
    tss_count_accumulator: defaultdict = defaultdict(lambda: defaultdict(float))

    # Bedgraph accumulators: {condition: {strand: {chrom: {pos: count}}}}
    _bedgraph_acc: dict = defaultdict(lambda: {'+': defaultdict(lambda: defaultdict(int)),
                                               '-': defaultdict(lambda: defaultdict(int))})
    _condition_totals: dict = defaultdict(int)

    # Condition mapping for bedgraph (built lazily per sample_id)
    def _sample_to_condition(sid: str) -> str:
        if 'condition' in manifest_df.columns:
            _row = manifest_df[manifest_df['sample_id'] == sid]
            if not _row.empty:
                return str(_row.iloc[0]['condition'])
        return extract_condition_from_sample(sid)

    for s in samples:
        sample_id = s['sample_id']
        tsv_path = s['path']
        _condition = _sample_to_condition(sample_id)

        # Try index first (CPA counts only; no bedgraph/TSS columns available)
        idx_df = load_position_index(tsv_path, sample_id, normalize_chroms=True, chrom_format=chrom_format)
        if idx_df is not None:
            _idx_total = float(idx_df['count'].sum())
            _idx_assigned = 0.0
            for row in idx_df.itertuples(index=False):
                cid = lookup(row.chrom, row.strand, int(row.corrected_position))
                if cid is not None:
                    count_accumulator[cid][sample_id] += float(row.count)
                    _idx_assigned += float(row.count)
                # Bedgraph from index
                _bedgraph_acc[_condition][row.strand][row.chrom][int(row.corrected_position)] += float(row.count)
                _condition_totals[_condition] += float(row.count)
            print(f"  {sample_id}: {_idx_assigned:,.0f}/{_idx_total:,.0f} reads assigned to clusters (index)")
            # TSS: not available from index — will be handled by full-TSV path below if no index
            continue

        # Stream full TSV — load all needed columns in one pass
        _header = pd.read_csv(tsv_path, sep='\t', nrows=0)
        _pos_col = None
        for _c in ['corrected_position', 'corrected_3prime', 'position']:
            if _c in _header.columns:
                _pos_col = _c
                break
        if _pos_col is None:
            print(f"  WARNING: No position column in {tsv_path}, skipping.")
            continue

        _has_fraction = 'fraction' in _header.columns
        _has_tss = 'five_prime_position' in _header.columns
        _usecols = ['chrom', 'strand', _pos_col]
        if _has_fraction:
            _usecols.append('fraction')
        if _has_tss:
            _usecols.append('five_prime_position')

        # AG_RICH reads are KEPT in cluster counts (v0.9.2+). The WT-vs-Ysh1AA
        # contrast shows condition-dependent AG-richness reflecting real
        # CPA-readthrough biology rather than internal-priming artefact;
        # downstream consumers can filter on the cluster-level `cluster_ag_flag`
        # column in cpa_clusters.tsv (or use DRS orthogonal voting).
        n_assigned = 0
        for _chunk in pd.read_csv(tsv_path, sep='\t', chunksize=100_000, usecols=_usecols):
            _chunk['chrom'] = _chunk['chrom'].map(lambda x: normalize_chromosome(x, chrom_format))
            weights = _chunk['fraction'].values if _has_fraction else None

            # CPA cluster assignment
            cids = [lookup(c, st, int(p)) for c, st, p in
                    zip(_chunk['chrom'], _chunk['strand'], _chunk[_pos_col])]
            for i, cid in enumerate(cids):
                if cid is not None:
                    count_accumulator[cid][sample_id] += (weights[i] if weights is not None else 1.0)
                    n_assigned += 1

            # Bedgraph accumulation (3' corrected positions, per-condition)
            for chrom, strand, pos in zip(_chunk['chrom'], _chunk['strand'], _chunk[_pos_col]):
                _bedgraph_acc[_condition][strand][chrom][int(pos)] += 1
                _condition_totals[_condition] += 1

            # TSS cluster assignment
            if _has_tss and tss_lookup is not None:
                _tss_chunk = _chunk.dropna(subset=['five_prime_position'])
                for chrom, strand, pos in zip(
                    _tss_chunk['chrom'], _tss_chunk['strand'], _tss_chunk['five_prime_position']
                ):
                    tcid = tss_lookup(chrom, strand, int(pos))
                    if tcid is not None:
                        tss_count_accumulator[tcid][sample_id] += 1.0

        print(f"  {sample_id}: {n_assigned:,} reads assigned to clusters")

    # Convert accumulator to DataFrame (clusters × samples)
    if not count_accumulator:
        print("ERROR: No reads could be assigned to clusters.", flush=True)
        return 1

    all_sample_ids = [s['sample_id'] for s in samples]
    all_cluster_ids = clusters_df['cluster_id'].tolist()

    # Ensure every sample is present for every cluster, even with zero counts.
    # A sample that had no corrected positions mapping to any cluster would
    # otherwise be absent from count_accumulator entirely, causing DESeq2 to
    # receive a matrix with wrong dimensions.
    for cluster_id in all_cluster_ids:
        for sample_id in all_sample_ids:
            count_accumulator[cluster_id][sample_id] = count_accumulator[cluster_id].get(sample_id, 0)

    count_matrix = pd.DataFrame(count_accumulator).T.fillna(0)
    # Ensure all samples present as columns (belt-and-suspenders guard)
    for sid in all_sample_ids:
        if sid not in count_matrix.columns:
            count_matrix[sid] = 0.0
    count_matrix = count_matrix[all_sample_ids]
    count_matrix.index.name = 'cluster_id'
    print(f"  Matrix shape: {count_matrix.shape[0]:,} clusters × {count_matrix.shape[1]} samples")

    counts_path = output_dir / 'cluster_counts.tsv'
    count_matrix.to_csv(counts_path, sep='\t')

    # ─────────────────────────────────────────────────────────────────────────
    # TSS count matrix
    # ─────────────────────────────────────────────────────────────────────────
    if tss_clusters_df is not None and tss_count_accumulator:
        _tss_all_cluster_ids = tss_clusters_df['cluster_id'].tolist()
        for _cid in _tss_all_cluster_ids:
            for _sid in all_sample_ids:
                tss_count_accumulator[_cid][_sid] = tss_count_accumulator[_cid].get(_sid, 0)
        tss_count_matrix = pd.DataFrame(tss_count_accumulator).T.fillna(0)
        for _sid in all_sample_ids:
            if _sid not in tss_count_matrix.columns:
                tss_count_matrix[_sid] = 0.0
        tss_count_matrix = tss_count_matrix[all_sample_ids]
        tss_count_matrix.index.name = 'cluster_id'
        tss_counts_path = output_dir / 'tss_cluster_counts.tsv'
        tss_count_matrix.to_csv(tss_counts_path, sep='\t')
        print(f"\n  TSS count matrix: {tss_count_matrix.shape[0]:,} clusters × {tss_count_matrix.shape[1]} samples → {tss_counts_path}")

    # ─────────────────────────────────────────────────────────────────────────
    # Bedgraph output (from accumulated per-condition counts)
    # ─────────────────────────────────────────────────────────────────────────
    if not getattr(args, 'no_bedgraph', False) and _bedgraph_acc:
        import os as _os
        bedgraph_dir = Path(args.bedgraph_dir) if getattr(args, 'bedgraph_dir', None) else output_dir / 'bedgraph'
        bedgraph_dir.mkdir(parents=True, exist_ok=True)
        CHROM_ORDER = [
            'chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII',
            'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrmt',
        ]
        print(f"\n[Bedgraph] Writing strand-specific bedgraphs for {len(_bedgraph_acc)} conditions...")
        for _cond, _strand_dict in _bedgraph_acc.items():
            _total = _condition_totals.get(_cond, 0)
            _rpm = 1e6 / _total if _total > 0 else 1.0
            for _strand_name, _chrom_dict in _strand_dict.items():
                _out = bedgraph_dir / f"{_cond}_{'+' if _strand_name == '+' else 'minus' if _strand_name == '-' else _strand_name}.bedgraph"
                _strand_label = 'plus' if _strand_name == '+' else 'minus'
                _out = bedgraph_dir / f"{_cond}_{_strand_label}.bedgraph"
                _tmp = _out.with_suffix('.tmp')
                try:
                    with open(_tmp, 'w') as _f:
                        _f.write(f'track type=bedGraph name="{_cond}_{_strand_label}" '
                                 f'description="RECTIFY 3\' ends ({_strand_label} strand)"\n')
                        for _chrom in CHROM_ORDER:
                            if _chrom not in _chrom_dict:
                                continue
                            for _pos in sorted(_chrom_dict[_chrom]):
                                _val = _chrom_dict[_chrom][_pos] * _rpm
                                _f.write(f"{_chrom}\t{int(_pos)-1}\t{int(_pos)}\t{_val:.4f}\n")
                    _os.rename(_tmp, _out)
                except Exception as _bg_exc:
                    if _tmp.exists():
                        _tmp.unlink()
                    print(f"  Warning: bedgraph write failed for {_cond}/{_strand_label}: {_bg_exc}")
            print(f"  {_cond}: {int(_total):,} reads → {bedgraph_dir}")

    # ─────────────────────────────────────────────────────────────────────────
    # Genomic distribution (Pass 3 — stream needed columns only)
    # ─────────────────────────────────────────────────────────────────────────
    if annotation_df is not None and not getattr(args, 'no_genomic_distribution', False):
        print(f"\n[Genomic Distribution] Running per-condition analysis (Pass 3)...")
        _dist_cols_needed = set()
        _dist_cols_needed.update(['chrom', 'strand'])
        # We'll load corrected_3prime, five_prime_position, alignment_start, alignment_end
        _GD_LOAD_COLS = ['chrom', 'strand', 'corrected_3prime', 'corrected_position',
                         'five_prime_position', 'alignment_start', 'alignment_end', 'fraction']

        # Accumulate per-condition lightweight DataFrames
        _gd_chunks: dict = defaultdict(list)  # condition → list of dfs
        for s in samples:
            sample_id = s['sample_id']
            tsv_path = s['path']
            _condition = _sample_to_condition(sample_id)
            _hdr = pd.read_csv(tsv_path, sep='\t', nrows=0)
            _load = [c for c in _GD_LOAD_COLS if c in _hdr.columns]
            if not _load or 'chrom' not in _load:
                continue
            try:
                _df = pd.read_csv(tsv_path, sep='\t', usecols=_load)
                _df['chrom'] = _df['chrom'].map(lambda x: normalize_chromosome(x, chrom_format))
                _df['_sample'] = sample_id
                _gd_chunks[_condition].append(_df)
            except Exception as _gd_exc:
                print(f"  Warning: could not load distribution data for {sample_id}: {_gd_exc}")

        for _cond, _dfs in _gd_chunks.items():
            if not _dfs:
                continue
            _cond_df = pd.concat(_dfs, ignore_index=True)
            _pos_col_gd = 'corrected_3prime' if 'corrected_3prime' in _cond_df.columns else 'corrected_position'
            _count_col_gd = 'fraction' if 'fraction' in _cond_df.columns else None
            _cond_df['sample'] = _cond

            try:
                run_3prime_distribution_analysis(
                    _cond_df,
                    annotation_df,
                    output_dir=str(plots_dir),
                    sample_column='sample',
                    position_column=_pos_col_gd,
                    count_column=_count_col_gd,
                )
            except Exception as _d3_exc:
                print(f"  Warning: 3' distribution failed for {_cond}: {_d3_exc}")

            if 'alignment_start' in _cond_df.columns and 'alignment_end' in _cond_df.columns:
                try:
                    run_transcript_body_distribution_analysis(
                        _cond_df,
                        annotation_df,
                        output_dir=str(plots_dir),
                        sample_column='sample',
                        start_column='alignment_start',
                        end_column='alignment_end',
                    )
                except Exception as _db_exc:
                    print(f"  Warning: body distribution failed for {_cond}: {_db_exc}")

            if 'five_prime_position' in _cond_df.columns:
                try:
                    run_5prime_distribution_analysis(
                        _cond_df,
                        annotation_df,
                        output_dir=str(plots_dir),
                        sample_column='sample',
                        position_column='five_prime_position',
                        count_column=_count_col_gd,
                    )
                except Exception as _d5_exc:
                    print(f"  Warning: 5' distribution failed for {_cond}: {_d5_exc}")

    # ─────────────────────────────────────────────────────────────────────────
    # From here: same as run_analyze() starting from sample metadata
    # ─────────────────────────────────────────────────────────────────────────

    sample_names = count_matrix.columns.tolist()

    if args.reference:
        reference_condition = args.reference
    else:
        control_samples = detect_control_samples(sample_names)
        if control_samples:
            reference_condition = extract_condition_from_sample(control_samples[0])
            print(f"  Auto-detected reference condition: {reference_condition}")
        else:
            print("  Warning: Could not auto-detect reference. Using first condition.")
            reference_condition = extract_condition_from_sample(sample_names[0])

    sample_metadata = create_sample_metadata(sample_names, control_samples if not args.reference else None)
    sample_metadata.to_csv(output_dir / 'sample_metadata.tsv', sep='\t')

    # If manifest has a condition column, use it to override auto-detected conditions
    if 'condition' in manifest_df.columns:
        cond_map = dict(zip(manifest_df['sample_id'], manifest_df['condition']))
        sample_metadata['condition'] = sample_metadata.index.map(lambda x: cond_map.get(x, sample_metadata.loc[x, 'condition']))
        _manifest_conditions = manifest_df['condition'].unique().tolist()
        if args.reference:
            # Case-insensitive match of --reference against actual manifest conditions
            _ref_lower = args.reference.lower()
            _matched = [c for c in _manifest_conditions if c.lower() == _ref_lower]
            if _matched:
                reference_condition = _matched[0]
            else:
                print(f"  Warning: --reference '{args.reference}' not found in manifest conditions {_manifest_conditions}. Using as-is.")
        else:
            # Auto-detect reference from manifest conditions (wt/ctrl/control keywords)
            _control_conditions = [c for c in _manifest_conditions if 'wt' in c.lower() or 'ctrl' in c.lower() or 'control' in c.lower()]
            if _control_conditions:
                reference_condition = _control_conditions[0]
                print(f"  Reference condition from manifest: {reference_condition}")
            else:
                reference_condition = _manifest_conditions[0]

    # PCA
    print(f"\n[4/9] Running PCA analysis...")
    sample_sets = None
    if args.sample_sets:
        import json
        if args.sample_sets.startswith('{'):
            sample_sets = json.loads(args.sample_sets)
        elif Path(args.sample_sets).exists():
            with open(args.sample_sets) as f:
                sample_sets = json.load(f)

    if sample_sets:
        for set_name, conditions in sample_sets.items():
            set_samples = sample_metadata[sample_metadata['condition'].isin(conditions)].index.tolist()
            set_samples = [s for s in set_samples if s in count_matrix.columns]
            if len(set_samples) < 3:
                print(f"  Skipping PCA for {set_name} (only {len(set_samples)} samples)")
                continue
            set_matrix = count_matrix[set_samples]
            set_meta = sample_metadata.loc[sample_metadata.index.isin(set_samples)]
            pca_results = run_pca_analysis(set_matrix)
            if pca_results['pca_coords'] is not None and not pca_results['pca_coords'].empty:
                pca_path = plots_dir / f'pca_{set_name}.png'
                plot_pca(pca_results, sample_metadata=set_meta, color_by='condition',
                         output_path=str(pca_path), title=f'PCA: {set_name}')
                print(f"  Saved PCA plot for {set_name}: {pca_path}")
    else:
        pca_results = run_pca_analysis(count_matrix)
        if pca_results['pca_coords'] is not None and not pca_results['pca_coords'].empty:
            pca_path = plots_dir / 'pca_samples.png'
            plot_pca(pca_results, sample_metadata=sample_metadata, color_by='condition',
                     output_path=str(pca_path), title='Sample PCA (CPA Clusters)')
            print(f"  Saved PCA plot to {pca_path}")

    # Sample heatmap
    print(f"\n[5/9] Creating sample clustering heatmap...")
    if count_matrix.shape[1] < 2:
        print(f"  Skipping heatmap — requires ≥2 samples")
    else:
        heatmap_path = plots_dir / 'sample_heatmap.png'
        plot_sample_heatmap(count_matrix, sample_metadata=sample_metadata,
                            color_by='condition', output_path=str(heatmap_path))
        print(f"  Saved heatmap to {heatmap_path}")

    # DESeq2
    deseq2_gene_results = {}
    deseq2_cluster_results = {}

    if args.run_deseq2:
        print(f"\n[6/9] Running DESeq2 differential expression...")

        # Guard: DESeq2 requires at least 2 distinct conditions
        _n_conditions = sample_metadata['condition'].nunique()
        if _n_conditions < 2:
            print(f"  Warning: DESeq2 requires ≥2 conditions but only {_n_conditions} found. Skipping.")
        else:
            print("  Gene-level analysis...")
            deseq2_gene_results = run_deseq2_gene_level(
                count_matrix, clusters_df, sample_metadata, reference_condition, n_cpus=args.threads)
            for condition, result_df in deseq2_gene_results.items():
                result_df.to_csv(tables_dir / f'deseq2_genes_{condition}.tsv', sep='\t')
                n_sig = (result_df['padj'] < 0.05).sum()
                print(f"    {condition}: {n_sig:,} significant genes")

            import gc
            gc.collect()

            print("  Cluster-level analysis...")
            deseq2_cluster_results = run_deseq2_cluster_level(
                count_matrix, clusters_df, sample_metadata, reference_condition, n_cpus=args.threads)
            for condition, result_df in deseq2_cluster_results.items():
                result_df.to_csv(tables_dir / f'deseq2_clusters_{condition}.tsv', sep='\t')
                n_sig = (result_df['padj'] < 0.05).sum()
                print(f"    {condition}: {n_sig:,} significant clusters")

            gc.collect()
    else:
        print(f"\n[6/9] Skipping DESeq2 (use --run-deseq2 to enable)")

    # GO enrichment
    if args.go_annotations and deseq2_gene_results:
        print(f"\n[7/9] Running GO enrichment analysis...")
        from .go_enrichment import load_go_annotations
        go_annotations = load_go_annotations(
            args.go_annotations,
            gene_col='gene_name',
            go_term_col='description',
            category_col='go_category',
        )
        # Build systematic → common name mapping so DESeq2 indices match GO file
        _sys2common = {}
        if 'gene_id' in clusters_df.columns and 'gene_name' in clusters_df.columns:
            for _, _r in clusters_df[['gene_id', 'gene_name']].dropna().drop_duplicates().iterrows():
                _sys2common[_r['gene_id']] = _r['gene_name']

        def _to_common_names(genes):
            return [_sys2common.get(g, g) for g in genes]

        for condition, result_df in deseq2_gene_results.items():
            up_genes = _to_common_names(result_df[(result_df['padj'] < 0.05) & (result_df['log2FoldChange'] > 1)].index.tolist())
            if len(up_genes) >= 10:
                go_up = run_go_enrichment(up_genes, go_annotations)
                if not go_up.empty:
                    go_up.to_csv(tables_dir / f'go_enrichment_up_{condition}.tsv', sep='\t', index=False)
                    plot_go_enrichment(go_up, output_path=str(plots_dir / f'go_enrichment_up_{condition}.png'),
                                       title=f'GO Enrichment: Upregulated in {condition}')
            down_genes = _to_common_names(result_df[(result_df['padj'] < 0.05) & (result_df['log2FoldChange'] < -1)].index.tolist())
            if len(down_genes) >= 10:
                go_down = run_go_enrichment(down_genes, go_annotations)
                if not go_down.empty:
                    go_down.to_csv(tables_dir / f'go_enrichment_down_{condition}.tsv', sep='\t', index=False)
                    plot_go_enrichment(go_down, output_path=str(plots_dir / f'go_enrichment_down_{condition}.png'),
                                       title=f'GO Enrichment: Downregulated in {condition}')
    else:
        print(f"\n[7/9] Skipping GO enrichment (provide --go-annotations)")

    # Motif discovery
    if args.genome and args.run_motif and deseq2_cluster_results:
        print(f"\n[8/9] Running de novo motif discovery...")
        for condition, result_df in deseq2_cluster_results.items():
            enriched = result_df[(result_df['padj'] < 0.05) & (result_df['log2FoldChange'] > 1)]
            depleted = result_df[(result_df['padj'] < 0.05) & (result_df['log2FoldChange'] < -1)]
            if len(enriched) >= 20 and len(depleted) >= 20:
                motif_results = run_differential_motif_analysis(
                    enriched, depleted, args.genome,
                    str(output_dir / 'motifs' / condition),
                    upstream_window=args.motif_upstream,
                    downstream_window=args.motif_downstream,
                )
                summary = summarize_motif_results(motif_results)
                summary.to_csv(tables_dir / f'motif_summary_{condition}.tsv', sep='\t', index=False)
                print(f"  {condition}: {summary['motifs_found'].sum()} motifs found")
    else:
        print(f"\n[8/9] Skipping motif discovery (provide --genome and --run-motif)")

    # Shift analysis
    shift_results = None
    if len(sample_metadata['condition'].unique()) >= 2:
        print(f"\n[9/9] Running cluster shift analysis...")
        conditions = [c for c in sample_metadata['condition'].unique() if c != reference_condition]
        for condition in conditions:
            shift_df = analyze_cluster_shifts(
                count_matrix, clusters_df, reference_condition, condition, sample_metadata)
            if not shift_df.empty:
                shift_df.to_csv(tables_dir / f'shift_analysis_{condition}.tsv', sep='\t', index=False)
                plot_shift_summary(shift_df, output_path=str(plots_dir / f'shift_summary_{condition}.png'))
                top_shifted = get_top_shifted_genes(shift_df, n_top=20)
                for _, row in top_shifted.head(5).iterrows():
                    plot_gene_browser(
                        row['gene_id'], count_matrix, clusters_df, sample_metadata,
                        [reference_condition, condition],
                        output_path=str(plots_dir / f'browser_{row["gene_name"]}_{condition}.png'),
                    )
                print(f"  {condition}: {len(shift_df)} genes analyzed, "
                      f"{(shift_df['distribution_divergence'] > 0.2).sum()} with large shifts")
                shift_results = shift_df
    else:
        print(f"\n[9/9] Skipping shift analysis (need >=2 conditions)")

    # Summary
    print(f"\n[Summary] Generating report...")
    n_genes = clusters_df['gene_id'].nunique() if 'gene_id' in clusters_df.columns else 0
    summary_df = generate_analysis_summary(
        n_samples=len(sample_names),
        n_clusters=len(clusters_df),
        n_genes=n_genes,
        deseq2_gene_results=deseq2_gene_results,
        deseq2_cluster_results=deseq2_cluster_results,
        reference_condition=reference_condition,
    )
    summary_df.to_csv(output_dir / 'analysis_summary.tsv', sep='\t', index=False)

    plots_dict = {p.stem: str(p) for p in plots_dir.glob('*.png')}
    tables_dict = {t.stem: str(t) for t in tables_dir.glob('*.tsv')}
    html_path = output_dir / 'report.html'
    generate_html_report(summary_df, plots_dict, tables_dict, str(html_path),
                         title='RECTIFY Analysis Report (Manifest Mode)')

    print(f"\n" + "=" * 70)
    print(f"Manifest-mode analysis complete!")
    print(f"  Output directory: {output_dir}")
    print(f"  HTML report: {html_path}")
    print("=" * 70)

    return 0
