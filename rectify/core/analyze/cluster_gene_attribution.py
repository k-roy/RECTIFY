"""
Weighted CPA-cluster to gene attribution helpers.

The cluster table can only hold one ``gene_id`` per CPA cluster, but real 3'
clusters can be shared by reads from multiple upstream transcription origins.
This module builds a separate weighted mapping table:

    cluster_id  gene_id  gene_name  attribution_weight  source

Weights are normalized within each cluster and can be used by gene-level DESeq2
and shift analysis without collapsing shared clusters to an arbitrary nearest
annotated TES.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Callable, Iterable, Optional, Union

import pandas as pd


AttributionLookup = Callable[[str, str, int], Optional[str]]


def annotation_cluster_attributions(
    clusters_df: pd.DataFrame,
    source: str = "annotation",
) -> pd.DataFrame:
    """Build one-row-per-annotated-cluster attribution from ``clusters_df``."""
    if clusters_df.empty or "gene_id" not in clusters_df.columns:
        return _empty_attribution_df()

    rows = []
    for row in clusters_df.itertuples(index=False):
        gene_id = getattr(row, "gene_id", None)
        if _is_missing_gene(gene_id):
            continue
        rows.append({
            "cluster_id": row.cluster_id,
            "gene_id": str(gene_id),
            "gene_name": str(getattr(row, "gene_name", gene_id) or gene_id),
            "raw_attributed_count": float(getattr(row, "n_reads", 1.0) or 1.0),
            "attribution_weight": 1.0,
            "source": source,
        })
    return pd.DataFrame(rows) if rows else _empty_attribution_df()


def load_reference_position_attributions(
    paths: Iterable[Union[str, Path]],
    *,
    chrom_normalizer: Optional[Callable[[str], str]] = None,
    chunk_size: int = 1_000_000,
    include_intergenic: bool = False,
) -> pd.DataFrame:
    """Load one or more per-position DRS attribution TSVs.

    Expected columns are those emitted by
    ``analyses/.../scripts/80_attribute_drs_reads.py``:
    ``chrom, position, strand, gene_id, attributed_count``.
    """
    required = {"chrom", "position", "strand", "gene_id", "attributed_count"}
    frames = []
    for path in paths:
        p = Path(path)
        if p.is_dir():
            files = sorted(p.glob("*.tsv"))
        else:
            files = [p]
        for file_path in files:
            header = pd.read_csv(file_path, sep="\t", nrows=0)
            missing = required - set(header.columns)
            if missing:
                raise ValueError(f"{file_path} missing attribution columns: {sorted(missing)}")
            for chunk in pd.read_csv(file_path, sep="\t", chunksize=chunk_size, usecols=list(required)):
                if not include_intergenic:
                    chunk = chunk[chunk["gene_id"].notna() & (chunk["gene_id"] != "intergenic")]
                if chunk.empty:
                    continue
                if chrom_normalizer is not None:
                    chunk["chrom"] = chunk["chrom"].map(chrom_normalizer)
                chunk["position"] = chunk["position"].astype(int)
                chunk["attributed_count"] = chunk["attributed_count"].astype(float)
                grouped = (
                    chunk.groupby(["chrom", "strand", "position", "gene_id"], sort=False)["attributed_count"]
                    .sum()
                    .reset_index()
                )
                frames.append(grouped)

    if not frames:
        return pd.DataFrame(columns=["chrom", "strand", "position", "gene_id", "attributed_count"])

    combined = pd.concat(frames, ignore_index=True)
    return (
        combined.groupby(["chrom", "strand", "position", "gene_id"], sort=False)["attributed_count"]
        .sum()
        .reset_index()
    )


def reference_position_attributions_to_clusters(
    position_attributions: pd.DataFrame,
    cluster_lookup: AttributionLookup,
    *,
    source: str = "reference",
) -> pd.DataFrame:
    """Map per-position gene attributions onto the current cluster set."""
    if position_attributions.empty:
        return _empty_attribution_df()

    agg: defaultdict[tuple[str, str], float] = defaultdict(float)
    for row in position_attributions.itertuples(index=False):
        cluster_id = cluster_lookup(row.chrom, row.strand, int(row.position))
        if cluster_id is None:
            continue
        gene_id = str(row.gene_id)
        if _is_missing_gene(gene_id):
            continue
        agg[(cluster_id, gene_id)] += float(row.attributed_count)

    rows = [
        {
            "cluster_id": cluster_id,
            "gene_id": gene_id,
            "gene_name": gene_id,
            "raw_attributed_count": count,
            "source": source,
        }
        for (cluster_id, gene_id), count in agg.items()
        if count > 0
    ]
    return normalize_cluster_gene_attributions(pd.DataFrame(rows)) if rows else _empty_attribution_df()


def build_reference_cluster_lookup(
    clusters_df: pd.DataFrame,
    *,
    max_distance: int = 0,
) -> AttributionLookup:
    """Build a cluster lookup with optional nearest-cluster tolerance.

    Exact cluster overlap is always preferred. When ``max_distance`` is > 0,
    reference positions just outside a current cluster can still assign to the
    nearest same-strand cluster within that many bases.
    """
    max_distance = max(0, int(max_distance or 0))

    try:
        from intervaltree import IntervalTree

        trees = {}
        for _, row in clusters_df.iterrows():
            key = (row["chrom"], row["strand"])
            trees.setdefault(key, IntervalTree())[
                int(row["start"]) - max_distance:int(row["end"]) + max_distance + 1
            ] = {
                "cluster_id": row["cluster_id"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "modal_position": int(row.get("modal_position", row["start"])),
            }

        def lookup(chrom, strand, pos):
            hits = trees.get((chrom, strand), IntervalTree())[int(pos)]
            if not hits:
                return None
            best = min(
                (iv.data for iv in hits),
                key=lambda d: (
                    _distance_to_interval(int(pos), d["start"], d["end"]),
                    abs(int(pos) - d["modal_position"]),
                    d["start"],
                    str(d["cluster_id"]),
                ),
            )
            return best["cluster_id"]

        return lookup
    except ImportError:
        lists = {}
        for _, row in clusters_df.iterrows():
            lists.setdefault((row["chrom"], row["strand"]), []).append({
                "cluster_id": row["cluster_id"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "modal_position": int(row.get("modal_position", row["start"])),
            })

        def lookup(chrom, strand, pos):
            candidates = []
            for cluster in lists.get((chrom, strand), []):
                distance = _distance_to_interval(int(pos), cluster["start"], cluster["end"])
                if distance <= max_distance:
                    candidates.append(cluster)
            if not candidates:
                return None
            best = min(
                candidates,
                key=lambda d: (
                    _distance_to_interval(int(pos), d["start"], d["end"]),
                    abs(int(pos) - d["modal_position"]),
                    d["start"],
                    str(d["cluster_id"]),
                ),
            )
            return best["cluster_id"]

        return lookup


def body_attributions_from_corrected_tsvs(
    samples: Iterable[dict],
    cluster_lookup: AttributionLookup,
    annotation_df: pd.DataFrame,
    *,
    chrom_normalizer: Optional[Callable[[str], str]] = None,
    mode: str = "origin",
    chunk_size: int = 250_000,
    source: str = "body",
) -> pd.DataFrame:
    """Build cluster-gene attributions by read-body overlap.

    ``mode='origin'`` assigns each read to the transcription-upstream overlapped
    CDS/gene. ``mode='proportional'`` splits each read by CDS overlap length.
    """
    from intervaltree import IntervalTree

    if mode not in {"origin", "proportional"}:
        raise ValueError("mode must be 'origin' or 'proportional'")

    cds_trees = _build_body_feature_trees(annotation_df, chrom_normalizer)
    agg: defaultdict[tuple[str, str], float] = defaultdict(float)
    required = {"chrom", "strand", "alignment_start", "alignment_end"}

    for sample in samples:
        tsv_path = Path(sample["path"])
        header = pd.read_csv(tsv_path, sep="\t", nrows=0)
        pos_col = "corrected_3prime" if "corrected_3prime" in header.columns else "corrected_position"
        missing = required | {pos_col}
        missing = missing - set(header.columns)
        if missing:
            continue
        usecols = ["chrom", "strand", pos_col, "alignment_start", "alignment_end"]
        has_fraction = "fraction" in header.columns
        if has_fraction:
            usecols.append("fraction")
        for chunk in pd.read_csv(tsv_path, sep="\t", chunksize=chunk_size, usecols=usecols):
            chunk = chunk.dropna(subset=["alignment_start", "alignment_end", pos_col])
            if chunk.empty:
                continue
            if chrom_normalizer is not None:
                chunk["chrom"] = chunk["chrom"].map(chrom_normalizer)
            for row in chunk.itertuples(index=False):
                chrom = getattr(row, "chrom")
                strand = getattr(row, "strand")
                pos = int(getattr(row, pos_col))
                start = int(getattr(row, "alignment_start"))
                end = int(getattr(row, "alignment_end"))
                weight = float(getattr(row, "fraction", 1.0)) if has_fraction else 1.0
                cluster_id = cluster_lookup(chrom, strand, pos)
                if cluster_id is None:
                    continue
                overlaps = list(cds_trees.get((chrom, strand), IntervalTree())[start:end])
                if not overlaps:
                    continue
                if mode == "origin":
                    gene_id = _choose_origin_gene(overlaps, strand)
                    if gene_id:
                        agg[(cluster_id, gene_id)] += weight
                else:
                    gene_overlap: defaultdict[str, int] = defaultdict(int)
                    for iv in overlaps:
                        ov = max(0, min(end, iv.end) - max(start, iv.begin))
                        if ov > 0:
                            gene_overlap[iv.data["gene_id"]] += ov
                    total = sum(gene_overlap.values())
                    if total <= 0:
                        continue
                    for gene_id, ov in gene_overlap.items():
                        agg[(cluster_id, gene_id)] += weight * (ov / total)

    rows = [
        {
            "cluster_id": cluster_id,
            "gene_id": gene_id,
            "gene_name": gene_id,
            "raw_attributed_count": count,
            "source": source,
        }
        for (cluster_id, gene_id), count in agg.items()
        if count > 0
    ]
    return normalize_cluster_gene_attributions(pd.DataFrame(rows)) if rows else _empty_attribution_df()


def apply_annotation_fallback(
    primary: pd.DataFrame,
    annotation_fallback: pd.DataFrame,
    clusters_df: pd.DataFrame,
) -> pd.DataFrame:
    """Fill clusters missing primary attribution from annotation fallback."""
    if annotation_fallback.empty:
        return normalize_cluster_gene_attributions(primary)
    if primary.empty:
        return annotation_fallback.copy()
    covered = set(primary["cluster_id"].astype(str))
    fallback = annotation_fallback[
        ~annotation_fallback["cluster_id"].astype(str).isin(covered)
    ].copy()
    if fallback.empty:
        return normalize_cluster_gene_attributions(primary)
    return normalize_cluster_gene_attributions(pd.concat([primary, fallback], ignore_index=True))


def normalize_cluster_gene_attributions(attr_df: pd.DataFrame) -> pd.DataFrame:
    """Normalize attribution weights within each cluster."""
    if attr_df is None or attr_df.empty:
        return _empty_attribution_df()
    df = attr_df.copy()
    for col in ("gene_name", "source"):
        if col not in df.columns:
            df[col] = df["gene_id"] if col == "gene_name" else "unknown"
    if "raw_attributed_count" not in df.columns:
        if "attributed_count" in df.columns:
            df["raw_attributed_count"] = df["attributed_count"]
        else:
            df["raw_attributed_count"] = 1.0
    df = df[df["gene_id"].notna() & (df["gene_id"] != "intergenic")].copy()
    if df.empty:
        return _empty_attribution_df()
    grouped = (
        df.groupby(["cluster_id", "gene_id", "gene_name", "source"], sort=False)["raw_attributed_count"]
        .sum()
        .reset_index()
    )
    totals = grouped.groupby("cluster_id")["raw_attributed_count"].transform("sum")
    grouped = grouped[totals > 0].copy()
    grouped["attribution_weight"] = grouped["raw_attributed_count"] / totals[totals > 0]
    return grouped[
        ["cluster_id", "gene_id", "gene_name", "raw_attributed_count", "attribution_weight", "source"]
    ]


def apply_primary_gene_to_clusters(
    clusters_df: pd.DataFrame,
    attr_df: pd.DataFrame,
) -> pd.DataFrame:
    """Set cluster ``gene_id`` to the highest-weight attribution for compatibility."""
    clusters = clusters_df.copy()
    if attr_df is None or attr_df.empty:
        return clusters
    ranked = attr_df.sort_values(
        ["cluster_id", "attribution_weight", "raw_attributed_count"],
        ascending=[True, False, False],
    ).drop_duplicates("cluster_id")
    mapping = ranked.set_index("cluster_id")
    clusters["gene_id"] = clusters["cluster_id"].map(mapping["gene_id"]).where(
        clusters["cluster_id"].isin(mapping.index),
        clusters.get("gene_id"),
    )
    clusters["gene_name"] = clusters["cluster_id"].map(mapping["gene_name"]).where(
        clusters["cluster_id"].isin(mapping.index),
        clusters.get("gene_name"),
    )
    clusters["gene_attribution_source"] = clusters["cluster_id"].map(mapping["source"])
    clusters["gene_attribution_weight"] = clusters["cluster_id"].map(mapping["attribution_weight"])
    return clusters


def apply_gene_names_from_annotation(
    attr_df: pd.DataFrame,
    annotation_df: Optional[pd.DataFrame],
) -> pd.DataFrame:
    """Fill attribution gene_name values from loaded annotation when possible."""
    if attr_df is None or attr_df.empty or annotation_df is None:
        return attr_df
    if "gene_id" not in annotation_df.columns or "gene_name" not in annotation_df.columns:
        return attr_df

    name_map = (
        annotation_df[["gene_id", "gene_name"]]
        .dropna()
        .drop_duplicates("gene_id")
        .set_index("gene_id")["gene_name"]
        .to_dict()
    )
    if not name_map:
        return attr_df

    df = attr_df.copy()
    mapped = df["gene_id"].map(name_map)
    if "gene_name" not in df.columns:
        df["gene_name"] = mapped
    else:
        needs_name = (
            df["gene_name"].isna()
            | (df["gene_name"].astype(str) == "")
            | (df["gene_name"].astype(str) == df["gene_id"].astype(str))
        )
        df.loc[needs_name, "gene_name"] = mapped.loc[needs_name]
    df["gene_name"] = df["gene_name"].fillna(df["gene_id"])
    return df


def _build_body_feature_trees(annotation_df: pd.DataFrame, chrom_normalizer=None) -> dict:
    from intervaltree import IntervalTree

    feature_types = {"CDS", "gene", "mRNA", "ncRNA", "tRNA", "snoRNA", "snRNA"}
    trees: defaultdict[tuple[str, str], IntervalTree] = defaultdict(IntervalTree)
    for row in annotation_df.itertuples(index=False):
        feature_type = getattr(row, "feature_type", "gene")
        if feature_type not in feature_types:
            continue
        chrom = getattr(row, "chrom")
        if chrom_normalizer is not None:
            chrom = chrom_normalizer(chrom)
        strand = getattr(row, "strand", "+")
        start = int(getattr(row, "start"))
        end = int(getattr(row, "end"))
        if end <= start:
            continue
        gene_id = getattr(row, "gene_id", None)
        if _is_missing_gene(gene_id):
            continue
        gene_name = getattr(row, "gene_name", gene_id)
        trees[(chrom, strand)][start:end] = {
            "gene_id": str(gene_id),
            "gene_name": str(gene_name or gene_id),
            "start": start,
            "end": end,
        }
    return dict(trees)


def _choose_origin_gene(overlaps, strand: str) -> Optional[str]:
    gene_bounds = {}
    for iv in overlaps:
        gene_id = iv.data["gene_id"]
        if gene_id not in gene_bounds:
            gene_bounds[gene_id] = [iv.begin, iv.end]
        else:
            gene_bounds[gene_id][0] = min(gene_bounds[gene_id][0], iv.begin)
            gene_bounds[gene_id][1] = max(gene_bounds[gene_id][1], iv.end)
    if not gene_bounds:
        return None
    if strand == "+":
        return min(gene_bounds.items(), key=lambda kv: (kv[1][0], kv[1][1], kv[0]))[0]
    return max(gene_bounds.items(), key=lambda kv: (kv[1][1], kv[1][0], kv[0]))[0]


def _distance_to_interval(pos: int, start: int, end: int) -> int:
    if start <= pos <= end:
        return 0
    if pos < start:
        return start - pos
    return pos - end


def _is_missing_gene(gene_id) -> bool:
    return gene_id is None or pd.isna(gene_id) or str(gene_id) == "" or str(gene_id) == "None"


def _empty_attribution_df() -> pd.DataFrame:
    return pd.DataFrame(columns=[
        "cluster_id",
        "gene_id",
        "gene_name",
        "raw_attributed_count",
        "attribution_weight",
        "source",
    ])
