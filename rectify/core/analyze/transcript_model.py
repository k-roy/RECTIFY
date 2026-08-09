"""
Transcript-model-aware CPA-site classifier (planning/167 + 169).

RECTIFY historically attributed each CPA cluster to the *nearest gene 3' end*
within a fixed ``-500..+100`` window (``clustering.annotate_clusters_with_genes``).
That worldview is correct for canonical 3'UTR ends but structurally blind to
CPA sites that fire *upstream* of the canonical site — intronic pA, CDS-internal
pA (NSD), 5'UTR pA, and cryptic-3'UTR-intron pA — which all fall out as
``gene_id=None`` / intergenic.

This module reconstructs a per-gene transcript model from the loaded annotation
DataFrame (``loaders.load_annotation`` → flat one-row-per-feature table) plus an
optional pluggable set of supplementary ncRNA atlas tracks (CUT/SUT/XUT/…), and
exposes:

* :class:`TranscriptModel` — interval trees over gene bodies, per-feature spans
  (CDS / 5'UTR / 3'UTR / intron), an opposite-strand feature index (antisense),
  and an ncRNA index; annotated introns read directly from ``intron`` rows
  (with a CDS-gap fallback for annotations that lack them); cryptic introns
  registered from per-read splice junctions.
* :meth:`TranscriptModel.classify_site` — position → a ``region_class`` plus the
  continuous transcript coordinates and antisense/also-within context that make
  the premature-vs-canonical continuum *data* rather than a frozen boundary.
* :func:`containment_attributions_from_clusters` — a containment-first cluster→gene
  attribution table that plugs into the existing weighted-attribution machinery
  (``cluster_gene_attribution.apply_primary_gene_to_clusters``).
* :func:`annotate_clusters_with_transcript_model` — adds the per-cluster
  ``region_class`` + continuous-coordinate columns onto ``cpa_clusters.tsv``,
  running the cryptic-intron read pass when corrected TSVs are available.

Coordinate convention: everything is **0-based half-open** (the loader subtracts
1 from GFF ``start``). Per-read ``junctions`` (0-based, end-exclusive from CIGAR
``N``) live in the same space, so the annotated-intron set MUST be built from the
loader-converted annotation, not from raw GFF coordinates.

Design provenance: planning/167, 169, 170a–d (Chanfreau_Lab).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# region_class vocabulary (169 D4: antisense is NEVER a region_class)
# ─────────────────────────────────────────────────────────────────────────────
REGION_5UTR = "5primeUTR"
REGION_CDS_INTERNAL = "CDS_internal"
REGION_INTRONIC = "intronic"
REGION_CRYPTIC_INTRONIC = "cryptic_intronic"
REGION_3UTR_PROXIMAL = "3primeUTR_proximal"
REGION_3UTR_DISTAL = "3primeUTR_distal"
REGION_DOWNSTREAM = "downstream_readthrough"
REGION_INTERGENIC = "intergenic"
# ncRNA sense hits are emitted as ``ncRNA:<ncrna_class>`` (e.g. ``ncRNA:CUT``).

INTRON_EVIDENCE_ANNOTATED = "annotated"
INTRON_EVIDENCE_CRYPTIC = "cryptic-read-junction"
INTRON_EVIDENCE_NONE = "none"

# Feature types that denote a non-coding-RNA transcript. Protein-coding is
# clean (container=``gene``, transcript=``mRNA``); ncRNA container + transcript
# rows collapse to the same normalized type (170a §3/§5).
NCRNA_CLASSES = {"tRNA", "snoRNA", "snRNA", "rRNA", "ncRNA", "CUT", "SUT", "XUT"}

# Default PRIMARY-region precedence for sense-overlapping features (169 D5).
# Lower index wins. ``mRNA`` == protein-coding gene body.
DEFAULT_PRECEDENCE: Tuple[str, ...] = (
    "mRNA",
    "snoRNA",
    "snRNA",
    "CUT",
    "SUT",
    "XUT",
    "ncRNA",
    "tRNA",
    "rRNA",
)

# The full ordered column list the classifier emits (superset of 169).
CLASSIFIER_COLUMNS: Tuple[str, ...] = (
    "region_class",
    "attributed_gene",
    "attributed_feature_id",
    "attribution_source",
    "distance_to_stop_codon",
    "orf_fraction",
    "distance_to_canonical_CPA",
    "is_premature",
    "within_intron_id",
    "intron_rank",
    "pos_within_intron_frac",
    "intron_evidence",
    "spliced_support",
    "antisense_context",
    "also_within",
    "distance_to_ncrna_3prime",
    "pos_within_ncrna_frac",
)


def _opp(strand: str) -> str:
    return "-" if strand == "+" else "+"


def _exon_blocks_from(
    body_start: int, body_end: int, introns_asc: Sequence[Tuple[int, int]]
) -> Tuple[Tuple[int, int], ...]:
    """Transcript exon blocks (genomic ascending, 0-based half-open) = the body
    span minus the introns lying within it (F1).

    ``introns_asc`` is an ascending list of ``(start, end)`` intron spans. Introns
    are clamped to the body and coalesced so overlapping/degenerate rows cannot
    produce negative-length blocks; the result partitions the transcript into the
    exonic segments used for spliced-mRNA coordinate projection.
    """
    inner = sorted(
        (max(s, body_start), min(e, body_end))
        for (s, e) in introns_asc
        if e > body_start and s < body_end
    )
    blocks: List[Tuple[int, int]] = []
    cur = body_start
    for (s, e) in inner:
        if s > cur:
            blocks.append((cur, s))
        cur = max(cur, e)
    if cur < body_end:
        blocks.append((cur, body_end))
    return tuple(blocks)


def _is_missing(value) -> bool:
    return value is None or (isinstance(value, float) and pd.isna(value)) or str(value) in ("", "None", "nan")


# ─────────────────────────────────────────────────────────────────────────────
# Per-gene model
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class GeneModel:
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    body_start: int          # widest transcript span, 0-based
    body_end: int            # widest transcript span, half-open
    is_coding: bool
    ncrna_class: Optional[str]
    source: str
    cds_start: Optional[int] = None   # union CDS extent (coding only)
    cds_end: Optional[int] = None
    canonical_cpa: int = 0            # outermost (3'-most) transcript end, genomic 0-based
    mrna_ends: Tuple[int, ...] = ()   # every isoform 3' end (genomic 0-based)
    # introns as (start, end, rank) sorted 5'->3'
    introns: Tuple[Tuple[int, int, int], ...] = ()
    # transcript exon blocks (genomic ascending, 0-based half-open) = body minus
    # introns; the substrate for the SPLICED-mRNA projection (F1). Cached at build.
    exon_blocks: Tuple[Tuple[int, int], ...] = ()

    # ---- derived geometry (all genomic 0-based) ----
    @property
    def stop_boundary(self) -> Optional[int]:
        """Genomic coord of the last CDS base (3' end of the ORF)."""
        if not self.is_coding or self.cds_start is None or self.cds_end is None:
            return None
        return (self.cds_end - 1) if self.strand == "+" else self.cds_start

    def utr5_span(self) -> Optional[Tuple[int, int]]:
        if not self.is_coding or self.cds_start is None:
            return None
        if self.strand == "+":
            return (self.body_start, self.cds_start) if self.cds_start > self.body_start else None
        return (self.cds_end, self.body_end) if self.body_end > self.cds_end else None

    def utr3_span(self) -> Optional[Tuple[int, int]]:
        if not self.is_coding or self.cds_end is None:
            return None
        if self.strand == "+":
            return (self.cds_end, self.body_end) if self.body_end > self.cds_end else None
        return (self.body_start, self.cds_start) if self.cds_start > self.body_start else None


# ─────────────────────────────────────────────────────────────────────────────
# ncRNA atlas track (169 D1/D2)
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class NcrnaTrack:
    """A supplementary ncRNA annotation source with provenance tags."""
    df: pd.DataFrame          # annotation-style rows (chrom,start,end,strand,gene_id,gene_name,feature_type,...)
    source: str               # e.g. Xu2009, vanDijk2011
    ncrna_class: str          # CUT / SUT / XUT / ...
    genome_build: str = "R64"


def _assert_build_compatible(atlas_build: str, genome_build: Optional[str]) -> None:
    """Assert an atlas is coordinate-compatible with the working genome (169 D3).

    Tolerant of the SGD ``#!assembly R64-4-1`` label lag: assert on the shared
    ``R64`` / ``sacCer3`` nuclear frame, NOT an exact ``-5-1`` substring (170b §5e),
    since a hard exact-string check would false-fail on every R64 point release.
    """
    if not genome_build:
        return
    ab = str(atlas_build).lower()
    gb = str(genome_build).lower()
    ab_frame = "r64" in ab or "saccer3" in ab
    gb_frame = "r64" in gb or "saccer3" in gb
    if ab_frame and gb_frame:
        return
    if ab == gb:
        return
    raise ValueError(
        f"ncRNA atlas build {atlas_build!r} is not coordinate-compatible with genome "
        f"build {genome_build!r}; re-lift the atlas and update the registry (169 D3)."
    )


def load_ncrna_track_from_gff(
    path: str,
    *,
    source: str,
    ncrna_class: str,
    genome_build: str = "R64",
    chrom_format: str = "passthrough",
    normalize_chroms: bool = True,
) -> NcrnaTrack:
    """Load one supplementary ncRNA track from a GFF/BED file.

    The file's own feature-type column is IGNORED for classification purposes —
    every feature is force-tagged ``ncrna_class`` and ``source`` (169 D1 ad-hoc
    ``file:source:class`` form), so a single-class GFF or a headerless BED both
    work. Non-feature lines (``NF<8``, stray headers like the SUT file's leading
    'CUTs' artifact, 170b §2c) are skipped by the underlying loader.
    """
    from .loaders import load_annotation

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"ncRNA atlas track not found: {path}")
    is_bed = p.name.endswith((".bed", ".bed.gz"))
    if is_bed:
        df = _load_bed_as_annotation(str(p))
        if normalize_chroms:
            from ...utils.chromosome import normalize_dataframe_chromosomes
            df = normalize_dataframe_chromosomes(df, "chrom", chrom_format)
    else:
        df = load_annotation(str(p), normalize_chroms=normalize_chroms, chrom_format=chrom_format)
    # force-tag class regardless of the file's own col-3 type
    df = df.copy()
    df["feature_type"] = ncrna_class
    return NcrnaTrack(df=df, source=source, ncrna_class=ncrna_class, genome_build=genome_build)


def _load_bed_as_annotation(path: str) -> pd.DataFrame:
    """Read a BED6 into the annotation-style schema (0-based half-open kept)."""
    import gzip as _gz
    _open = _gz.open if str(path).endswith(".gz") else open
    rows = []
    with _open(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                start, end = int(f[1]), int(f[2])   # BED is already 0-based half-open
            except ValueError:
                continue
            name = f[3] if len(f) > 3 and f[3] not in (".", "") else f"{f[0]}:{start}-{end}"
            strand = f[5] if len(f) > 5 and f[5] in ("+", "-") else "+"
            rows.append({
                "chrom": f[0], "start": start, "end": end, "strand": strand,
                "gene_id": name, "gene_name": name, "feature_type": "ncRNA",
                "raw_feature": "ncRNA", "source": "bed", "transcript_id": None, "parent_id": None,
            })
    return pd.DataFrame(rows, columns=[
        "chrom", "start", "end", "strand", "gene_id", "gene_name",
        "feature_type", "raw_feature", "source", "transcript_id", "parent_id",
    ])


def parse_ncrna_annotation_spec(spec: str) -> Tuple[str, str, str]:
    """Parse an ad-hoc ``path:source:class`` spec (169 D1).

    Splits from the RIGHT so a Windows drive / absolute path with colons in the
    path itself still yields ``(path, source, class)``.
    """
    parts = spec.rsplit(":", 2)
    if len(parts) != 3 or not all(parts):
        raise ValueError(
            f"--ncrna-annotations spec must be 'file.gff:source:class', got {spec!r}"
        )
    return parts[0], parts[1], parts[2]


def _default_registry_path() -> Path:
    # rectify/core/analyze/transcript_model.py -> rectify/data/ncrna_atlases/atlases.yaml
    return Path(__file__).resolve().parents[2] / "data" / "ncrna_atlases" / "atlases.yaml"


def load_ncrna_atlas(
    name: Optional[str] = None,
    *,
    registry_path: Optional[str] = None,
    genome_build: Optional[str] = None,
    chrom_format: str = "passthrough",
) -> List[NcrnaTrack]:
    """Load a named atlas from the registry (169 D2).

    ``name`` is a registry key (or None → the registry ``default``). Each track's
    build is asserted compatible with ``genome_build`` (tolerant of the R64 label
    lag). Track file paths in the registry are resolved relative to the registry
    file's directory. Returns a list of :class:`NcrnaTrack`.
    """
    import yaml

    reg_path = Path(registry_path) if registry_path else _default_registry_path()
    if not reg_path.exists():
        raise FileNotFoundError(f"ncRNA atlas registry not found: {reg_path}")
    with open(reg_path) as fh:
        registry = yaml.safe_load(fh) or {}
    atlases = registry.get("atlases", {})
    key = name or registry.get("default")
    if key not in atlases:
        raise ValueError(
            f"ncRNA atlas {key!r} not in registry {reg_path} (available: {sorted(atlases)})"
        )
    entry = atlases[key]
    atlas_build = entry.get("genome_build", "R64")
    _assert_build_compatible(atlas_build, genome_build)
    tracks: List[NcrnaTrack] = []
    for t in entry.get("tracks", []):
        track_path = (reg_path.parent / t["file"]).resolve()
        tracks.append(load_ncrna_track_from_gff(
            str(track_path),
            source=t.get("source", key),
            ncrna_class=t.get("ncrna_class", "ncRNA"),
            genome_build=atlas_build,
            chrom_format=chrom_format,
        ))
    return tracks


def build_ncrna_tracks_from_args(
    ncrna_atlas: Optional[str],
    ncrna_annotations: Optional[Sequence[str]],
    *,
    genome_build: Optional[str] = None,
    chrom_format: str = "passthrough",
) -> List[NcrnaTrack]:
    """Assemble ncRNA tracks from the two CLI forms (registry key + ad-hoc specs)."""
    tracks: List[NcrnaTrack] = []
    if ncrna_atlas:
        tracks.extend(load_ncrna_atlas(
            ncrna_atlas, genome_build=genome_build, chrom_format=chrom_format,
        ))
    for spec in (ncrna_annotations or []):
        path, source, cls = parse_ncrna_annotation_spec(spec)
        tracks.append(load_ncrna_track_from_gff(
            path, source=source, ncrna_class=cls,
            genome_build=genome_build or "R64", chrom_format=chrom_format,
        ))
    return tracks


def build_transcript_model_for_analyze(args, annotation_df, chrom_format="passthrough"):
    """Build a :class:`TranscriptModel` from analyze ``args`` + the loaded annotation.

    Reads ``--ncrna-atlas`` / ``--ncrna-annotations`` / ``--gene-attribution-window`` /
    ``--utr3-proximal-distal-split`` off ``args``. Returns ``None`` (and prints a
    warning) if the annotation is absent or the model cannot be built — the caller
    then falls back to legacy attribution and skips region classification, so the
    default-flip never turns a working run into a hard failure.
    """
    if annotation_df is None:
        return None
    from ...utils.chromosome import normalize_chromosome
    normalizer = lambda x: normalize_chromosome(x, chrom_format)
    try:
        genome_build = getattr(args, "assembly", None) or getattr(args, "genome_build", None) or "R64"
        tracks = build_ncrna_tracks_from_args(
            getattr(args, "ncrna_atlas", None),
            getattr(args, "ncrna_annotations", None),
            genome_build=genome_build,
            chrom_format=chrom_format,
        )
        model = TranscriptModel(
            annotation_df,
            ncrna_tracks=tracks,
            utr3_proximal_distal_split=int(getattr(args, "utr3_proximal_distal_split", 150) or 150),
            readthrough_window=int(getattr(args, "gene_attribution_window", 100) or 100),
            junction_match_tol=int(getattr(args, "junction_match_tol", 8) or 8),
            min_gap_intron_len=int(getattr(args, "min_gap_intron_len", 20) or 20),
            chrom_normalizer=normalizer,
        )
        for w in model.warnings:
            print(f"  transcript-model: {w}")
        return model
    except Exception as exc:
        print(f"  WARNING: could not build transcript model ({exc}); "
              f"falling back to legacy gene attribution, region_class not emitted")
        return None


class TranscriptModel:
    """A transcript-relative feature model built once from the annotation.

    Parameters
    ----------
    annotation_df:
        The flat feature table from ``loaders.load_annotation``. Must be already
        chromosome-normalized to the same convention the clusters use.
    ncrna_tracks:
        Optional supplementary ncRNA sources (169). Each is tagged with a
        ``source`` and ``ncrna_class``; features are keyed on the *sense* strand.
    precedence:
        Ordered region-class precedence for sense-overlapping features (169 D5).
    utr3_proximal_distal_split:
        3'UTR proximal/distal boundary in nt past the stop codon (default 150).
    readthrough_window:
        Downstream window (nt) past a gene 3' end within which an uncontained
        site is ``downstream_readthrough`` rather than ``intergenic``. Kept equal
        to the attribution proximity window so ``gene_id`` and ``region_class``
        stay coherent for readthrough clusters.
    """

    def __init__(
        self,
        annotation_df: pd.DataFrame,
        *,
        ncrna_tracks: Optional[Sequence[NcrnaTrack]] = None,
        precedence: Sequence[str] = DEFAULT_PRECEDENCE,
        utr3_proximal_distal_split: int = 150,
        readthrough_window: int = 100,
        junction_match_tol: int = 8,
        min_gap_intron_len: int = 20,
        chrom_normalizer: Optional[Callable[[str], str]] = None,
    ) -> None:
        self.precedence = {name: i for i, name in enumerate(precedence)}
        self.utr3_split = int(utr3_proximal_distal_split)
        self.readthrough_window = int(readthrough_window)
        # F2: ± bp tolerance for treating a read junction as an ANNOTATED intron
        # (ONT/DRS splice-site imprecision). ``0`` restores exact tuple matching.
        self.junction_tol = int(junction_match_tol)
        # F7: minimum CDS-gap length to derive an annotated intron. Gaps below this
        # are +1 frameshift sites (Ty ORFs), not spliceosomal introns.
        self.min_gap_intron_len = int(min_gap_intron_len)
        self._norm = chrom_normalizer
        self.warnings: List[str] = []

        self.genes: Dict[str, GeneModel] = {}
        # interval trees keyed by (chrom, strand)
        self._gene_body_trees: Dict[Tuple[str, str], "object"] = {}
        self._ncrna_trees: Dict[Tuple[str, str], "object"] = {}
        self._feature_trees: Dict[Tuple[str, str], "object"] = {}   # typed spans, both strands
        self._annotated_intron_trees: Dict[Tuple[str, str], "object"] = {}
        self._cryptic_intron_trees: Dict[Tuple[str, str], "object"] = {}
        # annotated intron (start,end) set per (chrom,strand) for cryptic detection
        self._annotated_intron_set: Dict[Tuple[str, str], set] = defaultdict(set)
        # sorted-by-donor index of annotated introns for fuzzy junction match (F2)
        self._annotated_introns_sorted: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
        self._annotated_donors: Dict[Tuple[str, str], List[int]] = {}

        self._build(annotation_df, ncrna_tracks or [])

    # ── construction ─────────────────────────────────────────────────────────
    def _c(self, chrom: str) -> str:
        return self._norm(chrom) if self._norm is not None else chrom

    def _build(self, annotation_df: pd.DataFrame, ncrna_tracks: Sequence[NcrnaTrack]) -> None:
        from intervaltree import IntervalTree

        if annotation_df is None or annotation_df.empty:
            self.warnings.append("annotation is empty; classifier will report intergenic for all sites")
            return

        df = annotation_df
        cols = set(df.columns)
        for required in ("chrom", "start", "end", "strand", "feature_type"):
            if required not in cols:
                raise ValueError(f"annotation_df missing required column '{required}'")

        # transcript token -> systematic gene (170a §2). Sub-feature gene_id is
        # the transcript token; only this join is reliable.
        tx_types = {"mRNA"} | NCRNA_CLASSES
        tx = df[df["feature_type"].isin(tx_types)]
        if "parent_id" in cols:
            tx = tx[tx["parent_id"].notna()]
        tok2gene: Dict[str, str] = {}
        tok2name: Dict[str, str] = {}
        for row in tx.itertuples(index=False):
            token = getattr(row, "gene_name", None)
            gid = getattr(row, "gene_id", None)
            if not _is_missing(token) and not _is_missing(gid):
                tok2gene.setdefault(str(token), str(gid))

        def resolve_gene(row) -> Optional[str]:
            """Map a sub-feature row to its systematic gene via the token join."""
            gid = getattr(row, "gene_id", None)
            if _is_missing(gid):
                return None
            gid = str(gid)
            return tok2gene.get(gid, gid)

        # F7: whether a gene uses CDS-gap intron derivation is decided PER-GENE
        # (does *this* gene have explicit 'intron' rows?), not by a single global
        # flag — one explicit intron row elsewhere must not silently disable the
        # gap fallback for every other gene.
        genes_with_explicit_introns: set = set()

        # ---- gather per-gene geometry ----
        # gene container rows carry the clean systematic id + common name (170a §5)
        gene_meta: Dict[str, Tuple[str, str]] = {}   # gene_id -> (gene_name, source)
        raw_has = "raw_feature" in cols
        parent_has = "parent_id" in cols
        for row in df.itertuples(index=False):
            ft = getattr(row, "feature_type", None)
            raw = str(getattr(row, "raw_feature", "")) if raw_has else ""
            parent = getattr(row, "parent_id", None) if parent_has else None
            is_container = (raw.endswith("gene")) or (ft == "gene")
            if is_container and _is_missing(parent):
                gid = getattr(row, "gene_id", None)
                if not _is_missing(gid):
                    gene_meta[str(gid)] = (
                        str(getattr(row, "gene_name", gid) or gid),
                        str(getattr(row, "source", "SGD") or "SGD"),
                    )

        mrna_spans: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        cds_extent: Dict[str, List[int]] = {}          # gene -> [min_start, max_end]
        gene_strand: Dict[str, str] = {}
        gene_chrom: Dict[str, str] = {}
        gene_is_coding: Dict[str, bool] = {}
        ncrna_span: Dict[str, List[int]] = {}          # ncRNA gene -> [min_start,max_end]
        ncrna_class_of: Dict[str, str] = {}
        intron_rows_by_gene: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        cds_blocks_by_tx: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

        for row in df.itertuples(index=False):
            ft = getattr(row, "feature_type", None)
            chrom = self._c(getattr(row, "chrom"))
            strand = getattr(row, "strand", "+")
            try:
                start = int(getattr(row, "start"))
                end = int(getattr(row, "end"))
            except (TypeError, ValueError):
                continue
            if end <= start:
                continue

            if ft == "mRNA":
                gid = getattr(row, "gene_id", None)   # mRNA gene_id == systematic (170a §2)
                if _is_missing(gid):
                    continue
                gid = str(gid)
                mrna_spans[gid].append((start, end))
                gene_strand[gid] = strand
                gene_chrom[gid] = chrom
                gene_is_coding[gid] = True
            elif ft == "CDS":
                gene = resolve_gene(row)
                if gene is None:
                    continue
                gene_is_coding[gene] = True
                gene_strand.setdefault(gene, strand)
                gene_chrom.setdefault(gene, chrom)
                if gene not in cds_extent:
                    cds_extent[gene] = [start, end]
                else:
                    cds_extent[gene][0] = min(cds_extent[gene][0], start)
                    cds_extent[gene][1] = max(cds_extent[gene][1], end)
                token = getattr(row, "parent_id", None) if parent_has else None
                if not _is_missing(token):
                    cds_blocks_by_tx[str(token)].append((start, end))
            elif ft == "intron":
                gene = resolve_gene(row)
                self._annotated_intron_set[(chrom, strand)].add((start, end))
                if gene is not None:
                    intron_rows_by_gene[gene].append((start, end))
                    genes_with_explicit_introns.add(gene)
            elif ft in NCRNA_CLASSES:
                # ncRNA container or transcript row; take widest span per gene id
                gene = resolve_gene(row)
                if gene is None:
                    gid = getattr(row, "gene_id", None)
                    gene = str(gid) if not _is_missing(gid) else None
                if gene is None:
                    continue
                gene_strand.setdefault(gene, strand)
                gene_chrom.setdefault(gene, chrom)
                gene_is_coding.setdefault(gene, False)
                ncrna_class_of.setdefault(gene, ft)
                if gene not in ncrna_span:
                    ncrna_span[gene] = [start, end]
                else:
                    ncrna_span[gene][0] = min(ncrna_span[gene][0], start)
                    ncrna_span[gene][1] = max(ncrna_span[gene][1], end)

        # supplementary ncRNA atlas tracks (169): tag source + class explicitly
        atlas_source: Dict[str, str] = {}
        for track in ncrna_tracks:
            for row in track.df.itertuples(index=False):
                chrom = self._c(getattr(row, "chrom"))
                strand = getattr(row, "strand", "+")
                try:
                    start = int(getattr(row, "start"))
                    end = int(getattr(row, "end"))
                except (TypeError, ValueError):
                    continue
                if end <= start:
                    continue
                gid = getattr(row, "gene_id", None)
                name = getattr(row, "gene_name", None)
                feature_id = None
                for cand in (name, gid):
                    if not _is_missing(cand):
                        feature_id = str(cand)
                        break
                if feature_id is None:
                    feature_id = f"{track.ncrna_class}:{chrom}:{start}-{end}"
                gene_strand[feature_id] = strand
                gene_chrom[feature_id] = chrom
                gene_is_coding[feature_id] = False
                ncrna_class_of[feature_id] = track.ncrna_class
                atlas_source[feature_id] = track.source
                ncrna_span[feature_id] = [start, end]

        # ---- derive introns via CDS-block gaps, PER-GENE, for genes lacking
        # explicit 'intron' rows (170a §6b + F7). A gap shorter than
        # ``min_gap_intron_len`` is a +1 frameshift site (Ty ORFs, e.g. YBL005W-B
        # has a 1 bp gap), NOT a spliceosomal intron, so it is skipped. ----
        n_gap_derived = 0
        for token, blocks in cds_blocks_by_tx.items():
            gene = tok2gene.get(token, token)
            if gene in genes_with_explicit_introns:
                continue
            blocks = sorted(blocks)
            strand = gene_strand.get(gene, "+")
            chrom = gene_chrom.get(gene)
            for (s0, e0), (s1, e1) in zip(blocks, blocks[1:]):
                if s1 - e0 >= self.min_gap_intron_len:
                    intron_rows_by_gene[gene].append((e0, s1))
                    n_gap_derived += 1
                    if chrom is not None:
                        self._annotated_intron_set[(chrom, strand)].add((e0, s1))
        if n_gap_derived:
            self.warnings.append(
                f"derived {n_gap_derived} annotated intron(s) from CDS-block gaps for "
                "gene(s) lacking explicit 'intron' rows (per-gene fallback)"
            )

        # ---- assemble GeneModel objects ----
        for gene in set(list(mrna_spans) + list(ncrna_span)):
            strand = gene_strand.get(gene, "+")
            chrom = gene_chrom.get(gene)
            if chrom is None:
                continue
            is_coding = gene_is_coding.get(gene, gene in mrna_spans)
            name, source = gene_meta.get(gene, (gene, atlas_source.get(gene, "SGD")))
            if gene in atlas_source:
                source = atlas_source[gene]

            if is_coding and gene in mrna_spans:
                spans = mrna_spans[gene]
                body_start = min(s for s, _ in spans)
                body_end = max(e for _, e in spans)
                ends = tuple(e - 1 if strand == "+" else s for s, e in spans)
                canonical = max(ends) if strand == "+" else min(ends)
                cds = cds_extent.get(gene)
                cds_s = cds[0] if cds else None
                cds_e = cds[1] if cds else None
            else:
                span = ncrna_span.get(gene)
                if span is None:
                    continue
                body_start, body_end = span[0], span[1]
                canonical = (body_end - 1) if strand == "+" else body_start
                ends = (canonical,)
                cds_s = cds_e = None
                is_coding = False

            # introns sorted 5'->3' (rank), plus ascending exon blocks for the
            # spliced-mRNA projection (F1)
            introns_asc = sorted(set(intron_rows_by_gene.get(gene, [])))
            exon_blocks = _exon_blocks_from(body_start, body_end, introns_asc)
            ordered = list(reversed(introns_asc)) if strand == "-" else introns_asc
            introns = tuple((s, e, i + 1) for i, (s, e) in enumerate(ordered))

            self.genes[gene] = GeneModel(
                gene_id=gene,
                gene_name=name if not _is_missing(name) else gene,
                chrom=chrom,
                strand=strand,
                body_start=body_start,
                body_end=body_end,
                is_coding=is_coding,
                ncrna_class=ncrna_class_of.get(gene),
                source=source if not _is_missing(source) else "SGD",
                cds_start=cds_s,
                cds_end=cds_e,
                canonical_cpa=canonical,
                mrna_ends=ends,
                introns=introns,
                exon_blocks=exon_blocks,
            )

        self._build_trees(IntervalTree)
        self._finalize_annotated_index()
        self._sanity_warn(df)

    def _finalize_annotated_index(self) -> None:
        """Build the sorted-by-donor annotated-intron index for fuzzy junction
        matching (F2). Called once after the annotated-intron set is complete."""
        for key, spans in self._annotated_intron_set.items():
            arr = sorted(spans)
            self._annotated_introns_sorted[key] = arr
            self._annotated_donors[key] = [d for (d, _a) in arr]

    def _build_trees(self, IntervalTree) -> None:
        gene_body: Dict[Tuple[str, str], object] = defaultdict(IntervalTree)
        ncrna: Dict[Tuple[str, str], object] = defaultdict(IntervalTree)
        feat: Dict[Tuple[str, str], object] = defaultdict(IntervalTree)
        introns: Dict[Tuple[str, str], object] = defaultdict(IntervalTree)

        for gene, gm in self.genes.items():
            key = (gm.chrom, gm.strand)
            if gm.body_end <= gm.body_start:
                continue
            if gm.is_coding:
                gene_body[key][gm.body_start:gm.body_end] = gene
                # typed feature spans (for antisense context + also_within)
                if gm.cds_start is not None and gm.cds_end is not None and gm.cds_end > gm.cds_start:
                    feat[key][gm.cds_start:gm.cds_end] = {"kind": "CDS", "gene_id": gene, "gene_name": gm.gene_name}
                u5 = gm.utr5_span()
                if u5:
                    feat[key][u5[0]:u5[1]] = {"kind": "5primeUTR", "gene_id": gene, "gene_name": gm.gene_name}
                u3 = gm.utr3_span()
                if u3:
                    feat[key][u3[0]:u3[1]] = {"kind": "3primeUTR", "gene_id": gene, "gene_name": gm.gene_name}
            else:
                ncrna[key][gm.body_start:gm.body_end] = gene
                feat[key][gm.body_start:gm.body_end] = {
                    "kind": gm.ncrna_class or "ncRNA", "gene_id": gene, "gene_name": gm.gene_name,
                }
            for (s, e, rank) in gm.introns:
                if e > s:
                    introns[key][s:e] = {"gene_id": gene, "rank": rank, "start": s, "end": e}

        self._gene_body_trees = dict(gene_body)
        self._ncrna_trees = dict(ncrna)
        self._feature_trees = dict(feat)
        self._annotated_intron_trees = dict(introns)

    def _sanity_warn(self, df: pd.DataFrame) -> None:
        n_coding = sum(1 for g in self.genes.values() if g.is_coding)
        if n_coding == 0:
            self.warnings.append(
                "no protein-coding genes with mRNA+CDS found; ORF/UTR classes unavailable "
                "(is this an exon-only GTF?)"
            )
        n_utr3 = sum(1 for g in self.genes.values() if g.is_coding and g.utr3_span())
        if n_coding and n_utr3 == 0:
            self.warnings.append(
                "no derivable 3'UTRs (mRNA co-extensive with CDS?); canonical CPA defaults to the stop codon"
            )

    # ── cryptic-intron registration (167 Part B.2) ───────────────────────────
    def is_annotated_junction(
        self, chrom: str, strand: str, start: int, end: int, *, tol: Optional[int] = None
    ) -> bool:
        """True if ``(start, end)`` corresponds to an annotated intron.

        Fuzzy by default (F2): a read junction whose donor AND acceptor are each
        within ``tol`` bp of an annotated intron is treated as that annotated
        intron — ONT/DRS splice-site imprecision on a REAL intron must not
        register a false cryptic intron. ``tol`` defaults to ``self.junction_tol``;
        pass ``tol=0`` to force an exact tuple match.
        """
        import bisect

        chrom = self._c(chrom)
        start, end = int(start), int(end)
        key = (chrom, strand)
        if (start, end) in self._annotated_intron_set.get(key, ()):   # fast exact path
            return True
        t = self.junction_tol if tol is None else int(tol)
        if t <= 0:
            return False
        donors = self._annotated_donors.get(key)
        if not donors:
            return False
        arr = self._annotated_introns_sorted[key]
        lo = bisect.bisect_left(donors, start - t)
        hi = bisect.bisect_right(donors, start + t)
        for i in range(lo, hi):
            _d, a = arr[i]
            if abs(a - end) <= t:
                return True
        return False

    def _within_transcribed_body(self, chrom: str, strand: str, start: int, end: int) -> bool:
        """True if a (start,end) span overlaps any sense-strand gene/ncRNA body."""
        for trees in (self._gene_body_trees, self._ncrna_trees):
            tree = trees.get((chrom, strand))
            if tree and tree.overlap(start, end):
                return True
        return False

    def register_cryptic_introns(
        self,
        junctions: Iterable[Tuple[str, str, int, int]],
        *,
        min_support: int = 1,
        require_gene_body: bool = True,
    ) -> int:
        """Register read-derived splice junctions that are NOT annotated.

        ``junctions`` is an iterable of ``(chrom, strand, start, end)`` (0-based
        half-open). A junction is registered as a cryptic intron only if it (a)
        matches no annotated intron, (b) is supported by ``>= min_support`` reads,
        and (c) (when ``require_gene_body``) lies within a transcribed body on its
        strand. (b)+(c) guard against DRS junction noise — single misaligned
        N-CIGAR reads registering spurious genome-wide cryptic introns. Returns
        the number of distinct cryptic introns registered.
        """
        from intervaltree import IntervalTree

        support: Dict[Tuple[str, str, int, int], int] = defaultdict(int)
        for chrom, strand, start, end in junctions:
            chrom = self._c(chrom)
            start, end = int(start), int(end)
            if end <= start:
                continue
            # F2: fuzzy — an imprecise placement of a REAL annotated intron (within
            # self.junction_tol on both donor and acceptor) is NOT a cryptic intron.
            if self.is_annotated_junction(chrom, strand, start, end):
                continue
            support[(chrom, strand, start, end)] += 1

        trees: Dict[Tuple[str, str], object] = {
            k: (self._cryptic_intron_trees.get(k) or IntervalTree())
            for k in {(c, s) for (c, s, _, _) in support}
        }
        n = 0
        for (chrom, strand, start, end), count in support.items():
            if count < min_support:
                continue
            if require_gene_body and not self._within_transcribed_body(chrom, strand, start, end):
                continue
            key = (chrom, strand)
            trees.setdefault(key, IntervalTree())[start:end] = {
                "start": start, "end": end, "support": count,
            }
            n += 1
        self._cryptic_intron_trees.update(trees)
        return n

    def has_cryptic_introns(self) -> bool:
        return any(len(t) for t in self._cryptic_intron_trees.values())

    # ── primary sense feature (single source of truth) ───────────────────────
    def primary_sense_feature(self, chrom: str, strand: str, pos: int):
        """Return the highest-precedence sense-strand feature containing ``pos``.

        Returns a dict ``{gene_id, gene_name, kind, source, is_coding}`` or None.
        Used by BOTH :meth:`classify_site` and
        :func:`containment_attributions_from_clusters` so ``attributed_gene`` and
        the attribution ``gene_id`` never diverge for contained clusters.
        """
        chrom = self._c(chrom)
        pos = int(pos)
        hits = []
        for gene in self._sense_hits(self._gene_body_trees, chrom, strand, pos):
            gm = self.genes.get(gene)
            if gm is not None:
                hits.append((self._prec("mRNA"), gm))
        for gene in self._sense_hits(self._ncrna_trees, chrom, strand, pos):
            gm = self.genes.get(gene)
            if gm is not None:
                hits.append((self._prec(gm.ncrna_class or "ncRNA"), gm))
        if not hits:
            return None
        best_prec = min(p for p, _ in hits)
        tied = [gm for p, gm in hits if p == best_prec]
        gm = self._choose_upstream(tied, strand)
        kind = "mRNA" if gm.is_coding else (gm.ncrna_class or "ncRNA")
        return {
            "gene_id": gm.gene_id,
            "gene_name": gm.gene_name,
            "kind": kind,
            "source": gm.source,
            "is_coding": gm.is_coding,
        }

    def sense_overlapping_genes(self, chrom: str, strand: str, pos: int) -> List[Dict[str, object]]:
        """All sense-strand gene / ncRNA features whose body contains ``pos`` (167
        Part A), the PRIMARY first (identical to :meth:`primary_sense_feature`)
        followed by the rest ordered by precedence then transcription-upstream.

        Used by :func:`containment_attributions_from_clusters` so multi-gene
        (nested / overlapping) loci emit ALL overlapped genes with weights rather
        than silently collapsing to the primary.
        """
        primary = self.primary_sense_feature(chrom, strand, pos)
        if primary is None:
            return []
        chrom_n = self._c(chrom)
        pos = int(pos)
        best: Dict[str, Tuple[int, GeneModel]] = {}
        for gene in self._sense_hits(self._gene_body_trees, chrom_n, strand, pos):
            gm = self.genes.get(gene)
            if gm is not None:
                best[gm.gene_id] = (self._prec("mRNA"), gm)
        for gene in self._sense_hits(self._ncrna_trees, chrom_n, strand, pos):
            gm = self.genes.get(gene)
            if gm is not None and gm.gene_id not in best:
                best[gm.gene_id] = (self._prec(gm.ncrna_class or "ncRNA"), gm)

        def _key(pgm: Tuple[int, GeneModel]):
            p, gm = pgm
            return (p, gm.body_start, gm.body_end, gm.gene_id) if strand == "+" \
                else (p, -gm.body_end, -gm.body_start, gm.gene_id)

        others = [
            gm for (_p, gm) in sorted(
                (v for k, v in best.items() if k != primary["gene_id"]), key=_key
            )
        ]
        return [primary] + [{
            "gene_id": gm.gene_id,
            "gene_name": gm.gene_name,
            "kind": "mRNA" if gm.is_coding else (gm.ncrna_class or "ncRNA"),
            "source": gm.source,
            "is_coding": gm.is_coding,
        } for gm in others]

    def _prec(self, name: str) -> int:
        return self.precedence.get(name, len(self.precedence) + 1)

    @staticmethod
    def _sense_hits(trees, chrom, strand, pos) -> List[str]:
        tree = trees.get((chrom, strand))
        if not tree:
            return []
        return [iv.data for iv in tree[pos]]

    @staticmethod
    def _choose_upstream(models: List[GeneModel], strand: str) -> GeneModel:
        """Transcription-upstream tiebreak (mirrors _choose_origin_gene)."""
        if strand == "+":
            return min(models, key=lambda g: (g.body_start, g.body_end, g.gene_id))
        return max(models, key=lambda g: (g.body_end, g.body_start, g.gene_id))

    # ── the classifier ───────────────────────────────────────────────────────
    def classify_site(
        self,
        chrom: str,
        strand: str,
        pos: int,
        *,
        junctions: Optional[Sequence[Tuple[int, int]]] = None,
        spliced_support: bool = False,
        fallback_gene: Optional[str] = None,
    ) -> Dict[str, object]:
        """Classify a single 3'-end position on the transcript-relative axis.

        Parameters
        ----------
        chrom, strand, pos:
            The corrected 3'-end (0-based). ``strand`` is the gene/RNA strand.
        junctions:
            This site's own supporting-read splice junctions (0-based half-open).
            Used to set ``spliced_support`` and to opportunistically register a
            cryptic intron. Population-level cryptic introns should be registered
            via :meth:`register_cryptic_introns` before classification.
        fallback_gene:
            The mode-selected ``gene_id`` already on the cluster (e.g. from the
            proximity-window attribution). When the site is not body-contained,
            it lets ``downstream_readthrough`` align with the attribution window
            so ``gene_id`` and ``region_class`` stay coherent.
        """
        chrom = self._c(chrom)
        pos = int(pos)
        out = {c: None for c in CLASSIFIER_COLUMNS}
        out["intron_evidence"] = INTRON_EVIDENCE_NONE
        out["spliced_support"] = bool(spliced_support)

        # cryptic-intron evidence from this site's own reads (opportunistic)
        if junctions:
            for (js, je) in junctions:
                if not self.is_annotated_junction(chrom, strand, js, je):
                    out["spliced_support"] = True

        out["antisense_context"] = self._antisense_context(chrom, strand, pos)

        primary = self.primary_sense_feature(chrom, strand, pos)
        also: List[str] = []

        if primary is not None and primary["is_coding"]:
            gm = self.genes[primary["gene_id"]]
            # F3: a pA may sit in an INNER/overlapping sense gene's ANNOTATED intron
            # while the primary (widest/upstream) gene is exonic there. Attribute the
            # intronic call to the intron's owner so within_intron_id / intron_rank
            # are not lost (nested/overlapping genes; the previous gate demanded the
            # intron belong to the primary and dropped inner-gene introns).
            ann_owner = self._annotated_intron_owner(chrom, strand, pos)
            if ann_owner is not None:
                gm = ann_owner
            self._classify_within_coding_gene(gm, pos, out)
            out["attributed_gene"] = gm.gene_id
            out["attributed_feature_id"] = gm.gene_id
            out["attribution_source"] = gm.source
            also.extend(self._other_sense_overlaps(chrom, strand, pos, exclude=gm.gene_id))
        elif primary is not None:
            gm = self.genes[primary["gene_id"]]
            out["region_class"] = f"ncRNA:{gm.ncrna_class or 'ncRNA'}"
            out["attributed_gene"] = gm.gene_id
            out["attributed_feature_id"] = gm.gene_id
            out["attribution_source"] = gm.source
            # ncRNA continuous coords (169 D6): frac 0 = 5' end, 1 = 3' end
            denom = (gm.body_end - 1) - gm.body_start   # bases between 5' and 3' ends
            if gm.strand == "+":
                out["distance_to_ncrna_3prime"] = pos - (gm.body_end - 1)
                if denom > 0:
                    out["pos_within_ncrna_frac"] = (pos - gm.body_start) / denom
            else:
                out["distance_to_ncrna_3prime"] = gm.body_start - pos
                if denom > 0:
                    out["pos_within_ncrna_frac"] = ((gm.body_end - 1) - pos) / denom
            also.extend(self._other_sense_overlaps(chrom, strand, pos, exclude=gm.gene_id))
        else:
            self._classify_uncontained(chrom, strand, pos, fallback_gene, out)

        out["also_within"] = ";".join(also) if also else None
        return out

    def _classify_within_coding_gene(self, gm: GeneModel, pos: int, out: Dict[str, object]) -> None:
        # canonical-CPA distance (signed; negative = upstream)
        out["distance_to_canonical_CPA"] = (
            pos - gm.canonical_cpa if gm.strand == "+" else gm.canonical_cpa - pos
        )
        # premature = 5' of the stop codon along the transcript
        stop_b = gm.stop_boundary
        if stop_b is not None:
            out["is_premature"] = (pos < gm.cds_end) if gm.strand == "+" else (pos >= gm.cds_start)

        # 1) annotated intron? (F3: of THIS gene OR an inner/overlapping sense gene
        # — the intron identity is attributed to the intron's OWNER, not dropped
        # just because the primary gene is the wider one)
        ann = self._intron_at(self._annotated_intron_trees, gm.chrom, gm.strand, pos)
        if ann is not None:
            owner = self.genes.get(ann["gene_id"], gm)
            self._fill_intron(out, owner, ann, pos, REGION_INTRONIC, INTRON_EVIDENCE_ANNOTATED)
            return
        # 2) cryptic intron? (registered from population junctions)
        cry = self._intron_at(self._cryptic_intron_trees, gm.chrom, gm.strand, pos)
        if cry is not None:
            self._fill_cryptic(out, gm, cry, pos)
            return

        # 3) on the mature transcript — continuous coords are SPLICED (F1): introns
        # 5' of the pA (orf_fraction) or between the pA and the stop (distance) are
        # NOT counted, so intron-bearing genes get transcript-relative coordinates.
        if gm.cds_start is not None and gm.cds_end is not None and gm.cds_start <= pos < gm.cds_end:
            out["region_class"] = REGION_CDS_INTERNAL
            out["orf_fraction"] = self._orf_fraction(gm, pos)
            out["distance_to_stop_codon"] = self._spliced_distance_to_stop(gm, pos)
            return

        u5 = gm.utr5_span()
        if u5 and u5[0] <= pos < u5[1]:
            out["region_class"] = REGION_5UTR
            out["distance_to_stop_codon"] = self._spliced_distance_to_stop(gm, pos)
            return

        u3 = gm.utr3_span()
        if u3 and u3[0] <= pos < u3[1]:
            d = self._spliced_distance_to_stop(gm, pos)
            if d is not None:
                out["distance_to_stop_codon"] = d
                out["region_class"] = (
                    REGION_3UTR_PROXIMAL if 0 < d <= self.utr3_split else REGION_3UTR_DISTAL
                )
            else:
                out["region_class"] = REGION_3UTR_DISTAL
            return

        # within body but no CDS (non-canonical) — treat as CDS_internal-ish fallback
        out["region_class"] = REGION_CDS_INTERNAL if stop_b is not None else REGION_INTERGENIC

    @staticmethod
    def _fill_intron(out, gm, ann, pos, region, evidence) -> None:
        out["region_class"] = region
        out["within_intron_id"] = f"intron{ann['rank']}"
        out["intron_rank"] = ann["rank"]
        out["intron_evidence"] = evidence
        s, e = ann["start"], ann["end"]
        if e > s:
            out["pos_within_intron_frac"] = (pos - s) / (e - s) if gm.strand == "+" else (e - 1 - pos) / (e - s)

    @staticmethod
    def _fill_cryptic(out, gm, cry, pos) -> None:
        out["region_class"] = REGION_CRYPTIC_INTRONIC
        out["within_intron_id"] = f"cryptic:{gm.chrom}:{cry['start']}-{cry['end']}"
        out["intron_evidence"] = INTRON_EVIDENCE_CRYPTIC
        out["spliced_support"] = True
        s, e = cry["start"], cry["end"]
        if e > s:
            out["pos_within_intron_frac"] = (pos - s) / (e - s) if gm.strand == "+" else (e - 1 - pos) / (e - s)

    def _classify_uncontained(self, chrom, strand, pos, fallback_gene, out) -> None:
        # a cryptic intron in an intergenic locus (no body) still counts
        cry = self._intron_at(self._cryptic_intron_trees, chrom, strand, pos)
        if cry is not None:
            out["region_class"] = REGION_CRYPTIC_INTRONIC
            out["within_intron_id"] = f"cryptic:{chrom}:{cry['start']}-{cry['end']}"
            out["intron_evidence"] = INTRON_EVIDENCE_CRYPTIC
            out["spliced_support"] = True
            return
        # readthrough fallback aligned with the attribution proximity window
        if fallback_gene and fallback_gene in self.genes:
            gm = self.genes[fallback_gene]
            d = (pos - gm.canonical_cpa) if gm.strand == "+" else (gm.canonical_cpa - pos)
            if 0 <= d <= self.readthrough_window:
                out["region_class"] = REGION_DOWNSTREAM
                out["attributed_gene"] = gm.gene_id
                out["attributed_feature_id"] = gm.gene_id
                out["attribution_source"] = gm.source
                out["distance_to_canonical_CPA"] = (pos - gm.canonical_cpa) if gm.strand == "+" else (gm.canonical_cpa - pos)
                out["is_premature"] = False
                return
        out["region_class"] = REGION_INTERGENIC

    # ── helpers ──────────────────────────────────────────────────────────────
    @staticmethod
    def _intron_at(trees, chrom, strand, pos):
        tree = trees.get((chrom, strand))
        if not tree:
            return None
        hits = tree[pos]
        if not hits:
            return None
        # most-specific: smallest span
        return min((iv.data for iv in hits), key=lambda d: d["end"] - d["start"])

    def _annotated_intron_owner(self, chrom, strand, pos) -> Optional[GeneModel]:
        """GeneModel owning the most-specific ANNOTATED intron containing ``pos``,
        or None (F3). Selects the same smallest-span intron ``_classify_within_
        coding_gene`` uses, so the gm-switch and the in-method gate stay coherent."""
        ann = self._intron_at(self._annotated_intron_trees, self._c(chrom), strand, pos)
        if ann is None:
            return None
        return self.genes.get(ann["gene_id"])

    # ── spliced-mRNA projection (F1) ─────────────────────────────────────────
    @staticmethod
    def _spliced_coord_5p(gm: GeneModel, pos: int) -> Optional[int]:
        """Spliced-mRNA coordinate of ``pos`` measured from the transcript 5' end
        (0-based nt), or None if ``pos`` is intronic / outside every exon block.
        Increases monotonically 5'->3' along the transcript on BOTH strands."""
        blocks = gm.exon_blocks
        if not blocks:
            return None
        if gm.strand == "+":
            c = 0
            for (s, e) in blocks:                 # ascending == 5'->3' on +
                if s <= pos < e:
                    return c + (pos - s)
                c += (e - s)
        else:
            c = 0
            for (s, e) in reversed(blocks):       # descending == 5'->3' on -
                if s <= pos < e:
                    return c + (e - 1 - pos)
                c += (e - s)
        return None

    def _spliced_distance_to_stop(self, gm: GeneModel, pos: int) -> Optional[int]:
        """Signed SPLICED nt from the stop codon to ``pos`` (negative = upstream /
        premature). Introns between ``pos`` and the stop are not counted (F1).
        Falls back to the genomic signed distance if either endpoint is not
        exonically projectable (e.g. ``pos`` inside an intron)."""
        stop_b = gm.stop_boundary
        if stop_b is None:
            return None
        sp_pos = self._spliced_coord_5p(gm, pos)
        sp_stop = self._spliced_coord_5p(gm, stop_b)
        if sp_pos is None or sp_stop is None:
            return (pos - stop_b) if gm.strand == "+" else (stop_b - pos)
        return sp_pos - sp_stop

    def _orf_fraction(self, gm: GeneModel, pos: int) -> Optional[float]:
        """Spliced-ORF fraction (0 = start codon, 1 = stop codon / last CDS base),
        computed over the SPLICED CDS so a 5'-proximal intron is not counted (F1).
        None if ``pos`` is not projectable onto the spliced ORF (e.g. intronic)."""
        if gm.cds_start is None or gm.cds_end is None:
            return None
        start_codon = gm.cds_start if gm.strand == "+" else (gm.cds_end - 1)
        stop_codon = (gm.cds_end - 1) if gm.strand == "+" else gm.cds_start
        sp_pos = self._spliced_coord_5p(gm, pos)
        sp_start = self._spliced_coord_5p(gm, start_codon)
        sp_stop = self._spliced_coord_5p(gm, stop_codon)
        if sp_pos is None or sp_start is None or sp_stop is None:
            return None
        denom = sp_stop - sp_start
        if denom <= 0:
            return None
        return (sp_pos - sp_start) / denom

    def _antisense_context(self, chrom, strand, pos) -> Optional[str]:
        tree = self._feature_trees.get((chrom, _opp(strand)))
        if not tree:
            return None
        hits = [iv.data for iv in tree[pos]]
        if not hits:
            return None
        order = {"CDS": 0, "3primeUTR": 1, "5primeUTR": 2}
        best = min(hits, key=lambda d: (order.get(d["kind"], 3), d["gene_id"]))
        return f"{best['kind']}:{best['gene_id']}"

    def _other_sense_overlaps(self, chrom, strand, pos, *, exclude: str) -> List[str]:
        out = []
        tree = self._feature_trees.get((chrom, strand))
        if not tree:
            return out
        for iv in tree[pos]:
            d = iv.data
            if d["gene_id"] == exclude:
                continue
            out.append(f"{d['kind']}:{d['gene_id']}")
        return sorted(set(out))


# ─────────────────────────────────────────────────────────────────────────────
# Containment-first cluster→gene attribution (167 Part A)
# ─────────────────────────────────────────────────────────────────────────────
def containment_attributions_from_clusters(
    clusters_df: pd.DataFrame,
    model: TranscriptModel,
    *,
    source: str = "containment",
    secondary_share: float = 0.25,
) -> pd.DataFrame:
    """Attribute each cluster to the gene(s) whose transcribed body contains its
    modal position (167 Part A). Emits the same schema as
    ``cluster_gene_attribution.body_attributions_from_corrected_tsvs`` so it
    plugs into the existing weighted-attribution machinery.

    Multi-gene (nested / overlapping) loci emit ALL sense-overlapping genes with a
    documented PRIMARY-WEIGHTED scheme (F8): the primary gene takes the full read
    count and each additional overlapped gene takes ``secondary_share`` of it,
    then weights are normalized within the cluster. This keeps the primary
    STRICTLY highest-weight (so ``apply_primary_gene_to_clusters`` deterministically
    picks it and canonical ``gene_id`` never coin-flips) while surfacing the other
    genes in the weight table that feeds gene-level DESeq2. A single-gene cluster
    (the common case) emits exactly one row at weight 1.0 — identical to before.
    """
    from .cluster_gene_attribution import _empty_attribution_df, normalize_cluster_gene_attributions

    if clusters_df is None or clusters_df.empty:
        return _empty_attribution_df()

    rows = []
    for row in clusters_df.itertuples(index=False):
        chrom = getattr(row, "chrom")
        strand = getattr(row, "strand")
        pos = int(getattr(row, "modal_position"))
        genes = model.sense_overlapping_genes(chrom, strand, pos)
        if not genes:
            continue
        n_reads = float(getattr(row, "n_reads", 1.0) or 1.0)
        cid = getattr(row, "cluster_id")
        for i, g in enumerate(genes):
            rows.append({
                "cluster_id": cid,
                "gene_id": g["gene_id"],
                "gene_name": g["gene_name"],
                "raw_attributed_count": n_reads if i == 0 else n_reads * secondary_share,
                "source": source,
            })
    if not rows:
        return _empty_attribution_df()
    return normalize_cluster_gene_attributions(pd.DataFrame(rows))


# ─────────────────────────────────────────────────────────────────────────────
# Cryptic-intron read pass + per-cluster classification (167 Part B/C)
# ─────────────────────────────────────────────────────────────────────────────
def collect_cryptic_junctions_from_tsvs(
    samples: Iterable[dict],
    model: TranscriptModel,
    cluster_lookup: Callable[[str, str, int], Optional[str]],
    *,
    chrom_normalizer: Optional[Callable[[str], str]] = None,
    chunk_size: int = 250_000,
) -> Tuple[List[Tuple[str, str, int, int]], set]:
    """Stream corrected TSV(s), returning (cryptic_junctions, clusters_with_splice).

    Reads the ``junctions`` / ``n_junctions`` columns (0-based half-open) that the
    main analyze loader drops. Junctions absent from the model's annotated-intron
    set are cryptic. ``clusters_with_splice`` = set of cluster_ids whose supporting
    reads carry ANY cryptic junction (for ``spliced_support``). Gracefully skips
    samples whose TSV lacks the junction columns (older corrected TSVs).
    """
    from ..unified_record import UnifiedReadRecord

    cryptic: List[Tuple[str, str, int, int]] = []
    spliced_clusters: set = set()
    for sample in samples:
        path = sample.get("path")
        if not path:
            continue
        try:
            header = pd.read_csv(path, sep="\t", nrows=0)
        except (OSError, ValueError):
            continue
        cols = set(header.columns)
        if "junctions" not in cols:
            continue
        pos_col = "corrected_3prime" if "corrected_3prime" in cols else "corrected_position"
        need = {"chrom", "strand", "junctions", pos_col}
        if need - cols:
            continue
        usecols = ["chrom", "strand", pos_col, "junctions"]
        for chunk in pd.read_csv(path, sep="\t", chunksize=chunk_size, usecols=usecols):
            chunk = chunk.dropna(subset=["junctions"])
            if chunk.empty:
                continue
            if chrom_normalizer is not None:
                chunk["chrom"] = chunk["chrom"].map(chrom_normalizer)
            for r in chunk.itertuples(index=False):
                js = getattr(r, "junctions")
                if _is_missing(js):
                    continue
                parsed = UnifiedReadRecord.junctions_from_string(str(js))
                if not parsed:
                    continue
                chrom = getattr(r, "chrom")
                strand = getattr(r, "strand")
                read_cryptic = False
                for (s, e) in parsed:
                    if not model.is_annotated_junction(chrom, strand, s, e):
                        cryptic.append((chrom, strand, int(s), int(e)))
                        read_cryptic = True
                if read_cryptic:
                    try:
                        cid = cluster_lookup(chrom, strand, int(getattr(r, pos_col)))
                    except (TypeError, ValueError):
                        cid = None
                    if cid is not None:
                        spliced_clusters.add(cid)
    return cryptic, spliced_clusters


def annotate_clusters_with_transcript_model(
    clusters_df: pd.DataFrame,
    model: TranscriptModel,
    *,
    samples: Optional[Iterable[dict]] = None,
    cluster_lookup: Optional[Callable[[str, str, int], Optional[str]]] = None,
    chrom_normalizer: Optional[Callable[[str], str]] = None,
    min_cryptic_support: int = 2,
) -> pd.DataFrame:
    """Add per-cluster ``region_class`` + continuous-coordinate columns.

    Owns ONLY the classification columns; never touches ``gene_id`` (that is the
    attribution dispatch's job). ``fallback_gene`` per cluster = its current
    ``gene_id`` so ``downstream_readthrough`` stays coherent with attribution.
    Runs the cryptic-intron read pass when ``samples`` are provided.
    """
    if clusters_df is None or clusters_df.empty:
        for c in CLASSIFIER_COLUMNS:
            if c not in clusters_df.columns:
                clusters_df[c] = None
        return clusters_df

    df = clusters_df.copy()

    spliced_clusters: set = set()
    if samples is not None and cluster_lookup is not None:
        try:
            cryptic, spliced_clusters = collect_cryptic_junctions_from_tsvs(
                samples, model, cluster_lookup, chrom_normalizer=chrom_normalizer,
            )
            if cryptic:
                model.register_cryptic_introns(cryptic, min_support=min_cryptic_support)
        except Exception as exc:  # never let the cryptic pass break the pipeline
            model.warnings.append(f"cryptic-intron read pass failed ({exc}); using annotated introns only")

    gene_col = df["gene_id"] if "gene_id" in df.columns else None
    records: List[Dict[str, object]] = []
    for i, row in enumerate(df.itertuples(index=False)):
        fallback = None
        if gene_col is not None:
            g = gene_col.iloc[i]
            fallback = None if _is_missing(g) else str(g)
        cid = getattr(row, "cluster_id", None)
        rec = model.classify_site(
            getattr(row, "chrom"),
            getattr(row, "strand"),
            int(getattr(row, "modal_position")),
            spliced_support=(cid in spliced_clusters),
            fallback_gene=fallback,
        )
        records.append(rec)

    class_df = pd.DataFrame.from_records(records, columns=list(CLASSIFIER_COLUMNS))
    for c in CLASSIFIER_COLUMNS:
        df[c] = class_df[c].values
    return df
