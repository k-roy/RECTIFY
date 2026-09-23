"""Per-library length-bias correction for ONT PCR-cDNA gene counts.

Replicate ONT PCR-cDNA libraries (SQK-PCB114) count long genes very differently:
relative to its replicates, a 3-kb gene can be counted a fraction of, or many
times, as often as a 0.5-kb gene, while direct RNA (DRS) libraries of the same
strains agree to within ~20 %. The bias is a property of the library
preparation, follows the library's read length, and can change sign between
libraries, so it has to be measured and removed library by library.

The correction measures each library's bias on a PANEL of genes whose true
relative abundance is known from direct RNA (genes that DRS shows to be stable
across the conditions the panel was built for), fits a smooth curve of the
bias against gene length, and divides it out:

    y_g   = log2(CPM_g + 0.5) - log2(REF_g + 0.5)       (panel genes only)
    f(L)  = natural cubic spline in log10 length, knots at panel quantiles,
            held FLAT outside the panel's length range (never extrapolated)
    B_g   = 2 ** f(L_g)                                  (bias factor)
    corrected CPM = CPM / B, renormalized to 1e6

``REF`` is the reference abundance: by default the cohort-mean CPM of direct
RNA libraries counted with the same rule. The curve's intercept is arbitrary
per library; corrected CPM (renormalized), DESeq2 normalization factors (each
gene's geometric mean scaled to 1) and ratios of two lengths are all immune to
it.

Counting rule (``count_cdna_clusters``): one UMI-deduplicated molecule per row
of ``cdna-analyze``'s ``clusters.tsv``; molecules with ``XF >= 1`` of BOTH read
types (Type 2 carries no UMI but its 3' end is a genuine poly(A) site; dropping
it under-counts long isoforms); 3' end = the walk-back ``anchor``; gene = the
gene whose stop codon is the nearest upstream one on the RNA strand, with the
3' end no more than 30 nt upstream of that stop and no more than 500 nt
downstream of it. Genes are CDS-defined (the 3'-most CDS base is the stop);
an isoform-union ``gene`` feature would move the stop into the 3' UTR.

Limits this module enforces rather than papers over:

* Gene counts use the curve at full strength (scale c = 1): it is fitted and
  applied between genes. Within-gene and short-molecule readouts (isoform
  ratios, fragments shorter than the panel's shortest gene) need a scale c
  that is not yet calibrated, so ``ratio_log2_correction`` takes c as a
  REQUIRED keyword and every output records the c used.
* Corrected values are on the reference's scale. A length bias common to every
  reference library is invisible.
* A panel is certified only for the strains, media and conditions in which its
  stability was measured. ``certify_libraries`` flags every library outside
  the panel's declared scope, including libraries with no declared metadata.
* Counts and curves must share one gene-length definition. Stored parameters
  carry a fingerprint of it, and ``apply_length_curves`` refuses a mismatch.
"""
from __future__ import annotations

import gzip
import hashlib
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

# ── Counting rule ─────────────────────────────────────────────────────────────
DEFAULT_MIN_XF = 1
ASSIGN_UPSTREAM_NT = 30      # a 3' end may lie up to 30 nt upstream of (inside) the stop codon
ASSIGN_DOWNSTREAM_NT = 500   # and up to 500 nt downstream of it
MIN_ASSIGNED_SHARE = 0.5     # below this, suspect a contig-naming mismatch, not biology

# ── Fit ───────────────────────────────────────────────────────────────────────
PSEUDO_CPM = 0.5
DEFAULT_MODEL = "ns"
DEFAULT_N_KNOTS = 4
DEFAULT_LENGTH_TYPE = "tx"
LENGTH_TYPES = ("tx", "cds")

# ── Panel construction (stable-gene panel from direct RNA) ─────────────────────
PANEL_MIN_MEAN = 100.0       # mean molecules per library, in BOTH assays
PANEL_SD_MAX = 0.5           # SD (log2) of log2 CPM across the DRS condition-group means
PANEL_STRICT_SD_MAX = 0.35
SLOPE_SET_MIN_MEAN = 20.0    # genes considered at all: mean >= 20 molecules in both assays

SCOPE_KEYS = ("background", "medium", "condition")

_MITO_ALIASES = {
    "chrm", "chrmt", "chrmito", "m", "mt", "mito", "ref|nc_001224|", "nc_001224",
    "nc_001224.1",
}

GENE_TABLE_COLUMNS = ("gene", "chrom", "strand", "cds_start0", "cds_end", "cds_len")

__all__ = [
    "LengthCurve", "FitResult", "StopIndex",
    "chrom_key", "read_gff_cds_models", "gene_table_from_gff", "load_gene_table",
    "resolve_lengths", "lengths_fingerprint", "genes_fingerprint",
    "read_cdna_clusters", "count_cdna_clusters", "count_three_prime_ends",
    "cpm", "reference_from_counts", "build_panel",
    "ns_basis", "library_log2_deviation", "fit_length_curve", "fit_length_bias",
    "bias_factors", "corrected_cpm", "apply_length_curves",
    "deseq2_normalization_factors", "ratio_log2_correction", "ratio_corrections",
    "certify_libraries", "replicate_slopes", "replicate_spread",
    "params_to_frame", "curves_from_params",
]


# ═════════════════════════════════════════════════════════════════════════════
# Gene table
# ═════════════════════════════════════════════════════════════════════════════

def chrom_key(name) -> str:
    """Contig key shared by gene tables and 3'-end tables.

    Strips a leading ``chr`` so UCSC-style (``chrI``, ``chr1``) and bare names
    (``I``, ``1``) meet, and maps every mitochondrial alias to ``MITO``.
    """
    s = str(name)
    if s.lower() in _MITO_ALIASES:
        return "MITO"
    if len(s) > 3 and s[:3].lower() == "chr":
        s = s[3:]
    return s


def _open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def _attrs(col9: str) -> Dict[str, str]:
    out = {}
    for a in col9.strip().split(";"):
        a = a.strip()
        if "=" in a:
            k, v = a.split("=", 1)
            out[k] = v
    return out


def _strip_prefix(s: str) -> str:
    for p in ("gene:", "transcript:"):
        if s.startswith(p):
            return s[len(p):]
    return s


def read_gff_cds_models(gff) -> Tuple[Dict[str, dict], Dict[str, set]]:
    """Read CDS intervals and introns per gene from a GFF3 (plain or .gz).

    Gene identity of a CDS row: SGD style ``Name=<ORF>_CDS`` first; otherwise
    the row's ``Parent`` transcript, mapped to its gene through the transcript
    row's own ``Parent`` (generic GFF3). Parsing stops at ``##FASTA``.

    Returns ``(cds, introns)``: ``cds[gene] = {'chrom', 'strand', 'intervals'}``
    with 1-based inclusive intervals; ``introns[gene]`` is a set of 1-based
    inclusive (start, end) pairs, keyed by the SGD ``Parent`` prefix
    (``YAL001C_mRNA`` -> ``YAL001C``).
    """
    tx_gene: Dict[str, str] = {}
    cds_rows = []
    intron_rows = []
    with _open_text(gff) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            if ftype == "CDS":
                cds_rows.append((f[0], f[6], int(f[3]), int(f[4]), f[8]))
            elif ftype == "intron":
                intron_rows.append((int(f[3]), int(f[4]), f[8]))
            elif ftype in ("mRNA", "transcript"):
                a = _attrs(f[8])
                if "ID" in a and "Parent" in a:
                    tx_gene[a["ID"]] = a["Parent"].split(",")[0]

    cds: Dict[str, dict] = {}
    for chrom, strand, s1, e1, col9 in cds_rows:
        a = _attrs(col9)
        name = a.get("Name", "")
        if name.endswith("_CDS"):
            gene = name[:-4]
        elif "Parent" in a:
            parent = a["Parent"].split(",")[0]
            gene = _strip_prefix(tx_gene.get(parent, parent))
        else:
            continue
        rec = cds.get(gene)
        if rec is None:
            cds[gene] = {"chrom": chrom, "strand": strand, "intervals": {(s1, e1)}}
        else:
            if rec["chrom"] != chrom or rec["strand"] != strand:
                raise ValueError(f"GFF: CDS of {gene} on two contigs or strands "
                                 f"({rec['chrom']}{rec['strand']} vs {chrom}{strand})")
            rec["intervals"].add((s1, e1))

    introns: Dict[str, set] = {}
    for s1, e1, col9 in intron_rows:
        par = _attrs(col9).get("Parent", "")
        for p in par.split(","):
            if p:
                introns.setdefault(p.split("_")[0], set()).add((s1, e1))
    return cds, introns


def _union_len(intervals) -> int:
    """Total length of the union of 1-based inclusive intervals (adjacent ones merge)."""
    iv = sorted(intervals)
    tot = 0
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s > ce + 1:
            tot += ce - cs + 1
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + ce - cs + 1


def gene_table_from_gff(
    gff,
    *,
    transcript_models: Optional[pd.DataFrame] = None,
    require_transcript_model: bool = False,
    exclude_genes: Iterable[str] = (),
    exclude_regions: Sequence[Tuple[str, int, int]] = (),
    drop_mito: bool = True,
) -> pd.DataFrame:
    """Build the CDS-defined gene table the counting rule and the curve share.

    Columns: gene, chrom, strand, cds_start0 (0-based), cds_end (half-open),
    cds_span, cds_len (spliced CDS: union of unique CDS intervals across
    isoforms), tx_len (spliced transcript length, NaN without a model).

    ``transcript_models`` (columns gene, tx_lo, tx_hi; 1-based inclusive) give
    ``tx_len = tx_hi - tx_lo + 1`` minus every annotated intron that lies
    inside it. ``exclude_regions`` entries are ``(contig, start0, end)``: a gene
    whose CDS span overlaps one is dropped. Exclusions REMOVE a gene from the
    candidate-stop set, so a 3' end near an excluded gene is not reassigned to
    it; this is part of the counting rule, not a post-filter.
    """
    cds, introns = read_gff_cds_models(gff)
    rows = []
    for gene, rec in cds.items():
        iv = rec["intervals"]
        s0 = min(s for s, _ in iv) - 1
        e = max(e for _, e in iv)
        rows.append((gene, rec["chrom"], rec["strand"], s0, e, e - s0, _union_len(iv)))
    gt = pd.DataFrame(rows, columns=["gene", "chrom", "strand", "cds_start0", "cds_end",
                                     "cds_span", "cds_len"])
    gt["tx_len"] = np.nan
    if transcript_models is not None:
        tm = transcript_models.drop_duplicates("gene").set_index("gene")
        have = gt.gene.isin(tm.index)
        if require_transcript_model:
            gt = gt[have].copy()
            have = pd.Series(True, index=gt.index)

        def _tx_len(g):
            r = tm.loc[g]
            L = int(r.tx_hi) - int(r.tx_lo) + 1
            for s, e in introns.get(g, ()):
                if s >= r.tx_lo and e <= r.tx_hi:
                    L -= e - s + 1
            return float(L)

        gt.loc[have, "tx_len"] = [_tx_len(g) for g in gt.gene[have]]
    keys = gt.chrom.map(chrom_key)
    drop = gt.gene.isin(set(exclude_genes))
    if drop_mito:
        drop |= keys.eq("MITO")
    for contig, start0, end in exclude_regions:
        drop |= keys.eq(chrom_key(contig)) & (gt.cds_end > start0) & (gt.cds_start0 < end)
    return gt[~drop].reset_index(drop=True)


def load_gene_table(path) -> pd.DataFrame:
    """Read a gene table TSV and check the columns the counting rule needs."""
    gt = pd.read_csv(path, sep="\t", dtype={"gene": str, "chrom": str, "strand": str})
    missing = [c for c in GENE_TABLE_COLUMNS if c not in gt.columns]
    if missing:
        raise ValueError(f"gene table {path} lacks column(s) {missing}; expected at least "
                         f"{list(GENE_TABLE_COLUMNS)} (+ tx_len for --length-type tx)")
    if gt.gene.duplicated().any():
        dup = gt.gene[gt.gene.duplicated()].iloc[0]
        raise ValueError(f"gene table {path}: gene {dup!r} appears twice")
    if "tx_len" not in gt.columns:
        gt["tx_len"] = np.nan
    return gt


def resolve_lengths(gene_table: pd.DataFrame, length_type: str = DEFAULT_LENGTH_TYPE
                    ) -> Tuple[pd.Series, dict]:
    """Per-gene length (nt) for the curve, indexed by gene.

    ``tx``: spliced transcript length; genes without one take their CDS length
    plus the median (transcript - CDS) difference of the genes that have both,
    and the number filled is reported. ``cds``: spliced CDS length.
    """
    if length_type not in LENGTH_TYPES:
        raise ValueError(f"length_type must be one of {LENGTH_TYPES}, not {length_type!r}")
    g = gene_table.set_index("gene")
    cds_len = g["cds_len"].astype(float)
    info = {"length_type": length_type, "n_filled": 0, "fill_utr_nt": None}
    if length_type == "cds":
        return cds_len.rename("length_nt"), info
    tx = g["tx_len"].astype(float)
    both = tx.notna() & cds_len.notna()
    if not both.any():
        raise ValueError("the gene table has no transcript lengths (tx_len); supply them or use "
                         "--length-type cds")
    missing = tx.isna()
    if missing.any():
        utr = float((tx[both] - cds_len[both]).median())
        tx = tx.where(~missing, cds_len + utr)
        info.update(n_filled=int(missing.sum()), fill_utr_nt=utr)
    return tx.rename("length_nt"), info


def lengths_fingerprint(lengths_nt: pd.Series) -> str:
    """sha256 of the gene-length definition (sorted gene\\tlength lines)."""
    s = lengths_nt.dropna().sort_index()
    h = hashlib.sha256()
    for g, v in s.items():
        h.update(f"{g}\t{float(v)!r}\n".encode())
    return h.hexdigest()


def genes_fingerprint(genes: Iterable[str]) -> str:
    """sha256 of a sorted gene list (panel identity)."""
    h = hashlib.sha256()
    for g in sorted(set(map(str, genes))):
        h.update(f"{g}\n".encode())
    return h.hexdigest()


# ═════════════════════════════════════════════════════════════════════════════
# Counting: 3' end -> gene
# ═════════════════════════════════════════════════════════════════════════════

class StopIndex:
    """Sorted stop codons per (contig key, strand) for nearest-upstream-stop assignment.

    The stop is the 3'-most CDS base on the RNA strand: ``cds_end - 1`` on
    ``+``, ``cds_start0`` on ``-``.
    """

    def __init__(self, gene_table: pd.DataFrame):
        self._idx: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]] = {}
        gt = gene_table.assign(_ck=gene_table.chrom.map(chrom_key))
        for (ck, st), g in gt.groupby(["_ck", "strand"], sort=False):
            if st == "+":
                stop = g.cds_end.values.astype(np.int64) - 1
            elif st == "-":
                stop = g.cds_start0.values.astype(np.int64)
            else:
                continue
            o = np.argsort(stop, kind="stable")
            self._idx[(ck, st)] = (stop[o], g.gene.values[o])
        self.contigs = sorted({k for k, _ in self._idx})

    def assign(self, contig: str, strand: str, pos, *, upstream: int = ASSIGN_UPSTREAM_NT,
               downstream: int = ASSIGN_DOWNSTREAM_NT) -> np.ndarray:
        """Gene for each 0-based 3' end on one contig and RNA strand (None = unassigned)."""
        x = np.asarray(pos, dtype=np.int64)
        out = np.full(len(x), None, dtype=object)
        k = self._idx.get((chrom_key(contig), strand))
        if k is None or len(x) == 0:
            return out
        stop, gn = k
        if strand == "+":
            j = np.searchsorted(stop, x + upstream, "right") - 1
            ok = j >= 0
            ok[ok] &= (x[ok] - stop[j[ok]]) <= downstream
        else:
            j = np.searchsorted(stop, x - upstream, "left")
            ok = j < len(stop)
            ok[ok] &= (stop[j[ok]] - x[ok]) <= downstream
        out[ok] = gn[j[ok]]
        return out


def read_cdna_clusters(path) -> pd.DataFrame:
    """The four ``clusters.tsv`` columns the counting rule reads (plain or .gz)."""
    need = ["chrom", "orient", "anchor", "xf"]
    head = pd.read_csv(path, sep="\t", nrows=0)
    missing = [c for c in need if c not in head.columns]
    if missing:
        raise ValueError(f"{path}: not a cdna-analyze clusters.tsv (missing {missing})")
    return pd.read_csv(path, sep="\t", usecols=need, dtype={"chrom": str, "orient": str})


def count_three_prime_ends(ends: pd.DataFrame, index: StopIndex, *,
                           upstream: int = ASSIGN_UPSTREAM_NT,
                           downstream: int = ASSIGN_DOWNSTREAM_NT) -> Tuple[pd.Series, dict]:
    """Count RNA 3' ends per gene.

    ``ends`` columns: chrom, strand ('+'/'-', the RNA strand), pos (0-based 3'
    end) and optional ``weight`` (default 1 per row). Use this for direct RNA
    3' ends (one row per read, or pre-aggregated with weights) so a reference
    is counted with exactly the rule used for cDNA.
    """
    w = ends["weight"].to_numpy(dtype=float) if "weight" in ends.columns else np.ones(len(ends))
    genes = np.full(len(ends), None, dtype=object)
    strand = ends["strand"].astype(str).to_numpy()
    chrom = ends["chrom"].astype(str).to_numpy()
    pos = ends["pos"].to_numpy(dtype=np.int64)
    bad_strand = ~np.isin(strand, ["+", "-"])
    for (c, st), idx in pd.Series(np.arange(len(ends))).groupby([chrom, strand]).groups.items():
        if st not in ("+", "-"):
            continue
        idx = np.asarray(idx)
        genes[idx] = index.assign(c, st, pos[idx], upstream=upstream, downstream=downstream)
    ok = pd.notna(genes)
    counts = pd.Series(w[ok], index=genes[ok]).groupby(level=0).sum()
    total = float(w.sum())
    stats = {
        "n_input": int(len(ends)),
        "weight_input": total,
        "weight_assigned": float(counts.sum()),
        "assigned_share": float(counts.sum() / total) if total else float("nan"),
        "n_bad_strand": int(bad_strand.sum()),
        "contigs_without_genes": sorted({chrom_key(c) for c in set(chrom)} - set(index.contigs)),
    }
    return counts, stats


def count_cdna_clusters(clusters: Union[pd.DataFrame, str, Path], index: StopIndex, *,
                        min_xf: int = DEFAULT_MIN_XF,
                        upstream: int = ASSIGN_UPSTREAM_NT,
                        downstream: int = ASSIGN_DOWNSTREAM_NT,
                        min_assigned_share: Optional[float] = MIN_ASSIGNED_SHARE,
                        ) -> Tuple[pd.Series, dict]:
    """Gene counts of one cDNA library from ``cdna-analyze``'s ``clusters.tsv``.

    One row = one UMI-deduplicated molecule (Type 1) or one UMI-less Type 2
    read; ``n_reads`` is NOT used. Rows with ``xf < min_xf`` are dropped (a
    Type-1 ``XF 0`` molecule has no poly(A) signature, so its 'anchor' is a
    truncation point inside the transcript); both read types are kept.
    ``orient`` fwd/rev is the RNA strand (+/-).

    Raises ``ValueError`` when fewer than ``min_assigned_share`` of the kept
    molecules reach a gene: that is the signature of a contig-naming mismatch
    between the BAM reference and the gene table, which would otherwise zero
    every count silently. Pass ``min_assigned_share=None`` to skip the check.
    """
    d = read_cdna_clusters(clusters) if not isinstance(clusters, pd.DataFrame) else clusters
    n_rows = len(d)
    d = d[d["xf"].astype(float) >= min_xf]
    strand = d["orient"].map({"fwd": "+", "rev": "-"})
    ends = pd.DataFrame({"chrom": d["chrom"].astype(str).to_numpy(),
                         "strand": strand.fillna("?").to_numpy(),
                         "pos": d["anchor"].to_numpy(dtype=np.int64)})
    counts, st = count_three_prime_ends(ends, index, upstream=upstream, downstream=downstream)
    stats = {
        "n_molecules": n_rows,
        "n_xf_pass": int(len(d)),
        "n_assigned": int(round(st["weight_assigned"])),
        "assigned_share": st["assigned_share"],
        "n_orient_unknown": st["n_bad_strand"],
        "contigs_without_genes": st["contigs_without_genes"],
        "min_xf": min_xf,
    }
    if (min_assigned_share is not None and len(d)
            and not (stats["assigned_share"] >= min_assigned_share)):
        raise ValueError(
            f"only {stats['assigned_share']:.1%} of {len(d):,} molecules were assigned to a gene "
            f"(< {min_assigned_share:.0%}). Suspect a contig-naming mismatch between the clusters "
            f"(contigs {sorted(set(d['chrom'].astype(str)))[:5]}...) and the gene table "
            f"(contig keys {index.contigs[:5]}...).")
    return counts.astype(np.int64), stats


# ═════════════════════════════════════════════════════════════════════════════
# Reference and panel
# ═════════════════════════════════════════════════════════════════════════════

def cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Counts per million over the genes (rows) of the table."""
    return counts / counts.sum() * 1e6


def reference_from_counts(ref_counts: pd.DataFrame, genes: Optional[Iterable[str]] = None
                          ) -> pd.Series:
    """Cohort-mean CPM of reference (e.g. direct RNA) libraries.

    Each library's CPM is computed over ``genes`` (default: every row), with
    genes absent from the table counted as 0, then averaged across libraries.
    """
    m = ref_counts if genes is None else ref_counts.reindex(list(genes)).fillna(0)
    return cpm(m).mean(axis=1).rename("ref_cpm")


def build_panel(cdna_counts: pd.DataFrame, ref_counts: pd.DataFrame,
                ref_groups: Mapping[str, str], *, genes: Optional[Iterable[str]] = None,
                min_mean: float = PANEL_MIN_MEAN, sd_max: float = PANEL_SD_MAX,
                strict_sd_max: float = PANEL_STRICT_SD_MAX,
                slope_set_min_mean: float = SLOPE_SET_MIN_MEAN) -> pd.DataFrame:
    """Stable-gene panel: high counts in both assays and stable reference abundance.

    A gene enters when its mean count is >= ``min_mean`` in the cDNA libraries
    AND in the reference libraries, and the SD (ddof 1) of log2 CPM across the
    reference condition-group means (``ref_groups``: library -> group) is
    <= ``sd_max``. Stability is therefore MEASURED across the conditions the
    reference covers, not assumed; the panel is valid only inside them.
    Considered genes: mean >= ``slope_set_min_mean`` in both assays.
    Returns one row per considered gene with ``panel`` and ``strict`` flags.
    """
    keep = cdna_counts.index.intersection(ref_counts.index)
    if genes is not None:
        keep = pd.Index(list(genes)).intersection(keep)
    keep = keep[(cdna_counts.loc[keep].mean(axis=1) >= slope_set_min_mean)
                & (ref_counts.loc[keep].mean(axis=1) >= slope_set_min_mean)]
    rc = cpm(ref_counts).loc[keep]
    groups: Dict[str, List[str]] = {}
    for lib, grp in ref_groups.items():
        if lib in rc.columns:
            groups.setdefault(grp, []).append(lib)
    if len(groups) < 2:
        raise ValueError("build_panel needs reference libraries from at least two condition groups")
    setmean = pd.DataFrame({g: np.log2(rc[libs].mean(axis=1) + PSEUDO_CPM)
                            for g, libs in groups.items()})
    out = pd.DataFrame({
        "cdna_mean": cdna_counts.loc[keep].mean(axis=1),
        "ref_mean": ref_counts.loc[keep].mean(axis=1),
        "ref_cpm": rc.mean(axis=1),
        "sd_sets": setmean.std(axis=1, ddof=1),
        "n_sets": len(groups),
    })
    high = (out.cdna_mean >= min_mean) & (out.ref_mean >= min_mean)
    out["panel"] = high & (out.sd_sets <= sd_max)
    out["strict"] = high & (out.sd_sets <= strict_sd_max)
    out.index.name = "gene"
    return out


# ═════════════════════════════════════════════════════════════════════════════
# Curve
# ═════════════════════════════════════════════════════════════════════════════

def ns_basis(x, knots) -> np.ndarray:
    """Natural cubic spline basis (Hastie, Tibshirani & Friedman, ESL eq. 5.4-5.5).

    Columns: 1, x and K-2 truncated-power terms; linear beyond the boundary
    knots. With K knots the basis has K columns.
    """
    x = np.asarray(x, float)
    xi = np.asarray(knots, float)
    K = len(xi)

    def d(k):
        return (np.clip(x - xi[k], 0, None) ** 3 - np.clip(x - xi[-1], 0, None) ** 3) / (xi[-1] - xi[k])

    return np.column_stack([np.ones_like(x), x] + [d(k) - d(K - 2) for k in range(K - 2)])


@dataclass
class LengthCurve:
    """One library's bias curve f(L), in log2, as a function of length in nt.

    ``x_lo``/``x_hi`` are the log10 length range of the panel genes that set
    it; lengths outside are clamped to the edge, so the curve is flat there
    and never extrapolated. Optional additive covariates (per-gene features
    such as GC content) are z-scored with the panel's mean/SD and clamped to
    the panel's z range, for the same reason.
    """
    model: str
    x_lo: float
    x_hi: float
    coef: np.ndarray
    knots: Optional[np.ndarray] = None
    length_type: str = DEFAULT_LENGTH_TYPE
    cov_names: Tuple[str, ...] = ()
    cov_center: np.ndarray = field(default_factory=lambda: np.zeros(0))
    cov_scale: np.ndarray = field(default_factory=lambda: np.ones(0))
    cov_zlo: np.ndarray = field(default_factory=lambda: np.zeros(0))
    cov_zhi: np.ndarray = field(default_factory=lambda: np.zeros(0))

    @property
    def len_lo_nt(self) -> float:
        return 10 ** self.x_lo

    @property
    def len_hi_nt(self) -> float:
        return 10 ** self.x_hi

    def length_design(self, x_log10) -> np.ndarray:
        xc = np.clip(np.asarray(x_log10, float), self.x_lo, self.x_hi)
        if self.model == "linear":
            return np.column_stack([np.ones(len(xc)), xc])
        return ns_basis(xc, self.knots)

    def cov_design(self, cov: Optional[np.ndarray], n: int) -> np.ndarray:
        """z-scored, clamped covariate columns; a missing value sits at the panel mean (z = 0)."""
        k = len(self.cov_names)
        if k == 0:
            return np.zeros((n, 0))
        if cov is None:
            raise ValueError(f"this curve uses covariates {list(self.cov_names)}; supply them")
        z = (np.asarray(cov, float).reshape(n, k) - self.cov_center) / self.cov_scale
        z = np.clip(z, self.cov_zlo, self.cov_zhi)
        return np.where(np.isfinite(z), z, 0.0)

    def design(self, x_log10, cov: Optional[np.ndarray] = None) -> np.ndarray:
        X = self.length_design(x_log10)
        return np.hstack([X, self.cov_design(cov, X.shape[0])])

    def log2_bias(self, lengths_nt, cov: Optional[np.ndarray] = None) -> np.ndarray:
        """f at each length (nt): the library's log2 over-count relative to the reference."""
        L = np.asarray(lengths_nt, float)
        return self.design(np.log10(L), cov) @ self.coef

    def clamped(self, lengths_nt) -> np.ndarray:
        """True where a length lies outside the fitted range (the edge value is used)."""
        x = np.log10(np.asarray(lengths_nt, float))
        return (x < self.x_lo) | (x > self.x_hi)

    def slope_equiv(self) -> float:
        """Average slope over the fitted range, log2 per tenfold length (exact for 'linear')."""
        f = self.length_design([self.x_lo, self.x_hi]) @ self.coef[: self.length_design([0.0]).shape[1]]
        return float((f[1] - f[0]) / (self.x_hi - self.x_lo))


def library_log2_deviation(counts: pd.DataFrame, reference_cpm: pd.Series,
                           pseudo: float = PSEUDO_CPM) -> pd.DataFrame:
    """y = log2(CPM + pseudo) - log2(REF + pseudo), CPM over the table's genes."""
    ref = reference_cpm.reindex(counts.index)
    return np.log2(cpm(counts) + pseudo).sub(np.log2(ref + pseudo), axis=0)


def fit_length_curve(x_log10: np.ndarray, y: np.ndarray, *, model: str = DEFAULT_MODEL,
                     n_knots: int = DEFAULT_N_KNOTS, length_type: str = DEFAULT_LENGTH_TYPE,
                     cov: Optional[np.ndarray] = None, cov_names: Sequence[str] = ()
                     ) -> Tuple[LengthCurve, dict]:
    """Least-squares fit of one library's curve on panel genes.

    Knots sit at equally spaced quantiles of the panel's log10 lengths
    (0, 1/3, 2/3, 1 for 4 knots) and are FIXED at fit time; apply never
    re-derives them. Returns the curve and its diagnostics (residual SD,
    equivalent slope, and a linear fit's slope and SE on the same genes).
    """
    x = np.asarray(x_log10, float)
    y = np.asarray(y, float)
    if model not in ("ns", "linear"):
        raise ValueError(f"model must be 'ns' or 'linear', not {model!r}")
    if model == "ns" and n_knots < 3:
        raise ValueError("a natural spline needs at least 3 knots")
    ok = np.isfinite(x) & np.isfinite(y)
    k = len(cov_names)
    C = None
    if k:
        C = np.asarray(cov, float).reshape(len(x), k)
    x, y = x[ok], y[ok]
    if C is not None:
        C = C[ok]
    lo, hi = float(x.min()), float(x.max())
    knots = np.quantile(x, np.linspace(0, 1, n_knots)) if model == "ns" else None
    curve = LengthCurve(model=model, x_lo=lo, x_hi=hi, coef=np.zeros(0), knots=knots,
                        length_type=length_type)
    if k:
        center = np.nanmean(C, axis=0)
        scale = np.nanstd(C, axis=0, ddof=1)
        if np.any(~(scale > 0)):
            raise ValueError(f"covariate(s) with no spread across the panel: "
                             f"{[n for n, s in zip(cov_names, scale) if not s > 0]}")
        z = (C - center) / scale
        curve.cov_names = tuple(cov_names)
        curve.cov_center, curve.cov_scale = center, scale
        curve.cov_zlo, curve.cov_zhi = np.nanmin(z, axis=0), np.nanmax(z, axis=0)
    X = curve.design(x, C)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    curve.coef = beta
    res = y - X @ beta
    dof = len(y) - X.shape[1]
    if dof <= 0:
        raise ValueError(f"{len(y)} panel genes cannot support a {X.shape[1]}-parameter curve")
    s2 = float((res ** 2).sum() / dof)
    Xl = np.column_stack([np.ones(len(x)), x])
    bl, *_ = np.linalg.lstsq(Xl, y, rcond=None)
    rl = y - Xl @ bl
    s2l = float((rl ** 2).sum() / (len(y) - 2))
    diag = {
        "n_panel_used": int(len(y)),
        "panel_resid_sd": math.sqrt(s2),
        "dof": int(dof),
        "slope_equiv": curve.slope_equiv(),
        "linear_slope": float(bl[1]),
        "linear_slope_se": math.sqrt(s2l * np.linalg.inv(Xl.T @ Xl)[1, 1]),
    }
    return curve, diag


@dataclass
class FitResult:
    curves: Dict[str, LengthCurve]
    params: pd.DataFrame
    qc: pd.DataFrame
    factors: pd.DataFrame
    corrected_cpm: pd.DataFrame
    genes: pd.Index
    lengths_nt: pd.Series
    y: pd.DataFrame


def _covariate_matrix(covariates: Optional[pd.DataFrame], genes: Sequence[str]
                      ) -> Tuple[Optional[np.ndarray], Tuple[str, ...]]:
    if covariates is None or covariates.shape[1] == 0:
        return None, ()
    c = covariates.reindex(list(genes)).astype(float)
    return c.to_numpy(), tuple(map(str, covariates.columns))


def fit_length_bias(counts: pd.DataFrame, reference_cpm: pd.Series, panel: Iterable[str],
                    lengths_nt: pd.Series, *, model: str = DEFAULT_MODEL,
                    n_knots: int = DEFAULT_N_KNOTS, length_type: str = DEFAULT_LENGTH_TYPE,
                    covariates: Optional[pd.DataFrame] = None,
                    pseudo: float = PSEUDO_CPM) -> FitResult:
    """Fit every library's curve on the panel and correct the whole table.

    ``counts``: genes x libraries (molecules, same counting rule as the
    reference). Genes without a length are dropped (they cannot be corrected)
    and CPM is computed over the rest. Panel genes need a finite reference.
    """
    genes = counts.index.intersection(lengths_nt.dropna().index)
    if len(genes) == 0:
        raise ValueError("no gene of the count table has a length in the gene table")
    C = counts.loc[genes].astype(float)
    L = lengths_nt.loc[genes].astype(float)
    x_all = np.log10(L.to_numpy())
    Y = library_log2_deviation(C, reference_cpm, pseudo)
    ref = reference_cpm.reindex(genes)
    pan = pd.Index([g for g in dict.fromkeys(map(str, panel)) if g in set(genes)])
    pan = pan[np.isfinite(ref.loc[pan].to_numpy(dtype=float))]
    if len(pan) == 0:
        raise ValueError("no panel gene is present in both the count table and the reference")
    x_pan = np.log10(L.loc[pan].to_numpy())
    cov_pan, cov_names = _covariate_matrix(covariates, pan)
    cov_all, _ = _covariate_matrix(covariates, genes)
    curves, prow, qrow, fac = {}, [], [], {}
    for lib in C.columns:
        curve, diag = fit_length_curve(x_pan, Y.loc[pan, lib].to_numpy(), model=model,
                                       n_knots=n_knots, length_type=length_type,
                                       cov=cov_pan, cov_names=cov_names)
        curves[lib] = curve
        f = curve.design(x_all, cov_all) @ curve.coef
        fac[lib] = 2.0 ** f
        clamp = (x_all < curve.x_lo) | (x_all > curve.x_hi)
        prow.append(_params_row(lib, curve))
        qrow.append(dict(library=lib, **diag, n_genes=int(len(genes)),
                         n_genes_clamped=int(clamp.sum()),
                         n_clamped_below=int((x_all < curve.x_lo).sum()),
                         n_clamped_above=int((x_all > curve.x_hi).sum())))
    B = pd.DataFrame(fac, index=genes)
    return FitResult(curves=curves, params=pd.DataFrame(prow), qc=pd.DataFrame(qrow),
                     factors=B, corrected_cpm=corrected_cpm(C, B), genes=genes,
                     lengths_nt=L, y=Y)


# ═════════════════════════════════════════════════════════════════════════════
# Stored parameters
# ═════════════════════════════════════════════════════════════════════════════

def _join(a) -> str:
    return ";".join(repr(float(v)) for v in np.asarray(a, float).ravel())


def _split(s) -> np.ndarray:
    if s is None or (isinstance(s, float) and math.isnan(s)) or str(s) == "":
        return np.zeros(0)
    return np.array([float(v) for v in str(s).split(";")])


def _params_row(lib: str, c: LengthCurve) -> dict:
    return {
        "library": lib,
        "model": c.model if c.model == "linear" else f"ns{len(c.knots)}",
        "length_type": c.length_type,
        "coefs": _join(c.coef),
        "knots": _join(c.knots) if c.knots is not None else "",
        "x_lo": repr(float(c.x_lo)),
        "x_hi": repr(float(c.x_hi)),
        "len_lo_nt": c.len_lo_nt,
        "len_hi_nt": c.len_hi_nt,
        "cov_names": ";".join(c.cov_names),
        "cov_center": _join(c.cov_center),
        "cov_scale": _join(c.cov_scale),
        "cov_zlo": _join(c.cov_zlo),
        "cov_zhi": _join(c.cov_zhi),
    }


def params_to_frame(curves: Mapping[str, LengthCurve], **constant_columns) -> pd.DataFrame:
    """One row per library; floats written with repr() so a round trip is exact."""
    df = pd.DataFrame([_params_row(lib, c) for lib, c in curves.items()])
    for k, v in constant_columns.items():
        df[k] = v
    return df


def curves_from_params(params: pd.DataFrame) -> Dict[str, LengthCurve]:
    """Rebuild curves from a params table (the inverse of ``params_to_frame``)."""
    out = {}
    for r in params.to_dict("records"):
        model = str(r["model"])
        names = tuple(n for n in str(r.get("cov_names") or "").split(";") if n and n != "nan")
        c = LengthCurve(
            model="linear" if model == "linear" else "ns",
            x_lo=float(r["x_lo"]), x_hi=float(r["x_hi"]), coef=_split(r["coefs"]),
            knots=None if model == "linear" else _split(r["knots"]),
            length_type=str(r.get("length_type", DEFAULT_LENGTH_TYPE)),
            cov_names=names,
            cov_center=_split(r.get("cov_center")), cov_scale=_split(r.get("cov_scale")),
            cov_zlo=_split(r.get("cov_zlo")), cov_zhi=_split(r.get("cov_zhi")),
        )
        out[str(r["library"])] = c
    return out


# ═════════════════════════════════════════════════════════════════════════════
# Apply
# ═════════════════════════════════════════════════════════════════════════════

def bias_factors(curves: Mapping[str, LengthCurve], lengths_nt: pd.Series,
                 covariates: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """B = 2^f per gene (rows, as ``lengths_nt``) and library (columns)."""
    L = lengths_nt.astype(float)
    cov_all, _ = _covariate_matrix(covariates, L.index)
    return pd.DataFrame({lib: 2.0 ** c.log2_bias(L.to_numpy(), cov_all)
                         for lib, c in curves.items()}, index=L.index)


def corrected_cpm(counts: pd.DataFrame, factors: pd.DataFrame) -> pd.DataFrame:
    """CPM / B, renormalized to 1e6 per library."""
    c = cpm(counts.astype(float)) / factors.loc[counts.index, counts.columns]
    return c / c.sum() * 1e6


def apply_length_curves(counts: pd.DataFrame, params: pd.DataFrame, lengths_nt: pd.Series, *,
                        covariates: Optional[pd.DataFrame] = None,
                        lengths_sha256: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Apply stored curves to a count table of the same libraries.

    Refuses when the gene-length definition differs from the one the curves
    were fitted with (``params['lengths_sha256']``), because a curve is a
    function of THAT length definition. Returns (factors, corrected CPM, qc).
    """
    if "lengths_sha256" in params.columns and params["lengths_sha256"].notna().any():
        want = set(params["lengths_sha256"].dropna().astype(str))
        have = lengths_sha256 or lengths_fingerprint(lengths_nt)
        if want != {have}:
            raise ValueError(
                "gene-length definition differs from the one the curves were fitted with "
                f"(params {sorted(want)[0][:12]}..., supplied {have[:12]}...). Use the same gene "
                "table and --length-type as the fit.")
    curves = curves_from_params(params)
    missing = [lib for lib in counts.columns if lib not in curves]
    if missing:
        raise ValueError(f"no stored curve for librar{'y' if len(missing) == 1 else 'ies'} "
                         f"{missing[:5]}")
    genes = counts.index.intersection(lengths_nt.dropna().index)
    C = counts.loc[genes].astype(float)
    B = bias_factors({lib: curves[lib] for lib in C.columns}, lengths_nt.loc[genes], covariates)
    qc = []
    for lib in C.columns:
        cl = curves[lib].clamped(lengths_nt.loc[genes].to_numpy(float))
        qc.append(dict(library=lib, n_genes=int(len(genes)), n_genes_clamped=int(cl.sum())))
    return B, corrected_cpm(C, B), pd.DataFrame(qc)


def deseq2_normalization_factors(counts: pd.DataFrame, factors: pd.DataFrame,
                                 size_factors: Optional[pd.Series] = None) -> pd.DataFrame:
    """DESeq2 ``normalizationFactors`` that carry the length bias.

    Mirrors DESeq2's ``estimateNormFactors`` with ``normMatrix = B``: B is
    scaled so each gene's geometric mean is 1, size factors are the
    median-of-ratios of ``counts / B`` over positive counts of genes with a
    finite log geometric mean (or ``size_factors`` if given), and the product
    is again scaled so each gene's geometric mean is 1.
    """
    B = factors.loc[counts.index, counts.columns].astype(float)
    nm = B.div(np.exp(np.log(B).mean(axis=1)), axis=0)
    if size_factors is None:
        with np.errstate(divide="ignore"):
            lnc = np.log(counts.astype(float)) - np.log(nm)
        lgm = lnc.mean(axis=1)
        fin = np.isfinite(lgm)
        if not fin.any():
            raise ValueError("every gene has a zero count in some library; cannot estimate size factors")
        sf = {}
        for lib in counts.columns:
            ok = fin & (counts[lib] > 0)
            sf[lib] = float(np.exp(np.median((lnc.loc[ok, lib] - lgm[ok]).to_numpy())))
        size_factors = pd.Series(sf)
    nf = nm.mul(size_factors.reindex(counts.columns).astype(float), axis=1)
    return nf.div(np.exp(np.log(nf).mean(axis=1)), axis=0)


def _require_scale(c) -> float:
    if c is None:
        raise TypeError("the scale c is required: a within-gene or short-molecule readout has no "
                        "calibrated default (gene counts use c = 1)")
    c = float(c)
    if not math.isfinite(c) or c < 0:
        raise ValueError(f"the scale c must be a finite number >= 0, not {c!r}")
    return c


def ratio_log2_correction(curve: LengthCurve, numerator_length_nt, denominator_length_nt, *, c
                          ) -> np.ndarray:
    """log2 correction to ADD to an observed log2(numerator / denominator) ratio.

    For two molecule classes of lengths L1 and L2 the library over-counts the
    ratio by 2^(c (f(L1) - f(L2))), so the correction is -c (f(L1) - f(L2)).
    ``c`` has NO default: within a gene the between-gene curve over-states the
    length effect (c near 0.5 on one reporter's isoforms), and below the
    panel's shortest gene the curve is held at its edge value.
    """
    c = _require_scale(c)
    f1 = curve.log2_bias(np.atleast_1d(np.asarray(numerator_length_nt, float)))
    f2 = curve.log2_bias(np.atleast_1d(np.asarray(denominator_length_nt, float)))
    return -c * (f1 - f2)


def ratio_corrections(curves: Mapping[str, LengthCurve], numerator_length_nt: float,
                      denominator_length_nt: float, *, c) -> pd.DataFrame:
    """Per-library correction of a two-class ratio; every row records the c used."""
    c = _require_scale(c)
    rows = []
    for lib, cv in curves.items():
        f1 = float(cv.log2_bias([numerator_length_nt])[0])
        f2 = float(cv.log2_bias([denominator_length_nt])[0])
        corr = -c * (f1 - f2)
        rows.append(dict(library=lib, numerator_length_nt=float(numerator_length_nt),
                         denominator_length_nt=float(denominator_length_nt), scale_c=c,
                         log2_bias_numerator=f1, log2_bias_denominator=f2,
                         log2_correction=corr, factor=2.0 ** corr,
                         numerator_clamped=bool(cv.clamped([numerator_length_nt])[0]),
                         denominator_clamped=bool(cv.clamped([denominator_length_nt])[0])))
    return pd.DataFrame(rows)


# ═════════════════════════════════════════════════════════════════════════════
# Certification and evaluation
# ═════════════════════════════════════════════════════════════════════════════

def _token(v) -> str:
    return " ".join(str(v).strip().lower().split())


def certify_libraries(libraries: Sequence[str], library_meta: Optional[pd.DataFrame],
                      scope: Optional[Mapping[str, Sequence[str]]]) -> pd.DataFrame:
    """Certified = the library's declared background, medium and condition all
    fall inside the panel's declared scope. Anything else, including missing
    metadata or a panel without a scope, is NOT certified, with the reason.
    """
    meta = None
    if library_meta is not None:
        meta = library_meta.set_index("library") if "library" in library_meta.columns else library_meta
    rows = []
    for lib in libraries:
        reasons = []
        if not scope:
            reasons.append("the panel declares no certification scope")
        else:
            for key in SCOPE_KEYS:
                allowed = {_token(a) for a in scope.get(key, ())}
                if not allowed:
                    continue
                val = None
                if meta is not None and lib in meta.index and key in meta.columns:
                    v = meta.at[lib, key]
                    val = None if (v is None or (isinstance(v, float) and math.isnan(v))
                                   or str(v).strip() == "") else str(v)
                if val is None:
                    reasons.append(f"{key} not declared")
                elif _token(val) not in allowed:
                    reasons.append(f"{key} '{val}' outside scope {sorted(scope.get(key, ()))}")
        rows.append(dict(library=lib, certified=not reasons,
                         certified_reason="; ".join(reasons) if reasons else "within the panel's scope"))
    return pd.DataFrame(rows)


def replicate_slopes(log2_values: pd.DataFrame, groups: Mapping[str, str], x_log10: pd.Series,
                     genes: Iterable[str]) -> pd.Series:
    """Each library's gene-length slope against the mean of its OTHER replicates.

    log2 per tenfold length, over ``genes`` (use genes that did not set the
    curve). Libraries without a replicate get NaN.
    """
    genes = list(genes)
    x = x_log10.loc[genes].to_numpy(float)
    grp = pd.Series(groups)
    out = {}
    for lib in log2_values.columns:
        g = grp.get(lib)
        others = [l for l in grp.index[grp == g] if l != lib and l in log2_values.columns]
        if g is None or not others:
            out[lib] = float("nan")
            continue
        d = log2_values.loc[genes, lib] - log2_values.loc[genes, others].mean(axis=1)
        out[lib] = float(np.polyfit(x, d.to_numpy(float), 1)[0])
    return pd.Series(out, name="rep_slope")


def replicate_spread(log2_values: pd.DataFrame, groups: Mapping[str, str],
                     genes: Iterable[str]) -> pd.Series:
    """Per replicate group: median over ``genes`` of the SD (ddof 1) of log2 values."""
    genes = list(genes)
    grp = pd.Series(groups)
    out = {}
    for g, libs in grp.groupby(grp):
        libs = [l for l in libs.index if l in log2_values.columns]
        if len(libs) < 2:
            continue
        out[g] = float(np.median(log2_values.loc[genes, libs].std(axis=1, ddof=1).to_numpy()))
    return pd.Series(out, name="rep_sd")
