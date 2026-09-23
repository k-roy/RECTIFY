#!/usr/bin/env python3
"""rectify cdna-length-correct — per-library length-bias correction of ONT PCR-cDNA gene counts.

PCR-cDNA libraries count long genes differently from one another, following
each library's read length (the bias belongs to the library preparation and
can change sign between replicates). This command measures the bias of each
library on a panel of stably expressed genes against a reference abundance
(by default the direct RNA cohort mean), fits a natural spline in log10 gene
length, holds it flat outside the panel's length range, and divides it out.
Algorithm and limits: ``rectify.core.cdna.length_bias``.

Modes:
  genes   CDS-defined gene table (stop codons + lengths) from a GFF3
  count   gene counts from cdna-analyze clusters.tsv (or 3'-end tables, e.g. direct RNA)
  panel   stable-gene panel + reference CPM from reference and cDNA counts
  fit     per-library curves: params, QC, bias factors, corrected CPM, DESeq2 factors
  apply   stored curves applied to another count table of the same libraries
  ratio   correction of a two-class ratio (isoforms, fragments); --scale is REQUIRED

With --Scer the bundled S. cerevisiae gene table, panel and direct RNA
reference are used unless overridden. The bundled panel is certified only for
W303 anchor-away strains in YPD with rapamycin; every other library is flagged.
"""
from __future__ import annotations

import hashlib
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from rectify.core.cdna import length_bias as lb

log = logging.getLogger("cdna-length-correct")


# ── small I/O helpers ─────────────────────────────────────────────────────────

def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_tsv(df: pd.DataFrame, path: Path, index: bool = True, float_format: Optional[str] = None):
    df.to_csv(path, sep="\t", index=index, float_format=float_format)
    log.info("wrote %s", path)


def _write_json(obj: dict, path: Path):
    path.write_text(json.dumps(obj, indent=2, default=str) + "\n")
    log.info("wrote %s", path)


def _provenance(mode: str, args, inputs: Dict[str, Optional[str]], **extra) -> dict:
    from rectify import __version__
    return {
        "command": f"rectify cdna-length-correct {mode}",
        "rectify_version": __version__,
        "created": datetime.now(timezone.utc).isoformat(),
        "argv": sys.argv[1:],
        "inputs": {k: ({"path": str(v), "sha256": _sha256(v)} if v else None) for k, v in inputs.items()},
        **extra,
    }


def _read_matrix(path) -> pd.DataFrame:
    m = pd.read_csv(path, sep="\t", index_col=0)
    m.index = m.index.astype(str)
    m.index.name = "gene"
    return m


def _read_gene_list(path) -> List[str]:
    """A panel: a TSV with a 'gene' column (rows with panel == False dropped) or one gene per line."""
    text = Path(path).read_text().splitlines()
    first = text[0].split("\t") if text else []
    if "gene" in first:
        df = pd.read_csv(path, sep="\t", dtype={"gene": str})
        if "panel" in df.columns:
            df = df[df["panel"].astype(str).str.lower().isin(["true", "1"])]
        return df["gene"].tolist()
    return [ln.strip() for ln in text if ln.strip() and not ln.startswith("#")]


def _read_reference(path) -> pd.Series:
    df = pd.read_csv(path, sep="\t", dtype={"gene": str})
    col = next((c for c in ("cpm", "ref_cpm", "drs_cpm") if c in df.columns), None)
    if "gene" not in df.columns or col is None:
        raise ValueError(f"reference {path} needs a 'gene' column and one of cpm / ref_cpm / drs_cpm")
    return df.set_index("gene")[col].astype(float).rename("ref_cpm")


def _split_list(value: Optional[str]) -> List[str]:
    """'A,B,C' or '@file' (one per line)."""
    if not value:
        return []
    if value.startswith("@"):
        return [ln.strip() for ln in Path(value[1:]).read_text().splitlines() if ln.strip()]
    return [v.strip() for v in value.split(",") if v.strip()]


def _parse_region(spec: str):
    try:
        contig, rng = spec.rsplit(":", 1)
        a, b = rng.replace(",", "").split("-")
        return contig, int(a), int(b)
    except ValueError:
        raise ValueError(f"--exclude-region {spec!r}: expected CONTIG:START-END (0-based, half-open)")


def _bundle(args) -> Optional[Dict[str, Path]]:
    org = getattr(args, "organism", None)
    if not org:
        return None
    from rectify.data import get_bundled_cdna_length_bias
    paths = get_bundled_cdna_length_bias(org)
    if paths is None:
        raise SystemExit(f"ERROR: no bundled cDNA length-bias calibration for organism {org!r}; "
                         "supply --gene-table, --panel and a reference")
    return paths


def _gene_table(args, bundle) -> Tuple[pd.DataFrame, str]:
    """--gene-table > --gff (built) > bundled. Returns (table, description)."""
    if getattr(args, "gene_table", None):
        return lb.load_gene_table(args.gene_table), str(args.gene_table)
    if getattr(args, "gff", None):
        tm = pd.read_csv(args.transcript_models, sep="\t") if getattr(args, "transcript_models", None) else None
        gt = lb.gene_table_from_gff(
            args.gff, transcript_models=tm,
            require_transcript_model=bool(getattr(args, "require_transcript_model", False)),
            exclude_genes=_split_list(getattr(args, "exclude_genes", None)),
            exclude_regions=[_parse_region(r) for r in (getattr(args, "exclude_region", None) or [])],
            drop_mito=not getattr(args, "keep_mito", False))
        return gt, f"built from {args.gff}"
    if bundle is not None:
        return lb.load_gene_table(bundle["gene_table"]), str(bundle["gene_table"])
    raise SystemExit("ERROR: no gene definition; use --gene-table, --gff or --Scer")


# ── modes ─────────────────────────────────────────────────────────────────────

def run_genes(args) -> int:
    gt, _ = _gene_table(args, None)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    _write_tsv(gt, out, index=False)
    log.info("%d genes (%d with a transcript length)", len(gt), int(gt.tx_len.notna().sum()))
    return 0


def _library_names(paths: List[Path], names: Optional[str]) -> List[str]:
    if names:
        out = _split_list(names)
        if len(out) != len(paths):
            raise SystemExit(f"ERROR: --names gives {len(out)} names for {len(paths)} inputs")
        return out
    out = []
    for p in paths:
        d = p.parent
        out.append(d.parent.name if d.name in ("analyze", "cdna_analyze") else d.name)
    if len(set(out)) != len(out):
        raise SystemExit("ERROR: library names derived from the paths are not unique; pass --names or --manifest")
    return out


def run_count(args) -> int:
    bundle = _bundle(args)
    gt, gt_desc = _gene_table(args, bundle)
    index = lb.StopIndex(gt)
    if args.manifest:
        man = pd.read_csv(args.manifest, sep="\t", dtype=str)
        paths, names = [Path(p) for p in man["path"]], man["library"].tolist()
    else:
        paths = [Path(p) for p in args.inputs]
        names = _library_names(paths, args.names)
    if not paths:
        raise SystemExit("ERROR: no inputs")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    counts, stats = {}, []
    for name, path in zip(names, paths):
        if args.input_format == "clusters":
            c, st = lb.count_cdna_clusters(path, index, min_xf=args.min_xf, upstream=args.upstream,
                                           downstream=args.downstream)
        else:
            usecols = ["chrom", "strand", args.position_column]
            head = pd.read_csv(path, sep="\t", nrows=0)
            if "weight" in head.columns:
                usecols.append("weight")
            ends = pd.read_csv(path, sep="\t", usecols=usecols, dtype={"chrom": str, "strand": str})
            ends = ends.rename(columns={args.position_column: "pos"})
            c, st0 = lb.count_three_prime_ends(ends, index, upstream=args.upstream, downstream=args.downstream)
            st = {"n_molecules": st0["n_input"], "n_assigned": int(round(st0["weight_assigned"])),
                  "assigned_share": st0["assigned_share"], "n_orient_unknown": st0["n_bad_strand"],
                  "contigs_without_genes": st0["contigs_without_genes"]}
            if not st["assigned_share"] >= lb.MIN_ASSIGNED_SHARE:
                raise SystemExit(f"ERROR: {name}: only {st['assigned_share']:.1%} of 3' ends reached a gene; "
                                 "suspect a contig-naming mismatch with the gene table")
        counts[name] = c
        stats.append(dict(library=name, path=str(path), **{k: v for k, v in st.items()
                                                            if k != "contigs_without_genes"},
                          contigs_without_genes=",".join(st["contigs_without_genes"])))
        log.info("%s: %s molecules, %.1f%% assigned", name, f"{st['n_molecules']:,}", 100 * st["assigned_share"])
    M = pd.DataFrame(counts).reindex(gt.gene).fillna(0).astype(np.int64)
    M.index.name = "gene"
    _write_tsv(M, out / "counts.tsv")
    _write_tsv(pd.DataFrame(stats), out / "count_stats.tsv", index=False)
    if getattr(args, "gff", None) and not getattr(args, "gene_table", None):
        _write_tsv(gt, out / "gene_table.tsv", index=False)
    _write_json(_provenance("count", args, {"gene_table": args.gene_table or (bundle or {}).get("gene_table")},
                            gene_table=gt_desc, input_format=args.input_format,
                            rule={"min_xf": args.min_xf if args.input_format == "clusters" else None,
                                  "upstream_nt": args.upstream, "downstream_nt": args.downstream,
                                  "one_row_is_one_molecule": args.input_format == "clusters"},
                            libraries=dict(zip(names, map(str, paths)))),
                out / "count_provenance.json")
    return 0


def run_panel(args) -> int:
    bundle = _bundle(args)
    gt, gt_desc = _gene_table(args, bundle)
    C = _read_matrix(args.cdna_counts)
    R = _read_matrix(args.reference_counts)
    grp = pd.read_csv(args.reference_groups, sep="\t", dtype=str)
    groups = dict(zip(grp["library"], grp["group"]))
    missing = [l for l in R.columns if l not in groups]
    if missing:
        raise SystemExit(f"ERROR: reference librar{'y' if len(missing) == 1 else 'ies'} without a group: {missing[:5]}")
    genes = pd.Index(gt.gene)
    P = lb.build_panel(C, R, groups, genes=genes, min_mean=args.min_mean, sd_max=args.sd_max,
                       strict_sd_max=args.strict_sd_max)
    universe = C.index.intersection(genes)
    ref = lb.reference_from_counts(R, universe)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    lengths = gt.set_index("gene")[["cds_len", "tx_len"]]
    _write_tsv(P.join(lengths), out / "panel.tsv")
    _write_tsv(pd.DataFrame({"gene": ref.index, "cpm": ref.values}), out / "reference_cpm.tsv", index=False,
               float_format="%.10g")
    _write_json(_provenance("panel", args, {"cdna_counts": args.cdna_counts, "reference_counts": args.reference_counts,
                                            "reference_groups": args.reference_groups},
                            gene_table=gt_desc, n_considered=int(len(P)), n_panel=int(P.panel.sum()),
                            n_strict=int(P.strict.sum()), n_reference_groups=int(len(set(groups.values()))),
                            criteria={"min_mean": args.min_mean, "sd_max": args.sd_max,
                                      "strict_sd_max": args.strict_sd_max,
                                      "considered_min_mean": lb.SLOPE_SET_MIN_MEAN}),
                out / "panel_provenance.json")
    log.info("panel: %d of %d considered genes (strict %d)", int(P.panel.sum()), len(P), int(P.strict.sum()))
    return 0


def _scope(args, bundle) -> Optional[dict]:
    if getattr(args, "panel_scope", None):
        return json.loads(Path(args.panel_scope).read_text())
    if bundle is not None and not getattr(args, "panel", None) and not getattr(args, "panel_high_count", None):
        return json.loads(Path(bundle["provenance"]).read_text()).get("certified_scope")
    return None


def _covariates(args) -> Optional[pd.DataFrame]:
    if not getattr(args, "covariates", None):
        return None
    cov = pd.read_csv(args.covariates, sep="\t", index_col=0)
    cov.index = cov.index.astype(str)
    cols = _split_list(args.covariate_columns) if getattr(args, "covariate_columns", None) else list(cov.columns)
    return cov[cols].astype(float)


def _write_correction(out: Path, C: pd.DataFrame, B: pd.DataFrame, Ccor: pd.DataFrame, deseq2: bool,
                      size_factors: Optional[pd.Series] = None):
    _write_tsv(B, out / "bias_factors.tsv.gz", float_format="%.8g")
    _write_tsv(Ccor, out / "corrected_cpm.tsv.gz", float_format="%.6g")
    if deseq2:
        nf = lb.deseq2_normalization_factors(C.loc[B.index, B.columns], B, size_factors)
        _write_tsv(nf, out / "deseq2_normalization_factors.tsv.gz", float_format="%.8g")


def run_fit(args) -> int:
    bundle = _bundle(args)
    gt, gt_desc = _gene_table(args, bundle)
    lengths, linfo = lb.resolve_lengths(gt, args.length_type)
    lsha = lb.lengths_fingerprint(lengths)
    C = _read_matrix(args.counts)
    genes = C.index.intersection(lengths.dropna().index)
    n_nolen = int(len(C.index) - len(genes))
    if n_nolen:
        log.warning("%d genes of the count table have no length in the gene table and are left out", n_nolen)

    # reference
    if args.reference_self:
        ref = lb.reference_from_counts(C.loc[genes])
        ref_kind = "self: cohort mean of the fitted libraries (fold changes valid; absolute scale not)"
        ref_path = None
    elif args.reference_counts:
        ref = lb.reference_from_counts(_read_matrix(args.reference_counts), genes)
        ref_kind, ref_path = "cohort mean of reference counts", args.reference_counts
    elif args.reference:
        ref, ref_kind, ref_path = _read_reference(args.reference), "reference CPM table", args.reference
    elif bundle is not None:
        ref, ref_kind, ref_path = _read_reference(bundle["reference"]), "bundled direct RNA cohort mean", bundle["reference"]
    else:
        raise SystemExit("ERROR: no reference; use --reference, --reference-counts, --reference-self or --Scer")

    # panel
    if args.panel_high_count is not None:
        m = C.loc[genes].mean(axis=1)
        panel = m.index[m >= args.panel_high_count].tolist()
        panel_kind, panel_path = f"every gene with mean count >= {args.panel_high_count} (no stability information)", None
    elif args.panel:
        panel, panel_kind, panel_path = _read_gene_list(args.panel), "user panel", args.panel
    elif bundle is not None:
        panel, panel_kind, panel_path = _read_gene_list(bundle["panel"]), "bundled stable-gene panel", bundle["panel"]
    else:
        raise SystemExit("ERROR: no panel; use --panel, --panel-high-count or --Scer")

    cov = _covariates(args)
    fr = lb.fit_length_bias(C.loc[genes], ref, panel, lengths, model=args.model, n_knots=args.knots,
                            length_type=args.length_type, covariates=cov)
    meta = pd.read_csv(args.library_meta, sep="\t", dtype=str) if args.library_meta else None
    scope = _scope(args, bundle)
    cert = lb.certify_libraries(list(C.columns), meta, scope)
    params = fr.params.assign(lengths_sha256=lsha, panel_sha256=lb.genes_fingerprint(panel),
                              n_panel=len(panel), reference=ref_kind)
    qc = fr.qc.merge(cert, on="library")
    for i, lib in enumerate(C.columns):
        cv = fr.curves[lib]
        for j, name in enumerate(cv.cov_names):
            qc.loc[qc.library == lib, f"cov_{name}_log2_per_sd"] = float(cv.coef[len(cv.coef) - len(cv.cov_names) + j])
    rep = None
    if meta is not None and "group" in meta.columns:
        groups = dict(zip(meta["library"], meta["group"]))
        m = C.loc[genes].mean(axis=1)
        eval_genes = [g for g in genes if g not in set(panel) and m[g] >= lb.SLOPE_SET_MIN_MEAN]
        x = np.log10(fr.lengths_nt)
        before = np.log2(lb.cpm(C.loc[genes]) + lb.PSEUDO_CPM)
        after = np.log2(fr.corrected_cpm + lb.PSEUDO_CPM)
        sb = lb.replicate_slopes(before, groups, x, eval_genes)
        sa = lb.replicate_slopes(after, groups, x, eval_genes)
        qc["rep_slope_before"] = qc.library.map(sb)
        qc["rep_slope_after"] = qc.library.map(sa)
        rep = pd.DataFrame({"sd_before": lb.replicate_spread(before, groups, eval_genes),
                            "sd_after": lb.replicate_spread(after, groups, eval_genes)})
        rep.index.name = "group"
        rep["n_eval_genes"] = len(eval_genes)

    cohort = None
    if not args.reference_self:
        cohort = lb.cohort_curve_summary(C.loc[genes], ref, panel, lengths, model=args.model, n_knots=args.knots)
        log.info("cohort-mean curve against the reference (shared by every library; cancels in fold changes): "
                 "%+.2f log2 per tenfold length, range %.2f log2", cohort["slope_equiv_log2_per_tenfold"],
                 cohort["range_log2"])

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    _write_tsv(params, out / "params.tsv", index=False)
    _write_tsv(qc, out / "qc.tsv", index=False)
    if rep is not None:
        _write_tsv(rep, out / "replicate_eval.tsv")
    _write_correction(out, C, fr.factors, fr.corrected_cpm, not args.no_deseq2)
    _write_json(_provenance(
        "fit", args, {"counts": args.counts, "reference": ref_path, "panel": panel_path,
                      "library_meta": args.library_meta, "covariates": getattr(args, "covariates", None)},
        gene_table=gt_desc, lengths=linfo, lengths_sha256=lsha, reference_kind=ref_kind, panel_kind=panel_kind,
        n_panel=len(panel), n_panel_used=int(fr.qc.n_panel_used.min()), n_genes=int(len(genes)),
        n_genes_without_length=n_nolen, model=args.model, knots=args.knots, certified_scope=scope,
        cohort_curve=cohort,
        scale_c=1.0, note="gene counts use the curve at full strength (c = 1); isoform or fragment "
                          "ratios need `ratio --scale`"),
        out / "fit_provenance.json")
    n_cert = int(cert.certified.sum())
    log.info("fitted %d libraries on %d panel genes; %d certified, %d flagged", len(C.columns),
             int(fr.qc.n_panel_used.min()), n_cert, len(C.columns) - n_cert)
    return 0


def run_apply(args) -> int:
    bundle = _bundle(args)
    params = lb.read_params(args.params)
    lt = set(params["length_type"].astype(str))
    if len(lt) != 1:
        raise SystemExit(f"ERROR: params mix length types {sorted(lt)}")
    gt, gt_desc = _gene_table(args, bundle)
    lengths, linfo = lb.resolve_lengths(gt, lt.pop())
    C = _read_matrix(args.counts)
    libs = [l for l in args.libraries.split(",")] if args.libraries else list(C.columns)
    try:
        B, Ccor, qc = lb.apply_length_curves(C[libs], params, lengths, covariates=_covariates(args))
    except ValueError as e:
        raise SystemExit(f"ERROR: {e}")
    sf = None
    if args.size_factors:
        s = pd.read_csv(args.size_factors, sep="\t", dtype={"library": str})
        sf = s.set_index("library")["size_factor"].astype(float)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    _write_tsv(qc, out / "apply_qc.tsv", index=False)
    _write_correction(out, C[libs], B, Ccor, not args.no_deseq2, sf)
    _write_json(_provenance("apply", args, {"params": args.params, "counts": args.counts,
                                            "size_factors": args.size_factors},
                            gene_table=gt_desc, lengths=linfo, scale_c=1.0),
                out / "apply_provenance.json")
    return 0


def run_ratio(args) -> int:
    params = lb.read_params(args.params)
    curves = lb.curves_from_params(params)
    if args.libraries:
        want = args.libraries.split(",")
        curves = {l: curves[l] for l in want}
    try:
        tab = lb.ratio_corrections(curves, args.numerator_length, args.denominator_length, c=args.scale)
    except (TypeError, ValueError) as e:
        raise SystemExit(f"ERROR: {e}")
    if args.out in (None, "-"):
        tab.to_csv(sys.stdout, sep="\t", index=False)
    else:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        _write_tsv(tab, Path(args.out), index=False)
    n_cl = int((tab.numerator_clamped | tab.denominator_clamped).sum())
    if n_cl:
        log.warning("%d of %d libraries: a length lies outside the fitted range, so the curve's edge value "
                    "was used (an extrapolation whose size depends on c)", n_cl, len(tab))
    return 0


def run(args) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    mode = getattr(args, "cdna_length_mode", None)
    handlers = {"genes": run_genes, "count": run_count, "panel": run_panel, "fit": run_fit,
                "apply": run_apply, "ratio": run_ratio}
    if mode not in handlers:
        print("usage: rectify cdna-length-correct {genes,count,panel,fit,apply,ratio} ...", file=sys.stderr)
        return 1
    try:
        return handlers[mode](args)
    except ValueError as e:
        log.error("%s", e)
        return 1


# ── parser ────────────────────────────────────────────────────────────────────

def _add_gene_args(p, organism=True):
    g = p.add_argument_group("gene definitions (one of)")
    g.add_argument("--gene-table", default=None,
                   help="Gene table TSV (gene, chrom, strand, cds_start0, cds_end, cds_len[, tx_len]), "
                        "e.g. from `cdna-length-correct genes`")
    g.add_argument("--gff", default=None, help="Build the gene table from this GFF3 (CDS-defined genes)")
    g.add_argument("--transcript-models", default=None,
                   help="With --gff: TSV gene, tx_lo, tx_hi (1-based inclusive) giving spliced transcript lengths")
    g.add_argument("--require-transcript-model", action="store_true",
                   help="With --gff: keep only genes that have a transcript model")
    g.add_argument("--exclude-genes", default=None,
                   help="With --gff: genes removed from the candidate-stop set (comma list or @file)")
    g.add_argument("--exclude-region", action="append", default=None, metavar="CONTIG:START-END",
                   help="With --gff: drop genes whose CDS overlaps this 0-based half-open region (repeatable), "
                        "e.g. the rDNA locus")
    g.add_argument("--keep-mito", action="store_true", help="With --gff: keep mitochondrial genes")
    if organism:
        from rectify.data import add_organism_args
        add_organism_args(g)


def create_cdna_length_correct_parser(subparsers):
    p = subparsers.add_parser(
        "cdna-length-correct",
        help="Per-library length-bias correction of ONT PCR-cDNA gene counts (stable-gene panel vs direct RNA)",
        description=__doc__, formatter_class=__import__("argparse").RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cdna_length_mode")

    g = sub.add_parser("genes", help="Build the CDS-defined gene table from a GFF3")
    _add_gene_args(g, organism=False)
    g.add_argument("-o", "--out", required=True, help="Output gene table TSV")

    c = sub.add_parser("count", help="Gene counts from cdna-analyze clusters.tsv (or 3'-end tables)")
    c.add_argument("inputs", nargs="*", help="clusters.tsv[.gz] files, one per library")
    c.add_argument("--names", default=None, help="Library names, comma list in input order (default: from the paths)")
    c.add_argument("--manifest", default=None, help="TSV with columns library, path (instead of positional inputs)")
    c.add_argument("--input-format", choices=("clusters", "three-prime"), default="clusters",
                   help="clusters: cdna-analyze clusters.tsv (one row = one molecule; XF filter). three-prime: a "
                        "table with chrom, strand (RNA strand), a 3'-end column and optional weight, e.g. a direct "
                        "RNA corrected_reads.tsv; it applies no alignment filter, so keep primary MAPQ >= 20 "
                        "alignments first to match the bundled reference")
    c.add_argument("--position-column", default="corrected_3prime",
                   help="three-prime input: the 0-based 3'-end column (default: corrected_3prime)")
    c.add_argument("--min-xf", type=int, default=lb.DEFAULT_MIN_XF,
                   help="Keep molecules with XF >= this (both read types; default 1)")
    c.add_argument("--upstream", type=int, default=lb.ASSIGN_UPSTREAM_NT,
                   help="A 3' end may lie this far upstream of (inside) the stop codon (default 30)")
    c.add_argument("--downstream", type=int, default=lb.ASSIGN_DOWNSTREAM_NT,
                   help="... and this far downstream of it (default 500)")
    _add_gene_args(c)
    c.add_argument("-o", "--out", required=True, help="Output directory (counts.tsv, count_stats.tsv)")

    pn = sub.add_parser("panel", help="Stable-gene panel and reference CPM from reference and cDNA counts")
    pn.add_argument("--cdna-counts", required=True, help="cDNA counts (genes x libraries)")
    pn.add_argument("--reference-counts", required=True,
                    help="Reference (e.g. direct RNA) counts made with the same rule (genes x libraries)")
    pn.add_argument("--reference-groups", required=True,
                    help="TSV library, group: the condition group of each reference library; stability is "
                         "measured across the group means")
    pn.add_argument("--min-mean", type=float, default=lb.PANEL_MIN_MEAN,
                    help="Minimum mean count per library in BOTH assays (default 100)")
    pn.add_argument("--sd-max", type=float, default=lb.PANEL_SD_MAX,
                    help="Maximum SD (log2) of log2 CPM across reference group means (default 0.5)")
    pn.add_argument("--strict-sd-max", type=float, default=lb.PANEL_STRICT_SD_MAX,
                    help="SD cut-off of the 'strict' flag (default 0.35)")
    _add_gene_args(pn)
    pn.add_argument("-o", "--out", required=True, help="Output directory (panel.tsv, reference_cpm.tsv)")

    f = sub.add_parser("fit", help="Fit each library's curve on the panel and correct the table")
    f.add_argument("--counts", required=True, help="cDNA counts (genes x libraries), e.g. from `count`")
    r = f.add_argument_group("reference (one of; --Scer = bundled direct RNA cohort mean)")
    r.add_argument("--reference", default=None, help="TSV gene, cpm")
    r.add_argument("--reference-counts", default=None, help="Reference counts (genes x libraries); cohort-mean CPM")
    r.add_argument("--reference-self", action="store_true",
                   help="Use the cohort mean of the fitted libraries: fold changes and DESeq2 factors are "
                        "unaffected by the reference choice, but absolute between-gene values keep the "
                        "cohort's shared bias")
    q = f.add_argument_group("panel (one of; --Scer = bundled 405-gene panel)")
    q.add_argument("--panel", default=None, help="Panel genes: TSV with a 'gene' column, or one gene per line")
    q.add_argument("--panel-high-count", type=float, default=None, metavar="N",
                   help="Use every gene with mean count >= N as the panel (no stability information; libraries "
                        "are not certified)")
    q.add_argument("--panel-scope", default=None,
                   help="JSON {background: [...], medium: [...], condition: [...]}: where a user panel is "
                        "certified")
    f.add_argument("--library-meta", default=None,
                   help="TSV library[, background, medium, condition, group]: certification against the panel's "
                        "scope; 'group' (replicates) adds before/after replicate slopes")
    f.add_argument("--length-type", choices=lb.LENGTH_TYPES, default=lb.DEFAULT_LENGTH_TYPE,
                   help="Gene length the curve uses: tx (transcript, default) or cds")
    f.add_argument("--model", choices=("ns", "linear"), default=lb.DEFAULT_MODEL,
                   help="Curve: natural cubic spline (default) or linear")
    f.add_argument("--knots", type=int, default=lb.DEFAULT_N_KNOTS,
                   help="Spline knots at panel length quantiles (default 4)")
    f.add_argument("--covariates", default=None,
                   help="OPTIONAL per-gene covariates TSV (gene + numeric columns), added to the curve; off by "
                        "default (in the calibration cohort, GC and UpA/CpA content changed held-out error by < 0.5%%)")
    f.add_argument("--covariate-columns", default=None, help="Columns of --covariates to use (default: all)")
    f.add_argument("--no-deseq2", action="store_true", help="Do not write DESeq2 normalization factors")
    _add_gene_args(f)
    f.add_argument("-o", "--out", required=True, help="Output directory")

    a = sub.add_parser("apply", help="Apply stored curves to another count table of the same libraries")
    a.add_argument("--params", required=True, help="params.tsv from `fit`")
    a.add_argument("--counts", required=True, help="Counts (genes x libraries) made with the same rule")
    a.add_argument("--libraries", default=None, help="Comma list of libraries to correct (default: all columns)")
    a.add_argument("--covariates", default=None, help="The covariates table used at fit time, if any")
    a.add_argument("--covariate-columns", default=None, help="Columns of --covariates (as at fit time)")
    a.add_argument("--size-factors", default=None, help="TSV library, size_factor for the DESeq2 factors")
    a.add_argument("--no-deseq2", action="store_true", help="Do not write DESeq2 normalization factors")
    _add_gene_args(a)
    a.add_argument("-o", "--out", required=True, help="Output directory")

    rt = sub.add_parser("ratio", help="Correction of a two-class ratio; --scale is required")
    rt.add_argument("--params", required=True, help="params.tsv from `fit`")
    rt.add_argument("--numerator-length", type=float, required=True, help="Length (nt) of the numerator class")
    rt.add_argument("--denominator-length", type=float, required=True, help="Length (nt) of the denominator class")
    rt.add_argument("--scale", type=float, required=True, metavar="C",
                    help="REQUIRED scale c: correction = -c (f(L1) - f(L2)). Gene counts use c = 1; within a "
                         "gene the calibrated value is uncertain (about 0.5 on one reporter's isoforms)")
    rt.add_argument("--libraries", default=None, help="Comma list of libraries (default: all in params)")
    rt.add_argument("-o", "--out", default=None, help="Output TSV (default: stdout)")
    return p
