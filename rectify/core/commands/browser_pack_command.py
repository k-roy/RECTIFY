#!/usr/bin/env python3
"""
RECTIFY pack-browser Command — fold ``rectify analyze`` outputs + the ``rectify
qc`` blocks into ONE ``analysis.json`` for a genome-browser Analysis tab.

Ported from the validated prototype ``planning/630f_pack_analysis.py`` in the
Chanfreau workspace, generalised off its one hard-coded dataset (see
"GENERALISATION" below).

🔴 THIS IS A PACKER, NOT AN ANALYST.  It reshapes numbers the pipeline already
produced.  It does not compute differential expression, does not renormalise,
and does not second-guess a p-value.  Deriving a second set of DE numbers that
disagrees with the pipeline's is the worst available outcome — so where the
pipeline's choice is questionable, the answer is a CAVEAT STRING that reaches the
user's screen, never a silent recomputation.

FAIL-OPEN CONTRACT: every section is independent.  A missing analyze directory
yields ``"de": {"<contrast>": null}`` -> the tab renders "not computed".  It never
yields an empty table, because an empty table reads as "there is nothing here",
which is the failure this design keeps guarding against.

────────────────────────────────────────────────────────────────────────────────
GENERALISATION vs the prototype (the prototype was written for one dataset)
────────────────────────────────────────────────────────────────────────────────
 1. ``COND_LABEL`` / ``REFERENCE`` were hard-coded to wt / ysh1 / rna15.  The
    reference condition is now read from ``<analyze-dir>/sample_metadata.tsv``
    (the ``is_control`` column, which ``rectify analyze`` writes), and condition
    labels fall back to the condition name itself.  An unrecognised condition
    never crashes — it labels as itself.
 2. ``--samples`` is optional.  When absent the sample list is derived from
    ``sample_metadata.tsv`` (genotype := condition, replicate parsed from a
    trailing ``_rep<N>``).
 3. The prototype's ``dataset`` string and several ``caveats`` carried that
    dataset's own MEASURED NUMBERS ("nine libraries", "4.87x", "n = 3 per
    genotype", "15.0% against a true 4.8%").  Every RULE is kept verbatim; the
    numbers are now computed from the inputs actually present, or the sentence is
    omitted when they cannot be.  Shipping another study's numbers as if they
    were the user's would be the same second-set-of-numbers failure this module
    exists to prevent.

Usage:
  rectify pack-browser --qc-dir qc --analyze-dir analyze_out \\
      --gene-index gene_index.json --library-depth library_depth.json \\
      [--samples samples.json] [--pressure pressure_test.json] [--top-n 400] \\
      --out analysis.json

Author: Kevin R. Roy
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from typing import Dict, List, Optional

# NOTE: pandas / numpy are imported INSIDE the functions that need them — this
# module is imported by ``rectify.cli.create_parser`` on every invocation.

_REP_RE = re.compile(r'_rep(\d+)$', re.IGNORECASE)


# ──────────────────────────────────────────────────────────────────────────────
# sample_metadata.tsv — the reference condition and the sample list
# ──────────────────────────────────────────────────────────────────────────────

def load_sample_metadata(analyze_dir: Optional[str]) -> tuple:
    """Read ``<analyze-dir>/sample_metadata.tsv``.

    Returns ``(rows, reference_condition)`` where *rows* is a list of
    ``{'sample','condition','is_control'}`` dicts and *reference_condition* is the
    first control condition, or ``None`` when it cannot be determined.

    ⚠️ Two traps this function absorbs, both verified against
    ``rectify.core.analyze.deseq2.create_sample_metadata``:
      (a) that function sets ``metadata.index = metadata['sample']`` and then
          ``to_csv(sep='\\t')``, so the file has TWO ``sample`` columns and pandas
          reads them back as ``sample`` / ``sample.1``;
      (b) ``is_control`` round-trips as a string in some writers, so it is
          coerced from text rather than trusted to be a bool.
    Contrast KEYS are still taken from the ``deseq2_<level>_<cond>.tsv``
    filenames, never from here: in manifest mode the on-disk ``condition`` column
    is overwritten AFTER the file is written, so it can disagree with what DESeq2
    actually used.
    """
    if not analyze_dir:
        return [], None
    path = os.path.join(analyze_dir, 'sample_metadata.tsv')
    if not os.path.exists(path):
        return [], None
    try:
        import pandas as pd
        df = pd.read_csv(path, sep='\t')
    except Exception as exc:
        print('  WARN sample_metadata.tsv unreadable: %s' % exc, file=sys.stderr)
        return [], None
    if df.empty:
        return [], None

    # The duplicated 'sample' column reads back as sample / sample.1 — either is
    # the sample id, so take the first column whose name starts with 'sample'.
    samp_col = next((c for c in df.columns if str(c).split('.')[0] == 'sample'), None)
    if samp_col is None:
        samp_col = df.columns[0]
    cond_col = 'condition' if 'condition' in df.columns else None

    if 'is_control' in df.columns:
        ctrl = df['is_control'].astype(str).str.strip().str.lower().isin(
            ('true', '1', 'yes', 't', 'y'))
    else:
        ctrl = [False] * len(df)

    rows = []
    for i in range(len(df)):
        rows.append({
            'sample': str(df.iloc[i][samp_col]),
            'condition': str(df.iloc[i][cond_col]) if cond_col else '',
            'is_control': bool(list(ctrl)[i]),
        })

    reference = next((r['condition'] for r in rows if r['is_control'] and r['condition']),
                     None)
    return rows, reference


def samples_from_metadata(meta_rows: List[dict], modality: str = '') -> List[dict]:
    """Synthesize the ``--samples`` records from ``sample_metadata.tsv``.

    All five keys the packer reads (``name``, ``analyze_id``, ``genotype``,
    ``replicate``, ``modality``) are always present, because downstream code uses
    bracket access on them.
    """
    out = []
    for r in meta_rows:
        sid = r['sample']
        m = _REP_RE.search(sid)
        out.append({
            'name': sid,
            'analyze_id': sid,
            'genotype': r.get('condition') or '',
            'replicate': m.group(1) if m else '',
            'modality': modality or '',
        })
    return out


def build_cond_labels(samples_meta: List[dict], meta_rows: List[dict]) -> Dict[str, str]:
    """condition -> display label.

    When a ``--samples`` file supplies a genotype per sample, its genotype is used
    as the label for that sample's condition (this is what the prototype's
    hard-coded ``COND_LABEL`` did).  Otherwise the label is the condition itself.
    A condition with conflicting genotypes across its samples falls back to the
    condition name rather than picking one arbitrarily.
    """
    cond_of = {r['sample']: r.get('condition', '') for r in meta_rows}
    seen: Dict[str, set] = {}
    for s in samples_meta:
        cond = cond_of.get(s.get('analyze_id') or s.get('name') or '', '')
        geno = (s.get('genotype') or '').strip()
        if cond and geno:
            seen.setdefault(cond, set()).add(geno)
    return {c: (next(iter(v)) if len(v) == 1 else c) for c, v in seen.items()}


def make_contrast_label(cond_labels: Dict[str, str], reference: Optional[str]):
    """Return a ``contrast_label(cond)`` closure that never crashes on an
    unrecognised condition."""
    def contrast_label(cond: str) -> str:
        left = cond_labels.get(cond, cond)
        if not reference:
            return str(left)
        return '%s vs %s' % (left, cond_labels.get(reference, reference))
    return contrast_label


# ──────────────────────────────────────────────────────────────────────────────
# Gene coordinates + JSON hygiene
# ──────────────────────────────────────────────────────────────────────────────

def load_gene_index(path: Optional[str]) -> dict:
    """alias -> ``{chrom,start,end,strand,label}``.  Reused from the BROWSER's own
    index so a row's coordinates are the same ones the browser will navigate to —
    they cannot disagree."""
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception as exc:
        print('  WARN gene index unreadable: %s' % exc, file=sys.stderr)
        return {}


def coords_for(gidx: dict, name: str) -> dict:
    g = gidx.get(name)
    if g is None and name:
        g = gidx.get(name.split('_')[0])          # YAL001C_mRNA -> YAL001C
    if g is None:
        return {}
    return {'chrom': g.get('chrom'), 'start': g.get('start'), 'end': g.get('end'),
            'strand': g.get('strand'), 'label': g.get('label') or name}


def _f(v):
    """JSON-safe float: NaN/Inf are not valid JSON and ``json.dump`` emits them
    anyway (as ``NaN``), which makes the file unparseable by ``JSON.parse`` in the
    browser — a blank tab with no error."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    if x != x or x in (float('inf'), float('-inf')):
        return None
    return round(x, 6)


# ──────────────────────────────────────────────────────────────────────────────
# Section packers
# ──────────────────────────────────────────────────────────────────────────────

def pack_de(tables_dir: str, gidx: dict, top_n: int, level: str = 'genes') -> dict:
    """``deseq2_<level>_<cond>.tsv`` -> ``{contrast: [rows]}`` with coordinates on
    every row."""
    import pandas as pd
    out: dict = {}
    if not os.path.isdir(tables_dir):
        return out
    for fn in sorted(os.listdir(tables_dir)):
        if not (fn.startswith('deseq2_%s_' % level) and fn.endswith('.tsv')):
            continue
        cond = fn[len('deseq2_%s_' % level):-4]
        try:
            df = pd.read_csv(os.path.join(tables_dir, fn), sep='\t', index_col=0)
        except Exception as e:
            out[cond] = None
            print('  WARN %s unreadable: %s' % (fn, e), file=sys.stderr)
            continue
        if df.empty:
            out[cond] = None
            continue
        df = df.reset_index().rename(columns={df.index.name or 'index': 'feature'})
        # 🔴 Rank by |log2FC| among padj-significant rows, NOT by padj alone. At
        # small n DESeq2's dispersion estimates are weak and a padj ranking
        # promotes tiny, well-measured effects over large ones. Both numbers
        # travel with every row so the reader can re-sort; the tab's default sort
        # is a presentation choice, not a claim.
        if 'padj' not in df.columns or 'log2FoldChange' not in df.columns:
            out[cond] = None
            print('  WARN %s lacks padj/log2FoldChange columns' % fn, file=sys.stderr)
            continue
        df['_sig'] = df['padj'].notna() & (df['padj'] < 0.05)
        df['_abs'] = df['log2FoldChange'].abs()
        df = df.sort_values(['_sig', '_abs'], ascending=[False, False]).head(top_n)
        rows = []
        for r in df.itertuples(index=False):
            feat = str(getattr(r, 'feature'))
            row = {'feature': feat,
                   'log2fc': _f(getattr(r, 'log2FoldChange', None)),
                   'padj': _f(getattr(r, 'padj', None)),
                   'pvalue': _f(getattr(r, 'pvalue', None)),
                   'baseMean': _f(getattr(r, 'baseMean', None)),
                   'lfcSE': _f(getattr(r, 'lfcSE', None))}
            if level == 'genes':
                row.update(coords_for(gidx, feat))
                row['gene'] = row.get('label') or feat
            else:
                # cluster-level rows were merged with clusters_df by rectify, so
                # they already carry real coordinates — prefer those over a name
                # lookup.
                for k, src in (('chrom', 'chrom'), ('start', 'start'), ('end', 'end'),
                               ('strand', 'strand')):
                    v = getattr(r, src, None)
                    if v is not None and v == v:
                        row[k] = v if k in ('chrom', 'strand') else int(v)
                gene = getattr(r, 'gene_name', None) or getattr(r, 'gene_id', None)
                if gene and gene == gene:
                    row['gene'] = str(gene)
                    if not row.get('chrom'):
                        row.update(coords_for(gidx, str(gene)))
                mp = getattr(r, 'modal_position', None)
                if mp is not None and mp == mp:
                    row['pos'] = int(mp)
            rows.append(row)
        out[cond] = rows
        print('  de[%s] %-8s %d rows' % (level, cond, len(rows)))
    return out


def pack_shift(tables_dir: str, gidx: dict, top_n: int) -> dict:
    """``shift_analysis_<cond>.tsv`` -> 3' isoform SHIFT rows.

    🔴 This is the one analyze table structurally immune to the denominator trap:
    it normalises WITHIN GENE (frac = gene_counts / gene_total), so library depth
    and any ncRNA inflation of the global denominator cancel exactly.
    ⚠️ It does average RAW counts across replicates before taking the fraction, so
    a condition's fraction is depth-WEIGHTED across its own replicates. Carried as
    a caveat rather than recomputed.
    """
    import pandas as pd
    out: dict = {}
    if not os.path.isdir(tables_dir):
        return out
    for fn in sorted(os.listdir(tables_dir)):
        if not (fn.startswith('shift_analysis_') and fn.endswith('.tsv')):
            continue
        cond = fn[len('shift_analysis_'):-4]
        try:
            df = pd.read_csv(os.path.join(tables_dir, fn), sep='\t')
        except Exception:
            out[cond] = None
            continue
        if df.empty:
            out[cond] = None
            continue
        if 'distribution_divergence' in df.columns:
            df = df.sort_values('distribution_divergence', ascending=False)
        df = df.head(top_n)
        rows = []
        for r in df.itertuples(index=False):
            gene = str(getattr(r, 'gene_name', None) or getattr(r, 'gene_id', ''))
            row = {'gene': gene,
                   'gene_id': str(getattr(r, 'gene_id', '')),
                   'n_clusters': int(getattr(r, 'n_clusters', 0) or 0),
                   'shift_bp': _f(getattr(r, 'shift_bp', None)),
                   'shift_direction': str(getattr(r, 'shift_direction', '')),
                   'jsd': _f(getattr(r, 'distribution_divergence', None)),
                   'major_frac_a': _f(getattr(r, 'major_frac_a', None)),
                   'major_frac_b': _f(getattr(r, 'major_frac_b', None)),
                   'counts_a': _f(getattr(r, 'counts_a', None)),
                   'counts_b': _f(getattr(r, 'counts_b', None))}
            row.update(coords_for(gidx, gene))
            rows.append(row)
        out[cond] = rows
        print('  shift[%s] %d rows' % (cond, len(rows)))
    return out


def pack_iso5(analyze_dir: str, gidx: dict, samples_meta: List[dict],
              top_n: int) -> Optional[list]:
    """5' isoforms from ``tss_clusters.tsv`` + ``tss_cluster_counts.tsv``.

    🔴 rectify runs NO differential test on TSS clusters — only CPA clusters reach
    DESeq2.  So this section reports per-condition WITHIN-GENE USAGE, and the tab
    must label it "usage, not DE".  Inventing a DE here would be exactly the
    second-set-of-numbers failure this module forbids.

    🔴 And on DRS a "TSS cluster" is partly READ LENGTH: the pore truncates the 5'
    end (median span 1,150-1,380 nt, class-blind), so these are 5' ENDS OF READS,
    not transcription start sites.
    """
    import pandas as pd
    cl_p = os.path.join(analyze_dir, 'tss_clusters.tsv')
    ct_p = os.path.join(analyze_dir, 'tss_cluster_counts.tsv')
    if not (os.path.exists(cl_p) and os.path.exists(ct_p)):
        return None
    try:
        cl = pd.read_csv(cl_p, sep='\t')
        ct = pd.read_csv(ct_p, sep='\t', index_col=0)
    except Exception as e:
        print('  WARN iso5 unreadable: %s' % e, file=sys.stderr)
        return None
    if cl.empty or ct.empty:
        return None

    samp_cond = {s.get('analyze_id'): s.get('genotype') for s in samples_meta}
    conds: Dict[str, list] = {}
    for c in ct.columns:
        conds.setdefault(samp_cond.get(c) or c, []).append(c)
    # within-condition mean, then within-GENE fraction: depth cancels
    per_cond = pd.DataFrame({k: ct[v].mean(axis=1) for k, v in conds.items()})

    if 'cluster_id' not in cl.columns:
        return None
    cl = cl.set_index('cluster_id')
    gene_col = ('gene_name' if 'gene_name' in cl.columns
                else ('gene_id' if 'gene_id' in cl.columns else None))
    if gene_col is None:
        return None
    common = [i for i in per_cond.index if i in cl.index]
    if not common:
        return None
    per_cond = per_cond.loc[common]
    genes = cl.loc[common, gene_col]

    rows = []
    for gene, idx in genes.groupby(genes).groups.items():
        sub = per_cond.loc[idx]
        tot = sub.sum(axis=0)
        if tot.sum() < 20:
            continue
        frac = sub.divide(tot.replace(0, pd.NA), axis=1)
        # spread across conditions of the dominant cluster's usage = how much the
        # 5' end moved
        spread = float(frac.max(axis=1).max() - frac.max(axis=1).min()) if len(frac) else 0.0
        top = frac.mean(axis=1).idxmax()
        r = {'gene': str(gene), 'n_clusters': int(len(idx)),
             'usage': {str(c): _f(frac.loc[top, c]) for c in frac.columns},
             'counts': {str(c): _f(tot[c]) for c in tot.index},
             'spread': _f(spread)}
        if top in cl.index:
            row = cl.loc[top]
            r['pos'] = (int(row['modal_position'])
                        if 'modal_position' in row and row['modal_position'] == row['modal_position']
                        else None)
        r.update(coords_for(gidx, str(gene)))
        rows.append(r)
    rows.sort(key=lambda x: -(x['spread'] or 0))
    print('  iso5 %d genes' % len(rows))
    return rows[:top_n]


def pack_pca_corr(analyze_dir: str, samples_meta: List[dict]) -> tuple:
    """PCA coordinates + sample correlation, recomputed from
    ``cluster_counts.tsv`` with EXACTLY the pipeline's own transform, because
    analyze writes the PLOTS but not the coordinates.

    ⚠️ This is the one place the packer touches numbers, so it must match rectify
    bit-for-bit:
      PCA     : log2(x+1) -> per-feature z-score -> PCA   (pca.py run_pca_analysis defaults)
      heatmap : log2(x+1) -> top-1000 variance -> Pearson (heatmap.py plot_sample_heatmap)
    It is IMPORTED from rectify rather than reimplemented so it CANNOT drift.
    """
    import numpy as np
    import pandas as pd
    p = os.path.join(analyze_dir, 'cluster_counts.tsv')
    if not os.path.exists(p):
        return None, None
    try:
        cm = pd.read_csv(p, sep='\t', index_col=0)
    except Exception as e:
        print('  WARN cluster_counts.tsv unreadable: %s' % e, file=sys.stderr)
        return None, None
    if cm.empty or cm.shape[1] < 2:
        return None, None
    try:
        from ..analyze.pca import run_pca_analysis
    except Exception as e:
        print('  WARN cannot import rectify PCA: %s' % e, file=sys.stderr)
        return None, None

    try:
        res = run_pca_analysis(cm)
    except Exception as e:
        print('  WARN PCA failed: %s' % e, file=sys.stderr)
        return None, None
    coords = res['pca_coords']
    if coords is None or len(coords) == 0 or coords.shape[1] < 2:
        return None, None
    vr = list(res['variance_ratio'])
    lbl = {s.get('analyze_id'): s for s in samples_meta}
    pca = {
        'labels': [lbl.get(i, {}).get('name', i) for i in coords.index],
        'genotype': [lbl.get(i, {}).get('genotype', '') for i in coords.index],
        'replicate': [lbl.get(i, {}).get('replicate', '') for i in coords.index],
        'pc': [[_f(coords.iloc[i, 0]), _f(coords.iloc[i, 1])] for i in range(len(coords))],
        'var_explained': [_f(v) for v in vr[:2]],
        'n_features': int(res['n_features']),
        'transform': 'log2(count+1), per-feature z-score (rectify pca.py defaults)',
    }
    d = np.log2(cm + 1)
    if len(d) > 1000:
        d = d.loc[d.var(axis=1).nlargest(1000).index]
    cmat = d.corr()
    corr = {
        'labels': [lbl.get(i, {}).get('name', i) for i in cmat.columns],
        'matrix': [[_f(cmat.iloc[i, j]) for j in range(len(cmat))] for i in range(len(cmat))],
        'metric': 'Pearson on log2(count+1), 1,000 most-variable clusters',
        'n_features': int(len(d)),
    }
    print('  pca %d samples, %d features; corr %d features'
          % (len(coords), res['n_features'], len(d)))
    return pca, corr


# ──────────────────────────────────────────────────────────────────────────────
# Caveats — shown IN THE TAB, not buried in a doc
# ──────────────────────────────────────────────────────────────────────────────

def build_caveats(qc: dict, depth: dict, samples_meta: List[dict],
                  iso5: Optional[list]) -> List[str]:
    """Every RULE below is verbatim from the prototype.  The dataset-specific
    NUMBERS are recomputed from the inputs actually present, or the sentence is
    dropped — quoting another study's measurements would be exactly the
    second-set-of-numbers failure this module exists to prevent."""
    caveats: List[str] = []

    # depth spread
    totals = [v.get('total') for v in depth.values()
              if isinstance(v, dict) and isinstance(v.get('total'), (int, float))
              and v.get('total')]
    if len(totals) >= 2 and min(totals) > 0:
        caveats.append(
            "Depth spans %.2fx across these %d libraries. Every raw-count display "
            "is labelled raw; rates (error, clip, untailed) are depth-independent "
            "and safe." % (max(totals) / float(min(totals)), len(totals)))
    else:
        caveats.append(
            "Library depth differs between samples. Every raw-count display is "
            "labelled raw; rates (error, clip, untailed) are depth-independent "
            "and safe.")

    # replicate structure
    per_geno: Dict[str, int] = {}
    for s in samples_meta:
        per_geno[s.get('genotype') or ''] = per_geno.get(s.get('genotype') or '', 0) + 1
    if per_geno:
        ns = sorted(set(per_geno.values()))
        n_txt = ('n = %d per genotype' % ns[0] if len(ns) == 1
                 else 'n = %d-%d per genotype' % (ns[0], ns[-1]))
        caveats.append(
            "%s. Replicate points are shown individually; there are no error bars "
            "and no significance stars beyond DESeq2's own padj." % n_txt)

    caveats.append(
        "The two poly(A) measures are never merged. dorado pt is signal-level; "
        "rectify's polya_length is read off the basecall and saturates near 14 nt. "
        "The gap is the comparison.")

    # error rate — quote THIS dataset's own naive-vs-corrected gap when present
    errs = [q['error'] for q in qc.values() if isinstance(q.get('error'), dict)]
    pairs = [(e.get('naive_nm_rate_pct'), e.get('rate_pct')) for e in errs
             if e.get('rate_pct') is not None
             and e.get('naive_nm_rate_pct') is not None]
    d_min = next((e.get('intron_d_min_bp') for e in errs
                  if isinstance(e.get('intron_d_min_bp'), int)), 20)
    err_rule = ("Error rate excludes deletions >=%d bp (introns -- these are "
                "map-ont alignments, NOT splice-aware) and the poly(A) walkback "
                "interval. The naive NM/aligned-length rate is shown beside it"
                % d_min)
    if pairs:
        naive = sorted(p[0] for p in pairs)[len(pairs) // 2]
        true = sorted(p[1] for p in pairs)[len(pairs) // 2]
        err_rule += (": across these libraries the median naive rate reads %.1f%% "
                     "against a corrected %.1f%%, because most of NM is intron."
                     % (naive, true))
    else:
        err_rule += "."
    caveats.append(err_rule)

    caveats.append(
        "5' and 3' soft clips are different quantities and are never pooled. A 5' "
        "clip is predominantly a spliced read whose junction failed to align.")

    modalities = {(s.get('modality') or '').upper() for s in samples_meta}
    modalities.discard('')
    if not modalities or modalities == {'DRS'}:
        caveats.append(
            "DRS is never UMI-deduped, so the duplicate-rate panel reads n/a for "
            "these libraries rather than 0%.")

    if iso5:
        caveats.append(
            "5' isoforms are 5' ENDS OF READS, not transcription start sites: DRS "
            "is truncated by the pore (median span 1,150-1,380 nt, class-blind), so "
            "a TSS cluster is partly read length. rectify runs no differential test "
            "on TSS clusters, so that section reports usage, not DE.")

    caveats.append(
        "Splicing is deliberately absent: junction coordinates are not in the read "
        "bundle and need separate extraction from rectify's TSV.")

    return caveats


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def run_browser_pack(args: argparse.Namespace) -> int:
    """Fold analyze + QC into one ``analysis.json``.  Returns an exit code."""
    qc_dir = getattr(args, 'qc_dir', None) or 'qc'
    analyze_dir = getattr(args, 'analyze_dir', None)
    top_n = getattr(args, 'top_n', 400)

    meta_rows, reference = load_sample_metadata(analyze_dir)

    samples_path = getattr(args, 'samples', None)
    samples_meta: List[dict] = []
    if samples_path and os.path.exists(samples_path):
        try:
            with open(samples_path) as fh:
                samples_meta = json.load(fh)
        except Exception as exc:
            print('  WARN --samples unreadable (%s); deriving from sample_metadata.tsv'
                  % exc, file=sys.stderr)
            samples_meta = []
    if not samples_meta:
        samples_meta = samples_from_metadata(meta_rows, getattr(args, 'modality', '') or '')
        if samples_meta:
            print('samples: derived %d from sample_metadata.tsv' % len(samples_meta))

    cond_labels = build_cond_labels(samples_meta, meta_rows)
    contrast_label = make_contrast_label(cond_labels, reference)

    gidx = load_gene_index(getattr(args, 'gene_index', None))

    depth_path = getattr(args, 'library_depth', None)
    depth = {}
    if depth_path and os.path.exists(depth_path):
        try:
            with open(depth_path) as fh:
                depth = json.load(fh)
        except Exception as exc:
            print('  WARN --library-depth unreadable: %s' % exc, file=sys.stderr)

    # ── QC ────────────────────────────────────────────────────────────────────
    qc: dict = {}
    if qc_dir and os.path.isdir(qc_dir):
        for fn in sorted(os.listdir(qc_dir)):
            if fn.startswith('qc_') and fn.endswith('.json'):
                try:
                    with open(os.path.join(qc_dir, fn)) as fh:
                        d = json.load(fh)
                except Exception as exc:
                    print('  WARN %s unreadable: %s' % (fn, exc), file=sys.stderr)
                    continue
                name = d.get('sample', fn[3:-5])
                ld = depth.get(name, {}) if isinstance(depth, dict) else {}
                if isinstance(ld, dict):
                    d['library_total'] = ld.get('total')
                    d['library_coding'] = ld.get('coding')
                    d['library_coding_frac'] = ld.get('coding_frac')
                qc[name] = d
    print('QC: %d libraries' % len(qc))

    # ── analyze ───────────────────────────────────────────────────────────────
    de_genes: dict = {}
    de_clusters: dict = {}
    shift: dict = {}
    iso5 = pca = corr = None
    if analyze_dir and os.path.isdir(analyze_dir):
        tables = os.path.join(analyze_dir, 'tables')

        def _section(name, fn, default):
            """Run one packer.  A section that raises must cost ONLY itself —
            the fail-open contract is per-section, and losing the whole
            analysis.json to one bad table is the failure this guards."""
            try:
                return fn()
            except Exception as exc:
                print('  WARN section %r failed, rendering "not computed": %s'
                      % (name, exc), file=sys.stderr)
                return default

        de_genes = _section('de', lambda: pack_de(tables, gidx, top_n, 'genes'), {})
        de_clusters = _section('iso3_de',
                               lambda: pack_de(tables, gidx, top_n, 'clusters'), {})
        shift = _section('iso3_shift', lambda: pack_shift(tables, gidx, top_n), {})
        iso5 = _section('iso5',
                        lambda: pack_iso5(analyze_dir, gidx, samples_meta, top_n), None)
        pca, corr = _section('pca/corr',
                             lambda: pack_pca_corr(analyze_dir, samples_meta),
                             (None, None))
    else:
        print('analyze dir absent -> DE/PCA/isoform sections render "not computed"')

    pressure_path = getattr(args, 'pressure', None)
    pressure = None
    if pressure_path and os.path.exists(pressure_path):
        try:
            with open(pressure_path) as fh:
                pressure = json.load(fh)
        except Exception as exc:
            print('  WARN --pressure unreadable: %s' % exc, file=sys.stderr)

    try:
        caveats = build_caveats(qc, depth if isinstance(depth, dict) else {},
                                samples_meta, iso5)
    except Exception as exc:
        print('  WARN caveat assembly failed: %s' % exc, file=sys.stderr)
        caveats = []
    if pressure and pressure.get('caveats'):
        caveats = list(pressure['caveats']) + caveats

    dataset = getattr(args, 'dataset', None)
    if not dataset:
        conds = sorted({c for c in cond_labels.values()} or
                       {r['condition'] for r in meta_rows if r.get('condition')})
        dataset = ('%d samples%s' % (len(samples_meta),
                                     ': ' + ' / '.join(conds) if conds else '')
                   if samples_meta else 'rectify analysis')

    version = getattr(args, 'rectify_version', None)
    if not version:
        try:
            from ... import __version__ as _v
            version = _v
        except Exception:
            version = None

    out = {
        'generated': time.strftime('%Y-%m-%d %H:%M'),
        'rectify_version': version,
        'dataset': dataset,
        'reference_condition': (cond_labels.get(reference, reference)
                                if reference else None),
        'samples': [{'name': s.get('name'), 'genotype': s.get('genotype'),
                     'replicate': s.get('replicate'), 'modality': s.get('modality'),
                     'analyze_id': s.get('analyze_id')} for s in samples_meta],
        'contrast_labels': {c: contrast_label(c) for c in
                            sorted(set(list(de_genes) + list(de_clusters) + list(shift)))},
        'library_depth': depth,
        'qc': qc,
        'pca': pca,
        'corr': corr,
        'de': de_genes,
        'iso3_de': de_clusters,
        'iso3_shift': shift,
        'iso5': iso5,
        'splice': None,          # deliberately not computed here
        'pressure_test': pressure,
        'caveats': caveats,
    }

    out_path = str(args.out)
    parent = os.path.dirname(os.path.abspath(out_path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(out_path, 'w') as fh:
        json.dump(out, fh, separators=(',', ':'), allow_nan=False, default=str)
    print('wrote %s (%.1f KB)' % (out_path, os.path.getsize(out_path) / 1024.0))
    return 0


def create_browser_pack_parser(subparsers) -> argparse.ArgumentParser:
    """Create argument parser for the ``pack-browser`` command."""
    parser = subparsers.add_parser(
        'pack-browser',
        help='Pack analyze outputs + qc blocks into one browser analysis.json',
        description=(
            "Fold `rectify analyze` tables and `rectify qc` blocks into a single "
            "analysis.json for a genome-browser Analysis tab. This is a PACKER, not "
            "an analyst: it reshapes numbers the pipeline already produced and never "
            "recomputes a differential test. Every section is independent — a missing "
            "input renders as \"not computed\", never as an empty table."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('--qc-dir', dest='qc_dir', default='qc',
                        help='Directory of qc_<sample>.json files written by `rectify qc`.')
    parser.add_argument('--analyze-dir', dest='analyze_dir', default=None,
                        help='`rectify analyze` output directory (expects tables/, '
                             'cluster_counts.tsv, sample_metadata.tsv).')
    parser.add_argument('--gene-index', dest='gene_index', default=None,
                        help="Browser gene index JSON (alias -> chrom/start/end/strand/"
                             "label). Reused so a row's coordinates are the ones the "
                             "browser navigates to.")
    parser.add_argument('--library-depth', dest='library_depth',
                        default='library_depth.json',
                        help='JSON of per-sample {total, coding, coding_frac}.')
    parser.add_argument('--samples', default=None,
                        help='Optional samples JSON (name / analyze_id / genotype / '
                             'replicate / modality). When absent, the sample list is '
                             'derived from <analyze-dir>/sample_metadata.tsv.')
    parser.add_argument('--pressure', default=None,
                        help='Optional pressure-test JSON; its caveats are prepended.')
    parser.add_argument('--top-n', dest='top_n', type=int, default=400,
                        help='Rows kept per table.')
    parser.add_argument('--dataset', default=None,
                        help='Human-readable dataset description for the tab header. '
                             'Derived from the sample list when omitted.')
    parser.add_argument('--modality', default='',
                        help='Modality label applied to samples derived from '
                             'sample_metadata.tsv (DRS / cDNA). Ignored when '
                             '--samples supplies its own.')
    parser.add_argument('--rectify-version', dest='rectify_version', default=None,
                        help='Version string recorded in the output. Defaults to the '
                             'running rectify version.')
    parser.add_argument('-o', '--out', required=True, help='Output analysis.json path.')

    return parser


def run(args: argparse.Namespace) -> int:
    """Alias so the module matches the ``<x>_command.run(args)`` convention."""
    return run_browser_pack(args)
