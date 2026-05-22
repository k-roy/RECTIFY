# rectify.visualize — Plot Skills Reference

Authoritative index of what's in `rectify.visualize`, when to use each tool,
and what to avoid. For runnable code recipes, see the README visualization
section and `docs/examples/viz_examples.ipynb` (forthcoming).

Install: `pip install rectify-rna[visualize]`

---

## Metagene Analysis

### When to use
Any time you want signal aggregated across genomic loci — CPA sites, TRT/T-tracts,
gene TES/TSS, motif matches, etc.

### Core objects

| Symbol | Purpose |
|--------|---------|
| `MetagenePipeline` | Orchestrates extraction + aggregation |
| `PositionIndex` | O(1) position lookup from a reads DataFrame |
| `MetageneConfig` | Window size, normalization, capping settings |
| `StrandOrientationError` | Raised when + and − strand peaks diverge (strand bug) |
| `verify_strand_balance()` | Standalone checker; call before any plot |

### Key method: `compute_center_profile()`

Use this for **center-based** metagene (TRT, CPA, motif matches — single position
per locus). Handles strand coordinate transformation and reversal internally.

```python
from rectify.visualize import MetagenePipeline, position_index_from_tsv, loci_from_pickle

index, total_reads = position_index_from_tsv(tsv_paths, position_col='corrected_3prime')
loci = loci_from_pickle(cache_path)

pipeline = MetagenePipeline()
result = pipeline.compute_center_profile(
    loci, index,
    window=(-50, 50),
    total_reads=total_reads,
    verify_strands=True,    # raises StrandOrientationError on strand bug
    cap_percentile=50,      # cap outlier loci before peak detection
)
# result keys: profile, profile_matrix, profiles_plus, profiles_minus,
#              x, window, n_loci, n_plus, n_minus, strand_verification
```

**Always set `cap_percentile=50`** when using window ≥ 100 bp — a single
high-expression gene can dominate the aggregate and mask the biological peak.

Use `compute_profile()` instead for **gene-body** metagenes (scaled body +
flanking windows).

### Locus loaders

| Function | Input | Use for |
|----------|-------|---------|
| `loci_from_pickle(path)` | `.pkl` cache | TRT/CPA caches from this project |
| `loci_from_gff(path, feature_type, center)` | GFF3 file | Gene TSS (`center='start'`), TES (`center='end'`) |
| `loci_from_motif_scan(sequences, motif)` | dict of sequences | Any regex motif, both strands |
| `loci_from_bed(path, center)` | BED file | Custom interval lists |
| `loci_from_tsv(path)` | TSV with chrom/strand/center | Pre-computed loci |
| `position_index_from_tsv(paths)` | RECTIFY output TSV | WT/mutant 3′ end signal |
| `position_index_from_bigwig(path)` | bigWig | Any coverage signal |

### Strand verification rules
- **Always** run `verify_strands=True` (default on `compute_center_profile`).
- If verification fails: check center coordinate convention, not just the reversal.
  Common bug: using `trt_end` instead of `trt_end - 1` for minus strand center.
- If peaks agree but at wrong position: window may be too small to see the signal.
- High-expression outliers (e.g. TDH3) can shift the aggregate peak to the window
  edge — use `cap_percentile=50` to suppress them.

---

## Stacked Read Browser

### When to use
Inspecting individual nanopore reads at a gene locus — splice patterns, read
length distributions, internal pA sites, condition comparisons.

### Key functions

| Function | Purpose |
|----------|---------|
| `assign_rows(starts, ends, gap)` | Row packing for non-overlapping layout |
| `parse_junction_strings(col)` | Convert RECTIFY junction strings to lists |
| `draw_stacked_reads(ax, starts, ends, rows, colors, junction_lists)` | Low-level, LineCollection-based |
| `plot_stacked_read_panel(ax, reads_df, color_col, ...)` | High-level convenience wrapper |

### Performance note
Use `plot_stacked_read_panel` / `draw_stacked_reads` rather than per-read
`ax.plot()` calls. LineCollection reduces 400 draw calls to ~6 (one per color
category). For 400 reads × 5 conditions the difference is ~10-50× faster render.

### Layout tip
Use `constrained_layout=True` in `plt.subplots()` instead of `plt.tight_layout()`
when mixing axes that have `axis('off')` (legend rows, gene track rows). This
avoids the `UserWarning: Axes not compatible with tight_layout` warning.

```python
fig, axes = plt.subplots(n_rows, 1, figsize=(14, h),
                          gridspec_kw={'height_ratios': [...]},
                          constrained_layout=True)
# Do NOT call plt.tight_layout()
```

### Color assignment is your responsibility
The library draws what you give it. Classification (CDS-internal, canonical,
readthrough, etc.) is project-specific and must happen before calling the library.
Assign a `'color'` column to your DataFrame before passing it to
`plot_stacked_read_panel`.

---

## Gene Track

### When to use
Genome browser-style gene structure panels showing CDS exons, UTRs, and strand
direction.

| Function | Purpose |
|----------|---------|
| `assign_feature_levels(features)` | Prevent gene overlaps by staggering y-levels |
| `draw_gene_track(ax, gene_name, gff_features, ...)` | Full gene track from GFF |
| `get_genes_in_region(gff_features, chrom, start, end)` | Subset features to a window |
| `draw_gene_arrow(ax, ...)` | Single gene arrow shape |

---

## Coverage Tracks

| Function | Purpose |
|----------|---------|
| `extract_coverage_from_bam(bam, chrom, start, end)` | Per-base coverage array from BAM |
| `extract_coverage_from_array(arr, ...)` | Smooth / subsample a coverage array |
| `draw_coverage_track(ax, coverage, ...)` | Fill-under coverage plot |
| `draw_strand_coverage(ax, plus, minus, ...)` | Dual-strand coverage (+ above, − below) |
| `compare_coverage_tracks(ax, tracks, ...)` | Overlay multiple conditions |

---

## Multi-Track Browser

### When to use
Composite figures with gene track + coverage + reads in a vertically stacked layout.

```python
from rectify.visualize import MultiTrackFigure

fig = (MultiTrackFigure()
    .add_gene_track(gff_features, highlight_gene="ENA1")
    .add_coverage_track("NET-seq", netseq_coverage)
    .add_coverage_track("RECTIFY WT", rectify_coverage)
)
fig.save("browser.png", region_start=530000, region_end=535000)
```

---

## Figure Utilities

| Function | Purpose |
|----------|---------|
| `set_publication_style()` | Apply consistent fonts, tick sizes, rc params |
| `despine(ax)` | Remove top and right spines |
| `save_multi_format(fig, path)` | Save as `.png` + `.pdf` + `.svg` in one call |
| `trimmed_mean(arr, q)` | Mean after excluding top/bottom q fraction |
| `apply_window_sum_capping(profiles, percentile)` | Cap outlier loci by total signal |
| `plot_metagene_line(ax, x, mean, sem, color, label)` | One-liner for mean ± SEM line |
| `add_metagene_annotations(ax, positions, labels)` | Add labeled vertical reference lines |

---

## Color Palettes

| Name | Contents |
|------|----------|
| `WONG_COLORS` | Color-blind-safe 8-color palette (blue, orange, green, red, sky, vermillion, bluegreen, yellow) |
| `CODON_VARIANT_COLORS` | syn / missense / nonsense / splice |
| `EFFECT_DIRECTION_COLORS` | up / down / neutral |
| `GENE_TYPE_COLORS` | essential / nonessential / etc. |

**Always prefer `WONG_COLORS`** for new figures. Use the other palettes only for
their specific data types.

---

## What NOT to do

- **Don't call `ax.plot()` per read** in stacked read plots — use `draw_stacked_reads`.
- **Don't skip `verify_strands`** in metagene — the visual is the last place to
  catch a strand bug.
- **Don't use `plt.tight_layout()`** when any panel has `axis('off')` — use
  `constrained_layout=True` in `plt.subplots()`.
- **Don't set `cap_percentile=None`** with windows ≥ 100 bp — one high-expression
  gene will dominate the aggregate peak.
- **Don't load deprecated caches** — use `load_loci_from_cache()` in
  `metagene_workflow.py` which checks `DEPRECATED_CACHES` and raises on bad inputs.
