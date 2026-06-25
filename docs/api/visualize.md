# rectify.visualize

Visualization subpackage. Requires `pip install rectify-rna[visualize]` (adds `matplotlib` and `seaborn`).

```python
# Visualization deps are optional; guard the import
try:
    from rectify.visualize.metagene import MetagenePipeline
    HAVE_VISUALIZE = True
except ImportError:
    HAVE_VISUALIZE = False
```

---

## metagene

Metagene analysis — aggregate signal around genomic regions and plot.

::: rectify.visualize.metagene
    options:
      members:
        - MetageneConfig
        - MetagenePipeline

---

## metagene_loaders

Data loaders for metagene input formats.

::: rectify.visualize.metagene_loaders
    options:
      members:
        - loci_from_bed
        - loci_from_gff
        - loci_from_tsv
        - loci_from_motif_scan
        - loci_from_pickle
        - position_index_from_tsv
        - position_index_from_bigwig

---

## multi_track

Multi-track genome browser figure generation (coverage, gene annotation, reads).

::: rectify.visualize.multi_track
    options:
      members:
        - MultiTrackFigure
        - create_gene_browser

---

## coverage

Per-base coverage track plotting.

::: rectify.visualize.coverage
    options:
      members:
        - draw_coverage_track
        - draw_strand_coverage
        - extract_coverage_from_bam
        - extract_coverage_from_array
        - compare_coverage_tracks

---

## gene_track

Gene annotation track drawing.

::: rectify.visualize.gene_track
    options:
      members:
        - draw_gene_track
        - draw_gene_arrow
        - assign_feature_levels
        - get_genes_in_region
        - GENE_TYPE_COLORS

---

## read_browser

Read-level alignment visualization.

::: rectify.visualize.read_browser
    options:
      members:
        - plot_stacked_read_panel
        - draw_stacked_reads
        - assign_rows
        - parse_junction_strings

---

## figure_utils

Plotting utilities and style helpers.

::: rectify.visualize.figure_utils
    options:
      members:
        - set_publication_style
        - save_multi_format
        - despine
        - format_genomic_axis
        - plot_metagene_line
        - add_metagene_annotations

---

## Example: metagene around CPA sites

```python
import pandas as pd
import matplotlib.pyplot as plt
from rectify.visualize.metagene import MetagenePipeline, MetageneConfig
from rectify.visualize.metagene_loaders import (
    loci_from_bed,
    position_index_from_bigwig,
)

config = MetageneConfig(window_upstream=200, window_downstream=200)
pipeline = MetagenePipeline(config)

# Loci (regions) as a DataFrame, plus a position index of the signal.
# position_index_from_bigwig returns (PositionIndex, total_signal) per strand.
loci = pd.DataFrame(loci_from_bed('cpa_clusters.bed', center='start'))
position_index, _total = position_index_from_bigwig('wt.plus.bw', strand='+')

result = pipeline.compute_profile(loci, position_index)

fig, ax = plt.subplots()
pipeline.plot_profile(ax, result, label="WT 3' ends")
ax.set_title("WT 3' end signal around CPA sites")
fig.savefig('metagene_wt.png', dpi=200)
```

The exact loader and `position_index_from_*` arguments depend on your input
format — see the `metagene_loaders` reference above.
