# rectify.visualize

Visualization subpackage. Requires `pip install rectify-rna[visualize]` (adds `matplotlib` and `seaborn`).

```python
# Check if visualization is available
import rectify
if rectify.visualize:
    from rectify.visualize.metagene import MetagenePipeline
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
        - load_bam_signal
        - load_bigwig_signal
        - load_corrected_tsv_signal

---

## multi_track

Multi-track genome browser figure generation (coverage, gene annotation, reads).

::: rectify.visualize.multi_track
    options:
      members:
        - MultiTrackFigure
        - add_coverage_track
        - add_gene_track
        - add_read_browser_track

---

## coverage

Per-base coverage track plotting.

::: rectify.visualize.coverage
    options:
      members:
        - plot_coverage
        - compute_coverage

---

## gene_track

Gene annotation track drawing.

::: rectify.visualize.gene_track
    options:
      members:
        - plot_gene_track
        - GeneTrackConfig

---

## read_browser

Read-level alignment visualization.

::: rectify.visualize.read_browser
    options:
      members:
        - ReadBrowser
        - plot_read_alignment

---

## figure_utils

Plotting utilities and style helpers.

::: rectify.visualize.figure_utils
    options:
      members:
        - set_style
        - save_figure
        - add_scale_bar

---

## Example: metagene around CPA sites

```python
from rectify.visualize.metagene import MetagenePipeline, MetageneConfig

config = MetageneConfig(
    flank_upstream=200,
    flank_downstream=200,
    bin_width=1,
    normalize_by_transcript_length=False,
)

pipeline = MetagenePipeline(config)
pipeline.load_regions(
    'cpa_clusters.bed',      # BED file with CPA cluster positions
    filter=lambda r: r['count'] > 50,  # Only clusters with > 50 reads
)

# Aggregate signal from corrected BAM
signal_matrix = pipeline.aggregate(
    bam_path='corrected.bam',
    output_prefix='metagene_wt',
)

# Plot
pipeline.plot(
    signal_matrix,
    output_path='metagene_wt.png',
    title='WT 3\' end signal around CPA sites',
    colormap='Blues',
)
```
