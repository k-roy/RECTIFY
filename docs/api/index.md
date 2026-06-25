# API Reference

RECTIFY's Python API is organized into four subpackages. All modules can be imported directly:

```python
import rectify
from rectify.core.bam.parallel import process_bam_streaming
from rectify.core.analyze.clustering import cluster_cpa_sites_adaptive
from rectify.utils.genome import load_genome, reverse_complement
```

---

## Package overview

| Subpackage | Description |
|-----------|-------------|
| [`rectify.core`](core.md) | Core correction pipeline: BAM processing, indel correction, alignment consensus |
| [`rectify.core.analyze`](analyze.md) | Downstream analysis: clustering, DESeq2, APA detection, GO enrichment, motifs |
| [`rectify.utils`](utils.md) | Shared utilities: genome I/O, CIGAR parsing, chromosome names, statistics |
| [`rectify.visualize`](visualize.md) | Visualization: metagene plots, genome browser figures, heatmaps |

---

## Quick reference

### Correct a BAM file programmatically

```python
from rectify.core.bam.parallel import process_bam_streaming

stats = process_bam_streaming(
    bam_path='reads.bam',
    genome_path='genome.fa.gz',   # a path, not a loaded genome
    output_path='corrected.tsv',
    chunk_size=10000,
)
print(f"Reads in BAM: {stats.total_reads_in_bam}, processed: {stats.reads_processed}")
```

### Cluster CPA positions

```python
from rectify.core.analyze.clustering import cluster_cpa_sites_adaptive
import pandas as pd

# A DataFrame with chrom / strand / corrected_position columns
# (count_col is optional; pass it to weight clusters by read count).
positions_df = pd.DataFrame({
    'chrom': ['chrI', 'chrI'],
    'strand': ['+', '+'],
    'corrected_position': [34521, 34530],
})
clusters = cluster_cpa_sites_adaptive(
    positions_df, max_cluster_radius=10, min_peak_separation=5, min_reads=5,
)
```

### Load genome

```python
from rectify.utils.genome import load_genome, fetch_genomic_sequence

genome = load_genome('genome.fa.gz')  # Returns Dict[str, str]
seq = fetch_genomic_sequence(genome, 'chrI', 34500, 34550, strand='+')
```

### Run metagene analysis

```python
import matplotlib.pyplot as plt
from rectify.visualize.metagene import MetagenePipeline, MetageneConfig

config = MetageneConfig(window_upstream=200, window_downstream=200)
pipeline = MetagenePipeline(config)

# loci (DataFrame of regions) + position_index of the signal — see the
# metagene_loaders reference and rectify.visualize docs for builders.
result = pipeline.compute_profile(loci, position_index)

fig, ax = plt.subplots()
pipeline.plot_profile(ax, result, label="3' ends")
fig.savefig('metagene.png', dpi=200)
```

---

## Version

```python
import rectify
print(rectify.__version__)  # '0.9.0'
```
