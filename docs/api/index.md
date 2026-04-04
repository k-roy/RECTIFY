# API Reference

RECTIFY's Python API is organized into four subpackages. All modules can be imported directly:

```python
import rectify
from rectify.core.bam_processor import process_bam_streaming
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
from rectify.core.bam_processor import process_bam_streaming
from rectify.utils.genome import load_genome

genome = load_genome('genome.fa.gz')

stats = process_bam_streaming(
    bam_path='reads.bam',
    genome=genome,
    output_path='corrected.tsv',
    num_processes=8,
    streaming=True,
    chunk_size=10000,
)
print(f"Corrected {stats.total_reads} reads, mean shift: {stats.mean_shift:.1f} bp")
```

### Cluster CPA positions

```python
from rectify.core.analyze.clustering import cluster_cpa_sites_adaptive
import pandas as pd

# positions: list of dicts with 'chrom', 'strand', 'pos', 'count'
positions = {'chrI': [{'pos': 34521, 'strand': '+', 'count': 150}, ...]}
clusters = cluster_cpa_sites_adaptive(positions, min_cluster_size=10)
```

### Load genome

```python
from rectify.utils.genome import load_genome, fetch_genomic_sequence

genome = load_genome('genome.fa.gz')  # Returns Dict[str, str]
seq = fetch_genomic_sequence(genome, 'chrI', 34500, 34550, strand='+')
```

### Run metagene analysis

```python
from rectify.visualize.metagene import MetagenePipeline, MetageneConfig

config = MetageneConfig(flank_upstream=200, flank_downstream=200, bin_width=1)
pipeline = MetagenePipeline(config)
pipeline.load_regions('genes.bed')
signal = pipeline.aggregate('corrected.bam', 'output_prefix')
pipeline.plot(signal, 'metagene.png', title='CPA metagene')
```

---

## Version

```python
import rectify
print(rectify.__version__)  # '2.7.6'
```
