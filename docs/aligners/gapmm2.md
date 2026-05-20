# gapmm2: DRS aux tags + duplicate UUIDs + PAF→BAM sequence injection

`run_gapmm2()` in `multi_aligner.py` always writes a cleaned,
deduplicated FASTQ before calling gapmm2. Two separate issues require
this.

## Issue 1: DRS auxiliary-tag headers

`trim_drs_bam_polya()` + `samtools fastq -T pt` produces FASTQ with
tab-separated auxiliary tags embedded in the read header, e.g.:

```
@6c606d1b-2310-4292-a285-d519fbd52502	pt:i:25
```

gapmm2 runs minimap2 to produce a PAF, then looks up each aligned read's
sequence via `query_idx.seq(paf[0])`. minimap2 strips everything after
the first whitespace → `paf[0]` = bare UUID. mappy 2.30 also strips tabs
when indexing, so in practice both names match. We strip tags anyway
for robustness.

## Issue 2: Duplicate UUIDs with empty sequences (TypeError crash)

DRS-trimmed FASTQs from Dorado contain a small number of reads (~1–2%)
with **duplicate UUIDs**: one entry has an empty sequence (a Dorado
placeholder) and one has the real sequence. When mappy builds its index
from a FASTQ containing duplicate read names, `seq(name)` returns `None`
for **both** entries, regardless of whether they have real sequences.
This causes a `TypeError` in gapmm2's refinement loop (`len(None)` at
`align.py:883`), which crashes gapmm2 after processing only the reads
before the first duplicate — typically ~22k out of 695k.

**Diagnosis:** crash is deterministic (same UUID every run); output PAF
has exactly N lines, where N = number of aligned reads before the first
duplicate UUID appears in minimap2's output stream.

**Fix (in `run_gapmm2`):** always write a cleaned FASTQ that:
1. Strips DRS auxiliary tags from read names (`@UUID\tpt:i:N` → `@UUID`).
2. Skips reads with empty sequences.
3. Skips subsequent occurrences of the same UUID (deduplicates).

```python
seen_uuids: set = set()
with opener(reads_path, 'rt') as src, open(tmp_fastq, 'w') as dst:
    while True:
        header = src.readline()
        if not header:
            break
        seq  = src.readline()
        plus = src.readline()
        qual = src.readline()
        clean_name = header[1:].split()[0]  # strip DRS tags
        seq_stripped = seq.rstrip()
        if not seq_stripped:
            continue  # skip empty-seq placeholders
        if clean_name in seen_uuids:
            continue  # skip duplicate UUIDs
        seen_uuids.add(clean_name)
        dst.write(f'@{clean_name}\n{seq}{plus}{qual}')
```

## gapmm2 PAF → BAM: sequence injection required

gapmm2 outputs PAF format only (no SAM/BAM). PAF does not carry query
sequences. `_paf_to_bam()` must inject sequences from the FASTQ into
each BAM record:

- Load `{uuid: (fwd_seq, fwd_qual)}` dict from the ORIGINAL FASTQ (not
  the cleaned temp FASTQ) using `_load_fastq_sequences()`, which strips
  DRS tags for key matching.
- For minus-strand reads: `query_sequence = _reverse_complement(fwd_seq)`,
  `query_qualities = array(fwd_qual)[::-1]` (pysam expects alignment
  orientation).
- For reads where `cigar_qlen != len(fwd_seq)`: skip the record (gapmm2
  cs-overrun bug: ~0.02% of reads have cs tags that over-consume query
  bases).

## mapPacBio: pt:i tag embedded in READ NAME

**mapPacBio** embeds `pt:i:N` into the READ NAME of the aligned BAM.
Separator depends on processing stage:

- **Direct mapPacBio output** (before `samtools sort`): space-separated,
  `UUID pt:i:25`. This is the format in live production data handled by
  `consensus.py` and `corrected_consensus.py`.
- **After `samtools sort/merge`**: samtools converts the space to an
  underscore (BAM spec forbids spaces in QNAME), yielding `UUID_pt:i:25`.
  This is the format in merged dev-run BAMs and the sorted validation
  aligner BAMs.

Use whichever form is present:

```python
# Direct mapPacBio output (pre-sort):
if " pt:i:" in name:
    clean = name.split(" pt:i:")[0]

# Merged / sorted BAMs:
if "_pt:i:" in name:
    clean = name.split("_pt:i:")[0]
```

This has been encountered multiple times when extracting validation
reads from DRS-trimmed dev-run outputs.
