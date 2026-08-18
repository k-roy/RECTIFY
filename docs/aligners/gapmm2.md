# gapmm2: DRS aux tags + duplicate UUIDs + PAF→BAM sequence injection

> ⚠ **PINNED to `gapmm2==25.4.5`** (`pyproject.toml` + `install_aligners_command.py`).
> **Do NOT upgrade.** gapmm2 25.4.13+ added a `splice_aligner_minimap2` code path
> (used whenever the `minimap2` binary is on PATH) whose tag loop does
> `if ts is None: continue` — and minimap2 emits the `ts:A:` transcript-strand tag
> **only for spliced alignments**, so 25.8.12 **silently drops every single-exon
> (unspliced) read** (10/35 validation reads vs 25.4.5's 35/35; the drop is invisible
> in `--debug`). The faster 25.8.12 binary-path timings quoted in the Performance
> section below are real, but that version is **unusable** until upstream restores
> tolerance for non-spliced alignments. BUGS_TO_FIX NEW-082 (root-caused 2026-06-24).

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

## gapmm2 CLI: PAF-only output (no SAM mode)

The gapmm2 CLI (v25.4.5) only outputs PAF or GFF3:

```
-f , --out-format   output format [paf,gff3] (default: paf)
```

There is no `-a` / SAM flag. Piping `gapmm2 ... | samtools sort` fails with
`samtools sort: failed to read header from "-"`. To get a BAM from the CLI,
write PAF to a file and parse the `cs:Z:` diff-string tags (introns encoded
as `~XYlenXY`). RECTIFY's `run_gapmm2()` uses the Python API (`mappy` +
`edlib`) and handles PAF→BAM internally.

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

(For the full mapPacBio parameter/suitability story, see
[mapPacBio.md](mapPacBio.md).)

---

## Performance: the terminal-refinement loop is the bottleneck (not mapping, not `-i`)

gapmm2 = minimap2 mapping + a single-threaded Python edlib loop that rescues
terminal exons minimap2 soft-clipped. On large mammalian samples that refinement
loop dominates and blows past the 6 h `ALIGNER_TIMEOUT` (multi_aligner.py:232).

**Empirical (human chr5 ONT DRS, CNTL_HB53 = 883,350 reads, 2026-05-25):**

- `-i` / `--max-intron` (the edlib terminal-window size) has **no effect on wall
  time**: 582 s → 606 s across `-i` ∈ {10k…500k} on a 3,000-read subset (~6 %
  spread). It is NOT the bottleneck — and a larger window rescues slightly more
  terminal exons for free, so do **not** cap it for speed.
- Decomposition: chr5 mappy index build is only ~5 s; per-read cost ≈ **0.19
  s/read** (v25.4.5) → full sample ≈ **46–50 h** single-threaded.
- v25.8.12 (uses the minimap2 **binary** when on PATH) with `--no-refine` maps
  1,000 reads in **7 s** (~6 min full sample) — mapping is essentially free. With
  refinement on it is **0.116 s/read** (~28 h). The refinement loop is the entire
  bottleneck; the binary upgrade only buys ~1.6×. gapmm2 `--no-refine` ≈ plain
  minimap2 (redundant with the panel's minimap2).

**Node architecture note (2026-05-25):** SG-NEx A549 chr5, 10k-read benchmark on sh02-03n47
(old Sherlock Intel node) took **4,644 s → 0.46 s/read**: ~2× slower than AMD Milan nodes
nodes. A separate Intel owners-queue node (sh04-05n15) measured **0.21 s/read** vs 0.19 s/read
on AMD Milan — ~10% variance between similar-generation Intel nodes. Treat 0.19–0.21 s/read
as the range for modern nodes (AMD Milan / recent Intel); old sh02 Intel is ~0.45 s/read.

**Scaling fix (if gapmm2 is needed):** the cost is per-read, so chunk the FASTQ
into N parts and run gapmm2 per part in parallel (no built-in chunking — wrap it:
split → parallel `rectify align --aligners gapmm2` → `samtools merge`), keeping
`-i` large. The per-chunk index build (~5 s) is negligible.

### Minimum-overhang gate: precision improvement + potential speedup

gapmm2 triggers terminal refinement for **any** terminal soft-clip: `h.q_st > 0`
(5′ clip) or `len(seq) > h.q_en` (3′ clip) in `align.py`. It pays the full
per-attempt cost — window fetch + `find_all_splice` + edlib — for a 2 bp clip
exactly as for a 200 bp one.

A two-line in-place patch gates those triggers on a minimum overhang (same
pattern as the prior `--min-mapq` argparse fix — apply per-install):

```python
MIN_OVERHANG = 10
if h.q_st >= MIN_OVERHANG: ...               # was: if h.q_st > 0
if (len(seq) - h.q_en) >= MIN_OVERHANG: ...  # was: if len(seq) > h.q_en
```

**Measured** (v25.8.12 binary path, 3,000-read CNTL_HB53 subset, `GATE=10`,
2026-05-25 — `env GAPMM2_MIN_OVERHANG` gate):

- **Speed: only ~10 %** (220 s → 199 s). The gate skipped 28.6 % of refinement
  *calls*, but those are the **cheap** ones — short clips where `find_all_splice`
  finds nothing and edlib barely runs. The expensive calls are the long overhangs
  the gate deliberately keeps. So the gate is **not** a meaningful speedup; the
  small-overhang attempts are not where the time goes. **Chunking remains the only
  real speed fix.**
- **Precision: clean win.** The gate dropped 98 introns (315 → 217, keeps 69 %),
  and **100 % of the dropped introns were novel/unannotated — 0 annotated (real)
  junctions lost.** Confirms short clips are the forced-GT-AG artifact source
  (≥85 % of gapmm2's unique introns; see next section). Genuine short terminal
  exons are rare and already covered by uLTRA/deSALT.

**Bottom line:** apply the ≥10 bp gate for **quality** (a near-free artifact filter
on gapmm2's terminal rescues), not for speed. It stacks with chunking; it does not
replace it. (The initial "this should be a big speedup" hypothesis was wrong — a
good reminder to measure: the synthetic `-i` microbenchmark and the speedup
intuition both missed where the cost actually is.)

## Does gapmm2 recover real junctions the other aligners miss?

Per-read intron comparison vs minimap2/uLTRA/deSALT/mapPacBio (CNTL_21.8 and
SMA_GSB2394, 2026-05-25): only **1.3–1.8 %** of gapmm2 introns are unique to
gapmm2 for that read, and **~85 % of those uniques are novel/unannotated** with an
artifact signature (clustered forced GT-AG terminal rescues to the same acceptor
with scattered donors). Only ~550–840 unique introns/sample match an annotated
GENCODE junction (~0.2 % of gapmm2's introns). Net: small genuine contribution,
bundled with ~5× as many likely artifacts, ~98 % redundant with uLTRA/deSALT.

**SG-NEx A549 chr5 10k-read benchmark (2026-05-25):** on a 10k-read subset gapmm2
achieved 97.2% GT-AG at **4,644 s** (sh02 Intel, 0.46 s/read). After correcting a strand
bug in `compute_junction_stats.py` (the BAM-based script omitted `aln.is_reverse` handling
— reverse-strand reads had donor/acceptor swapped, producing ~50% spurious non-canonical),
the corrected figures are: minimap2 91.6%, minisplice_mm2 96.0%, **winnowmap2 99.0%**,
gapmm2 97.2%. Winnowmap2 matches or exceeds gapmm2's GT-AG at 1/516th the wall time; the
forced GT-AG terminal rescues are not the differentiator that initially appeared.

**Decision (Sumner ONT-DRS chr5, 2026-05-25):** gapmm2 dropped from the human-DRS
panel — marginal real-junction yield does not justify ~200 chunked array tasks
(hourly billing), and uLTRA/deSALT already cover splice-aware terminal exons.

---

## Primary-alignment & duplicate handling

gapmm2 emits PAF, which rectify converts to BAM in `_paf_to_bam()`. Secondary
alignments are handled at conversion: any PAF record with `tp:A` ≠ `P` is marked
`FLAG |= 0x100` (in `_paf_to_bam()`, `multi_aligner.py`), so `rectify correct`'s `is_secondary`
filter drops them. Separately, `run_gapmm2()` **deduplicates the input FASTQ by UUID
before alignment** (Issue 2 above) — it skips empty-sequence placeholders and
subsequent occurrences of the same read name. So a doubled input FASTQ that would
2×-inflate minimap2 is *collapsed* by gapmm2's pre-align dedup — a useful contrast,
and the reason gapmm2's per-read counts can legitimately differ from the other
aligners on the same (duplicated) input.

That pre-align dedup lives inside `run_gapmm2()` and does not protect an external BAM
fed to `rectify correct`. The panel-wide unguarded hazard remains duplicate
**primary** records in a BAM `rectify correct` reads directly. See the canonical
writeup and cross-aligner table in
[minimap2.md](minimap2.md#-duplicate-primary-alignments--2-double-counted-3-ends-external-bam-hazard).
