# Minimap2: junction annotation + the `-G` / `--max-intron` human pitfall

The `rectify align` command generates a junction BED from the GFF
annotation and passes it to minimap2 via `--junc-bed`. This improves
splice junction accuracy for long-read RNA-seq.

The junction BED is cached as `annotation.junc.bed` in the sample output
directory. Exact minimap2 command:

```
minimap2 -ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD -y
    -t <threads>
    --junc-bed <sample_dir>/annotation.junc.bed
    --junc-bonus 9
    <genome.fsa.gz> <reads.fastq.gz>
```

Key flags:
- `-uf`: forward-strand-only (correct for direct RNA / cDNA sense reads).
- `-k14`: smaller k-mer for sensitivity on noisy nanopore reads.
- `-G 5000`: max intron (chainable gap) size. **This is the yeast default and
  the #1 misconfiguration on mammalian data — see below.**
- `--splice-flank=no`: disabled for compatibility (per the `run_minimap2`
  source comment).
- `--MD`: required for indel artifact correction downstream.
- `-y`: copies FASTQ comment fields (SAM-format aux tags, e.g. the DRS
  `pt:i` poly(A) tag) through into the BAM aux records.

---

## ⚠️ `-G` comes from `--max-intron` (default 5000 = yeast) and `--organism` does NOT raise it

`-G` is minimap2's maximum intron length: minimap2 will not chain a spliced
alignment across a gap larger than `-G`. RECTIFY sets `-G` from the
`rectify align --max-intron` flag, which **defaults to 5000**
(`align_command.py`, help text: *"appropriate for S. cerevisiae. Use 500000 or
larger for human data."*). Critically, **`--organism homo_sapiens` selects the
bundled genome/annotation but does NOT raise `--max-intron`** — so unless you
pass it explicitly, minimap2 runs with `-G 5000` on human data.

**Consequence on mammalian RNA:** human introns are routinely 10 kb–1 Mb. With
`-G 5000`, a read spanning a larger intron cannot be chained as one spliced
alignment — minimap2 soft-clips the distal exon or maps a fragment. minimap2 is
the consensus workhorse (~49 % of per-read wins on A549 chr5), so a wrong `-G`
silently degrades the *primary* alignment for large-intron transcripts.

**How bad, measured (A549 chr5 ONT DRS, 2026-05-25, panel run WITHOUT
`--max-intron` → `-G 5000`):** the impact on minimap2 is real but **bounded**,
not catastrophic. The consensus combo `deSALT+minimap2+uLTRA` still agreed on
312,796 reads, i.e. minimap2 tracked the splice-aware aligners on most reads —
the chr5 read set is dominated by transcripts whose introns are all <5 kb.
Damage is concentrated in the minority of reads spanning a >5 kb intron. (Contrast
mapPacBio, which was badly broken by its own intron misconfiguration — ~1k introns
total, a disjoint read set — see [mapPacBio.md](mapPacBio.md).)

**Caveat for junction analysis:** at large-intron loci a `-G 5000` minimap2
clipping at a 5 kb boundary while uLTRA/deSALT place the real 50 kb intron looks
like a genuine algorithmic *disagreement* (Cat6/Cat9-style) when it is only a
parameter artifact. Do not trust minimap2-vs-{uLTRA,deSALT} junction
disagreements at large-intron genes on human data until `-G` is set correctly.

**Fix — always pass `--max-intron 500000` (or larger) for human/mammalian:**

```
rectify align <reads> --genome <GRCh38> --organism homo_sapiens \
    --max-intron 500000  ...
```

**`--max-intron` is overloaded across the panel** — know what it does to each:
- **minimap2** `-G` ← this is the one it controls correctly; set it to 500000 for human.
- **gapmm2** `-i` (edlib terminal-rescue window) — no wall-time effect, larger is
  free; see [gapmm2.md](gapmm2.md).
- **mapPacBio** `intronlen` — historically (and incorrectly) wired here; mapPacBio
  needs `intronlen=10` + `maxindel=200000` instead; see [mapPacBio.md](mapPacBio.md).

---

## Quick check

```bash
rectify align <reads> --genome <GRCh38> --organism homo_sapiens \
    --max-intron 500000 --aligners minimap2 -o /tmp/mm2_smoke -t 8
# spliced reads should show N ops well over 5 kb on a mammalian sample:
samtools view /tmp/mm2_smoke/*.minimap2.bam \
  | awk '$6 ~ /[0-9]{4,}N/' | head
```

If the largest `N` op tops out near 5,000 on human data, `-G` is still at the
yeast default — `--max-intron` was not raised.

---

## ⚠️ Duplicate **primary** alignments → 2× double-counted 3′ ends (external-BAM hazard)

> **This section is the canonical writeup of secondary/supplementary/duplicate
> handling for the whole panel.** Each aligner doc carries a short
> "Primary-alignment & duplicate handling" note pointing back here; the
> cross-aligner table below is the single source of truth.

**Incident (by4742 in-house DRS, 2026-05-29).** The pre-existing minimap2 BAMs in
`projects/TRT/intermediate_data/nanopore/inhouse_by4742_dst1_4nqo/` contain **every
read as TWO byte-identical PRIMARY alignment records** (same QNAME, POS, CIGAR,
both FLAG 0/16 — neither secondary nor supplementary). Measured on `wt_by4742_rep1`:

```
chrI:    300,010 records / 149,752 distinct QNAMEs        ≈ 2.0×
chrMito: 1,656 primary QNAMEs each appearing exactly 2×   (-F 0x900)
```

The source **Dorado BAMs are clean** (`samtools view -F 0x900 | cut -f1 | sort |
uniq -d` → 0 duplicates). The duplication was introduced by the *external*
alignment command that produced these intermediate BAMs —
`minimap2 -ax splice -G 3000` (no `--secondary=no`), fed a FASTQ that contained
each read twice (doubled/concatenated input or the file passed twice). It is **not
a rectify-produced artifact** and **not** a secondary/supplementary-alignment
problem.

**Why it matters:** `rectify correct` emits one `corrected_3ends.tsv` row per
processed read and counts each `corrected_3prime` once. Two identical primaries →
two rows → the per-position pileup and every emitted bedgraph are inflated **~2×**.
Verified empirically: a chrI smoke processed 281,088 primary records for ~149,752
distinct reads — exactly the doubling.

### "Restrict to primary alignments" — already done, but it does NOT fix this case

`rectify correct` already restricts to primary alignments. Every BAM-iteration loop
skips unmapped/secondary/supplementary and nothing else
(`rectify/core/bam/parallel.py:952-958`; non-parallel fallback
`rectify/core/bam/bam_processor.py:1184-1188`):

```python
if read.is_unmapped:      continue   # 0x4
if read.is_secondary:     continue   # 0x100
if read.is_supplementary: continue   # 0x800
```

So secondary/supplementary records from **any** aligner are dropped before 3′-end
counting — that protection is real and panel-wide. **But the by4742 duplicates are
PRIMARY records**, so this filter keeps both. Two gaps in `rectify correct` let them
through (verified by code audit, 2026-05-29):

- It does **not** inspect the duplicate flag (`0x400` / `read.is_duplicate`) — a read
  marked by `samtools markdup` is still processed as a normal primary.
- It does **not** deduplicate by `query_name`. The `seen_read_ids`/`is_primary_result`
  mechanism (`_parse_tsv_result`, `parallel.py:368-386`) gates *summary-stat counters only* (NET-seq
  proportional row-splitting); it does not suppress TSV rows or position counts.

The two existing one-primary-per-read protections both live in `rectify align` and
**never touch an external BAM handed to `rectify correct`**:
`_dedup_desalt_bam()` (deSALT-only) and the consensus winner-promotion
(`rectify/core/consensus/consensus.py:617`, `flag &= ~0x900`) — and per-aligner
`rectify correct` runs *before* the consensus merge ([[per_aligner_rescue_runs_first]]).
`tests/test_no_duplicate_primaries.py` exercises only the consensus path; the
per-aligner `correct` → TSV path has no such test or guard.

### Cross-aligner secondary-suppression status (`rectify align`)

What flag (if any) `rectify align` passes to keep each aligner from emitting
secondary alignments. `rectify correct` filters secondary/supplementary records
regardless, so these flags are belt-and-suspenders — the real, unguarded hazard is
duplicate **primaries**, which only an upstream/data fault produces.

| Aligner | rectify-align secondary control | Source | Clean primaries? |
|---|---|---|---|
| **minimap2** | `--secondary=no` | `run_minimap2`, `multi_aligner.py:391` | ✓ |
| **winnowmap2** | `--secondary=no` | `run_winnowmap2`, `multi_aligner.py:2001` | ✓ |
| **minisplice_mm2** | `--secondary=no` | `run_minisplice_mm2`, `multi_aligner.py:2143` | ✓ |
| **gapmm2** | PAF output; PAF→BAM marks `tp:A`≠`P` as `0x100`; also dedups duplicate UUIDs in the input FASTQ pre-align | `_paf_to_bam`, `multi_aligner.py:1735-1741`; [gapmm2.md](gapmm2.md) Issue 2 | ✓ |
| **bwa mem** | `-M` (marks split hits secondary → filtered) | `run_bwa_mem`, `multi_aligner.py:1136` | ✓ |
| **deSALT** | no flag; post-filtered by `_dedup_desalt_bam()` on `(name,flag,chrom,pos,cigar)` (deSALT emits each aln `-N`× by default) | `_dedup_desalt_bam`, `multi_aligner.py:2601`; [deSALT.md](deSALT.md) | ✓ (after dedup) |
| **mapPacBio** | **no secondary-control flag passed**; BBMap default behavior not determinable from rectify source | `run_map_pacbio`, `multi_aligner.py:693` | secondaries (if any) dropped by `correct`'s `is_secondary` filter |
| **bbmap** (short-read) | `ambiguous=best` (picks one best site; not a secondary-record control) | `run_bbmap`, `multi_aligner.py:1015` | as above |
| **uLTRA** | **no secondary-control flag passed** | `run_ultra`, `multi_aligner.py:2387` | secondaries (if any) dropped by `correct`'s `is_secondary` filter |

`rectify align` does **no** `samtools markdup`/`rmdup` and **no** `view -F 0x900`
post-filter (audit: zero hits in `rectify/`); deSALT's `_dedup_desalt_bam` is the
only physical dedup, and it's deSALT-specific.

### Remediation

1. **Preferred — re-align through rectify.** A `rectify align` / `rectify run-all`
   pass produces clean BAMs (`--secondary=no`) *and* applies the DRS poly-A/adapter
   pre-trim, so externally-doubled BAMs become moot. This is the correct fix when
   the upstream FASTQ provenance is suspect.
2. **If reusing an external BAM you cannot re-align:** dedup *before* `rectify correct`.
   `samtools view -F 0x900` will **not** help (the dups are primary). Dedup by
   coordinate+name, e.g. `samtools markdup -r` after `collate`+`fixmate`, or keep one
   record per `(QNAME,FLAG,RNAME,POS,CIGAR)` tuple (the `_dedup_desalt_bam` pattern).
3. **Last resort — dedup the output.** `corrected_3ends.tsv` is dedup-able by
   `read_id`; this also collapses any region-shard boundary dups. Do this only if
   (1)/(2) are impossible, since the wasted ~2× correct-stage compute remains.

**Recommended code hardening (not yet implemented):** give `rectify correct` an
opt-in `--dedup-primary` (skip `is_duplicate`, and/or keep first record per
`query_name`) so the per-aligner path is robust to malformed external BAMs the way
the consensus path already is. Tracked for AGENT_FIXES.md if adopted.
