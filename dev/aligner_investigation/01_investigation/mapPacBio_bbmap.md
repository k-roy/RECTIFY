# mapPacBio / BBMap (BBTools) — Source-Level Investigation

> **Build note:** code-level claims here were verified against `master`; see `../CORRECTIONS_vs_DRS_BUILD.md` for re-verification vs `origin/drs-validation-rebuild`. **Corrected on the build:** mapPacBio uses `intronlen=10` (not 50) and now sets an explicit `maxindel=max(200000, max_intron)` (`multi_aligner.py:749,754`); the redteam_denovo B5 "no maxindel" concern is RESOLVED. These corrections are applied inline below.

> **Scope.** This report covers **BBMap's `mapPacBio.sh`** (Java class `align2.BBMapPacBio`),
> a component of Brian Bushnell's **BBTools** suite. RECTIFY uses it as a Tier-1 splice-aware
> long-read aligner (`run_map_pacbio()` in `rectify/core/align/multi_aligner.py`, ~line 418),
> where it **wins 18.2 %** of reads in the correct-first consensus — second only to deSALT.

---

## ⚠️ This is BBMap mapPacBio, NOT minimap2's `map-pb` preset

This is a recurring point of confusion, so it is stated first and explicitly:

- **mapPacBio = BBTools / BBMap.** The shell script `mapPacBio.sh` is a thin wrapper that
  launches the Java class **`align2.BBMapPacBio`**. **[FACT — confirmed from the script source:
  the command line invokes `"align2.BBMapPacBio"`.]** It is a *distinct aligner program* with
  its own k-mer index, its own affine-transform scoring matrix, and its own banded global
  alignment engine.
- **minimap2 `-x map-pb`** is a **preset of an entirely different aligner** (Heng Li's minimap2),
  using minimizer seeding + chaining + SIMD extension. It shares nothing with BBMap except the
  marketing word "PacBio."
- The two have **opposite design philosophies**: minimap2 is a heuristic seed-chain-align
  aligner optimized for speed; BBMap mapPacBio is a **slower, more exhaustive seed-and-extend
  aligner that runs a full affine-scored alignment** of every candidate site. That difference is
  the central reason BBMap is RECTIFY's second-best aligner for exact 3′-end / junction placement
  (see *Why It's Second-Best*).

RECTIFY also wraps **vanilla `bbmap.sh`** (`align2.BBMap`) separately in `run_bbmap()` (~line 616)
for short reads — same engine family, different long-read tuning. This report focuses on
mapPacBio but notes the differences.

---

## BBMap Index & Seeding

### K-mer index structure
**[FACT]** BBMap indexes the **reference genome** (not the reads) into a hash-style structure
that maps each k-mer **key → a list of genomic sites** (positions) where that k-mer occurs.
It indexes *all* k-mers of the reference (an "all-kmers" index), so any seed in a read can be
looked up directly against every reference position carrying that exact k-mer.

- **Default k = 13.** **[FACT — the official guide uses `k=13` as the canonical example for
  memory estimation; lower values (`k=12`, `k=8`) are documented as "high"/"super-high"
  sensitivity, higher (`k=14`) for high-precision contaminant removal.]** The BBMap site
  documents k=13 as the standard reference build.
- **Memory ≈ 6 bytes per reference base** in normal mode; a low-memory `usemodulo` mode drops
  to **≈ 3 bytes/base** at a slight sensitivity cost. **[FACT — from the BBMap guide.]**
- The index is built once and **cached on disk** under a `path=` directory. RECTIFY exploits
  this: both `run_map_pacbio()` and `run_bbmap()` point `path=<genome_dir>/bbmap_index` so all
  parallel jobs **share one prebuilt index** (`mpb_index_dir` in the source) rather than
  re-indexing per chunk. **[FACT — RECTIFY source, lines 502–512.]**

### Multi-kmer seed-and-extend
**[FACT]** Bushnell describes BBMap's strategy as a **"multi-kmer-seed-and-extend approach,
analogous to growing polycrystalline silicon."** Rather than relying on a single seed, BBMap
collects many k-mer hits from a read, tallies which reference sites accumulate the most hits
(`minhits` controls the floor), and nucleates alignments outward from those high-density seed
clusters. For long, error-prone reads this multi-seed voting is robust: a single seed may be
destroyed by an error, but dozens of surviving seeds still localize the read.

- **`minratio`** is the master sensitivity knob: a read is "mapped" only if its best alignment
  score reaches at least `minratio × (perfect-match score)`. **mapPacBio defaults to
  `minratio=0.40`** (vs ~0.56 short-read default), i.e. it accepts alignments only ~40 % as good
  as perfect — appropriate for high-error long reads. **[FACT — mapPacBio.sh default is
  `minratio=0.40`; RECTIFY re-asserts `minratio=0.4` explicitly, line 515.]**
- `minid` is an alternative identity-based way to specify the same threshold
  (e.g. `minid=0.9 → minratio≈0.816`). **[FACT — BBMap guide.]**

---

## Base-Level Alignment (why it's accurate)

### Full affine-transform scoring, not a heuristic
**[FACT]** After seeding, BBMap performs a **slow, accurate, banded global-style alignment**
(a Smith-Waterman/Needleman-Wunsch-class dynamic program) of the read against each candidate
reference window, scored with a **custom affine-transform matrix** whose weights were
**"determined empirically through extensive testing."** Documented base weights:

| Event | Score |
|---|---|
| Match | **+100** |
| Single mismatch | **−127** |
| Second consecutive mismatch | additional **−51** |

**[FACT — affine weights quoted from the BBMap guide.]** Penalties differ by mutation class
(substitution vs insertion vs deletion), and consecutive errors are penalized progressively —
an affine gap model. This is a genuine alignment-matrix DP, not the heuristic gap-filling that
fast aligners use.

### Why this yields precise end/indel placement
**[INFERENCE — grounded in the documented affine DP + RECTIFY's empirical win rates]**

- A full affine DP finds the **globally optimal** placement of every indel within the aligned
  window, given the scoring matrix. minimap2's extension is heuristic and chain-anchored: it can
  place an indel a few bases off when a homopolymer or low-complexity run admits several
  near-equal placements, because it stops at the first good-enough chain. BBMap's DP scores the
  whole window and snaps the indel to the single best-scoring column.
- For RECTIFY this matters at the **3′ cleavage-and-polyadenylation site (CPA)**: the exact
  reference base where the terminal exon alignment ends is decided by the DP boundary. A
  globally-scored boundary lands closer to the true CPA than a heuristically-extended one,
  which is why BBMap recovers exact 3′ ends well.
- The **PacBio/Nanopore error model is "dominated by short indels"** and mapPacBio uses a
  **"different error weight profile designed for long reads with a high error rate."**
  **[FACT — BBMap guide; it explicitly calls mapPacBio "the recommended version for Nanopore
  data."]** Tuning the affine penalties toward short indels (rather than substitutions) matches
  the ONT error spectrum, improving boundary precision over a substitution-tuned matrix.

### Modern CIGAR output
**[FACT]** BBMap emits **SAM 1.4 CIGARs with `=` (match) and `X` (mismatch)** rather than
collapsing both to `M`. This explicit match/mismatch encoding is directly useful to RECTIFY's
downstream CIGAR surgery (the correction code reads `=`/`X`/`D`/`N` ops). RECTIFY additionally
runs `samtools calmd -e` (`_apply_calmd_eq`, line 611) after alignment to normalize MD/`=`
encoding. A `sam=1.3` flag exists to force legacy `M` for tools that need it.

---

## Splice / Intron Handling

**[FACT]** BBMap is splice-aware via two cooperating parameters:

- **`maxindel`** — the maximum indel length explicitly searched. **Default 16000** for vanilla
  BBMap. It is a **soft limit**: "indels longer than maxindel … may still be found" unless
  `strictmaxindel` is set. For organisms with long introns (human) the guide recommends
  `maxindel=200k`. RECTIFY's short-read `run_bbmap()` sets `maxindel=100000` to allow
  yeast-scale gaps. **[FACT — BBMap guide + RECTIFY source line 674.]**

> **⚠️ BUILD CORRECTION (vs `drs-validation-rebuild`):** On the build, `run_map_pacbio()` **does**
> set `maxindel` explicitly — `multi_aligner.py:754` `f'maxindel={max(200000, max_intron)}'` (≥200 kb
> for human RNA, scales with `max_intron`). The "no `maxindel` cap on the long-read path" statement
> below (and the redteam_denovo B5 "possibly load-bearing gap") is **stale and now RESOLVED**:
> mapPacBio no longer relies on BBMap's soft ~16000 default for long introns.

- **`intronlen`** — any reference gap (deletion) **≥ `intronlen` is re-encoded as an `N`
  (intron-skip) CIGAR op** instead of `D`. The guide's RNA-seq example uses `intronlen=20`.
  RECTIFY sets **`intronlen=10`** for mapPacBio and **`intronlen=20`** for vanilla
  bbmap (line 673). **[FACT — build `multi_aligner.py:749`; vs master, which used `intronlen=50`.]**

> **⚠️ BUILD CORRECTION (vs `drs-validation-rebuild`):** mapPacBio's D→N threshold on the build is
> **`intronlen=10`** (`multi_aligner.py:749`), not `intronlen=50` as recorded against master. Every
> sentence below that cites `intronlen=50` is corrected to `intronlen=10` — yeast introns (typically
> >50 bp) still reliably become `N` ops; the lower threshold also reclassifies short introns and any
> ≥10 bp deletion-gap into `N`.

The mechanism: BBMap's affine DP can open a long reference gap when the alignment score favors
it (subject to `maxindel`); `intronlen` is purely the **D→N reclassification threshold** applied
afterward. There is no canonical-splice-site (GT-AG) model in the core scorer — splices are
discovered as **long scored deletions**, then relabeled. `xstag` (e.g. `xstag=us`) controls
emission of the XS strand tag for downstream splice-strand tools. RECTIFY does not rely on XS;
its own `junction_refiner` (Module 2H) scans the `N` ops mapPacBio produces and refines the
boundaries against the genome.

**Implication for RECTIFY:** because mapPacBio finds introns as *scored* gaps (full DP), its
junction boundaries are sequence-optimal within the window, complementing deSALT's
splice-graph approach. The `intronlen=10` choice means yeast introns (typically >50 bp,
often 100s of bp) reliably become `N` ops; with the threshold at 10 bp, short introns and any
≥10 bp deletion-gap are also reclassified to `N`, while shorter homopolymer/error deletions stay `D`.

---

## mapPacBio-specific behavior

**[FACT — from `mapPacBio.sh` source + BBMap guide]**

mapPacBio = `align2.BBMapPacBio` with these baked-in defaults (from the script's command line):
```
build=1 overwrite=true minratio=0.40 fastareadlen=6000 ambiguous=best
minscaf=100 startpad=10000 stoppad=10000 midpad=6000
```
plus dynamic heap sizing via `calcXmx` (`freeRam 3200m 84`; initial `-Xmx1g -Xms1g`).

Key differences from vanilla BBMap:

1. **Read-length ceiling = 6 kbp.** Vanilla `align2.BBMap` is tuned for ≤ ~700 bp; mapPacBio
   raises the supported length but the script default **`fastareadlen=6000`** still caps reads
   at ~6 kb. Reads longer than this trigger an internal **AssertionError around 6019 bp**.
   **[FACT — documented in RECTIFY's own docstring, line 430, and consistent with the
   6 kbp ceiling in the guide.]**
2. **Different error-weight profile** for high-error long reads "dominated by short indels"
   — the substantive long-read tuning. **[FACT — BBMap guide.]**
3. **Larger per-thread alignment matrices.** The guide notes the per-thread memory for the
   alignment DP is "relatively small in normal mode, but **bigger in PacBio mode**" — because
   the banded DP window scales with read length. **[FACT.]**
4. **`ambiguous=best`** default — pick the single best site for multi-mappers (RECTIFY relies
   on this for deterministic single-record-per-read behavior).

### Why reads > 6 kb must be split
**[FACT — RECTIFY source]** RECTIFY does **not** trust the internal 6 kb ceiling. It both:
- patches the script default to `fastareadlen=100000` and re-passes `fastareadlen=100000`
  on the command line (line 513, "belt-and-suspenders"), AND
- **pre-splits any read > `MAX_MPB_READ_LENGTH = 6000` bp** into ≤6 kb chunks
  (`split_long_reads`, `mpb_split_reads.py`), aligns the chunks, then **stitches the chunk
  alignments back** into one record (`stitch_split_bam`, line 604). **[FACT — `MAX_MPB_READ_LENGTH
  = 6000` in `mpb_split_reads.py`; split/stitch in `run_map_pacbio`, lines 522–609.]**

This belt-and-suspenders design exists because the 6 kb assertion is a **hard crash**, not a
soft truncation — a single >6 kb read can abort the whole job.

---

## RECTIFY Integration Quirks

These are RECTIFY-side observations (CLAUDE.md + source), not BBMap-intrinsic behavior. The
distinction matters: **the `pt:i:N` tag is NOT produced by BBMap.**

### The `pt:i:N` tag — origin and propagation
**[FACT]** `pt:i:N` is **Dorado's poly-A tail-length tag** (poly-tail; emitted by the ONT
basecaller `dorado` with `--estimate-poly-a`). It is **not** a BBMap/mapPacBio output.

The path by which it pollutes QNAMEs:
1. RECTIFY's DRS pre-trim step extracts reads to FASTQ; if the poly-A length is embedded in the
   FASTQ **header** (e.g. `@UUID\tpt:i:25`, the form `samtools fastq -T pt` produces),
2. **mapPacBio reads the whole header verbatim** and, after `samtools sort` converts the
   forbidden whitespace to `_`, the **BAM QNAME becomes `UUID_pt:i:25`** (or `UUID pt:i:25`
   pre-sort, space-separated). **[FACT — CLAUDE.md DRS notes + RECTIFY source.]**
3. RECTIFY therefore **strips the `_pt:i:N` suffix** from every QNAME after alignment
   (lines 578–598, the "Read ID Purity Policy") so downstream parquet lookups keyed by bare UUID
   work. This is why CLAUDE.md repeatedly warns: **do NOT use `samtools fastq -T pt`** — keep
   the FASTQ read name a bare UUID; the poly-A length lives in the parquet metadata instead.

So mapPacBio's only "fault" here is **faithfully copying the header into QNAME** — a documented
BBMap behavior (`trd`/`trimreaddescriptions` exists precisely to trim names at whitespace).

### Unmapped duplicate records
**[FACT — CLAUDE.md v3.1.8 note]** mapPacBio **emits unmapped (FLAG 0x4) copies of reads
alongside the primary alignment.** These records pass the usual
`not is_secondary and not is_supplementary` filter (because unmapped reads are neither), so
RECTIFY must add an explicit **`read.is_unmapped` guard** to skip them. If ignored, the same
read appears twice (once aligned, once unmapped), corrupting per-read consensus.

**[INFERENCE]** This is consistent with BBMap's `out=` stream behavior: by default BBMap writes
*all* reads (mapped + unmapped) to the main output unless `outm=`/`outu=` are used to split them.
RECTIFY pipes everything to one SAM, so unmapped reads are interleaved and must be filtered.

### Other quirks
- **Name-sorted output** (`samtools sort -n`) so consensus selection can stream across aligners
  without a second sort (lines 543–559). The space→underscore QNAME conversion happens *here*.
- **`-Xmx32g` heap cap** (line 516) to stop BBMap's JVM from OOM-killing neighboring SLURM tasks
  — BBMap will grab a large fraction of free RAM by default via `calcXmx`.
- **gapmm2 contrast:** unlike gapmm2, mapPacBio does *not* need a cleaned/deduplicated FASTQ for
  the alignment itself, but it *does* import duplicate-UUID and empty-seq problems from DRS data
  in the same way; RECTIFY's pt-strip + unmapped-guard handle the mapPacBio-specific fallout.

---

## Why It's Second-Best (accuracy mechanisms)

RECTIFY correct-first win rates: **deSALT 78.9 %, mapPacBio 18.2 %**, uLTRA 2 %, gapmm2 0.8 %,
minimap2 0.1 %. **[FACT — CLAUDE.md.]** mapPacBio is a clear #2. Mechanistic reasons:

1. **Full affine-scored global alignment of every candidate window.** mapPacBio does not stop at
   a heuristic chain; it runs a real DP and reports the column-optimal boundary. For 3′-end (CPA)
   recovery — where RECTIFY's whole purpose is the exact terminal base — a globally optimal
   alignment boundary beats a heuristically extended one. **[INFERENCE from documented DP +
   empirical win rate.]**
2. **Long-read-tuned, indel-dominated error model** matched to ONT's error spectrum, so indel
   placement near the 3′ end is biased correctly. **[FACT — guide; INFERENCE for CPA impact.]**
3. **Multi-kmer voting** localizes noisy reads robustly, so it finds *a* good seed cluster even
   where a single-minimizer scheme would miss, then refines with DP. **[FACT seeding + INFERENCE.]**
4. **Explicit `=`/`X` CIGARs + `N` intron reclassification** give RECTIFY's `junction_refiner`
   and CPA walk-back clean, unambiguous ops to operate on. **[FACT.]**

Why it loses to deSALT: deSALT uses a **reference splice graph / de Bruijn index** that models
known and concordant junctions explicitly, so on spliced reads it places junctions with graph
support that BBMap (which discovers splices only as scored gaps, with no GT-AG model) cannot
match. BBMap's edge is on **unspliced terminal exons and exact 3′ ends**, where its DP boundary
is hard to beat — hence a strong but second-place 18.2 %.

---

## Strengths

- **[FACT]** Full affine-transform DP alignment → **column-optimal indel and end placement**;
  the single biggest reason for exact CPA recovery.
- **[FACT]** Splice-aware with tunable D→N threshold (`intronlen`) and large `maxindel`; no
  reliance on a fixed splice motif, so it finds non-canonical introns as scored gaps.
- **[FACT]** Modern `=`/`X` CIGAR output — directly consumable by RECTIFY's CIGAR surgery.
- **[FACT]** Long-read error profile explicitly tuned for ONT/PacBio ("recommended for Nanopore").
- **[FACT]** Disk-cached, shareable reference index (one build, many parallel jobs).
- **[FACT]** Fast indexing relative to FM-index aligners; `ambiguous=best` gives deterministic
  single-record output.

## Weaknesses

- **[FACT]** **6 kbp read ceiling** (`fastareadlen=6000`, hard AssertionError ~6019 bp) →
  RECTIFY must pre-split and re-stitch long reads. A genuine operational burden.
- **[FACT]** **RAM-hungry**: ~6 bytes/ref base for the index plus larger per-thread DP matrices
  in PacBio mode; JVM grabs most free RAM unless `-Xmx` capped → can OOM-kill neighbors on HPC.
- **[FACT]** **Java/JVM startup + heap overhead** and a slower, exhaustive DP make it
  substantially slower than minimap2; acceptable only because RECTIFY chunks and parallelizes.
- **[INFERENCE]** **No GT-AG / splice-graph model** — junctions are discovered as scored
  deletions, so on heavily spliced isoforms it underperforms graph-aware deSALT/uLTRA.
- **[FACT/INFERENCE]** **Homopolymers & very long introns**: long-read homopolymer indel noise
  can still admit multiple near-equal DP placements (RECTIFY adds HP-aware correction precisely
  because no aligner, including BBMap, nails every homopolymer boundary). Very long introns
  approach/exceed `maxindel` and may be missed or mis-scored — but on the build this is mitigated:
  RECTIFY sets `maxindel=max(200000, max_intron)` explicitly (`multi_aligner.py:754`), so the cap
  scales to the organism rather than relying on BBMap's soft ~16000 default.
- **[FACT]** **Operational quirks imported from DRS data**: propagates `pt:i:N` into QNAME and
  emits unmapped duplicate records — both require RECTIFY-side post-processing.

---

## Source References

**Primary (BBMap / Bushnell):**
- Bushnell B. *BBMap: A Fast, Accurate, Splice-Aware Aligner.* 9th Annual Genomics of Energy &
  Environment Meeting, March 2014. OSTI 1241166 / eScholarship `1h3515gn`. (Abstract: "superior
  speed, sensitivity, and specificity to … bowtie2, bwa, smalt, GSNAP, and BLASR.")
  https://www.osti.gov/biblio/1241166 · https://escholarship.org/uc/item/1h3515gn
- BBMap Guide (`docs/guides/BBMapGuide.txt`): k-mer index, affine weights (+100/−127/−51),
  `maxindel` default 16000, `intronlen`, `minratio`, `usemodulo` ~3 B/base, PacBio mode notes.
  https://github.com/BioInfoTools/BBMap/blob/master/docs/guides/BBMapGuide.txt
- `mapPacBio.sh` source — class `align2.BBMapPacBio`; defaults `minratio=0.40 fastareadlen=6000
  ambiguous=best minscaf=100 startpad/stoppad=10000 midpad=6000`; `calcXmx`/`freeRam`.
  https://github.com/BioInfoTools/BBMap/blob/master/sh/mapPacBio.sh
- Official BBTools site: https://bbmap.org/ · https://bbmap.org/tools/bbmap
- Ubuntu/Debian manpages for `bbmap.sh` (k-mer, splice-aware, parameters).

**`pt:i:` tag origin (NOT BBMap):**
- Dorado basecaller — `pt:i` poly-A tail tag emitted with `--estimate-poly-a`.
  https://github.com/nanoporetech/dorado (CHANGELOG, issues #593/#1090/#1562).

**Long-read splice-alignment benchmarking context:**
- "Evaluation of tools for long read RNA-seq splice-aware alignment," PMC6192213.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6192213/

**RECTIFY usage (ground truth for integration):**
- `/home/user/RECTIFY/rectify/core/align/multi_aligner.py` — `run_map_pacbio()` (≈418–613),
  `run_bbmap()` (≈616+), `_apply_calmd_eq()` (≈400).
- `/home/user/RECTIFY/rectify/core/align/mpb_split_reads.py` — `MAX_MPB_READ_LENGTH = 6000`,
  `split_long_reads`, `stitch_split_bam`.
- `/home/user/RECTIFY/CLAUDE.md` — `pt:i:N` QNAME notes, unmapped-duplicate guard (v3.1.8),
  correct-first win rates (deSALT 78.9 %, mapPacBio 18.2 %).

---

*FACT vs INFERENCE is labeled inline. "FACT" = stated in BBMap source/docs or RECTIFY source.
"INFERENCE" = mechanistic reasoning connecting documented behavior to RECTIFY's observed results.*
