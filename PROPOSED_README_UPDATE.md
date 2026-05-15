# Proposed README Update — Handoff Doc

**Filed:** 2026-05-14  
**For:** agent implementing README + ARCHITECTURE doc updates  
**Context:** RECTIFY now supports three sequencing protocol tracks. The current README describes only the DRS (direct RNA-seq) track. This doc specifies what to add.

---

## Background

RECTIFY has gained two new protocol tracks since the last README pass:

1. **ONT PCR-cDNA (PCB114.24)** — `rectify correct-cdna` subcommand. Full UMI-aware per-molecule consensus pipeline for SSP-primed PCB114 chemistry. Produces `stage1_consensus.bam` with a rich BAM tag schema, plus `clusters.tsv`, `isoforms.tsv`, and `t1t2_pairs.tsv`.

2. **QuantSeq 3' REV** — `rectify correct --dT-primed-cDNA --short-read`. Antisense short-read chemistry; gene/RNA strand is the inverse of BAM strand. The read-vs-reference walkback (`walkback.py`) replaces the old reference-only A-tract approach, fixing internal-priming misassignment at convergent loci.

All three tracks share the same pipeline skeleton, though the point at which the two aligners compete differs:

- **Poly-A pre-trim** — strips basecalled poly-A tail before alignment
- **Per-chromosome chunking** — fan-out to `--region chrX` jobs; merge after
- **Multi-aligner consensus** — two aligners per track; consensus selection picks the best junction set per read
  - *DRS / QuantSeq REV:* both aligners run on the **raw input reads** in parallel
  - *ONT cDNA:* input BAM is pre-aligned (minimap2); UMI deduplication + abPOA runs first; both aligners then compete on the **per-cluster consensus reads**
- **Read-vs-reference walkback** — protocol-agnostic core; only the 3'-end side and strand mapping differ

### Current vs target aligner sets

| Track | Current aligners | Target aligners | Optional |
|---|---|---|---|
| DRS | minimap2 + mapPacBio + gapmm2 | minimap2 + mapPacBio + gapmm2 | deSALT, uLTRA |
| ONT cDNA (PCB114.24) | minimap2 only | minimap2 + **mapPacBio (mpb)** | deSALT, uLTRA |
| QuantSeq REV | single aligner | **BWA + BBMap** | — |

**Gap — ONT cDNA:** `correct-cdna` currently runs minimap2 only. Adding mapPacBio (mpb) would bring it to parity with DRS and enable the same consensus quality improvement. The multi-aligner consensus module (`core/chimeric_consensus.py`) is already protocol-agnostic — the cDNA path needs to call it.

**Gap — QuantSeq REV:** No multi-aligner consensus in the current short-read path. BWA-MEM and BBMap both handle short spliced reads well; BBMap's splicing sensitivity is particularly strong for yeast introns. Both aligners produce SAM/BAM that the existing consensus selector can consume.

---

## Changes to make in README.md

### 1. Replace the pipeline prose (after "Bundled data for yeast") with a three-track overview table

```markdown
## Protocol Tracks

RECTIFY supports three sequencing chemistries. All share poly-A pre-trim and
per-chromosome chunking; they differ in how the 3' end is identified and
which strand convention applies.

| Track | Chemistry | Subcommand | Aligners | 3'-end side | Gene strand |
|---|---|---|---|---|---|
| **DRS** | ONT direct RNA-seq | `rectify correct --direct-rna` | minimap2 + mapPacBio + gapmm2 (± deSALT, uLTRA) | right | = BAM strand |
| **ONT cDNA** | PCB114.24 SSP+UMI | `rectify correct-cdna` | minimap2 + mapPacBio* | right | = BAM strand |
| **QuantSeq REV** | dT-primed short-read | `rectify correct --dT-primed-cDNA --short-read` | BWA + BBMap* | left | opposite of BAM strand |

\* = target state; not yet implemented (see Implementation gaps below)
```

### 2. Update the Key Features table — add these rows

| Feature | Description |
|---|---|
| **ONT cDNA UMI Consensus** | PCB114.24 SSP+UMI architecture: three read classes (Type-1a SSP-first, Type-1b pA-first UMI-captured, Type-2 pA-first truncated); directional UMI clustering (Levenshtein ≤ 3, 2× rule); two-stage strand-split abPOA consensus per molecule |
| **5'/3' End Isoform Clustering** | 5' UMI-anchored reads cluster by both TSS and CPA; 3' pA-tail-anchored reads cluster by CPA only; same-molecule strand pairing links unpaired Type-2 clusters to their Type-1 partner via `XL:Z` tag |
| **Read-vs-Reference Walkback** | Compares read bases to reference at each position from the 3' end inward; walks past every A (genomic or not) until a non-A match is found — fixes internal-priming into genomic A-tracts that the old reference-only approach missed |
| **NMD-AS Splice Audit** | Per-read splice classification (unspliced / annotated / alternative / novel) against a gene panel; quantifies NMD-sensitive isoform accumulation across conditions |

### 3. Add a new "How It Works" section after the existing ones

````markdown
### ONT PCR-cDNA UMI Consensus (PCB114.24)

PCB114.24 library prep produces a dsDNA cDNA molecule where one physical strand
carries `SSP → UMI → GGG bridge → mRNA body → poly-A` and the complementary
strand carries the reverse complement. Because either strand can be sequenced,
three read classes arise:

| Class | How identified | UMI available | orient tag |
|-------|---------------|---------------|------------|
| **Type-1a** (5' UMI-anchored, SSP-first) | SSP_FWD found in read | yes — upstream of SSP | `fwd` |
| **Type-1b** (5' UMI-anchored, pA-first, full-length) | SSP_RC found near read end | yes — downstream of SSP_RC (rev-comp) | `rev` |
| **Type-2** (3' pA-tail-anchored, truncated) | No SSP detected; poly-A/adapter anchor only | no | `fwd` or `rev` |

Type-1a and Type-1b both carry the UMI and land in the same UMI cluster.
Type-2 reads did not travel far enough to capture the UMI and are carried
forward as independent observations grouped by CPA position for isoform
counting only — not deduplicated.

> **Note for implementing agents:** Type-1a and Type-1b must be added as
> explicit code-level distinctions so the code matches the README terminology.
> The implementation is a new `read_subtype: str` field on `ReadInfo`, set in
> `parse_read_info()`:
> - `read_type=1, orient="fwd"` → `read_subtype="1a"`
> - `read_type=1, orient="rev"` → `read_subtype="1b"`
> - `read_type=2` → `read_subtype="2"`
>
> The existing `read_type` integer (1 or 2) is unchanged — clustering logic
> uses it as before (both subtypes route identically through UMI clustering).
> `read_subtype` adds visibility for isoform analysis and should be emitted
> as a BAM tag (`XY:Z`) on output records. Update `clusters.tsv` and
> `isoforms.tsv` to include a `read_subtype` column.

**Key asymmetry:** all reads in a UMI cluster (Type-1a + Type-1b) can be
merged regardless of 3' truncation status — the UMI identifies the molecule.
Type-2 reads are not deduplicated at all: each read is treated as an
independent molecule observation. Multiple molecules commonly share the same
CPA site, and singleton UMI counts in Type-1 data confirm that PCR duplication
rates are low — so collapsing Type-2 reads by position would discard real
signal, not duplicates.

<p align="center">
  <img src="docs/figures/cdna_umi_architecture.png" alt="PCB114.24 read classes" width="700">
</p>

#### Two-stage abPOA consensus within a UMI cluster

Within each UMI cluster the consensus is built in two sequential steps to
cancel strand-specific error profiles:

**Step 1 — Within-strand sub-consensus:** reads are split by `is_reverse`
(BAM flag). Top-strand reads (is_reverse=True, typically Type-1b) form one
sub-pool; bottom-strand reads (is_reverse=False, typically Type-1a) form
another. abPOA builds a consensus from each sub-pool independently.

**Step 2 — Cross-strand merge:** the two per-strand sub-consensuses are
fed into a second abPOA pass to produce the final per-cluster consensus.
Merging at this stage instead of pooling all reads first prevents
systematic strand-specific basecaller biases from skewing the consensus.

Falls back to single-pass POA when one strand has fewer than two reads.

<p align="center">
  <img src="docs/figures/cdna_poa_consensus.png" alt="Two-stage strand-split abPOA" width="680">
</p>

#### Same-molecule strand pairing (Type-1 ↔ Type-2)

A molecule that produced a Type-1 cluster may also have generated Type-2
reads (pA-first reads that were truncated before reaching the UMI). These
are independent observations with no Stage-1 deduplication; they are simply
carried forward as individual records until this pairing step.
RECTIFY links them post-hoc: within each (gene, orientation) group, each
Type-2 cluster is matched to the nearest Type-1 cluster whose 5' and 3'
termini are both within tolerance (`|Δ5'| ≤ 5 bp`, `|Δ3'| ≤ 5 bp`).
Matched pairs receive an `XL:Z` partner tag on both records and are
reported in `t1t2_pairs.tsv`. A paired CPA call (confirmed from both
strands of the same molecule) carries higher confidence than a single-strand
call alone.

<p align="center">
  <img src="docs/figures/cdna_isoform_clustering.png" alt="Isoform clustering and strand pairing" width="700">
</p>

Key output tags on `stage1_consensus.bam`:

| Tag | Meaning |
|-----|---------|
| `XU:Z` | UMI sequence |
| `XC:i` | UMI cluster ID |
| `XR:i` | Number of raw reads merged |
| `XF:i` | POA consensus flag (1=pileup, 2=abPOA) |
| `XA:i` | Poly-A tail length (basecalled A's after alignment end) |
| `XG:Z` | Gene/RNA strand (`+`/`-`) |
| `XS:Z` | Sense/antisense classification |
| `XL:Z` | Partner cluster ID (same-molecule strand pairing) |

### Read-vs-Reference Poly-A Walkback

The core algorithmic improvement over earlier RECTIFY versions. Three-case
terminal gate at the 3' end of each read:

1. **Anchored non-A base matching genome** → already correctly placed, no correction.
2. **Read base is A, genome base is also A** → this is the internal-priming scenario; walk inward (old approach incorrectly stopped here).
3. **Mismatch** → walk inward until `read_base == ref_base` and `ref_base != 'A'`.

This protocol-agnostic logic is implemented in `rectify/core/correct/walkback.py`
and wrapped by per-protocol modules:

| Module | Chemistry | is_reverse → gene strand | 3'-end side |
|--------|-----------|--------------------------|-------------|
| `walkback_drs` | ONT DRS | False→+, True→− | right |
| `walkback_quantseq_rev` | QuantSeq REV | False→−, True→+ | left |

<p align="center">
  <img src="docs/figures/walkback_readvsref.png" alt="Read-vs-reference walkback" width="680">
</p>

### Splice Classification

`rectify analyze splice` classifies each read's N-cigar ops against a set of
annotated introns using strict 0-bp matching (with optional snap tolerance for
basecaller jitter):

| Class | Condition |
|-------|-----------|
| `no_intron_span` | Read does not overlap any annotated intron |
| `unspliced` | Read spans intron but has no N-cigar op |
| `annotated` | All N-ops match annotated donor+acceptor exactly |
| `alternative` | One end of an N-op matches an annotated site |
| `novel` | Neither end of any N-op matches any annotated site |

<p align="center">
  <img src="docs/figures/splice_classification.png" alt="Splice classification decision tree" width="620">
</p>
````

### 4. Update the Commands table — add

| Command | Purpose |
|---------|---------|
| `correct-cdna` | PCB114.24 ONT cDNA: UMI extraction → directional clustering → abPOA consensus → 5'/3' isoform calling |
| `analyze splice` | Per-read splice classification against a gene panel (unspliced/annotated/alternative/novel) |

### 5. Update "Supported Technologies"

Add to the existing list:
- ONT PCR-cDNA (PCB114.24, SSP+UMI chemistry)
- QuantSeq 3' REV (dT-primed antisense short-read)

---

## Changes to make in ARCHITECTURE.md

### 1. Replace the Pipeline Overview (Step 0–5) with a three-track diagram

```
Shared pre-processing
  └── Poly-A pre-trim (all tracks)
  └── Per-chromosome chunking (all tracks)

DRS track
  Step 0: Poly-A pre-trim
  Step 1: Multi-aligner alignment (minimap2 / mapPacBio / gapmm2)
  Step 2: Per-aligner correction (walkback_drs)
  Step 3: Consensus selection
  Step 4: Restore poly-A tail
  Step 5: Multi-sample downstream analysis

ONT cDNA track (PCB114.24)
  Input: pre-aligned BAM (minimap2, supplied by user)
  Step 0: SSP/UMI extraction + poly-A pre-trim
  Step 1: Per-chrom UMI directional clustering
  Step 2: abPOA consensus per UMI cluster  (one read per molecule)
  Step 3a: re-align consensus with minimap2
  Step 3b: re-align consensus with mapPacBio (mpb)  [target; not yet in code]
  Step 3c: chimeric_consensus.select_best(minimap2, mpb)
  Step 4: Type-1/Type-2 isoform clustering + T1↔T2 reconciliation
  Step 5: Multi-sample downstream analysis (shared with DRS)

QuantSeq REV track
  Step 0: Poly-A pre-trim (same module as DRS)
  Step 1: Per-chrom alignment
  Step 2: walkback_quantseq_rev (inverted strand convention)
  Step 3: 3'-end clustering (shared with DRS analyze)
  Step 4: Multi-sample downstream analysis (shared)
```

### 2. Layer 2 — add `cdna_correct_command.py`

### 3. Layer 4 — add new correction modules

```
core/correct/
  walkback.py          — protocol-agnostic walkback core + walkback_drs()
                         wrapper (both live in this single file; there is
                         no separate protocols/drs.py)
  protocols/
    quantseq_rev.py    — walkback_quantseq_rev wrapper (strand inverted;
                         planned — see HANDOFF_QUANTSEQ_REV_WALKBACK.md)
```

### 4. Layer 5 — add splice analysis

```
core/analyze/
  splice_summary.py    — GeneIntronSet, classify_read, classify_bam_for_genes
```

### 5. Update Mermaid call graph

Add a parallel branch from `cli.py` to `cdna_correct_command.py` alongside the existing `correct_command.py` branch. Show `walkback.py` called from both DRS and QuantSeq REV correction paths.

---

## New figures to generate (docs/figures/)

All figures should follow `FIGURE_STYLE_GUIDE.md`. Use `generate_figures.py` as the template.

**Existence status:** `cdna_umi_architecture`, `cdna_isoform_clustering`,
`walkback_readvsref`, and `splice_classification` already exist in
`docs/figures/` (created in a prior session). Only `cdna_poa_consensus`
needs to be built from scratch. `cdna_isoform_clustering` requires a full
redraw — the existing file has Type-1 and Type-2 biology inverted.

| Filename | Content | Notes |
|----------|---------|-------|
| `cdna_umi_architecture.png/.svg` | Show the PCB114.24 dsDNA molecule once (horizontal, with SSP/UMI/GGG/mRNA/pA labelled on one strand and the RC on the other). Below it, show three read classes as arrows: Type-1a (SSP-first, full-length, orient=fwd), Type-1b (pA-first, reaches UMI at far end, orient=rev), Type-2 (pA-first, truncated within exon, no UMI). Annotate which XU/XA/XF/XG/XS tags each class emits. | Three-row read diagram below a single dsDNA schematic, ~700px wide |
| `cdna_poa_consensus.png/.svg` | Two-stage abPOA consensus within a UMI cluster. Input: a UMI cluster containing Type-1a reads (blue, is_reverse=False) and Type-1b reads (orange, is_reverse=True). Step 1: split by strand into two sub-pools, abPOA each → sub-consensus A (blue) and sub-consensus B (orange). Step 2: feed both sub-consensuses into second abPOA pass → final consensus (green). Show singleton/fallback path (fewer than 2 reads per strand → single-pass POA). | Vertical flow diagram, ~600px wide. NEW figure — not in current figure set. |
| `cdna_isoform_clustering.png/.svg` | **Complete redraw — current figure has Type-1 and Type-2 biology inverted.** Panel 1 (5' UMI-anchored, Type-1): show multiple reads sharing the same 5' TSS (UMI anchor) but with variable 3' ends (some truncated before CPA). Label: all merged by UMI regardless of 3' truncation. Panel 2 (3' pA-tail-anchored, Type-2): show reads with variable 5' truncation points converging on the same 3' CPA. Label: cannot merge across different 5' starts — may be different molecules. Panel 3 (same-molecule strand pairing): a Type-1 cluster and a Type-2 cluster with matching TSS + CPA (within 5 bp tolerance) are linked by XL:Z. | Three-panel layout. Replace the existing figure entirely. |
| `walkback_readvsref.png/.svg` | Side-by-side: old reference-only A-tract (stops at genomic A — wrong) vs new read-vs-reference three-case gate (walks through — correct). Highlight the internal-priming scenario. | Two-panel layout |
| `splice_classification.png/.svg` | Decision tree flowchart: no_span → unspliced → annotated → alternative → novel | Vertical flowchart |

---

## Implementation gaps (items for a follow-up agent)

These are not doc changes — they are code changes that the docs will reflect once done.

### 1. Add consensus pre-trimming before multi-aligner step (cDNA-specific, NOT yet implemented)

The abPOA consensus sequence currently passes directly to the aligner with the SSP/UMI/GGG bridge and poly-A tail still attached. Aligners soft-clip these, but that leaves junction calling at the mercy of aligner-specific soft-clip heuristics. The DRS track avoids this with an explicit poly-A pre-trim before alignment. The cDNA track needs two analogous steps:

**Step A — 5' SSP/UMI/GGG strip (cDNA-specific):**
The SSP + UMI (27 nt) + GGG bridge at the 5' of every Type-1a consensus (and SSP_RC + UMI_RC + CCC at the 3' end of Type-1b) is non-genomic. The exact trim boundary is already known from UMI extraction — it is the GGG/CCC bridge end position in the query sequence. Strip this prefix/suffix from the consensus string before passing to any aligner.

**Step B — 3' poly-A pre-trim (same module as DRS):**
Strip the basecalled poly-A tail and any trailing adapter sequence from the consensus 3' end before alignment. Record the trimmed length so it can be restored as the `XA:i` tag after alignment. Use the same pre-trim module (`core/polya_trim.py` or equivalent) already used in the DRS track.

```
abPOA consensus (raw: SSP-UMI-GGG-mRNA-pA-adapter)
  ↓
Step A: strip SSP/UMI/GGG from 5' (and adapter_RC/pT from 5' for orient=rev)
Step B: strip poly-A + adapter from 3' (record trim length → XA)
  ↓ (clean mRNA sequence, ready for multi-aligner)
  ├── minimap2      (already present, applied to trimmed consensus)
  ├── mapPacBio     ← TO ADD
  └── gapmm2        ← TO ADD
  ↓
chimeric_consensus.select_best(...)
  ↓
Restore poly-A tail annotation (XA tag from Step B trim length)
  ↓
walkback + 5' TSS walk-forward + output tags
```

**File:** `rectify/core/cdna_correct_command.py` — insert trimming before `realign_consensus()` call (around line 1229). The UMI/GGG boundary position is available from the `ReadInfo.pos5_corrected` field computed during extraction.  
**Why:** without pre-trim, all three aligners see non-genomic sequence at both ends. This degrades junction calls near the TSS and CPA site — exactly where accuracy matters most.  
**Test:** add a fixture where the raw consensus has a 5-nt SSP stub + 8-nt poly-A; verify the trimmed sequence passed to the aligner matches the expected mRNA-only string.

### 2. Add mapPacBio + gapmm2 to `correct-cdna` multi-aligner step (NOT yet implemented)

Once pre-trimming is in place (gap 1 above), add the two additional aligners so the cDNA track matches DRS. Re-alignment applies to **all** consensus reads — singletons (XR=1) benefit from the multi-aligner splice model just as much as polished clusters.

**File:** `rectify/core/cdna_correct_command.py` — after the pre-trim step, replace the single `realign_consensus(poa_seq, mp_aligner)` call with parallel calls to minimap2 (`realign_consensus`), mapPacBio (`multi_aligner.run_map_pacbio`), and gapmm2 (`multi_aligner.run_gapmm2`), then feed all three to `chimeric_consensus.select_best`.  
**Why:** mapPacBio's splice model handles GC-poor yeast splice sites; gapmm2 adds terminal-exon homopolymer refinement. Same rationale as the DRS triple-aligner default.  
**Test:** existing `tests/test_cdna_correct_command.py` pattern; add a synthetic consensus read where minimap2 miscalls a junction and mapPacBio or gapmm2 rescues it.

### 2. Add BWA + BBMap to the QuantSeq REV short-read path

**File:** `rectify/core/correct_command.py` (the `--short-read --dT-primed-cDNA` branch, around lines 354–369)  
**What:** run BWA-MEM and BBMap in parallel (already have `core/multi_aligner.py`); feed both alignments to `chimeric_consensus.select_best`; then apply `walkback_quantseq_rev`.  
**Why:** short-read splice sites are harder to call from a single aligner; two-aligner consensus catches the ~15% of reads where one aligner miscalls the junction.  
**Aligners bundled?** BWA binary is in `rectify/data/bin/linux_x86_64/` (check). BBMap is a JVM tool — may need a `subprocess` call to the user's BBMap install, similar to how deSALT is handled.  
**Test:** mirror `tests/test_walkback_readvsref.py::TestQuantSeqRevWrapper` — add a two-aligner consensus fixture.

### 3. Update `ARCHITECTURE.md` pipeline diagram

Once both gaps above are closed, update the Mermaid call graph to show the three parallel aligner branches converging at `chimeric_consensus.select_best` before the walkback step. The three-track overview table in the README should then drop the asterisk footnote.

### 4. Fix stale docstring in `cdna_correct_command.py`

The v1.15 changelog block (around line 27) describes Type-2 reads as using
"position-based deduplication." This is the same stale framing corrected
above: Type-2 reads are not deduplicated. Update to read "independent
observations grouped by CPA position for isoform counting."

---

## Files NOT to touch

- Existing 9 algorithm figures (DRS-specific, still accurate)
- `CHANGELOG.md`, `TODO.md`, `PROJECT_STATUS.md` (internal)
- `mkdocs.yml` nav (separate pass once content is settled)

---

## Reference implementations to read before editing

- `rectify/core/correct/walkback.py` — the walkback core
- `rectify/core/cdna_correct_command.py` — full cDNA pipeline orchestrator
- `rectify/core/analyze/splice_summary.py` — splice classifier
- `tools/nmd_as_splice_audit.py` — splice audit tool using the classifier
- `tests/test_walkback_readvsref.py` — walkback test patterns
- `tests/test_splice_summary.py` — splice classifier test patterns
- `HANDOFF_QUANTSEQ_REV_WALKBACK.md` — QuantSeq REV wrapper design (not yet wired into CLI)
