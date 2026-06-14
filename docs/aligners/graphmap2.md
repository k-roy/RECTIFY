# Graphmap2 — candidate aligner (spaced-seed graph region-voting)

**Status (2026-06-14): UNDER EVALUATION, smoke-test gated.** Algorithmically
orthogonal in its **seeding**, splice-usable, but **STALE (last release Feb 2020)**
— adoption is gated on a compile + RNA004-DRS smoke-test confirming it still emits
valid N-op alignments on current chemistry (Sherlock job 29546795).

## Why it's interesting (and where its orthogonality actually is)

Graphmap2's orthogonality is in its **region selection / seeding**, NOT its splice
model:

- **Seeding/region selection (NOVEL):** gapped-spaced-seed **graph region voting**
  → LCSk anchor chaining → **L1 linear regression (Hough-transform-like) anchor
  fit** to a diagonal → **knapsack** anchor filtering. None of the other panel
  members do spaced-seed graph voting + L1-regression anchor fitting. This is the
  orthogonality case.
- **Splice DP (NOT novel):** Graphmap2 **vendors minimap2's `ksw_exts2_sse`
  spliced kernel verbatim**. So its intron-scoring is *identical* to minimap2's,
  not a distinct splice model. Splice-model diversity must come from elsewhere
  (GMAP/uLTRA/deSALT/BBMap/Magic-BLAST).

## GT-AG motif handling — verified at source (2026-06-14)

This was scrutinized because RECTIFY is a **motif-agnostic non-canonical-junction
discovery pipeline** ([[feedback_no_motifs_unbiased_discovery]]) and the word
"motif" in a splice aligner is a red flag. Source-level finding:

- **Intron EXISTENCE is decided motif-agnostically** — from anchor-gap geometry
  (region-selection → knapsack exon assembly), then any deletion ≥10 bp is
  reclassified `D→N` with **no motif check** (`process_read.cc:928,1036`,
  `MIN_INTRON_LEN=10`).
- **GT-AG is only a finite SOFT bonus** inside the vendored DP: non-canonical
  boundaries cost `noncan = 9` points (additive, **never `-INF`/hard block**;
  `ksw2_exts2_sse.cc:97-112`). Non-canonical introns ARE emitted given enough
  flanking evidence.
- **`noncan = 9` is the IDENTICAL constant minimap2 `-x splice` uses** (same
  vendored kernel). So Graphmap2 is **exactly as motif-biased as the minimap2 we
  already run and accept — not more.** It is *not* disqualified, and it is *not*
  uniquely unbiased.
- **Caveat — less tunable:** minimap2 has `-u n` to switch GT-AG off entirely;
  Graphmap2 hardcodes `noncan=9` (and `flag=1600`) with no off-switch
  (`aligner_ksw2.cc:202-207`). The only splice CLI knob is the boolean `--spliced`.

Conclusion: usable in our motif-agnostic pipeline at the same canonical-bias level
as minimap2. "Annotation-agnostic (no GTF) ≠ motif-agnostic" — Graphmap2 is the
former, and carries minimap2's mild soft GT-AG bonus. (An earlier note calling it
the "best motif-agnostic fit" was wrong and has been corrected.)

## ONT suitability (evidence)

- Built specifically for ONT (Marić, Sović, Šikić et al.; bioRxiv 720458). `-x
  rnaseq` mode. Authors report higher correct-alignment / exon-boundary recall
  than minimap2 and GMAP on their ONT datasets; slower than minimap2, faster than
  GMAP. Numbers predate RNA004 — no current-chemistry DRS numbers exist.
- The paper PDF returned 403 to automated fetch; %-non-canonical-junction numbers
  could not be retrieved (the mechanism is structurally capable of emitting them).

## Invocation

```bash
# rnaseq mode = spliced; index auto-builds next to the (writable) reference
graphmap2 align -r GRCh38_chr5.fa -d reads.fastq -o out.sam -x rnaseq -t <threads>
```

- Output: SAM with **N-op introns**. Index files are written next to the
  reference → use a **writable copy** of the reference (the shared
  `sumner_lab/references/` dir may be read-only).

## Install / license / maintenance / RISK

- bioconda: `conda install -c bioconda graphmap` (v0.6.4 line = the graphmap2-era
  build; original `isovic/graphmap` is deprecated → `lbcb-sci/graphmap2`).
- License: MIT. Linux-64 + aarch64; macOS build older (0.5.2).
- **RISK: last release Feb 2020 (~6 yr stale).** This is the headline risk —
  the smoke-test must confirm it (a) compiles/runs in the current env and (b)
  produces sane N-op alignments on RNA004 reads before any integration.

## Sources

- Marić, Sović, Šikić et al., Graphmap2, bioRxiv [720458](https://www.biorxiv.org/content/10.1101/720458v1.full).
- Source: [github.com/lbcb-sci/graphmap2](https://github.com/lbcb-sci/graphmap2) — GT-AG
  mechanism verified in `src/ksw2/ksw2_exts2_sse.cc`, `src/aligner/aligner_ksw2.cc`,
  `src/graphmap/process_read.cc`, `src/program_parameters.cc` (2026-06-14 source read).
- 2026-06-14 orthogonality panel: memory `project-aligner-orthogonality-panel`.
