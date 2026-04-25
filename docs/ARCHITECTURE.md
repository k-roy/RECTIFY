# RECTIFY Architecture

This document describes the full processing pipeline from raw DRS BAM to
corrected 3' and 5' end coordinates, with emphasis on the pre-consensus
scoring pass and the correction modules in `rectify correct`.

---

## End-to-End Pipeline

```
Dorado BAM (DRS)
      │
      ▼  rectify trim-polya  (Step 0 — DRS only)
      │  ┌──────────────────────────────────────────────┐
      │  │ Three-pass poly(A) + adapter removal          │
      │  │ Output: unaligned BAM  +  metadata.parquet   │
      │  └──────────────────────────────────────────────┘
      │
      ▼  samtools fastq  →  trimmed.fastq.gz
      │
      ▼  rectify align  (Step 1)
      │  ┌──────────────────────────────────────────────┐
      │  │ Phase 1: mapPacBio alone (all threads)        │
      │  │ Phase 2: minimap2 + gapmm2 in parallel        │
      │  │                                               │
      │  │ Per aligner BAM → extract AlignmentInfo:      │
      │  │   • effective 5' clip  (MD-free, genome-only) │
      │  │   • 5' clip rescue     (single-pass, scoring) │
      │  │   • A-tract 3' depth   (genome-only estimate) │
      │  │   • 3' terminal errors (non-poly-A clip)      │
      │  │   • junction-proximity errors                 │
      │  │                                               │
      │  │ Composite score → select best aligner per read│
      │  │ Tiebreakers: canonical GT-AG → annotated →   │
      │  │              majority 3' vote → wider span    │
      │  │                                               │
      │  │ Output: consensus.bam  (winning aligner BAM   │
      │  │         record + MD tags via samtools calmd)  │
      │  └──────────────────────────────────────────────┘
      │
      ▼  rectify correct  (Step 2)
      │  ┌──────────────────────────────────────────────┐
      │  │  Per-read correction modules (in order):      │
      │  │                                               │
      │  │  Module 2H — N-op junction refinement         │
      │  │  Module 2G — 3' soft-clip HP rescue           │
      │  │  Module 2E — A-tract ambiguity (walk-back)    │
      │  │  Module 2B — Poly(A) trimming                 │
      │  │  Module 2C — Indel artifact correction        │
      │  │  Module 2F — 5' junction rescue (Cat3/Cat4)   │
      │  │  Module FJF — False junction filter           │
      │  │  Module AG  — AG-mispriming (dT-cDNA only)    │
      │  │  Module NET — NET-seq refinement (optional)   │
      │  │                                               │
      │  │  Output: corrected_3ends.tsv                  │
      │  │          rectified_corrected_3end.bam         │
      │  │          rectified_pA_tail_trimmed.bam        │
      │  └──────────────────────────────────────────────┘
      │
      ▼  rectify restore-softclip  (Step 3 — DRS only)
         ┌──────────────────────────────────────────────┐
         │ Re-attaches trimmed poly(A) tail as soft-clip │
         │ using per-read tail lengths from metadata.parquet│
         │ Output: rectified_pA_tail_soft_clipped.bam   │
         └──────────────────────────────────────────────┘
```

---

## Step 1: Pre-Consensus Scoring (`rectify align` / `rectify consensus`)

### Motivation

The key insight driving this design is: **it is cheaper to select the right
alignment before correction than to correct a bad alignment after the fact**.
Choosing the aligner that already has the cleanest junction boundaries and the
3' end closest to the true CPA site means that `rectify correct` — especially
the expensive junction refinement step (Module 2H) — encounters fewer reads
that need intervention, and the interventions it does make are less dramatic.

### What happens to each read before selection

`extract_alignment_info()` in `consensus.py` is called for every read from
every aligner's BAM. It computes five signals from the raw alignment using
genome sequence alone (no MD tags required at this stage):

#### Signal 1 — Effective 5' clip (`_get_effective_5prime_clip`)

Scans up to 20 bp from the 5' end of the alignment for terminal imperfections:

- **Explicit soft-clip** (minimap2, gapmm2): the `S` op length in the CIGAR.
- **mapPacBio forced-mismatch region**: mapPacBio deliberately forces
  mismatches rather than soft-clipping at splice junction boundaries.
  `_get_effective_5prime_clip` scans the terminal aligned bases using a
  greedy sliding-window scan and extends the effective clip length to include
  any contiguous terminal mismatch/indel region.

Both representations are treated identically in scoring. This levels the
playing field between aligners that soft-clip and those that force mismatches.

The extracted terminal sequence (`five_prime_softclip_seq` for explicit clips,
`effective_five_prime_clip_seq` for MPB mismatch regions) is used as the
rescue sequence in Signal 2.

#### Signal 2 — 5' soft-clip rescue (`_rescue_5prime_softclip`)

For each read with a non-zero effective 5' clip, a lightweight sequence-based
rescue is attempted:

1. Collect all candidate junctions for this read (annotated + all aligners'
   observed junctions from this read group).
2. For each candidate junction within `junction_proximity_bp` of the
   alignment's 5' end: extract the upstream exon-1 sequence window and compute
   the edit distance against the rescue sequence.
3. If `edit_distance / clip_length ≤ max_edit_frac` (default 20%): **rescue
   fires** — the 5' clip penalty is waived entirely for this aligner.

This is a simplified, single-pass version of the full rescue in
`splice_aware_5prime.py`. It does not perform the shift-aware, HP-weighted,
ambiguity-resolving rescue that runs later in `rectify correct`. Its sole
purpose here is to avoid penalising aligners that correctly identified an
intron but could not align the upstream exon fragment — a penalty that would
cause the selection to unfairly prefer aligners that avoided the junction
entirely.

**Important:** the rescue result here does **not** change the read's position
or CIGAR. It only affects the alignment's score for the purposes of aligner
selection.

#### Signal 3 — A-tract 3' depth (`calculate_atract_ambiguity`)

Using genome sequence, estimates how far the aligner's 3' end is from the
true CPA site in A-tract regions. The `downstream_a_count` (number of A's
downstream of the raw 3' end, within 10 bp) quantifies how deep into the
A-tract the aligner landed. Penalty: −1 per downstream A, capped at 10.

This pre-estimate uses genome sequence only. The full MD-tag-dependent
walk-back correction runs in `rectify correct` (Module 2E/2C).

A pre-corrected `corrected_3prime` position is also computed (genome-only,
best-guess) and used as a tiebreaker in `select_best_alignment` to prefer
aligners that agree with the majority 3' position.

#### Signal 4 — 3' non-poly(A) terminal errors (`_get_effective_3prime_clip`)

Scans the 3' terminal region for non-A/T errors — mismatches or indels that
are NOT part of a poly(A) tail. These indicate the aligner stopped before the
true transcript end. Penalty: −2 per base, capped at 10.

#### Signal 5 — Junction-proximity errors (`_count_junction_proximity_errors`)

For each N-op (splice junction) in the read, counts mismatches and indels
within 5 bp of each junction boundary. Penalty: −1 per error, capped at 10.

This signal favours aligners that produce clean junction handling (mapPacBio
is typically best here). An aligner that places junction boundaries with
forced mismatches immediately flanking the N-op will score worse than one
that produces clean M-ops up to the splice site.

### Composite score and selection

```
score = − 2 × effective_5prime_clip     (0 if rescued)
        − 1 × atract_depth              (capped at 10)
        − 2 × effective_3prime_clip     (capped at 10)
        − 1 × junction_proximity_errors (capped at 10)
```

**Tiebreakers** (applied in order when scores are equal):

1. Most canonical GT-AG junctions
2. Most annotated junction matches
3. Majority vote on `corrected_3prime` (genome-only estimate)
4. Wider reference span (more of the transcript covered)

The winning aligner's raw BAM record is written to `consensus.bam`. MD tags
are then added via `samtools calmd` to enable the MD-dependent modules in
`rectify correct`.

---

## Step 2: Correction Modules (`rectify correct`)

`rectify correct` processes the consensus BAM read by read, applying modules
in a fixed order. The order matters: some modules update the read's coordinate
and CIGAR, and later modules operate on the updated state.

### Module 2H — N-Op Junction Refinement (`junction_refiner.py`)

**When:** First, before any 3' or 5' end corrections.

**What:** For every N-op (splice junction) in every read, tests all candidate
junctions within `search_radius` (default 5000 bp) and replaces imprecise
N-op boundaries with the best sequence-supported junction.

**Scoring:** HP-aware split-alignment. The rescue sequence (bases downstream
of the current N-op split point in transcript orientation) is split at every
candidate `k ∈ [0, L)` and both halves scored against the flanking exon
genomic context:

```
t1(k) = hp_score(rescue[k:],       genome[intron_end  : intron_end + buf])
t2(k) = hp_score(rescue[:k][::-1], genome[intron_end - buf : intron_end][::-1])
total(k) = t1(k) + t2(k)
```

HP-aware edit distance: indels within a homopolymer run cost 0.5; non-HP
indels and substitutions cost 1.0. When the current N-op is non-canonical
(tier ≥ 4), canonical-tier alternatives receive a 0.5 discount
(`_CANONICAL_HP_PRIOR`) — equal to the expected noise floor for one Nanopore
HP deletion.

**Priority:** sequence score first → current-junction stability (equal-scoring
candidates never displace an already-correct junction) → canonical GT-AG →
annotated → smallest boundary shift. **Annotation never overrides a
better-scoring junction.**

**Fast path:** Reads already at an annotated canonical-tier-0 junction skip
scoring entirely (255× speedup).

**Requires:** `--aligner-bams` to supply per-aligner BAMs for novel junction
discovery.

### Module 2G — 3' Soft-Clip Rescue at Homopolymer Boundaries (`rescue_softclip_at_homopolymer`)

Homopolymer under-calling causes aligners to terminate early and soft-clip
non-T/non-A bases that belong to the transcript body. Extends the 3' end
through remaining homopolymer reference bases to match downstream sequence.

### Module 2E — A-Tract Ambiguity Detection (`calculate_atract_ambiguity`)

Full MD-tag-dependent A-tract correction. Refines the genome-only pre-estimate
from consensus scoring using per-base mismatch information from the MD tag.

### Module 2B — Poly(A) Trimming

Trims poly(A) tail bases from the alignment. For DRS data these were stripped
before alignment; this module handles residual cases.

### Module 2C — Indel Artifact Correction (`find_polya_boundary`)

Walk-back algorithm: scans backward from the soft-clip boundary position,
skipping A's, deletions, T sequencing errors, and N-ops, until the first
unambiguous non-A/T genome–read agreement. This recovers the true CPA site
even when the aligner has spread poly(A) signal across genomic A-runs or
introduced spurious junctions to reach downstream A-tracts.

Key guards:
- **Large-deletion pre-scan:** detects over-calling artifacts where a large
  deletion bridges a poly-A over-extension back to the exon body.
- **N-op boundary guard (minus strand):** stops the forward scan at the first
  N-op boundary, preventing the scan from crossing into an earlier exon.
- **Trailing-base false-stop guard:** if the last base of a poly-A tail
  coincidentally matches a genomic base, the scan continues rather than
  stopping prematurely.

### Module 2F — 5' Junction Rescue (`rescue_3ss_truncation` in `splice_aware_5prime.py`)

Full Cat3 / Cat4 rescue, using the HP-weighted, shift-aware algorithm with
ambiguity-window detection. Four cases handled in priority order:

| Case | Trigger | Action |
|------|---------|--------|
| **1** | 5' soft-clip whose bases match upstream exon within `max_edit_frac` | Snap 5' end to intron boundary; emit M/I/D exon CIGAR via local NW aligner |
| **2** | mapPacBio forced-mismatch terminal region matching upstream exon | Same as Case 1 |
| **3** | 5' end within `junction_proximity_bp + five_clip` of a known 3'SS but no sequence evidence | Record junction hit; do not move 5' end |
| **4** | 5' end strictly inside an annotated intron, no existing N-op covers it | Snap to exon-1-side boundary (intronic snap) — only when exon sequence is a strictly better match than intron sequence (guards against misclassifying unspliced pre-mRNA reads) |

The proximity threshold for Cases 1/2/3 is extended by the soft-clip length
(`dist > junction_proximity_bp + five_clip`), so reads whose soft-clipped
bases reach across the junction boundary are not filtered out by the distance
check.

### Module FJF — False Junction Filter (`false_junction_filter.py`)

Detects and removes poly(A)-artifact junctions: N-ops where the bases flanking
the junction boundary are A-rich (characteristic of poly(A) tail misalignment).

### Module AG — AG-Mispriming Detection (dT-primed cDNA only)

Detects reverse-transcriptase slippage artifacts at internal AG-rich sequences.
Disabled for DRS data (enabled with `--dT-primed-cDNA`).

### Module NET — NET-seq Refinement (optional)

NNLS deconvolution using a reference NET-seq signal to refine 3' end positions
at multi-peak loci. Enabled with `--netseq`.

---

## Key Design Decisions

### Why score before correcting?

Consensus scoring is fast (genome-only, no MD tags) and run in parallel across
aligners. Running the full correction pipeline on all three aligner BAMs and
comparing outputs would be 3× slower and would not necessarily produce a better
winner selection — the lightweight pre-signals are sufficient to identify the
best-aligned read. The correction modules then operate on a cleaner starting
point, with fewer errors to fix.

### Why does Module 2H run first?

Junction refinement (2H) needs the N-op boundaries to be as accurate as
possible before the 5' end rescue (2F) attempts to snap reads to the exon-1
boundary. If 2F ran first with a poorly-placed N-op, it might snap to the
wrong exon boundary. Running 2H first ensures the N-op positions are as
accurate as the sequence evidence allows before 2F uses them as targets.

### Why does the 5' rescue in consensus scoring NOT change positions?

The rescue in `score_alignment` (`_rescue_5prime_softclip`) is a scoring
heuristic only — it waives the clip penalty for aligners that found a correct
junction but could not fully align the upstream exon end. Position changes
are deferred to `rectify correct` (Module 2F), which runs the full
shift-aware, HP-weighted rescue with proper CIGAR surgery.

### Unspliced pre-mRNA guard (Case 4)

Case 4 (intronic snap) fires when a read's 5' end is inside an annotated
intron. Without a guard, this would incorrectly snap the 5' end of unspliced
pre-mRNA reads to the exon-1 boundary. The guard compares the intronic query
bases against both the intron reference and the exon-1 reference using HP-aware
edit distance: snap only fires when the exon-1 match is strictly better (ties
favour the unspliced interpretation, keeping the read in the intron).

---

## Implementation Files

| File | Responsibility |
|------|---------------|
| `consensus.py` | `extract_alignment_info`, `score_alignment`, `select_best_alignment`, `_rescue_5prime_softclip`, `_get_effective_5prime_clip`, `_get_effective_3prime_clip`, `_count_junction_proximity_errors` |
| `multi_aligner.py` | Per-aligner subprocess management, two-phase scheduling |
| `bam_processor.py` | Per-read correction orchestration; calls all modules in order |
| `splice_aware_5prime.py` | Module 2F: full Cat3/Cat4 5' junction rescue |
| `junction_refiner.py` | Module 2H: N-op refinement, HP-aware scoring |
| `indel_corrector.py` | Module 2C/2E: walk-back, A-tract, poly(A) boundary |
| `false_junction_filter.py` | Module FJF: poly(A) artifact junction detection |
| `bam_writer.py` | CIGAR surgery for Cat3 extension, intronic tail clipping, 3' soft-clip rescue |
| `local_aligner.py` | Semi-global NW (Gotoh affine gap) for Cat3 exon CIGAR |
| `atract_detector.py` | A-tract ambiguity calculation (genome-only, used pre-consensus) |
