# SPEC — Round 2: cDNA "discovery → assignment" alignment (MASTER)

**Status:** scoped, GO-WITH-CONDITIONS (test-first). **Date:** 2026-06-14.
**Authoritative master** consolidating two designer drafts + adversarial review.
Designer detail (consume for full algorithms):
- `SPEC_round2_isoform_library_designer1_20260614.md` — isoform library + Round-2 aligners
- `SPEC_round2_cdna_scoring_and_liftover_DESIGNER2_20260614.md` — scoring + lift-over + integration

This master records the architecture (Kevin's framing), the decisions that supersede
the drafts, the adversary's blocker + must-fixes, and the **phased build plan with an
empirical kill-gate first**.

---

## 1. Motivation & framing — DISCOVERY → ASSIGNMENT

Separate junction **discovery** from alignment-quality **maximization over known junctions**.
Analogy (Kevin): MAGeSTIC barcoding — `bc0` guide-donor *linkage/discovery* vs `bc1`
*quantitation/assignment*. Here:

- **Round 1 = DISCOVERY** (exists): genome consensus discovers junctions, motif-agnostic
  (no GT-AG gating), with the k≥10 perfect-match anchor gate vetting them.
- **Round 2 = ASSIGNMENT** (new): align a *subset* of reads to a library of full-length
  cDNA isoform sequences (introns pre-removed, 5'/3' genomic-flank padded) using
  **non-spliced contiguous aligners**; the read picks its most-likely isoform; the winning
  cDNA alignment is lifted back to genome coords with N-ops at the isoform's junctions.

Hypothesis: a contiguous aligner freed from intron discovery places a non-trivial subset
of reads BETTER than any genome aligner (recovers exonic bases the genome aligner
soft-clipped at a junction it couldn't seed across; places micro-exons seed-aligners skip).

## 2. Why this is NOT "gaming" (Kevin's correction — symmetric scoring)

An earlier draft called the cDNA's intron-free reference a "trap" (it would trivially win).
**That was wrong.** HP-ED scores N-ops as **free** in BOTH directions:
- genome alignment: `…= N =…` → junction costs 0;
- cDNA alignment after lift-over: `…= N =…` (identical) → 0; equivalently pre-lift `…= = =…`
  seamless match also 0.

So the cDNA gets **no bonus** for the missing intron — the genome never paid for it. A cDNA
can win ONLY by reducing real **exonic** error (recovered soft-clip / micro-exon). Comparison =
**HP-ED on the cDNA's lifted genome CIGAR vs the genome alignment, under the same anchor gate.**
Adversary verdict: SOUND **with conditions** — symmetry holds in the clean case, but two
conditions are required (BLOCKER-1 fix + the no-shrink guard for noisy intron-retention).

## 3. Library admission (Kevin, authoritative)

cDNA isoforms are built ONLY from **(a) literature GTF/GFF isoforms ∪ (b) de-novo junctions
that PASS the Round-1 k≥10 perfect-match anchor gate**. Consequences:
- Mutant/pathological NOVEL junctions enter (real + gate-passing); spurious ones don't
  (fail k≥10) → unbiased-but-clean, resolving the WT-bias-vs-pollution tension.
- Every junction the cDNA round can transfer was vetted once in Round 1 (construction-time gate).
- **Adversary caveat:** construction-time admission does NOT protect the *transfer* step —
  a transferred N must RE-clear the gate at selection (see BLOCKER-1).

Build the library as **read-evidence-bounded isoform chains** (one cDNA per distinct
read-supported intron chain over annotation ∪ gate-passed junctions), NOT naive splice-graph
path enumeration (2^N blow-up at high-exon genes like TTN/363 exons). Prefer an **in-house
collapse** (~150 lines) over FLAIR defaults — FLAIR `correct`/`--check_splice` silently
re-imposes GT-AG + ≥3-read gating (FLAIR #371) and would erase the non-canonical junctions we
must keep; **mandatory spike test**: a read with a known non-canonical gate-passed Round-1
junction must survive into the library. Junction-wobble snap W=3 to highest-support junction
(validate it collapses, not fragments, on a noisy sample).

## 4. Padding & block-map (designer-1 interface)

- Pad **1–2 kb** each side (NOT flat 10 kb — swallows neighbors in the compact human genome),
  adaptive up to ~10 kb but **hard-clipped at the nearest intergenic boundary** (neighbor exon
  −50 bp guard). 3' pad ≥ the existing poly-A walkback look-ahead; pad is REAL genomic sequence.
- Per-cDNA **block map** (parquet/JSONL sidecar): ordered `blocks` of
  `{t_start, g_start, g_end, length, strand, is_pad}` (1:1 transcript→genome) + explicit
  `junctions` array. The transcript-contiguous / genome-jump boundary between adjacent blocks
  IS the N-op site for lift-over. **Minus-strand:** cDNA is RNA-sense (revcomp of plus genome);
  consecutive blocks have decreasing g_start; lift-over MUST revcomp + reverse CIGAR.
- PROVENANCE.json sidecar (GENCODE version, Round-1 junction-set hash, collapse method, params).

## 5. Round-2 aligners (non-spliced, contiguous to transcript)

| Aligner | Role | Invocation | Note |
|---|---|---|---|
| **minimap2** | PRIMARY | `-ax map-ont -N 100 -p 0.8 -Y`, NO splice/`-G` | installed |
| **mapPacBio** | secondary | long-read mode, `intronlen=0 maxindel=100` | installed (compute-node) |
| **GMAP→cDNA** | tie-breaker | `-n 1`, NO `--nofails` | installed (works; smoke 2026-06-14) |
| ~~BLAT, Magic-BLAST~~ | DROP | — | BLAT >95%-id (unfit raw ONT); Magic-BLAST emits pathological mega-gapped CIGARs on ONT DRS (smoke 2026-06-14) |

`-Y` preserves soft-clips (clips into pad must be visible to lift-over). **`--for-only` is
PER-PROTOCOL**: correct for DRS; for PCR-cDNA it silently drops the antisense half unless
`ont_cdna` pre-oriented the reads — VERIFY (grep XS/XR), branch by protocol.

## 6. Read-subset gate (cost control)

Run Round 2 only on reads that are **(A) spliced (span ≥1 junction)** AND **(B) weak in
Round 1** (large soft-clip OR high HP-ED OR failed anchor gate OR low-MAPQ / multi-aligner
disagreement — OR'd). Single-exon reads gain nothing (Round 2 ≡ Round 1). Round 2 may LOSE;
scoring keeps Round 1 when Round 2 isn't strictly better.

## 7. Scoring & winner integration

cDNA-derived genome alignment enters the EXISTING `merge_corrected_tsvs` as another candidate
(`_aligner = "cdna_round2"`) — reuse = the whole safety story (same gate, same sort). It wins
the read iff ALL hold (Designer 2 win-guard):
1. passes the perfect-match anchor gate on its transferred N-ops,
2. strictly lower HP-ED than the best genome candidate by ≥ ε,
3. **no-shrink**: aligned_bases not lower than the genome candidate's (over the read's mapped extent).

### ⚠ BLOCKER-1 (must fix before build) — the gate is bypassable via the Cat3 exemption
The selection engine computes `_effective_chimera_ok = _chimera_ok & (_five_rescued == 0)`
(`corrected_consensus.py` ~L1375). A `cdna_round2` row carrying `five_prime_rescued=1` gets
`_effective_chimera_ok = 0` **vacuously — anchor gate switched off** — and Round-2's proudest
wins (5' clip recovery / TSS-into-pad) are exactly the rows the genome path tags
`five_prime_rescued=1`. **Fix (adopt one, prove with a test):**
- (b, preferred) win-guard condition (1) gates on the RAW `min_junction_anchor ≥
  min_junction_anchor_bp`, NOT on `_effective_chimera_ok`; OR
- (a) force `five_prime_rescued = 0` on all `cdna_round2` rows.
**Test:** a lifted candidate with a sub-K-anchor N-op AND `five_prime_rescued=1` must be REJECTED.

### Safety precondition
Round 2 active ⇒ anchor gate ON (`min_junction_anchor_bp > 0`). On yeast (gate off by default)
a bad transferred junction would smuggle through — Round 2 must refuse to run, or force the gate
on, when `min_junction_anchor_bp == 0`.

## 8. Lift-over (transcript→genome + N-op insertion)

Designer 2 §2 algorithm + 12-row edge-case table. Highest-bug-density step = **minus-strand**:
- revcomp + CIGAR-reversal for `-` genes; assert the **orientation invariant** (first/last
  matched read base == genome at reference_start / reference_end−1) — length/span invariants are
  revcomp-blind and won't catch a missed revcomp.
- **Add test (adversary):** insertion/deletion abutting a junction ON a minus-strand gene — the
  "I attaches to the genome-5' (donor) exon" convention must survive CIGAR reversal (E6×E11 interaction).
- reads ending in 5'/3' pad = TSS/CPA variation → matches, not clips (then poly-A walkback runs).

## 9. Validation / falsification

- annotated-junction fraction = VALIDATION LABEL ONLY (never a selection input); novel
  gate-passed junctions must survive.
- **win-vs-uLTRA-specifically (adversary MUST-FIX):** uLTRA is already annotation-guided in
  Round 1. Report the fraction of `cdna_round2` wins that beat **uLTRA's own Round-1 record for
  that read**, not just the consensus winner — else the feature may be redundant.
- per-read win attribution: every cDNA win must have an attributable cause (recovered-clip /
  micro-exon / fragmented-seed). Trivial-win leak >5% with no attributable cause ⇒ metric is
  reverse-gaming ⇒ reject.
- noisy intron-retention test (the no-shrink guard, not HP-ED, is the SOLE IR protector in the
  noisy case — test with an SNV-bearing/noisy retained intron, not a clean one).

## 10. PHASED BUILD PLAN (test-first)

- **Phase 0 — EMPIRICAL KILL-GATE (do this FIRST, before any library machinery).** On a handful
  of loci with known micro-exons / multi-junction reads, hand-build a minimal cDNA library, align
  the Round-1-weak subset (minimap2-no-splice + mapPacBio), lift, score. GO criteria: a
  non-trivial read subset beats the Round-1 consensus **and uLTRA specifically**, >95% wins
  attributable, measurable net HP-ED reduction. **If it fails → NO-GO, stop.** (Riskiest
  assumption, cheapest test.) Substrate: lab human/yeast DRS with the anchor gate ON.
- **Phase 1 — Library construction**: in-house collapse (annotation ∪ gate-passed chains) +
  spike test + padding + block-map + PROVENANCE.
- **Phase 2 — Scoring + lift-over + integration**: implement the win-guard WITH the BLOCKER-1 fix;
  lift-over with the minus-strand tests; `cdna_round2` candidate into `merge`; gate-on precondition.
- **Phase 3 — Validation**: §9 metrics + falsification on full A549 chr5 (+ a mutant sample).

## 11. Risk register (adversary, prioritized)

| Sev | Item | Fix/test |
|---|---|---|
| BLOCKER | Cat3 (`five_prime_rescued`) exemption bypasses anchor gate (§7) | gate on raw anchor; test sub-K + rescued ⇒ reject |
| MUST-FIX | Empirical subset may not exist | Phase 0 kill-gate FIRST |
| MUST-FIX | IR protection rests on no-shrink guard, not HP-ED (noisy case) | implement + test with noisy IR |
| MUST-FIX | Redundancy with uLTRA | win-vs-uLTRA-specifically metric |
| WATCH | minus-strand lift-over (I/D at boundary × reversal) | dedicated unit tests + orientation invariant |
| WATCH | `--for-only` drops antisense PCR-cDNA half | verify ont_cdna orientation, branch by protocol |
| WATCH | library scale on noisy mutants | validate W=3 collapse-not-fragment vs depth |
| WATCH (benign) | uncombined co-junction coverage gap | document recall ceiling; Round-1 fallback guaranteed |

**Verdict:** GO-WITH-CONDITIONS. Single riskiest assumption to test first: that the rescuable
subset exists (beats the genome consensus incl. uLTRA). Run Phase 0 before building anything else.
