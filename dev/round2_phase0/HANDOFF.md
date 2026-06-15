# HANDOFF — Round-2 cDNA Phase-0 empirical kill-gate — **COMPLETE: NO-GO**

**Session:** 2026-06-14 · **Branch:** `walkback-guard-refactor-todo` (harness committed `5a40cbc`;
`cluster/` scripts + verdict docs UNCOMMITTED — see §4) · **Repo:** `~/work/rectify`

**One-line state:** Phase 0 RAN end-to-end and returned **NO-GO**. The cDNA "discovery→assignment"
round provides no meaningful spliced-read rescue on the A549 RNA002 substrate, robust against three
controls. **Per spec §10/§11, NO-GO ⇒ do NOT build the production library/lift-over (Phases 1–3)**
unless the substrate is changed — and that change is **Kevin's decision** (§5 below), not an auto-launch.

**Read:** `PHASE0_VERDICT.md` (full writeup) · `PROVENANCE.json` (params + table hashes).

---

## 1. Verdict

**NO-GO**, robust against three independent controls:
- **Corrected-genome baseline** (advisor: a raw baseline would be a false-GO generator).
- **Single-exon null control:** Round-2 ≡ Round-1 on single-exon reads, yet they "won" at ~10–14%
  from alignment jitter — the win-guard's false-win floor already exceeds the spec's 5% leak bar.
  Genuine structural spliced wins ≈ **7 / 4,579 (0.15%)**.
- **Selection-bias check (Kevin):** of 844,262 genome-UNMAPPED reads, 642 hit the cDNA library but
  **0 confident spliced rescues** — all confident placements single-exon; the 17 spliced hits are
  68–83%-identity force-fits of junk reads.

## 2. What was verified
- Harness 15/15 unit tests; **F.0 real-read lift round-trip: 833 reads, 660(+)/173(−), 0 mismatches.**
- Step D: 3-aligner Round-1 + per-aligner `rectify correct` on human DRS works (walkback fired).
- Step C: 42-cDNA library (7 MANE ∪ 35 gate-passed de-novo), per-base 1:1 identity validated.
- Bias check: 0/844k genome-unmapped reads get a confident spliced cDNA rescue.

## 3. Two lift/projection bugs found+fixed by F.0 (don't regress)
- **Projection** (`cluster/cdna_models.py`): must enforce transcript-coordinate contiguity, else a
  read using a NON-MANE junction is silently remapped onto MANE coords. Fixed; rejects such reads.
- **F.0 comparison** (`cluster/liftover_roundtrip.py`): compare semantic matched-(query,ref) pairs,
  NOT raw CIGAR — an insertion exactly at a junction is `M N I M` (minimap2) vs `M I N M` (lift),
  same alignment, different representation. `test_projection_local.py` 7/7 guards this.

## 4. Open items
- **`cluster/` scripts + `PHASE0_VERDICT.md` + `PROVENANCE.json` are UNCOMMITTED.** Stage with
  `git add dev/round2_phase0/` (NOT `-A`). Not committed because Kevin hasn't said to (CLAUDE.md).
- **mapPacBio dropped** from Round-2: BBMapPacBio caps reads at 6019 bp (throws on full-length DRS)
  AND games HP-ED. minimap2-map-ont is the sole Round-2 aligner; verdict is minimap2-cDNA only.
- **CPA-regression gate is confounded** (raw-cDNA-3' vs walked-genome-3'); reported as a diagnostic
  only (the leak/null-floor already decide the verdict). Fix = raw-vs-raw if ever revisited.

## 5. The decision that is Kevin's (do NOT auto-launch)
Honest caveat: A549 is **SRRM4-low** (micro-exon-poor) + **RNA002**. The **multi-junction** case was
well-powered (TRIO 56 junc, XPO5 31) and is **tested-and-weak** → that part of the NO-GO is solid.
The **micro-exon** sub-case is **under-tested**. Whether to retest on a **neuronal / SRRM4-high
RNA004** substrate before fully killing the feature is **Kevin's call**.

## 6. Files
`dev/round2_phase0/`: harness (`liftover.py`, `score_phase0.py`, `test_phase0.py`, `select_loci.py`,
`selected_loci.tsv`) + `cluster/` (`stepD_*.sbatch`, `build_library.py`, `cdna_models.py`,
`liftover_roundtrip.py`, `phase0_score.py`, `stepCF_*.sbatch`, `stepBias_*.sbatch`,
`rescue_unmapped.py`, `test_projection_local.py`, `debug_one_read.py`) + `PHASE0_VERDICT.md` +
`PROVENANCE.json` + this file.
Cluster workdir: Sherlock `/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0/round2/`.
