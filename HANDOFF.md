# Handoff — anchor gate + aligner panel + Round-2 cDNA scope (session 2026-06-13 → 14)

**Date:** 2026-06-14
**Branch:** `walkback-guard-refactor-todo`  ⚠ NOT `drs-validation-rebuild` (see Open items #1)
**HEAD:** `5ff7a6a`

---

## 1. What was done

- **HP-ED junction-anchor gate** — closes the mapPacBio consensus-gaming vector (N-ops
  scored free → noisy soft-clips relabeled as spurious introns win winner-selection).
  `_cigar_min_junction_anchor()` + a gate branch in `_add_chimera_flag` (the primary
  winner-sort key ahead of HP-ED). Default `_MIN_JUNCTION_ANCHOR=0` → gate OFF →
  byte-identical behavior (yeast safe by construction); human sets 10. 12 new tests. (`b448b4b`)
- **Aligner algorithmic-orthogonality survey** (4 Opus agents) + reject registry, written as
  the versioned primary source. (`ba00e5b`)
- **Round-2 cDNA "discovery → assignment" feature SCOPED** (2 designers + 1 adversary;
  Kevin's idea; MAGeSTIC bc0/bc1 analogy). Master spec + 2 designer drafts. (`878c3a2`, `5ff7a6a`)
- **Empirical aligner smoke-tests on Sherlock** (env `apanel`): GMAP viable; Graphmap2 +
  Magic-BLAST disqualified (see Verified).

## 2. What's verified

- `pytest` consensus/merge blast radius: **75 passed + 4 skipped**; **12 new** anchor-gate
  tests pass. Default-off change is non-regressing.
- **End-to-end A/B on the REAL 4-aligner consensus** (Sherlock job 29498330, COMPLETED): K=0 vs
  K=10 over identical per-aligner inputs — **novel (spurious) junctions −25%** (38,385→28,767),
  **annotated (real) +0.3%** (905,081→908,004); mapPacBio wins −5,266, flips go toward
  minimap2/deSALT. Anchor separation: real introns ~13 bp min-flank perfect-match anchor vs
  ~2 bp spurious (A549 chr5, 102,769 N-ops).
- **GMAP smoke** (job 29549741): fresh `gmap_build` + align exit 0, **7,198/10,000 spliced**,
  clean N-op CIGARs → VIABLE.
- **Graphmap2**: `-x rnaseq` CRASHES in postprocessing (cluster 16/52 → 0-byte SAM); default
  mode runs (no splice). 2020 staleness → NOT viable for RNA. **Magic-BLAST** (genome): emits
  pathological mega-gapped CIGARs (20–35 kb D) → unsuitable on genome (may be OK on contiguous
  cDNA — see Round-2).
- **NOT VERIFIED:** anchor gate under a SHARED prescan pool (the A/B used self-pools — Open #4);
  GMAP at full depth / the expanded 5-aligner consensus; the Round-2 feature (scoped, zero code).

## 3. Open items (with why)

1. **Branch.** Every commit this session is on `walkback-guard-refactor-todo`, not the canonical
   `drs-validation-rebuild`. Decide: cherry-pick `b448b4b ba00e5b 878c3a2 5ff7a6a` onto
   `drs-validation-rebuild`, or keep here. *Why deferred:* Kevin's call; working tree had heavy
   unrelated WIP, didn't switch branches.
2. **Anchor-gate productionization (READY, validated).** Add CLI `--min-junction-anchor-bp` +
   per-organism default (human=10 / yeast=0); wire ALL 3 `merge_corrected_tsvs` call sites incl.
   the `split_command.py` CHUNK TEMPLATE (the two `single_sample.py` sites have `args` in scope);
   re-admit mapPacBio to the human panel; reconcile docs. *Why deferred:* validated end-to-end, but
   Kevin pivoted to the aligner survey before flipping the default. Gate currently default-OFF.
3. **Sherlock overlay.** `corrected_consensus.py` on Sherlock (`/oak/.../software/rectify`, editable
   install) is overlaid AHEAD of M1 with the gate (md5 `a056f33`; original backed up
   `.bak_pregate_20260613`). Default-off → safe. Resync from M1 properly when convenient.
4. **Gate A/B used SELF-POOLS** (Kevin's catch — not production-faithful). Production-faithful
   expanded-panel measurement: `rectify prescan --aligner-bams <all panel BAMs>` → per-aligner
   `correct --junction-pool-cache/--variant-scan-cache` (shared pool) → gated `merge` → per-aligner
   sole-win. This SUPERSEDES the quick A/B (whose −25% relative effect still holds). Guard pool
   composition (deSALT/Graphmap2/mapPacBio artifacts; sane `--max-junction-size`).
5. **Expanded 5-aligner panel (GMAP integration), not started.** GMAP full-depth align A549 chr5
   (env `apanel`) → prescan(5) → shared-pool correct → gated merge → per-aligner sole-win quality.
   Target human panel = minimap2 · uLTRA · deSALT · mapPacBio(gated) · GMAP.
6. **Round-2 cDNA feature — scoped only.** MUST run **Phase 0 kill-gate FIRST** (does a non-trivial
   read subset align better to cDNA than Round-1, and beat **uLTRA specifically**?). **BLOCKER to fix
   in code before build:** the Cat3 `five_prime_rescued` exemption
   (`_effective_chimera_ok = _chimera_ok & (_five_rescued==0)`) lets a lifted cDNA row bypass the
   anchor gate → gate the win-guard on RAW `min_junction_anchor`. Non-splice bake-off panel:
   minimap2-no-splice · mapPacBio · GMAP · **LAST** (new lead) · **Magic-BLAST `-reftype
   transcriptome`** (gate on HSP-fragmentation). Full plan: `dev/specs/SPEC_round2_cdna_discovery_assignment_20260614.md`.
7. **⭐ ORIGINAL DEFERRED DELIVERABLE — human validation read-set vetting.** The session STARTED here
   and it is still unfinished. The committed 9-category human DRS validation set (public SG-NEx A549,
   chr1/5/11/17/19) is bundled at `~/igv_data/a549_validation/combined/` (72 reads, `VETTING_all.tsv`)
   awaiting Kevin's IGV sign-off on ALL 9 categories before building the committed artifact. Deferred
   when Kevin asked about aligner combining. See `dev/handoffs/HANDOFF_2026-06-13_human_validation_morechrom.md`
   + `dev/handoffs/STATUS_human_validation_readset.md`.

## 4. Resume

Pick the thread by intent:
- **Ship the anchor gate** → Open #2 (CLI flag + organism default + chunk template + re-admit
  mapPacBio); prove no-yeast-regression via the yeast validation suite; flip human default to 10.
- **Measure GMAP's orthogonal value** → Open #5 (env `apanel` on Sherlock;
  `W=/scratch/users/kevinroy/rectify_human_validation/sgnex_a549`; `$W/orthopanel/` has the index +
  smoke scripts; `$W/gate_ab/` has the A/B harness to clone).
- **Build Round-2** → Open #6; fix the Cat3 BLOCKER, then Phase 0 kill-gate on a minimal
  hand-built cDNA library.
- **Finish the original task** → Open #7: `bash ~/igv_data/a549_validation/combined/load_igv.sh`,
  step Kevin through `VETTING_all.tsv` category by category.
- **Before ANY further commit:** resolve Open #1 (branch).
- **M1 swap near-full** (12.3/13.3 GB used, ~0.1 GB RAM free) — keep all heavy work on Sherlock.

## 5. Files touched

Committed (branch `walkback-guard-refactor-todo`):
- `rectify/core/consensus/corrected_consensus.py` — anchor gate (`b448b4b`)
- `tests/test_junction_anchor_gate.py` — new, 12 tests (`b448b4b`)
- `docs/ALIGNER_RECOMMENDATIONS.md`; `docs/aligners/{gmap,graphmap2,magicblast,EVALUATED_AND_DISQUALIFIED}.md` (new); `docs/aligners/mapPacBio.md` (re-admission note) (`ba00e5b`)
- `dev/specs/SPEC_round2_cdna_discovery_assignment_20260614.md` (master) + `..._isoform_library_designer1_...` + `..._scoring_and_liftover_DESIGNER2_...` (`878c3a2`, `5ff7a6a`)
- `HANDOFF.md` archived → `dev/handoffs/HANDOFF_2026-05-25_93d99b8.md`

Uncommitted, pre-existing (NOT this session's gate/docs work — from the earlier human-validation
thread): `tests/test_calibrate_overhang_header.py`, `tests/test_data_bundling_human.py`,
`dev/specs/TODO_cdna_corrected_reads_tsv.md`, ` M dev/specs/SPEC_overcall_rescue_and_ed_metric_20260529.md`.

Off-repo:
- **Sherlock:** `corrected_consensus.py` overlaid (Open #3); conda env `apanel`
  (gmap/graphmap2/magicblast/blast/samtools); `$W/orthopanel/` (smoke + GMAP index),
  `$W/gate_ab/` (A/B outputs + `merge_ab.py`/`compare_gate_ab.py`).
- **Memory:** `project_hped_anchor_gate.md`, `project_aligner_orthogonality_panel.md`, the
  Round-2 spec pointer, + `MEMORY.md` index lines.
