# HANDOFF — RECTIFY cDNA Pipeline Refactor (Task 4)

**Date:** 2026-05-14
**Working dir:** `/Users/kevinroy/work/rectify`
**Authoritative spec:** `OPUS_IMPL_BRIEF.md` in this repo (Task 4 section)
**Project memory:** `/Users/kevinroy/.claude/projects/-Users-kevinroy-Library-CloudStorage-GoogleDrive-kevinrjroy-gmail-com-My-Drive-Work-Chanfreau-Lab/memory/project_cdna_refactor.md` (read this first — has the locked-in design decisions)

---

## Architecture (current target)

```
[user-pre-aligned input BAM]
  ↓
rectify correct-cdna
  - UMI extract → directional cluster → POA consensus per cluster
  - emits stage1_consensus.fastq.gz with per-cluster SAM tags in the FASTQ comment line
  - emits clusters.tsv / isoforms.tsv / t1t2_pairs.tsv  ← STILL HERE, NEEDS REMOVAL (see Task 4 work item 4)
  ↓
rectify align (done by another agent in a parallel session)
  - triple-aligner (minimap2 + mapPacBio + gapmm2) + chimeric-consensus integration
  - uses `minimap2 -y` to propagate FASTQ comments → BAM tags
  - output: one BAM record per cluster, tagged with XU/XC/XR/XO/XT/XY/XF/XB
  ↓
rectify cdna-analyze  ← we're building this
  - reads the post-align BAM
  - recomputes XA (walkback), pos5_corrected (walk-forward), XG/XS (gene/strand),
    XI (isoform), XL (T1↔T2 pair) on POST-ALIGN positions
  - writes the FINAL clusters.tsv / isoforms.tsv / t1t2_pairs.tsv
```

---

## Task status (from OPUS_IMPL_BRIEF.md)

| Task | Status | Notes |
|------|--------|-------|
| 1 — `read_subtype` field + XY tag | ✅ done | values `"umi_captured_fwd"`/`"umi_captured_rev"`/`"umi_not_captured"` (Kevin's preferred naming, not the brief's `"1a"`/`"1b"`/`"2"` — see memory) |
| 2 — `pretrim_consensus()` | ✅ done | `TestPretrimConsensus` passes; function kept as utility, no longer called inside correct-cdna |
| 3 — triple-aligner in correct-cdna | ⚠️ **DROPPED** | Architectural pivot: triple aligner moved to `rectify align`. Do NOT re-add re-alignment to correct-cdna. See memory `project_cdna_refactor.md`. |
| 4 — Move manifest analysis to downstream | 🚧 **IN FLIGHT** | This handoff is for finishing Task 4. |
| 5 — Stale Type-2 docstring | ✅ done |

---

## Task 4 progress (what's in flight)

### ✅ Done in this session

1. **Created [rectify/core/cdna_analyze_command.py](rectify/core/cdna_analyze_command.py)** — the new `cdna-analyze` module.
   - `_read_info_from_bam_record(rec, chrom_seq)` builds a synthetic `ReadInfo` per cluster, recomputing tail_len + pos5_corrected on post-align coords
   - `run(args)` orchestrates: load reference + GFF → stream BAM → classify XS/XG → assign_isoforms → reconcile_t1_t2_pairs → write 3 TSVs
   - Each cluster's `ReadInfo` is inflated by `XC` cluster_size so existing `assign_isoforms` / `reconcile_t1_t2_pairs` weight=`len(clusters[cid])` logic works unchanged
   - Reuses (imported from `cdna_correct_command`): `ReadInfo`, `walk_back_anchor_and_tail`, `walk_forward_tss`, `load_gff_genes`, `classify_sense_antisense`, `assign_isoforms`, `reconcile_t1_t2_pairs`

2. **Registered `cdna-analyze` subcommand in [rectify/cli.py](rectify/cli.py)** — both the parser block (just after the `correct-cdna` parser) and the dispatch `elif args.command == 'cdna-analyze':` (just before `restore-softclip`).
   - Required args: `bam` (positional), `--out`, `--gff`, `--reference`
   - Optional args: `--min-gene-frac`, `--min-read-frac`, `--isoform-tol-5`, `--isoform-tol-3`, `--t1t2-tol-5`, `--t1t2-tol-3` (all default 5 except gene-frac=0.3 / read-frac=0.8)

### 🚧 Still to do (in priority order)

1. **Verify the module imports and CLI help renders.** Run:
   ```bash
   cd /Users/kevinroy/work/rectify
   python -c "from rectify.core import cdna_analyze_command; print('ok')"
   python -m rectify cdna-analyze --help
   ```

2. **Strip pre-align tags from `write_stage1_fastq` in `cdna_correct_command.py`.**
   - The FASTQ comment currently still includes `XA`/`XS`/`XG`/`XI`/`XL` from the pre-align computation done inside `correct-cdna`. Those values will be wrong / misleading once `cdna-analyze` recomputes them post-align.
   - Keep only the alignment-independent tags: `XU`, `XO`, `XC`, `XR`, `XM`, `XF`, `XT`, `XY`, `XB`, and `XP` if/when paired_partner is used.
   - File: [rectify/core/cdna_correct_command.py](rectify/core/cdna_correct_command.py), look for `tag_parts = [...]` inside `write_stage1_fastq` (around line 1709 — line may have shifted; grep for `tag_parts = [`).

3. **Remove manifest-writing from `correct-cdna`'s `run()`.**
   - File: [rectify/core/cdna_correct_command.py](rectify/core/cdna_correct_command.py), look for the lines that write `clusters.tsv`, `isoforms.tsv`, `t1t2_pairs.tsv` (currently around lines 1842–1897). Also drop the calls to `assign_isoforms()` and `reconcile_t1_t2_pairs()` upstream (find with `grep -n "assign_isoforms\|reconcile_t1_t2_pairs" rectify/core/cdna_correct_command.py`).
   - Also remove the pre-align `cluster_xs` / `cluster_xg` computation (`classify_sense_antisense` call) if it's only used for the manifest now.
   - **Caveat:** UMI clustering still uses 3' anchor from input BAM as a bucketing key — that's fine, it doesn't need gene assignment.

4. **Update `correct-cdna`'s output-summary print block.** Remove references to `cluster_xs` counts and isoform/pair counts (they're not computed here anymore). Keep input-read / cluster-count / read-type stats.

5. **Update [rectify/cli.py](rectify/cli.py)** `correct-cdna` help text to drop the isoform/t1t2 description (it now only emits the FASTQ + a barebones cluster manifest — or no manifest at all; see decision below).

6. **Decide: does `correct-cdna` still emit a `clusters_pre.tsv`?**
   - Argument for: a downstream consumer that wants to skip `rectify align` (e.g., quick debugging) can pull cluster membership without parsing the FASTQ.
   - Argument against: redundant with the FASTQ-comment tags + read_ids in `XR`. Removing keeps the boundary clean.
   - **Recommendation:** drop it. The cluster_id ↔ input-read mapping is in the FASTQ's `XR` tag. If someone needs a TSV, they can `samtools view` the BAM and parse tags.

7. **Add a smoke test for `cdna-analyze`** in [tests/test_cdna_correct.py](tests/test_cdna_correct.py) (or a new `tests/test_cdna_analyze.py`). The test should:
   - Skip if the test BAM isn't on disk
   - Run `correct-cdna` → `rectify align` → `cdna-analyze` end-to-end on chrI, or use a pre-built consensus BAM if one exists
   - Assert that `clusters.tsv`, `isoforms.tsv`, `t1t2_pairs.tsv` exist and have row counts in expected bands (24k–30k clusters, 3k–4.5k isoforms, etc. — match the bands already in `test_correct_cdna_chri_smoke`)

8. **(Optional, follow-up)** Rewrite the cdna-analyze output to include a tagged BAM (one record per cluster with XU/XC/XR/XO/XT/XY/XF/XB/**XA**/**XS**/**XG**/**XI**/**XL** — i.e. the input BAM with the analyzer's new tags added). The brief doesn't mandate this; the user may want it for downstream consumers that prefer BAM-tag access over TSV joins.

---

## Critical context (do not violate)

- **`read_subtype` values**: `"umi_captured_fwd"` / `"umi_captured_rev"` / `"umi_not_captured"` — NOT `"1a"`/`"1b"`/`"2"`. See memory.
- **`read_type` (XT, integer 1/2)** is UNCHANGED — ~20 clustering sites depend on it. Don't refactor it away.
- **Triple aligner does NOT belong in `correct-cdna`.** If you see a future agent reaching for `realign_consensus*` or `mapPacBio`/`gapmm2` subprocess calls in `cdna_correct_command.py`, stop them — that work belongs in `rectify align`. The triple-aligner machinery was deleted from this file deliberately.
- **`correct-cdna`'s output is now `stage1_consensus.fastq.gz`** (gzipped FASTQ with per-cluster SAM tags in comment line), NOT a BAM. There is no `pysam.sort` / `pysam.index` in `run()`.
- **rDNA upstream masking is on the to-do list** but DEFERRED (see memory) — don't conflate with Task 4.

---

## How to verify what's there

```bash
cd /Users/kevinroy/work/rectify

# Imports clean?
python -c "from rectify.core import cdna_analyze_command, cdna_correct_command; print('ok')"

# CLI help works?
python -m rectify --help                # should list correct-cdna + cdna-analyze
python -m rectify cdna-analyze --help

# Unit tests (5s on M1)
python -m pytest tests/test_cdna_correct.py --deselect tests/test_cdna_correct.py::test_correct_cdna_chri_smoke --override-ini="addopts=" -v

# Smoke test (~3.5 min — only run after edits, not for quick iteration)
python -m pytest "tests/test_cdna_correct.py::test_correct_cdna_chri_smoke" --override-ini="addopts=" -v
```

Expected unit-test result: **23 passed** (TestUmi*, TestPosition*, TestBoundary*, TestRevComp, TestPretrimConsensus, TestSelectBestChimericSplice). Smoke test should pass in ~200–250s producing `stage1_consensus.fastq.gz` + the 3 TSV files (until Task-4 work item 3 is done, in which case the TSVs come from `cdna-analyze` instead).

---

## Files touched in this session

- `rectify/core/cdna_analyze_command.py` — **new module**, full content in place
- `rectify/core/cdna_correct_command.py` — already pivoted to FASTQ output earlier this session (see memory); manifest-writing logic still present and is item #3 to remove
- `rectify/cli.py` — `cdna-analyze` parser + dispatch added
- `rectify/__main__.py` — untouched
- `tests/test_cdna_correct.py` — smoke-test assertion updated to expect `stage1_consensus.fastq.gz`

---

## Test data paths (M1 dev only — skipped on CI)

- `/Users/kevinroy/work/ont_cdna/test_data/wt_rep1.chrI.bam` — chrI smoke BAM
- `/Users/kevinroy/work/ont_cdna/test_data/saccharomyces_cerevisiae_R64-5-1_20240529.gff`
- `/Users/kevinroy/work/ont_cdna/test_data/S288C_reference_sequence_R64-5-1_20240529.fsa`

---

## Open question for Kevin (don't proceed without confirming)

The other agent wired in `rectify align` to consume `stage1_consensus.fastq.gz`. Before running `cdna-analyze` end-to-end, **confirm that the post-align BAM has the per-cluster SAM tags** (specifically `XU`, `XC`, `XR`, `XO`, `XT`, `XY`, `XF`, `XB`):

```bash
samtools view <post-align-bam> | head -5 | tr '\t' '\n' | grep -E '^X[UCROTYFB]:'
```

If any tag is missing, that's a `rectify align` issue — fix that before `cdna-analyze` will work. (Most likely culprit: missing `-y` on the `minimap2` invocation, or the mapPacBio/gapmm2 wrapper not propagating `RG`/comment-format strings.)
