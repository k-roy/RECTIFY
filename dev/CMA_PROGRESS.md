# CMA (compressed-multialign) — footprint reclaim + format-core build

**Session:** 2026-07-23 · branch `drs-validation-rebuild` · rectify repo (M1/"the Mac" local + Hoffman2 measurement)
**Design of record:** `~/work/UCLA/Chanfreau_Lab/planning/254_rectify_multialign_compression_plan.md`
(M0 measurement gate = GO, see planning/256). This doc is the build + measurement record for the two
tracks Kevin approved ("both in parallel", target = Hoffman2 Chanfreau panels).

Milestone names (to avoid the M0/M1/M2 ↔ "M1 laptop" collision):
- **measurement gate** (doc M0) — DONE, GO (planning/256).
- **format-core build** (doc M1) — **DONE this session** (Track 2 below).
- **backlog reclaim** (doc M2) — measured + verified-safe, **proposal below awaiting Kevin's go** (Track 1).

---

## Track 1 — footprint reclaim on Hoffman2 (measured, read-only)

### What's on disk (du sweep — H2 `~/m0_scratch/footprint_sweep_20260723_1123.txt`)

| Run dir (`/u/scratch/k/kevinroy/`) | Size | Notes |
|---|--:|---|
| `drs_decay_mutants_rectify_20260702` | **187 GB** | 11 samples; **2-aligner panel (minimap2 + uLTRA)** |
| `408d_plants` | **127 GB** | not yet broken down |
| `398l_quantseq_epap` | 29 GB | short-read (QuantSeq) |
| `pare_invitro_381` | 11 GB | |
| `rectify_*` array dirs (×407) | 31 GB | small per-dir chunk/array scratch |

Per-sample split of the decay run: **`aligner_chunks/` ≈ 104 GB total**, `merged_bams/` ≈ 31 GB, rest =
raw chunk fastqs (~33 GB, regenerable) + sidecar parquet + consensus + pkls.

### Layout of one sample (`Dcp2_AID_25C_repA`, verified)
```
<sample>/
  <sample>_chunk_NNN_of_032.fastq.gz     # 32 raw chunks (~3 GB) — REGENERABLE via `rectify split`
  <sample>.read_num_sidecar.parquet      # 174 MB
  aligner_chunks/{minimap2,uLTRA}/chunk_NNN/
      *.<aligner>.bam                    # raw per-chunk per-aligner  → merged into merged_bams/
      *.<aligner>.junction_refined_*.bam # Module-2H intermediates    → pool already in junction_pool.pkl
      *.rectified.bam                    # per-chunk corrected         → superseded by consensus/
  merged_bams/<sample>.<aligner>.bam     # merged per-aligner (the durable pre-correct store)
  consensus/corrected_consensus.bam + corrected_reads.tsv   # FINAL corrected output
  junction_pool.pkl, rescue_scan.pkl, logs/
```

### Completeness verification (the deletion-safety gate — H2 `reclaim_verify2_20260723_1136.txt`)
Read-count equality, **raw-per-chunk-only vs merged**, primary reads (`-F 0x900`):

| Sample | aligner | raw Σ | merged | COMPLETE |
|---|---|--:|--:|:--:|
| Dcp2_AID_25C_repA | minimap2 | 3,158,533 | 3,158,533 | ✅ |
| Dcp2_AID_25C_repA | uLTRA | 514,857 | 514,857 | ✅ |
| Xrn1_AID_25C_repB | minimap2 | 2,030,361 | 2,030,361 | ✅ |
| Xrn1_AID_25C_repB | uLTRA | 1,037,369 | 1,037,369 | ✅ |

`consensus/corrected_consensus.bam` present for both. (An earlier substring-glob check reported a spurious
2× because `*minimap2*.bam` also matched the per-chunk `.junction_refined_*.bam`; the precise `*.minimap2.bam`
matcher above resolves it.)

### Reclaim proposal (NOT executed — needs Kevin's go + explicit file list, never a glob)
**Lever B (no new format, verified lossless-redundant):** the entire `aligner_chunks/` tree per sample is
reconstructible from durable products that all exist:
- raw per-chunk `*.<aligner>.bam` → **`merged_bams/` (verified complete)**
- `*.junction_refined_*.bam` → **`junction_pool.pkl` (already built)**
- per-chunk `*.rectified.bam` → **`consensus/corrected_consensus.bam` (final, present)**

⇒ **~104 GB reclaimable from the decay run alone**, no information loss. Likely similar on `408d_plants`
(127 GB) — not yet inventoried.

**Before any delete (pre-delete checklist):**
1. Run the same read-count completeness check on the remaining 9 decay samples (and inventory `408d_plants`).
2. Confirm no resume path consumes per-chunk BAMs (doc §3.3: `correct` resumes from the consensus BAM +
   `--junction-pool-cache`, both present).
3. **Read-count equality proves no reads were *dropped* — it does NOT prove the merged records are
   field-identical to the chunk records.** For strict safety, spot-check a sample's merged vs concatenated-chunk
   records on the normalized field view before deleting (cheap on one sample).
4. **"`junction_pool.pkl` makes the `.junction_refined_*.bam` redundant" is asserted, not verified.** Confirm
   the pool pkl is the sole downstream consumer of the refined per-chunk BAMs before deleting them.
5. Build an **explicit file list**, show Kevin, then remove named paths only (never a glob, never `rm -rf`).
   Optional extra: raw chunk fastqs (~33 GB) are regenerable but keep unless space-critical.

---

## Track 2 — CMA format core (BUILT + proven lossless)

New package `rectify/core/multialign/` + `tests/test_cma_drs_expand.py` (7/7 pass):
- `cma_schema.py` — tags `Za/Zp/Zk/Zq/Zm`, `collapse_key = (reference_name, is_reverse, reference_start,
  cigarstring)`, payload-donor choice, `hclip_bounds`, `revcomp`, **`decode_eq_seq`**.
- `cma_writer.py` — `build_cma(stream, header, out, panel, genome)` + `load_aligner_records` (prototype grouper).
- `cma_reader.py` — `expand(cma, genome)` — rebuilds per-aligner dict **from `Za`** (never re-feeds the
  grouper, doc §3.1 hard guard); reverse-frame SEQ/QUAL reconstruction (doc §2.4).
- `validate.py` — `validate_cma()` structural invariants.

**Test suite: `tests/test_cma_drs_expand.py` — 8/8 pass.** Two gates (planning/254 §3.3):

- **PRIMARY gate (the one that licenses deletion) — PASSES.** The pipeline's own consumer,
  `extract_alignment_info` + `select_best_alignment`, picks the **identical winner + corrected 3′ end on
  `expand(CMA)` vs the original per-aligner records, for all 36 reads** — compared on explicit-SEQ (the
  production case; see finding 1). This proves the pipeline *consumes* the expanded records identically, not
  merely that fields reproduce.
- **SECONDARY gate — PASSES.** `expand()` reproduces the load-bearing field view (SEQ nucleotides / POS /
  CIGAR / MAPQ / NM / MD / strand) on the fixture (36 reads, **18 reverse-strand frames**) + synthesized
  cases absent from the fixture: reverse-strand **+ hard-clip** worked-check (doc §2.4:
  `revcomp(F)[4:16] == "GGGGTTTTTTTT"`), gapmm2 SEQ=None restore, per-aligner MAPQ, RN-key join, no-donor
  fallback. Byte equality is never asserted (BGZF/@PG/aux-order noise).

**Footprint:** SEQ/QUAL payload stored **36× not 180× (5.0× fewer copies)**, losslessly.

**Scope of "proven lossless":** proven on the 36-read fixture (format core + both gates). NOT yet proven at
scale on a real run — that is the backlog-reclaim / DAG-wiring milestones, each with the deletion gate on a
real sample before any per-aligner BAM is deleted.

### Three findings the fixture surfaced (real, worth acting on)
1. **`=`-encoded SEQ is placement-relative — and the current pipeline mis-scores it.** The fixture SEQ uses
   SAM `=` (match-to-ref, from `calmd -e`). `=` cannot transfer to a variant at another locus, so the writer
   **decodes `=`→explicit nucleotides against the genome** before storing the payload (no-op on production
   explicit-SEQ BAMs). `get_aligned_pairs(with_seq=True)` is useless here (returns `=` for matches) — decode
   needs the actual reference. **Consequence measured:** feeding raw `=`-SEQ to `select_best_alignment`
   changes the winner for **14/36** fixture reads vs the decoded-explicit version; the CMA path is
   reference-correct. If any real BAMs are `=`-encoded, the current selection mis-scores them — worth checking
   independently of CMA.
2. **Per-aligner QUAL is inconsistent; CMA canonicalizes it to the donor.** In the fixture **deSALT stores
   reverse-strand QUAL in the opposite orientation** to the other aligners and **uLTRA drops QUAL entirely
   (`QUAL=*`)**. Whether deSALT's is *wrong* or a fixture-generation artifact is **not verified** — what is
   certain is it is inconsistent. QUAL is read-intrinsic and consumed by **no** algorithm (only MAPQ is, doc
   §2.6), so the CMA stores it once from the payload donor and reconstructs it for all placements;
   **per-aligner QUAL divergence is not preserved.** Real consequence to sign off on: for a read where a
   deSALT/uLTRA record wins consensus, the **final consensus-output QUAL will change** after a CMA round-trip.
   That is a genuine output diff (not test-only) for Kevin + the deletion gate, not just a fidelity note.

### Note on this run's panel
`drs_decay_mutants_rectify_20260702` is a **2-aligner** run (minimap2 + uLTRA), so its CMA SEQ/QUAL dedup is
~2× on `merged_bams` (~31 GB → ~17–20 GB) — smaller than Lever B here. The full ~3–5× CMA win lands on
5-aligner runs (mex67aa-style, 85 GB merged) and is what enables deleting the merged_bams double-copy too.

---

## Next steps
- **Track 1:** on Kevin's go — completeness-check the other 9 decay samples, inventory `408d_plants`, build the
  explicit delete list. (Reclaim ≈104 GB decay + TBD plants.)
- **Track 2:** next doc-milestones are **backlog reclaim ingest** (`rectify cma build --from-aligner-bams`, doc
  M2) and **DRS DAG wiring** (doc M4). The format core de-risks both.
- **Unrelated but open (inbox):** 2 CRITICAL Sumner-filed bugs in `AGENT_FIXES.md` — (a) `align_command.py`
  single-aligner output path is sample-keyed (per-aligner fan-out clobbers itself; `.multialigned.bam` rename
  didn't fix it — needs the aligner component); (b) `corrected_consensus.py` on the shared Oak env matches no
  commit. Both untouched by this session's work.
