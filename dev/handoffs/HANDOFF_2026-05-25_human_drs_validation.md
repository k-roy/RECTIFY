# Handoff — Human DRS Validation Set (SG-NEx A549) + RECTIFY human-readiness

**Date:** 2026-05-25
**Branch:** `drs-validation-rebuild`
**HEAD at handoff:** `961c844` — `perf(splice): dual-site interval-tree pool index for 3'SS-rescue fetch`
**Session commit:** `99ff250` — `fix(genome): organism-aware chrom standardization` (in history)
**Scope:** Build a committable 9-category human DRS validation set from **public** data,
modeled on the yeast validation set; surfaced and (partly) fixed several RECTIFY
human-readiness issues along the way.

> ⚠️ **The current A549 validation selection is NOT final** — it was built on
> alignments run with the yeast intron ceiling (`-G 5000`). Junction categories
> must be re-derived after re-aligning with `--max-intron 500000`. See §5.

---

## 1. Goal & hard constraints

- Produce **8 reads/category × 9 categories** (4 per strand) matching the yeast
  RECTIFY validation set at `rectify/rectify/data/validation`, but on **human** DRS.
- **Privacy (load-bearing):** the committed set must contain **no patient reads**
  (human fingerprint/variant exposure). Patient Sumner data (CNTL_21.8 etc.) is
  prototyping/calibration only. Committed set is sourced from **public, consented**
  human DRS.
- Cat8 (yeast NET-seq analog) substitute = **PolyASite 2.0** (bp-resolution CPA
  atlas with per-cluster TPM).

## 2. Decisions made (with rationale)

| Decision | Rationale |
|---|---|
| **Source = SG-NEx A549** (4 directRNA runs, 319k chr5 reads) | GM12878 rejected: not in standard open buckets, older R9.4. SG-NEx fully accessible from Sherlock login via `wget`/`samtools`-over-HTTPS (curl is broken on Sherlock — error 43). Bucket `s3://sg-nex-data` region `ap-southeast-1`, Ensembl chrom naming (`5`). |
| **4-aligner panel** (minimap2, mapPacBio, uLTRA, deSALT) | gapmm2 dropped on human — see §3. |
| **Subsample to 25%** (~80k reads, max depth 16.8k→4.2k) | A549 cancer-line pileups choke per-region `correct`; name-hash seed (`samtools -s 7.25`) keeps the *same* reads across all 4 aligner BAMs. |
| **`correct` run WITHOUT `--annotation`** | `--annotation` is intractably slow on human (§3); the classifier is self-sufficient for Cat3/7/8 + gene attribution. |
| **Cat8 reframed** to "A-tract correction lands on a TPM≥20 PolyASite CPA" | Genuine comparable multi-CPA usage is ~absent in A549 (0/0 strict). |
| **Cat3 uses NO aligner-agreement quorum** (CAT3_MIN=1) | Kevin: RECTIFY selects the single best alignment; single-aligner rescue is the *informative* case, not a thing to gate out. |
| **Cat5 paralog-multimap guard** (intron <500kb, pair span <400kb, rank by compactness) | Without it, picks were SMN-locus segmental-dup multimaps, not genuine 2-intron reads. |

## 3. Key RECTIFY human-readiness findings

1. **chrom-naming `chr5`→`chrV`** (FIXED, commit `99ff250`). `standardize_chrom_name`
   mapped `chr1`..`chr16` to yeast roman numerals → emptied `gene_id`, wrote `chrV`.
   Fixed via genome-derived `_KNOWN_CONTIGS` registry. See `AGENT_FIXES.md` 2026-05-24.

2. **`--max-intron` defaults to 5000 (yeast); `--organism homo_sapiens` does NOT
   raise it.** Feeds minimap2 `-G`, gapmm2 `-i`, mapPacBio `intronlen`. Both the
   A549 *and* patient Sumner aligns ran with `-G 5000` → reads spanning human
   introns >5 kb mishandled. **minimap2 impact bounded** (still agreed with
   uLTRA/deSALT on 312k reads — chr5 set dominated by <5kb-intron transcripts);
   **mapPacBio badly broken** (~1k introns, disjoint set — though its root issue is
   BBMap `intronlen`/`maxindel` semantics). **Fix: `--max-intron 500000` for human.**
   Docs: `docs/aligners/{minimap2,mapPacBio}.md`.

3. **`rectify correct --annotation` intractably slow on human** (~20 s/region ×
   ~1000 regions → 2 h timeout even at 80k reads; ~475× slower than without).
   Cause is NOT the 3'SS-rescue flag (enabled in both fast & slow runs) — it's the
   **annotated junctions in the junction pool**: long/numerous human introns fan out
   the per-read candidate-junction lookup. The patient full-chr5 run was fast only
   because the GTF-intron bug left 0 annotated junctions. **Being actively addressed
   by HEAD `961c844` (dual-site interval-tree pool index) + `dev/verify_dualsite_correct.sbatch`.**
   Until that's validated, run `correct` without `--annotation` on human.

4. **gapmm2 dropped from human panel.** Profiled: the single-threaded per-read
   terminal-refinement (edlib) loop dominates (~0.12–0.19 s/read → ~28–50 h full
   sample); `-i` has <6% wall-time effect. Marginal real-junction yield (~1.3–1.8%
   unique, ~85% artifacts, ~98% redundant with uLTRA/deSALT). Full detail:
   `docs/aligners/gapmm2.md`. Hit the 6 h `ALIGNER_TIMEOUT`.

5. **mapPacBio is a weak ONT-DRS splice aligner** and was mis-parameterized
   (`intronlen` should be 10, `maxindel` 200000; RECTIFY historically wired
   `intronlen={max_intron}`). Full detail: `docs/aligners/mapPacBio.md`.

## 4. Artifacts & where things live

**Sherlock — `/scratch/users/kevinroy/rectify_human_validation/sgnex_a549/`:**
- `chr5_bams/`, `a549_chr5.primary.bam` (319k primary reads), `a549_chr5.reads.fastq.gz`
- `trim/a549_chr5_trimmed.fastq.gz` (+ `..._polya_trim_metadata.tsv`) — **re-align input**
- `alignments/` (full-depth, `-G 5000`), `alignments_sub25/` (25% subset, `-G 5000`)
- `correct_sub25_noann/a549_chr5_corrected.tsv` — corrected, no annotation
- `a549_candidates.tsv` — classifier inventory (all categories)
- `a549_validation_selection.tsv` — **72-read down-selection** (provisional, see §5)
- Scripts: `classify_candidates.py`, `select_validation.py`, `01_subset_chr5.sh`,
  `02_trim_align.sbatch`, `03c_correct_noann.sbatch`
- Refs reused from `/scratch/users/kevinroy/sumner_lab/references/` (GRCh38_chr5.fa,
  gencode.v44.basic.chr5.gtf, gencode.v44.basic.chr5.introns.gtf.gz, polyasite_chr5.bed)

**M1:** `~/sumner_scratch/{classify_candidates,select_validation}.py` + pulled TSVs
(local scratch, NOT Drive-synced — these are analysis scripts, not committable artifacts).

**Patient IGV bundle (local-only, never commit):** `~/igv_data/sumner_lab/`.

## 5. Remaining goals (in order)

> ### ⭐ HIGH PRIORITY — train a human empirical error model (DRS + cDNA)
>
> Mirror the yeast empirical HP/STR penalty tables for human, using human DRS **and**
> cDNA data. The yeast model lives at
> `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/`
> (`penalty_scores.tsv` = DRS, `penalty_scores_cdna.tsv` + `_umi1/_umi2/_umi3plus`,
> `penalty_scores_qsrev.tsv`, `str_penalty_scores.tsv`, `junction_overhang_table.tsv`).
> Produce the `homo_sapiens/penalty_tables/` analogs (at minimum DRS + cDNA, plus the
> per-UMI-depth cDNA variants if UMI cDNA is used).
>
> - **Method:** `scripts/calibration/empirical_cigar_error_profiler.py` derives per-`(op,
>   base_class AT/CG, HP-length)` error rates at positions where **all aligners agree**
>   (multi-aligner consensus), normalized to the HP=1 substitution baseline per
>   `docs/EMPIRICAL_HP_PENALTY_SCORING.md`. Consumer/loader is
>   `rectify/core/splice/hp_penalty.py` — confirm it resolves per-organism table paths
>   (and that `correct`/refinement actually load the human tables once present).
> - **Inputs:** human DRS + human cDNA. SG-NEx carries both directRNA and PCR-/direct-cDNA
>   for the same cell lines (A549/Hct116/…); the Sumner cohort has DRS + PCR-cDNA. Use
>   public cell-line data for anything shipped; patient data prototyping-only.
> - **Hard dependency on the re-align:** the profiler trusts multi-aligner *consensus*, so it
>   must train on alignments produced with `--max-intron 500000` + the chrom fix — **not** the
>   current `-G 5000` BAMs, whose consensus is intron-cap-contaminated (§3.2). Sequence this
>   after goal 1.
> - **Watch-outs:** (a) chr5 alone is likely too thin for rare HP/STR contexts — prefer
>   multi-chromosome/full-genome (the yeast model used the whole genome, 9.7M reads);
>   (b) **error models are chemistry-specific** — SG-NEx is **R9.4**, whereas the yeast table
>   and the Sumner patient data are **R10.4.1**. An SG-NEx-trained table will not transfer
>   cleanly to R10.4.1 Sumner production; train the production DRS table on R10.4.1 human data
>   and record basecaller + chemistry in the table provenance.

1. **Re-align A549 with `--max-intron 500000`** from `trim/a549_chr5_trimmed.fastq.gz`
   (4-aligner panel; trim step + read set unaffected). Then re-run `correct`
   (without `--annotation`, or *with* it once `961c844` is validated) → `classify`
   → `select`. **This supersedes the current `a549_validation_selection.tsv`** for
   the junction categories (Cat3/5/6/7/9). Cat1/2/8 are intron-independent and
   largely unaffected.
2. **Decide: re-align patient Sumner samples with `--max-intron 500000`** — their
   BAMs (`01_align_array.sh`) have the same 5 kb cap, affecting downstream
   isoform/splicing analysis, not just this validation effort.
3. **Kevin's visual IGV vetting** of the re-derived 72-read set — especially Cat7
   (non-canonical junctions, artifact-prone) and Cat5. Needs an A549 IGV bundle
   (model on `~/igv_data/sumner_lab/` patient bundle).
4. **Build the committed artifact**: extract the vetted reads × 4 aligners, apply
   XV/XG (category/gene) tags, add to the RECTIFY validation set, commit publicly.
   No patient reads.
5. **Commit the doc/code changes from this session** (currently uncommitted working
   tree — see below).

## 6. Working-tree state (uncommitted at handoff)

This session edited (not yet committed): `docs/aligners/minimap2.md` (expanded from
stub — `-G`/`--max-intron` pitfall), `docs/ALIGNER_RECOMMENDATIONS.md` (reconciled
gapmm2 "Why" to the profiled terminal-refinement-loop finding), `AGENT_FIXES.md`
(chrom-fix entry, committed at `99ff250`; re-check for drift). Other working-tree
changes (`HANDOFF.md`, `dev/PERF_AUDIT*.md`, `junction_scoring.py`, `docs/aligners/{mapPacBio,gapmm2,deSALT}.md`,
`verify_dualsite_correct.sbatch`, test files) are **concurrent work by other
sessions** (esp. the `961c844` 3'SS-rescue perf fix) — do not assume they're mine;
stage surgically by explicit path.

---

## 7. Addendum — validate the pool-filter / 3'SS-rescue changes on human (perf/splice workstream, 2026-05-25)

Added by the RECTIFY perf/splice workstream (the `961c844` dual-site lineage in §3.3).
**Pick this up after commit `ed3df74`** (`perf(splice): pool anchor floor +
cross-family concordance relaxation`) is on `drs-validation-rebuild` and pulled
to Sherlock.

### Context — what changed, and why human is the real test

3'SS-rescue candidate generation is now controlled at two stages:

- **Pool-build (hard, low floor + support relaxation):** an observed junction enters the
  rescue candidate pool only if some read crosses it with a clean **≥10 bp exon anchor on
  both flanks** (after adjacent D/N-op collapse) **OR** the junction has sufficient
  cross-aligner / cross-read support. This drops gapmm2-style tiny-anchor splits
  (e.g. `4M 4250N 223M`) without losing real-but-short-anchor junctions. Lives in
  `junction_scoring.py`. Its only job is performance (keep junk out of the candidate set).
- **Consensus (soft, calibrated, never disqualifying):** existing `_add_chimera_flag` —
  length-dependent overhang (`junction_overhang_table.tsv`, calibrated 22 bp for yeast) +
  `hp_ed / (left_ov + right_ov)` ratio + short-intron support relaxation.

A **periodicity/complexity dimension** (motif-agnostic; flags anchors that are
homopolymer / 2–3 bp repeats — *not* splice motifs) was designed for the consensus soft
gate but **deferred** because it is a **no-op on yeast**:

> Yeast 5-aligner derivation (10k reads): worst-flank periodicity ≥0.9 fired on **0/210**
> gapmm2-only artifacts and **0/139** real concordant junctions. The gene-dense yeast
> genome has essentially no long low-complexity anchors. The yeast discriminator is the
> 10 bp clean-anchor floor (caught **98%** of gapmm2-only artifacts) — but it also clips
> **20%** of real concordant junctions, which is exactly why the pool floor carries the
> support-relaxation.

Human (repeats, segmental dups, long introns) is where periodicity should earn its keep.
Helpers are ready: `_seq_periodicity`, `_junction_worst_flank_periodicity` in
`junction_scoring.py`.

### Tasks (in order)

1. **Re-enable `--annotation` on human** — this is the direct test of §3.3. With the
   dual-site interval-tree pool index, confirm `correct --annotation` is now tractable on
   A549 (was ~475× slower / 2 h timeout). If so, retire §2's "run without `--annotation`"
   workaround for human.
2. **Output-equivalence decider for the pool floor:** run A549 `correct` (re-aligned at
   `--max-intron 500000`, §5.1) **with vs without** the pool floor. Confirm (a) it drops
   the chimeric long-intron gapmm2/deSALT artifacts, and (b) it causes **no real-junction
   correction loss** (concordant/annotated junctions unchanged). This is the agreed
   criterion for keeping the floor hard-with-relaxation vs softening further.
3. **Periodicity earn-its-keep test on human:** **extend the committed
   `rectify/core/splice/calibrate_junction_overhang.py`** (the canonical concordance-based
   derivation — do *not* re-derive from scratch) to also emit `_junction_worst_flank_periodicity`
   per junction. On the human 4-aligner panel (no gapmm2), compare the distribution for
   real (concordant) vs spurious junctions. If it separates → derive a human **T2** and wire
   `_junction_worst_flank_periodicity` into `_add_chimera_flag` (genome is already available
   in `merge_corrected_tsvs`; `rep_df` carries `chrom`; thread through both call sites,
   ~lines 1060 & 1244). If it's a no-op on human too → drop the dimension entirely.

### Guardrail (load-bearing principle)

**No splice-motif gating in rectify.** Prune candidates by anchor quality / read support,
never by canonical dinucleotide (GT-AG), to preserve unbiased discovery of non-canonical
junctions. The periodicity check is complexity-based (repeat structure), not motif-based —
keep it that way.
