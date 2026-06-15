# STATUS — Human DRS validation read set (A549) — LIVING DOC

**Keep this current.** Supersedes the point-in-time `HANDOFF_2026-05-25_human_drs_validation.md`
for the *validation-read-set* deliverable (that handoff's §1–4 context is still good; its §5
"remaining goals" is what this doc tracks). Last updated **2026-06-12**.

## What this deliverable is

The human analog of the yeast committed validation set (`rectify/data/validation/`): ~8
reads/category × 9 RECTIFY correction categories (4/strand), built from **public** SG-NEx
**A549** chr5 direct-RNA (no patient reads). It is a **correction-behavior regression
fixture**, not an error-rate model — so chemistry (A549 = R9.4) is a provenance caveat, not a
disqualifier; native poly-A (A549 is native) is what matters for the 3′-end categories.

## Decisions locked (Kevin, 2026-06-12)

- **Substrate = finish A549** (not rebuild on GM12878 RNA004 — the only public GM12878 RNA004
  is IVT with synthetic ~15 bp poly-A, weaker for cat1/2/8; native RNA004 GM12878 is the
  un-emailed Miten-Jain item).
- **Kevin visually vets cat5 (chimeric) + cat7 (alt-splice) in IGV** before commit (artifact-prone;
  not self-certifiable). Deliverable handed back this session = re-derived candidate set + IGV bundle.

## The bug that made the prior selection invalid (now understood)

The provisional 72-read selection (`a549_validation_selection.tsv`, May 25) was built on
alignments run at the **yeast intron ceiling `-G 5000`** — confirmed: the `alignments/` BAM
`prov.json` shows `rectify align` was called with **no `--max-intron`** (and `--organism
homo_sapiens` does NOT raise it — that was the bug). Junction categories (cat3/5/6/7/9) are
intron-cap-contaminated. **Transcript search (2026-06-12) confirmed the re-align was deferred to
"the Human DRS agent" and never executed** — so this is genuinely outstanding work, not a redo.

**Key simplification:** `--max-intron` only feeds minimap2 (`-G`) and gapmm2 (dropped). uLTRA and
deSALT are intron-model-independent → their existing BAMs are reused unchanged. So the re-align is
**minimap2-only** (minutes), not a 5–6 h 3-aligner redo.

## Code-state fix (load-bearing for a regression fixture)

Sherlock's rectify (`/oak/.../software/rectify`, 0.9.0 rsync copy, `git_sha: null`) was **behind
M1 HEAD `9f613a6` on `walkback.py`** — missing `1b1db38 fix(walkback): bypass lone-A early-exit
guard for +strand DRS reads`, which directly governs cat1 +strand poly-A walkback. The stale copy
dated **May 21** (~3 weeks behind). **Fixed 2026-06-12:** synced `walkback.py` to M1 HEAD
(md5 `50ce85c…`); backup at `walkback.py.bak_pre_head_sync_20260612`. Because cat1/2/8 picks
depend on walkback, **all 9 categories are re-derived from the new corrected run** (not just the
junction cats).
- ⚠ `correct_command.py` on Sherlock is *ahead* of M1 HEAD with an uncommitted `--netseq` arm
  (every delta is `is_netseq`-gated → DRS path byte-identical, so left untouched). **Needs
  reconciliation to M1** (M1 also has WIP on this file). Not this deliverable's job — flagged.

## Pipeline (all on Sherlock; A549 data lives there)

Driver script: `/scratch/users/kevinroy/rectify_human_validation/sgnex_a549/a549val_realign_correct.sbatch`
(owners partition + AVX-512 constraint, idempotent per-stage skip-checks).

1. **Stage A** — minimap2 re-align @ `--max-intron 500000` → `alignments_500k/` (reuse existing
   uLTRA + deSALT BAMs via symlink).
2. **Stage A2** — subsample all 3 aligners to 25% with name-hash seed `samtools -s 7.25` (same
   reads across aligners; A549 cancer-line pileups choke per-region `correct` at full depth) →
   `alignments_500k_sub25/`.
3. **Stage B** — `rectify correct` 3-aligner consensus, **no `--annotation`** (the working human
   path; `--annotation` is ~475× slower / 2 h timeout on human) → `correct_500k_noann/a549_chr5_corrected.tsv`.
4. **Stage C** — `classify_candidates.py` (→ `a549_candidates_500k.tsv`) then `select_validation.py`
   (→ `a549_validation_selection_500k.tsv`). Re-derives all 9 categories. Fast (pysam cross-walk).
5. **Stage D** — build A549 IGV bundle for Kevin's cat5/cat7 vetting.

Stage C invocations (run after Stage B):
```
python classify_candidates.py --tsv correct_500k_noann/a549_chr5_corrected.tsv \
  --genome $REF/GRCh38_chr5.fa --gencode $REF/gencode.v44.basic.chr5.gtf \
  --polyasite $REF/polyasite_chr5.bed --aligner-dir alignments_500k_sub25 \
  --prefix a549_chr5_trimmed --out a549_candidates_500k.tsv
python select_validation.py --inventory a549_candidates_500k.tsv \
  --gtf $REF/gencode.v44.basic.chr5.gtf --out a549_validation_selection_500k.tsv
```

## Current state (2026-06-12)

- Walkback synced ✅. **SLURM job 29321539 COMPLETED** (7.5 min): minimap2 re-align @ -G 500000 →
  25% subsample (mm2 79722 / uLTRA 79713 / deSALT 79722) → correct (79,716 reads, no `--annotation`).
  Correction signal healthy: A-tract walking 46.1%, indel 10.0%, 5′ rescue 0.1%.
- **Stage C done.** classify → 9 categories all well above 4/strand; select → 72 reads (8/cat, 4/strand,
  8 loci each): `a549_candidates_500k.tsv`, `a549_validation_selection_500k.tsv`.
- **cat7 selector fixed (length→support).** Was ranking by longest intron, which enshrined
  single-aligner 277–453 kb mismapping artifacts (cat7 was the only junction category with no
  multi-aligner check — the others, cat5/6/9, cross-check via `same_intron` and were clean). Now
  ranks by cross-aligner support (`same_intron` agreement count; guardrail-compliant — support not
  motif). New cat7 picks: all **3/3-aligner-supported**, 48 bp–99 kb, real genes (TNPO1/BTF3/RPL37/
  ERCC8/SSBP2). Fillability checked honestly: 16+/17− at full support, 113+/120− at ≥2 — ample, no
  force-fill. Edited scratch scripts `classify_candidates.py` (cat7 block) + `select_validation.py`
  (cat7 strength); pulled to M1 `/tmp` copies. **Latent same flaw in the yeast selector** (length
  bound masked it) — note only, not fixed this session.
- cat4 eyeballed (advisor): legit — it captures false N-junctions over poly-A (term=AAAAAAAAA/TTTTTT),
  single-aligner by design. cat5 verified: compact 2-intron reads, real genes.
- **Stage D DONE — vetting bundle on M1 at `~/igv_data/a549_validation/`** (local-only, public data):
  `aligners/val.{minimap2,uLTRA,deSALT}.bam` (72 reads each), `VETTING_cat5_cat7.tsv` (16 reads +
  padded IGV loci), `selection_all72.tsv`, `load_igv.sh`, `README_VETTING.md`, `scripts/` (the
  updated classify/select + sbatch). Sherlock copy at `sgnex_a549/vetting_bundle/`.
### IGV vetting round 1 (2026-06-12, live with Kevin) — found a SECOND selector bug

IGV setup gotchas resolved: IGV had the **SGD/yeast genome** loaded (human chr5 BAMs showed
nothing); switching to hosted hg38 gave gene track but **no sequence bases** (reads rendered `=`).
Fix = load the **local `~/igv_data/sumner_lab/reference/GRCh38_chr5.fa`** (exact alignment
reference) + the matching `gencode.v44.basic.chr5.sorted.gtf.gz` track. Port has no `genome`
command → genome load is a manual IGV step.

**cat5 selector bug (analogous to cat7).** cat5 #1 (UQCRQ) and #2 (IK) both REJECTED on sight:
minimap2 silently resolved BOTH introns, while the `same_intron` (±4 bp) holder check credited
intron-2 to deSALT alone (deSALT placed it a few bp off) — a boundary disagreement (cat9 flavour)
masquerading as chimeric partitioning. **Fixed:** the "no single aligner spans both" exclusion now
uses `overlaps` (loose) not `same_intron`. After the fix, genuine cat5 = **3 reads, all +strand,
0 −strand** on the 25% subsample → true chimeric (no aligner resolves both introns) is genuinely
rare because minimap2 resolves most multi-intron reads alone.

**Kevin decisions (2026-06-12):** (1) cat5 → **re-search at FULL depth** (un-subsample; 4× the
rare-category yield, may surface −strand cat5); (2) **spot-vet EVERY category** before commit (2/2
vetted categories had selector bugs → don't trust the rest).

### Full-depth re-derivation IN FLIGHT

- **SLURM job 29352347** (owners): `rectify correct` on full-depth `alignments_500k/` (minimap2
  @500k + uLTRA/deSALT symlinks, 319k reads), no `--annotation`, `--checkpoint-dir` (per-region
  resume = preemption-resilient) → `correct_500k_fulldepth/`. Background waiter on M1.
- Full-depth correct DONE (job 29352347, 319,266 rows). Re-classified/selected; spot-vet began.
- ⏸ The sub25 selection (`a549_validation_selection_500k.tsv`) is SUPERSEDED by the full-depth pass.

### IGV vetting (live, 2026-06-12/13) — TWO more selector bugs found + a quality issue

Spot-vetting EVERY category (not just cat5/cat7) paid off. Total **4 selector fixes**, all
Kevin-validated in IGV (all in the scratch `classify_candidates.py`/`select_validation.py`):
1. **cat7 length→support** — was ranking by longest intron (450 kb mismapping artifacts); now by
   cross-aligner `same_intron` support (3/3-supported picks).
2. **cat5 same_intron→overlaps** — "no aligner spans both" used ±4 bp, so minimap2 silently
   resolving both introns (UQCRQ/IK rejected on sight) looked like clean partitioning; now loose
   overlap. Genuine cat5 is sparse (24 on chr5, +strand low-quality) → drove Round 2.
3. **cat8 poly-A-run preference** — picks were degraded no-poly-A reads. `polya_len` is UNRELIABLE
   (ligation-based DRS selects for poly-A+, but the basecaller mangles the tail into low-Q G-poor
   garble → real tail reads as polya_len=0; Kevin's insight, confirmed: read `ca4e8975` clip
   *starts* with a 23 bp poly-A run then degrades). Now prefers a recognizable A-run in the 3′
   soft-clip, then TPM. Picks now pA 0.50–0.80 at high-TPM CPAs.
4. **Quality floor (all cats)** — read quality is continuous/unimodal (median 87.8% identity, NOT a
   junk subset). Selection prefers ≥90% identity with relaxation; added per-read identity (NM/aln)
   to classify. 8/9 cats now ≥90%; cat5 the exception (genuine chimeric = inherently low-quality).

Confirmed in IGV: cat5 genuine chimeric (fix works), cat8 real-CPA+poly-A. cat4 legit (false
N-junction over poly-A, single-aligner by design). **Decided: committed set preserves untrimmed
dorado_source + trim metadata** (poly-A recoverable, mirrors yeast set).

### ROUND 2: multi-chromosome expansion (IN FLIGHT) — Kevin: "align a few more chromosomes"

cat5 +strand can't fill 4 clean reads on chr5 alone. Enlarging the pool with **chr1, chr11, chr17,
chr19** from the same 4 A549 runs (only those 4 exist). 2.42M reads streamed → **subsampled to ~1M**.
- **classify/select are now CHROMOSOME-AWARE** (were chr5-hardcoded): per-chrom gencode/PolyASite/
  genome.fetch, chrom-keyed gene_at. **Regression-tested: reproduces chr5 exactly** (cat5=24, 8/8).
- Resources on Sherlock `sgnex_a549/morechrom/`: `morechrom_ref.fa` (chr1/11/17/19 from
  GRCh38.primary_assembly), `gencode.morechrom.gtf`, `polyasite_morechrom.bed` (atlas normalized
  1→chr1, repPos from col4). Streamed BAM `morechrom.merged.bam` (2.42M; the .fastq.gz was
  login-truncated to 531k — regenerate from BAM in-job).
- **SLURM job 29463887** (owners 16c/96G/10h): subsample→trim-polya→fastq→align (minimap2@500k +
  uLTRA + deSALT), per-aligner skip-checks for preemption resilience.
- **NEXT when align lands:** correct (full depth, no-annotation, checkpoint-dir) → classify morechrom
  (chrom-aware) → **merge chr5 + morechrom inventories** → select (quality floor) → extract →
  re-vet cat5 (should now fill 4 clean +strand) → finalize 9 cats → tag (XV/XG) + commit. Committed
  set becomes **multi-chromosome** (chr1/5/11/17/19), no longer chr5-only.

## Remaining after this session (still open)

1. Kevin IGV-vets cat5/cat7 on the re-derived set.
2. Extract vetted reads × 3 aligners, apply **XV (category) / XG (gene)** tags, assemble the
   committed artifact mirroring `rectify/data/validation/` layout (per-aligner BAMs in `aligners/`,
   `corrected_reads.tsv`, PROVENANCE.json **recording M1 HEAD sha `9f613a6`** + R9.4 caveat,
   VALIDATION_READS.md). Commit publicly — **no patient reads**.
3. Reconcile the Sherlock `--netseq` `correct_command.py` divergence back to M1.

## Provenance notes for the committed artifact

- Source: SG-NEx A549 directRNA chr5, public. Chemistry **R9.4** (caveat: fixture exercises
  correction *behavior*, not error rates). Trimmed input `trim/a549_chr5_trimmed.fastq.gz`.
- Aligners: minimap2 @ `--max-intron 500000` (re-run 2026-06-12) + uLTRA + deSALT (reused, May 25,
  intron-model-independent of `--max-intron`). 25% name-hash subsample (seed 7).
- Code: rectify M1 HEAD **`9f613a6`** (walkback synced to this on Sherlock). Record the sha
  manually — the rsync copy can't stamp `git_sha`.
