## [2026-07-25] 🔴🔴 TRAP (reference) — a FASTA/GFF **contig-naming mismatch** makes the whole panel annotation-blind, silently

**Found by:** Chanfreau 5′-rescue session, 2026-07-25 (`planning/460`). **Not a code bug — a reference-staging
trap, and the most dangerous one yet, because every job exits 0 and every BAM looks normal.**

**Symptom.** Nothing fails. The one and only visible signal is a single line in minimap2's stderr:

```
[mm_idx_bed_read_merge] read 0 introns, 0 of which are non-redundant
```

immediately after rectify logs `Extracted 402 introns` from the GFF. Junction hints are silently discarded.

**Root cause.** The genome FASTA and the annotation use different chromosome vocabularies:

| file | seqids |
|---|---|
| `S288C_reference_simple.fasta` (and every BAM built from it) | `I, II, … XVI, ref\|NC_001224\|` |
| SGD `saccharomyces_cerevisiae_R64-5-1_*.gff` / `.gtf` | `chrI, chrII, … chrXVI, chrmt` |

Both ship **side by side in the same lab reference directory**, so the pairing looks correct. minimap2 builds the
junc-BED from the GFF (`chrI…`), then matches it against index contigs (`I…`) and finds nothing. uLTRA's
`ultra_norm.gtf` inherits the same wrong names.

**Why it matters far beyond alignment quality:** rectify's 5′ soft-clip rescue draws its candidate pool from
*annotated junctions ∪ cross-aligner novel junctions*. With the annotated half empty, the pool is **novel-only** —
i.e. the unfiltered short-anchor substrate that `min_junction_anchor_bp` exists to gate. Any rescue or
junction-pool measurement taken on a mismatched reference is measuring the wrong thing.

**Detection (do this once per reference dir, before any run):**
```bash
diff <(cut -f1 "$REF.fai" | sort -u) \
     <(awk '/^##FASTA/{exit} !/^#/{print $1}' "$GFF" | sort -u)   # must be empty
```
Or, after a first minimap2 chunk: `grep mm_idx_bed_read_merge <log>` — a **0** there means annotation-blind.

**Fix.** Normalise the annotation to the FASTA's vocabulary (whole-token replacement only — a naive
`s/chrX/X/` will corrupt `chrXI`), covering feature column 1, `##sequence-region` pragmas, and `>id` headers in
any embedded `##FASTA` block; then delete/regenerate `*.ultra_norm.gtf` so uLTRA rebuilds it. Keep the originals.
Chanfreau copy of the converter: `planning/460` (`460_fix_ref_contigs.py`).

---

## [2026-07-25] 🔴 TRAP (align) — `minimap2 -y` on a non-SAM FASTQ comment aborts samtools, and the error is MASKED as SIGPIPE

**Found by:** Chanfreau 5′-rescue session, 2026-07-25 (`planning/460`).

**Symptom.** `rectify align` reports `minimap2 failed (exit -13)` (SIGPIPE) a second or two after
`loaded/built the index`, with retries failing identically. The stderr shown is minimap2's, and it looks healthy.

**Root cause.** `run_minimap2()` (`rectify/core/align/multi_aligner.py`) passes **`-y`**, which copies the FASTQ
**comment** verbatim into the SAM record as aux fields. That is correct for the cDNA/UMI pipeline, whose comments
are tab-separated SAM aux tags (`RN:i:5`, `pt:i:25`) — and wrong for any FASTQ whose comment is free text, e.g. a
raw SRA header `@SRR32518273.1 e464c9e6-290a-…/1`. `samtools sort` then aborts:

```
[E::aux_parse] unrecognized type '4'
samtools sort: truncated file. Aborting
```

killing the pipe, so minimap2 dies of SIGPIPE. **The real error is invisible** because `run_minimap2` checks
`sam_output.returncode` (minimap2) *before* `sort_proc.returncode` (samtools) and raises on the former, discarding
samtools' stderr.

**Reproduce / confirm** — run the pipeline by hand and capture BOTH stderrs; `samtools sort` will name the problem:
```bash
minimap2 -ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD -y -t 2 \
    --junc-bed $BED $GENOME $FASTQ 2>mm2.err | samtools sort -n -@2 -o /dev/null 2>sort.err
```

**Workaround today.** Feed minimap2 only FASTQs that went through `rectify split` (it rewrites headers to
`RN:i:<n>`), or strip the comment. **Chunked runs are safe; a run pointed straight at an SRA/dorado FASTQ is not.**

**Suggested fixes (not yet applied).** (a) report samtools' stderr when minimap2's rc is `-13`, since a broken
pipe means the *consumer* is the real failure; (b) sniff the first FASTQ comment and only pass `-y` when it parses
as SAM aux (`^[A-Za-z][A-Za-z0-9]:[AifZHB]:`).

---

## [2026-07-25] 🔴 TRAP (align) — concurrent bbmap index build: `java.io.EOFException` on ONE array task

**Found by:** Chanfreau 5′-rescue session, 2026-07-25 (`planning/460`).

**Symptom.** In a mapPacBio array, exactly one task (usually the first) dies while its siblings pass:
```
NOTE:  Ignoring reference file because it already appears to have been processed.
Exception in thread "main" java.lang.RuntimeException: java.io.EOFException
    at dna.ChromosomeArray.read(ChromosomeArray.java:63)
    at align2.BBMapPacBio.loadIndex(BBMapPacBio.java:392)
```

**Root cause.** All array tasks share one `path=<ref_dir>/bbmap_index`. The first task to arrive starts building
it; a sibling sees the freshly-written `summary.txt`, concludes the index is ready, and reads a **still-being-
written `chr1.chrom.gz`** → EOF. Later tasks read the finished file and succeed, so the failure looks
input-specific ("chunk_000 must be corrupt") when it is purely a startup race.

**Fix.** **Pre-build every shared index serially before submitting any array** — bbmap's `path=` index, deSALT's
`desalt_index`, and uLTRA's `*.ultra_norm.gtf` (which uLTRA writes next to the annotation and which races the same
way). One throwaway single-task alignment on a few thousand reads warms all of them. At 21-sample fan-out scale
this is a blocker, not a nuisance.

---

## [2026-07-25] ✅ FIXED (run-all) — single-sample `--chunked-alignment` SKIPS the entire consensus chain (so 5′ rescue never runs)

**Found & fixed by:** Chanfreau 5′-rescue session, 2026-07-25 (`planning/460`).
**Fix:** `rectify/core/commands/run/chunked_batch.py`, single-sample branch of `_generate_chunked_pipeline`.

**Symptom.** The chain runs to completion, then `run_correct_analyze.sh` dies with
`ERROR: consensus BAM not found: <out>/consensus/<sample>.consensus.bam`.

**Root cause — two independent defects, both required for a working chain.**
1. **The middle of the pipeline was never wired.** The generated `submit_pipeline.sh` chained only
   `split → mapPacBio/others → merge_aligners → correct_analyze`. `generate_alignment_scripts()` already returns
   and writes `prescan`, `correct` (per aligner), `chunk_merge`, `final_merge`, `consensus_per_chunk` and
   `merge_consensus_chunks` — none were submitted. **Consequence: no consensus BAM, and therefore no 5′
   soft-clip rescue at all**, because the rescue happens inside consensus scoring (`consensus/scoring.py`
   `_rescue_5prime_softclip`), not in `correct`.
2. **The consensus path did not match.** `run_correct_analyze.sh` looked in `<out>/consensus/`, while
   `run_merge_consensus_chunks.sh` writes to `<out>/chunks/consensus/` (its `OUTDIR` is the chunks dir). Fixing
   (1) alone still fails with the identical error.

Same family as the `KeyError: 'merge'` entry below: the single-sample branch is systematically behind the
multi-sample one. ⚠️ **Still open:** the multi-sample branch stops at `final_merge` and likewise never produces a
consensus BAM — so multi-sample `--chunked-alignment` also cannot rescue.

**Validated:** `pytest -m "not slow" -k "chunk or split or run_command"` → **74 passed, 1 skipped**; regenerated
pipeline inspected (15-stage chain, consensus path matching on both sides).

---

## [2026-07-25] 🔴 TRAP (align) — a reference genome in ANOTHER USER's directory silently fails EVERY uLTRA and deSALT task

**Found by:** Chanfreau lariat/rescue session, 2026-07-25 (`planning/455`, `planning/456`). **Not a code bug —
an environment trap that looks like an aligner bug, and it silently reduces a 4-aligner panel to 2.**

**Symptom.** A `run-all --chunked-alignment` fan-out completes with `exit_status 0` for every minimap2 and gapmm2
task and `exit_status 1` for **every** uLTRA and deSALT task (here: tasks 1-8 OK, 9-16 all failed, 3 retries each).
No BAMs for those aligners. The consensus then runs on a **2-aligner** pool, so the cross-aligner junction pool —
the whole basis of 5' soft-clip rescue — is silently crippled. `qstat` is empty and the top-level job looks fine.

**Root cause (ONE cause, two messages).** Both aligners must WRITE derived artifacts *adjacent to the reference*:
- uLTRA: `[Errno 13] Permission denied: '<ref_dir>/<genome>.ultra_norm.gtf'` — it normalises the GTF
  (mRNA→transcript, deriving exons) into the reference directory.
- deSALT: `deSALT index not found (or empty) adjacent to <genome.fa>. Looked at: <ref_dir>/desalt_index, ...`

If the reference lives in a **read-only directory owned by another lab member** (ours was
`/u/project/guillom/boli0204/S288C_reference_genome_R64-5-1_20240529/`), uLTRA cannot write its normalised GTF and
deSALT's index can never be built there. minimap2/gapmm2 need no such sidecar, so they pass and mask the problem.

**Fix.** Stage the reference into a **writable** location you own and build the index once:
```bash
REF=/u/project/guillom/kevinroy/refs/R64-5-1        # writable
cp -p <src>/S288C_reference_simple.fasta <src>/*.gff <src>/*.gtf $REF/
samtools faidx $REF/S288C_reference_simple.fasta
deSALT index $REF/S288C_reference_simple.fasta $REF/desalt_index
```
then pass `--genome $REF/... --annotation $REF/...`. **Never point a multi-aligner run at a reference you cannot
write to** — and per `shared/README.md` hygiene, do not write into another user's dir to work around it.

**Detection.** After any multi-aligner run, assert one BAM per aligner per chunk before trusting the consensus:
`for a in minimap2 gapmm2 uLTRA deSALT; do echo -n "$a "; ls <chunks>/aligner_chunks/$a/chunk_*/*.$a.bam | wc -l; done`
A zero for any aligner means the pool is degraded even though the pipeline "succeeded".

---

## [2026-07-25] ✅ FIXED (run-all) — `--chunked-alignment` dies instantly with `KeyError: 'merge'`

**Filed & fixed by:** Chanfreau lariat/rescue session (Kevin Roy), 2026-07-25.
**Fix:** `rectify/core/commands/run/chunked_batch.py` — two stale dict keys
`align_scripts['merge']` → `align_scripts['merge_aligners']` (L267 and the summary print at L365).

**Symptom.** ANY `rectify run-all --chunked-alignment` run aborts in <60 s, after the split/chunk-count step:

```
File ".../run/chunked_batch.py", line 267, in _generate_chunked_pipeline
    merge_script  = align_scripts['merge']
KeyError: 'merge'
```

**Root cause.** `split_command.generate_alignment_scripts()` returns the key **`merge_aligners`**
(alongside `mpb`, `others`, `prescan`, `correct`, `chunk_merge`, `final_merge`,
`consensus_per_chunk`, `merge_consensus_chunks`, `submit`). The single-sample branch of
`_generate_chunked_pipeline` still read the pre-rename name `merge`. Key-set diff confirms `'merge'`
was the ONLY consumed key with no producer, and the **multi-sample branch (~L508) already used the
correct `merge_aligners`** — so this was a rename that was propagated to one call path and not the other.
Scheduler- and organism-independent: it fails before any scheduler branch is reached.

**Validated:** key-set diff (consumed ⊆ returned) now closes; `pytest -m "not slow" -k "chunk or split or
run_command"` → **74 passed, 1 skipped**. Chanfreau `planning/455`.

**Note for callers:** `--chunked-alignment` is the mandatory submission mode for cluster rectify work in the
Chanfreau workspace, so this bug blocked *all* compliant fan-outs. Also heed the runtime warning it prints on
yeast: mapPacBio is in the panel with `min_junction_anchor_bp=0`, which re-opens the spurious-intron gaming
vector — pass `--min-junction-anchor-bp 10` (or `--skip-map-pacbio`) for any junction-sensitive analysis.

---

## [2026-07-22] ✅ FIXED (correct) — `get_reference_sequence()` heap corruption / SIGABRT on spurious mega-intron reads

**Filed & fixed by:** Chanfreau Session-B (Kevin Roy), 2026-07-22. **Fix:** `rectify/utils/alignment.py`
`extract_deletions()` — guard the `read.get_reference_sequence()` call (skip when any N-op > `_MAX_PLAUSIBLE_INTRON_BP`
= 50 kb). One-line-ish, additive; the existing `ref_seq=None` fallback already handles the skip.

**Symptom.** `rectify correct` on a genomic BAM aborts with glibc `*** double free or corruption ***` /
`corrupted size vs prev_size` inside `pysam/libcalignedsegment`, then HANGS. rc=134, **uncatchable** by the
surrounding Python try/except (C-level SIGABRT). Reproduced with `-j 1` AND multithreaded; `--threads 1` and
`--legacy-single-threaded` do NOT help.

**Root cause.** `extract_deletions()` calls `read.get_reference_sequence()` (reconstructs ref bases from the MD tag
over the read's full reference span). On reads carrying a spurious multi-kb `N`-op — 142–159 kb "introns" emitted by
long-read splice aligners (deSALT/uLTRA) as long-range splice artifacts; real *S. cerevisiae* introns are <3 kb —
pysam 0.23.3's C implementation corrupts the heap. Construct/intron-less reads never hit it; ~10% of genomic reads
in the test locus did.

**Diagnosis path (for future ref):** module bisect isolated `--skip-indel-correction`; per-read isolation found the
mega-N-op CIGARs; keeping huge-N (skip only big indels) still crashed → the N-op is the trigger; `extract_deletions`
is the only pysam-C call in that path.

**Validated:** the exact original crash case (`correct` on a construct+chrXIII ROI, all modules) now completes
rc=0, 25,110 reads corrected. **Note:** deeper upstream fix would cap the aligner's max-intron so these spurious
alignments never win consensus; this correction-layer guard is the robust backstop. Full writeup: Chanfreau
`planning/281`.

---

## [2026-07-21] 🔴 CRITICAL BUG (align) — **OPEN, NOT FIXED**: single-aligner output path is sample-keyed, so a per-aligner fan-out silently clobbers itself

**Filed by:** Sumner (JHU) session-C, 2026-07-21 01:40 PT. **Reported by Kevin Roy**, who stated the
invariant: *"There should always be the aligner and sample ID name in the prefix to avoid the
race/clobbering issue."* **Not patched here — see "Why not fixed in this commit" below.**

**Symptom.** Running the 3-aligner panel as one job per (sample × aligner) — the correct pattern, since
uLTRA rebuilds a whole-genome index per invocation and must not be chunked — leaves a
`<sample>_trimmed.rectified.bam` that is **not a consensus**. On the Sumner pilot (`WT_21.8_rep3`) it was
byte-identical to the uLTRA BAM plus **one `@PG` record**: `rectified.bam` 3,692,200,234 B vs
`uLTRA.bam` 3,692,200,211 B; `@PG` chain `samtools.4/samtools.5` vs `samtools.3/samtools.4`. **All three
job logs contain** `Single-aligner output (sorted+indexed): …_trimmed.rectified.bam` — three writers,
one path. Last writer wins (uLTRA, finishing last).

**Cause.** `rectify/core/commands/align_command.py` (deployed Oak copy ~L810; M1 `drs-validation-rebuild`
~L838):
```python
if len(successful_aligners) == 1 or getattr(args, 'no_consensus', False):
    logger.info("Skipping consensus selection (single aligner or --no-consensus)")
    if len(successful_aligners) == 1:            # ← fires REGARDLESS of --no-consensus
        rectified_bam = args.output_dir / f"{prefix}.rectified.bam"
        samtools sort -@ threads -o rectified_bam single_bam
        samtools index rectified_bam
    return 0
```
`prefix` is **sample-level** (`:416`, derived from the reads stem) with **no aligner component** ⇒ every
single-aligner job for a sample targets the same path.

**🔴 `--no-consensus` does NOT prevent it.** The outer guard is `len==1 OR no_consensus`, but the inner
guard is `len==1` alone. A per-aligner fan-out always has `len==1`, so the block fires on every job.
A caller passing `--no-consensus` to opt out of this behaviour **still gets it.**

**Two costs.** (1) *Silent wrong data*: a consumer of `rectified.bam` gets **single-aligner** output while
appearing to use the panel — and it fails in the worst way, since the file exists, is a plausible size,
and passes `samtools quickcheck`. (2) *Pure waste*: each job pays an extra `samtools sort` + `index` of a
3.4–4.1 GB BAM to produce a file that is immediately clobbered. On the Sumner 15-sample × 3-aligner
panel that is **45 redundant multi-GB sorts.**

**FIX (one line).** `f"{prefix}.rectified.bam"` → **`f"{prefix}.{aligner_name}.rectified.bam"`**
(the single surviving key of `successful_aligners`). Optionally also gate the inner block on
`not args.no_consensus` so the flag means what callers expect.
⚠️ **Renaming alone does NOT fix it:** `drs-validation-rebuild` already renamed this to
`.multialigned.bam` (~L843) and it is **still sample-level**, so the identical race exists there under a
new filename. **The aligner component is the load-bearing part.**

**Why not fixed in this commit.** At filing time **23 of 45 jobs of a live Sumner panel were RUNNING**;
changing the deployed copy mid-panel would make early and late jobs behave differently. `align_command.py`
is also traversed by COMPASS on this shared branch. **Apply deliberately once the panel drains, and
re-run the two `--short-read` PE tests plus `tests/test_compass_aligner_cmds.py`.**

**Interim handling for any consumer:** treat `*_trimmed.rectified.bam` / `*.multialigned.bam` produced by
a per-aligner fan-out as **garbage**; never glob it into a downstream stage. A real consensus exists only
where `rectify consensus` was run over the per-aligner BAMs (those ARE correctly named).

**Generalizable rule:** any output path written by a parallel fan-out must be keyed on **every dimension
the fan-out varies**. A path keyed on only some of them is a silent race.
(Sumner cross-ref: `planning/sma_drs_program/30_CORRECTIONS_REGISTRY.md` §0.18, §0.25.)

## [2026-07-17] CHANGE (align): bbmap `intronlen` 10/20 → 40 — stop sub-40bp deletions being mislabeled as introns

bbmap relabels any reference gap ≥ `intronlen` bp as an N-op (fake intron). It was too low
(`intronlen=10` in `run_map_pacbio`, `=20` in `run_bbmap`, `rectify/core/align/multi_aligner.py`), manufacturing
fake introns from short deletions: **55.6% of the COMPASS panel's junction set was bbmap artifact** (97.4%
motif-less, median 51 bp; Chanfreau `planning/303`). Raised both to **40** (S. cerevisiae min intron ~40 bp;
"allow down to 40 to be safe"). HONEST CAVEAT — this does NOT fully de-contaminate: `intronlen=40` stops only the
SUB-40bp deletion noise; ≥40bp fake introns (incl. the median-51bp artifacts) still pass bbmap and are handled
DOWNSTREAM by junction QC (100%-persistent gap + zero intron-body coverage + no canonical motif ⇒ genomic-deletion
allele, not splicing). Supersedes the 2026-05-25 mapPacBio note below (which set `intronlen=10`–`20` so mammalian
introns emit as N-ops); for the yeast COMPASS panel, 40 is the correct floor. (this commit)

## [2026-07-17] BUG: cdna-analyze crashes on multi-aligner BAMs — walkback chrom_seq[rp] IndexError

`rectify cdna-analyze` on a multi-aligner `multialigned.bam` (minimap2 + deSALT + gapmm2 + mapPacBio, from the
cDNA full-pool run-all, planning/259) crashes: `IndexError: string index out of range` at
`rectify/core/correct/walkback.py` `walkback_3prime_with_qpos` (`ref_base = chrom_seq[rp].upper()`). Cause:
deSALT/gapmm2/mapPacBio can report a reference position AT or PAST the chromosome end, and the two walkback
functions indexed `chrom_seq[first_rp]` / `chrom_seq[rp]` without bound-checking (unlike the guarded scan at
line ~238/556 which already does `0 <= x < len(chrom_seq)`). FIX: bound-guard `first_rp` (only fire the terminal
gate when in-range) and `continue` on any out-of-range/None `rp` in the scan loop, in BOTH walkback functions.
Defensive-only (skips positions that would crash; cannot change correct-path behavior). H2 overlay `rectify_patched_250` patched identically to unblock the run-all analyze step.
SECOND bug exposed after that fix: `get_aligned_pairs(matches_only=False)` on a pathological deSALT/mapPacBio
alignment with a huge reference span → `MemoryError` (pysam materialises every intron/deletion gap). FIX:
`matches_only=True` (behaviour-identical — the code already filtered to both-not-None positions) so pysam
never builds the giant list. Both walkback functions.

## [2026-07-20] OPS: H2's shared rectify is a NON-GIT 0.9.0 tree missing BOTH cDNA crash fixes — check `rectify.__file__`, not `__version__`

Found while staging the prp5/prp28/snp1 ONT PCR-cDNA run-all (Chanfreau `planning/377`).
**`/u/home/k/kevinroy/software/rectify` is NOT a git repository** (`fatal: Not a git repository`) — it is a
plain copy of the frozen **0.9.0** release, and the shared env
`/u/project/guillom/shared/envs/rectify` imports rectify FROM IT. So on H2 you cannot `git pull`, and the
HEAD is unknowable. Verified it lacks BOTH documented cDNA-path fixes: `walkback.py` has
`matches_only=False` x3 / `matches_only=True` x0 (⇒ the pysam MemoryError + `chrom_seq[rp]` IndexError are
LIVE), and `hp_penalty._hp_run_length` has no bounds guard (⇒ the whole-chunk_merge IndexError is LIVE).
**Any cDNA multi-aligner run-all on the stock H2 env will hit these.**

**Trap:** `rectify.__version__` reports `0.9.0` in the patched overlays too (the string was never bumped) —
**identify a tree by `rectify.__file__`, never by the version string.**

**Known-good source trees on H2** (carry the walkback fix, `matches_only=True` x4):
`/u/project/guillom/kevinroy/rectify_patched_250/rectify`, `…/rectify_umi_run/rectify`.
**Stale (0.9.0):** `/u/home/k/kevinroy/software/rectify`, `/u/project/guillom/shared/software/rectify`.

**Workaround used:** isolated overlay `/u/project/guillom/kevinroy/rectify_377` (copy of patched_250 + the
hp_penalty guard + the corrected_consensus per-read try/except), selected via
`export PYTHONPATH=/u/project/guillom/kevinroy/rectify_377:$PYTHONPATH` on top of the shared env.
⚠️ Note those last two fixes live as **uncommitted WIP on M1**, so a `git archive HEAD` deploy would MISS
them — deploying "HEAD" is not sufficient to get a working cDNA pipeline on H2.

**Env gotchas that cost 3 round-trips:** `module load conda/23.11.0` must precede sourcing `conda.sh`, and it
needs `ssh … 'bash -lc'` — without a login shell there is no `module`, the activate silently fails, and you
land on **/usr/bin/python 2.7.5**, which emits SyntaxErrors that look like corrupt source. H2's system
`python3` is **3.6/ASCII**: patch scripts must pass `encoding="utf-8"` (rectify sources contain em-dashes)
and must not `py_compile` a 3.7+ codebase.

# AGENT_FIXES.md

Fast coordination log for active debugging sessions across M1 / H2 / Sherlock.
**Read this before touching pipeline code. Update it when you find a bug.**
Archive entries into CHANGELOG.md when the session wave is done.

> **[2026-07-06] TODO — PROPAGATE the split-continuation fix (9efd4f0): H2/Sherlock rectify installs are STALE (0.9.0, pre-fix).**
> `rectify split --generate-slurm --scheduler uge` in rectify **0.9.0** (the H2 env at
> `/u/project/guillom/shared/envs/rectify`) emits `run_array_correct_*.sh` with **blank lines inside the
> backslash-continued `rectify correct` command** (after `--annotation "..." \`). A blank line has no
> trailing `\`, so bash TERMINATES the command at `--annotation`, dropping `-o`, `--aligner-bams`, etc. →
> `rectify correct: error: the following arguments are required: -o/--output` on EVERY chunk → silent empty
> aggregates. This bit the **xrn1-AID `drs_decay_mutants` fan-out (2026-07-06)** — all 9 conditions produced
> zero 3'/5' output (diagnosed + patched-in-place; see Chanfreau `planning/111_xrn1_resubmit.md`).
> **The generator bug is ALREADY FIXED in code by `9efd4f0` (2026-07-03, split_command.py + a 270-line test)
> — but NOT propagated:** `drs-validation-rebuild` is **[ahead 2] of origin** (unpushed), and the H2 +
> Sherlock envs still run **0.9.0**. **TODO: (1) push `9efd4f0` to GitHub; (2) reinstall/update the rectify
> env on H2 (`/u/project/guillom/shared/envs/rectify`) AND Sherlock so future `rectify split` emits clean
> scripts.** Interim workaround: `perl -i -0pe 's/\\\n\n+/\\\n/g'` on the generated correct scripts.

> **[2026-05-29] OPS NOTICE — M1 is sluggish. Offload moderately-heavy work to a login node.**
> All agents: do NOT run moderately-heavy in-process work on M1 right now (anything
> loading multiple bedgraphs/BAMs, bootstraps/permutations, or peaking toward ~1 GB).
> Push it to an **H2/Sherlock login node** — `ssh hoffman2 'bash -lc "..."'` (no qsub,
> no queue, far more headroom than M1's crowded desktop). Reserve M1 for editing, git,
> reading, and single small scripts. Escalate to a compute allocation (sdev/qrsh/sbatch)
> only for OOM-risk or subprocess fan-out (numpy on a Sherlock *compute* node needs the
> AVX-512 constraint; H2 login is friction-free). See memory `feedback-m1-memory-discipline`.

**Perf work:** see `dev/PERF_AUDIT.md` — playbook for finding/fixing per-read
over-computation (profile-don't-guess, py-spy recipe, anti-patterns, suspect
hotspots). The 3'SS-rescue full-pool stall (entry below) is the worked example;
it is likely not the only such bottleneck.

---

## [2026-07-09] BUG: one pathological read crashes the WHOLE chunk_merge (Lazy HP scoring IndexError)

**Status:** FIXED on M1 (uncommitted, drs-validation-rebuild) + SCG rectify_scg_src. Chanfreau
`planning/180_bug3_hp_scoring_indexerror_fix.md`.

`_read_hp_edit_distances_from_raw_bam` (core/consensus/corrected_consensus.py) scores each corrected read
under `strict=True`; a single edge-case read (`IndexError: string index out of range` during
`apply_corrected_edits_to_read`/`_cigar_hp_edit_distance` — e.g. a CIGAR/coord near a contig boundary)
**re-raised and failed the entire chunk_merge task** (35k+ good reads lost), `afterok`-cascade-cancelling
final_merge→consensus→merge_consensus → sample `.done_rc`=1. Data-dependent (hit 1 of wt_rep1's 109 chunks
on the SCG panel), so it recurs across samples. **Fix (2 additive guards, no behavior change for valid
input):** (1) wrap the per-read scoring loop in try/except → log + score worst-case
`(inf, 0, _NO_JUNCTION_ANCHOR)` (read is already in `seen_corrected_ids`, so strict completeness still
passes); (2) `core/splice/hp_penalty.py::_hp_run_length` bounds-guard `pos` before `seq[pos]`.

## [2026-06-14] FEAT: strip terminal (AAG/GAA)n triplet-repeat basecaller artifact in DRS pre-trim

**Status:** DONE. Implements `dev/specs/TODO_aag_pretrim_strip_20260614.md`. Trim-level fix +
validation only (a full Sumner re-align/re-correct of all 11 samples is a separate decision,
NOT done here).

**Problem:** Dorado v5.2.0 / RNA004 mis-basecalls the poly-A homopolymer tail as a low-period
multi-base repeat `(AAG/GAA)ₙ` (read-side artifact, period-3, k=3). The DRS pre-trim's poly-A
scan finds no terminal A-run, reports `polya_len≈1`, and the artifact is carried into alignment
and force-aligned (the cat2_minus_2 11D-bridge mess).

**Fix (`rectify/core/commands/drs_trim_command.py`):** new `_find_terminal_repeat_block` peels a
terminal low-period **multi-base** repeat block BEFORE the poly-A scan, exposing any genuine
poly-A behind it. Reuses `repeat_expansion.dominant_repeat_period` (no new motif scanner). Key
correctness point the spec under-specified: **poly-A is period-consistent for *every* k**, so a
naive period walk merges the genuine poly-A run into the artifact block and destroys recovery.
Two-stage algorithm: (1) coarse extent = longest terminal phase-walk over candidate periods, to
identify the *true* period via `dominant_repeat_period`; (2) re-walk at the true period for the
precise 5' boundary, anchored on the leftmost on-phase non-A base. Off-phase basecall errors
tolerated (break at 2 consecutive, mirroring `_scan_polya`); purity gate (`min_frac=0.8`) bounds
body absorption; multi-base-motif rule = strict no-op on pure poly-A.

`find_polya_and_adapter` now returns a **6-tuple** (added `repeat_len`, `repeat_motif`); all 5
return points + 3 callers + parallel path + `_trim_unmapped_task` metadata updated. **Critical:**
the trim-decision condition is now `(polya_len >= 1 or repeat_len > 0)` — the dominant artifact
case has `polya_len=0` (the (GAA)ₙ replaced the tail) yet must still be peeled.

**CLI:** `--strip-repeat-expansion` (BooleanOptionalAction, **default ON**), `--repeat-min-len`
(default **15**), `--repeat-min-frac` (default 0.8). New metadata columns `repeat_len`/`repeat_motif`.

**Validation (2026-06-14, Sherlock):**
- **Yeast no-op PROVEN:** 0/36 Cat1–9 validation reads trigger the strip; 0 `polya_len` changes
  vs strip-off. (This decided default-ON being safe for the validation suite.)
- **CNTL_21.8 (Sumner SMA), 20k primary, measured at the shipped `min_len=15`:** terminal-GAA
  fraction on trimmed 3' ends **42.6% → 3.5%**; strip fires on 40.5% of reads, motif 99.9% `AAG`;
  14.2% of reads recover poly-A≥10 behind the artifact (median `polya_len` 1→4). Residual 3.5% =
  artifacts >150bp (adapter_window cap) + sub-15bp short blocks now intentionally spared + low-purity
  edges. (At `min_len=9` the residual was 2.3% / fired 41.7% / recovery 14.6%, p5=17 — that
  characterization run is what informed the min_len choice below; the 15-numbers above are the
  shipped config.)
- **min_len=15** chosen from the CNTL_21.8 block-length distribution (median 59bp, p5=17; only
  2.9% of fired reads in 9–14bp) — keeps ~97% of true artifacts while sparing genuine short
  GAA/AAG genomic tracts (SMA = repeat-disease context). Confirmed with Kevin.
- Unit tests: `tests/test_drs_trim_repeat_strip.py` (15 tests, M1+Sherlock green). NOTE: the spec
  pointed at `tests/test_polya_trimming.py` but that tests a *different* module
  (`rectify.core.polya.polya_trimmer`); a dedicated file was created instead.

**NOT exercised by `test_validation_reads*.py`** — that suite runs `rectify correct` on pre-built
per-aligner BAMs, never `trim-polya`, so it is orthogonal to this change. Run anyway as a
no-regression / import check: **124 passed, 85 skipped, EXIT=0** (Sherlock compute-node sbatch job
29577944; the login node throttles the `--annotation` correct fan-out → use a compute node, not
login, for this suite).

---

## [2026-06-13] FIX: walkback guard-refactor (structural) + softclip_rescue stop-base bug + CI 100%-non-A gate

**Status:** DONE. Implements `dev/specs/TODO_walkback_guard_refactor_20260613.md` (+ ADDENDUM:
target is **100% non-A**, not ≤2%). Cat1–9 e2e green on H2 (108 passed / 8 skipped); 61 M1
walkback/indel unit tests green.

**What changed (3 files):**

1. **`rectify/core/correct/walkback.py` — Option A: removed `early_exit_homopolymer_check`.**
   The genomic-A-tract proxy is invalid under the DRS 100%-non-A policy (pre-trim removes the
   tail; absence of a genomic A×4 is NOT a reason to skip the walk) and was the recurring lone-A
   suppressor (`30d2280`/`77ced6e`/`1b1db38`/`a1728eb` whack-a-mole). A "provably-safe terminal
   short-circuit" (skip when 3' end is a non-stop read==genome match) was **tried and rejected** —
   it reopened the cat1 regression, because a force-aligned-past-pA read can have a *coincidental*
   single terminal match after a pA-mismatch run that only the main scan's tail-context guard
   catches. So: **never short-circuit; always scan.** Guards 1–3 only *restrict* the anchor.

2. **`rectify/core/correct/walkback.py` — Option B: `_enforce_non_stop_anchor` post-condition.**
   When the guarded walk returns None and would leave the end on a stop base **and no real-N-op /
   large-del guard clipped the scan**, re-run the unguarded clean core (never `max_scan_depth`-
   capped) and adopt its non-stop anchor. Guard-respecting → defers to genuine artifact rejections
   and the deliberate left-side real-N-op fallback (the addendum's permitted all-A-span residual);
   also rescues the `max_scan_depth` truncation edge. No-op on Cat1–9 (anchors well within depth).

3. **`rectify/core/correct/indel_corrector.py` — removed the cat2_minus_2 "TTGC motif" 2-bp
   deletion-extension block** in `rescue_softclip_at_homopolymer` (−strand path). It hard-coded
   the 2026-05-18 directive and re-extended *past* the strip-trailing-T's guard onto genome
   chrI:128096 (a T = RNA A) — a 100%-non-A violation. Kevin (2026-06-13): the directive came from
   a **+strand** case and was mis-applied to this −strand read after the validation bundle was
   regenerated → **superseded, removed.** cat2_minus_2 now lands at 128102 (first non-A;
   genome[128102]=A → RNA T). NOTE: 128102 is the outward-side non-A boundary of the genomic
   T-run; inward/gene-body-side non-A is 128117 — final CPA boundary pending Kevin confirmation.

4. **`tests/test_validation_reads.py` — CI metric `TestCorrectedEndsAreNonA`** (the §5 gate):
   asserts **0** corrected DRS 3' ends on a gene-strand stop base, enumerating any residual with
   strand (the +strand skew is the canary). Updated the pinned cat2_minus_2 exact position
   128096→128102.

**The `=`-SEQ trap (don't repeat my detour):** `rectify/data/validation/validation_reads.bam`
stores SEQ with SAM `=` match-encoding, so `pysam.query_sequence` returns literal `=` for match
positions. Calling `walkback_drs_full` directly on those reads compares `ord('=')` to the genome →
every match looks like a mismatch and the scan can't anchor. **Do NOT byte-identity-probe walkback
via the raw bundled BAM.** The authoritative guardrail is `test_validation_reads.py` (it realigns
to real SEQ; cat1→10611 is impossible under `=`-pathology, so green proves real-SEQ operation).

**Sumner human-DRS verification (DONE 2026-06-14, job 29551915+29552532, `correct_20260614_walkbackfix/`):**
pooled on-A **7.09% → 0.0749%** (2532/3.38M reads), strand-balanced (+0.073% / −0.076%; the 105:1
+strand skew is GONE = lone-A bug dead). All 11 samples 0.046–0.265%. Matches the prior `1b1db38`
point-patch (0.08%) but achieved structurally + skew-free.

**Universal non-A post-condition (DONE 2026-06-14, Kevin approved option A; `bam_processor.py`
before `result['corrected_3prime'] = current_position`):** the 0.0749% walkbackfix residual was NOT
the lone-A bug — it was rescue/correction modules landing on a stop base by design (74% `homopolymer_rescue`
= variant-aware `rescue_polya_prefix_in_softclip`, whose internal-priming premise "a genomic-A CPA is
real" is exactly what the 100%-non-A policy overrides; 11% `indel_correction`; 7% `atract_ambiguity`).
**Fix = one chokepoint:** after all modules settle the corrected 3' end, if it is on a gene-strand
stop base, re-anchor via the **guard-respecting** `walkback_drs_full` (which defers to real introns /
large-dels — NOT a fresh genomic scan, which would re-cross them), or revert to the raw 3' end when
walkback can't anchor off the stop base. Reverts rescue CIGAR-surgery metadata (`sc_*`/`oc_*`/
`homopolymer_rescue_bases`) + retags `polya_walkback` so the output BAM CIGAR stays coherent, and
collapses ambiguity so NET-seq (NEW-065) can't reintroduce a poly-A base. Generalizes the walkback's
own `_enforce_non_stop_anchor` from "walkback returned None" to "ANY path left the end on a stop base".

**Result (job 29554774+29554900, `correct_20260614_postcond/`):** pooled Sumner on-A
**0.0749% → 0.0056%** (190/3.38M), strand-balanced (+0.0014% / −0.0096%); the homopolymer_rescue /
indel / softclip on-A classes are eliminated. The 190-read residual is almost entirely
`atract_ambiguity` ± `polya_walkback` — the addendum's permitted "no inward non-A anchor in the
artifact-permitted range" class (≈0). Cat1–9 e2e green (108 passed/8 skipped) — the post-condition is
a no-op there (no Cat read sits on a stop base), so its efficacy is verified on real human DRS, not in
the bundled CI; `TestCorrectedEndsAreNonA` remains the regression gate.

**Still open — cat2_minus_2 ACCURACY (not policy):** lands at 128102 (non-A, satisfies policy + CI),
but that's the OUTWARD side of the poly-A tail; the true CPA is 128117 (toward the SYN8 body). The
universal post-condition does NOT move it (128102 is already non-A). Fixing it needs a `softclip_rescue`
read-evidence change (don't extend the 3' end across a poly-A tract) — deferred as a focused follow-up
per Kevin. PNG for eye-verification: `/tmp/cat2m2/cat2_minus_2_CPA_boundary.png`. See memory
`feedback_drs_100pct_nonA_policy`.

---

## [2026-06-12] FINDING: lone-A walkback bug is ~6× worse on human DRS; recurrence is structural

**Context:** Sumner human chr5 DRS (corrected 2026-05-26, PREDATES `1b1db38`). Measured corrected
3' ends on a gene-strand A with CORRECT 0-based coords (corrected_3prime is 0-based-inclusive —
do NOT treat as 1-based; that bug produced a phantom "77% on A / walkback over-walks 1-2 bp" that
was retracted). Real numbers (CNTL_21.8): **9.94%** of corrected ends on A (vs yeast 1.6%), signature
**+strand 105:1 / 98% no-move / 79% atract_ambiguity** — the SAME [2026-06-04] lone-A guard bug, just
heavier on A-rich human 3'UTRs + v5.2.0-basecaller reads. When `polya_walkback` DOES fire it is correct
(0.01% on A) — the residual is reads the GUARD skips, not the walk failing.

**Why it keeps recurring (the user's lament — "fixed definitively 2-3×"):** production uses
`walkback_3prime_guarded` (line 319+) with 3 accreted early-exit guards (large-deletion, intron-boundary,
`early_exit_homopolymer_check`). The clean core `walkback_3prime` (43–186) is never buggy. Each guard
suppresses a real poly-A-noise artifact but ALSO skips legitimate walkback for short/lone-A cases; fixes
are per-guard patches (`30d2280` remove Case-2, `77ced6e` bridging guard, `1b1db38` lone-A bypass,
`a1728eb` cat1 regression) → whack-a-mole. **Structural fix:** always run the core walk; guards may only
VETO positively-identified artifacts, never skip the walk. **Add a CI regression metric:** assert +strand
corrected-on-A rate ≈ biological ~2% (the strand-skew is the canary). Repro: Sumner
`apa_3prime/{onA_check,onA_decomp}.py`.

**VERIFIED 2026-06-12 (eve):** re-corrected all 11 Sumner samples at HEAD (job 29376801) → pooled
corrected-on-A **7.09% → 0.08% (88×)**, every sample ≤0.27%. So fix `1b1db38` IS sufficient on the
6×-heavier human case. The Sumner production outputs should come from the re-correction
(`correct_HEAD_20260612/`), not the May-26 pre-fix run. Recurrence-prevention (structural guard
refactor + CI on-A-strand-skew metric) still stands — this just confirms the current patch holds here.
**→ Scoped as a fresh-agent TODO: `dev/specs/TODO_walkback_guard_refactor_20260613.md`** (self-contained:
problem, line numbers, Option A guard-inversion + Option B post-condition fallback, CI metric, guardrails,
verification scripts). Pick this up to end the "fixed it for the final time" cycle.

---

## [2026-06-04] BUG: DRS +strand lone-A terminus not corrected by walkback — FIXED

**Status:** FIXED in `walkback_3prime_guarded` (early-exit bypass when 3' base == stop_base).

**Symptom:** ~1.60% of DRS reads (146,826 in wt_by4742_rep1, 9.18 M) had `corrected_3prime`
on a genomic A. cDNA corrected 3' ends are 0.000% on A (walkback correctly walks all the way
back). Strand-skewed 43:1 (+). 96% had `original_3prime == corrected_3prime` (no movement);
86% flagged `correction_applied = atract_ambiguity`. Discovered during cDNA isoform atlas
analysis: `analyses/cdna_isoform_atlas_20260604/results/DRS_walkback_bug.md`.

**Root cause:** `walkback_3prime_guarded` (DRS +strand path) has an early-exit guard at
`rectify/core/correct/walkback.py:437–459` that returns `None` when the 22bp window around
the 3' end contains no run of ≥4 A's ("AAAA").  A read ending on a **lone** genomic A (one
A flanked by non-A on both sides) fails this check → walkback never fires → position stays on
the genomic A.  The −strand path deliberately has NO such guard, which directly explains the
43:1 strand skew.

**Fix** (`walkback_3prime_guarded`, lines ~453–465): added a `_base_at_3p` bypass — skip the
AAAA-run check if the 3' base itself IS the stop base (A for +strand, T for −strand):

    _base_at_3p = chrom_seq[_raw_3p].upper() if 0 <= _raw_3p < len(chrom_seq) else ""
    if _base_at_3p != stop_base and stop_base * early_exit_min_homopolymer_len not in _hp_window:
        return None  # original guard, applied only when 3' end is NOT on A

**Tests:** `TestLoneATerminationFix` added to `tests/test_walkback_readvsref.py`; all 20
walkback tests pass, including the existing parity check vs `find_polya_boundary`.

**Residual / known limitation:** the bypass operates on `read.reference_end - 1` (the raw
aligner-placed 3' end).  If Module 2C (indel correction) in `bam_processor.py` already moved
`current_position` off the raw end (which occurs for ~4% of the buggy reads), the bypass
won't fire for those reads. Impact is minor; no code change needed now — note for future
`bam_processor` refactor (thread `current_position` into the early-exit).

**Note from DRS_walkback_bug.md:** the bug doc also flagged that Module 1 (atract ambiguity,
`bam_processor.py:618-629`) annotates `atract_ambiguity` but never updates `current_position`.
That is a separate concern: Module 2E (walkback) handles position update correctly; Module 1
is a pure annotator. The ~1.6% A-on-CPA cases are explained by Module 2E early-exiting, NOT
by Module 1 missing a position update. The strand-skew (43:1 +/-) is the decisive evidence.

---

## [2026-05-30] FOLLOW-UP (deferred, gated): bundled DRS penalty_scores.tsv lacks empirical insertion rows

**Status:** OPEN — deferred by Kevin to the penalty-table maintainers; NOT a wind-down item.

**What:** the bundled yeast **DRS** `penalty_scores.tsv` has only 5 trivial I-rows (counts 1–3)
because it was generated *before* the profiler's insertion-counting code existed (insertion path
`39458d3` landed 2026-05-17 13:43, **7 min after** the table commit `ba602c3` 13:36; DRS table
never regenerated). So HP-ED scoring on DRS uses the flat `--default-ins 1.25` fallback instead of
empirical per-HP-length insertion penalties. **Verified 2026-05-29** (Sherlock jobs 26829646 /
26834583): HP over-calls are real (~1%/HP position, peak hp6 ~1.35%) and the gap is NOT the
`--isolation-flank` filter (controlled flank10-vs-0 refuted that). Only DRS is stale — the cDNA
(33 I-rows) and qsrev (25) tables were regenerated 05-20 (`ed331e7`) after the insertion code.

**Why deferred (gated):** regenerating the DRS table SHIFTS production HP-ED winner-selection. Doing
it right needs (a) the original provenance — WT by4742 full-read **5-aligner** panel (the upf1Δ
3-aligner 500k subset used for grounding is the wrong genotype + wrong panel for a production WT
table); (b) re-staging those inputs (the `wt_by4742_rep1_chunked_20260412` dev-run was cleaned);
(c) a pre/post HP-ED winner diff on real DRS data + advisor review before swapping the bundled
table. Full context + grounding numbers: `dev/specs/SPEC_overcall_rescue_and_ed_metric_20260529.md`
(INVESTIGATION FINDINGS) and memory `feedback_bundled_drs_table_lacks_insertions`.

---

## [2026-05-27] BUG: validation cat3_minus_2 (4 tests) regressed by 5ceb243 offset early-exit — FIXED

**Status:** FIXED at `_rescue_3ss_truncation_body` (both strand offset loops). Bisected to the
breaking commit; fix gates the early-exit on the current `_ed`, not the cumulative
`_best_local_ed`. Independent of the `_RESCUE_DP_CAP` perf change.

**Symptom:** `tests/test_validation_reads.py::TestCategory3JunctionRescue` failed 4 `cat3_minus_2`
tests (chrII, − strand, 5' junction rescue): the minus-strand canonical-AC rescue landed at
intron_end **366598 with N(76)D(6)** instead of the expected **366584 with clean N(82)** (+14bp
donor mis-placement, equivalence-extension never absorbs the D).

**Bisect (2026-05-27, M1 worktree):** 8e8dc8c / bd20f9e / 25d7a30 PASS (366584); **5ceb243
FAILS (366598)**; HEAD same. So `5ceb243` ("perf(splice): offset-loop ED=0 early exit") is the
breaking commit. Also confirmed independent of `_RESCUE_DP_CAP` (cap-toggle identical; clip 25bp).

**Root cause:** 5ceb243 added `if _best_local_ed == 0: break` to the offset loop. But
`_best_local_ed` is initialized ONCE before the `for _shift` loop and accumulates across shifts.
So once *any* shift reaches ED=0, every *subsequent* shift breaks after only its first offset —
skipping the offset a later (e.g. canonical-donor) shift needs to win. The commit's safety claim
("tiebreakers are shift-only, later offsets can't change the winner") is true per-shift but the
break was keyed on the cross-shift cumulative value.

**Fix:** gate the break on `_ed == 0` (the current offset's distance), not `_best_local_ed`. This
breaks only the current shift's remaining offsets (the legitimate perf intent) and reproduces the
pre-5ceb243 (25d7a30) selection exactly, while keeping the per-shift early-exit. Verified: the 4
cat3_minus_2 tests + the cat3 class green; full `pytest -m "not slow"` clean. NOTE: this also
restores `test_plus_offset_junction_rescued` to its correct `xfail` state (5ceb243 had
accidentally flipped it to XPASS — the underlying `_off`-tracking design item it tracks is
unaffected by this fix).

---

## [2026-05-26] BUG: winnowmap2 align path broken on any from-scratch genome (TWO bugs) — #2 FIXED (d89b10b), #1 OPEN

**Status:** Bug #2 **FIXED** by human-track commit `d89b10b` (2026-05-27, `-k14`→`-k15`,
deployed to Sherlock). Bug #1 (meryl `output=` syntax) **empirically confirmed real, left OPEN
per Kevin's 2026-05-27 instruction** — see below. The winnowmap2/meryl code in `multi_aligner.py`
is owned by the concurrent human-data track (see `HANDOFF_2026-05-26_yeast_perf_aligner.md` §3).
Found during the yeast aligner-swap-in eval (Sherlock job 26231454).

**Bug #1 — meryl `output=<db>` syntax — CONFIRMED REAL, OPEN.**
`multi_aligner.py:1469` builds `[meryl, 'count', 'k=15', f'output={meryl_db}', genome]`.
Tested directly on Sherlock (2026-05-27, env `rectify` meryl): the installed meryl **rejects**
`output=<db>` with `Can't interpret 'output=...': not a meryl command, option, or recognized
input file`. It wants `output <db>` as two separate tokens. So this DOES break the FIRST
from-scratch winnowmap2 build on any genome. ("meryl built fine" in the eval refers to the
*workaround* `winnowmap2_fix.sbatch`, which used the correct `output <db>` syntax — NOT the
in-code path.) Code left unchanged per Kevin's instruction; until fixed, from-scratch winnowmap2
builds need the rep-kmers pre-built externally (workaround in `winnowmap2_fix.sbatch`). Fix when
ready: `f'output={meryl_db}'` → `'output', str(meryl_db)`.

**Bug #2 — meryl-k vs winnowmap-k mismatch — FIXED (`d89b10b`).**
Was: `multi_aligner.py:1469` runs `meryl count k=15` (→ `_repetitive_k15.txt`) but :1502 called
`winnowmap ... -k14`; winnowmap requires meryl-count-k == winnowmap-`-k`, so it failed with
`input list of k-mers and winnowmap parameter k are inconsistent` (exit -6). Fixed to `-k15`
(matches the meryl build and winnowmap's documented map-ont default).

**a549 implication (for the human track to validate):** the handoff says a549 winnowmap2
"only works via a cached DB." Given bug #2, that cached rep-kmers DB must have been
hand-built at k=14 to match the `-k14` call, OR a549 winnowmap2 silently produced no/garbage
output. Worth confirming the a549 winnowmap2 BAM is actually valid.

**Eval impact:** winnowmap2 produced NO BAM, so the yeast swap-in eval has data for 4
aligners only (minimap2, gapmm2, minisplice_mm2, deSALT). winnowmap2 conclusion is pending
these fixes. Note winnowmap2 is itself an mm2-derivative, so prior expectation (like minisplice)
is it collapses onto minimap2 on the small, low-repeat yeast genome.

---

## [2026-05-26] PERF: `_rescue_3ss_truncation_body` O(N) junction iteration on dense splice loci — FIXED

**Status:** Fixed at `splice_aware_5prime.py:_rescue_3ss_truncation_body`. Tests green (1318 passed). Deploy: commit + push + rsync + resubmit.

**Symptom:** After the `_terminal_peel_rescue` K=20 peel-path cap (entry below), py-spy on WM2/MSP/A549 jobs showed stacks STILL pinned at `rescue_3ss_truncation:1110` → `_rescue_3ss_truncation_body` via the DIRECT base-rescue path (not through `_terminal_peel_rescue`). Workers at 99% CPU with near-zero checkpoint progress on human chr5.

**Root cause:** `_rescue_3ss_truncation_body` collects `_nearby_junctions` via the proximity filter at lines ~1260–1285. On junction-dense human chr5 (SMN1/SMN2 locus), the filter can admit 50–200+ junctions even in the direct-rescue call path (not just the peel path). The inner loops then iterate all junctions × N_shifts × N_offsets × `_hp_edit_distance` calls — O(N × 40K ops) per read, with N=200+.

**Fix:** Partitioned K=25 cap immediately after `_nearby_junctions` collection (lines ~1287–1316):
- When `_n_intervals` is populated (mapPacBio leading-D path): preserve all N-op-matched junctions first (they can sit far from `align_5prime` and would be wrongly dropped by edge-distance sort), then apply the K budget to the remaining proximity-only junctions.
- When `_n_intervals` is empty (MSP, WM2, and most non-mapPacBio reads): simple sort-by-edge-distance + slice.

This is correctness-safe: N-op-matched junctions are the *evidence basis* for the mapPacBio leading-D hypothesis; dropping them would silently regress rescue rates. The partition collapses to a no-op when N ≤ 25.

---

## [2026-05-26] PERF: `_rescue_3ss_truncation_body` offset-loop early exit on ED=0 — FIXED

**Status:** Fixed at `splice_aware_5prime.py:_rescue_3ss_truncation_body` (plus and minus strand offset loops). Tests green (1318 passed, 1 xfailed). Deployed to Sherlock.

**Root cause:** Within each (junction, shift) pair, the offset loop iterates up to `_off_limit` (≤110) offsets and calls `_hp_edit_distance` for each. Once a perfect match (ED=0) is found at a given offset, the remaining offsets for that shift cannot improve the result — the tiebreaker tuple `(not _in_amb, not _donor_ok, _shift_abs)` depends only on the shift, not on the offset. Iterating further offsets wastes CPU without any possibility of changing the winner.

**Fix:** Add `if _best_local_ed == 0: break` immediately after the per-offset update block (both plus-strand ~line 1524 and minus-strand ~line 1655). This exits the offset loop early as soon as a perfect match is found at the current shift.

**Correctness argument:** The tiebreaker has three components: `_in_amb` (junction ambiguity status — fixed per junction), `_donor_ok` (junction validity — fixed per junction), and `_shift_abs` (absolute shift — fixed per shift). None depend on `_off`. A different offset at the same shift would produce the same or worse ED and the same tiebreaker, so skipping it cannot change which (junction, shift, offset) triplet wins overall. NOTE: shift-loop and junction-loop breaks on ED=0 are UNSAFE — different shifts/junctions can tie at ED=0 but differ on tiebreakers (e.g., a later junction may have `_in_amb=True` which is preferred). Only the inner offset loop is safe to break early.

---

## [2026-05-26] PERF: global `junction_pool.pkl` inflated by deSALT nosec artifact junctions — DIAGNOSED

**Status:** Diagnosed. deSALT nopool resubmit running (job 26131450). Pool rebuild strategy TBD.

**Symptom:** deSALT solo smoke jobs (26128808, 26130266) stuck with 0 region checkpoints for
25+ min in Module 2H (`_score_hp_anchored` / `_precompute_del_costs`). Stack moving but rate
is ~0 regions/hr. Used the `--junction-pool-cache` pre-built pool.

**Root cause:** `chunks/junction_pool.pkl` contains 1,685,525 junctions vs 385 annotated.
Scoring 1.68M junctions per read × reads per region in pure Python (base Python 3.12 lacks
Numba) is computationally infeasible regardless of stack movement. Intron length distribution:
median 87 bp, p99 110 KB, max 448 KB — impossible for yeast (max known intron ~1 KB).
46,724 junctions are >10 KB. The pool was almost certainly built from a prior run that included
deSALT nosec BAMs; the htslib fetch() corruption produces reads with artifact N-ops at random
genomic positions, all of which get collected into the global pool.

**Immediate fix:** Resubmit without `--junction-pool-cache` (job 26131450: `run_desalt_nopool.sh`).
Per-region discovery will use only junctions seen in local reads, expected ~1K–10K vs 1.68M.

**Longer-term options:**
1. Rebuild pool from non-deSALT aligners only (minimap2, mapPacBio, gapmm2, minisplice)
2. Filter pool by max intron length (<10 KB for yeast, <1 MB for human) and min support (≥3 reads)
3. Don't use global pool for deSALT single-aligner runs at all

**Do NOT reuse `chunks/junction_pool.pkl`** until it is rebuilt without deSALT-nosec contributions.

---

## [2026-05-26] PERF: `align_clip_to_exon` O(clip_len²) hang for deeply-intronic reads — FIXED

**Status:** Fixed at `local_aligner.py:align_clip_to_exon`. Deployed when M1→Sherlock rsync runs.

**Symptom:** winnowmap2 A549 correction (job 26124625) hung from spawn: 8 workers at 99%
CPU for 80+ min, 0 new region checkpoints. deSALT solo-seq smoke (26123608) would hit the
same code path. py-spy on worker PID 204961 (sh02-17n28) showed stack pinned for 79+ min at
`_align_left_anchored:317` ← `align_clip_to_exon` ← `_rescue_3ss_truncation_body:1939`
← `_terminal_peel_rescue` ← `rescue_3ss_truncation` ← `correct_read_3prime`.

**Root cause:** `align_clip_to_exon` dispatches to `_align_right_anchored` (+ strand) or
`_align_left_anchored` (− strand) — both are O(Q×R) pure-Python Gotoh DP where Q = clip_len
and R ≈ clip_len + max_indel. When `_rescue_3ss_truncation_body` Case 4 computes
`_intronic_seq4 = _get_intronic_query_bases(read, _clip_bd4, strand)` for a read whose 5'
end is deeply inside an intron (e.g., 70k+ bp from the exon boundary), `clip_len` can be
tens of thousands of bases → Q×R ≈ 5B cells → hours of runtime in Python.

**Fix:** 500 bp guard at the top of `align_clip_to_exon` (after the `clip_len == 0` check).
Clips longer than 500 bp return an all-M fallback CIGAR identical to the existing
`not ref_region` early-return semantics. Biological 5' clips at exon boundaries are <100 bp;
clips >500 bp are artifact/pathology (read mapped deeply inside an intron).

**Confirmed callers of `align_clip_to_exon`:**
- `splice_aware_5prime.py:1817` (Case 1/2 `_exon_cigar_str`)
- `splice_aware_5prime.py:1939` (Case 4 `_exon_cigar_str4`)  ← where the hang fires
- `junction_refiner.py:717`
All three are covered by the entry-point guard.

---

## [2026-05-26] PERF: `_terminal_peel_rescue` O(D×N) junction scan on dense splice loci — FIXED

**Status:** Fixed at `splice_aware_5prime.py:_terminal_peel_rescue`. Deploy: commit + push + rsync + resubmit.

**Symptom:** WM2 (26130147) slowed to 12-min per checkpoint batch after the `align_clip_to_exon`
fix removed the prior bottleneck. MSP (26129522) showed 57-min zero-checkpoint gap with all 8
workers at 99.7% CPU. A549 (26104022) also stuck for 37+ min. py-spy: all workers in
`_hp_edit_distance` ← `_rescue_3ss_truncation_body:1596` ← `_terminal_peel_rescue:951`.

**Root cause:** `_terminal_peel_rescue` collects `nearby_junctions` (junction pool entries within
`reach = junction_proximity_bp + peel_max_bp`). On densely-spliced human chr5 (SMN1/SMN2 locus),
there can be 50–200+ junctions within `reach` even after the proximity filter. For each peel depth
`d`, `_rescue_3ss_truncation_body` iterates the full nearby set × N_shifts × N_offsets ×
`_hp_edit_distance` calls → O(D × N × 40K ops) per read. With 100+ depths × 200 junctions × many
reads, this dominates runtime.

**Fix:** Cap `nearby_junctions` to the K=20 closest junctions (by edge distance from the read's
5' end) before passing to the body. Sort by distance, keep top-20. Recompute `max_edge_dist` from
capped set so `effective_max_peel` also tightens. Provably non-regressing: a peel aligns the 5'
boundary to a 3'SS edge; the closest edge is always the most likely geometric match. The body's
own proximity window further narrows the set anyway, so distant junctions rarely win.
Fix location: `splice_aware_5prime.py` lines ~898–926 (between `reach` computation and acceptance
baseline).

---

## [2026-05-26] PERF: A549 `_hp_edit_distance` slow (not hung) on junction-dense region

**Status:** Confirmed slow progress, not a hang. Root cause now fixed by `_terminal_peel_rescue` K-cap above.

**Symptom:** A549 workers (PIDs 59510, 59527, 59544 on sh02-10n14, job 26104022) at 99%
CPU; appeared to "stop" at 451 checkpoints. Multi-sample py-spy showed stacks MOVING
(different `_hp_edit_distance` lines and call sites across 5 samples, including a direct
`rescue_3ss_truncation:1097` call with no `_terminal_peel_rescue`). By end of diagnosis
checkpoint count had grown to 455.

**Root cause:** Dense alternative-splice region in human genome — many candidate junctions
within 110 bp of read 5' end triggers O(N_junctions × N_depths × N_shifts × N_offsets)
`_hp_edit_distance` calls. Each call is bounded to 200×200 = 40K ops (by `_HP_ED_MAX_LEN`)
but the outer multiplication (100 depths × 50 junctions × 341 shift/offset combos) is
expensive. The `effective_max_peel` cap and candidate narrowing (both on Sherlock) limit
but don't eliminate the cost on pathologically dense loci.

**Action:** Left running. A549 job will complete eventually. If future runs hit the same
dense locus as a bottleneck, consider capping `len(nearby_junctions)` in
`_terminal_peel_rescue` to the K closest junctions (K=10–20).

---

## [2026-05-25] BUG (CORRECTNESS): mapPacBio `intronlen`/`maxindel` misconfigured for mammalian RNA — emits ~no N introns

**Status:** Correctness bug in `run_map_pacbio` (multi_aligner.py ~L723). On human
chr5 ONT DRS, mapPacBio emitted only ~1,046 introns vs 200k–418k for
minimap2/uLTRA/deSALT on the same reads. Two compounding param errors, confirmed
against the installed `bbmap.sh --help`:

- `intronlen=500000` (set to `max_intron`) is **backwards**. BBMap help:
  *"intronlen=999999999 — Set to a lower number like 10 to change 'D' to 'N' in
  cigar strings for deletions of at least that length."* It is a MINIMUM D→N
  relabel threshold, not a max. At 500000, only deletions ≥500 kb become `N`;
  all human introns (1–50 kb) stay `D`. The code comment ("Max intron length;
  use >=100k for mammalian RNA-seq") is wrong. **Correct: `intronlen=10`–`20`.**
- `maxindel` is never set → BBMap default **16000**. BBMap help: *"maxindel=16000
  — Don't look for indels longer than this... Set to >=100k for RNAseq with long
  introns like mammals."* At 16 kb, mapPacBio cannot even SEARCH across most
  human introns (soft-clips/fragments instead). **Correct: `maxindel=200000`.**

NOTE: the 2026-05-24 "intronlen fix" (`intronlen=50` → `intronlen=max_intron`)
was based on a misreading of the param and made intron labeling worse, not
better. For yeast (introns mostly <1 kb, default maxindel fine) the symptom was
masked.

**Proper fix (TODO):** in `run_map_pacbio`, set `intronlen=10` (or ~20) and
`maxindel=max(200000, max_intron)`; fix the comment. Decouple from the
genome-wide `max_intron` that legitimately feeds minimap2 `-G`.

**Caveat (do not over-invest):** even correctly parameterized, BBMap/mapPacBio
is a poor fit for ONT DRS spliced alignment (Křižanović 2018: 26.8% correct
exon-junction spanning on ONT vs GMAP 87.1%; PacBio short-indel error model ≠
ONT). minimap2 `-ax splice -uf -k14` is the ONT-DRS standard; uLTRA/deSALT are
the strong splice-aware alternates. Consider whether mapPacBio belongs in an
ONT splice panel at all before spending effort re-running it.

---

## [2026-05-25] BUG (salvageable, output OK): `rectify align --mapPacBio-chunk-idx K` exits 1 on success

**Status:** Orchestration bug in chunk-idx mode. The chunk alignment **succeeds**
and writes a valid `{prefix}.mapPacBio.chunk_K_of_N.bam`, but the SLURM task
**exits 1**, so an `afterok` merge dependency never fires.

**Root cause:** In chunk-idx mode `run_map_pacbio` (multi_aligner.py ~L694)
redirects its output to the `chunk_K_of_N.bam` path and **returns** that path.
But `align_command._run_one_aligner` (align_command.py L442-450) does **not
capture** the return value — it then coord-sorts/indexes (L504-512) and validates
(L590-599) the original `output_bam` = standard `{prefix}.mapPacBio.bam`, which
chunk-idx mode never writes → "output BAM is missing or empty; treating as
failed" → "No aligners succeeded" → exit 1.

**Evidence:** Sumner SMA_GSB2713 chunk array (job 26040807, 2026-05-25): all 8
`chunk_K_of_8.bam` produced (~53 MB each) yet tasks 0-7 FAILED with exit 1.

**Workaround (used):** chunk BAMs are valid. Merge mode
(`rectify align --aligners mapPacBio --mapPacBio-chunks N` with **no**
`--chunk-idx`) is unaffected — it writes the standard path, so the post-step
succeeds. So ignore the spurious chunk-task exit-1 and run merge directly (guard
on all N chunk BAMs existing). Do NOT use `--dependency=afterok` on the chunk
array (it will never be satisfied); gate the merge on chunk-BAM presence instead.

**Proper fix (TODO, not yet done):** in `_run_one_aligner`, capture
`run_map_pacbio(...)`'s returned path and use it for coord-sort + validation
(or skip coord-sort/validation when `mapPacBio_chunk_idx is not None`). Same
pattern would let chunk tasks exit 0 and `afterok` merges work as intended.

---

## [2026-05-24] AUDIT: redundant outputs + pipeline over-computation candidates logged in PERF_AUDIT

**Status:** Static audit only; no pipeline behavior changed in this entry.

**Finding:** `dev/PERF_AUDIT.md` now has a comprehensive run-all/align/correct/
manifest-analysis pass with prioritized findings and an output-necessity ledger.
Highest-priority items:
- correct-first `run-all` currently pays for raw pre-correction consensus outputs
  that the normal corrected-consensus path does not use;
- manifest analysis defaults still emit bedgraph/genomic-distribution convenience
  outputs despite the documented manifest memory contract;
- per-aligner correction can rebuild constant per-sample junction context instead
  of reusing `prescan`-style caches;
- streaming correction hand-writes a reduced position-index schema and drops
  `count_ag_rich`.

**Next action:** use `dev/PERF_AUDIT.md` as the fix queue. Any code fix must get
focused tests plus a scale rerun/output-equivalence check before being marked done.

---

## [2026-05-24] PARTIAL FIX (OOM gone; baseline-rescue perf DEFERRED): 3'SS rescue O(reads x full-pool) on big junction pools

**Status:** Perf/scaling bug, **NOT a correctness bug** (output is correct, just
slow). Partial fix committed: the **OOM is eliminated**; the full-pool config no
longer crashes but is still **CPU-bound-slow**. The remaining 87% (baseline
rescue) is DEFERRED per Kevin's "profile then pause" call (2026-05-24). RECTIFY
is v0.9.x; manifest mode is experimental — this is a documented 0.9.x limitation,
NOT a release blocker.

**Symptom:** `rectify run-all --manifest` (DRS, 5% wt+upf1 subset, 556K-read
minimap2 BAM) hung at `Processing 126 regions across 16 workers / 0%` for 6.5h
then OOM-killed at the 64GB mem cap (Sherlock job 25846844). Align + poly-A trim
+ Module-2H refine all succeeded first; the stall is in the inline `correct`
parallel region step. Run-all builds a **16,911-junction pool** from the 5
aligner BAMs and feeds it into per-read 3'SS rescue; standalone/GFF paths fed
only ~385 annotated junctions, so this never surfaced before.

**Root cause (py-spy profile, 29,999 samples, single-thread full-pool subset =
job DIAGE 25923520):** `_hp_edit_distance` (the HP-aware edit-distance DP,
splice_aware_5prime.py:693-721) is **~86% of all CPU**. Split by caller:
- **Baseline rescue `_rescue_3ss_truncation_body`: 87.1%** of the hp-ed cost.
  Two loops scan the FULL candidate set and call `_hp_edit_distance` per
  candidate, per read: the sequence-rescue loop (**line 1409 = 10,914 samples**)
  and the Case-4 intronic-snap loop (**lines 1877/1878 = 10,856 samples**).
  On a big pool in dense alt-splice loci, many introns contain the read's 5'
  boundary -> many DP calls/read.
- Terminal peel `_terminal_peel_rescue`: only **12.9%** (already narrowed below).
- `_align_right_anchored` (local_aligner): ~5%.

The OOM (vs slowness) came from the PEEL: up to `peel_max_bp=100` depths/read,
each re-scanning the full pool + a DP; `imap` is ordered, so while an early
high-coverage region stalled, workers buffered the other 125 regions' results
-> 64GB. Narrowing the peel removed the buffering blow-up (RSS 64GB -> ~14GB).

**Fix SO FAR (committed; localized to `_terminal_peel_rescue`, baseline path
untouched -> correctness-clean, 238 focused tests + 8 peel tests green):**
1. Narrow `candidate_junctions` to the nearby subset once/read (Gate B scan),
   pass that small set to the per-depth body calls -> kills the peel's O(pool)
   per-depth scan (THIS eliminated the OOM).
2. Cap peel depth at `max_edge_dist + 15` (farthest reachable 3'SS). NOTE: in
   dense regions `max_edge_dist` approaches `reach=110`, so this cap is often a
   no-op — the narrowing did the memory work, not the cap.

**DEFERRED fix (next session — the remaining 87% CPU):** apply the SAME nearby
narrowing to the BASELINE `_rescue_3ss_truncation_body` candidate loops (lines
~1409 and ~1877/1878), OR narrow once at `rescue_3ss_truncation` entry, so
`_hp_edit_distance` is computed only for nearby candidates. Same provably-safe
argument: a 5'-terminal 3'SS rescue can only win on a junction near the read's
5' boundary. Best: narrow at the rescue entry so both baseline + peel share it.
Verify: (a) focused splice/validation tests green; (b) `pytest -m "not slow"`;
(c) re-run full-pool standalone correct (`--aligner-bams x5`) on the refined BAM
at `$SCRATCH/rectify_wt_by4742_rep1_25846844_0/...minimap2.junction_refined_*.bam`
— must COMPLETE (DIAGA GFF-only baseline = 7:56) AND peak <32GB; (d) rescue
counts within DIAGA ballpark.

**Tooling notes:** py-spy installed in Sherlock rectify env (use
`py-spy record -f raw --pid <PID>` then aggregate by leaf frame). gdb `py-bt` is
UNUSABLE on this conda python3.9 (DWARF v5 vs gdb's v2-4) — use py-spy.
Diagnostic jobs + subset BAM live in `$SCRATCH/rectify_runall_diag_20260524/`.

---

## [2026-05-24] FIXED: standardize_chrom_name mangled non-yeast chroms (chr5 -> chrV)

**Status:** FIXED at this commit. Validated on human DRS (Sumner SMA chr5).

**Bug:** `standardize_chrom_name()` (rectify/utils/genome.py) unconditionally
maps `chr1`..`chr16` to yeast roman numerals, so human `chr5` became `chrV`.
Because gene interval trees are keyed with `normalize_chroms=False` (`chr5`) but
per-read lookup standardizes to `chrV`, **every** gene attribution missed ->
empty `gene_id` across all reads. The corrected BAM/TSV chrom column was also
written as `chrV`. (Latent: `chrV` is a real yeast chrom in CHROM_SIZES, so
`clamp_position`/`validate_coordinates` silently used the wrong size; on this
path clamp is a no-op for human via the `inf` default and validate is unused.)

**Fix:** module-level `_KNOWN_CONTIGS` registry in genome.py, populated by
`load_genome()` (both FASTA + pickle paths) and by
`register_genome_contigs_from_fasta()` called early in `correct_command` before
annotation/junction loading. `standardize_chrom_name` returns a registered
contig verbatim (checked AFTER CHROM_SIZES/GENOME_TO_CHROM, BEFORE the
arabic->roman fallback). Empty registry => identical legacy behaviour, so yeast
`chr1`->`chrI` and NCBI->canonical are unchanged (157 targeted tests pass).

**Verified:** patched smoke on human SMN locus -> chrom column `chr5` (was
`chrV`); `gene_id` populated on 454/547 reads with correct ENSG IDs (SMN1/2,
SERF1A). Empty `gene_id` only on intergenic/antisense reads.

**Still open (deferred, separate issue):** rectify's `load_annotated_junctions`
needs explicit `intron` feature rows (SGD convention); gencode GTFs have none,
so annotated junctions load 0 and native 3'SS rescue is disabled on human data.
Workaround for human runs: preprocess the GTF to add intron rows derived from
per-transcript exons. Only affects native Cat3 rescue / junction-guided scoring.

---

## [2026-05-21] BAM-first 10k correct smoke + RN sidecar check

**Status:** FIRST PASS SUCCEEDED on H2 for TSV correction; follow-up 2H-enabled
pass needed because the first command omitted `--aligner-bams`.

**Why this was run:** the old chunk directories often have BAMs but zero-byte
`corrected_reads.tsv` files, so smoke-test candidates must be identified from
the BAMs themselves rather than from successful `rectify correct` outputs.

**Input panel:** H2
`/u/project/guillom/shared/processed/alignments/wt_tfiiib_rep2/chunks/aligner_chunks`,
chunk `chunk_000`, all five aligners.

**BAM-derived candidate scan highlights:**
- Existing per-aligner chunk BAMs were present for `mapPacBio`, `minimap2`,
  `gapmm2`, `uLTRA`, and `deSALT`; their `corrected_reads.tsv` files were all
  0 bytes.
- The scan selected 10,000 normalized read IDs from BAM CIGAR/read-pattern
  triggers, including junction/N-op reads, terminal clips, terminal
  clip+junction combinations, short terminal exons, and terminal indel-rich
  reads.
- Reads shared by all five aligners in the source BAM panel: 346,618.

**Subset BAMs created:**
`/u/project/guillom/shared/processed/alignments/wt_tfiiib_rep2/chunks/correct_smoke_10k_from_bam/aligner_chunks/{aligner}/chunk_000/*.subset_10k.bam`

Primary read IDs retained:
- `mapPacBio`: 9,491 primary read IDs / 10,143 records
- `minimap2`: 8,704 primary read IDs / 10,168 records
- `gapmm2`: 8,366 primary read IDs / 8,372 records
- `uLTRA`: 8,697 primary read IDs / 8,706 records
- `deSALT`: 9,788 primary read IDs / 12,397 records

**First `rectify correct` smoke run:** H2 SGE array `13465761.1-5` succeeded
for all five aligners with `--streaming --emit-merged-tsv`, no corrected BAM
write, and no heap corruption. Outputs are under:
`/u/project/guillom/shared/processed/alignments/wt_tfiiib_rep2/chunks/correct_smoke_10k_from_bam/correct_outputs/{aligner}/chunk_000/`

Timings:
- `mapPacBio`: 9,504 reads processed; BAM processing 25.2s; total 27.4s
- `minimap2`: 8,707 reads processed; BAM processing 23.7s
- `gapmm2`: 8,372 reads processed; total completed successfully
- `uLTRA`: 8,700 reads processed; total completed successfully
- `deSALT`: 9,788 reads processed; total 33.3s

**Important caveat:** this first smoke run exercised 2F 3'SS rescue and core
correction, but Module 2H junction refinement was logged as disabled:
`Junction refinement: DISABLED (pass --aligner-bams to enable)`. A follow-up
10k run with the same subset BAM panel passed as repeated `--aligner-bams`
is needed before treating this as a full 2F+2H correction smoke.

**RN sidecar test:** wrote a post-hoc BAM-panel sidecar:
`/u/project/guillom/shared/processed/alignments/wt_tfiiib_rep2/chunks/correct_smoke_10k_from_bam/wt_tfiiib_rep2_chunk000_bam_panel_10k.read_num_sidecar.parquet`

Provenance:
`/u/project/guillom/shared/processed/alignments/wt_tfiiib_rep2/chunks/correct_smoke_10k_from_bam/wt_tfiiib_rep2_chunk000_bam_panel_10k.read_num_sidecar.POSTHOC_PROVENANCE.json`

Validation:
- Sidecar rows: 10,000; missing selected reads: 0.
- Normalized QNAME lookup covered every primary record in all five subset
  BAMs: `mapPacBio` 9,504/9,504, `minimap2` 8,707/8,707, `gapmm2` 8,372/8,372,
  `uLTRA` 8,700/8,700, `deSALT` 9,788/9,788.
- Fingerprint verification is intentionally not a pass/fail criterion for this
  BAM-derived sidecar because aligners can hard-clip/trim/reorient query
  sequence and qualities. For production provenance, sidecars reconstructed
  from chunk FASTQs remain the correct fingerprint source.
- Local focused module tests passed:
  `tests/test_read_num_sidecar.py tests/test_consensus_tag_restoration.py`
  → `7 passed, 1 skipped`. H2 production env lacks `pytest`, so the cluster
  validation was the live sidecar open/lookup test above.

## [2026-05-21] Module 2H junction pool explosion — CANDIDATE LOOKUP + YEAST SIZE CAP IN WORKING TREE

**Status:** FIXED in working tree, uncommitted. Needs deployment to Sherlock/H2
RECTIFY checkouts before rerunning production chunk-correct jobs.

**Scope:** DRS/CPA `rectify correct` Module 2H junction refinement, especially
all-5-aligner per-chunk correction on Sherlock where `refine_bam_junctions`
was spending tens of minutes after loading a 53k min-support-filtered pool and
where the unfiltered all-5 pool reached ~5.0M junctions for `ysh1_rep2`.

**Root cause clarified:**
- Candidate lookup was already restricted to the same chromosome.
- There was a `--junction-search-radius` default of 5000 bp, but the active
  lookup call also used `start_radius=end_radius=--junction-max-boundary-shift`
  (default 50 bp). That means candidates were endpoint-bounded to ±50 bp, not
  a free 5 kb window.
- `_candidates_near()` still scanned the per-chromosome sorted list from the
  chromosome start for every N-op until `js > ns + radius`. With large noisy
  pools, this made the cost grow with chromosome-position × N-op count.
- Module 2H had no explicit max junction/intron-size cap in pool construction
  or candidate scoring. Yeast should use an organism-tuned cap; 10 kb is a
  conservative S. cerevisiae value.

**Fixes landed:**
- `junction_scoring.py`: `_candidates_near()` now uses `bisect_left` /
  `bisect_right` to scan only the intron-start window around the current N-op
  instead of repeatedly walking from chromosome start.
- `junction_scoring.py`: added optional `max_junction_size` filtering to
  `collect_junctions_from_bam()`, `collect_junction_counts_from_bam()`,
  `build_junction_pool()`, and `_candidates_near()`.
- `junction_scoring.py`: `build_junction_pool()` now supports
  `min_observed_support` and counts N-op support across aligner BAMs; annotated
  junctions are retained regardless of support or max-size cap.
- `junction_scoring.py`: pool scans skip secondary and supplementary records.
- `junction_refiner.py`: `max_junction_size` is propagated through
  `refine_bam_junctions()`, sequential/parallel workers, and
  `refine_read_junctions()`, so it applies both to newly built pools and to
  pre-built pools loaded from `rectify prescan`.
- `correct_command.py`: added `--junction-max-size BP`.
- `prescan_command.py`: added `--junction-min-support N` and
  `--junction-max-size BP`; both are recorded in `junction_pool.pkl` metadata.

**Recommended yeast command policy:**
- Prefer preserving all aligners, but pass `--junction-max-size 10000` to both
  `rectify prescan` and per-chunk `rectify correct`.
- Keep `--junction-min-support` available as an optional guardrail, but the
  10 kb cap plus bounded lookup may make min-support unnecessary for yeast.
- If using a pre-built pool, still pass `--junction-max-size 10000` to
  `rectify correct`; it filters oversized candidates at scoring time even if
  the pickle was built before the cap existed.

**Tests run:** focused splice/refiner suite passed with single-threaded BLAS
environment variables:
```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
pytest tests/test_junction_scoring_parallel.py tests/test_junction_refiner.py
```
Result: `41 passed, 17 skipped, 1 warning in 4.41s`.

**Caveat:** running the same pytest command without BLAS thread caps aborted
during NumPy's macOS import/runtime check before test collection. This appears
environmental and not specific to the junction-refiner patch.

**Operational note:** Sherlock test array `25600177_[0-4]` was still RUNNING at
25 min on the old code path after loading `junction_pool.min10.all5.pkl`
(53,392 total junctions). Replace/rerun the single-chunk all-5-aligner test
after syncing this patch, using `--junction-max-size 10000`.

**[2026-05-21 09:38 PDT] Cluster sync + first speed check:**
- M1 working-tree code/docs/tests were synced to Sherlock and H2 so the three
  checkouts are consistent. Large `rectify/data/` payloads were excluded.
- Backups were written before sync:
  - Sherlock: `.codex_backups/pre_m1_code_sync_20260521_092607`
  - H2: `.codex_backups/pre_m1_code_sync_20260521_092608`
- A first broad rsync briefly copied package/test/docs contents into the
  checkout root. Those accidental top-level copies were moved, not deleted:
  - Sherlock: `.codex_backups/accidental_top_level_20260521_093324`
  - H2: `.codex_backups/accidental_top_level_20260521_093330`
- Verification:
  - Sherlock production Python: `tests/test_junction_scoring_parallel.py`
    + `tests/test_junction_refiner.py` → `41 passed, 17 skipped`.
  - H2 production Python (`/u/project/guillom/shared/envs/rectify/bin/python`):
    import smoke passed and `rectify correct --help` shows
    `--junction-max-size`.
- Patched Sherlock benchmark `25601448_[0-4]` is running with the previous
  min10 all-aligner pool and `--junction-max-size 10000`.
  At 09:38 PDT, tasks were still in Module 2H:
  `mapPacBio` 12:15 elapsed, other aligners 7:35 elapsed, no outputs yet.
  That is not convincing speedup evidence. Since the 10 kb cap removes only
  1,342 / 53,392 junctions from the already-min10 pool (~2.5%), this benchmark
  mostly tests the bisect lookup. If it remains slow, the next bottleneck is
  per-candidate scoring and we should add scoring/memoization instrumentation
  rather than relying on lookup changes alone.

**[2026-05-21 11:27 PDT] Speed check result — not enough:**
- Patched benchmark `25601448_[0-4]` was stopped after proving the current
  optimization is insufficient. No final corrected TSV/BAM outputs were
  produced.
- Module 2H did complete for several aligners, but still took about an hour:
  - mapPacBio: 402,070 reads, 205,619 N-op reads, 101,183 refined;
    Module 2H timing `3915.4s` (~65.3 min), then timed out in downstream BAM
    parallel streaming before final output.
  - gapmm2: 330,634 reads, 211,002 N-op reads, 81,214 refined;
    Module 2H timing `3607.7s` (~60.1 min), then remained in downstream BAM
    streaming.
  - uLTRA: 340,584 reads, 209,726 N-op reads, 85,200 refined;
    Module 2H timing `3547.4s` (~59.1 min), then remained in downstream BAM
    streaming.
  - minimap2 was still in Module 2H at ~54 min when the benchmark was stopped.
  - deSALT log ended with a Python memory-map dump while still running; treat
    this test as failed/indeterminate.
- Interpretation: bisect lookup and a 10 kb max-size cap are necessary cleanup,
  but do not make all-5-aligner Module 2H tractable on this chunk. The dominant
  cost is now per-read/per-candidate scoring and/or too many N-op reads being
  refined. Next plan should instrument candidate counts and add a hard
  candidate cap/support-ranked candidate selection or a much narrower yeast
  junction pool before running production correction.

**[2026-05-21] Plan C implementation — profiling + candidate cap:**
- Added opt-in Module 2H profiling:
  - `rectify correct --junction-profile profile.json`
  - `--junction-profile-sample-rate N` (default 1)
- Profile JSON records aggregate timings, counts, raw/after-cap candidate
  histograms, p50/p90/p99/max candidate counts, `_score_junction` calls,
  `_score_hp_anchored` tier-1/tier-2 calls, boundary-filter counts, and final
  Module 2H stats.
- Added deterministic candidate cap:
  - `--junction-max-candidates-per-nop N`
  - ranking order: current junction first, annotated junctions next, then
    smallest boundary delta, smallest intron-length delta, coordinates.
  - Recommended first speed test value: `32`.
- Normal runs are unchanged unless either new flag is supplied.
- Focused local validation:
  ```bash
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  pytest tests/test_junction_refiner.py tests/test_junction_scoring_parallel.py -q
  ```
  Result: `43 passed, 17 skipped`.
- Synced to Sherlock and H2. Sherlock focused tests passed
  (`43 passed, 17 skipped`); H2 production Python import smoke passed.
- Sherlock Plan C benchmark submitted:
  - job array: `25624678_[0-4]`
  - sample/chunk: `ysh1_rep2`, chunk `000`
  - pool: previous `junction_pool.min10.all5.pkl`
  - flags: `--junction-max-size 10000`,
    `--junction-max-candidates-per-nop 32`,
    `--junction-profile <aligner>/chunk_000/junction_profile.json`
  - output root:
    `/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/ysh1_rep2/chunks/planc_test`

**[2026-05-21] Bypass test hooks:**
- Added `rectify correct --skip-junction-refinement` to bypass Module 2H
  N-op boundary realignment/rewrite while still allowing an optional junction
  pool to be indexed for downstream Module 2F.
- Added `rectify correct --skip-3ss-rescue` to bypass Module 2F 3'SS
  truncation rescue entirely, including annotation, pool, and read-own N-op
  rescue triggers. Annotation-dependent gene attribution remains available.
- Planned benchmark: ysh1_rep2 chunk 000 all five aligners with both flags,
  to measure baseline correction runtime when both junction-specific rescue
  layers are disabled.

**[2026-05-21] ysh1 chunk 000 TSV-only vs corrected-BAM timing estimate:**
- In the no-2F/no-2H bypass benchmark, TSV correction finished quickly for
  the four usable aligners:
  - mapPacBio: BAM processing `93.0s`
  - minimap2: BAM processing `86.3s`
  - gapmm2: BAM processing `92.7s`
  - uLTRA: BAM processing `90.6s`
- The only fully completed corrected-BAM write was minimap2:
  `Corrected BAM write: 111.6s`, `Correction total: 202.6s`.
- mapPacBio/gapmm2/uLTRA had already emitted their corrected TSVs but were
  still in corrected-BAM writing when canceled at `14:35` elapsed, implying a
  lower bound of roughly `>12.6 min` extra BAM-writing time for each of those
  aligners on this chunk.
- Operational conclusion: corrected TSV manifests are sufficient for
  `rectify analyze`, so production correction should omit corrected BAMs by
  default. Generate corrected BAMs only for targeted IGV/debug chunks until
  the corrected-BAM writer is profiled/optimized.

**[2026-05-21] Lazy corrected consensus design note:**
- Dedicated plan written at `dev/specs/lazy_corrected_consensus_plan.md`.
- Core idea: keep homopolymer-aware corrected-CIGAR scoring and still emit a
  final corrected consensus BAM, but stop materializing full per-aligner
  corrected BAMs as an intermediate. Instead, score raw BAM records in memory
  using the per-aligner corrected TSV rows, select `winning_aligner`, then
  write only the winning corrected records into `corrected_consensus.bam`.
- Advisor review confirmed the main coupling points:
  `merge_corrected_tsvs()` uses per-aligner corrected BAMs only for
  `hp_edit_distance`/`aligned_bases`; generated split scripts use corrected
  BAMs both for that scoring and as the source for `corrected_consensus.bam`;
  non-chunked `run/stages.py` also unconditionally wires corrected-BAM output.
- Implementation should first target `realign_exon_blocks()` so it realigns
  only exon blocks containing homopolymer-position `X` ops, then extract the
  corrected-read edit sequence into a shared helper used by sequential writer,
  parallel writer, lazy HP scoring, and final consensus BAM generation.
- Important legacy-run caveat added to the plan: CPA DRS and H2 mex67aa
  predate RN-tagged chunk FASTQs / `*.read_num_sidecar.parquet`. Before lazy
  raw-BAM scoring is trusted on those outputs, reconstruct post-hoc read_num
  sidecars from existing chunk FASTQs using old round-robin chunk order
  (`read_num = local_index * n_chunks + chunk_index`) and write explicit
  `POSTHOC_PROVENANCE.json`. Use sidecar-backed normalized QNAME lookup first;
  inject `RN:i` into existing BAMs only if QNAME fallback coverage is
  insufficient.

**[2026-05-21] Lazy corrected consensus implementation checkpoint:**
- `realign_exon_blocks()` now targets only exon blocks that contain an `X` op
  at a homopolymer reference position, instead of realigning every eligible
  short exon block after a read-level trigger.
- Added shared `apply_corrected_edits_to_read()` in `bam_writer.py`; sequential
  and parallel corrected-BAM writers now use the same hard-clipped correction
  edit sequence that lazy scoring/final consensus BAM writing uses.
- `merge_corrected_tsvs()` now supports `per_aligner_raw_bams=...` for lazy
  HP-aware edit-distance scoring from raw BAM + corrected TSV, with strict
  identity checks by default. The old `per_aligner_corrected_bams` path remains
  supported and takes precedence when supplied.
- Added `write_corrected_consensus_bam()` to write the final corrected consensus
  BAM directly from raw BAMs + per-aligner corrected TSVs + merged
  `winning_aligner`, without materializing full per-aligner corrected BAMs.
- Generated chunk correction scripts now omit per-aligner
  `--write-corrected-bam` by default; chunk merge prefers lazy HP scoring and
  writes one `corrected_consensus.bam`.
- Added `reconstruct_posthoc_sidecar_from_chunks()` in
  `core/chunking/sidecar.py` for legacy round-robin chunk FASTQs, emitting
  `<sample>.POSTHOC_PROVENANCE.json`.
- Focused local validation:
  - `tests/test_bam_writer_parallel_smoke.py`,
    `tests/test_corrected_consensus_tiebreaker.py`,
    `tests/test_qname_sanitizer_and_validator.py`,
    `tests/test_read_num_sidecar.py`, `tests/test_splice_junction.py`,
    `tests/test_parallel_processing.py` →
    `145 passed, 1 skipped, 1 xfailed`.
  - Full `tests/test_validation_reads.py` did not run in the local sandbox:
    subprocess `rectify correct` calls aborted at OpenMP startup with
    `OMP: Error #179: Function Can't open SHM failed`, before exercising the
    assertions.
- Synced this implementation checkpoint to Sherlock and H2. Compile smoke
  passed on both. Sherlock focused test pass
  (`tests/test_bam_writer_parallel_smoke.py tests/test_read_num_sidecar.py`)
  passed: `11 passed`.

---

## [2026-05-20] CODEX audit bug sweep — CHECKPOINT/IO/PRESCAN FIXES IN WORKING TREE

**Status:** FIXED in working tree, uncommitted.

**Scope:** Follow-up fixes from `dev/BUGS_TO_FIX.md` NEW-067 through NEW-074 and
related `CODEX_AUDIT.md` durability/performance findings.

**Fixes landed:**
- `consensus.py`: checkpoint batches now fsync/validate BAMs before atomic
  `.done` sentinels; resume trusts only a contiguous valid batch prefix.
- `consensus.py`: scratch-to-output BAM and `.bai` copies now fsync destination
  files and parent dirs.
- `single_sample.py`: DRS `samtools fastq` failures preserve stderr/stdout in the
  raised error instead of discarding context.
- `single_sample.py`: DRS restored-polya BAM sort swaps now use `Path.replace()`
  instead of `unlink()` + `rename()`.
- `multi_aligner.py`: if minimap2 was requested but failed, the pipeline raises
  instead of returning a partial/empty aligner map into consensus.
- `align_command.py`: MD-tagged `.md.bam` is indexed before replacing
  `rectified.bam`; the `.bai` is atomically swapped only after index success.
- `bam/parallel.py`: parallel correction prescans the BAM once for read stats,
  coverage-region planning, and optional variant-aware rescue, removing the
  previous stats scan + per-chromosome coverage scans + variant scan fan-out.

**Tests run:** focused suite passed:
`tests/test_consensus_checkpoint_safety.py`, `tests/test_align_command_index_commit.py`,
`tests/test_multi_aligner_failures.py`, `tests/test_run_single_sample_safety.py`,
`tests/test_parallel_processing.py`, `tests/test_corrected_consensus_tiebreaker.py`,
`tests/test_analyze.py`, `tests/test_bam_writer.py`,
`tests/test_bam_writer_parallel_smoke.py` (`93 passed`). Earlier related slices
also passed for GTF feature expansion, checkpoint resume, and BAM writer smoke.

**Caveat:** full repository pytest has not been run in this session.

---

## [2026-05-20] mapPacBio QNAME sanitizer — ONT Dorado FASTQs

**Status:** FIXED (`e8c8070`)

**Affects:** any dataset sequenced with ONT Dorado (mex67aa, wtaa, and any
future DRS/cDNA deposits from the Nanopore). Does NOT affect cDNA pipeline
FASTQs produced by `rectify trim-polya` (those already have bare UUID QNAMEs).

**Symptom 1 — before fix 1:**
```
samtools view exit code 1: [E::sam_parse1] query name too long
mapPacBio failed after ~14400s
```
Dorado embeds full run metadata in the FASTQ header
(`@<uuid> runid=... ch=... flow_cell_id=... basecall_gpu=...`).
mapPacBio.sh copies the full header verbatim into SAM QNAME (346 chars);
SAM spec limit is 254. `samtools view -bS` exits 1 on every read.

**Symptom 2 — introduced by fix 1 (`838293c`), fixed in `382fcc7`:**
```
BBMap AssertionError: Error in mpb_san.fastq, line N, sequence line is blank
BBMap terminated in an error state
```
`split(' ', 1)[0]` on a no-space header (e.g. `__mpbsplit_*` sub-reads)
retains the trailing `\n`, so appending `\n` produced a double-newline
(blank sequence line).

**Symptom 3 — introduced by fix 2 (`382fcc7`), fixed in `e8c8070`:**
```
BBMap validateQualityLength: quality string length != sequence length
This can be bypassed with the flag 'tossbrokenreads' or 'nullifybrokenquality'
```
Line-by-line `startswith('@')` mis-identified quality score lines beginning
with `@` (valid FASTQ — Phred Q31) as headers and truncated them.

**Root fix (`e8c8070`):** parse sanitized FASTQ in 4-line blocks; only the
header line (block line 1) is ever touched. Quality lines are written as-is
regardless of content.

**Safe to pull:** strict correctness improvement; any FASTQ that worked
before still works.

**Downstream symptom in `rectify correct` (pre-fix):** if you didn't run
alignment yourself but inherited mapPacBio BAMs from a pre-`e8c8070` run,
`rectify correct` workers will hang/crash in the parallel-correction phase
even with the fix in place — because the bad BAMs already exist on disk.
Telltale:
```
INFO - rectify.core.bam.parallel - X regions across 8 workers
*** Error in `python': free(): invalid next size (fast): ... ***
*** Error in `python': free(): corrupted unsorted chunks: ... ***
```
followed by silent hang (Slurm reports RUNNING until walltime kill).
pysam/htslib's C code corrupts its buffer state when fetch() hits a
malformed QNAME or quality-length mismatch in the BAM. The bug LOOKS like
it's in `rectify.core.bam.parallel` (post-`0024fa3` refactor module) but
the rot is upstream in the BAM file itself. **Fix: re-align (re-run
`run_array_mapPacBio.sh`) with `e8c8070` in place, then re-run correct.**
Set2 alignments produced 2026-05-07 through 2026-05-19 all need re-alignment.

**Faster recovery for an existing pile of tainted mapPacBio BAMs (2026-05-20):**
Re-alignment costs ~60 min per chunk × 11 samples × ~10 chunks avg = ~6+ hours
of compute. The actual malformed records in the BAMs are reads with **empty
SEQ** (column 10 = `*`) that BBmap emitted when its FASTQ input had the
double-`\n` blank-sequence-line bug. Pysam scan of `wt_rep1_trimmed.mapPacBio.bam`
on 2026-05-20 found 366,728 / 4,347,067 reads (8.4%) with `query_sequence is
None`. Zero QNAME issues, zero seq/qual length mismatches — only empty SEQ.

Recovery without re-alignment: filter the empty-SEQ reads from the **merged
mapPacBio BAM only** (the one used by `rectify correct` via `--aligner-bams`
for cross-aligner consensus lookups; per-chunk input BAMs can be left alone —
they're iterated, not random-access fetched, so empty-SEQ reads get skipped):

```bash
samtools view -h --threads 4 INPUT.bam \
  | awk 'BEGIN{OFS="\t"} /^@/{print; next} $10!="*" && length($10)>0' \
  | samtools view --threads 4 -bS - > INPUT.bam.filtering
mv INPUT.bam INPUT.bam.pre_empty_seq_filter
mv INPUT.bam.filtering INPUT.bam
samtools index --threads 4 INPUT.bam
# Then invalidate stale prescan PKLs — they were built from corrupted merged BAM:
rm -f junction_pool.pkl rescue_scan.pkl
```

After this, re-submit `run_prescan.sh` → 5 `run_array_correct_*.sh` →
`run_array_chunk_merge.sh` → `run_final_merge.sh`. **Verified 2026-05-20:**
all 11 set2 samples cleared prescan in 1-7 min after filter (prescan would
crash within seconds on the unfiltered merged BAM); correct stage started
firing without `free():` errors.

**Do NOT try to also filter the per-chunk input BAMs** — they're not
coordinate-sorted (BBmap output order), so `samtools index` fails after
filter (`[E::hts_idx_push] Unsorted positions on sequence #4`). Either
sort them first, or skip — empty-SEQ reads in the input BAM are
iterated-not-fetched, so they don't trigger the htslib corruption.

---

## [2026-05-20] `rectify correct` pysam heap corruption on Han 2023 wt_R1 — STILL OPEN

**Status:** UNDER INVESTIGATION. v0.9.2 (`cb2fe6c`) only triggers the symptom because it correctly de-duplicates consensus rows (6.7M unique reads vs the previous 12.95M-row buggy output that masked it). The crash predates v0.9.2 — was always there, only now reachable.

**Affects:** `rectify correct` on a post-consensus rectified.bam with ≥ ~4M reads on H2 (16 threads, 16-core pod_smp.q node). Sherlock chunked pipeline (~290k reads per chunk in Stage D) is unaffected; Stage F.5 (full ~6.7M-read merged consensus) is the risk point.

**Symptom (matches the mapPacBio downstream symptom above):**
```
INFO - rectify.core.bam.parallel - 107 regions across 16 workers
*** Error in `python3.11`: free(): invalid next size (fast): 0x00002b... ***
*** Error in `python3.11`: free(): invalid next size (fast): 0x00002b... ***
======= Backtrace: =========
libcalignedsegment.cpython-311-x86_64-linux-gnu.so(+0x2d4b9)
libcalignedsegment.cpython-311-x86_64-linux-gnu.so(+0x4fc67)
libcalignedsegment.cpython-311-x86_64-linux-gnu.so(+0x4fb34)
...
Processing regions:  16%|██ | 17/107 [02:05<10:37, 7.09s/it]
```
Multiple workers crash, then silent hang (SGE reports `r` until walltime kill). Same pattern as the mapPacBio symptom.

**What's been ruled out (no fix needed for these — verified clean on the H2 BAM):**
- `samtools quickcheck` → exit 0
- SEQ length != QUAL length scan over full 6.7M reads → 0 hits
- QNAME length > 50 chars → 0 hits (BBmap ran with `trd=t`, qnames are bare SRR accessions)
- The `382fcc7` / `e8c8070` FASTQ sanitizer bugs do NOT apply — BBmap reads FASTQ directly, no rectify sanitizer in the BBmap path

**Scale evidence:**
| reads | threads | result |
|---|---|---|
| 1,072 | 1 | ✓ 31s |
| 1,072 | 8 | ✓ 9s |
| 102,157 (chrI) | 16 | ✓ 3:14 |
| 2,009,670 (chrI-V) | 16 | ✓ full run 61:40 wall (parallel + BAM write); exit 0 |
| 6,708,292 (full) | 16 | CRASH at heap corruption, ~17 of 107 regions completed before workers died |

**Working hypothesis:** memory pressure or specific read pattern only present on chrVI+ at scale. Need to test chrXII (rDNA, high coverage) or chrVI-XVI subset to localize.

**Mitigation in place for Sherlock Stage F.5 (job 25510786, after F's benign exit-2 detached from afterok dep):**
- Cut `--cpus-per-task=16 → 8` (halves worker count, gives each worker more memory budget)
- Added `--checkpoint-dir $OUTDIR/correct_merged_checkpoint` (per-region `.done` sentinels + rescue_scan.pkl on /oak; resume preserves work-done on crash)
- Bumped `--time=1:00:00 → 3:00:00`
- Template + deployed script both patched.

**[2026-05-20 09:47] Initial mitigation attempt: F.5 (25510786) running, no explicit crash.**
At 20+ min elapsed on sh03-07n10: 9 worker processes at 97–99% CPU each, all in `R` state. No `free()` heap-corruption log signature.

**[2026-05-20 10:10] Mitigation FAILED — silent-hang variant of the same bug.**
At 41+ min elapsed: workers STILL pinned at 99% CPU but **zero output** in scratch dir, **zero checkpoint sentinels**, **strace shows zero syscalls** — workers stuck in tight userspace loops in corrupted memory state. Same pysam heap-corruption bug, just not crashing immediately at threads=8 — instead silent-hanging. Cancelled 25510786.

**Working alternative: per-chunk consensus correct (Stage E.5, job 25514775).**
Bypass the merged-BAM crash entirely. Stage E already produced 24 per-chunk consensus.bam files (~290k reads each — the proven-safe Stage D scale). New sbatch array runs `rectify correct` on each per-chunk consensus.bam (24 small jobs at threads=8), then concat the 24 corrected_reads.tsv files into the sample-level TSV. Same scale as Stage D's per-aligner correct calls, all 48 of which succeeded. ETA ~4 hr.

**Updated guidance:** for any sample with merged BAM ≥ ~4M reads, run `rectify correct` PER CHUNK and concatenate TSVs. Do not run on the merged BAM regardless of thread count.

**Cross-cluster usage:** RECTIFY ships `rectify split --scheduler {slurm,uge,pbs}` that emits the per-chunk array scripts natively. H2 (SGE) uses `--scheduler uge --uge-queue pod_smp.q --uge-pe shared`. Sherlock (SLURM) uses `--scheduler slurm --slurm-partition larsms`. Do not use `rectify run-all` for samples >4M reads on either cluster — use `rectify split` and submit the generated array scripts instead.

**Path forward:** if 25509811 succeeds with threads=8, ship that as the workaround config in `05b_correct_merged_consensus.sh.tmpl` and re-run H2 with the same config. If it still crashes, resume from checkpoint to make incremental progress while we localize the malformed read pattern.

**Do not pull yet** — v0.9.2 is correct; this is a downstream pysam-level issue surfaced BY the correctness fix, not caused by it. Pulling `cb2fe6c` is safe; just be aware that whole-sample `rectify correct` at ≥4M-read scale may need the `--checkpoint-dir` + lower thread workaround.

**[2026-05-20] Commit B hypothesis — tempered.** Commit B's `write_corrected_bam_parallel` PRE-PARTITIONS the input BAM into per-region BAMs (≤100k ref bp each) BEFORE region workers start — each worker opens its own small region BAM via pysam, not the full 6.7M-read merged BAM. This is meaningfully different from the `--checkpoint-dir` mitigation above (which still fetched from the full merged BAM, just with sentinel-based resume). Hypothesis: workers operating on ≤100k-bp region BAMs may not trigger the pysam C-level state corruption that the merged-BAM workers do at scale.

**But this is unproven and the prior threads=8 mitigation already silent-hung.** The structural resolution test (6.7M-read Han wt_R1 full-scale run on H2 16-core) was **DEFERRED**: H2 pod_smp.q had 25,600 waiting jobs at submission time; Sherlock ControlMaster was not available (Duo re-auth needed). The chrI-V subset (2,009,670 reads, 16 threads) from prior Sherlock runs passed (61:40 wall, exit 0), validating Commit B's architecture at the proven-safe scale. The Han wt_R1 6.7M-read test remains as Outcome A/B/C to be run in a coordinated follow-up session when queues clear.

**Current working policy stands:** for samples >4M reads, use `rectify split` per-chunk array (290k reads per chunk, Stage E.5 pattern) — do NOT run `rectify correct` on merged BAMs regardless of which cluster, thread count, or commit. Commit B's region-pre-partitioning may or may not change this; until the deferred test runs, the chunk-first policy is the recommendation.

See briefing: `dev/specs/briefings/commit_b_briefing.md` §4 for the smoke deferral context.

---

## [2026-05-20] parallel BAM writer process aborts — MITIGATED IN WORKING TREE

**Status:** MITIGATED (uncommitted Codex working tree, 2026-05-20). The
parallel BAM writer no longer uses worker processes by default. It preserves
the region-planning code path but executes region workers sequentially unless
`RECTIFY_ENABLE_PARALLEL_BAM_WRITER=1` is explicitly set.

**Affects:** `rectify.core.bam.bam_writer_parallel.write_corrected_bam_parallel()`
and tests that call it with `n_threads > 1`.

**Symptom:** intermittent fatal Python aborts in/around `multiprocessing.Pool`
or process launch while exercising `write_corrected_bam_parallel()`. In the M1
sandbox, attempts to switch to spawn/forkserver/subprocess workers exposed
related low-level failures:

- `spawn` stalled in pytest on the 36-read smoke fixture.
- `forkserver` failed with sandbox `PermissionError` binding the Unix socket.
- clean subprocess workers hit OpenMP shared-memory startup failures and could
  still emit fatal abort traces during `subprocess.Popen`.

**Adjacent import bug found:** plain `import rectify` eagerly imported
`rectify.visualize`, which imported Matplotlib and ran font-manager subprocess
checks during test collection. This violates the repo rule that optional
numpy/plotting stacks must not load before thread limits are set. The working
tree now makes visualization lazy by default; `import rectify.visualize` still
works when plotting is needed, and `RECTIFY_IMPORT_VISUALIZE=1` restores the
legacy eager probe.

**Fix/mitigation in working tree:**

- Region BAM commits are now atomic: write unsorted temp → sort to temp BAM →
  fsync → `os.replace()` final region BAM → atomic `.ok` sentinel.
- Resume only trusts a region when both the `.ok` sentinel and region BAM exist.
- Unsafe worker-process execution is opt-in via
  `RECTIFY_ENABLE_PARALLEL_BAM_WRITER=1`; default `n_threads > 1` logs a warning
  and runs region workers sequentially.

**Validation:** `pytest tests/test_bam_writer.py tests/test_bam_writer_parallel_smoke.py -q`
passes on M1 (`18 passed`). This is intentionally a correctness/safety
mitigation, not a throughput fix. Re-enable true process parallelism only after
cluster-specific smoke tests prove the local Python/pysam/OpenMP runtime can
launch workers without aborts.

---

## [2026-05-20] `rectify analyze` per-condition bedgraph: 1-bp left shift

**Status:** FIXED on M1 working tree (uncommitted at time of entry); awaiting commit + push.

**Affects:** every per-condition strand-specific bedgraph emitted by
`rectify analyze` (filenames like `<condition>_plus.bedgraph` /
`<condition>_minus.bedgraph`, header
`description="RECTIFY 3' ends (<strand> strand)"`). Hits DRS and QuantSeq REV
equally — anything routed through `rectify/core/analyze/bedgraph.py::generate_bedgraphs`.

**Symptom:** the per-position 3'-end pileup peak in a `rectify analyze`
bedgraph appears 1 bp to the LEFT of the true position. Empirical:
CST6 ysh1 peak at chrIX in 2026-05-12 wbfix Han run wrote
`chrIX 287748 287749 70.66` (IGV 1-based 287,749); the matching DRS
3'-end signal from `bedtools genomecov -3 -bg` on the corrected
strand-split DRS BAM wrote `chrIX 287749 287750 18` (IGV 287,750). Gap
1 bp, consistent.

**Root cause:** `rectify/core/analyze/bedgraph.py:99-100` wrote
```python
start = int(pos) - 1
end   = int(pos)
```
as if `pos` (a `corrected_3prime` value) were 1-based. But
`corrected_3prime` is **0-based-inclusive** everywhere it is computed —
`reference_end - 1` for `is_reverse=True` (with explicit comment
`# 0-based inclusive` at `walkback.py:142,214,471`, `indel_corrector.py:1661,2027`)
and `reference_start` for `is_reverse=False`. BED is 0-based half-open,
so a single base at 0-based `pos` is the interval `[pos, pos+1)`. The
subtract-1 spelling was a leftover from an earlier 1-based convention
that was never updated when corrected_3prime semantics settled at 0-based.

**Fix:**
```python
start = int(pos)
end   = int(pos) + 1
```
Plus three regression tests in `tests/test_analyze.py::TestBedgraphCoordinates`
covering single-base emission, strand routing, and same-position aggregation.
Tests confirmed to FAIL under the pre-fix code and PASS under the fix.

**Safe to pull:** strict correctness improvement. Existing pipeline runs
that produced left-shifted bedgraphs were never used as the source of
truth for clustering or shift analysis (those work from `corrected_reads.tsv`
position columns directly, not from the bedgraphs), so the bug only
affected downstream visualization and any analysis that cross-referenced
the per-condition bedgraph against an independently-generated track. The
He et al manuscript-anchored TRT analysis (`analyses/cross_modality_trt_20260519`)
caught it during a CST6 cross-modality diagnostic; v2 classification was
unaffected because it worked from per-position attributed counts via the
DRS bedtools-derived bedgraphs, not from the `rectify analyze` output.

**Audit completed 2026-05-20** — see `dev/audits/bedgraph_coordinate_audit_20260520.md` for the full per-file table. **Three instances of the same bug found and fixed:**

- `rectify/core/analyze/manifest.py:598` — manifest-mode equivalent of `generate_bedgraphs` (streams per-sample TSVs without holding the combined frame in memory). Same `int(pos) - 1` / `int(pos)` spelling against the same 0-based `corrected_position` / `corrected_3prime` column. Corrected to `int(pos)` / `int(pos) + 1`.
- `scripts/generate_bedgraph_from_polished.py:113` — standalone CLI utility. Same off-by-1 with an explicitly wrong inline comment claiming the position was 1-based. Caught by a parallel multi-agent audit run after my initial close; my recipe had been scoped to `rectify/core/` and `rectify/data/` and missed `scripts/`. Fixed and noted in the audit doc as a scope-expansion miss.

All other bedgraph and bigwig emitters verified to already use the correct 0-based half-open convention (`netseq/netseq_output.py:122,191`, `commands/export_command.py:88,128`, and the delegate chains through `bam/bedgraph_writers.py`, `commands/{consensus,correct,analyze}_command.py`).

**Adjacent coordinate-convention findings (not bedgraph emission)** surfaced by the same multi-agent audit — tracked as separate BUGS_TO_FIX entries NEW-077 through NEW-081:

- NEW-077 (HIGH) — `bam_processor.py:826-836` minus-strand artifact N snap can land on a SKIPPED base.
- NEW-078 (HIGH/MEDIUM) — `bam_processor.py:847` ambiguity clip leaves NET-seq window in N span.
- NEW-079 (MEDIUM) — `clustering.py:532` plus-strand `distance_to_gene_3prime` is off by 1 (uses `gene['end']` against half-open annotation).
- NEW-080 (MEDIUM) — `false_junction_filter.py:282 vs 300` minus-strand wrong-flank fetch.
- NEW-081 (LOW) — `analyze_command.py:130` + `manifest.py:217` rDNA exclusion uses `<= end` against half-open regions.

**Action items still open:**
- Commit all three bedgraph fixes (`analyze/bedgraph.py` + `analyze/manifest.py` + `scripts/generate_bedgraph_from_polished.py` + `tests/test_analyze.py`) to `drs-validation-rebuild`, push, pull on H2 + Sherlock.
- Run the fast test suite on a cluster login node (M1 is memory-constrained; pytest there was killing M1 with the multi-agent audit also running).
- Triage and fix NEW-077 / NEW-078 / NEW-080 (HIGH-severity, biology-affecting); NEW-079 / NEW-081 are easy clean-ups.

---

## [2026-05-20] Read-number sidecar + RN aux-tag hybrid

**Status:** IMPLEMENTED (working tree; uncommitted at time of entry).

**Problem solved:** cDNA FASTQ comment tags (`XA/XC/XF/XU/XR/XI/XK/XS`) were
only propagated by minimap2 `-y`; non-minimap2 consensus winners lost that
metadata. `rectify split` now assigns a per-sample integer `RN:i:<read_num>`
to each derived FASTQ header and writes `<sample>.read_num_sidecar.parquet`
with the original QNAME, full FASTQ comment, chunk id, and seq/qual MD5s.

**Architecture:** Option B from `HANDOFF.md` landed. Original QNAMEs stay in
FASTQ/BAMs for IGV and grep ergonomics. minimap2 carries RN natively via `-y`;
mapPacBio, gapmm2, deSALT, uLTRA, bbmap, and bwa build QNAME→RN maps before
comment stripping and stream-inject `RN:i` into their BAMs after the existing
QNAME validator. Consensus uses RN as the K-way merge key only when every input
BAM has RN; otherwise it falls back to the existing normalized-QNAME path.

**Restoration:** `rectify consensus --read-num-sidecar` (or autodetection beside
input BAMs) restores SAM-format tags from the sidecar FASTQ comment onto each
winning consensus record without overwriting tags already present on the read.
Old RN-less BAMs remain valid and require no migration.

**Tests:** New focused tests cover sidecar round-trip/lookups/fingerprint checks,
chunker RN emission, RN BAM injection, RN-keyed merge fallback behavior, and
consensus tag restoration. Local PyArrow is unavailable, so the 1M-row streaming
RSS test is skipped locally; deployments with PyArrow exercise true parquet
row-group writing.

## [2026-05-20] QNAME pipeline hardening — sanitizer + validator + cross-aligner audit

**Status:** FIXED (working tree; uncommitted at time of entry).

**Affects:** any RECTIFY pipeline run. The active-panel mutation
class was latent in production DRS (set1 spot-check showed 0 leakage),
but multiple structural gaps were silently active for the chimeric path
and the TSV merge.

**Three connected fixes shipped today (one session):**

1. **mapPacBio sanitizer: always per-record, handle tab**
   (`multi_aligner.py: _sanitize_mpb_fastq`).
   The previous gate `_need_san = _first.startswith('@') and ' ' in _first`
   evaluated only the first FASTQ record AND only checked space (not tab),
   so bare-first files or tab-aux FASTQs (minimap2 `-y`, cDNA-pipeline shape)
   silently bypassed the sanitizer. Always rewrites per-record now; strips
   both space and tab; 254-char cap.

2. **Extended `_normalize_bam_read_name` regex** (`consensus.py:171`).
   Now covers generic SAM aux tag suffixes (`_<2c>:[AZifHB]:`) and
   enumerated Dorado metadata keys (`_runid=`, `_ch=`, `_start_time=`,
   `_flow_cell_id=`, ...). Defense in depth for any aligner that emits
   underscore-encoded comment leaks.

3. **`validate_post_alignment_qnames` + auto-sanitize**
   (`rectify/core/align/qname_validator.py`, NEW).
   Wired into ALL aligner wrappers (minimap2, mapPacBio, gapmm2, deSALT,
   bbmap, bwa, uLTRA). Samples 200 primary BAM reads; tracks 4 mutation
   classes (whitespace, overlength, cosmetic-needs-normalize,
   non-round-tripping). Default behavior: if a recoverable mutation is
   detected, stream-rewrite the BAM with normalized QNAMEs and re-validate.
   Raises only on unrecoverable mutation. Bypass with
   `RECTIFY_NO_AUTO_QNAME_SANITIZE=1` (detect-only) or
   `RECTIFY_SKIP_POST_ALIGN_VALIDATION=1` (skip entirely).

**Three cross-aligner audit follow-ups (audit:
`scripts/diagnostics/qname_mutation_survey/EDGE_CASES.md`):**

- **Chimeric path QNAME normalize** (`consensus.py:459`) — mirrored the
  non-chimeric path's line-489 normalize. Without this, mapPacBio /
  uLTRA-templated chimeric reads carried mutated QNAMEs into all
  downstream joins.
- **TSV-side normalizers reuse canonical regex**
  (`corrected_consensus.py:48`, `:518`) — deleted `_bare_uuid`, replaced
  body of `_normalize_read_id` with vectorized form of
  `_normalize_bam_read_name`; imports `_UNDERSCORE_COMMENT_RE` so the
  two cannot drift.
- **`_restore_sequence_from_aligner_reads` strand guard**
  (`consensus.py:362`) — skip donor when `donor_read.is_reverse !=
  best_read.is_reverse` to prevent RC'd SEQ injection.

**Open follow-ups (NOT done in this session):**
- Issue #3 from EDGE_CASES.md (cDNA FASTQ tags lost for non-minimap2
  winners) — structural; needs sidecar architecture or per-wrapper `-y`
  equivalent. Held for separate design session.
- Issues #4 (cosmetic tiebreaker comment), #6 (validator sequential
  sampling) — low severity.

**Survey artifacts:**
- `scripts/diagnostics/qname_mutation_survey/results.tsv` — 4-aligner ×
  53-read adversarial synthetic + 3 real-data sanity checks
- `scripts/diagnostics/qname_mutation_survey/REPORT.md` — full audit
  narrative
- `scripts/diagnostics/qname_mutation_survey/EDGE_CASES.md` —
  cross-aligner punch list (6 entries, 3 of which fixed today)

**Tests:** 113 pass in the QNAME + consensus + restore suite (47 new
covering the sanitizer, validator, vectorized TSV normalizer, and
strand guard). Pre-existing unrelated breakage in
`test_validation_reads.py` and `test_quantseq_rev_integration.py`
(Commit B manifest-output WIP) is not addressed by this work.

**Safe to pull:** strict correctness. Per-aligner BAMs will have an
added <1s validator pass at the end; nothing else slows. Existing BAMs
on disk are unchanged unless the auto-sanitizer fires on a re-run.

---

## [2026-05-20] set2 cascade scale: merged-BAM `--aligner-bams` lookups crash at chunk scale — STILL OPEN

**Status:** UNDER INVESTIGATION. Same root cause as entry 2 (Han 2023 6.7M
heap corruption) but reproduced at set2's smaller scale, demonstrating the
problem is **not** Han 2023-specific. Empty-SEQ filter (entry 1's "Faster
recovery" section) is necessary but **not sufficient**.

**Affects:** `rectify correct` invoked with `--aligner-bams` pointing at
merged per-aligner BAMs ≥ ~4M reads, when the per-chunk input BAM is just
~280k-466k reads. This is the standard set2 wrapper invocation generated
by `rectify split --generate-slurm` — each per-chunk correct does
cross-aligner consensus lookups against the merged BAMs.

**Symptom (mixed silent-hang + explicit free()):**

Same pattern as entries 1 ("Downstream symptom") and 2:
```
INFO - rectify.core.bam.parallel -   X regions across 8 workers
*** Error in `python': free(): invalid next size (fast): 0x... ***
```
… followed by memory-map dump and either:
- Silent hang (Slurm reports RUNNING until walltime kill — most cases)
- The explicit `free():` + map dump in the .err log (minority cases)

Only ~2 of dozens of logs at first inspection showed the verbatim
`free()` line; the rest had logs frozen at the "regions across N workers"
line with no further output. Workers may or may not show CPU activity —
the silent-hang variant can keep workers in tight loops in corrupted
memory state.

**Scale at which it manifests for set2 (verified 2026-05-20):**
- Per-chunk INPUT BAM: ~280k-466k reads (well under the 4M threshold)
- `--aligner-bams` MERGED BAMs: 4-15M reads per sample (above threshold)
- Cross-aligner fetch() calls into the large merged BAMs trigger the
  same htslib state corruption that entry 2 documents on Han 2023

**What the empty-SEQ filter (entry 1) fixed vs left open:**
- ✅ Fixed: merged-BAM crashes caused by the mapPacBio QNAME bug's empty-SEQ
  reads. Pysam scan showed 8.4% (366k/4.3M) of mapPacBio merged BAM reads
  had `query_sequence is None`. Filtering them allowed PRESCAN to complete
  cleanly (it would crash on the unfiltered merged BAM within seconds).
- ❌ Still open: the broader scale-related heap corruption at merged BAM
  ≥ ~4M reads, INDEPENDENT of the empty-SEQ pattern. Filtering removes
  one trigger but not all.

**Verification path that exposed it:** post-empty-SEQ-filter cascade
launched 99 jobs (11 samples × prescan + 5 correct + chunk_merge +
final_merge). Prescan succeeded on all 11. Correct stage ran but only
~19 / 530+ chunk tasks COMPLETED — and most of those 19 were the
empty-BAM-guard fallbacks (deSALT SIGSEGV chunks that exit-0 via the
6c8f5a6 guard in seconds, never doing real correction). Real-correct
tasks at 1:24+ elapsed showed the silent-hang / explicit-crash pattern.

**Working theory (from entry 2):** the issue is htslib's per-region
fetch() against multi-million-read indexed BAMs under multi-worker fan-out.
Workers' shared file handles or memory state corrupt in a way that
manifests as `free(): invalid next size`. The empty-SEQ filter removed
~8% of mapPacBio reads (an outlier-removal that brought some samples
under threshold) but didn't change the fundamental dynamic for samples
still over 4M reads in merged BAMs (every set2 sample is).

**Practical workarounds for set2 today (in order of preference):**

1. **Single-aligner correct (no `--aligner-bams`).** Bypass cross-aligner
   lookups entirely. Each per-chunk correct uses ONLY its own input BAM
   (e.g., minimap2-only). Removes the merged-BAM fetch path that's
   crashing. 3'-end accuracy stays ~95-99% of multi-aligner consensus
   per Kevin's earlier analysis (DRS yeast data; cross-aligner consensus
   mostly improves splice junction edge cases, less the 3'-end call
   itself). Used by this session as the salvage path after the 5-aligner
   cascade collapsed.

2. **Stage E.5 per-chunk-consensus correct.** Documented in entry 2 as the
   Han 2023 workaround. Requires producing per-chunk consensus.bam files
   FIRST (Stage E), then running rectify correct on each per-chunk
   consensus BAM (~290k reads each, proven-safe scale, no merged-BAM
   lookup at all). More compute than (1) but retains multi-aligner
   consensus value.

3. **Commit B `write_corrected_bam_parallel` pre-partitioning.**
   Architectural fix in flight (entry 2 references). When merged with M1,
   correct workers operate on per-region BAM slices ≤100k bp ref-space
   rather than fetching against the full merged BAM. Should resolve at
   the code level. Smoke-test pending.

**For the next agent debugging this:**

- Don't trust Slurm RUNNING — most "running" tasks at the second hour
  are hung in userspace memory-corrupted state, not making progress.
- The signature in logs is line: `rectify.core.bam.parallel -   X regions
  across N workers` followed by either silence OR `free():` errors.
- Check log progression with `tail -1 <log>.err` — if frozen at the
  "regions" line for >10 min, the task is hung.
- `grep -l "free()" *.err` undercounts because some logs go to silent
  hang without printing the free() line at all.
- Workers may show high CPU (50-99%) while hung — they're stuck in
  userspace loops in corrupted memory, not stuck in syscall waits.
- See sssion 5 handoff in `project_status_markdowns/DRS_CPA_PROJECT_STATUS.md`
  for the empty-SEQ-filter recovery context and the full set2 cascade
  timeline.

**Do not pull this entry for guidance** — pull the actual fix (entry 2's
Commit B branch or whatever supersedes it) once it lands. This entry is a
field report, not a fix.

---

*To add an entry: symptom (exact error string), root cause (one sentence),
fix commit, safe-to-pull verdict.*

---

## [2026-05-21] Analyze gene attribution for APA/shift analysis — SHERLOCK HAN VALIDATED

**Status:** Implemented in working tree, uncommitted. Local focused tests pass
and Sherlock Han BWA-only analyze completed successfully with DRS-origin
reference attribution.

**Affects:** `rectify analyze` gene-level DESeq2 and cluster shift analysis,
especially short-read QuantSeq REV runs and any long-read run where a CPA
cluster may derive from multiple upstream genes.

**Problem:** `rectify analyze` currently annotates CPA clusters with a single
nearest same-strand annotated TES using `annotate_clusters_with_genes()`
(`-500/+100` bp default). There is no protocol/read-length fork: manifest
analysis uses this heuristic for both DRS and QSrev once it has
`chrom/strand/corrected_position`. This is incorrect for:

- long reads, where read-body/TSS evidence should drive gene attribution;
- shared CPA clusters, where one cluster can legitimately have weighted
  contributions from multiple upstream genes;
- short reads, where body overlap is weak and the best attribution source is
  an imported long-read-derived cluster/position gene map, with annotation as
  an explicit fallback rather than the primary rule.

**Implementation landed in working tree:**

- Add a first-class weighted `cluster_gene_attributions.tsv` table:
  `cluster_id, gene_id, gene_name, raw_attributed_count,
  attribution_weight, source`.
- Add attribution modes for analyze:
  - `annotation`: legacy nearest-TES behavior;
  - `body`: compute attribution from corrected TSV read-body spans
    (`alignment_start/end`) against CDS/gene features;
  - `reference`: map external per-position DRS attribution TSVs
    (`chrom, position, strand, gene_id, attributed_count`) onto the current CPA
    clusters;
  - `body-then-annotation` / `reference-then-annotation`: fill unattributed
    clusters with the legacy annotation fallback.
- Thread weighted attributions into gene-level DESeq2 and shift analysis so a
  shared CPA cluster contributes fractionally to multiple genes instead of
  being forced into one arbitrary `gene_id`.
- Fill reference/body attribution `gene_name` values from loaded annotation
  when possible, so systematic IDs still map to common names for GO/plots.
- Explicit non-default attribution modes fail the analyze run if attribution
  cannot be built; the legacy `annotation` mode remains backward compatible.
- Add `--gene-attribution-reference-window BP` for reference modes. Exact
  DRS-position-to-current-cluster overlap is still preferred; if absent, a
  DRS attributed position can map to the nearest same-strand current CPA
  cluster within BP bases. This handles small DRS/QSrev peak offsets without
  reopening the broad `-500/+100` annotation window.

**Files touched for this entry:**

- `rectify/core/analyze/cluster_gene_attribution.py` (new)
- `rectify/core/analyze/__init__.py`
- `rectify/core/analyze/deseq2.py`
- `rectify/core/analyze/manifest.py`
- `rectify/core/analyze/shift_analysis.py`
- `rectify/core/commands/analyze_command.py`
- `tests/test_analyze.py`

**Local tests run:**

```bash
pytest tests/test_analyze.py -q
```

Result after the reference-window addition: `53 passed, 2 warnings in 35.93s`.
The new tests cover weighted
gene-level aggregation, weighted shift analysis, and mapping external
per-position DRS attributions onto current CPA clusters, including the bounded
nearest-cluster reference window.

**Validation:** reran Han BWA-only analyze on Sherlock with
`--gene-attribution-mode reference-then-annotation` using staged DRS origin
attribution tables under
`analyses/cross_modality_trt_20260519/inputs/*/attribution_origin/`.
Submitted full rerun as Sherlock SLURM job `25601498` using:
`/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/han2023_bwa_fast/run_han_bwa_analyze_drs_origin_ref.sbatch`.
Output target:
`/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/han2023_bwa_fast/analyze_drs_origin_ref_20260521`.
Job `25601498` failed before analysis started because the sbatch script exited
on `git log/status` from a compute-node path that did not resolve as a git
checkout. Fix the script to make git reporting best-effort, add
`--gene-attribution-reference-window`, then resubmit.
Resubmitted as Sherlock SLURM job `25602373` with
`--gene-attribution-reference-window 25`.
Result: `COMPLETED`, exit `0:0`, elapsed `01:11:16`.

Output:
`/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/han2023_bwa_fast/analyze_drs_origin_ref_20260521`

Key output metrics:
- CPA clusters retained: 231,182
- `cluster_gene_attributions.tsv`: 214,114 rows, 170,819 clusters, 8,442 genes
- Attribution row sources: 194,778 reference; 19,336 annotation fallback
- Shared multi-gene clusters: 36,513
- Gene-level DESeq2 tested 7,238 genes per contrast
- Significant genes: Ysh1AA 4,403; Rat1AA 3,727
- Weighted shift genes analyzed: Ysh1AA 6,859; Rat1AA 6,858

**Safe-to-pull verdict:** candidate fix is validated on the Han BWA-only data
and ready for review/commit. Still inspect the dirty tree carefully and stage
only the attribution-related paths; the repo contains unrelated WIP.

---

## [2026-05-21] H2 mex67aa_rep3 DRS continuation — MERGE SUBMITTED

**Status:** H2 merge-aligners job submitted and under monitoring.

**Dataset:** `/u/project/guillom/shared/processed/mex67aa_vs_wtaa_rectify_v0.9.0`

**Observed state before submission:**
- `mex67aa_rep1` and `mex67aa_rep2` each had all five merged aligner BAMs.
- `mex67aa_rep3` had completed chunk outputs for all five aligners but had no
  `chunks/merged_bams/` outputs yet.
- `mex67aa_rep3` deSALT chunks `003`, `005`, and `006` produced tiny empty BAMs
  after deSALT exited `-11`; logs identify this as likely SIGSEGV in
  `Loop-ProcessReads` from the upstream deSALT bug and intentionally emit empty
  BAMs so those chunks proceed with 4-aligner consensus.

**Submission:**
- First `qsub` failed because H2 JSV treated propagated `LC_ALL=C.UTF-8`
  locale warnings as fatal.
- Resubmitted with locale variables cleared:
  `env -u LC_ALL -u LANG qsub .../mex67aa_rep3/chunks/run_merge_aligners.sh`
- Job: `13463288` (`mex67aa_rep3_merge_aln`)
- Initial queue state at 2026-05-21 15:21 PDT: `qw`, waiting for 16 shared
  slots on `campus2.q`.

**Next steps when merge completes:**
Run the manual UGE chain from `mex67aa_rep3/chunks/submit_pipeline.sh`:
`run_prescan.sh`, per-aligner correction arrays, chunk merge, final merge,
per-chunk consensus, then consensus-chunk merge. Confirm whether correction
uses the current safer single-aligner/per-chunk strategy before launching the
large merged-BAM cross-aligner correct stage.

---

## [2026-05-21] Post-hoc read-number sidecars — CPA SHERLOCK COMPLETE

**Status:** Sherlock CPA set2 sidecars reconstructed and validated. H2 pending
environment/data prerequisites.

**Why:** Production CPA and mex67aa chunk FASTQs predate the RN/read-number
sidecar split path. They had `chunks_manifest.json`, but no
`*.read_num_sidecar.parquet`; FASTQ headers also lacked `RN:i`.

**Method:** used the legacy round-robin inverse:
`read_num = local_read_index * n_chunks + chunk_index`.
The helper writes:
- `<sample>.read_num_sidecar.parquet`
- `<sample>.read_num_sidecar.PROVENANCE.json`
- `<sample>.POSTHOC_PROVENANCE.json`

**Sherlock CPA output root:**
`/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery`

**Validated row counts:** Parquet metadata row counts exactly match each
sample's `chunks_manifest.json` `n_reads`.

| Sample | Rows |
| --- | ---: |
| wt_rep1 | 4,347,602 |
| wt_rep2 | 4,198,935 |
| wt_rep3 | 2,897,972 |
| rna15_rep1 | 1,090,278 |
| rna15_rep2 | 1,910,111 |
| rna15_rep3 | 945,995 |
| ysh1_rep1 | 1,264,614 |
| ysh1_rep2 | 2,010,907 |
| ysh1_rep3 | 1,859,080 |
| wt_tfiiib_rep1 | 11,543,753 |
| wt_tfiiib_rep2 | 8,234,198 |
| wt_tfiiib_rep3 | 7,657,191 |

**Important operational notes:**
- The first batch process was killed after writing `wt_rep1`/`wt_rep2` because
  validation used `ReadNumSidecar.open()`, which loads all rows into memory.
  The completed sidecars were valid. The rerun used one sample per Python
  process and `pyarrow.parquet.ParquetFile(...).metadata.num_rows`.
- Do not validate multi-million-row sidecars by opening them through
  `ReadNumSidecar.open()` on login nodes.
- These are post-hoc sidecars only. Existing BAMs do **not** have `RN:i` tags.
  Consensus can still use the sidecar for tag restoration via normalized-QNAME
  fallback, but RN-keyed consensus requires either BAM RN injection or
  re-alignment from RN-tagged chunk FASTQs.

**H2 status:**
- `/u/project/guillom/shared/processed/mex67aa_vs_wtaa_rectify_v0.9.0` has
  38 chunk FASTQs and 3 `chunks_manifest.json` files, so post-hoc
  reconstruction is possible.
- Installed `pyarrow==16.1.0` into
  `/u/project/guillom/shared/envs/rectify` using the env's own `pip` with
  `--only-binary=:all:`. The first unpinned `pip install pyarrow` attempted a
  source build and failed before installing anything; the pinned wheel install
  succeeded. Verified a Parquet write/read smoke test with the env Python.
- H2 mex67aa post-hoc sidecars were reconstructed and validated by Parquet
  metadata against each sample's `chunks_manifest.json`:
  - `mex67aa_rep1`: 7,925,277 rows
  - `mex67aa_rep2`: 6,311,767 rows
  - `mex67aa_rep3`: 4,111,368 rows
- `/u/project/guillom/shared/processed/alignments/wt_tfiiib_rep*` has no
  chunk FASTQs and no `chunks_manifest.json`; it cannot be reconstructed
  in-place from that copied H2 artifact tree.

---

## [2026-05-21] Lazy Corrected Consensus — ysh1 Pilot Timings and Fast-Path Test

**Status:** implemented, synced to Sherlock for testing; H2 sync pending this
checkpoint.

**Code changes:**
- `merge_corrected_tsvs()` now accepts `lazy_scoring_workers` and can score
  aligners in parallel from raw BAMs + corrected TSVs.
- Lazy raw-BAM HP scoring skips the transient corrected-BAM edit stack when a
  TSV row has no correction surgery fields that can alter the CIGAR.
- Generated split and single-sample merge paths pass their thread counts into
  lazy HP scoring.

**Full ysh1 chunk 000 baseline, Sherlock job `25635672`:**
- Aligners: `mapPacBio`, `minimap2`, `gapmm2`, `uLTRA`; deSALT excluded.
- Lazy merge + whole-read HP scoring: `1055.0s`.
- Final corrected consensus BAM writing/sort/index: `227.4s`.
- Output:
  `/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/ysh1_rep2/chunks/bypass_2f2h_lazy_consensus_region_test/chunk_000/corrected_consensus.bam`

**10k stratified fast-path pilot, Sherlock job `25637055`:**
- `N_READS_TO_TEST=10000` builds temporary per-aligner TSV subsets before
  merge/consensus.
- Trigger coverage selected from ysh1 chunk 000:
  - softclip HP rescue: 2,127 selected / 2,790 available
  - overcall rescue: 2,087 / 2,492
  - changed 3' end: 6,426 / 41,431
  - junction reads: 6,301 / 223,577
  - no module trigger: 3,957 / 140,323
  - no available examples in this bypass chunk for 2F five-prime rescue,
    intronic-tail clip, or reanchor/upstream-trim.
- Lazy merge + HP scoring: `18.8s`.
- Final corrected consensus BAM: `14.1s`.
- Per-aligner scoring stats:
  - mapPacBio: 9,944 no-edit fast path, 0 transient-edit, 9,813 scored
  - minimap2: 5,992 no-edit, 4,000 transient-edit, 9,876 scored
  - uLTRA: 6,483 no-edit, 3,520 transient-edit, 9,857 scored
  - gapmm2: 5,542 no-edit, 3,903 transient-edit, 9,040 scored

**Important next step:** implement differential HP scoring over only the read
regions where candidate aligners disagree. Current whole-read scoring is
correct but still recomputes shared blocks that cancel out when aligners have
identical or near-identical CIGAR representations.

**2F-rich targeted panel:** RPL19B/RPL20B short-exon-1 genes are good real
examples for minimap2 2F rescue. In current TSVs these appear as systematic IDs
`YBL027W` and `YOR312C`. Scanning full minimap2 corrected TSVs found dense
RPL19B 2F examples; `rna15_rep3/chunk_001` was selected because all four stable
aligners had nonempty corrected TSVs and raw BAMs.

Sherlock job `25637632`:
- `N_READS_TO_TEST=10000`
- `MIN_READS_PER_TRIGGER=1000`
- `MIN_READS_TARGET_GENES=1000`
- `TARGET_GENE_IDS=YBL027W,YOR312C`
- selected all 134 target-gene reads and all 93 target-gene 2F reads
- selected 7,552 / 9,590 available 2F reads overall
- also included intronic-tail, reanchor/upstream-trim, junction, 3'-changed,
  softclip HP, overcall, and no-trigger examples
- lazy merge + HP scoring: `18.3s`
- final corrected consensus BAM: `7.9s`
- output:
  `/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/rna15_rep3/chunks/lazy_consensus_target_2f_test/chunk_001`

---

## [2026-05-22] Lazy corrected consensus P0–P3 bug fixes — committed in `1ab71f0`

**Status:** Fixed and committed. Full suite 1219 passed, 35 skipped.

### P0 — DELETED `_correction_requires_transient_edit` fast path (DO NOT RE-ADD)

**Symptom:** 22/36 reads (61%) showed divergent HP-ED values between the lazy
path and the reference `rectify correct` path in the audit fixture.

**Root cause:** `_correction_requires_transient_edit` returned `False` for any
read where the TSV row had no surgery fields that could alter the CIGAR. The
function then skipped `realign_exon_blocks` entirely for that read. But
`realign_exon_blocks` is NOT just a CIGAR-surgery function — it also computes
the HP-edit-distance score. Skipping it set those reads' HP-ED to 0 (the
aligner default), which silently changed winner selection in the consensus step.

**Fix:** Deleted `_correction_requires_transient_edit` in full. `realign_exon_blocks`
is now called unconditionally — that function has its own cheap internal pre-check
that returns early when no HP scoring work is required. The old fast path was
not just premature optimization; it was incorrect.

**Lesson for future agents:** If you see `_correction_requires_transient_edit`
referenced anywhere (old comments, stale branches, spec drafts), treat it as a
**deleted error**. Do not re-introduce any "skip HP scoring if no CIGAR changes"
guard in the lazy consensus path.

### P1 — `_chunk_index_from_path` now raises on unparseable filenames

**Symptom:** Corrupt or non-standard chunk filenames silently returned `None`
for the chunk index, causing downstream KeyErrors with confusing messages.

**Fix:** `_chunk_index_from_path` in `sidecar.py` now raises `ValueError` on
unparseable filenames. Callers must handle the exception explicitly.

### P2 — `_load_tsv` drops pre-existing `winning_aligner` column

**Symptom:** When a corrected TSV was written by a previous run that had already
added a `winning_aligner` column, reloading it caused a duplicate-column error
in the merge step.

**Fix:** `_load_tsv` in `corrected_consensus.py` drops any pre-existing
`winning_aligner` column before returning the DataFrame.

### P3 — Removed dead `mv corrected.bam` line; wired `--write-per-aligner-corrected-bams`

**Symptom:** A dead `mv corrected.bam` call at `split_command.py:942-946` ran
after the new lazy consensus path and tried to rename a file that no longer
existed, causing a spurious error at the end of successful split runs.

**Fix:** Removed the dead line. Also wired the `--write-per-aligner-corrected-bams`
CLI flag through `split_command.py` so it propagates to the lazy consensus path.

### Sherlock timing (job 25662145, committed in timing spec)

- Dataset: `wt_tfiiib_rep3/chunk_008`, 3 aligners, ~15k corrected reads.
- Lazy path: 495 s wall, 286 MB peak RSS.
- Full path (`rectify correct --write-corrected-bam` × 3): TIMEOUT >81 min.
- Lower-bound speedup: **>9.9×** (true speedup likely 20–30×).
