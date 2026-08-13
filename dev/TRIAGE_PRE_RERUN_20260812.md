# Pre-re-run triage — every known issue, sorted by whether it changes the re-run

**Compiled 2026-08-12 by the Rectify Agent**, from `HANDOFF.md`, `HANDOFF_REALIGNER_NEXT.md`,
the three untracked session handoffs (`HANDOFF_681/682/684_*.md`), `dev/TODO.md`,
`dev/BUGS_TO_FIX.md`, the agent inbox, and `planning/691`.

Kevin's question was **"fix all known issues and bugs before re-doing Rectify"** — so the sort key
here is *does the fix change bytes in the re-run's output*, **not** severity. An item can be a real
bug and still belong in bucket 2.

**Verification tags.** `✅VERIFIED` = I read the code this session and cite `file:line`, or ran it.
`📄REPORTED` = another session's measurement, carried forward unverified. `⛔SUPERSEDED` = already
fixed; listed so nobody re-opens it. Presenting both confidences as one list is how the
"27 % fabrication" story survived eight reviewers — hence the tags.

**Stage** tells you what a fix costs: `align` = re-align (expensive), `consensus` = re-run
`rectify consensus` off the arm BAMs already on disk (cheap, no realignment), `correct` = stage 4
only, `pool-gate` = Station C only (~34 s/sample), `analyze` = re-pack.

---

## Bucket 1 — BLOCKS the re-run (fix it, or you re-run twice)

| # | item | tag | stage | failure mode | fix cost |
|---|---|---|---|---|---|
| 1.1 | **No corrected BAM is written** unless `--write-corrected-bam` is passed. The 682 production stage 4 omitted it (`correct.command.provenance.json`: `corrected_bam = None`), discarding a 121,215-read 5′ rescue per sample. | ✅VERIFIED (`planning/691`) | `correct` | Every BAM consumer sees uncorrected alignments; rbrowse imported `multialigned.bam` and drew raw 5′ soft clips at RPL22B under a "rectify-corrected" label. | **Add the flag.** Guard already landed (`2d77575`) so a future omission is loud. |
| 1.2 | **Cohort was built on tree `277c708`**, which predates `d0e3a0f` (impossible-intron truncation), `ae69e79` (write-time contig invariant) and `0128840` (`Xn`/`Xi`). | ✅VERIFIED (`align.command.provenance.json`) | `consensus` | The delivered consensus still carries the >10 kb / past-contig junction population — 3,789/3,789 survive into corrected 3′ ends. | **Consensus-only re-run on current master.** Arm BAMs on disk; no realignment. |
| 1.3 | **Station C ignores `flagged` on the canonical branch.** `flagged` is computed at `station_c.py:379` and read only inside the non-canonical `else` at `:389`. A canonical-in-class junction with `repeat_flag` **and** `selfhom_flag` lit is admitted on support+q alone. | ✅VERIFIED (`station_c.py:379-393`) | `pool-gate` | A 111 kb LTR-to-LTR artifact got `admit_candidate` with both flags lit, into the table rbrowse consumes as reviewable. | Code fix is small. **⚠️ Remedy needs Kevin — it changes admissions** (bucket 3c). |
| 1.4 | **The module docstring contradicts its own verdict table.** `station_c.py:24-25` promises unconditionally *"A flag DEMOTES a candidate…; it never discards"*; the table at `:29-37` restricts demotion to the non-canonical track. **The code follows the table.** | ✅VERIFIED (`station_c.py:24-25` vs `:29-37`) | docs | This contradiction is *why* 1.3 survived review — a reader checking the prose concludes the guard is unconditional. | One-line docstring fix, no behaviour change. **Do this regardless of 1.3.** |
| 1.5 | **2–10 kb impossible-intron residual.** The 10 kb bound catches >10 kb only; ~15 % of the >2 kb population survives. **Cannot** be closed by lowering the bound — minimap2's clean max N-op is 5,064 bp. | 📄REPORTED (668-drs-arm, 400k reads/arm) | `consensus` | Spurious 2–10 kb junctions persist into the re-run's consensus and into Station C's `review` pile. | **Design work, not a patch** — needs a per-read minimap2-disagreement discriminator. Decide whether to re-run now and accept the band, or hold. |
| 1.6 | **`--junction-pool-max-intron-len` is unreachable in production.** Exposed only on `rectify consensus`; absent from `align` and `run-all`. | ✅VERIFIED (`consensus_command.py:113`, `:424` — the only two hits) | `align`/`run-all` plumbing | "The guard did not fail, it was never reachable." Any non-default bound cannot be applied on the production path. | Add the arg + thread it. Small. |
| 1.7 | **`correct` does not propagate `Xa`/`Xc`/`Xn`/`Xt`** into `corrected_3ends.tsv` — no aligner column at all (production header is 34 columns, none of them aligner). | ✅VERIFIED (no `Xa` handling in `bam_processor.py`; header dumped from `ski2d_rep1`) | `correct` | The concordance page cannot be built from `correct` output alone. **Blocks only if the concordance page ships off this re-run** — rbrowse can otherwise join on `read_id` against the multialigned BAM. | Add columns. 🔴 **SEPARATE columns, never a pre-computed concordance class** — the collapse is a viewer decision and must stay reversible (rbrowse's constraint, and they are right). |
| 1.8 | **Re-run recipe, not a bug:** production align used `--aligners minimap2 --junction-aligners uLTRA deSALT overhang_resolver` with **no `--require-aligners`**. | ✅VERIFIED (`align.command.provenance.json`) | `align` | A missing binary silently drops an arm; the run still exits 0. That is how a "3-aligner" run shipped with 2. | Pass `--require-aligners` (on master, opt-in, `4533de5`). |

---

## Bucket 2 — real, but does NOT block the re-run

| # | item | tag | source |
|---|---|---|---|
| 2.1 | Splice-classification terminology diverges from the literature ("alternative" for one-side-novel; "novel" narrower than standard). README + figure + code. | 📄REPORTED | `dev/TODO.md:70` |
| 2.2 | `analyze` does not emit per-cluster read-weighted **center-of-mass**; consumers recompute it by streaming the full TSV. | 📄REPORTED | `dev/TODO.md:269` |
| 2.3 | DRS validation-bundle regen (stale `cat2_plus_1` minimap2 clip). Listed as blocked by gapmm2 NEW-082 — **that bug is fixed** (pinned to `25.4.5`), so this is now unblocked. | 📄REPORTED | `dev/TODO.md:308`, `dev/BUGS_TO_FIX.md` |
| 2.4 | `rectify/data/validation/corrected_3ends.tsv` not re-emitted alongside the regenerated aligner BAMs. | 📄REPORTED | `dev/TODO.md:333` |
| 2.5 | cat7 validation-set selector ranks by **intron length** instead of cross-aligner support. Latent: lives in the down-selection script, not the shipped package; yeast unaffected. Only matters *if* the yeast validation set is regenerated (i.e. 2.3/2.4). | 📄REPORTED | `dev/TODO.md:364` |
| 2.6 | RN-metadata sidecar feature parked on `wip/sidecar-rn-before-shim-repair` — revive or leave archived. | 📄REPORTED | `dev/TODO.md:29` |
| 2.7 | mapPacBio index caching — status note only, no action outstanding. | 📄REPORTED | `dev/TODO.md:228` |
| 2.8 | **The three session handoffs and Kevin's two notes files are UNTRACKED** (`HANDOFF_681/682/684_*.md`, `TODO_KR_NOTES.txt`, `Rectify_readme_changes_KR_proposed.txt`). The sessions that wrote them are stood down; their open items exist only in files no commit protects. | ✅VERIFIED (`git status`) | this session |

---

## Bucket 3 — needs Kevin: a decision, not a fix

These change **scientific output**, and rbrowse consumes the schema live. None can be worked around.

| # | decision | evidence | recommendation |
|---|---|---|---|
| 3.a | **SUS / `I_GR` repetitiveness gate.** *An overhang cannot justify a placement at a locus whose SUS exceeds the overhang length used to score it.* Parameter-free. | `planning/690`, 110,219 Station C rows: long novel `admit_candidate` 0.875 frac SUS>32 vs 0.049 for real spliceosomal introns and 0.049 for what Station C *rightly* admits — **18× separation**. Mechanistic close: `overhang_cap = 60` but 72 % of long-admitted termini need >64 bases to be unique. | My predecessor's recommendation, and it **supersedes 3.b and 3.c**. Caveat: n=72 termini, cDNA only, DRS untested. |
| 3.b | **Length pre-gate on the verdict.** Nuclear max annotated intron is 1,002 bp; all longer ones are mitochondrial (self-splicing). | Measured cost: 99 rows, 294 reads. | Works, but needs an organism constant — 3.a does not. |
| 3.c | **Consult `flagged` on both verdict branches** (the 1.3 hole). | Verified this session. | Do it, or explicitly decide the canonical track is exempt and fix the docstring to say so. |
| 3.d | **D/N motif filter.** A cliff at exactly 40 bp shows these are deletions relabelled as introns by the length rule — **panel-wide, minimap2 included**. | `dev/DN_THRESHOLD_ARTEFACT_20260812.md`. Reporting-side cost is *structurally zero* (Station C skips annotated introns before writing rows). | 🔴 Note the trap: *"switch to minimap2-only" fixes the long-junction artefact and leaves this one untouched.* |
| 3.e | **Quarantine `ysh1_rep1`/`rep2`?** The corpus's only `pod5_skip` libraries: 178k/197k junctions vs 24k for their sibling; 14.7 %/13.9 % untailed vs 1.9 %. | 📄REPORTED (668-drs-arm) | Excluding them takes ysh1-AA to n=1. |
| 3.f | **676 uBAM deletion** — 313 GB, 12 files on project. | Data safe; `pt` tags confirmed preserved; 320 GB scratch already reclaimed. | 🔴 **Needs your sign-off on the file list. md5 is the WRONG test** (27-byte `@PG` delta on 9 files; the 3 `wtaa` are genuine shard merges). |
| 3.g | **pod5 preservation email** — drafted, recipient unconfirmed. Raw signal for 3 un-rebasecalled replicates exists ONLY on two external drives. | 📄REPORTED | Time-sensitive in the "single point of failure" sense. |
| 3.h | **683 G-arbiter patch** — parked, correct, 17 tests, **zero measured impact** (0 of 23,663 suspect clips within 10 bp of a 3′SS). | `dev/683_G_ARBITER_EXPLORED_NOT_NEEDED.md` | Land as a latent-defect guard, or drop. Either way stop carrying it. |

---

## ⛔ Closed — listed only so they are not re-opened

| item | resolution |
|---|---|
| **Station C `ambiguity_window` IndexError** (668-drs-arm's crash report) | **Already fixed on master** — the bounds guard is at `overhang_informativeness.py:392-410` with an explicit comment that it must never be the only fix. Probed this session: `ambiguity_window(g[:200], 150, 260)` → `(0, 0)`, no `IndexError`. **Their uncommitted H2 patch in `/u/scratch/k/kevinroy/tree_master_6490803` is SUPERSEDED — do not retrieve it.** |
| **`Xn` tag collision** | Fixed `0128840`. Verified: `consensus.py:817` writes `Xi` (offending bp), `:964` writes `Xn` (`n_aligners_agree`). `682_drs1m` is clean (values only 1/2/3). |
| **"The 5′ soft-clip rescue is broken"** | Not a defect. Rescues 1 bp → 57 bp, annotated **and** novel-pool, in every execution mode, verified end-to-end including full `run-all`. `planning/691`. |
| **`I`-padded exon CIGAR on acceptor overrun** | **Withdrawn.** Artefact of a `MockRead` stub whose `get_aligned_pairs()` reports no mismatches. On real records: clean `<K>M`, position exact, 0 mismatched and 0 indel bases across 36 reads at overrun 0–10. `planning/691 §6.3`. |
| **gapmm2 NEW-082 read-dropping** | Fixed — installer pinned to `25.4.5`. `dev/BUGS_TO_FIX.md` "Open: None". |
| **cDNA consensus frame-trim (679)** | Fixed `b3a8c35`; independently confirmed in 684 production (`XQ==0` 52.5 % → 12.4 %). |

---

## Scope, not a bug — the largest gap between here and "done"

`TODO_KR_NOTES.txt` (Kevin's own paragraph, untracked) defines **the final deliverable of RECTIFY**
and none of it is in the buckets above because none of it is broken — it is unbuilt:

> CPA clusters → PCA → sample-clustering heatmaps → DESeq2 at **both** the individual-CPA-cluster
> and the gene level (clusters summed per gene) → reference samples auto-detected from
> `WT`/`wt`/`control`/`ctrl`/`untreated` → summary tables (genes up/down, clusters up/down, top-20
> genes by cluster-distribution shift with browser-style plots) → GO enrichment on gene-level
> changes → **four** de-novo motif scans (100 bp upstream / 50 bp downstream × enriched / depleted
> vs background) + plots.

Kevin's note says much of this exists in `roadblocks` and needs porting for production readiness.
**This is a programme, not a bug fix** — it should be scoped separately from the re-run.
