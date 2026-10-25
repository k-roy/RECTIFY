# RECTIFY TODO

Known gaps, planned improvements, and deferred work.
Add items here rather than inline `# TODO` comments where possible.

---

## Housekeeping

### ✅ DONE 2026-06-17 — Retired the duplicate `~/repos/rectify` clone (reclaimed ~2.9 G)
**Done by:** CLEANUP_UCLA agent. **Outcome:** the dup is removed; everything unique is preserved off-M1 on
H2 `/u/project/guillom/shared/archives/repos_rectify_offload_20260617/` (md5-verified):
`repos_rectify_allrefs_20260617.bundle` (`git bundle --all` — all history incl. 56b00b5/190e82f),
`repos_rectify_dirty_untracked_20260617.tgz` (the 17 dirty + untracked, incl. the 148K `stage/`), plus
`live_rectify_precious_datadirs_20260617.tgz` (live-tree align_out_4/_h2_rebaseline/handoffs, the only
non-git non-reproducible data, <1 MB). **Did NOT push to GitHub:** 56b00b5/190e82f are already present in the
live `~/work/rectify` on `wip/sidecar-rn-before-shim-repair`, and the active `drs-validation-rebuild` is 0/0
with origin → the commits are **SUPERSEDED, not pending** (the RN-metadata feature was reworked). They were
NOT merged onto the active branch; pending Kevin's call on whether to revive the RN-metadata feature. Bundle
restores them if so: `git clone <archive>/repos_rectify_allrefs_20260617.bundle`.
_(original task below, for reference)_

### ✅ DECIDED 2026-06-25 — RN-metadata sidecar: leave archived (do not revive)
**Decision (Kevin):** leave the RN-metadata sidecar feature archived (off-branch bundle on H2 +
the `wip/sidecar-rn-before-shim-repair` local branch). It was deliberately excluded from
`drs-validation-rebuild` as superseded/reworked; no revival planned. Restorable from the bundle
if ever wanted. _Original task below._

### TODO — decide whether the RN-metadata sidecar feature is worth reviving
**Priority:** Low · **Added:** 2026-06-17 (CLEANUP_UCLA, during dup retirement) · **Owner:** next RECTIFY agent
The `feat(sidecar): add read-number RN metadata path` work (commit `56b00b5`; plus `190e82f docs(agent): Codex
shims`) was developed in the now-retired `~/repos/rectify` dup and is also parked in `~/work/rectify` on
`wip/sidecar-rn-before-shim-repair`. The active branch `drs-validation-rebuild` **deliberately excludes it**
(superseded/reworked — it's 0/0 with origin). **Action: evaluate whether the RN-metadata sidecar path is still
wanted.** If yes → cherry-pick from the local wip branch, or restore from the off-M1 bundle
`/u/project/guillom/shared/archives/repos_rectify_offload_20260617/repos_rectify_allrefs_20260617.bundle`,
onto `drs-validation-rebuild`. If no → leave it archived; no action.

### ~~Retire the duplicate `~/repos/rectify` clone (reclaim ~2.9 G on the M1)~~
**Priority:** Low (housekeeping) · **Added:** 2026-06-16 (M1-cleanup pass)
**Context:** A *second* full RECTIFY clone lives at `~/repos/rectify` (2.9 G), separate from this
active working copy (`~/work/rectify`). It holds local-only work NOT on `k-roy/RECTIFY`:
- 2 commits not in any remote: `56b00b5 feat(sidecar): add read-number RN metadata path`,
  `190e82f docs(agent): add Codex startup shims`
- 17 dirty: modified `dev/{BUGS_TO_FIX,PLOT_SKILLS,TODO}.md`,
  `dev/specs/SPEC_overcall_rescue_and_ed_metric_20260529.md`,
  `dev/validation_review/validation_read_review/cat2_softclip_findings.md`,
  `rectify/core/commands/correct_command.py`, `scripts/validation_data/render_read_alignment.py`;
  untracked `dev/_*.py` experiment scripts, `dev/specs/TODO_cdna_corrected_reads_tsv.md`,
  `scripts/validation_data/upf1d_2026_05/stage/`
**Action:**
1. Check whether `56b00b5` / `190e82f` (+ the WIP) are already present here in `~/work/rectify` —
   if so, the other clone is fully redundant.
2. Cherry-pick/merge any *unique* commits; decide the untracked `dev/_*.py` (commit vs discard —
   likely throwaway). Explicit `git add <paths>` only (never `-A`); target the active branch
   `drs-validation-rebuild` (NOT `master` — that's frozen at the 0.9.0 release).
3. Push to `k-roy/RECTIFY` (`git push origin drs-validation-rebuild`).
4. Then `rm -rf ~/repos/rectify` to reclaim ~2.9 G.

---

## Documentation

### ✅ DONE 2026-06-24/25 — Splice-classification terminology aligned with the literature
**Resolved:** README class table + prose, `docs/ARCHITECTURE.md`, and both light+dark
`splice_classification` figures now use **one-side novel / both-side novel** (output labels
`alternative`/`novel` kept, mapping + "novel = not in the annotation" caveat documented).
Internal labels intentionally unchanged (non-breaking). _Original task below._

### Splice-classification terminology — align README/docs/code with the literature
**Priority:** Medium
**Context:** The splice-junction classes are named `unspliced` / `annotated` /
`alternative` / `novel` in the README ("(b)/(c)" sections), the
`splice_classification` figure, and (where surfaced) code/output. Two problems vs
standard usage:
1. **"alternative"** for a one-side-novel junction (one annotated splice site, one not)
   is non-standard and collides with "alternative splicing" (a broader phenomenon). The
   literature term is a **novel donor/acceptor** (one-side-novel).
2. **"novel"** for a both-sides-unannotated junction is narrower than common usage: in
   the literature "novel" means **not in the annotation** generally (which covers the
   one-side-novel class too), and crucially ≠ "never observed before" — many
   junctions absent from a given GENCODE *basic* build are catalogued elsewhere
   (recount3/intropolis) or present in the *comprehensive* set.
**Done so far (2026-06-15):** the `splice_classification` figure was relabeled to
`unspliced` / `annotated` / `one-side novel` / `both-side novel` with subtitle
"'novel' = not in the annotation" (`docs/figures/generate_splice_classification_v3.py`;
light PNG/SVG regenerated).
**Done (2026-06-24):**
- Regenerated the **dark** `splice_classification_dark.{png,svg}` to match the
  relabeled light variant (one-side novel / both-side novel + "novel = not in the
  annotation" subtitle). Regenerated only this figure to avoid churning all 12 dark figs.
- Updated README (class table + step-(c) prose) and `docs/ARCHITECTURE.md` (added the
  output-label → literature-term mapping + the "novel = not annotated, not
  never-observed; basic vs comprehensive" caveat).
- **Decision (made, non-breaking):** KEEP the internal/output labels
  `alternative`/`novel` (renaming would break downstream consumers' TSV schema);
  documented the display-label mapping instead. Re-surfaced as a deferred user
  decision in §"Algorithm Semantics" / Bucket C — Kevin can override to a hard rename.

**Remaining (optional):** none required for v1.0.0.
- Cross-ref: Sumner chr5 summary §4 already uses the corrected terms + caveat.

---

## Alignment / Consensus

### ✅ DONE 2026-06-25 — Fix deSALT on Sherlock (binary already replaced + verified)
**Priority:** High → RESOLVED.
**Resolution:** The crashing binary has already been replaced everywhere with the
working bioconda build. The repo's vendored `rectify/data/bin/linux_x86_64/deSALT`,
H2's `~/.conda/envs/rectify/bin/deSALT`, and Sherlock's `~/.conda/envs/rectify` +
oak-deployed binary are now ALL **md5-identical** (`e923d866…`, `desalt-1.5.6-h577a1d6_7`).
The original SIGSEGV was a DIFFERENT (older) binary. `run_desalt()` also carries three
runtime SIGSEGV mitigations: `-f` tmp file on local fs (not NFS mmap), `-G` annotation
skip, and `LD_LIBRARY_PATH` strip.
**Empirical confirmation (Sherlock compute node, job 31178641, COMPLETED):** rectify's
`run_desalt()` on the 36 DRS validation reads → **33,567-byte BAM, 34 mapped primaries,
no SIGSEGV** (VERDICT: PASS). vs the old 274B empty fallback.
**Note:** the original `project_status_markdowns/TASK_desalt_conda_sherlock_test.md`
crash-window validation (14.9k/30k/40k reads) targeted the OLD broken binary and is moot
now that the binary is the proven-working H2 build; a full-depth run can confirm at scale
if desired but is not required to close this.

---

### Evaluate Newer Aligners as Panel Additions (Minisplice, GLASS, Winnowmap2)
**Priority:** Medium · **Status (2026-06-25): mostly DONE.**

| Candidate | Status |
| --- | --- |
| **Winnowmap2** | ✅ wrapped — `run_winnowmap2()` in `multi_aligner.py` (meryl k15 repetitive-kmer build + cache). Smoke script `dev/_smoke_winnowmap2.py`. |
| **Minisplice** | ✅ wrapped — `run_minisplice_mm2()` in `multi_aligner.py` (minimap2 + DL splice signals). |
| **GMAP** | ✅ added this branch (commit `7b32fa0`, opt-in splice-aware junction aligner) — `run_gmap()`. |
| **GLASS** (bioRxiv 2025.04) | ✅ EVALUATED → **DISQUALIFIED**. GLASS is a *post-alignment BAM filter* (graph-ML "Read-AS Map" that removes falsely-spliced reads from minimap2 output), **NOT a standalone aligner** — and its **code is not publicly released**, so it cannot be installed or wrapped. See `docs/aligners/EVALUATED_AND_DISQUALIFIED.md`. (Conceptually adjacent to RECTIFY's own anchor gate; revisit as a *post-consensus* false-splice filter only IF the code is ever released.) |

**Status: COMPLETE for v1.0.0.** All three candidates resolved — Winnowmap2 + Minisplice
wrapped, GMAP added, GLASS disqualified (no code / not an aligner). No open aligner-eval work.

---

### ~~`merge_corrected_tsvs` — Winner-Selection Tie-Breaker Bug — Fixed, verified 2026-05-23~~
**Status:** RESOLVED. `_n_agree` now groups on `(read_id, chrom, corrected_3prime)` (`corrected_consensus.py:1207-1220`), with an explicit comment about the paralog same-chrom case, so a 3-aligner same-chrom position consensus beats a single-aligner outlier on a different chromosome. Historical detail below.
**File:** `rectify/core/consensus/corrected_consensus.py` — `merge_corrected_tsvs()`
**Priority:** Medium (correctness; only affects paralog-ambiguous reads)

When a read's `hp_edit_distance` is N/A across all aligners (e.g. when the
per-aligner correction wasn't given a `--junction-penalty-table`), the sort
falls through to `_chimera_ok` / `_n_agree` / `_span` / `_n_junc` — but the
fallback can still let a single-aligner outlier win over a multi-aligner
position consensus.

**Observed:** cat1_plus_1 (read 0cb5a111) is a paralog read whose sequence
matches both chrXIV:10610 and chrVI:8703 nearly identically. minimap2,
gapmm2, uLTRA agree on chrXIV:10610. deSALT picks chrVI:8703. With no
penalty table active, the merge selects deSALT — despite a 3-aligner
position consensus on chrXIV.

The user has flagged repeatedly that **mapq is NOT a cross-aligner
metric**, so any sort key derived from mapq is fragile. The replacement
should be:

1. When `hp_edit_distance` is unavailable, weight by `_n_agree` on the
   exact `(chrom, corrected_3prime)` tuple — currently `_n_agree` only
   groups on `corrected_3prime` (the integer position), so paralog
   alignments at the same numeric position on different chromosomes
   spuriously combine, and different positions on the same chromosome
   stay separate.
2. Make `_n_agree` require chrom match. A 3-aligner chrXIV consensus
   should beat a 1-aligner chrVI alignment regardless of any per-aligner
   self-reported confidence.

**Repro:** Run `python -m rectify.cli correct
rectify/data/validation/validation_reads.bam --Scer --aligner-bams ...`
on the bundle's per-aligner BAMs WITHOUT `--junction-penalty-table`.
cat1_plus_1's `corrected_3prime` lands at chrVI:8703 instead of
chrXIV:10610.

**Note:** With the bundled penalty tables (now auto-resolved from
`--Scer`), `hp_edit_distance` should be populated and this fallback path
shouldn't be reached on the validation set. But the fallback is still
the wrong shape and will bite users who omit the table or run on
non-S. cerevisiae data.

---

### ~~5' Soft-Clip Rescue — Sequence-Based Matching — Fixed, verified 2026-05-23~~
**Status:** RESOLVED. Sequence-based rescue is implemented in `rectify/core/splice/splice_aware_5prime.py` (`rescue_3ss_truncation` / `_rescue_3ss_truncation_body`): clipped 5' query bases are compared by HP-aware edit distance against the upstream-exon reference (`_hp_edit_distance(_rseq, exon_seq)`, lines ~1310/1324), and `consensus._rescue_5prime_softclip` rewards an explained clip instead of the old blind `* 2` length penalty. The more ambitious multi-hypothesis terminal-peel refinement (`dev/specs/terminal_junction_peel_2f_plan.md`) remains an opt-in, conditional FUTURE enhancement — not a v1.0.0 blocker. Historical detail below.
**File:** `rectify/core/consensus/consensus.py` — `score_alignment()`
**Priority:** High

The current 5' rescue is a blind soft-clip length penalty (`score -= five_prime_softclip * 2`).
This causes aligners that *missed* a real intron (and therefore start after it with 0 soft clip)
to beat aligners that *found* the intron but have unrescued 5' bases.

Observed example: read SRR32518284.448567 at SNC1 (chrI:87,446-87,878).
- minimap2 + gapmm2: found the intron, reference start = 87,329, 25bp 5' soft clip → score penalised
- mapPacBio: missed the intron, reference start = 87,445, 0bp soft clip → wins incorrectly

**Intended implementation:**
1. When a read has a 5' soft clip, take those clipped bases (from `read.query_sequence`)
2. Find annotated junctions within range of the read's 5' alignment start (using the
   junction BED already loaded during consensus)
3. Fetch the genomic sequence at the end of the upstream exon (same length as the clip)
4. Compute edit distance between soft-clipped query bases and exon-end reference sequence
5. If edit distance ≤ threshold (~20% of clip length), the clip is "explained" — reward
   the aligner rather than penalising it

`edit_distance()` already exists in `rectify/core/spikein_filter.py` and can be imported.

---

## Statistics / Observability

### ✅ DONE 2026-06-24 — Tie-Break Rate + Aligner Combo Breakdown surfaced
**File:** `rectify/core/consensus/consensus.py` (tracking), `rectify/core/bam/processing_stats.py`
(`write_consensus_stats_tsv`), `rectify/core/analyze/summary.py` (`generate_consensus_html_report`).
**Status:** `tied_score` already in the stats TSV; `by_aligner_combo` (the panel
available/compared per read, frozenset → count) now surfaced as `aligner_combo_<panel>`
rows in `*.consensus_aligner_stats.tsv` AND as an "Aligner Combination Breakdown"
table in the consensus HTML report. Both degrade gracefully when the key is absent.
Tests: `tests/test_consensus_stats_surfacing.py`. (Subsumes the "Expose Aligner Stats
in HTML Report" item under Analysis/Downstream.)

---

## Performance

### ✅ DONE 2026-09-05 — Resolver numba kernel default-on + loud checkpoint-schema guard (a28714a, on master 269ebe7)

`overhang_informativeness._load_numba_kernel()` imports numba lazily on the first large DP;
`RECTIFY_HP_ED_NUMBA=0` opts out. The correct-stage kernel in `splice_aware_5prime.py` stays
opt-in (spawn workers; documented OOM). Why: ISSUE-010 lifted the yeast-calibrated candidate cap on
human windows and the pure-Python DP is then 26× slower (2,259 s vs ~135 s per 1,000 reads). Tester
(145k, 8 shards): pure-Python 21,542 s vs numba 10,621 s, XJ 317 both, bit-identical counts.
`_rebuild_output_from_region_files` now refuses (naming the dir) a `correct` checkpoint whose row
width ≠ the current corrected-TSV header instead of silently dropping every row (TSV 39 → 44 cols).
**Operational:** discard pre-D3 `correct` checkpoint dirs before resuming.

### ✅ LANDED 2026-09-05 — AT-AC intron class default-on for all organisms (todo/integration 56abf3a → `origin/master` 6485226)

2H side f307311 (`_canonical_tier` treats the paired AT-AC class as canonical at GC-AG rank,
`_ATAC_DINUCS` is the single source of truth); resolver side d20ddd7 (`ResolverConfig.atac=True`,
`--no-resolver-atac` on `align` + `run-all`, `--resolver-atac` kept as a no-op; the AT-AC pass has
its OWN candidate budget — the shared budget was starving the canonical pass). Price: yeast ×1.15
on the validation suite; **human ×2.93 wall on holdout1k (candidates ×1.67)** — Kevin's decision;
surface the human number when the tester times the next tip. The `test_validation_reads.py`
148 → 593 s scare was machine load (BUGS_TO_FIX A15, CLOSED). Still open: A14 `run_multi_aligner`
drops the resolver knobs (dead code path today).

### 🔄 IN PROGRESS 2026-09-05 — `correct` per-module timing counter + resolver-verdict tags → 2F skip (Kevin's NEW item)

(a) DONE on `todo/correct-timing-counter` 950f74d (pushed): `correct_read_3prime` times 2E pre-pass /
2F / 2C / 2G / 2E walkback / NET-seq + `read_total` per read; `ProcessingStats.module_seconds` sums
them (primary rows only; merges; region-checkpoint JSON round-trip); `correct_command` adds the stage
walls (2H, pool setup, per-read loop, writers); surfaced as `module_seconds_*` rows in the stats TSV,
a "Module wall time" report block, one `[TIMING] correct modules:` log line, and `stats.module_seconds`
in the correct sidecar. 9 new tests + a sidecar assertion; 103 neighbors green.
(b) DONE 567f30d (same branch, pushed): the resolver writes `XW:Z:<side>:<token>:<bits>:<window>:<clip_len>`
for every assessed clip (`splice/resolver_verdict.py`); 2F SKIPS the sequence search on `low_info` (same 200 bp
slice + information bound on both sides; structural paths live; applied in the ≥30 bp short-circuit, the body
floor and the peel re-entry), still TRIES after `rejected_edit` / `ambiguous` / `no_candidates` / `blowup` /
`repeat`, and COUNTS a sequence rescue there as `ends_2f_rescued_over_resolver_refusal` beside
`ends_2f_skipped_by_resolver`. 22 tests; regression set 476 passed / 0 failed. NOT landed on master (Kevin).
(c) Only if (a) shows 2F above a few percent: cross-arm cache keyed on read name + clip sequence.
Measure on the 748 fixture + holdout1k (tester asked for one human run's `module_seconds_*` rows).

### mapPacBio Index Caching
**File:** `rectify/core/align/multi_aligner.py` — `run_map_pacbio()`
**Status:** `nodisk` removed, `path=bbmap_index/` added (2026-03-29). Index pre-built for
bundled S. cerevisiae genome. New genomes will build on first run and cache thereafter.

---

## Analysis / Downstream

### ✅ DONE — Bedgraph and Genomic Distribution in Manifest Mode
**File:** `rectify/core/analyze/manifest.py` — `_run_analyze_manifest()`
**Resolved (option 1):** Manifest mode now generates BOTH. Bedgraphs accumulate
per-condition from either the position index (CPA counts) or the full-TSV streaming
path (`manifest.py:515-741`); genomic distribution runs as a lightweight **Pass 3**
(`manifest.py:744-813`) that streams only the needed columns
(`alignment_start`/`alignment_end`/`five_prime_position`/...) one sample at a time and
gracefully skips the transcript-body distribution when those columns are absent (rather
than the old blanket skip). Gated by `--no-bedgraph` / `--no-genomic-distribution`. The
function docstring documents the Pass-3 design. Bedgraph emitters are unit-tested
(`tests/test_analyze.py` — `generate_bedgraphs`, 0-based-shift regression).
**Optional follow-up:** a dedicated end-to-end `_run_analyze_manifest` integration test
asserting bedgraph files appear (currently only the emitter functions are unit-tested).

### ✅ DONE 2026-06-24 — Expose Aligner Stats in HTML Report
Resolved together with the Statistics/Observability items above (`tied_score` +
`by_aligner_combo` now in `*.consensus_aligner_stats.tsv` and the consensus HTML
report). Tests: `tests/test_consensus_stats_surfacing.py`.

### ✅ DONE 2026-06-24 — Emit cluster COM (read-weighted center-of-mass) from `analyze`
**Resolved:** `cluster_com` column added to both `cluster_cpa_sites()` and
`cluster_cpa_sites_adaptive()` in `rectify/core/analyze/clustering.py` (floor of
the count-weighted mean position, same coordinate convention as
`modal_position`/`start`/`end`). Flows to `cpa_clusters.tsv` in both standard and
manifest mode (both call the adaptive clusterer). Tests:
`tests/test_analyze.py::TestClusterCenterOfMass`. **Convention: 0-based, matching
`modal_position`/`start`/`end`/`corrected_position`** (NOT the "1-based" the original
spec mentioned — internal row consistency wins). ⚠ Before retiring the post-hoc
`projects/TRT/.../62_compute_cluster_com.py`, verify it computes COM in the same
(0-based) convention — if it emits 1-based, values differ by 1 and the script can't be
dropped as-is. Original task below.

### Emit cluster COM (read-weighted center-of-mass) from `analyze`
**File:** `rectify/core/commands/analyze_command.py`
**Priority:** Medium

`rectify analyze` currently emits per-cluster bedgraphs and a cluster table, but the
**read-weighted center-of-mass position** within each cluster is not exported. This
forces downstream consumers (motif scans, IGV navigation, cross-sample COM comparison)
to recompute it post-hoc by streaming the full `corrected_3ends.tsv` and grouping
by cluster ID — see `projects/TRT/scripts/rectify/han2023/62_compute_cluster_com.py`
as the canonical workaround.

**Intended behaviour:** add a `cluster_com` column (1-based genomic coordinate, integer
floor of read-weighted mean) to the existing per-cluster summary TSV that `analyze`
emits. Read weight = count of reads whose corrected 3' end falls within the cluster.
For motif-window extraction, this is the anchor point — the existing cluster start/end
coordinates do not capture it (they're peak boundaries, not signal centroid).

Acceptance: a manifest-mode analyze run produces a cluster table with `cluster_id`,
`chrom`, `strand`, `start`, `end`, `cluster_com`, `n_reads`, ... and the post-hoc
script `62_compute_cluster_com.py` can be retired.

---

## Validation

### ✅ ALREADY DONE (verified 2026-06-25) — DRS validation bundle is current, NOT stale
**Verification:** The committed `validation_reads.minimap2.bam` for cat2_plus_1 (`61b0c014`)
already ends in `…49M 9D 39M 8S` (3′≈23759) — the **through-alignment** that was the *fix*,
NOT the old stale `49M 47H` clip. gapmm2 BAM = 36/36 primaries (built with working 25.4.5).
`tests/test_validation_reads.py` expects cat2_plus_1 → 23754 (mapPacBio winner) and **passes**
(green in the full `not slow` suite). The bundle was regenerated in commit `c69524a`
(PROVENANCE run #35) after this TODO was written. **Recommendation: do NOT regenerate** — a
re-run through a non-identical align env would only risk injecting noise into a
currently-correct, test-passing regression backbone. The NEW-082 gapmm2 fix unblocks a
*future* regen if ever needed, but none is needed now. `corrected_3ends.tsv` is absent but is
NOT a required test/runtime input (only referenced in comments). The cat7 length→support
yeast port applies only *if* the yeast set is regenerated — N/A here.
_(Original task below for history.)_

### Regenerate the DRS validation bundle — stale cat2_plus_1 minimap2 clip
**Files:** `rectify/data/validation/aligners/*.bam`, `rectify/data/validation/rectified/`,
`dev/validation_review/.../cat2_softclip_review_pngs/cat2_plus_1.png`
**Priority:** Medium · **Blocked by:** the gapmm2 read-dropping bug (BUGS_TO_FIX NEW-082)

The committed bundle's `minimap2.trimmed.bam` for cat2_plus_1 (read `61b0c014`, chrI A-tract)
carries a **stale soft-clip** (corrected 3′=23711, EER ED 56.6, `…49M 47H`). Root-caused
2026-06-01 as a **one-off artifact of the original bundle build** — NOT walkback/trim/version/
flags/junc-bed/arch. Every current minimap2 (incl. production 2.30-r1287 on H2/Sherlock) aligns
this read **through** with a clean single 9-bp deletion (`…49M 9D 39M 8S`, 3′=23759). Full
diagnosis: `dev/validation_review/validation_read_review/cat2_softclip_findings.md` (RESOLUTION
section) + memory `project_cat2_plus1_minimap2_clip`.

**To do (focused pass, once NEW-082 is fixed):** re-align the validation reads with all 5 aligners
(needs working gapmm2 + a configured H2/prod align env), `update_validation_aligner_bams.py` the
fresh records, re-run `rectify correct`, re-render via `generate_review_report.py --arm drs`, and
**update the `tests/test_validation_reads.py` assertions** for this read (3′ 23711→~23758). Confirm
minimap2/gapmm2/uLTRA join deSALT/mapPacBio's through-cluster (EER 56.6 → ~15).

### ✅ RESOLVED 2026-06-25 — `corrected_3ends.tsv` is not a required artifact
**Verified:** `rectify/data/validation/corrected_3ends.tsv` is NOT git-tracked and NOT a test
or runtime input (referenced only in comments — `walkback.py`, a test docstring). The tracked
validation TSVs are `corrected_reads.tsv`, `corrected_reads_stats.tsv`, `corrected_3ends_stats.tsv`,
all current. No regeneration needed. _Original task below._

### Regenerate `rectify/data/validation/corrected_3ends.tsv`
**File:** `rectify/data/validation/corrected_3ends.tsv` (currently deleted in working tree)
**Priority:** Medium

PROVENANCE.json run #35 (2026-04-27) lists this TSV as an expected output of validation
BAM regeneration but it was not re-emitted alongside the new aligner BAMs. M1 dev work
will likely add additional validation reads as we debug the QuantSeq REV strand-flip
behaviour and the SLURM-script-gen bug in `rectify split --generate-slurm --short-read`,
so a full re-emit of the validation set (BAMs + TSV with consistent checksums) is
deferred until that work is more settled.

---

## Algorithm Semantics — Pending Review

### ✅ DONE 2026-06-24 — Inverted indel-correction tests verified
**File:** `tests/test_indel_correction.py`
**Resolved:** The 2026-05-10 inversion is VERIFIED correct against production.
`find_polya_boundary` delegates to `walkback_3prime_guarded` with `stop_base='A'`
(plus) / `'T'` (minus); the walkback stops at the first read==genome agreement whose
base is NOT the strand's poly-A `stop_base` (`walkback.py:169,181`). So T=T (plus) and
A=A (minus) ARE genuine exon boundaries. Justification: plus is probabilistic (~99% —
Nanopore HP errors are deletions, not substitutions, so a faked T=T is ~1%); minus is
deterministic (a minus tail base reverse-complements to 'T', so A=A can't come from the
tail). The stale class docstring (claimed `gb not in ('A','T')`) and the NOTE/TODO block
were corrected to document the verified contract. Both tests pass; full suite green.

---

## Validation-set selector — latent cat7 bug (for future yeast agents)

### cat7 (alt-splice) selection ranks by intron length, not aligner support
**Priority:** Low — LATENT, NOT actionable now. The flaw lives in the validation-set
**down-selection script** (used only when *regenerating* the committed yeast validation set),
not in the shipped rectify package. Yeast introns max ~800 bp so "longest" never pulled
artifacts → the committed yeast set is unaffected. Already fixed in the human scratch scripts.
**Action only IF a future agent regenerates the yeast validation set** — port the
length→cross-aligner-support ranking then. No package code change required. _(detail below)_

The validation-set down-selection ranks cat7 (non-canonical unannotated junction)
candidates by **longest intron**. On human (A549 chr5) this enshrined single-aligner
277–453 kb mismapping artifacts as "correct behavior" — exactly wrong for a regression
fixture. Fixed in the **human** scratch scripts (2026-06-12) by ranking on
**cross-aligner support** (`same_intron` agreement count across the aligner panel — the
same machinery cat5/6/9 already use; guardrail-compliant: prune by support, never by
motif). See `~/igv_data/a549_validation/scripts/{classify_candidates,select_validation}.py`
and `rectify/dev/handoffs/STATUS_human_validation_readset.md`.

**Yeast has the same latent flaw** — yeast introns max ~800 bp, so "longest" never
pulled artifacts, but the heuristic is still wrong in principle. **If a future yeast
agent regenerates the committed yeast validation set, port the length→support cat7 fix**
so the yeast cat7 exemplars are multi-aligner-confirmed, not just long.

### ✅ DONE 2026-06-25 — Reconcile Sherlock `correct_command.py` `--netseq` arm back to M1
**Resolved:** Sherlock's `/oak/.../software/rectify/.../correct_command.py` is now
**byte-identical to M1 HEAD** (`diff` = 0 lines, both 1974 lines, 41 netseq refs each). The
`--netseq` arm was already reconciled into M1's committed history (part of the 26 commits
ahead of origin), and the prior "M1 WIP" on this file has also been committed — M1 working
tree now differs from HEAD only by this session's Dorado pt-tag (C2) edits. No divergence
remains.

---

## Poly(A)

### ✅ DONE 2026-06-25 — Dorado pt:i poly(A) tail-length integration
**Files:** `rectify/core/polya/polya_trimmer.py` (`get_dorado_polya_length`,
`calculate_full_polya_length`, `trim_polya_from_read`), `bam_processor.py`,
`bam/parallel.py`, `commands/correct_command.py` (CLI `--use-dorado-polya`, default off).
Reads the Dorado `pt:i:` tail-length estimate; **always** records it as
`dorado_polya_length` in the per-read result dict; under `--use-dorado-polya` it becomes
the authoritative `polya_length`. Default-off → byte-identical behavior (the new dict key
is ignored by the positional `corrected_reads.tsv` writer, so the TSV schema and golden
hashes are unchanged — a dedicated TSV column was intentionally NOT added to avoid a
schema break; deferred/gated if ever wanted). Tests:
`tests/test_polya_trimming.py::TestDoradoPolyaIntegration` (7, incl. the bundled
`validation_reads_dorado_source.bam` pt-tagged reads). Replaces the long inline TODO at
`polya_trimmer.py` (now a concise docstring note).

## Known Non-Issues

- `minratio=0.4` in `run_map_pacbio()` cmd — duplicates BBTools default, harmless.
- Combined TSV (`corrected_3ends_combined.tsv`) in `combined/` is now obsolete for the
  standard pipeline (manifest mode replaces it) but harmless to leave in place.

---

## Human DRS `correct`-stage junction defects (Sumner panel, 2026-09-04)

Filed by the Sumner tester session from a 100-read misplacement panel
(`dev/sumner_misplaced_panel_20260904/`, see `SUBAGENT_BRIEF.md`). Each item has a full write-up in
`dev/sumner_misplaced_panel_20260904/ISSUES/`. Regression guards requested with every fix: human-DRS
fixture from the panel; FP side (no annotated canonical junction moved off annotation without ≥1.0
edit evidence; no non-canonical N-op created by 2F/2H); FN side (yeast cat3/cat7/refiner tests green).

- [x] **ISSUE-001** FIXED ee766c5 `todo/sumner-panel-issues` (on `todo/*` branch, NOT merged; tester verification pending — job 42115177 on fcafd1a) — `standardize_chrom_name` romanizes human chr names in spawned workers (chr5→chrV, chr10→chrX). CRITICAL. `ISSUES/ISSUE-001_*.md`
- [x] **ISSUE-002** FIXED 37a7685 `todo/sumner-2f-rescue` (extend gated on 'N lands on icp' + Fix 1b exon CIGAR from _intronic_seq; TSV↔BAM agreement on rescued rows 42/228 → 27/51 = all; refused-surgery TSV claim still OPEN → S2) — 2F rescue: BAM N-op ≠ TSV junction (extend-before-reroute); 42 % fire rate; junk exon CIGARs. HIGH. `ISSUES/ISSUE-002_*.md`
- [x] **ISSUE-003** FIXED 7cf6e97 `todo/sumner-panel-issues` (on `todo/*` branch, NOT merged; tester verification pending — job 42115177 on fcafd1a); docs now say 2H-before-2F — 2H silently skipped without `--aligner-bams`; docs/code disagree on 2F/2H order. MEDIUM. `ISSUES/ISSUE-003_*.md`
- [x] **ISSUE-004** FIXED 7d571f5 `todo/sumner-panel-issues` (on `todo/*` branch, NOT merged; tester verification pending — job 42115177 on fcafd1a); 81/81 → 0/81 wrong, pool phantom 85 % → 2 % — `_iter_n_ops`/pool walkers never advance past N → phantom coordinates for every intron after the first; drives R1/R2/R3 and 85 % pool contamination. **CRITICAL — fix first.** `ISSUES/ISSUE-004_*.md`
- [x] **ISSUE-005** (a) FIXED 21dbc57; (b) FIXED bb71473 as a NEW signal `move_dist` (L1 distance from the aligner's placement, 5th ordering key; byte-identical on panel + hold-out); PLUS d8b27e8 D/N-twin 2H fix found by the hold-out (`_dn_run_extent`, annotation decides) PLUS e85c0d0 `_CANONICAL_TIER_MAX=1` (a GT donor with no AG acceptor is not canonical — three sites; panel FN 3→1 / TP 8→10, 0 FP cost; the remaining FN 1c48253d is a U12 AT-AC alternative the read's own sequence rejects by 1.32 ED — left by decision) — formerly: `todo/sumner-panel-issues` (on `todo/*` branch, NOT merged; tester verification pending — job 42115177 on fcafd1a); (b) `abs_delta` tie-breaker BLOCKED-BY-DESIGN (`test_junction_refiner.py:225` pins delta == 0 — needs a new signal, Kevin); (c) move-guard defaults not invented; R1 evidence gate 250a608 landed instead (`RECTIFY_ANNOT_CANON_HOLD=0` to A/B) — `correct` loads no human penalty/STR tables; `abs_delta` tie-breaker dead; move guards 0.0. HIGH, land after 004. `ISSUES/ISSUE-005_*.md`
- [x] **ISSUE-006** FIXED e533ec4 `todo/sumner-2f-rescue` (informative-clip floor from the information bound = 10 nt, matching the shortest yeast cat3 clip exactly; Case-4 guard reachable; writer refuses self-invented non-canonical introns; hold-out non-canonical created 158 → 0, rescued 228 → 51, annotated share 11 % → 96 %; FN 10 → 14 flagged) — 2F rescues 1–3 nt clips (`ed_exon==0`, no min informative length); Case-4 guard unreachable on softclip branch → the 42 % fire rate. HIGH. `ISSUES/ISSUE-006_*.md`
- [x] **ISSUE-007** FIXED 5cb2b90 `todo/sumner-2f-rescue` (refuse when the trim CUT a junction; pre-existing N-ops destroyed 3 → 0; L1 rows verifiable only on the integrated tree with 004) — `reroute` merge fuses a PARTIALLY trimmed N with the new N → 97 kb "intron" replaces six real junctions (L1 rows). HIGH. `ISSUES/ISSUE-007_*.md`
- [x] **ISSUE-008** FIXED faf93d1 `todo/sumner-station-c` (+ census bound → parameter default 30, follow-up in progress) (on `todo/*` branch, NOT merged; tester verification pending — job 42115177 on fcafd1a); created-junction visibility 16/121 → 80/121 (97/121 at 30 bp) — Station C blind to 86 % of created junctions (I/D planted beside N by 2H defeats the anchor rule); GTF parser yields 0 introns; human tracks dead. HIGH. Proposed hermetic test in the panel dir fails 2/3 on HEAD. `ISSUES/ISSUE-008_*.md`
- [x] **ISSUE-009** FIXED f95d582 `todo/sumner-station-c` (on `todo/*` branch, NOT merged; tester verification pending — job 42115177 on fcafd1a); opt-in `rectify triage --original-bams`, panel reverts 7/94 → 80/94 (upper bound), C1 5/5 unchanged; run-all wiring DONE 1ef63f5 (`--triage-original-bams` nargs='*': absent = inert, bare = auto-discover `<sample>.multialigned.bam`, paths = override; Station C over the triaged BAM: admit 14→0, report 39→18 rows) — Station B proposer = 2H itself; arbiter would accept the original back (16/16 R1) but it is never proposed. MEDIUM. `ISSUES/ISSUE-009_*.md`
- [x] **ISSUE-010** FIXED 8427583 `todo/resolver-a7-a9-a8` (cap derives from each clip's information-bounded W, floor 2000; yeast byte-identical; human hold-out 2 → 8 resolved / 31 → 2 abandoned; `_warn_blowup` organism-neutral; NEW `_warn_slow_window`. COST: pure-Python DP 87 s → 2,259 s per 1k reads — Station A on mammalian genomes REQUIRES RECTIFY_HP_ED_NUMBA=1, ~135 s) — Station A candidate cap + yeast-calibrated skip-region advice on human chromosomes. LOW/MED. `ISSUES/ISSUE-010_*.md`
- ISSUE-002 mechanism CORRECTED (see file): reorder alone does not fix; softclip-only exon CIGAR + icp published for all rescue types; gate `extend` on "N lands on icp".
- [x] **ISSUE-011** FIXED 26263d6 `todo/sumner-panel-issues` (on `todo/*` branch, NOT merged; tester verification pending — job 42115177 on fcafd1a); cache stamp (+ prescan writer), self-pool + missing-chrom WARNINGs, `--junction-min-observed-support` (default 1; A/B at 2) — pool discards support (flat Set); singletons full-weight; self-pool warning; prescan caches contaminated (delete + stamp). MED-HIGH, after 004. `ISSUES/ISSUE-011_*.md`
- ISSUE-001 mechanism CORRECTED: parent-side ordering (registration at :840 AFTER the 2H block at :641), not spawn; hoist + roman-map guard. ISSUE-003: 2H-is-a-pre-pass CONFIRMED with writer line cites.
- [x] **ISSUE-012** FIXED 68d35fb `todo/sumner-station-c` (relative bound gap > 1.0 × hp_ed(original): fires on 1.2 % of bypassed hold-out reads, catches both panel reads; the information bound would have fired on 84.7 %) — Station B classifier bypasses huge-Δ hp_ed reads as high_confidence (fixer S3 finding). MED. `ISSUES/ISSUE-012_*.md`
- [x] **ISSUE-013** FIXED 4a230df `todo/sumner-station-c` (pre-gate = p99.5 × 2 of annotated intron length: yeast 5,000 exactly, GENCODE chr5 310,100; run-all `--pool-gate` wiring DONE 1ef63f5 — logged as `[station-c] length pre-gate: N bp (2x p99.5 …)`) — human max_intron clamps to 500,000 → Station C length pre-gate dead; ties to 010. LOW-MED. `ISSUES/ISSUE-013_*.md`
- [x] **ISSUE-014** FIXED 96f1c06 `todo/resolver-a7-a9-a8` (`--junction-aligners` help: "99% of" = the % o conversion; `tests/test_cli_help_all_subcommands.py` sweeps every subcommand, 75 checks) — `rectify align --help` crashes (unescaped `%` in a help string; pre-existing on 837990c). LOW. `ISSUES/ISSUE-014_*.md`
- [x] **ISSUE-015** (2026-09-05: ALL THREE LAYERS DONE — L1-2 2a88cdc `todo/sumner-2f-rescue`: admission = contig bounds (`clip_exceeds_contig`), placement invariant once in `apply_5prime_rescue_surgery` (`placement_off_contig`; the miss was `reroute`'s unbounded `new_ref_start`), 25 kb backstop cap only — GENCODE first exons reach 12,358 nt so an exon-1 cap could not have caught this read; L3 07523af `todo/sumner-panel-issues`: guarded sort+index at all 5 sites, offenders → <out>.unindexable_reads.tsv, run continues) 2F rescue writes reference_start < 0 (chrM 8 kb clip → pos −5843) → index fails, run aborts. HIGH. `ISSUES/ISSUE-015_*.md`
- [x] **ISSUE-016** FIXED `todo/issue-016-station-c-ledger` (2026-09-05 ~15:30 PT; VERIFIED by the tester ~16:10 PT on the 145k: 797 created → annotated 512 / reported 146 / refused 139 / not_seen 0 at ops2 — 83 % censused vs 13 %; residual 139 = read_end 69 / indel_ops 41 / indel_bp 29, a design call for Kevin): the 693 were NOT anchor refusals — the table lists non-annotated junctions only and keys them at the LEFTMOST ambiguity-equivalent coordinate; created junctions are mostly annotated and sit at motif coordinates. Now `<prefix>.pool_gate.annotated.tsv` + `<prefix>.census_refusals.tsv` (reason per side per raw N-op) + `start_raw`/`end_raw` columns + JSON `census` block + `rectify pool-gate --attribute FILE` (panel: 10 created → 9 annotated + 1 annotated_not_seen, 0 refused) — Station C on real 145k output: 693/797 created junctions ABSENT at any op budget / min-support → need per-N-op census-refusal reasons. MED. `ISSUES/ISSUE-016_*.md`
- [x] **ISSUE-017** FIXED `todo/2f-novel-evidence-gate` (2026-09-05 ~17:30 PT; tester verification pending) — 2F places 1–8 matched bases across a NEW intron onto UNANNOTATED pool candidates (added_nov = 83 % of cohort FP; Station C `read_end` residue). Provenance-aware evidence gate: annotated landing sites keep today's acceptance; novel ones need ≥10 matched bases, no junction-side gap, bounded indel burden (`five_prime_rescue_refused` = `novel_exon_*`). HIGH. `ISSUES/ISSUE-017_*.md`
- [x] **ISSUE-018** FIXED `todo/issue-018-2h-noncanon-preponderance` 58ffabb + 64bed7a (2026-09-05 ~17:45 PT; tester verification pending — prediction: class A 47 → 34 gated, to_noncanon 74 → ~30 as labeled today / ~3 if annotated destinations stop being FP-labeled, rescue_annot Δ 0, controls 0) — 2H moves a canonical junction onto an UNANNOTATED non-canonical site with a compensating indel (47 class-A rows, anchor 0 on a side 47/47): the pool is built through `_merge_del_into_intron`, so read A's D-merged form is read B's "alternative" and escapes the d8b27e8 own-twin rule; the vanished deletion cost (1.07–1.41 units) clears the R1 hold. Mechanism replayed op-for-op on 10/47 rows (real hg38 slices). RULING 3 rule (b) as a preponderance gate read from a dry run of the surgery (no alternative within δ, no I/D beside the new N, no burden rise, both anchors ≥ informative floor); counter `noncanon_destination_refused` + reason tokens; `RECTIFY_2H_NONCANON_PREPONDERANCE=0` = pre-fix. Class B (annotated destination) untouched by design — Kevin call. Leading-soft-clip q_split window shift found here, out of scope, filed separately as ISSUE-021. HIGH. `ISSUES/ISSUE-018_*.md`
- [x] **ISSUE-020** FIXED `todo/issue-020-anchored-ranking` (2026-09-05 evening PT; commit f53d770 on top of d6acaaa, base e3906ce; tester verification pending) — 2F RANKED 5' rescue candidates with an unanchored hp-ED shift x offset sweep (a free junction-side gap of up to 10 nt; the GTRAGT +4 decoy won only through it, `_off` = 2/10/7/3/3 on the six residual reads) but PLACED them with the anchored affine aligner. Now the anchored DP (`score_right_anchored` / `score_left_anchored`, same core as `align_clip_to_exon`) is the FINAL rank, hp-ED only prunes shifts; the compared segment always ends at the junction (clip + intron-mapped bases when the 5' end is inside the intron; dist > 0 trims exon-2 bases from the READ, never slides the genome); consistency invariant asserted (`tests/test_2f_anchored_ranking.py`, `RECTIFY_2F_CHECK_CONSISTENCY=1`); counter `five_prime_reranked_between_annotated`. All 10 bundle reads annotated, gap token 0, the three ISSUE-023 within-1 reads now draw the annotated N-op. HIGH. `ISSUES/ISSUE-020_2F_anchored_ranking.md`
- [x] **ISSUE-021** FIXED `todo/issue-018-2h-noncanon-preponderance` ee811b6 (2026-09-05 ~18:04 PT, separate commit after 018; tester verification pending — on the 74 real to_noncanon records the scored window equals the true exon-2 start 6/74 → 74/74; 4/68 soft-clipped decisions change (2 to_noncanon → annotated, 1 stock → annotated, 1 new non-canonical → non-canonical move), 62 change scores; 018-gate refusals 28 → 10; validation reads changed 0) — 2H indexed `read.query_sequence` (clip INCLUDED) with `_iter_n_ops`' `q_split` (leading clip EXCLUDED), so every N-op of a soft-clipped read was scored on a 30-nt rescue window shifted len(clip) bases into exon 1 (68/74 to_noncanon records carry a leading S; DRS reads almost always do). Fix: `q_split += read.query_alignment_start` once at the loop head; `_iter_n_ops`, the scorer and the surgery untouched; unclipped reads byte-identical. `tests/test_2h_softclip_window.py` (10). Bundle `holdout/events/ee811b6/rbrowse/issue021_softclip_10/`. HIGH. `ISSUES/ISSUE-021_*.md`
- [x] **ISSUE-025** FIXED `todo/issue-018-2h-noncanon-preponderance` 7174c47 (2026-09-05 ~18:27 PT; verified on the M1 at `--threads 4`: all five upf1Δ `correct` children `refine_bam_junctions done`, mapPacBio 20 N-op reads / 7 refined, merged cat9_plus_1 = 439093-439323; the 017 branch's darwin xfail on cat9_plus_1 should now be removed) — 2H's `_run_parallel` relied on fork-inherited `_WORKER_POOL_STATE`; on spawn (macOS) every worker raised KeyError('header') and `correct` skipped 2H as "non-fatal" — 2H did nothing at ≥ 2 threads with no marker. Fix: `ctx.Pool(initializer=_init_worker_pool_state, initargs=(state,))`, platform-default start method (`RECTIFY_2H_MP_START_METHOD` overrides), pool failure → RuntimeError; `correct` logs at ERROR, writes `module_2h_failed` as the first `_stats.tsv` row, `RECTIFY_2H_FAILURE_FATAL=1` aborts. `tests/test_2h_parallel_driver_spawn.py` (9). Linux production unaffected (fork). HIGH (M1 T0 validity). `ISSUES/ISSUE-025_*.md`
