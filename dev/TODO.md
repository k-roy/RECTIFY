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
