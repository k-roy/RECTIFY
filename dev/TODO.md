# RECTIFY TODO

Known gaps, planned improvements, and deferred work.
Add items here rather than inline `# TODO` comments where possible.

---

## Alignment / Consensus

### Fix deSALT on Sherlock — Install Conda Binary
**Priority:** High
**Context:** deSALT's vendored binary at `rectify/data/bin/linux_x86_64/deSALT`
(v1.5.6) crashes with SIGSEGV on every chunk on Sherlock, producing 274B empty
fallback BAMs. The same tool works correctly on H2 where it is conda-installed
(`~/.conda/envs/rectify/bin/deSALT`, 856 KB, built 2025-09-29): 5065 alignments,
99.5% mapped, real CIGAR strings with splice N-ops. See
`docs/ALIGNER_RECOMMENDATIONS.md` for full details.

**Fix:** install deSALT on Sherlock via conda (same channel/version as H2).
After fixing, re-run set2 and set3 merges and correct with deSALT enabled.

**Also investigate:** the `rectify correct` hang observed in smoke test 25410215
when empty (274B) deSALT BAMs were used as primary input — confirm whether the
hang is triggered by the empty BAM itself or is a separate bug in
`rectify.core.bam.parallel`.

**Agent task:** full step-by-step instructions (install, crash-window validation at
14,900 / 30,000 / 40,000 reads, result interpretation, and rollout if passing) are
in `project_status_markdowns/TASK_desalt_conda_sherlock_test.md`.

---

### Evaluate Newer Aligners as Panel Additions (Minisplice, GLASS, Winnowmap2)
**Priority:** Medium
**Context:** Regardless of the deSALT Sherlock fix, newer aligners have emerged
that may improve consensus quality or provide algorithmic diversity.

Candidates:
- **Minisplice** (~June 2025) — minimap2 + deep-learning splice signals
- **GLASS** (April 2025) — graph-learning splice-aware alignment
- **Winnowmap2** — optimized for repetitive regions; lower priority for yeast

Integration path: add a wrapper in `rectify/core/align/multi_aligner.py`, register
in `SUPPORTED_ALIGNERS`, validate with `pytest tests/test_consensus_selection.py`.

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

### Tie-Break Rate per Sample
**File:** `rectify/core/consensus/consensus.py` — now tracks `stats['tied_score']`
**Status:** Counter added (2026-03-29). Not yet surfaced in HTML report or corrected_3ends_stats.tsv.

### Aligner Combo Breakdown per Sample
**File:** `rectify/core/consensus/consensus.py` — now tracks `stats['by_aligner_combo']`
**Status:** Logged at INFO level (2026-03-29). Not yet surfaced in HTML report or stats TSV.

---

## Performance

### mapPacBio Index Caching
**File:** `rectify/core/align/multi_aligner.py` — `run_map_pacbio()`
**Status:** `nodisk` removed, `path=bbmap_index/` added (2026-03-29). Index pre-built for
bundled S. cerevisiae genome. New genomes will build on first run and cache thereafter.

---

## Analysis / Downstream

### Bedgraph and Genomic Distribution in Manifest Mode
**File:** `rectify/core/commands/analyze_command.py` — `_run_analyze_manifest()`
**Priority:** Medium

Manifest mode currently skips bedgraph generation and genomic distribution analysis
because both require per-read alignment coordinates (`alignment_start`, `alignment_end`),
which are not stored in the position index. Two options:
1. During Pass 2, fall back to streaming the full TSV for samples that need these columns
2. Add a `rectify export --bedgraph --manifest` subcommand (cleaner — bedgraphs are
   large and users may not always want them)

### Expose Aligner Stats in HTML Report
**Files:** `rectify/core/consensus/consensus.py`, HTML report template
**Priority:** Low

`stats['tied_score']` and `stats['by_aligner_combo']` are tracked but not shown in
the HTML report or `corrected_3ends_stats.tsv`. Surfacing these would help diagnose
aligner quality per sample.

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

### Inverted indel-correction tests (test_indel_correction.py)
**File:** `tests/test_indel_correction.py`
**Priority:** Medium

Two tests were renamed and had their assertions inverted (commit pending push, 2026-05-10):
- `test_plus_strand_T_match_is_genuine_exon_boundary` (was `_skips_T_match_before_true_CPA`)
- `test_minus_strand_A_match_is_genuine_exon_no_correction` (was `_skips_A_match_*`)

The new names assert that T=T (plus) and A=A (minus) matches are GENUINE exon boundaries
rather than poly-A false stops to skip past. **Author intent and biological justification
need explicit verification** before treating this as a stable contract; see the inline
NOTE block at the top of the regression section in the test file.

---

## Validation-set selector — latent cat7 bug (for future yeast agents)

### cat7 (alt-splice) selection ranks by intron length, not aligner support
**Priority:** Medium — latent; masked on yeast, surfaced on human 2026-06-12.

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

### Reconcile Sherlock `correct_command.py` `--netseq` arm back to M1
**Priority:** Medium — cross-cluster divergence.

Sherlock's rsync rectify (`/oak/.../software/rectify`) carries an **uncommitted `--netseq`
arm** in `rectify/core/commands/correct_command.py` that is *ahead* of M1 HEAD `9f613a6`
(every delta is `is_netseq`-gated → DRS/cDNA paths byte-identical, which is why it was
safe to leave during the 2026-06-12 human-validation work). M1 ALSO has its own
uncommitted WIP on this file. These need reconciling into one committed state on M1 (M1
authoritative; do NOT discard the Sherlock netseq work — cherry-pick it back per
CLAUDE.md). Until then the file differs across M1-WIP / M1-HEAD / Sherlock.

---

## Known Non-Issues

- `minratio=0.4` in `run_map_pacbio()` cmd — duplicates BBTools default, harmless.
- Combined TSV (`corrected_3ends_combined.tsv`) in `combined/` is now obsolete for the
  standard pipeline (manifest mode replaces it) but harmless to leave in place.
