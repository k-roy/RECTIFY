# Session Handoff — mpb 5'-edge reanchor wired via path (b)

**Date:** 2026-05-18 (late evening, continuing from `845a573`)
**Branch:** `drs-validation-rebuild`
**Top of branch:** `560f82c` — reanchor wired via TSV persist (path b)
**Predecessor handoffs:**
- `handoffs/HANDOFF_2026-05-18_6943450_cat3p2_cat2m2_reanchor.md` — shipped
  cat3_plus_2 + cat2_minus_2 fixes; deferred the reanchor hooks.

---

## 1. What was done

Wired the deferred `reanchor_5prime_for_rescue` function via "path (b)"
from the predecessor handoff (persist as TSV column). The five
production-code edits are scoped to a single commit (no test changes —
the existing `tests/test_validation_reads.py` covers the surface).

### `rescue_3ss_truncation` (`splice_aware_5prime.py`)

Added a reanchor pre-pass at the top of the function body. The function
saves `(read.cigartuples, read.reference_start)`, runs
`reanchor_5prime_for_rescue` in-place, compares the resulting CIGAR
against the saved one, and only treats the reanchor as material when
the cigartuples differ. On a material reanchor, the new leading-S
length is captured as `_reanchor_clip_len`. The rest of the function
body — now extracted into `_rescue_3ss_truncation_body` — runs against
the (possibly reanchored) read so that
`align_clip_to_exon` sizes `five_prime_exon_cigar` against the
reanchored geometry. A `try/finally` restores the read state so callers
see no mutation. On the rescued return path, the result dict gains
`reanchor_clip_len`.

The non-material gate is load-bearing: uLTRA-style reads with a
pre-existing leading S that already matches the first ≥10 match run
cause `reanchor_5prime_for_rescue` to return `True` while emitting a
cigartuples list identical to the input. Without the comparison the
TSV would silently carry a phantom clip_len equal to the original S
length, and `bam_writer` would call reanchor on the raw read for no
reason (still a no-op, but a noisier control path).

### TSV plumbing (`bam_processor.py`, `output.py`)

`bam_processor.py` reads `reanchor_clip_len` from the rescue result and
threads it into the per-read result dict. `output.py` adds the column
to the schema (after `five_prime_upstream_trim`) and writes it on
every row (default `0`).

### `bam_writer.py` — three write paths

`write_corrected_bam`, `write_softclipped_bam`, and `write_dual_bam`
each gained a reanchor pre-pass *before* `realign_exon_blocks`,
gated on `correction['reanchor_clip_len'] > 0`. Same
`anchor_min_run=10` as `bam_processor`, so the function is
deterministic and both call sites produce the same leading S. This
inverts the previous session's failed wiring (reanchor happened
*after* realign in bam_writer), which was the documented cause of the
six regressions.

`_load_corrections_from_tsv` parses the new column with the same
"optional column" pattern as the surrounding fields (graceful default
when the column is absent in legacy TSVs).

### Docstring (`read_edits.py`) + queue entry (`debugger_queue.md`)

`reanchor_5prime_for_rescue`'s DEFERRED docstring replaced with a wired
description. `docs/handoffs/debugger_queue.md` "Deferred: mapPacBio
5'-edge reanchor" entry replaced with "RESOLVED — see commit X" plus
implementation summary and verification notes.

---

## 2. What's verified

- `pytest tests/test_validation_reads.py --no-header -q` →
  **107 passed, 8 skipped** — matches the predecessor baseline.
  Targeted run `-k "cat3 or cat4 or cat7"` (the prior-session
  regression hot list): 45 passed, 1 skipped.

- Targeted CIGAR-rewrite check on the bundled per-aligner BAMs
  (`rectify/data/validation/aligners/`) — script run inline, not
  committed:
  - mpb cat3_plus_1: pre `1X2=7I86=2I65=…`, post `10S86=2I65=…`,
    `reference_start` 168805 → 168808 (the canonical 5' marker for the
    YBL027W intron).
  - mpb cat3_minus_1: pre `…10=5I1=2I2=1I5=`, post `…10=16S` (trailing
    cluster collapsed; reference_start unchanged for − strand).
  - mpb cat3_minus_2 + every deSALT/minimap2/gapmm2 cat3 read: correctly
    no-op (function returns False or the cigar-equality gate suppresses).
  - uLTRA cat3 reads: function returns True but `_cigar_after ==
    _cigar_before`, so the material-change gate emits `reanchor_clip_len
    = 0`. Verifies the gate is doing real work.

### Bundle regenerated in-session — end-to-end results

`scripts/validation_data/regen_pa_rest_bundle.py` ran cleanly
(36 reads, all spliced; winner counts:
minimap2=0, gapmm2=3, mapPacBio=8, deSALT=19, uLTRA=6).

**`reanchor_clip_len > 0` in merged `corrected_reads.tsv`:** 2 reads.

- `8f86cb34-…` — `clip_len=6`. mpb is the per-read winner (all 5
  aligners tie at HP-ED 4.12 after rescue; mpb wins on tiebreak). The
  reanchored mpb alignment was clean enough to match the others.
- `d3357db5-…` — `clip_len=19`. uLTRA wins (HP-ED 12.37); mpb's HP-ED
  26.96. mpb's `_five_rescued=0` despite the reanchor — rescue failed
  for a different reason (probably the 5' bases past the reanchor
  still don't match exon-1 closely enough).

**Cat3 read consensus after reanchor:**

| read | mpb HP-ED | winner HP-ED | mpb `_five_rescued` | winner |
| --- | --- | --- | --- | --- |
| cat3_plus_1 (0a28167d) | 18.91 | 10.16 | **1** (was 0) | deSALT |
| cat3_minus_1 (ac4db6da) | 24.90 | 9.90 | **1** (was 0) | deSALT |
| cat3_minus_2 (28ea9379) | 20.83 | 20.63 | 0 (correctly no-op) | deSALT |
| cat3_plus_2 (79f61403) | 27.60 | 26.21 | 0 (mpb didn't need rescue) | deSALT |

**The reanchor is doing exactly its job** — it enables mpb's 5'-rescue
to succeed on cat3_plus_1 and cat3_minus_1, where it previously
failed. mpb's HP-ED still trails the winner cluster on those reads
(probably from body mismatches the reanchor doesn't touch), so the
*winner* didn't flip on cat3_plus_1 / cat3_minus_1.

**cat3_plus_2 is a separate winner-selection finding (user, on plot
review):** mpb's raw alignment for cat3_plus_2 is *already* correct —
exon-1 tail is `…14= 1D 9=` with both terminal G's properly in exon 1
and a clean 366-bp intron. mpb's `_five_rescued=0 /
correction_applied=none` is the *right* answer there: no rescue
needed. The winner cluster (deSALT/gapmm2/minimap2/uLTRA) needed
rescue (raw alignments had the 365-bp off-by-one) and ended up with
`…14= 1D 7= 1D 1=` — same ref end position, but with one terminal G
pushed off into a `1D 1=` split instead of staying inside the clean
`9=`. HP-ED scores deSALT 26.21 vs mpb 27.60 — only 1.39 points
apart, and mpb loses even though its junction is the more
parsimonious one (per the plotter image
`validation_read_review/cat3_junction_review_pngs/cat3_plus_2.png`).

So the actionable finding for the *next* session is:
**HP-ED is under-penalizing the winner's `1D 1=` tail relative to
mpb's clean `9=` tail.** This is the same family of complaint as the
existing "Cat1 cluster (HP-ED metric)" entry in
`docs/handoffs/debugger_queue.md` — terminal `1D ... 1=` splits at
junction boundaries should cost more than a single equivalent `D` in
the middle of a clean match run, but currently they don't. Worth
tracing the per-op HP-ED contributions for both alignments to confirm
exactly which weights are mis-calibrated.

The earlier "mpb wins on body-quality HP-ED after reanchor"
prediction from the predecessor handoff does not hold here — but
that prediction was about cat3_plus_1 / cat3_minus_1, where the
mpb HP-ED gap *is* a body-quality issue, not a junction-tail one.

- **Pre-existing tiebreaker test failures in
  `test_corrected_consensus_tiebreaker.py`.** Unchanged from the
  predecessor handoff. Not in scope here.

---

## 3. Open items

### Diagnose why mpb's body HP-ED is still high on cat3_plus_1 / cat3_minus_1

The reanchor produced clean 5'-rescue for mpb on these two reads, but
mpb's overall HP-ED is still 8–15 points above deSALT/minimap2/uLTRA.
That suggests body mismatches the reanchor does not address. Pull up
the mpb BAM CIGAR for both reads in the regen bundle and compare to
the deSALT/uLTRA CIGAR for the same read; the rest of the alignment
(post-reanchor) is the source of the gap.

If the diagnosis points at a real-data correctness issue in mpb's
output (e.g. forced mismatches in the body), it overlaps with the open
"Cat1 cluster" item where mpb force-aligns 4+ mismatches into the body
to anchor at a non-A past a poly-A tail. May want to revisit the
HP-ED scoring weights or mpb's penalty configuration.

### cat3_plus_2: HP-ED winner-selection should pick mpb

mpb's raw alignment is already canonical (no rescue needed). The
winner cluster's rescued alignment has a less-parsimonious `1D 1=` at
the exon-1 tail. HP-ED scores them within 1.4 points and picks the
winner cluster anyway. Trace per-op HP-ED contributions for
cat3_plus_2 deSALT vs mpb and find which weight (likely either the
match-after-D bonus or the per-base soft-clip penalty) under-counts
the winner's split-D tail. Shares root cause with the open
"Cat1 cluster (HP-ED metric)" entry. The user's prior call ("don't
need a hack; once CIGARs converge after canonical-intron fixes,
HP-ED naturally decides") needs revisiting — HP-ED is in fact
deciding *wrong* here.

### Defensive belt-and-suspenders (optional)

Bam_writer's reanchor gate is `correction['reanchor_clip_len'] > 0`,
relying on the invariant that `bam_processor` only emits non-zero
clip_len on `five_prime_rescued=True` paths. If a future change breaks
that invariant the gate would fire on un-rescued reads (with no
extend_read_5prime to consume the new S, so likely a no-op but
strictly wasteful). Adding `and correction.get('five_prime_rescued',
False)` to the gate would harden it. Not done because (a) the
invariant is clear in `bam_processor.py:411-413`, (b) the wasted call
on a hypothetical regression is still a no-op.

### Items inherited from the predecessor handoff (still open)

- **Cat3_plus_2 consensus pick** — same as above; bundle regen
  required to confirm mpb wins.
- **Cat3_plus_1 chimera exemption** — per user 2026-05-18 don't need
  a hack; HP-ED decides naturally once CIGARs converge.
- **Pre-existing tiebreaker test failures** in
  `test_corrected_consensus_tiebreaker.py` (NaN in `_eff_key`).
- **Ambiguity-window + motif-strength tiebreaker** for consensus.
- **Cat1 HP-mode metric design** (architectural).
- **Phase C 5' rescue calibration application**.

---

## 4. Resume command

```bash
cd /Users/kevinroy/work/rectify
git log --oneline HEAD~2..HEAD  # see what shipped this session
.venv/bin/python scripts/validation_data/regen_pa_rest_bundle.py
# then: inspect rectified/per_aligner_summary.tsv winning_aligner for cat3 reads
.venv/bin/pytest tests/test_validation_reads.py -q --no-header  # sanity (expect 107/8)
```

---

## 5. Files touched

**CODE (to be committed):**

- `rectify/core/splice/splice_aware_5prime.py` — reanchor pre-pass +
  body-split refactor.
- `rectify/core/bam/bam_processor.py` — pick up + thread
  `reanchor_clip_len` into the per-read result dict.
- `rectify/core/bam/output.py` — new TSV column.
- `rectify/core/bam/bam_writer.py` — parse the new column; apply
  reanchor before `realign_exon_blocks` in three write paths.
- `rectify/core/bam/read_edits.py` — docstring update only (function
  body unchanged from prior session).

**DOCS (to be committed):**

- `docs/handoffs/debugger_queue.md` — close out the DEFERRED entry as
  RESOLVED with verification notes.
- This file (`handoffs/HANDOFF_2026-05-18_reanchor_wired.md`).

**`[uncommitted]` — pre-existing WIP, NOT touched this session:**

- `rectify/core/commands/correct_command.py`,
  `rectify/core/splice/junction_refiner.py`,
  `docs/figures/*.png/.svg`, various `*.pre_regen` TSVs, and the
  other ~155 entries `git status` shows. Surgical staging only.
