# Session Handoff — Verification of cat3 equivalence-extension (0653172); + strand path disabled

**Date:** 2026-05-18 (evening, continuing from `0653172_cat3_equivalence`)
**Branch:** `drs-validation-rebuild`
**Top of branch:** `01c2a18` (disable buggy + strand equivalence-extension)
**Predecessor handoff:** `HANDOFF_2026-05-18_0653172_cat3_equivalence.md` — read
that first; this one only records the verification outcome.

---

## 1. What was done

Ran the 5 verification checks from §4 of the 0653172 handoff. Four pass;
one (the + strand mirror geometry trace) caught a bug.

Detailed walk through each:

1. **+ strand mirror geometry trace** — FAIL. The previous commit's
   `_overshoot = read.reference_start - _intron_end > 0` trigger fires on
   body M starting PAST canonical intron_end (a gap between intron and
   body M). The structural mirror of - strand `ref_end > intron_start`
   (body M extending INTO intron) is the OPPOSITE sign:
   `_intron_end - read.reference_start > 0`. The case the buggy code
   handled produces a non-canonical intron 2k bases longer than the
   annotated junction (`canonical + 2k` from the
   `extend_read_5prime_for_junction_rescue` arithmetic — see commit
   message of `01c2a18` for the full derivation). There's no D op to
   absorb in the buggy case either: body M's first k bases are pure
   matches inside the downstream exon.

   **Action:** Disabled the + strand `elif` branch with an explanatory
   inline comment. The proper fix (the mirror, structurally identical
   to the deferred cat3_plus_2 "off-by-1 acceptor" undershoot case)
   is documented in `docs/handoffs/debugger_queue.md` as the first
   top-level note.

2. **k-sweep edge cases** — PASS by reading code.
   - `_overshoot == 0` → `if 0 < _overshoot:` False → no-op ✓
   - `_overshoot > 10` → `_max_k = min(overshoot, _MAX_K=10) = 10`, loop
     iterates `range(10, 0, -1)` ✓
   - `_overshoot == 2, k=2 fails, k=1 succeeds` → loop hits k=2 first,
     fails equality, continues; next iter k=1 succeeds → sets
     `_upstream_trim=1` and breaks ✓

3. **Joint `_borrowed == _ref_old AND _borrowed == _ref_new` check** —
   PASS. `_borrowed == _ref_old` excludes minimap2-style encoded
   mismatches (M op with non-matching bases at that ref position).
   `_borrowed == _ref_new` ensures the alternate landing position is
   also a real match. Both are necessary; neither is sufficient.

4. **Interaction with shift loop (lines 1191-1218)** — PASS, with caveat.
   The equivalence-extension runs AFTER `best_junction` is committed
   and uses whatever `(_intron_start, _intron_end)` the rescue chose.
   When that's slid by the shift loop, the equivalence check fires
   relative to the slid coordinates — geometrically valid. The
   "wrong placement" concern raised in §4.5 of the 0653172 handoff
   is about consensus selection (winner-takes-all on which aligner's
   slid junction to keep), not about the equivalence-extension
   transformation itself. That concern is covered by the separate
   ambiguity-window + motif-strength tiebreaker work queued at the
   top of `debugger_queue.md`.

5. **Synthetic + strand mirror test** — DELETED from the queue. Writing
   a test against the buggy code would have locked in non-canonical
   intron coordinates. The proper + strand fix needs its own test once
   the correct trigger sign is implemented.

---

## 2. What's verified

- `pytest tests/test_validation_reads.py --no-header -q` → **106 passed,
  8 skipped** (unchanged from 0653172 baseline). The disabled + strand
  branch was never exercised by the validation bundle, so no test
  flipped.
- cat3_minus_2 - strand equivalence-extension still fires and produces
  canonical CIGAR (all 5 aligners share `M256 N(366502-366584) M15`),
  confirmed by the existing
  `test_cat3_minus_2_rescued_aligners_have_clean_intron_cigar` test
  passing after the change.

**Not changed:**

- Behavior on - strand 5'-rescue equivalence-extension (the
  cat3_minus_2 path).
- All other rescue paths (mpb_mismatch, forced-snap, etc.).
- TSV column `five_prime_upstream_trim` is still written; remains 0 for
  + strand reads (was previously 0 in practice too, since the buggy
  branch never fired on validation data).

---

## 3. Open items

- **Proper + strand equivalence-extension** (replaces the disabled
  branch). Full sketch at the top of `docs/handoffs/debugger_queue.md`.
  Trigger should be `_intron_end - read.reference_start > 0`, and the
  borrowed-bases destination may differ from the - strand case (the
  surgery needs careful re-derivation because - strand routes borrowed
  bases through the rescue M, which is the DOWNSTREAM-exon side on
  - strand but the UPSTREAM-exon side on + strand). May need a
  separate code path in `extend_read_5prime_for_junction_rescue`.

- **cat3_plus_2 "off-by-1 acceptor" undershoot pattern** (4 of 5
  aligners). This is the same structural case as the proper + strand
  mirror — fixing one would naturally cover the other.

- All other items from §3 of the 0653172 handoff still open
  (ambiguity-window + motif-strength tiebreaker, Cat1 HP-mode metric
  design, Phase C 5' rescue calibration application, uncommitted
  raw aligner BAM WIP).

---

## 4. Resume command

If the user wants to ship the proper + strand fix:

```bash
cd /Users/kevinroy/work/rectify
git log --oneline 01c2a18^..01c2a18  # see what was disabled
sed -n '1230,1310p' rectify/core/splice/splice_aware_5prime.py  # disabled site
head -100 docs/handoffs/debugger_queue.md  # full sketch
```

The sketch in `debugger_queue.md` flags one tricky point: on - strand,
"borrowed bases" travel from body M's end INTO the rescue M's start
(both on the downstream-exon side); on + strand, the structurally
equivalent move is INTO the rescue M's end (upstream-exon side). The
asymmetry — borrowed bases moving to opposite sides of the intron —
means the BAM-writer surgery may need a + strand-specific code path,
not just a sign-flipped trigger. Re-derive the `intron_len` formula on
paper before coding to confirm it collapses to canonical.

If the user instead wants the ambiguity-window + motif tiebreaker
(the consensus-selection concern that originally motivated the
session), that's the second top-level note in `debugger_queue.md`
and is independent of the + strand fix.

---

## 5. Files touched

**CODE:**

- `rectify/core/splice/splice_aware_5prime.py` — replaced the + strand
  `elif` block in `rescue_3ss_truncation` with an inline explanatory
  comment. Lines roughly 1374-1404 (where the elif was). The
  - strand block immediately above (lines 1336-1372) is unchanged.

**DOCS:**

- `docs/handoffs/debugger_queue.md` — new top-level "Bug note" entry
  documenting the geometry inversion. Includes the smoking-gun
  arithmetic and a sketch of the proper fix.
- This file
  (`handoffs/HANDOFF_2026-05-18_01c2a18_plus_strand_disabled.md`).

**`[uncommitted]` — working tree items present BEFORE this session and
not changed here:**

- `CLAUDE.md`, `HANDOFF.md`, `README.md`, `docs/figures/*` — pre-existing
  WIP.
- `rectify/core/commands/correct_command.py`,
  `rectify/core/splice/junction_refiner.py` — pre-existing WIP from the
  parallel UMI-bin Phase D work track (see in-tree `HANDOFF.md`).
- `rectify/data/validation/aligners/validation_reads.{minimap2,uLTRA}.bam`
  — pre-existing WIP from before the 0653172 session (mtime
  2026-05-18 10:07).
