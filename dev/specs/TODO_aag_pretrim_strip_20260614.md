# TODO (fresh agent): strip the (AAG)ₙ/(GAA)ₙ triplet-repeat basecaller artifact in DRS pre-trim

**Status:** scoped, not started · **Drafted:** 2026-06-14 · **Owner:** unassigned
**Read first:** memory `project_sumner_lab` (the "(AAG)n/(GAA)n re-confirmed" 2026-06-14 bullet) and
`triplet-repeat-expansion-artifact-readside`. This spec is self-contained otherwise.

---

## 0. TL;DR for a context-less agent

ONT DRS on the **Dorado v5.2.0 / RNA004** modern cohort mis-basecalls the poly-A **homopolymer tail as
`(AAG)ₙ/(GAA)ₙ` triplet repeats** (a read-side basecaller artifact, period-3, k=3 enriched ~278×). RECTIFY's
DRS pre-trim (`rectify trim-polya` → `drs_trim_command.py`) scans for an **A-run** from the 3′ end, hits the
GAA block immediately, and reports `polya_len ≈ 1`. So the artifact tail is **NOT trimmed** — it's carried
into alignment, where the aligner force-aligns it (the cat2_minus_2 11D-bridge mess) and the genome-aware
walkback has to clean it up.

**Your job:** add a triplet-repeat **strip** to the pre-trim — detect a terminal `(AAG/GAA/AGA)ₙ`-type
low-period repeat block, peel it off (like an adapter), then run the existing poly-A scan on what's behind
it. Reuse the existing motif-agnostic detector. Validate it recovers/cleans the 3′ end on Sumner data and
is a **strict no-op on clean data** (old cohort, yeast validation, R9.4). This is a **trim-level fix +
validation only** — a full Sumner re-align/re-correct of all 11 samples is a SEPARATE downstream decision,
not in this scope.

**This is a mitigation, not the root fix.** The upstream fix is re-basecalling with a model that handles
the poly-A homopolymer (Stephen Brown's domain). Don't block on that; ship the trim-level strip.

---

## 1. Evidence (don't re-derive — measured 2026-06-14, CNTL_21.8 chr5, 30k primary reads)

- **43.4%** of reads have a terminal `(GAA/AAG)ₙ` repeat vs only **19%** clean poly-A(≥8/10).
- **99.1%** of those GAAs are in the **unaligned 3′ soft-clip** → read-side artifact, not a genomic GAA locus.
- **4.8%** have a genuine poly-A **blocked behind** a terminal GAA (the clearest "recoverable poly-A" case).
- Median `polya_len = 1`; **`--max-error-rate` does NOTHING** (0.0 vs 0.1 on 20k reads: median 1→1, ≤5
  80.7%→80.7%, 6% of reads shift). The issue is "no A-run exists," not error tolerance. **Do NOT re-try the
  error-rate angle.**
- The "43% of trimmed fastq reads are empty" was a FALSE ALARM = the 78,602 SECONDARY alignments (no SEQ via
  `samtools fastq`); trim metadata has 0% emptied.

## 2. The code (verified line numbers, `rectify/core/commands/drs_trim_command.py`, HEAD 2026-06-14)

- **`_scan_polya`** (line 54) — right-to-left A-run scanner. NOT the place to change (it correctly finds the
  A-run; the problem is the GAA block sits between the 3′ end and the poly-A).
- **`find_polya_and_adapter`** (line 120) — the 3-pass detector. **This is where the strip goes.**
  - `_ADAPTER_RE = re.compile(r'[TU][CTU]{0,10}$')` (line 44) — Pass-1 adapter stub. Does NOT match GAA.
  - Pass 0/1 strip the adapter (line ~158-166), then `_scan_polya` runs (line 173), then Pass 2 iterative
    peel (line ~178-194). **Add the triplet-repeat strip so the poly-A scan runs on the post-repeat window.**
  - Returns `(polya_len, adapter_seq, last_base, adapter_pass)`.
- **Three trim paths** each do `total_trim = polya_len + len(adapter_seq)` then build the removed-bases record
  `trimmed_3prime_seq`: `trim_drs_bam_polya` (428/435), `trim_drs_fastq_polya` (561/564),
  `_process_mapped_read` (628). The parallel path is `_trim_drs_bam_polya_parallel` (737) via
  `_trim_region_task` (663) / `_trim_unmapped_task` (698). All consume `find_polya_and_adapter`'s return —
  so if you account for the stripped repeat **inside that function's return**, all callers stay correct.

## 3. Reuse the existing detector — do NOT write a new motif scanner

`rectify/core/splice/repeat_expansion.py` already does exactly the structure detection you need
(motif-agnostic, period-1 homopolymers excluded so it never fires on a clean poly-A):
- `dominant_repeat_period(seq, min_frac, k_range)` (line 40) → `(k, motif, frac)` or `None`.
- `is_repeat_expansion(seq, ...)` (line 76) → `True` if `seq` is a multi-base low-period repeat; **rejects
  pure poly-A/poly-T** and keeps A-rich `(AAG)ₙ` like `"AAAAAGAAGAAG"` (line 87-88 docstring). Exactly the
  AAG-vs-polyA discriminator you want.
- `_canonical("GAA") == _canonical("AAG") == "AAG"` (line 35), `DEFAULT_K_RANGE=(2,3,4,5)`.

This is also the **no-op guarantee**: on clean poly-A tails `is_repeat_expansion` returns `False`, so the
strip never fires → clean data (old cohort / yeast validation / R9.4) is untouched.

## 4. Implementation

**Add a triplet-repeat peel to `find_polya_and_adapter`** (after the Pass-1 adapter strip, before/around the
poly-A scan). Algorithm:

1. From the 3′ end of the (adapter-stripped) window, identify the **maximal terminal low-period repeat
   block** — walk left while the trailing window stays a period-2..4 repeat (use
   `dominant_repeat_period` on a growing/shrinking suffix, or scan period-3 phase directly). Require a
   minimum block length (e.g. ≥9 nt / ≥3 tandem copies) and high purity (`min_frac` ~0.8) so a couple of
   incidental `…GAA` codons in real sequence don't trigger it. Gate with `is_repeat_expansion` so pure
   poly-A is never peeled.
2. **Peel that block**, then run `_scan_polya` on the window **behind** it → recovers the poly-A for the
   4.8% blocked case; for the (majority) GAA-replaced-tail case there may be no recoverable poly-A, and
   that's fine — the value is removing the artifact junk before alignment.
3. **Return the peeled length.** Cleanest: add a 5th return value `repeat_len` (extend the tuple) and have
   the 3 callers fold it into `total_trim` (`total_trim = repeat_len + polya_len + len(adapter_seq)`) +
   record motif/len in metadata. Minimal-diff shortcut (acceptable, but muddies semantics): fold the repeat
   block into `adapter_seq` so `total_trim` already covers it and **no caller changes** — if you do this,
   tag `adapter_pass`/metadata so the artifact is still distinguishable from a real adapter.

**Metadata (`_write_metadata`, schema ~line 18-19):** add `repeat_len` and `repeat_motif` columns so the
artifact prevalence is queryable post-trim. (Parquet preferred; TSV fallback already handled.)

**CLI:** expose a flag to enable/disable + tune (e.g. `--strip-repeat-expansion` default ON for DRS,
`--repeat-min-len`, `--repeat-min-frac`) via `create_trim_polya_parser` (line 977) and thread through `run`
(line 836). Default-ON is reasonable since it's a strict no-op on clean tails, but confirm with Kevin.

## 5. Validation (REQUIRED)

1. **Unit tests** (extend `tests/test_polya_trimming.py` — verify it's the right file): synthetic RNA-oriented
   reads — `[body][polyA20][(GAA)10]` must strip the GAA AND recover ~20 polyA; `[body][(GAA)10]` strips the
   GAA, polya_len≈0; a **clean** `[body][polyA30]` is UNCHANGED (no-op); a real sequence with an incidental
   `…GAAGAA` (2 copies) is NOT stripped (purity/length gate).
2. **No-op proof on clean data:** re-trim the OLD cohort (`CNTL_4.2`/`SMA_7.12`, clean poly-A) and a yeast
   validation BAM → `polya_len` distribution UNCHANGED. If clean data shifts, the gate is too loose.
3. **Sumner recovery:** re-trim `CNTL_21.8` (trim-only, cheap — reuse the §7 harness) → terminal-GAA
   fraction should drop from 43% to ~0; the 4.8% blocked-poly-A reads recover realistic `polya_len`.
4. **Full validation suite** (`tests/test_validation_reads*.py`) still green — the yeast Cat1–9 path uses
   this trimmer; the strip must be a no-op there.

## 6. Guardrails (do NOT regress)

- **Strict no-op on clean poly-A** (period-1 homopolymer) — `is_repeat_expansion` guarantees this; keep it.
- **No motif gating** beyond "is it a low-period repeat" — stay motif-agnostic (lab principle
  `no-motifs-unbiased-discovery`); don't hard-code `GAA`. `repeat_expansion.py` is already motif-agnostic.
- **Don't strip genuine genomic 3′ ends:** gate on length + purity. drs_trim is pre-alignment so you can't
  cheaply use soft-clip status, but 99.1% of the artifact is soft-clipped and genuine 3′UTRs rarely end in a
  pure multi-copy (GAA)ₙ block — the length/purity gate is the protection. Flag borderline cases in metadata.
- **All three trim paths + the parallel path** must stay consistent (they all flow through
  `find_polya_and_adapter` — account for the strip in its return, not per-caller).
- **0-based-inclusive `corrected_3prime`** convention downstream is unaffected (this is pre-alignment).

## 7. How to verify on real data (reuse this session's harness)

- Input: `/scratch/users/kevinroy/sumner_lab/chr5_bams/CNTL_21.8_chr5.bam` (Sherlock).
- Re-trim demo pattern (subsample 20k PRIMARY → trim → compare `polya_len` parquet), used 2026-06-14:
  ```bash
  samtools view -h -F 0x900 $IN | awk '/^@/{print;next}{if(c<20000){print;c++}}' | samtools view -b - > sub.bam
  rectify trim-polya sub.bam -o out_new --strip-repeat-expansion   # vs out_old without
  # compare polya_len distributions from out_*/sub_polya_trim_metadata.parquet (pandas)
  ```
- Measure terminal-GAA fraction before/after with the period-3 + soft-clip probe from the
  `project_sumner_lab` 2026-06-14 investigation (re.compile(r'(GAA|AAG|AGA){3,}') on RNA-oriented 3′ ends).
- Sherlock rectify is a non-git rsync copy; deploy changed `.py` from M1 via
  `rsync -a -f'+ */' -f'+ *.py' -f'- *' rectify/ sherlock:<oak>/rectify/rectify/` (the validation env needs
  `MAX_SS_SHIFT` etc. → sync ALL `.py`, not just the one file — learned the hard way 2026-06-14).

## 8. Scope boundary (read before expanding)

- **IN scope:** the trim-level strip + unit tests + the CNTL_21.8 recovery measurement + the clean-data
  no-op proof. Land it on `drs-validation-rebuild` (or the active feature branch) with CHANGELOG/AGENT_FIXES.
- **OUT of scope (separate decision, needs Kevin):** re-running the full Sumner pipeline (re-align + 3-aligner
  consensus + re-correct all 11) with the new trim — that's significant compute and re-stamps every downstream
  result. Quantify the expected benefit from the CNTL_21.8 trim measurement first, then decide.
- **Also out of scope:** the cat2_minus_2 softclip_rescue accuracy fix (128102→128117) — tracked separately
  in AGENT_FIXES [2026-06-13/14]; it's a different (post-alignment) layer.

## 9. Deliverables checklist

- [ ] Triplet-repeat strip added to `find_polya_and_adapter`, accounted for in its return (all callers correct).
- [ ] Reuses `repeat_expansion.is_repeat_expansion` / `dominant_repeat_period` (no new motif scanner).
- [ ] `repeat_len`/`repeat_motif` metadata columns + CLI flag(s).
- [ ] Unit tests (recover-behind, strip-replaced, clean-no-op, incidental-codon-not-stripped) green.
- [ ] No-op proven on old-cohort + yeast validation; full validation suite green.
- [ ] CNTL_21.8 re-trim: terminal-GAA 43%→~0; blocked-poly-A reads recover realistic `polya_len`. Numbers in AGENT_FIXES.
- [ ] CHANGELOG + AGENT_FIXES updated; this TODO marked done; memory `project_sumner_lab` updated.
