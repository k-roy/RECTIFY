# Known issues

Defects that are **currently open, or fixed only on an unmerged branch**. If you are about to
trust RECTIFY output, read this first.

This file is **tracked**, so it travels with every clone, worktree and cluster checkout — unlike
the maintainers' internal notes, which are deliberately not published. Entries here are written to
be useful without access to any lab infrastructure.

> ⚠️ **SCOPE — this list is NOT yet a complete audit.** It was created 2026-08-29 and seeded with
> the two most recent issues. The maintainers' internal record holds roughly ten further open
> entries (alignment fan-out, reference-naming traps, `run-all` paths, cost scaling) that have not
> yet been triaged into public form. **Absence from this file is not evidence that something
> works.** Triage is tracked as an open task.

**Conventions.** Each entry carries a status, the affected versions/commits, a workaround, and
where the fix lives. **Closing an entry is part of the fix**: when a fix merges to `master`, delete
its entry here and record it in `CHANGELOG.md`. A stale issues list that claims fixed things are
broken is worse than no list at all.

---

## 🟡 Type-2 (no-UMI) cDNA reads are deduplicated by coordinate — fixed on a branch, not on `master`

- **Status:** fixed on `feat/cdna-stage1-qc` (`599260c`); **`master` still has it**
- **Affects:** `rectify correct-cdna` (ONT PCR-cDNA Stage 1), all versions up to and including `47e3b39`
- **Impact:** Type-2 record counts understated ~2×, **depth-dependently**

`correct-cdna` routed every Type-2 (SSP-less) anchor bucket through a coordinate collapse that
grouped reads on exact `(aln_start, aln_end)` and treated each group as PCR duplicates of one
molecule. **Type-2 reads carry no UMI**, so there is no evidence by which to call two of them the
same molecule; the collapse measured *positional concentration*, not amplification.

Measured on one 18-library cohort: **13,292,754 Type-2 reads → 6,450,950 records, 51.5 % removed.**
The rate scaled with sequencing depth — 4–6 % on ~50 k-read libraries versus 44–57 % on
multi-million-read libraries — while true UMI-measured PCR duplication on the *same* libraries was
24–41 %. PCR duplication is a property of a library, not of how deeply it was sequenced; positional
crowding is what scales that way. The excess was therefore genuinely distinct molecules merged away.

🔴 **Because the bias tracks depth, it does not cancel in a between-sample ratio.** Any comparison
of Type-2 abundance across libraries of differing depth is confounded.

- **Type-1 records are unaffected** — UMI-anchored deduplication was never in question, and Type-1
  is typically 82–88 % of reads.
- **Workaround on `master`:** treat Type-2 record counts as unusable for abundance or cross-sample
  comparison; use Type-1 only. Or run the branch above, where the default is `--type2-collapse none`.
- **After the fix:** each Type-2 read is one observation. Grouping Type-2 reads by 3′ end is still
  correct as **isoform / CPA-site clustering in `cdna-analyze` (Stage 3)**, where it is labelled as
  such — it is simply not deduplication and no longer happens in Stage 1.

Note: `docs/quickstart_cdna.md` already specified the correct behaviour ("Deduplication: None
(each read is one observation)"), so the code was violating its own documented contract.

---

## 🟡 Stage-1 cDNA QC is missing whenever `--workers > 1` — fixed on a branch, not on `master`

- **Status:** fixed on `feat/cdna-stage1-qc` (`7ff8f5c`); **`master` still has it**
- **Affects:** `rectify correct-cdna` with `--workers > 1` — i.e. effectively every production run
- **Impact:** no correctness effect on the output FASTQ; QC reporting only

The region-parallel path computed the read-type, XF-tier and tail-length metrics per region and then
discarded them, and the parent process aggregated only a fixed key list that excluded them. The
serial path printed a full QC block while the parallel path printed almost none of it, so runs
shipped with no read-type breakdown and the numbers had to be reconstructed by hand from the output
FASTQ.

- **Workaround on `master`:** rerun a subset with `--workers 1` to see the QC block, or derive the
  metrics from the `XT` / `XY` / `XC` tags in `stage1_consensus.fastq.gz`.
- **After the fix:** one shared implementation serves both paths (verified to produce identical
  output), adds UMI duplication rate and `XY` sub-type breakdown, and writes a machine-readable
  `stage1_qc.json` beside the FASTQ.

🔴 **Interpretation trap, independent of this bug:** the read-level and molecule-level Type-1
fractions are **different quantities**. The ~82 % figure documented in
`docs/algorithms/cdna_correct.md` is the **read-level** one; comparing a molecule-level fraction
against it reads as a false failure. Both are now reported and explicitly labelled.

---

## 🔴 Junction machinery is minor-intron (U12) blind — OPEN

- **Status:** open, not fixed
- **Affects:** junction scoring and splice-site indexing on U12-type introns
- **Impact:** AT–AC introns are unrepresentable and score at the worst tier

Four sites are involved: the splice-site index has no plus-strand `AC` acceptor kind
(`splice_site_index.py`); the canonical dinucleotide set omits AT–AC and AT–AG
(`overhang_informativeness.py`); plus-strand AT–AC scores at tier 8 and the treatment is
strand-asymmetric (`junction_scoring.py`); and the canonical homopolymer prior applies a 0.5-unit
handicap against true U12 junctions (`junction_scoring.py`).

Measured consequence: **92.5 % of STX10 long reads land on a phantom unannotated `AG` 5–6 nt off**
the true junction.

- **Workaround:** treat U12/minor-intron junction calls as unreliable; do not use RECTIFY junction
  tiers to adjudicate AT–AC introns.
- **Fix guidance:** address AT–AC and AT–AG together and bump `_FORMAT_VERSION` 2 → 3.

⚠️ The internal record points at a `TODO_MINOR_INTRON_GRAMMAR.md` patch spec at the repo root, but
**that file is not present in the tree** — the pointer is currently dangling.
