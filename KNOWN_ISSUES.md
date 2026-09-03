# Known issues

Defects that are **currently open, or fixed only on an unmerged branch**. If you are about to
trust RECTIFY output, read this first.

This file is **tracked**, so it travels with every clone, worktree and cluster checkout — unlike
the maintainers' internal notes, which are deliberately not published. Entries here are written to
be useful without access to any lab infrastructure.

> ⚠️ **SCOPE.** Created 2026-08-29; a triage pass over the maintainers' internal record was
> completed the same day. **Every entry below was verified against the current `master` source**,
> not taken from a stale label — two internally-recorded issues turned out to be already fixed and
> were deliberately excluded (see *Verified fixed* at the end). A handful of narrow
> environment-specific operational notes remain unpublished because they concern one lab's cluster
> layout rather than the tool. **Absence from this file is weak evidence, not proof, that something
> works.**

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

---

## 🔴 A per-aligner fan-out silently clobbers its own output — OPEN

- **Status:** open on `master`, verified 2026-08-29
- **Affects:** `rectify align` run as one job per (sample × aligner) with a shared `--prefix` and `-o`
- **Impact:** the "consensus" BAM is not a consensus — it is whichever aligner finished last

The single-aligner output path is built as `<output_dir>/<prefix>.multialigned.bam` — keyed on the
**sample prefix only, not on the aligner**. Running the panel as separate per-aligner jobs (the
correct pattern, since uLTRA rebuilds a whole-genome index per invocation and must not be chunked)
gives three writers and one destination. **Last writer wins, and nothing errors.**

The symptom is quiet: the resulting file is byte-identical to one aligner's BAM apart from an extra
`@PG` record, and every job log cheerfully reports `Single-aligner output (sorted+indexed): …`.

- **Workaround:** give every per-aligner job its own `--prefix` (or its own `-o` directory), then
  merge deliberately. Never let two aligner jobs share both.
- **Detection:** compare the size of the "consensus" BAM against each single-aligner BAM. If it
  matches one of them within a few dozen bytes, you have this bug.
- **Fix guidance:** the aligner name belongs in the output prefix.

---

## 🔴 `correct` cost scales with junction-pool DENSITY, and the candidate cap does not gate it — OPEN

- **Status:** open on `master`, verified 2026-08-29
- **Affects:** `rectify correct` on reads crossing junction-dense loci
- **Impact:** runtime, not correctness — but it can look like a hang

The 3′SS rescue runs an O(n·m) dynamic program **once per candidate junction**, and while the DP
bounds its *sequence length*, **nothing bounds candidate multiplicity**. Cost therefore tracks how
junction-dense the locus is rather than how many reads there are, so a modest BAM over a dense
region can take far longer than a larger one over a sparse region.

`--junction-max-candidates-per-nop` looks like the relevant control but **is not**: it is consumed
only by junction *refinement* (`junction_refiner.py`); the 3′SS rescue path never receives it.
Verified — that parameter appears in no other module.

- **Workaround:** none via that flag. Budget wall-time by locus density, not read count, and chunk
  junction-dense regions more aggressively.

---

## 🔴🔴 TRAP — a FASTA/GFF contig-naming mismatch makes the whole panel annotation-blind, silently

- **Status:** user/environment trap, not a code defect — but it produces confidently wrong output
- **Impact:** annotated junctions silently absent; every downstream junction measurement invalid

Reference FASTA and annotation frequently disagree on contig naming — bare (`I, II, … XVI`) versus
`chr`-prefixed (`chrI, chrII, … chrmt`) — and the two files often ship **side by side in the same
directory**, so the pairing looks correct. minimap2 builds its junction BED from the annotation,
matches those names against the index contigs, finds nothing, and **proceeds without error**.

**This matters well beyond alignment quality.** The 5′ soft-clip rescue draws its candidate pool
from *annotated junctions ∪ cross-aligner novel junctions*. With the annotated half empty, the pool
is novel-only — precisely the unfiltered short-anchor substrate that the anchor-length gate exists
to control. Any rescue or junction-pool number taken on a mismatched reference is measuring
something other than what it claims.

- **Detection — do this once per reference directory, before any run:** compare the contig names in
  the FASTA (or its `.fai`) against column 1 of the annotation. They must match exactly, including
  the mitochondrial contig, which is the one most often named differently in an otherwise matching
  pair.
- **Workaround:** keep matched FASTA/annotation pairs together and verify, rather than assuming
  co-location implies compatibility.

---

## 🔴 TRAP — `minimap2 -y` with a non-SAM FASTQ comment aborts samtools, and the error is masked as SIGPIPE

- **Status:** usage trap with a diagnosis hazard
- **Impact:** the run dies with a misleading error that points at the wrong tool

`minimap2 -y` copies the FASTQ comment field into the SAM record as auxiliary tags. If the comment
is not valid `TAG:TYPE:VALUE` SAM syntax, the downstream `samtools` fails with
`[E::aux_parse] unrecognized type` and aborts, which closes the pipe, so **minimap2 dies of
SIGPIPE and that is the error you see.** The real cause is invisible because the pipeline checks
minimap2's return code before samtools'.

- **Workaround:** ensure every FASTQ comment is well-formed SAM aux syntax before using `-y`, or
  drop `-y`.
- **Diagnosis:** run the alignment and the sort by hand and capture **both** stderr streams;
  samtools names the problem precisely.

---

## Verified fixed — deliberately NOT listed as open

Both were recorded internally as outstanding but were checked against the current `master` source on
2026-08-29 and found already fixed. Recorded here so they are not re-filed:

- **`run-all` per-region manifest handling** (single-sample producing no analysis; multi-sample
  crashing on a missing `corrected_reads.tsv`) — fixed; the sample column now goes to the region
  tables and the manifest stays an index.
- **ONT PCR-cDNA via `run-all`** — a Path A route (pre-alignment → UMI collapse) is implemented, so
  the previous advice to use only the explicit `trim-cdna-polya → minimap2 → correct` chain no
  longer applies.

⚠️ Historical-data caveat, not a live bug: ONT-cDNA outputs produced **before** `gene_id` was
assigned on the gene strand carry corrupt `gene_id` for antisense reads. Regenerate rather than
trust `gene_id` in old cDNA outputs.
