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

## 🔴 Chimeric consensus can emit a record assembled from TWO aligners — fixed on a branch, not on `master`

- **Status:** fixed on `fix/consensus-contig-swap-864` (which supersedes and includes the
  narrower `fix/runall-quantseq-862`); **`master` (`fd2e2d2`) still has it**
- **Affects:** any multi-aligner run with `--chimeric-consensus`, which is **the `run-all`
  default** (`run_command.py:1021-1034`). `rectify align` / `rectify consensus` default it OFF.
  The mechanism is datatype-independent — it lives in the disagreement fallback, not in anything
  short-read-specific — and it has now been measured on **two** short-read panels (below).
  Long-read panels are not yet measured; how often minimap2 + the resolver arms, which share a
  seed chain, disagree on contig is unknown.
- **Impact:** an output record can carry one aligner's **chromosome, strand, MAPQ and aux tags**
  with another aligner's **position and CIGAR**, while the `Xa` tag credits the second one.
  Sometimes it crashes; otherwise it is written silently to a locus no aligner reported.
- **Workaround on `master`:** `run-all --no-chimeric-consensus`, or a single-aligner panel.

`select_best_chimeric` refuses to assemble unless every candidate is on one contig AND one strand,
and otherwise falls back to `_fallback_simple_selection` — which returns the **winner's**
`chimeric_ref_start` and CIGAR **but no chromosome and no strand**.
`consensus._process_and_write_batch` then chooses a "template" read by **iteration order over
`aligner_reads.values()`**, gated only on sequence length, and
`chimeric_consensus.build_chimeric_read` builds a fresh record taking `reference_id`, the strand
bits, `mapping_quality` and the aux tags from that template. When the candidates disagree, the
record is `<template's chrom/strand>` × `<winner's position/CIGAR>`.

Measured on `fd2e2d2`:

| panel | rate |
|---|---|
| QuantSeq `--short-read --dT-primed-cDNA`, 200k reads, bbmap + bwa | 2,532 / 179,062 reads took the cross-contig fallback |
| TruSeq short-read panel, 200k reads SE, 5 arms | **115 / 157,475 = 0.073 %** written to the wrong chromosome |
| TruSeq short-read panel, 200k read pairs PE, 7 arms | **25 / 20,000 = 0.125 %** (first two checkpoint batches) |

Only ~8 % of the swapped records overran the borrowed contig and were caught by
`_validate_bam_sample`; the rest passed every guard. One QuantSeq example: read
`D00689:118:C890GANXX:8:2204:16881:55011` was placed by bwa at `chrIX:300228` (`2S48M`) and by
bbmap at `chrXIII:480619` (`1=1X47=1X`), and the consensus wrote `chrIX:480619` — past the end of
chrIX (439,888 nt). One TruSeq example: bbmap won with `chrI:165474 49=1X25=`, and the first arm
in the dict supplied chrX, MAPQ 60, FLAG 89 and its entire tag set.

**Fix:** `ChimericResult` now carries an **anchor** — the candidate whose placement
`chimeric_ref_start`/`chimeric_cigar` describe — and `build_chimeric_read` takes RNAME, strand and
MAPQ from it, with the template preferred to be the anchor itself. An **emit-time invariant**
refuses to write any record whose RNAME/POS/CIGAR/strand are not all the selector's, and, for a
pass-through (unstitched) result, whose triple is not exactly the winning candidate's own
placement. A no-op when all arms share a contig and strand, which includes every true-chimeric
assembly.

---

## 🟡 Type-2 (no-UMI) cDNA reads are deduplicated by coordinate — fixed on a branch, not on `master`

- **Status:** fixed on `feat/cdna-stage1-qc` (`599260c`), also merged into
  `feat/netseq-junction-rescue-836`; **`master` still has it**
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

- **Status:** fixed on `feat/cdna-stage1-qc` (`7ff8f5c`), also merged into
  `feat/netseq-junction-rescue-836`; **`master` still has it**
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

## ✅ NET-seq / QuantSeq defects — ALL SIX FIXED ON `master`

- **Status:** **fixed and merged** (`feat/netseq-junction-rescue-836` → `master` @ `457fe00`,
  2026-09-02). This entry is kept as the record of what the defects were and what the numbers
  looked like on either side of them; the "workaround" lines apply only to a checkout OLDER than
  `457fe00`.
- **Affected:** `rectify netseq`, `rectify correct --netseq` / `--dT-primed-cDNA`,
  `rectify run-all --short-read`, `rectify align --help`, `rectify correct --netseq-dir`

Six independent defects, each verified by execution before the fix:

1. **`rectify netseq` tracks are built on the RAW alignment terminus.** `aggregate_positions`
   defaults to `three_prime_raw` and the command never overrides it, so every bedgraph, BigWig and
   the deconvolution input ignore the module's own corrections; they survive only in a parquet
   column. *Workaround on `master`:* read `three_prime_corrected` from the parquet yourself.
2. **The minus-strand terminal-oligo(A) trim sign is inverted.** The correction must move the 3'
   end 5'-ward IN RNA, which is `-1` on a `+` gene and `+1` on a `-` gene; `master` subtracts on
   both, doubling the error on minus-strand reads (synthetic read: 97 returned, 103 correct). Latent
   for the tracks because of (1), but the exclusion-region test and the UMI molecule track consume
   the wrong value.
3. **Module 2F reads the wrong soft clip under any antisense protocol.** `rescue_3ss_truncation`
   takes coordinates from the gene strand but the clip from `read.is_reverse`
   (`splice_aware_5prime.py:363-372`, `:420`). Under `--netseq` / `--dT-primed-cDNA` those are
   opposite ends of the read, so it compares the RNA-3'-end clip — which holds the randomer and the
   poly(A) tail — against exon-1 sequence, and can write a false `five_prime_position`.
4. **`run-all --short-read --dT-primed-cDNA` exits 1.** `--short-read` is listed as a mutually
   exclusive protocol, so the documented QuantSeq invocation is unreachable. And `run-all`'s
   hand-built `correct` Namespace omits `short_read`, so once that is fixed, poly(A) trimming and
   the position-shifting indel module run on 150-bp Illumina reads with nothing erroring.
   *Workaround on `master`:* run `align` / `correct` / `split` stepwise instead of `run-all`.
5. **`rectify align --help` raises `TypeError: %o format`** — a literal `%` in a help string
   (`align_command.py:195-196`) that argparse tries to %-format. The help is unreachable.
   *Workaround:* read `docs/user_guide/commands/align.md`.
6. **`correct --netseq-dir` sums every loaded BigWig** regardless of strand or raw-vs-deconv
   (`netseq_refiner.py:225-229`). With the standard `<sample>.<kind>.<strand>.bw` layout a `+`
   query gets the minus strand added AND the deconvolved track added to the raw one — and the
   regularised NNLS does not conserve mass, so the inflation varies with local A-content.
   *Workaround on `master`:* point `--netseq-dir` at a directory holding one flavour only, and
   accept that both strands are still summed.

Landed alongside them, and missing capability rather than defects: the donor-side junction rescue
for NET-seq geometry (nothing in RECTIFY re-placed a 3' soft clip across a junction — Module 2F is
transcript-5'-only and the three 3'-side modules are contiguous-reference) and the randomer-aware
non-templated 3'-tail call. See `docs/algorithms/netseq_refinement.md` § *Donor-side junction
rescue*.

🔴 **Still open: there is no donor-side rescue inside `rectify correct`.** The rescue lives in
`rectify netseq`, which is where the NET-seq geometry is. `correct --netseq` on its own still
calls a spliced read with a soft-clipped exon-2 overhang at the 5' splice site — a false splicing
intermediate. For NET-seq use `run-all --netseq` (or `rectify netseq`), not `correct --netseq`.

---

## ✅ Mitochondrial group I/II introns were in the NET-seq junction-rescue pool

- **Status:** fixed 2026-09-03 (`--pool-include-organellar` restores the old behaviour).
  **A checkout between `b806e87` and that fix fabricates rescues on the mitochondrial genome.**
- **Affects:** `rectify netseq`, `rectify run-all --netseq` — any run given a BAM that retains the
  mitochondrial contig.

`JunctionPool.from_annotation` dropped **tRNA** introns (not spliceosomal) by parent feature type,
but yeast **mitochondrial** introns are annotated with parent type `mRNA`, exactly like a nuclear
intron, so the parent-type filter could never see them. The SGD GFF names the contig `chrmt`
(which standardizes to `chrMito`) and carries **32 intron features** on COX1 / COB / 21S — **group
I and group II self-splicing** introns, with no spliceosomal donor/acceptor grammar, on a genome
Pol II does not transcribe. A nascent Pol II 3' end cannot legitimately sit at one.

**Measured on wt_rep3 (2.36 M reads): 94 of 580 rescues — 16 % — were on chrMito.** Isolating the
47,023 chrMito reads reproduced the entire discrepancy against a mito-free reference run
(near-donor 210, rescued 94, exon1_end 45, ambiguous 12, intronic_end 59).

*Workaround on an older checkout:* hand `netseq` a BAM with the organellar contig removed
(`samtools view -b in.bam chrI … chrXVI`), which is what a hand-built NET-seq pipeline usually did
anyway. Nuclear results are unaffected either way — the pool entries are organellar only.

---

## 🟠 The NET-seq molecule (`--dedup`) track carries a randomer correction the read track does not

- **Status:** narrowed on `feat/netseq-junction-rescue-836`, not eliminated; **worse on `master`**
- **Affects:** `rectify netseq --dedup --umi-length N`

`netseq_umi.iter_netseq_fragments` used to key molecules on the terminal-oligo(A)-trimmed position
while the read track keyed on the raw alignment terminus, so the two disagreed on every read with an
aligned tail base — despite a docstring promising they "sit on the SAME positions". Both now run the
identical correction chain (trim → walkback → junction rescue).

**One deliberate divergence remains.** `randomer_overshoot(L, N)` returns `N - L` for every clip
shorter than the randomer, `L == 0` included, because on a library where every read carries an N-nt
randomer a zero-length clip means all N randomer bases aligned by chance. On a **mixed** library that
is false for the randomer-free class: measured on PRJNA1521488 wt_rep3 (chrI+chrII, `--umi-length 6`),
**50.1 % of reads are shifted, 52 % of reads have a zero-length clip**, and 11,152 of the 47,647
molecule positions do not exist in the read track.

- **Workaround:** do not use `--dedup` unless the randomer is UNIVERSAL in the library. Verify from
  the aligned 5'-clip histogram, not from the kit — for PRJNA1521488 the standing conclusion is no
  `--umi-length`, no `--dedup` (a barcode on ~30 % of reads cannot deduplicate the rest).
- The command now prints this warning whenever `--dedup` runs.

---

## 🟠 Two NET-seq correction defaults are calibrated, not obvious — and the summary JSON reports both sides

- **Status:** OPEN by design on `feat/netseq-junction-rescue-836`; not a bug, a calibration
- **Affects:** `rectify netseq --junction-rescue` and the poly(A) walkback

Both were calibrated on measurement and the defaults now reflect it; this entry records what the
defaults are and what reverting them costs.

1. **The rescue floor is remainder-aware, because the chance channel is the remainder.**
   `--rescue-min-k` (default 1) applies when the clip is exon-2 sequence and nothing else
   (`r == 0`); `--rescue-min-k-with-remainder` (default 4) applies when a randomer is invoked to
   explain the rest. Measured on a 194 k-read yeast slice, 504 candidate reads: the decoy-acceptor
   null produced `k=1` 70 times against 67 observed in the randomer channel, `k=2` 24 against 19,
   `k=3` 10 against 6 — at or above the observed rate — and **never reached `k >= 4`**. With the
   defaults the run rescues **160 reads at `decoy_rescued` = 0**; with a flat `min_k = 1` it
   rescued 227 at a chance floor of 54. **Read `decoy_rescued`, `rescued_by_k_clean` and
   `near_donor_k_randomer` in `<sample>.netseq_summary.json` before trusting a rescue count.**
   The floor is tested BEFORE the remainder, so a chance randomer match stays an `exon1_end`
   rather than being hidden in `ambiguous`; and tRNA introns are dropped from the pool by default
   (`--pool-include-trna` to keep them) because they are not spliceosomal.
2. **The poly(A) walkback is GATED on clip evidence by default in `rectify netseq`.** Invariant 7
   (a terminal read A over a genomic A is not skipped) is right when every read has a tail; most
   nascent 3' ends have none, and ~25 % of reads end on an A over a genomic A by chance. Measured
   unconditionally: 42,644 walkbacks, **41,711 with no clip evidence**, 22.06 % of all ends moved,
   and at *RPL32* (exon 1 ends `...AAAA`) 24 of the 33 reads on the exon-1 3' end were walked 4 nt
   off it, taking the splicing-intermediate peak from 33 to 1. Gated, the same run moves 0.56 % of
   ends and that peak reads 33 → 29, while real CPA tails are unaffected (*TEF2* `tail >= 3` is
   57.9 % either way). **`--walkback-unconditional` restores invariant-7 behaviour — use it for a
   poly(A)-SELECTED input, not for nascent RNA.**

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
- **`run-all --short-read` collected only the bbmap arm** — `_ALIGNER_NAMES` (`run/helpers.py`)
  listed no COMPASS aligner, so `_collect_per_aligner_bams` never saw the STAR / HISAT2 /
  magicblast / gsnap BAMs and `corrected_reads.tsv` was bbmap-only while the 5- or 7-way consensus
  sat unused in `multialigned.bam`. Fixed 2026-09-03; the list is now a superset of
  `align_command`'s `--aligners` choices and `tests/test_run_all_wp0_remnants.py` parametrizes over
  `COMPASS_{SE,PE}_ALIGNERS` so a new arm cannot be added without being collectable.
- **`run-all --ONT-cDNA` renamed every sample to `stage1_consensus`** (and `--drs` appended
  `_trimmed`) — the sample id was re-derived from `input_path.stem` at four sites AFTER Step 0 had
  replaced `input_path` with a derived FASTQ. Fixed 2026-09-03: one `_canonical_sample_id()` call
  on the ORIGINAL input, before Step 0.

⚠️ Historical-data caveat, not a live bug: ONT-cDNA outputs produced **before** `gene_id` was
assigned on the gene strand carry corrupt `gene_id` for antisense reads. Regenerate rather than
trust `gene_id` in old cDNA outputs.

## `subprocess.run(capture_output=True)` can deadlock a whole `run-all` for 6 hours

**Status: OPEN (mitigated, not fixed).** Every aligner arm in
`rectify/core/align/multi_aligner.py` is launched as
`subprocess.run(cmd, capture_output=True, text=True, timeout=ALIGNER_TIMEOUT)`
with `ALIGNER_TIMEOUT = 21600` (6 h). If the aligner is a **wrapper script** that
re-spawns itself and the re-spawned process is reparented to init while still
holding the inherited stdout/stderr pipe write ends, `communicate()` never sees
EOF and the run blocks at 0 % CPU until the timeout — i.e. past most cluster
walltimes, so the job is killed with no diagnostic in the log.

Measured on Hoffman2, `rectify run-all --short-read --Scer` (planning 861 §S3),
where the COMPASS indices do not exist inside the package data dir:

```
python -m rectify run-all …            State S, fds 5,6,11,13 -> pipe:[…]
                                       2 task wchans = poll_schedule_timeout
[perl] <defunct>  PPID <python>        rectify's direct hisat2 children, exited
perl …/hisat2 …   PPID 1   0.0 %CPU    re-spawned wrapper, orphaned, holds the pipes
```

The same `hisat2 -x <missing index>` run from a shell exits **2** in under a
second — the deadlock requires the parent to be holding the pipes.

**Mitigation shipped (`fix/runall-truseq-861`):** `_require_compass_index()`
pre-flights the derived STAR / HISAT2 / BLAST / GSNAP index paths in
`run_star` / `run_hisat2` / `run_magicblast` / `run_gsnap` and raises
`FileNotFoundError` *before* the binary is launched, so a missing index becomes
the `DROPPED-ALIGNER` the caller already handles. **This does not fix the
deadlock** — any other fast-failing wrapper-script aligner (corrupt index,
unreadable path, a wrapper that forks for another reason) reaches it. The real
fix is to stop using `capture_output=True` for aligner subprocesses: redirect
stdout/stderr to files, use `start_new_session=True`, and kill the whole process
group on timeout.

## `rectify run-all --chunked-alignment` was a hard `TypeError` on every datatype

**Status: FIXED on `fix/runall-truseq-861`.** `run/chunked_batch.py:221` has
passed `write_per_aligner_corrected_bams=True` since `9d35ae1`, but the
parameter was never added to `split_command.generate_alignment_scripts`, so the
call raised `TypeError: generate_alignment_scripts() got an unexpected keyword
argument 'write_per_aligner_corrected_bams'` and only `00_split.sh` was written.
The call site is unconditional, so this affected DRS, ONT-cDNA, QuantSeq and
TruSeq alike. `9d35ae1` is an ancestor of `429f1c9`, so the path was already
broken when the 830–834 audits were written — planning `[[833]]` G-3's
"generates a long-read pipeline with exit 0" was never reachable.

## `hisat2` accepts a non-existent `--known-splicesite-infile` silently

**Status: OPEN (warning added).** Index present + splice-site file absent →
exit 0, normal alignment rate, no message. A half-built
`HISAT2_annotated_index/` therefore yields a green run with **zero annotated
junctions**. `run_hisat2` now logs a WARNING; it does not fail, because novel
splicing is still on and the arm is usable.

## RECTIFY's bundled S. cerevisiae GTF cannot build any "annotated" COMPASS index

**Status: OPEN.** `rectify/data/genomes/saccharomyces_cerevisiae/
saccharomyces_cerevisiae_R64-5-1_20240529.gtf` has feature types
`mRNA / CDS / gene / intron` and **no `exon` lines**. Consequently
`STAR --sjdbGTFfile` aborts (`Fatal INPUT FILE error, no exon lines in the GTF
file`), and both `hisat2_extract_splice_sites.py` and `gtf_splicesites` emit
**0 sites** — so none of `STAR_annotated_<L>_bp_SJDB_index`,
`HISAT2_annotated_index` or the GSNAP `--use-splicing` map can be produced from
the annotation RECTIFY ships. Workaround: derive exons as *mRNA span minus its
`intron` children* from the GFF (planning `861_scripts/861_make_exon_gtf.py`,
11,562 transcripts / 12,175 exons / 613 introns).
