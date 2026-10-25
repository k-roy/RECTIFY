# RECTIFY Bugs to Fix

## 2026-09-13 audit reconciliation — Codex / Fable

Initial audit LANDED on master by Fable at `d7d92f7`, with independent 517/1-skip
validation. Follow-up `audit/consensus-provenance-20260914` @ `4ec966a` is published,
315 passed / 4 skipped, awaiting Fable integration. Details and review inbox:
`dev/audits/codex_fable_20260913/REVIEW.md`; portable report
`docs/development/2026-09-13-junction-audit.md` on that branch.

Ancestry on master `af2f788` confirms old branch-only narratives below are stale
for A4, A7, A8 Case-B support, A9, C3, D2, D3, D5 and both D8 repairs. A10/A12/A13/A14
plumbing/performance work is already landed too. Keep AT-AC enabled and the 2000
ceiling, as the existing A/B evidence supports.

New repairs: triage original-arm order FN; incomplete micro-exon candidate sets;
junction-local mismatch / reference-compressed SEQ evidence; single-arm unspliced
count doubling and lost fallback signals; primary-only cDNA molecule counting;
paired AT-AC aggregate filtering; pileup runtime error; undeclared PyYAML.

**A11 follow-up / CFX-03 ✅ FIXED 2026-09-16 (`f34cb21`):** B2/B3 now validate EVERY
relocated base (`_whole_block_improves`: the whole block past / before the split must
beat the current placement by `arb_margin`; refusals counted `arb_mm_whole_block_refused`);
the 30→94 witness is refused and the real spliced tail/head controls still splice
(`tests/test_resolver_whole_block_cfx03.py`). **2H realizability ✅ FIXED 2026-09-16
(`ae5cd5b`, supersedes `1597a3b`)** — see the section below. Consensus canonical credit
and micro-exon provenance are fixed on the follow-up branch above. Versioned Xb
records successful calls and credits exact selected Ns; legacy metadata remains
explicitly unverified. Fable has runnable witnesses and integration packets.

Full non-slow 3699 passed, with named existing data/environment exceptions in the
report; targeted exception checks and final triage tests passed. No cohort-level
FN/FP-rate claim. Root WIP and active cluster production sources preserved.


## Last Updated: 2026-09-16 (CFX-03 + 2H per-candidate walk + D11 XO-frame LANDED; ISSUE-047 open) — previously 2026-09-10 evening (A12 ask 2 + A13 ask 1/3 knob + A14 BUILT on `todo/resolver-a12-a13`, not merged) — previously 2026-09-10 (A13 filed: ceiling-abandonment is a CORRECTNESS bug — 5.2-5.4% of clips abandoned, 5' rescues silently fail; A12 filed: overhang_resolver single-threaded + per-contig fixed cost; suspect the 2026-09-05 AT-AC default flip) — previously 2026-09-05 (overnight todo run: A4/A7/A8ii/A9/C3/D2/D5 resolved on `todo/*` branches, NOT merged; new A10/A11; Sumner human-DRS ISSUE-001..011 in flight on `todo/sumner-*` — `dev/todo_run_20260904/README.md`) — previously 2026-08-21 19:40Z (D1/E1/E2 actioned; D7 XA-clobber fixed in tree) — previously 2026-08-20 (reconciled against HANDOFF.md, HANDOFF_REALIGNER_NEXT.md, dev/TRIAGE_PRE_RERUN_20260812.md, Chanfreau planning/689a–722b, and the rbrowse inbox; code anchors re-verified on `95bc039`)

---

## Open — filed 2026-10-07 (Chanfreau 1090u): per-read 5′ truncation correction for read-level panels (FEATURE, Kevin)

**Ask (Kevin, 2026-10-07, verbatim):** "The main issue with (B) is that it shows only coverage, so we lose the 3' cluster grouping that
is critical for showing the IPA. I see that our truncation correction only applies to coverage. We'd need to develop a per-read level
truncation correction that would stochastically "add" that coverage to individual reads."

**Context.** `rectify.visualize.pileup` has the survival correction (`Survival`: S(d) = 2^(−d / t½), weight 1/S capped at 8;
`spliced_distance`; `corrected_coverage`, Horvitz–Thompson; `five_prime_profile`), on master since d7d92f7 and in c287748. It corrects
COVERAGE only. Chanfreau planning/1090u drew IPA manuscript Fig. 2g both ways (PSP2, CQD1; WT-AA ×2, Prp28-AA, Snp1-AA): the
coverage version loses the 3′-sorted molecule piles, i.e. the grouping by 3′-end cluster that shows the intronic poly(A) isoforms,
so the paper keeps the uncorrected molecules (comparison page https://claude.ai/artifact/3rpDZPHA4RcSf1H9Ei2FsX; script
`Chanfreau planning/1090u_fig2g_truncation_versions.py`; numbers `1090_figs/1090u_fig2g_measure.json`).

**Wanted.** A per-read correction that keeps every read's own 3′ end and cluster membership and stochastically extends truncated reads
toward their 5′ end, so a dense pile (`dense_pile`, `read_panel`, `pile_stack`) can be drawn corrected:
- Sample the missing 5′ part per read: the hazard is memoryless, so given an observed read end at y the true 5′ end is drawn from the
  locus's 5′-end distribution upstream of y (template: reads of the same 3′ cluster with `clip5 == 0` that reach an atlas TSS, or the
  frozen cDNA isoform atlas 5′ reference). An extension that crosses an intron must also draw its path (spliced or unspliced) from
  the reads of the same 3′ cluster that span that intron.
- Draw imputed bases in their own style (lighter tint or hatch), never as observed sequence; fixed seed; print the imputed fraction
  per pile and the half-life source in the legend (R6, as `Survival` does).
- Acceptance: (1) imputed coverage reproduces `corrected_coverage` in expectation (test over seeds); (2) truncate full-length reads
  with the fitted hazard and check that the imputed 5′ ends recover the held-out truth; (3) report loci where the library-wide t½
  under-restores: in 1090u, PSP2 Barr WT-AA reached only 0.33 of molecules over exon 1 after correction (0.13 as sequenced), PSP2
  Snp1-AA 0.45 (0.25), CQD1 WT-AA Sept 2025 0.68 (0.42), so a locus- or length-dependent hazard may be needed.
- Suggested home: `pileup.impute_five_prime(reads, survival, tss_model, rng) -> List[Read]` carrying an `imputed5` segment, and an
  `imputed=` style switch on `dense_pile`.

## Open — filed 2026-10-12 (Chanfreau 1090 inbox note, 2026-10-11: cdna-analyze aborted on records past a contig end)

### ISSUE-094 🧭 FILED 2026-10-12 — deSALT aligns across contig boundaries (concatenated index); consensus kept 875 such placements in a July cDNA BAM, and master still keeps the wrong contig
Record: `dev/fixes_20261012/HANDOFF.md` §1 and `oob_origin.json` (same folder; `oob_origin.py`, H2 job 15139832).
- **What.** July Path-A wtaa_rep1 `align2_fullpool/wtaa_rep1_consensus.multialigned.bam` (H2 `shared/processed/cdna_wt_rna15_ysh1_20260711`,
  2026-07-16; consensus from minimap2, gapmm2, deSALT) holds 875 of 8,103,114 records whose reference span runs past the end of their
  contig (invalid SAM), on all 16 nuclear contigs: chrXV 360, chrIV 174, chrV 73, chrXIII 64, chrXII 36, chrI 28, … (1090's "all on
  chrI" came from the first 20 records of a coordinate-sorted BAM). All 875 are primary, `Xa=deSALT`, `Xn=1`, `Xc=low`. The deSALT arm
  BAM holds 1,094 (958 primary, 136 secondary); gapmm2 and minimap2 hold none.
- **Mechanism (verified).** The bases aligned past each contig's end match the NEXT contig in @SQ order at (position − contig
  length): identity 0.990 over 280,190 bases, vs 0.261 at a random shift (chance). deSALT indexes the genome as one concatenated
  sequence and reports a boundary-crossing alignment on the earlier contig. The stage-1 cluster contigs agree: the 346 molecules
  clustered on chrXVI are placed on chrXV, the contig before it. Typical shape: a short spurious 5′ block on the earlier contig, then
  a 100–270 kb N into the next contig (`chrII_cluster_30810`: chrI 31,483–290,962 on a 230,218-nt chrI, `169S30M162433N19M…`).
- **Which step.** deSALT writes the placement; RECTIFY's consensus selection (July: no contig check) chose it as the lone,
  low-confidence winner over the minimap2/gapmm2 placements.
- **Master (`f15b523`).** `intron_sanity` (d0e3a0f, 2026-08-12) soft-clips at the first N that runs past the contig end, so new
  multialigned BAMs no longer overrun through an N (and fa63c12's max intron now stops deSALT making the jump at all; below). Two gaps remain: (a) what it keeps is the spurious fragment on the WRONG contig,
  while the molecule's real placement (another arm's) is discarded; (b) an overrun through M/D ops is only tagged `Xi` and written as
  is, and the checkpoint validator samples 1,000 records per batch, so such a record can still reach a written BAM.
- **Downstream.** fa7f315's cdna-analyze aborted on these records (rc 1). The ISSUE-093 follow-up skips and counts them in
  cdna-analyze and `correct` (intake and writers emit them as stored) via `alignment_tags.span_fits_reference`.
- **Suggested fix.** In consensus selection, reject an arm placement whose span leaves its contig before scoring, so another arm
  wins; reject or truncate M/D overruns too; make the past-contig check total, not sampled. At the arm, post-filter deSALT records
  that cross a contig end (or pad its index).
- **Exposure.** Since fa63c12 (2026-08-17) every arm, deSALT included, runs with the annotation-derived max intron (`-I`), which
  rules out the 100–270 kb jumps: two master-era deSALT arms (1090z42 DRS chunks barrwt_rep1 / snp1aa_rep1 s000, 194,918 and
  89,251 records; H2 job 15142227, `arm_oob_scan.py`) hold none, where July's rate (1,094 of 6,497,653) predicts ~33 and ~15.
  Residual: a crossing within one max intron of a contig end. Old data: multialigned BAMs made before d0e3a0f (2026-08-12) and
  deSALT arm BAMs made before fa63c12 (consensus.py's comment cites planning/684c at 36 per 400k DRS reads; intron_sanity's
  docstring counts 3 of 400,001 running off the end: the two are not reconciled here). Re-correcting such an arm with master's
  correct-first `run-all` fails on the first `=`-encoded one, and `run-all` then continues WITHOUT that arm (a WARNING only).
  1090's scan found them in 1 of 33 July/933 cDNA BAMs. Pre-08-12 DRS BAMs: not scanned.

**EXPOSURE NARROWED 2026-10-25 (Chanfreau `[[1190]]` atlas, 194 inventory rows / 93 admitted
libraries).** The line above — *"Pre-08-12 DRS BAMs: not scanned"* — was the live risk for the
junction atlas, because a boundary-crossing record carries a 100–270 kb `N` and a junction caller
emits that as a novel junction. Checked, and **it does not reach the atlas**, on two independent
grounds:

| check | result |
| --- | --- |
| `n_junctions_out_of_span_only` over the **93 admitted** libraries | **93 of 93 measured, 0 with any** |
| admitted libraries whose `aligner_cmd` involves **deSALT** | **0** |

The second is the structural reason for the first: deSALT is the only arm that produces these (this
entry records gapmm2 and minimap2 holding none), and no admitted atlas library runs it.

🔴 **Two precautions that made this answer trustworthy rather than reassuring:**
1. The naive read was *"0 across all 194 rows"*, which is **wrong** — 102 rows carry a literal `0`
   and **92 carry an empty string**, i.e. were never called. `n_junctions`, `annotated_recall` and
   `n_spliced_reads` are non-empty in exactly the same 102, confirming the 92 are uncalled rows
   rather than a missing measurement. Summing a column that mixes measured zeros with blanks reports
   a clean zero either way. (Chanfreau memory `feedback-max-gate-passes-on-all-zeros`.)
2. All **93 admitted** libraries fall inside the measured 102, so the atlas's own working set is
   fully covered — the blanks are rows it does not use.

**What this does NOT establish:** whether pre-08-12 DRS **BAMs contain such records**. This tests
whether an out-of-span *junction* reaches a result, which is the binding question for the atlas but a
weaker claim than a direct record-level scan of those BAMs. That scan is still owed for any consumer
that reads those BAMs at record level rather than through the junction caller.

Independently, Chanfreau `planning/1189_pressure_test.py` gate 3 (*"placed by minimap2"*) already
discards any junction whose placement was won only by deSALT or uLTRA, so the deck-facing novel
junctions carry a third layer of protection against this class.

## Open — filed 2026-10-10 (Kevin: "dive deeper and sort this out", the cDNA `=` bug)

### ISSUE-093 🧭 FILED 2026-10-10 — cdna-analyze never walked back cDNA 3′ ends (calmd `=` SEQ); #6 then broke cDNA tail lengths
Record: `dev/prod_20261009/CDNA_EQ_BUG.md`. Since `7e408b1` (2026-05-15) cdna-analyze ran the shared walkback on `=`-encoded
records, which can never anchor, so `corrected_3prime` / cluster anchors stayed the raw alignment end: 25.6 % of molecules sit
1–50 nt downstream of the walked-back end (smoke, verified per record). Gene counts and isoforms are unaffected (isoforms use
the raw end by design). Every one of the 276 cdna-analyze outputs on H2 (2026-06-25 … 10-06) carries it. #6 (`83af3d1`) made
the multialigned SEQ literal, which fixes the anchors but lets the walk's residual A count override the carried stage-1 `XA`
(median tail 33 → 5 for 26 % of molecules): master writes wrong cDNA tails. Fix on `fix/cdna-analyze-eq-decode` (`a72721d`
decode + `=`-safe shared walkback; `dd9f681` carried tail on `XK > 0`), NOT pushed: Kevin's call, with the remediation options.
Chanfreau planning/578 found this root cause on 2026-08-24 and fixed only stage 1.

**🔴 EXPOSURE READ 2026-10-25 — the IPA work is the most exposed consumer, and the bug's worst class
is exactly its class.** No new measurement; this is a connection across §4, §10 and §3 that the entry
had not drawn.

§10 splits the unwalked rate by **whether the walked-back end lies in an annotated intron on the
molecule's strand**, and the intronic subset is the worse one in two of three libraries:

| library | other 3′ ends unwalked | **intronic** 3′ ends unwalked |
| --- | ---: | ---: |
| 1076 tun_poly_bc03 | 40.3 % | **56.6 %** |
| 1088 rrp6aa_rep2 | 28.7 % | 23.4 % |
| 1082 wt_30min_rep1 | 27.8 % | **32.8 %** |

§4's consumer list includes **907c/d/e (intronic and internal PAS)**, **1075/1076 (IPA QC, polysome
IPA)** and the **1090 manuscript drafts (IPA)**. Intronic polyadenylation is a claim about **where in
an intron a 3′ end sits**, at nucleotide resolution — precisely what §3 says the bug changes, and
nothing else in the consumer list depends on the intronic subset specifically.

Two properties make it worth a decision rather than a shrug:
- **The shift is systematic, not noise.** ~75 % of moved ends move 1–2 nt, ~16 % 3–4 nt, ~8 % 5–10 nt,
  0.5–1 % more than 10 nt — and **never downstream**. A uniformly upstream shift can move a called
  site relative to a hexamer or relative to an intron boundary; random error would not.
- **Gene-level and isoform-level results stand** (§3; isoforms use the raw end by design), so this is
  narrow rather than sweeping. 1095b/1077's length correction is gene-level and unaffected.

**What would settle it cheaply:** §7 option 2 — re-run **only** `cdna-analyze` on the existing
multialigned BAMs (~1.5 min on the 79 k-molecule smoke) on a `c287748` + fix pin, which keeps every
consensus pick and changes only the anchors. Then compare the IPA site calls before and after. That
is a far smaller job than re-deriving anything, and it answers whether the IPA conclusions move at
all rather than leaving it as a known unknown in a manuscript draft.

**Not exposed:** the `[[1190]]` junction atlas and the Mex67 splicing deck. Checked — no document of
that campaign consumes `corrected_3prime`, `cdna-analyze` output or `tail_len`; it is junction work,
and junctions come from `N` ops rather than 3′ anchors. §9 also records rbrowse exposure as none.

## Open — filed 2026-10-09 (Chanfreau 1099 inbox note of 2026-10-07)

### ISSUE-092 🧭 FILED 2026-10-09 — the DRS panel's consensus pick is not run-to-run deterministic on repeat-derived reads
= GitHub #8 (filed 2026-10-09 by the queue-audit session: 46 of 159,051 primary records, all inside R64_selfhomology.bed). Same chunk (rrp6aa_rep2 c000_chrI_1, 159,051 reads), same pin `c287748`, same inputs and `1066d_drs_chain_v2.sh` logic: two
identical 4-thread runs differ on 50 reads (0.031 %); 8 vs 4 threads differ on 33 and 48, so thread count is not the cause.
Every differing read with a gene_id lies in duplicated sequence (Ty1/Ty2 copies, subtelomeric YAR060C, FLO1); the rest are
intergenic, LTR-like. 9 of the 50 had per-aligner MAPQ > 20 in one run, and the consensus picks a different copy (a different
chromosome). So the tie-break among near-identical copies varies between runs. Kevin's #6 smoke (83af3d1 message) saw one such
read in 3,000. Paths: H2 `/u/project/guillom/kevinroy/1099/{drs,drs_t4test,drs_t4test2}/rrp6aa_rep2/c000_chrI_1/`; record
`~/work/UCLA/Chanfreau_Lab/planning/1099_anchoraway_drs_rectify_c287748.md`. Consequence for any A-vs-B smoke: ~0.03 % of reads
differ with no code change (≈ 13 of the 42,627-read route smoke), so a diff at that level on repeat copies is noise. Fix
direction: find the unordered step (dict/set order across aligner arms, or a worker-order merge) and break ties by a stable key.

## Open — filed 2026-10-03 (independent audit of the consolidated queue, `dev/prod_20260928/audit_20261002/`)

### ISSUE-087 addendum (2026-10-03) — the loss class is error-dense WRITTEN tails; the chance test cannot see them
Two independent reviewers (verified in `audit_20261002/AUDIT.md`): `chance_3prime_ext.py` fits the clip only within clip + 12 nt
(`EXTRA = 12`), so its p never scores the tail the landing wrote (P3-3: written span 34 nt; P3-20 passes the test while its
written tail scores 15.0 bits, under the null q99 of 21.2). Tails carrying `CCATC` = 35 (composition-preserving shuffle
expectation 7.5), pin end better on 33 of 35; 20 of them at chrI + 143,367. 646a397f: a 16-nt poly(A)-like tail written as
`15=/5875D/1S` moves the corrected end +5,878 nt; 7 of 925 extended tails carry D ≥ 20 or move the end ≥ 100 nt. Corrected count:
925 alignments extend; the corrected 3′ end moves downstream for 851, upstream 15, not at all 59. **Gate options costed on the
clean-read metric (2026-10-03, `audit_20261002/AUDIT.md`; reads fixed − broken of 925):** the filed gate (≥ 16 nt, chance-fail on
both statistics) +25, catches 32/35 CCATC tails, reverts P3-3; written-tail score < null q99 restricted to ≥ 16 nt +27, catches
33/35, reverts the positive control P3-20; written tail at ANY length −202 (rejected); error-dense written tail ((X+I+D)/len > 0.3,
≥ 8 nt) −47 (rejected); total tail D ≥ 20 −4 (rejected: 5 of its 6 reads are landing-better). Add to either surviving gate: a
deletion ≥ 100 nt inside a tail never beats a soft clip (only 646a397f; fixes 1, breaks 0). Kevin picks (queue v17 D1).

### ISSUE-088 🧭 FILED 2026-10-03 — TSV `junctions` keeps a junction the 3′ walkback removed (CF-1, 5cac4ae6)
After the walkback drops a terminal exon and its intron, the BAM record (ISSUE-085-corrected start 40,834,509) carries one N,
while the TSV row still lists `40832463-40834509`. The writer audit tags it (`Xh: jx-:…`) instead of making the TSV match.
read-review §5 ("TSV junctions must equal the written N ops; the BAM is the record") already settles it: fix, not discuss
(both reviewers, P11). 3 of 62,602 chr5 cohort reads, all walkback reads. Open question for Kevin (B): should the walkback
remove a whole terminal exon and its intron at all?

### ISSUE-089 🧭 FILED 2026-10-03 — the 5′ rescue destination guard accepts canonical motifs from EITHER strand
`bam_writer._revert_selfinflicted_noncanonical_n` judges a writer-drawn N with
`overhang_informativeness.is_canonical_junction`, which accepts all six forward-genome dinucleotide pairs (GT-AG, GC-AG, AT-AC
and their reverse complements CT-AC, CT-GC, GT-AT). For stranded RNA (DRS, RNA-sense cDNA) only the read's own strand's motifs
are canonical, so a wrong-strand pair passes: Reviewer B found novel CT-GC (RS-28, a skipped tandem-repeat unit) and CT-AC (RS-29,
a chimera) drawn on + strand DRS reads in the re-split census. **Not a landing regression:** master `c287748` has the identical
guard (checked 2026-10-03). Fix: pass the RNA strand and require that strand's motif when the protocol is stranded.

### ISSUE-090 🧭 FILED 2026-10-03 — 2F writes Case-4 5′ snaps under the 12-bit evidence floor (master too)
Queue v17 ruling R1 (083-26…29, chr5 cohort): the written 5′ block scores −1.5, −0.5, 1.5 and 4.5 bits on 2–4 matched bases,
under the floor every landing must meet (read-review §5). 083-29 creates a NOVEL 56-nt AT-AC junction inside an annotated RPS23
exon; 083-26's clip fits the annotated RPS14 exon-1 end (150,449,702) 8/9, so the snap took the wrong donor. Master
`base_3957b6c` writes the same four records (same junction, same bits; Sherlock census
`ultracode_census_20260921_b94b9a1/out/base_3957b6c/cohort/records.tsv`), so it is NOT a landing regression; the landing only
slides each insertion one base off the N. Both reviewers rejected "draw the tip" on all four. Fix direction: no 5′ write under
the floor; never a novel junction inside an annotated exon from a few matched bases. Kevin's ruling: queue v17 R1.

### ISSUE-091 🧭 FILED 2026-10-03 — a REFUSED 5′ rescue keeps the re-anchored longer clip, dropping bases the aligner placed (master too)
Both reviewers (re-split census cards RS-2, 8, 11, 14, 15, 16, 18, 25, 27, 33; Reviewer A: 20 of the 342 census reads start 4–57 nt
short of the aligner's placement). Master writes identical records for RS-2/8/11/15/16/25 and the same 5′ clip for RS-14/18/27/33;
across all 347 census reads the landing never clips more at the 5′ end than master (checked 2026-10-03). Old behaviour, not a
landing regression. Fix direction: when the rescue is refused, restore the aligner's own 5′ placement.

Moved here from the queue on 2026-10-03 (both reviewers agreed with each; no pill needed): P10 = ISSUE-086 report-only flag
(N from a shuffled-SSP null; Kevin picks N), and P11's carried items (A13 ask 2 max_intron 2×, D10 XP vs XA, the ISSUE-048
open list, the deliverable program's multi-condition run).

## Open — filed 2026-09-28 (production-route smoke of the landing, `dev/prod_20260928/`)

### ISSUE-087 🧭 FILED 2026-09-28 — DRS consensus: a forced 3′ tail of ≥ 16 nt that fits the genome no better than chance beats the soft clip; a forced tail can lift a MAPQ-0 multimapper to 60
**Status 2026-10-09 23:15 PDT:** the landing is on master as `f15b523` (D1); this gate is D1's condition for any Chanfreau
reprocessing at the new sha. Re-measured at 83af3d1 → f15b523: long chance-fail 21 landing / 46 pin / 2 tie, unchanged from P1.
`score_alignment` (`rectify/core/consensus/scoring.py`) charges an RNA-3′ soft clip −2/nt (cap 10) and charges nothing for
an aligner that aligns the same bases with mismatches or deletions, because `_get_effective_3prime_clip` counts only
non-A mismatches/insertions in 8-nt windows at ≥ 40 %. Master hid this by mis-scoring uLTRA's `=`-encoded SEQ; the
landing's decode (ISSUE-066) exposes it. Route smoke (wt_by4742_rep1, chrI, 30 %, H2 `/u/scratch/k/kevinroy/prod_20260928/
route_smoke`): 925 reads get a longer 3′ alignment in the landing. Against the library's clean reads (37,392: same pick in
both versions, no 3′ clip) the landing's corrected end is better supported (clean ends within 2 nt) in 619, the pin's in
203 (103 ties); outside genomic A runs 441 to 189. Tails ≤ 15 nt are mostly real (recurrent at their sites; per-read
chance power is too low to show it). The losing class: tails ≥ 16 nt failing the per-read chance test
(`dev/prod_20260928/chance_3prime_ext.py`, p ≥ 0.01 on the whole-tail and the end fit): 69 reads, pin better 46, landing
21. Candidate fix: treat such a tail as soft clip for scoring (net +25 of 925 here); per-length bits thresholds from the
model's null (median q99 ≈ 16.5 bits at 16–20 nt, 21 at 21–30, 26.5 at 31–50). **Rejected:** gating every tail by the
chance test (fixes 198, breaks 574). Separately 2 of 216 MAPQ-0 reads become MAPQ 60 (220b3a40: identical body on chrI
and chrXII, 128-nt tail forced onto chrI at NM 74). Review cards: consolidated queue section P3. **Kevin's call:** the
gate and its threshold. Not a blocker for landing (the flux pin is unaffected until re-pinned).

## Open — filed 2026-09-28 (Chanfreau 1021 inbox note of 2026-09-24, Kevin's review card C09)

### ISSUE-086 🧭 FILED 2026-09-28 — cDNA Path A: SSP strand-invasion molecules pile their 5′ ends at genomic SSP matches (CDC48)
Kevin asked why CDC48 cDNA reads "have a strange cut off half way through the CDS". Measured by chanfreau-1021 on H2
(`~/work/UCLA/Chanfreau_Lab/planning/1036_rbrowse_correction_pills.md`, "Orchestrator findings on C09 and C12/D08"):
Path A consensus molecules pile their 5′ ends at chrIV:237,413 inside the CDC48 CDS (YDL126C, − strand) — ysh1aa_rep1
233 reads at 237,405–237,425 (182 exactly at 237,413), wtaa_rep2 ~225 — and direct RNA of the same strains has no pile.
225/233 carry a ≥ 20-nt 5′ soft clip (29 ≥ 100 nt) of adapter-like sequence with UMI-like variable 5′ parts; only 2 have
an SA tag; all Type 1 (`XT:i:1`). The genome there reads GGTGCTGATATTGC (RNA sense) = SSP positions 9–22
(`rectify/core/cdna/_constants.py`); the exact 14-mer occurs twice in R64-5-1 (chrIV:237,400 and chrXIII:446,405). This
is the template-switching-oligo strand-invasion artifact: the SSP anneals internally, RT jumps onto it, and the
molecule's aligned 5′ end is the motif. **Build (safe, report-only):** flag molecules whose aligned 5′ end sits at a
genomic match to the SSP's 3′ region starting at the 5′ end with ≥ N identical nt (a BAM tag + a Path A stats count
with N reported); a genome-wide scan of SSP-like k-mers against cDNA 5′-end piles sizes it (shorter partial matches are
common and may make weaker piles). **Kevin's call:** dropping or down-weighting flagged molecules (result-changing for
every cDNA 5′-end analysis) and N. Parked behind the 2026-09-28 production landing.

## Open — filed 2026-09-23 (Chanfreau 1021, Kevin: "I think we need to add this cDNA correction approach to the Rectify package")

### ISSUE-048 🧭 FILED 2026-09-23 — FEATURE: per-library length-bias correction for ONT PCR-cDNA gene counts, calibrated on a stable-gene panel against direct RNA
Replicate PCR-cDNA libraries count a 3-kb gene 0.16–12× as often as a 0.5-kb gene relative to their replicates (direct RNA 0.82–1.14×).
Chanfreau 1021x fits, per library, a natural spline in log10 transcript length on 405 genes whose direct RNA abundance is stable across
11 DRS genotype sets, and divides it out. On held-out genes: replicate slope median |1.20| → 0.10 log2/decade (DRS 0.05); replicate SD of
genes > 2 kb 0.75 → 0.26 (DRS 0.18); Rrp6-AA ÷ WT-AA fold changes vs DRS r −0.16 → 0.47. Limits to carry: within-gene / short-molecule
readouts need an explicit scale c (≈ 0.5 on the PSP2 isoform ratio; unresolved below 314 nt); values are on the DRS scale; the panel is
certified only for anchor-away W303 in YPD. Full spec, reference code and acceptance numbers:
`dev/specs/FEAT_cdna_length_bias_correction_1021.md`. Result-changing defaults are Kevin's call.
**✅ LANDED 2026-09-23 on master `3d55963` (fast-forward; origin + h2 mirror; Sherlock tree not updated)** — branch `feat/cdna-length-correction`
(`.worktrees/cdna-length-correction`), 8 commits `05391b6`…`3d55963` on `3957b6c`: `rectify/core/cdna/length_bias.py`, command
`rectify cdna-length-correct {genes,count,panel,fit,apply,ratio}`, bundled Scer panel/reference/gene lengths with PROVENANCE, docs page +
README figure. Kevin's rulings: that command name with cdna-analyze unchanged; bundle all three; c has NO default; the cohort tables stay
out of git (`RECTIFY_CDNA_LENGTH_FIXTURES`). Verified on H2: reproduces 1021c/1021v/1021x exactly (coefficients < 1e-9); 34/34 targeted
tests incl. the cohort regression; 12/12 planted defects caught; full suite failure set identical to base. Open for Kevin: keep/drop the
report-only h(L) cohort line and the opt-in `--reference-self` / `--panel-high-count`; scope tokens; a gene table restoring the 8
reporter-host loci; ISSUE-049 (DRS truncation fit/apply/deconvolve) vs a 716e redraw. The PSP2 isoform scale was refit after the D11 recount to about 0.6 (0.54–0.60
leave-one-out; wording only, `3d55963`). Handoff:
`handoffs/HANDOFF_issue048_cdna_length_correction.md`; work files (fixture, gates, figure scripts): `dev/issue048/`.

## Open — filed 2026-09-10 (Chanfreau 907d full-depth run)

### ISSUE-047 🧭 FILED 2026-09-14 — acceptor-only 2H moves are realized as `N` + glued `kD`, never a block shift; with 031 + realizability the cryptic-acceptor discovery path is dead at surgery
Witness: `tests/test_microhom_drift_guard.py` PG fixture (`40M90N16M`, exon-2 from +6): surgery says "adds indel 0→6, refusing" on every
tree since 031; honest CIGAR `40M96N16M`. Two tests xfail(strict). Fix = block-shift realization gated on an anchored positional
score + a review class; result-changing → Kevin. Full text: `dev/sumner_misplaced_panel_20260904/ISSUES/ISSUE-047_*.md`.

### realizability-before-ranking ✅ LANDED 2026-09-16 (`ae5cd5b`; 1597a3b's claim repaired per Codex's HOLD)

Every decision is now PER CANDIDATE inside the walk: incumbent stops it; the move gates
(margins/cap/positional spare, annotated-canonical hold) run per candidate BEFORE the dry
run (`gate_vetoed_candidate_skipped`); the surgery dry run skips what it cannot write
(`unrealizable_winner_skipped`); N-ops are walked RIGHT-TO-LEFT on an evolving trial copy so
each is judged on the read the writer will see (Codex CFX-12 (a)+(b) both pinned in
`tests/test_2h_realizable_ranking.py`; joint result `21M100N1M100N23M`). The 19 policy
tests are explicitly policy-only (`policy_only_surgery` fixture, documented) and every file
carries writable controls with final-CIGAR assertions (pure slides with identical flipping
bases; the D/N twin where annotation permits). The two "inherited" failures were fixture
conditions: `test_correct_no_bam_warning` now straddles `min_informative_clip_bp()` (1, 8
stay clipped; 12, 57 rescue) and `test_restore_polya_from_parquet` skips when the gitignored
script/metadata is absent. Dry-run cost is timed separately (`realizability_dry_run`) and is
paid only by N-ops with a non-incumbent head. Codex's efficiency note (a writable winner is
surgered twice in refine + once in the writer) stands as a known, bounded cost.

### D11. ✅ FIXED 2026-09-16 (`ba856b5`) — `cdna-analyze` took the poly(A) side from the carried `XO`, which names the PRE-alignment frame — filed by Chanfreau 936 (RPL20B reporter)
On an `XN:i:1` record the flag is the frame (flag 0 → 'fwd', 16 → 'rev'); disagreements counted +
logged (`XO frame taken from the alignment flag on N oriented molecules`); pre-XN records keep XO.
Measured before: 43 % of reverse-aligned RPL20B molecules, 13.6 % of forward-aligned construct ones,
concentrated on 5'-truncated molecules; every wrong `corrected_3prime` sat at the 5' boundary.
Chanfreau must re-run `cdna-analyze` (not the consensus) on 936 bp_{1,2,3} and 907d.

### ISSUE-031 ✅ LANDED 2026-09-13 (`af2f788` on master) — no I/D adjacent to an N; HP-insertion exemption per Kevin's 09-09 pills
19/19 cards reviewed by Kevin 2026-09-09 (9 approve / 3 disapprove / 7 discuss; the disapproves are the HP exemption, a
third placement, and a micro-exon ask — none object to the rule). `todo/issue-018` now carries only the preponderance gate.
Environment correction from Codex: `test_cluster_detect.py::test_detect_m1_on_darwin`
failed because the sandbox denied `sysctl`; the real query and targeted test passed
outside the sandbox. Order dependence has not been established by that evidence.


### D10. ✅ BUILT 2026-09-12 (on master, `7e7d691`) — dorado's `pt:i` rides Path A into the consensus as `XP:f` (MEAN over reads with pt > 0) + `XD:i` (n) + `XW:f` (sample SD, n ≥ 2); stage 1 reports within-molecule pt agreement (SD median/p90, CV median); `XA` documented as the SEQUENCE-level A-count (per-cluster median) — filed by rbrowse-50 from Kevin ("I thought Dorado puts out pt tags for cDNA")

Measured by rbrowse on WW1: uBAM `pt` on 5,000/5,000 reads (91 % positive, mean 40 nt); consensus carried only `XA`
(`pt` median 41.5 vs `XA` 25.5, `XA − pt` median −8.5 nt, IQR −37..0); the aligned consensus has NO tail in SEQ by
construction (pretrim strips it). Two drops fixed: `extract_read_info` never read `pt`, and the Step-0 pre-trim rewrote
the FASTQ header as UUID + ro + pl, discarding `samtools fastq -T pt`'s comment before `minimap2 -y`. Now `ReadInfo.pt`
→ `cluster_pt_summary` → `XP`/`XD` on the FASTQ comment, in `_CDNA_COMMENT_TAGS` and the CMA whitelist; the pre-trim
keeps well-formed SAM tags. **The reader is `read.get_tag('pt')` on the input BAM record**, so a Stage-A BAM that
carries `pt` as a real aux field (Chanfreau 907: `samtools fastq -T pt,pa,ns,ts | minimap2 -y` upstream — pt on
2000/2000 records) is covered directly; the `-T pt` requirement applies to `run-all --ONT-cDNA`'s FASTQ input only.
A pt-less input is LOUD (WARNING + `dorado pt carried -> XP/XD 0` in the pretrim-health block; `pt_reads` /
`pt_clusters` in stats) — never a silent `XD:i:0`. 907 status: their 03eaead stage-1 rerun is DONE (WW1 05:58–06:04Z,
15 libraries + figures derived); XP/XD needs its OWN correct-cdna pass (~5 min/library at 1M reads, Stage-A BAMs on
disk) — timing is Kevin's call. NOT changed:
`corrected_reads.tsv` `polya_length`/`polya_source` on Path A still derive from `XA` — switching them to `XP` is a
🧭DECISION (which number is "the" tail for a cDNA molecule). `tests/test_cdna_pt_tag.py` (13).

### D9. ✅ FIXED 2026-09-11 (landed on master the same night) — cDNA stage 1 classified full-length molecules with an error-bearing SSP as Type 2 (no UMI, no dedup, adapter left as a 5′ soft clip) — filed by rbrowse-psp2 from Kevin's r120_4585 / WT-AA_WT_rep2

`read_info.extract_read_info` matched the 23-nt SSP EXACTLY (`seq.find`) while `detect_full_length_tier` and the
stage-1 trimmer's `_find_ssp` were already edlib-tolerant (planning/681). One miscalled base (r120_4585: the last T
read as C) → Type 2, `XU=''`, never merged, 143-nt adapter kept as a 5′ clip. rbrowse census WW2: 79/502 construct
reads (16 %) in this class, 67 with an intact SSP+UMI. **Fix:** one shared finder `read_info.find_ssp_span` (exact
first, unwindowed; fuzzy fallback ≤ `SSP_MAX_EDIT` window-gated to the adapter end; one span rule for every read),
used by Type-1 detection, the trimmer and `_detect_frame`. Repro on the reported clip: Type 1, UMI
`TTAAGATTCGCATTGCCCTTGGGATTC`, `trim_5p` 143. `tests/test_cdna_ssp_fuzzy.py` (13). **Consequences to re-measure
(Chanfreau 907c/907d): the Type-2 share (~17 %) and the XF ≥ 1 denominator; the A13 "35 of 35 intron-start reads
with a median 147-nt 5′ clip" is partly THIS class, not resolver abandonment.** Every Path-A cDNA output written
before `03eaead` carries the old classification — re-run stage 1, do not re-score.

### A17. 🔴 minimap2 splice mode aligns an mRNA 3′ end + poly(A) tail as a JUNCTION onto a downstream genomic A-run — filed 2026-09-11 (Chanfreau 907, `planning/922_junctions_beyond_3p.md` addendum)

RPS26B (chrV): donor 32 nt inside the annotated mRNA end; "acceptor exon" = `TATTGGTTTTC` + A₁₁ at 425,041–425,051;
~half of the junction reads end exactly at 425,052; exon-2 alignment 10–29 nt; **93 % of 57,393 spanning molecules**
in plasmid-free WT-AA cDNA carry it. It survives a ≥ 20 nt double-anchor gate because the A-run IS the anchor. Not
a resolver defect — the aligner's own record. Fix shape: in `rectify align` (or station C / the 3′-end walkback), veto
an N op whose downstream exon is < ~30 nt AND ends in a genomic A-run (or whose read ends at the run) — the same
"tail fragment aligns to genome" class the pA walkback already handles at the 3′ end, now presented as a junction.
Needs a read-level review class (10–20 RPS26B reads on the queue) before any gate: the veto must not eat a real
short terminal exon. Owner: rectify. Tag `📄DOC` until the class is built.

### A13. 🔴🔴 CORRECTNESS, not just speed: the resolver's candidate ceiling SILENTLY ABANDONS 5.2–5.4 % of clips, so 5′ junction rescues fail — measured 2026-09-10

**STATUS 2026-09-10 (rectify, branch `todo/resolver-a12-a13`, worktree `~/work/rectify_worktrees/devq-resolver`):**
ask 1 ✅ BUILT — `ResolverStats` now carries `blowup_by_contig` + `blowup_first` (ceiling, n_cand, window
of the first refusal per contig) and `as_dict()` adds `abandoned_frac`; all three land in
`<prefix>.overhang_resolver.stats.json`. The end-of-run WARNING names the consequence ("N junction rescues
DID NOT OCCUR — read ends pile up at a splice site instead of splicing across it") and the per-contig split.
ask 3 ✅ KNOB — `--resolver-candidate-ceiling N` on `align` and `run-all` (`ResolverConfig.max_candidates_per_clip`,
default unchanged at 2000). **A/B DONE (Chanfreau 907, 2026-09-11, branch @1f01064, `-t 8`, 37,219-read chunk):** ceiling 2000 → 7,732 junction
reads, abandoned 4.81 % (3,268/67,960), 546 s; ceiling 20000 → 7,732 junction reads, abandoned 0, 812 s; the two BAM
bodies are BYTE-IDENTICAL (md5), XE moves 857 vs 857. ⇒ on yeast cDNA the abandoned clips are NOT where the rescues
are; **default stays 2000**; the warning text was softened to "not assessed" (`db7c796`+, next landing). Production
907d runs at 2000: 100k-read chunks in ~24 min at `-t 8` vs TIMEOUT at 1:25 on 78df859. ask 2
(`max_intron` per organism) is unchanged: note `derive_max_intron` ALREADY derives it from the annotation as
2× the longest intron (yeast → 5000, closed-this-month table); the question is the 2× factor, Kevin's call.
ask 4 → A12. Tests: `tests/test_resolver_parallel.py::TestAbandonmentIsAMetric`.


**This is A12's ceiling seen from the science side, and it is the more serious half.**

`overhang_resolver` exists to re-place 5′-soft-clipped reads whose junction failed to align
(Chanfreau invariant #4 — those reads otherwise manufacture a phantom 5′-end pile at splice sites).
When clip enumeration exceeds the ceiling the log says *"was ABANDONED (the read passes through
untouched)"* — i.e. **the rescue the module exists to perform silently does not happen.**

**Measured across nine PSP2 cDNA libraries (rectify `78df859`, yeast R64 + one plasmid contig):**

```
resolver: 10531 of 200214 clip(s) ABANDONED on the candidate ceiling (2000 per 5000 bp)
          11085 of 221083 · 9170 of 180155 · 6738 of 136967 · 6666 of 127505
           6302 of 114057 ·  6257 of 108813 · 6185 of 128612 · 2748 of  50579
```
⇒ **5.2–5.4 % of all clips abandoned, in every library**, and `candidate blow-up` fires on
**every chromosome** (chrI…chrXVI) in all nine runs — not a repetitive-contig edge case.

**Consequence at a locus we looked at by eye** (endogenous PSP2 5′UTR intron, chrXIII:236,592–236,953):
35 of 35 reads that start *inside* the intron carry a 5′ soft clip with **median 147 nt** at Stage A —
the textbook failed-junction signature — and only 1 read in the window is spliced. A reader looking
at the browser sees reads apparently initiating mid-intron; those are mis-placed spliced reads.

**A/B on a fixed 37,219-read chunk (same input, one variable at a time):**

| arm | wall | junction reads |
|---|---|---|
| **default** (AT-AC on, `max_intron` 5000) | **3016 s** | 7766 |
| `--no-resolver-atac` (5000) | 1590 s | 7765 |
| **AT-AC on + `--max-intron 2500`** | **1566 s** | 7724 |

⇒ **The AT-AC pass costs 47 % of resolver runtime and returned +1 junction read here.** That is the
apportionment A12 asked for — reported as a *cost measurement only*: AT-AC stays ON (Kevin
2026-09-10, real yeast AT-AC introns, kept for discovery), so this says *optimise that pass*, not
disable it — ideally by sharing candidate enumeration with the canonical pass rather than repeating it.

⇒ **`--max-intron 2500` is a 1.93× speedup for a 0.54 % junction-call loss with AT-AC retained.** No
S. cerevisiae intron approaches 5 kb (longest ≈ 2.5 kb), so the default window is ~2× wider than the
biology and the candidate count grows with it (`near_sites × far_sites`).

**Asks:**
1. **Treat ceiling-abandonment as a correctness metric, not a log line.** Emit abandoned-clip count
   and rate in the run's provenance JSON, and make a high rate a WARNING that names the consequence
   ("N 5′ rescues did not occur"). Today it is one INFO line among thousands.
2. **Derive `max_intron` per organism / from `--annotation` when supplied** (yeast ≈ 2,500). 🧭DECISION:
   result-changing, so it needs Kevin — but note it likely *improves* recall by keeping clips under
   the ceiling, on top of the 1.93× speedup.
3. Raise the ceiling, or make it adaptive, when the numba kernel is available
   (`RECTIFY_HP_ED_NUMBA`) — the log itself names this as lever (2).
4. Parallelise by contig (A12 ask 2) — still the structural fix.

**Filed by:** Chanfreau `chanfreau-lab-65`, from `planning/907c` + `922`.

### A12. 🔴 `overhang_resolver` is SINGLE-THREADED and its cost is dominated by a PER-CONTIG fixed term — it makes read-chunking nearly useless and times out real jobs (measured 2026-09-10)

**STATUS 2026-09-10 (same branch):** ask 2 ✅ BUILT — `run_overhang_resolver(threads=N)` now scores the stream
on a process pool of N workers (batches of 256 records, ONE writer keeps input order, back-pressure at 3N
batches in flight, dead-worker detection per the D8 lesson). Output is BYTE-IDENTICAL to `threads=1` and the
tally merges field-for-field (`tests/test_resolver_parallel.py::TestParallelDriver`). Start method follows
`RECTIFY_BAM_MP_START_METHOD` (default spawn — each worker loads the genome; on a mammalian genome set
`fork` on Linux to share it copy-on-write). ask 4 ✅ the SINGLE-THREADED warning is gone (it was true; it is
not any more). ask 1 (AT-AC pass apportionment) — measured by the filer (47 % of resolver wall on one chunk);
sharing enumeration between the passes is NOT possible as such (different site classes), the lever is the
per-clip DP cost — see the profile note below once recorded. ask 3 → A13 (Kevin). 🔴 The "1.8 min FIXED
per-contig" fit is from three points with different contig progress; NOT confirmed — a cProfile of the
serial resolver on 7,022 real yeast DRS reads is in this session's record (`prof_resolver.py`).
**PROFILE RESULT (M1, load ~15, cProfile-inflated):** 7,022 reads / 7,221 clips / **0 resolved** / 154 refused
low-info / 1 blow-up → 2,053 s = **0.28 s per clip**. 15.4 M `hp_edit_distance_bounded` calls = **2,133 candidates per
clip** (both passes) at ~115 µs each = 87 % of wall; index build+save 20 s once; `sites_in` 31 s; `is_repeat_expansion`
21 s. **There is NO per-contig fixed term** — the cost is clips × candidates × one DP, and on this DRS locus every one of
the 7,221 clips was enumerated to ~2,000 candidates and none placed: the bill is paid on clips that can never resolve
and that `assess_overhang` does not refuse. Census of that BAM (`550f_tracks/WT_rep1.YAL003W.clean.bam`): 7,009 LEFT
clips on plus-strand reads, median 138 nt, p90 146; 212 right clips median 13 nt; and the alignments carry 30–106 kb
N-ops with 145-nt insertions — it is the 550f impossible-intron ARTIFACT track (C2 class), so "0 resolved" is the
expected answer for junk, not an accuracy finding, and the 7,043 `no_candidates` clips each cost ~2,000 DPs to reach
that answer. Levers, in order: parallelism (DONE); a cheap exact refusal for clips whose candidate set cannot contain
a placement (none exists under the 20 % edit rule — the cutoff is too loose for q-gram filtering; a seeded
"best-first" order is exact but only helps clips that DO resolve); the ceiling/window policy (A13, Kevin). The AT-AC pass cannot share enumeration (different site classes); it
is simply a second ~1,000-candidate pass, hence the filer's 47 %.


**Anchor:** `core/align/overhang_resolver.py` (the post-pass run by default via
`junction_aligners`); the warning it emits itself: *"threads=N requested but the resolver is
SINGLE-THREADED (parameter accepted for API compatibility, not implemented). Size jobs for 1 core on
this stage."*

**Measured on Sherlock (Chanfreau PSP2 cDNA, yeast R64 + a 9.5 kb plasmid contig, 18 contigs,
`--aligners minimap2`, `-t 4`/`-t 8`, rectify @ `78df859`):**

| input | minimap2 | overhang_resolver | contig progress |
|---|---|---|---|
| 2,836,207 reads, 8 threads | **483 s** | **>43 min, not finished** (cancelled) | — |
| 500,000 reads, 4 threads | **95.5 s** | **>48 min, not finished** | reached chrII in 22 min ⇒ **~11 min/contig ⇒ ~176 min projected** |
| **37,219 reads**, 4 threads | ~10 s | **>25 min, not finished** | reached chrX in 25 min ⇒ **~2.5 min/contig ⇒ ~40 min projected** |

🔑 **A 13.4× smaller input ran only ~4.4× faster.** Fitting `t_per_contig ≈ a + b·N` gives
**a ≈ 1.8 min of FIXED per-contig cost** and `b ≈ 1.8e-5 min/read`. Consequences:

1. **minimap2 is now ~1 % of `rectify align` wall time**; the resolver is ~99 %.
2. **Read-chunking does not amortise** — splitting a library into *k* chunks multiplies the fixed
   term by *k*. Splitting 98.8 M reads into 205 chunks was **slower in aggregate** than not
   splitting, while still exceeding a 1:30 task limit per chunk (151 of 205 chunks would have
   TIMED OUT; the array was cancelled).
3. The irreducible read-term alone is **~485 CPU-hours for 98.8 M reads** on this genome.

**Suspected COST contributor (not a defect): AT-AC enumeration was flipped ON by default on
2026-09-05** (`--resolver-atac`, `default=True`). Its own docstring calls it *"a second
clip-enumeration pass"*, i.e. a ~2× candidate cost, and it landed **five days before** these
measurements; resolver wall time does not appear to have been re-benchmarked after the flip.

🔴 **To be explicit, because an earlier draft of this ticket got it wrong: the AT-AC default is
CORRECT and must stay on.** Kevin, 2026-09-10: *"This motif combination has been shown to operate on
several yeast introns so we activated it for discovery potential."* AT-AC introns **do** occur in
S. cerevisiae. **Do not make this opt-in for yeast and do not disable it to buy speed** — the
discovery capability is the point of the flip. (Related: `docs/todo-minor-intron-grammar-825j`, the
minor-intron/U12 grammar request — AT-AC is currently unrepresentable in the splice-site index and
scored worst-tier, which is its own open item.)

⇒ The ask here is **make AT-AC affordable, not optional**: a timed `--no-resolver-atac` vs default
A/B on one fixed chunk (running now as Sherlock `42767198`) is diagnostic only — it measures how much
of the cost the second pass carries, so the optimisation can be aimed at the right pass. Whatever it
shows, the fix is parallelisation and enumeration efficiency (asks 2–4), not switching it off.

**Second contributor — the candidate ceiling is being hit constantly.** Every contig logs
`candidate blow-up … 2,001 candidates against a 5,000 bp search window, over its ceiling of 2,000,
and was ABANDONED (the read passes through untouched)`. So on this data the resolver is **paying
near-worst-case enumeration cost and then discarding the result** — the expensive path and the
useless path are the same path. `max_intron=5000` is the yeast-calibrated default, but **no real
S. cerevisiae intron approaches 5 kb** (longest ≈ 2.5 kb), so the window is ~2× larger than the
biology requires and the cost is quadratic-ish in it (`candidates = near_sites × far_sites`).

**Asks, in priority order:**
1. **Re-benchmark the 2026-09-05 AT-AC default** to apportion the cost between the two
   clip-enumeration passes — as a *diagnostic*. AT-AC stays ON (see above); the goal is to know which
   pass to optimise, and whether the AT-AC pass can share candidate enumeration with the canonical
   pass instead of repeating it.
2. **Parallelise the resolver by contig** — it is already contig-sequential and per-read, and its
   `SpliceSiteIndex` is genome-derived and cached (`*.splice_sites.npz`), so a process pool over
   contigs is embarrassingly parallel and would recover the `-t N` the API already accepts.
3. **Lower the default `max_intron` per organism** (yeast ≈ 2,500), or derive it from the annotation
   when `--annotation` is supplied. Note this is **result-changing**, so it needs a 🧭DECISION, not a
   silent default change.
4. Make the single-threaded warning **name the wall-clock consequence** ("this stage will not use
   the other N−1 cores and may dominate runtime"), since the current wording reads as an API note.

**Filed by:** Chanfreau session `chanfreau-lab-65`, from `planning/907d_full_depth.md`. The run that
hit it needs 98.8 M reads through Stage A and is currently **blocked** on this.

---

## Open — reconciled 2026-08-20 (Re-aligner emphasis: Resolver + Stations A–C)

Tags: `✅CODE` = anchor re-verified in the tree today · `📄DOC` = carried from a planning/handoff
measurement · `🧭DECISION` = needs Kevin, changes scientific output. Station A = overhang
resolver (`core/align/overhang_resolver.py`), B = triage (`core/consensus/triage.py`),
C = pool-gate (`core/consensus/station_c.py`).

### A0. ✅ FIXED `3457ecc` + SMOKE PASS (fix C, Kevin-approved 2026-08-21; WT rep1 1M: W7 decoys 37→0, W17 3→0, annotated up everywhere, XN 100 %, XO↔flag 99.91 %; H2 `730_triage/730p_fixC_smoke/`) — minimap2 `-uf` on UNORIENTED PCR-cDNA manufactured decoy junctions (found 2026-08-21, planning/730)

**Fix:** `write_stage1_fastq` emits every Path-A molecule RNA-sense (rev-frame → reverse-complement
after the pretrim) and tags `XN:i:1`; XO unchanged; XN on the CMA whitelist; four-case test
`tests/test_cdna_emission_orientation.py`. **Residual:** every Path-A output WITHOUT `XN:i:1` is
pre-fix (684 cohort, the MS2 1M runs) — re-run, do not re-score. `-uf` itself is left as is (correct
now that reads are sense); Path B still carries both orientations and its consumers use XO/ro.

`rectify align` passes `-uf` unconditionally (`multi_aligner.py:601/2249/2393`), but Path-A
cDNA molecules are emitted in basecalled orientation (`ont_cdna.py` docstring; `XO` is only a
label). For the pA-first half, splice motifs are scored on the wrong strand → decoy junctions
at loci with an opposite-orientation canonical pair near a real intron (RPS24B: 4.6 % of blind
spliced reads in EVERY sample, 0 in the aware orientation; RPS4A GT..GG "Prp18" event likewise
0/33,452 aware). The resolver re-arbitrates and amplifies them. `ts:A:+` on every read = the
fingerprint. **Fix before the MS2 run:** orientation-normalise molecules to RNA-sense at
emission (`core/cdna/io.py:216`, XO + S1 flag) so `-uf` is valid for all reads; fallback `-ub`
for cDNA. Discriminator for existing data: the sequence-deficit test in the motif-aware
orientation (planning/730 §REVERSAL). DRS unaffected (reads are RNA-sense).

### A. Resolver / Station A

| # | item | tag | where | fix shape |
|---|---|---|---|---|
| A1 | **Non-canonical discovery gap** — recovers 99.1 % of mapPacBio-unique ANNOTATED junctions but only **34.5 %** of non-canonical ones (planning/721); structural: `SpliceSiteIndex` enumerates GT/AG-class sites only. Documented in the module docstring + `tests/test_resolver_noncanonical_gap.py`; 722 REFUTED an index extension; 722b shipped opt-in `--resolver-acceptor-classes prp18`. | ✅CODE (`overhang_resolver.py:11-20`) | Station A | Not a code fix. Open follow-ups: (i) the pre-registered upf1Δ flag-on/off test vs the 512 S1 junctions (722b §Recommended); (ii) **re-score 721 with ±6–8 bp boundary tolerance** — rbrowse's 722 eyeball shows mapPacBio carries a systematic 1–6 bp boundary error even on real introns (SUS1 dominant call 3 bp off, n=807), so the TOL=5 recovery table under-counts recovery; the "PARTIAL" verdict may soften. |
| A2 | **W_max calibration defect** (planning/700) — the search bound stops binding at ~15 bases (exponent saturates at an arbitrary 64, not `log2(G)`), and the uniform-genome null ignores the repeat tail. Worse on human. | 📄DOC | Station A | Saturate at `log2(G)`; SUS-based per-locus bound (= 690's rule, reached independently). **Not implemented, not costed** (cost = how many legitimate placements it refuses). |
| A3 | **`rectify split` chunked-HPC path NOT de-paneled** — still `OTHER_ALIGNERS_DEFAULT = ['minimap2','gapmm2','uLTRA','deSALT']` + a mapPacBio array; the resolver cannot be an "others" leg (needs the finished minimap2 arm). Deliberate, but the chunked path now diverges from `run-all`/`align` defaults. | ✅CODE (`split_command.py:106,405-408`) | Station A | Design decision: resolver leg after the minimap2 chunk task, or resolver at merge. |
| A4 | ✅ RESOLVED `todo/docs-a4` cb1469e+17a8e9d (landed on a `todo/*` branch overnight 2026-09-05 — NOT merged to master; see `dev/todo_run_20260904/README.md`) — **Docs frame the resolver as a clip re-placer; its measured value is v2 re-arbitration** (`arb_dmerge` 24,686 vs `XJ` 2,491 per 897k reads, planning/720). README:64 / ARCHITECTURE Step 1 still lead with "re-places terminal soft-clipped overhangs". | ✅CODE (README.md:64) | docs | One paragraph: "junction re-arbitrator first, clip resolver second". Also record the tag trap: `XB` is on EVERY read — "has XJ or XB" ≠ changed; use `XJ` + typed `XB` kinds (`dmerge/shift/dop/mm/mmL`). |
| A5 | planning/708 three-way triage is **unvalidated against real resolutions** (do the 2,491 `699b` resolutions survive it?). | 📄DOC | Station A/B | Run it on `699b`. |
| A6 | `_two_sided_candidate_hook` is a documented no-op (`dev/CLIP_LEG_WIRING_PLAN` §4b) — the two-sided donor-held/acceptor-held enumeration from 644h integration is not plugged in. | 📄DOC | Station A↔B seam | Post-landing resolver-library work; not a regression. |
| A7 | ✅ RESOLVED `todo/resolver-a7-a9-a8` c85da29+8cdd509 (landed on a `todo/*` branch overnight 2026-09-05 — NOT merged to master; see `dev/todo_run_20260904/README.md`): `_block_cigar` re-aligns the placed block (affine Gotoh, anchored at the junction); gap-free blocks emit the old flat M exactly; allowance `max(5, ceil(max_edit_frac*m))`; ANY gap touching the junction falls back to flat M (`emit_fallback_flat`) so the junction stays arbitrable; 23+2 tests — **Resolver accepts a clip on an indel-tolerant score but emits an ungapped CIGAR** (found 2026-09-03 on P06/RPS24B `r060_2601`, Chanfreau planning/862 §8b). `resolve_clip` scores candidates with `hp_edit_distance_bounded` (indels allowed), then `_rewrite_cigar` writes the matched segment as ONE flat `(0, m)` M op — no path for an indel — and drops MD/NM so nothing downstream recomputes. On P06 the 24 bp block has an exact 14-mer core with both flanks needing a 1-bp shift; the flat `24M` forces them into mismatches (13/24 apparent identity vs 16/24 best ungapped, vs the true gapped fit). Applies to EVERY resolved clip; inflates the read's mismatch/indel burden, which feeds the arbiter's own ED test and violates the spirit of the "a junction edit must not raise the read's indel burden" invariant. | ✅CODE (`overhang_resolver.py` `_rewrite_cigar` ~:540, scoring ~:470) | Station A | Route the accepted segment through `align_clip_to_exon()` (`local_aligner.py`, affine Gotoh) and emit its M/I/D CIGAR — exactly what `rescue_3ss_truncation()` already does for Cat3 5′ rescue. Then re-derive MD/NM or leave them absent consistently. |
| A8 | ✅ (i) DEFAULT-ON decided by Kevin 2026-09-05 for all organisms — 2H side landed f307311 `todo/sumner-panel-issues` (`_canonical_tier` recognizes the paired class at GC-AG rank via `_ATAC_DINUCS`; `atac=True` default in `canonical_in_class`; yeast byte-identical; +0.2 µs/call), resolver-side d20ddd7 `todo/resolver-a7-a9-a8` (`ResolverConfig.atac=True`, `--no-resolver-atac`, AT-AC pass gets its OWN candidate budget — the shared budget discarded canonical placements; yeast ×1.15 suite, human re-timed by the tester on 6485226 holdout1k: **×1.15 wall (210 s vs 182 s), candidates ×1.67, resolved 8 both ways** — R's ×2.93 superseded) — **MERGED into `todo/integration` 56abf3a (2026-09-05 ~12:50 PT; the 4-file conflict resolved keeping both sides; resolver gate 98/98; pushed to origin)** — **LANDED on `origin/master` 6485226 (2026-09-05 ~13:35 PT; full suite on the land branch = master's known state: 3,167 passed, the same 2 netseq flakes (pass standalone) + 6 worktree env errors; bisect A15 clean)** · ✅ follow-up (ii) DONE `todo/resolver-a7-a9-a8` c33d693 (landed on a `todo/*` branch overnight 2026-09-05 — NOT merged to master; see `dev/todo_run_20260904/README.md`): Case B1/B2/B3 + mm-rescue now enumerate the AT-AC pair as a separate paired pass; (i) default-on and (iii) SUT635 validation still open — **AT-AC introns now enumerable, opt-in** (2026-09-04, Kevin: Talkish et al. 2019 PLoS Genet 15:e1008249 — yeast splices AT-AC via the MAJOR spliceosome, SUT635; human U12-type needs it). `ResolverConfig.atac` / `--resolver-atac`: `SpliceSiteIndex` format v3 (+4 paired arrays), a second PAIRED clip-enumeration pass (never AT..AG / GT..AC — not the 722-refuted motif-free extension), rank GT 0 < GC 1 < AT-AC 2, `canonical_in_class(atac=)` at all arbiter sites, class-preserving Case A shifts. NOT wired: Case B D-op snap + mm-flagged linear rescue (GT/GC-only). Default off (722b precedent). Uncommitted; `tests/test_resolver_atac_introns.py` 12/12; vocabulary pins re-pinned to 11 kinds. Softens A1's "GT/AG-class only" — the acceptor (prp18) and AT-AC pair classes now both exist as flags. | ✅CODE (`overhang_resolver.py` `ResolverConfig.atac`, `_site_kinds_atac`; `splice_site_index.py` v3) | Station A | Open follow-ups: (i) decide default-on for yeast (needs a re-price of `alpha` vs candidate density); (ii) wire the pair into Case B / mm-flagged rescue; (iii) validate on SUT635 in a real library. |
| A9 | ✅ RESOLVED `todo/resolver-a7-a9-a8` 1b8031d (landed on a `todo/*` branch overnight 2026-09-05 — NOT merged to master; see `dev/todo_run_20260904/README.md`): SEQ decoded against the INPUT alignment on the first rewrite, qualities preserved, passthrough byte-identical; 9 tests — **Resolver preserves `=`-compressed SEQ from its input while rewriting the CIGAR** (found 2026-09-04 on P07, Chanfreau 748 doc §P07 + `862_eq_compressed_seq_survey.py`). The 748 fixture's minimap2 BAM is `=`-compressed on every read (42,409/42,409); the resolver decodes `=` for SCORING (`_decoded_query`, `inside_seq`) but writes the record back with the original `=` characters and a NEW alignment, so the `k_inside` bases it re-assigns across the junction decode against the wrong reference position (they read as matches at the new placement). Production minimap2 BAMs are not `=`-compressed, so the blast radius is `calmd -e` / CMA-expanded / fixture inputs — but the output is silently inconsistent whenever it happens. | ✅CODE (`overhang_resolver.py` `resolve_read`, `_rewrite_cigar`) | Station A | On any rewrite of a read whose SEQ contains `=`: replace SEQ with the decoded letters (decode against the INPUT alignment before rewriting; `_decoded_query` already does this), preserving `query_qualities` (pysam clears them on SEQ assignment). One test: `=`-compressed input read with k_inside > 0 → output SEQ has no `=` and the moved bases equal the input's decoded bases. Also worth asking why the fixture's `rectify align` output is `=`-compressed at all. |
| A10 | ✅ FIXED 2026-09-10 on `todo/resolver-a12-a13` — resolver move tag renamed `XB` → `XE` (5 call sites, tests, ARCHITECTURE Tags note; `tests/test_resolver_parallel.py::TestMoveTagIsNotShared` pins XE/XJ against `_CDNA_COMMENT_TAGS`, the CMA whitelist and every other `set_tag` in the tree). Only two throwaway debug scripts read the old tag (Chanfreau 644e/862). — **`XB` tag-name COLLISION** (found by the docs agent 2026-09-05, verified): the ONT cDNA pipeline writes `XB` = `n_top/n_bot` strand-split counts (`core/cdna/io.py:250`) and the resolver writes `XB` = its move family (`dmerge`/`shift:`/`dop:`/`mm:`/`mmL:`); `consensus._restore_sidecar_tags` (`consensus.py:408` lists XB in `_CDNA_COMMENT_TAGS`, overwrite at `:457-483`, call sites `:923`, `:980`) restores the cDNA value unconditionally, so on cDNA reads the resolver's XB is invisible in the final consensus BAM (read the resolver arm BAM or the `arb_*` stats instead). Also: the `0/1` seen on P07's read was this cDNA tag, not a resolver sentinel — the resolver writes XB/XJ ONLY on rewritten reads. | ✅CODE | cDNA ↔ resolver seam | Rename one side (resolver → e.g. `XR`) or make the sidecar restore skip a tag the resolver rewrote; documented as a Warning in ARCHITECTURE's Tags note (17a8e9d). |
| A11 | **Arbiter Case B2/B3 default arm splices to a chance GT..AG decoy** ~11–16 bp from the true site at WORSE ED (measured by agent R on the AT-AC synthetic genome; pre-existing — the planning/721 residual in miniature; `atac=True` fixes only the AT-AC instance). | 📄DOC | Station A | Needs the margin/ambiguity discipline of `_score_alts` applied to B1's best-vs-current rule, or an ED-improvement floor for mm-rescue proposals. Not a regression. |
| A12 | ✅ DONE `todo/integration` a28714a (2026-09-05, Kevin): **resolver numba kernel ON by default, imported lazily** (`overhang_informativeness._load_numba_kernel`, first large DP only; `RECTIFY_HP_ED_NUMBA=0` opts out; the correct-stage kernel in `splice_aware_5prime.py` stays opt-in because it runs in spawn workers where the documented OOM happened and human 2H/2F need 84 s/145k without it). Why: ISSUE-010 lifted the yeast-calibrated candidate cap on human windows and the pure-Python DP is then 26× slower (2,259 s vs ~135 s per 1,000 reads). Yeast: same outputs, same speed, workers never import numba. **Plus a loud checkpoint-schema guard:** `_rebuild_output_from_region_files` refuses (naming the dir) a checkpoint whose row width ≠ the current corrected-TSV header instead of silently dropping every row (D3 + ISSUE-002 took the TSV 39 → 44 columns). | ✅CODE | Station A / correct resume | Docs DONE 2026-09-05: README:68 + ARCHITECTURE § resolver now say the resolver kernel is default-on and the correct-stage kernel stays opt-in (R d20ddd7, merged 56abf3a; verified against `_HP_ED_NUMBA_REQUESTED` defaults '1' vs '0'). |
| A13 | **Per-pass candidate budget** (semantics change, R d20ddd7, merged 56abf3a): the optional AT-AC pass used to share `max_candidates_per_clip` with the canonical GT/GC pass and could starve it — cost 3 real yeast junctions on the AT-AC synthetic genome. Each paired pass now has its own budget, so worst-case per-clip work rises in proportion to the number of classes searched (2 today; 3 with `prp18`). | ✅CODE (`overhang_resolver.py` ~:655-670) | Station A | Tracked so a future "why did per-clip work go up" question has an answer; re-price if a third class is added. |
| A14 | ✅ FIXED 2026-09-10 on `todo/resolver-a12-a13` — `run_multi_aligner(resolver_acceptor_classes, resolver_atac, resolver_candidate_ceiling)` forwards them plus `max_intron` to the post-pass (`tests/test_resolver_parallel.py::TestKnobPlumbing`) — **`run_multi_aligner`'s resolver call ignores every resolver knob** (`acceptor_classes`, `atac`, `max_intron`): `multi_aligner.py:3504` calls `run_overhang_resolver(base_bam, genome_path, output_bam, threads)` only. Harmless today — the only callers are `tests/test_multi_aligner_failures.py` — but a live trap: anyone wiring it into a command would silently drop `--no-resolver-atac`. (Found by agent R 2026-09-05; verified in the merged tree.) | ✅CODE | Station A / align seam | Thread a `ResolverConfig` (or the three knobs) through `run_multi_aligner` → `run_overhang_resolver`, or delete the dead branch. |
| A16 | 🟡 OPEN 2026-09-05 (arbiter rectify-f5, events-tier walls, for the T3 run-all sizing) — `correct` on ~6k event reads costs ≈ 0.6 CPU-s/read at 4 workers but ≈ 1.3 CPU-s/read at 8: 8 workers do NOT scale on small inputs (per-worker startup / genome-cache load, or 2F candidate-loop contention). Not for now; when the T3 run-all is configured, pick worker counts from data (measure 2/4/8 on one events BAM) rather than `nproc`. | 📄PERF | correct/parallel | Measure once on Sherlock with the timing counter (567f30d) before T3; record the curve in `docs/architecture/pipeline_and_io.md`. |
| A15 | ✅ CLOSED 2026-09-05 — **machine load, not code.** Serial bisect of `tests/test_validation_reads.py` (one pytest at a time, same M1): 850388e 726 s at load-1 **85**; a28714a 266 s at load 15; 550fa24 228 s at load 10; 56abf3a (AT-AC default-on merge) **220 s at load 7** — the tip is the fastest, and R's 148 → 593 s was measured with 12+ pytest procs running. Neither 8427583 (window-scaled ceiling; yeast windows are byte-identical) nor f307311 costs anything measurable. | 📄DOC | tests | Nothing to fix. If the number ever matters, time it alone at load < 5; the 850388e row above is not a baseline. |

### B. Station B (triage) — opt-in, so non-blocking

| # | item | tag | fix shape |
|---|---|---|---|
| B1 | **Cost ~30× Station C** (751–1,062 s vs 34 s/sample; plain-DRS clip legs 48 min/100k reads at a 94.9 % triage rate — linear candidate scan in rescue). | 📄DOC | Entry preconditions + bisect-indexed lookup + memoization (planned, not done). |
| B2 | Fail-soft kept (`except Exception → WARNING`), now **machine-readable** via `<sample>.station_failures.json`. | ✅CODE (`run/stages.py:1106-1252`) | Closed as "loud enough" unless Kevin wants fatal. Smoke checklist: the JSON must be ABSENT. |

### C. Station C (pool-gate)

| # | item | tag | where | fix shape |
|---|---|---|---|---|
| C1 | ✅ **684c trio LANDED** (`fa7342a`): flags on BOTH tracks, annotation-derived length pre-gate + `over_max_intron` column, docstring contradiction fixed. Do not re-open. | ✅CODE (`station_c.py:402-447`) | — | — |
| C2 | **2–10 kb impossible-intron residual** (planning/689a): the consensus truncation is destructive so it cannot go below ~5.1 kb (minimap2's clean max N-op = 5,064 bp); ~15 % of the >2 kb population (DRS) / 8.7 % of gaps (cDNA) survive. Station C's new ≤`max_intron` pre-gate (5,000 yeast) now MARKS 5–10 kb non-destructively, so the reporting residual is **2–5 kb**. | 📄DOC | consensus + C | Needs a per-read **minimap2-disagreement** discriminator — design work, not a threshold. |
| C3 | ✅ RESOLVED `todo/plumbing-d5-d2-c3` c431fdc+854db78+ee5062a (landed on a `todo/*` branch overnight 2026-09-05 — NOT merged to master; see `dev/todo_run_20260904/README.md`): flag (and sibling `--junction-pool-min-anchor-bp`) on `run-all` + `align`, threaded single_sample → stages → `run_align` → `run_consensus_selection` — NOTE it governs the ALIGN-stage consensus selection, not the merge; 16 tests incl. a getsource count check on both call sites — **`--junction-pool-max-intron-len` unreachable in production** — exposed only on `rectify consensus`; absent from `align`/`run-all` (`run-all --max-intron` is the DIFFERENT aligner/Station-C knob). | ✅CODE (`consensus_command.py:113,424` only) | plumbing | Add the arg to `run-all` + thread it. Small. |
| C4 | **D/N `intronlen=40` threshold artefact** — deletions relabelled as introns, panel-wide incl. minimap2 (cliff at exactly 40; `dev/DN_THRESHOLD_ARTEFACT_20260812.md`). Raising `intronlen` is WRONG (kills ~1/3 of real introns). | 🧭DECISION | C (reporting-side, structurally zero cost to annotated introns) or D/N decision (alignment-side) | Kevin's call; nothing implemented. |
| C5 | **SUS / I_GR gate** (planning/690; 18× separation on long admits). 692 shows it adds nothing over `selfhom_flag` for Station C's long-admit verdicts; it survives as the RESOLVER's W_max calibration (A2). | 🧭DECISION | C / A | Probably moot at C now that C1 landed; decide explicitly. |
| C6 | Cross-sample recurrence + orthogonal-evidence columns are recorded but **not computed** (v0 single-sample); the 644h phase-2 **placement-uniqueness leg** (T6 decoy generalised) is still the missing axis for the non-canonical track. | 📄DOC | C | Build when the MS2 cohort exists (recurrence needs the cohort). |
| C7 | **Boundary-tolerant matching vs mapPacBio** — when mapPacBio is used as a scout, its junctions carry a systematic 1–6 bp placement error; exact-coordinate matching double-counts every real intron as jittered neighbours (rbrowse 722). | 📄DOC (rbrowse, 9 loci) | C / 721 | ±6 bp collapse before any recurrence or recovery count. |
| C8 | `ysh1_rep1`/`ysh1_rep2` (the two `pod5_skip` libraries; 178k/197k junctions vs 24k sibling) — quarantine from junction-level analysis? | 🧭DECISION | cohort | Excluding takes ysh1-AA to n=1. Junction count doubles as a library-QC metric. |

### D. Consensus / correct — bugs adjacent to the stations

| # | item | tag | fix shape |
|---|---|---|---|
| D1 | ✅ **719 negative-POS FIXED** (`17d1c35`: 5′ rescue refuses off-edge extensions; invariant two-sided at `consensus.py:172-173`). **Residual:** `4nqo_rep1`'s −80,234 bp overshoot implies a SECOND path (long N-op × rescue) — symptom now refused, mechanism uninvestigated. **2026-08-21 19:20Z: both samples RE-RUNNING on H2 job `14500215` (`r710fail2`) against the new git-clone tree `tree_master_95bc039`** (driver `710_correct_rerun_v2.sh`: skip requires `.bam` AND `.bai`; `--force-all`). rc=1 root cause confirmed from the logs = `samtools index` failure on the corrected BAM (off-edge record). | ✅CODE | Investigate the 80 kb path (still open). |
| D2 | ✅ RESOLVED `todo/plumbing-d5-d2-c3` 6fc4ef5 (landed on a `todo/*` branch overnight 2026-09-05 — NOT merged to master; see `dev/todo_run_20260904/README.md`): ERROR + `RuntimeError` naming the BAM; fatal to `refine_bam_junctions` itself — `correct`'s 2H step stays fail-soft (warning, continues on the UNREFINED BAM) but the broken BAM is no longer promoted — **`refine_bam_junctions` logs a `samtools sort`/index failure as WARNING and continues**, dying ~90 s later with an unrelated error. | ✅CODE (`junction_refiner.py:1858-1860`) | Make it fatal (raise). Small. |
| D3 | ✅ RESOLVED `todo/d3-x-tags` 5ca607e+45d0081+4ebfbd7+f266293 — the production gap closed too: `merge_corrected_tsvs(consensus_bam=)` joins the REAL Xa/Xc/Xn/Xt by normalized read name from `<sample>.multialigned.bam` (auto-discovered in run-all, `--consensus-bam` to override, path always logged; empty cells only; 0-match WARNING); `split` chunk merge deliberately unwired (no per-chunk consensus artifact) (landed on a `todo/*` branch overnight 2026-09-05 — NOT merged to master; see `dev/todo_run_20260904/README.md`): four SEPARATE columns appended at the END of the row — `consensus_aligner` / `consensus_confidence` / `consensus_n_agree` / `consensus_tied` — read once at intake by `output.consensus_tag_fields`; header 39→43; every reader resolves by name (audited: zero positional readers); golden re-recorded (pure addition). CAVEATS: populated only on the consensus-first path (`multialigned.bam` input) — empty by design in the correct-first per-aligner order, where `winning_aligner` is the merge verdict; discard stale `correct` checkpoint dirs across this commit (`parallel.py:498` drops rows whose width ≠ header). `Xn` is now per-read in the TSV as `consensus_n_agree` (what D6 needs). — **`correct` does not propagate `Xa`/`Xc`/`Xn`/`Xt`** into `corrected_3ends.tsv` (no hits in `core/bam`, `core/correct`). | ✅CODE | Add SEPARATE columns — never a pre-computed concordance class (rbrowse's constraint). |
| D4 | xdist flake: `test_validation_reads_upf1d cat7_minus_2 test_no_five_prime_rescue` fails under `-n 4`, passes serially. | 📄DOC | Chase corrected-fixture ordering; don't blame the change du jour. |
| ✅ RESOLVED `todo/plumbing-d5-d2-c3` 4dfbc67 (landed on a `todo/*` branch overnight 2026-09-05 — NOT merged to master; see `dev/todo_run_20260904/README.md`): test chdirs into tmp_path; stray file removed from the main tree 2026-09-04 — D5 | `tests/test_universal_command_provenance.py:66` exercises the uncreatable-`output_dir` fallback, which writes `export.command.provenance.json` into **cwd = repo root** (stray file dated 08-17). | ✅CODE | `monkeypatch.chdir(tmp_path)` in that test; delete the stray file. |
| D7 | ✅ **FIXED in tree (uncommitted, 2026-08-21)** — `cdna-analyze` recomputed `tail_len` by walkback on the PRETRIMMED post-align record and `set_tag("XA", 0)` over the carried value → 100 % `tail_len==0` on a 4.63 M-molecule library (Chanfreau 736). Now: walkback wins only when it finds a tail in SEQ, else the carried `XA:i` is honoured; provenance INFO line; WARNING at >95 % zeros. Tests `test_cdna_analyze.py::test_cdna_analyze_honours_carried_xa_on_pretrimmed_records`, `::test_cdna_analyze_warns_when_tail_len_is_dead`. See `AGENT_FIXES.md [2026-08-21]`. | ✅CODE | Commit (Kevin). |
| D8 | 🔴 **`correct` HANGS FOREVER when a BAM-region worker dies (found 2026-08-22, 682 `correct_v2`):** the 4 "slow" DRS samples (`dis3rrp6_rep1`, `rna15_rep1`, `ysh1_rep1`, `ysh1_rep2`) each lose a spawned worker to glibc `free(): invalid next size` (heap corruption in a native extension, 2–3 min into region processing, reproducible — same samples left `core.8482` on 08-12), and `ctx.Pool(...).imap` (`bam/parallel.py:862`) then waits forever: Pool respawns the worker but never re-queues the lost task. 9.7 h of idle wallclock × 4 × 8 slots before it was noticed by a flat CPU counter. TWO bugs: (a) the corruption itself — data-dependent, kernel unknown (numba hp-ED kernel looks bounds-safe: inputs truncated to 200); 746 diagnostics on H2 (`746_hang/A_faulthandler` = `PYTHONFAULTHANDLER=1` for the Python frame at SIGABRT; `B_nonumba` = `RECTIFY_HP_ED_NUMBA=0`); (b) the hang — replace the bare `imap` with a dead-worker-aware loop (`concurrent.futures.ProcessPoolExecutor` raises `BrokenProcessPool`; or `apply_async` + per-region timeout + `pool._pool` exitcode check) so a worker death FAILS the run loudly with the region name. **(b) FIXED in tree (local commit, not pushed): `_iter_pool_results` + `RegionWorkerDied`, both imap sites; `tests/test_parallel_dead_worker.py` aborts a real spawned worker.** **(a) ROOT-CAUSED + FIXED in tree (uncommitted): pysam `get_reference_sequence()` on an MD/CIGAR-inconsistent read (spurious multi-N-op reads that evaded the 50 kb max-N guard: ysh1_rep1 has 2, ysh1_rep2 1, 4nqo_rep1 0) → `extract_deletions` now gates on MD consistency + total N and takes deleted bases from the genome dict when available (`utils/alignment.py`, `indel_corrector.py:1837`; `tests/test_extract_deletions_md_guard.py`). **VERIFIED 2026-08-22 05:45 PT: ysh1_rep1 completes in 24:41 on the patched tree (rc=0, indexed corrected BAM, 0 error signatures). Commits `5ed4bd5` (hang) + `5e688a2` (corruption), local.** | ✅CODE | push (Kevin) → re-run the 4 samples → 48/48. |
| D6 | Agreement analyses: matched-support control (696 confounds support with class) and the single-arm exemption for `Xn>=2` are untested — do before any `Xn` threshold. | 📄DOC | Analysis, not code. |

### E. Production line / ops (the MS2 mission)

| # | item | state 2026-08-20 (live probe) |
|---|---|---|
| E1 | 682 `correct_v2` re-run | **46/48 at 2026-08-22 05:30 PT; the 4 "slow" samples are the D8 hang — `14382552` qdel'd after 9.7 h idle; 746 diagnostics running.** Earlier: `14382552` (4 slow, tree 0128840) sat `qw` on highp 08-17→08-21; `qalter -l h_data=8G,h_rt=86400` (highp dropped) at 19:17 → all 4 tasks running within seconds (pod_smp.q). `4nqo_rep1`/`dis3rrp6_rep2` → `14500215` on tree `95bc039` (D1). The 48 will carry TWO `rectify_git_sha` values (46 × 0128840, 2 × 95bc039); rbrowse informed (inbox 19:35Z). `sen1_rep1` is an ORPHAN 49th dir (08-12 trim+align only, never in `710_samples.txt`). rbrowse's nsd re-slice is blocked on 48/48 clean. |
| E2 | ✅ **`/u/scratch/k/kevinroy/tree_master_95bc039/rectify` exists (2026-08-21): a real `git clone --depth 1` of origin/master, so sidecars read `git:95bc039…`, not `pkghash`.** Older trees remain: `tree_master_7c6a536` — **6 commits behind** `origin/master 95bc039`, including `e532ec0` (base long-read default `['minimap2','mapPacBio','gapmm2']` → `['minimap2']`; bare-`align` junction default += uLTRA/deSALT; protocol-keyed short-read panels). A run-all on that tree still launches mapPacBio+gapmm2 unless `--aligners minimap2` is explicit. | Advance the tree (suite-gated 2,379/0) or pin the flags in the driver. |
| E3 | MS2 full-depth: driver + smoke + 33-sample fan-out **not started**. prp5 staging DONE — `_pre.bam` on **scratch** `/u/scratch/k/kevinroy/662_cdna/prp5_rep{2,3}/` (17.2 GB / 1.2 GB; NOT `/u/project`, not backed up). | — |
| E4 | RTD dead since 2026-04-05 (webhook secret) — 1-min Kevin action (Resync webhook). Stale Zenodo hook with a token in its URL. | Kevin |
| E5 | git-history purge (filter-repo + force-push) for handoffs/dev docs — decision. | Kevin |
| E6 | Local `master` ref stale at `0128840` (origin/master = HEAD = `95bc039`): `git branch -f master origin/master`. | housekeeping |
| E7 | 676 uBAM deletion (313 GB, file list needs sign-off; md5 is the wrong test) · pod5 preservation (raw signal on two external drives only) · 683 G-arbiter patch parked (land as latent guard or drop). | Kevin |

### F. Carried from `dev/TODO.md` (bucket 2, unchanged)

DRS validation-bundle regen (unblocked since NEW-082) · `rectify/data/validation/corrected_3ends.tsv`
re-emit · cat7 selector ranks by intron length · per-cluster centre-of-mass emission · RN-metadata
sidecar stays archived (decided).

### Not bugs — scope (Kevin's own notes, untracked)

`TODO_KR_NOTES.txt` = the final-deliverable programme (CPA clusters → PCA/heatmaps → DESeq2 at cluster
AND gene level with auto-detected reference samples → summary tables → GO → 4 de-novo motif scans).
**INVENTORIED 2026-09-11 (rectify-85): every step is ALREADY IMPLEMENTED and wired in `rectify analyze`
(`analyze_command.py` steps 1–9; manifest mode mirrors it): `detect_control_samples` (wt / wildtype / control / ctrl /
untreated / mock / dmso / reference patterns, `deseq2.CONTROL_PATTERNS`), `cluster_cpa_sites_adaptive`,
`run_pca_analysis` + `plot_pca`, `plot_sample_heatmap` / `plot_cluster_heatmap`, `run_deseq2_gene_level` +
`run_deseq2_cluster_level` + `summarize_deseq2_results`, GO up/down (`run_go_enrichment` + `plot_go_enrichment`,
padj < 0.05 & |log2FC| > 1, ≥ 10 genes), `run_differential_motif_analysis` = exactly the 4 scans (enriched/depleted ×
upstream 100 / downstream 50, STREME/MEME), `analyze_cluster_shifts` + `get_top_shifted_genes(n_top=20)` +
`plot_gene_browser`. "Port for production readiness" therefore means EXERCISE it on a real multi-condition manifest,
not build it. Two gaps against the note: (i) motif discovery runs only when enriched ≥ 20 AND depleted ≥ 20 — one thin
side silences all four scans; gate each side independently; (ii) no matplotlib plot of the discovered motifs (only the
STREME/MEME HTML) — add logo plots (`logomaker`) of the top motif per scan into `plots/`. Neither built tonight (no
multi-condition data on the M1 to smoke-test the path on).**
`Rectify_readme_changes_KR_proposed.txt` = README items (species-agnostic framing, per-library
quickstarts with run-all + sequential + chunked-array examples, `[visualize]` extra rationale, penalty
tables beyond 8 / human tables / pan-species question, configurable isoform-cluster tolerance,
walkback figure pA tail, adaptive-valley figure realism, chimeric-reconstruction evaluation). RE-DIFFED
2026-09-11 (rectify-85): every item was already in README/figures except four, landed in `db7c796`; only the
chimeric-reconstruction EVALUATION remains (an experiment).

---

## Closed this month — listed so nobody re-opens them

| item | resolution |
|---|---|
| Station C `flagged` unreachable on the canonical track · docstring contradiction · no length pre-gate (684c) | `fa7342a` |
| 719 negative POS (one-sided invariant) | `17d1c35` |
| Resolver not in the default panel · `max_intron` hard-coded 5000 · per-aligner `correct` ran on the RAW minimap2 arm | `fa63c12` (`derive_max_intron` = 2× longest annotated intron, yeast derives exactly 5000; `_collect_per_aligner_bams` substitutes the resolver BAM) |
| mapPacBio/gapmm2 in code defaults (align + run-all) | `e532ec0` — NOT the `split` path (A3) |
| Station C `ambiguity_window` IndexError · `Xn` tag collision (`Xi`) · `pool-gate -o` silent overwrite · cDNA frame-trim (679) · gapmm2 NEW-082 | master (`0128840` era); see `dev/TRIAGE_PRE_RERUN_20260812.md` §Closed |
| `cdna-analyze` zeroes carried `XA` (D7) | working tree 2026-08-21 (uncommitted) — `cdna_analyze_command.py` honours carried `XA:i` when walkback finds no tail |
| 722 GT/AG index extension | REFUTED by measurement (flat motif spectrum, 86 % singletons, ×191 candidates); rbrowse's 722 eyeball of all 9 recurrent/singleton loci: **none real** (poly(A)-tail snaps + mapPacBio boundary error) |

---

## Previously open (pre-2026-08)

**None.** (NEW-082 resolved 2026-06-24 — see below.)

---

### ~~NEW-082 (MEDIUM) — gapmm2 v25.8.12 drops most reads — ROOT-CAUSED + FIXED, verified 2026-06-24~~
**File:** `rectify/core/commands/install_aligners_command.py` (the unpinned installer).
**Found:** 2026-06-01 (M1); **root-caused + fixed:** 2026-06-24.

**Root cause (proven):** gapmm2 **25.4.13+** added a `splice_aligner_minimap2()` code path
(used whenever the `minimap2` binary is on PATH) whose tag-validation loop does
`if ts is None: continue` (align.py ~855). minimap2 emits the `ts:A:` transcript-strand tag
**only for spliced alignments** (those with an N/intron op), so gapmm2 **silently drops every
single-exon (unspliced) read**. The drop is invisible in `--debug` (only `low-mapq` is printed;
the `missing-ts` count is tracked but never surfaced). v25.4.5 (H2/Sherlock prod) predates this
path entirely (no `splice_aligner_minimap2`, no `ts` filter).

**Empirical proof (M1, 36 DRS validation reads → 35 after the empty-seq placeholder is cleaned):**
- gapmm2 **25.8.12**: minimap2 aligns all 35; gapmm2 emits **10/35** — and the 10 survivors are
  exactly minimap2's `ts:`-tagged (spliced) reads (`comm` diff empty). The 25 dropped are all
  single-exon.
- gapmm2 **25.4.5**: emits **35/35**. Verified end-to-end through rectify's own `run_gapmm2`
  (clean FASTQ 35 → PAF 35 → BAM 35 mapped primaries, one per input read).

**Fix:** pin gapmm2 to `==25.4.5` in `install_aligners_command.py` (line ~234) to match the
existing `pyproject.toml [aligners]` pin. The stale `gapmm2>=0.2` in the installer is what let
25.8.12 onto the M1 in the first place (the pyproject dep was already correct). M1 base env
downgraded to 25.4.5 to clear the drift. **Unblocks** the DRS validation-bundle regen (TODO).

---

**Previously: None.** As of 2026-05-23 every tracked NEW-### bug has been verified resolved
in the current `drs-validation-rebuild` code: NEW-066–074 + NEW-076 landed in
prior commits; NEW-077/078 fixed by the `1ab71f0` artifact-N rework; NEW-080
fixed in `cc86cc0`; NEW-079/081 already correct in code; NEW-075 (dT-primed
walkback) integrated with gene-strand emission. Resolved entries (with
verification pointers) are retained below and in the strikethrough archive.

---

### ~~NEW-077 (HIGH) — `bam_processor.py` minus-strand artifact N snap can emit an intronic / skipped coordinate — Fixed by `1ab71f0`, verified 2026-05-23~~

**Verified:** the artifact-N rework in `1ab71f0` (2026-05-22) snaps minus-strand to the first ALIGNED base after the N. `rectify/core/bam/bam_processor.py:841-844`: `if strand == '+': current_position = _art.junction_start - 1` else `current_position = _art.junction_end` — exactly the prescribed fix, with an explicit comment (837-840) that `junction_start` would land inside the artifact gap.

**Original symptom (for history):** the routine used to set `current_position = junction_start` on the minus strand, landing the 3' end on a skipped base.

**Source:** caught by the user's multi-agent coordinate-convention audit, 2026-05-20.

---

### ~~NEW-078 (HIGH/MEDIUM) — `bam_processor.py` minus-strand artifact ambiguity clipping leaves NET-seq window touching the N span — Fixed by `1ab71f0`, verified 2026-05-23~~

**Verified:** `rectify/core/bam/bam_processor.py:855-859` (Case B, minus) now clips `ambiguity_min` to `_art.junction_end` (first aligned base on the 3' side of the N), not `junction_start`, keeping NET-seq refinement off the skipped span.

**Original symptom (for history):** minus-strand clipping set `ambiguity_min = junction_start`, leaving the refinement window starting at the first skipped base.

**Source:** caught by the user's multi-agent coordinate-convention audit, 2026-05-20.

---

### ~~NEW-079 (MEDIUM) — `clustering.py` plus-strand `distance_to_gene_3prime` is off by 1 — Already fixed, verified in code 2026-05-23~~

**Verified:** `rectify/core/analyze/clustering.py:530-533` already computes `gene_3prime = gene['end'] - 1` for `+` strand and `gene['start']` for `-` strand — exactly the prescribed fix. `loaders.py:469` confirms the 0-based half-open convention (`start = int(fields[3]) - 1`). No code change needed.

**Original symptom (for history):** the bug report claimed `gene['end']` was used as the plus-strand gene 3' coordinate, making `distance_to_gene_3prime` 1 bp too large. The current code does not have this bug.

**Source:** caught by the user's multi-agent coordinate-convention audit, 2026-05-20.

---

### ~~NEW-080 (MEDIUM) — `false_junction_filter.py` minus-strand wrong-flank lookup against artifact N target — Fixed 2026-05-23~~

**Fix:** the genomic *target* flank now mirrors the (already strand-branched) read-side flank in `analyze_junction_for_artifact`. Plus strand fetches `[end, end+20)` (A-richness); minus strand fetches `[max(0, start-20), start)` (T-richness) — the genomic T-run the read's poly-T aligns over, upstream of the N. Previously `[end, end+20)` was fetched for both strands, so minus-strand genomic-A-tract artifacts were missed unless the `not has_canonical_motif` fallback also fired.

**Verified:** regression test `tests/test_false_junction_filter.py::test_minus_strand_artifact_uses_upstream_genomic_flank` fails pre-fix (`assert 0 == 1`) and passes post-fix; splice blast-radius slice (150 tests) green.

**Source:** caught by the user's multi-agent coordinate-convention audit, 2026-05-20.

---

### ~~NEW-081 (LOW) — `analyze_command.py` / `analyze/manifest.py` rDNA exclusion uses inclusive `<= end` — Already fixed, verified in code 2026-05-23~~

**Verified:** both rDNA exclusion sites already use half-open `>= start) & (< end)` — `rectify/core/commands/analyze_command.py:212-213` and `rectify/core/analyze/manifest.py:308-311`. (The `start <= pos <= end` at `analyze_command.py:95` is an unrelated *cluster*-membership lookup with inclusive cluster boundaries, built via `IntervalTree[start:end+1]` at `:80` — correct as-is.) No code change needed.

**Source:** caught by the user's multi-agent coordinate-convention audit, 2026-05-20.

---

### ~~NEW-075 (HIGH) — `correct`: dT-primed cDNA walkback integrated; emits gene/RNA strand — Fixed, verified 2026-05-23~~

**Verified:** the read-vs-genome walkback is integrated into `rectify correct` and emits the **gene/RNA strand**, not the BAM strand. `rectify/core/correct/walkback.py` holds the protocol-agnostic core (`walkback_3prime`, commit `986a19d`) plus `walkback_drs`/`walkback_drs_full`; `rectify/core/correct/protocols/quantseq_rev.py` is the `--dT-primed-cDNA` antisense wrapper. `bam_processor.py:775` calls `walkback_quantseq_rev(...)` returning `_gene_strand_wb` and `:798` calls `walkback_drs_full`; `correct_command.py:1374` defines `--dT-primed-cDNA` and the antisense strand flag (361-363). Regression covered by `tests/test_quantseq_rev_walkback.py` (gene-strand assertions). See `dev/specs/ISSUE_walkback_integration_retrospective.md` — all acceptance criteria checked off, validated on Han 2023 wt_R1.

**Remaining as deferred follow-ups (out of scope per the retrospective):** dedicated ONT `ont_drs.py` / `ont_cdna.py` protocol wrappers — file separately; DRS is already handled by `walkback_drs_full`.

**Original symptom (for history):** walkback lived only in the post-hoc `11_polya_walkback_recompute.py`, which emitted the BAM strand and corrupted cluster taxonomy at convergent loci.

---

### ~~NEW-074 (MEDIUM) — `consensus.py`: `shutil.copy2()` to NFS without `os.fsync()` risks partial writes — Fixed (landed on drs-validation-rebuild), verified in code 2026-05-23~~

**Verified:** `consensus.py:131` `_copy2_and_fsync`.

**File:** `rectify/core/consensus/consensus.py` (~line 1761)

**Symptom:** Sorted BAM copied from local scratch to NFS Oak with `shutil.copy2()`. Under NFS load, OS buffers the write; if the job crashes or the NFS connection drops before the buffer flushes, the destination BAM is silently truncated. Downstream readers crash with "truncated file".

**Fix:**
```python
shutil.copy2(src, dst)
with open(dst, 'r+b') as f:
    os.fsync(f.fileno())
```
Apply to both the sorted BAM copy and the `.bai` copy at the same location. Alternatively, switch to `rsync` (already used in `sync_to_oak()`) which performs checksum-verified transfer.

**Discovered:** 2026-04-22 post-mortem audit (job 22314279)

---

### ~~NEW-073 (MEDIUM) — `bam_processor.py`: 17× redundant full BAM scans before correction begins — Fixed (landed on drs-validation-rebuild), verified in code 2026-05-23~~

**Verified:** `bam/parallel.py:262` `_prescan_bam_once` (one `bam.fetch(until_eof=True)`).

**File:** `rectify/core/bam/bam_processor.py` (lines 1335, 1352, 1768)

**Symptom:** Three independent full BAM scans run before any correction worker starts:
1. Line 1335: full linear scan for read-flag stats (`ProcessingStats`)
2. Line 1352: calls `get_processing_regions()` → `find_coverage_gaps()` → one `bam.fetch(chrom)` scan per chromosome (16 scans for yeast)
3. Line 1768: another full linear scan for the variant-aware prescan

For a 200k-read yeast BAM this consumes ~47 min of wall time before correction starts, contributing to the observed 6h timeout.

**Fix:** Add `prescan_bam_once(bam_path)` that performs a single `bam.fetch(until_eof=True)` pass, simultaneously accumulating flag stats, building per-chromosome coverage intervals for gap-finding, and optionally running `VariantAwareHomopolymerRescue.scan_read()` when variant-aware mode is enabled. Call this once before spawning workers; remove the redundant scans.

**Discovered:** 2026-04-22 post-mortem audit (job 22314279); estimated ~47 min wall time savings per sample

---

### ~~NEW-072 (MEDIUM) — `multi_aligner.py`: deSALT/uLTRA sort subprocess deadlock — Not a bug, verified 2026-05-23~~

**Verified (NO-OP):** in both `run_ultra` (view_proc/sort_proc Popens ~1748-1755) and `run_desalt` (~1958-1965), the `view_proc` sets only `stdout=subprocess.PIPE` (stderr **inherited**) and `sort_proc` sets neither `stdout` nor `stderr` (both inherited). There is no `stderr=PIPE` left undrained, so `sort_proc.communicate()` cannot deadlock. The deadlock pattern requires a producer with `stderr=PIPE` that nothing reads while blocking on the consumer — that is the minimap2 case (~388/399, which correctly drains via a daemon thread), not deSALT/uLTRA. No code change needed.

**File:** `rectify/core/align/multi_aligner.py` (`run_ultra` ~1566, `run_desalt` ~1818)

**Symptom:** `sort_proc` spawned with `stderr=PIPE` (or inheriting parent stderr). If samtools sort writes ≥64 KB of warnings, the OS pipe buffer fills; since no thread drains it, `sort_proc.communicate()` deadlocks. (The minimap2 path already has a drain thread at lines 209–216; deSALT does not.)

**Note:** mapPacBio was already fixed (`stderr=subprocess.DEVNULL` with explanatory comment). This affects deSALT and uLTRA only.

**Fix:** Add a daemon stderr-drain thread matching the minimap2 pattern (lines 209–216):
```python
desalt_stderr_chunks = []
def _drain_stderr():
    desalt_stderr_chunks.append(sort_proc.stderr.read())
drain_thread = threading.Thread(target=_drain_stderr, daemon=True)
drain_thread.start()
sort_proc.wait()
drain_thread.join()
```

**Discovered:** 2026-04-22 post-mortem audit (job 22314279); closed as not-a-bug 2026-05-23.

---

### ~~NEW-071 (MEDIUM) — `run_command.py`: `stderr=DEVNULL` on samtools fastq hides error context — Fixed (landed on drs-validation-rebuild), verified in code 2026-05-23~~

**Verified:** `single_sample.py:36-55` `_run_samtools_fastq` (no DEVNULL; stderr included in RuntimeError).

**File:** `rectify/core/commands/run_command.py` (line 751, DRS trim → FASTQ conversion)

**Symptom:** `subprocess.run(['samtools', 'fastq', ...], check=True, stderr=subprocess.DEVNULL)`. If samtools fails, `check=True` raises `CalledProcessError` (error is caught), but all samtools stderr is discarded, making the failure impossible to diagnose. Line 1157 (same operation in single-sample path) correctly omits `DEVNULL` — only line 751 is affected.

**Fix:** Remove `stderr=subprocess.DEVNULL`. Capture or let stderr propagate; `check=True` already raises on non-zero exit.

**Discovered:** 2026-04-22 post-mortem audit (job 22314279)

---

### ~~NEW-070 (MEDIUM) — `run_command.py`: non-atomic `unlink()` + `rename()` instead of `Path.replace()` — Fixed (landed on drs-validation-rebuild), verified in code 2026-05-23~~

**Verified:** atomic `Path.replace` at `single_sample.py:343` and `:790`.

**File:** `rectify/core/commands/run_command.py` (lines 878–879, 1403–1404)

**Symptom:** Two locations use:
```python
_restored_bam.unlink(missing_ok=True)
Path(_restored_tmp).rename(_restored_bam)
```
If a crash or SIGKILL occurs between the two operations, `_restored_bam` is permanently absent (data loss). On NFS, `rename()` across directories can also fail.

**Fix:** Replace with the atomic POSIX `replace()` call:
```python
Path(_restored_tmp).replace(_restored_bam)
```
`Path.replace()` is atomic on POSIX filesystems and handles the unlink+rename as a single syscall.

**Discovered:** 2026-04-22 post-mortem audit (job 22314279)

---

### ~~NEW-069 (MEDIUM) — `consensus.py`: `.done` checkpoint sentinel written without BAM validity check — Fixed (landed on drs-validation-rebuild), verified in code 2026-05-23~~

**Verified:** `consensus.py:137` `_validate_bam_sample`, called at `:863` before trusting checkpoint batches.

**File:** `rectify/core/consensus/consensus.py` (lines 1667, 1698)

**Symptom:** After each batch BAM is written and `out_bam.close()` is called, the `.done` sentinel file is written immediately. If the batch BAM is corrupt (CIGAR/seq mismatch, partial write), the sentinel marks it complete. On resume, `pysam.cat()` or `pysam.sort()` reads the corrupt batch BAM and crashes.

**Fix:** Before writing the sentinel, call `_validate_bam_cigar(_cur_batch_path, sample_size=1000)` (the same function added for pre-sort validation). If validation fails, log the error and do NOT write the sentinel — the batch will be re-written on resume. A simple file-size check (`os.path.getsize > 0`) is insufficient: CIGAR-corrupted BAMs have non-zero size.

**Discovered:** 2026-04-22 post-mortem audit (job 22314279)

---

### ~~NEW-068 (MEDIUM) — `multi_aligner.py`: no assertion that minimap2 succeeded after aligner loop — Fixed (landed on drs-validation-rebuild), verified in code 2026-05-23~~

**Verified:** `multi_aligner.py:2103`.

**File:** `rectify/core/align/multi_aligner.py` (after aligner loop, ~line 1267)

**Symptom:** When all aligners are tried in a loop, exceptions are caught and logged per-aligner (`logger.error(f"Aligner {aligner} failed: {e}")`), but the function returns `results` with however many aligners succeeded — including zero. The caller receives an empty dict with no indication that minimap2 (the required baseline aligner) failed. Downstream consensus selection runs with no aligners and silently produces no output.

**Fix:** After the aligner loop, assert that minimap2 succeeded:
```python
if 'minimap2' not in results:
    raise RuntimeError(
        f"minimap2 (required baseline aligner) failed — "
        f"cannot proceed with consensus selection. "
        f"Succeeded: {list(results.keys())}"
    )
```

**Discovered:** 2026-04-22 post-mortem audit (job 22314279)

---

### ~~NEW-067 (MEDIUM) — `align_command.py`: `calmd_bam.replace(rectified_bam)` before index creation leaves BAM without `.bai` on index failure — Fixed (landed on drs-validation-rebuild), verified in code 2026-05-23~~

**Verified:** `align_command.py:34-36` (index temp BAM, then atomic replace of bam + bai).

**File:** `rectify/core/commands/align_command.py` (lines 576–579)

**Symptom:** Current code:
```python
calmd_bam.replace(rectified_bam)          # destructive rename
_sp.run(['samtools', 'index', str(rectified_bam)], check=True)   # might fail
```
If `samtools index` fails (disk full, permissions, samtools crash), `rectified_bam` exists without a `.bai`. Downstream tools requiring an index fail with "index not found", and the skip-logic in `run_command.py` (now using `_validate_bam_integrity()`) will correctly detect the missing `.bai` and re-run alignment — but the original `calmd_bam` has been destroyed, so alignment must restart from scratch.

**Fix:** Create the index first on the temp BAM, then do two atomic renames:
```python
if result.returncode == 0 and calmd_bam.stat().st_size > 0:
    _sp.run(['samtools', 'index', str(calmd_bam)], check=True)
    calmd_bam.replace(rectified_bam)
    Path(str(calmd_bam) + '.bai').replace(Path(str(rectified_bam) + '.bai'))
    logger.info("  MD tags added successfully")
```

**Discovered:** 2026-04-22 post-mortem audit (job 22314279); not present in squishy-waddling-stream.md

---

### ~~NEW-066 (MEDIUM) — uLTRA reuses stale cached database with empty genome, silently produces zero alignments — Fixed (landed on drs-validation-rebuild), verified in code 2026-05-23~~

**Verified:** fixed in `multi_aligner.py` `run_ultra`.

**File:** `rectify/core/align/multi_aligner.py`, `run_ultra()` (~line 848)

**Symptom:** uLTRA exits 0 but produces a SAM with only `@SQ` headers and no alignment records. Observed for chunks 0–6 of `wt_by4742_rep1_chunked_20260412` across three separate SLURM runs.

**Root cause:** uLTRA caches its genome index in `{output_bam.stem}_ultra_tmp/database.db`. When the first run fails mid-way (after the `ultra_tmp/` dir is created but before the genome is indexed), the `_tmp_dir` decompressed genome is cleaned up. On retry, uLTRA prints "Database found in directory — using this one" and skips genome indexing. `refs_sequences.fa` is 0 bytes → "Number of ref seqs in fasta: 0" → all chromosomes excluded → zero alignments. `run_ultra()` only checks `sam_path.exists() and sam_path.stat().st_size > 0`, but uLTRA writes a valid-looking SAM with just headers even when it aligned nothing.

**Fix:** After `ultra_out_dir.mkdir()`, check if `database.db` exists but `refs_sequences.fa` is empty or missing — this is the fingerprint of a failed prior indexing run. Remove the stale directory and recreate it so uLTRA re-indexes from scratch. A valid cache (non-empty `refs_sequences.fa`) should be preserved — the index is genome/GTF-derived and is intentionally shared across chunks. Fixed in `multi_aligner.py` `run_ultra()`.

---

### ~~NEW-076 (HIGH) — `rectify analyze` per-condition bedgraph: 1-bp left shift in `start`/`end` BED columns — Fixed (landed on drs-validation-rebuild), verified in code 2026-05-23~~

**Symptom:** Per-condition strand-specific bedgraphs from `rectify analyze` showed their 3'-end pileup peak 1 bp to the LEFT of the true position. Caught 2026-05-20 in the He et al cross-modality TRT analysis: CST6 ysh1 peak in `han2023_analyze_wbfix_20260512/analyze/bedgraph/Ysh1AA_plus.bedgraph` at `chrIX 287748 287749` (IGV 287,749) vs the matching DRS bedtools-derived 3'-end peak at `chrIX 287749 287750` (IGV 287,750). Affected DRS, QuantSeq REV, and PCR-cDNA equally — the emitter is protocol-agnostic.

**Root cause:** the writer treated `corrected_3prime` (which is 0-based-inclusive, from `reference_end - 1` / `reference_start`) as 1-based, doing `start = pos - 1; end = pos`. Correct: `start = pos; end = pos + 1`.

**Fix:** corrected in 3 emitters — `rectify/core/analyze/bedgraph.py`, `rectify/core/analyze/manifest.py`, and `scripts/generate_bedgraph_from_polished.py`. Audit complete (no further bedgraph emitters affected). Regression tests added in `tests/test_analyze.py::TestBedgraphCoordinates`. See `dev/audits/bedgraph_coordinate_audit_20260520.md` for the convention statement and per-file audit table.

**Blast radius:** every prior `rectify analyze` output's `bedgraph/` subdir has 1-bp-left-shifted per-condition tracks and needs regeneration. Clustering, shift analysis, and per-position attribution work from the `corrected_reads.tsv` / `corrected_3ends.tsv` position columns directly and are unaffected.

---

### ~~NEW-061 (HIGH) — `clip_read_to_corrected_3prime` shrinks N ops instead of eliminating them — Fixed 2026-04-09 (v2.7.9)~~

When `polya_walkback` produces a corrected_3prime that falls *inside* a near-3' N op (intron artifact), `clip_read_to_corrected_3prime` partially trims the N rather than eliminating it. The walkback position is numerically correct (validated shifts of −45, −72 for cat4_plus_1/2) but a stub N op remains in the rectified CIGAR (1520→1501 bp, 100→40 bp). Any downstream junction caller will report this residual N as a spurious splice junction.

**Root cause:** The clipping loop walks backward consuming `n_ref_clip` reference bases; when it runs out of bases mid-N, it shrinks the N and stops. There is no special handling to detect "target is inside N → snap to junction_start−1."

**Fix (two options):**
- Option A: In `clip_read_to_corrected_3prime`, detect when walking lands inside an N op and snap the clip to the N's start (fully clip the N and all trailing ops). Update corrected_3prime in the TSV to `junction_start − 1`.
- Option B: In `bam_processor._run_correction`, when FJF detects a near-3' artifact junction, propagate FJF's `corrected_3prime` (junction_start−1) and do not override it with the walkback position when the walkback lands inside the artifact N span.

**Files:** `rectify/core/bam/bam_processor.py` (`clip_read_to_corrected_3prime`), `rectify/core/splice/false_junction_filter.py`
**Discovered:** 2026-04-09 validation audit (cat4_plus_1, cat4_plus_2)

---

### ~~NEW-062 (MEDIUM) — `five_prime_rescued` not written to TSV; `correction_applied` omits 5' rescue — Fixed 2026-04-09 (v2.7.9)~~

When `rescue_3ss_truncation` corrects a read's 5' end, `five_prime_rescued=True` is stored in the result dict (bam_processor.py:277) but never written to any TSV column. The `correction_applied` field for cat3 reads shows only `atract_ambiguity,polya_walkback` — the 5' junction rescue is invisible to downstream consumers. Users cannot identify rescued reads from the TSV without comparing `five_prime_position` against raw alignment coordinates.

**Fix:** Add `five_prime_rescued` (bool) to the TSV header and row in `write_results_tsv`, and append `five_prime_rescued` to `correction_applied` when True.

**Files:** `rectify/core/bam/bam_processor.py:567–613` (`write_results_tsv`)
**Discovered:** 2026-04-09 validation audit (cat3 reads)

---

### ~~NEW-063 (MEDIUM) — Chimeric reads with 3'-hard-clipped alignment not exempt from poly-A walkback — Fixed 2026-04-09 (v2.7.9)~~

Reads with XK=1 (chimeric reconstruction) that have a hard-clip at the 3' end are subject to the full polya_walkback pipeline. For cat5_plus_3aligner (100H hard-clip), the unresolved 3' end causes an ambiguity window of 211 bp (positions 9753–9964), producing 5 TSV rows with walkback shifts of 6–185 bp. The correction pipeline has no short-circuit for chimeric reads whose 3' sequence is not present in the query.

**Fix:** In `bam_processor.correct_read_3prime()`, skip polya_walkback (or cap ambiguity window to 0) when the read has a hard-clip at its 3' end (indicates missing sequence; walkback is unanchored).

**Files:** `rectify/core/bam/bam_processor.py` (`correct_read_3prime`), `rectify/core/polya/atract_detector.py`
**Discovered:** 2026-04-09 validation audit (cat5_plus_3aligner)

---

### ~~NEW-064 (LOW) — `netseq_refinement` listed in `correction_applied` when ambiguity_range=1 — Fixed 2026-04-09 (v2.7.9)~~

When `ambiguity_range == 1`, there is only one candidate position; NET-seq refinement runs but is a no-op. Listing `netseq_refinement` in `correction_applied` is misleading — it implies the NET-seq signal was used to resolve ambiguity when no ambiguity existed. Seen in cat6_minus_single (ambiguity_range=1, fraction=1.0, correction_applied includes `netseq_refinement`).

**Fix:** In `bam_processor.correct_read_3prime()`, only add `netseq_refinement` to `correction_applied` when `ambiguity_range > 1` (i.e., when NET-seq signal was actually consulted to break a tie).

**Files:** `rectify/core/bam/bam_processor.py` (correction_applied assembly)
**Discovered:** 2026-04-09 validation audit (cat6_minus_single)

---

### ~~Bug 37 (HIGH) — Zero unit tests for `terminal_exon_refiner.py` — Fixed 2026-04-11 (v0.9.0)~~

The module (1690 lines, multiple coordinate-sensitive code paths including Bugs 29, 33, 35) had no dedicated test file. Highest-priority test gap.

**Fix:** Added `tests/test_terminal_exon_refiner.py` with 51 tests covering: `SpliceSiteIndex` add/lookup/nearby, `load_splice_sites_from_gff` (plus and minus strand coordinate conversion, GFF position conflicts), `merge_splice_indices`, `detect_junction_truncated_reads` (both strands), `detect_partial_junction_crossings`, `get_soft_clip_info`, `simple_align`, `detect_mismatch_clusters`. Real-data class (`TestLoadSpliceSitesFromBundledGff`, `TestWithRealBam`) validates known S. cerevisiae intron positions (YAL030W, YAL001C) against the bundled R64-5-1 GFF and wt_by4742_rep1 BAM.
**File:** `tests/test_terminal_exon_refiner.py`

---

### ~~Bug 38 (HIGH) — `consensus.py` alignment selection only indirectly tested — Fixed 2026-04-11 (v0.9.0)~~

Core per-read consensus selection logic was exercised only through two peripheral tests (`test_xr_flag.py`, `test_gapmm2_seq_restore.py`). Tie-breaking, penalty scoring, and multi-aligner merging had no direct unit tests.

**Fix:** Added `tests/test_consensus_selection.py` with 40 tests covering: `extract_junctions_from_cigar` (multi-intron, soft-clip vs N op), `check_canonical_splice_sites` (GT/AG, GC/AG, non-canonical, real YAL030W intron), `score_alignment` (5' clip penalty −2/base, A-tract depth cap at 10, 3' clip penalty), `select_best_alignment` (winner selection, `was_5prime_rescued` flag, tiebreaker by annotated junction count / 3' agreement, `confidence` levels). Real-data class validates against wt_by4742_rep1 BAM at the YAL030W locus.
**File:** `tests/test_consensus_selection.py`

---

### ~~Bug 41 (MEDIUM) — Trained `--polya-model` is never used — Fixed 2026-04-11 (v0.9.0)~~

`--polya-model` flag, model training pipeline, and `load_model()` infrastructure existed, but `polya_model_path` was captured in config and then discarded. The hardcoded 80% A-richness threshold was always applied regardless.

**Fix:** Wired `polya_model_path` through `bam_processor.correct_read_3prime()` (both streaming and parallel modes). Added `pt_tag`, `polya_score`, and `polya_source` columns to `corrected_3ends.tsv`. Added `rectify tag-polya` subcommand for retroactive annotation of aligned BAMs. Added unaligned dorado BAM auto-detection and tag-preserving alignment pipeline in `preprocess.py`.
**Files:** `rectify/core/bam/bam_processor.py`, `rectify/core/commands/correct_command.py`, `rectify/core/commands/tag_polya_command.py`, `rectify/core/align/preprocess.py`, `rectify/cli.py`

---

### ~~Bug 55 (MEDIUM) — Several APA clustering parameters not CLI-configurable — Fixed 2026-04-11 (v0.9.0)~~

`DEFAULT_MIN_PEAK_SEPARATION = 5`, `DEFAULT_MAX_CLUSTER_RADIUS = 10`, and `DEFAULT_MIN_SAMPLES = 2` had no corresponding CLI arguments.

**Fix:** Added `--min-peak-sep`, `--max-cluster-radius`, `--min-cluster-samples` to `create_analyze_parser()`. Both clustering call sites (single-sample and manifest mode) dispatch to `cluster_cpa_sites_adaptive()` when any non-default value is provided; otherwise fall through to the existing fixed-distance `cluster_cpa_sites()`.
**Files:** `rectify/core/commands/analyze_command.py`

---

## Historical bug fixes

Resolved bugs are no longer tracked here. See `../CHANGELOG.md` for the
release-by-release history.
