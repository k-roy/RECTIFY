# Junction accuracy and efficiency audit — 2026-09-13

The priority is to make placement evidence reliable before changing scientific
thresholds. This audit found order-dependent loss of original candidates,
incomplete micro-exon enumeration, corrupted junction evidence, and inconsistent
canonical-motif handling. Small, reproducible fixes are separated from proposed
changes that still need read-level validation.

## Revision and scope

Initial source baseline: `origin/master` `0df276a`. Integrated Fable's ISSUE-031
commit `af2f788` before the final audit tests. The original checkout and its
pre-existing documentation edits were preserved; work uses
`audit/dev-queue-stations-20260913` in an isolated worktree.

All 504 tracked Python files received a Python 3.8 grammar check and the runtime
error subset of Ruff checks. Manual review concentrated on align/resolver,
consensus/triage, splice/refinement, BAM read/write, cDNA analysis, aggregation,
provenance, visualization, and benchmark evidence. This is not a claim that every
line was manually inspected or that cohort FN/FP rates have been measured.

Fable leads development and landing. Codex supplies independent reproductions,
review verdicts and isolated fixes. No additional agents were spawned. Test runs
are serial, with one pytest process tree at a time and numerical-library threads
capped at one on the M1.

## Implemented corrections

| ID | Defect and effect | Correction / evidence |
|---|---|---|
| CFX-01 | A high-confidence corrected read was compared with only the first original arm. A better later arm could never restore it. | Compare the minimum eligible original HP-ED across arms before the guard. A same-sequence fixture with corrected 43.75, originals 43.75 and 10.00 now restores the same CIGAR in either arm order. |
| CFX-02 | Micro-exon search returned a prefix when its candidate limit was reached; an omitted candidate could change the winning split or ambiguity. | Search one beyond the limit, refuse incomplete enumeration, preserve exactly-at-limit sets. Two bisects also replace chromosome-tail list copying. Overflow refusal is conservative; its cohort sensitivity cost remains unmeasured. |
| CFX-04 | Squished pileup rendering referenced an undefined poly(A)-source variable. | Propagate clip/pt choice; both-strand rendering controls. |
| CFX-05 | SEQ `=` counted as mismatch; every intron inherited the read's combined near-junction errors. | Respect reference encoding, use local near windows per junction, retain shared read-body background. An error at one intron no longer contaminates its clean sibling. Monotone edge lookup avoids scanning all introns at each base. |
| CFX-06 | Single-arm pool construction doubled unspliced support; process-pool fallback omitted optional signal arguments. | Separate worker and aggregate counters; keep evidence arguments in fallback. Cache format 3 rejects previous evidence. |
| CFX-07 | `cdna-analyze` counted secondary/supplementary records as extra molecules and assigned them primary-placement tags. | Primary-only molecule intake; retain other alignments with their existing tags; report skipped non-primary counts. Integration control verifies three molecules from six alignment records. |
| CFX-08 | Aggregate canonical filtering omitted AT-AC despite promising it. | Paired motifs on both strands; AT-AG, GT-AC and GC-AC negative controls. Removed an unused read-name list per junction observation. |
| CFX-09 | Atlas loading imports PyYAML, but neither pip nor conda runtime metadata declared it. | Declare PyYAML in both manifests; the failing atlas test passes in the existing environment that supplies it. |

Two developer tools were also repaired: the indel-correction figure script's
missing `Path` import and the validation HTML generator's Python-before-3.12
f-string syntax error.

## Remaining junction priorities

### 1. Select among candidates that can actually be emitted

The ISSUE-031 guard prevents harmful new I/D beside an N, with the approved
homopolymer insertion exception. However, `refine_read_junctions` retains only
its best proposal. A later surgery refusal does not try the next candidate.

A controlled-score witness on `af2f788` proves the mechanism:

```
input:          20M100N5D20M, intron [120,220)
best proposal:  acceptor +2 → forbidden remaining 3D → refused
runner-up:      acceptor +5 → 20M105N20M → realizable
actual output:  unchanged input; runner-up never attempted
```

This is a real loss-of-candidate mechanism, not just inaccurate proposal counters.
Its frequency in biological reads is not established. Prefer a non-mutating
realizability probe before ranking, plus a final output check. Test interactions
between multiple edits and distinguish proposed, realizable, applied and refused
counts. Do not weaken the new I/D invariant to make application counts larger.

The HP exception needs independent literal-sequence controls: run length 1 versus
2, mixed-base insertions, both genomic flanks, reverse-strand geometry, and
literal versus `=`-encoded input with qualities preserved. An assertion that
uses the production helper to excuse its own outputs is insufficient alone.

For a fixed inserted base and independent uniform flanks, the chance that either
flank has a touching run of length at least 2 is `31/256 ≈ 12.11%`; at least 3 is
`127/4096 ≈ 3.10%`. These toy nulls are not FP rates among selected proposals.
Raising the floor to 3 also excludes approved AA→AAA cases. Calibrate on matched
negative controls before changing it. Either flank can explain an insertion in
the spliced sequence: I and N consume different axes, so their order in CIGAR
alone does not identify a biological side.

### 2. Make canonical credit consistent throughout consensus

The aggregate-table fix does not repair upstream selection:

- `consensus/extract.py::check_canonical_splice_sites` omits AT-AC.
- `chimeric_consensus.py::_canonical_within_window` checks only forward GT/GC-AG;
  `score_segment` receives no strand. Minus-strand canonical junctions can receive
  -3 instead of +5 per junction, as can minor-class junctions.

Propagate strand and use paired transcript-oriented motifs. Verify whole
selection on both strands, including biologically noncanonical negative controls
and ambiguity-equivalent placements. Retain default-on AT-AC in the resolver and
2H scorer, where support already exists.

### 3. Keep population support independent of the correction

Use `scripts/recount_junction_support.py` with explicit move targets and ORIGINAL
BAMs. It streams each arm once, unions primary QNAMEs, excludes the moving QNAME,
and reports exact-coordinate support and a separate contiguous M/= anchor count.
It does not infer molecule independence, sequence identity or biological truth.
Memory scales with requested-destination supporters, not all BAM records.

Do not count corrected co-movers as corroboration. Also fix the separate Station-B
provenance issue: micro-exon `XB` collides with cDNA strand-split `XB`, and aggregation
credits all introns on an XB-tagged read, including distant original introns.
Use a distinct tag with an explicit migration policy and attribute each created
junction to its actual recovered exon configuration.

Older Station C scores only the first 50 eligible reads per junction for `q_max`
and `q_2nd`. These are censored, order-dependent statistics. Report saturation and
compare deterministic or complete evidence collection before treating the values
as a calibrated gate. A collector that returns partial counters after a BAM read
error also needs an explicit incomplete-state/cache policy.

### 4. Finish the resolver's evidence and ambiguity checks

A11's description is stale: Case B2 and B3 already call `_score_alts`. The
remaining B1 D-to-N selection takes its first best-scoring candidate without an
ambiguity check. Distinguish harmless equivalent D/N representations from
inequivalent tied splice placements before changing this path.

**CFX-03 is reproduced in B2, without scoring mocks.** A synthetic 360M read
contains 200 bases matching its left flank, 40 bases matching a distant canonical
acceptor, then 120 bases matching its original linear continuation. Production
`resolve_read` emits `200M300N160M`: the local 40-base window wins, but the entire
160-base block moves. Literal base-by-base mismatches rise from **30 to 94**, as
does emitted-CIGAR HP-ED. A separate CIGAR walker confirms the mismatch counts.

Validate the entire block whose coordinates change before accepting a B2/B3 move.
A linear-time emitted-CIGAR non-regression check is one candidate; if gapped
rescue is intended, explicitly score and emit that whole alignment. Preserve
both-strand genuine-rescue controls. Cross-donor/acceptor candidates also need
comparable evidence windows and ambiguity checks, not just the lowest local raw
score. The B3 mirror remains source-reviewed; biological incidence is unmeasured.

The two-sided candidate hook (A6) is still intentionally unconnected. The new
micro-exon pass is annotation-driven; de-novo rescue needs a repeat-aware null and
read-level controls, not merely more candidate motifs.

### 5. Measure operating points before changing defaults

Keep the resolver ceiling at 2000. The existing yeast A/B admitted the previously
unassessed candidates at 20000 with byte-identical BAM output and +49% wall time;
that is evidence about that dataset, not a universal claim about abandoned clips.
Human and yeast require separate abandonment/rescue curves.

A2's literal exponent-cap replacement is not repeat calibration: below the
configured maximum window, the arbitrary exponent ceiling is not the active
bound. Measure per-locus repeat/sequence informativeness and the number of true
placements excluded at each operating point.

A17's terminal poly(A)-to-genomic-A false junction needs identity and anchored
sequence evidence, not an arbitrary minimum exon length. Include genuine short
terminal exons and intron-retention controls before adding a veto.

Use read-level truth restricted to introns each read actually spans. Report exact
coordinates and independently derived sequence-equivalent matches separately;
±6–8 bp is a sensitivity analysis, not permission to merge real nearby isoforms.
Stratify precision/recall by motif class, annotation, repeat context, exon length,
strand, protocol and aligner. Verify both CIGAR and reconstructed sequence/quality
invariants after every stage.

## Queue and validation state

Git ancestry confirms these old “branch only” fixes are already on master:
A4 (`cb1469e`), A7 (`c85da29`), A9 (`1b8031d`), A8 Case-B support (`c33d693`),
C3 (`c431fdc`), D2 (`6fc4ef5`), D3 (`f266293`), D5 (`4dfbc67`), and both D8
repairs (`5ed4bd5`, `5e688a2`). Chimeric placement identity, Type-2 no-collapse,
parallel cDNA QC and Path-A carried tails are also landed; stale open entries were
removed from `KNOWN_ISSUES.md`.

The unrelated motif-analysis queue still needs independently gated enriched and
depleted scans; one small group currently silences the other. Motif logo plots
and chimeric-reconstruction evaluation remain separate work.

Validation: 94 earlier focused tests passed; 54 further tests passed with one
fixture skip after integrating `af2f788`. The controlled runner-up witness passed.
The full non-slow suite completed in **486.90 seconds: 3,699 passed, 50 skipped,
3 xfailed, 4 deselected; two failures and two setup errors**. The failures were
the blocked hardware query and missing PyYAML; both targeted tests subsequently
passed in the existing conda environment outside the hardware-query sandbox.
The two pre-existing setup errors require an absent uncompressed yeast FASTA and
validation trim metadata. A final constant-time header-membership optimization
in triage was followed by **11 passing** triage tests. No new pipeline regression
was found; this is not an assertion that a fresh clone has an all-green suite.

The full suite needs five pre-existing untracked benchmark/calibration helpers;
these were copied into the isolated worktree for testing and were not included in
the changes. Their hashes are recorded in the local coordination log. A clean
clone's test reproducibility remains an existing repository issue.

The live M1-detection test fails when the sandbox denies `sysctl`; the same query
succeeds outside it. This directly explains an environment failure without
asserting leaked test state. Prefer hermetic arm64/Rosetta controls and keep any
live hardware test optional. The remaining 14 Ruff findings comprise 13 type-only
names and one unused test helper; the syntax sweep reports zero errors.
