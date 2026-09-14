# Consensus and micro-exon provenance follow-up — 2026-09-14

Fable reviewed and fast-forwarded the initial nine Codex audit commits to master
at `d7d92f7`, with an independent 517 passed / 1 skipped check. This follow-up
implements the consensus and micro-exon provenance items assigned back to Codex.
Fable retains ownership of resolver B2/B3 validation and 2H candidate realizability.

## CFX-10: canonical credit reaches actual selection

`consensus/extract.py` and chimeric segment scoring now share a paired motif
check. The genomic pairs are GT-AG / GC-AG / AT-AC on plus and CT-AC / CT-GC /
GT-AT on minus. Chimeric selection passes the read strand into the segment
scorer; sequence-equivalent placements retain the same canonical credit.
Independent donor and acceptor sets would admit AT-AG and GT-AC incorrectly.

The synthetic whole-selection regression has two candidate CIGARs:
`40M100N60M` at a canonical junction and `44M100N56M` at a noncanonical junction.
Before the fix, plus AT-AC and all three minus classes lose in both input orders.
After the fix, the emitted CIGAR preserves the canonical junction. Tests also
cover the ordinary consensus selector, mixed-pair refusals, both orientations,
ambiguous coordinates, and invalid bounds. These are mechanism tests, not an
estimate of biological FN/FP rates. The fix is commit `44d5ee9`.

## CFX-11: provenance records successful edits at specific junctions

The old writer copied planned micro-exon calls from the TSV into XB/XV, even
when live surgery refused a call. XB already carries cDNA strand-split counts;
XV labels validation reads. Aggregation treated any nonempty XB as evidence
that every intron on the read was drawn by Station B.

New output uses **Xb:Z**, a versioned compact JSON object. XB and XV are
preserved. Example, in genomic zero-based half-open coordinates:

```json
{"v":1,"chrom":"chrT","strand":"+","calls":[{"intron":[60,460],"exons":[[200,206]],"alternatives":"chrT:300-306"}]}
```

This call creates exactly `[60,200)` and `[206,460)`. It gives no credit to a
distant original intron. `alternatives` uses the existing format: commas join
exons within a configuration; semicolons separate equally good configurations.
The call is appended only after its live CIGAR rewrite succeeds. If just one of
several calls succeeds, only that call and its alternatives are recorded.
Repeating the writer preserves previous successful calls without duplicating them.

A call may also contain **`junctions`**, a nonempty subset of the new introns:

```json
{"intron":[60,460],"exons":[[200,206]],"alternatives":"chrT:300-306","junctions":[[60,200]]}
```

This is necessary when chimeric selection retains only part of a configuration.
The writer takes provenance from the source of each selected N-op, intersects
it with the final CIGAR, and records that subset. Copying Xb from the sequence
donor would attach the wrong history even when donor and winner have identical
coordinates. Fallback selection uses the winning placement's provenance.
Later clipping/refinement also cannot transfer credit to different coordinates:
aggregation matches only the surviving exact N-ops.

`build_chimeric_read` accepts the source `aligner_reads` map for this purpose;
the production consensus caller supplies it. External callers that construct a
multi-source result must supply source events and reads to retain this provenance.
Xb is **placement-specific**, never a read-intrinsic tag to broadcast across
aligner arms. Fable owns the corresponding CMA whitelist/documentation review.
Untagged reads take a fast path before any extra segment/CIGAR traversal.

## Migration and interpretation

- `station_b_reads` counts verified Xb draws for that exact current junction.
- `station_b_alternatives` reports alternatives belonging to those calls only.
- `station_b_unverified_reads` counts supporting reads whose legacy or invalid
  provenance leaves their micro-exon history unresolved. It is a read-level
  uncertainty flag repeated on the read's junction rows, not a drawn-junction count.
- Legacy XB coordinate strings receive the unverified flag and a warning to
  reprocess original alignments. cDNA `n_top/n_bottom` XB values do not count.
  Unknown Xb versions, malformed calls, or contig/strand mismatches receive no
  verified credit. No automatic legacy-to-verified conversion is possible:
  the old tags did not record whether the writer actually succeeded.
- Missing Xb alone does not prove that an aligner supplied independent support.
  Other correction stages and legacy files may lack these tags. Recount exact
  destination support on original BAMs with `scripts/recount_junction_support.py`.

## Validation

- Initial red reproduction: 10 failed / 8 passed among extraction and actual
  whole-selection cases, before the canonical fix.
- Canonical change and neighboring tests: 100 passed / 4 existing data skips.
- Provenance/writer/aggregation/chimeric handoff: 61 passed; existing matplotlib
  deprecation warnings only.
- Broader consensus and BAM-writer checks: **315 passed / 4 existing data skips**,
  18.71 seconds; existing matplotlib deprecation warnings only. Log:
  `dev/audits/codex_fable_20260913/followup_consensus_provenance_tests.log`.
- New helper/tests pass Ruff's runtime-error and unused-import checks; diff
  whitespace checks pass. One local pytest process tree at a time, numeric
  worker counts capped at one, no subagents or cohort jobs.

No annotation threshold, homopolymer exemption floor, or resolver search ceiling
was changed. Cohort-level sensitivity/specificity still needs Fable's replay work.
