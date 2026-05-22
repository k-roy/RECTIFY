# Terminal Junction Peel Plan for Module 2F

Date: 2026-05-21

Status: design note, not implemented.

## Context

Module 2H performs global N-op boundary refinement by scoring candidate
junctions for existing splice `N` operations and rewriting an intermediate BAM.
This is expensive on all-aligner chunk correction.

Module 2F is narrower: it rescues terminal 5' truncation or terminal
misalignment at a candidate 3' splice site. The goal of this plan is to improve
2F's terminal rescue without recreating 2H inside the per-read correction path.

## Core Insight

For terminally misaligned reads, a single inferred rescue boundary is too brittle.
Instead, once a read passes cheap gates, consider each possible peel depth from
the 5' end as a potential boundary. Score all compatible candidate junction
hypotheses plus the current mapping, then pick the best hypothesis.

This handles cases where an early annotated junction gives a mediocre score but
a deeper alternative pooled junction explains the peeled terminal sequence much
better.

## Cheap Gates

Do not peel every noisy 5' end.

Enter terminal peel mode only when both conditions are met:

1. A candidate 3' splice site from annotation or the junction pool is within a
   bounded window of the read's 5' terminal region.
2. The terminal alignment has evidence of a boundary problem, such as a soft
   clip, mismatch, insertion, or deletion within `X` bp of the 5' end.

Read-own N-ops should not be used as primary candidate generators for this
feature. They may be used only as stop boundaries while peeling.

## Peel Generation

Generate candidate peel depths inward from the 5' end until the first of:

- an existing `N` operation is reached
- at least `10` consecutive clean matches are observed
- a hard maximum peel cap is reached

For aligners that emit `M` rather than `=`/`X`, clean matches must be inferred
by comparing aligned query bases to the reference genome or by using equivalent
MD-tag information if available.

Initial conservative defaults to test:

- terminal evidence window: `25-50` bp
- clean anchor: `10` consecutive matches
- maximum peel: `75-100` bp
- candidate 3'SS proximity window: organism-tuned; yeast can start at `10 kb`
  for pooled/annotated candidate discovery, followed by tighter per-candidate
  compatibility checks

## Hypothesis Scoring

For each generated peel depth, score:

- current mapping / no rescue
- each compatible annotation or pool candidate junction
- local ambiguity-shifted versions of each compatible junction

Candidate scoring should compare the peeled terminal read sequence against the
candidate upstream exon sequence. When the bases are currently aligned by the
aligner, also score the current genomic interpretation so rescue must beat the
status quo rather than merely look plausible.

Annotation, pool support, canonical splice-site status, and ambiguity-window
membership should act as priors or tie-breakers, not replacements for sequence
fit.

Require an acceptance margin before overriding the current mapping. Ties should
favor no rescue.

## Non-Goals

This feature must not become Module 2H in disguise.

Do not:

- refine arbitrary internal `N` operations
- use read-own approximate N-ops to search for better N-op boundaries
- rewrite the BAM as a global splice-junction preprocessing step
- score all junctions for all N-ops

This feature should answer only:

> Is the terminal 5' segment better explained as upstream exon sequence across a
> nearby annotated or pooled 3' splice site than by the current terminal mapping?

## Implementation Sketch

Likely touch points:

- `rectify/core/splice/splice_aware_5prime.py`
  - add terminal peel generation
  - add candidate-specific hypothesis scoring
  - preserve current rescue behavior unless the new mode is enabled
- `rectify/core/bam/bam_processor.py`
  - add/query a candidate 3'SS index from annotation plus optional pool
  - avoid passing full candidate sets into the expensive rescue path
- focused tests in `tests/test_splice_junction.py` or validation-read tests
  - early annotated candidate loses to deeper pooled candidate
  - no-rescue tie wins
  - clean terminal alignment does not peel
  - existing `N` stops peel generation
  - `M`-only CIGAR uses genome comparison to identify clean anchors

## Suggested Rollout

1. First run the cheaper 2H bypass tests:
   - no Module 2H, no pool
   - no Module 2H, pool-only for Module 2F candidate discovery
2. Compare ysh1 chunk outputs and runtime against the Plan C profiled run.
3. Implement terminal peel as an opt-in flag only if bypass output loses
   important terminal 3'SS rescue cases.
4. Keep counters/profiling for:
   - reads passing candidate gate
   - reads passing terminal evidence gate
   - peel depths evaluated
   - candidates evaluated
   - rescues accepted
   - no-rescue wins

