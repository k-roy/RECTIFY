# bundle — `validation_reads_dorado_source.bam` is stale

**Status**: RESOLVED — commit `c3202f6` on branch `drs-validation-rebuild`
(2026-05-19). Verification passes: 0 missing reads. Both the immediate
staleness and the structural cause (regen_pa_rest_bundle.py skipping the
dorado-source step) are fixed. See resolution notes below.

**Filed**: 2026-05-19 plotter session.

---

## Problem

`rectify/data/validation/validation_reads_dorado_source.bam` is dated
2026-05-10, while the rest of the bundle (`validation_reads.bam`,
`rectified/per_aligner/*.bam`, etc.) was regenerated on 2026-05-18 at
21:50. The dorado-source BAM was not updated when the bundle was rebuilt,
and as a result it is missing 7 of the 36 validation reads:

| XV tag         | qname (first 8) |
| -------------- | --------------- |
| cat1_minus_1   | 77b392d9        |
| cat1_minus_2   | 34ba198b        |
| cat1_plus_1    | 0cb5a111        |
| cat1_plus_2    | a146838d        |
| cat2_minus_1   | b313b50d        |
| cat2_plus_1    | 61b0c014        |
| cat2_plus_2    | 88953e9c        |

(All four cat1 reads and three of the four cat2 reads; cat2_minus_2 is
the only cat2 read still present.)

## Why this matters

The "minimap2 (unrectified)" row in `render_read_alignment.py` is the
visual answer to the question "what was wrong before RECTIFY?". The
dorado-source BAM is minimap2's alignment of the **full untrimmed read**,
so the body typically extends 1–3 bp into the genomic A-tract before
soft-clipping the rest of the poly-A tail. That extension — and the
resulting drift of the naive 3'-end position vs. the corrected one — is
the smoking gun.

When the dorado BAM is missing a read, the renderer falls back to
`rectified/per_aligner/minimap2_unrectified.bam` (the spliced version,
which aligns the *trimmed* body and re-adds the pA as a soft-clip). The
body never extends into the A-tract there, so the boundary-extension
behaviour is invisible. We see "minimap2's body looks the same as the
winner" instead of "minimap2 over-extends into the tail."

For the 7 reads above, the per-base view currently shows the spliced
fallback. Once the dorado BAM is regenerated to cover these reads, the
visualization will switch back to the faithful source automatically (no
plotter change needed).

## What the debugger needs to do

Regenerate `validation_reads_dorado_source.bam` so it includes the 36
reads currently in `validation_reads.bam`. The builder lives at
`scripts/validation_data/rebuild_2026_05/build_combined_dorado_source.py`
— that path matches the file that produced earlier dorado_source BAMs.
A plain regen-from-current-FASTQ should be sufficient; the prior staling
appears to be because the bundle regen step in
`regen_pa_rest_bundle.py` does not touch the dorado-source BAM.

Two options for keeping this from going stale again:

1. **Add a dorado-source regen step** to `regen_pa_rest_bundle.py`, so
   the dorado BAM is rebuilt every time the bundle is. Lowest-friction.
2. **Add a freshness check** at the top of the bundle regen: if
   `validation_reads_dorado_source.bam`'s mtime is older than
   `validation_reads.bam`'s, fail loudly with a remediation hint.

Either is fine. Option 1 is the safer default; option 2 is fine if
regenerating dorado is expensive (it needs the original untrimmed
FASTQs).

## Plotter workaround (already shipped, WIP)

`scripts/validation_data/render_read_alignment.py` now probes the dorado
BAM per-read and falls back to `minimap2_unrectified.bam` when the read
is missing. This guarantees the "minimap2 (unrectified)" row always
renders, but for the 7 affected reads it shows the less-faithful spliced
alignment. No further plotter change needed once the dorado BAM is
regenerated.

## Verification

After regen, this should print nothing missing:

```python
import pysam
qmap = {r.get_tag('XV'): r.query_name
        for r in pysam.AlignmentFile(
            'rectify/data/validation/validation_reads.bam', 'rb')
        if r.has_tag('XV')}
dorado_qnames = {r.query_name for r in pysam.AlignmentFile(
    'rectify/data/validation/validation_reads_dorado_source.bam', 'rb')}
missing = [(xv, q) for xv, q in qmap.items() if q not in dorado_qnames]
assert not missing, missing
```

**Result (2026-05-19):** `Total: 36  In dorado: 36  Missing: 0` ✓

---

## Resolution notes (2026-05-19)

**Root cause:** `regen_pa_rest_bundle.py` never touched the dorado-source
BAM. When v3.2.5 replaced Cat1/Cat2 reads, the bundle was regenerated but
`validation_reads_dorado_source.bam` was left at its v3.2.4 state (dated
Apr 24), containing the old Cat1/Cat2 qnames. The 7 new qnames were absent.

**Fix:**

1. **`scripts/validation_data/build_dorado_source.py`** (new) — rebuilds
   `validation_reads_dorado_source.bam` from scratch: reads XV/XG from
   `validation_reads.bam`, extracts each read from `wt_by4742_rep1.bam`
   (primary mapped only), falls back to the DRS-trim minimap2 merged BAM
   for reads unmapped in the Dorado source (cat1_plus_1 / cat1_minus_1).
   Verifies count, XV tags, and qname coverage.

2. **`regen_pa_rest_bundle.py` Step 5** (added) — calls
   `build_dorado_source.py` at the end of every bundle regen, so the
   dorado-source BAM is always regenerated alongside `validation_reads.bam`.

**CIGAR note:** Dorado and the minimap2 merged BAM both emit `M` ops
(not `=/X`). The "M-conversion" trap from the v3.2.5 CAUTION does not
apply here — the sources are consistent.

**Commit:** `c3202f6` on branch `drs-validation-rebuild`.
