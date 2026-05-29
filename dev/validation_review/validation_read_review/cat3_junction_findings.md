# cat3_junction — findings surfaced by the End-correction panel

Written by the **plotter session** for the **debugger session** to consume.

**Source bundle**: `rectify/data/validation/rectified/` regenerated 2026-05-18
with `=`-SEQ decoded + HP-ED winner selection.

---

## cat3_minus_1 — clean (positive control)

Reviewed by user 2026-05-18: looks correct as rendered. Use as a diff
target when investigating the cat3_minus_2 / cat3_plus_2 issues below.

---

## Common issue across cat3_minus_2 + cat3_plus_2: 2D-Nop CIGAR instead of clean N-op

**User's diagnosis (verbatim, 2026-05-18, on cat3_minus_2)**:

> "We should be capturing the correct intron with a clean N-op, rather
> than 2D-Nop. All we need to do is move the AG aligned bases to exon 1.
> Also, it seems our 3' pileup signal is stale."

> "Same issues for cat3 plus2."

> "gapmm2 and mpb actually get the intron mapped cleanly, without the
> 2D80N issue. They get the proper splice junction (GT-AG) with a clean
> 82N."

**This is the key debugger insight from cat3_minus_2**: two of the five
aligners (gapmm2, mapPacBio) emit the canonical `82N` CIGAR — donor and
acceptor land on the proper GT/AG motif. The other three (deSALT,
minimap2, uLTRA) emit `2D 80N` — the 2 bp of "AG" at the start of the
downstream exon get charged as a deletion instead of being absorbed into
exon 1, shifting the donor 2 bp upstream of the canonical site.

**Hypothesis for the debugger**: the junction refiner (Module 2H) IS
producing the right answer for gapmm2 and mapPacBio because those
aligners' native CIGAR around this junction already aligns the AG bases
into exon 1 — the refiner's "no change needed" path triggers. For
deSALT/minimap2/uLTRA, the native CIGAR puts the AG bases at the start
of exon 2 (right after the intron), and the refiner's "junction slide"
detection isn't firing for this 2-bp case.

The clean-N-op aligners' result is the **correct biological answer**;
the debugger should look at why the slide isn't firing on the 2-bp case
for the other three aligners. Candidate: the slide fast-path in
`junction_refiner.py` may only handle 1-bp equivalences (per the
`feedback_rectify_junction_slide` memory) — extending it to k=2,3,...
should resolve.

The Phase 2 painted-CIGAR grouping (rendered 2026-05-18) makes this
visible at a glance: gapmm2+mapPacBio collapse into one group with a
clean grey N-block; the other three form a second group with an orange
2D segment immediately before the N-block.

### Why don't gapmm2/mpb win this read?

Per `corrected_consensus.merge_corrected_tsvs`, the HP-ED scoring DOES
correctly penalize `2D80N` over `82N`:

| aligner   | CIGAR around junction | HP-ED  | _chimera_ok | _five_rescued | _effective_chimera_ok |
| ---       | ---                   | ---    | ---         | ---           | ---                   |
| gapmm2    | `… 82N …`             | **20.63** (best) | 1 | 0 | **1** (penalized) |
| mapPacBio | `… 82N …`             | 20.83  | 1           | 0             | **1** (penalized)     |
| deSALT (winner) | `… 2D 80N …`    | 21.49  | 1           | 1             | 0 (exempt)            |
| minimap2  | `… 2D 80N …`          | 21.49  | 1           | 1             | 0 (exempt)            |
| uLTRA     | `… 2D 80N …`          | 21.49  | 1           | 1             | 0 (exempt)            |

The 0.86 ED gap between `82N` and `2D80N` is exactly the cost of the 2
extra deletions (2 × `del_cost(hp_context, base)` ≈ 2 × 0.43). The
scorer is doing the right thing — `82N` is correctly ranked above
`2D80N`.

**The problem is the chimera-filter exemption.** All 5 aligners trip
`_chimera_ok=1` (short 5' overhang at the junction). The exemption
rule:

```python
_effective_chimera_ok = _chimera_ok & (_five_rescued == 0)
```

…protects only 5'-rescued aligners from the chimera penalty. gapmm2
and mapPacBio aligned **cleanly through the intron without needing
rescue** (their native CIGAR already had the 82N), so the exemption
doesn't apply. They remain in tier 2 even though they have the lowest
HP-ED and the cleanest CIGAR.

**Suggested debugger fix**: extend the chimera exemption to apply to
any aligner whose junction matches an annotated splice site with a
canonical motif (GT-AG or GC-AG), regardless of whether 5' rescue
fired. With that exemption, the sort would land gapmm2 (lowest HP-ED
with proper junction) as the winner — the biologically correct outcome.

A simpler partial fix: if ALL 5 aligners share the same junction
coords and the junction is in the annotated set, exempt all of them.
That removes the 5'-rescue bias for the well-supported-junction case
without changing behavior elsewhere.

This finding also resolves the upstream "2D80N CIGAR" symptom: the
debugger doesn't need to extend the slide fast-path to handle 2-bp
equivalences if the chimera-exemption fix simply lets the correct-CIGAR
aligners win. The slide-fast-path extension is still desirable
(produces the cleaner CIGAR in the OTHER three aligners' corrected
output), but it's not on the critical path to fixing this read's
winner.

**What the rendered alignment rows show:**

In both cat3_minus_2 and cat3_plus_2, the corrected track displays:
- The end of upstream exon 1
- `--` (orange) representing a 2-bp **deletion** in the alignment
- Then `||||...||||` (the N-op intron block)
- Then the downstream exon starting with `AG...`

So the CIGAR pattern is `...M 2D N{intron_len} AG...` instead of the
cleaner `...M+2 N{intron_len-2} ...` where the `AG` bases are absorbed
into the upstream exon and the N-op shifts left by 2.

**This is the "junction slide" pattern** noted in
`feedback_rectify_junction_slide.md`:

> When `genome[old_donor] == genome[old_acceptor]`, the canonical 1-bp
> slide is a trivial = op length adjustment (e.g. `8=N51=` → `9=N50=`),
> NOT a local realignment with new I/D ops. Refiner currently does the
> wrong thing for this case.

Here the slide is 2 bp instead of 1 bp, but the same principle applies:
the genomic donor and acceptor positions are equivalent (both flanked by
the same bases), so the CIGAR should be allowed to slide AG into the
upstream exon and shorten the N-op accordingly. Module 2H's refiner is
producing `2D N` instead of the parsimonious `=+2 N(short)` form.

**Hypothesis for the debugger**: extend the junction-slide fast-path in
the refiner (`rectify/core/splice/junction_refiner.py`) to handle 2+ bp
equivalence cases, not just 1-bp. The decision criterion is
`genome[old_donor:old_donor+k] == genome[old_acceptor:old_acceptor+k]`
for k=1, 2, 3, …; when the equivalence holds, emit a length-adjustment
CIGAR fix rather than a local realignment that produces `kD N`.

---

## Stale 3' pileup bedgraph

**User's diagnosis**: "Our 3' pileup signal is stale."

**Confirmed**: `regen_pa_rest_bundle.py` regenerates the corrected BAMs
and `corrected_reads.tsv`, but it does NOT regenerate the bedgraphs at
`rectify/data/validation/rectified/corrected_3ends.{plus,minus}.bedgraph`.
After winner-selection switched aligners during the HP-ED regen, the
per-read corrected_3prime values changed → bedgraph counts at the OLD
positions are now stale.

**Plotter-side fix candidate**: add a bedgraph regen step to
`regen_pa_rest_bundle.py`. The corrected_reads.tsv has the per-read
corrected_3prime and strand columns, which is everything needed to
recompute the bedgraphs. The existing `rectify analyze` command, or
`rectify/core/bam/bedgraph_writers.py`, has the writer logic.

This is borderline plotter/debugger. The fix is mechanical (add another
step to a regen script), but it touches the bundle's data contract.
**Awaiting user direction** on whether I should ship the regen-script
fix or queue it for the debugger.

---

## Panel display bug fixed in this session

While investigating, noticed the panel was rendering wrong intron coords
because the summary TSV's `junctions` column has a related bug:

```
deSALT     ...  five_prime_position=142252  junctions=1   ...   five_prime_rescued
gapmm2     ...  five_prime_position=142252  junctions=1   ...   five_prime_rescued
mapPacBio  ...  five_prime_position=142229  junctions=142253-142619  ...  none
minimap2   ...  five_prime_position=142252  junctions=1   ...   five_prime_rescued
uLTRA      ...  five_prime_position=142252  junctions=1   ...   five_prime_rescued
```

For aligners where 5' rescue fired, the summary TSV's `junctions` field
contains literal `1` (presumably a count of rescued junctions) instead
of the proper `donor-acceptor` coord pair. mapPacBio, which didn't
rescue, carries the proper coords.

**Plotter workaround applied**: when the winner has empty corr_juncs
but other aligners have non-empty, the panel falls back to the longest
non-empty list for the intron column structure. This keeps the columns
populated correctly even with the upstream data bug.

**Debugger ticket**: the summary TSV writer (in `merge_corrected_tsvs`
or its helper) should always write the `donor-acceptor` coord pair for
the `junctions` column, regardless of how the junction was arrived at
(direct alignment vs. 5' rescue). The current behavior loses
information — for these 4 aligners we now know "a junction was
rescued" but not WHICH one.

---

## cat3_plus_1

Not flagged by user this turn. PNG in `cat3_junction_review_pngs/`.
Per the existing PLOTTING.md note: cat3_plus_1 has a known issue where
the corrected BAM has no N-op despite the TSV claiming rescued=1 — that
may also be related to the same TSV-serialization bug.

---

## Update — 2026-05-18, formal cat1–3 review (post Phase A/E/cat3 cleanup)

The freshly regenerated bundle resolves the cat3_minus_2 / cat3_plus_2
2D-N issue (debugger commits `6d2cf59` + `0653172`). The remaining cat3
finding from this review is **mapPacBio's 5'-end rescue failure**.

### Finding: mapPacBio's 5'-end mismatch/indel cluster blocks 5'-rescue

**Reads affected**: cat3_minus_1 and cat3_plus_1 (both with mapPacBio in
its own non-winning group despite being structurally aligned through
the junction).

**What the overview panels show**:

- cat3_minus_1: mapPacBio's overview row starts (RNA-5' side) with a
  tight cluster of `+1 +2 +5` (purple insertions) right before the
  long body alignment. The winner group (deSALT + gapmm2 + minimap2 +
  uLTRA) doesn't have this cluster — they start clean and the body
  flows into the canonical junction.
- cat3_plus_1: mapPacBio shows `+1 +1` at the 5' edge where the winner
  group shows a single `+1`. Looks like a 1-bp encoding / off-by-one
  difference inside the rescue logic — small but enough to push
  mapPacBio out of the winner cluster.

**User's diagnosis (verbatim, 2026-05-18)**:

> "For cat3 minus1, I would have thought that mapPacBio's mismatch
> heavy 5' end would have been rescued. We should be showing insertions
> in the 5' end/junction base-level panels as well."

> "For cat3 plus1 it is clear that the rescue for mapPacBio has a
> simple-off-by-one issue that prevents it from getting into the
> winner's section. In general, it seems that the rescue step is
> having a hard time identifying what section of mapPacBio's 5' ends
> to attempt rescue with."

**User's proposed fix (verbatim)**:

> "I thought it would be straightforward to essentially turn mpb's
> alignment into a 5' soft clip by identifying where the mismatch/indel
> heavy portion ends, where a long string of = starts, and treat
> everything upstream as a soft-clip for the rescue step."

### Suggested debugger algorithm

In the 5'-rescue pre-pass (before `extend_read_5prime_for_junction_rescue`
/ `reroute_intronic_tail_5prime_via_junction` decides what to do with
mapPacBio's read), insert a **5'-anchor-finding** step:

1. Walk the CIGAR from the 5' edge (= start of CIGAR in BAM order for
   `+` strand reads, or end of CIGAR for `−` strand reads).
2. Identify the **first sustained = run** of length ≥ N (suggest N=10
   as a starting point). This is the "clean anchor" — the point where
   mapPacBio's alignment becomes well-supported by the reference.
3. Reinterpret everything upstream of that anchor as a leading
   soft-clip. The mismatch/indel/= ops in that messy region become a
   single `kS` op where k = total query bases consumed upstream.
4. Hand off to the existing 5'-rescue machinery, which now sees a
   clean `kS body_align…` and can do its junction-rescue normally
   (matching mapPacBio's geometry to the other aligners that already
   produced a clean leading soft-clip).

Pseudo-code:

```python
def reanchor_mpb_5p_for_rescue(read, anchor_min_run=10):
    """If the 5' edge has a mismatch/indel cluster followed by a sustained
    = run, collapse the cluster into a leading soft-clip so the 5'-rescue
    machinery sees the same geometry as the other aligners."""
    cigar = read.cigartuples
    # ... walk from 5' edge of CIGAR (strand-aware), find first =/M run of
    # length >= anchor_min_run, collapse everything upstream into S, then
    # re-write read.cigartuples + read.reference_start.
```

This would also help any other aligner whose native CIGAR happens to
emit a messy 5' edge (gapmm2 occasionally does the same on Cat3-style
reads). The cost is one CIGAR walk per read in the rescue path.

**Bonus**: with mpb in the winner cluster, the overall HP-ED ranking
gets cleaner, and the `effective_group` column would consolidate even
further (likely all 5 aligners in cluster A on these reads).

### Plotter-side related observation

The base-level junction-zoom panels (donor-side + acceptor-side
stub views) do NOT currently render insertions. The painted-CIGAR
overview shows the `+1 +2 +5` purple ticks, but when the user
zooms in to the junction the insertions disappear. Tracking as a
plotter-scope todo: add I-op rendering to `render_zoom_alignment_row`
(probably the same purple-tick treatment used in the per-base
alignment row, just scaled to the zoom panel's narrower window).
