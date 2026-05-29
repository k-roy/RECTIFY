# cat1_indel — walkback findings surfaced by the End-correction panel

Written by the **plotter session** for the **debugger session** to consume.
Plotter does not own walkback / indel_corrector / scoring code; these are
observations only, with hypotheses for the debugger to validate.

**Source bundle**: `rectify/data/validation/rectified/` regenerated 2026-05-18
with `regen_pa_rest_bundle.py` (HP-ED winner selection, `=`-SEQ decoded at
bam_writer.py).

**Source observations**: rendered PNGs at
`validation_read_review/cat1_indel_review_pngs/`. Panel data driven by
`rectify/data/validation/rectified/per_aligner_summary.tsv`.

---

## Policy reminder (per user, 2026-05-18)

> Reads can NEVER end in A in RECTIFY, as per policy we always walkback to
> the first non-A due to the inherent ambiguity at the genomic-A /
> poly(A)-tail boundary.

This is stronger than "leftmost-possible-CPA" — it's an *absolute terminal
base constraint* on the corrected 3' end. A corrected position must point
at a non-A on the RNA strand.

---

## cat1_minus_1 (77b392d9, chrII, − strand) — under-shoots by 1 bp

**Genome at the corrected position:**
```
chrII pos:   9825      9831       9840
plus seq:    AAAAAA A T  CTCTTTCACT
RNA  seq:    UUUUUU U A  GAGAAAGUGA
(RNA = revcomp of plus)
```

For minus-strand reads, RNA-3' is at the LOWER genomic coordinate.

- `corrected_3prime = 9831`. chrII[9831] = `T` on plus strand → `A` on RNA.
- **The corrected RNA-3' base IS A** — violates the "never end in A" policy.
- The next non-A walking in RNA-3' direction (= increasing genomic coord)
  is chrII[9832] = `C` plus → `G` on RNA. Walkback should land there.

**Per-aligner summary (HP-ED panel)**: all 5 aligners agree on corr_3p =
9831. So the bug isn't aligner-specific — it's in the walkback termination
condition that lets the algorithm stop on an A-on-RNA position.

**Hypothesis**: the walkback's "is this position non-A" check is operating
on the genome-plus-strand base rather than the RNA-strand base for minus
reads. Or the check uses `read_base == ref_base AND ref_base != 'A'` but
the ref slice is post-revcomp / pre-revcomp inconsistently for minus strand.

---

## cat1_minus_2 (34ba198b, chrXII, − strand) — CORRECT, use as positive control

The GG was properly rescued by invoking a 1-bp del in the A4 tract. Module
2C (indel_correction) fired and the panel shows the corrected position
landing at a non-A on the RNA strand. Comparing this read against the three
buggy ones may help isolate the walkback divergence.

---

## cat1_plus_1 (0cb5a111, chrXIV, + strand) — over-extends past poly-A tail

**Genome at the corrected position:**
```
chrXIV pos:  10600       10610      10617    10620
plus seq:    AAAAAAAAAAA  T  A  G  C  T  C  T  A  T
                          11    13    14    15    16    17    18
```

Long genomic A-tract chrXIV[10600..10610], then `TAGCTCT...`.

- 4 of 5 aligners stop at corr_3p = 10610 (inside the A-tract).
- **mapPacBio extends to 10617** (`C` position), and HP-ED rewards it for
  +12 aligned bases despite forcing 4 mismatches at chrXIV[10611..10614] =
  `TAGC`: the read bases there are all A (poly-A tail), genome is non-A.
- HP-ED awards mapPacBio the win because the H/S penalty doesn't catch
  mismatches inside `=`/`X` ops as expensively as soft-clip extension would.

**User's diagnosis** (verbatim):
> "We are matching terminal non-A residues at the expense of too many
> mismatches. Over/under-calls of A-tracts are allowed, with perhaps a
> mismatch or two, but in this case we need to compare what's more likely:
> (a) two A → T nanopore sequencing errors that would mean the pA tail is
> simply longer, or (b) four consecutive errors GCTC → ATAA."

(a) is clearly more parsimonious than (b).

**Hypothesis**: the walkback / 3'-extension heuristic needs a guard on the
number of consecutive mismatches accepted when extending toward a "non-A
anchor". Threshold of ~2–3 consecutive A→non-A mismatches should switch
the interpretation from "real body extension" to "poly-A tail
continuation", and stop the alignment shorter.

**Secondary observation**: this is the case where the panel surfaced a
*winner-selection bug* downstream of the walkback bug. Even if walkback
is fixed for the losing aligners, mapPacBio's over-extension would still
be a real-data correctness issue. So this is two problems:
1. HP-ED scoring shouldn't reward mapPacBio's over-extension.
2. mapPacBio's extension itself shouldn't be produced.

---

## cat1_plus_2 (a146838d, chrI, + strand) — ambiguous; HP-aware ins should resolve

**Genome at the corrected position:**
```
chrI pos:    31535       31545       31555
plus seq:    GTCACCGAAAAGAAAAGGTAAAAAG
                    |    |
                  31543 31548
```

Two A-tracts: `GAAAA` (31541..31544), one `G`, `AAAA` (31546..31549), then `GGT`, then `AAAAA` (31553..31557).

- corr_3p = 31545 across all 5 aligners. chrI[31545] = `G` on plus = `C` on
  RNA. That's a non-A — **policy is satisfied**.
- But the user notes this is genuinely ambiguous between two paths:

**User's diagnosis** (verbatim):
> "If we let the A be a relatively free overcall of the AAAA inside GAAAAG
> to GAAAA(A)G, then we only need to invoke a single mismatch of G → A
> plus an undercall of the A-tract for the genomic match, vs GAG → AAA for
> the A-tail being longer. So... I actually think plus2 is a crabshoot."

User's preferred resolution:
> "Regardless, it would be optimal to align the first A mismatch as an
> insertion GAAAA(A)G, as the hp-aware scoring should favor that anyways."

**Hypothesis**: the HP-aware insertion cost `ins_cost(hp_len=4, base='A')`
should be very cheap, so a CIGAR like `…GAAAA 1I G…` should outscore the
current solution which probably has a terminal mismatch + flat penalty.
If the indel corrector isn't producing this form here, the HP-aware ins
lookup may not be firing on terminal positions, or the indel corrector
isn't being invoked on this read's terminal alignment.

This may already be covered by the debugger's Phase D ("Per-aligner
investigation of consensus drift").

---

## Suggested debugger workflow

1. Use **cat1_minus_2** as the positive control — its 3' end is correctly
   landed at a non-A on RNA. Diff the walkback execution path against
   cat1_minus_1 to isolate the minus-strand termination bug.
2. **cat1_minus_1** is the cleanest "obvious bug" — terminal RNA base IS
   an A, should be one position further. Likely a single-bp off-by-one in
   the walkback termination guard for minus-strand.
3. **cat1_plus_1** is the harder one — needs the mismatch-cap heuristic
   for the "extend-to-non-A" path, plus a re-evaluation of whether HP-ED
   is the right metric when the winning alignment forces 4+ mismatches.
4. **cat1_plus_2** ties into Phase D's consensus-drift investigation.
   Verify the HP-aware ins lookup fires on terminal A-tract positions.

The panel data lives at:
`rectify/data/validation/rectified/per_aligner_summary.tsv` (180 rows × 12
cols, all 36 validation reads × 5 aligners). Filterable by read_id.
