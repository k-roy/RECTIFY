# FEATURE REQUEST — teach the junction machinery minor-intron (U12) grammar

**Filed 2026-08-29 by the Chanfreau m6A/AT–AC session (Kevin Roy). Priority: HIGH — blocks
all long-read minor-intron work, and the fix is small. Full audit with evidence:
`~/work/UCLA/Chanfreau_Lab/planning/825j_panelB_rectify_scan.md` (read it before patching —
it carries the measured failure, the patch spec, and the verification design).**

## Problem (verified in source, four places)

RECTIFY as shipped cannot represent or fairly score AT–AC (minor/U12-type) introns:

1. `rectify/core/splice/splice_site_index.py:87` `SpliceSiteIndex.build` — **no kind key
   resolves to a plus-strand AC acceptor.** The AC and GT position arrays ARE already
   computed (lines 93–94) but bound to `don_minus`/`don_gt_plus` — wrong semantics for U12.
2. `overhang_informativeness.py:437` `_CANONICAL_DINUCS = {(GT,AG),(GC,AG),(CT,AC),(CT,GC)}`
   — a plus-strand U12 intron `(AT…AC)` and a minus-strand one `(GT…AT)` are **both absent**.
3. `junction_scoring.py:396/415` — plus-strand AT–AC lands at **tier 8 (worst)**; minus-strand
   at tier 3 (strand-asymmetric; arithmetic in the audit §0.2).
4. `junction_scoring.py:368` `_CANONICAL_HP_PRIOR = 0.5` — a true U12 junction at tier 8 must
   beat a canonical/phantom AG competitor at tier ≤4 by >0.5 HP-edit units to win.

## Measured consequence (real data, on disk)

At STX10 (twin AT–AC acceptors 11 nt apart, HeLa long reads in
`/scratch/users/kevinroy/rectify_human_validation/human_polyAminus_344/`): minimap2 places
**92.5% of reads (247/267 TLDR; 9/20 TERA) on a phantom AG acceptor at chr19:13149599 that
is annotated in NO GENCODE v46 transcript**, 5–6 nt from the true acceptors; neither real AC
acceptor receives more than 1 read. Chimera artifacts excluded (SA:Z 13.0% at phantom vs
12.2% window-wide). Any vanilla long-read acceptor quantification at such loci is ~100% wrong
— and rescue is the only route: the strand-aware `clip == 0` gate keeps 0/151 TLDR and 0/25
TERA reads here.

## Patch spec (do AT–AC and AT–AG TOGETHER — patching one alone creates the mirror bias)

1. `splice_site_index.py:87` `build()` — add plus-strand AC-acceptor and minus-strand
   AT-donor kinds using the arrays it already computes.
2. `splice_site_index.py:62` — **bump `_FORMAT_VERSION` 2 → 3** (non-optional:
   `load_or_build` will otherwise serve stale indexes).
3. `overhang_informativeness.py:437` — add `('AT','AC')` and `('GT','AT')` to
   `_CANONICAL_DINUCS`.
4. `junction_scoring.py:396/415` — give the minor classes their own tiers (not the NAT
   fallback), strand-symmetric; the `_CANONICAL_HP_PRIOR` interaction is the point of the
   patch (audit §"the 0.5-HP discount").

## Verification (design in the audit; zero new data)

Before/after recovery measurement on the 271 GENCODE v46 AT–AC introns against the
`human_polyAminus_344` BAMs (TERA DRS = ratios; TLDR = discovery-only, no UMI), stratified
on presence of a competing canonical AG within the slide window; `--junction-window 15`
covers the observed 5–6 nt displacement. Positive control: VEZT (+, recovered pre-patch);
failure exemplar: STX10 (−). Publishable as a methods result.

## Constraints for the implementing session

- Repo rules apply (`rectify/CLAUDE.md`): work on a branch off `master`, surgical staging,
  M1 is the source of truth; add tests for the new kinds/tiers; check `AGENT_FIXES.md` first.
- Do not regress GT–AG scoring: the recovery measurement above doubles as the regression
  gate (GT–AG placements must be unchanged).
