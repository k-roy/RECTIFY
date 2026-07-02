# RECTIFY Validation Read Case Studies

Mechanistic deep-dives of individual validation reads — *why* a particular
aligner won, what RECTIFY's correction did (or correctly declined to do), and
the genomic context that makes the call defensible. These are the "show your
work" companion to the bare per-read metadata in
[`VALIDATION_READS.md`](VALIDATION_READS.md) and the figure renderer in
[`../../scripts/validation_data/PLOTTING.md`](../../scripts/validation_data/PLOTTING.md).

## Purpose & audience

Two audiences, both of which keep re-deriving the same facts from scratch:

1. **Agents reviewing / extending the validation panels.** When you add reads or
   re-render figures, this doc tells you what "correct" looks like for each
   correction scenario and how to confirm it from the BAMs — not just from the
   pretty picture.
2. **Agents building the RECTIFY de-novo aligner.** These cases are a catalogue
   of *real ONT-DRS failure modes* (homopolymer undercalls, force-align
   over-extension, terminal soft-clip burial, cross-chromosome mis-placement)
   and of what a good aligner must do to win them. Treat them as acceptance
   criteria / regression targets for the new aligner.

Each case states only what the BAMs and reference actually show; reviewer
interpretation is flagged as such.

## How to reproduce the analysis (inspection toolkit)

All commands run from `rectify/data/validation/`. Genome:
`../genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa`
(contigs are `chrI`…`chrXVI`, `chrMito`). Use `conda run -n base` for samtools.

```bash
RID=<read-uuid>

# 1. Per-aligner RAW placement + terminal CIGAR (the 5 aligners RECTIFY compares)
for al in minimap2 gapmm2 mapPacBio deSALT uLTRA; do
  echo -n "$al: "
  samtools view aligners/validation_reads.$al.bam 2>/dev/null \
    | awk -v r="$RID" '$1~r{print $3":"$4" flag="$2" cig="$6}'
done

# 2. Per-aligner CORRECTED placement (after walkback / pA restore)
for f in rectified/per_aligner/*.bam; do
  out=$(samtools view "$f" 2>/dev/null | awk -v r="$RID" '$1~r{print $4" "$6}')
  [ -n "$out" ] && echo "$(basename $f): $out"
done

# 3. The final winner row (one row per read)
#    cols: 4 original_3prime, 5 corrected_3prime, 22 correction_applied,
#          26 gene_id, 39 winning_aligner  (see header: head -1 corrected_reads.tsv)
awk -F'\t' -v r="$RID" '$1~r{print $4,$5,$22,$26,$39}' corrected_reads.tsv

# 4. Reference context. For a MINUS-strand read, add -i to get the RNA (5'->3')
#    sequence so poly-A reads as poly-A, not poly-T.
samtools faidx <genome> chr<...>:<start>-<end>        # plus / genomic
samtools faidx -i <genome> chr<...>:<start>-<end>     # minus / RNA orientation
```

## Coordinate conventions (recap — see root `CLAUDE.md`)

All coordinates **0-based, half-open**. The **3′ end (CPA)** is:

| strand | 3′ end (CPA) | poly-A tail lives at | figure side |
|---|---|---|---|
| `+` | `reference_end - 1` (high coord) | BAM 3′ soft-clip | right |
| `-` | `reference_start` (low coord) | BAM 5′ soft-clip (read is RC) | right (axis inverted) |

`samtools view` `POS` is **1-based** → subtract 1 to compare with the
0-based `corrected_3prime`. For a minus-strand read, `reference_start` =
`POS-1` = the 3′ end.

## Cross-cutting aligner-behavior principles

Distilled from the cases below — the generalizable lessons for both reviewers
and the de-novo aligner:

1. **A non-A genomic landmark at/near the CPA is the telltale of true genomic
   origin.** When the true 3′ end abuts a non-homopolymer base (e.g. a `GG`
   dinucleotide), an aligner that *threads* it has proven the read is templated
   from the locus rather than floating on its poly-A tail. Aligners that
   soft-clip the landmark have thrown away the evidence and will mis-place the
   3′ end.

2. **ONT under-calls homopolymer length by ~1; the right model is a single `1D`,
   not a soft-clip and not a `1X`.** A winning aligner represents "one missing A"
   as a clean one-base deletion *inside* the alignment, keeping the anchor on
   the genomic landmark. Modeling it as a terminal mismatch (`1X`) lands one
   base short; soft-clipping it discards the anchor entirely.

3. **Force-alignment (no terminal soft-clip) is double-edged.** `mapPacBio`
   force-aligns the read ends. That *wins* when the terminus is real but
   slightly noisy (it keeps the anchor instead of clipping — see `cat1_plus_1`),
   and *loses* when it drags the end into non-homologous sequence as a wall of
   mismatches (see `cat1_minus_1`, `9X`). The de-novo aligner needs the upside
   without the runaway downside.

4. **RECTIFY's walkback is not always a rescue — it faithfully reports a bad
   alignment as bad.** Walkback trims terminal poly-A that an aligner
   force-aligned. If the aligner had *also* soft-clipped the true terminus
   (the landmark), walkback trims the few A's it did align and lands far short —
   correctly producing a *worse* EER-ED so that aligner loses. The win comes
   from a *different* aligner whose native alignment was right (see
   `cat1_minus_2`: walkback "fired" on minimap2/gapmm2/deSALT and they still
   lost to uLTRA's untouched alignment).

5. **The 3′ end can be perfect and the read still loses on the 5′ end.** A
   leading soft-clip (`5S`, `19S`) buries an otherwise-correct 3′ alignment in
   the EER-ED comparison. Winner selection is whole-read, not 3′-only.

6. **Cross-chromosome placement happens and must be surfaced, not hidden.** A
   single aligner can map the read to an entirely different locus
   (`cat1_plus_1`: deSALT → chrVI). It is excluded from the per-base 3′ view and
   shown as its own mini-panel (now rendered above the overview).

7. **Zero-shift is a first-class outcome.** "No correction" and "module fired
   but net-zero shift" are *successes* — RECTIFY must not over-correct an
   already-clean terminus (`cat1_plus_2`, `cat1_minus_2`).

## Case studies — `cat1_indel` (A-tract indel scenario)

Summary of all four reads (0-based 3′; `shift = corr − orig`):

| read | chrom | strand | orig 3′ | corr 3′ | shift | correction_applied | winner | gene |
|---|---|---|---|---|---|---|---|---|
| cat1_plus_2  | chrI   | + | 31551 | 31551 | 0  | `none` | mapPacBio | YAL062W |
| cat1_minus_2 | chrXII | − | 15345 | 15345 | 0  | `indel_correction` | uLTRA | YLL063C |
| cat1_plus_1  | chrXIV | + | 10617 | 10611 | −6 | `polya_walkback` | mapPacBio | — |
| cat1_minus_1 | chrII  | − | 9831  | 9834  | +3 | `atract_ambiguity,polya_walkback` | deSALT | YBL107C |

### cat1_plus_2 — zero-shift, nothing fired (clean negative control)

`a146838d` · chrI · **+** · span [31118, 31552) · orig = corr = **31551** ·
`correction_applied=none` · winner mapPacBio.

The terminal alignment was already at the leftmost-possible-CPA (first non-A
read=ref match), so **no module fired**. This is the control that proves RECTIFY
does not invent a correction when none is warranted. Verify: the base at 31551
is a non-A with an A-tract one base 3′ (right).

### cat1_minus_2 — module fired, net-zero, won by a *different* aligner's native CIGAR

`34ba198b` · chrXII · **−** · orig = corr = **15345** ·
`correction_applied=indel_correction` · **winner uLTRA**.

Reference (RNA 5′→3′, high→low coord; 0-based):

```
 ...15348 15347   15346 15345   15344 15343
      A     A   |   G     G   |   A     A
       A-tract     G G (landmark)  poly-A tail →
                       ▲ CPA = 15345 (the 3'-most G; a non-A)
```

Per-aligner, minus-strand 3′ end = `reference_start` (0-based):

| aligner | raw 3′ | corrected 3′ | walked back? | terminal CIGAR |
|---|---|---|---|---|
| **uLTRA** | 15345 | **15345** | no (byte-identical) | `2=1D17=…` |
| mapPacBio | 15346 | 15346 | no | `1=1X17=…` |
| minimap2 | 15348 | 15351 | **yes** (−3) | `2S17M…` → `5H14M…` |
| gapmm2 | 15348 | 15351 | yes | `2S…` → `5H…` |
| deSALT | 15348 | 15351 | yes | `2S…` → `5H…` |

- **uLTRA wins on its native alignment.** Its `2=` covers the `GG` landmark
  (anchoring the true 3′ end), it models the single missing homopolymer base as
  a clean `1D`, then matches into the body. Its raw and corrected CIGARs are
  identical — RECTIFY did **not** edit uLTRA. `correction_applied=indel_correction`
  is the read-level winning-path label, not an edit uLTRA received.
- **minimap2 / gapmm2 / deSALT soft-clipped the `GG`** (`2S`), force-aligned a
  few downstream A's, and walkback trimmed those A's → landed at 15351, ~6 bp
  short. Walkback "fired" on them and they **still lost** (principle #4).
- **mapPacBio modeled the missing base as `1X`** (mismatch, no indel) → landed
  at 15346, one base short. No terminal indel ⇒ indel_correction has nothing to
  act on ⇒ mapPacBio left untouched (principle #2).

### cat1_plus_1 — −6 bp poly-A walkback; force-align wins, 5′ clip buries the rest; one aligner elopes

`0cb5a111` · chrXIV · **+** · orig **10617** → corr **10611** (**−6 bp**) ·
`correction_applied=polya_walkback` · **winner mapPacBio** · 5′ pos 10429.

Reference (plus strand, 1-based) `chrXIV:10605-10625 = AAAAAAATAGCTCTATTCCGT`
— a 7-A run (10605-10611) ending at the corrected CPA `10611` (0-based) = `T`
(non-A), with the aligner's naive end 6 bp further into the tail at 10617.

Per-aligner raw placement:

| aligner | placement | terminal/leading CIGAR | note |
|---|---|---|---|
| **mapPacBio** | chrXIV:10430 | force-aligned, **no leading clip** | winner |
| minimap2 | chrXIV:10436 | `5S29M…` | 5′ soft-clip |
| gapmm2 | chrXIV:10436 | `5S29M…` | 5′ soft-clip |
| uLTRA | chrXIV:10436 | `5S29=…` | 5′ soft-clip |
| deSALT | **chrVI:8527** | — | **cross-chrom (wrong locus)** |

- **mapPacBio is the only aligner that both anchors the 5′ end** (force-aligns
  from 10429, no `5S`) **and handles the 3′ semantically** — placing the
  terminal A/T-rich bases against downstream genomic sequence rather than
  cramming a terminal insertion (reviewer interpretation of the CIGAR). RECTIFY's
  `polya_walkback` then recognizes the trailing T's as miscalled poly-A and pulls
  the end 10617 → 10611.
- **minimap2 / gapmm2 / uLTRA all carry a `5S` leading clip** that buries them in
  the whole-read EER-ED (principle #5), despite landing on the right locus.
- **deSALT mapped to chrVI:8527 entirely** — a cross-chrom digression, now shown
  in a mini-panel above the overview.

### cat1_minus_1 — A-tract ambiguity + walkback; force-align *fails* here

`77b392d9` · chrII · **−** · orig **9831** → corr **9834** (**+3 bp**; minus
strand, so the 3′ end / `reference_start` moved up 3, trimming 3 bp of tail) ·
`correction_applied=atract_ambiguity,polya_walkback` · **winner deSALT** ·
gene YBL107C.

Per-aligner raw placement:

| aligner | placement | leading/terminal CIGAR | note |
|---|---|---|---|
| minimap2 | chrII:9832 | `19S11M3I2M1I2M2D…` | 19-bp 5′ clip |
| gapmm2 | chrII:9832 | `19S11M3I2M1I2M2D…` | 19-bp 5′ clip |
| **deSALT** | chrII:9832 | `19S11M3I2M1I2M2D…` | winner |
| uLTRA | chrII:9832 | `19S11=3I2=1I2=2D…` | 19-bp 5′ clip |
| mapPacBio | **chrII:9810** | `2=2D4=9X2=1X1=1D…` | **force-align ran 22 bp out into `9X` garbage** |

- Four aligners agree at chrII:9832 (1-based) with a shared 19-bp leading clip;
  walkback + a-tract-ambiguity resolution trims 3 bp and **deSALT** wins the
  EER-ED tiebreak among the equivalent group.
- **mapPacBio's force-alignment is the cautionary opposite of `cat1_plus_1`**:
  here it dragged the terminus 22 bp out to 9810 as a wall of mismatches (`9X`),
  i.e. force-aligning into non-homologous sequence. Same mechanism, opposite
  outcome — the key tension the de-novo aligner must navigate (principle #3).

## Correction-module glossary (as observed in cat1)

| label (`correction_applied`) | what it does |
|---|---|
| `none` | terminus already at leftmost-possible-CPA; verified, not moved |
| `indel_correction` | resolve terminal A-tract indel artifact; anchor at first non-A read=ref match (read-level label; the winner may be an untouched native alignment) |
| `polya_walkback` | walk back through poly-A stop-base matches to the first non-stop read=ref agreement (= leftmost-possible-CPA) |
| `atract_ambiguity` | resolve which of several equivalent A-tract endpoints is the CPA when the homopolymer makes the boundary ambiguous |

## Case studies — `cat2_softclip` (3′ soft-clip rescue scenario)

### cat2_plus_1 — genomic homopolymer undercall the walkback used to eat (FIXED 2026-06-29)

`61b0c014` · chrI · **+** · winner **deSALT** @ corr 3′ **23754** (`correction_applied=none`).

Reference (plus strand, 1-based): a **24-bp genomic poly-A homopolymer** `chrI:23713-23736`, then an
AT-rich 3′UTR `TTAAATAAATAAATAAAATAAAT` (23737-23759), then ordinary sequence `ACAATGAT…` (23760+).

Three aligners' RAW alignments independently model the DRS undercall of the 24-A run as a **`9D` + `39=/39M`**
(matching through the AT-rich UTR): `…2I49(M|=)9D39(M|=)8S`. deSALT/mapPacBio instead use an insertion-heavy
model with no `9D`. The `9D` is real: uLTRA emits `39=` (exact matches) over a region containing **T's** — a
pure poly-A tail cannot match genomic T's, so those are genuine templated residues.

**The bug (now fixed):** the walkback's large-deletion pre-scan
(`walkback.py`, `_MIN_GENOMIC_ANCHOR_3P` guard) skipped past the `9D` (treating it as inside a force-aligned
poly-A tail) and over-walked through the `39=`, clipping the genomic match and landing ~23711. So the `9D`
appeared **only** in the unrectified row; every rectified per-aligner row had it clipped. The fix: do NOT skip
a ≥`large_del_min_bp` deletion when ≥`_MIN_GENOMIC_ANCHOR_3P` (=5) non-stop-base read=ref matches sit 3′ of
it (a genomic-homopolymer undercall flanked by real sequence, not a poly-A tail). Post-fix, minimap2/gapmm2/
uLTRA keep `…49(M|=)9D38(=|M)9S`, anchored at 23758. **Same bug class found in the upf1d cat1_minus_1 read**
(37 bp perfect C/G-rich match over-walked past a `6D`; test expectation corrected 1162861→1162817). Net
RECTIFY-output blast radius across all 36 reads: **0 winner/3′ changes** — the fix corrects per-aligner
diagnostic alignments only.

**Open EER-ED hypothesis (deferred Option B — NOT yet implemented):** even with the `9D` preserved, deSALT
still WINS this read because the HP-aware edit distance charges uLTRA a flat **9.0** clip penalty for honestly
soft-clipping the 9-base poly-A tail (`del_cost(hp=24)=0.0087`/base makes the `9D` itself ~free at 0.078; the
clip is the whole gap: uLTRA ED 12.99 vs deSALT 9.90). So the scorer **rewards force-aligning the poly-A tail
over honestly clipping it** — a likely systematic flaw (see `_cigar_hp_edit_distance`, corrected_consensus.py).
Proposed fix: make the 3′-clip penalty **graded by the probability the clipped tail is poly-A-derived**
(strand-aware A/T-richness; e.g. `9S=AAATAAAAT`=7A/9 → discount toward 0; a low-A genomic clip keeps full
penalty, auto-preserving the `cat1_minus_2`-class protection). Must be validated on reads with **unambiguous
CPA** (this read's AT-rich tail makes its CPA a genuine toss-up: 23754 vs 23758) + a full before/after winner
audit + a regression gate that clip-to-win still loses. Two independent advisors + reviewer agreed: ship the
walkback fix now, pursue the graded penalty as a separate truth-anchored effort. cat2_plus_1→uLTRA should be
a *consequence* of that effort, not its calibration target.

**OPTION B SHELVED — the graded clip penalty re-breaks clip-to-win (2026-07-01).** Investigation arc:
- First audit (jobs 32183378/32193027, `graded_clip_audit.py`) looked clean: 3/1929 flips, 0 SUSPICIOUS.
- **But that audit had a STRAND BUG** (advisor-caught, then verified): it only graded the *trailing* soft-clip.
  On minus-strand reads the 3′ poly-A tail is the *leading* clip (BAM SEQ is reference-forward; a minus gene's
  CPA sits at reference_start; the tail — poly-T on the forward strand — extends to lower coord = leading).
  So the grader was a **no-op on minus strand** (and could spuriously discount a T-rich genomic 5′ fragment —
  2 of the 3 "clean" flips, 902999d7/b164de40, were exactly this artifact).
- Strand-fixed the estimator (grade the tail-side clip: leading for minus / trailing for plus, scanning from
  the RNA 3′ terminus). Re-ran on the 36-read bundle → it now flips **cat1_minus_2 (the load-bearing
  regression case): uLTRA @ CPA 15345 (CORRECT) → minimap2 @ 15351 (6 bp WRONG).**
- Root cause: minimap2's leading clip `CCTTT` = `AAAGG` on the RNA — it contains the **GG genomic landmark
  uLTRA correctly threads**, but reads T-rich on the forward strand. A pure A/T-richness estimator CANNOT
  distinguish a mis-clipped genomic homopolymer-adjacent landmark from a real poly-A tail, so it forgives the
  clip and lets the wrong aligner win at the wrong CPA — exactly the clip-to-win failure the flat penalty
  exists to prevent.
- **Verdict: the simple tail-richness graded penalty is NOT safe to wire.** The flat clip penalty stays. A
  future Option-B would need a fundamentally stronger tail-vs-genomic discriminator (e.g. cross-aligner 3′
  consensus, or requiring the forgiven clip to NOT abut a genomic non-stop landmark), validated on BOTH
  strands with cat1_minus_2 as a hard regression gate. cat2_plus_1 keeps deSALT as winner for now.
- The **walkback fix (fc44ee2)** — restoring the 9D to diagnostic rows — is independent and unaffected.

### cat2 CATEGORY DRIFT — 3 of 4 reads no longer exercise `softclip_rescue` (found 2026-07-01)

**Validation-bundle-integrity finding.** cat2 reads were deliberately curated to demonstrate `softclip_rescue`
(VALIDATION_READS.md: `replace_cat2_reads.py` "replaces 3 Cat2 reads with DRS minimap2 reads showing
softclip_rescue"; category def: "aligner ended inside a genomic homopolymer... left matching downstream bases
as a 3' soft-clip... corrected position shifts outward"). But current `corrected_reads.tsv` shows the module
firing on **only 1 of 4**:

| read | correction_applied | winner |
|---|---|---|
| cat2_plus_1 (61b0c014)  | `none` | deSALT |
| cat2_plus_2 (88953e9c)  | `none` | uLTRA |
| cat2_minus_1 (b313b50d) | `none` | deSALT |
| cat2_minus_2 (9dbd37bf) | `atract_ambiguity, indel_correction, **softclip_rescue**` | deSALT |

**Root cause (NOT winner-selection picking a lucky aligner — an earlier mis-diagnosis).** For cat2_minus_1,
ALL FIVE aligners thread the homopolymer undercall as an inline **`3D`** and reach the correct CPA natively —
e.g. minimap2/deSALT/uLTRA `1S 6M 3D 24M…`, mapPacBio `1X 6= 3D 24=…` (minus strand, 3' end = reference_start
~186). The only clip is a trivial `1S` (1 base), NOT the matching-downstream-bases soft-clip the rescue module
recovers. So there is nothing to rescue → `correction_applied=none` is *correct*, but the read does not test
the module it's filed under. The reads were presumably selected when the aligner (older minimap2 params/
version) *soft-clipped* the post-homopolymer bases; the 2026-05-16 re-trim + re-align (and aligner-version
drift) now emit the `3D` inline instead. The `3D`-vs-soft-clip choice is aligner/param/version-dependent, so
the curated examples decayed from "demonstrates soft-clip rescue" into "aligner already threaded it."

**Why the tests didn't catch it:** `TestCategory2SoftClipRescue` asserts only the final 3' *position*
(`test_3prime_exact_position`), not that `softclip_rescue` fired — and the docstrings were updated to expect
`correction_applied=none`. So it's a *known-but-unguarded* state, not a silent regression, but the category no
longer guards the module it's named for.

**Action (bundle refresh / cDNA-agent + validation re-curation):** re-curate cat2 so it genuinely exercises
`softclip_rescue` — pick reads where CURRENT aligners leave a rescue-requiring soft-clip (not a `3D`), or pin
the alignment condition (aligner+params) that produces the soft-clip. Add a positive assertion that
`softclip_rescue` appears in `correction_applied` for cat2 reads (mirroring cat3's
`test_correction_applied_includes_five_prime_rescued`). Until then, only cat2_minus_2 is a valid cat2 example.

## Extending this doc

Add a `## Case studies — cat<N>_<name>` section as each category is reviewed.
Keep the discipline: a one-line metadata header, a small per-aligner table built
from the toolkit commands above, BAM-verifiable facts first, reviewer
interpretation explicitly flagged. Promote any new generalizable lesson into the
**Cross-cutting principles** list so the de-novo aligner agents get it without
reading every case.

## Case studies — human non-canonical junction discovery (A549, on-cluster)

A *junction-discovery* counterpart to the yeast CPA cases above — the de-novo
aligner acceptance criteria for **5′/3′ splice-site placement**, not 3′-CPA.
These are NOT in the bundled validation panel (human A549, data on Sherlock);
the full multi-platform + WGS derivation lives in
[`dev/COMPASS_2corroborated_CROSSPLATFORM.md`](../../../dev/COMPASS_2corroborated_CROSSPLATFORM.md)
(verdict locked in `adjudication_111_v3.json`). Summarized here as acceptance targets.

**Context.** Of GMAP's 111 chr5 "gmap-only recurrent novel GT-AG" candidates, rigorous
adjudication (ambiguity-tolerant short-read COMPASS panel + cross-aligner long-read
placement + A549 WGS) found ~107 artifacts and **3 real intragenic novel RNA splice
junctions** (DNA-confirmed genomic — A549 is near-triploid but these loci are flat-copy,
not deletions/rearrangements/circRNA/fusions). The diagnostic wrinkle:

| junction (host gene, +) | dominant SUPPORTED motif | canonical motif nearby | support |
|---|---|---|---|
| chr5:179824400-179832205 (SQSTM1) | `GT..GA` (non-canon) | GT-AG **4 bp** away | 2959 short + deSALT/uLTRA/GMAP |
| chr5:177592500-177593474 (TMED9)  | `CG..CA` (non-canon) | GC-AG **1 bp** away | 323 short + 4 LR aligners |
| chr5:140564954-140565547 (SLC35A4)| `AG..CA` (non-canon) | GT-AG **3 bp** away | 168 short + 3 LR aligners |

Every aligner (incl. accurate short reads) places these at the NON-canonical coordinate,
1–4 bp off a canonical GT-AG/GC-AG, and the ambiguity window is 0 (so the reads support the
non-canonical placement as a *distinct* junction, not the canonical neighbor).

**Acceptance criteria for the de-novo aligner (the open test):** given the supporting
reads at each locus, does motif-aware empirical-penalty realignment (a) **snap** the
junction to the nearby canonical motif (→ the truth was canonical, aligners mis-placed by
breakpoint read-errors), or (b) **hold** the non-canonical placement (→ genuinely
non-canonical splicing)? Either answer is a pass *if defensible from the reads*; the
failure mode is silently inheriting whichever the seed-and-chain panel preferred. RT-PCR /
Sanger across each junction is the wet-lab arbiter.

---

_Provenance: yeast cat1 reads from `wt_by4742_rep1` ONT-DRS, S. cerevisiae S288C
(R64-5-1); per-aligner BAMs under `aligners/`, corrected per-aligner outputs
under `rectified/per_aligner/`, final calls in `corrected_reads.tsv`. Authored
during the cat1 read-by-read figure review (branch `drs-validation-rebuild`)._
