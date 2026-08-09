# ★★ DNA TEST RESULT (2026-06-28) — the "structural deletion" call is REFUTED; these are RNA splicing events

A549 WGS (ENCODE ENCSR521ELB, Illumina; downloaded via python urllib — Sherlock login curl is broken;
subsampled ~1.3x, `$W/wgs/a549_wgs_partial.bam`) DNA coverage at the 3 SUPPORTED_NONCANONICAL candidate
"deletions":

| locus | interior meandepth | interior %cov | flank / genome avg | call |
| --- | --- | --- | --- | --- |
| 7805bp (chr5:179824400-179832205) | 1.63 | 73% | 1.49 / 1.52 | **NOT deleted (present)** |
| 974bp (chr5:177592500-177593474) | 0.95 | 68% | 1.57 | **NOT deleted (present)** |
| J1 593bp (chr5:140564954-140565547) | 1.50 | 62% | 2.86 | **NOT deleted (present)** |
| random-genome control | 1.52 | 71% | — | (normal-copy reference) |

Every interior has NORMAL DNA coverage = its flanks = genome average (~1.3x → ~73% Poisson-expected %cov).
A homozygous deletion would read ~0% covered. **The genomic sequence is PRESENT at all three loci → NOT
structural deletions.** Sequence present in DNA + removed in RNA (the N-gap) = **RNA processing (splicing)**,
not a genomic deletion and not an alignment artifact. So the 3 SUPPORTED_NONCANONICAL junctions are **real
non-canonical RNA junctions** in A549, corroborated across short + long reads AND DNA-confirmed as genomic.
The aligners place them at the non-canonical breakpoint because that IS the RNA junction (not a coincidental
nearby GT-AG). **This REVERSES the "structural deletion" verdict in the table below and in adjudication_111_v2.json.**

MECHANISM follow-up (2026-06-28): **NOT circRNA** (`dev/circ_check.py`): forward N-op orientation (direct read
dump), rolling-circle=0, head-to-tail=0 across minimap2/deSALT/uLTRA/GMAP. **NOT a distant fusion**
(`dev/chimera_invest.py`): the 7805 locus's 6446 minimap2-chimeric reads map their SA segments back into
SQSTM1 itself (chr5:179.8-179.9Mb, 6579) — local multi-exon read splitting, no distant partner. **All 3 are
INTRAGENIC, sense strand:** 7805∈SQSTM1(+), 974∈TMED9(+), J1-593∈SLC35A4(+) → candidate alt-splicing/
exon-skipping isoforms (SQSTM1 hosts TWO: the 7805 + J2's alt-acceptor — recurrent novel splicing in an
autophagy/cancer gene). Each dominant (non-canonical) coord has a canonical GT-AG/GC-AG 1-4bp away (ambiguity
window 0 → a DISTINCT placement, so reads support the non-canonical junction, not the canonical neighbor).
OPEN: genuinely non-canonical splice vs a true canonical junction mis-placed a few bp by read errors near the
breakpoint — resolvable by motif-aware RE-ALIGNMENT (the C1 empirical-penalty realigner) or RT-PCR across each junction.
A549 PLOIDY/CN CONTEXT (2026-06-28): A549 is NEAR-TRIPLOID (modal 66 chromosomes, Cellosaurus; has i(5p) +
high genome-wide SV burden) — NOT diploid. From our own ~10x WGS, broad 1Mb CN: all 3 junction loci + the rest
of chr5 sit uniformly at ~9.5x = the chr5 triploid baseline (≈3 copies); contrast chr1@100Mb 6.4x (~2n) and
chr17@40Mb 13x (~4n), so the assay DOES see CN variation — chr5q31.1/q35.3 are just at flat baseline with NO
focal amp/del. (Triploid reframes the het-test: a 3→2 single-copy loss = 0.67x, not 0.5x; the 974 locus's
small-window "0.78" was flank noise — its 1Mb CN 9.63 = chr5 baseline → no partial loss either.) No documented
SV at these exact 5q loci in literature; our WGS shows them structurally clean.
HET-SV EXCLUDED + NO REARRANGEMENT (2026-06-28, deep WGS ~9-10x, `$W/wgs/a549_wgs_deep.bam`,
`dev/dna_split.py`): all 3 interiors at NORMAL copy (interior/flank meandepth 0.99 / 0.78 / 1.0; 100%
covered at ~9-10x — a het deletion would read ~0.5x≈5x). DNA reads cross ALL breakpoints CONTINUOUSLY
(0 soft-clip≥15bp, 0 split/SA, 0 SA-to-partner at every donor+acceptor) → NO balanced rearrangement either.
**Every genomic alternative is ruled out: not homozygous-del, not het-del, not balanced-rearrangement, not
circRNA, not fusion. The genomic DNA is intact, continuous, normal-copy. The RNA gap is purely RNA-level =
SPLICING.** FINAL: the 3 are real novel intragenic RNA splice junctions (SQSTM1/TMED9/SLC35A4); the ONLY
open question is genuinely-non-canonical vs canonical-placed-a-few-bp-off (motif-aware realign / RT-PCR).
CAVEATS: (1) [resolved above — het-SV excluded at ~10x]. (2) Three non-canonical junctions is surprising → MECHANISM is
OPEN: non-canonical cis-splicing vs **circRNA back-splicing** (produces exactly this signature: present in
DNA, gapped in RNA, non-canonical in linear analysis) vs trans-splicing. NEXT: deeper WGS (parallel-stream
download — single-stream ENCODE/S3 was ~570 KB/s) + a circRNA/back-splice analysis (out-of-order / head-to-tail
junction orientation) would resolve the mechanism. The ONT A549 WGS (ENA PRJEB90580) would confirm at the DNA SV level.

---

# ★ 111 RE-DERIVATION (2026-06-26 PM) — the locked "109 artifacts / 2 corroborated" is superseded

The locked `adjudication_111.json` used EXACT ambiguity-normalized coordinate matching, which (a) UNDERCOUNTS
corroboration (strict base-equality normalize() leaves same-length placements a few bp apart unmerged) and
(b) never checked the splice MOTIF at the supported coordinate. Re-derived with ambiguity-tolerant
(same-intron-length, |Δstart|≤5bp) matching + motif-check at the supported coord:
`dev/rederive_111.py` → `$W/rederive_111.json`.

**EXACT-supported = 2 (reproduces the locked verdict) → RELAXED-supported = 4.** Classes: **107 ARTIFACT**
(no independent short-read support even relaxed — these are the solid GMAP artifacts), **4 SUPPORTED**.
The 4 supported, with motif checked at the dominant short-read coordinate + an offset scan to the nearest
canonical site:

| junction | len | gmap | short-rd (relaxed) | dominant SR coord motif | nearest canonical | verdict |
| --- | --- | --- | --- | --- | --- | --- |
CROSS-ALIGNER long-read placement (the decisive, non-bespoke discriminator; `dev/lr_probe_4loci.py` →
`$W/lr_probe_4loci_out.txt`). Each cell = anchored long-read count placing the same-length gap near the locus:

| junction | len | short-rd | minimap2 | deSALT | uLTRA | GMAP | mapPacBio | nearest canon | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 179824400-179832205 | **7805** | **2959** | 0 | 137 | 107 | 149 | 0 | GT-AG 4bp away | **STRUCTURAL DELETION** (recurrent ~7.8kb del in A549; INSIDE SQSTM1) |
| 177592500-177593474 | 974 | 323 | **628** | **626** | **485** | 491 | 0 | GC-AG 1bp away (gets ~0) | **STRUCTURAL DELETION** (~974bp; all aligners place 1bp off the GC-AG, splice-aware ones do NOT snap to it) |
| 140564954-140565547 (J1) | 593 | 168 | 102 | 103 | 81 | 56 | 0 | GT-AG 3bp away (12 reads) | **STRUCTURAL/genomic** (non-canonical, multi-aligner) |
| 179823051-179823857 (J2) | 806 | 3 | 1 | 20 | 8 | 37 | 0 | at coord (GT-AG) | INCONCLUSIVE minor +2 alt-acceptor of annotated SQSTM1 intron (within ONT acceptor smear) |

**THE DECISIVE PATTERN (one lens, not four stories):** for 7805 / 974 / J1, MULTIPLE INDEPENDENT long-read
aligners (deSALT, uLTRA, minimap2) AND Illumina short reads all UNIFORMLY place the gap at the NON-canonical
coordinate, and splice-aware aligners do NOT snap to the canonical GT-AG/GC-AG sitting 1–4bp away (974's
canonical GC-AG gets ~0 reads vs 600+ at the non-canonical). That is the signature of **recurrent structural
DELETIONS in A549** (a CNV-rich cancer line), NOT novel splice junctions — the nearby canonical dinucleotides
are coincidental. J2 is a minor alt-acceptor within the long-read acceptor smear (3 short reads) → inconclusive.

**FINAL VERDICT ON THE 111:** ~107 are clean GMAP artifacts (no independent support — the locked json is right
for these). The 4 with independent short-read support are NOT validated novel splice junctions: **3 are
structural-deletion signatures** (real A549 DNA features, the standout being the recurrent ~7.8kb deletion
inside SQSTM1, 2959 short + ~400 long reads) and **1 is an inconclusive alt-acceptor**. **NONE of the 111 is
rescued as a clean canonical novel splice junction.** The locked `adjudication_111.json` "2 corroborated →
likely real novel junctions" is WRONG in interpretation (those 2 are an SV-signature + an inconclusive
alt-acceptor); and its "109 artifacts" both over-counts (107) and mis-frames the 4 supported (they ARE
independently corroborated — as structural/alt-SS features, not as splicing).

RESIDUAL (noted, not chased): the same-intron-LENGTH relaxed match cannot see a real alt-SS that CHANGES
intron length (e.g. J2's 804→806) — so "107 artifacts" may still slightly over-count. Same-length is the
correct conservative choice (|Δstart|≤5 alone false-matches different-acceptor alt-SS).
CLASSIFIER CAVEATS (auto-`class` field — fix before reuse): 974 mislabeled SUPPORTED_NONCANONICAL vs the
1bp-off canonical; J2 should be ALT_SS_OF_ANNOTATED (XOR-on-normalized-coord missed it). Both moot for the
verdict (cross-aligner lens overrides), but fix before the classifier is reused on other chromosomes/samples.

NEXT (confirmation): IGV the 7805bp + 974bp loci (deletion-breakpoint signatures: soft-clips/microhomology vs
clean exon/exon) and overlay an A549 genomic-DNA/WGS/CNV track to confirm deletion vs splicing. The locked
`adjudication_111.json` (chmod 444) now actively contradicts this re-derivation → ask Kevin whether to drop a
SUPERSEDED sidecar / re-lock a corrected json.

---

# COMPASS 2-corroborated junctions — cross-platform verdict (2026-06-26)

Follow-up to the COMPASS short-read 111-adjudication (`adjudication_111.json`,
109/111 = artifacts, **2 corroborated**). Characterizes the 2 corroborated
junctions with the specifics the adjudication booleans hide — splice motif +
strand, overlapping-gene strand, short-read splice-in, the locus splice
landscape, AND a **cross-platform (long-read vs short-read) usage-fraction**
check (per Kevin's point: a real low-rate alt-acceptor should show a *similar
usage frequency* across platforms; a seq/alignment error would not).

Data: short-read = locked Oak COMPASS consensus
`$OAK/A549_rep1.consensus.bam`; long-read = SG-NEx A549 ONT panel
`…/sgnex_a549/alignments/a549_chr5_trimmed.{minimap2,deSALT,mapPacBio,GMAP,uLTRA,rectified}.bam`.
All counts = anchored N-op junctions, min_anchor=10, via the live rectify
`collect_junction_counts_from_bam`. Scripts: `dev/compass_corroborated_junctions_characterize.py`,
`dev/longread_probe.py`. Self-check: recomputed short-read depths = 12 and 3 (match locked json). ✅

---

## J1 — chr5:140564954-140565547 (593 bp intron, SLC35A4 +) → **REAL novel junction (high confidence)**

Canonical **GT..AG**, plus strand, strand-AGREE with SLC35A4 (+). **No annotated
junction anywhere in the locus** (genuinely novel). Intron 593 bp.

The biological junction appears at two placements 3 bp apart (140564954-547 and
140564957-550, same 593 bp intron) that strict base-equality normalization left
UNMERGED — so the adjudication's "depth 12" UNDERCOUNTS it.

| Platform / aligner | reads at the 593 bp novel intron |
| --- | --- |
| COMPASS short-read consensus (Illumina) | **168** (12 @…954-547 + 156 @…957-550) |
| minimap2 (ONT) | 101 |
| deSALT (ONT) | 102 |
| uLTRA (ONT) | 81 |
| GMAP (ONT) | ~44 (8 + 34 + smear) |
| mapPacBio (ONT) | ~0 |
| long-read consensus | 105 |

**Conclusion:** corroborated by FOUR independent long-read aligners AND the
independent short-read panel → a real novel 593 bp intron in SLC35A4. **It is NOT
gmap-only** — the "111 gmap-only" classification for this junction was an artifact
of a 3 bp ambiguity-normalization gap (minimap2/deSALT/uLTRA place it 3 bp shifted
from the gmap coordinate, so the strict-normalized sets didn't intersect). IGV
should confirm a clean canonical novel junction; recommend ACCEPT.

---

## J2 — chr5:179823051-179823857 (806 bp, SQSTM1 +) → **INCONCLUSIVE (leans not-confirmed)**

Canonical **GT..AG**, plus strand, strand-AGREE with SQSTM1 (+). It is a **+2 bp
alternative acceptor** off the **annotated** SQSTM1 intron: same donor
(179823053), acceptor 179823859 vs the annotated 179823857 (804 bp, annotated).

Acceptor usage at the shared donor 179823053:

| Platform | acceptor 857 (annotated) | acceptor 859 (+2 novel) | acceptor smear 860-862 |
| --- | --- | --- | --- |
| COMPASS short reads (precise 3′) | **2258 (sharp)** | 3 (0.13%) | ~0 |
| GMAP (ONT) | 6 | 8 | 860×11, 861×3, 862×3 |
| deSALT (ONT) | 9 | — (049-859×1) | 860×7, 861×3 |
| uLTRA (ONT) | 2 | — | 860×5, 861×1 |
| long-read consensus | 5 | — | 860×6, 861×1 |

**Why inconclusive (Kevin's test, applied):** short reads have precise 3′ ends and
place the acceptor overwhelmingly at the annotated 857 (2258), with 859 at only
0.13% — i.e. splice-noise level. Long reads (ONT) **cannot resolve the exact
acceptor**: at this donor they smear across 855–862 with no clean 859 peak (859 is
no more populated than its neighbors 857/860/861/862, and GMAP is the only aligner
landing exactly on 859 at all). So the two platforms do **NOT** agree on a specific
~0.13% alt-acceptor frequency at 859 — the long-read "support" is 3′-end
imprecision around the dominant junction, not a frequency-matched corroboration.

So J2 is most consistent with low-rate short-read noise + ONT acceptor smear rather
than a robustly distinguishable novel alt-acceptor. It is **canonical and not
provably a seq error**, but the cross-platform frequency agreement that would argue
"real low-rate alt-acceptor" is absent. Recommend: treat as unconfirmed; if pursued,
the deciding evidence would be precise-3′ short-read reads landing *sharply* and
*reproducibly* at 859 above background (here only 3, vs 2258 at 857).

> Note: at this donor BOTH platforms are dominated by entirely different SQSTM1
> isoforms (the len-340 acceptor family at 179823721-824061 ×47-58, and the
> len-7805 junction at 179823506-831311 ×28-42); the 857/859 micro-distinction is a
> minor sub-population in long reads.

---

## Process note (worth a fix upstream)

The "gmap-only recurrent novels" set (`gmap_only_recurrent_novels_chr5.tsv`) is
built on STRICT base-equality normalization (`normalize_junction`), which leaves
same-intron placements a few bp apart UNMERGED when an intervening base breaks the
equality slide. J1 shows this directly: minimap2/deSALT/uLTRA support the SAME 593
bp intron 3 bp from the gmap coordinate, yet it was labeled "gmap-only." A
length-aware / windowed merge (group placements sharing intron length within the
ambiguity window) before the gmap-only test would prevent over-counting gmap-only
junctions. Recommend auditing the 111 for other 3 bp-shift co-supported cases.
