# Phase-0 cDNA-alignment stress loci (A549 direct-RNA, chr5-first)

Hand-built cDNA-library seed set for testing whether a **contiguous aligner vs
a pre-spliced cDNA reference** beats a splice-aware genome aligner on
micro-exon-crossing and dense-multi-junction reads.

- Coordinates: **GENCODE v46 / Ensembl GRCh38, MANE Select transcript**, Ensembl
  chr naming (`5`, not `chr5`).
- Exon lengths/counts computed from **raw GTF coordinates** (`end-start+1`),
  not from any web prose (Ensembl REST page summaries mis-stated DIAPH1 as
  "8 bp/26 bp" — the coordinates are 9 bp/27 bp; TRIO prose said 46 exons,
  raw MANE = 57). Trust the table below.
- Expression breadth from Human Protein Atlas cell-line RNA (A549 included in
  the n=232 lung-cancer panel). HPA reports gene-level detection, **not exon
  inclusion**. Microexon *inclusion* in A549 must be confirmed in the SG-NEx
  A549 BAM before building the library (flagged per locus).
- Render date: 2026-06-14. Source GTF:
  ftp.ebi.ac.uk/.../gencode.v46.basic.annotation.gtf.gz

## Ranked set

| rank | gene | Ensembl ID | chr | strand | gene body (GRCh38) | MANE tx | exons | stress feature (MANE) | A549 expr (HPA) | inclusion status |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | DIAPH1 | ENSG00000131504 | 5 | - | 141,515,021-141,619,000 | ENST00000389054 | 28 | **9 bp** internal exon @141,571,428-141,571,436 + 27 bp @141,588,224-141,588,250 | Detected in ALL (Tau 0.15) | likely constitutive (not a documented neural-skip); confirm in BAM |
| 2 | TRIO | ENSG00000038382 | 5 | + | 14,143,342-14,510,204 | ENST00000344204 | **57** | 28 bp @14,497,847-14,497,874 in a 57-exon gene; seed chain fragments | Detected in ALL (Tau 0.28) | gene-expression sufficient (junction-count case) |
| 3 | MYO10 | ENSG00000145555 | 5 | - | 16,661,907-16,936,288 | ENST00000513610 | **41** | **14 bp** internal exon @16,780,728-16,780,741 in a 41-exon gene | Detected in ALL (Tau 0.27) | 14 bp exon is coding/internal, NOT the known headless 5'UTR neural event; confirm in BAM |
| 4 | ABLIM3 | ENSG00000173210 | 5 | + | 149,141,493-149,260,542 | ENST00000309868 | 24 | **8 bp** @149,252,201-149,252,208 + two 30 bp (@149,239,248; @149,249,815); dense | Detected in MANY (kidney-enhanced, Tau 0.45) | moderate A549; confirm in BAM |
| 5 | KIF3A | ENSG00000131437 | 5 | - | 132,692,631-132,737,546 | ENST00000403231 | 19 | **9 bp** internal exon @132,706,451-132,706,459 | broadly expressed kinesin | confirm in BAM |
| 6 | MPC1 | ENSG00000060762 | **6** | - | 166,364,919-166,382,940 | ENST00000360961 | 5 | **4 bp constitutive** microexon (exon 2) @166,370,218-166,370,221 | ubiquitous metabolic gene | constitutive (documented) — best-anchored microexon |

### Off-chr5 budget
- MPC1 (chr6) is the **one off-chr5 pick** — taken because its 4 bp exon is a
  *documented constitutive* microexon (Lin et al. 2023, Mol Ther Nucleic
  Acids), so inclusion in A549 is guaranteed, unlike the alternatively-spliced
  microexons above. Highest-confidence true-microexon stress case in the set.
- Budget allows a 2nd off-chr5 locus if desired; not used.

### Bench / backup (not in core 6)
- **COL23A1** (ENSG00000050767, chr5 -, 178,237,618-178,590,393, MANE
  ENST00000390654, 29 exons): the structural champion — **9 internal exons
  <=30 bp** (one 8 bp @178,290,362 + eight 27 bp exons). Dropped from core only
  on A549-expression doubt (HPA: cancer-enhanced bone/leukemia, Tau 0.80,
  "detected in some"). **Promote to core if the A549 BAM shows coverage** — it
  is by far the strongest multi-microexon stressor on chr5.
- **MATR3** (ENSG00000280987/...394805): dropped — MANE transcript's shortest
  internal exon is 42 bp (the <=30 bp exons were non-canonical). Not a stress
  case in the isoform A549 expresses.
- **ANXA6** (chr5, 18 bp @151,110,627): solid A549 expression, single short
  exon; weaker than the above. Available as a backup.

## Why these stress a splice-aware genome aligner

minimap2 (`-x splice`, k=15) seeds anchors of length k and chains them across
introns. A **9-14 bp internal exon is shorter than the seed length**, so it
cannot nucleate an anchor; the aligner either (a) skips it — folding it into a
flanking intron and emitting a wrong junction, or (b) soft-clips/mismaps the
read. Dense multi-junction genes (TRIO 57 exons, MYO10 41, ABLIM3 24) fragment
the seed chain so co-linear chaining drops short exons even when individually
seedable. A **contiguous aligner against the pre-spliced cDNA** has no introns
to chain across: the read is one ungapped/near-ungapped match, so the
micro-exon bases are placed correctly by construction. That is the win this
Phase-0 test is built to demonstrate. (Confirmed framing: uLTRA / ExceS-A
papers note exons <30 nt are below seed length and are the dominant split-
alignment failure mode.)

## Microexon-benchmark literature note

The canonical neuronal-microexon benchmark set is **Irimia et al. 2014 Cell**
(SRRM4/nSR100-regulated, 3-27 nt). Those microexons are **SRRM4-dependent and
SRRM4 is low in A549** (and SRRM4 silencing across cancers suppresses microexon
inclusion — Gonatopoulos-Pournatzis/Torres-Mendez line of work), so most
classic neuronal microexons are **skipped** in A549 and would test nothing.
The 2 bp human microexons in the literature are **GRK6** (ENSG00000198055,
chr5:177,426,659-177,442,901) and **SEPTIN7** (chr7) — GRK6 is on chr5 and
expressed in lung, but its 2 bp exon is in a minor/alternatively-spliced
transcript (canonical GRK6 has no <=30 bp exon), so include only if the A549
BAM shows reads crossing it. Net: the neuronal-microexon set is mostly off-chr5
and inclusion-suppressed in A549 — hence this set leans on **constitutive
microexons (MPC1)** and **canonical short exons in broadly-expressed,
dense-junction genes (DIAPH1, TRIO, MYO10, ABLIM3, KIF3A)**.
