# Sumner SMA leads — SWEEP B: differential splicing (native re-placer discovery, recurrence-filtered)
25% genome-wide run (13/15 samples; 6 SMA + 7 WT — NOTE 2 largest SMA missing → SMA under-represented).
SMA-specific recurrent (>=3 SMA samples, 0 WT) revealed non-canonical junctions, annotated to genes.
★ DUAL PURPOSE: this IS the rigor test of the native re-placer on real disease data.

## Biology-anchored standouts (hard to dismiss as fabrication)
| gene | junction | in N SMA | why compelling |
|---|---|---|---|
| UBA1 | chrX:47211252-47212397 | 5 | ★★ mutations cause X-linked infantile SMA (SMAX2); ubiquitin-activating enzyme (proteostasis) |
| SNRPN | chr15:24977922-24978198 | 5 | ★★ spliceosomal snRNP (SmN) — snRNP assembly is SMN's CORE function; directly on-pathway |
| CACNA2D3 | chr3:54837229-54838566 | 5 | neuronal Ca channel subunit |
| PCBP2 | chr12:53468808-53471613 | 5 | poly(C)-binding, splicing regulator |
| L1CAM | chrX:153875935-... | 4 | neuronal adhesion / neurodevelopment |
| CLU | chr8:27598622-... | 4 | clusterin, neurodegeneration |
| HSPA5/BiP | chr9:125240306-... | 4 | ER stress / UPR (relevant to MN degeneration) |

## RED FLAG — expression-correlated fabrication pattern (the rigor concern)
Ribosomal-protein / high-expression genes dominate the rest: RPL7A, RPL12, RPL8, EEF1D, PABPC1, ADRM1, PSMC6.
These are the highest-expressed genes → most reads → most fabrication opportunity for the motif-blind re-placer.
Likely artifact, NOT biology. This is exactly why the SPIKE-IN precision number is needed before trusting discovery.

## RIGOR-TEST VERDICT (the point of the SMA panel)
The re-placer surfaces BOTH real SMA-pathway junctions (UBA1, SNRPN) AND likely expression-correlated fabrication
(ribosomal cluster). Encouraging that it FINDS real biology; concerning that fabrication co-occurs. => the native
aligner is NOT yet trustworthy for unguarded genome-wide discovery on real data without the precision calibration
(spike-in track) + short-read corroboration. The SMA panel did its job: it exposed the re-placer's real-data
precision ceiling.

## CAVEATS: fabrication-dominated raw; 2 SMA samples missing (SMA under-powered); highly-expressed-gene bias;
SMN positive control junctions appear WT-only single-sample (consistent with SMN1-null but the re-placer moves
some SMN paralog junctions — paralog ambiguity). NOT significance-tested; recurrence+specificity+annotation only.
VALIDATION PATH: spike-in precision → short-read (Illumina) corroboration of UBA1/SNRPN junctions → RT-PCR.

---
## ★★★ SHORT-READ CORROBORATION (Snaptron/recount3, ~316K samples) — the precision verdict (2026-07-09)
Method VALIDATED: SMN1 exon7-8 control = EXACT match, samples_count=216,950 (coord conversion js+1,je correct).
Per-lead exact/near junction lookup (samples_count = # short-read samples with that junction):
| lead | nearest short-read junction | verdict |
| SMN1 ctrl | EXACT, 216950 | method OK |
| CACNA2D3 | off-by-2 (1bp/side), 46386 | REAL junction (coord convention) — corroborated |
| SNRPN | off-by-12 (6bp shift), 205554 | REAL locus, re-placer DRIFTED junction ~6bp -> false non-canonical |
| UBA1 | off-by-26, 1-63 samples | NO support -> fabrication / gross mis-placement |
| PCBP2 | off-by-28, 1-2 samples | NOT corroborated -> fabrication |

## ★ THE RIGOR-TEST VERDICT (native re-placer, real data, orthogonal short-read truth)
The re-placer's genome-wide NON-CANONICAL discovery is dominated by DRIFT/MIS-PLACEMENT of REAL junctions, not
genuine novel biology: it shoves real (often highly-expressed) junctions a few-to-many bp off their true short-read-
supported coordinate, manufacturing spurious "non-canonical" calls (SNRPN: real 205K-sample junction drifted 6bp;
UBA1/PCBP2: no real junction within ~26bp -> pure fabrication). CRITICAL: this drift is at NON-HOMOPOLYMER positions
(6-28bp) -> the HP-drift guard does NOT catch it. => the native aligner needs a GENERAL drift/precision control
(short-read-anchored or expression-aware), not just the HP-specific guard, before genome-wide novel discovery is
trustworthy. The Snaptron short-read route is a POWERFUL, immediate, zero-download validation (per-lead definitive).
NEXT native-aligner build: use Snaptron short-read support as the precision oracle to characterize + fix the
non-HP drift (the sweep-B leads become the test set). COMPASS-mode (GSE108094 iPSC-MN FASTQ) = deeper confirmation.
