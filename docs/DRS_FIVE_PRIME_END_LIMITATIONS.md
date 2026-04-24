# 5' End Resolution in Native Poly(A) DRS: Limitations and Alternatives

> **Author:** Kevin R. Roy  
> **Context:** Assessment of whether 5' end clustering from standard poly(A) DRS data
> is scientifically defensible. Conclusion: it is not, without orthogonal validation or
> cap-selection library preparation. This document explains why and describes the
> correct approaches.

---

## Short Answer

Standard native poly(A) DRS performs poorly for resolving true 5' ends. This is not
a minor edge case — it is an intrinsic, multi-layered limitation of the sequencing
chemistry and library architecture. There are two distinct problems:

1. Most reads are truncated and never reach the 5' end at all.
2. Even reads that do reach the physical 5' end are missing the terminal ~11 nt.

Together, these make TSS calling from standard DRS unreliable without orthogonal
methods.

---

## 1. The Directionality Problem: DRS Reads 3'→5'

The poly(A) DRS workflow loads the RNA 3' end first. The poly(T)/RTA adapter
hybridizes to the poly(A) tail, the motor protein engages, and the RNA is threaded
through the pore **3'→5'**. This means the 5' end of the transcript is the **last**
thing to translocate — if it translocates at all.

Any interruption — pore stalling, RNA secondary structure, degradation, or motor
dissociation — terminates the read before reaching the 5' end. All such truncations
appear identical in the data to genuine 5' ends. There is no signal in the BAM
distinguishing a truncated read from a full-length read.

---

## 2. The Read Truncation Rate

Nanopore DRS reads are frequently truncated prior to annotated transcription start
sites, resulting in a strong **3' coverage bias**. This is not a small effect:

- Read coverage declines towards 5' ends genome-wide.
- Annotation tools frequently misclassify truncated reads as independent isoforms
  with alternative short 5' exons.
- Approximately 20% of reads are truncated during sequencing due to signal noise
  alone; real-world truncation rates from RNA degradation, secondary structure, and
  motor dissociation are considerably higher on long transcripts.
- With standard ONT nanopore libraries, only a fraction of reads start at the
  annotated 5' end; a large fraction start within annotated exons.

**The critical ambiguity:** standard DRS is unable to differentiate full-length
transcripts from RNA fragments resulting from degradation or incomplete sequencing.
Every internal 5' soft-clip looks the same regardless of whether it represents a
decay intermediate, a truncated-read artifact, or a real alternative TSS.

---

## 3. Even "Full-Length" Reads Miss ~11 nt of the True 5' End

Even when a read traverses the entire transcript body, there is a systematic,
chemistry-intrinsic loss of terminal 5' sequence. Studies using synthetic ERCC RNA
Spike-Ins and *in vitro* transcribed RNAs have shown that reads lack approximately
**11 nt** of authentic 5' sequence.

- The precise length of missing 5' sequence varies across individual molecules,
  suggesting sequence- or context-specific effects on accuracy.
- This offset reflects the fact that when the 5' terminal nucleotides exit the pore,
  insufficient downstream RNA remains to maintain motor engagement and generate
  interpretable signal — the last bases in the pore do not produce a reliable ionic
  current event before the molecule dissociates.
- Nanopore DRS can only sequence to within **5–10 nt of the 5' cap**.
- Distance-to-nearest-TSS analyses show a consistent downstream offset even for
  the best-aligned reads; distance-to-CAGE analyses confirm the bias is real and
  genome-wide, with a substantial fraction of called TSSs falling downstream of CAGE
  TSSs.

---

## 4. Downstream Consequences for TSS Annotation

The combination of:

- **(a)** pervasive truncation
- **(b)** inability to distinguish truncation from biological 5' ends
- **(c)** the intrinsic ~11 nt cap-proximal gap

means that 5' end calls from standard DRS are **systematically shifted 3' of the
true TSS**. Bioinformatic tools that attempt isoform reconstruction can partially
compensate using reference annotations, but annotation tools frequently produce
transcript models shorter than the aligned reads, with truncated isoforms arising
from DRS limitations, and a substantial fraction of predictions lack their 5' exons
entirely.

---

## 5. Why 3' Ends Are More Reliable (Contrast)

It is worth noting what DRS does well: **3' end and poly(A) tail resolution**.

- The poly(A) tail enters the pore first.
- Cleavage/polyadenylation site mapping from DRS is relatively accurate because the
  terminal 3' nucleotides are directly adjacent to the adapter junction, not subject
  to motor drop-off.
- Internal priming artifacts that plague dT-primed cDNA libraries are largely absent.

This asymmetry is why rectify focuses on 3' end correction as its primary output.
The `five_prime_position` column is retained in `corrected_3ends.tsv` and can be
used for rough 5' coverage analysis, but should **not** be interpreted as a reliable
TSS measurement.

---

## 6. Solutions: Cap-Selection Approaches

The community has converged on cap-based selection as the correct approach for TSS
mapping from long-read data. These methods restrict full-length calling to reads that
carry confirmed evidence of the m7G cap:

| Method | Strategy | Key advantage |
|--------|----------|---------------|
| **NRCeq** (Nanopore ReCappable-seq) | Enzymatic cap replacement with azido-modified cap; click chemistry to ligate RNA adapter | Directly sequenced on nanopore; avoids ligation bias at 5'-truncated ends |
| **CapTrap-seq** | Cap-trapping (bisulfite oxidation of non-capped 5' ends → degradation) + oligo-dT; produces cDNA for long-read seq | Platform-agnostic; compatible with ONT and PacBio |
| **TLDR-seq** | Cap-trapping + 3' adapter ligation; works on pA+ and pA− molecules | Accurately determines 5' and 3' ends without ribo-depletion |
| **Biotinylated 5' adapter ligation** (e.g., Parker lab method) | Ligation of a biotinylated RNA oligo to the cap; streptavidin pulldown | Simple enrichment for capped, full-length molecules |

NRCeq avoids the common shortfall of ligation methods such as sequence bias and
ligation at 5'-truncated ends, enabling sequencing of individual RNA molecules
containing the poly(A) tail through the 5' cap (defining full-length RNA isoform
scaffolds). Using NRCeq-derived high-confidence scaffolds, overlap with orthogonal
TSS markers was **92.2%** for high-confidence scaffolds, compared to only ~36–38%
for untreated nanopore reads.

TLDR-seq shows high overlap with CAGE TSSs and substantially improves TSS detection
compared with the standard ONT nanopore protocol.

---

## 7. Summary

| Capability | Standard poly(A) DRS | Cap-selected DRS |
|------------|---------------------|-----------------|
| True 5' end resolution | ❌ Poor (~11 nt short even for longest reads; most reads truncated well upstream) | ✅ Good (TSS overlap with CAGE ~92%) |
| Distinguishing full-length from degradation/truncation | ❌ Impossible | ✅ By design |
| 3' end / poly(A) site | ✅ Reliable | ✅ Reliable |
| Splicing & isoform structure | ✅ Good (body of transcript) | ✅ Good |
| Native modifications (m6A, etc.) | ✅ Yes | ✅ Yes (if using DRS-based method like NRCeq) |

---

## Implications for Rectify

Rectify's `analyze` subcommand includes 5' end clustering (`tss_clusters.tsv`) as
a convenience output. **This output should be interpreted with these caveats in mind:**

- The 75 bp clustering window (vs 25 bp for 3' ends) reflects the known noisiness
  of DRS 5' ends — it is not a sign that 5' clustering is working well, it is an
  acknowledgment that 5' calls are imprecise.
- `five_prime_position` in `corrected_3ends.tsv` reflects where the aligner found
  the 5' end of the read, which for most reads is NOT the true TSS.
- The 5' cluster counts are useful for **relative comparisons** between conditions
  (e.g., detecting shifts in read 5' coverage between WT and upf1Δ), but should not
  be used to make claims about the positions of TSSs without CAGE or NRCeq
  cross-validation.
- For the NMD analysis specifically: upf1Δ-specific 5' extensions in the TSV may
  reflect genuine 5' UTR extensions accumulating due to NMD protection, but the
  exact nucleotide positions cannot be trusted without cap-seq validation.

---

*Last updated: 2026-04-22*
