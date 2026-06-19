# SG-NEx A549 orthogonal validation data — CONFIRMED (2026-06-18)

Verified by live anonymous S3 listing + GitHub repo + AWS Open Data registry + the SG-NEx paper
(Chen et al., *Nat. Methods* 22, 801–812, 2025). For validating the 111 GMAP-only recurrent novel
junctions (`dev/gmap_only_recurrent_novels_chr5.tsv`) and the 609 corroborated novels with ORTHOGONAL,
high-base-accuracy evidence (the "error-model realism" validity claim, distinct from the sim benchmark).

## Access
- Bucket `s3://sg-nex-data`, **public, `--no-sign-request`, region `ap-southeast-1`**. No GEO/SRA/ENA
  accession — AWS S3 is the canonical path. Pull ON SHERLOCK (not the M1).
- Reference build: **GRCh38, Ensembl-style** (`Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa`) under
  `data/annotations/genome_fasta/`. ⚠ **CHROM-NAME MISMATCH:** SG-NEx uses Ensembl names (`5`), our A549
  junctions are on GENCODE/UCSC names (`chr5`). Harmonize (`chr5`↔`5`) before intersecting, OR confirm
  the SG-NEx BAMs' actual `@SQ` names first (some SG-NEx BAMs may be UCSC-named — check `samtools view -H`).

## Illumina short-read RNA-seq for A549 — YES (the gold standard for junction validation)
3 paired-end 150bp replicates: **replicate1, replicate3, replicate5** (run1 each).
- Genome BAMs (use directly): `s3://sg-nex-data/data/sequencing_data_illumina/bam/genome/SGNex_A549_Illumina_replicate{1,3,5}_run1/` (~5 GB each, ~15 GB total).
- Raw FASTQ (only if STAR 2-pass): `…/sequencing_data_illumina/fastq/…_R1/_R2.fastq.gz` (~25.5 GB).

## PacBio Iso-Seq for A549 — YES, single SMRTcell (orthogonal isoform support, NOT the gold standard)
`s3://sg-nex-data/data/sequencing_data_pacbio/bam/genome/SGNex_A549_PacBio-SMRTcell_replicate7_run1/` (3.27 GB).
HiFi = high base accuracy but still long-read (independent error profile). Good for isoform-level / phasing
support, complementary to Illumina split-read evidence.

## Validating a NOVEL junction with short reads — YES
A 150bp PE read spanning an exon-exon boundary maps as a SPLIT read across the intron gap; annotation
membership is irrelevant to detection. Two paths:
- **regtools `junctions extract`** directly on the provided Illumina genome BAM → BED of junctions +
  read support. FASTEST (no re-alignment); the recommended first pass.
- **STAR 2-pass** on the FASTQ → `SJ.out.tab` (col-6 annotated/novel flag + unique-read counts). More
  sensitive for novel junctions; needed only if regtools support is borderline.
Require ≥N uniquely-mapped split reads with adequate anchor on both flanks across ≥2 of 3 replicates.
**Caveats:** STAR `--alignIntronMax` (~590kb default) can miss ultra-long introns (raise it); low-expression
junctions may lack short-read coverage (absence ≠ refutation); short reads validate SINGLE junctions, not
multi-junction isoform phasing (that's where PacBio helps).

## TURN-KEY NEXT ACTION (on Sherlock scratch)
```
aws s3 sync --no-sign-request --region ap-southeast-1 \
  s3://sg-nex-data/data/sequencing_data_illumina/bam/genome/SGNex_A549_Illumina_replicate1_run1/ <dest>/
# repeat replicate3, replicate5 (~15 GB total)
regtools junctions extract -s XS <bam> > rep1.junc.bed   # per replicate
# harmonize chrom names (chr5<->5) then intersect against the 111 (and 609, and a noncanonical FP sample)
```
This is the orthogonal validation that turns "GMAP-only recurrent candidate" → "validated novel junction"
(or refutes it). Also lets us measure GMAP's ~198k non-canonical junctions' true FP rate against short-read
truth — quantifying exactly what the scoring fences must suppress.

Sources: github.com/GoekeLab/sg-nex-data · registry.opendata.aws/sgnex · nature.com/articles/s41592-025-02623-4
