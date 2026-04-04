# rectify.core

Core correction pipeline modules.

---

## bam_processor

The main BAM processing orchestrator. Handles parallel and streaming correction modes.

::: rectify.core.bam_processor
    options:
      members:
        - process_bam_file_parallel
        - process_bam_streaming
        - correct_read_3prime
        - write_position_index
        - get_read_3prime_position
        - get_read_5prime_position

---

## indel_corrector

Walk-back algorithm for poly(A)/A-tract indel artifacts.

::: rectify.core.indel_corrector
    options:
      members:
        - find_polya_boundary
        - rescue_softclip_at_homopolymer
        - detect_indel_artifacts
        - correct_position_for_indels

---

## consensus

Multi-aligner consensus selection.

::: rectify.core.consensus
    options:
      members:
        - AlignmentInfo
        - ConsensusResult
        - score_alignment
        - select_best_alignment
        - run_consensus_selection

---

## multi_aligner

Aligner wrappers for minimap2, mapPacBio, gapmm2, uLTRA, and deSALT.

::: rectify.core.multi_aligner
    options:
      members:
        - run_minimap2
        - run_mapPacBio
        - run_gapmm2
        - run_ultra
        - run_desalt

---

## atract_detector

A-tract ambiguity detection and classification.

::: rectify.core.atract_detector
    options:
      members:
        - detect_atract_ambiguity
        - count_downstream_as
        - classify_ambiguity_level

---

## netseq_refiner

NET-seq-guided position refinement for A-tract ambiguous reads.

::: rectify.core.netseq_refiner
    options:
      members:
        - NetseqRefiner
        - NetseqLoader

---

## unified_record

Per-read data structure capturing all correction results.

::: rectify.core.unified_record
    options:
      members:
        - UnifiedReadRecord

---

## processing_stats

Correction statistics accumulation.

::: rectify.core.processing_stats
    options:
      members:
        - ProcessingStats

---

## splice_aware_5prime

5' splice junction rescue.

::: rectify.core.splice_aware_5prime
    options:
      members:
        - rescue_5prime_junction
        - find_upstream_exon_match

---

## polya_trimmer

Poly(A) tail detection and length measurement.

::: rectify.core.polya_trimmer
    options:
      members:
        - trim_polya_tail
        - measure_polya_length
        - PolyAModel

---

## ag_mispriming

AG-mispriming detection for oligo-dT libraries.

::: rectify.core.ag_mispriming
    options:
      members:
        - screen_ag_mispriming
        - compute_ag_richness
