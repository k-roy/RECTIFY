# rectify.core

Core correction pipeline modules.

---

## bam_processor

Per-read correction orchestrator. Owns the `correct_read_3prime` workhorse and
the small helpers it needs (3'/5' position extraction).

::: rectify.core.bam.bam_processor
    options:
      members:
        - correct_read_3prime
        - get_read_3prime_position
        - get_read_5prime_position
        - process_bam_file

---

## parallel

Parallel and streaming dispatch over a BAM. Wraps the per-read
`correct_read_3prime` worker for multi-process and chunked execution.

::: rectify.core.bam.parallel
    options:
      members:
        - process_bam_file_parallel
        - process_bam_streaming
        - process_bam_streaming_parallel

---

## indel_corrector

Walk-back algorithm for poly(A)/A-tract indel artifacts.

::: rectify.core.correct.indel_corrector
    options:
      members:
        - find_polya_boundary
        - rescue_softclip_at_homopolymer
        - detect_indel_artifacts
        - correct_position_for_indels

---

## consensus

Multi-aligner consensus selection.

::: rectify.core.consensus.consensus
    options:
      members:
        - AlignmentInfo
        - ConsensusResult
        - score_alignment
        - select_best_alignment
        - run_consensus_selection

---

## multi_aligner

Aligner wrappers. The documented long-read panel is minimap2, mapPacBio,
gapmm2, uLTRA, and deSALT; additional opt-in wrappers (GMAP, Winnowmap2,
Minisplice, plus short-read bbmap/bwa) also live in this module.

::: rectify.core.align.multi_aligner
    options:
      members:
        - run_minimap2
        - run_map_pacbio
        - run_gapmm2
        - run_ultra
        - run_desalt

---

## atract_detector

A-tract ambiguity detection and classification.

::: rectify.core.polya.atract_detector
    options:
      members:
        - calculate_atract_ambiguity
        - get_ambiguity_category

---

## netseq_refiner

NET-seq-guided position refinement for A-tract ambiguous reads.

::: rectify.core.netseq.netseq_refiner
    options:
      members:
        - NetseqLoader
        - refine_with_netseq
        - refine_batch

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

::: rectify.core.bam.processing_stats
    options:
      members:
        - ProcessingStats

---

## splice_aware_5prime

5' splice junction rescue.

::: rectify.core.splice.splice_aware_5prime
    options:
      members:
        - correct_5prime_for_splicing
        - rescue_3ss_truncation
        - build_intron_interval_tree

---

## polya_trimmer

Poly(A) tail detection and length measurement.

::: rectify.core.polya.polya_trimmer
    options:
      members:
        - trim_polya_from_read
        - trim_polya_from_bam_read
        - calculate_full_polya_length

---

## ag_mispriming

AG-mispriming detection for oligo-dT libraries.

::: rectify.core.polya.ag_mispriming
    options:
      members:
        - screen_ag_mispriming
        - calculate_ag_content
