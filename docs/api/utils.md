# rectify.utils

Shared utility modules used throughout the RECTIFY pipeline.

---

## genome

Genome FASTA loading, sequence fetching, and A-tract detection.

::: rectify.utils.genome
    options:
      members:
        - load_genome
        - fetch_genomic_sequence
        - is_atract
        - count_contiguous_a_tract
        - complement
        - reverse_complement
        - standardize_chrom_name

---

## alignment

CIGAR string parsing, soft-clip extraction, and indel detection.

::: rectify.utils.alignment
    options:
      members:
        - parse_cigar
        - extract_soft_clips
        - extract_deletions
        - extract_insertions
        - get_query_length

---

## chromosome

Chromosome name normalization and size lookup.

::: rectify.utils.chromosome
    options:
      members:
        - normalize_chromosome
        - normalize_dataframe_chromosomes
        - detect_chromosome_format
        - build_chromosome_map
        - YEAST_NCBI_TO_CHR
        - YEAST_CHR_TO_NCBI

---

## junction_bed

BED file generation for splice junction annotation.

::: rectify.utils.junction_bed
    options:
      members:
        - generate_junction_bed
        - parse_gff_introns
        - parse_gtf_introns
        - write_junction_bed
        - get_minimap2_junc_args

---

## splice_motif

Canonical splice site detection.

::: rectify.utils.splice_motif
    options:
      members:
        - SpliceMotifScorer
        - NonCanonicalSSDetector
        - get_splice_site_dinucleotides
        - get_splice_site_sequences

---

## stats

Statistics aggregation for QC reporting.

::: rectify.utils.stats
    options:
      members:
        - calculate_summary_stats
        - format_summary_report
        - calculate_acount_distribution
        - calculate_shift_by_acount
        - assign_qc_flag

---

## provenance

Tool version tracking for reproducibility.

::: rectify.utils.provenance
    options:
      members:
        - ProvenanceTracker
        - init_provenance
