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
        - cigar_query_length

---

## chromosome

Chromosome name normalization and size lookup.

::: rectify.utils.chromosome
    options:
      members:
        - normalize_chrom_name
        - get_chrom_size
        - CANONICAL_CHROMS
        - NCBI_TO_CHROM
        - CHROM_SIZES

---

## junction_bed

BED file generation for splice junction annotation.

::: rectify.utils.junction_bed
    options:
      members:
        - gff_to_junction_bed
        - gtf_to_junction_bed

---

## splice_motif

Canonical splice site detection.

::: rectify.utils.splice_motif
    options:
      members:
        - is_canonical_donor
        - is_canonical_acceptor
        - score_splice_motif

---

## stats

Statistics aggregation for QC reporting.

::: rectify.utils.stats
    options:
      members:
        - ProcessingStats
        - update_from_result
        - generate_summary_report

---

## provenance

Tool version tracking for reproducibility.

::: rectify.utils.provenance
    options:
      members:
        - get_tool_versions
        - record_provenance
